= Bundled Typst packages

The renderer embeds these packages for offline use; preparation never downloads
packages. They precede the local Typst package path and cache, which supply
every other package.

- oxifmt 1.0.0 (`preview/oxifmt/1.0.0`), MIT OR Apache-2.0
  - Source: `https://packages.typst.org/preview/oxifmt-1.0.0.tar.gz`
  - Nix recursive hash: `sha256-RtGKdyiX2kJbUjChPohSGNeYOKVlI2VM0k1uFaEqDC8=`
- MiTeX 0.2.6 (`preview/mitex/0.2.6`), Apache-2.0
  - Source: `https://packages.typst.org/preview/mitex-0.2.6.tar.gz`

oxifmt and MiTeX retain their licenses and manifests in their package directories.

CeTZ (LGPL-3.0-or-later) is not bundled. Linnest draws with stock CeTZ 0.5.2
from Typst Universe, which the Typst route reads from the local package
directories; native drawings need no Typst package.
