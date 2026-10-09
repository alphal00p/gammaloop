= Bundled Typst packages

The renderer embeds these packages for offline use; preparation never downloads
packages. External stores supply additional packages, but bundled packages win.

- CeTZ 0.5.2, LGPL-3.0-or-later
  - Pristine archive: `archives/cetz-0.5.2.tar.gz`
  - Source: `https://packages.typst.org/preview/cetz-0.5.2.tar.gz`
  - Archive SHA-256: `77cf8490114ae04c6e665a11efa691d284a0cadb9719771b5708c1197292f23f`
  - Upstream release: `v0.5.2`, commit `40714e2f939be64206eaf86d7f5e45e054b05333`
  - License: `archives/LICENSE.cetz`, identical to the archive's `LICENSE`.
  - Preparation stages the original archive files under `preview/cetz/0.5.2`.
    No local patches or source fork are maintained. CeTZ's external painter
    remains LGPL; its code is not ported into MIT components.
- oxifmt 1.0.0 (`preview/oxifmt/1.0.0`), MIT OR Apache-2.0
  - Source: `https://packages.typst.org/preview/oxifmt-1.0.0.tar.gz`
  - Nix recursive hash: `sha256-RtGKdyiX2kJbUjChPohSGNeYOKVlI2VM0k1uFaEqDC8=`
- MiTeX 0.2.6 (`preview/mitex/0.2.6`), Apache-2.0
  - Source: `https://packages.typst.org/preview/mitex-0.2.6.tar.gz`

oxifmt and MiTeX retain their licenses and manifests in their package directories.
The archive, standalone CeTZ license and this notice are distribution assets,
not additional Typst packages.
