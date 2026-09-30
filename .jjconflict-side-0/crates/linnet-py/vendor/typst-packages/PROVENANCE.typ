= Bundled Typst packages

These packages are embedded in the `linnet` wheel so rendering
does not depend on a network connection or an installed Typst package cache.

- CeTZ 0.5.1 (`preview/cetz/0.5.1`), LGPL-3.0-or-later
  - Source: `https://packages.typst.org/preview/cetz-0.5.1.tar.gz`
  - Nix recursive hash: `sha256-qqYcK/qrjU0T1X2+1g0CZHbzB3L5QVefmiL+XSyuCrY=`
- oxifmt 1.0.0 (`preview/oxifmt/1.0.0`), MIT OR Apache-2.0
  - Source: `https://packages.typst.org/preview/oxifmt-1.0.0.tar.gz`
  - Nix recursive hash: `sha256-RtGKdyiX2kJbUjChPohSGNeYOKVlI2VM0k1uFaEqDC8=`

The corresponding license files and package manifests remain in each package
directory.

== Local performance patch

CeTZ's `merge-path` skips per-element bounding-box calculations when debug
rendering is disabled. The merged path still computes its final bounds; debug
behavior and path geometry are unchanged. The hashes above identify the
unmodified upstream archives.

CeTZ content separates preparation from placement. Repeated content can reuse
resolved styles, measurements and local geometry while preserving per-position
coordinate resolution, transformed bounds, frames and named anchors. Linnest
uses this preparation for consecutive identical inspection targets.
