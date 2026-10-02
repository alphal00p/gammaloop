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
rendering is disabled, and binds its processing helper directly instead of
recreating a closure. The merged path still computes its final bounds; debug
behavior and path geometry are unchanged. The hashes above identify the
unmodified upstream archives.

CeTZ content separates preparation from placement. Repeated content can reuse
resolved styles, measurements and local geometry while preserving per-position
coordinate resolution, transformed bounds, frames and named anchors. Linnest
uses this preparation for consecutive identical inspection targets. Numeric
centers with rectangular frames are placed in one native batch, retaining the
same transform arithmetic, line order and final coordinate context. Other
coordinate systems, resolver hooks, frames and nonfinite geometry use scalar
placement from the same prepared plan. Preparation receives only the style,
length and transform it uses; coordinate resolution and placement retain the
full canvas context. Prepared plans store geometry and drawable templates;
a shared placement function consumes them, avoiding per-plan closures and their
repeated hashing during batch processing.

Content computes its bounds from the four frame corners and evaluates secondary
anchors when requested. `line-strip` normalizes its closing and zero-length
segments while constructing the path, avoiding a second traversal.

Path normalization filters zero-length lines in one pass. Drawable tag helpers
avoid recursive dispatch for individual drawables, and the canvas inlines cubic
coordinate transforms while retaining their arithmetic order.
After bounds and anchors have been computed, the canvas omits paths whose stroke
and fill are both `none`. These paths cannot paint, so skipping their final bounds
queries and curve construction removes redundant SVG elements. Content drawables,
including transparent inspection links and their hit rectangles, remain unchanged.

CeTZ core batches each path's cubic extrema and AABB reduction in one plugin
call. The scalar extrema and AABB routines are unchanged. Finite paths reuse
their two bounds corners during aggregation; nonfinite paths retain the full
point stream.

Marker placement batches resolved reference contacts and carrier trimming in
CeTZ core, using the existing fixed-sample cubic measurements. Shape evaluation,
style resolution, and custom markers remain in Typst. Painting and annotation
collision bounds consume the same placed geometry. The numeric core preserves
contact and transform arithmetic, validates finite geometry and dimensions,
and returns floating-point coordinates. Typst's `libm` trigonometric dependency
is pinned in the marker patch.

Candidate marker footprints reuse shared canonical templates and resolve each
carrier's stations in one native batch. The marker painter owns contacts and
trimming, and the existing bounds owner supplies cubic extrema. Collision boxes
retain per-command shafts, whole marker paths, stroke radii, and content
position/dimension behavior. Shape callbacks that require their full carrier
context continue through the ordinary drawing owner.
Canonical carrier geometry and resulting footprint bounds can remain in CBOR
between native calls, avoiding repeated Typst path reconstruction. Encoded and
expanded inputs enter the same numeric painter after shaft-style resolution.

The source patches and reproducible build instructions are in
`preview/cetz/0.5.1/cetz-core/batched-bounds.patch`,
`preview/cetz/0.5.1/cetz-core/content-placement.patch`,
`preview/cetz/0.5.1/cetz-core/sampled-curves.patch`,
`preview/cetz/0.5.1/cetz-core/marker-geometry.patch`,
`preview/cetz/0.5.1/cetz-core/marker-footprints.patch`, and
`preview/cetz/0.5.1/cetz-core/REBUILD.typ`. They apply to upstream CeTZ `v0.5.1`,
commit `100ea3c52342c8bf866e8f154589ce8cd309265c`, under the same LGPL license. The core uses release optimization level 3 for
speed; the exact compiler, linker, dependency lockfile, and override are
recorded in the build instructions.
