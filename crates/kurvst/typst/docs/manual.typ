#import "@preview/tidy:0.4.3"
#import "../src/lib.typ" as kurvst
#import "@preview/cetz:0.5.1" as cetz

#set document(title: "kurvst Typst API")
#set page(margin: 22mm)
#set text(size: 10pt)

#let demo-scale = 36pt
#let demo-segment = (
  start: (0.0, 1.0),
  control-start: (0.85, 0.1),
  control-end: (2.15, 1.9),
  end: (3.0, 1.0),
)
#let demo-start = (0.0, 1.2)
#let demo-through = (1.5, 0.2)
#let demo-end = (3.0, 1.2)

#let native-point(point, dx: 0, dy: 0) = (
  (point.at(0) + dx) * demo-scale,
  (point.at(1) + dy) * demo-scale,
)

#let native-cubic-items(segment, dx: 0, dy: 0) = (
  curve.move(native-point(segment.start, dx: dx, dy: dy)),
  curve.cubic(
    native-point(segment.control-start, dx: dx, dy: dy),
    native-point(segment.control-end, dx: dx, dy: dy),
    native-point(segment.end, dx: dx, dy: dy),
  ),
)

#let native-scene(body, width: 128pt, height: 82pt) = block(
  width: width,
  height: height,
  inset: 0pt,
  body,
)

#let native-cubic(segment, stroke: black + 0.6pt, dx: 0, dy: 0) = place(
  dx: 0pt,
  dy: 0pt,
  curve(stroke: stroke, ..native-cubic-items(segment, dx: dx, dy: dy)),
)

#let native-cubics(segments, stroke: black + 0.6pt, dx: 0, dy: 0) = place(
  dx: 0pt,
  dy: 0pt,
  curve(
    stroke: stroke,
    ..segments.map(segment => native-cubic-items(segment, dx: dx, dy: dy)).flatten(),
  ),
)

#let native-polyline(points, stroke: black + 0.6pt, dx: 0, dy: 0) = {
  if points.len() < 2 {
    ()
  } else {
    place(
      dx: 0pt,
      dy: 0pt,
      curve(
        stroke: stroke,
        curve.move(native-point(points.at(0), dx: dx, dy: dy)),
        ..points.slice(1).map(point => curve.line(native-point(point, dx: dx, dy: dy))),
      ),
    )
  }
}

#let native-dot(point, fill: black, radius: 2.2pt, dx: 0, dy: 0) = {
  let pos = native-point(point, dx: dx, dy: dy)
  place(
    dx: pos.at(0) - radius,
    dy: pos.at(1) - radius,
    circle(radius: radius, fill: fill),
  )
}

= kurvst Typst API

`kurvst` is a Typst package backed by `kurvst.wasm`, a small Kurbo-based
plugin. Its public API is path-based: construct a path dictionary, pass that
path into transforms, and draw the returned path.

It exposes Bezier utilities that are awkward or unavailable in native Typst:

- constructing cubic and Hobby paths,
- trimming paths by arc length,
- generating sampled path patterns,
- fitting Kurbo offset/parallel paths,
- converting returned path geometry to CeTZ drawing commands.

Path geometry points are two-item tuples:

```typ
#let p = (1.2, -0.4)
```

Cubic segments use `start`, `control-start`, `control-end`, and `end`:

```typ
#let segment = (
  start: (0, 0),
  control-start: (1, 0.5),
  control-end: (2, -0.5),
  end: (3, 0),
)
```

Use `from-cubic(segment)` once to turn that cubic dictionary into a path
dictionary. Returned paths contain a curve-command wire path in `path`. Rust
converts that wire path to Kurbo's `BezPath` internally, and every path-level
Kurvst function accepts either a returned dictionary or the raw wire path.

```typ
#let path = kurvst.from-cubic(segment)
#let trimmed = kurvst.trim(path, start-outset: 0.2, end-outset: 0.1)
#let shifted = kurvst.parallel(trimmed, distance: 0.15)
```

== Path Wire Format

A Kurvst wire path is a dictionary with one `elements` array. The elements mirror
Typst's native `curve` commands, but points are plain numeric tuples:

```typ
#let path = (
  elements: (
    (kind: "move", start: (0, 0)),
    (kind: "line", end: (0.6, 0.3)),
    (kind: "quad", control: (1.0, 1.0), end: (1.4, 0.3)),
    (
      kind: "cubic",
      control-start: (1.8, -0.4),
      control-end: (2.4, 1.0),
      end: (3.0, 0),
    ),
    (kind: "close", mode: "straight"),
  ),
)
```

The supported element kinds are:

- `move`: starts a subpath at `start`.
- `line`: adds a straight segment ending at `end`.
- `quad`: adds a quadratic segment using `control` and `end`.
- `cubic`: adds a cubic segment using `control-start`, `control-end`, and `end`.
- `close`: closes the current subpath. `mode` is optional and defaults to
  `"straight"` when emitted by Kurvst.

Prefer the path fragment helpers for new code. `line`, `quad`, and `cubic`
include their start points, so fragments can be built independently:

```typ
#let custom = kurvst.path(
  kurvst.line((0, 0), (0.6, 0.3)),
  kurvst.quad((0.6, 0.3), (1.0, 1.0), (1.4, 0.3)),
  kurvst.cubic((1.4, 0.3), (1.8, -0.4), (2.4, 1.0), (3.0, 0)),
)
```

`path` and `append` flatten fragments. If an appended fragment starts where
the current path ends, its leading `move` is skipped; if it starts elsewhere, the
`move` begins a new subpath. The edited path can go straight back into Kurvst
transforms:

```typ
#let base = kurvst.from-cubic(segment)
#let extended = kurvst.append(base, kurvst.line(segment.end, (3.4, 0.2)))

#let trimmed = kurvst.trim(extended, start-outset: 0.15)
#let shifted = kurvst.parallel(trimmed, distance: 0.1)

#kurvst.to-native(shifted, unit: 36pt, stroke: rgb("#1b7f4c") + 0.7pt)
```

Use `elements` when you need the raw command stream, `points` for the
visited endpoints, `segments` for cubic segment dictionaries, `to-native`
for native Typst drawing, and `to-cetz-data` or `to-cetz` for CeTZ.

A typical chain keeps the returned dictionaries intact until the final drawing
step:

```typ
#let base = kurvst.hobby-through(demo-start, demo-through, demo-end)
#let trimmed = kurvst.trim(base, start-outset: 0.2, end-outset: 0.1)
#let shifted = kurvst.parallel(trimmed, distance: 0.15)

#cetz.canvas({
  kurvst.to-cetz(base, stroke: rgb("#bbbbbb") + 0.4pt)
  kurvst.to-cetz(trimmed, stroke: rgb("#d72638") + 0.6pt)
  kurvst.to-cetz(shifted, stroke: rgb("#1b7f4c") + 0.6pt)
})
```

== Paths And Trimming

`cubic` constructs a path from one cubic Bezier. `trim` trims by curve
distance from the start and/or end and returns another path dictionary.

```typ
#let path = kurvst.from-cubic(segment)
#let trimmed = kurvst.trim(path, start-outset: 0.15, end-outset: 0.1)
```

`hobby-through(start, through, end)` builds a smooth open curve through three
points; `segments` returns the two cubic halves split at `through`.
`hobby-spline(points)` does the same for any open point sequence with at least
two points. Both return path dictionaries and can be passed directly to
path-level helpers.

```typ
#let spline = kurvst.hobby-spline((
  (0.0, 0.0),
  (0.9, 0.8),
  (1.8, -0.3),
  (2.8, 0.4),
))
#let shifted = kurvst.parallel(spline, distance: 0.14)
#let wiggle = kurvst.pattern(
  spline,
  pattern: kurvst.wave(samples-per-period: 24),
  amplitude: 0.1,
  wavelength: 0.55,
)
#cetz.canvas({
  kurvst.to-cetz(spline, stroke: rgb("#bbbbbb") + 0.45pt)
  kurvst.to-cetz(shifted, stroke: rgb("#1b7f4c") + 0.6pt)
  kurvst.to-cetz(wiggle, stroke: rgb("#355c9a") + 0.55pt)
})
```

`split-through` combines Hobby spline construction with endpoint trimming and
returns one path part between each consecutive point. Consumers such as Linnest
can wrap it for graph-edge geometry.

== Path Patterns

`pattern` generates a one-dimensional pattern along any path dictionary.
Built-ins are regular Typst dictionaries:

```typ
#let wave = kurvst.wave()
#let zigzag = kurvst.zigzag()
#let coil = kurvst.coil(longitudinal-scale: 1.6)
```

String names `"wave"`, `"zigzag"`, and `"coil"` are accepted for convenience,
but they are resolved in Typst before calling wasm.

Custom patterns use the same object shape:

```typ
#let hook = (
  kind: "points",
  name: "hook",
  interpolation: "linear",
  points: (
    (at: 0, x: 0, y: 0),
    (at: 0.25, x: 0.15, y: 1),
    (at: 0.5, x: 0, y: 0),
    (at: 0.75, x: -0.15, y: -1),
    (at: 1, x: 0, y: 0),
  ),
)
```

`at` is the phase position inside one wavelength, from `0` to `1`. The `x`
coordinate offsets along the path tangent and `y` offsets along the path normal;
both are scaled by `amplitude`. If `at` is omitted, points are spaced evenly.
Use `interpolation: "linear"` for corners and `"smooth"` for a spline through
sampled points.

`endpoint-ramp: true` tapers amplitude at anchored endpoints. This is useful
for coils, whose longitudinal offset would otherwise put the first and last
visible points inside the turn near nodes.

```typ
#let base = kurvst.from-cubic(segment)
#let path = kurvst.pattern(
  base,
  pattern: hook,
  amplitude: 0.15,
  wavelength: 0.7,
)
#native-scene({
  native-cubics(kurvst.segments(base), stroke: rgb("#c8c8c8") + 0.45pt)
  native-polyline(kurvst.points(path), stroke: rgb("#1b7f4c") + 0.8pt)
})
```

== Parallel Paths

`parallel` uses Kurbo's offset curve fitter to produce a path at a fixed
normal distance from the source path. Positive distances follow the left normal
of the path direction; negative distances follow the right normal.
`start-outset` and `end-outset` trim the fitted path by arc length after the
offset is computed.

```typ
#let base = kurvst.from-cubic(segment)
#let left = kurvst.parallel(base, distance: 0.18, start-outset: 0.35, end-outset: 0.35)
#let right = kurvst.parallel(base, distance: -0.18)
#native-scene({
  native-cubics(kurvst.segments(base), stroke: rgb("#c8c8c8") + 0.45pt)
  native-cubics(kurvst.segments(left), stroke: rgb("#1b7f4c") + 0.8pt)
  native-cubics(kurvst.segments(right), stroke: rgb("#355c9a") + 0.8pt)
})
```

== Path Layers

`layer` is a convenience wrapper for drawing derived visible paths. It
combines side-aware offsetting, endpoint outsets, and centered shortening into
one path-in/path-out operation. `length` is a fixed visible arc length, `ratio`
is a fraction of the input path length, and `resolve-length` decides how to
combine them. The default `"min"` keeps whichever limit is shorter.

```typ
#let base = kurvst.hobby-spline((
  (0.0, 0.0),
  (0.9, 0.8),
  (1.8, -0.3),
  (2.8, 0.4),
))
#let layer = kurvst.layer(
  base,
  offset: 0.16,
  length: 1.6,
  ratio: 0.5,
  resolve-length: "min",
)
#cetz.canvas({
  kurvst.to-cetz(base, stroke: rgb("#bbbbbb") + 0.45pt)
  kurvst.to-cetz(layer, stroke: rgb("#1b7f4c") + 0.8pt)
})
```

Use `side-point` to choose the sign of the offset from a point on the desired
side of the path. Drawing packages can use this to place derived layers on the
same side as an edge label without knowing anything about the physics or graph
style that requested the layer.

== Native Drawing Primitives

Kurvst returns dictionaries and arrays. The core geometry can be drawn without
CeTZ by emitting native `curve` content:

#tidy.show-example.show-example(
  ```typ
  #kurvst.to-native(kurvst.from-cubic(segment), unit: 36pt, stroke: black + 0.7pt)
  ```,
  scope: (kurvst: kurvst, segment: demo-segment),
)

The generated reference below keeps drawing code inline. It only uses the small
common helpers above for coordinate scaling, native `curve` construction, point
markers, and fixed-size preview blocks.

== CeTZ Interoperability

The CeTZ helpers are thin adapters over the same returned geometry. Use them
when the surrounding document already lives in a CeTZ canvas, or when you want
CeTZ path merging and styling:

#tidy.show-example.show-example(
  ```typ
  #cetz.canvas({
    let base = kurvst.from-cubic(segment)
    kurvst.to-cetz(base, stroke: rgb("#d72638") + 0.6pt)

    let path = kurvst.pattern(
      base,
      pattern: kurvst.coil(longitudinal-scale: 1.6),
      amplitude: 0.12,
      wavelength: 0.55,
    )
    kurvst.to-cetz(path, stroke: rgb("#355c9a") + 0.55pt)
  })
  ```,
  scope: (kurvst: kurvst, cetz: cetz, segment: demo-segment),
)

== Generated Reference

#let tidy-style = dictionary(tidy.styles.default)
#let _ = tidy-style.insert("show-example", tidy-style.show-example.with(scale-preview: 100%))

#let docs = tidy.parse-module(
  read("../src/lib.typ"),
  name: "kurvst",
  scope: (
    cetz: cetz,
    kurvst: kurvst,
    demo-segment: demo-segment,
    demo-start: demo-start,
    demo-through: demo-through,
    demo-end: demo-end,
    native-scene: native-scene,
    native-cubic: native-cubic,
    native-cubics: native-cubics,
    native-polyline: native-polyline,
    native-dot: native-dot,
  ),
)
#tidy.show-module(docs, style: tidy-style)
