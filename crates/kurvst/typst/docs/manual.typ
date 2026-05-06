#import "@preview/tidy:0.4.3"
#import "../src/lib.typ" as kurvst
#import "@preview/cetz:0.3.4" as cetz

#set document(title: "kurvst Typst API")
#set page(margin: 22mm)
#set text(size: 10pt)

#let demo-scale = 36pt
#let demo-segment = (
  start: (x: 0.0, y: 1.0),
  ctrl-a: (x: 0.85, y: 0.1),
  ctrl-b: (x: 2.15, y: 1.9),
  end: (x: 3.0, y: 1.0),
)
#let demo-start = (x: 0.0, y: 1.2)
#let demo-through = (x: 1.5, y: 0.2)
#let demo-end = (x: 3.0, y: 1.2)
#let demo-nodes = (
  (pos: demo-start),
  (pos: demo-end),
)
#let demo-edge = (
  source: (node: 0),
  sink: (node: 1),
  pos: demo-through,
)

#let native-point(point, dx: 0, dy: 0) = (
  (point.x + dx) * demo-scale,
  (point.y + dy) * demo-scale,
)

#let native-cubic-items(segment, dx: 0, dy: 0) = (
  curve.move(native-point(segment.start, dx: dx, dy: dy)),
  curve.cubic(
    native-point(segment.ctrl-a, dx: dx, dy: dy),
    native-point(segment.ctrl-b, dx: dx, dy: dy),
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

All points are dictionaries with numeric `x` and `y` fields:

```typ
#let p = (x: 1.2, y: -0.4)
```

Cubic segments use `start`, `ctrl-a`, `ctrl-b`, and `end`:

```typ
#let segment = (
  start: (x: 0, y: 0),
  ctrl-a: (x: 1, y: 0.5),
  ctrl-b: (x: 2, y: -0.5),
  end: (x: 3, y: 0),
)
```

Use `cubic-path(..segment)` once to turn that cubic dictionary into a path
dictionary. Returned paths carry Kurvst's explicit wire path in `path`, plus
Typst-friendly `points`, cubic `curves`, line `segments`, and `length` fields
for drawing or inspection. Rust converts the wire path to Kurbo's `BezPath`
internally.

```typ
#let path = kurvst.cubic-path(..segment)
#let trimmed = kurvst.trim-path(path, start-outset: 0.2, end-outset: 0.1)
#let shifted = kurvst.parallel-path(trimmed, distance: 0.15)
```

A typical chain keeps the returned dictionaries intact until the final drawing
step:

```typ
#let base = kurvst.hobby-through(demo-start, demo-through, demo-end)
#let trimmed = kurvst.trim-path(base, start-outset: 0.2, end-outset: 0.1)
#let shifted = kurvst.parallel-path(trimmed, distance: 0.15)

#cetz.canvas({
  kurvst.cetz-path(base, stroke: rgb("#bbbbbb") + 0.4pt)
  kurvst.cetz-path(trimmed, stroke: rgb("#d72638") + 0.6pt)
  kurvst.cetz-path(shifted, stroke: rgb("#1b7f4c") + 0.6pt)
})
```

== Paths And Trimming

`cubic-path` constructs a path from one cubic Bezier. `trim-path` trims by curve
distance from the start and/or end and returns another path dictionary.

```typ
#let path = kurvst.cubic-path(..segment)
#let trimmed = kurvst.trim-path(path, start-outset: 0.15, end-outset: 0.1)
```

`hobby-through(start, through, end)` builds a smooth open curve through three
points; its `curves` field contains the two cubic halves split at `through`.
`hobby-spline(points)` does the same for any open point sequence with at least
two points. Both return path dictionaries and can be passed directly to
path-level helpers.

```typ
#let spline = kurvst.hobby-spline((
  (x: 0.0, y: 0.0),
  (x: 0.9, y: 0.8),
  (x: 1.8, y: -0.3),
  (x: 2.8, y: 0.4),
))
#let shifted = kurvst.parallel-path(spline, distance: 0.14)
#let wiggle = kurvst.pattern-path(
  spline,
  pattern: kurvst.wave(samples-per-period: 24),
  amplitude: 0.1,
  wavelength: 0.55,
)
#cetz.canvas({
  kurvst.cetz-path(spline, stroke: rgb("#bbbbbb") + 0.45pt)
  kurvst.cetz-path(shifted, stroke: rgb("#1b7f4c") + 0.6pt)
  kurvst.cetz-pattern(wiggle, stroke: rgb("#355c9a") + 0.55pt)
})
```

`edge-halves` combines the three-point Hobby helper with endpoint trimming for
graph-edge style geometry.

== Path Patterns

`pattern-path` generates a one-dimensional pattern along any path dictionary.
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
#let base = kurvst.cubic-path(..segment)
#let path = kurvst.pattern-path(
  base,
  pattern: hook,
  amplitude: 0.15,
  wavelength: 0.7,
)
#native-scene({
  native-cubics(base.curves, stroke: rgb("#c8c8c8") + 0.45pt)
  native-polyline(path.points, stroke: rgb("#1b7f4c") + 0.8pt)
})
```

== Parallel Paths

`parallel-path` uses Kurbo's offset curve fitter to produce a path at a fixed
normal distance from the source path. Positive distances follow the left normal
of the path direction; negative distances follow the right normal. The output
has the same drawable shape as pattern paths: `points`, fitted cubic `curves`,
straight `segments`, and `length`.
`start-outset` and `end-outset` trim the fitted path by arc length after the
offset is computed.

```typ
#let base = kurvst.cubic-path(..segment)
#let left = kurvst.parallel-path(base, distance: 0.18, start-outset: 0.35, end-outset: 0.35)
#let right = kurvst.parallel-path(base, distance: -0.18)
#native-scene({
  native-cubics(base.curves, stroke: rgb("#c8c8c8") + 0.45pt)
  native-cubics(left.curves, stroke: rgb("#1b7f4c") + 0.8pt)
  native-cubics(right.curves, stroke: rgb("#355c9a") + 0.8pt)
})
```

== Path Layers

`layer-path` is a convenience wrapper for drawing derived visible paths. It
combines side-aware offsetting, endpoint outsets, and centered shortening into
one path-in/path-out operation. `length` is a fixed visible arc length, `ratio`
is a fraction of the input path length, and `resolve-length` decides how to
combine them. The default `"min"` keeps whichever limit is shorter.

```typ
#let base = kurvst.hobby-spline((
  (x: 0.0, y: 0.0),
  (x: 0.9, y: 0.8),
  (x: 1.8, y: -0.3),
  (x: 2.8, y: 0.4),
))
#let layer = kurvst.layer-path(
  base,
  offset: 0.16,
  length: 1.6,
  ratio: 0.5,
  resolve-length: "min",
)
#cetz.canvas({
  kurvst.cetz-path(base, stroke: rgb("#bbbbbb") + 0.45pt)
  kurvst.cetz-path(layer, stroke: rgb("#1b7f4c") + 0.8pt)
})
```

Use `side-point` to choose the sign of the offset from a point on the desired
side of the path. Drawing packages can use this to place derived layers on the
same side as an edge label without knowing anything about the physics or graph
style that requested the layer.

== Native Drawing Primitives

Kurvst returns dictionaries and arrays. The core geometry can be drawn without
CeTZ by mapping points to Typst lengths and feeding them to native `curve`
components:

#tidy.show-example.show-example(
  ```typ
  #let point(p) = (p.x * 36pt, p.y * 36pt)
  #let draw-path(path) = curve(
    stroke: black + 0.7pt,
    ..path.curves.map(segment => (
      curve.move(point(segment.start)),
      curve.cubic(point(segment.ctrl-a), point(segment.ctrl-b), point(segment.end)),
    )).flatten(),
  )

  #draw-path(kurvst.cubic-path(..segment))
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
    let base = kurvst.cubic-path(..segment)
    kurvst.cetz-path(base, stroke: rgb("#d72638") + 0.6pt)

    let path = kurvst.pattern-path(
      base,
      pattern: kurvst.coil(longitudinal-scale: 1.6),
      amplitude: 0.12,
      wavelength: 0.55,
    )
    kurvst.cetz-pattern(path, stroke: rgb("#355c9a") + 0.55pt)
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
    demo-nodes: demo-nodes,
    demo-edge: demo-edge,
    native-scene: native-scene,
    native-cubic: native-cubic,
    native-cubics: native-cubics,
    native-polyline: native-polyline,
    native-dot: native-dot,
  ),
)
#tidy.show-module(docs, style: tidy-style)
