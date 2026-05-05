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

#let native-pattern-object(pattern, stroke: rgb("#1b7f4c") + 0.75pt) = {
  let points = pattern.points.map(point => (
    x: point.at * 3.1 + point.x * 0.18,
    y: 1 + point.y * 0.32,
  ))
  native-scene({
    native-polyline(((x: 0.0, y: 1.0), (x: 3.15, y: 1.0)), stroke: rgb("#bbbbbb") + 0.35pt)
    native-polyline(points, stroke: stroke)
  })
}

#let native-split-demo(split) = native-scene({
  if split.keys().contains("curve") {
    native-cubic(split.curve, stroke: rgb("#c8c8c8") + 0.45pt)
  } else {
    native-cubics(split.segments, stroke: rgb("#c8c8c8") + 0.45pt)
  }
  native-cubic(split.segments.at(0), stroke: rgb("#d72638") + 0.95pt)
  native-cubic(split.segments.at(1), stroke: rgb("#355c9a") + 0.95pt)
  native-dot(split.segments.at(0).end, fill: black)
})

#let native-trim-demo(original, trimmed) = native-scene({
  native-cubic(original, stroke: rgb("#c8c8c8") + 0.45pt)
  native-cubic(trimmed, stroke: rgb("#d72638") + 1pt)
  native-dot(trimmed.start, fill: rgb("#d72638"))
  native-dot(trimmed.end, fill: rgb("#d72638"))
})

#let native-pattern-path(path, base: demo-segment) = native-scene({
  native-cubic(base, stroke: rgb("#c8c8c8") + 0.45pt)
  if path.curves.len() > 0 {
    native-cubics(path.curves, stroke: rgb("#1b7f4c") + 0.8pt)
  } else {
    native-polyline(path.points, stroke: rgb("#1b7f4c") + 0.8pt)
  }
})

#let native-outset-demo(from, toward, moved) = native-scene({
  native-polyline((from, toward), stroke: rgb("#c8c8c8") + 0.65pt)
  native-dot(from, fill: rgb("#355c9a"))
  native-dot(moved, fill: rgb("#d72638"))
  native-dot(toward, fill: black)
})

#let native-edge-halves-demo(halves) = native-scene({
  native-cubic(halves.source, stroke: rgb("#d72638") + 0.9pt)
  native-cubic(halves.sink, stroke: rgb("#355c9a") + 0.9pt)
  native-dot(demo-through, fill: black)
})

= kurvst Typst API

`kurvst` is a Typst package backed by `kurvst.wasm`, a small Kurbo-based
plugin. It exposes Bezier utilities that are awkward or unavailable in native
Typst:

- splitting cubic Beziers,
- trimming cubic Beziers by arc length,
- constructing smooth Hobby curves through a control point,
- generating sampled path patterns,
- converting returned geometry to CeTZ drawing commands.

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

== Splitting And Trimming

`split-cubic` splits one cubic at a parameter `t`. `trim-cubic` and
`trim-segment` trim by curve distance from the start and/or end. These helpers
return cubic segment dictionaries, so the output can be passed back into
`kurvst`, drawn with native Typst `curve`, or forwarded to CeTZ.

```typ
#let split = kurvst.split-cubic(
  segment.start,
  segment.end,
  segment.ctrl-a,
  segment.ctrl-b,
  t: 0.5,
)
#let trimmed = kurvst.trim-segment(split.segments.at(0), start-outset: 0.15)
```

`hobby-through(start, through, end)` builds a smooth open curve through all
three points and returns two cubic halves split at `through`. `edge-halves`
combines that with endpoint trimming for graph-edge style geometry.

== Path Patterns

`pattern-cubic` and `pattern-segment` generate a one-dimensional pattern along a
cubic path. Built-ins are regular Typst dictionaries:

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
#let path = kurvst.pattern-segment(
  segment,
  pattern: hook,
  amplitude: 0.15,
  wavelength: 0.7,
)
#native-pattern-path(path, base: segment)
```

== Native Drawing Primitives

Kurvst returns dictionaries and arrays. The core geometry can be drawn without
CeTZ by mapping points to Typst lengths and feeding them to native `curve`
components:

#tidy.show-example.show-example(
  ```typ
  #let point(p) = (p.x * 36pt, p.y * 36pt)
  #let draw-segment(segment) = curve(
    stroke: black + 0.7pt,
    curve.move(point(segment.start)),
    curve.cubic(point(segment.ctrl-a), point(segment.ctrl-b), point(segment.end)),
  )

  #draw-segment(segment)
  ```,
  scope: (segment: demo-segment),
)

The generated reference below uses the same idea for each core function: native
`curve` for cubics and polylines, `circle` for point markers, and `place` inside
a fixed-size `block` when several primitives need to be overlaid.

== CeTZ Interoperability

The CeTZ helpers are thin adapters over the same returned geometry. Use them
when the surrounding document already lives in a CeTZ canvas, or when you want
CeTZ path merging and styling:

#tidy.show-example.show-example(
  ```typ
  #cetz.canvas({
    kurvst.cetz-bezier(segment, stroke: rgb("#d72638") + 0.6pt)

    let path = kurvst.pattern-segment(
      segment,
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
    native-pattern-object: native-pattern-object,
    native-split-demo: native-split-demo,
    native-trim-demo: native-trim-demo,
    native-pattern-path: native-pattern-path,
    native-outset-demo: native-outset-demo,
    native-edge-halves-demo: native-edge-halves-demo,
  ),
)
#tidy.show-module(docs, style: tidy-style)
