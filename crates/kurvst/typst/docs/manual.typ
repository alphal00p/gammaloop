#import "@preview/tidy:0.4.3"
#import "../src/lib.typ" as kurvst
#import "@preview/cetz:0.3.4" as cetz

#set document(title: "kurvst Typst API")
#set page(margin: 22mm)
#set text(size: 10pt)

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

Cubic segments use `start`, `ctrl_a`, `ctrl_b`, and `end`:

```typ
#let segment = (
  start: (x: 0, y: 0),
  ctrl_a: (x: 1, y: 0.5),
  ctrl_b: (x: 2, y: -0.5),
  end: (x: 3, y: 0),
)
```

== Splitting And Trimming

`split-cubic` splits one cubic at a parameter `t`. `trim-cubic` and
`trim-segment` trim by curve distance from the start and/or end. These helpers
return cubic segment dictionaries, so the output can be passed back into
`kurvst` or drawn with CeTZ.

```typ
#let split = kurvst.split-cubic(
  segment.start,
  segment.end,
  segment.ctrl_a,
  segment.ctrl_b,
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

`endpoint_ramp: true` tapers amplitude at anchored endpoints. This is useful
for coils, whose longitudinal offset would otherwise put the first and last
visible points inside the turn near nodes.

```typ
#let path = kurvst.pattern-segment(
  segment,
  pattern: hook,
  amplitude: 0.15,
  wavelength: 0.7,
)
#cetz.canvas({
  kurvst.cetz-pattern(path, stroke: black + 0.5pt)
})
```

== Generated Reference

#let tidy-style = dictionary(tidy.styles.default)
#let _ = tidy-style.insert("show-example", tidy-style.show-example.with(scale-preview: 100%))

#let docs = tidy.parse-module(
  read("../src/lib.typ"),
  name: "kurvst",
  scope: (cetz: cetz, kurvst: kurvst),
)
#tidy.show-module(docs, style: tidy-style)
