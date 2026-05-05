# kurvst

`kurvst` is a Typst package backed by a small Kurbo WebAssembly plugin. It
provides Bezier splitting, arc-length trimming, Hobby-style edge curves,
parallel paths, and one-dimensional path patterns that can be emitted as CeTZ
drawing commands.

```typ
#import "src/lib.typ" as kurvst

#let segment = (
  start: (x: 0, y: 0),
  ctrl-a: (x: 1, y: 0.5),
  ctrl-b: (x: 2, y: -0.5),
  end: (x: 3, y: 0),
)

#let path = kurvst.pattern-segment(
  segment,
  pattern: kurvst.coil(longitudinal-scale: 1.6),
  amplitude: 0.15,
  wavelength: 0.7,
)

#let parallel = kurvst.parallel-segment(segment, distance: 0.18)
```

See `docs/manual.typ` for the full API manual.
