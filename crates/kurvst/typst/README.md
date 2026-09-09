# kurvst

`kurvst` is a Typst package backed by a small Kurbo WebAssembly plugin. It
provides Bezier splitting, arc-length trimming, Hobby-style edge curves,
parallel paths, and one-dimensional path patterns that can be emitted as CeTZ
drawing commands.

```typ
#import "src/lib.typ" as kurvst

#let segment = (
  start: (0, 0),
  control-start: (1, 0.5),
  control-end: (2, -0.5),
  end: (3, 0),
)

#let base = kurvst.from-cubic(segment)
#let path = kurvst.pattern(
  base,
  pattern: kurvst.coil(longitudinal-scale: 1.6),
  amplitude: 0.15,
  wavelength: 0.7,
)

#let parallel = kurvst.parallel(base, distance: 0.18)
```

Kurvst requires Typst 0.15.0 or newer. See `docs/manual.typ` for the full API manual and
`examples/` for compile-checked drawings.
