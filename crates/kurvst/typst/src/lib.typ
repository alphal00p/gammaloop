// Public Kurvst API.
//
// Keep path constructors, geometry helpers, and emitters here. The implementation
// lives in `impl.typ`.

#import "impl.typ" as _impl

/// Build a numeric point tuple.
/// -> array
#let point(
  /// X coordinate. -> int | float
  x,
  /// Y coordinate. -> int | float
  y,
) = _impl.point(x, y)

/// Default geometry options for derived path layers.
/// -> dictionary
#let layer-defaults = _impl.layer-defaults

/// A smooth sinusoidal path pattern.
/// -> dictionary
#let wave(
  /// Number of samples used to approximate one wave period. -> int
  samples-per-period: 16,
) = _impl.wave(samples-per-period: samples-per-period)

/// A straight-segment triangular path pattern.
/// -> dictionary
#let zigzag() = _impl.zigzag()

/// A smooth coil path pattern.
/// -> dictionary
#let coil(
  /// Number of samples used to approximate one coil period. -> int
  samples-per-period: 16,
  /// Horizontal scale of the coil before it is mapped onto a path. -> int | float
  longitudinal-scale: 1.25,
) = {
  _impl.coil(samples-per-period: samples-per-period, longitudinal-scale: longitudinal-scale)
}

/// Return `from` moved toward `toward` by `distance`.
/// -> array
#let outset-point(
  /// Point to move. -> array
  from,
  /// Target point that defines the direction. -> array
  toward,
  /// Distance to move from `from` toward `toward`. -> int | float
  distance: 0,
) = _impl.outset-point(from, toward, distance: distance)

/// Build a `move` path element.
/// -> dictionary
#let move-to(
  /// New current point and subpath start. -> array
  start,
) = _impl.move-to(start)

/// Build a `line` path element.
/// -> dictionary
#let line-to(
  /// Line endpoint. -> array
  end,
) = _impl.line-to(end)

/// Build a `quad` path element.
/// -> dictionary
#let quad-to(
  /// Quadratic control point. -> array
  control,
  /// Quadratic endpoint. -> array
  end,
) = _impl.quad-to(control, end)

/// Build a `cubic` path element.
/// -> dictionary
#let cubic-to(
  /// Cubic control point near the start point. -> array
  control-start,
  /// Cubic control point near the endpoint. -> array
  control-end,
  /// Cubic endpoint. -> array
  end,
) = _impl.cubic-to(control-start, control-end, end)

/// Build a `close` path element.
/// -> dictionary
#let close(
  /// Native Typst curve close mode. -> string
  mode: "straight",
) = _impl.close(mode: mode)

/// Build a path dictionary from an existing element array.
/// -> dictionary
#let from-elements(
  /// Array of Kurvst path elements. -> array
  elements,
) = _impl.from-elements(elements)

/// Build a path dictionary from path fragments or elements.
/// -> dictionary
#let path(
  /// Path fragments, path elements, or element arrays to concatenate. -> any
  ..parts,
) = _impl.path(..parts)

/// Return a path with additional fragments or elements appended.
/// -> dictionary
#let append(
  /// Base Kurvst path dictionary. -> dictionary
  path,
  /// Path fragments, path elements, or element arrays to append. -> any
  ..parts,
) = _impl.append(path, ..parts)

/// Build a straight-line path fragment.
/// -> dictionary
#let line(
  /// Start point. -> array
  start,
  /// Endpoint. -> array
  end,
) = _impl.line(start, end)

/// Build a quadratic path fragment.
/// -> dictionary
#let quad(
  /// Start point. -> array
  start,
  /// Quadratic control point. -> array
  control,
  /// Endpoint. -> array
  end,
) = _impl.quad(start, control, end)

/// Build a cubic path fragment.
/// -> dictionary
#let cubic(
  /// Start point. -> array
  start,
  /// Cubic control point near `start`. -> array
  control-start,
  /// Cubic control point near `end`. -> array
  control-end,
  /// Endpoint. -> array
  end,
) = _impl.cubic(start, control-start, control-end, end)

/// Build a path fragment from a cubic segment dictionary.
/// -> dictionary
#let from-cubic(
  /// Segment with `start`, `control-start`, `control-end`, and `end`. -> dictionary
  segment,
) = _impl.from-cubic(segment)

/// Build a cubic segment dictionary for a straight line.
/// -> dictionary
#let line-segment(
  /// Start point. -> array
  start,
  /// Endpoint. -> array
  end,
) = _impl.line-segment(start, end)

/// Return the command elements that make up a Kurvst path.
/// -> array
#let elements(
  /// Kurvst path dictionary to inspect. -> dictionary
  path,
) = _impl.elements(path)

/// Return the points visited by a Kurvst path's command stream.
/// -> array
#let points(
  /// Kurvst path dictionary to inspect. -> dictionary
  path,
) = _impl.points(path)

/// Evaluate a cubic segment at parameter `t`.
/// -> array
#let cubic-point(
  /// Segment with `start`, `control-start`, `control-end`, and `end`. -> dictionary
  segment,
  /// Segment parameter in the range `[0, 1]`. -> int | float
  t,
) = _impl.cubic-point(segment, t)

/// Evaluate the tangent of a cubic segment at parameter `t`.
/// -> array
#let cubic-tangent(
  /// Segment with `start`, `control-start`, `control-end`, and `end`. -> dictionary
  segment,
  /// Segment parameter in the range `[0, 1]`. -> int | float
  t,
) = _impl.cubic-tangent(segment, t)

/// Return drawable cubic segments for any Kurvst path dictionary.
/// -> array
#let segments(
  /// Kurvst path dictionary to convert. -> dictionary
  path,
) = _impl.segments(path)

/// Compute the arc length of a path dictionary.
/// -> int | float
#let length(
  /// Kurvst path dictionary to measure. -> dictionary
  path,
  /// Arc-length approximation accuracy passed to the Rust geometry engine. -> float
  accuracy: 0.001,
) = _impl.length(path, accuracy: accuracy)

/// Resolve a fixed and relative visible path length.
/// -> none | int | float
#let resolve-length(
  /// Full base path arc length. -> int | float
  base-length,
  /// Fixed target arc length. -> none | int | float
  length: none,
  /// Relative target length as a fraction of `base-length`. -> none | int | float
  ratio: none,
  /// Resolution strategy for fixed and relative targets. -> string | function
  method: "min",
) = {
  _impl.resolve-length(base-length, length: length, ratio: ratio, method: method)
}

/// Compute the symmetric trim needed to center a shorter path layer.
/// -> int | float
#let center-outset(
  /// Full base path arc length. -> int | float
  base-length,
  /// Fixed target visible length. -> none | int | float
  length: none,
  /// Relative target visible length as a fraction of `base-length`. -> none | int | float
  ratio: none,
  /// Resolution strategy for fixed and relative targets. -> string | function
  resolve-length: "min",
  /// Already-applied trim at the start of the path. -> int | float
  start-outset: 0,
  /// Already-applied trim at the end of the path. -> int | float
  end-outset: 0,
) = _impl.center-outset(
  base-length,
  length: length,
  ratio: ratio,
  resolve-length: resolve-length,
  start-outset: start-outset,
  end-outset: end-outset,
)

/// Trim a path by arc length from each end.
/// -> dictionary
#let trim(
  /// Kurvst path dictionary to trim. -> dictionary
  path,
  /// Arc length removed from the start. -> int | float
  start-outset: 0,
  /// Arc length removed from the end. -> int | float
  end-outset: 0,
  /// Arc-length approximation accuracy passed to the Rust geometry engine. -> float
  accuracy: 0.001,
) = {
  _impl.trim(path, start-outset: start-outset, end-outset: end-outset, accuracy: accuracy)
}

/// Construct a cubic Hobby path through three points.
/// -> dictionary
#let hobby-through(
  /// Start point. -> array
  start,
  /// Intermediate point that the curve passes through. -> array
  through,
  /// Endpoint. -> array
  end,
  /// Hobby curl/tension parameter. -> float
  omega: 1.0,
  /// Geometry approximation accuracy passed to the Rust geometry engine. -> float
  accuracy: 0.001,
) = {
  _impl.hobby-through(start, through, end, omega: omega, accuracy: accuracy)
}

/// Construct a Hobby spline through an arbitrary point sequence.
/// -> dictionary
#let hobby-spline(
  /// Two or more points for the open spline to pass through. -> array
  points,
  /// Hobby curl/tension parameter. -> float
  omega: 1.0,
  /// Geometry approximation accuracy passed to the Rust geometry engine. -> float
  accuracy: 0.001,
) = {
  _impl.hobby-spline(points, omega: omega, accuracy: accuracy)
}

/// Apply a repeated path pattern to a base path.
/// -> dictionary
#let pattern(
  /// Base Kurvst path dictionary. -> dictionary
  path,
  /// Pattern dictionary or built-in pattern name. -> string | dictionary
  pattern: "wave",
  /// Normal amplitude of the pattern in path units. -> int | float
  amplitude: 0.1,
  /// Arc length of one pattern period. -> int | float
  wavelength: 1.0,
  /// Initial phase offset in pattern periods. -> int | float
  phase: 0,
  /// Samples per period for string-resolved smooth patterns. -> int
  samples-per-period: 16,
  /// Longitudinal scale used when resolving the built-in coil pattern. -> int | float
  coil-longitudinal-scale: 1.25,
  /// Force the generated path to start on the base path. -> bool
  anchor-start: true,
  /// Force the generated path to end on the base path. -> bool
  anchor-end: true,
  /// Geometry approximation accuracy passed to the Rust geometry engine. -> float
  accuracy: 0.001,
) = _impl.pattern(
  path,
  pattern: pattern,
  amplitude: amplitude,
  wavelength: wavelength,
  phase: phase,
  samples-per-period: samples-per-period,
  coil-longitudinal-scale: coil-longitudinal-scale,
  anchor-start: anchor-start,
  anchor-end: anchor-end,
  accuracy: accuracy,
)

/// Generate a parallel path for a path.
/// -> dictionary
#let parallel(
  /// Base Kurvst path dictionary. -> dictionary
  path,
  /// Signed normal offset distance. -> int | float
  distance: 0,
  /// Arc length removed from the start before offsetting. -> int | float
  start-outset: 0,
  /// Arc length removed from the end before offsetting. -> int | float
  end-outset: 0,
  /// Geometry approximation accuracy passed to the Rust geometry engine. -> float
  accuracy: 0.001,
  /// Let Kurbo simplify/optimize the fitted path. -> bool
  optimize: true,
) = {
  _impl.parallel(
    path,
    distance: distance,
    start-outset: start-outset,
    end-outset: end-outset,
    accuracy: accuracy,
    optimize: optimize,
  )
}

/// Build a derived visible path layer.
/// -> dictionary
#let layer(
  /// Base Kurvst path dictionary. -> dictionary
  path,
  /// Signed normal offset distance. -> int | float
  offset: 0,
  /// Fixed target visible length. -> none | int | float
  length: none,
  /// Relative target visible length as a fraction of the base length. -> none | int | float
  ratio: none,
  /// Resolution strategy for fixed and relative targets. -> string | function
  resolve-length: "min",
  /// Arc length removed from the start. -> int | float
  start-outset: 0,
  /// Arc length removed from the end. -> int | float
  end-outset: 0,
  /// Optional point used to choose the sign of `offset`. -> none | array
  side-point: none,
  /// Geometry approximation accuracy passed to the Rust geometry engine. -> float
  accuracy: 0.001,
  /// Let Kurbo simplify/optimize fitted parallel paths. -> bool
  optimize: true,
) = _impl.layer(
  path,
  offset: offset,
  length: length,
  ratio: ratio,
  resolve-length: resolve-length,
  start-outset: start-outset,
  end-outset: end-outset,
  side-point: side-point,
  accuracy: accuracy,
  optimize: optimize,
)

/// Emit a Kurvst path as native Typst `curve` content.
/// -> content
#let to-native(
  /// Kurvst path dictionary to emit. -> dictionary
  path,
  /// Coordinate multiplier for emitted Typst curve points. -> int | float | length | ratio
  unit: 1,
  /// Native `curve` style arguments. -> any
  ..style,
) = _impl.to-native(path, unit: unit, ..style)

/// Emit a Kurvst path as CeTZ path data.
/// -> array
#let to-cetz-data(
  /// Kurvst path dictionary to convert. -> dictionary
  path,
  /// Coordinate multiplier for emitted CeTZ points. -> int | float | length | ratio
  unit: 1,
) = _impl.to-cetz-data(path, unit: unit)

/// Draw a path dictionary through CeTZ.
/// -> content
#let to-cetz(
  /// Kurvst path dictionary to draw. -> dictionary
  path,
  /// Coordinate multiplier for emitted CeTZ points. -> int | float | length | ratio
  unit: 1,
  /// CeTZ draw style arguments forwarded to `merge-path`. -> any
  ..style,
) = _impl.to-cetz(path, unit: unit, ..style)

/// Split a path through a point sequence into per-span paths.
/// -> dictionary
#let split-through(
  /// Two or more points for the curve to pass through. -> array
  points,
  /// Hobby curl/tension parameter. -> float
  omega: 1.0,
  /// Arc length removed from the first span start. -> int | float
  start-outset: 0,
  /// Arc length removed from the last span end. -> int | float
  end-outset: 0,
  /// Geometry approximation accuracy passed to the Rust geometry engine. -> float
  accuracy: 0.001,
) = {
  _impl.split-through(
    points,
    omega: omega,
    start-outset: start-outset,
    end-outset: end-outset,
    accuracy: accuracy,
  )
}
