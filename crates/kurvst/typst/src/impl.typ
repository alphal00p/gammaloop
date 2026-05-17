// Internal Kurvst implementation. Public users should import `lib.typ`.

#import "@preview/cetz:0.5.1" as cetz

#let _plugin = plugin("../kurvst.wasm")

#let _point-x(p) = p.at(0)
#let _point-y(p) = p.at(1)
#let _point-pair(p) = (_point-x(p), _point-y(p))
#let _point(p, unit: 1) = (_point-x(p) * unit, _point-y(p) * unit)

/// Build a numeric point tuple.
///
/// -> array
#let point(x, y) = (x, y)

/// Default geometry options for derived path layers.
///
/// These defaults are shared by @layer and drawing packages that build on
/// Kurvst. The dictionary shape intentionally mirrors the public function
/// arguments so callers can merge user styles into it and unpack the result.
/// -> dictionary
#let layer-defaults = (
  offset: 0,
  length: none,
  ratio: none,
  resolve-length: "min",
  start-outset: 0,
  end-outset: 0,
  side-point: none,
  accuracy: 0.001,
  optimize: true,
)

#let _sampled-pattern(samples-per-period, offset) = {
  let samples = calc.max(1, samples-per-period)
  range(0, samples + 1).map(index => {
    let at = index / samples
    let theta = 2 * calc.pi * at
    let point = offset(theta)
    (at: at, x: _point-x(point), y: _point-y(point))
  })
}

/// A smooth sinusoidal path pattern.
///
/// ```example
/// #let pattern = kurvst.wave(samples-per-period: 24)
/// #let points = pattern.points.map(point => (
///   point.at * 3.1 + point.x * 0.18,
///   1 + point.y * 0.32,
/// ))
/// #native-scene({
///   native-polyline(((0.0, 1.0), (3.15, 1.0)), stroke: rgb("#bbbbbb") + 0.35pt)
///   native-polyline(points, stroke: rgb("#1b7f4c") + 0.75pt)
/// })
/// ```
///
/// -> dictionary
#let wave(samples-per-period: 16) = (
  kind: "points",
  name: "wave",
  interpolation: "smooth",
  points: _sampled-pattern(samples-per-period, theta => (0, calc.sin(theta))),
)

/// A straight-segment triangular path pattern.
///
/// ```example
/// #let pattern = kurvst.zigzag()
/// #let points = pattern.points.map(point => (
///   point.at * 3.1 + point.x * 0.18,
///   1 + point.y * 0.32,
/// ))
/// #native-scene({
///   native-polyline(((0.0, 1.0), (3.15, 1.0)), stroke: rgb("#bbbbbb") + 0.35pt)
///   native-polyline(points, stroke: rgb("#1b7f4c") + 0.75pt)
/// })
/// ```
///
/// -> dictionary
#let zigzag() = (
  kind: "points",
  name: "zigzag",
  interpolation: "linear",
  points: (
    (at: 0, x: 0, y: 0),
    (at: 0.25, x: 0, y: 1),
    (at: 0.75, x: 0, y: -1),
    (at: 1, x: 0, y: 0),
  ),
)

/// A smooth coil path pattern.
///
/// ```example
/// #let pattern = kurvst.coil(samples-per-period: 24, longitudinal-scale: 1.6)
/// #let points = pattern.points.map(point => (
///   point.at * 3.1 + point.x * 0.18,
///   1 + point.y * 0.32,
/// ))
/// #native-scene({
///   native-polyline(((0.0, 1.0), (3.15, 1.0)), stroke: rgb("#bbbbbb") + 0.35pt)
///   native-polyline(points, stroke: rgb("#1b7f4c") + 0.75pt)
/// })
/// ```
///
/// -> dictionary
#let coil(samples-per-period: 16, longitudinal-scale: 1.25) = (
  kind: "points",
  name: "coil",
  interpolation: "smooth",
  endpoint-ramp: true,
  points: _sampled-pattern(
    samples-per-period,
    theta => (longitudinal-scale * calc.cos(theta), calc.sin(theta)),
  ),
)

#let _resolve-pattern(pattern, samples-per-period: 16, coil-longitudinal-scale: 1.25) = {
  if type(pattern) != str {
    pattern
  } else {
    let name = pattern.trim()
    if name == "wave" or name == "sine" or name == "sin" {
      wave(samples-per-period: samples-per-period)
    } else if name == "zigzag" or name == "zig-zag" or name == "triangle" {
      zigzag()
    } else if name == "coil" or name == "helix" or name == "spring" {
      coil(samples-per-period: samples-per-period, longitudinal-scale: coil-longitudinal-scale)
    } else {
      panic("Unsupported path pattern: " + pattern)
    }
  }
}

#let outset-point(from, toward, distance: 0) = {
  let dx = _point-x(toward) - _point-x(from)
  let dy = _point-y(toward) - _point-y(from)
  let length = calc.sqrt(dx * dx + dy * dy)
  if distance == 0 or length == 0 {
    _point-pair(from)
  } else {
    let applied = calc.min(distance, length * 0.45)
    (
      _point-x(from) + dx / length * applied,
      _point-y(from) + dy / length * applied,
    )
  }
}

#let _path-value(path) = {
  if type(path) == dictionary and path.keys().contains("path") {
    path.path
  } else {
    path
  }
}

#let _same-point(a, b) = _point-pair(a) == _point-pair(b)

#let _origin = (0, 0)

/// Build a `move` path element.
///
/// -> dictionary
#let move-to(start) = (kind: "move", start: _point-pair(start))

/// Build a `line` path element.
///
/// -> dictionary
#let line-to(end) = (kind: "line", end: _point-pair(end))

/// Build a `quad` path element.
///
/// -> dictionary
#let quad-to(control, end) = (
  kind: "quad",
  control: _point-pair(control),
  end: _point-pair(end),
)

/// Build a `cubic` path element.
///
/// -> dictionary
#let cubic-to(control-start, control-end, end) = (
  kind: "cubic",
  control-start: _point-pair(control-start),
  control-end: _point-pair(control-end),
  end: _point-pair(end),
)

/// Build a `close` path element.
///
/// -> dictionary
#let close(mode: "straight") = (kind: "close", mode: mode)

#let _normalized-element(element) = {
  if element.kind == "move" {
    move-to(element.start)
  } else if element.kind == "line" {
    line-to(element.end)
  } else if element.kind == "quad" {
    quad-to(element.control, element.end)
  } else if element.kind == "cubic" {
    cubic-to(element.control-start, element.control-end, element.end)
  } else if element.kind == "close" {
    close(mode: element.at("mode", default: "straight"))
  } else {
    panic("Unsupported Kurvst path element kind: " + str(element.kind))
  }
}

#let _part-elements(part) = {
  let value = _path-value(part)
  let elements = if type(value) == dictionary and value.keys().contains("elements") {
    value.elements
  } else if type(value) == dictionary and value.keys().contains("kind") {
    (value,)
  } else if type(value) == array {
    value
  } else {
    panic("Expected a Kurvst path, path element, or element array")
  }
  elements.map(_normalized-element)
}

#let _path-state-after(state, element) = {
  if element.kind == "move" {
    (current: element.start, subpath: element.start)
  } else if element.kind == "line" or element.kind == "quad" or element.kind == "cubic" {
    (current: element.end, subpath: state.subpath)
  } else if element.kind == "close" {
    (current: if state.subpath == none { _origin } else { state.subpath }, subpath: state.subpath)
  } else {
    state
  }
}

#let _path-state(elements) = {
  let state = (current: none, subpath: none)
  for element in elements {
    state = _path-state-after(state, element)
  }
  state
}

#let _append-elements(base-elements, parts) = {
  let result = base-elements.map(_normalized-element)
  let state = _path-state(result)
  for part in parts {
    for element in _part-elements(part) {
      if element.kind == "move" and state.current != none and _same-point(state.current, element.start) {
        ()
      } else {
        if state.current == none and (element.kind == "line" or element.kind == "quad" or element.kind == "cubic") {
          panic("Path segment elements need a current point; use line(start, end), quad(start, control, end), cubic(start, control-start, control-end, end), or move-to(start) first")
        }
        result.push(element)
        state = _path-state-after(state, element)
      }
    }
  }
  result
}

/// Build a path dictionary from an existing element array.
///
/// -> dictionary
#let from-elements(elements) = (elements: elements.map(_normalized-element))

/// Build a path dictionary from path fragments or elements.
///
/// Leading `move` elements are kept. Later `move` elements are skipped when they
/// match the current endpoint, so `path(line(a, b), line(b, c))` produces one
/// continuous subpath.
///
/// -> dictionary
#let path(..parts) = from-elements(_append-elements((), parts.pos()))

/// Return a path with additional fragments or elements appended.
///
/// If an appended fragment starts with a `move` at the current endpoint, the
/// `move` is skipped. If it starts elsewhere, the `move` begins a new subpath.
///
/// -> dictionary
#let append(path, ..parts) = from-elements(_append-elements(_elements(path), parts.pos()))

/// Build a straight-line path fragment.
///
/// -> dictionary
#let line(start, end) = path(move-to(start), line-to(end))

/// Build a quadratic path fragment.
///
/// -> dictionary
#let quad(start, control, end) = path(move-to(start), quad-to(control, end))

/// Build a cubic path fragment.
///
/// -> dictionary
#let cubic(start, control-start, control-end, end) = path(
  move-to(start),
  cubic-to(control-start, control-end, end),
)

/// Build a path fragment from a cubic segment dictionary.
///
/// -> dictionary
#let from-cubic(segment) = cubic(
  segment.start,
  segment.control-start,
  segment.control-end,
  segment.end,
)

/// Build a cubic segment dictionary for a straight line.
///
/// -> dictionary
#let line-segment(start, end) = {
  let control-start = (
    _point-x(start) + (_point-x(end) - _point-x(start)) / 3,
    _point-y(start) + (_point-y(end) - _point-y(start)) / 3,
  )
  let control-end = (
    _point-x(start) + (_point-x(end) - _point-x(start)) * 2 / 3,
    _point-y(start) + (_point-y(end) - _point-y(start)) * 2 / 3,
  )
  (
    start: _point-pair(start),
    control-start: control-start,
    control-end: control-end,
    end: _point-pair(end),
  )
}

#let _lerp(a, b, t) = a + (b - a) * t

#let _elements(path) = {
  let path = _path-value(path)
  if type(path) == dictionary and path.keys().contains("elements") {
    path.elements
  } else {
    panic("Expected a Kurvst path dictionary with an `elements` field")
  }
}

#let _quad-cubic-segment(start, control, end) = (
  start: _point-pair(start),
  control-start: (
    _point-x(start) + (_point-x(control) - _point-x(start)) * 2 / 3,
    _point-y(start) + (_point-y(control) - _point-y(start)) * 2 / 3,
  ),
  control-end: (
    _point-x(end) + (_point-x(control) - _point-x(end)) * 2 / 3,
    _point-y(end) + (_point-y(control) - _point-y(end)) * 2 / 3,
  ),
  end: _point-pair(end),
)

#let _path-cursor(current) = if current == none { _origin } else { current }

#let _path-element-end(element, current, subpath-start) = {
  if element.kind == "move" {
    element.start
  } else if element.kind == "line" or element.kind == "quad" or element.kind == "cubic" {
    element.end
  } else if element.kind == "close" {
    if subpath-start == none { _origin } else { subpath-start }
  } else {
    current
  }
}

/// Return the command elements that make up a Kurvst path.
///
/// The elements mirror Typst's native `curve` components with numeric point tuples:
/// `move`, `line`, `quad`, `cubic`, and `close`.
/// -> array
#let elements(path) = _elements(path)

/// Return the points visited by a Kurvst path's command stream.
///
/// -> array
#let points(path) = {
  let points = ()
  let current = none
  let subpath-start = none
  for element in _elements(path) {
    if element.kind == "move" {
      points.push(element.start)
      current = element.start
      subpath-start = element.start
    } else if element.kind == "line" or element.kind == "quad" or element.kind == "cubic" {
      points.push(element.end)
      current = element.end
    } else if element.kind == "close" {
      current = _path-element-end(element, current, subpath-start)
      points.push(current)
    }
  }
  points
}

/// Evaluate a cubic segment at parameter `t`.
///
/// -> dictionary
#let cubic-point(segment, t) = {
  let ab = (
    _lerp(_point-x(segment.start), _point-x(segment.control-start), t),
    _lerp(_point-y(segment.start), _point-y(segment.control-start), t),
  )
  let bc = (
    _lerp(_point-x(segment.control-start), _point-x(segment.control-end), t),
    _lerp(_point-y(segment.control-start), _point-y(segment.control-end), t),
  )
  let cd = (
    _lerp(_point-x(segment.control-end), _point-x(segment.end), t),
    _lerp(_point-y(segment.control-end), _point-y(segment.end), t),
  )
  let abc = (
    _lerp(_point-x(ab), _point-x(bc), t),
    _lerp(_point-y(ab), _point-y(bc), t),
  )
  let bcd = (
    _lerp(_point-x(bc), _point-x(cd), t),
    _lerp(_point-y(bc), _point-y(cd), t),
  )
  (
    _lerp(_point-x(abc), _point-x(bcd), t),
    _lerp(_point-y(abc), _point-y(bcd), t),
  )
}

/// Evaluate the tangent of a cubic segment at parameter `t`.
///
/// -> dictionary
#let cubic-tangent(segment, t) = {
  let ab = (
    _lerp(_point-x(segment.start), _point-x(segment.control-start), t),
    _lerp(_point-y(segment.start), _point-y(segment.control-start), t),
  )
  let bc = (
    _lerp(_point-x(segment.control-start), _point-x(segment.control-end), t),
    _lerp(_point-y(segment.control-start), _point-y(segment.control-end), t),
  )
  let cd = (
    _lerp(_point-x(segment.control-end), _point-x(segment.end), t),
    _lerp(_point-y(segment.control-end), _point-y(segment.end), t),
  )
  let abc = (
    _lerp(_point-x(ab), _point-x(bc), t),
    _lerp(_point-y(ab), _point-y(bc), t),
  )
  let bcd = (
    _lerp(_point-x(bc), _point-x(cd), t),
    _lerp(_point-y(bc), _point-y(cd), t),
  )
  (
    _point-x(bcd) - _point-x(abc),
    _point-y(bcd) - _point-y(abc),
  )
}

/// Return drawable cubic segments for any Kurvst path dictionary.
///
/// Lines and quadratics are converted to equivalent cubic segments so
/// downstream renderers can use one code path.
/// -> array
#let segments(path) = {
  let segments = ()
  let current = none
  let subpath-start = none
  for element in _elements(path) {
    if element.kind == "move" {
      current = element.start
      subpath-start = element.start
    } else if element.kind == "line" {
      let start = _path-cursor(current)
      segments.push(line-segment(start, element.end))
      current = element.end
    } else if element.kind == "quad" {
      let start = _path-cursor(current)
      segments.push(_quad-cubic-segment(start, element.control, element.end))
      current = element.end
    } else if element.kind == "cubic" {
      let start = _path-cursor(current)
      segments.push((
        start: _point-pair(start),
        control-start: _point-pair(element.control-start),
        control-end: _point-pair(element.control-end),
        end: _point-pair(element.end),
      ))
      current = element.end
    } else if element.kind == "close" and subpath-start != none {
      let start = _path-cursor(current)
      if start != subpath-start {
        segments.push(line-segment(start, subpath-start))
      }
      current = subpath-start
    }
  }
  segments
}

/// Compute the arc length of a path dictionary.
///
/// -> int | float
#let _length(path, accuracy: 0.001) = {
  cbor(_plugin.curve_path_length(cbor.encode((
    path: _path-value(path),
    accuracy: accuracy,
  ))))
}

#let length(path, accuracy: 0.001) = _length(path, accuracy: accuracy)

/// Resolve a fixed and relative visible path length.
///
/// `length` is interpreted as a fixed arc length. `ratio` is interpreted as a
/// fraction of `base-length`. `method` may be `"min"`/`"shorter"`,
/// `"max"`/`"longer"`, `"length"`/`"fixed"`, `"ratio"`/`"relative"`,
/// `"none"`/`"full"`, or a function receiving
/// `(base-length, length, ratio)`.
/// -> none | int | float
#let _resolve-length-target(base-length, length: none, ratio: none, method: "min") = {
  let fixed = if length != none and length > 0 { length } else { none }
  let relative = if ratio != none and ratio > 0 { base-length * ratio } else { none }
  if type(method) == function {
    method((base-length: base-length, length: fixed, ratio: relative))
  } else if fixed == none {
    relative
  } else if relative == none {
    fixed
  } else if method in ("max", "longer") {
    calc.max(fixed, relative)
  } else if method in ("length", "fixed") {
    fixed
  } else if method in ("ratio", "relative") {
    relative
  } else if method in ("none", "full") {
    none
  } else {
    calc.min(fixed, relative)
  }
}

#let resolve-length(base-length, length: none, ratio: none, method: "min") = {
  _resolve-length-target(base-length, length: length, ratio: ratio, method: method)
}

/// Compute the symmetric trim needed to center a shorter path layer.
///
/// -> int | float
#let center-outset(
  base-length,
  length: none,
  ratio: none,
  resolve-length: "min",
  start-outset: 0,
  end-outset: 0,
) = {
  let target = _resolve-length-target(base-length, length: length, ratio: ratio, method: resolve-length)
  if target == none {
    0
  } else {
    let visible-length = calc.max(0, base-length - start-outset - end-outset)
    calc.max(0, (visible-length - target) / 2)
  }
}

#let _side-segment(path) = {
  let segments = segments(path)
  if segments.len() == 0 { none } else { segments.at(calc.quo(segments.len(), 2)) }
}

#let _offset-toward-side-point(path, offset, side-point) = {
  if side-point == none or offset == none or offset == 0 {
    offset
  } else {
    let segment = _side-segment(path)
    if segment == none {
      offset
    } else {
      let mid = cubic-point(segment, 0.5)
      let tangent = cubic-tangent(segment, 0.5)
      let cross = (
        _point-x(tangent) * (_point-y(side-point) - _point-y(mid))
          - _point-y(tangent) * (_point-x(side-point) - _point-x(mid))
      )
      let sign = if cross < 0 { -1 } else { 1 }
      sign * calc.abs(offset)
    }
  }
}

/// Trim a path by curve distance from either end.
///
/// The returned dictionary is another path dictionary with a curve-command
/// wire path.
///
/// ```example
/// #let base = kurvst.from-cubic(demo-segment)
/// #let trimmed = kurvst.trim(base, start-outset: 0.35, end-outset: 0.3)
/// #native-scene({
///   native-cubics(kurvst.segments(base), stroke: rgb("#c8c8c8") + 0.45pt)
///   native-cubics(kurvst.segments(trimmed), stroke: rgb("#d72638") + 1pt)
/// })
/// ```
///
/// -> dictionary
#let trim(path, start-outset: 0, end-outset: 0, accuracy: 0.001) = {
  cbor(_plugin.curve_trim_path(cbor.encode((
    path: _path-value(path),
    start-outset: start-outset,
    end-outset: end-outset,
    accuracy: accuracy,
  ))))
}

/// Build an open Hobby curve through `start`, `through`, and `end`.
///
/// The returned dictionary is a path dictionary. For this three-point variant,
/// @segments returns the two cubic halves meeting at `through`.
///
/// ```example
/// #let path = kurvst.hobby-through(demo-start, demo-through, demo-end, omega: 1.2)
/// #native-scene({
///   native-cubics(kurvst.segments(path), stroke: rgb("#1b7f4c") + 0.8pt)
///   native-dot(demo-through, fill: black)
/// })
/// ```
///
/// -> dictionary
#let hobby-through(start, through, end, omega: 1.0, accuracy: 0.001) = {
  cbor(_plugin.curve_hobby_through(cbor.encode((
    start: start,
    through: through,
    end: end,
    omega: omega,
    accuracy: accuracy,
  ))))
}

/// Build an open Hobby spline through two or more points.
///
/// The returned dictionary is a path dictionary backed by a curve-command wire
/// path. Rust converts that wire path to Kurbo's `BezPath` internally, so it
/// can be passed directly to @parallel, @pattern, @trim, or
/// @to-native.
///
/// ```example
/// #let spline = kurvst.hobby-spline((
///   (0.0, 0.0),
///   (0.9, 0.8),
///   (1.8, -0.3),
///   (2.8, 0.4),
/// ))
/// #let shifted = kurvst.parallel(spline, distance: 0.14)
/// #native-scene({
///   native-cubics(kurvst.segments(spline), stroke: rgb("#c8c8c8") + 0.45pt)
///   native-cubics(kurvst.segments(shifted), stroke: rgb("#1b7f4c") + 0.8pt)
/// })
/// ```
///
/// -> dictionary
#let hobby-spline(points, omega: 1.0, accuracy: 0.001) = {
  cbor(_plugin.curve_hobby_spline(cbor.encode((
    points: points,
    omega: omega,
    accuracy: accuracy,
  ))))
}

/// Generate a sampled 1D pattern along a path.
///
/// `pattern` may be `wave()`, `zigzag()`, `coil()`, a compatible string name,
/// or a point pattern:
/// `(kind: "points", interpolation: "linear" or "smooth", points: ((at: 0, x: 0, y: 0), ...))`.
/// The whole input path is sampled continuously, so pattern phase does not
/// restart at cubic segment boundaries.
///
/// ```example
/// #let spline = kurvst.hobby-spline((
///   (0.0, 0.0),
///   (0.9, 0.8),
///   (1.8, -0.3),
///   (2.8, 0.4),
/// ))
/// #let path = kurvst.pattern(
///   spline,
///   pattern: kurvst.wave(samples-per-period: 24),
///   amplitude: 0.12,
///   wavelength: 0.55,
/// )
/// #native-scene({
///   native-cubics(kurvst.segments(spline), stroke: rgb("#c8c8c8") + 0.45pt)
///   native-cubics(kurvst.segments(path), stroke: rgb("#1b7f4c") + 0.8pt)
/// })
/// ```
///
/// -> dictionary
#let pattern(
  path,
  pattern: "wave",
  amplitude: 0.1,
  wavelength: 1.0,
  phase: 0,
  samples-per-period: 16,
  coil-longitudinal-scale: 1.25,
  anchor-start: true,
  anchor-end: true,
  accuracy: 0.001,
) = {
  let pattern = _resolve-pattern(
    pattern,
    samples-per-period: samples-per-period,
    coil-longitudinal-scale: coil-longitudinal-scale,
  )
  cbor(_plugin.curve_pattern_path(cbor.encode((
    path: _path-value(path),
    pattern: pattern,
    amplitude: amplitude,
    wavelength: wavelength,
    phase: phase,
    samples-per-period: samples-per-period,
    coil-longitudinal-scale: coil-longitudinal-scale,
    anchor-start: anchor-start,
    anchor-end: anchor-end,
    accuracy: accuracy,
  ))))
}

/// Generate a parallel path for a path.
///
/// Pass any path returned by this module, such as @hobby-spline. The result is a
/// drawable path dictionary with the fitted parallel path and its curve-command
/// wire path.
///
/// ```example
/// #let spline = kurvst.hobby-spline((
///   (0.0, 0.0),
///   (0.9, 0.8),
///   (1.8, -0.3),
///   (2.8, 0.4),
/// ))
/// #let parallel = kurvst.parallel(spline, distance: 0.14)
/// #native-scene({
///   native-cubics(kurvst.segments(spline), stroke: rgb("#c8c8c8") + 0.45pt)
///   native-cubics(kurvst.segments(parallel), stroke: rgb("#1b7f4c") + 0.8pt)
/// })
/// ```
///
/// -> dictionary
#let parallel(
  path,
  distance: 0,
  start-outset: 0,
  end-outset: 0,
  accuracy: 0.001,
  optimize: true,
) = {
  cbor(_plugin.curve_parallel_path(cbor.encode((
    path: _path-value(path),
    distance: distance,
    start-outset: start-outset,
    end-outset: end-outset,
    accuracy: accuracy,
    optimize: optimize,
  ))))
}

/// Build a derived visible path layer.
///
/// `layer` combines the common operations needed by drawing packages:
/// optional side-aware offsetting, endpoint trimming, and centered shortening by
/// a fixed `length`, a relative `ratio`, or both. The return value is a normal
/// Kurvst path dictionary that can be passed to @pattern, @parallel,
/// @trim, or @to-cetz.
///
/// ```example
/// #let base = kurvst.hobby-spline((
///   (0.0, 0.0),
///   (0.9, 0.8),
///   (1.8, -0.3),
///   (2.8, 0.4),
/// ))
/// #let layer = kurvst.layer(
///   base,
///   offset: 0.16,
///   length: 1.6,
///   ratio: 0.5,
///   resolve-length: "min",
/// )
/// #native-scene({
///   native-cubics(kurvst.segments(base), stroke: rgb("#c8c8c8") + 0.45pt)
///   native-cubics(kurvst.segments(layer), stroke: rgb("#1b7f4c") + 0.8pt)
/// })
/// ```
///
/// -> dictionary
#let layer(
  path,
  offset: 0,
  length: none,
  ratio: none,
  resolve-length: "min",
  start-outset: 0,
  end-outset: 0,
  side-point: none,
  accuracy: 0.001,
  optimize: true,
) = {
  let distance = _offset-toward-side-point(path, offset, side-point)
  let center-trim = center-outset(
    _length(path, accuracy: accuracy),
    length: length,
    ratio: ratio,
    resolve-length: resolve-length,
    start-outset: start-outset,
    end-outset: end-outset,
  )
  let start-outset = start-outset + center-trim
  let end-outset = end-outset + center-trim
  if distance == none or distance == 0 {
    trim(path, start-outset: start-outset, end-outset: end-outset, accuracy: accuracy)
  } else {
    parallel(
      path,
      distance: distance,
      start-outset: start-outset,
      end-outset: end-outset,
      accuracy: accuracy,
      optimize: optimize,
    )
  }
}

/// Emit a Kurvst path as native Typst `curve` content.
///
/// -> content
#let to-native(path, unit: 1, ..style) = {
  let components = ()
  for element in _elements(path) {
    if element.kind == "move" {
      components.push(curve.move(_point(element.start, unit: unit)))
    } else if element.kind == "line" {
      components.push(curve.line(_point(element.end, unit: unit)))
    } else if element.kind == "quad" {
      components.push(curve.quad(_point(element.control, unit: unit), _point(element.end, unit: unit)))
    } else if element.kind == "cubic" {
      components.push(curve.cubic(
        _point(element.control-start, unit: unit),
        _point(element.control-end, unit: unit),
        _point(element.end, unit: unit),
      ))
    } else if element.kind == "close" {
      components.push(curve.close(mode: element.at("mode", default: "straight")))
    }
  }

  if components.len() == 0 {
    ()
  } else {
    curve(..style.named(), ..components)
  }
}

/// Emit a Kurvst path as CeTZ path data.
///
/// The result has CeTZ's subpath shape: `(origin, closed, segments)`, where
/// line segments are `("l", end)` and cubics are `("c", control-start,
/// control-end, end)`.
/// -> array
#let to-cetz-data(path, unit: 1) = {
  let subpaths = ()
  let origin = none
  let current = none
  let segments = ()
  let closed = false

  for element in _elements(path) {
    if element.kind == "move" {
      if origin != none {
        subpaths.push((_point(origin, unit: unit), closed, segments))
      }
      origin = element.start
      current = element.start
      segments = ()
      closed = false
    } else {
      if origin == none {
        origin = _origin
        current = _origin
        segments = ()
        closed = false
      }

      if element.kind == "line" {
        segments.push(("l", _point(element.end, unit: unit)))
        current = element.end
      } else if element.kind == "quad" {
        let cubic = _quad-cubic-segment(_path-cursor(current), element.control, element.end)
        segments.push((
          "c",
          _point(cubic.control-start, unit: unit),
          _point(cubic.control-end, unit: unit),
          _point(cubic.end, unit: unit),
        ))
        current = element.end
      } else if element.kind == "cubic" {
        segments.push((
          "c",
          _point(element.control-start, unit: unit),
          _point(element.control-end, unit: unit),
          _point(element.end, unit: unit),
        ))
        current = element.end
      } else if element.kind == "close" {
        closed = true
        subpaths.push((_point(origin, unit: unit), closed, segments))
        current = origin
        origin = none
        segments = ()
        closed = false
      }
    }
  }

  if origin != none {
    subpaths.push((_point(origin, unit: unit), closed, segments))
  }
  subpaths
}

/// Draw a path dictionary through CeTZ.
///
/// ```example
/// #let path = kurvst.parallel(kurvst.from-cubic(demo-segment), distance: 0.18)
/// #cetz.canvas({
///   kurvst.to-cetz(path, stroke: rgb("#355c9a") + 0.55pt)
/// })
/// ```
///
/// -> content
#let to-cetz(path, unit: 1, ..style) = {
  let style = style.named()
  let subpaths = to-cetz-data(path, unit: unit)
  if subpaths.len() == 0 {
    ()
  } else {
    cetz.draw.merge-path({
      for (origin, closed, segments) in subpaths {
        let current = origin
        for (kind, ..args) in segments {
          if kind == "l" {
            cetz.draw.line(current, args.last())
          } else if kind == "c" {
            cetz.draw.bezier(current, args.last(), args.at(0), args.at(1))
          }
          current = args.last()
        }
        if closed and current != origin {
          cetz.draw.line(current, origin)
        }
      }
    }, ..style)
  }
}

/// Split a path through a point sequence into per-span paths.
///
/// The returned dictionary has `parts` and `curve`. Each entry in `parts` is
/// the path between two consecutive points, and all parts join smoothly through
/// the input points.
///
/// ```example
/// #let split = kurvst.split-through(
///   (demo-start, demo-through, demo-end),
///   start-outset: 0.2,
///   end-outset: 0.2,
/// )
/// #native-scene({
///   native-cubics(kurvst.segments(split.parts.at(0)), stroke: rgb("#d72638") + 0.9pt)
///   native-cubics(kurvst.segments(split.parts.at(1)), stroke: rgb("#355c9a") + 0.9pt)
///   native-dot(demo-through, fill: black)
/// })
/// ```
///
/// -> dictionary
#let split-through(points, omega: 1.0, start-outset: 0, end-outset: 0, accuracy: 0.001) = {
  if points.len() < 2 {
    panic("split-through expects at least two points")
  }

  let curve = hobby-spline(points, omega: omega, accuracy: accuracy)
  let split-segments = segments(curve)
  let parts = ()
  let last-index = split-segments.len() - 1
  for (index, segment) in split-segments.enumerate() {
    let part-start-outset = if index == 0 { start-outset } else { 0 }
    let part-end-outset = if index == last-index { end-outset } else { 0 }
    parts.push(trim(
      from-cubic(segment),
      start-outset: part-start-outset,
      end-outset: part-end-outset,
      accuracy: accuracy,
    ))
  }

  (
    parts: parts,
    curve: curve,
  )
}
