// Shape vocabulary and parameter definitions follow tiptoe 0.4 (MIT licensed,
// shape API). Independently implemented; no CeTZ code.

#let _plugin = plugin("../kurvst.wasm")

#let fields = (
  triangle: ("length", "width", "inset", "fill", "stroke", "rev"),
  stealth: ("length", "width", "inset", "fill", "stroke", "rev"),
  round: ("length", "width", "inset", "fill", "stroke", "rev"),
  straight: ("length", "width", "stroke", "rev"),
  tikz: ("width", "stroke"),
  barb: ("width", "arc", "stroke", "rev"),
  hooks: ("width", "arc", "stroke", "rev"),
  bar: ("width", "stroke", "align"),
  bracket: ("length", "width", "stroke", "rev"),
  circle: ("length", "width", "fill", "stroke", "align"),
  square: ("length", "width", "fill", "stroke", "align"),
  diamond: ("length", "width", "fill", "stroke", "align"),
  rays: ("length", "n", "phase", "stroke", "align"),
  combine: ("parts",),
)

#let finite(value, name) = {
  assert(
    type(value) in (int, float)
      and value == value
      and calc.abs(value) < float("inf"),
    message: "kurvst.mark: " + name + " must be finite",
  )
  float(value)
}

#let dimension(value) = {
  let result = if type(value) == length {
    (points: value / 1pt, ratio: 0.0)
  } else if type(value) == ratio {
    (points: 0.0, ratio: value / 100%)
  } else if type(value) == relative {
    (points: value.length / 1pt, ratio: value.ratio / 100%)
  } else {
    panic(
      "kurvst.mark: dimensions must be fixed lengths, ratios, or mixed relative lengths",
    )
  }
  (
    points: finite(result.points, "dimension points"),
    ratio: finite(result.ratio, "dimension ratio"),
  )
}

// Validation and normalization share one boundary, also for hand-written specs.
#let normalize(spec, child: false) = {
  assert(
    type(spec) == dictionary
      and spec.at("kind", default: none) == "kurvst-mark",
    message: "kurvst.mark: replace CeTZ mark dictionaries with kurvst.mark constructors",
  )
  let shape = spec.at("shape", default: none)
  assert(shape in fields, message: "kurvst.mark: unknown shape")
  let result = (shape: shape)
  for (key, value) in spec {
    if key in ("kind", "shape") { continue }
    assert(
      key in fields.at(shape) or key in ("fit", "shorten"),
      message: "kurvst.mark."
        + shape
        + ": unsupported option "
        + key
        + "; use kurvst.mark options instead of CeTZ symbol/scale/anchor/shorten-to",
    )
    assert(
      not child or key not in ("fit", "shorten"),
      message: "kurvst.mark.combine: fit and shorten belong to the composite, not its children",
    )
    if (
      value == auto
        and (
          key in ("fill", "stroke", "width")
            or (key == "phase" and shape == "rays")
            or (key == "length" and shape == "bracket")
        )
    ) { continue }
    result.insert(key, if key in ("length", "width") {
      dimension(value)
    } else if key in ("inset", "shorten") {
      assert(
        type(value) == ratio,
        message: "kurvst.mark: " + key + " must be a ratio",
      )
      finite(value / 100%, key)
    } else if key in ("arc", "phase") {
      assert(
        type(value) == angle,
        message: "kurvst.mark: " + key + " must be an angle",
      )
      finite(value / 1rad, key)
    } else if key == "rev" {
      assert(type(value) == bool, message: "kurvst.mark: rev must be boolean")
      value
    } else if key == "n" {
      assert(
        type(value) == int and value > 0,
        message: "kurvst.mark: n must be a positive integer",
      )
      value
    } else if key == "align" {
      assert(
        value in ("center", "end"),
        message: "kurvst.mark: align must be center or end",
      )
      value
    } else if key == "fit" {
      assert(
        value in ("chord", "bend"),
        message: "kurvst.mark: fit must be chord or bend",
      )
      value
    } else if key in ("fill", "stroke") {
      if key == "stroke" and type(value) == stroke {
        assert(
          value.thickness == auto,
          message: "kurvst.mark: set line thickness on the line, not the mark stroke",
        )
        assert(
          value.cap == auto
            and value.join == auto
            and value.dash == auto
            and value.miter-limit == auto,
          message: "kurvst.mark: head stroke is paint-only; geometry styles belong to the shape",
        )
        assert(
          value.paint == auto or type(value.paint) in (color, gradient, tiling),
          message: "kurvst.mark: stroke must contain valid paint",
        )
      } else {
        assert(
          value == none or type(value) in (bool, color, gradient, tiling),
          message: "kurvst.mark: "
            + key
            + " must be paint, none, auto, or boolean",
        )
      }
      if type(value) == bool { value } else { value != none }
    } else if key == "parts" {
      assert(
        type(value) == array and value.len() > 0,
        message: "kurvst.mark.combine: parts must be a nonempty array",
      )
      value.map(part => if type(part) == dictionary {
        normalize(part, child: true)
      } else {
        (gap: dimension(part))
      })
    })
  }
  assert(
    shape != "combine" or "parts" in result,
    message: "kurvst.mark.combine: parts is required",
  )
  result
}

#let construct(shape, args) = {
  assert(
    args.pos().len() == 0,
    message: "kurvst.mark." + shape + ": options must be named",
  )
  for key in args.named().keys() {
    assert(
      key in fields.at(shape) or key in ("fit", "shorten"),
      message: "kurvst.mark."
        + shape
        + ": unsupported option "
        + key
        + "; replace CeTZ options with kurvst.mark options",
    )
  }
  let spec = (kind: "kurvst-mark", shape: shape, ..args.named())
  let numeric = normalize(spec)
  for key in args.named().keys() {
    if not (key in numeric) { let _ = spec.remove(key) }
  }
  spec
}

#let triangle(..args) = construct("triangle", args)
#let straight(..args) = construct("straight", args)
#let stealth(..args) = construct("stealth", args)
#let round(..args) = construct("round", args)
#let tikz(..args) = construct("tikz", args)
#let barb(..args) = construct("barb", args)
#let hooks(..args) = construct("hooks", args)
#let bar(..args) = construct("bar", args)
#let bracket(..args) = construct("bracket", args)
#let circle(..args) = construct("circle", args)
#let square(..args) = construct("square", args)
#let diamond(..args) = construct("diamond", args)
#let rays(..args) = construct("rays", args)
#let combine(..args) = {
  for key in args.named().keys() {
    assert(
      key in ("fit", "shorten"),
      message: "kurvst.mark.combine: only fit and shorten are named; supply parts positionally",
    )
  }
  construct("combine", arguments(parts: args.pos(), ..args.named()))
}

/// Prepare CBOR-safe engine data without resolving shared Rust defaults.
/// Drawing context belongs to the batch template, not the shape specification.
#let prepare(spec) = normalize(spec)

/// Prepare and place a batch of marks using the shared native/Wasm engine.
///
/// Each template is `(mark: spec, context: (units-per-pt: ..., line-thickness: ...))`.
/// Context may include `shaft-stroke` with cap, join, and miter-limit.
/// Carriers are Kurvst paths; placements refer to template/carrier indices.
/// Paint stays in the original specs and is matched through the returned indices.
#let geometry(templates, carriers, placements, mode: "candidates", format: "native") = {
  assert(
    format in ("native", "cbor"),
    message: "kurvst.mark.geometry: format must be native or cbor",
  )
  assert(
    mode in ("candidates", "selected"),
    message: "kurvst.mark.geometry: mode must be candidates or selected",
  )
  assert(
    type(templates) == array
      and type(carriers) == array
      and type(placements) == array,
    message: "kurvst.mark.geometry: templates, carriers, and placements must be arrays",
  )
  let templates = templates.map(template => {
    assert(
      type(template) == dictionary
        and template.keys().sorted() == ("context", "mark"),
      message: "kurvst.mark.geometry: a template needs exactly mark and context",
    )
    let ctx = template.at("context")
    assert(
      type(ctx) == dictionary,
      message: "kurvst.mark.geometry: context must be a dictionary",
    )
    assert(
      finite(ctx.at("units-per-pt"), "units-per-pt") > 0,
      message: "kurvst.mark.geometry: units-per-pt must be positive",
    )
    assert(
      finite(ctx.at("line-thickness"), "line-thickness") >= 0,
      message: "kurvst.mark.geometry: line-thickness must be nonnegative",
    )
    (
      mark: prepare(template.mark),
      "context": ctx,
    )
  })
  let carriers = carriers.map(carrier => {
    assert(
      type(carrier) == dictionary,
      message: "kurvst.mark.geometry: carriers must be Kurvst paths",
    )
    (path: carrier.at("path", default: carrier))
  })
  let request = cbor.encode((
    mode: mode,
    templates: templates,
    carriers: carriers,
    placements: placements,
  ))
  if format == "cbor" {
    _plugin.mark_geometry_packed(request)
  } else {
    cbor(_plugin.mark_geometry(request))
  }
}
