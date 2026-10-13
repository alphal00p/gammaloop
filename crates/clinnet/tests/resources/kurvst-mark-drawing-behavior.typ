#import "/crates/linnest/typst/src/curve.typ" as curve
#import "/crates/linnest/typst/src/impl/draw.typ" as drawing
#import "/crates/linnest/typst/src/lib.typ": draw, graph
#import "@preview/cetz:0.5.2" as cetz

#set page(width: auto, height: auto, margin: 0pt)

#let empty-batch = (
  templates: (), template-keys: (), carriers: (), placements: (), plans: (),
)
#let path = curve.cubic((0, 0), (0, 2), (3, -1), (3, 0))

// Coordinate normalization must preserve closure origins and local previous
// positions. Hooks observe synthetic cubic coordinates in reverse hook order.
#cetz.canvas({
  (ctx => {
    let style = (mark: curve.mark.triangle(), stroke: black + 1pt)
    let original = curve.from-elements((
      curve.move-to((-0.0, 0)), curve.quad-to((1, 2), (3, 0)),
      curve.move-to((3, 0)), curve.line-to((4, 1)),
      curve.cubic-to((5, 2), (6, 2), (7, 0)), curve.close(),
    ))
    let transform = (
      (2, 0.5, 4, 10), (0, -3, 2, 20),
      (0, 0, 1, 0), (0, 0, 0, 1),
    )
    let transformed = ctx + (transform: transform)
    for unit in (2, 0, -2, 1pt) {
      let packet = drawing._mark-carrier-input(original, style + (unit: unit))
      let plan = drawing._prepare-pending-mark(transformed, packet)
      // An identity hook forces the general resolver without changing geometry.
      let reference = drawing._prepare-pending-mark(transformed + (
        resolve-coordinate: (((ctx, point) => point),),
      ), packet)
      assert.eq(plan.carrier, reference.carrier)
      assert.eq(plan.template, reference.template)
      assert.eq(plan.placement, reference.placement)
      assert.eq(plan.drawable, reference.drawable)
      assert.eq(plan.after-ctx.prev, reference.after-ctx.prev)
      reference.after-ctx.resolve-coordinate = plan.after-ctx.resolve-coordinate
      assert.eq(plan.after-ctx, reference.after-ctx)
      let expected = drawing._map-path(original, point =>
        drawing._canvas-point(transform, (
          ..point.map(value => cetz.util.resolve-number(ctx, value * unit)), 0.0,
        )))
      if unit != 0 {
        assert.eq(curve.segments(plan.carrier).len(), curve.segments(expected).len())
        for (actual, wanted) in curve.segments(plan.carrier).zip(curve.segments(expected)) {
          for key in ("start", "control-start", "control-end", "end") {
            assert(drawing._point-distance(actual.at(key), wanted.at(key)) < 1e-12)
          }
        }
      }
      assert.eq(plan.after-ctx.prev.pt,
        (cetz.util.resolve-number(ctx, 3 * unit), 0.0, 0.0))
      assert.eq(plan.preparation.key.root, ())
    }
    let derived = drawing._derived-path-elements(original, style, 0.7, false, true)
    assert.eq(derived, (drawing._mark-carrier-input(
      drawing._segments-path(curve.segments(original)), style,
      phase: 0.7, anchor-start: false, anchor-end: true),))
    let packet = derived.first()
    assert.eq(curve.segments(packet.path), curve.segments(original))
    assert.eq((packet.phase, packet.anchor-start, packet.anchor-end), (0.7, false, true))
    assert.eq(drawing._segments-elements(curve.segments(original),
      style + (label-only: true), auto, true, true), ())
    assert.eq(drawing._derived-segments-elements((), style, auto, true, true), ())
    for empty in (curve.from-elements(()), curve.path(curve.move-to((2, 3)))) {
      assert.eq(drawing._derived-path-elements(empty, style, auto, true, true), ())
      let plan = drawing._prepare-pending-mark(ctx,
        drawing._mark-carrier-input(empty, style))
      let point = curve.points(empty).at(0, default: (0, 0))
      assert.eq(curve.points(plan.carrier), (point,))
      assert.eq(plan.after-ctx.prev.pt, (..point, 0.0))
    }
    // A real gap must remain a gap for the engine to reject, not a new line.
    let disconnected = curve.path(
      curve.line((0, 0), (1, 0)), curve.line((2, 0), (3, 0)))
    let gap = drawing._prepare-pending-mark(ctx,
      drawing._mark-carrier-input(disconnected, style))
    assert.eq(curve.segments(gap.carrier), curve.segments(disconnected))
    let hook-context = transformed + (
      prev: ctx.prev + (pt: (7.0, 8.0, 9.0)),
      resolve-coordinate: (
        ((ctx, point) => {
          ctx.insert("hook-local-mutation", true)
          (point.at(0) + ctx.prev.pt.at(0), point.at(1), 3)
        }),
        ((ctx, point) => (2 * point.at(0), point.at(1))),
      ),
    )
    let hooked = drawing._prepare-pending-mark(hook-context,
      drawing._mark-carrier-input(curve.line((0, 0), (3, 0)), style))
    assert.eq(curve.segments(hooked.carrier).first(), (
      start: (36.0, 26.0), control-start: (40.0, 26.0),
      control-end: (48.0, 26.0), end: (60.0, 26.0),
    ))
    assert.eq(hooked.after-ctx.prev.pt, (19.0, 0.0, 3.0))
    assert(not ("hook-local-mutation" in hooked.after-ctx))
    // Root selection is based on drawable cubic count, including implicit
    // closure, not on the retained command or control representation.
    let inherited = ctx + (style: ctx.style + (bezier: ctx.style.bezier + (
      stroke: blue + 2pt, fill: orange,
    )))
    for (path, close, root) in (
      (curve.line((0, 0), (3, 0)), false, "bezier"),
      (curve.quad((0, 0), (1, 2), (3, 0)), true, "bezier"),
      (curve.path(curve.line((0, 0), (3, 0)), curve.close()), false, ()),
    ) {
      let plan = drawing._prepare-pending-mark(inherited,
        drawing._mark-carrier-input(path, (mark: style.mark, close: close)))
      let expected = cetz.styles.resolve(inherited.style,
        merge: (mark: none, close: close), root: root)
      assert.eq(plan.preparation.key.root, root)
      assert.eq(plan.drawable.stroke, expected.stroke)
      assert.eq(plan.drawable.fill, expected.fill)
    }
    (ctx: ctx, drawables: ())
  },)
})

#let legacy = sys.inputs.at("legacy", default: none)
#if legacy != none {
  let obsolete = (
    symbol: (symbol: ">"),
    ends: (start: "<", end: ">"),
    scale: (symbol: "triangle", scale: 0.5),
    callback: (symbol: (ctx => none)),
  ).at(legacy)
  let _ = drawing._mark-carrier-input(path, (stroke: black + 1pt, mark: obsolete))
}

#cetz.canvas({
  (ctx => {
    let packets = (
      drawing._mark-carrier-input(path, (
        stroke: red + 1pt,
        mark: curve.mark.triangle(length: 8pt, width: 6pt),
      )),
      drawing._mark-carrier-input(path, (
        stroke: blue + 1pt,
        mark: curve.mark.triangle(length: 8pt, width: 6pt),
      )),
      drawing._mark-carrier-input(path, (
        stroke: black + 1pt,
        mark: curve.mark.combine(
          curve.mark.bar(stroke: green),
          2pt,
          curve.mark.circle(fill: orange),
          fit: "bend",
        ),
        pattern: "wave", pattern-amplitude: 0.05, pattern-wavelength: 0.4,
      )),
    )
    let registered = drawing._register-mark-elements(ctx, packets, empty-batch)
    let batch = registered.batch
    // Paint does not participate in numeric template identity.
    assert.eq(batch.templates.len(), 2)
    assert.eq(batch.placements.len(), 3)
    let candidates = curve.mark.geometry(batch.templates, batch.carriers, batch.placements)
    let selected = curve.mark.geometry(batch.templates, batch.carriers, batch.placements, mode: "selected")
    assert.eq(selected.shafts.len(), 3)
    for mark in selected.marks {
      assert.eq(curve.elements(mark.shaft), ())
    }
    let finished = drawing._materialize-mark-groups(
      registered.elements, batch, selected, selected: true,
    ).flatten()
    let painted = cetz.process.many(ctx, finished, compute-bounds: false)
    let heads = painted.drawables.filter(drawable => cetz.drawable.TAG.mark in drawable.tags)
    assert.eq(heads.len(), 4)
    assert.eq(heads.at(0).fill, red)
    assert.eq(heads.at(1).fill, blue)
    assert.eq(heads.at(2).stroke.paint, green)
    assert.eq(heads.at(3).fill, orange)
    assert.eq(heads.at(2).stroke.cap, selected.marks.at(2).paths.at(0).cap)
    let shafts = painted.drawables.filter(drawable =>
      cetz.drawable.TAG.mark not in drawable.tags and cetz.drawable.TAG.hidden not in drawable.tags)
    assert.eq(shafts.len(), 3)
    // The wave follows the bent authoritative shaft, including its end contact.
    let painted-end = shafts.last().segments.last().last().last().last().slice(0, 2)
    let engine-end = curve.points(selected.shafts.last().shaft).last()
    assert(drawing._point-distance(painted-end, engine-end) < 1e-9)
    for (plan, mark) in batch.plans.zip(candidates.marks) {
      let bounds = drawing._engine-mark-bounds(plan, mark)
      assert(bounds.len() > 0)
      assert(bounds.all(box => box.left <= box.right and box.bottom <= box.top))
    }
    (ctx: ctx, drawables: painted.drawables)
  },)
})

// Invisible paint and sizing thickness are independent. Neither a plain nor
// patterned invisible shaft may repel labels, even when its head is visible.
#cetz.canvas(length: 1pt, {
  (ctx => {
    for pattern in (none, "wave") {
      let style = (
        stroke: (paint: none, thickness: 20pt),
        mark: curve.mark.triangle(length: 10pt, width: 8pt, fill: red),
      )
      if pattern != none {
        style += (pattern: pattern, pattern-amplitude: 2, pattern-wavelength: 8)
      }
      let packet = drawing._mark-carrier-input(curve.line((0, 0), (100, 0)), style)
      let registered = drawing._register-mark-elements(ctx, (packet,), empty-batch)
      let batch = registered.batch
      let geometry = curve.mark.geometry(batch.templates, batch.carriers, batch.placements)
      let boxes = drawing._engine-mark-bounds(batch.plans.first(), geometry.marks.first())
      assert(boxes.len() > 0 and boxes.all(box => box.left >= 90 - 1e-9),
        message: "unpainted shaft thickness is sizing data, not a candidate obstacle")
      let painted = cetz.process.many(ctx,
        drawing._materialize-mark-groups(registered.elements, batch, geometry).flatten(),
        compute-bounds: false)
      assert.eq(drawing._edge-collision-lines(ctx, painted.drawables, painted: true), (),
        message: "unpainted shafts produce no collision lines")
      assert.eq(drawing._annotation-drawable-bounds(ctx,
        (cetz.drawable.path(curve.to-cetz-data(batch.carriers.first()),
          stroke: (paint: none, thickness: 20pt), fill: none),)), (),
        message: "ordinary unpainted drawables also produce no obstacles")
    }
    // Stroke expansion must reach root/draw-after bounds, not just candidate
    // footprints. The bar centerline is at x=100; its painted edge is x=110.
    let packet = drawing._mark-carrier-input(curve.line((0, 0), (100, 0)), (
      stroke: black + 20pt, mark: curve.mark.bar(width: 20pt),
    ))
    let registered = drawing._register-mark-elements(ctx, (packet,), empty-batch)
    let batch = registered.batch
    let geometry = curve.mark.geometry(
      batch.templates, batch.carriers, batch.placements, mode: "selected")
    let painted = cetz.process.many(ctx,
      drawing._materialize-mark-groups(registered.elements, batch, geometry, selected: true).flatten())
    assert(calc.abs(painted.bounds.high.at(0) - 110) < 1e-9,
      message: "root bounds include the engine's painted head outline")
    (ctx: ctx, drawables: painted.drawables)
  },)
})

#let constructors = (
  curve.mark.triangle, curve.mark.straight, curve.mark.stealth,
  curve.mark.round, curve.mark.tikz, curve.mark.barb, curve.mark.hooks,
  curve.mark.bar, curve.mark.bracket, curve.mark.circle, curve.mark.square,
  curve.mark.diamond, curve.mark.rays,
)
#let painted-center(mark) = {
  // This test needs exact extrema, not the production collision adapter's
  // conservative control hulls (stroked circles may exceed those extrema).
  let points = mark.paths.fold((), (all, path) =>
    all + cetz.path-util.bounds(curve.to-cetz-data(path.outline)))
  let left = calc.min(..points.map(point => point.at(0)))
  let right = calc.max(..points.map(point => point.at(0)))
  ((left + right) / 2, 0)
}
// Numeric placements from the actual Typst drawing preparation must center
// painted outlines, not merely the origin/tip or a triangle-specific size.
#cetz.canvas({
  (ctx => {
    let specs = (
      curve.mark.triangle(), curve.mark.circle(), curve.mark.square(),
      curve.mark.diamond(), curve.mark.bar(),
      curve.mark.triangle(length: 3pt + 200%, width: 2pt + 150%, rev: true),
      curve.mark.combine(curve.mark.circle(), 2pt, curve.mark.bar()),
    )
    let packets = ()
    for spec in specs {
      for direction in ("forward", "backward") {
        packets.push(drawing._mark-carrier-input(curve.line((0, 0), (3, 0)), (
          stroke: black + 0.8pt, mark: spec,
          mark-position: "center", mark-direction: direction, mark-shift: -0.07,
        )))
      }
    }
    let registered = drawing._register-mark-elements(ctx, packets, empty-batch)
    let batch = registered.batch
    assert.eq(batch.placements.len(), 14)
    assert(batch.placements.all(placement =>
      placement.station == (kind: "ratio", value: 0.5) and placement.shift == -0.07))
    let selected = curve.mark.geometry(batch.templates, batch.carriers, batch.placements, mode: "selected")
    for head in selected.marks {
      assert(drawing._point-distance(painted-center(head), (1.43, 0)) < 1e-9)
      assert.eq(curve.elements(head.shaft), ())
    }
    for shaft in selected.shafts {
      assert(drawing._point-distance(curve.points(shaft.shaft).first(), (0, 0)) < 1e-9)
      assert(drawing._point-distance(curve.points(shaft.shaft).last(), (3, 0)) < 1e-9)
    }
    let painted = cetz.process.many(ctx, drawing._materialize-mark-groups(
      registered.elements, batch, selected, selected: true,
    ).flatten(), compute-bounds: false)
    (ctx: ctx, drawables: painted.drawables)
  },)
})

#cetz.canvas({
  (ctx => {
    let mark-style = (stroke: black + 1pt,
      mark: curve.mark.triangle(length: 4pt, width: 3pt))
    let packets = (
      drawing._mark-carrier-input(curve.line((0, 0), (3, 0)),
        mark-style + (mark-direction: "backward")),
      drawing._mark-carrier-input(curve.line((0, 0), (3, 0)), mark-style),
    )
    let registered = drawing._register-mark-elements(ctx, packets, empty-batch)
    let batch = registered.batch
    // Exercise grouped selected shafts explicitly, not merely the one-mark
    // carrier case used by ordinary single-spec drawing styles.
    batch.carriers = (batch.carriers.first(),)
    batch.placements = batch.placements.map(placement => placement + (carrier: 0))
    let selected = curve.mark.geometry(batch.templates, batch.carriers, batch.placements, mode: "selected")
    assert.eq(selected.shafts.len(), 1)
    let painted = cetz.process.many(ctx, drawing._materialize-mark-groups(
      registered.elements, batch, selected, selected: true,
    ).flatten(), compute-bounds: false)
    assert.eq(painted.drawables.filter(drawable => cetz.drawable.TAG.hidden not in drawable.tags).len(), 3)
    assert.eq(painted.drawables.filter(drawable => cetz.drawable.TAG.mark in drawable.tags).len(), 2)
    (ctx: ctx, drawables: painted.drawables)
  },)
})

#cetz.canvas({
  (ctx => {
    let packets = constructors.enumerate().map(((index, constructor)) =>
      drawing._mark-carrier-input(curve.line((0, index), (3, index)), (
        stroke: black + 0.8pt,
        mark: constructor(),
      )))
    let registered = drawing._register-mark-elements(ctx, packets, empty-batch)
    let batch = registered.batch
    let selected = curve.mark.geometry(batch.templates, batch.carriers, batch.placements, mode: "selected")
    let painted = cetz.process.many(ctx, drawing._materialize-mark-groups(
      registered.elements, batch, selected, selected: true,
    ).flatten(), compute-bounds: false)
    assert.eq(selected.shafts.len(), 13)
    assert.eq(painted.drawables.filter(drawable => cetz.drawable.TAG.hidden not in drawable.tags).len(), 26)
    for (head, mark) in painted.drawables.filter(drawable =>
      cetz.drawable.TAG.mark in drawable.tags).zip(selected.marks) {
      if head.stroke != none {
        assert.eq(head.stroke.join, mark.paths.first().join)
        assert.eq(head.stroke.cap, mark.paths.first().cap)
        assert.eq(head.stroke.miter-limit, mark.paths.first().miter-limit)
      }
    }
    (ctx: ctx, drawables: painted.drawables)
  },)
})

#cetz.canvas({
  (ctx => {
    let line = curve.line((0, 0), (10, 0))
    let style = (
      stroke: black + 0.6pt, crossing-gap: 2,
      pattern: "wave", pattern-amplitude: 0.05,
      pattern-wavelength: 1, pattern-phase: 0.25,
    )
    let mark-style = style + (
      mark: curve.mark.triangle(length: 4pt, width: 3pt),
      mark-position: 0.5,
    )
    let pieces = drawing._cut-path-elements(line, style, (4,), mark-style: mark-style)
    let packet = pieces.find(drawing._pending-mark)
    // Visible lengths are 3 + 5: their halfway station is x=6, not x=5.
    assert.eq(curve.points(packet.path).first(), (5.0, 0.0))
    assert.eq(curve.points(packet.path).last(), (10.0, 0.0))
    assert.eq(packet.style.mark-position, 0.2)
    assert.eq(packet.phase, 0.25 + 10 * calc.pi)
    assert.eq(packet.anchor-start, false)
    assert.eq(packet.anchor-end, true)
    let registered = drawing._register-mark-elements(ctx, pieces, empty-batch)
    let batch = registered.batch
    let selected = curve.mark.geometry(batch.templates, batch.carriers, batch.placements, mode: "selected")
    // The public layer path is cubicized before the shared arc-length solve.
    assert(drawing._point-distance(painted-center(selected.marks.first()), (6, 0)) < 1e-9)
    let ends = drawing._cut-path-elements(line, style, (4,),
      mark-style: mark-style + (mark-position: "end"))
    let end-packet = ends.find(drawing._pending-mark)
    assert.eq(curve.points(end-packet.path).first(), (5.0, 0.0))
    let finished = drawing._materialize-mark-groups(registered.elements, batch, selected, selected: true).flatten()
    let painted = cetz.process.many(ctx, finished, compute-bounds: false)
    assert.eq(painted.drawables.filter(drawable => cetz.drawable.TAG.hidden not in drawable.tags).len(), 3)
    let halves = (
      whole: curve.segments(line),
      source: curve.segments(curve.line((0, 0), (5, 0))),
      sink: curve.segments(curve.line((5, 0), (10, 0))),
    )
    let same = drawing._pattern-edge-halves(halves, mark-style, mark-style)
    assert.eq(same.filter(drawing._pending-mark).len(), 1)
    assert(drawing._point-distance(curve.points(same.first().path).last(), (10, 0)) < 1e-9)
    let different = drawing._pattern-edge-halves(halves,
      mark-style + (stroke: red + 0.6pt),
      mark-style + (stroke: blue + 0.6pt, mark: curve.mark.circle(fill: green)))
    // Distinct half styles retain both marks and their separate shaft paints.
    assert.eq(different.filter(drawing._pending-mark).len(), 2)
    let halves-prepared = drawing._register-mark-elements(ctx, different, empty-batch)
    let halves-batch = halves-prepared.batch
    let halves-selected = curve.mark.geometry(halves-batch.templates,
      halves-batch.carriers, halves-batch.placements, mode: "selected")
    assert.eq(halves-selected.shafts.len(), 2)
    (ctx: ctx, drawables: painted.drawables)
  },)
})

#cetz.canvas({
  (ctx => {
    let transformed = ctx + (transform: cetz.matrix.mul-mat(
      cetz.matrix.transform-translate(2, 1, 0),
      cetz.matrix.transform-rotate-z(25deg),
      cetz.matrix.transform-scale((2, 0.7, 1)),
    ))
    let named = cetz.draw.line((0, 0), (1, 1), name: "custom-line")
    let callback = (ctx => (ctx: ctx + (custom-calls: ctx.at("custom-calls", default: 0) + 1), drawables: ()))
    let packet = drawing._mark-carrier-input(curve.line((0, 0), (3, 0)), (
      stroke: (paint: purple, thickness: 1pt, cap: "round", join: "bevel", miter-limit: 9),
      mark: curve.mark.bar(width: 0pt, stroke: false),
      mark-position: "center",
    ))
    let registered = drawing._register-mark-elements(transformed,
      (named, callback, packet, drawing._mark-carrier-input(path, (
        name: "marked", stroke: purple + 1pt,
        mark: curve.mark.triangle(length: 8pt, width: 6pt),
      )), cetz.draw.line("marked.end", (4, 0), name: "after-mark")), empty-batch)
    let batch = registered.batch
    let candidates = curve.mark.geometry(batch.templates, batch.carriers, batch.placements)
    assert.eq(candidates.marks.first().shaft-style, (cap: "round", join: "bevel", miter-limit: 9.0))
    assert.eq(batch.templates.first().context.line-thickness, cetz.util.resolve-number(ctx, 1pt))
    let bounds = drawing._engine-mark-bounds(batch.plans.first(), candidates.marks.first())
    // Footprints come from engine stroke outlines, not a line-radius estimate.
    for point in curve.points(candidates.marks.first().footprint) {
      assert(bounds.any(box => point.at(0) >= box.left and point.at(0) <= box.right
        and point.at(1) >= box.bottom and point.at(1) <= box.top))
    }
    let selected = curve.mark.geometry(batch.templates, batch.carriers, batch.placements, mode: "selected")
    let painted = cetz.process.many(transformed, drawing._materialize-mark-groups(
      registered.elements, batch, selected, selected: true,
    ).flatten(), compute-bounds: false)
    assert.eq(painted.ctx.custom-calls, 1)
    assert("custom-line" in painted.ctx.nodes)
    assert("after-mark" in painted.ctx.nodes)
    // Named marked paths retain unshortened start/end/control anchors.
    let original = curve.segments(batch.plans.last().carrier).first()
    for (anchor, expected) in (
      ("start", original.start), ("end", original.end),
      ("ctrl-0", original.control-start), ("ctrl-1", original.control-end),
    ) {
      let actual = (painted.ctx.nodes.at("marked").anchors)(anchor).slice(0, 2)
      assert(drawing._point-distance(actual, expected) < 1e-9)
    }
    assert(drawing._point-distance(curve.points(selected.shafts.last().shaft).last(),
      original.end) > 0.01)
    (ctx: ctx, drawables: painted.drawables)
  },)
})

#cetz.canvas({
  (ctx => {
    let packet = drawing._mark-carrier-input(curve.line((0, 0), (3, 0)), (
      stroke: black + 1pt, mark: curve.mark.triangle(length: 4pt, width: 3pt),
    ))
    let shifted = ctx + (
      transform: cetz.matrix.transform-translate(20, 0, 0),
      resolve-coordinate: ((ctx, point) => (point.at(0) + 0.4, point.at(1) - 0.2),),
    )
    let first = drawing._register-mark-elements(ctx, (packet,), empty-batch)
    let second = drawing._register-mark-elements(shifted, (packet,), first.batch)
    let batch = second.batch
    assert.eq(batch.templates.len(), 1)
    let candidates = curve.mark.geometry(batch.templates, batch.carriers, batch.placements)
    assert(drawing._point-distance(candidates.marks.first().tip, (3, 0)) < 1e-9)
    assert(drawing._point-distance(candidates.marks.last().tip, (23.4, -0.2)) < 1e-9)
    // The resolver affected only the carrier, not the mark's physical size.
    assert.eq(candidates.marks.first().end, candidates.marks.last().end)
    let selected = curve.mark.geometry(batch.templates, batch.carriers, batch.placements, mode: "selected")
    let finished = drawing._materialize-mark-groups(
      first.elements + second.elements, batch, selected, selected: true).flatten()
    let painted = cetz.process.many(ctx, finished, compute-bounds: false)
    (ctx: ctx, drawables: painted.drawables)
  },)
})

#cetz.canvas({
  (ctx => {
    let collapsed = ctx + (transform: cetz.matrix.transform-scale((0, 1, 1)))
    let packet = drawing._mark-carrier-input(curve.line((0, 0), (0, 3)), (
      stroke: black + 1pt, mark: curve.mark.stealth(fit: "bend"),
      pattern: "wave", pattern-amplitude: 0.1, pattern-wavelength: 0.5,
    ))
    let registered = drawing._register-mark-elements(collapsed, (packet,), empty-batch)
    let batch = registered.batch
    let selected = curve.mark.geometry(batch.templates, batch.carriers, batch.placements, mode: "selected")
    let painted = cetz.process.many(collapsed, drawing._materialize-mark-groups(
      registered.elements, batch, selected, selected: true,
    ).flatten(), compute-bounds: false)
    assert.eq(painted.drawables.filter(drawable => cetz.drawable.TAG.hidden not in drawable.tags).len(), 2)
    for drawable in painted.drawables {
      for (origin, _, commands) in drawable.segments {
        assert(origin.all(value => value == value))
        assert(commands.all(command => command.slice(1).flatten().all(value => value == value)))
      }
    }
    (ctx: ctx, drawables: painted.drawables)
  },)
})

// Bulk registration must match independently prepared candidate contexts,
// including callbacks, named-anchor ordering, resolver updates, and paints.
#cetz.canvas({
  (ctx => {
    let red-packet = drawing._mark-carrier-input(path, (
      name: "candidate", stroke: red + 1pt,
      mark: curve.mark.combine(curve.mark.circle(fill: orange), curve.mark.bar(stroke: green)),
    ))
    let blue-packet = red-packet
    blue-packet.style.stroke = blue + 1pt
    let update = (ctx => (ctx: ctx + (
      custom-calls: ctx.at("custom-calls", default: 0) + 1,
    ), drawables: ()))
    let shifted = ctx + (
      length: ctx.length * 2,
      transform: cetz.matrix.transform-translate(20, 0, 0),
      resolve-coordinate: ((ctx, point) =>
        (point.at(0) + ctx.at("custom-calls", default: 0), point.at(1)),),
    )
    let groups = (
      (),
      (red-packet, update, cetz.draw.line("candidate.end", (4, 0))),
      (red-packet, update, cetz.draw.line("candidate.end", (4, 0))),
      (blue-packet, update, cetz.draw.line("candidate.end", (4, 0))),
      (update, red-packet, update, blue-packet,
        cetz.draw.hide(cetz.draw.line((0, 0), (1, 1)))),
    )
    let contexts = (ctx, ctx, ctx, ctx, shifted)
    let bulk = drawing._register-mark-elements(ctx, groups, empty-batch,
      group-contexts: contexts)
    let reference = empty-batch
    let reference-elements = ()
    for (group, independent) in groups.zip(contexts) {
      let registered = drawing._register-mark-elements(independent, (group,), reference)
      reference = registered.batch
      reference-elements += registered.elements
    }
    assert.eq(bulk.batch, reference)
    assert.eq(bulk.batch.carriers.len(), 5)
    assert.eq(bulk.batch.templates.len(), 2)
    assert.eq(bulk.batch.plans.at(0).line-paint, red)
    assert.eq(bulk.batch.plans.at(2).line-paint, blue)
    assert.eq(bulk.batch.plans.at(3).ctx.custom-calls, 1)
    assert.eq(bulk.batch.plans.at(4).ctx.custom-calls, 2)
    assert.eq(bulk.ctx.custom-calls, 2)
    // Exact native request and painted-structure equality, not just footprints.
    for mode in ("candidates", "selected") {
      let geometry = curve.mark.geometry(bulk.batch.templates, bulk.batch.carriers,
        bulk.batch.placements, mode: mode)
      let expected = curve.mark.geometry(reference.templates, reference.carriers,
        reference.placements, mode: mode)
      assert.eq(geometry, expected)
      let finished = drawing._materialize-mark-groups(
        bulk.elements, bulk.batch, geometry, selected: mode == "selected")
      assert.eq(finished.len(), groups.len())
      assert.eq(finished.first(), ())
      let boxes = if mode == "candidates" {
        drawing._candidate-groups-bounds(ctx, bulk.elements, bulk.batch, geometry,
          group-contexts: contexts)
      }
      for (index, independent) in contexts.enumerate() {
        let per-group = drawing._materialize-mark-groups(
          (bulk.elements.at(index),), bulk.batch, geometry, selected: mode == "selected").first()
        let actual = cetz.process.many(independent, finished.at(index), compute-bounds: false)
        let painted = cetz.process.many(independent, per-group, compute-bounds: false)
        assert.eq(actual.drawables, painted.drawables)
        assert.eq(actual.ctx, painted.ctx)
        if mode == "candidates" {
          assert.eq(boxes.at(index), drawing._candidate-groups-bounds(
            independent, (bulk.elements.at(index),), bulk.batch, geometry).first())
        }
      }
      let actual = cetz.process.many(ctx, drawing._materialize-mark-groups(
        bulk.elements, bulk.batch, geometry, selected: mode == "selected").flatten(),
        compute-bounds: false)
      let painted = cetz.process.many(ctx, drawing._materialize-mark-groups(
        reference-elements, reference, expected, selected: mode == "selected").flatten(),
        compute-bounds: false)
      assert.eq(actual.drawables, painted.drawables)
      assert.eq(actual.ctx, painted.ctx)
    }
    // CeTZ's inherited bezier root applies only to single-segment carriers.
    // Sharing the context and mark is not sufficient to share resolved style.
    let inherited = ctx + (style: ctx.style + (bezier: (
      fill: none, stroke: purple + 2pt,
    )))
    let single = drawing._mark-carrier-input(path, (mark: curve.mark.triangle()))
    let multiple = drawing._mark-carrier-input(curve.path(
      ..curve.elements(path), ..curve.elements(curve.line((3, 0), (4, 0))),
    ), single.style)
    let roots = drawing._register-mark-elements(ctx, ((single,), (multiple,), (single,)),
      empty-batch, group-contexts: (inherited, inherited, inherited))
    let root-reference = empty-batch
    for packet in (single, multiple, single) {
      root-reference = drawing._register-mark-elements(inherited, (packet,), root-reference).batch
    }
    assert.eq(roots.batch, root-reference)
    assert.eq(roots.batch.plans.first().line-paint, purple)
    assert.eq(roots.batch.plans.last().line-paint, purple)
    assert(roots.batch.plans.at(1).line-paint != purple)
    (ctx: ctx, drawables: ())
  },)
})

// Exercise the actual drawing owner, including movable label candidates.
#let g = graph.build({
  graph.node(<a>, pos: graph.pos(x: graph.pin(0), y: graph.pin(0)))
  graph.node(<b>, pos: graph.pos(x: graph.pin(4), y: graph.pin(0)))
  graph.edge(graph.source(<a>), graph.sink(<b>), <ab>, pos: graph.pos(x: graph.pin(2), y: graph.pin(1)))
})
#draw(g, node-label: none, node-style: (radius: 0.04), edge-style: (
  stroke: black + 0.6pt,
  mark: curve.mark.stealth(length: 5pt, width: 4pt, fit: "bend"),
  mark-position: "center-if-dangling", mark-orientation: "edge",
  label: [$p$], label-slide: sys.inputs.at("slide", default: "true") == "true",
  label-path: (
    stroke: blue + 0.5pt, offset: 0.3, length: 1.1,
    mark: curve.mark.straight(length: 4pt, width: 3pt),
  ),
))
