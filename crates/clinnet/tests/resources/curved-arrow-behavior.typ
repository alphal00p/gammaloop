#import "crates/linnest/typst/src/curve.typ" as curve
#import "crates/linnest/typst/src/impl/draw.typ" as drawing
#import "@preview/cetz:0.5.1" as cetz
#import "map-style.typ" as feynman

// Include this file, or import and render this value: canvas assertions need layout.
#let curved-arrow-behavior = {
  let accuracy = 1e-6
  let epsilon = 1e-8
  let source = curve.cubic((0, 0), (0, 2), (1, 3), (3, 3))
  let sink = curve.cubic((3, 3), (5, 3), (6, 2), (6, 0))
  let path = curve.path(source, sink)
  let geometry-checks = ()

  // Length, ratio and shift are measured on the offset path, not its centerline.
  for offset in (-0.35, 0.35) {
    let parallel = curve.parallel(path, distance: offset, accuracy: accuracy)
    let total = curve.length(parallel, accuracy: accuracy)
    for ratio in (none, 0.6) {
      for shift in (-0.2, 0.2) {
        let style = (
          offset: offset,
          length: if ratio == none { 2 } else { none },
          ratio: ratio,
          shift: shift,
          accuracy: accuracy,
        )
        let target = if ratio == none { 2 } else { ratio * total }
        let margin = (total - 0.1 - 0.2 - target) / 2
        let expected = curve.trim(
          parallel,
          start-outset: 0.1 + margin + shift,
          end-outset: 0.2 + margin - shift,
          accuracy: accuracy,
        )
        let case = repr((offset: offset, ratio: ratio, shift: shift))
        geometry-checks.push((
          "offset layer " + case,
          drawing._path-layer(path, style, 0.1, 0.2, none, auto),
          expected,
          target,
        ))
        for mode in ("whole", "gap", "opposite") {
          let gap = if mode == "gap" { 0.2 } else { 0 }
          let source-style = style + (split-gap: gap)
          let sink-style = (
            source-style
              + (
                offset: if mode == "opposite" { -offset } else { offset },
              )
          )
          let paired = drawing._split-edge-geometry(
            source,
            sink,
            path,
            source-style,
            sink-style,
            0.1,
            0.2,
            none,
            accuracy,
          )
          assert.eq(paired.split-gap, gap, message: "paired gap " + case)
          if mode == "whole" {
            assert.ne(paired.whole, none, message: "missing whole path " + case)
            geometry-checks.push((
              "paired whole " + case,
              curve.path(..paired.whole.map(curve.from-cubic)),
              expected,
              target,
            ))
          } else {
            assert.eq(
              paired.whole,
              none,
              message: "unexpected whole path " + case,
            )
          }
          for (index, half-style) in (source-style, sink-style).enumerate() {
            let halves = (source, sink).map(part => curve.parallel(
              part,
              distance: half-style.offset,
              accuracy: accuracy,
            ))
            let lengths = halves.map(part => curve.length(
              part,
              accuracy: accuracy,
            ))
            let target = if ratio == none { 2 } else { ratio * lengths.sum() }
            let margin = (lengths.sum() - 0.1 - 0.2 - target) / 2
            let origin = if index == 0 { 0 } else { lengths.first() }
            let start = calc.max(
              if index == 1 { gap / 2 } else { 0 },
              0.1 + margin + shift - origin,
            )
            let end = calc.min(
              lengths.at(index) - if index == 0 { gap / 2 } else { 0 },
              lengths.sum() - 0.2 - margin + shift - origin,
            )
            assert(
              end > start,
              message: "fixture must retain both paired halves",
            )
            geometry-checks.push((
              "paired " + repr((mode, index)) + " " + case,
              curve.path(
                ..(paired.source, paired.sink).at(index).map(curve.from-cubic),
              ),
              curve.trim(
                halves.at(index),
                start-outset: start,
                end-outset: lengths.at(index) - end,
                accuracy: accuracy,
              ),
              end - start,
            ))
          }
        }
      }
    }
  }
  for (case, actual, expected, target) in geometry-checks {
    assert(
      calc.abs(curve.length(actual, accuracy: accuracy) - target)
        < 8 * accuracy,
      message: case + ": visible arc length",
    )
    let actual = curve.segments(actual)
    let expected = curve.segments(expected)
    for (a, b) in (
      (actual.first().start, expected.first().start),
      (actual.last().end, expected.last().end),
    ) {
      assert(
        cetz.vector.dist(a, b) < 8 * accuracy,
        message: case + ": shifted endpoint",
      )
    }
  }

  cetz.canvas({
    cetz.draw.get-ctx(ctx => {
      let segment = (
        start: (0, 0),
        control-start: (0, 2),
        control-end: (2, -2),
        end: (2, 0),
      )
      let curved = curve.from-cubic(segment)
      let straight = curve.from-cubic(curve.line-segment((0, 0), (2, 0)))
      let head = (
        anchor: "center",
        fill: black,
        stroke: black + 0.8pt,
        length: 0.6,
        width: 0.3,
        inset: 0,
        shorten-to: auto,
      )
      let identity = cetz.matrix.ident(4)
      let cases = ()
      let comparisons = ()

      // Native physical sizes, reversal and named anchors survive transforms;
      // boundary heads now move inward instead of extending beyond the carrier.
      for transform in (
        identity,
        cetz.matrix.mul-mat(
          cetz.matrix.transform-translate(3, 2, 0),
          cetz.matrix.transform-rotate-z(37deg),
        ),
        cetz.matrix.transform-scale((-1.5, 0.6, 1)),
        cetz.matrix.transform-shear-x(0.6),
        cetz.matrix.transform-rotate-xyz(40deg, 30deg, 0deg),
      ) {
        for shape in (false, true) {
          for path in (straight, curved) {
            cases.push((path: path, transform: transform, shape: shape))
          }
        }
      }
      for entry in (
        (scale: 1.5),
        (length: 700%, width: 400%, stroke: black + 0.5pt),
      ) {
        cases.push((path: curved, ratio: 0.5, head: entry))
      }
      // Interior heads overlay the entire shaft, even at repeated controls.
      // Exercise the full path, with a head spanning a cubic join, not a terminal half.
      for path in (
        curve.from-cubic(segment + (control-start: segment.start)),
        curve.from-cubic(segment + (control-end: segment.end)),
        curve.path(curved, curve.from-cubic(curve.line-segment(
          (2, 0),
          (2, 0.02),
        ))),
        curve.path(
          curve.cubic((0, 0), (0, 0.5), (0.5, 1), (1, 1)),
          curve.cubic((1, 1), (1.5, 1), (2, 0.5), (2, 0)),
        ),
      ) {
        for ratio in (none, 0, 0.2, 0.5, 0.8, 1) {
          for shift in (-0.1, 0.1) {
            cases.push((path: path, ratio: ratio, shift: shift))
          }
        }
      }
      // A fully consumed shaft still has a compressed mark; only an empty path has none.
      for length in (0, 1e-9, 1e-6, 0.01, 0.1) {
        for ratio in (none, 0.5) {
          cases.push((
            path: curve.from-cubic(curve.line-segment((0, 0), (length, 0))),
            ratio: ratio,
          ))
        }
      }

      // End heads and centered heads have two actual geometric contacts, not tangent alignment.
      for (index, config) in cases.enumerate() {
        let path = config.path
        let head = head + config.at("head", default: (:))
        let transform = config.at("transform", default: identity)
        let shape = config.at("shape", default: false)
        let ratio = config.at("ratio", default: none)
        let shift = config.at("shift", default: 0)
        let local = ctx + (transform: transform)
        let space = ctx + (transform: if shape { identity } else { transform })
        let carrier = curve
          .to-cetz(path, mark: none)
          .first()(space)
          .drawables
          .first()
        if not shape {
          carrier = cetz
            .drawable
            .apply-transform(cetz.matrix.transform-scale((1, 1, 0)), carrier)
            .first()
        }
        let pieces = curve.segments(path)
        let native = if pieces.len() == 1 {
          let piece = pieces.first()
          cetz
            .draw
            .bezier(
              piece.start,
              piece.end,
              piece.control-start,
              piece.control-end,
              mark: none,
            )
            .first()(local)
        } else { curve.to-cetz(path, mark: none).first()(local) }
        let total = cetz.path-util.length(carrier.segments)
        for root in ("start", "end") {
          for symbol in (">", "<", "straight") {
            for reverse in (false, true) {
              for anchor in if ratio == none {
                ("tip", "center", "base")
              } else { ("center",) } {
                let entry = (
                  head + (symbol: symbol, reverse: reverse, anchor: anchor)
                )
                let mark = (transform-shape: shape)
                mark.insert(root, entry)
                let style = (
                  name: "named",
                  stroke: black + 0.8pt,
                  mark: mark,
                  mark-position: if ratio == none { "end" } else { ratio },
                  mark-direction: if root == "start" { "backward" } else {
                    "forward"
                  },
                  mark-shift: shift,
                  accuracy: accuracy,
                )
                let processed = cetz
                  .mark
                  .process-style(
                    ctx,
                    cetz
                      .styles
                      .resolve(
                        ctx.style,
                        merge: drawing._draw-style(style),
                        root: "bezier",
                      )
                      .mark,
                    root,
                    total,
                  )
                  .first()
                let result = cetz.process.many(
                  local,
                  drawing
                    ._derived-path-elements(path, style, auto, true, true)
                    .elements
                    .flatten(),
                  compute-bounds: false,
                )
                let actual = result.drawables

                let marks = actual.filter(item => (
                  cetz.drawable.TAG.mark in item.tags
                ))
                let case = repr((index, root, symbol, reverse, anchor))
                assert.eq(
                  marks.len(),
                  if total == 0 { 0 } else { 1 },
                  message: case + ": head count",
                )
                if total == 0 { continue }
                let (origin, _, commands) = marks.first().segments.first()
                let tip = if symbol == "straight" {
                  commands.first().last()
                } else { origin }
                let wings = if symbol == "straight" {
                  (origin, commands.at(1).last())
                } else { (commands.at(0).last(), commands.at(1).last()) }
                let back = cetz.vector.lerp(..wings, 0.5)
                let reversed = reverse != (symbol == "<")
                let span = calc.min(processed.length, total)
                let tip-at = (tip: 0, center: -span / 2, base: -span).at(anchor)
                let back-at = tip-at + span
                if reversed { (tip-at, back-at) = (back-at, tip-at) }
                let at = if ratio == none { 0 } else {
                  if root == "start" { ratio * total + shift } else {
                    (1 - ratio) * total - shift
                  }
                }
                let low = calc.min(tip-at, back-at) + at
                let inward = calc.clamp(low, 0, total - span) - low
                let stations = (tip-at, back-at).map(s => s + at + inward)
                let contacts = stations.map(s => {
                  cetz
                    .path-util
                    .point-at(
                      carrier.segments,
                      s,
                      reverse: root == "end",
                    )
                    .point
                })
                let axis = cetz.vector.norm(cetz.vector.sub(
                  contacts.last(),
                  contacts.first(),
                ))
                let width = processed.width
                if shape {
                  let normal = (-axis.at(1), axis.at(0), 0)
                  width *= cetz.vector.dist(
                    cetz.matrix.mul4x4-vec3(transform, normal),
                    cetz.matrix.mul4x4-vec3(transform, (0, 0, 0)),
                  )
                  contacts = contacts.map(p => cetz.matrix.mul4x4-vec3(
                    transform,
                    p,
                  ))
                }
                for (actual, expected) in (tip, back).zip(contacts) {
                  assert(
                    cetz.vector.dist(actual, expected) < epsilon,
                    message: case
                      + ": actual tip/back contact "
                      + repr((actual, expected)),
                  )
                }
                assert(
                  calc.abs(cetz.vector.dist(..wings) - width) < epsilon,
                  message: case + ": width must not compress",
                )
                assert.eq(
                  marks.first().stroke.thickness,
                  processed.stroke.thickness,
                  message: case + ": stroke size",
                )
                let shaft = actual.filter(item => (
                  cetz.drawable.TAG.mark not in item.tags
                ))
                if ratio != none {
                  comparisons.push((
                    case + ": continuous interior shaft",
                    shaft.map(d => d.segments),
                    cetz
                      .drawable
                      .apply-transform(
                        if shape { transform } else { identity },
                        carrier,
                      )
                      .map(d => d.segments),
                  ))
                } else if total > processed.length {
                  let contact = if symbol == "straight" or reversed {
                    tip
                  } else { back }
                  let endpoint = if root == "start" {
                    cetz.path-util.first-subpath-start(shaft.first().segments)
                  } else {
                    cetz.path-util.last-subpath-end(shaft.last().segments)
                  }
                  assert(
                    cetz.vector.dist(endpoint, contact) < epsilon,
                    message: case + ": painted shaft contact",
                  )
                }
                // Named anchors continue to describe the unshortened carrier.
                for anchor in (
                  ("start", "end")
                    + if pieces.len() == 1 { ("ctrl-0", "ctrl-1") } else { () }
                ) {
                  let point = (result.ctx.nodes.at("named").anchors)(anchor)
                  let expected = (native.anchors)(anchor)
                  assert(
                    cetz.vector.dist(point, expected) < epsilon,
                    message: case + ": named anchor",
                  )
                }
              }
            }
          }
        }
      }

      // Both ends can shorten the same cubic, or consume the entire short shaft.
      for path in (
        curved,
        curve.from-cubic(curve.line-segment((0, 0), (0.1, 0))),
      ) {
        let style = (mark: (symbol: ">", ..head))
        let actual = cetz
          .process
          .many(
            ctx,
            drawing
              ._derived-path-elements(path, style, auto, true, true)
              .elements
              .flatten(),
            compute-bounds: false,
          )
          .drawables
        assert.eq(actual.len(), 3, message: "two heads retained")
        if actual.first().segments.len() > 0 {
          for (mark, endpoint) in actual
            .slice(1)
            .zip((
              cetz.path-util.first-subpath-start(actual.first().segments),
              cetz.path-util.last-subpath-end(actual.first().segments),
            )) {
            let commands = mark.segments.first().last()
            assert(
              cetz.vector.dist(endpoint, cetz.vector.lerp(
                commands.at(0).last(),
                commands.at(1).last(),
                0.5,
              ))
                < epsilon,
              message: "two-ended shaft contacts",
            )
          }
        } else {
          assert(
            curve.length(path, accuracy: accuracy) < 2 * head.length,
            message: "only overlapping heads consume the shaft",
          )
        }
      }

      // Mnemonics and custom overrides resolve before built-in contact handling.
      let custom = ctx
      custom.marks.marks.insert("straight", cetz.mark-shapes.marks.diamond)
      custom.marks.mnemonics.insert("arrow-alias", "<")
      for (local, symbol, target) in (
        (custom, "straight", "diamond"),
        (custom, "arrow-alias", "<"),
      ) {
        for reverse in (false, true) {
          let variants = (symbol, target).map(symbol => {
            let style = (
              mark: (
                end: head
                  + (symbol: symbol, reverse: reverse, flip: true, slant: 30%),
              ),
              mark-position: "center",
            )
            cetz
              .process
              .many(
                local,
                drawing
                  ._derived-path-elements(curved, style, auto, true, true)
                  .elements
                  .flatten(),
                compute-bounds: false,
              )
              .drawables
          })
          assert.eq(..variants, message: "resolved symbol " + symbol)
        }
      }

      // Patterns and crossing gaps only change paint, never the single full-path mark.
      for position in ("end", "center") {
        for symbol in (">", "straight") {
          let style = (
            mark: (end: head + (symbol: symbol)),
            mark-position: position,
          )
          let reference = cetz
            .process
            .many(
              ctx,
              drawing
                ._derived-path-elements(curved, style, auto, true, true)
                .elements
                .flatten(),
              compute-bounds: false,
            )
            .drawables
            .filter(d => cetz.drawable.TAG.mark in d.tags)
          for pattern in (none, "wave", "coil") {
            let paint = style + (pattern: pattern, crossing-gap: 0.4)
            let total = curve.length(curved, accuracy: accuracy)
            for elements in (
              drawing
                ._derived-path-elements(curved, paint, auto, true, true)
                .elements,
              drawing._cut-path-elements(
                curved,
                paint,
                (total / 2,),
                mark-style: style,
              ),
            ) {
              let marks = cetz
                .process
                .many(ctx, elements.flatten(), compute-bounds: false)
                .drawables
                .filter(d => cetz.drawable.TAG.mark in d.tags)
              assert.eq(
                marks,
                reference,
                message: repr((position, symbol, pattern))
                  + ": undistorted carrier",
              )
            }
          }
        }
      }

      // Coordinate resolvers apply to the curve, never to the mark's private anchors.
      let local = cetz
        .draw
        .register-coordinate-resolver((ctx, point) => {
          if type(point) == array { point.map(value => value * 2) } else {
            point
          }
        })
        .first()(ctx)
        .ctx
      let style = (mark: (end: head + (symbol: ">")), mark-position: "center")
      let resolved = cetz
        .process
        .many(
          local,
          drawing
            ._derived-path-elements(straight, style, auto, true, true)
            .elements
            .flatten(),
          compute-bounds: false,
        )
        .drawables
      let doubled = curve.from-cubic(curve.line-segment((0, 0), (4, 0)))
      let expected = cetz
        .process
        .many(
          ctx,
          drawing
            ._derived-path-elements(doubled, style, auto, true, true)
            .elements
            .flatten(),
          compute-bounds: false,
        )
        .drawables
      comparisons.push((
        "coordinate resolver",
        resolved.map(d => d.segments),
        expected.map(d => d.segments),
      ))

      // Endpoint clamps remain exact even below the arc-length accuracy.
      let label-accuracy = drawing._style-value((:), "accuracy")
      for length in (0, label-accuracy / 2) {
        let segment = curve.line-segment((0, 0), (length, 0))
        let path = curve.from-cubic(segment)
        for (shift, t) in ((-1, 0), (1, 1)) {
          let frame = drawing._path-mid-frame(path, label-accuracy, shift: shift)
          assert.eq(frame.point, curve.cubic-point(segment, t))
          assert.eq(frame.tangent, curve.cubic-tangent(segment, t))
        }
      }

      // Measure the rendered frame for automatic clearance and the actual CeTZ
      // anchor for explicit placement, never plain Typst text metrics.
      let label-checks = ()
      let halves = (
        source: curve.segments(source),
        sink: curve.segments(sink),
        whole: curve.segments(path),
      )
      for (path-index, (path, paired)) in (
        (straight, none),
        (curved, none),
        (path, (halves: halves, source: true, sink: true)),
        (path, (halves: halves + (whole: none), source: true, sink: true)),
        (source, (halves: halves, source: true, sink: false)),
        (sink, (halves: halves, source: false, sink: true)),
      ).enumerate() {
        for side in ("left", "right") {
          for (label-index, label) in (
            [$k-p_2$],
            [Hgyp],
            [Hg\ yp],
            box(width: 6pt, height: 6pt),
            box(width: 36pt, height: 6pt),
          ).enumerate() {
            for (style-index, label-style) in (
              (:),
              (anchor: auto),
              (anchor: "auto"),
              (anchor: "center"),
              (anchor: "east"),
              (
                wrap: body => box(width: 18pt, text(size: 9pt, body)),
                padding: (left: 0.15, right: 0.25, top: 0.1, bottom: 0.2),
              ),
              (anchor: "east", angle: 13deg, padding: (left: 0.15, top: 0.2)),
              (angle: -23deg, auto-scale: true),
              (angle: 37deg, anchor: "north-east"),
              (angle: -23deg, anchor: "south-west", auto-scale: true),
            ).enumerate() {
              label-checks.push((
                repr((path-index, side, label-index, style-index)),
                path,
                (
                  label: text(size: 6pt, label),
                  label-style: label-style,
                  label-side: side,
                  label-gap: 0.2,
                ),
                paired,
              ))
            }
          }
        }
      }

      // Fixed xbox-opened2 sink geometry: label shifts are absolute on the full
      // offset path. Equal and independent shifts must not move the arrow, and
      // arrow length must not change label placement or its endpoint clamps.
      let endpoints = (
        (3.7141557137182426, -3.510279312439297),
        (1.5082746017729896, -2.213011602860704),
      )
      let native = ctx + (length: 2.6mm)
      for (side, shift, kind) in (
        ("right", -1, "incoming"), ("right", -0.4, "incoming"),
        ("right", 1, "incoming"), ("left", -0.4, "incoming"),
        ("right", -0.4, "outgoing"), ("right", 0.4, "paired"),
        ("left", -0.4, "paired"),
      ) {
        let path = if kind == "paired" { path } else {
          drawing._dangling-path(
            endpoints.at(if kind == "incoming" { 0 } else { 1 }),
            endpoints.at(if kind == "incoming" { 1 } else { 0 }),
            0,
            (dangling-tangent: "horizontal"),
            dangling-at-start: kind == "incoming",
          )
        }
        let fields = (
          momentum-arrow-side: side,
          momentum-arrow-offset: 0.4,
          momentum-arrow-length: 0.7,
          momentum-arrow-shift: shift,
        )
        let reference = feynman.edge-style((momentum: [], fields: fields)).at(1)
        let reference-path = drawing._path-layer(path, reference, 0, 0, none, auto)
        let reference-paint = cetz
          .process
          .many(
            native,
            drawing
              ._derived-path-elements(
                reference-path,
                reference,
                auto,
                true,
                true,
              )
              .elements
              .flatten(),
            compute-bounds: false,
          )
          .drawables
        let full-style = reference + (
          length: none, ratio: none, resolve-length: "none", shift: 0,
        )
        let full-path = drawing._path-layer(path, full-style, 0, 0, none, auto)
        let full-length = curve.length(
          full-path, accuracy: drawing._style-value(full-style, "accuracy"),
        )
        let paired = if kind == "paired" {
          (
            halves: drawing._split-edge-geometry(
              source, sink, path, full-style, full-style, 0, 0, none,
              drawing._style-value(full-style, "accuracy"),
            ),
            source: true,
            sink: true,
          )
        } else { none }
        for (gap, label, anchor) in (
          (0.2, [$k-p_2$], auto),
          (0.6, [Hg\ yp], "auto"),
          (0, box(width: 6pt, height: 6pt), auto),
          (-0.2, [$k-p_2$], auto),
          (0.2, box(width: 36pt, height: 6pt), "center"),
          (0.2, box(width: 36pt, height: 6pt), "east"),
          (0.2, box(width: 36pt, height: 6pt), "north-east"),
          (0.2, box(width: 36pt, height: 6pt), "south-west"),
          (-0.2, box(width: 36pt, height: 6pt), "center"),
        ) {
          let momentum = text(size: 6pt, label)
          let centers = ()
          for requested in (
            none, shift, shift - 1e-6, shift + 1e-6, -0.1, 0,
            -full-length, full-length,
          ) {
            let label-shift = if requested == none { shift } else { requested }
            let case = "momentum " + repr((side, shift, requested, gap, anchor, kind))
            let fields = fields + (
              momentum-label-gap: gap,
              momentum-label-anchor: anchor,
            ) + if requested == none { (:) } else {
              (momentum-label-shift: requested)
            }
            let layers = feynman.edge-style((momentum: momentum, fields: fields))
            let arrow = layers.at(1)
            let carrier = layers.last()
            assert.eq(layers.len(), 3, message: case + ": carrier count")
            assert.eq(
              (carrier.length, carrier.ratio, carrier.resolve-length, carrier.shift),
              (none, none, "none", 0),
              message: case + ": full unshifted label carrier",
            )
            assert.eq(
              carrier.label-shift,
              label-shift,
              message: case + ": requested absolute label shift",
            )
            assert.eq(carrier.stroke, none, message: case + ": invisible carrier")
            assert.eq(carrier.mark, none, message: case + ": no duplicate head")
            assert.eq(
              arrow.at("label", default: none), none,
              message: case + ": no label attached to arrow",
            )
            let arrow-path = drawing._path-layer(path, arrow, 0, 0, none, auto)
            let label-path = drawing._path-layer(path, carrier, 0, 0, none, auto)
            assert.eq(
              arrow-path,
              reference-path,
              message: case + ": unchanged arrow path",
            )
            assert.eq(
              cetz
                .process
                .many(
                  native,
                  drawing
                    ._derived-path-elements(arrow-path, arrow, auto, true, true)
                    .elements
                    .flatten(),
                  compute-bounds: false,
                )
                .drawables,
              reference-paint,
              message: case + ": unchanged shaft and chord head",
            )
            assert.eq(label-path, full-path, message: case + ": full offset path")
            let frame = drawing._path-mid-frame(
              label-path, drawing._style-value(carrier, "accuracy"),
              shift: carrier.label-shift,
            )
            let rendered = cetz
              .process
              .many(
                native,
                drawing
                  ._layer-label-element(native, label-path, carrier, none, (:))
                  .flatten(),
                compute-bounds: false,
              )
              .drawables
            let labels = rendered.filter(d => d.type == "content")
            assert.eq(labels.len(), 1, message: case + ": one momentum label")
            centers.push(labels.first().pos)
            label-checks.push((case, label-path, carrier, paired))
            for length in (0.2, 0.8 * full-length, 2 * full-length) {
              let case = case + " length " + repr(length)
              let layers = feynman.edge-style((
                momentum: momentum,
                fields: fields + (momentum-arrow-length: length),
              ))
              assert.eq(
                layers.at(1), arrow + (length: length),
                message: case + ": only arrow length changes",
              )
              let carrier = layers.last()
              let label-path = drawing._path-layer(path, carrier, 0, 0, none, auto)
              assert.eq(
                label-path, full-path,
                message: case + ": length-independent path",
              )
              let shifted = drawing._path-mid-frame(
                label-path, drawing._style-value(carrier, "accuracy"),
                shift: carrier.label-shift,
              )
              assert.eq(
                shifted.point, frame.point,
                message: case + ": unchanged label point",
              )
              assert.eq(
                cetz.vector.norm(shifted.tangent), cetz.vector.norm(frame.tangent),
                message: case + ": unchanged label normal",
              )
              assert.eq(
                cetz.process.many(
                  native,
                  drawing
                    ._layer-label-element(native, label-path, carrier, none, (:))
                    .flatten(),
                  compute-bounds: false,
                ).drawables,
                rendered,
                message: case + ": unchanged final label placement",
              )
            }
          }
          for center in centers.slice(1, 4) {
            assert(
              cetz.vector.dist(centers.first(), center) < 8 * accuracy,
              message: "equal/adjacent label continuity "
                + repr((side, shift, gap, anchor, kind, centers)),
            )
          }
        }
      }
      for (transform-index, transform) in (
        identity,
        cetz.matrix.mul-mat(
          cetz.matrix.transform-translate(3, 2, 0),
          cetz.matrix.transform-rotate-z(37deg),
        ),
        cetz.matrix.transform-scale((-1.5, 0.6, 1)),
        cetz.matrix.transform-shear-x(0.6),
        cetz.matrix.transform-rotate-xyz(40deg, 30deg, 0deg),
      ).enumerate() {
        let local = native + (transform: transform)
        for (case, path, style, paired) in label-checks {
          let case = "label placement " + repr(transform-index) + " " + case
          let label-style = style.label-style + (name: "checked-label")
          let style = style + (label-style: label-style)
          let anchor = label-style.at("anchor", default: auto)
          let automatic = anchor in (auto, "auto")
          let gap = calc.max(0, style.label-gap)
          let accuracy = drawing._style-value(style, "accuracy")
          let label-shift = style.at("label-shift", default: 0)
          let frame = drawing._path-mid-frame(path, accuracy, shift: label-shift)
          let total = curve.length(path, accuracy: accuracy)
          let at = calc.clamp(total / 2 + label-shift, 0, total)
          let (segment, t) = if at == 0 {
            (curve.segments(path).first(), 0)
          } else if at == total {
            (curve.segments(path).last(), 1)
          } else {
            (curve.segments(curve.trim(
              path, end-outset: total - at, accuracy: accuracy,
            )).last(), 1)
          }
          let tolerance = if at in (0, total) { epsilon } else { 8 * accuracy }
          assert(
            cetz.vector.dist(frame.point, curve.cubic-point(segment, t)) < tolerance,
            message: case + ": shifted arc-length point and full-path endpoint clamp",
          )
          let tangent = cetz.vector.norm(frame.tangent)
          assert(
            cetz.vector.dist(
              tangent, cetz.vector.norm(curve.cubic-tangent(segment, t)),
            ) < tolerance,
            message: case + ": local tangent at shifted label point",
          )
          let normal = cetz.vector.scale(
            (-tangent.at(1), tangent.at(0)),
            if style.label-side == "right" { -1 } else { 1 },
          )
          let origin = cetz.matrix.mul4x4-vec3(transform, (..frame.point, 0))
          let outward = cetz.vector.sub(
            cetz.matrix.mul4x4-vec3(transform, (
              ..cetz.vector.add(frame.point, normal),
              0,
            )),
            origin,
          )
          let element = if paired == none {
            drawing._layer-label-element(local, path, style, none, (:))
          } else {
            drawing._paired-layer-label-element(
              local,
              paired.halves,
              if paired.source { style } else { none },
              if paired.sink { style } else { none },
              none,
              (:),
            )
          }
          let rendered = cetz.process.many(
            local, element.flatten(), compute-bounds: false,
          )
          let point = (rendered.ctx.nodes.at("checked-label").anchors)(
            if automatic { "center" } else { anchor },
          )
          let distance = if automatic {
            (
              cetz.vector.dot(cetz.vector.sub(point, origin), outward)
                / cetz.vector.dot(outward, outward)
            )
          } else { gap }
          assert(
            cetz.vector.dist(
              point,
              cetz.vector.add(origin, cetz.vector.scale(outward, distance)),
            ) < epsilon,
            message: case + if automatic { ": centered legacy position" } else {
              ": explicit anchor at shifted arc-length point plus signed normal gap"
            },
          )
          let frames = rendered.drawables.filter(d => "content-frame" in d.tags)
          assert.eq(frames.len(), 1, message: case + ": one rendered box")
          if automatic {
            let (start, _, commands) = frames.first().segments.first()
            let nearest = calc.min(..(start, ..commands.map(c => c.last())).map(
              point => (
                cetz.vector.dot(cetz.vector.sub(point, origin), outward)
                  / cetz.vector.dot(outward, outward)
              ),
            ))
            assert(
              calc.abs(nearest - gap) < epsilon,
              message: case + ": rendered clearance " + repr((nearest, gap)),
            )
          }
        }
      }

      // Compare drawable structure and points in CeTZ's resolved 3D coordinates.
      for (case, actual, expected) in comparisons {
        assert.eq(
          actual,
          expected,
          message: case + ": unchanged shaft geometry",
        )
      }
      ()
    })
  })
}

#curved-arrow-behavior
