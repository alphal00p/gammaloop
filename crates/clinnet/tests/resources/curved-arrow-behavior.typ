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

      // Measure the rendered frame, not plain Typst text metrics or label anchors.
      let label-checks = ()
      for (path-index, path) in (straight, curved).enumerate() {
        for side in ("left", "right") {
          for (label-index, label) in (
            [$k-p_2$],
            [Hgyp],
            [Hg\ yp],
            box(width: 6pt, height: 6pt),
          ).enumerate() {
            for (style-index, label-style) in (
              (:),
              (
                wrap: body => box(width: 18pt, text(size: 9pt, body)),
                padding: (left: 0.15, right: 0.25, top: 0.1, bottom: 0.2),
              ),
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
              ))
            }
          }
        }
      }

      // Fixed xbox-opened2 sink geometry: label shifts are absolute, and switching
      // carrier modes must not move the arrow or introduce different endpoint clamps.
      let incoming = drawing._dangling-path(
        (3.7141557137182426, -3.510279312439297),
        (1.5082746017729896, -2.213011602860704),
        0,
        (dangling-tangent: "horizontal"),
        dangling-at-start: true,
      )
      let native = ctx + (length: 2.6mm)
      for shift in (-1, -0.4, 1) {
        let fields = (
          momentum-arrow-side: "right",
          momentum-arrow-offset: 0.4,
          momentum-arrow-length: 0.7,
          momentum-arrow-shift: shift,
        )
        let reference = feynman.edge-style((momentum: [], fields: fields)).at(1)
        let reference-path = drawing._path-layer(
          incoming,
          reference,
          0,
          0,
          none,
          auto,
        )
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
        for (gap, label) in (
          (0.2, [$k-p_2$]),
          (0.6, [Hg\ yp]),
          (0, box(width: 6pt, height: 6pt)),
        ) {
          let centers = ()
          for label-shift in (shift, shift - 1e-6, shift + 1e-6, -0.1) {
            let case = "momentum " + repr((shift, label-shift, gap))
            let layers = feynman.edge-style((
              momentum: text(size: 6pt, label),
              fields: fields
                + (momentum-label-shift: label-shift, momentum-label-gap: gap),
            ))
            let arrow = layers.at(1)
            let carrier = layers.last()
            assert.eq(
              layers.len(),
              if label-shift == shift { 2 } else { 3 },
              message: case + ": carrier count",
            )
            assert.eq(
              carrier.length,
              arrow.length,
              message: case + ": shared clamp length",
            )
            assert.eq(
              carrier.shift,
              label-shift,
              message: case + ": absolute label shift",
            )
            if label-shift != shift {
              assert.eq(
                carrier.stroke,
                none,
                message: case + ": invisible carrier",
              )
              assert.eq(
                carrier.mark,
                none,
                message: case + ": no duplicate head",
              )
            }
            let arrow-path = drawing._path-layer(
              incoming,
              arrow,
              0,
              0,
              none,
              auto,
            )
            let label-path = drawing._path-layer(
              incoming,
              carrier,
              0,
              0,
              none,
              auto,
            )
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
            assert.eq(
              label-path,
              drawing._path-layer(
                incoming,
                arrow + (shift: label-shift),
                0,
                0,
                none,
                auto,
              ),
              message: case + ": independently shifted copy",
            )
            if calc.abs(shift) == 1 and calc.abs(label-shift - shift) < 2e-6 {
              assert.eq(
                label-path,
                arrow-path,
                message: case + ": identical endpoint clamp",
              )
            }
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
              .filter(d => d.type == "content")
            assert.eq(rendered.len(), 1, message: case + ": one momentum label")
            centers.push(rendered.first().pos)
            label-checks.push((case, label-path, carrier))
          }
          for center in centers.slice(1, 3) {
            assert(
              cetz.vector.dist(centers.first(), center) < 8 * accuracy,
              message: "equal/adjacent label continuity "
                + repr((shift, gap, centers)),
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
        for (case, path, style) in label-checks {
          let case = "label box " + repr(transform-index) + " " + case
          let frame = drawing._path-mid-frame(path, drawing._style-value(
            style,
            "accuracy",
          ))
          let tangent = cetz.vector.norm(frame.tangent)
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
          let rendered = cetz
            .process
            .many(
              local,
              drawing
                ._layer-label-element(local, path, style, none, (:))
                .flatten(),
              compute-bounds: false,
            )
            .drawables
          let frames = rendered.filter(d => "content-frame" in d.tags)
          assert.eq(frames.len(), 1, message: case + ": one rendered box")
          let (start, _, commands) = frames.first().segments.first()
          let nearest = calc.min(..(start, ..commands.map(c => c.last())).map(
            point => (
              cetz.vector.dot(cetz.vector.sub(point, origin), outward)
                / cetz.vector.dot(outward, outward)
            ),
          ))
          assert(
            calc.abs(nearest - style.label-gap) < epsilon,
            message: case
              + ": rendered clearance "
              + repr((nearest, style.label-gap)),
          )
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
