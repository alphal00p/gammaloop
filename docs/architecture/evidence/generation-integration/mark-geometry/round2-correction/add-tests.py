from pathlib import Path
p=Path('crates/kurvst/src/marks.rs');s=p.read_text();anchor='    #[test]\n    fn actual_zero_area_polygons_keep_the_as_given_exception()'
assert anchor in s
new='''    #[test]
    fn painted_miters_preserve_complete_grown_polygon_bounds() {
        // Independent lobe inradius: area / semiperimeter. No offset-normal
        // or stroker calculation is used to predict these outer vertices.
        for scale in [1e-6, 1.0, 1e6] {
            for (shape, width, inset) in [
                ("triangle", 1e-5, 0.0),
                ("triangle", 1e-7, 0.0),
                ("triangle", 1e-9, 0.0),
                ("stealth", 8.0, 0.999999),
                ("stealth", 8.0, 0.99999999),
            ] {
                let length: f64 = 10.0;
                let half = width / 2.0;
                let radius = half * length * (1.0 - inset)
                    / (length.hypot(half) + (length * inset).hypot(half));
                let grow = (1.01 * 0.5 / radius).max(1.0);
                let expected = [-grow * length * scale, 0.0, grow * half * scale];
                let mut spec = sized(shape, length * scale, width * scale);
                spec.inset = (shape == "stealth").then_some(inset);
                spec.fill = Some(false);
                let prepared = spec.prepare(MarkContext {
                    units_per_pt: 1.0,
                    line_thickness: scale,
                    shaft_stroke: MarkStrokeStyle::default(),
                }).unwrap();
                let paint = &prepared.paths[0].outline.path;
                let bounds = paint.bounding_box();
                let label = format!("{shape} width={width} inset={inset} scale={scale}");
                // Subtracting the grown tip's coordinates has absolute error
                // proportional to its extent. The wing expectation additionally
                // allows 1e-9 relative conditioning error in the emitted edges.
                let tip_error = 128.0 * f64::EPSILON * expected[0].abs();
                assert!(bounds.x1.abs() <= tip_error, "{label}: {bounds:?}");
                for (actual, wanted) in [(bounds.x0, expected[0]), (bounds.y1, expected[2]), (bounds.y0, -expected[2])] {
                    assert!((actual - wanted).abs() <= 1e-9 * wanted.abs() + tip_error,
                        "{label}: {actual} versus {wanted}");
                }
                // A real interior point near the grown wing must be painted;
                // matching a metadata rectangle alone cannot satisfy this.
                let near_wing = Point::new(expected[0] * 0.999999, expected[2] * 0.999998);
                assert_ne!(paint.winding(near_wing), 0, "{label}: missing wing");
                assert!(prepared.end().is_finite() && prepared.end() > 0.0);
            }
        }
    }

'''
s=s.replace(anchor,new+anchor);p.write_text(s)
p=Path('crates/kurvst/typst/tests/mark-geometry.typ');s=p.read_text();s+='''
// Complete paint extents, independently predicted from each dart lobe's
// area/semiperimeter inradius. Exercise the public WASM owner, not metadata.
#for scale in (0.000001, 1.0, 1000000.0) {
  for (shape, width, inset) in (("triangle", 1e-7, 0.0), ("triangle", 1e-9, 0.0), ("stealth", 8.0, float("0.99999999"))) {
    let length = 10.0
    let half = width / 2
    let radius = half * length * (1 - inset) / (calc.sqrt(length * length + half * half) + calc.sqrt(length * length * inset * inset + half * half))
    let grow = calc.max(1.0, 1.01 * 0.5 / radius)
    let spec = if shape == "triangle" { mark.triangle(length: length * scale * 1pt, width: width * scale * 1pt, fill: false, stroke: true) } else { mark.stealth(length: length * scale * 1pt, width: width * scale * 1pt, inset: inset * 100%, fill: false, stroke: true) }
    let head = mark.geometry(((mark: spec, "context": (units-per-pt: 1, line-thickness: scale)),), (line,), ((template: 0, carrier: 0, station: (kind: "end")),)).marks.first()
    let boxes = head.paths.first().outline-bounds
    let extent = grow * length * scale
    let tip-error = 128 * float("2.220446049250313e-16") * extent + 1e-12
    assert(calc.abs(calc.max(..boxes.map(box => box.right)) - 100) <= tip-error)
    assert(calc.abs(calc.min(..boxes.map(box => box.left)) - (100 - extent)) <= 1e-9 * extent + tip-error)
    assert(calc.abs(calc.max(..boxes.map(box => box.top)) - grow * half * scale) <= 1e-9 * grow * half * scale + tip-error)
  }
}
''';p.write_text(s)
