from pathlib import Path
p=Path('vendor/kurbo/src/stroke.rs');s=p.read_text();anchor='    // A degenerate stroke with a cusp at the endpoint.';assert anchor in s;s=s.replace(anchor,'''    #[test]
    fn finite_near_half_turn_miters_survive_scaling_and_translation() {
        for scale in [1e-6, 1.0, 1e6] {
            for angle in [1e-4, 1e-8, 1e-10] {
                let origin = 1e6 * scale;
                let mut path = BezPath::new();
                path.move_to((origin - scale, angle * scale));
                path.line_to((origin, 0.0));
                path.line_to((origin - scale, -angle * scale));
                let pen = Stroke::new(scale).with_join(Miter).with_miter_limit(3.0 / angle);
                let paint = stroke(path.elements().iter().copied(), &pen, &StrokeOpts::default(), 1e-6);
                // Intersect the two offset lines analytically for this
                // symmetric V, using the actual represented edge length.
                let dx = origin - (origin - scale);
                let height = angle * scale;
                let expected = origin + 0.5 * scale * dx.hypot(height) / height;
                assert!((paint.bounding_box().x1 - expected).abs() < 1e-12 * expected.abs());
                let clipped = stroke(path, &pen.with_miter_limit(2.0), &StrokeOpts::default(), 1e-6);
                assert!(clipped.bounding_box().x1 <= origin + scale);
            }
        }
        let mut line = BezPath::new();
        line.move_to((-1.0, 0.0));
        line.line_to((0.0, 0.0));
        line.line_to((-1.0, 0.0));
        let paint = stroke(line, &Stroke::new(1.0).with_join(Miter), &StrokeOpts::default(), 1e-6);
        assert!(paint.is_finite());
        assert!(paint.bounding_box().x1 <= 0.5);
    }

'''+anchor,1);p.write_text(s)
