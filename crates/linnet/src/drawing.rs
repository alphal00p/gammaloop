//! # Graph Drawing and Visualization
//!
//! This module is intended to provide utilities for visualizing graph structures,
//! primarily by generating Scalable Vector Graphics (SVG) outputs. It aims to
//! facilitate the visual inspection and analysis of graphs represented within
//! the `linnet` library.
//!
//! Currently, the module's direct public API for drawing `HedgeGraph` instances
//! is not yet fully exposed or implemented. The existing code includes examples
//! and tests demonstrating the use of the [`kurbo`](https://crates.io/crates/kurbo)
//! crate for creating and manipulating 2D Bézier paths, and serializing these
//! paths into SVG format.
//!
//! Future development in this module will likely focus on:
//! - Providing functions to take a `HedgeGraph` and produce an SVG representation.
//! - Integrating layout algorithms to position nodes and edges aesthetically.
//! - Offering customization options for styling the output (colors, stroke widths, etc.).
//!
//! Users interested in path drawing can look at the internal tests for examples of
//! using `kurbo::BezPath` and its SVG export capabilities. For direct `HedgeGraph`
//! visualization, the DOT output functionalities (see `dot_parser` module and
//! `HedgeGraph::dot_...` methods) can be used in conjunction with external
//! Graphviz tools.

#[cfg(test)]
pub mod test {

    use kurbo::{
        simplify::{simplify_bezpath, SimplifyOptLevel, SimplifyOptions},
        BezPath, Point, Rect, Shape, Size,
    };
    use std::{f64::consts::PI, fs::File, io::Write};

    #[test]
    fn helix() {
        // Generate helix points
        let alpha = 0.35;
        let radius = 10.0;
        let steps = 1000;
        let t_range = (0.0, 46.0 * PI);
        let step_size = (t_range.1 - t_range.0) / steps as f64;

        let mut points = Vec::new();
        for i in 0..=steps {
            let t = t_range.0 + i as f64 * step_size;
            let x = (alpha * t).sin() * radius + t;
            let y = (alpha * t).cos() * radius;
            points.push(Point::new(x, y));
        }

        // Create path
        let mut path = BezPath::new();
        path.move_to(points[0]);
        for point in points.iter().skip(1) {
            path.line_to(*point);
        }
        let options = SimplifyOptions::default().opt_level(SimplifyOptLevel::Optimize);
        let simplified = simplify_bezpath(path.clone(), 1.0, &options);

        // Get the bounding box of the path
        let bbox: Rect = path.bounding_box();

        // Define target SVG dimensions and padding
        let svg_width = 800.0;
        let svg_height = 600.0;
        let padding = 50.0;

        // Expand the bounding box by padding (to keep a margin)
        let padded_bbox = Rect::from_origin_size(
            Point::new(bbox.min_x() - padding, bbox.min_y() - padding),
            Size::new(bbox.width() + 2.0 * padding, bbox.height() + 2.0 * padding),
        );

        // Create the SVG with a viewBox set to the padded bounding box.
        let svg_string = format!(
        "<svg width='{width}' height='{height}' viewBox='{min_x} {min_y} {vb_width} {vb_height}' xmlns='http://www.w3.org/2000/svg'>\n  \
         <path d='{path}' stroke='#000' fill='none' stroke-width='2'/>\n</svg>",
        width = svg_width,
        height = svg_height,
        min_x = padded_bbox.min_x(),
        min_y = padded_bbox.min_y(),
        vb_width = padded_bbox.width(),
        vb_height = padded_bbox.height(),
        path = simplified.to_svg(),
    );

        // Write to file
        let mut file = File::create("output.svg").unwrap();
        write!(file, "{svg_string}").unwrap();
    }
}
