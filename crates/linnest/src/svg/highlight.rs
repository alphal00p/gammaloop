//! Ordered highlights resolved against the drawing's node and edge order.

use std::collections::BTreeSet;

use super::{Dash, NodeDrawing, Scene, Stroke};

/// A named highlight treatment. See [`Scene::apply_highlights`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum HighlightStyle {
    /// Keep the owner's traditional amber stroke and node treatment.
    #[default]
    Amber,
    /// Blue for the region an operator acts on.
    ActedOn,
    /// Violet for an intermediate nested region.
    Middle,
    /// Red for the inner, already-treated region.
    Treated,
}

impl HighlightStyle {
    /// Theme role, light fallback paint and dark paint for the layered palette.
    /// Amber retains the palette of the graph's renderer. See [`Scene::apply_highlights`].
    pub const THEME: [(&'static str, &'static str, &'static str); 3] = [
        ("acted-on", "#3b82f6", "#75a7ff"),
        ("middle", "#8b5cf6", "#b197fc"),
        ("treated", "#e4504f", "#ff7775"),
    ];

    /// Named stroke paint and preferred width in points. Amber instead derives
    /// its stroke from the graph's renderer. See [`Scene::apply_highlights`].
    pub fn stroke(self) -> Option<(&'static str, f64)> {
        match self {
            Self::Amber => None,
            Self::ActedOn => Some((Self::THEME[0].1, 1.45)),
            Self::Middle => Some((Self::THEME[1].1, 1.55)),
            Self::Treated => Some((Self::THEME[2].1, 1.65)),
        }
    }
}

/// One styled region, indexed in scene order. See [`Scene::apply_highlights`].
pub struct HighlightRegion {
    pub edges: BTreeSet<usize>,
    pub nodes: BTreeSet<usize>,
    pub style: HighlightStyle,
}

/// The owner-specific amber and outside paints. See [`Scene::apply_highlights`].
pub struct HighlightPalette {
    pub amber_paint: &'static str,
    /// Override the base edge width and dash for amber if provided.
    pub amber_edge_width: Option<f64>,
    pub amber_node_fill: Option<&'static str>,
    pub outside_node_fill: Option<&'static str>,
}

impl HighlightPalette {
    /// The muted complement stroke. See [`Scene::apply_highlights`].
    pub const OUTSIDE_PAINT: &'static str = "#77777773";

    fn outside_stroke() -> Stroke {
        Stroke {
            paint: Self::OUTSIDE_PAINT.into(),
            width: 0.6,
            dash: Dash::Dotted,
            round_cap: false,
        }
    }

    fn paint_outside_node(&self, node: &mut NodeDrawing) {
        node.stroke = Self::outside_stroke();
        if let Some(fill) = self.outside_node_fill {
            node.fill = fill.into();
        }
    }
}

impl Scene {
    /// Minimum difference in stroke width between adjacent casings, in points.
    /// The extra width is split evenly across the two sides of the centerline.
    /// See [`Self::apply_highlights`].
    pub const CASING_WIDTH_STEP: f64 = 1.1;

    /// Paint regions in input order. The innermost stroke keeps its requested
    /// width and earlier strokes grow only as needed to remain visible.
    /// Named layered styles use solid strokes on the existing edge geometry
    /// and keep the base node appearance. Amber uses the
    /// owner's treatment. Outside paint applies only when at least one region
    /// is supplied and `mute_outside` is true.
    ///
    /// ```
    /// use linnest::svg::{HighlightPalette, HighlightRegion, HighlightStyle, Scene};
    /// # fn highlight(scene: &mut Scene) {
    /// let palette = HighlightPalette {
    ///     amber_paint: "#ffd166",
    ///     amber_edge_width: None,
    ///     amber_node_fill: None,
    ///     outside_node_fill: None,
    /// };
    /// let region = HighlightRegion {
    ///     edges: [0, 1].into(),
    ///     nodes: [0, 1].into(),
    ///     style: HighlightStyle::ActedOn,
    /// };
    /// scene.apply_highlights(&[region], &palette, true);
    /// # }
    /// ```
    pub fn apply_highlights(
        &mut self,
        regions: &[HighlightRegion],
        palette: &HighlightPalette,
        mute_outside: bool,
    ) {
        if regions.is_empty() {
            return;
        }
        for (index, node) in self.nodes.iter_mut().enumerate() {
            match regions.iter().rfind(|region| region.nodes.contains(&index)) {
                Some(region) if region.style == HighlightStyle::Amber => {
                    node.stroke = Stroke {
                        paint: palette.amber_paint.into(),
                        width: 1.2,
                        dash: Dash::Solid,
                        round_cap: false,
                    };
                    if let Some(fill) = palette.amber_node_fill {
                        node.fill = fill.into();
                    }
                }
                None if mute_outside => palette.paint_outside_node(node),
                _ => {}
            }
        }
        for (index, edge) in self.edges.iter_mut().enumerate() {
            let mut strokes = regions
                .iter()
                .enumerate()
                .filter(|(_, region)| region.edges.contains(&index))
                .rev()
                .map(|(layer, region)| {
                    let named = region.style.stroke();
                    let (paint, width) = named.unwrap_or((
                        palette.amber_paint,
                        palette.amber_edge_width.unwrap_or(edge.stroke.width),
                    ));
                    let solid = named.is_some() || palette.amber_edge_width.is_some();
                    let stroke = Stroke {
                        paint: paint.into(),
                        width,
                        dash: if solid { Dash::Solid } else { edge.stroke.dash },
                        round_cap: named.is_some() || (!solid && edge.stroke.round_cap),
                    };
                    (layer, stroke)
                });
            edge.highlight_strokes.clear();
            if let Some((layer, top)) = strokes.next() {
                if regions.len() > 1 {
                    let mut width = top.width;
                    edge.highlight_strokes.push((layer, top));
                    for (layer, mut stroke) in strokes {
                        width = stroke.width.max(width + Self::CASING_WIDTH_STEP);
                        stroke.width = width;
                        edge.highlight_strokes.push((layer, stroke));
                    }
                    edge.highlight_strokes.reverse();
                } else {
                    edge.stroke = top;
                }
            } else if mute_outside {
                edge.stroke = HighlightPalette::outside_stroke();
            }
        }
    }
}
