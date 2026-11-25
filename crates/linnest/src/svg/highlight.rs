//! Ordered highlights resolved against the drawing's node and edge order.

use std::borrow::Cow;
use std::collections::BTreeSet;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::str::FromStr;

use thiserror::Error;

use super::{Dash, NodeDrawing, Scene, Stroke};

/// Geometry used where several highlight regions meet. See [`Scene::apply_highlights_with_overlap`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum HighlightOverlap {
    /// Nested strokes, with the last layer on the centerline.
    #[default]
    Casings,
    /// Parallel 2.4pt strokes on 3pt centers and segmented vertex rings.
    Bands,
    /// Parallel strokes on 4pt centers for small diagrams.
    BandsMiniature,
}

impl HighlightOverlap {
    /// Name used by the Python and Typst rendering options.
    /// See [`Scene::apply_highlights_with_overlap`].
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Casings => "casings",
            Self::Bands => "bands",
            Self::BandsMiniature => "bands-miniature",
        }
    }

    /// Band center-to-center distance in points, if applicable.
    /// See [`Scene::apply_highlights_with_overlap`].
    pub fn band_spacing(self) -> Option<f64> {
        match self {
            Self::Casings => None,
            Self::Bands => Some(3.0),
            Self::BandsMiniature => Some(4.0),
        }
    }
}

impl FromStr for HighlightOverlap {
    type Err = UnknownHighlightOverlap;

    fn from_str(mode: &str) -> Result<Self, Self::Err> {
        [Self::Casings, Self::Bands, Self::BandsMiniature]
            .into_iter()
            .find(|candidate| candidate.as_str() == mode)
            .ok_or(UnknownHighlightOverlap)
    }
}

/// An unknown overlap mode. See [`Scene::apply_highlights_with_overlap`].
#[derive(Clone, Copy, Debug, Error, PartialEq, Eq)]
#[error("overlap must be 'casings', 'bands', or 'bands-miniature'")]
pub struct UnknownHighlightOverlap;

/// Theme paints for an automatic palette slot.
/// See [`Scene::apply_highlights_with_overlap`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SlotPaint {
    /// SVG fallback and light-theme paint.
    pub light: String,
    /// Dark-theme paint.
    pub dark: String,
}

impl SlotPaint {
    /// The shared CFF and highlight palette. Slots are one-based. The first
    /// seven use Tol qualitative hues, then golden-angle hues extend the palette.
    ///
    /// ```
    /// use linnest::svg::SlotPaint;
    /// assert_eq!(SlotPaint::new(1).light, "#0077bb");
    /// assert_ne!(SlotPaint::new(8).light, SlotPaint::new(8).dark);
    /// ```
    pub fn new(slot: usize) -> Self {
        const FIXED: [(&str, &str); 7] = [
            ("#0077bb", "#77bbee"),
            ("#ee7733", "#ffaa66"),
            ("#009988", "#44ccbb"),
            ("#ee3377", "#ff77aa"),
            ("#33bbee", "#66ddff"),
            ("#cc3311", "#ff8866"),
            ("#aa3377", "#cc88cc"),
        ];
        if let Some(&(light, dark)) = slot.checked_sub(1).and_then(|i| FIXED.get(i)) {
            return Self {
                light: light.into(),
                dark: dark.into(),
            };
        }
        let hue = (slot as f64 * 137.508) % 360.0;
        let hsl = |saturation: f64, lightness: f64| {
            let amplitude = saturation * lightness.min(1.0 - lightness);
            let channel = |shift: f64| {
                let k = (shift + hue / 30.0) % 12.0;
                let value = lightness - amplitude * (-1.0_f64).max((k - 3.0).min(9.0 - k).min(1.0));
                (value * 255.0).round() as u8
            };
            format!(
                "#{:02x}{:02x}{:02x}",
                channel(0.0),
                channel(8.0),
                channel(4.0)
            )
        };
        Self {
            light: hsl(0.55, 0.38),
            dark: hsl(0.65, 0.73),
        }
    }
}

/// A named highlight treatment. See [`Scene::apply_highlights`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
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
    /// An automatic one-based palette slot.
    Slot(usize),
}

impl HighlightStyle {
    /// Named styles available in both themes. Amber is renderer-specific.
    /// See [`Scene::apply_highlights`].
    pub const LAYERED: [Self; 3] = [Self::ActedOn, Self::Middle, Self::Treated];

    /// Stable role name used by Python and CSS. See [`Scene::apply_highlights`].
    pub fn role(self) -> &'static str {
        match self {
            Self::Amber => "amber",
            Self::ActedOn => "acted-on",
            Self::Middle => "middle",
            Self::Treated => "treated",
            Self::Slot(_) => "slot",
        }
    }

    /// Light paint. Amber uses the supplied palette. See [`Scene::apply_highlights`].
    pub fn light_paint(self) -> Option<&'static str> {
        match self {
            Self::Amber | Self::Slot(_) => None,
            Self::ActedOn => Some("#3b82f6"),
            Self::Middle => Some("#8b5cf6"),
            Self::Treated => Some("#e4504f"),
        }
    }

    /// Resolve the light paint, including automatic slots and owner-specific amber.
    /// See [`Scene::apply_highlights_with_overlap`].
    pub fn paint(self, amber: &str) -> Cow<'_, str> {
        match self {
            Self::Slot(slot) => Cow::Owned(SlotPaint::new(slot).light),
            style => Cow::Borrowed(style.light_paint().unwrap_or(amber)),
        }
    }

    /// Dark paint. Amber uses the supplied palette. See [`Scene::apply_highlights`].
    pub fn dark_paint(self) -> Option<&'static str> {
        match self {
            Self::Amber | Self::Slot(_) => None,
            Self::ActedOn => Some("#75a7ff"),
            Self::Middle => Some("#b197fc"),
            Self::Treated => Some("#ff7775"),
        }
    }

    /// Preferred width in points. Amber can inherit the particle's width.
    /// See [`Scene::apply_highlights`].
    pub fn width(self) -> Option<f64> {
        match self {
            Self::Amber | Self::Slot(_) => None,
            Self::ActedOn => Some(1.45),
            Self::Middle => Some(1.55),
            Self::Treated => Some(1.65),
        }
    }
}

/// The input did not name a highlight style.
/// See [`Scene::apply_highlights`].
#[derive(Clone, Copy, Debug, Error, PartialEq, Eq)]
#[error("highlight style must be 'amber', 'acted-on', 'middle', or 'treated'")]
pub struct UnknownHighlightStyle;

impl FromStr for HighlightStyle {
    type Err = UnknownHighlightStyle;

    fn from_str(role: &str) -> Result<Self, Self::Err> {
        [Self::Amber]
            .into_iter()
            .chain(Self::LAYERED)
            .find(|style| style.role() == role)
            .ok_or(UnknownHighlightStyle)
    }
}

impl Display for HighlightStyle {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> FmtResult {
        match self {
            Self::Slot(slot) => write!(formatter, "slot-{slot}"),
            style => formatter.write_str(style.role()),
        }
    }
}

/// One styled region, indexed in scene order. See [`Scene::apply_highlights`].
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub struct HighlightRegion {
    /// Selected indices in the scene's edge order.
    pub edges: BTreeSet<usize>,
    /// Selected indices in the scene's node order.
    pub nodes: BTreeSet<usize>,
    /// Treatment of this region.
    pub style: HighlightStyle,
}

impl HighlightRegion {
    /// Bind a style to scene edge and node indices. See [`Scene::apply_highlights`].
    pub fn new(edges: BTreeSet<usize>, nodes: BTreeSet<usize>, style: HighlightStyle) -> Self {
        Self {
            edges,
            nodes,
            style,
        }
    }
}

/// The owner-specific amber and outside paints. See [`Scene::apply_highlights`].
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct HighlightPalette<'a> {
    /// Amber paint for edge and node outlines.
    pub amber_paint: Cow<'a, str>,
    /// Solid amber width in points. Without an override, inherit edge width and dash.
    pub amber_solid_width: Option<f64>,
    /// Replacement fill for amber nodes, if any.
    pub amber_node_fill: Option<Cow<'a, str>>,
    /// Replacement fill for outside nodes, if any.
    pub outside_node_fill: Option<Cow<'a, str>>,
}

impl<'a> HighlightPalette<'a> {
    /// Set the amber outline and otherwise inherit the drawing's styles.
    /// See [`Scene::apply_highlights`].
    pub fn new(amber_paint: impl Into<Cow<'a, str>>) -> Self {
        Self {
            amber_paint: amber_paint.into(),
            amber_solid_width: None,
            amber_node_fill: None,
            outside_node_fill: None,
        }
    }
}

/// An edge stroke and its global paint rank. See [`Scene::apply_highlights`].
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct HighlightStroke {
    /// Earlier layers paint underneath later layers.
    pub layer: usize,
    /// Resolved paint and casing width.
    pub stroke: Stroke,
    /// Signed band offset in points. Absent for casings.
    pub offset: Option<f64>,
}

impl HighlightPalette<'_> {
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
        if let Some(fill) = &self.outside_node_fill {
            node.fill = fill.as_ref().into();
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
    /// Named layer styles keep the particle's dash and base node appearance. Amber uses the
    /// owner's treatment. Outside paint applies only when at least one region
    /// is supplied and `mute_outside` is true.
    ///
    /// ```
    /// use linnest::svg::{HighlightPalette, HighlightRegion, HighlightStyle, Scene};
    /// # fn highlight(scene: &mut Scene) {
    /// let palette = HighlightPalette::new("#ffd166");
    /// let region = HighlightRegion::new([0, 1].into(), [0, 1].into(), HighlightStyle::ActedOn);
    /// scene.apply_highlights(&[region], &palette, true);
    /// # }
    /// ```
    pub fn apply_highlights(
        &mut self,
        regions: &[HighlightRegion],
        palette: &HighlightPalette<'_>,
        mute_outside: bool,
    ) {
        self.apply_highlights_with_overlap(
            regions,
            palette,
            mute_outside,
            HighlightOverlap::Casings,
        );
    }

    /// Paint the selected regions using the requested overlap geometry.
    /// The legacy single amber selection keeps its original appearance.
    /// See [`Self::apply_highlights`] for a scene and selection example.
    pub fn apply_highlights_with_overlap(
        &mut self,
        regions: &[HighlightRegion],
        palette: &HighlightPalette<'_>,
        mute_outside: bool,
        overlap: HighlightOverlap,
    ) {
        if regions.is_empty() {
            return;
        }
        if let Some(spacing) = overlap
            .band_spacing()
            .filter(|_| !matches!(regions, [region] if region.style == HighlightStyle::Amber))
        {
            self.apply_bands(regions, palette, mute_outside, spacing);
            return;
        }
        for (index, node) in self.nodes.iter_mut().enumerate() {
            node.highlight_ring.clear();
            match regions.iter().rfind(|region| region.nodes.contains(&index)) {
                Some(region) if region.style == HighlightStyle::Amber => {
                    node.stroke = Stroke {
                        paint: palette.amber_paint.as_ref().into(),
                        width: 1.2,
                        dash: Dash::Solid,
                        round_cap: false,
                    };
                    if let Some(fill) = &palette.amber_node_fill {
                        node.fill = fill.as_ref().into();
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
                    let style = region.style;
                    let amber_solid_width = palette
                        .amber_solid_width
                        .filter(|_| style == HighlightStyle::Amber);
                    HighlightStroke {
                        layer,
                        offset: None,
                        stroke: Stroke {
                            paint: style.paint(palette.amber_paint.as_ref()).into_owned(),
                            width: style
                                .width()
                                .or(amber_solid_width)
                                .unwrap_or(edge.stroke.width),
                            dash: if amber_solid_width.is_some() {
                                Dash::Solid
                            } else {
                                edge.stroke.dash
                            },
                            round_cap: amber_solid_width.is_none() && edge.stroke.round_cap,
                        },
                    }
                });
            edge.highlight_strokes.clear();
            if let Some(top) = strokes.next() {
                if regions.len() > 1 {
                    let mut width = top.stroke.width;
                    edge.highlight_strokes.push(top);
                    for mut casing in strokes {
                        width = casing.stroke.width.max(width + Self::CASING_WIDTH_STEP);
                        casing.stroke.width = width;
                        casing.stroke.dash = Dash::Solid;
                        edge.highlight_strokes.push(casing);
                    }
                    edge.highlight_strokes.reverse();
                } else {
                    edge.stroke = top.stroke;
                }
            } else if mute_outside {
                edge.stroke = HighlightPalette::outside_stroke();
            }
        }
    }
    fn apply_bands(
        &mut self,
        regions: &[HighlightRegion],
        palette: &HighlightPalette<'_>,
        mute_outside: bool,
        spacing: f64,
    ) {
        let colors: Vec<_> = regions
            .iter()
            .map(|region| region.style.paint(&palette.amber_paint))
            .collect();
        for (index, node) in self.nodes.iter_mut().enumerate() {
            node.highlight_ring = regions
                .iter()
                .zip(&colors)
                .filter(|(region, _)| region.nodes.contains(&index))
                .map(|(_, paint)| paint.to_string())
                .collect();
            if node.highlight_ring.is_empty() && mute_outside {
                palette.paint_outside_node(node);
            }
        }
        for (index, edge) in self.edges.iter_mut().enumerate() {
            let members: Vec<_> = regions
                .iter()
                .zip(&colors)
                .enumerate()
                .filter(|(_, (region, _))| region.edges.contains(&index))
                .collect();
            let center = (members.len() as f64 - 1.0) / 2.0;
            edge.highlight_strokes = members
                .into_iter()
                .enumerate()
                .map(|(slot, (layer, (_, paint)))| HighlightStroke {
                    layer,
                    offset: Some((slot as f64 - center) * spacing),
                    stroke: Stroke {
                        paint: paint.to_string(),
                        width: 2.4,
                        dash: Dash::Solid,
                        round_cap: true,
                    },
                })
                .collect();
            if edge.highlight_strokes.is_empty() && mute_outside {
                edge.stroke = HighlightPalette::outside_stroke();
            }
        }
    }
}
