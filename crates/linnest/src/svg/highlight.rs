//! Ordered highlights resolved against the drawing's node and edge order.

use std::borrow::Cow;
use std::cmp::Reverse;
use std::collections::BTreeSet;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::mem;
use std::str::FromStr;

use thiserror::Error;

use super::{Dash, EdgeDrawing, NodeDrawing, Scene, Stroke};

/// Geometry used where several highlight regions meet. See [`Scene::apply_highlights_with_overlap`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum HighlightOverlap {
    /// Use the original stroke treatment for one layer and hulls for several.
    #[default]
    Auto,
    /// Nested strokes, with the last layer on the centerline.
    Casings,
    /// Translucent unions of selected edge tubes and vertex discs.
    Hulls,
    /// Recolour edges by the last matching layer.
    Recolour,
    /// Parallel 2.4pt strokes on 3pt centers and segmented vertex rings.
    Bands,
    /// Parallel strokes on 4pt centers for small diagrams.
    BandsMiniature,
    /// Global parallel slots joined by cubic turns around vertices.
    Lanes,
}

impl HighlightOverlap {
    /// Mode name. Python selects `Auto` by omitting the option.
    /// See [`Scene::apply_highlights_with_overlap`].
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Casings => "casings",
            Self::Hulls => "hulls",
            Self::Recolour => "recolour",
            Self::Bands => "bands",
            Self::BandsMiniature => "bands-miniature",
            Self::Lanes => "lanes",
        }
    }

    /// Band center-to-center distance in points, if applicable.
    /// See [`Scene::apply_highlights_with_overlap`].
    pub fn band_spacing(self) -> Option<f64> {
        match self {
            Self::Auto | Self::Casings | Self::Hulls | Self::Recolour | Self::Lanes => None,
            Self::Bands => Some(3.0),
            Self::BandsMiniature => Some(4.0),
        }
    }
}

impl FromStr for HighlightOverlap {
    type Err = UnknownHighlightOverlap;

    fn from_str(mode: &str) -> Result<Self, Self::Err> {
        [
            Self::Casings,
            Self::Bands,
            Self::BandsMiniature,
            Self::Hulls,
            Self::Recolour,
            Self::Lanes,
        ]
        .into_iter()
        .find(|candidate| candidate.as_str() == mode)
        .ok_or(UnknownHighlightOverlap)
    }
}

/// An unknown overlap mode. See [`Scene::apply_highlights_with_overlap`].
#[derive(Clone, Copy, Debug, Error, PartialEq, Eq)]
#[error("overlap must be 'hulls', 'recolour', 'casings', 'bands', 'bands-miniature', or 'lanes'")]
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

/// Subtraction scheme cue on a ghost counterterm. See [`HighlightDetails`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum CountertermScheme {
    /// Plain open circled cross, with no scheme in the hover text.
    #[default]
    Unspecified,
    /// Open circled cross.
    MsBar,
    /// Solid disc with a paper-coloured cross.
    OnShell,
}

impl CountertermScheme {
    /// The scheme as it appears on hover, if specified.
    /// See [`HighlightDetails`].
    pub fn name(self) -> Option<&'static str> {
        match self {
            Self::Unspecified => None,
            Self::MsBar => Some("MS-bar"),
            Self::OnShell => Some("on-shell"),
        }
    }
}

/// Unrecognized scheme keyword. See [`CountertermScheme`].
#[derive(Clone, Copy, Debug, Error, PartialEq, Eq)]
#[error("counterterm scheme must be 'MS-bar', 'on-shell', or 'unspecified'")]
pub struct UnknownCountertermScheme;

impl FromStr for CountertermScheme {
    type Err = UnknownCountertermScheme;

    fn from_str(name: &str) -> Result<Self, Self::Err> {
        if name.eq_ignore_ascii_case("ms-bar") {
            Ok(Self::MsBar)
        } else if name.eq_ignore_ascii_case("on-shell") {
            Ok(Self::OnShell)
        } else if name.eq_ignore_ascii_case("unspecified") {
            Ok(Self::Unspecified)
        } else {
            Err(UnknownCountertermScheme)
        }
    }
}

/// Region hover text and ghost-in-place treatment.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct HighlightDetails {
    /// Short name shown on hover.
    pub label: String,
    /// Longer description shown on hover.
    pub description: String,
    /// Show a ghost counterterm.
    pub contracted: bool,
    /// Optional subtraction scheme cue.
    pub scheme: CountertermScheme,
}

/// Hull radius and increment per inclusion level, in points.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HullOptions {
    /// Radius for an innermost region.
    pub radius: f64,
    /// Extra radius per nesting level.
    pub step: f64,
}

impl Default for HullOptions {
    fn default() -> Self {
        Self {
            radius: 4.2,
            step: 2.6,
        }
    }
}

/// Region inputs to resolve after layout.
#[derive(Clone, Debug)]
pub struct HullRegion {
    pub(super) edges: BTreeSet<usize>,
    pub(super) nodes: BTreeSet<usize>,
    pub(super) paint: String,
    pub(super) radius: f64,
    pub(super) details: HighlightDetails,
}

/// Lane membership resolved against the routed edge carriers.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LaneRegion {
    pub(super) edges: BTreeSet<usize>,
    pub(super) paint: String,
    pub(super) details: HighlightDetails,
}

/// Derived highlight geometry and the paint it temporarily replaces.
/// Initialize with `HighlightDrawing::default()` when constructing a scene.
#[derive(Default)]
pub struct HighlightDrawing {
    pub(super) hulls: Vec<HullRegion>,
    pub(super) lanes: Vec<LaneRegion>,
    pub(super) recolour_stamp: bool,
    original: Option<OriginalPaint>,
}

struct OriginalPaint {
    nodes: Vec<NodeDrawing>,
    edges: Vec<EdgeDrawing>,
}

impl HullOptions {
    /// Validate the base radius and nesting step before resolving highlights.
    ///
    /// # Errors
    /// Returns a message if the base radius is outside (0, 1000] points or
    /// the nesting step is negative or non-finite.
    pub fn validate(self) -> Result<(), &'static str> {
        if !valid_radius(self.radius) || !self.step.is_finite() || self.step < 0.0 {
            return Err(
                "hull_radius must be in (0, 1000] pt and hull_step must be finite and nonnegative",
            );
        }
        Ok(())
    }
}

fn valid_radius(radius: f64) -> bool {
    radius > 0.0 && radius <= 1000.0
}

impl HullRegion {
    fn resolve(
        regions: &[HighlightRegion],
        palette: &HighlightPalette<'_>,
        overlap: HighlightOverlap,
        hulls: HullOptions,
    ) -> Vec<Self> {
        let mut heights = vec![0_u32; regions.len()];
        let mut order: Vec<_> = regions
            .iter()
            .enumerate()
            .filter(|(_, region)| !region.edges.is_empty() || !region.nodes.is_empty())
            .map(|(index, _)| index)
            .collect();
        order.sort_by_key(|&index| regions[index].edges.len() + regions[index].nodes.len());
        for (position, &outer) in order.iter().enumerate() {
            let region = &regions[outer];
            heights[outer] = order[..position]
                .iter()
                .filter_map(|&inner| {
                    let child = &regions[inner];
                    (child.edges.is_subset(&region.edges)
                        && child.nodes.is_subset(&region.nodes)
                        && (child.edges != region.edges || child.nodes != region.nodes))
                        .then_some(heights[inner] + 1)
                })
                .max()
                .unwrap_or(0);
        }
        order.extend(regions.iter().enumerate().filter_map(|(index, region)| {
            (region.details.contracted && region.edges.is_empty() && region.nodes.is_empty())
                .then_some(index)
        }));
        order.sort_by_key(|&index| (Reverse(heights[index]), index));
        order
            .into_iter()
            .filter_map(|index| {
                let region = &regions[index];
                if overlap != HighlightOverlap::Hulls && !region.details.contracted {
                    return None;
                }
                let style = if region.details.contracted {
                    HighlightStyle::Treated
                } else {
                    region.style
                };
                let radius = hulls.radius
                    + hulls.step * f64::from(heights[index])
                    + if region.details.contracted { 1.0 } else { 0.0 };
                Some(HullRegion {
                    edges: region.edges.clone(),
                    nodes: region.nodes.clone(),
                    paint: style.paint(&palette.amber_paint).into_owned(),
                    radius,
                    details: region.details.clone(),
                })
            })
            .collect()
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
    /// Optional hover text and counterterm treatment.
    pub details: HighlightDetails,
}

impl HighlightRegion {
    /// Bind a style to scene edge and node indices. See [`Scene::apply_highlights`].
    pub fn new(edges: BTreeSet<usize>, nodes: BTreeSet<usize>, style: HighlightStyle) -> Self {
        Self {
            edges,
            nodes,
            style,
            details: HighlightDetails::default(),
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

    /// Draw one ordinary layer with strokes, or several with reading-size hulls.
    /// A contracted layer uses a hull and a counterterm glyph.
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
        self.apply_highlights_with_overlap(regions, palette, mute_outside, HighlightOverlap::Auto);
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
        self.apply_highlights_with_hulls(
            regions,
            palette,
            mute_outside,
            overlap,
            HullOptions::default(),
        );
    }

    /// Check effective radii, including nesting and the extra counterterm radius.
    ///
    /// # Errors
    /// Returns a message when a resolved radius is outside (0, 1000] points.
    pub fn validate_hull_radii(&self) -> Result<(), &'static str> {
        if self
            .highlights
            .hulls
            .iter()
            .any(|region| !valid_radius(region.radius))
        {
            return Err("effective hull radius must be in (0, 1000] pt after nesting and counterterm padding; reduce hull_radius or hull_step");
        }
        Ok(())
    }

    /// Resolve highlight geometry with configurable hull dimensions.
    /// See [`Self::apply_highlights`] for a selection example.
    pub fn apply_highlights_with_hulls(
        &mut self,
        regions: &[HighlightRegion],
        palette: &HighlightPalette<'_>,
        mute_outside: bool,
        overlap: HighlightOverlap,
        hulls: HullOptions,
    ) {
        if let Some(original) = mem::take(&mut self.highlights).original {
            self.nodes = original.nodes;
            self.edges = original.edges;
        }
        if regions.is_empty() {
            return;
        }
        let overlap = match overlap {
            HighlightOverlap::Auto if regions.len() < 2 => HighlightOverlap::Casings,
            HighlightOverlap::Auto => HighlightOverlap::Hulls,
            overlap => overlap,
        };
        let legacy = matches!(regions, [region] if region.style == HighlightStyle::Amber && !region.details.contracted && region.details.label.is_empty() && region.details.description.is_empty());
        let has_counterterm = regions.iter().any(|region| region.details.contracted);
        if !legacy && (overlap == HighlightOverlap::Hulls || has_counterterm) {
            self.highlights.hulls = HullRegion::resolve(regions, palette, overlap, hulls);
        }
        if !legacy && overlap == HighlightOverlap::Hulls {
            return;
        }
        let contracted_edges: BTreeSet<_> = regions
            .iter()
            .filter(|region| region.details.contracted)
            .flat_map(|region| region.edges.iter().copied())
            .collect();
        let contracted_nodes: BTreeSet<_> = regions
            .iter()
            .filter(|region| region.details.contracted)
            .flat_map(|region| region.nodes.iter().copied())
            .collect();
        let ordinary;
        let regions = if has_counterterm {
            ordinary = regions
                .iter()
                .filter(|region| !region.details.contracted)
                .cloned()
                .collect::<Vec<_>>();
            ordinary.as_slice()
        } else {
            regions
        };
        self.highlights.recolour_stamp = !legacy && overlap == HighlightOverlap::Recolour;
        if regions.is_empty() && !self.highlights.recolour_stamp {
            return;
        }
        if overlap == HighlightOverlap::Lanes && !legacy {
            self.highlights.lanes = regions
                .iter()
                .map(|region| LaneRegion {
                    edges: region.edges.clone(),
                    paint: region.style.paint(&palette.amber_paint).into_owned(),
                    details: region.details.clone(),
                })
                .collect();
            return;
        }
        self.highlights.original = Some(OriginalPaint {
            nodes: self.nodes.clone(),
            edges: self.edges.clone(),
        });
        if self.highlights.recolour_stamp {
            // At stamp sizes double lines, coils and dashes merge with the
            // coloured membership. Use one plain carrier per particle.
            for edge in &mut self.edges {
                edge.pattern = None;
                edge.flow = None;
                edge.stroke.dash = Dash::Solid;
            }
            for (index, edge) in self.edges.iter_mut().enumerate() {
                if contracted_edges.contains(&index) {
                    continue;
                }
                if let Some(region) = regions.iter().rfind(|region| region.edges.contains(&index)) {
                    edge.stroke.paint = region.style.paint(&palette.amber_paint).into_owned();
                    edge.stroke.width = edge.stroke.width.max(2.2);
                } else if mute_outside && !has_counterterm {
                    edge.stroke.paint = HighlightPalette::OUTSIDE_PAINT.into();
                    edge.stroke.width = 1.0;
                }
            }
            return;
        }
        if let Some(spacing) = overlap.band_spacing().filter(|_| {
            !legacy
                && (has_counterterm
                    || !matches!(regions, [region] if region.style == HighlightStyle::Amber))
        }) {
            self.apply_bands(
                regions,
                palette,
                mute_outside,
                spacing,
                &contracted_edges,
                &contracted_nodes,
            );
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
                None if mute_outside && !contracted_nodes.contains(&index) => {
                    palette.paint_outside_node(node)
                }
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
                if regions.len() > 1 || contracted_edges.contains(&index) {
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
            } else if mute_outside && !contracted_edges.contains(&index) {
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
        contracted_edges: &BTreeSet<usize>,
        contracted_nodes: &BTreeSet<usize>,
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
            if node.highlight_ring.is_empty() && mute_outside && !contracted_nodes.contains(&index)
            {
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
            if edge.highlight_strokes.is_empty()
                && mute_outside
                && !contracted_edges.contains(&index)
            {
                edge.stroke = HighlightPalette::outside_stroke();
            }
        }
    }
}
