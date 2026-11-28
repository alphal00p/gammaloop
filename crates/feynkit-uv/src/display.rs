//! Native drawings of spinneys, forests, counterterms, woods and forest
//! operation graphs of FeynKit diagrams.
//!
//! Single drawings of the parent diagram are returned as highlight
//! [`UvRegion`]s for [`FeynmanDiagram::to_scene_with_highlights`]: a spinney is
//! one region, a forest stacks its spinneys as nested regions and a
//! counterterm contracts them in place. Woods and operation graphs are
//! complete SVG documents whose nodes are miniature stamps of the diagram,
//! recoloured with the innermost region winning: the spinney being
//! approximated is blue and the one already treated is red.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::slice;

use feynkit_graph::{
    CountertermScheme as Glyph, DiagramError, FeynmanDiagram, HighlightDetails, HighlightLayer,
    HighlightOverlap, HighlightStyle, LoopMomentumBasis, MINIATURE_STAMP_WEIGHT, SceneOptions,
};
use linnest::svg::{
    CANVAS_INSET, Config, Dash, Details, EdgeDrawing, LabelMetrics, LabelPage, LabelSource,
    NodeDrawing, Scene, Stroke, Typeset,
};
use linnest::{TypstEdgeSpec, TypstEndpointSpec, TypstGraphSpec, TypstNodeSpec};
use linnet::half_edge::NodeIndex;
use linnet::half_edge::subgraph::{Inclusion, SuBitGraph, SubSetLike, SubSetOps};
use serde_json::{Map, Value, json};
use thiserror::Error;

use crate::{
    Counterterm, CountertermScheme, Expansion, ForestTerm, Forests, Spinney, TermFactor, Wood,
};

/// Ink of the drawing, matching FeynKit's physics palette.
const INK: &str = "#3d2645";
/// Box outline, connector and fills of the `tree-auto-root` operation graphs.
const BOX_STROKE: &str = "#60656f";
const CONNECTOR: &str = "#7890b4";
/// Fill of the unsubtracted root.
pub const ROOT_FILL: &str = "#fff2b8";
/// Fill of forests and spinneys that cover every line of every spinney.
pub const FULL_FILL: &str = "#dff0ff";
/// Default stamp size, in points: the larger side of the bare diagram stamp.
pub const STAMP_SIZE: f64 = 28.0;

/// Size and line weight of the miniature stamps in woods, operation graphs and
/// expansions.
///
/// ```
/// use feynkit_uv::{STAMP_SIZE, StampOptions};
///
/// let heavier = StampOptions { weight: 1.6, ..StampOptions::default() };
/// assert_eq!(heavier.size, STAMP_SIZE);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StampOptions {
    /// The larger side of the bare stamp, in points.
    pub size: f64,
    /// Line weight multiplier; [`MINIATURE_STAMP_WEIGHT`] by default.
    pub weight: f64,
}

impl Default for StampOptions {
    fn default() -> Self {
        Self {
            size: STAMP_SIZE,
            weight: MINIATURE_STAMP_WEIGHT,
        }
    }
}

/// Failure to draw a UV object.
#[derive(Debug, Error)]
pub enum UvDisplayError {
    /// The diagram could not build a drawing scene.
    #[error(transparent)]
    Diagram(#[from] DiagramError),
    /// Layout or SVG output failed.
    #[error("failed to draw: {0}")]
    Render(String),
    /// The stamp size is not a positive finite number of points.
    #[error("stamp size must be positive and finite, got {0}")]
    StampSize(f64),
    /// The stamp line weight is not a positive finite multiplier.
    #[error("stamp weight must be positive and finite, got {0}")]
    StampWeight(f64),
}

/// One highlighted region of a UV drawing: half-edges, paint and hover details.
///
/// Build the diagram's scene from regions with [`Self::layers`]:
///
/// ```
/// use feynkit_graph::{FeynmanDiagram, SceneOptions};
/// use feynkit_model::Model;
/// use feynkit_uv::{RenormalizationPrescriptionSettings, UvGraph, UvRegion};
///
/// let sunrise = FeynmanDiagram::from_dot(
///     Model::phi4(),
///     "digraph sunrise { edge [particle=phi]; ext [style=invis]; \
///          ext -> a; a -> b; a -> b; a -> b; b -> ext; }",
/// )?;
/// let spinneys = sunrise.classify_spinneys(
///     &sunrise.hedge_graph().full_filter(),
///     sunrise.model(),
///     &RenormalizationPrescriptionSettings::default(),
///     sunrise.loop_momentum_basis(),
/// )?;
/// let full = spinneys.iter().max_by_key(|spinney| spinney.loop_count()).unwrap();
/// let regions = [full.region(&sunrise)];
/// assert!(regions[0].details.description.contains("2 loops"));
/// let scene = sunrise.to_scene_with_highlights(
///     &UvRegion::layers(&regions),
///     None,
///     &SceneOptions::default(),
/// )?;
/// assert_eq!(scene.edges.len(), 5);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UvRegion {
    /// The selected half-edges.
    pub hedges: SuBitGraph,
    /// Paint of the region.
    pub style: HighlightStyle,
    /// Hover text, and the counterterm treatment of contracted regions.
    pub details: HighlightDetails,
}

/// No isolated vertices: UV regions select lines only.
static NO_ISOLATED: BTreeSet<usize> = BTreeSet::new();

impl UvRegion {
    /// Borrow `regions` as highlight layers, in order.
    ///
    /// See [`UvRegion`] for an example.
    pub fn layers(regions: &[Self]) -> Vec<HighlightLayer<'_>> {
        regions
            .iter()
            .map(|region| {
                let mut layer = HighlightLayer::new(&region.hedges, &NO_ISOLATED, region.style);
                layer.details = region.details.clone();
                layer
            })
            .collect()
    }
}

impl From<CountertermScheme> for Glyph {
    fn from(scheme: CountertermScheme) -> Self {
        match scheme {
            CountertermScheme::Unspecified => Self::Unspecified,
            CountertermScheme::MsBar => Self::MsBar,
            CountertermScheme::OnShell => Self::OnShell,
        }
    }
}

impl Spinney<LoopMomentumBasis> {
    /// Loop count, degrees and scheme, for hover text.
    ///
    /// See [`UvRegion`] for an example.
    pub fn description(&self, diagram: &FeynmanDiagram) -> String {
        let loops = self.loop_count();
        let mut text = format!(
            "{loops} loop{}, degree {}",
            if loops == 1 { "" } else { "s" },
            self.dod()
        );
        if let Ok(superficial) = self.superficial_degree(diagram) {
            let _ = write!(text, " (superficial {superficial})");
        }
        if self.n_components() > 1 {
            let _ = write!(text, ", {} components", self.n_components());
        }
        let _ = write!(text, ", {} approximation", self.renormalization_scheme());
        text
    }

    /// The spinney as one region: blue, with hover details.
    ///
    /// See [`UvRegion`] for an example.
    pub fn region(&self, diagram: &FeynmanDiagram) -> UvRegion {
        UvRegion {
            hedges: self.filter().clone(),
            style: HighlightStyle::ActedOn,
            details: HighlightDetails {
                label: self.label(diagram),
                description: self.description(diagram),
                ..HighlightDetails::default()
            },
        }
    }
}

/// Stack `spinneys` outermost first, painted by nesting depth: the outermost
/// blue, the innermost red and those in between violet. Disjoint spinneys at
/// the same depth share a paint.
fn nested_regions(
    spinneys: &[&Spinney<LoopMomentumBasis>],
    diagram: &FeynmanDiagram,
) -> Vec<UvRegion> {
    let depths: Vec<usize> = spinneys
        .iter()
        .map(|spinney| {
            spinneys
                .iter()
                .filter(|other| {
                    other.filter() != spinney.filter()
                        && other.subgraph().includes(spinney.subgraph())
                })
                .count()
        })
        .collect();
    let deepest = depths.iter().copied().max().unwrap_or(0);
    let mut order: Vec<usize> = (0..spinneys.len()).collect();
    order.sort_by_key(|&index| (depths[index], spinneys[index].filter().clone()));
    order
        .into_iter()
        .map(|index| {
            let mut region = spinneys[index].region(diagram);
            region.style = match depths[index] {
                0 => HighlightStyle::ActedOn,
                depth if depth == deepest => HighlightStyle::Treated,
                _ => HighlightStyle::Middle,
            };
            region
        })
        .collect()
}

impl Forests<LoopMomentumBasis> {
    /// The forest at `node` as nested regions, outermost first: blue outside,
    /// red innermost and violet in between. The empty forest has none.
    ///
    /// ```
    /// use feynkit_graph::FeynmanDiagram;
    /// use feynkit_model::Model;
    /// use feynkit_uv::{RenormalizationPrescriptionSettings, UvGraph, Wood};
    /// use feynkit_graph::HighlightStyle;
    ///
    /// let eye = FeynmanDiagram::from_dot(
    ///     Model::phi4(),
    ///     "digraph eye { edge [particle=phi]; ext [style=invis]; ext -> a; ext -> a; \
    ///          a -> b; a -> c; b -> c; b -> c; b -> ext; c -> ext; }",
    /// )?;
    /// let spinneys = eye.classify_spinneys(
    ///     &eye.hedge_graph().full_filter(),
    ///     eye.model(),
    ///     &RenormalizationPrescriptionSettings::default(),
    ///     eye.loop_momentum_basis(),
    /// )?;
    /// let forests = Wood::from_spinneys(spinneys, &eye, eye.loop_momentum_basis())?.unfold();
    /// let (node, _) = forests.operations().find(|(_, forest)| forest.depth() == 2).unwrap();
    /// let styles: Vec<_> = forests.regions(node, &eye).into_iter().map(|r| r.style).collect();
    /// assert_eq!(styles, [HighlightStyle::ActedOn, HighlightStyle::Treated]);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if `node` is not a forest of these forests.
    pub fn regions(&self, node: NodeIndex, diagram: &FeynmanDiagram) -> Vec<UvRegion> {
        nested_regions(&self.forest_spinneys(node), diagram)
    }

    /// Draw the forests as an operation graph in the style of Linnest's
    /// `tree-auto-root` example, each box holding its forest's term.
    ///
    /// A term is written `Op[m(S|G) · inner]`: the stamp `m(S|G)` shows the
    /// approximated spinney `S` in blue over the treated spinneys `G` in red, and
    /// `inner` is the nested term for `G`; disjoint factors multiply. The empty
    /// forest is `𝟙`. Boxes carry no cofactor; see [`Self::expansion_svg`] for
    /// the full formula. Boxes are ranked left to right by the size of their
    /// cover; the unsubtracted root is yellow and forests covering every line of
    /// every spinney are light blue. Hover shows each term as `K_{S|G}`, its
    /// stored quantity `k` and its sign.
    ///
    /// `expansion` is [`Forests::expansion`] of these forests, possibly with
    /// [renamed](Expansion::rename) spinneys. `config` applies to every stamp,
    /// so stamps share the layout of the diagram's reading-size drawing. `op`
    /// is the operator symbol, such as `K`.
    ///
    /// # Errors
    ///
    /// Fails when a stamp cannot be drawn, when layout fails or when the stamp
    /// options are not positive and finite.
    ///
    /// ```
    /// use feynkit_graph::FeynmanDiagram;
    /// use feynkit_model::Model;
    /// use feynkit_uv::{RenormalizationPrescriptionSettings, StampOptions, UvGraph, Wood};
    /// use linnest::svg::Config;
    ///
    /// let sunrise = FeynmanDiagram::from_dot(
    ///     Model::phi3(),
    ///     "digraph sunrise { edge [particle=phi]; a -> b; a -> b; a -> b; }",
    /// )?;
    /// let spinneys = sunrise.classify_spinneys(
    ///     &sunrise.hedge_graph().full_filter(),
    ///     sunrise.model(),
    ///     &RenormalizationPrescriptionSettings::default(),
    ///     sunrise.loop_momentum_basis(),
    /// )?;
    /// let forests = Wood::from_spinneys(spinneys, &sunrise, sunrise.loop_momentum_basis())?.unfold();
    /// let expansion = forests.expansion(&sunrise);
    /// let svg = forests.operation_graph_svg(
    ///     &sunrise,
    ///     &expansion,
    ///     &Config::default(),
    ///     StampOptions::default(),
    ///     "K",
    /// )?;
    /// assert_eq!(svg.matches("data-feynkit-forest=").count(), forests.n_forests());
    /// assert!(svg.contains("<title>K_{Γ|a} = k5 · depth 2, sign +</title>"));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn operation_graph_svg(
        &self,
        diagram: &FeynmanDiagram,
        expansion: &Expansion,
        config: &Config,
        stamps: StampOptions,
        op: &str,
    ) -> Result<String, UvDisplayError> {
        let mut sheet = StampSheet::new(diagram, config, stamps)?;
        let domain = self.wood().divergent_cover(diagram);
        let nodes: Vec<NodeIndex> = self.operations().map(|(node, _)| node).collect();
        let index: BTreeMap<NodeIndex, usize> = nodes
            .iter()
            .enumerate()
            .map(|(position, node)| (*node, position))
            .collect();
        let covers: Vec<SuBitGraph> = nodes
            .iter()
            .map(|node| {
                self.operation(*node)
                    .covers()
                    .unwrap_or_else(|| diagram.underlying().empty_subgraph())
            })
            .collect();
        let ranks = ranks_by_size(covers.iter().map(SubSetLike::n_included));
        let mut boxes = Vec::with_capacity(nodes.len());
        for (position, node) in nodes.iter().enumerate() {
            let term = expansion
                .term(*node)
                .expect("every forest has a term in its expansion");
            let block = if term.factors.is_empty() {
                Block::text("𝟙")
            } else {
                sheet.factors(self.wood(), &term.factors, op)?
            };
            let fill = if term.factors.is_empty() {
                ROOT_FILL
            } else if covers[position] == domain {
                FULL_FILL
            } else {
                "none"
            };
            boxes.push(GraphBox {
                block,
                title: term_title(expansion, term, op),
                rank: ranks[position] + 1,
                fill,
                attribute: format!(
                    r#"data-feynkit-forest="{position}" data-feynkit-k="{}""#,
                    term.k
                ),
            });
        }
        let links = self
            .edges()
            .map(|edge| (index[&edge.parent], index[&edge.child]))
            .collect::<Vec<_>>();
        sheet.render(boxes, &links, Direction::Right)
    }

    /// Draw the forest formula term by term, one depth per group of lines: a
    /// depth that does not fit one line wraps by width, its continuation lines
    /// indented, with more space between depths than between wrapped lines.
    /// The formula reads `m(Γ|∅) − m(Γ|γ)·Op[m(γ|∅)] ⋯`, where each term carries its
    /// cofactor, the diagram `Γ` in blue with the term's outermost spinneys in red.
    ///
    /// `expansion`, `config`, `stamps` and `op` are as in
    /// [`Self::operation_graph_svg`].
    /// Hover shows each term as `K_{S|G}` and its stored quantity `k`.
    ///
    /// # Errors
    ///
    /// Fails as [`Self::operation_graph_svg`] does.
    ///
    /// ```
    /// use feynkit_graph::FeynmanDiagram;
    /// use feynkit_model::Model;
    /// use feynkit_uv::{RenormalizationPrescriptionSettings, StampOptions, UvGraph, Wood};
    /// use linnest::svg::Config;
    ///
    /// let sunrise = FeynmanDiagram::from_dot(
    ///     Model::phi3(),
    ///     "digraph sunrise { edge [particle=phi]; a -> b; a -> b; a -> b; }",
    /// )?;
    /// let spinneys = sunrise.classify_spinneys(
    ///     &sunrise.hedge_graph().full_filter(),
    ///     sunrise.model(),
    ///     &RenormalizationPrescriptionSettings::default(),
    ///     sunrise.loop_momentum_basis(),
    /// )?;
    /// let forests = Wood::from_spinneys(spinneys, &sunrise, sunrise.loop_momentum_basis())?.unfold();
    /// let expansion = forests.expansion(&sunrise);
    /// let svg = forests.expansion_svg(
    ///     &sunrise,
    ///     &expansion,
    ///     &Config::default(),
    ///     StampOptions::default(),
    ///     "K",
    /// )?;
    /// assert_eq!(svg.matches("data-feynkit-term=").count(), 8);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn expansion_svg(
        &self,
        diagram: &FeynmanDiagram,
        expansion: &Expansion,
        config: &Config,
        stamps: StampOptions,
        op: &str,
    ) -> Result<String, UvDisplayError> {
        let mut sheet = StampSheet::new(diagram, config, stamps)?;
        // Width before a line wraps, in multiples of the stamp size.
        let max_width = LINE_STAMPS * sheet.options.size;
        let mut depths: Vec<(usize, Vec<Block>)> = Vec::new();
        for term in expansion.terms() {
            let depth = self.operation(term.forest).depth();
            let block = sheet.term(self.wood(), expansion, term, op, term.k != 0)?;
            match depths.last_mut() {
                Some((last, blocks)) if *last == depth => blocks.push(block),
                _ => depths.push((depth, vec![block])),
            }
        }
        let mut lines = Vec::new();
        for (depth, blocks) in depths {
            let mut first = true;
            let mut pieces: Vec<Block> = Vec::new();
            let mut width = 0.0;
            for block in blocks {
                let indent = if first { 0.0 } else { CONTINUATION_INDENT };
                if !pieces.is_empty() && indent + width + GAP + block.width > max_width {
                    lines.push(Line::new(Block::row(pieces, None), depth, indent, first));
                    first = false;
                    pieces = Vec::new();
                    width = 0.0;
                }
                width += if pieces.is_empty() { 0.0 } else { GAP } + block.width;
                pieces.push(block);
            }
            let indent = if first { 0.0 } else { CONTINUATION_INDENT };
            lines.push(Line::new(Block::row(pieces, None), depth, indent, first));
        }
        sheet.document(lines)
    }

    /// Draw the term of the forest at `node`: its sign, cofactor and factors,
    /// as one term of [`Self::expansion_svg`].
    ///
    /// # Errors
    ///
    /// Fails as [`Self::operation_graph_svg`] does.
    ///
    /// # Panics
    ///
    /// Panics if `node` is not a forest of these forests or `expansion` is not
    /// theirs.
    ///
    /// See [`Self::expansion_svg`] for an example of the expansion.
    pub fn term_svg(
        &self,
        node: NodeIndex,
        diagram: &FeynmanDiagram,
        expansion: &Expansion,
        config: &Config,
        stamps: StampOptions,
        op: &str,
    ) -> Result<String, UvDisplayError> {
        let mut sheet = StampSheet::new(diagram, config, stamps)?;
        let term = expansion
            .term(node)
            .expect("the forest belongs to these forests");
        let row = sheet.term(self.wood(), expansion, term, op, true)?;
        let depth = self.operation(node).depth();
        sheet.document(vec![Line::new(row, depth, 0.0, true)])
    }
}

/// Hover text of a term: `K_{S|G} = k5 · depth 2, sign +`.
fn term_title(expansion: &Expansion, term: &ForestTerm, op: &str) -> String {
    format!(
        "{} = k{} · depth {}, sign {}",
        expansion.notation(term, op),
        term.k,
        term.factors.iter().map(factor_depth).sum::<usize>(),
        if term.sign > 0 { "+" } else { "−" },
    )
}

/// Number of operations in a factor and its nested factors.
fn factor_depth(factor: &TermFactor) -> usize {
    1 + factor.inner.iter().map(factor_depth).sum::<usize>()
}

impl Wood<LoopMomentumBasis> {
    /// Draw the wood as its Hasse diagram: one stamp per spinney with the
    /// spinney in blue, ranked by loop count with the largest spinneys on top and
    /// the bare diagram (the empty root) at the bottom. Lines join each spinney
    /// to the spinneys it covers. Each box carries the spinney's hover text.
    ///
    /// `config` and `stamps` are as in [`Forests::operation_graph_svg`].
    ///
    /// # Errors
    ///
    /// Fails as [`Forests::operation_graph_svg`] does.
    ///
    /// ```
    /// use feynkit_graph::FeynmanDiagram;
    /// use feynkit_model::Model;
    /// use feynkit_uv::{RenormalizationPrescriptionSettings, StampOptions, UvGraph, Wood};
    /// use linnest::svg::Config;
    ///
    /// let basketball = FeynmanDiagram::from_dot(
    ///     Model::phi4(),
    ///     "digraph basketball { edge [particle=phi]; a -> b; a -> b; a -> b; a -> b; }",
    /// )?;
    /// let spinneys = basketball.classify_spinneys(
    ///     &basketball.hedge_graph().full_filter(),
    ///     basketball.model(),
    ///     &RenormalizationPrescriptionSettings::default(),
    ///     basketball.loop_momentum_basis(),
    /// )?;
    /// let wood = Wood::from_spinneys(spinneys, &basketball, basketball.loop_momentum_basis())?;
    /// let svg = wood.hasse_svg(&basketball, &Config::default(), StampOptions::default())?;
    /// assert_eq!(svg.matches("data-feynkit-spinney=").count(), 12);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn hasse_svg(
        &self,
        diagram: &FeynmanDiagram,
        config: &Config,
        stamps: StampOptions,
    ) -> Result<String, UvDisplayError> {
        let mut sheet = StampSheet::new(diagram, config, stamps)?;
        let domain = self.divergent_cover(diagram);
        let nodes = self.nodes_by_edges(diagram);
        let index: BTreeMap<NodeIndex, usize> = nodes
            .iter()
            .enumerate()
            .map(|(position, node)| (*node, position))
            .collect();
        let most = self.max_loop_count();
        let mut boxes = Vec::with_capacity(nodes.len());
        for (position, node) in nodes.iter().enumerate() {
            let spinney = self.spinney(*node);
            let block = if spinney.is_empty() {
                sheet.stamp(None, None)?
            } else {
                sheet.stamp(Some(spinney.filter()), None)?
            };
            let (title, fill) = if spinney.is_empty() {
                ("unsubtracted diagram".to_owned(), ROOT_FILL)
            } else {
                (
                    format!(
                        "{} · {}",
                        spinney.label(diagram),
                        spinney.description(diagram)
                    ),
                    if *spinney.filter() == domain {
                        FULL_FILL
                    } else {
                        "none"
                    },
                )
            };
            boxes.push(GraphBox {
                block,
                title,
                rank: most - spinney.loop_count() + 1,
                fill,
                attribute: format!("data-feynkit-spinney=\"{position}\""),
            });
        }
        let links = self
            .cover_relation()
            .into_iter()
            .map(|(smaller, larger)| (index[&larger], index[&smaller]))
            .collect::<Vec<_>>();
        sheet.render(boxes, &links, Direction::Down)
    }
}

impl Counterterm<LoopMomentumBasis> {
    /// The contracted spinneys as regions, outermost first. Each outermost
    /// spinney keeps the layout, fades its lines inside a red hull and carries
    /// one circled cross: open for MS-bar, filled on shell and in a dashed ring
    /// while the scheme is unspecified. Its hover lists the contracted chain,
    /// inner to outer. Nested spinneys show only as hulls inside it, violet,
    /// then blue, then amber by nesting depth.
    ///
    /// ```
    /// use feynkit_graph::{FeynmanDiagram, HighlightStyle};
    /// use feynkit_model::Model;
    /// use feynkit_uv::{Counterterm, CountertermScheme, RenormalizationPrescriptionSettings, UvGraph};
    ///
    /// let sunrise = FeynmanDiagram::from_dot(
    ///     Model::phi4(),
    ///     "digraph sunrise { edge [particle=phi]; ext [style=invis]; \
    ///          ext -> a; a -> b; a -> b; a -> b; b -> ext; }",
    /// )?;
    /// let mut spinneys = sunrise.classify_spinneys(
    ///     &sunrise.hedge_graph().full_filter(),
    ///     sunrise.model(),
    ///     &RenormalizationPrescriptionSettings::default(),
    ///     sunrise.loop_momentum_basis(),
    /// )?;
    /// spinneys.retain(|spinney| !spinney.is_empty());
    /// spinneys.sort_by_key(|spinney| spinney.loop_count());
    /// let chain = [spinneys[3].clone(), spinneys[0].clone()];
    /// let counterterm =
    ///     Counterterm::new(chain, CountertermScheme::OnShell, &sunrise, sunrise.loop_momentum_basis())?;
    /// let regions = counterterm.regions(&sunrise);
    /// // One crossed region for the full graph, the bubble only as a nested hull.
    /// assert!(regions[0].details.contracted && !regions[1].details.contracted);
    /// assert!(regions[0].details.label.ends_with(" ⊂ γ{0,1,2}"));
    /// assert_eq!(regions[1].style, HighlightStyle::Middle);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn regions(&self, diagram: &FeynmanDiagram) -> Vec<UvRegion> {
        let spinneys: Vec<_> = self.contracted().iter().collect();
        let inside = |inner: &Spinney<LoopMomentumBasis>, outer: &Spinney<LoopMomentumBasis>| {
            inner.filter() != outer.filter() && outer.subgraph().includes(inner.subgraph())
        };
        let depth = |spinney: &Spinney<LoopMomentumBasis>| {
            spinneys
                .iter()
                .filter(|other| inside(spinney, other))
                .count()
        };
        let mut order: Vec<_> = spinneys.to_vec();
        order.sort_by_key(|spinney| (depth(spinney), spinney.filter().clone()));
        order
            .into_iter()
            .map(|spinney| {
                let mut region = spinney.region(diagram);
                match depth(spinney) {
                    0 => {
                        // The nested spinneys, inner to outer; siblings inside
                        // one parent are grouped in parentheses.
                        fn tree(
                            spinney: &Spinney<LoopMomentumBasis>,
                            spinneys: &[&Spinney<LoopMomentumBasis>],
                            diagram: &FeynmanDiagram,
                            inside: &dyn Fn(
                                &Spinney<LoopMomentumBasis>,
                                &Spinney<LoopMomentumBasis>,
                            ) -> bool,
                        ) -> String {
                            let children: Vec<String> = spinneys
                                .iter()
                                .filter(|child| {
                                    inside(child, spinney)
                                        && !spinneys.iter().any(|between| {
                                            inside(child, between) && inside(between, spinney)
                                        })
                                })
                                .map(|child| tree(child, spinneys, diagram, inside))
                                .collect();
                            let label = spinney.label(diagram);
                            match children.as_slice() {
                                [] => label,
                                [only] => format!("{only} ⊂ {label}"),
                                many => format!("({}) ⊂ {label}", many.join(", ")),
                            }
                        }
                        region.details.label = tree(spinney, &spinneys, diagram, &inside);
                        region.style = HighlightStyle::Treated;
                        region.details.contracted = true;
                        region.details.scheme = self.scheme().into();
                        // Linnest appends a specified scheme to the hover text itself.
                        region.details.description.push_str("; contracted");
                    }
                    1 => region.style = HighlightStyle::Middle,
                    2 => region.style = HighlightStyle::ActedOn,
                    _ => region.style = HighlightStyle::Amber,
                }
                region
            })
            .collect()
    }
}

impl Wood<LoopMomentumBasis> {
    /// The union of every spinney: the lines any forest can cover.
    fn divergent_cover(&self, diagram: &FeynmanDiagram) -> SuBitGraph {
        let mut cover = diagram.underlying().empty_subgraph::<SuBitGraph>();
        for (_, spinney) in self.spinneys() {
            cover.union_with(spinney.filter());
        }
        cover
    }
}

/// The union of the spinneys at `nodes`, or `None` for no spinney.
fn union_of(wood: &Wood<LoopMomentumBasis>, nodes: &[NodeIndex]) -> Option<SuBitGraph> {
    let (first, rest) = nodes.split_first()?;
    let mut union = wood.spinney(*first).filter().clone();
    for node in rest {
        union.union_with(wood.spinney(*node).filter());
    }
    Some(union)
}

/// Rank values by size: equal sizes share a rank, smaller sizes come first.
fn ranks_by_size(sizes: impl Iterator<Item = usize>) -> Vec<usize> {
    let sizes: Vec<usize> = sizes.collect();
    let distinct: BTreeSet<usize> = sizes.iter().copied().collect();
    let rank: BTreeMap<usize, usize> = distinct
        .into_iter()
        .enumerate()
        .map(|(rank, size)| (size, rank))
        .collect();
    sizes.iter().map(|size| rank[size]).collect()
}

/// A piece of a drawn formula, around its centre line `y = 0`: it spans
/// `[0, width] × [-ascent, descent]`, in points.
#[derive(Clone, Debug)]
struct Block {
    width: f64,
    ascent: f64,
    descent: f64,
    svg: String,
}

/// Formula text size, the gap between juxtaposed pieces and the padding inside
/// brackets, in points.
const TEXT_SIZE: f64 = 9.0;
const GAP: f64 = 1.6;
const BRACKET_PAD: f64 = 0.8;
/// Extra space in front of an operator symbol, in points.
const OP_SPACE: f64 = 2.4;
/// Average advance of one character, as a fraction of the text size.
const ADVANCE: f64 = 0.56;
/// Math font stack of the formula text.
const MATH_FONT: &str = "STIX Two Math, Cambria Math, Times New Roman, serif";

impl Block {
    /// A stamp of size `width × height`, centred on the line.
    fn stamp(id: &str, width: f64, height: f64) -> Self {
        Self {
            width,
            ascent: height / 2.0,
            descent: height / 2.0,
            svg: format!(
                r##"<use href="#{id}" y="{}" width="{}" height="{}"/>"##,
                number(-height / 2.0),
                number(width),
                number(height)
            ),
        }
    }

    /// Math text centred on the line. `^{…}` and `_{…}` become superscripts and
    /// subscripts; everything else is drawn as written.
    fn text(text: &str) -> Self {
        let mut spans = String::new();
        let mut width = 0.0;
        let mut rest = text;
        while !rest.is_empty() {
            let script = ["^{", "_{"]
                .into_iter()
                .filter_map(|marker| rest.find(marker).map(|at| (at, marker)))
                .min();
            let (plain, tail) = match script {
                Some((at, _)) => rest.split_at(at),
                None => (rest, ""),
            };
            if !plain.is_empty() {
                let _ = write!(spans, "<tspan>{}</tspan>", xml_text(plain));
                width += ADVANCE * TEXT_SIZE * plain.chars().count() as f64;
            }
            let Some((_, marker)) = script else { break };
            let body = &tail[2..];
            let end = body.find('}').unwrap_or(body.len());
            let shift = if marker == "^{" { "0.45em" } else { "-0.3em" };
            let _ = write!(
                spans,
                r#"<tspan font-size="70%" baseline-shift="{shift}">{}</tspan>"#,
                xml_text(&body[..end])
            );
            width += 0.7 * ADVANCE * TEXT_SIZE * body[..end].chars().count() as f64;
            rest = body.get(end + 1..).unwrap_or("");
        }
        Self {
            width,
            ascent: 0.7 * TEXT_SIZE,
            descent: 0.3 * TEXT_SIZE,
            svg: format!(
                r#"<text y="{}" font-family="{MATH_FONT}" font-size="{TEXT_SIZE}" fill="{INK}">{spans}</text>"#,
                number(0.35 * TEXT_SIZE)
            ),
        }
    }

    /// `pieces` side by side, separated by `separator` (or a gap).
    fn row(pieces: Vec<Block>, separator: Option<&str>) -> Self {
        let mut row = Self {
            width: 0.0,
            ascent: 0.0,
            descent: 0.0,
            svg: String::new(),
        };
        for (index, piece) in pieces.into_iter().enumerate() {
            if index > 0 {
                row.width += GAP;
                if let Some(separator) = separator {
                    row.push(Self::text(separator));
                    row.width += GAP;
                }
            }
            row.push(piece);
        }
        row
    }

    fn push(&mut self, piece: Block) {
        let _ = write!(
            self.svg,
            r#"<g transform="translate({} 0)">{}</g>"#,
            number(self.width),
            piece.svg
        );
        self.width += piece.width;
        self.ascent = self.ascent.max(piece.ascent);
        self.descent = self.descent.max(piece.descent);
    }

    /// `content` in square brackets drawn to its height.
    fn bracketed(content: Block) -> Self {
        const ARM: f64 = 1.5;
        let (top, bottom) = (content.ascent + BRACKET_PAD, content.descent + BRACKET_PAD);
        let bracket = |x: f64, arm: f64| {
            format!(
                r#"<path d="M{} {}H{}V{}H{}" fill="none" stroke="{INK}" stroke-width="0.7"/>"#,
                number(x + arm),
                number(-top),
                number(x),
                number(bottom),
                number(x + arm),
            )
        };
        let mut block = Self {
            width: 0.0,
            ascent: top,
            descent: bottom,
            svg: bracket(0.0, ARM),
        };
        block.width = ARM + BRACKET_PAD;
        let width = content.width;
        let _ = write!(
            block.svg,
            r#"<g transform="translate({} 0)">{}</g>"#,
            number(block.width),
            content.svg
        );
        block.width += width + BRACKET_PAD + ARM;
        block.svg.push_str(&bracket(block.width, -ARM));
        block
    }

    /// `Op[stamp · inner₁ · inner₂ ⋯]`.
    fn factor(op: &str, stamp: Block, inner: Vec<Block>) -> Self {
        let mut content = vec![stamp];
        content.extend(inner);
        let mut factor = Self::row(
            vec![
                Self::text(op),
                Self::bracketed(Self::row(content, Some("·"))),
            ],
            None,
        );
        factor.svg = format!(
            r#"<g transform="translate({} 0)">{}</g>"#,
            number(OP_SPACE),
            factor.svg
        );
        factor.width += OP_SPACE;
        factor
    }
}

fn number(value: f64) -> String {
    let rounded = (value * 1e4).round() / 1e4;
    if rounded == 0.0 {
        "0".into()
    } else {
        format!("{rounded}")
    }
}

/// Width of an expansion line before it wraps, in stamp sizes.
const LINE_STAMPS: f64 = 20.0;
/// Indent of a wrapped line under the first line of its depth, in points.
const CONTINUATION_INDENT: f64 = 18.0;
/// Vertical space between depths, and between wrapped lines of one depth.
const DEPTH_GAP: f64 = 12.0;
const LINE_GAP: f64 = 3.0;

/// One line of the expansion: terms of one depth, the first line of the depth
/// or a wrapped continuation indented under it.
struct Line {
    block: Block,
    depth: usize,
    indent: f64,
    starts_depth: bool,
}

impl Line {
    fn new(block: Block, depth: usize, indent: f64, starts_depth: bool) -> Self {
        Self {
            block,
            depth,
            indent,
            starts_depth,
        }
    }
}

/// Stamps of one diagram, rendered once each and shared through `<use>`.
struct StampSheet<'a> {
    diagram: &'a FeynmanDiagram,
    config: &'a Config,
    scale: Option<f64>,
    options: StampOptions,
    stamps: BTreeMap<(Option<SuBitGraph>, Option<SuBitGraph>), Block>,
    defs: String,
}

/// A box of a stamp graph: its label, hover title, layout rank and fill.
struct GraphBox {
    block: Block,
    title: String,
    rank: usize,
    fill: &'static str,
    attribute: String,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Direction {
    Right,
    Down,
}

impl<'a> StampSheet<'a> {
    fn new(
        diagram: &'a FeynmanDiagram,
        config: &'a Config,
        options: StampOptions,
    ) -> Result<Self, UvDisplayError> {
        if !(options.size.is_finite() && options.size > 0.0) {
            return Err(UvDisplayError::StampSize(options.size));
        }
        if !(options.weight.is_finite() && options.weight > 0.0) {
            return Err(UvDisplayError::StampWeight(options.weight));
        }
        Ok(Self {
            diagram,
            config,
            scale: None,
            options,
            stamps: BTreeMap::new(),
            defs: String::new(),
        })
    }

    /// The stamp with `current` in blue and `given` in red, drawn over the
    /// muted diagram; neither gives the bare stamp.
    fn stamp(
        &mut self,
        current: Option<&SuBitGraph>,
        given: Option<&SuBitGraph>,
    ) -> Result<Block, UvDisplayError> {
        if self.scale.is_none() && (current.is_some() || given.is_some()) {
            self.stamp(None, None)?;
        }
        let key = (current.cloned(), given.cloned());
        if let Some(stamp) = self.stamps.get(&key) {
            return Ok(stamp.clone());
        }
        let region = |hedges: &SuBitGraph, style| UvRegion {
            hedges: hedges.clone(),
            style,
            details: HighlightDetails::default(),
        };
        let mut regions: Vec<UvRegion> = [
            (current, HighlightStyle::ActedOn),
            (given, HighlightStyle::Treated),
        ]
        .into_iter()
        .filter_map(|(hedges, style)| hedges.map(|hedges| region(hedges, style)))
        .collect();
        if regions.is_empty() {
            // An empty region still selects recolouring, which mutes every line
            // to the plain carriers of the other stamps.
            regions.push(region(
                &self.diagram.underlying().empty_subgraph(),
                HighlightStyle::Treated,
            ));
        }
        let options = SceneOptions {
            show_particle: false,
            overlap: HighlightOverlap::Recolour,
            stamp_weight: self.options.weight,
            ..SceneOptions::default()
        };
        let mut scene =
            self.diagram
                .to_scene_with_highlights(&UvRegion::layers(&regions), None, &options)?;
        self.config
            .apply(&mut scene)
            .map_err(UvDisplayError::Render)?;
        scene.title = None;
        let svg = scene
            .render_static(&Typeset::default())
            .map_err(UvDisplayError::Render)?;
        let (width, height, inner) = svg_parts(&svg)?;
        let scale = *self
            .scale
            .get_or_insert_with(|| self.options.size / width.max(height));
        // Crop Linnest's blank border so brackets and neighbours sit on the ink;
        // the scale still follows the uncropped bare stamp.
        let inset = CANVAS_INSET.min(width / 4.0).min(height / 4.0);
        let (width, height) = (width - 2.0 * inset, height - 2.0 * inset);
        let mut hasher = DefaultHasher::new();
        inner.hash(&mut hasher);
        let id = format!("feynkit-uv-stamp-{:016x}", hasher.finish());
        let _ = write!(
            self.defs,
            r#"<symbol id="{id}" viewBox="{} {} {} {}">{inner}</symbol>"#,
            number(inset),
            number(inset),
            number(width),
            number(height)
        );
        let stamp = Block::stamp(&id, width * scale, height * scale);
        self.stamps.insert(key, stamp.clone());
        Ok(stamp)
    }

    /// The factors `Op[m(S|G) · inner]` of a term, side by side with `·`.
    fn factors(
        &mut self,
        wood: &Wood<LoopMomentumBasis>,
        factors: &[TermFactor],
        op: &str,
    ) -> Result<Block, UvDisplayError> {
        let mut blocks = Vec::with_capacity(factors.len());
        for factor in factors {
            let given = union_of(wood, &factor.given);
            let stamp = self.stamp(Some(wood.spinney(factor.spinney).filter()), given.as_ref())?;
            let inner = factor
                .inner
                .iter()
                .map(|inner| self.factors(wood, slice::from_ref(inner), op))
                .collect::<Result<Vec<_>, _>>()?;
            blocks.push(Block::factor(op, stamp, inner));
        }
        Ok(Block::row(blocks, Some("·")))
    }

    /// `± m(Γ|γ) · factors`, with the hover text and `k` of the term. The sign is
    /// omitted for a positive term unless `signed`.
    fn term(
        &mut self,
        wood: &Wood<LoopMomentumBasis>,
        expansion: &Expansion,
        term: &ForestTerm,
        op: &str,
        signed: bool,
    ) -> Result<Block, UvDisplayError> {
        let whole = self.diagram.internal_subgraph();
        let treated = union_of(wood, &term.outermost);
        let mut pieces = Vec::new();
        if signed || term.sign < 0 {
            pieces.push(Block::text(if term.sign > 0 { "+" } else { "−" }));
        }
        let cofactor = self.stamp(Some(&whole), treated.as_ref())?;
        if term.factors.is_empty() {
            pieces.push(cofactor);
        } else {
            pieces.push(Block::row(
                vec![cofactor, self.factors(wood, &term.factors, op)?],
                Some("·"),
            ));
        }
        let mut block = Block::row(pieces, None);
        block.svg = format!(
            r#"<g data-feynkit-term="{}"><title>{}</title>{}</g>"#,
            term.k,
            xml_text(&term_title(expansion, term, op)),
            block.svg
        );
        block.width += 2.0 * GAP;
        Ok(block)
    }

    /// A standalone SVG document of `rows`, stacked and left aligned.
    fn document(self, lines: Vec<Line>) -> Result<String, UvDisplayError> {
        const MARGIN: f64 = 6.0;
        let width = lines
            .iter()
            .map(|line| line.indent + line.block.width)
            .fold(0.0, f64::max)
            + 2.0 * MARGIN;
        let mut body = String::new();
        let mut y = MARGIN;
        for (index, line) in lines.iter().enumerate() {
            if index > 0 {
                y += if line.starts_depth {
                    DEPTH_GAP
                } else {
                    LINE_GAP
                };
            }
            y += line.block.ascent;
            let _ = write!(
                body,
                r#"<g data-feynkit-depth="{}" data-feynkit-continued="{}" transform="translate({} {})">{}</g>"#,
                line.depth,
                !line.starts_depth,
                number(MARGIN + line.indent),
                number(y),
                line.block.svg
            );
            y += line.block.descent;
        }
        let height = y + MARGIN;
        Ok(format!(
            r#"<svg viewBox="0 0 {w} {h}" width="{w}pt" height="{h}pt" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">{body}<defs>{}</defs></svg>"#,
            self.defs,
            w = number(width),
            h = number(height),
        ))
    }

    /// Lay the boxes out by rank with Linnest's layered layout and draw them.
    fn render(
        self,
        boxes: Vec<GraphBox>,
        links: &[(usize, usize)],
        direction: Direction,
    ) -> Result<String, UvDisplayError> {
        // Compass ports make the layered router join connectors at one point of
        // a box side; Hasse lines keep their own attachment points instead.
        let ((exit, exit_compass), (entry, entry_compass)) = match direction {
            Direction::Right => (("east", Some("e")), ("west", Some("w"))),
            Direction::Down => (("south", None), ("north", None)),
        };
        let mut scene = Scene {
            highlights: Default::default(),
            graph: TypstGraphSpec {
                name: None,
                data: None,
                statements: BTreeMap::new(),
                default_edge_statements: BTreeMap::new(),
                default_node_statements: BTreeMap::new(),
                nodes: Vec::with_capacity(boxes.len()),
                edges: Vec::with_capacity(links.len()),
            },
            nodes: Vec::with_capacity(boxes.len()),
            edges: Vec::with_capacity(links.len()),
            preamble: String::new(),
            title: None,
            pages: Vec::with_capacity(boxes.len()),
            layout: layered_layout(direction),
            layout_passes: Vec::new(),
            layout_edges: None,
            label_feedback: false,
        };
        let mut pages = Vec::with_capacity(boxes.len());
        for (position, graph_box) in boxes.into_iter().enumerate() {
            let title = xml_text(&graph_box.title);
            let block = &graph_box.block;
            let height = block.ascent + block.descent;
            pages.push(LabelPage {
                metrics: LabelMetrics {
                    width: block.width,
                    height,
                    baseline: height,
                },
                body: format!(
                    r#"<g {} transform="translate(0 {})"><title>{title}</title>{}</g>"#,
                    graph_box.attribute,
                    number(block.ascent),
                    block.svg
                ),
            });
            scene.pages.push(LabelSource::text(&graph_box.title));
            scene.graph.nodes.push(TypstNodeSpec {
                name: Some(format!("n{position}")),
                index: Some(position),
                data: None,
                pos: None,
                statements: BTreeMap::from([(
                    "layout-rank".to_owned(),
                    graph_box.rank.to_string(),
                )]),
            });
            scene.nodes.push(NodeDrawing {
                highlight_ring: Vec::new(),
                radius: 0.3,
                label_padding: Some(0.3),
                label: Some(position),
                rectangular: true,
                fill: graph_box.fill.to_owned(),
                stroke: Stroke {
                    paint: BOX_STROKE.to_owned(),
                    width: 0.4,
                    dash: Dash::Solid,
                    round_cap: false,
                },
                details: Details::default(),
            });
        }
        for (position, &(source, sink)) in links.iter().enumerate() {
            let endpoint = |node: usize, compass: Option<&str>| TypstEndpointSpec {
                node,
                statement: None,
                id: None,
                data: None,
                port_label: None,
                compass: compass.map(str::to_owned),
                in_subgraph: false,
                route_points: Vec::new(),
            };
            scene.graph.edges.push(TypstEdgeSpec {
                name: None,
                source: Some(endpoint(source, exit_compass)),
                sink: Some(endpoint(sink, entry_compass)),
                data: None,
                orientation: None,
                flow: None,
                id: Some(position),
                pos: None,
                statements: BTreeMap::from([
                    ("source-route-exit".to_owned(), exit.to_owned()),
                    ("sink-route-exit".to_owned(), entry.to_owned()),
                ]),
            });
            scene.edges.push(EdgeDrawing {
                stroke: Stroke {
                    paint: CONNECTOR.to_owned(),
                    width: 0.6,
                    dash: Dash::Solid,
                    round_cap: true,
                },
                highlight_strokes: Vec::new(),
                pattern: None,
                flow: None,
                momentum: false,
                momentum_offset: EdgeDrawing::DEFAULT_MOMENTUM_OFFSET,
                label: None,
                details: Details::default(),
            });
        }
        let typeset = Typeset {
            title: None,
            pages,
            defs: self.defs,
        };
        scene
            .render_static(&typeset)
            .map_err(UvDisplayError::Render)
    }
}

/// Layered layout options of the `tree-auto-root` operation graphs: wide gaps
/// between ranks for the connectors, narrow gaps between boxes of one rank.
/// `tree-dx` is horizontal and `tree-dy` vertical, in drawing units.
fn layered_layout(direction: Direction) -> Map<String, Value> {
    let (name, dx, dy) = match direction {
        Direction::Right => ("right", 3.0, 0.42),
        Direction::Down => ("down", 0.42, 2.2),
    };
    let options = json!({
        "layout-algo": "dot",
        "layout-direction": name,
        "rank-align": "center",
        "tree-dx": dx,
        "tree-dy": dy,
        "route-label-width-cap": 0.0,
        "label-steps": 0,
    });
    match options {
        Value::Object(map) => map,
        _ => unreachable!("the layered layout options are an object"),
    }
}

/// Width, height and content of a complete SVG document from Linnest.
fn svg_parts(svg: &str) -> Result<(f64, f64, &str), UvDisplayError> {
    let malformed = || UvDisplayError::Render("stamp is not a complete SVG document".into());
    let open = svg.find('>').ok_or_else(malformed)?;
    let close = svg.rfind("</svg>").ok_or_else(malformed)?;
    let header = &svg[..open];
    let view_box = header
        .split("viewBox=\"")
        .nth(1)
        .and_then(|rest| rest.split('"').next())
        .ok_or_else(malformed)?;
    let values: Vec<f64> = view_box
        .split_whitespace()
        .map(str::parse)
        .collect::<Result<_, _>>()
        .map_err(|_| malformed())?;
    let [_, _, width, height] = values[..] else {
        return Err(malformed());
    };
    Ok((width, height, &svg[open + 1..close]))
}

/// Escape plain text for SVG character data. Control characters, which XML
/// cannot carry even escaped, become U+FFFD.
fn xml_text(text: &str) -> String {
    let mut escaped = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '&' => escaped.push_str("&amp;"),
            '<' => escaped.push_str("&lt;"),
            '>' => escaped.push_str("&gt;"),
            '"' => escaped.push_str("&quot;"),
            '\t' | '\n' | '\r' => escaped.push(c),
            c if c.is_control() || c == '\u{FFFE}' || c == '\u{FFFF}' => escaped.push('\u{FFFD}'),
            c => escaped.push(c),
        }
    }
    escaped
}
