//! Label pages typeset by the host, their placement candidates, and the
//! annotation search `draw.typ` runs. Candidates and the search exchange the
//! same records as the Typst renderer, so both route labels identically.
use std::collections::BTreeMap;

use kurbo::{BezPath, Point};
use serde::{de::DeserializeOwned, Deserialize, Serialize};

use super::MOMENTUM_ARROW;
use super::{marks, UNIT};
use crate::annotation::MomentumCandidates;
pub(super) use crate::annotation::{
    carrier_candidates, endpoint_candidates, Anchor, Bounds, Candidate, Places,
};
pub(super) use crate::label_placement::Leader;
use crate::label_placement::{
    label_overlaps, PlacedLabel, Point as PlanePoint, StrokeLine, LABEL_PADDING, OBSTACLE_PADDING,
    PAIR_PADDING, STROKE_ACCURACY,
};

/// One typeset page: its size and SVG body in points. Label pages also carry
/// the content measurements CeTZ uses (width, height, line baseline, line
/// bounds) and the surrounding text size.
#[derive(Clone, Debug)]
pub struct LabelPage {
    pub width: f64,
    pub height: f64,
    pub body: String,
    pub metrics: Option<[f64; 4]>,
    pub text_size: Option<f64>,
}

/// A scene's typeset title and label pages, with their merged glyph definitions.
#[derive(Clone, Debug, Default)]
pub struct Typeset {
    pub title: Option<LabelPage>,
    pub pages: Vec<LabelPage>,
    pub defs: String,
}

const METRICS_LINK: &str = "linnest-metrics:";

impl LabelPage {
    /// A page holding `content` in CeTZ's content box: cap height down to its
    /// glyph bounds. An invisible link reports the measurements CeTZ uses.
    pub(super) fn measured(content: &str) -> String {
        format!(
            "#context {{ let body = {content}; \
             let content = text(top-edge: \"cap-height\", bottom-edge: \"baseline\", body); \
             let size = measure(content); \
             let line(edge) = measure(text(top-edge: \"cap-height\", bottom-edge: edge, [ #show linebreak: [ ]; #body])).height; \
             let (base, bounds) = (line(\"baseline\"), line(\"bounds\")); \
             block(width: size.width, height: size.height + bounds - base, {{ place(top + left, content); \
             place(top + left, link(\"{METRICS_LINK}\" + (size.width, size.height, base, bounds, text.size).map(v => repr(v.pt())).join(\",\"), box(width: 1pt, height: 1pt))) }}) }}"
        )
    }

    /// Read one SVG page, collecting its glyph symbols into `defs`.
    fn read(svg: &str, defs: &mut BTreeMap<String, String>) -> Result<Self, String> {
        let invalid = || "invalid typeset label page".to_owned();
        let view_box = svg.split("viewBox=\"0 0 ").nth(1).ok_or_else(invalid)?;
        let mut size = view_box
            .split('"')
            .next()
            .ok_or_else(invalid)?
            .split(' ')
            .map(str::parse::<f64>);
        let (Some(Ok(width)), Some(Ok(height))) = (size.next(), size.next()) else {
            return Err(invalid());
        };
        let open = svg.find('>').ok_or_else(invalid)? + 1;
        let close = svg
            .find("<defs")
            .or_else(|| svg.rfind("</svg>"))
            .ok_or_else(invalid)?;
        let mut body = svg[open..close].to_owned();
        let metrics = match body.find(&format!("<a href=\"{METRICS_LINK}")) {
            None => None,
            Some(start) => {
                let end = start + body[start..].find("</a>").ok_or_else(invalid)? + "</a>".len();
                let link: String = body.drain(start..end).collect();
                let mut values = link
                    .split(METRICS_LINK)
                    .nth(1)
                    .and_then(|values| values.split('"').next())
                    .ok_or_else(invalid)?
                    .split(',')
                    .map(str::parse::<f64>);
                let mut metrics = [0.0; 5];
                for metric in &mut metrics {
                    *metric = values.next().and_then(Result::ok).ok_or_else(invalid)?;
                }
                Some(metrics)
            }
        };
        let (metrics, text_size) = match metrics {
            Some([width, height, base, bounds, size]) => {
                (Some([width, height, base, bounds]), Some(size))
            }
            None => (None, None),
        };
        if let (Some(start), Some(end)) = (svg.find("<defs"), svg.rfind("</defs>")) {
            let inner = &svg[start + svg[start..].find('>').ok_or_else(invalid)? + 1..end];
            for symbol in inner.split("<symbol").filter(|symbol| !symbol.is_empty()) {
                let id = symbol
                    .split("id=\"")
                    .nth(1)
                    .and_then(|id| id.split('"').next())
                    .ok_or_else(invalid)?;
                defs.entry(id.to_owned())
                    .or_insert_with(|| format!("<symbol{symbol}"));
            }
        }
        Ok(Self {
            width,
            height,
            body,
            metrics,
            text_size,
        })
    }
}

impl Typeset {
    /// Read the SVG pages of a scene's label document, title first when present.
    pub(super) fn read(svgs: &[String], title: bool) -> Result<Self, String> {
        let mut defs = BTreeMap::new();
        let mut pages = svgs
            .iter()
            .map(|svg| LabelPage::read(svg, &mut defs))
            .collect::<Result<Vec<_>, _>>()?;
        let title = if title {
            if pages.is_empty() {
                return Err("typeset label document has no title page".to_owned());
            }
            Some(pages.remove(0))
        } else {
            None
        };
        Ok(Self {
            title,
            pages,
            defs: defs.into_values().collect(),
        })
    }

    /// The labels' line height in points: their largest surrounding text size.
    pub(super) fn line_height(&self) -> f64 {
        self.pages
            .iter()
            .filter_map(|page| page.text_size)
            .fold(0.0, f64::max)
    }
}

fn decode<T: DeserializeOwned>(bytes: &[u8]) -> Result<T, String> {
    ciborium::from_reader(bytes).map_err(|error| error.to_string())
}

fn cbor<T: Serialize>(value: &T) -> Result<Vec<u8>, String> {
    crate::graph_api::encode_cbor(value)
}

/// One label's candidates, with its momentum arrow carriers per side.
pub(super) struct Placement {
    pub edge: usize,
    pub page: usize,
    pub candidates: Vec<Candidate>,
    pub href: String,
}

/// Momentum candidates whose arrow windows are mark placements of the
/// drawing's batch, resolved with the marks' painted footprints.
pub(super) struct MarkedMomentum {
    proposals: MomentumCandidates,
    marks: Vec<usize>,
}

impl MarkedMomentum {
    pub(super) fn new(
        path: &BezPath,
        metrics: [f64; 4],
        side: f64,
        template: usize,
        batch: &mut marks::Batch,
    ) -> Result<Self, String> {
        let proposals = MomentumCandidates::new(
            path,
            content_corners(metrics, Anchor::Center),
            side,
            MOMENTUM_ARROW,
            Places::Every,
        )?;
        let marks = proposals
            .windows()
            .iter()
            .map(|window| batch.push(template, window, kurvst::marks::MarkStation::End, true))
            .collect();
        Ok(Self { proposals, marks })
    }

    pub(super) fn resolve(
        self,
        geometry: &[kurvst::marks::MarkGeometryOutput],
    ) -> Result<Vec<Candidate>, String> {
        // As `draw.typ` places them, labels clear the arrow's shaft by their
        // gap, and the search weighs the painted arrow, heads included.
        let arrows: Vec<_> = self
            .marks
            .iter()
            .zip(self.proposals.windows())
            .map(|(&index, window)| {
                (
                    window.clone(),
                    Bounds::painted(&geometry[index].footprint.path),
                )
            })
            .collect();
        let mut candidates = self.proposals.resolve(&arrows)?;
        for (candidate, &index) in candidates.iter_mut().zip(&self.marks) {
            candidate.arrow_mark = Some(index);
        }
        Ok(candidates)
    }
}

/// CeTZ `content` corners (NW, NE, SW, SE) of a typeset page relative to the
/// placement point.
pub(super) fn content_corners(
    [width, height, baseline, bounds]: [f64; 4],
    anchor: Anchor,
) -> [[f64; 3]; 4] {
    let measure = |pt: f64| (pt / UNIT).abs();
    anchor.corners([
        measure(width) / 2.0,
        (measure(height) + measure(bounds) - measure(baseline)).abs() / 2.0,
    ])
}

/// A painted stroke for collisions: its radius, its path, and the edge
/// that paints it.
pub(super) struct Stroke {
    pub radius: f64,
    pub path: BezPath,
    pub edge: usize,
}

#[derive(Serialize)]
struct PaintedStroke {
    radius: f64,
    segments: Vec<(Point2, bool, Vec<serde_json::Value>)>,
}

type Point2 = [f64; 2];

impl From<&Stroke> for PaintedStroke {
    /// CeTZ subpaths `(origin, closed, commands)` of the painted path; quadratics
    /// are raised to cubics.
    fn from(stroke: &Stroke) -> Self {
        use kurbo::{PathEl, QuadBez};
        use serde_json::json;
        let xy = |p: Point| [p.x, p.y];
        let mut segments: Vec<(Point2, bool, Vec<serde_json::Value>)> = Vec::new();
        let mut current = Point::ZERO;
        for &element in stroke.path.elements() {
            let command = match element {
                PathEl::MoveTo(p) => {
                    segments.push((xy(p), false, Vec::new()));
                    current = p;
                    continue;
                }
                PathEl::ClosePath => {
                    if let Some(last) = segments.last_mut() {
                        last.1 = true;
                    }
                    continue;
                }
                PathEl::LineTo(p) => json!(["l", xy(p)]),
                PathEl::QuadTo(c, p) => {
                    let raised = QuadBez::new(current, c, p).raise();
                    json!(["c", xy(raised.p1), xy(raised.p2), xy(p)])
                }
                PathEl::CurveTo(c1, c2, p) => json!(["c", xy(c1), xy(c2), xy(p)]),
            };
            current = element.end_point().unwrap_or(current);
            if let Some(last) = segments.last_mut() {
                last.2.push(command);
            }
        }
        Self {
            radius: stroke.radius,
            segments,
        }
    }
}

#[derive(Serialize)]
struct StrokeLinesSpec {
    strokes: Vec<PaintedStroke>,
    accuracy: f64,
}

#[derive(Serialize)]
#[serde(rename_all = "kebab-case")]
struct SearchCandidate<'a> {
    bounds: Bounds,
    cost: f64,
    relative: f64,
    corners: Vec<[f64; 2]>,
    arrow_bounds: &'a [Bounds],
    anchor: [f64; 2],
    normal: Option<[f64; 2]>,
    side: f64,
}

/// Leaders as `draw.typ` requests them: one line of text keeps leadered
/// labels apart, each painted line names its edge, and the drawing has no
/// canvas transform.
#[derive(Serialize)]
#[serde(rename_all = "kebab-case")]
struct LeaderSpec<'a> {
    line_height: f64,
    line_edges: &'a [usize],
}

#[derive(Serialize)]
struct SearchPlacement<'a> {
    edge: usize,
    candidates: Vec<SearchCandidate<'a>>,
}

#[derive(Serialize)]
#[serde(rename_all = "kebab-case")]
struct SearchGeometry<'a> {
    placements: Vec<SearchPlacement<'a>>,
    obstacles: &'a [Bounds],
    label_padding: f64,
    fixed_arrows: &'a [Bounds],
    coordinated: bool,
    edge_lines: &'a [StrokeLine],
    obstacle_padding: f64,
    leaders: Option<LeaderSpec<'a>>,
}

#[derive(Serialize)]
#[serde(rename_all = "kebab-case")]
struct SearchMath {
    temperatures: Vec<f64>,
    exponentials: Vec<(f64, f64)>,
    pair_padding: f64,
}

#[derive(Deserialize)]
struct ExponentialCheck {
    argument: f64,
    lower: Option<f64>,
    upper: Option<f64>,
}

#[derive(Deserialize)]
struct SearchResult {
    choices: Vec<usize>,
    missing: Vec<ExponentialCheck>,
    leaders: Vec<Option<Leader>>,
}

/// The chosen candidate per placement, the leaders of labels set off their
/// carrier, and the painted collision lines.
pub(super) struct Searched {
    pub choices: Vec<usize>,
    pub leaders: Vec<Option<Leader>>,
    pub lines: Vec<StrokeLine>,
}

impl Searched {
    /// The candidate drawn for a placement: its leader's, if it has one.
    pub fn drawn(&self, index: usize) -> usize {
        self.leaders
            .get(index)
            .and_then(Option::as_ref)
            .map_or(self.choices[index], |leader| leader.choice)
    }
}

/// `_relax-label-placements`: one annealing search over all labels, with its
/// exponential comparisons certified by host values as Typst does.
pub(super) fn search(
    placements: &[Placement],
    obstacles: &[Bounds],
    strokes: &[Stroke],
    fixed_arrows: &[Bounds],
    leaders: Option<f64>,
) -> Result<Searched, String> {
    // Each stroke flattens on its own to the lines it adds to all of them.
    let mut edge_lines = Vec::new();
    let mut line_edges = Vec::new();
    for stroke in strokes {
        let spec = StrokeLinesSpec {
            strokes: vec![PaintedStroke::from(stroke)],
            accuracy: STROKE_ACCURACY,
        };
        let lines: Vec<StrokeLine> =
            decode(&crate::label_placement::stroke_lines_bytes(&cbor(&spec)?)?)?;
        line_edges.extend(std::iter::repeat_n(stroke.edge, lines.len()));
        edge_lines.extend(lines);
    }
    if placements
        .iter()
        .all(|placement| placement.candidates.len() == 1)
    {
        return Ok(Searched {
            choices: vec![0; placements.len()],
            leaders: vec![None; placements.len()],
            lines: edge_lines,
        });
    }
    let geometry = cbor(&SearchGeometry {
        placements: placements
            .iter()
            .map(|placement| SearchPlacement {
                edge: placement.edge,
                candidates: placement
                    .candidates
                    .iter()
                    .map(|candidate| SearchCandidate {
                        bounds: candidate.bounds,
                        cost: candidate.cost,
                        relative: candidate.relative,
                        corners: candidate.corners.iter().map(|c| [c[0], c[1]]).collect(),
                        arrow_bounds: &candidate.arrow_bounds,
                        anchor: candidate.anchor,
                        normal: candidate.normal,
                        side: candidate.side,
                    })
                    .collect(),
            })
            .collect(),
        obstacles,
        label_padding: LABEL_PADDING,
        fixed_arrows,
        coordinated: true,
        edge_lines: &edge_lines,
        obstacle_padding: OBSTACLE_PADDING,
        leaders: leaders.map(|line_height| LeaderSpec {
            line_height,
            line_edges: &line_edges,
        }),
    })?;
    let mut math = SearchMath {
        temperatures: (0..84).map(|sweep| 0.15 * 0.9f64.powi(sweep)).collect(),
        exponentials: Vec::new(),
        pair_padding: PAIR_PADDING,
    };
    loop {
        let result: SearchResult = decode(&crate::label_placement::search_bytes(
            &geometry,
            &cbor(&math)?,
        )?)?;
        let mut verified = true;
        for check in &result.missing {
            let probability = check.argument.exp();
            math.exponentials.push((check.argument, probability));
            verified &= check.lower.is_none_or(|lower| lower < probability)
                && check.upper.is_none_or(|upper| upper >= probability);
        }
        if result.missing.is_empty() || verified {
            return Ok(Searched {
                choices: result.choices,
                leaders: if result.leaders.is_empty() {
                    vec![None; placements.len()]
                } else {
                    result.leaders
                },
                lines: edge_lines,
            });
        }
    }
}

/// Real (unpadded) overlaps of the chosen labels, fixed ones included, with
/// painted lines and other labels' text: (overlapping labels, line hits,
/// overlapping label pairs). The shared placement clips each line exactly
/// against the text box, stroke radius included, so a line merely touching
/// it, such as a zero-radius mark outline, counts as a hit. Labels set off
/// their carrier with a leader count where the search placed them on it, so
/// label feedback keeps making room for them there.
pub(super) fn overlaps(placements: &[Placement], searched: &Searched) -> (usize, usize, usize) {
    let labels: Vec<_> = placements
        .iter()
        .zip(&searched.choices)
        .map(|(placement, &choice)| {
            let candidate = &placement.candidates[choice];
            PlacedLabel {
                bounds: candidate.bounds,
                corners: Some(candidate.corners.map(|[x, y, _]| PlanePoint(x, y))),
                arrow_bounds: Vec::new(),
                own_lines: 0..0,
                leader: None,
                edge: None,
            }
        })
        .collect();
    let mut overlapping = 0;
    let mut hits = 0;
    let mut pairs = 0;
    for found in label_overlaps(&labels, &searched.lines, &[], 0.0, &[], &[]) {
        hits += found.lines.len();
        pairs += found.labels.len();
        overlapping += usize::from(!found.lines.is_empty() || !found.labels.is_empty());
    }
    (overlapping, hits, pairs / 2)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn straight(from: (f64, f64), to: (f64, f64)) -> BezPath {
        let mut path = BezPath::new();
        path.move_to(from);
        path.line_to(to);
        path
    }

    /// Dotted edges are dashed only when painted; their collision stroke is the
    /// whole path. A narrow label prefers the side a long dotted edge crosses
    /// between the 2-unit spacing of points sampled along it.
    #[test]
    fn labels_avoid_and_report_a_dotted_edge_crossing_their_preferred_side() {
        let carrier = straight((-4.0, 0.0), (4.0, 0.0));
        // 8.1pt wide text: 0.6 drawing units.
        let candidates = carrier_candidates(
            &carrier,
            content_corners([8.1, 13.5, 10.0, 13.5], Anchor::Center),
            0.15,
            1.0,
            Places::Every,
        )
        .unwrap();
        // The preferred text lies above the carrier, strictly between the
        // dotted edge's sample points x = -0.5 and x = 1.5.
        let preferred = candidates[0].bounds;
        assert!(preferred.left > -0.4 && preferred.right < 0.4 && preferred.bottom > 0.0);
        let middle = (preferred.bottom + preferred.top) / 2.0;
        let strokes = [
            Stroke {
                radius: 0.05,
                path: carrier.clone(),
                edge: 0,
            },
            Stroke {
                radius: 0.05,
                path: straight((-8.5, middle), (7.5, middle)),
                edge: 1,
            },
        ];
        let placements = [Placement {
            edge: 0,
            page: 0,
            candidates,
            href: String::new(),
        }];
        let searched = search(&placements, &[], &strokes, &[], Some(9.0 / UNIT)).unwrap();
        let chosen_index = searched.choices[0];
        let chosen = &placements[0].candidates[chosen_index];
        assert!(chosen.bounds.top < 0.0, "the label takes the clear side");
        assert!(searched.leaders.iter().all(Option::is_none));
        assert_eq!(overlaps(&placements, &searched), (0, 0, 0));
        // On its preferred side, the text crosses the dotted edge.
        let preferred = Searched {
            choices: vec![0],
            leaders: vec![None],
            lines: searched.lines,
        };
        assert_eq!(overlaps(&placements, &preferred), (1, 1, 0));
        // However short, a stroke reaching the chosen text is a hit: here a
        // zero-radius cut across its corner, about 0.14 units inside.
        let (right, top) = (chosen.bounds.right, chosen.bounds.top);
        let cut = Stroke {
            radius: 0.0,
            path: straight((right - 0.2, top + 0.1), (right + 0.1, top - 0.2)),
            edge: 1,
        };
        let cut = Searched {
            choices: vec![chosen_index],
            leaders: vec![None],
            lines: search(&placements, &[], &[cut], &[], Some(9.0 / UNIT))
                .unwrap()
                .lines,
        };
        assert_eq!(overlaps(&placements, &cut), (1, 1, 0));
    }
}
