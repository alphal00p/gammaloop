//! Label pages typeset by the host, their placement candidates, and the
//! annotation search `draw.typ` runs. Candidates and the search exchange the
//! same records as the Typst renderer, so both route labels identically.
use std::collections::BTreeMap;

use kurbo::{BezPath, CubicBez, Point};
use serde::{de::DeserializeOwned, Deserialize, Serialize};

use super::{curves, marks, UNIT};

/// Clearance between an external label and its free endpoint.
const EXTERNAL_GAP: f64 = 0.25;
/// Momentum arrows: offset beside the edge, visible length and ratio cap, label gap.
pub(super) const ARROW_OFFSET: f64 = 0.35;
pub(super) const ARROW_WINDOW: (f64, f64) = (1.4, 0.5);
const MOMENTUM_GAP: f64 = 0.2;

#[derive(Serialize, Deserialize, Clone, Copy, Debug, PartialEq)]
pub(super) struct Bounds {
    pub left: f64,
    pub right: f64,
    pub bottom: f64,
    pub top: f64,
}

impl Bounds {
    fn overlap(&self, other: &Self) -> f64 {
        let width = self.right.min(other.right) - self.left.max(other.left);
        let height = self.top.min(other.top) - self.bottom.max(other.bottom);
        width.max(0.0) * height.max(0.0)
    }

    fn distance(&self, point: [f64; 2]) -> f64 {
        let dx = (self.left - point[0]).max(0.0).max(point[0] - self.right);
        let dy = (self.bottom - point[1]).max(0.0).max(point[1] - self.top);
        dx.hypot(dy)
    }
}

/// One typeset page: its size and SVG body in points. Label pages also carry
/// the content measurements CeTZ uses (width, height, line baseline, line bounds).
#[derive(Clone, Debug)]
pub struct LabelPage {
    pub width: f64,
    pub height: f64,
    pub body: String,
    pub metrics: Option<[f64; 4]>,
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
             place(top + left, link(\"{METRICS_LINK}\" + (size.width, size.height, base, bounds).map(v => repr(v.pt())).join(\",\"), box(width: 1pt, height: 1pt))) }}) }}"
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
                let mut metrics = [0.0; 4];
                for metric in &mut metrics {
                    *metric = values.next().and_then(Result::ok).ok_or_else(invalid)?;
                }
                Some(metrics)
            }
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
}

fn cbor<T: Serialize>(value: &T) -> Result<Vec<u8>, String> {
    crate::graph_api::encode_cbor(value)
}

fn decode<T: DeserializeOwned>(bytes: &[u8]) -> Result<T, String> {
    ciborium::from_reader(bytes).map_err(|error| error.to_string())
}

#[derive(Serialize)]
#[serde(rename_all = "kebab-case")]
struct Cubic {
    start: [f64; 2],
    control_start: [f64; 2],
    control_end: [f64; 2],
    end: [f64; 2],
}

impl From<&CubicBez> for Cubic {
    fn from(cubic: &CubicBez) -> Self {
        let xy = |p: Point| [p.x, p.y];
        Self {
            start: xy(cubic.p0),
            control_start: xy(cubic.p1),
            control_end: xy(cubic.p2),
            end: xy(cubic.p3),
        }
    }
}

fn wire(path: &BezPath) -> Vec<Cubic> {
    curves::cubics(path).iter().map(Cubic::from).collect()
}

#[derive(Serialize)]
struct Attachment {
    carrier: Vec<Cubic>,
    paths: Vec<Vec<Cubic>>,
}

#[derive(Serialize)]
struct Frame {
    point: [f64; 2],
    tangent: [f64; 2],
}

#[derive(Serialize)]
#[serde(rename_all = "kebab-case")]
struct CandidateSpec {
    frames: Vec<Frame>,
    positions: Vec<f64>,
    transform: [[f64; 4]; 4],
    origin: [f64; 3],
    corners: [[f64; 3]; 4],
    total: f64,
    preferred: f64,
    accuracy: f64,
    side: f64,
    preferred_side: f64,
    path_index: usize,
    path_shift: f64,
    gap: f64,
    clear_box: bool,
    fixed: bool,
    endpoint_direction: Option<[f64; 2]>,
    attachment: Option<Attachment>,
}

const IDENTITY: [[f64; 4]; 4] = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
];

impl CandidateSpec {
    fn new(frames: Vec<curves::Frame>, positions: Vec<f64>, corners: [[f64; 3]; 4]) -> Self {
        Self {
            frames: frames
                .into_iter()
                .map(|frame| Frame {
                    point: frame.point,
                    tangent: frame.tangent,
                })
                .collect(),
            positions,
            transform: IDENTITY,
            origin: [0.0; 3],
            corners,
            total: 0.0,
            preferred: 0.0,
            accuracy: curves::ACCURACY,
            side: 1.0,
            preferred_side: 1.0,
            path_index: 0,
            path_shift: 0.0,
            gap: 0.0,
            clear_box: true,
            fixed: false,
            endpoint_direction: None,
            attachment: None,
        }
    }

    fn candidates(&self) -> Result<Vec<Candidate>, String> {
        decode(&crate::label_placement::candidates_bytes(&cbor(self)?)?)
    }
}

#[derive(Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub(super) struct Candidate {
    pub bounds: Bounds,
    pub cost: f64,
    pub corners: Vec<[f64; 3]>,
    #[serde(skip)]
    pub arrow_bounds: Vec<Bounds>,
    #[serde(skip)]
    pub arrow_mark: Option<usize>,
}

/// One label's candidates, with its momentum arrow carriers per side.
pub(super) struct Placement {
    pub edge: usize,
    pub page: usize,
    pub candidates: Vec<Candidate>,
    pub href: String,
}

/// CeTZ `content` corners (NW, NE, SW, SE) relative to the placement point.
fn content_corners([width, height, baseline, bounds]: [f64; 4], anchor: Anchor) -> [[f64; 3]; 4] {
    let measure = |pt: f64| (pt / UNIT).abs();
    let w = measure(width) / 2.0;
    let h = (measure(height) + measure(bounds) - measure(baseline)).abs() / 2.0;
    let offset = match anchor {
        Anchor::East => [-w, 0.0],
        Anchor::West => [w, 0.0],
        Anchor::North => [0.0, -h],
        Anchor::South => [0.0, h],
        Anchor::Center => [0.0, 0.0],
    };
    [[-w, h], [w, h], [-w, -h], [w, -h]].map(|[x, y]| [x + offset[0], y + offset[1], 0.0])
}

#[derive(Clone, Copy)]
pub(super) enum Anchor {
    Center,
    East,
    West,
    North,
    South,
}

impl Anchor {
    /// The side of an external label facing away from its free endpoint.
    pub(super) fn outward(dx: f64, dy: f64) -> Self {
        if dx < -1e-9 {
            Self::East
        } else if dx > 1e-9 {
            Self::West
        } else if dy > 0.0 {
            Self::South
        } else if dy < 0.0 {
            Self::North
        } else {
            Self::Center
        }
    }
}

/// The preferred position and the 3/32..29/32 arc positions draw.typ tries,
/// clamped into `[low, high]`, deduplicated and away from the preferred one.
fn positions(total: f64, preferred: f64, low: f64, high: f64) -> Vec<f64> {
    let mut positions = vec![preferred];
    let mut seen = Vec::new();
    for i in 3..30 {
        let at = (total * f64::from(i) / 32.0).clamp(low, high);
        if seen.contains(&at) {
            continue;
        }
        seen.push(at);
        if (at - preferred).abs() > curves::ACCURACY {
            positions.push(at);
        }
    }
    positions
}

/// A particle label sliding along its carrier at `gap`, both sides,
/// interleaved by position.
pub(super) fn carrier_candidates(
    path: &BezPath,
    metrics: [f64; 4],
    gap: f64,
) -> Result<Vec<Candidate>, String> {
    let total = curves::length(path);
    let preferred = (total / 2.0).clamp(0.0, total);
    let positions = positions(total, preferred, 0.0, total);
    let frames = curves::frames(path, positions.clone())?;
    let [left, right] = [1.0, -1.0].map(|side: f64| {
        CandidateSpec {
            total,
            preferred,
            side,
            path_index: usize::from(side < 0.0),
            gap,
            attachment: Some(Attachment {
                carrier: wire(path),
                paths: positions.iter().map(|_| Vec::new()).collect(),
            }),
            ..CandidateSpec::new(
                frames.clone(),
                positions.clone(),
                content_corners(metrics, Anchor::Center),
            )
        }
        .candidates()
    });
    let (left, right) = (left?, right?);
    Ok(left
        .into_iter()
        .zip(right)
        .flat_map(|(a, b)| [a, b])
        .collect())
}

/// A momentum label riding its arrow: per side, an offset carrier, arrow
/// windows at each position, and their footprints. Sides are not
/// interleaved. Returns the candidates and the offset carriers.
pub(super) struct MomentumCandidates {
    specs: Vec<CandidateSpec>,
    marks: Vec<usize>,
}

impl MomentumCandidates {
    pub(super) fn new(
        path: &BezPath,
        metrics: [f64; 4],
        side: f64,
        template: usize,
        batch: &mut marks::Batch,
    ) -> Result<Self, String> {
        let carriers = [side, -side]
            .iter()
            .map(|&s| curves::layer(path, ARROW_OFFSET * s, None, 0.0))
            .collect::<Result<Vec<_>, _>>()?;
        let mut specs = Vec::new();
        let mut marks = Vec::new();
        for (index, carrier) in carriers.iter().enumerate() {
            let total = curves::length(carrier);
            let centered = curves::layer(carrier, 0.0, Some(ARROW_WINDOW), 0.0)?;
            let half = curves::length(&centered) / 2.0;
            let low = half.clamp(0.0, total);
            let high = (total - half).clamp(low, total);
            let preferred = (total / 2.0).clamp(low, high);
            let positions = positions(total, preferred, low, high);
            let frames = curves::frames(carrier, positions.clone())?;
            let shifts: Vec<f64> = positions.iter().map(|at| at - total / 2.0).collect();
            let windows = curves::layers(carrier, &shifts, 0.0, Some(ARROW_WINDOW))?;
            for window in &windows {
                marks.push(batch.push(template, window, kurvst::marks::MarkStation::End, true));
            }
            specs.push(CandidateSpec {
                total,
                preferred,
                side: if index == 0 { side } else { -side },
                preferred_side: side,
                path_index: index,
                gap: MOMENTUM_GAP,
                attachment: Some(Attachment {
                    carrier: wire(path),
                    paths: Vec::new(),
                }),
                ..CandidateSpec::new(frames, positions, content_corners(metrics, Anchor::Center))
            });
        }
        Ok(Self { specs, marks })
    }

    pub(super) fn resolve(
        self,
        geometry: &[kurvst::marks::MarkGeometryOutput],
    ) -> Result<Vec<Candidate>, String> {
        let mut candidates = Vec::new();
        let mut marks = self.marks.into_iter();
        for mut spec in self.specs {
            let indices: Vec<_> = marks.by_ref().take(spec.frames.len()).collect();
            let geometry: Vec<_> = indices.iter().map(|&index| &geometry[index]).collect();
            // Attachment clearance and search obstacles use the same geometry that
            // selected painting will use, including a bent shaft and every head.
            spec.attachment.as_mut().expect("momentum attachment").paths = geometry
                .iter()
                .map(|mark| wire(&mark.footprint.path))
                .collect();
            let mut resolved = spec.candidates()?;
            for ((candidate, mark), index) in resolved.iter_mut().zip(geometry).zip(indices) {
                candidate.arrow_bounds = vec![Bounds::painted(&mark.footprint.path)];
                candidate.arrow_mark = Some(index);
            }
            candidates.extend(resolved);
        }
        Ok(candidates)
    }
}

/// A fixed external label, cleared outward from its free endpoint.
pub(super) fn endpoint_candidates(
    at: [f64; 2],
    direction: [f64; 2],
    anchor: Anchor,
    metrics: [f64; 4],
) -> Result<Vec<Candidate>, String> {
    CandidateSpec {
        gap: EXTERNAL_GAP,
        clear_box: false,
        fixed: true,
        endpoint_direction: Some(direction),
        ..CandidateSpec::new(
            vec![curves::Frame {
                point: at,
                tangent: [1.0, 0.0],
            }],
            vec![0.0],
            content_corners(metrics, anchor),
        )
    }
    .candidates()
}

/// A painted stroke for collisions: its radius and path.
pub(super) struct Stroke {
    pub radius: f64,
    pub path: BezPath,
}

#[derive(Deserialize, Clone, Copy)]
pub(super) struct Line {
    pub start: [f64; 2],
    pub end: [f64; 2],
    pub radius: f64,
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
    corners: Vec<[f64; 2]>,
    arrow_bounds: &'a [Bounds],
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
    edge_lines: &'a serde_json::Value,
    obstacle_padding: f64,
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
}

/// The chosen candidate per placement and the painted collision lines.
pub(super) struct Searched {
    pub choices: Vec<usize>,
    pub lines: Vec<Line>,
}

/// `_relax-label-placements`: one annealing search over all labels, with its
/// exponential comparisons certified by host values as Typst does.
pub(super) fn search(
    placements: &[Placement],
    obstacles: &[Bounds],
    strokes: &[Stroke],
    fixed_arrows: &[Bounds],
) -> Result<Searched, String> {
    let spec = StrokeLinesSpec {
        strokes: strokes.iter().map(PaintedStroke::from).collect(),
        accuracy: 0.005,
    };
    let edge_lines: serde_json::Value =
        decode(&crate::label_placement::stroke_lines_bytes(&cbor(&spec)?)?)?;
    let lines: Vec<Line> =
        serde_json::from_value(edge_lines.clone()).map_err(|error| error.to_string())?;
    if placements
        .iter()
        .all(|placement| placement.candidates.len() == 1)
    {
        return Ok(Searched {
            choices: vec![0; placements.len()],
            lines,
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
                        corners: candidate.corners.iter().map(|c| [c[0], c[1]]).collect(),
                        arrow_bounds: &candidate.arrow_bounds,
                    })
                    .collect(),
            })
            .collect(),
        obstacles,
        label_padding: 0.6,
        fixed_arrows,
        coordinated: true,
        edge_lines: &edge_lines,
        obstacle_padding: 0.08,
    })?;
    let mut math = SearchMath {
        temperatures: (0..84).map(|sweep| 0.15 * 0.9f64.powi(sweep)).collect(),
        exponentials: Vec::new(),
        pair_padding: 0.35,
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
                lines,
            });
        }
    }
}

/// Real (unpadded) overlaps of the chosen movable labels with painted lines
/// and other labels: (overlapping labels, line hits, overlapping label pairs).
pub(super) fn overlaps(
    placements: &[Placement],
    chosen: &[Bounds],
    lines: &[Line],
) -> (usize, usize, usize) {
    let mut overlapping = 0;
    let mut hits = 0;
    let mut pairs = 0;
    for (i, placement) in placements.iter().enumerate() {
        if placement.candidates.len() == 1 {
            continue;
        }
        let bounds = &chosen[i];
        let line_hits = lines
            .iter()
            .filter(|line| {
                (0..=8)
                    .map(|k| {
                        let t = f64::from(k) / 8.0;
                        bounds.distance([
                            line.start[0] + (line.end[0] - line.start[0]) * t,
                            line.start[1] + (line.end[1] - line.start[1]) * t,
                        ])
                    })
                    .fold(f64::INFINITY, f64::min)
                    < line.radius
            })
            .count();
        let label_hits = chosen
            .iter()
            .enumerate()
            .filter(|(j, other)| *j != i && bounds.overlap(other) > 1e-9)
            .count();
        hits += line_hits;
        pairs += label_hits;
        if line_hits + label_hits > 0 {
            overlapping += 1;
        }
    }
    (overlapping, hits, pairs / 2)
}
