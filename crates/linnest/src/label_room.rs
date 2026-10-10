//! Label room: whether a laid-out drawing leaves each measured edge label a
//! clear place, and the layout changes that make room where it does not.
//!
//! The check predicts the drawing's own annotation placement: Kurvst's Hobby
//! curves through the solved routes, label candidates built as `draw.typ`
//! builds them, and the annotation search's costs and final descent. Painted
//! bands grow by a small margin, so a label that would only just clear a line
//! counts as meeting it. A layout whose labels all have a clear place is left
//! unchanged. Otherwise the edges whose paint crosses a label's best place
//! bend around it while every other edge, vertex and external leg holds
//! still, eased back to the least bend that gives the labels room and never
//! bringing painted bands into contact; failing that, the drawing scales up
//! uniformly, which keeps its shape exactly. A remedy is kept only when the
//! drawing's own search leaves its labels less severely crossed.
use std::{
    cell::{OnceCell, RefCell},
    collections::BTreeMap,
    ops::Range,
};

use kurbo::{BezPath, Point};
use linnet::half_edge::{
    involution::EdgeIndex,
    layout::impred::{
        movement::{AxisConstraint, MovementDofs, PointRecord},
        EdgeLabel, ImpredConfig, ImpredLayout,
    },
};
use serde_json::Value;

use crate::{
    annotation::{
        carrier_candidates, endpoint_candidates, Anchor, ArrowStyle, Bounds, Candidate,
        MomentumCandidates, Places, NODE_RADIUS,
    },
    curves,
    label_placement::{predict_placement, LabelOverlaps, LineGrid, StrokeLine, Thoroughness},
    NodeIndex, TypstGraph,
};

/// Largest uniform scale of a drawing, its default band margin, and the
/// bisection steps that find the smallest sufficient scale.
const MAX_SCALE: f64 = 1.5;
const MARGIN: f64 = 0.03;
const BISECTIONS: usize = 4;
/// The bending pass: cooling progress resumed from and its steps.
const BEND_WARM: f64 = 0.8;
const BEND_STEPS: usize = 100;
/// What a bend may change where two edges' painted bands touch: no new
/// place, and at most this growth of the length within reach, as a factor
/// and in canvas units beyond it, so a crossing may become shallower but no
/// band comes to lie along another.
const CONTACT_GROWTH: f64 = 2.0;
const CONTACT_SLACK: f64 = 0.05;

/// One edge's measured annotation style, as `graph.style` records it.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct RoomEdge {
    /// Half extents of its label's text.
    label: Option<[f64; 2]>,
    /// Half width of its painted band; `None` for an unpainted edge.
    band: Option<f64>,
    /// Clearance of its label from the carrier, or from the arrow it rides.
    gap: f64,
    /// Its momentum arrow, carrying the label on an internal edge.
    arrow: Option<RoomArrow>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct RoomArrow {
    style: ArrowStyle,
    /// Arrowhead length and width.
    head: [f64; 2],
}

impl RoomArrow {
    /// The arrow painted over `window`: its shaft and a triangular head at
    /// the window's end.
    fn footprint(self, window: &BezPath) -> BezPath {
        let mut path = window.clone();
        let Some(last) = curves::cubics(window).pop() else {
            return path;
        };
        let [tx, ty] = curves::cubic_tangent(&last, 1.0);
        let length = tx.hypot(ty);
        let [head_length, head_width] = self.head;
        if length > 1e-12 && head_length > 0.0 {
            let (tx, ty) = (tx / length, ty / length);
            let tip = last.p3;
            let base = Point::new(tip.x - tx * head_length, tip.y - ty * head_length);
            let side = |sign: f64| {
                Point::new(
                    base.x - sign * ty * head_width / 2.0,
                    base.y + sign * tx * head_width / 2.0,
                )
            };
            path.move_to(side(1.0));
            path.line_to(tip);
            path.line_to(side(-1.0));
            path.close_path();
        }
        path
    }
}

/// What a predicted placement leaves unresolved: labels whose text meets a
/// painted band, another label, a node or an arrow, how many such contacts
/// there are, and the band length inside text.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct Defects {
    pub labels: usize,
    pub contacts: usize,
    pub crossing: f64,
}

impl Defects {
    fn of(overlaps: &[LabelOverlaps]) -> Self {
        let mut defects = Self::default();
        for (index, found) in overlaps.iter().enumerate() {
            let pairs = found.labels.iter().filter(|&&other| other > index).count();
            let contacts = found.lines.len() + pairs + found.obstacles.len() + found.arrows.len();
            defects.labels += usize::from(contacts > 0 || !found.labels.is_empty());
            defects.contacts += contacts;
            defects.crossing += found.lines.iter().map(|(_, inside)| inside).sum::<f64>();
        }
        defects
    }

    pub(crate) fn is_clear(self) -> bool {
        self.labels == 0
    }

    /// How much the defects show: one per contact, and one per canvas unit
    /// of band inside text, so a band running through a label weighs more
    /// than a few that only touch.
    fn severity(self) -> f64 {
        self.contacts as f64 + self.crossing
    }

    /// Less severe first, then fewer defective labels.
    fn better_than(self, other: Self) -> bool {
        let (severity, others) = (self.severity(), other.severity());
        severity < others - 1e-9 || (severity <= others + 1e-9 && self.labels < other.labels)
    }
}

/// How a bending pass holds the labels it bends around: each in its
/// predicted slot, or free to slide along its carrier.
#[derive(Clone, Copy, PartialEq)]
enum Hold {
    Slots,
    Carriers,
}

/// The best remedy so far: its defects and its layout, if any changed.
struct Remedy {
    defects: Defects,
    layout: Option<ImpredLayout>,
}

impl Remedy {
    /// Keep `layout` if it leaves less severe defects; whether it leaves
    /// none.
    fn keep(&mut self, defects: Defects, layout: ImpredLayout) -> bool {
        if defects.better_than(self.defects) {
            *self = Self {
                defects,
                layout: Some(layout),
            };
        }
        defects.is_clear()
    }
}

/// A layout's predicted paint, which all its predictions share: each edge's
/// visible carrier, the painted bands grown by the margin with each edge's
/// range of them, the fixed arrows and the node boxes, and each internal
/// label's candidates at every place once a prediction has built them.
struct Painted {
    paths: Vec<BezPath>,
    lines: Vec<StrokeLine>,
    own: Vec<Range<usize>>,
    /// The edge of each painted line, and, once needed, where they lie.
    line_edges: Vec<usize>,
    grid: OnceCell<LineGrid>,
    fixed_arrows: Vec<Bounds>,
    nodes: Vec<Bounds>,
    every: RefCell<Vec<Option<Vec<Candidate>>>>,
}

/// Where two edges' painted bands touch: in how many separate places, each
/// a run of consecutive lines along one of them, and the length of either's
/// band within reach of the other's, halved.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct Touch {
    places: usize,
    length: f64,
}

impl Painted {
    fn grid(&self) -> &LineGrid {
        self.grid.get_or_init(|| LineGrid::new(&self.lines))
    }

    /// Painted contact between edges, by pair of edges: bands closer than
    /// `clearance`.
    fn contact(&self, clearance: f64) -> BTreeMap<(usize, usize), Touch> {
        let mut contact: BTreeMap<(usize, usize), Touch> = BTreeMap::new();
        // Per edge and edge it touches: the last line touching it, and the
        // runs of lines that do.
        let mut runs: BTreeMap<(usize, usize), (usize, usize)> = BTreeMap::new();
        let mut met = Vec::new();
        for (index, line) in self.lines.iter().enumerate() {
            let edge = self.line_edges[index];
            met.clear();
            for other in self.grid().near(line.extent().expanded(clearance)) {
                let other_edge = self.line_edges[other];
                if other_edge != edge
                    && !met.contains(&other_edge)
                    && line.touches(&self.lines[other], clearance)
                {
                    met.push(other_edge);
                }
            }
            for &other in &met {
                let (last, count) = runs.entry((edge, other)).or_insert((usize::MAX, 0));
                if last.wrapping_add(1) != index {
                    *count += 1;
                }
                *last = index;
                let touch = contact
                    .entry((edge.min(other), edge.max(other)))
                    .or_default();
                touch.length += line.length() / 2.0;
                touch.places = touch.places.max(*count);
            }
        }
        contact
    }
}

/// A layout under evaluation: its paint, and what has been predicted of it.
struct Probe {
    layout: ImpredLayout,
    painted: Painted,
    /// Each label at its preferred place, those meeting something cleared
    /// among their other places: the first stages of [`LabelRoom::staged`].
    cleared: OnceCell<Prediction>,
    contact: OnceCell<BTreeMap<(usize, usize), Touch>>,
    /// How close bands may come before they count as touching.
    clearance: f64,
}

impl Probe {
    /// The painted contact between edges, measured once.
    fn contact(&self) -> &BTreeMap<(usize, usize), Touch> {
        self.contact
            .get_or_init(|| self.painted.contact(self.clearance))
    }

    /// Whether the probe's bands touch where those of `before` do not: in a
    /// new place, or along much more of their length.
    fn adds_contact(&self, before: &Self) -> bool {
        let before = before.contact();
        self.contact().iter().any(|(pair, touch)| {
            let was = before.get(pair).copied().unwrap_or_default();
            touch.places > was.places || touch.length > CONTACT_GROWTH * was.length + CONTACT_SLACK
        })
    }
}

/// A predicted drawing's labels: each label's own carrier lines and
/// candidates and, for internal labels, its edge.
struct PredictedScene {
    labels: Vec<(Range<usize>, Vec<Candidate>)>,
    owners: Vec<Option<usize>>,
}

/// A layout's predicted placement: defects, and each internal label's chosen
/// text and arrow box.
#[derive(Clone)]
struct Prediction {
    defects: Defects,
    slots: Vec<Option<Bounds>>,
    /// The edges whose internal labels meet something.
    meeting: Vec<usize>,
}

/// The label room of a drawing: its measured annotation styles and node
/// boxes, and the room options.
#[derive(Clone, Debug)]
pub(crate) struct LabelRoom {
    edges: Vec<RoomEdge>,
    nodes: Vec<[f64; 2]>,
    margin: f64,
    max_scale: f64,
}

impl LabelRoom {
    /// The measured label room of a layout graph under the `impred-label-room`
    /// options; `None` when room is off or no edge has a measured label.
    pub(crate) fn new(graph: &TypstGraph, options: &Value) -> Result<Option<Self>, String> {
        let flag = |key: &str| -> Result<bool, String> {
            options
                .get(key)
                .map_or(Some(true), Value::as_bool)
                .ok_or_else(|| format!("{key} must be a boolean"))
        };
        let number = |key: &str, default: f64| -> Result<f64, String> {
            options
                .get(key)
                .map_or(Some(default), Value::as_f64)
                .filter(|value| value.is_finite() && *value >= 0.0)
                .ok_or_else(|| format!("{key} must be a nonnegative number"))
        };
        let margin = number("impred-label-room-margin", MARGIN)?;
        let max_scale = number("impred-label-room-max-scale", MAX_SCALE)?;
        if max_scale < 1.0 {
            return Err("impred-label-room-max-scale must be at least 1".to_owned());
        }
        if !flag("impred-label-room")? {
            return Ok(None);
        }
        let edges: Vec<RoomEdge> = (0..graph.n_edges())
            .map(|edge| RoomEdge::read(&graph[EdgeIndex(edge)].statements))
            .collect::<Result<_, _>>()?;
        if edges.iter().all(|edge| edge.label.is_none()) {
            return Ok(None);
        }
        let nodes = (0..graph.node_bound())
            .map(|node| {
                let statements = &graph[NodeIndex(node)].statements;
                ["layout-width", "layout-height"]
                    .map(|key| positive(statements, key).map_or(0.0, |size| size / 2.0))
            })
            .collect();
        Ok(Some(Self {
            edges,
            nodes,
            margin,
            max_scale,
        }))
    }

    /// The layout that makes room, if `layout` needs and can take one: none
    /// when every label has a clear place. Otherwise the first remedy the
    /// drawing's own search clears: bend the edges crossing labels, which
    /// keeps every vertex; else scale the drawing up by the smallest
    /// sufficient factor, which keeps its shape exactly; else scale it by
    /// the smallest factor after which bending suffices. Without a clear
    /// remedy, the one leaving the least severe defects, if it improves on
    /// the layout. No bend brings painted bands into contact. Authored fixed
    /// coordinates never move: such drawings only bend.
    pub(crate) fn remedy(
        &self,
        layout: &ImpredLayout,
        config: ImpredConfig,
    ) -> Result<Option<ImpredLayout>, String> {
        let original = self.probe(layout.clone())?;
        let predicted = self.staged(&original, Thoroughness::Hopeful)?;
        if predicted.defects.is_clear() {
            return Ok(None);
        }
        // A remedy must improve on the drawing's own placement.
        let mut best = Remedy {
            defects: self.predict(&original, Thoroughness::Exact)?.defects,
            layout: None,
        };
        if best.defects.is_clear() {
            return Ok(None);
        }
        if let Some((bent, defects)) = self.bend(&original, &predicted, Hold::Slots, config)? {
            if best.keep(defects, bent.layout) {
                return Ok(best.layout);
            }
        }
        let Some(scaling) = Scaling::new(layout) else {
            return Ok(best.layout);
        };
        let Some(widest) = scaling.at(self.max_scale) else {
            return Ok(best.layout);
        };
        let widest = self.probe(widest)?;
        // Scaling alone: the smallest factor whose settled placement is as
        // good as the widest's, growing, when the widest's is clear, until
        // the drawing's search agrees.
        let widest_prediction = self.staged(&widest, Thoroughness::Hopeful)?;
        let target = widest_prediction.defects;
        if target.better_than(best.defects) || target.is_clear() {
            let (factor, mut passed) = self.smallest(&scaling, |scaled| {
                Ok(!target.better_than(self.staged(scaled, Thoroughness::Settled)?.defects))
            })?;
            let growing = if target.is_clear() { 3 } else { 1 };
            for factor in [factor, (factor + self.max_scale) / 2.0, self.max_scale]
                .into_iter()
                .take(growing)
            {
                // The bisection's last success already holds the first factor.
                let scaled = match passed.take() {
                    Some(scaled) => scaled,
                    None => match scaling.at(factor) {
                        Some(scaled) => self.probe(scaled)?,
                        None => continue,
                    },
                };
                if best.keep(self.confirmed(&scaled)?, scaled.layout) {
                    return Ok(best.layout);
                }
                if factor >= self.max_scale {
                    break;
                }
            }
        }
        // Scaling, then bending: the smallest factor whose bent drawing the
        // search clears, if the widest's is, holding labels in their slots
        // or, failing that, letting them slide along their carriers.
        let mut cleared = None;
        for hold in [Hold::Slots, Hold::Carriers] {
            let Some((bent, defects)) = self.bend(&widest, &widest_prediction, hold, config)?
            else {
                continue;
            };
            if best.keep(defects, bent.layout) {
                cleared = Some(hold);
                break;
            }
        }
        let Some(hold) = cleared else {
            return Ok(best.layout);
        };
        self.smallest(&scaling, |scaled| {
            let slots = self.predict(scaled, Thoroughness::Settled)?;
            let Some((bent, defects)) = self.bend(scaled, &slots, hold, config)? else {
                return Ok(false);
            };
            if defects.is_clear() {
                best = Remedy {
                    defects,
                    layout: Some(bent.layout),
                };
            }
            Ok(defects.is_clear())
        })?;
        Ok(best.layout)
    }

    /// The smallest factor, by bisection between 1 and the maximum scale,
    /// whose scaled layout passes `test`, with that layout's probe; the
    /// maximum, untested, when none does.
    fn smallest(
        &self,
        scaling: &Scaling,
        mut test: impl FnMut(&Probe) -> Result<bool, String>,
    ) -> Result<(f64, Option<Probe>), String> {
        let (mut low, mut high) = (1.0, self.max_scale);
        let mut passed = None;
        for _ in 0..BISECTIONS {
            let middle = (low + high) / 2.0;
            let scaled = scaling
                .at(middle)
                .map(|scaled| self.probe(scaled))
                .transpose()?;
            match scaled {
                Some(scaled) if test(&scaled)? => {
                    high = middle;
                    passed = Some(scaled);
                }
                _ => low = middle,
            }
        }
        Ok((high, passed))
    }

    /// Defects of a remedy: none when a clear placement is sure; otherwise
    /// quickly settled, then, if that is clear, as the drawing's own
    /// annealing search leaves them.
    fn confirmed(&self, probe: &Probe) -> Result<Defects, String> {
        Ok(self.staged(probe, Thoroughness::Confirmed)?.defects)
    }

    /// A clear placement the drawing's search is sure to reach, else the
    /// placement `thoroughness` predicts. Labels clear at their preferred
    /// places, or clear once those that are not move, need nothing more:
    /// the search finishes no worse than every label at its preferred place,
    /// cleared of what it can avoid. Only labels meeting something at their
    /// preferred places take their other candidates, so most layouts that
    /// leave labels their places are quick to tell; that cleared placement
    /// also stands for a settled one.
    fn staged(&self, probe: &Probe, thoroughness: Thoroughness) -> Result<Prediction, String> {
        let cleared = match probe.cleared.get() {
            Some(cleared) => cleared.clone(),
            None => {
                let preferred =
                    self.predict_in(probe, |_| Places::Preferred, Thoroughness::Cleared)?;
                let cleared = if preferred.defects.is_clear() {
                    preferred
                } else {
                    let moving = |edge| {
                        if preferred.meeting.contains(&edge) {
                            Places::Every
                        } else {
                            Places::Preferred
                        }
                    };
                    self.predict_in(probe, moving, Thoroughness::Cleared)?
                };
                probe.cleared.get_or_init(|| cleared).clone()
            }
        };
        if cleared.defects.is_clear() || thoroughness == Thoroughness::Settled {
            return Ok(cleared);
        }
        self.predict_in(probe, |_| Places::Every, thoroughness)
    }

    /// The layout bent around each meeting internal label's predicted
    /// place, with its confirmed defects, if a bend can help and brings no
    /// painted bands into contact: a short warm pass that moves only the
    /// edges whose bands cross those places, holding every other route
    /// point, vertex and external leg, with each label's text and arrow box,
    /// grown by the margin, held in its slot or beside its carrier; then
    /// eased back as far as the labels allow.
    fn bend(
        &self,
        probe: &Probe,
        prediction: &Prediction,
        hold: Hold,
        config: ImpredConfig,
    ) -> Result<Option<(Probe, Defects)>, String> {
        let layout = &probe.layout;
        let painted = &probe.painted;
        // The edges crossing a meeting label's place, other than its own.
        let mut free = vec![false; self.edges.len()];
        for &edge in &prediction.meeting {
            let Some(slot) = prediction.slots[edge] else {
                continue;
            };
            for line in painted.grid().near(slot.expanded(1.0)) {
                let crossing = painted.line_edges[line];
                if crossing != edge
                    && layout.external[crossing].is_none()
                    && painted.lines[line].reaches(slot, self.margin)
                {
                    free[crossing] = true;
                }
            }
        }
        if !free.contains(&true) {
            return Ok(None);
        }
        let mut bent = layout.clone();
        let held = layout.node_points.iter().copied().chain(
            layout
                .routes
                .iter()
                .zip(&free)
                .filter(|(_, free)| !**free)
                .flat_map(|(route, _)| route.iter().copied()),
        );
        for point in held {
            let record = &mut bent.constraints[point];
            record.fixed_axes = [true, true];
            record.reference = bent.positions[point];
        }
        bent.labels = prediction
            .slots
            .iter()
            .map(|slot| {
                slot.map(|slot| EdgeLabel {
                    extents: [
                        (slot.right - slot.left) / 2.0 + self.margin,
                        (slot.top - slot.bottom) / 2.0 + self.margin,
                    ],
                    center: Some([
                        (slot.left + slot.right) / 2.0,
                        (slot.bottom + slot.top) / 2.0,
                    ]),
                    pinned: hold == Hold::Slots,
                    pin: None,
                })
            })
            .collect();
        bent.solve(ImpredConfig {
            labels: true,
            warm_start: BEND_WARM,
            steps: BEND_STEPS,
            external_max_points: 0,
            ..config
        })?;
        // The held coordinates keep their authored constraints from here on.
        bent.constraints.clone_from(&layout.constraints);
        // Label-aware layouts keep their labels, now beside the bent edges.
        if config.labels {
            for label in bent.labels.iter_mut().flatten() {
                label.pinned = false;
                label.pin = None;
            }
        } else {
            bent.labels.clear();
        }
        let bent = self.probe(bent)?;
        if bent.adds_contact(probe) {
            // The most of the bend that brings no bands into contact.
            let Some(partial) = self.partial(probe, &bent, true, |partial| {
                Ok(!partial.adds_contact(probe))
            })?
            else {
                return Ok(None);
            };
            let defects = self.confirmed(&partial)?;
            return Ok(Some((partial, defects)));
        }
        // The least of it that settles the labels as well as the whole.
        let full = self.confirmed(&bent)?;
        let target = self.staged(&bent, Thoroughness::Settled)?.defects;
        let eased = self.partial(probe, &bent, false, |partial| {
            Ok(!partial.adds_contact(probe)
                && !target.better_than(self.staged(partial, Thoroughness::Settled)?.defects))
        })?;
        if let Some(eased) = eased {
            let defects = self.confirmed(&eased)?;
            if !full.better_than(defects) {
                return Ok(Some((eased, defects)));
            }
        }
        Ok(Some((bent, full)))
    }

    /// Part of a bend: by bisection of the way from the probed layout to
    /// the bent one, the most of it that passes `test` if `most`, otherwise
    /// the least; `None` when no part short of the whole and beyond none
    /// passes.
    fn partial(
        &self,
        probe: &Probe,
        bent: &Probe,
        most: bool,
        mut test: impl FnMut(&Probe) -> Result<bool, String>,
    ) -> Result<Option<Probe>, String> {
        let (mut low, mut high) = (0.0, 1.0);
        let mut found = None;
        for _ in 0..BISECTIONS {
            let middle = (low + high) / 2.0;
            let along = |from: [f64; 2], to: [f64; 2]| {
                [0, 1].map(|axis| from[axis] + middle * (to[axis] - from[axis]))
            };
            let mut layout = bent.layout.clone();
            for (point, &from) in layout.positions.iter_mut().zip(&probe.layout.positions) {
                *point = along(from, *point);
            }
            for (label, from) in layout.labels.iter_mut().zip(&probe.layout.labels) {
                if let (Some(label), Some(from)) = (label, from) {
                    if let (Some(center), Some(from)) = (label.center.as_mut(), from.center) {
                        *center = along(from, *center);
                    }
                }
            }
            let partial = if constraints_hold(&layout) {
                Some(self.probe(layout)?)
            } else {
                None
            };
            let passed = match &partial {
                Some(partial) => test(partial)?,
                None => false,
            };
            // The passing side moves toward the bend's other end.
            if passed == most {
                low = middle;
            } else {
                high = middle;
            }
            if passed {
                found = partial;
            }
        }
        Ok(found)
    }

    /// The drawing's predicted placement of the layout's labels.
    fn predict(&self, probe: &Probe, thoroughness: Thoroughness) -> Result<Prediction, String> {
        self.predict_in(probe, |_| Places::Every, thoroughness)
    }

    /// The predicted placement of the layout's labels among the `places` of
    /// each internal label's edge.
    fn predict_in(
        &self,
        probe: &Probe,
        places: impl Fn(usize) -> Places,
        thoroughness: Thoroughness,
    ) -> Result<Prediction, String> {
        let painted = &probe.painted;
        let scene = self.scene(&probe.layout, painted, places)?;
        let boxes: Vec<Vec<Bounds>> = scene
            .labels
            .iter()
            .map(|(_, candidates)| {
                candidates
                    .iter()
                    .map(|candidate| {
                        candidate
                            .arrow_bounds
                            .iter()
                            .fold(candidate.bounds, |all, &arrow| all.union(arrow))
                    })
                    .collect()
            })
            .collect();
        let (choices, overlaps) = predict_placement(
            scene.labels,
            painted.lines.clone(),
            self.margin,
            &painted.nodes,
            &painted.fixed_arrows,
            thoroughness,
        )?;
        // Each internal label's chosen text with its arrow, by edge.
        let mut slots = vec![None; self.edges.len()];
        for ((owner, boxes), choice) in scene.owners.iter().zip(&boxes).zip(choices) {
            if let Some(edge) = *owner {
                slots[edge] = Some(boxes[choice]);
            }
        }
        let meeting = scene
            .owners
            .iter()
            .zip(&overlaps)
            .filter_map(|(owner, found)| owner.filter(|_| found.contacts() > 0))
            .collect();
        Ok(Prediction {
            defects: Defects::of(&overlaps),
            slots,
            meeting,
        })
    }

    /// Whether `layout`'s painted bands touch where those of `before` do
    /// not, as no bend may make them.
    #[cfg(all(test, feature = "svg"))]
    pub(crate) fn adds_contact(&self, layout: &ImpredLayout, before: &ImpredLayout) -> bool {
        let probe = |layout: &ImpredLayout| self.probe(layout.clone()).unwrap();
        probe(layout).adds_contact(&probe(before))
    }

    /// A layout with its predicted paint, for predictions to share.
    fn probe(&self, layout: ImpredLayout) -> Result<Probe, String> {
        let painted = self.painted(&layout)?;
        Ok(Probe {
            layout,
            painted,
            cleared: OnceCell::new(),
            contact: OnceCell::new(),
            clearance: self.margin,
        })
    }

    /// The paint a layout predicts: each edge's visible carrier, its painted
    /// band, an arrow beside an edge without a label to ride it, and the node
    /// boxes.
    fn painted(&self, layout: &ImpredLayout) -> Result<Painted, String> {
        if layout.routes.len() != self.edges.len() || layout.node_points.len() != self.nodes.len() {
            return Err("label room does not match the layout graph".to_owned());
        }
        let nodes: BTreeMap<usize, usize> = layout
            .node_points
            .iter()
            .enumerate()
            .map(|(node, &point)| (point, node))
            .collect();
        // Edges stop at the vertex they meet, at least the smallest radius.
        let outset = |point: &usize| {
            nodes.get(point).map_or(0.0, |&node| {
                NODE_RADIUS
                    .max(self.nodes[node][0])
                    .max(self.nodes[node][1])
            })
        };
        let paths: Vec<BezPath> = layout
            .routes
            .iter()
            .map(|route| {
                let points: Vec<[f64; 2]> =
                    route.iter().map(|&point| layout.positions[point]).collect();
                curves::trim_routed(
                    &curves::routed_curve(&points)?,
                    outset(&route[0]),
                    outset(route.last().expect("routes have two points")),
                )
            })
            .collect::<Result<_, _>>()?;
        let painted: Vec<Vec<StrokeLine>> = paths
            .iter()
            .zip(&self.edges)
            .map(|(path, style)| {
                style.band.map_or_else(Vec::new, |band| {
                    StrokeLine::flattened(
                        curves::cubics(path).into_iter().map(|cubic| {
                            [cubic.p0, cubic.p1, cubic.p2, cubic.p3].map(|p| [p.x, p.y])
                        }),
                        band,
                    )
                })
            })
            .collect();
        // Each edge's lines, which its own label clears by construction.
        let own: Vec<Range<usize>> = painted
            .iter()
            .scan(0, |start, lines| {
                let range = *start..*start + lines.len();
                *start = range.end;
                Some(range)
            })
            .collect();
        let line_edges = own
            .iter()
            .enumerate()
            .flat_map(|(edge, lines)| lines.clone().map(move |_| edge))
            .collect();
        let lines: Vec<StrokeLine> = painted.into_iter().flatten().collect();
        // External legs and edges without a label carry fixed arrows.
        let mut fixed_arrows = Vec::new();
        for ((path, style), external) in paths.iter().zip(&self.edges).zip(&layout.external) {
            if let Some(arrow) = style.arrow.filter(|_| {
                !path.elements().is_empty() && (external.is_some() || style.label.is_none())
            }) {
                let window =
                    curves::layer(path, arrow.style.offset, Some(arrow.style.window), 0.0)?;
                fixed_arrows.push(Bounds::painted(&arrow.footprint(&window)));
            }
        }
        let nodes = layout
            .node_points
            .iter()
            .zip(&self.nodes)
            .map(|(&point, &[w, h])| {
                let [x, y] = layout.positions[point];
                Bounds {
                    left: x - w,
                    right: x + w,
                    bottom: y - h,
                    top: y + h,
                }
            })
            .collect();
        Ok(Painted {
            every: RefCell::new(vec![None; paths.len()]),
            paths,
            grid: OnceCell::new(),
            lines,
            own,
            line_edges,
            fixed_arrows,
            nodes,
        })
    }

    /// The labels a layout's paint predicts: each label's candidates among
    /// the `places` of its edge.
    fn scene(
        &self,
        layout: &ImpredLayout,
        painted: &Painted,
        places: impl Fn(usize) -> Places,
    ) -> Result<PredictedScene, String> {
        let mut labels = Vec::new();
        // The edge of each internal label, by label.
        let mut owners = Vec::new();
        for (edge, (path, style)) in painted.paths.iter().zip(&self.edges).enumerate() {
            if path.elements().is_empty() {
                continue;
            }
            let route = &layout.routes[edge];
            let own = painted.own[edge].clone();
            if layout.external[edge].is_some() {
                // External labels stay beyond their free end.
                let Some(half) = style.label else { continue };
                let tip = layout.positions[layout.anchor_points[edge]];
                let starts = route[0] == layout.anchor_points[edge];
                let node = layout.positions[if starts {
                    *route.last().expect("routes have two points")
                } else {
                    route[0]
                }];
                let cubics = curves::cubics(path);
                let mut direction = match (starts, cubics.first(), cubics.last()) {
                    (true, Some(first), _) => curves::cubic_tangent(first, 0.0).map(|v| -v),
                    (false, _, Some(last)) => curves::cubic_tangent(last, 1.0),
                    _ => [0.0; 2],
                };
                let (dx, dy) = (tip[0] - node[0], tip[1] - node[1]);
                if direction[0].hypot(direction[1]) <= 1e-9 {
                    direction = [dx, dy];
                }
                if direction[0].hypot(direction[1]) <= 1e-9 {
                    direction = [1.0, 0.0];
                }
                let corners = Anchor::outward(dx, dy).corners(half);
                labels.push((own, endpoint_candidates(tip, direction, corners)?));
                owners.push(None);
                continue;
            }
            let Some(half) = style.label else { continue };
            let places = places(edge);
            owners.push(Some(edge));
            if let Some(every) = painted.every.borrow()[edge]
                .as_ref()
                .filter(|_| places == Places::Every)
            {
                labels.push((own, every.clone()));
                continue;
            }
            let corners = Anchor::Center.corners(half);
            let anchor = layout.positions[layout.anchor_points[edge]];
            let candidates = if let Some(arrow) = style.arrow {
                // The bend relative to the chord sets the preferred side.
                let [a, b] = [route[0], *route.last().expect("routes have two points")]
                    .map(|point| layout.positions[point]);
                let cross = (b[0] - a[0]) * (anchor[1] - a[1]) - (b[1] - a[1]) * (anchor[0] - a[0]);
                let side = if cross.abs() > 1e-9 && cross < 0.0 {
                    -1.0
                } else {
                    1.0
                };
                let proposals = MomentumCandidates::new(path, corners, side, arrow.style, places)?;
                // The text clears the arrow's shaft, as `draw.typ` attaches it;
                // the search weighs the arrow's box with its head.
                let arrows: Vec<_> = proposals
                    .windows()
                    .iter()
                    .map(|window| (window.clone(), Bounds::painted(&arrow.footprint(window))))
                    .collect();
                proposals.resolve(&arrows)?
            } else {
                // A label-aware layout's label point sets the preferred side.
                let toward = layout
                    .labels
                    .get(edge)
                    .copied()
                    .flatten()
                    .and_then(|label| label.center)
                    .map(|[x, y]| [x - anchor[0], y - anchor[1]]);
                let side = match (
                    toward,
                    curves::frames(path, vec![curves::length(path) / 2.0])?.first(),
                ) {
                    (Some([dx, dy]), Some(frame)) if dx.hypot(dy) > 1e-9 => {
                        let [tx, ty] = frame.tangent;
                        if tx * dy - ty * dx < 0.0 {
                            -1.0
                        } else {
                            1.0
                        }
                    }
                    _ => 1.0,
                };
                carrier_candidates(path, corners, style.gap, side, places)?
            };
            if places == Places::Every {
                painted.every.borrow_mut()[edge] = Some(candidates.clone());
            }
            labels.push((own, candidates));
        }
        Ok(PredictedScene { labels, owners })
    }
}

impl RoomEdge {
    /// The `label-*` statements `graph.style` measures.
    fn read(statements: &BTreeMap<String, String>) -> Result<Self, String> {
        let label = match (
            positive(statements, "label-width"),
            positive(statements, "label-height"),
        ) {
            (Some(width), Some(height)) => Some([width / 2.0, height / 2.0]),
            _ => None,
        };
        let band = statements
            .get("label-room-band")
            .map(|value| number(value, "label-room-band"))
            .transpose()?;
        let gap = statements
            .get("label-room-gap")
            .map(|value| number(value, "label-room-gap"))
            .transpose()?
            .unwrap_or(0.0);
        let arrow = statements
            .get("label-room-arrow")
            .map(|value| {
                let values: Vec<&str> = value.split_whitespace().collect();
                let [offset, length, ratio, head_length, head_width] = values[..] else {
                    return Err(format!(
                        "label-room-arrow must be `offset length ratio head-length head-width`, not {value:?}"
                    ));
                };
                Ok(RoomArrow {
                    style: ArrowStyle {
                        offset: number(offset, "label-room-arrow offset")?,
                        window: (
                            number(length, "label-room-arrow length")?,
                            if ratio == "none" {
                                f64::INFINITY
                            } else {
                                number(ratio, "label-room-arrow ratio")?
                            },
                        ),
                        gap,
                    },
                    head: [
                        number(head_length, "label-room-arrow head length")?,
                        number(head_width, "label-room-arrow head width")?,
                    ],
                })
            })
            .transpose()?;
        Ok(Self {
            label,
            band,
            gap,
            arrow,
        })
    }
}

fn number(value: &str, name: &str) -> Result<f64, String> {
    value
        .trim()
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite() && *value >= 0.0)
        .ok_or_else(|| format!("{name} must be a nonnegative number, not {value:?}"))
}

fn positive(statements: &BTreeMap<String, String>, key: &str) -> Option<f64> {
    statements
        .get(key)
        .and_then(|value| value.trim().parse::<f64>().ok())
        .filter(|value| value.is_finite() && *value > 0.0)
}

/// Uniform scalings of a layout about a center that keeps every authored
/// fixed coordinate: on an axis with fixed coordinates, their common value;
/// otherwise the vertex centroid. Labels of a label-aware layout scale with it.
struct Scaling<'a> {
    layout: &'a ImpredLayout,
    center: [f64; 2],
}

impl<'a> Scaling<'a> {
    /// `None` when fixed coordinates of one axis differ, or the layout holds
    /// its vertices in place: no scaling keeps them.
    fn new(layout: &'a ImpredLayout) -> Option<Self> {
        if layout.nodes_fixed {
            return None;
        }
        let count = layout.node_points.len().max(1) as f64;
        let mut center = layout.node_points.iter().fold([0.0; 2], |sum, &point| {
            let [x, y] = layout.positions[point];
            [sum[0] + x / count, sum[1] + y / count]
        });
        for (axis, center) in center.iter_mut().enumerate() {
            let mut fixed = layout
                .constraints
                .iter()
                .zip(&layout.positions)
                .filter(|(record, _)| fixed_axis(record, axis))
                .map(|(_, position)| position[axis]);
            if let Some(value) = fixed.next() {
                if fixed.any(|other| (other - value).abs() > 1e-9) {
                    return None;
                }
                *center = value;
            }
        }
        Some(Self { layout, center })
    }

    /// The layout scaled by `factor`, if its constraints still hold.
    fn at(&self, factor: f64) -> Option<ImpredLayout> {
        let [cx, cy] = self.center;
        let map = |[x, y]: [f64; 2]| [cx + factor * (x - cx), cy + factor * (y - cy)];
        let mut scaled = self.layout.clone();
        for (point, record) in scaled.positions.iter_mut().zip(&self.layout.constraints) {
            *point = map(*point);
            // Fixed coordinates lie on the center: keep them bit for bit.
            for (axis, (coordinate, &center)) in point.iter_mut().zip(&self.center).enumerate() {
                if fixed_axis(record, axis) {
                    *coordinate = center;
                }
            }
        }
        for label in scaled.labels.iter_mut().flatten() {
            label.center = label.center.map(map);
        }
        constraints_hold(&scaled).then_some(scaled)
    }
}

/// Whether an authored coordinate is fixed on `axis`.
fn fixed_axis(record: &PointRecord, axis: usize) -> bool {
    record.fixed_axes[axis] || matches!(record.constraints[axis], AxisConstraint::Fixed)
}

/// Whether the native affine constraints and directed signs hold.
fn constraints_hold(layout: &ImpredLayout) -> bool {
    MovementDofs::new(&layout.constraints, layout.nodes_fixed).is_ok_and(|dofs| {
        dofs.validate(&layout.positions).is_ok()
            && dofs
                .sign_rows(&layout.positions)
                .iter()
                .all(|row| row.bound >= 0.0)
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    /// A theta graph between vertices at `(0, 0)` and `(3, 0)`: a straight
    /// edge 0 through `(1.5, 0)`, labelled `label` wide, between edges 1
    /// and 2 through `(1.5, ±arc)`. Vertices are held unless `free`.
    fn theta(arc: f64, label: f64, free: bool) -> (LabelRoom, ImpredLayout) {
        let point = |node: bool| {
            json!({"constraints": ["Free", "Free"], "shift": [0.0, 0.0],
                "reference": [0.0, 0.0], "fixed_axes": [false, false], "is_node": node})
        };
        let layout: ImpredLayout = serde_json::from_value(json!({
            "positions": [[0.0, 0.0], [3.0, 0.0], [1.5, 0.0], [1.5, arc], [1.5, -arc]],
            "routes": [[0, 2, 1], [0, 3, 1], [0, 4, 1]],
            "node_points": [0, 1],
            "anchor_points": [2, 3, 4],
            "external": [null, null, null],
            "constraints": [point(true), point(true), point(false), point(false), point(false)],
            "pull_scales": [1.0, 1.0, 1.0],
            "samples": [false, false, false, false, false],
            "nodes_fixed": !free,
            "labels": [],
        }))
        .unwrap();
        let edge = RoomEdge {
            label: None,
            band: Some(0.05),
            gap: 0.15,
            arrow: None,
        };
        let room = LabelRoom {
            edges: vec![
                RoomEdge {
                    label: Some([label / 2.0, 0.2]),
                    ..edge
                },
                edge,
                edge,
            ],
            nodes: vec![[NODE_RADIUS; 2]; 2],
            margin: MARGIN,
            max_scale: MAX_SCALE,
        };
        (room, layout)
    }

    #[test]
    fn a_band_through_text_outweighs_touches() {
        let crossed = Defects {
            labels: 1,
            contacts: 1,
            crossing: 5.9,
        };
        let touched = Defects {
            labels: 1,
            contacts: 2,
            crossing: 0.6,
        };
        assert!(touched.better_than(crossed));
        assert!(!crossed.better_than(touched));
        // Equally severe, fewer defective labels win; nothing beats itself.
        let spread = Defects {
            labels: 2,
            ..touched
        };
        assert!(touched.better_than(spread) && !spread.better_than(touched));
        assert!(!touched.better_than(touched));
    }

    #[test]
    fn remedies_never_replace_a_less_severe_layout() {
        let (_, layout) = theta(0.45, 1.2, false);
        let defects = |contacts| Defects {
            labels: 1,
            contacts,
            crossing: 0.0,
        };
        let mut best = Remedy {
            defects: defects(2),
            layout: None,
        };
        assert!(!best.keep(defects(3), layout.clone()));
        assert!(best.layout.is_none());
        assert!(!best.keep(defects(1), layout.clone()));
        assert_eq!(best.defects, defects(1));
        assert!(best.keep(Defects::default(), layout));
        assert!(best.defects.is_clear());
    }

    #[test]
    fn touches_count_their_places_and_bends_may_not_add_any() {
        let (room, layout) = theta(0.45, 1.2, false);
        // A polyline band, as lines at most 0.1 long.
        let band = |points: &[[f64; 2]], radius: f64| {
            let pieces = points.windows(2).flat_map(|w| {
                let steps = ((w[1][0] - w[0][0]).hypot(w[1][1] - w[0][1]) / 0.1).ceil() as usize;
                (0..steps).map(move |step| {
                    let at = |s: usize| {
                        let t = s as f64 / steps as f64;
                        [0, 1].map(|axis| w[0][axis] + t * (w[1][axis] - w[0][axis]))
                    };
                    let (a, b) = (at(step), at(step + 1));
                    [a, a, b, b]
                })
            });
            StrokeLine::flattened(pieces.collect::<Vec<_>>(), radius)
        };
        let painted = |edges: Vec<Vec<StrokeLine>>| {
            let mut probe = room.probe(layout.clone()).unwrap();
            let line_edges = edges
                .iter()
                .enumerate()
                .flat_map(|(edge, lines)| lines.iter().map(move |_| edge))
                .collect();
            probe.painted.lines = edges.into_iter().flatten().collect();
            probe.painted.line_edges = line_edges;
            probe.painted.grid = OnceCell::new();
            probe
        };
        // Edge 1 crosses edge 0 once; edge 2 runs along it at both ends.
        let along = |gap: f64| {
            band(
                &[[0.0, gap], [1.0, gap], [2.0, 1.0], [3.0, gap], [4.0, gap]],
                0.05,
            )
        };
        let before = painted(vec![
            band(&[[0.0, 0.0], [4.0, 0.0]], 0.05),
            band(&[[2.0, -1.0], [2.0, 0.5]], 0.05),
            along(0.08),
        ]);
        let contact = before.contact();
        assert_eq!(contact[&(0, 1)].places, 1);
        assert_eq!(contact[&(0, 2)].places, 2);
        assert!(!contact.contains_key(&(1, 2)));
        assert!(!before.adds_contact(&before));
        // Lifting edge 2 off edge 0 removes contact; laying edge 1 along it
        // adds a place.
        let lifted = painted(vec![
            band(&[[0.0, 0.0], [4.0, 0.0]], 0.05),
            band(&[[2.0, -1.0], [2.0, 0.5]], 0.05),
            along(0.3),
        ]);
        assert!(!lifted.adds_contact(&before));
        assert!(before.adds_contact(&lifted));
        let laid = painted(vec![
            band(&[[0.0, 0.0], [4.0, 0.0]], 0.05),
            band(&[[2.0, -1.0], [2.0, 0.05], [3.5, 0.05]], 0.05),
            along(0.08),
        ]);
        assert!(laid.adds_contact(&before));
    }

    #[test]
    fn room_bends_only_the_edges_crossing_a_label() {
        // Wide arcs leave the label a place: nothing changes.
        let (room, layout) = theta(1.2, 1.2, false);
        let config = ImpredConfig::default();
        assert!(room.remedy(&layout, config).unwrap().is_none());
        // Close arcs cross the label on either side. With the vertices held,
        // only the arc crossing its best place bends, away from it, and just
        // far enough to give the label its place.
        let (room, layout) = theta(0.45, 1.2, false);
        let before = room.probe(layout.clone()).unwrap();
        let crowded = room.predict(&before, Thoroughness::Exact).unwrap();
        assert!(!crowded.defects.is_clear());
        let bent = room.remedy(&layout, config).unwrap().expect("a bend");
        for point in [0, 1, 2, 4] {
            assert_eq!(bent.positions[point], layout.positions[point]);
        }
        let rise = bent.positions[3][1] - layout.positions[3][1];
        assert!(rise > 0.0 && rise < 1.0, "{:?}", bent.positions);
        let after = room.probe(bent).unwrap();
        assert!(room
            .predict(&after, Thoroughness::Exact)
            .unwrap()
            .defects
            .is_clear());
        assert!(!after.adds_contact(&before));
    }
}
