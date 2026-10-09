//! Native paint adapter for Kurvst's data-only mark batches.
use kurbo::{BezPath, Cap, Shape};
use kurvst::{
    CurvePathOutput,
    marks::{
        MarkContext, MarkDirection, MarkGeometryBatchOutput, MarkGeometryMode, MarkGeometrySpec,
        MarkPlacement, MarkSpec, MarkStation, MarkStrokeStyle, MarkTemplateSpec,
    },
};
use serde::Deserialize;

use super::{Stroke, UNIT, labels::Bounds};

/// Paint overrides stay outside the shared geometry flags, in flattened leaf order.
#[derive(Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MarkPaint {
    pub fill: Option<String>,
    pub stroke: Option<String>,
}

impl super::EdgeDrawing {
    /// Preserve the notebook flow size rather than the catalogue's larger default.
    pub fn default_flow_arrow() -> MarkSpec {
        serde_json::from_value(serde_json::json!({
            "shape": "triangle",
            "length": {"points": 0.21 * 72.0 / 2.54, "ratio": 0.0},
            "width": {"points": 0.1575 * 72.0 / 2.54, "ratio": 0.0}
        }))
        .expect("explicit flow dimensions are valid")
    }

    /// Preserve the notebook momentum size.
    pub fn default_momentum_arrow() -> MarkSpec {
        serde_json::from_value(serde_json::json!({
            "shape": "straight",
            "length": {"points": 0.16 * 72.0 / 2.54, "ratio": 0.0},
            "width": {"points": 0.12 * 72.0 / 2.54, "ratio": 0.0}
        }))
        .expect("explicit momentum dimensions are valid")
    }
}

pub(super) struct Batch {
    spec: MarkGeometrySpec,
    templates: std::collections::BTreeMap<String, usize>,
}

#[cfg(test)]
std::thread_local! {
    static GEOMETRY_CALLS: std::cell::RefCell<Vec<MarkGeometryMode>> = const { std::cell::RefCell::new(Vec::new()) };
}

impl Batch {
    pub(super) fn new() -> Self {
        Self {
            spec: MarkGeometrySpec {
                templates: Vec::new(),
                carriers: Vec::new(),
                placements: Vec::new(),
                mode: MarkGeometryMode::Candidates,
            },
            templates: std::collections::BTreeMap::new(),
        }
    }

    pub(super) fn template(&mut self, mark: &MarkSpec, stroke: &Stroke) -> Result<usize, String> {
        let context = Self::context(stroke);
        let key = serde_json::to_string(&(mark, context)).map_err(|e| e.to_string())?;
        if let Some(&index) = self.templates.get(&key) {
            return Ok(index);
        }
        let index = self.spec.templates.len();
        self.spec.templates.push(MarkTemplateSpec {
            mark: mark.clone(),
            context,
        });
        self.templates.insert(key, index);
        Ok(index)
    }

    pub(super) fn context(stroke: &Stroke) -> MarkContext {
        MarkContext {
            units_per_pt: 1.0 / UNIT,
            line_thickness: stroke.width / UNIT,
            shaft_stroke: MarkStrokeStyle {
                cap: if stroke.round_cap {
                    Cap::Round
                } else {
                    Cap::Butt
                },
                ..MarkStrokeStyle::default()
            },
        }
    }

    pub(super) fn push(
        &mut self,
        template: usize,
        path: &BezPath,
        station: MarkStation,
        forward: bool,
    ) -> usize {
        let carrier = self.spec.carriers.len();
        self.spec
            .carriers
            .push(CurvePathOutput { path: path.clone() });
        self.place(template, carrier, station, forward)
    }

    /// Numeric engine stations center the painted head, including composites.
    pub(super) fn centered(
        &mut self,
        template: usize,
        path: &BezPath,
        ratio: f64,
        forward: bool,
    ) -> usize {
        self.push(template, path, MarkStation::Ratio { value: ratio }, forward)
    }

    pub(super) fn place(
        &mut self,
        template: usize,
        carrier: usize,
        station: MarkStation,
        forward: bool,
    ) -> usize {
        let index = self.spec.placements.len();
        self.spec.placements.push(MarkPlacement {
            template,
            carrier,
            station,
            shift: 0.0,
            direction: if forward {
                MarkDirection::Forward
            } else {
                MarkDirection::Backward
            },
        });
        index
    }

    pub(super) fn geometry(
        &self,
        mode: MarkGeometryMode,
    ) -> Result<MarkGeometryBatchOutput, String> {
        let mut spec = self.spec.clone();
        spec.mode = mode;
        Self::execute(spec)
    }

    pub(super) fn selected(&self, choices: &[usize]) -> Result<MarkGeometryBatchOutput, String> {
        let mut spec = self.spec.clone();
        spec.mode = MarkGeometryMode::Selected;
        spec.placements = choices
            .iter()
            .map(|&choice| self.spec.placements[choice].clone())
            .collect();
        Self::execute(spec)
    }

    fn execute(spec: MarkGeometrySpec) -> Result<MarkGeometryBatchOutput, String> {
        #[cfg(test)]
        GEOMETRY_CALLS.with(|calls| calls.borrow_mut().push(spec.mode));
        spec.geometry()
    }

    #[cfg(test)]
    pub(super) fn take_calls() -> Vec<MarkGeometryMode> {
        GEOMETRY_CALLS.with(|calls| std::mem::take(&mut *calls.borrow_mut()))
    }
}

impl Bounds {
    pub(super) fn painted(path: &BezPath) -> Self {
        let rect = path.bounding_box();
        Self {
            left: rect.x0,
            right: rect.x1,
            bottom: rect.y0,
            top: rect.y1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use kurbo::{PathEl, Point};
    use kurvst::marks::MarkFit;

    fn mark(shape: &str) -> MarkSpec {
        serde_json::from_value(serde_json::json!({"shape": shape})).unwrap()
    }

    fn curve() -> BezPath {
        let mut path = BezPath::new();
        path.move_to((0.0, 0.0));
        path.curve_to((0.0, 2.0), (3.0, 2.0), (3.0, 0.0));
        path
    }

    #[test]
    fn explicit_defaults_resolve_in_one_drawing_unit_conversion() {
        let context = Batch::context(&super::super::arrow_stroke());
        let flow = super::super::EdgeDrawing::default_flow_arrow();
        let momentum = super::super::EdgeDrawing::default_momentum_arrow();
        assert!(
            (flow.length.unwrap().resolve(context).unwrap() - 0.4409448818897638).abs() < 1e-14
        );
        assert!(
            (flow.width.unwrap().resolve(context).unwrap() - 0.33070866141732286).abs() < 1e-14
        );
        assert!(
            (momentum.length.unwrap().resolve(context).unwrap() - 0.3359580052493438).abs() < 1e-14
        );
        assert!(
            (momentum.width.unwrap().resolve(context).unwrap() - 0.25196850393700787).abs() < 1e-14
        );
        let prepared = flow.prepare(context).unwrap();
        assert!(prepared.paths()[0].fill);
        assert!(!prepared.paths()[0].stroke);
    }

    #[test]
    fn templates_are_deduplicated_by_geometry_and_line_context() {
        let mut batch = Batch::new();
        let mark = super::super::EdgeDrawing::default_momentum_arrow();
        let stroke = super::super::arrow_stroke();
        let first = batch.template(&mark, &stroke).unwrap();
        assert_eq!(batch.template(&mark.clone(), &stroke).unwrap(), first);
        let mut other_paint = stroke.clone();
        other_paint.paint = "red".into();
        assert_eq!(batch.template(&mark, &other_paint).unwrap(), first);
        other_paint.width = 2.0;
        assert_ne!(batch.template(&mark, &other_paint).unwrap(), first);
        assert_eq!(batch.spec.templates.len(), 2);
    }

    #[test]
    fn selected_batch_has_head_only_marks_and_one_authoritative_shaft() {
        let mut path = BezPath::new();
        path.move_to((0.0, 0.0));
        path.line_to((3.0, 0.0));
        let mut batch = Batch::new();
        let template = batch
            .template(&mark("triangle"), &super::super::arrow_stroke())
            .unwrap();
        batch.push(template, &path, MarkStation::End, true);
        batch.place(template, 0, MarkStation::Start, false);
        let result = batch.geometry(MarkGeometryMode::Selected).unwrap();
        assert_eq!(result.shafts.len(), 1);
        assert_eq!(result.marks.len(), 2);
        for head in &result.marks {
            assert!(head.shaft.path.is_empty() && head.shaft_outline.path.is_empty());
        }
        assert_eq!(
            result.shafts[0].shaft.path.elements()[0],
            PathEl::MoveTo(Point::from(result.marks[1].shaft_contact))
        );
        assert_eq!(
            result.shafts[0].shaft.path.elements()[1],
            PathEl::LineTo(Point::from(result.marks[0].shaft_contact))
        );
    }

    #[test]
    fn centered_flow_uses_the_visible_carrier_and_preserves_direction() {
        let path = BezPath::from_vec(vec![
            PathEl::MoveTo(Point::new(10.0, 0.0)),
            PathEl::LineTo(Point::new(13.0, 0.0)),
        ]);
        let mark = super::super::EdgeDrawing::default_flow_arrow();
        for forward in [true, false] {
            let mut batch = Batch::new();
            let template = batch
                .template(&mark, &super::super::arrow_stroke())
                .unwrap();
            batch.centered(template, &path, 0.5, forward);
            let result = batch.geometry(MarkGeometryMode::Selected).unwrap();
            let head = &result.marks[0];
            assert!(((head.tip.x + head.back.x) / 2.0 - 11.5).abs() < 1e-12);
            assert_eq!(head.tip.x > head.back.x, forward);
            assert_eq!(result.shafts[0].shaft.path, path);
        }
    }

    #[test]
    fn centered_catalogue_and_composites_use_their_actual_painted_extent() {
        let path = BezPath::from_vec(vec![
            PathEl::MoveTo(Point::new(10.0, 0.0)),
            PathEl::LineTo(Point::new(30.0, 0.0)),
        ]);
        let mut specs: Vec<_> = [
            "triangle", "straight", "stealth", "round", "tikz", "barb", "hooks", "bar", "bracket",
            "circle", "square", "diamond", "rays",
        ]
        .into_iter()
        .map(mark)
        .collect();
        specs.push(serde_json::from_value(serde_json::json!({
            "shape":"combine","parts":[{"shape":"circle"},{"gap":{"points":2,"ratio":0}},{"shape":"bar"}]
        })).unwrap());
        for mark in specs {
            for forward in [true, false] {
                let mut batch = Batch::new();
                let template = batch
                    .template(&mark, &super::super::arrow_stroke())
                    .unwrap();
                let index = batch.centered(template, &path, 0.5, forward);
                let geometry = batch.selected(&[index]).unwrap();
                let bounds = geometry.marks[0].footprint.path.bounding_box();
                assert!(
                    ((bounds.x0 + bounds.x1) / 2.0 - 20.0).abs() < 1e-12,
                    "{}",
                    mark.shape
                );
            }
        }
    }

    #[test]
    fn catalogue_and_combine_candidate_bounds_match_selected_paint_for_both_fits() {
        let mut specs: Vec<_> = [
            "triangle", "straight", "stealth", "round", "tikz", "barb", "hooks", "bar", "bracket",
            "circle", "square", "diamond", "rays",
        ]
        .into_iter()
        .map(mark)
        .collect();
        specs.push(serde_json::from_value(serde_json::json!({
            "shape": "combine", "parts": [{"shape":"bar"},{"gap":{"points":2,"ratio":0}},{"shape":"stealth"}]
        })).unwrap());
        for mut mark in specs {
            for fit in [MarkFit::Chord, MarkFit::Bend] {
                mark.fit = Some(fit);
                for path in [
                    curve(),
                    BezPath::from_vec(vec![
                        PathEl::MoveTo(Point::ZERO),
                        PathEl::LineTo(Point::new(0.001, 0.0)),
                    ]),
                ] {
                    let request = || {
                        let mut batch = Batch::new();
                        let template = batch
                            .template(&mark, &super::super::arrow_stroke())
                            .unwrap();
                        batch.push(template, &path, MarkStation::End, true);
                        batch
                    };
                    let candidate = request().geometry(MarkGeometryMode::Candidates).unwrap();
                    let selected = request().geometry(MarkGeometryMode::Selected).unwrap();
                    let bounds = Bounds::painted(&candidate.marks[0].footprint.path);
                    assert_eq!(candidate.marks[0].shaft.path, selected.shafts[0].shaft.path);
                    for (a, b) in candidate.marks[0]
                        .paths
                        .iter()
                        .zip(&selected.marks[0].paths)
                    {
                        assert_eq!(a.path, b.path);
                        let painted = Bounds::painted(&b.outline.path);
                        assert!(bounds.left <= painted.left && bounds.right >= painted.right);
                        assert!(bounds.bottom <= painted.bottom && bounds.top >= painted.top);
                    }
                }
            }
        }
    }
}
