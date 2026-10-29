use kurvst::marks::{MarkContext, MarkSize, MarkSpec, MarkStrokeStyle};
use std::f64::consts::TAU;
fn spec(shape: &str) -> MarkSpec {
    MarkSpec {shape: shape.into(), fit: None, shorten: None, length: None, width: None,
        inset: None, arc: None, phase: None, n: None, rev: None, align: None,
        fill: None, stroke: Some(true), parts: None}
}
fn size(points: f64) -> Option<MarkSize> { Some(MarkSize {points, ratio: 0.0}) }
fn context() -> MarkContext {
    MarkContext {units_per_pt: 1.0, line_thickness: 1.0, shaft_stroke: MarkStrokeStyle::default()}
}
fn right(spec: &MarkSpec) -> f64 {
    let p = spec.prepare(context()).unwrap();
    p.paths().iter().flat_map(|p| &p.outline_bounds).map(|b| b.right).fold(f64::NEG_INFINITY, f64::max)
}
#[test]
fn positive_area_deep_notch_must_not_paint_past_tip() {
    let mut s = spec("stealth"); s.length=size(10.0); s.width=size(8.0); s.inset=Some(0.999999);
    let front=right(&s); assert!(front <= 1e-8, "front={front}; a positive-area dart must grow/inset, not fall back to nominal vertices");
}
#[test]
fn tiny_positive_triangle_must_grow_to_hold_stroke() {
    let mut s = spec("triangle"); s.length=size(1e-20); s.width=size(1e-20);
    let front=right(&s); assert!(front <= 1e-8, "front={front}; finite bisection exhaustion must not masquerade as zero area");
}
#[test]
fn full_turn_barb_must_stop_at_rear_axis_crossing() {
    let mut s=spec("barb"); s.width=size(10.0); s.arc=Some(TAU);
    let p=s.prepare(context()).unwrap();
    assert!((p.end()-9.5).abs() < 1e-12,"end={} but the drawn circle centered at -5 with radius 4.5 meets the shaft first at -9.5",p.end());
}
#[test]
fn ordinary_ninety_percent_notch_is_a_control() {
    let mut s=spec("stealth"); s.length=size(10.0); s.width=size(8.0); s.inset=Some(0.9);
    assert!(right(&s).abs() < 1e-8);
}
