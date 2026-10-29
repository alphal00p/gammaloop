use kurvst::marks::{MarkContext, MarkSize, MarkSpec, MarkStrokeStyle};
use std::f64::consts::TAU;

fn spec(shape: &str, length: f64, width: f64, inset: Option<f64>) -> MarkSpec {
    let size = |points| Some(MarkSize { points, ratio: 0.0 });
    MarkSpec { shape: shape.into(), fit: None, shorten: None, length: size(length),
        width: size(width), inset, arc: None, phase: None, n: None, rev: None,
        align: None, fill: Some(false), stroke: Some(true), parts: None }
}
fn context(stroke: f64) -> MarkContext {
    MarkContext { units_per_pt: 1.0, line_thickness: stroke, shaft_stroke: MarkStrokeStyle::default() }
}
fn bounds(spec: &MarkSpec, stroke: f64) -> (f64, f64, f64) {
    let p = spec.prepare(context(stroke)).unwrap();
    let mut b = (f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for part in p.paths() {
        assert!(part.path.path.is_finite());
        assert!(part.outline.path.is_finite());
        for q in &part.outline_bounds { b.0 = b.0.min(q.left); b.1 = b.1.max(q.right); b.2 = b.2.max(q.top); }
    }
    b
}
#[test]
fn thin_positive_triangle_keeps_the_miter_tip() {
    let b = bounds(&spec("triangle", 10.0, 1e-7, None), 1.0);
    assert!(b.1.abs() < 1e-6, "nominal tip is 0, actual paint bounds={b:?}");
}
#[test]
fn deeper_positive_stealth_preserves_uniform_grown_wings() {
    let (l,w,i,s) = (10.0_f64,8.0_f64,0.99999999_f64,0.5_f64);
    let half=w/2.0;
    // Independent triangle inradius of one dart lobe, from area / semiperimeter.
    // All four inset edges collapse at this distance; uniform growth must retain
    // the outer wing vertices at (-grow*l, +/-grow*w/2).
    let limit=half*l*(1.0-i)/(l.hypot(half)+(l*i).hypot(half));
    let grow=(1.01*s/limit).max(1.0);
    let b=bounds(&spec("stealth",l,w,Some(i)),2.0*s);
    assert!((b.0+grow*l).abs() < 1e-4*grow*l, "grown left expected {}, actual bounds={b:?}",-grow*l);
    assert!((b.2-grow*half).abs() < 1e-4*grow*half, "grown top expected {}, actual bounds={b:?}",grow*half);
}
#[test]
fn nearby_triangle_control_still_reaches_tip() {
    let b=bounds(&spec("triangle",10.0,1e-5,None),1.0);
    assert!(b.1.abs()<1e-6,"{b:?}");
}
#[test]
fn nearby_stealth_control_preserves_grown_wings() {
    let (l,w,i,s)=(10.0_f64,8.0_f64,0.999999_f64,0.5_f64);
    let h=w/2.0; let limit=h*l*(1.0-i)/(l.hypot(h)+(l*i).hypot(h)); let g=(1.01*s/limit).max(1.0);
    let b=bounds(&spec("stealth",l,w,Some(i)),2.0*s);
    assert!((b.0+g*l).abs()<1e-4*g*l,"{b:?}");
}
#[test]
fn signed_exact_and_adjacent_barb_contacts_remain_correct() {
    for arc in [-TAU, -f64::from_bits(TAU.to_bits()-1), -0.0, 0.0, f64::from_bits(1),f64::from_bits(TAU.to_bits()-1),TAU] {
        for rev in [false,true] {
            let mut s=spec("barb",0.0,10.0,None);s.length=None;s.fill=None;s.arc=Some(arc);s.rev=Some(rev);
            let p=s.prepare(context(1.0)).unwrap();let expected=if rev{4.5}else if arc.abs()==TAU{9.5}else{0.5};
            assert!((p.end()-expected).abs()<1e-12,"{arc} {rev}: {}",p.end());
        }
    }
    for arc in [f64::from_bits(TAU.to_bits()+1),-f64::from_bits(TAU.to_bits()+1)] {
        let mut s=spec("barb",0.0,10.0,None);s.length=None;s.fill=None;s.arc=Some(arc);assert!(s.prepare(context(1.0)).is_err());
    }
}
