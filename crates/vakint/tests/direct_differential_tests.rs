//! Differential check of the PRE-EXISTING entries against the base they were
//! built on: a generator of v1 direct records (one and two loops; bubble,
//! sunset and inserted-bubble topologies; line powers 1 and 2; open slots
//! carried by the external momentum (ranks 0..8), by loop momenta (ranks
//! 0..4) or both; with and without a projector that closes every port;
//! scalar, gluon-like and closed-quark-loop numerators), single-scale vacuum
//! records with open slots and external vectors, and the in-tree B -> B DOT
//! graphs through the graph entry. Every case goes through
//! `FactorizedNumerator::from_text` (+ `validate`), `vakint_input`,
//! `evaluate_factorized` (or `evaluate_factorized_vacuum`, `evaluate_graph`),
//! and the text round trip `parse_text` -> `to_text`, and is reduced to rows
//! `(case, stage, ok | err, canonical text)`. The reviewer's bubble with 0..7
//! external spectators (a vendored v1 text) and the audit's witnesses (vacuum
//! values with a projector, unsorted ports, a port only in the projector,
//! amputated fermion legs, compact slashes, absent spinor ports) are cases.
//!
//! This file uses only API the base already had, so it compiles on both
//! sides: with `VAKINT_DIFF_RECORD=<path>` it writes the rows (run on the
//! base), otherwise it compares the current build with the vendored rows of
//! the base (tests/fixtures/direct_differential/base.tsv). A value without a
//! new-feature marker is validated and evaluated by exactly the base's rules,
//! so every row must be IDENTICAL: accept / refuse per stage, the refusal
//! message (variant and text), the canonical text of the value, of the vakint
//! input (raw indices) and of the written record. Any difference fails.
#![cfg(feature = "feynkit-ingress")]

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;

use feynkit_graph::FeynmanDiagram;
use symbolica::atom::{Atom, AtomCore};
use vakint::Vakint;
use vakint::p_integral::factorized::{
    FromFeynkitOptions, ProjectorChoice, evaluate_factorized, evaluate_factorized_vacuum,
    evaluate_graph, vakint_input,
};
use vakint::p_integral::{FactorizedNumerator, FactorizedParts};
use vakint::refpath;

const DIM: &str = "gammalooprs::{}::dim";
const MINK: &str = "spenso::{spenso::representation,spenso::self_dual}::mink";
const BIS: &str = "spenso::{spenso::representation,spenso::self_dual}::bis";
const HEDGE: &str = "gammalooprs::{spenso::index}::hedge";
const EDGE: &str = "gammalooprs::{spenso::index,spenso::index-label:e}::edge";
const Q: &str = "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q";
const G: &str = "spenso::{symmetric,linear,real,spenso::tensor}::g";
const GAMMA: &str = "spenso::{linear,spenso::tensor}::gamma";

fn mink(index: &str) -> String {
    format!("{MINK}({DIM},{index})")
}

fn hedge(h: usize, k: usize) -> String {
    format!("{HEDGE}({h},{k})")
}

fn q(edge: usize, slot: &str) -> String {
    format!("{Q}({edge},{slot})")
}

fn product(factors: &[String]) -> String {
    if factors.is_empty() {
        "1".into()
    } else {
        factors.join("*")
    }
}

/// One line of a v1 record: `edge` with its factor.
struct Edge {
    id: usize,
    source: String,
    target: String,
    hedges: (Option<usize>, Option<usize>),
    particle: &'static str,
    line: &'static str,
    power: i64,
    momentum: Vec<i64>,
    mass: String,
    factor: String,
}

struct Record {
    label: String,
    loops: usize,
    externals: usize,
    loop_edges: Vec<usize>,
    external_edge: Option<usize>,
    edges: Vec<Edge>,
    vertices: Vec<(usize, &'static str, String)>,
    ports: Vec<String>,
    projector: Option<(i64, String)>,
    loop_sign: i64,
    nf_power: u32,
}

impl Record {
    fn text(&self) -> String {
        let mut l = vec![
            "schema\tpint6.factorized-numerator.v1".to_string(),
            format!("label\t{}", self.label),
            "kind\tdirect".into(),
            format!("dimension\t{DIM}"),
            format!(
                "momenta\tloops={}\texternals={}",
                self.loops, self.externals
            ),
        ];
        for (i, e) in self.loop_edges.iter().enumerate() {
            l.push(format!("loop\t{}\tedge={e}", i + 1));
        }
        for j in 0..self.externals {
            l.push(format!(
                "external\t{}\tedge={}",
                j + 1,
                self.external_edge
                    .map_or("none".to_string(), |e| e.to_string())
            ));
        }
        l.push(format!("edges\t{}", self.edges.len()));
        for e in &self.edges {
            let h = |h: Option<usize>| h.map_or("-".to_string(), |h| h.to_string());
            l.push(format!(
                "edge\t{}\tsource={}\ttarget={}\thedges={},{}\tparticle={}\tline={}\tpower={}\tmomentum={}\t{}",
                e.id,
                e.source,
                e.target,
                h(e.hedges.0),
                h(e.hedges.1),
                e.particle,
                e.line,
                e.power,
                e.momentum
                    .iter()
                    .map(|c| c.to_string())
                    .collect::<Vec<_>>()
                    .join(","),
                e.mass
            ));
            l.push(format!("edge_factor\t{}\t{}", e.id, e.factor));
        }
        l.push(format!("vertices\t{}", self.vertices.len()));
        for (id, rule, factor) in &self.vertices {
            l.push(format!("vertex\t{id}\trule={rule}\t{factor}"));
        }
        l.push(format!("ports\t{}", self.ports.len()));
        for p in &self.ports {
            l.push(format!("port\t{p}"));
        }
        match &self.projector {
            None => l.push("projector\tnone".into()),
            Some((p2, t)) => l.push(format!("projector\tp2_power={p2}\t{t}")),
        }
        l.push(format!(
            "prefactor\tloop_sign={}\tnf_power={}",
            self.loop_sign, self.nf_power
        ));
        l.push("prefactor_overall\t1".into());
        l.push("prefactor_scalars\t1".into());
        l.push("prefactor_colour\t1".into());
        l.push("colour_branch\tnone".into());
        l.push("insertions\t0".into());
        l.push(format!("lines\t{}", l.len() + 1));
        l.push("end".into());
        let mut t = l.join("\n");
        t.push('\n');
        t
    }
}

fn scalar_edge(
    id: usize,
    s: &str,
    t: &str,
    hs: (Option<usize>, Option<usize>),
    m: Vec<i64>,
    power: i64,
) -> Edge {
    Edge {
        id,
        source: s.into(),
        target: t.into(),
        hedges: hs,
        particle: "phi",
        line: "scalar",
        power,
        momentum: m,
        mass: "0".into(),
        factor: "1".into(),
    }
}

/// A propagator topology: the edges (scalar lines), the vertex of the
/// incoming leg (edge 0, half-edge 0) and the loop edges.
struct Topology {
    name: &'static str,
    loops: usize,
    edges: Vec<Edge>,
    vertices: Vec<usize>,
    loop_edges: Vec<usize>,
    /// The internal edge whose (loop) momentum carries loop slots.
    loop_carrier: usize,
}

fn topologies(power: i64) -> Vec<Topology> {
    let s = |v: usize| v.to_string();
    vec![
        // one loop: k, -k-p (as the reviewer's witness)
        Topology {
            name: "bubble",
            loops: 1,
            edges: vec![
                scalar_edge(0, "ext", &s(1), (None, Some(0)), vec![0, 1], 1),
                scalar_edge(1, &s(0), "ext", (Some(1), None), vec![0, 1], 1),
                scalar_edge(2, &s(0), &s(1), (Some(2), Some(3)), vec![1, 0], power),
                scalar_edge(3, &s(0), &s(1), (Some(4), Some(5)), vec![-1, -1], 1),
            ],
            vertices: vec![0, 1],
            loop_edges: vec![2],
            loop_carrier: 2,
        },
        // two loops: k1, k2, -k1-k2-p
        Topology {
            name: "sunset",
            loops: 2,
            edges: vec![
                scalar_edge(0, "ext", &s(1), (None, Some(0)), vec![0, 0, 1], 1),
                scalar_edge(1, &s(0), "ext", (Some(1), None), vec![0, 0, 1], 1),
                scalar_edge(2, &s(0), &s(1), (Some(2), Some(3)), vec![1, 0, 0], power),
                scalar_edge(3, &s(0), &s(1), (Some(4), Some(5)), vec![0, 1, 0], 1),
                scalar_edge(4, &s(0), &s(1), (Some(6), Some(7)), vec![-1, -1, -1], 1),
            ],
            vertices: vec![0, 1],
            loop_edges: vec![2, 3],
            loop_carrier: 2,
        },
        // two loops: a one-loop bubble inserted on the k1 line
        Topology {
            name: "inserted",
            loops: 2,
            edges: vec![
                scalar_edge(0, "ext", &s(0), (None, Some(0)), vec![0, 0, 1], 1),
                scalar_edge(1, &s(1), "ext", (Some(1), None), vec![0, 0, 1], 1),
                scalar_edge(2, &s(0), &s(2), (Some(2), Some(3)), vec![1, 0, 0], power),
                scalar_edge(3, &s(2), &s(3), (Some(4), Some(5)), vec![0, 1, 0], 1),
                scalar_edge(4, &s(2), &s(3), (Some(6), Some(7)), vec![1, -1, 0], 1),
                scalar_edge(5, &s(3), &s(1), (Some(8), Some(9)), vec![1, 0, 0], 1),
                scalar_edge(6, &s(0), &s(1), (Some(10), Some(11)), vec![-1, 0, 1], 1),
            ],
            vertices: vec![0, 1, 2, 3],
            loop_edges: vec![2, 3],
            loop_carrier: 3,
        },
    ]
}

/// The vertex that holds half-edge 0 (the incoming leg's).
fn leg_vertex(t: &Topology) -> usize {
    t.edges[0].target.parse().unwrap()
}

/// Scalar-line records with `ext` external-momentum slots and `lp` loop
/// slots at the incoming leg's vertex, optionally closed by a projector.
fn scalar_cases(out: &mut Vec<(String, String)>) {
    for power in [1, 2] {
        for t in topologies(power) {
            let max_loop = if t.loops == 1 { 4 } else { 3 };
            for ext in 0..=8usize {
                for lp in 0..=max_loop {
                    if ext + lp > 8 || (ext > 4 && lp > 0) {
                        continue;
                    }
                    for projected in [false, true] {
                        if projected && (ext + lp == 0 || ext + lp > 4) {
                            continue;
                        }
                        let v = leg_vertex(&t);
                        let mut ports = Vec::new();
                        let mut factors = Vec::new();
                        for k in 1..=ext {
                            let slot = mink(&hedge(0, k));
                            factors.push(q(0, &slot));
                            ports.push(slot);
                        }
                        for k in 1..=lp {
                            let slot = mink(&hedge(0, 100 + k));
                            factors.push(q(t.loop_carrier, &slot));
                            ports.push(slot);
                        }
                        ports.sort();
                        // a projector closing every port: metrics on pairs,
                        // the external momentum on a last odd port
                        let projector = projected.then(|| {
                            let mut tensor = Vec::new();
                            let mut i = 0;
                            while i + 1 < ports.len() {
                                tensor.push(format!("{G}({},{})", ports[i], ports[i + 1]));
                                i += 2;
                            }
                            if i < ports.len() {
                                tensor.push(q(0, &ports[i]));
                            }
                            (-1, product(&tensor))
                        });
                        let mut vertices = t
                            .vertices
                            .iter()
                            .map(|id| (*id, "V", "1".to_string()))
                            .collect::<Vec<_>>();
                        let at = vertices.iter().position(|x| x.0 == v).unwrap();
                        vertices[at].2 = product(&factors);
                        let edges = topologies(power)
                            .into_iter()
                            .find(|x| x.name == t.name)
                            .unwrap()
                            .edges;
                        let record = Record {
                            label: format!(
                                "{}-p{power}-ext{ext}-loop{lp}-{}",
                                t.name,
                                if projected { "closed" } else { "open" }
                            ),
                            loops: t.loops,
                            externals: 1,
                            loop_edges: t.loop_edges.clone(),
                            external_edge: Some(0),
                            edges,
                            vertices,
                            ports,
                            projector,
                            loop_sign: 1,
                            nf_power: 0,
                        };
                        out.push((record.label.clone(), record.text()));
                    }
                }
            }
        }
    }
}

/// One-loop quark bubble (closed fermion loop, two open Lorentz ports) and
/// gluon-like bubble (metrics on the lines); with and without projector.
fn line_cases(out: &mut Vec<(String, String)>) {
    for power in [1, 2] {
        for projected in [false, true] {
            let ports = vec![mink(&hedge(0, 1)), mink(&hedge(1, 1))];
            let projector = projected.then(|| (-1, format!("{G}({},{})", ports[0], ports[1])));
            let b = |h: usize| format!("{BIS}(4,{})", hedge(h, 1));
            let m = |e: usize| mink(&format!("{EDGE}({e},1)"));
            // quark loop: e2 v0 -> v1 (k), e3 v1 -> v0 (k + p)
            let edges = vec![
                scalar_edge(0, "ext", "1", (None, Some(0)), vec![0, 1], 1),
                scalar_edge(1, "0", "ext", (Some(1), None), vec![0, 1], 1),
                Edge {
                    line: "fermion",
                    particle: "q",
                    factor: format!("{}*{GAMMA}({},{},{})", q(2, &m(2)), b(3), b(2), m(2)),
                    ..scalar_edge(2, "0", "1", (Some(2), Some(3)), vec![1, 0], power)
                },
                Edge {
                    line: "fermion",
                    particle: "q",
                    factor: format!("{}*{GAMMA}({},{},{})", q(3, &m(3)), b(5), b(4), m(3)),
                    ..scalar_edge(3, "1", "0", (Some(4), Some(5)), vec![1, 1], 1)
                },
            ];
            let vertices = vec![
                (
                    0,
                    "V_qqB",
                    format!("{GAMMA}({},{},{})", b(2), b(5), ports[1]),
                ),
                (
                    1,
                    "V_qqB",
                    format!("{GAMMA}({},{},{})", b(4), b(3), ports[0]),
                ),
            ];
            let record = Record {
                label: format!(
                    "quark-bubble-p{power}-{}",
                    if projected { "closed" } else { "open" }
                ),
                loops: 1,
                externals: 1,
                loop_edges: vec![2],
                external_edge: Some(0),
                edges,
                vertices,
                ports: ports.clone(),
                projector: projector.clone(),
                loop_sign: -1,
                nf_power: 1,
            };
            out.push((record.label.clone(), record.text()));
            // gluon-like: vector lines with metrics, vertices g x momentum
            let mh = |h: usize| mink(&hedge(h, 1));
            let edges = vec![
                scalar_edge(0, "ext", "1", (None, Some(0)), vec![0, 1], 1),
                scalar_edge(1, "0", "ext", (Some(1), None), vec![0, 1], 1),
                Edge {
                    line: "vector",
                    particle: "g",
                    factor: format!("{G}({},{})", mh(2), mh(3)),
                    ..scalar_edge(2, "0", "1", (Some(2), Some(3)), vec![1, 0], power)
                },
                Edge {
                    line: "vector",
                    particle: "g",
                    factor: format!("{G}({},{})", mh(4), mh(5)),
                    ..scalar_edge(3, "0", "1", (Some(4), Some(5)), vec![-1, -1], 1)
                },
            ];
            let vertices = vec![
                (
                    0,
                    "V3",
                    format!("{}*{G}({},{})", q(2, &mh(1)), mh(2), mh(4)),
                ),
                (
                    1,
                    "V3",
                    format!("{}*{G}({},{})", q(3, &mh(0)), mh(3), mh(5)),
                ),
            ];
            let record = Record {
                label: format!(
                    "gluon-bubble-p{power}-{}",
                    if projected { "closed" } else { "open" }
                ),
                loops: 1,
                externals: 1,
                loop_edges: vec![2],
                external_edge: Some(0),
                edges,
                vertices,
                ports,
                projector,
                loop_sign: 1,
                nf_power: 0,
            };
            out.push((record.label.clone(), record.text()));
        }
    }
}

/// Edge cases at the boundary of what the base accepted.
fn edge_cases(out: &mut Vec<(String, String)>) {
    let base =
        |label: &str, factor: String, ports: Vec<String>, projector: Option<(i64, String)>| {
            Record {
                label: label.into(),
                loops: 1,
                externals: 1,
                loop_edges: vec![2],
                external_edge: Some(0),
                edges: topologies(1).into_iter().next().unwrap().edges,
                vertices: vec![(0, "V", "1".into()), (1, "V", factor)],
                ports,
                projector,
                loop_sign: 1,
                nf_power: 0,
            }
        };
    // an integer port index (the base mapped it)
    let slot = format!("{MINK}({DIM},7)");
    let r = base("edge-integer-index", q(2, &slot), vec![slot], None);
    out.push((r.label.clone(), r.text()));
    // spinor slots on a vertex without a fermion line (refused on both)
    let b1 = format!("{BIS}(4,{})", hedge(0, 1));
    let b2 = format!("{BIS}(4,{})", hedge(0, 2));
    let r = base(
        "edge-spinor-slot-no-line",
        format!(
            "{GAMMA}({b1},{b2},{})*{}",
            mink(&format!("{EDGE}(2,1)")),
            q(2, &mink(&format!("{EDGE}(2,1)")))
        ),
        vec![b1.clone(), b2.clone()],
        None,
    );
    out.push((r.label.clone(), r.text()));
    // a one-loop quark self-energy with amputated quark legs whose spinor
    // ports a gamma projector closes (the base refused it)
    let bh = |h: usize| format!("{BIS}(4,{})", hedge(h, 1));
    let mh = |h: usize| mink(&hedge(h, 1));
    let me = mink(&format!("{EDGE}(2,1)"));
    let closing = mink(&hedge(0, 9));
    let quark = |id, s: &str, t: &str, hs, m, factor: String| Edge {
        line: "fermion",
        particle: "q",
        factor,
        ..scalar_edge(id, s, t, hs, m, 1)
    };
    let self_energy = |label: &str, projector: Option<(i64, String)>| Record {
        label: label.into(),
        loops: 1,
        externals: 1,
        loop_edges: vec![2],
        external_edge: Some(0),
        edges: vec![
            quark(0, "ext", "1", (None, Some(0)), vec![0, 1], "1".into()),
            quark(1, "0", "ext", (Some(1), None), vec![0, 1], "1".into()),
            quark(
                2,
                "1",
                "0",
                (Some(2), Some(3)),
                vec![1, 0],
                format!("{}*{GAMMA}({},{},{me})", q(2, &me), bh(3), bh(2)),
            ),
            Edge {
                line: "vector",
                particle: "g",
                factor: format!("{G}({},{})", mh(4), mh(5)),
                ..scalar_edge(3, "1", "0", (Some(4), Some(5)), vec![-1, 1], 1)
            },
        ],
        vertices: vec![
            (
                0,
                "V_qqg",
                format!("{GAMMA}({},{},{})", bh(1), bh(3), mh(5)),
            ),
            (
                1,
                "V_qqg",
                format!("{GAMMA}({},{},{})", bh(2), bh(0), mh(4)),
            ),
        ],
        ports: vec![bh(0), bh(1)],
        projector,
        loop_sign: 1,
        nf_power: 0,
    };
    let r = self_energy(
        "v1-gamma-projector",
        Some((
            0,
            format!("{GAMMA}({},{},{closing})*{}", bh(0), bh(1), q(0, &closing)),
        )),
    );
    out.push((r.label.clone(), r.text()));
    // the same with a COMPACT slash in the projector (audit A-2 (c), A-4: the
    // base kept the text as written and refused the projector)
    let r = self_energy(
        "v1-compact-slash-projector",
        Some((
            0,
            format!(
                "{GAMMA}({},{},{})",
                bh(0),
                bh(1),
                q(0, &format!("{MINK}({DIM})"))
            ),
        )),
    );
    out.push((r.label.clone(), r.text()));
    // amputated fermion legs without a projector (audit A-3)
    let r = self_energy("v1-amputated-legs-open", None);
    out.push((r.label.clone(), r.text()));
    // a projector closing one of two ports (the base refused it)
    let (p1, p2) = (mink(&hedge(0, 1)), mink(&hedge(0, 2)));
    let r = base(
        "v1-partial-projector",
        format!("{}*{}", q(2, &p1), q(2, &p2)),
        vec![p1.clone(), p2.clone()],
        Some((0, q(0, &p1))),
    );
    out.push((r.label.clone(), r.text()));
    // a soft leg with the smallest id (the base refused it): edge 0 at zero
    // momentum, p in through a new edge 9
    let mut r = base(
        "v1-soft-first-leg",
        q(2, &mink(&hedge(0, 1))),
        vec![mink(&hedge(0, 1))],
        None,
    );
    r.edges[0].momentum = vec![0, 0];
    r.edges
        .push(scalar_edge(9, "ext", "1", (None, Some(9)), vec![0, 1], 1));
    r.edges.sort_by_key(|e| e.id);
    r.external_edge = Some(9);
    out.push((r.label.clone(), r.text()));
    // a scalar product of two compact momenta in a factor (the base took it
    // for tensor content)
    let dot = format!(
        "spenso::{{symmetric,linear,scalar}}::dot({},{})",
        q(2, &format!("{MINK}({DIM})")),
        q(3, &format!("{MINK}({DIM})"))
    );
    let r = base("v1-dot-in-factor", dot, Vec::new(), None);
    out.push((r.label.clone(), r.text()));
    // ports listed out of canonical order
    let (pa, pb) = (mink(&hedge(0, 2)), mink(&hedge(0, 1)));
    let r = base(
        "edge-unsorted-ports",
        format!("{}*{}", q(2, &pa), q(0, &pb)),
        vec![pa.clone(), pb.clone()],
        None,
    );
    out.push((r.label.clone(), r.text()));
    // a port that occurs only in the projector, never in a factor (audit A-1)
    let absent = mink(&hedge(99, 0));
    let r = base(
        "edge-port-only-in-projector",
        q(2, &pb),
        vec![pb.clone(), absent.clone()],
        Some((0, format!("{}*{}", q(0, &pb), q(0, &absent)))),
    );
    out.push((r.label.clone(), r.text()));
    // spinor ports on a half-edge the graph lacks, no projector (audit A-6)
    let (s0, s1) = (
        format!("{BIS}(4,{})", hedge(99, 0)),
        format!("{BIS}(4,{})", hedge(99, 1)),
    );
    let r = base(
        "edge-absent-spinor-ports",
        format!("{GAMMA}({s0},{s1},{pb})"),
        vec![pb.clone(), s0, s1],
        None,
    );
    out.push((r.label.clone(), r.text()));
    // the reviewer's bubble with 0..7 external spectators (a vendored v1
    // text; ranks below 7 drop the last factors and ports)
    for rank in 0..=7 {
        out.push((format!("spectators-rank{rank}"), spectators(rank)));
    }
}

/// The reviewer's v1 bubble with `rank` external spectators
/// `Q(0, mink(D, hedge(0, i)))`, i = 1..rank (vendored texts; 0, 6 and 7 are
/// the reviewer's files, the others drop the last factors and ports of 7).
fn spectators(rank: usize) -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!(
        "tests/fixtures/direct_differential/rank{rank}-v1.txt"
    ));
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// Single-scale vacuum records: a massive one-loop bubble with two open
/// slots (in the given and in the reversed order), and the same closed by a
/// projector of two constant external vectors.
fn vacuum_cases(out: &mut Vec<(String, String)>) {
    let m2 = "vakint_diff::{}::msq";
    for (label, order, projected) in [
        ("vacuum-open", false, false),
        ("vacuum-open-reversed", true, false),
        ("vacuum-closed-vectors", false, true),
    ] {
        let (a, b) = (mink(&hedge(0, 1)), mink(&hedge(0, 2)));
        let mut ports = vec![a.clone(), b.clone()];
        if order {
            ports.reverse();
        }
        let massive = |id, s: &str, t: &str, hs, m| Edge {
            mass: m2.into(),
            ..scalar_edge(id, s, t, hs, m, 1)
        };
        let record = Record {
            label: label.into(),
            loops: 1,
            externals: 0,
            loop_edges: vec![0],
            external_edge: None,
            edges: vec![
                massive(0, "0", "1", (Some(0), Some(1)), vec![1]),
                massive(1, "1", "0", (Some(2), Some(3)), vec![1]),
            ],
            vertices: vec![
                (0, "V", format!("{}*{}", q(0, &a), q(0, &b))),
                (1, "V", "1".into()),
            ],
            ports,
            projector: projected.then(|| (0, format!("{}*{}", q(98, &a), q(99, &b)))),
            loop_sign: 1,
            nf_power: 0,
        };
        out.push((label.to_string(), record.text()));
    }
}

/// Vacuum values the vacuum entry takes WITHOUT the contract's checks (read
/// by `parse_text`, as a caller hands a value over): the audit's A-9
/// witnesses: (a) a projector of the loop momentum on both ports, (b) the
/// same with the edges in descending id order, (c) a spinor metric in the
/// projector closing a trace of the factors; and two ports on distinct
/// constant vectors in both orders (A-10).
fn vacuum_raw_cases(out: &mut Vec<(String, String)>) {
    let m2 = "vakint_diff::{}::msq";
    let (a, b) = (mink(&hedge(0, 1)), mink(&hedge(0, 2)));
    let massive = |id, s: &str, t: &str, hs, m| Edge {
        mass: m2.into(),
        ..scalar_edge(id, s, t, hs, m, 1)
    };
    let record = |label: &str,
                  factor: String,
                  ports: Vec<String>,
                  projector: Option<String>,
                  descending: bool| {
        let mut edges = vec![
            massive(0, "0", "1", (Some(0), Some(1)), vec![1]),
            massive(1, "1", "0", (Some(2), Some(3)), vec![1]),
        ];
        if descending {
            edges.reverse();
        }
        Record {
            label: label.into(),
            loops: 1,
            externals: 0,
            loop_edges: vec![0],
            external_edge: None,
            edges,
            vertices: vec![(0, "V", factor), (1, "V", "1".into())],
            ports,
            projector: projector.map(|t| (0, t)),
            loop_sign: 1,
            nf_power: 0,
        }
    };
    let loop_projector = format!("{}*{}", q(0, &a), q(0, &b));
    for (label, descending) in [
        ("vacuum-raw-loop-projector", false),
        ("vacuum-raw-loop-projector-descending", true),
    ] {
        let r = record(
            label,
            format!("{}*{}", q(0, &a), q(0, &b)),
            vec![a.clone(), b.clone()],
            Some(loop_projector.clone()),
            descending,
        );
        out.push((label.to_string(), r.text()));
    }
    let s = |k: usize| format!("{BIS}(4,{})", hedge(0, 10 + k));
    let (mu, nu) = (mink(&hedge(0, 21)), mink(&hedge(0, 22)));
    let trace = format!(
        "{GAMMA}({},{},{mu})*{GAMMA}({},{},{nu})*{}*{}",
        s(1),
        s(2),
        s(2),
        s(3),
        q(0, &mu),
        q(0, &nu)
    );
    let r = record(
        "vacuum-raw-spinor-metric-projector",
        trace.clone(),
        Vec::new(),
        Some(format!("{G}({},{})", s(3), s(1))),
        false,
    );
    out.push((r.label.clone(), r.text()));
    // the same trace opened at a bispinor port pair the projector closes
    let r = record(
        "vacuum-raw-bispinor-ports-closed",
        trace,
        vec![s(1), s(3)],
        Some(format!("{G}({},{})", s(3), s(1))),
        false,
    );
    out.push((r.label.clone(), r.text()));
    // two open ports carried by two different constant vectors, no projector,
    // listed in and against canonical order (audit A-10: the raw index
    // PORT + i follows the stored order)
    for (label, reversed) in [
        ("vacuum-raw-two-vectors", false),
        ("vacuum-raw-two-vectors-reversed", true),
    ] {
        let mut ports = vec![a.clone(), b.clone()];
        if reversed {
            ports.reverse();
        }
        let r = record(
            label,
            format!("{}*{}", q(98, &a), q(99, &b)),
            ports,
            None,
            false,
        );
        out.push((label.to_string(), r.text()));
    }
}

fn settings(loops: usize) -> vakint::VakintSettings {
    refpath::vakint_settings(loops as i64 + 1)
}

fn row(case: &str, stage: &str, r: std::result::Result<String, String>) -> String {
    let one = |s: &str| s.replace(['\n', '\t'], " ");
    match r {
        Ok(v) => format!("{case}\t{stage}\tok\t{}", one(&v)),
        Err(e) => format!("{case}\t{stage}\terr\t{}", one(&e)),
    }
}

/// Every case through the old entries: the rows.
fn run_all() -> Vec<String> {
    Vakint::initialize_vakint_symbols();
    let vakint = Vakint::new().unwrap();
    let mut cases = Vec::new();
    scalar_cases(&mut cases);
    line_cases(&mut cases);
    edge_cases(&mut cases);
    let mut vacuum = Vec::new();
    vacuum_cases(&mut vacuum);
    let mut rows = Vec::new();
    for (case, text) in cases.iter().chain(&vacuum) {
        rows.push(row(
            case,
            "roundtrip",
            FactorizedNumerator::parse_text(text)
                .and_then(|v| v.to_text())
                .map_err(|e| e.to_string()),
        ));
    }
    // from_parts of the value as read (the base built it unchanged)
    for (case, text) in &cases {
        let Ok(v) = FactorizedNumerator::parse_text(text) else {
            continue;
        };
        let built = FactorizedNumerator::from_parts(FactorizedParts {
            label: v.label.clone(),
            dimension: v.dimension,
            loop_edges: v.loop_edges.clone(),
            external_edges: v.external_edges.clone(),
            edges: v.edges.clone(),
            vertices: v.vertices.clone(),
            ports: v.ports.clone(),
            projector: v.projector.clone(),
            prefactor: v.prefactor.clone(),
        });
        rows.push(row(
            case,
            "from-parts",
            built.and_then(|b| b.to_text()).map_err(|e| e.to_string()),
        ));
    }
    for (case, text) in &cases {
        let value = FactorizedNumerator::from_text(text);
        rows.push(row(
            case,
            "validate",
            value
                .as_ref()
                .map(|_| "valid".into())
                .map_err(|e| e.to_string()),
        ));
        let Ok(value) = value else { continue };
        let s = settings(value.loops());
        let prefactor = value.prefactor.total();
        rows.push(row(
            case,
            "input",
            vakint_input(&s, &value, &prefactor)
                .map(|(i, _)| i.to_canonical_string())
                .map_err(|e| e.to_string()),
        ));
        rows.push(row(
            case,
            "value",
            evaluate_factorized(&vakint, &s, &value, &prefactor)
                .map(|v| v.value.to_canonical_string())
                .map_err(|e| e.to_string()),
        ));
    }
    let mut raw = Vec::new();
    vacuum_raw_cases(&mut raw);
    let vacuum = vacuum
        .iter()
        .map(|(case, text)| (case, text, true))
        .chain(raw.iter().map(|(case, text)| (case, text, false)));
    for (case, text, contract) in vacuum {
        // the contract's checks for the plain records; the raw values go to
        // the vacuum entry as read
        let value = if contract {
            FactorizedNumerator::from_text(text)
        } else {
            FactorizedNumerator::parse_text(text)
        };
        let stage = if contract { "validate" } else { "parse" };
        rows.push(row(
            case,
            stage,
            value
                .as_ref()
                .map(|_| "valid".into())
                .map_err(|e| e.to_string()),
        ));
        let Ok(value) = value else { continue };
        let s = settings(value.loops());
        let got = evaluate_factorized_vacuum(&vakint, &s, &value, &Atom::num(1));
        rows.push(row(
            case,
            "vacuum-value",
            got.as_ref()
                .map(|v| v.value.to_canonical_string())
                .map_err(|e| e.to_string()),
        ));
        rows.push(row(
            case,
            "vacuum-vakint-value",
            got.map(|v| v.vakint_value.to_canonical_string())
                .map_err(|e| e.to_string()),
        ));
    }
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/bb");
    for level in ["L1", "L2"] {
        let directory = root.join(level);
        let model = refpath::load_model(&directory).unwrap();
        let mut files = std::fs::read_dir(&directory)
            .unwrap()
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.extension().is_some_and(|x| x == "dot"))
            .collect::<Vec<_>>();
        files.sort();
        for file in files {
            let dot = std::fs::read_to_string(&file).unwrap();
            let diagram = FeynmanDiagram::from_dot(Arc::clone(&model), &dot).unwrap();
            let s = refpath::vakint_settings(diagram.loop_count() as i64 + 1);
            for (name, projector) in [
                ("metric", ProjectorChoice::Metric),
                ("open", ProjectorChoice::None),
            ] {
                let case = format!(
                    "graph-{level}-{}-{name}",
                    file.file_stem().unwrap().to_string_lossy()
                );
                let options = FromFeynkitOptions {
                    projector,
                    ..FromFeynkitOptions::default()
                };
                rows.push(row(
                    &case,
                    "graph-value",
                    evaluate_graph(&vakint, &s, &diagram, &model, &options)
                        .map(|v| v.value.to_canonical_string())
                        .map_err(|e| e.to_string()),
                ));
            }
        }
    }
    rows
}

#[derive(Debug, PartialEq)]
enum Difference {
    /// the base accepted, the current build refuses
    Narrower(String),
    /// the base refused, the current build accepts
    Wider(String),
    /// both accept, the texts differ (`algebraic`: the difference expands
    /// to zero, i.e. only the form differs)
    Value { name: String, algebraic: bool },
    /// both refuse with other messages (another variant or text)
    Message(String),
    /// a row of one side is missing on the other
    Missing(String),
}

fn parse(text: &str) -> Atom {
    symbolica::try_parse!(text).unwrap_or_else(|e| panic!("{text}: {e}"))
}

/// The differences between the base rows and the current rows.
fn compare(base: &[String], current: &[String]) -> Vec<Difference> {
    let key = |r: &String| {
        let mut f = r.splitn(4, '\t');
        let (c, s, st, v) = (
            f.next().unwrap(),
            f.next().unwrap(),
            f.next().unwrap(),
            f.next().unwrap_or(""),
        );
        (
            (c.to_string(), s.to_string()),
            (st.to_string(), v.to_string()),
        )
    };
    let b = base.iter().map(key).collect::<BTreeMap<_, _>>();
    let c = current.iter().map(key).collect::<BTreeMap<_, _>>();
    let mut out = Vec::new();
    for (k, (bs, bv)) in &b {
        let name = format!("{} / {}", k.0, k.1);
        let Some((cs, cv)) = c.get(k) else {
            // a refusal at validate stops the later stages on that side
            out.push(Difference::Missing(name));
            continue;
        };
        match (bs.as_str(), cs.as_str()) {
            ("ok", "ok") => {
                if bv != cv {
                    // the written record is not an atom
                    let algebraic =
                        k.1 != "roundtrip" && (&parse(bv) - &parse(cv)).expand().is_zero();
                    out.push(Difference::Value {
                        name: format!("{name}: base {bv} / now {cv}"),
                        algebraic,
                    });
                }
            }
            ("ok", _) => out.push(Difference::Narrower(format!(
                "{name}: base {bv} / now refused: {cv}"
            ))),
            (_, "ok") => out.push(Difference::Wider(format!(
                "{name}: base refused: {bv} / now {cv}"
            ))),
            _ => {
                if bv != cv {
                    out.push(Difference::Message(format!("{name}: base {bv} / now {cv}")));
                }
            }
        }
    }
    for k in c.keys() {
        if !b.contains_key(k) {
            out.push(Difference::Missing(format!("{} / {} (only now)", k.0, k.1)));
        }
    }
    out
}

fn vendored() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/direct_differential/base.tsv")
}

/// The old entries against the base, case by case (module documentation).
#[test]
fn the_old_entries_behave_as_on_the_base() {
    let rows = run_all();
    if let Some(path) = std::env::var_os("VAKINT_DIFF_RECORD") {
        let mut text = rows.join("\n");
        text.push('\n');
        std::fs::write(&path, text).unwrap();
        eprintln!(
            "DIFF_RECORD {} rows -> {}",
            rows.len(),
            PathBuf::from(path).display()
        );
        return;
    }
    let base = std::fs::read_to_string(vendored())
        .unwrap_or_else(|e| panic!("the base rows {}: {e}", vendored().display()))
        .lines()
        .map(String::from)
        .collect::<Vec<_>>();
    assert!(base.len() > 100, "too few base rows: {}", base.len());
    let differences = compare(&base, &rows);
    for d in &differences {
        eprintln!("DIFF {d:?}");
    }
    let count = |f: fn(&Difference) -> bool| differences.iter().filter(|d| f(d)).count();
    let cases = base
        .iter()
        .map(|r| r.split('\t').next().unwrap())
        .collect::<std::collections::BTreeSet<_>>()
        .len();
    eprintln!(
        "DIFF_SUMMARY cases={cases} base_rows={} rows={} narrower={} wider={} value={} message={} missing={}",
        base.len(),
        rows.len(),
        count(|d| matches!(d, Difference::Narrower(_))),
        count(|d| matches!(d, Difference::Wider(_))),
        count(|d| matches!(d, Difference::Value { .. })),
        count(|d| matches!(d, Difference::Message(_))),
        count(|d| matches!(d, Difference::Missing(_))),
    );
    assert_eq!(
        rows.len(),
        base.len(),
        "the row count differs from the base's"
    );
    assert!(
        differences.is_empty(),
        "differences from the base: {differences:#?}"
    );
}

/// The comparator itself: a deliberately wrong oracle (a value doubled, an
/// acceptance turned into a refusal) must be reported at once.
#[test]
fn the_comparator_reports_a_wrong_oracle() {
    Vakint::initialize_vakint_symbols();
    let now = vec![
        "c\tvalue\tok\t2*vakint::{}::x".to_string(),
        "d\tvalue\tok\t1".to_string(),
    ];
    let doubled = vec![
        "c\tvalue\tok\t4*vakint::{}::x".to_string(),
        "d\tvalue\tok\t1".to_string(),
    ];
    let refused = vec![
        "c\tvalue\tok\t2*vakint::{}::x".to_string(),
        "d\tvalue\terr\tno".to_string(),
    ];
    let same_value = vec![
        "c\tvalue\tok\tvakint::{}::x+vakint::{}::x".to_string(),
        "d\tvalue\tok\t1".to_string(),
    ];
    assert!(matches!(
        compare(&doubled, &now)[..],
        [Difference::Value {
            algebraic: false,
            ..
        }]
    ));
    assert!(matches!(
        compare(&refused, &now)[..],
        [Difference::Wider(_)]
    ));
    assert!(matches!(
        compare(&now, &refused)[..],
        [Difference::Narrower(_)]
    ));
    // another form of the same value is still a difference (the base's text
    // is pinned), marked algebraic
    assert!(matches!(
        compare(&same_value, &now)[..],
        [Difference::Value {
            algebraic: true,
            ..
        }]
    ));
    let message = vec![
        "c\tvalue\tok\t2*vakint::{}::x".to_string(),
        "d\tvalue\terr\tother".to_string(),
    ];
    assert!(matches!(
        compare(&refused, &message)[..],
        [Difference::Message(_)]
    ));
    let fewer = vec!["c\tvalue\tok\t2*vakint::{}::x".to_string()];
    assert!(matches!(
        compare(&now, &fewer)[..],
        [Difference::Missing(_)]
    ));
}
