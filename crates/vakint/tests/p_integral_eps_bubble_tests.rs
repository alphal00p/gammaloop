//! ε-shifted families against explicit one-loop bubbles (skeleton-mode
//! acceptance (b)): an integral of RustRed's ε-SHIFTED family (a line at power
//! `n + m ε`, `eps_units = m`) equals the same integral written with explicit
//! one-loop bubbles (integer powers, `m` more loops) reduced in the integer
//! family, exactly.
//!
//! The identity (RustRed `docs/p_integrals.md`, conventions C1-C3). With the
//! Euclidean one-loop function
//! `G(a, b) = Γ(a+b-d/2) Γ(d/2-a) Γ(d/2-b) / (Γ(a) Γ(b) Γ(d-a-b))`, a bubble
//! with lines of powers `(a, b)` in series on a line `e` of power `n_e` is the
//! line at power `n_e + a + b - 2 + ε` times `G(a, b)`; for every numerator `N`
//! that does not contain the bubble momentum (the same Minkowski polynomial on
//! both sides)
//!
//! ```text
//! I_M^(L+1)[G + bubble(a, b) on e; N] = G(a, b) · I_M^(L)[G, e at n + ε; N]
//! ```
//!
//! with no sign (C1: the integer parts change by `a + b` and by `a + b - 2`).
//! A bubble inside a bubble multiplies its `G` with the inner ε parts in its
//! arguments (`G(1+ε, 1)` for a `(1, 1)` bubble on a line of the outer one).
//! `PLibrary::reduce` returns the G-scheme `R = (ε G(1,1))^-L I_M` (C3), so
//!
//! ```text
//! R_int = F' · R_eps,   F' = Π_bubbles G(α_i, β_i) / (ε G(1,1))^U
//! ```
//!
//! with `U` the number of bubbles = the total ε units. `F'` is built here from
//! `Γ(n + mε) = Γ(1 + mε) · (rational in d)`, independently of RustRed's
//! kernels.
//!
//! A case is an ε-graph (Forcer notation, one `*` per ε unit), its line powers
//! (dots) and a numerator (products of scalar products of its line momenta and
//! `q`), and one or more ARRANGEMENTS of explicit bubbles on its ε lines
//! (`1+(1,1)`: the line at power 1 and a bubble of two lines at power 1 in
//! series; `1+(1+(1,1),1)`: a bubble inside a bubble; `(2,1)`: a dotted bubble
//! that replaces the line; `-1+(1,1)`: a line at power -1, i.e. the numerator
//! `p_e²`, and a bubble: the ε line at `-1 + ε`). The requester's object
//! `ep(q, m) = (-q²/μ²)^(-mε)` is `m` bubbles `(1,1)` IN SERIES (each `1/ε` in
//! the G-scheme); nested arrangements are the other way r3's families arose.
//! An ε graph may split one line into pieces at two-valent vertices, each with
//! its own units and power (several tagged vertices on one line):
//! `normalise` merges them (units and powers add). Both sides are found through RustRed's canonical
//! graph key (`topology::{normalise, canonical}`, the DAG node's presentation),
//! lowered onto `(PresentationId, IntegralKey, Coefficient)` and reduced with
//! the ordinary entry `PLibrary::reduce`. They are compared
//!
//! 1. exactly on the terminals, when both land on the same ones;
//! 2. exactly in `d` on the masters of the package's catalogue
//!    (`p-master-catalogue.json` rows; `Γ(1 + mε)` tags as symbols);
//! 3. as Laurent series with the catalogue's master values (vakint's
//!    finaliser, `Vakint::refinalise_p_retained`, G-scheme at `Q² = μ² = 1`)
//!    through the deepest order both sides know, at most
//!    `VAKINT_TEST_EPS_BUBBLES_MAX_ORDER` (default 6).
//!
//! Every test is `#[ignore]`d: `just test_vakint` compiles this file but runs
//! none of it. Each needs the one-to-five-loop package (r3) named by
//! `VAKINT_TEST_P_LIBRARY_1_5` (a run without it fails naming the variable,
//! `support/wh_fixtures.rs`), and a vakint build whose RustRed compiles arity
//! 20 (feature `five-loop-p` with `RUSTRED_RUNTIME_ARITIES=1,...,16,20`). The
//! package loads once per process (about 25 s, 7 GiB): run them together,
//!
//! ```text
//! VAKINT_TEST_P_LIBRARY_1_5=<r3>/p_integral/p-library-1-5 \
//! RUSTRED_RUNTIME_ARITIES=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,20 \
//!   cargo test --profile ci-optim -p vakint --features five-loop-p --test p_integral_eps_bubble_tests \
//!   -- --ignored --test-threads 10 --skip four_loop_one_unit_diamond --skip four_loop_one_unit_custom_off_corner
//! ```
//!
//! (about 1 min after the build; the two skipped tests take about 5 and 12 min).
//! Reductions use `Exact` arithmetic, except the two slow tests, which take
//! vakint's engine arithmetic (modular, exact by reconstruction); set
//! `VAKINT_TEST_EPS_BUBBLES_ARITHMETIC=engine` (or `exact`) for every test.
//! `VAKINT_TEST_EPS_BUBBLES_TABLE=<file>` appends one TSV row per comparison.

#[path = "support/p_library4.rs"]
mod p_library4;
#[path = "support/p_rustred.rs"]
mod p_rustred;
#[path = "support/wh_fixtures.rs"]
mod wh_fixtures;

use std::collections::{BTreeMap, HashMap};
use std::io::Write as _;
use std::sync::{Arc, Mutex, OnceLock};
use std::time::Instant;

use rustred::algebra::{Coefficient, CoefficientContext, ExactAlgebraLimits};
use rustred::family::IntegralKey;
use rustred::propagator::PresentationId;
use rustred::propagator::library::{PBasis, PLibrary, PReductionRequest, PTerm, PUncoveredPolicy};
use rustred::propagator::presentation::Slot;
use rustred::propagator::tag::GammaTag;
use rustred::propagator::topology::{Fate, Normalised, PTopo, canonical, normalise};
use rustred::reduction::{OrderedSchedulerOptions, OrderedWeightArithmetic, ReductionLimits};
use serde_json::Value;
use symbolica::atom::{Atom, AtomCore};
use symbolica::domains::rational::Rational;
use vakint::p_integral::PMasterCatalogue;
use vakint::{
    EvaluationMethod, EvaluationOrder, RustRedEvaluationOptions, RustRedMasterMode, Vakint,
    VakintSettings,
};

const PACKAGE: &str = "VAKINT_TEST_P_LIBRARY_1_5";

// ---------------------------------------------------------------------------------------------
// The package, its catalogue and the engine (loaded once)
// ---------------------------------------------------------------------------------------------

/// A catalogue row: `[(master or None for a constant term, coefficient text)]`.
type CatalogueRow = Vec<(Option<String>, String)>;

struct Fixture {
    library: Arc<PLibrary>,
    engine: Vakint,
    /// Canonical graph key bytes -> presentations of that DAG node.
    nodes: HashMap<Vec<u8>, Vec<PresentationId>>,
    /// `(node, key)` -> `[(master or None, coefficient text)]`.
    rows: HashMap<(String, Vec<i64>), CatalogueRow>,
    dimension: Atom,
}

fn fixture() -> &'static Fixture {
    static FIXTURE: OnceLock<Fixture> = OnceLock::new();
    FIXTURE.get_or_init(|| {
        let directory = wh_fixtures::required_env(PACKAGE, "the 1-5 package (r3)");
        let started = Instant::now();
        let library = p_library4::load_p_library_directory(&directory);
        assert_eq!(*library.scope().end(), 5, "the 1-5 package");
        let text = std::fs::read_to_string(directory.join("p-master-catalogue.json")).unwrap();
        let mut engine = p_rustred::rustred_p_engine(Arc::clone(&library));
        engine.set_p_master_catalogue(
            PMasterCatalogue::from_json("p-master-catalogue.json", &text).unwrap(),
        );
        let mut nodes = HashMap::<Vec<u8>, Vec<PresentationId>>::new();
        for presentation in library.presentations() {
            if let Some(node) = &presentation.node {
                nodes
                    .entry(node.graph_key.to_vec())
                    .or_default()
                    .push(presentation.id);
            }
        }
        let catalogue: Value = serde_json::from_str(&text).unwrap();
        let mut rows = HashMap::new();
        for terminal in catalogue["terminals"].as_array().unwrap() {
            let node = terminal["node"].as_str().unwrap().to_owned();
            let key = (terminal["key"].as_array().unwrap().iter())
                .map(|k| k.as_i64().unwrap())
                .collect::<Vec<_>>();
            let combination = (terminal["combination"].as_array().into_iter().flatten())
                .map(|term| {
                    (
                        term["master"].as_str().map(str::to_owned),
                        term["coefficient"].as_str().unwrap().to_owned(),
                    )
                })
                .collect();
            rows.insert((node, key), combination);
        }
        let dimension = library
            .coefficient_context()
            .parameter("d")
            .unwrap()
            .to_expression();
        eprintln!(
            "EPS_BUBBLES package {} loaded in {:.1} s ({} presentations, {} catalogue rows)",
            directory.display(),
            started.elapsed().as_secs_f64(),
            library.presentations().len(),
            rows.len()
        );
        Fixture {
            library,
            engine,
            nodes,
            rows,
            dimension,
        }
    })
}

fn limits() -> ExactAlgebraLimits {
    ExactAlgebraLimits::default()
}

fn context() -> &'static CoefficientContext {
    fixture().library.coefficient_context()
}

// ---------------------------------------------------------------------------------------------
// Graphs, chains of explicit bubbles and their momenta
// ---------------------------------------------------------------------------------------------

#[derive(Clone, Debug)]
struct Line {
    u: u8,
    v: u8,
    eps: u8,
    power: i64,
}

#[derive(Clone, Debug)]
struct Graph {
    n: u8,
    s: u8,
    t: u8,
    lines: Vec<Line>,
}

impl Graph {
    /// `"0-1* 0-2 1-2"` (one `*` per ε unit) with the line powers `powers`.
    fn parse(edges: &str, s: u8, t: u8, powers: &[i64]) -> Self {
        let lines: Vec<Line> = edges
            .split_whitespace()
            .zip(powers.iter().chain(std::iter::repeat(&1)))
            .map(|(edge, &power)| {
                let eps = edge.chars().filter(|&c| c == '*').count() as u8;
                let (u, v) = edge.trim_end_matches('*').split_once('-').unwrap();
                Line {
                    u: u.parse().unwrap(),
                    v: v.parse().unwrap(),
                    eps,
                    power,
                }
            })
            .collect();
        assert_eq!(lines.len(), powers.len(), "{edges}: one power per line");
        let n = lines.iter().map(|l| l.u.max(l.v)).max().unwrap() + 1;
        Self { n, s, t, lines }
    }

    fn loops(&self) -> usize {
        self.lines.len() + 1 - usize::from(self.n)
    }

    /// Momentum rows over `(k_1, .., k_L, q)`, oriented `u -> v`: a
    /// breadth-first tree from `s`, chords carry the loop momenta, `q` enters
    /// at `s` and leaves at `t`.
    fn route(&self) -> Vec<Vec<i64>> {
        let n = usize::from(self.n);
        let width = self.loops() + 1;
        let mut parent: Vec<Option<(usize, usize)>> = vec![None; n]; // (vertex, line)
        let mut seen = vec![false; n];
        let mut order = vec![usize::from(self.s)];
        seen[usize::from(self.s)] = true;
        let mut tree = vec![false; self.lines.len()];
        let mut head = 0;
        while head < order.len() {
            let x = order[head];
            head += 1;
            for (index, line) in self.lines.iter().enumerate() {
                let (u, v) = (usize::from(line.u), usize::from(line.v));
                let other = if u == x {
                    v
                } else if v == x {
                    u
                } else {
                    continue;
                };
                if !seen[other] {
                    seen[other] = true;
                    tree[index] = true;
                    parent[other] = Some((x, index));
                    order.push(other);
                }
            }
        }
        assert!(seen.iter().all(|&s| s), "a connected graph");
        let mut rows = vec![vec![0i64; width]; self.lines.len()];
        // out(w): momentum leaving w through q and the chords
        let mut out = vec![vec![0i64; width]; n];
        out[usize::from(self.t)][width - 1] += 1;
        out[usize::from(self.s)][width - 1] -= 1;
        let mut chord = 0;
        for (index, line) in self.lines.iter().enumerate() {
            if tree[index] {
                continue;
            }
            rows[index][chord] = 1;
            out[usize::from(line.u)][chord] += 1;
            out[usize::from(line.v)][chord] -= 1;
            chord += 1;
        }
        assert_eq!(chord, width - 1, "chords = loops");
        // the tree line into w carries the sum of out() over the subtree of w
        let mut subtree = out.clone();
        for &w in order.iter().rev() {
            if let Some((p, index)) = parent[w] {
                let flow = subtree[w].clone();
                for (sum, x) in subtree[p].iter_mut().zip(&flow) {
                    *sum += x;
                }
                let line = &self.lines[index];
                let sign = if usize::from(line.u) == p { 1 } else { -1 };
                rows[index] = flow.iter().map(|x| sign * x).collect();
            }
        }
        rows
    }
}

/// A chain from one vertex to another: integer lines and bubbles in series.
#[derive(Clone, Debug)]
enum Part {
    Line(i64),
    Bubble(Chain, Chain),
}

#[derive(Clone, Debug)]
struct Chain(Vec<Part>);

/// A power `integer + eps ε`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Power {
    integer: i64,
    eps: i64,
}

impl Chain {
    /// `chain := part ('+' part)*`, `part := [-]INT | '(' chain ',' chain ')'`.
    fn parse(text: &str) -> Self {
        fn chain(s: &[u8], i: &mut usize) -> Chain {
            let mut parts = vec![part(s, i)];
            while *i < s.len() && s[*i] == b'+' {
                *i += 1;
                parts.push(part(s, i));
            }
            Chain(parts)
        }
        fn part(s: &[u8], i: &mut usize) -> Part {
            if s[*i] == b'(' {
                *i += 1;
                let a = chain(s, i);
                assert_eq!(s[*i], b',');
                *i += 1;
                let b = chain(s, i);
                assert_eq!(s[*i], b')');
                *i += 1;
                Part::Bubble(a, b)
            } else {
                // a line of any integer power: at n <= 0 it is a numerator
                // (the key leaves the line's sector), as on the ε side
                let start = *i;
                if s[*i] == b'-' {
                    *i += 1;
                }
                while *i < s.len() && s[*i].is_ascii_digit() {
                    *i += 1;
                }
                let n: i64 = std::str::from_utf8(&s[start..*i]).unwrap().parse().unwrap();
                Part::Line(n)
            }
        }
        let text: String = text.split_whitespace().collect();
        let mut i = 0;
        let parsed = chain(text.as_bytes(), &mut i);
        assert_eq!(i, text.len(), "{text}: trailing input");
        parsed
    }

    /// The Euclidean power of the chain once its bubbles are integrated.
    fn power(&self) -> Power {
        let mut total = Power { integer: 0, eps: 0 };
        for part in &self.0 {
            match part {
                Part::Line(n) => total.integer += n,
                Part::Bubble(a, b) => {
                    let (a, b) = (a.power(), b.power());
                    total.integer += a.integer + b.integer - 2;
                    total.eps += 1 + a.eps + b.eps;
                }
            }
        }
        total
    }

    /// Every bubble's `(α, β)`, inner ones first.
    fn bubbles(&self, out: &mut Vec<(Power, Power)>) {
        for part in &self.0 {
            if let Part::Bubble(a, b) = part {
                a.bubbles(out);
                b.bubbles(out);
                out.push((a.power(), b.power()));
            }
        }
    }

    /// Add the chain from `x` to `y` to `lines` (new vertices from `next`);
    /// returns the lines whose momenta sum to the chain's momentum `x -> y`.
    fn build(&self, x: u8, y: u8, next: &mut u8, lines: &mut Vec<Line>) -> Vec<usize> {
        let mut flow = Vec::new();
        let mut from = x;
        for (index, part) in self.0.iter().enumerate() {
            let to = if index + 1 == self.0.len() {
                y
            } else {
                *next += 1;
                *next - 1
            };
            let entering = match part {
                Part::Line(power) => {
                    lines.push(Line {
                        u: from,
                        v: to,
                        eps: 0,
                        power: *power,
                    });
                    vec![lines.len() - 1]
                }
                Part::Bubble(a, b) => {
                    let mut both = a.build(from, to, next, lines);
                    both.extend(b.build(from, to, next, lines));
                    both
                }
            };
            if index == 0 {
                flow = entering;
            }
            from = to;
        }
        flow
    }
}

/// The integer graph of `graph` with `chains` on its ε lines (in line order),
/// the lines whose momenta make each line of `graph`, and every bubble.
fn explicit(graph: &Graph, chains: &[Chain]) -> (Graph, Vec<Vec<usize>>, Vec<(Power, Power)>) {
    let mut next = graph.n;
    let mut lines = Vec::new();
    let mut flows = Vec::new();
    let mut bubbles = Vec::new();
    let mut chains = chains.iter();
    for line in &graph.lines {
        if line.eps == 0 {
            lines.push(Line { eps: 0, ..*line });
            flows.push(vec![lines.len() - 1]);
            continue;
        }
        let chain = chains.next().expect("one chain per ε line");
        assert_eq!(
            chain.power(),
            Power {
                integer: line.power,
                eps: i64::from(line.eps)
            },
            "the chain {chain:?} makes the line {line:?}"
        );
        chain.bubbles(&mut bubbles);
        flows.push(chain.build(line.u, line.v, &mut next, &mut lines));
    }
    assert!(chains.next().is_none(), "one chain per ε line");
    let int = Graph {
        n: next,
        s: graph.s,
        t: graph.t,
        lines,
    };
    (int, flows, bubbles)
}

// ---------------------------------------------------------------------------------------------
// Lowering a graph onto its DAG node's presentation
// ---------------------------------------------------------------------------------------------

/// `graph` with `numerator` (pairs of momentum rows over the graph's own
/// `(k, q)`) as inputs of the presentation of its DAG node.
struct Lowered {
    presentation: PresentationId,
    inputs: Vec<(IntegralKey, Rational)>,
}

fn rational_rows(rows: &[Vec<i64>]) -> Vec<Vec<Rational>> {
    rows.iter()
        .map(|row| row.iter().map(|&x| Rational::from(x)).collect())
        .collect()
}

fn rank(rows: &[Vec<Rational>]) -> usize {
    let mut m = rows.to_vec();
    let width = m.first().map_or(0, Vec::len);
    let mut rank = 0;
    for col in 0..width {
        let Some(pivot) = (rank..m.len()).find(|&r| !m[r][col].is_zero()) else {
            continue;
        };
        m.swap(rank, pivot);
        let pivot_row = m[rank].clone();
        for (r, row) in m.iter_mut().enumerate() {
            if r != rank && !row[col].is_zero() {
                let f = &row[col] / &pivot_row[col];
                for (x, p) in row.iter_mut().zip(&pivot_row) {
                    let sub = &f * p;
                    *x -= &sub;
                }
            }
        }
        rank += 1;
    }
    rank
}

/// The inverse of the square matrix whose COLUMNS are `columns`.
fn inverse(columns: &[Vec<i64>]) -> Option<Vec<Vec<Rational>>> {
    let n = columns.len();
    let mut a: Vec<Vec<Rational>> = (0..n)
        .map(|r| {
            let mut row: Vec<Rational> = (0..n).map(|c| Rational::from(columns[c][r])).collect();
            row.extend((0..n).map(|k| Rational::from(i64::from(k == r))));
            row
        })
        .collect();
    for col in 0..n {
        let pivot = (col..n).find(|&r| !a[r][col].is_zero())?;
        a.swap(col, pivot);
        let p = a[col][col].clone();
        for x in a[col].iter_mut() {
            *x = &*x / &p;
        }
        let pivot_row = a[col].clone();
        for (r, row) in a.iter_mut().enumerate() {
            if r != col && !row[col].is_zero() {
                let f = row[col].clone();
                for (x, p) in row.iter_mut().zip(&pivot_row) {
                    let sub = &f * p;
                    *x -= &sub;
                }
            }
        }
    }
    Some(a.into_iter().map(|row| row[n..].to_vec()).collect())
}

fn integer(x: &Rational) -> i64 {
    assert!(x.is_integer(), "{x} is not an integer");
    i64::try_from(x.numerator()).unwrap()
}

/// The integer map `T` (`T q = q`) that sends every `(caller, target)` row
/// onto `± target`, `|det T| = 1` (the requester's `eps_frame` search: over a
/// basis of lines, every sign).
fn frame_map(kept: &[(Vec<i64>, Vec<i64>)], width: usize) -> Option<Vec<Vec<i64>>> {
    let mut head = vec![0i64; width];
    head[width - 1] = 1;
    let mut basis = vec![(head.clone(), head)];
    for (caller, target) in kept {
        let mut rows: Vec<Vec<i64>> = basis.iter().map(|(c, _)| c.clone()).collect();
        rows.push(caller.clone());
        if rank(&rational_rows(&rows)) == rows.len() {
            basis.push((caller.clone(), target.clone()));
        }
        if basis.len() == width {
            break;
        }
    }
    if basis.len() != width {
        return None;
    }
    let callers: Vec<Vec<i64>> = basis.iter().map(|(c, _)| c.clone()).collect();
    let inv = inverse(&callers)?;
    'signs: for signs in 0u32..(1 << (width - 1)) {
        let target = |j: usize| -> Vec<i64> {
            let sign = if j == 0 || signs & (1 << (j - 1)) == 0 {
                1
            } else {
                -1
            };
            basis[j].1.iter().map(|x| sign * x).collect()
        };
        let mut t = vec![vec![0i64; width]; width];
        for (r, row) in t.iter_mut().enumerate() {
            for (k, entry) in row.iter_mut().enumerate() {
                let mut sum = Rational::zero();
                for (j, inv_row) in inv.iter().enumerate() {
                    sum += &(Rational::from(target(j)[r]) * &inv_row[k]);
                }
                if !sum.is_integer() {
                    continue 'signs;
                }
                *entry = integer(&sum);
            }
        }
        for (caller, target) in kept {
            let image = apply(&t, caller);
            let minus: Vec<i64> = target.iter().map(|x| -x).collect();
            if &image != target && image != minus {
                continue 'signs;
            }
        }
        let columns: Vec<Vec<i64>> = (0..width)
            .map(|c| (0..width).map(|r| t[r][c]).collect())
            .collect();
        let det = determinant(&columns);
        assert!(det.abs() == 1, "a momentum map of determinant {det}");
        return Some(t);
    }
    None
}

fn apply(t: &[Vec<i64>], v: &[i64]) -> Vec<i64> {
    t.iter()
        .map(|row| row.iter().zip(v).map(|(a, b)| a * b).sum())
        .collect()
}

fn determinant(columns: &[Vec<i64>]) -> i64 {
    let n = columns.len();
    let mut a = rational_rows(columns);
    let mut det = Rational::one();
    for col in 0..n {
        let Some(pivot) = (col..n).find(|&r| !a[r][col].is_zero()) else {
            return 0;
        };
        if pivot != col {
            a.swap(col, pivot);
            det = -det;
        }
        det = &det * &a[col][col];
        let pivot_row = a[col].clone();
        for row in a.iter_mut().skip(col + 1) {
            let f = &row[col] / &pivot_row[col];
            for (x, p) in row.iter_mut().zip(&pivot_row).skip(col) {
                let sub = &f * p;
                *x -= &sub;
            }
        }
    }
    integer(&det)
}

/// Exponents per slot -> coefficient.
type Polynomial = BTreeMap<Vec<u32>, Rational>;

/// `graph` with the numerator `Π_f (a_f · b_f)` (rows over the graph's
/// `(k, q)`, Minkowski scalar products, `q.q = -1` of the frame) on its node's
/// presentation; `Err` names a graph whose node the package does not hold.
fn lower(graph: &Graph, numerator: &[(Vec<i64>, Vec<i64>)]) -> Result<Lowered, String> {
    let fx = fixture();
    let topology = PTopo::new(
        graph.n,
        graph.s,
        graph.t,
        graph.lines.iter().map(|l| (l.u, l.v, l.eps)),
    );
    let Normalised::Topo(normalised, fates) = normalise(&topology) else {
        panic!("{graph:?} normalises to zero or to the trivial topology");
    };
    assert!(
        fates
            .iter()
            .all(|f| matches!(f, Fate::Kept(_) | Fate::SeriesMergedInto(_))),
        "{graph:?}: a contracted or absorbed line ({fates:?})"
    );
    let (key, canon) = canonical(&normalised);
    let presentations = fx
        .nodes
        .get(key.as_bytes())
        .ok_or_else(|| format!("no node {}", key.to_hex()))?;
    let presentation = fx.library.presentation(presentations[0]);
    assert_eq!(usize::from(presentation.loops), graph.loops());
    let width = graph.loops() + 1;
    let rows = graph.route();
    let mut powers = vec![0i64; presentation.slots.len()];
    let mut units = vec![0u8; presentation.slots.len()];
    let mut kept = Vec::new();
    for ((line, fate), row) in graph.lines.iter().zip(&fates).zip(&rows) {
        let slot = usize::from(canon.line_map[usize::from(fate.line().unwrap())]);
        let Slot::Line(target) = &presentation.slots[slot] else {
            panic!(
                "slot {slot} of presentation {} is not a line",
                presentation.id.0
            );
        };
        powers[slot] += line.power;
        units[slot] += line.eps;
        kept.push((row.clone(), target.momentum.clone()));
    }
    for (slot, entry) in presentation.slots.iter().enumerate() {
        if let Slot::Line(line) = entry {
            assert_eq!(line.eps_units, units[slot], "ε units on slot {slot}");
        }
    }
    let t = frame_map(&kept, width).ok_or("no momentum map onto the presentation")?;
    let mut polynomial = Polynomial::from([(vec![0u32; powers.len()], Rational::one())]);
    for (a, b) in numerator {
        let form = presentation
            .express_scalar_product(&apply(&t, a), &apply(&t, b))
            .map_err(|error| error.to_string())?;
        let mut product = Polynomial::new();
        for (exponents, value) in &polynomial {
            if !form.constant.is_zero() {
                *product
                    .entry(exponents.clone())
                    .or_insert_with(Rational::zero) += &(value * &form.constant);
            }
            for (slot, factor) in &form.terms {
                let mut raised = exponents.clone();
                raised[usize::from(*slot)] += 1;
                *product.entry(raised).or_insert_with(Rational::zero) += &(value * factor);
            }
        }
        product.retain(|_, value| !value.is_zero());
        polynomial = product;
    }
    let inputs = polynomial
        .into_iter()
        .map(|(exponents, value)| {
            let key = powers
                .iter()
                .zip(&exponents)
                .map(|(&power, &e)| power - i64::from(e))
                .collect::<Vec<_>>();
            (IntegralKey::try_new(key).unwrap(), value)
        })
        .collect();
    Ok(Lowered {
        presentation: presentation.id,
        inputs,
    })
}

// ---------------------------------------------------------------------------------------------
// Γ factors, reductions and the three comparisons
// ---------------------------------------------------------------------------------------------

/// `rational(d) · Π Γ(1 + mε)^e`.
#[derive(Clone, Debug)]
struct GammaFactor {
    rational: Coefficient,
    tag: GammaTag,
}

impl GammaFactor {
    fn one() -> Self {
        Self {
            rational: context().one(),
            tag: GammaTag::one(),
        }
    }

    /// `k + m ε` at `ε = (4 - d)/2`.
    fn linear(k: i64, m: i64) -> Coefficient {
        let c = context();
        let d = c.parameter("d").unwrap();
        let twice = c
            .try_sub(
                &c.integer(2 * k + 4 * m),
                &c.try_mul(&c.integer(m), &d, limits()).unwrap(),
                limits(),
            )
            .unwrap();
        c.try_div(&twice, &c.integer(2), limits()).unwrap()
    }

    /// `Γ(n + mε)^e`.
    fn gamma(n: i64, m: i64, e: i64) -> Self {
        let c = context();
        let mut value = c.one();
        if m == 0 {
            assert!(n >= 1, "Γ({n}) is a pole");
            for k in 1..n {
                value = c.try_mul(&value, &c.integer(k), limits()).unwrap();
            }
        } else if n >= 1 {
            for k in 1..n {
                value = c.try_mul(&value, &Self::linear(k, m), limits()).unwrap();
            }
        } else {
            for k in n..=0 {
                value = c.try_div(&value, &Self::linear(k, m), limits()).unwrap();
            }
        }
        let mut rational = c.one();
        for _ in 0..e.abs() {
            rational = if e > 0 {
                c.try_mul(&rational, &value, limits()).unwrap()
            } else {
                c.try_div(&rational, &value, limits()).unwrap()
            };
        }
        let tag = if m == 0 {
            GammaTag::one()
        } else {
            GammaTag::from_pairs([(m as i8, e as i16)]).unwrap()
        };
        Self { rational, tag }
    }

    fn times(&self, other: &Self) -> Self {
        Self {
            rational: context()
                .try_mul(&self.rational, &other.rational, limits())
                .unwrap(),
            tag: self.tag.add(&other.tag),
        }
    }

    /// The Euclidean `G(α, β)` at `d = 4 - 2ε`.
    fn g(alpha: Power, beta: Power) -> Self {
        let (a, ja, b, jb) = (alpha.integer, alpha.eps, beta.integer, beta.eps);
        Self::gamma(a + b - 2, 1 + ja + jb, 1)
            .times(&Self::gamma(2 - a, -1 - ja, 1))
            .times(&Self::gamma(2 - b, -1 - jb, 1))
            .times(&Self::gamma(a, ja, -1))
            .times(&Self::gamma(b, jb, -1))
            .times(&Self::gamma(4 - a - b, -2 - ja - jb, -1))
    }

    /// `1 / (ε G(1,1))`, checked against C3: `(d - 3) · tag{(1,-1), (-1,-2), (-2,1)}`.
    fn per_loop() -> Self {
        let c = context();
        let one = Power { integer: 1, eps: 0 };
        let g = Self::g(one, one);
        let eps = Self::linear(0, 1);
        let value = Self {
            rational: c
                .try_div(
                    &c.one(),
                    &c.try_mul(&eps, &g.rational, limits()).unwrap(),
                    limits(),
                )
                .unwrap(),
            tag: g.tag.neg(),
        };
        let d = c.parameter("d").unwrap();
        let expected = c.try_sub(&d, &c.integer(3), limits()).unwrap();
        assert!(
            c.try_sub(&value.rational, &expected, limits())
                .unwrap()
                .is_zero()
                && value.tag == GammaTag::from_pairs([(1, -1), (-1, -2), (-2, 1)]).unwrap(),
            "C3: 1/(ε G(1,1)) = (d-3) tag{{(1,-1),(-1,-2),(-2,1)}}"
        );
        value
    }

    /// `F' = Π G(α_i, β_i) / (ε G(1,1))^U`.
    fn of(bubbles: &[(Power, Power)]) -> Self {
        let per_loop = Self::per_loop();
        bubbles.iter().fold(Self::one(), |acc, &(alpha, beta)| {
            acc.times(&Self::g(alpha, beta)).times(&per_loop)
        })
    }

    /// The factor as a vakint expression at `d = 4 - 2ε`.
    fn atom(&self) -> Atom {
        tag_atom(&self.tag) * at_dimension(&self.rational)
    }
}

/// The arithmetic of the reductions: `Exact` (the default), or the one of
/// vakint's ordinary p requests (`modular_polynomial` at one to four loops,
/// `modular_sampled` at five; exact results by reconstruction) for a test
/// that sets it, or for every test with
/// `VAKINT_TEST_EPS_BUBBLES_ARITHMETIC=engine`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Arithmetic {
    Exact,
    Engine,
}

thread_local! {
    static ARITHMETIC: std::cell::Cell<Arithmetic> = const { std::cell::Cell::new(Arithmetic::Exact) };
}

fn arithmetic_choice() -> Arithmetic {
    match std::env::var("VAKINT_TEST_EPS_BUBBLES_ARITHMETIC").as_deref() {
        Ok("engine") => Arithmetic::Engine,
        Ok("exact") => Arithmetic::Exact,
        Err(_) => ARITHMETIC.with(|a| a.get()),
        Ok(other) => panic!("VAKINT_TEST_EPS_BUBBLES_ARITHMETIC={other}: exact or engine"),
    }
}

fn arithmetic(loops: usize) -> OrderedWeightArithmetic {
    match arithmetic_choice() {
        Arithmetic::Exact => OrderedWeightArithmetic::Exact,
        Arithmetic::Engine if loops >= 5 => {
            settings(3).rustred_five_loop_ordered_weight_arithmetic()
        }
        Arithmetic::Engine => settings(3).rustred_ordered_weight_arithmetic(),
    }
}

/// The reduction of `lowered`, nonzero terms only.
fn reduce(lowered: &Lowered) -> Vec<PTerm> {
    let c = context();
    let inputs = (lowered.inputs.iter())
        .map(|(key, value)| {
            let numerator = c.integer(i64::try_from(value.numerator()).unwrap());
            let denominator = c.integer(i64::try_from(value.denominator()).unwrap());
            (
                lowered.presentation,
                key.clone(),
                c.try_div(&numerator, &denominator, limits()).unwrap(),
            )
        })
        .collect();
    let mut reduction = fixture()
        .library
        .reduce(PReductionRequest {
            inputs,
            arithmetic: arithmetic(usize::from(
                fixture().library.presentation(lowered.presentation).loops,
            )),
            scheduler: OrderedSchedulerOptions::default(),
            step_cache: None,
            uncovered: PUncoveredPolicy::Reject,
            onto_masters: true,
            limits: ReductionLimits::default(),
            inspector: None,
        })
        .unwrap_or_else(|error| {
            panic!(
                "reduction of presentation {}: {error}",
                lowered.presentation.0
            )
        });
    reduction.terms.retain(|term| !term.coefficient.is_zero());
    reduction.terms
}

fn basis_name(basis: &PBasis) -> String {
    match basis {
        PBasis::Unit => "Unit".to_owned(),
        PBasis::Master(name) => format!("PMaster({name})"),
        PBasis::Residual(id) => {
            let residual = fixture().library.residual_terminal(*id).unwrap().unwrap();
            format!("{}{:?}", residual.node, residual.key.powers())
        }
    }
}

/// Terms by `(terminal, tag)`, each scaled by `factor`.
fn by_terminal(terms: &[PTerm], factor: &GammaFactor) -> BTreeMap<(String, GammaTag), Coefficient> {
    let c = context();
    let mut out = BTreeMap::<(String, GammaTag), Coefficient>::new();
    for term in terms {
        let coefficient = c
            .try_mul(&term.coefficient, &factor.rational, limits())
            .unwrap();
        let key = (basis_name(&term.basis), term.tag.add(&factor.tag));
        let sum = match out.remove(&key) {
            Some(previous) => c.try_add(&previous, &coefficient, limits()).unwrap(),
            None => coefficient,
        };
        out.insert(key, sum);
    }
    out.retain(|_, value| !value.is_zero());
    out
}

fn epsilon() -> Atom {
    Atom::var(symbolica::symbol!(
        &VakintSettings::default().epsilon_symbol
    ))
}

fn at_dimension(coefficient: &Coefficient) -> Atom {
    let dimension = Atom::num(4) - Atom::num(2) * epsilon();
    coefficient
        .to_expression()
        .replace(fixture().dimension.to_pattern())
        .with(dimension.to_pattern())
}

fn tag_atom(tag: &GammaTag) -> Atom {
    let gamma = symbolica::symbol!("vakint::PGamma");
    let mut atom = Atom::num(1);
    for &(m, e) in tag.exponents() {
        atom *= gamma.call((i64::from(m),)).pow(Atom::num(i64::from(e)));
    }
    atom
}

/// vakint's `Retain` form of a reduction (`Π PGamma(m)^e · basis ·
/// coefficient(d = 4 - 2ε)`, `PResidual(node, k..)` for a residual), and the
/// same with every residual replaced by its catalogue row (masters as symbols).
fn render(terms: &[PTerm]) -> (Atom, Atom) {
    let fx = fixture();
    let residual = symbolica::symbol!("vakint::PResidual");
    let mut retained = Atom::Zero;
    let mut on_masters = Atom::Zero;
    for term in terms {
        let weight = tag_atom(&term.tag) * at_dimension(&term.coefficient);
        let (basis, masters) = match &term.basis {
            PBasis::Unit => (Atom::num(1), Atom::num(1)),
            PBasis::Master(name) => panic!("a Forcer master {name} from a 1-5 package"),
            PBasis::Residual(id) => {
                let terminal = fx.library.residual_terminal(*id).unwrap().unwrap();
                let mut arguments = vec![Atom::var(symbolica::symbol!(&format!(
                    "vakint::{}",
                    terminal.node
                )))];
                arguments.extend(terminal.key.powers().iter().map(|&k| Atom::num(k)));
                let row = fx
                    .rows
                    .get(&(terminal.node.clone(), terminal.key.powers().to_vec()))
                    .unwrap_or_else(|| {
                        panic!(
                            "no catalogue row for {}{:?}",
                            terminal.node,
                            terminal.key.powers()
                        )
                    });
                let mut masters = Atom::Zero;
                for (master, coefficient) in row {
                    let master = match master {
                        Some(name) => Atom::var(symbolica::symbol!(&format!(
                            "vakint::EQ_{}",
                            name.replace(|c: char| !c.is_ascii_alphanumeric(), "_")
                        ))),
                        None => Atom::num(1),
                    };
                    masters += master * catalogue_coefficient(coefficient);
                }
                (residual.call(&arguments[..]), masters)
            }
        };
        retained += &weight * basis;
        on_masters += weight * masters;
    }
    (retained, on_masters)
}

/// `PGamma(m)` -> the symbol `EQG_m` (`|m| <= 8`): the Gamma tags as
/// independent variables of a rational function.
fn gamma_as_symbols(expression: &Atom) -> Atom {
    let gamma = symbolica::symbol!("vakint::PGamma");
    let mut out = expression.clone();
    for m in -8i64..=8 {
        if m == 0 {
            continue;
        }
        let name = format!("vakint::EQG_{}{}", if m < 0 { "n" } else { "p" }, m.abs());
        out = out
            .replace(gamma.call((m,)).to_pattern())
            .with(Atom::var(symbolica::symbol!(&name)).to_pattern());
    }
    assert!(
        !out.contains_symbol(gamma),
        "a Gamma tag beyond |m| = 8 in {out}"
    );
    out
}

/// A catalogue coefficient (`d`, `PGamma(m)`) at `d = 4 - 2ε`.
fn catalogue_coefficient(text: &str) -> Atom {
    let parsed = symbolica::try_parse!(text, default_namespace = "vakint")
        .unwrap_or_else(|error| panic!("{text}: {error}"));
    let d = symbolica::try_parse!("d", default_namespace = "vakint").unwrap();
    let dimension = Atom::num(4) - Atom::num(2) * epsilon();
    parsed.replace(d.to_pattern()).with(dimension.to_pattern())
}

fn settings(terms: i64) -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Symbolic,
            },
        )]),
        number_of_terms_in_epsilon_expansion: terms,
        use_dot_product_notation: true,
        n_cores: 1,
        ..VakintSettings::default()
    }
}

/// The Laurent series of `retained` (`loops` loops) through `ε^order`.
fn laurent(retained: &Atom, loops: usize, order: i64) -> Result<Atom, String> {
    fixture()
        .engine
        .refinalise_p_retained(
            &settings(order + loops as i64 + 1),
            retained.as_view(),
            loops,
        )
        .map_err(|error| error.to_string())
}

// ---------------------------------------------------------------------------------------------
// Cases
// ---------------------------------------------------------------------------------------------

/// One ε-family integral and its explicit-bubble arrangements.
struct Case {
    graph: &'static str,
    s: u8,
    t: u8,
    powers: &'static [i64],
    /// Factors `a.b` with `a, b` = `lN` (momentum of line N) or `q`.
    numerator: &'static str,
    /// Per arrangement, one chain per ε line (in line order).
    arrangements: &'static [&'static [&'static str]],
}

/// One comparison, as a row of the table.
struct Row {
    class: (usize, i64),
    case: String,
    arrangement: String,
    eps_presentation: u32,
    int_presentation: String,
    terminals: String,
    masters: String,
    laurent: String,
    arithmetic: Arithmetic,
}

fn table(row: &Row) {
    let line = format!(
        "({},{})\t{}\t{}\tp{}\t{}\t{}\t{}\t{}\t{:?}",
        row.class.0,
        row.class.1,
        row.case,
        row.arrangement,
        row.eps_presentation,
        row.int_presentation,
        row.terminals,
        row.masters,
        row.laurent,
        row.arithmetic
    );
    eprintln!("EPS_BUBBLES {line}");
    if let Some(path) = std::env::var_os("VAKINT_TEST_EPS_BUBBLES_TABLE") {
        static LOCK: Mutex<()> = Mutex::new(());
        let _guard = LOCK.lock().unwrap();
        let mut file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
            .unwrap();
        writeln!(file, "{line}").unwrap();
    }
}

fn max_order() -> i64 {
    std::env::var("VAKINT_TEST_EPS_BUBBLES_MAX_ORDER").map_or(6, |v| v.parse().unwrap())
}

fn momentum(
    token: &str,
    rows: &[Vec<i64>],
    flows: Option<&[Vec<usize>]>,
    width: usize,
) -> Vec<i64> {
    if token == "q" {
        let mut q = vec![0i64; width];
        q[width - 1] = 1;
        return q;
    }
    let line: usize = token.strip_prefix('l').unwrap().parse().unwrap();
    match flows {
        None => rows[line].clone(),
        Some(flows) => {
            let mut sum = vec![0i64; width];
            for &index in &flows[line] {
                for (s, x) in sum.iter_mut().zip(&rows[index]) {
                    *s += x;
                }
            }
            sum
        }
    }
}

fn numerator_rows(
    text: &str,
    graph: &Graph,
    flows: Option<&[Vec<usize>]>,
) -> Vec<(Vec<i64>, Vec<i64>)> {
    let rows = graph.route();
    let width = graph.loops() + 1;
    text.split_whitespace()
        .map(|factor| {
            let (a, b) = factor.split_once('.').unwrap();
            (
                momentum(a, &rows, flows, width),
                momentum(b, &rows, flows, width),
            )
        })
        .collect()
}

/// Every arrangement of `case` against its ε-family integral; panics on the
/// first inequality, returns the number of comparisons made.
fn check(case: &Case) -> usize {
    let eps_graph = Graph::parse(case.graph, case.s, case.t, case.powers);
    let loops = eps_graph.loops();
    let units: i64 = eps_graph.lines.iter().map(|l| i64::from(l.eps)).sum();
    let name = format!(
        "{} [{}]{}",
        case.graph,
        case.powers
            .iter()
            .map(i64::to_string)
            .collect::<Vec<_>>()
            .join(","),
        if case.numerator.is_empty() {
            String::new()
        } else {
            format!(" N={}", case.numerator)
        }
    );
    let started = Instant::now();
    let eps_lowered = lower(
        &eps_graph,
        &numerator_rows(case.numerator, &eps_graph, None),
    )
    .unwrap_or_else(|error| panic!("{name}: the ε family: {error}"));
    // A zero (e.g. a numerator odd under an automorphism) must be zero on
    // the explicit side too.
    let eps_terms = reduce(&eps_lowered);
    eprintln!(
        "EPS_BUBBLES_STEP {name}: ε side p{} ({} inputs) reduced in {:.1} s",
        eps_lowered.presentation.0,
        eps_lowered.inputs.len(),
        started.elapsed().as_secs_f64()
    );
    let (eps_retained, eps_masters) = render(&eps_terms);
    let mut compared = 0;
    let mut missing = Vec::new();
    for arrangement in case.arrangements {
        let chains: Vec<Chain> = arrangement.iter().map(|text| Chain::parse(text)).collect();
        let label = arrangement.join(" ; ");
        let (int_graph, flows, bubbles) = explicit(&eps_graph, &chains);
        assert_eq!(
            bubbles.len() as i64,
            units,
            "{name} {label}: one bubble per ε unit"
        );
        let numerator = numerator_rows(case.numerator, &int_graph, Some(&flows));
        let int_lowered = match lower(&int_graph, &numerator) {
            Ok(lowered) => lowered,
            Err(error) => {
                missing.push(format!("{label}: {error}"));
                continue;
            }
        };
        let step = Instant::now();
        let int_terms = reduce(&int_lowered);
        eprintln!(
            "EPS_BUBBLES_STEP {name} [{label}]: explicit side p{} ({} inputs) reduced in {:.1} s",
            int_lowered.presentation.0,
            int_lowered.inputs.len(),
            step.elapsed().as_secs_f64()
        );
        let factor = GammaFactor::of(&bubbles);
        // 1. the terminals
        let int_side = by_terminal(&int_terms, &GammaFactor::one());
        let eps_side = by_terminal(&eps_terms, &factor);
        let same = int_side.keys().eq(eps_side.keys());
        let terminals = if int_side.is_empty() && eps_side.is_empty() {
            "both zero".to_owned()
        } else if same {
            let c = context();
            let equal = int_side
                .iter()
                .zip(&eps_side)
                .all(|((_, a), (_, b))| c.try_sub(a, b, limits()).unwrap().is_zero());
            assert!(
                equal,
                "{name} {label}: same terminals, different coefficients"
            );
            format!("equal exactly on {} terminal terms", int_side.len())
        } else {
            format!(
                "different terminals ({} vs {})",
                int_side.len(),
                eps_side.len()
            )
        };
        // 2. the catalogue's masters, exactly in d
        let (int_retained, int_masters) = render(&int_terms);
        let difference =
            gamma_as_symbols(&(&int_masters - &factor.atom() * &eps_masters)).together();
        assert!(
            difference.is_zero(),
            "{name} {label}: the catalogue masters differ by {difference}"
        );
        // 3. Laurent series through the deepest order both sides know
        let step = Instant::now();
        let eps_scaled = (&factor.atom() * &eps_retained).expand();
        let mut reached = None;
        let mut limit = String::new();
        let orders = if int_retained.is_zero() && eps_retained.is_zero() {
            -(loops as i64 + units)..=-(loops as i64 + units) - 1
        } else {
            -(loops as i64 + units)..=max_order()
        };
        for order in orders.rev() {
            let int_series = laurent(&int_retained, loops + units as usize, order);
            let eps_series = laurent(&eps_scaled, loops, order);
            match (int_series, eps_series) {
                (Ok(a), Ok(b)) => {
                    let difference = (&a - &b).expand();
                    assert!(
                        difference.is_zero(),
                        "{name} {label}: the Laurent series through ε^{order} differ by {difference}"
                    );
                    reached = Some(order);
                    break;
                }
                (a, b) => {
                    if limit.is_empty() {
                        limit = [a.err(), b.err()]
                            .into_iter()
                            .flatten()
                            .map(|e| e.chars().take(160).collect::<String>())
                            .collect::<Vec<_>>()
                            .join(" | ");
                    }
                }
            }
        }
        let laurent = match reached {
            None if int_retained.is_zero() && eps_retained.is_zero() => {
                "both zero (all orders)".to_owned()
            }
            Some(order) if order == max_order() => format!("equal through eps^{order} (cap)"),
            Some(order) => format!("equal through eps^{order} (deeper: {limit})"),
            None => format!("no order known: {limit}"),
        };
        eprintln!(
            "EPS_BUBBLES_STEP {name} [{label}]: Laurent in {:.1} s",
            step.elapsed().as_secs_f64()
        );
        table(&Row {
            class: (loops, units),
            case: name.clone(),
            arrangement: label,
            eps_presentation: eps_lowered.presentation.0,
            int_presentation: format!(
                "p{} ({} inputs)",
                int_lowered.presentation.0,
                int_lowered.inputs.len()
            ),
            terminals,
            masters: "equal exactly in d".to_owned(),
            laurent,
            arithmetic: arithmetic_choice(),
        });
        compared += 1;
    }
    for entry in &missing {
        eprintln!("EPS_BUBBLES_MISSING ({loops},{units}) {name}: {entry}");
    }
    assert!(
        compared > 0,
        "{name}: no arrangement is a node of the package: {missing:?}"
    );
    eprintln!(
        "EPS_BUBBLES_CASE ({loops},{units}) {name}: {compared} arrangement(s) equal, {:.1} s",
        started.elapsed().as_secs_f64()
    );
    compared
}

fn check_all(class: (usize, i64), cases: &[Case]) {
    check_all_with(class, cases, Arithmetic::Exact);
}

fn check_all_with(class: (usize, i64), cases: &[Case], arithmetic: Arithmetic) {
    ARITHMETIC.with(|a| a.set(arithmetic));
    let mut total = 0;
    for case in cases {
        let graph = Graph::parse(case.graph, case.s, case.t, case.powers);
        let units: i64 = graph.lines.iter().map(|l| i64::from(l.eps)).sum();
        assert_eq!((graph.loops(), units), class, "{}: class", case.graph);
        total += check(case);
    }
    eprintln!(
        "EPS_BUBBLES_CLASS ({},{}): {total} comparisons equal",
        class.0, class.1
    );
}

// ---------------------------------------------------------------------------------------------
// The classes of the r3 package (loops, total ε units)
// ---------------------------------------------------------------------------------------------

/// (1,1): the one-loop bubble with one ε line (the only node of the class);
/// the ε line also at n = 0, -1, -2 (pinched: the key leaves its sector).
#[test]
#[ignore = "external data: VAKINT_TEST_P_LIBRARY_1_5 (the 1-5 package)"]
#[rustfmt::skip]
fn one_loop_one_unit() {
    check_all(
        (1, 1),
        &[
            Case { graph: "0-1 0-1*", s: 0, t: 1, powers: &[1, 1], numerator: "",
                arrangements: &[&["1+(1,1)"], &["(2,1)"], &["(1,2)"]] },
            Case { graph: "0-1 0-1*", s: 0, t: 1, powers: &[2, 2], numerator: "l0.q",
                arrangements: &[&["2+(1,1)"], &["1+(2,1)"], &["(2,2)"]] },
            Case { graph: "0-1 0-1*", s: 0, t: 1, powers: &[1, 1], numerator: "l1.q l1.q",
                arrangements: &[&["1+(1,1)"]] },
            Case { graph: "0-1 0-1*", s: 0, t: 1, powers: &[1, 0], numerator: "",
                arrangements: &[&["(1,1)"], &["-1+(2,1)"]] },
            Case { graph: "0-1 0-1*", s: 0, t: 1, powers: &[1, -1], numerator: "",
                arrangements: &[&["-1+(1,1)"], &["(1,1)+-1"]] },
            Case { graph: "0-1 0-1*", s: 0, t: 1, powers: &[2, -2], numerator: "l0.q",
                arrangements: &[&["-2+(1,1)"]] },
        ],
    );
}

/// (1,2): two units on one line (two bubbles in series, the requester's
/// chain; a bubble inside a bubble), one unit on each line, two tagged pieces
/// merged in series, and the ε line at n = 0, -1.
#[test]
#[ignore = "external data: VAKINT_TEST_P_LIBRARY_1_5 (the 1-5 package)"]
#[rustfmt::skip]
fn one_loop_two_units() {
    check_all(
        (1, 2),
        &[
            Case { graph: "0-1 0-1**", s: 0, t: 1, powers: &[1, 1], numerator: "",
                arrangements: &[&["1+(1,1)+(1,1)"], &["1+(1+(1,1),1)"], &["(1+(1,1),2)"]] },
            Case { graph: "0-1 0-1**", s: 0, t: 1, powers: &[1, 2], numerator: "l0.q",
                arrangements: &[&["1+(1,1)+(2,1)"], &["1+(2+(1,1),1)"]] },
            Case { graph: "0-1* 0-1*", s: 0, t: 1, powers: &[1, 1], numerator: "l0.l1",
                arrangements: &[&["1+(1,1)", "1+(1,1)"], &["(2,1)", "1+(1,1)"]] },
            Case { graph: "0-2* 2-1* 0-1", s: 0, t: 1, powers: &[1, 0, 1], numerator: "",
                arrangements: &[&["1+(1,1)", "(1,1)"]] },
            Case { graph: "0-1 0-1**", s: 0, t: 1, powers: &[1, 0], numerator: "",
                arrangements: &[&["(1,1)+(1,1)"]] },
            Case { graph: "0-1 0-1**", s: 0, t: 1, powers: &[2, -1], numerator: "l1.q",
                arrangements: &[&["-1+(1,1)+(1,1)"], &["(1,1)+-1+(1,1)"]] },
        ],
    );
}

/// (1,3): three units on one line (chain of three bubbles; nested forms),
/// units on two lines, three tagged pieces merged in series, n = -1.
#[test]
#[ignore = "external data: VAKINT_TEST_P_LIBRARY_1_5 (the 1-5 package)"]
#[rustfmt::skip]
fn one_loop_three_units() {
    check_all(
        (1, 3),
        &[
            Case { graph: "0-1 0-1***", s: 0, t: 1, powers: &[1, 1], numerator: "",
                arrangements: &[&["1+(1,1)+(1,1)+(1,1)"], &["1+(1+(1+(1,1),1),1)"], &["1+(1+(1,1),1+(1,1))"]] },
            Case { graph: "0-1 0-1***", s: 0, t: 1, powers: &[1, 2], numerator: "l0.q l1.q",
                arrangements: &[&["2+(1,1)+(1,1)+(1,1)"]] },
            Case { graph: "0-1* 0-1**", s: 0, t: 1, powers: &[1, 1], numerator: "l1.q",
                arrangements: &[&["1+(1,1)", "1+(1+(1,1),1)"], &["1+(1,1)", "1+(1,1)+(1,1)"]] },
            Case { graph: "0-2* 2-3* 3-1* 0-1", s: 0, t: 1, powers: &[1, 0, 0, 1], numerator: "",
                arrangements: &[&["1+(1,1)", "(1,1)", "(1,1)"]] },
            Case { graph: "0-1 0-1***", s: 0, t: 1, powers: &[1, -1], numerator: "",
                arrangements: &[&["-1+(1,1)+(1,1)+(1,1)"]] },
        ],
    );
}

/// (1,4): five-loop explicit graphs; four units on one line as a chain, two
/// tagged pieces of two units merged, n = 0 and -2.
#[test]
#[ignore = "external data: VAKINT_TEST_P_LIBRARY_1_5 (the 1-5 package)"]
#[rustfmt::skip]
fn one_loop_four_units() {
    check_all(
        (1, 4),
        &[
            Case { graph: "0-1 0-1****", s: 0, t: 1, powers: &[1, 1], numerator: "",
                arrangements: &[&["1+(1,1)+(1,1)+(1,1)+(1,1)"], &["1+(1+(1+(1+(1,1),1),1),1)"]] },
            Case { graph: "0-1 0-1****", s: 0, t: 1, powers: &[1, 2], numerator: "l0.q",
                arrangements: &[&["2+(1,1)+(1,1)+(1,1)+(1,1)"]] },
            Case { graph: "0-1* 0-1***", s: 0, t: 1, powers: &[1, 1], numerator: "l0.q",
                arrangements: &[&["1+(1,1)", "1+(1+(1,1),1)+(1,1)"], &["1+(1,1)", "1+(1,1)+(1,1)+(1,1)"]] },
            Case { graph: "0-1** 0-1**", s: 0, t: 1, powers: &[1, 1], numerator: "",
                arrangements: &[&["1+(1,1)+(1,1)", "1+(1+(1,1),1)"], &["1+(1,1)+(1,1)", "1+(1,1)+(1,1)"]] },
            Case { graph: "0-2** 2-1** 0-1", s: 0, t: 1, powers: &[1, 0, 1], numerator: "",
                arrangements: &[&["1+(1,1)+(1,1)", "(1,1)+(1,1)"]] },
            Case { graph: "0-1 0-1****", s: 0, t: 1, powers: &[1, 0], numerator: "",
                arrangements: &[&["(1,1)+(1,1)+(1,1)+(1,1)"]] },
            Case { graph: "0-1 0-1****", s: 0, t: 1, powers: &[2, -2], numerator: "l0.q",
                arrangements: &[&["-2+(1,1)+(1,1)+(1,1)+(1,1)"]] },
        ],
    );
}

/// (2,1): t1star05 (Custom node, master t105) and the Triangle node; the ε
/// line at n = 0, -1, -2; a tagged piece merged with a pinched neighbour.
#[test]
#[ignore = "external data: VAKINT_TEST_P_LIBRARY_1_5 (the 1-5 package)"]
#[rustfmt::skip]
fn two_loop_one_unit() {
    check_all(
        (2, 1),
        &[
            Case { graph: "0-1* 0-2 0-3 1-2 1-3", s: 2, t: 3, powers: &[1, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["1+(1,1)"], &["(2,1)"]] },
            Case { graph: "0-1* 0-2 0-3 1-2 1-3", s: 2, t: 3, powers: &[2, 1, 1, 1, 2], numerator: "l1.q",
                arrangements: &[&["1+(2,1)"], &["2+(1,1)"]] },
            Case { graph: "0-1* 0-2 0-3 1-2 1-3", s: 2, t: 3, powers: &[1, 1, 1, 1, 1], numerator: "l1.l4 l2.q",
                arrangements: &[&["1+(1,1)"]] },
            Case { graph: "0-1* 0-2 0-3 1-2 1-3", s: 2, t: 3, powers: &[0, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["(1,1)"]] },
            Case { graph: "0-1* 0-2 0-3 1-2 1-3", s: 2, t: 3, powers: &[-1, 1, 1, 1, 1], numerator: "l1.q",
                arrangements: &[&["-1+(1,1)"], &["(1,1)+-1"]] },
            Case { graph: "0-1* 0-2 0-3 1-2 1-3", s: 2, t: 3, powers: &[-2, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["-2+(1,1)"]] },
            Case { graph: "0-4* 4-1 0-2 0-3 1-2 1-3", s: 2, t: 3, powers: &[0, 1, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["(1,1)"]] },
            Case { graph: "0-1 0-2 0-3 1-2 1-3*", s: 2, t: 3, powers: &[1, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["1+(1,1)"], &["(1,2)"]] },
            Case { graph: "0-1 0-2 0-3 1-2 1-3*", s: 2, t: 3, powers: &[1, 2, 1, 1, 1], numerator: "l0.q",
                arrangements: &[&["1+(1,1)"]] },
            Case { graph: "0-1 0-2 0-3 1-2 1-3*", s: 2, t: 3, powers: &[1, 1, 1, 1, -1], numerator: "l1.q",
                arrangements: &[&["-1+(1,1)"]] },
        ],
    );
}

/// (2,2): t1star55 (two units on one line: chain and nested), t1star24 /
/// t1star45 (two lines), two tagged pieces merged, n = 0, -1.
#[test]
#[ignore = "external data: VAKINT_TEST_P_LIBRARY_1_5 (the 1-5 package)"]
#[rustfmt::skip]
fn two_loop_two_units() {
    check_all(
        (2, 2),
        &[
            Case { graph: "0-1** 0-2 0-3 1-2 1-3", s: 2, t: 3, powers: &[1, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["1+(1,1)+(1,1)"], &["1+(1+(1,1),1)"]] },
            Case { graph: "0-1** 0-2 0-3 1-2 1-3", s: 2, t: 3, powers: &[1, 1, 2, 1, 1], numerator: "l1.q",
                arrangements: &[&["1+(1,1)+(1,1)"], &["1+(1+(1,1),1)"]] },
            Case { graph: "0-1** 0-2 0-3 1-2 1-3", s: 2, t: 3, powers: &[0, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["(1,1)+(1,1)"]] },
            Case { graph: "0-1** 0-2 0-3 1-2 1-3", s: 2, t: 3, powers: &[-1, 1, 1, 1, 1], numerator: "l1.q",
                arrangements: &[&["-1+(1,1)+(1,1)"]] },
            Case { graph: "0-4* 4-1* 0-2 0-3 1-2 1-3", s: 2, t: 3, powers: &[1, 0, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["1+(1,1)", "(1,1)"]] },
            Case { graph: "0-1* 0-2 0-3 1-2 1-3*", s: 2, t: 3, powers: &[1, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["1+(1,1)", "1+(1,1)"], &["(2,1)", "1+(1,1)"]] },
            Case { graph: "0-1* 0-2 0-3 1-2 1-3*", s: 2, t: 3, powers: &[1, 1, 1, 2, 1], numerator: "l2.q l3.q",
                arrangements: &[&["1+(1,1)", "1+(1,1)"]] },
            Case { graph: "0-1* 0-2 0-3 1-2 1-3*", s: 2, t: 3, powers: &[-1, 1, 1, 1, 0], numerator: "",
                arrangements: &[&["-1+(1,1)", "(1,1)"]] },
            // l0.q is odd under the swap of vertices 0 and 1: zero on both sides
            Case { graph: "0-1 0-2 0-3* 1-2 1-3*", s: 2, t: 3, powers: &[1, 1, 1, 1, 1], numerator: "l0.q",
                arrangements: &[&["1+(1,1)", "1+(1,1)"]] },
            Case { graph: "0-1 0-2 0-3* 1-2 1-3*", s: 2, t: 3, powers: &[1, 1, 1, 1, 1], numerator: "l1.q",
                arrangements: &[&["1+(1,1)", "1+(1,1)"], &["(2,1)", "1+(1,1)"]] },
        ],
    );
}

/// (2,3): three units on one line (chain), on two lines, and one unit on each
/// of three lines (Custom callers, five-loop explicit graphs); a merged piece;
/// n = 0, -1.
#[test]
#[ignore = "external data: VAKINT_TEST_P_LIBRARY_1_5 (the 1-5 package)"]
#[rustfmt::skip]
fn two_loop_three_units() {
    check_all(
        (2, 3),
        &[
            Case { graph: "0-1*** 0-2 0-3 1-2 1-3", s: 2, t: 3, powers: &[1, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["1+(1,1)+(1,1)+(1,1)"], &["1+(1+(1,1),1)+(1,1)"]] },
            Case { graph: "0-1*** 0-2 0-3 1-2 1-3", s: 2, t: 3, powers: &[1, 1, 1, 1, 1], numerator: "l1.q",
                arrangements: &[&["1+(1,1)+(1,1)+(1,1)"]] },
            Case { graph: "0-1*** 0-2 0-3 1-2 1-3", s: 2, t: 3, powers: &[-1, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["-1+(1,1)+(1,1)+(1,1)"]] },
            Case { graph: "0-4* 4-1** 0-2 0-3 1-2 1-3", s: 2, t: 3, powers: &[1, 0, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["1+(1,1)", "(1,1)+(1,1)"]] },
            Case { graph: "0-1 0-2 0-3* 1-2* 1-3*", s: 2, t: 3, powers: &[1, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["1+(1,1)", "1+(1,1)", "1+(1,1)"]] },
            Case { graph: "0-1 0-2 0-3* 1-2* 1-3*", s: 2, t: 3, powers: &[1, 2, 1, 1, 1], numerator: "l0.q",
                arrangements: &[&["1+(1,1)", "1+(1,1)", "1+(1,1)"]] },
            Case { graph: "0-1 0-2 0-3* 1-2* 1-3*", s: 2, t: 3, powers: &[1, 1, 0, 1, -1], numerator: "",
                arrangements: &[&["(1,1)", "1+(1,1)", "-1+(1,1)"]] },
            Case { graph: "0-1* 0-2 0-3 1-2 1-3**", s: 2, t: 3, powers: &[1, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["1+(1,1)", "1+(1,1)+(1,1)"], &["1+(1,1)", "1+(1+(1,1),1)"]] },
            Case { graph: "0-1 0-2 0-3* 1-2 1-3**", s: 2, t: 3, powers: &[1, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["1+(1,1)", "1+(1,1)+(1,1)"]] },
        ],
    );
}

/// (3,1): Custom nodes (an ε node and a caller), a Triangle node; the ε line
/// at n = 0, -1, -2; a tagged piece merged with a pinched neighbour.
#[test]
#[ignore = "external data: VAKINT_TEST_P_LIBRARY_1_5 (the 1-5 package)"]
#[rustfmt::skip]
fn three_loop_one_unit() {
    check_all(
        (3, 1),
        &[
            Case { graph: "0-1 0-3 0-5 1-2 1-4 2-3* 2-5 3-4", s: 4, t: 5, powers: &[1, 1, 1, 1, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["1+(1,1)"], &["(2,1)"]] },
            Case { graph: "0-1 0-3 0-5 1-2 1-4 2-3* 2-5 3-4", s: 4, t: 5, powers: &[1, 2, 1, 1, 1, 1, 1, 1], numerator: "l0.q",
                arrangements: &[&["1+(1,1)"]] },
            Case { graph: "0-1 0-3 0-5 1-2 1-4 2-3* 2-5 3-4", s: 4, t: 5, powers: &[1, 1, 1, 1, 1, 0, 1, 1], numerator: "",
                arrangements: &[&["(1,1)"]] },
            Case { graph: "0-1 0-3 0-5 1-2 1-4 2-3* 2-5 3-4", s: 4, t: 5, powers: &[1, 1, 1, 1, 1, -1, 1, 1], numerator: "l0.q",
                arrangements: &[&["-1+(1,1)"]] },
            Case { graph: "0-1 0-3 0-5 1-2 1-4 2-3* 2-5 3-4", s: 4, t: 5, powers: &[1, 1, 1, 1, 1, -2, 1, 1], numerator: "",
                arrangements: &[&["-2+(1,1)"]] },
            Case { graph: "0-1 0-3 0-5 1-2 1-4 2-6* 6-3 2-5 3-4", s: 4, t: 5, powers: &[1, 1, 1, 1, 1, 0, 1, 1, 1], numerator: "",
                arrangements: &[&["(1,1)"]] },
            Case { graph: "0-1 0-3 0-5 1-2 1-4 2-3* 2-4 3-5", s: 4, t: 5, powers: &[1, 1, 1, 1, 1, 1, 1, 1], numerator: "l1.l6",
                arrangements: &[&["1+(1,1)"]] },
            Case { graph: "0-1 0-2 0-5 1-2 1-4 2-3* 3-4 3-5", s: 4, t: 5, powers: &[1, 1, 1, 1, 1, 2, 1, 1], numerator: "",
                arrangements: &[&["2+(1,1)"], &["1+(2,1)"]] },
            Case { graph: "0-1 0-2 0-5 1-2 1-4 2-3* 3-4 3-5", s: 4, t: 5, powers: &[1, 1, 1, 1, 1, -1, 1, 1], numerator: "",
                arrangements: &[&["-1+(1,1)"]] },
        ],
    );
}

/// (3,2): two units on one line (Custom caller: chain and nested), on two
/// lines (Custom, Triangle), two tagged pieces merged, n = 0, -1; five-loop
/// explicit graphs.
#[test]
#[ignore = "external data: VAKINT_TEST_P_LIBRARY_1_5 (the 1-5 package)"]
#[rustfmt::skip]
fn three_loop_two_units() {
    check_all(
        (3, 2),
        &[
            Case { graph: "0-1 0-3 0-5 1-2 1-4 2-3** 2-4 3-5", s: 4, t: 5, powers: &[1, 1, 1, 1, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["1+(1,1)+(1,1)"], &["1+(1+(1,1),1)"]] },
            Case { graph: "0-1 0-3 0-5 1-2 1-4 2-3** 2-4 3-5", s: 4, t: 5, powers: &[1, 1, 1, 1, 1, 1, 1, 1], numerator: "l0.q",
                arrangements: &[&["1+(1,1)+(1,1)"], &["1+(1+(1,1),1)"]] },
            Case { graph: "0-1 0-3 0-5 1-2 1-4 2-3** 2-4 3-5", s: 4, t: 5, powers: &[1, 1, 1, 1, 1, 0, 1, 1], numerator: "",
                arrangements: &[&["(1,1)+(1,1)"]] },
            Case { graph: "0-1 0-3 0-5 1-2 1-4 2-3** 2-4 3-5", s: 4, t: 5, powers: &[1, 1, 1, 1, 1, -1, 1, 1], numerator: "",
                arrangements: &[&["-1+(1,1)+(1,1)"]] },
            Case { graph: "0-1 0-3 0-5 1-2 1-4 2-6* 6-3* 2-4 3-5", s: 4, t: 5, powers: &[1, 1, 1, 1, 1, 1, 0, 1, 1], numerator: "",
                arrangements: &[&["1+(1,1)", "(1,1)"]] },
            Case { graph: "0-1 0-2 0-3* 1-2* 1-5 2-4 3-4 3-5", s: 4, t: 5, powers: &[1, 1, 1, 1, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["1+(1,1)", "1+(1,1)"]] },
            Case { graph: "0-1 0-2 0-5 1-2 1-4 2-3* 3-4 3-5*", s: 4, t: 5, powers: &[1, 1, 1, 1, 1, 1, 1, 1], numerator: "l0.q",
                arrangements: &[&["1+(1,1)", "1+(1,1)"]] },
            Case { graph: "0-1 0-2 0-5 1-2 1-4 2-3* 3-4 3-5*", s: 4, t: 5, powers: &[1, 1, 1, 1, 1, 0, 1, -1], numerator: "",
                arrangements: &[&["(1,1)", "-1+(1,1)"]] },
        ],
    );
}

/// (4,1): the corner of a Custom caller (five-loop explicit graphs, a GGPPSS
/// master behind the terminal), and the Insertion and SeriesParallel nodes
/// with a numerator and n = -1, 0. The Triangle, Diamond and off-corner Custom
/// keys are separate tests (the last two take minutes).
#[test]
#[ignore = "external data: VAKINT_TEST_P_LIBRARY_1_5 (the 1-5 package)"]
#[rustfmt::skip]
fn four_loop_one_unit() {
    check_all(
        (4, 1),
        &[
            Case { graph: "0-1 0-3 0-7 1-2 1-6 2-3 2-5 3-4 4-5* 4-6 5-7", s: 6, t: 7, powers: &[1; 11], numerator: "",
                arrangements: &[&["1+(1,1)"], &["(2,1)"]] },
            Case { graph: "0-1 0-2 0-3 1-2* 1-6 2-5 3-4 3-4 3-6 4-5", s: 5, t: 6, powers: &[1; 10], numerator: "l1.q",
                arrangements: &[&["1+(1,1)"]] },
            Case { graph: "0-1 0-2 0-3 1-2* 1-6 2-5 3-4 3-4 3-6 4-5", s: 5, t: 6, powers: &[1, 1, 1, -1, 1, 1, 1, 1, 1, 1], numerator: "",
                arrangements: &[&["-1+(1,1)"]] },
            Case { graph: "0-1 0-2 0-3 0-4 1-2 1-2 2-4 3-4*", s: 3, t: 4, powers: &[1; 8], numerator: "l0.q",
                arrangements: &[&["1+(1,1)"], &["(2,1)"]] },
            Case { graph: "0-1 0-2 0-3 0-4 1-2 1-2 2-4 3-4*", s: 3, t: 4, powers: &[1, 1, 1, 1, 1, 1, 1, 0], numerator: "",
                arrangements: &[&["(1,1)"]] },
        ],
    );
}

/// (4,1), the Triangle node: corner, numerator, n = -1 on the ε line.
#[test]
#[ignore = "external data: VAKINT_TEST_P_LIBRARY_1_5 (the 1-5 package)"]
#[rustfmt::skip]
fn four_loop_one_unit_triangle() {
    check_all(
        (4, 1),
        &[
            Case { graph: "0-1 0-3 0-7 1-2 1-6 2-3 2-4 3-4 4-5* 5-6 5-7", s: 6, t: 7, powers: &[1; 11], numerator: "",
                arrangements: &[&["1+(1,1)"]] },
            Case { graph: "0-1 0-3 0-7 1-2 1-6 2-3 2-4 3-4 4-5* 5-6 5-7", s: 6, t: 7, powers: &[1; 11], numerator: "l1.q",
                arrangements: &[&["1+(1,1)"]] },
            Case { graph: "0-1 0-3 0-7 1-2 1-6 2-3 2-4 3-4 4-5* 5-6 5-7", s: 6, t: 7, powers: &[1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1], numerator: "",
                arrangements: &[&["-1+(1,1)"]] },
        ],
    );
}

/// (4,1), the Diamond node, with vakint's engine arithmetic (the corner of
/// this node did not return in 50 min under `Exact`): about 5 min on one core.
#[test]
#[ignore = "slow, external data: VAKINT_TEST_P_LIBRARY_1_5 (the 1-5 package)"]
#[rustfmt::skip]
fn four_loop_one_unit_diamond() {
    check_all_with(
        (4, 1),
        &[
            Case { graph: "0-3 0-6 0-7 1-2 1-4 1-7 2-3 2-5 3-4 4-5* 5-6", s: 6, t: 7, powers: &[1; 11], numerator: "",
                arrangements: &[&["1+(1,1)"], &["(2,1)"]] },
            Case { graph: "0-3 0-6 0-7 1-2 1-4 1-7 2-3 2-5 3-4 4-5* 5-6", s: 6, t: 7, powers: &[1; 11], numerator: "l1.q",
                arrangements: &[&["1+(1,1)"]] },
            Case { graph: "0-3 0-6 0-7 1-2 1-4 1-7 2-3 2-5 3-4 4-5* 5-6", s: 6, t: 7, powers: &[1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1], numerator: "",
                arrangements: &[&["(1,1)"]] },
        ],
        Arithmetic::Engine,
    );
}

/// (4,1), the Custom caller off its corner: a numerator, n = 0 and 2 on the
/// ε line, a dot elsewhere with a numerator; vakint's engine arithmetic (these
/// keys did not return in 30-50 min under `Exact`): about 12 min on one core.
#[test]
#[ignore = "slow, external data: VAKINT_TEST_P_LIBRARY_1_5 (the 1-5 package)"]
#[rustfmt::skip]
fn four_loop_one_unit_custom_off_corner() {
    check_all_with(
        (4, 1),
        &[
            Case { graph: "0-1 0-3 0-7 1-2 1-6 2-3 2-5 3-4 4-5* 4-6 5-7", s: 6, t: 7, powers: &[1; 11], numerator: "l0.q",
                arrangements: &[&["1+(1,1)"]] },
            Case { graph: "0-1 0-3 0-7 1-2 1-6 2-3 2-5 3-4 4-5* 4-6 5-7", s: 6, t: 7, powers: &[1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1], numerator: "",
                arrangements: &[&["(1,1)"]] },
            Case { graph: "0-1 0-3 0-7 1-2 1-6 2-3 2-5 3-4 4-5* 4-6 5-7", s: 6, t: 7, powers: &[1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1], numerator: "",
                arrangements: &[&["2+(1,1)"], &["1+(2,1)"]] },
            Case { graph: "0-1 0-3 0-7 1-2 1-6 2-3 2-5 3-4 4-5* 4-6 5-7", s: 6, t: 7, powers: &[1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1], numerator: "l0.q",
                arrangements: &[&["1+(1,1)"]] },
        ],
        Arithmetic::Engine,
    );
}
