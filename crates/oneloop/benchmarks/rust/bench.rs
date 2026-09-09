//! Local-only stress + timing sweep of the full oneloop reducer (gitignored).
//!
//! Reduces every topology (tadpole..heptagon) at ranks 0-4, scalar and dotted, over
//! physical (4D-offset, degenerate-Gram for N>=6) kinematics.  For each reduction it
//! catches panics, checks every coefficient is finite (Symbolica renders a divide-by-zero
//! as the non-ASCII glyph), checks masters are valid, and times it.  Run with:
//!   cargo run --release --example bench -p oneloop

use std::panic::{AssertUnwindSafe, catch_unwind};
use std::time::Instant;

use oneloop::masters::MasterIntegral;
use oneloop::symbols::S;
use oneloop::{Integral, IntegralFamily, Kinematics, Propagator, reduce};
use symbolica::atom::Atom;
use symbolica::{function, symbol};

fn dot_ll() -> Atom {
    function!(S.dot, Atom::var(S.k), Atom::var(S.k))
}
fn dot_lq(j: usize) -> Atom {
    let q = symbol!(format!("oneloop::q{}", j + 1));
    function!(S.dot, Atom::var(S.k), Atom::var(q))
}

// Integer 4D offsets for the N-gon lines (r_1 = 0); 5+ of them are dependent in 4D, so the
// external Gram is rank-deficient for N >= 6 -- the degenerate regime we want to stress.
const OFFSETS: [[i64; 4]; 7] = [
    [0, 0, 0, 0],
    [2, 0, 0, 0],
    [1, 3, 0, 0],
    [0, 2, 4, 0],
    [1, 0, 2, 5],
    [3, 1, 4, 2],
    [2, 3, 1, 4],
];

// Lexicographic pairwise invariants s_ij = -(r_i - r_j)^2 for lines 0..n.
#[allow(clippy::needless_range_loop)]
fn invariants(n: usize) -> Vec<Atom> {
    let mut out = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            let d: i64 = (0..4).map(|k| (OFFSETS[i][k] - OFFSETS[j][k]).pow(2)).sum();
            out.push(Atom::num(-d));
        }
    }
    out
}

fn family(n: usize, exps: Vec<i32>, numerator: Atom) -> IntegralFamily {
    IntegralFamily {
        propagators: (1..=n)
            .map(|i| Propagator {
                momentum: Atom::Zero,
                mass_sq: Atom::num(i as i64),
            })
            .collect(),
        isps: vec![],
        kinematics: Kinematics {
            invariants: invariants(n),
        },
        targets: vec![Integral {
            propagator_exponents: exps,
            isp_exponents: vec![],
        }],
        numerator,
    }
}

// Scalar plus a spread of rank-1..4 numerators from the n_ext available external dots.
fn numerators(n_ext: usize) -> Vec<(String, Atom)> {
    let mut v = vec![("scalar".to_string(), Atom::num(1))];
    if n_ext == 0 {
        return v; // tadpole: a numerator with < 2 propagators is a todo! guard, skip.
    }
    let last = n_ext - 1;
    v.push(("r1 l.q1".into(), dot_lq(0)));
    v.push(("r2 l.l".into(), dot_ll()));
    v.push(("r2 l.q1 l.qN".into(), &dot_lq(0) * &dot_lq(last)));
    v.push(("r3 l.l l.q1".into(), &dot_ll() * &dot_lq(0)));
    v.push(("r4 (l.l)^2".into(), &dot_ll() * &dot_ll()));
    v.push((
        "r4 (l.q1)^4".into(),
        &(&dot_lq(0) * &dot_lq(0)) * &(&dot_lq(0) * &dot_lq(0)),
    ));
    v
}

struct Row {
    label: String,
    terms: usize,
    ms: f64,
    status: String,
}

fn run_one(label: String, fam: &IntegralFamily) -> Row {
    let t = Instant::now();
    let outcome = catch_unwind(AssertUnwindSafe(|| reduce(fam)));
    let ms = t.elapsed().as_secs_f64() * 1e3;
    match outcome {
        Err(e) => {
            let msg = e
                .downcast_ref::<&str>()
                .map(|s| s.to_string())
                .or_else(|| e.downcast_ref::<String>().cloned())
                .unwrap_or_else(|| "panic".into());
            Row {
                label,
                terms: 0,
                ms,
                status: format!("PANIC: {msg}"),
            }
        }
        Ok(r) => {
            let mut problems = Vec::new();
            if r.terms.is_empty() {
                problems.push("empty result".to_string());
            }
            for (c, m) in &r.terms {
                if !c.to_string().is_ascii() {
                    problems.push(format!("non-finite coeff on {m:?}"));
                    break;
                }
                if !matches!(
                    m,
                    MasterIntegral::Tadpole { .. }
                        | MasterIntegral::Bubble { .. }
                        | MasterIntegral::Triangle { .. }
                        | MasterIntegral::Box { .. }
                ) {
                    problems.push("unexpected master".into());
                    break;
                }
            }
            let status = if problems.is_empty() {
                "ok".into()
            } else {
                format!("ERROR: {}", problems.join("; "))
            };
            Row {
                label,
                terms: r.terms.len(),
                ms,
                status,
            }
        }
    }
}

fn emit(r: Row, rows: &mut Vec<Row>) {
    println!(
        "{:<34} {:>6} {:>9.2}ms  {}",
        r.label, r.terms, r.ms, r.status
    );
    let _ = std::io::Write::flush(&mut std::io::stdout());
    rows.push(r);
}

fn main() {
    gammalooprs::initialisation::initialise().expect("activate Symbolica license");
    std::panic::set_hook(Box::new(|_| {})); // catch_unwind reports; keep stderr clean

    let names = [
        "tadpole", "bubble", "triangle", "box", "pentagon", "hexagon", "heptagon",
    ];
    println!("{:<34} {:>6} {:>10}  status", "case", "terms", "time");
    println!("{}", "-".repeat(72));
    let mut rows = Vec::new();
    for n in 1..=7usize {
        let n_ext = n - 1;
        for (nlabel, num) in numerators(n_ext) {
            // Cap the expensive tail: skip rank 3-4 for N >= 6 (atom explosion).
            if n >= 6 && (nlabel.starts_with("r3") || nlabel.starts_with("r4")) {
                continue;
            }
            let label = format!("{} [{}]", names[n - 1], nlabel);
            emit(
                run_one(label, &family(n, vec![1; n], num.clone())),
                &mut rows,
            );
            // dotted (raised power on line 1); skip N>=7 (det(Y)=0 wall + minutes-long).
            if matches!(nlabel.as_str(), "scalar" | "r1 l.q1") && n <= 6 {
                let mut exps = vec![1; n];
                exps[0] = 2;
                let label = format!("{} [{} dotted]", names[n - 1], nlabel);
                emit(run_one(label, &family(n, exps, num)), &mut rows);
            }
        }
    }

    let _ = std::panic::take_hook();
    let errors: Vec<&Row> = rows.iter().filter(|r| r.status != "ok").collect();
    let mut slow: Vec<&Row> = rows.iter().collect();
    slow.sort_by(|a, b| b.ms.partial_cmp(&a.ms).unwrap());
    println!("\n{} reductions, {} error(s).", rows.len(), errors.len());
    println!("slowest:");
    for r in slow.iter().take(5) {
        println!("  {:>9.2}ms  {}", r.ms, r.label);
    }
    if !errors.is_empty() {
        println!("errors:");
        for r in &errors {
            println!("  {}  ->  {}", r.label, r.status);
        }
    }
}
