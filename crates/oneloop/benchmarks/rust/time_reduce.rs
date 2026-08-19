//! Per-topology timing of reduce() — absolute reduction speed for triangle/box/pentagon, scalar and
//! tensor. Gitignored, local-only.   cargo run --release --example time_reduce -p oneloop

use oneloop::symbols::S;
use oneloop::{Integral, IntegralFamily, Kinematics, Propagator, reduce};
use std::time::Instant;
use symbolica::atom::Atom;
use symbolica::function;

fn fam(n: usize, invs: &[i64], numerator: Atom) -> IntegralFamily {
    let m = |x: i64| Atom::num(x) / Atom::num(100);
    IntegralFamily {
        propagators: (0..n)
            .map(|i| Propagator {
                momentum: Atom::Zero,
                mass_sq: m(100 + 7 * i as i64),
            })
            .collect(),
        isps: vec![],
        kinematics: Kinematics {
            invariants: invs.iter().map(|&s| m(s)).collect(),
        },
        targets: vec![Integral {
            propagator_exponents: vec![1; n],
            isp_exponents: vec![],
        }],
        numerator,
    }
}

fn kq(j: usize) -> Atom {
    let q = symbolica::symbol!(format!("oneloop::q{}", j + 1));
    function!(S.dot, Atom::var(S.k), Atom::var(q))
}

fn time_it(label: &str, build: impl Fn() -> IntegralFamily, n: u32) {
    let _ = reduce(&build()); // warm up
    let t0 = Instant::now();
    let mut terms = 0;
    for _ in 0..n {
        terms = reduce(&build()).terms.len();
    }
    let ms = t0.elapsed().as_secs_f64() * 1e3 / f64::from(n);
    println!("  {label:34} {ms:8.3} ms/reduction   ({terms} terms)");
}

fn main() {
    gammalooprs::initialisation::initialise().expect("activate Symbolica license");
    // generic off-shell massive invariants (spacelike, non-degenerate)
    let tri = [-9, -13, -2];
    let bx = [-30, -25, -20, -35, -28, -22];
    let pent: Vec<i64> = [-9, -13, -2, -11, -7, -15, -6, -10, -8, -5].to_vec();

    println!("oneloop reduce() per-topology timing (generic massive kinematics):");
    time_it("triangle  scalar", || fam(3, &tri, Atom::num(1)), 500);
    time_it("triangle  rank-1 dot(k,q1)", || fam(3, &tri, kq(0)), 500);
    time_it(
        "triangle  rank-2 (k.q1)^2",
        || fam(3, &tri, &kq(0) * &kq(0)),
        500,
    );
    time_it("box       scalar", || fam(4, &bx, Atom::num(1)), 300);
    time_it("box       rank-1 dot(k,q1)", || fam(4, &bx, kq(0)), 300);
    time_it(
        "box       rank-2 (k.q1)^2",
        || fam(4, &bx, &kq(0) * &kq(0)),
        200,
    );
    time_it("pentagon  scalar", || fam(5, &pent, Atom::num(1)), 200);
    time_it("pentagon  rank-1 dot(k,q1)", || fam(5, &pent, kq(0)), 100);
}
