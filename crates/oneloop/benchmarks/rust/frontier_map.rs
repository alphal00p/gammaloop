//! Map the reducer's on-shell-massless-leg frontier: for a MASSIVE-internal triangle/box, vary the
//! number of on-shell (zero pairwise-invariant) legs and the numerator rank, catching panics so one
//! run reports the full OK/PANIC map. Gitignored, local-only.
//!   cargo run --release --example frontier_map -p oneloop

use oneloop::symbols::S;
use oneloop::{Integral, IntegralFamily, Kinematics, Propagator, reduce};
use symbolica::atom::Atom;
use symbolica::function;

fn fam(n: usize, exps: Vec<i32>, invs: &[(i64, i64)], numerator: Atom) -> IntegralFamily {
    IntegralFamily {
        propagators: (0..n)
            .map(|_| Propagator {
                momentum: Atom::Zero,
                mass_sq: Atom::num(1),
            }) // massive internal m^2=1
            .collect(),
        isps: vec![],
        kinematics: Kinematics {
            invariants: invs
                .iter()
                .map(|(a, b)| Atom::num(*a) / Atom::num(*b))
                .collect(),
        },
        targets: vec![Integral {
            propagator_exponents: exps,
            isp_exponents: vec![],
        }],
        numerator,
    }
}

fn kq(j: usize) -> Atom {
    let q = symbolica::symbol!(format!("oneloop::q{}", j + 1));
    function!(S.dot, Atom::var(S.k), Atom::var(q))
}
fn kk() -> Atom {
    function!(S.dot, Atom::var(S.k), Atom::var(S.k))
}

fn try_reduce(label: &str, f: IntegralFamily) {
    let res = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| reduce(&f)));
    match res {
        Ok(r) => println!("  {label:28} OK   ({} terms)", r.terms.len()),
        Err(e) => {
            let msg = e
                .downcast_ref::<String>()
                .map(String::as_str)
                .or_else(|| e.downcast_ref::<&str>().copied())
                .unwrap_or("panic");
            let short = msg.split(':').next().unwrap_or(msg);
            println!("  {label:28} PANIC ({short})");
        }
    }
}

fn main() {
    gammalooprs::initialisation::initialise().expect("activate Symbolica license");
    std::panic::set_hook(Box::new(|_| {})); // silence panic spew; we report via catch_unwind

    let q1q1 = &kq(0) * &kq(0);

    println!("== TRIANGLE (massive internal m^2=1): rank-2 (k.q1)^2 vs # of on-shell legs ==");
    // lex invariants [s01,s02,s12]; a zero invariant = an on-shell (massless) leg.
    try_reduce(
        "0 on-shell (all off)",
        fam(
            3,
            vec![1; 3],
            &[(30, 100), (20, 100), (40, 100)],
            q1q1.clone(),
        ),
    );
    try_reduce(
        "1 on-shell",
        fam(3, vec![1; 3], &[(0, 1), (20, 100), (40, 100)], q1q1.clone()),
    );
    try_reduce(
        "2 on-shell",
        fam(3, vec![1; 3], &[(0, 1), (0, 1), (40, 100)], q1q1.clone()),
    );
    try_reduce(
        "3 on-shell",
        fam(3, vec![1; 3], &[(0, 1), (0, 1), (0, 1)], q1q1.clone()),
    );

    println!("\n== TRIANGLE: rank by rank at 2 on-shell legs [0,0,0.4] ==");
    let onsh = [(0, 1), (0, 1), (40, 100)];
    try_reduce("scalar", fam(3, vec![1; 3], &onsh, Atom::num(1)));
    try_reduce("rank-1 dot(k,q1)", fam(3, vec![1; 3], &onsh, kq(0)));
    try_reduce("rank-1 k^2", fam(3, vec![1; 3], &onsh, kk()));
    try_reduce(
        "rank-2 (k.q1)^2",
        fam(3, vec![1; 3], &onsh, &kq(0) * &kq(0)),
    );
    try_reduce(
        "rank-2 (k.q1)(k.q2)",
        fam(3, vec![1; 3], &onsh, &kq(0) * &kq(1)),
    );
    try_reduce(
        "raised power [2,1,1]",
        fam(3, vec![2, 1, 1], &onsh, Atom::num(1)),
    );

    println!("\n== BOX (massive internal): rank-2 (k.q1)^2, off-shell vs 2 on-shell legs ==");
    let boxoff = [
        (30, 100),
        (25, 100),
        (20, 100),
        (35, 100),
        (28, 100),
        (22, 100),
    ];
    let boxon = [(0, 1), (25, 100), (20, 100), (0, 1), (28, 100), (22, 100)];
    try_reduce(
        "box all off-shell",
        fam(4, vec![1; 4], &boxoff, q1q1.clone()),
    );
    try_reduce("box 2 on-shell", fam(4, vec![1; 4], &boxon, q1q1.clone()));
}
