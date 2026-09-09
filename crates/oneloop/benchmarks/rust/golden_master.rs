//! Golden-master characterization of reduce() (gitignored).  Reduces a battery of families and
//! prints every (master, coefficient) evaluated at a fixed numeric point -- so two algebraically
//! different but equal reductions produce the SAME dump.  Capture on the committed code, then
//! require a byte-identical dump after a behaviour-preserving refactor:
//!   cargo run --release --example golden_master -p oneloop > /tmp/golden.txt
//!   ...refactor...
//!   cargo run --release --example golden_master -p oneloop | diff /tmp/golden.txt -

use oneloop::masters::MasterIntegral;
use oneloop::symbols::S;
use oneloop::{Integral, IntegralFamily, Kinematics, Propagator, reduce};
use symbolica::atom::{Atom, AtomCore};
use symbolica::{function, symbol};

fn v(name: &str) -> Atom {
    Atom::var(symbol!(format!("oneloop::{name}")))
}
fn dot_ll() -> Atom {
    function!(S.dot, Atom::var(S.k), Atom::var(S.k))
}
fn dot_lq(j: usize) -> Atom {
    function!(
        S.dot,
        Atom::var(S.k),
        Atom::var(symbol!(format!("oneloop::q{}", j + 1)))
    )
}

// A generic numeric point (distinct rationals, away from any pole): symbol name -> num/den.
const POINT: &[(&str, i64, i64)] = &[
    ("msq", 5, 7),
    ("m1sq", 5, 7),
    ("m2sq", 9, 11),
    ("m3sq", 13, 4),
    ("m4sq", 17, 6),
    ("m5sq", 19, 5),
    ("m6sq", 23, 8),
    ("psq", 3, 2),
    ("s1", 7, 5),
    ("s2", 11, 9),
    ("s3", 13, 6),
    ("p1", 3, 2),
    ("p2", 5, 3),
    ("p3", 7, 4),
    ("p4", 9, 5),
    ("sinv", 11, 7),
    ("tinv", 13, 8),
];

fn subst(a: &Atom) -> Atom {
    let mut out = a
        .replace(Atom::var(S.d).to_pattern())
        .with(Atom::num(11) / Atom::num(3));
    for (name, n, d) in POINT {
        out = out
            .replace(v(name).to_pattern())
            .with(Atom::num(*n) / Atom::num(*d));
    }
    out.cancel()
}

fn record(m: &MasterIntegral) -> String {
    let s = |x: &Atom| subst(x).to_string();
    match m {
        MasterIntegral::Tadpole { m_sq } => format!("A0[{}]", s(m_sq)),
        MasterIntegral::Bubble { p_sq, m1_sq, m2_sq } => {
            format!("B0[{};{},{}]", s(p_sq), s(m1_sq), s(m2_sq))
        }
        MasterIntegral::Triangle {
            p1_sq,
            p2_sq,
            p12_sq,
            m1_sq,
            m2_sq,
            m3_sq,
        } => format!(
            "C0[{},{},{};{},{},{}]",
            s(p1_sq),
            s(p2_sq),
            s(p12_sq),
            s(m1_sq),
            s(m2_sq),
            s(m3_sq)
        ),
        MasterIntegral::Box {
            p1_sq,
            p2_sq,
            p3_sq,
            p4_sq,
            s: si,
            t,
            m1_sq,
            m2_sq,
            m3_sq,
            m4_sq,
        } => {
            format!(
                "D0[{},{},{},{},{},{};{},{},{},{}]",
                s(p1_sq),
                s(p2_sq),
                s(p3_sq),
                s(p4_sq),
                s(si),
                s(t),
                s(m1_sq),
                s(m2_sq),
                s(m3_sq),
                s(m4_sq)
            )
        }
    }
}

fn emit(label: &str, fam: &IntegralFamily) {
    // Aggregate coefficients per master (numeric point) so the dump reflects the mathematical
    // result -- insensitive to term order and to how masters happen to merge/represent.
    use std::collections::BTreeMap;
    let r = reduce(fam);
    let mut agg: BTreeMap<String, Atom> = BTreeMap::new();
    for (c, m) in &r.terms {
        let e = agg.entry(record(m)).or_insert(Atom::Zero);
        *e = &*e + &subst(c);
    }
    println!("{label}");
    for (m, c) in agg {
        let c = c.cancel();
        if c != Atom::Zero {
            println!("    {c} * {m}");
        }
    }
}

fn fam(masses: &[&str], invariants: &[&str], exps: Vec<i32>, numerator: Atom) -> IntegralFamily {
    IntegralFamily {
        propagators: masses
            .iter()
            .map(|m| Propagator {
                momentum: Atom::Zero,
                mass_sq: v(m),
            })
            .collect(),
        isps: vec![],
        kinematics: Kinematics {
            invariants: invariants.iter().map(|s| v(s)).collect(),
        },
        targets: vec![Integral {
            propagator_exponents: exps,
            isp_exponents: vec![],
        }],
        numerator,
    }
}

fn main() {
    gammalooprs::initialisation::initialise().expect("activate Symbolica license");
    let one = || Atom::num(1);

    // Tadpole (N=1)
    emit("tad/scalar", &fam(&["msq"], &[], vec![1], one()));
    emit("tad/dot2", &fam(&["msq"], &[], vec![2], one()));
    emit("tad/dot3", &fam(&["msq"], &[], vec![3], one()));
    emit("tad/num_ll", &fam(&["msq"], &[], vec![1], dot_ll()));
    emit(
        "tad/num_ll2",
        &fam(&["msq"], &[], vec![1], &dot_ll() * &dot_ll()),
    );
    emit("tad/num_ll_dot2", &fam(&["msq"], &[], vec![2], dot_ll()));

    // Bubble (N=2)
    emit(
        "bub/scalar",
        &fam(&["m1sq", "m2sq"], &["psq"], vec![1, 1], one()),
    );
    emit(
        "bub/dot21",
        &fam(&["m1sq", "m2sq"], &["psq"], vec![2, 1], one()),
    );
    emit(
        "bub/dot31",
        &fam(&["m1sq", "m2sq"], &["psq"], vec![3, 1], one()),
    );
    emit(
        "bub/num_q1",
        &fam(&["m1sq", "m2sq"], &["psq"], vec![1, 1], dot_lq(0)),
    );
    emit(
        "bub/num_q1q1",
        &fam(
            &["m1sq", "m2sq"],
            &["psq"],
            vec![1, 1],
            &dot_lq(0) * &dot_lq(0),
        ),
    );
    emit(
        "bub/num_ll",
        &fam(&["m1sq", "m2sq"], &["psq"], vec![1, 1], dot_ll()),
    );
    emit(
        "bub/dot21_num_q1",
        &fam(&["m1sq", "m2sq"], &["psq"], vec![2, 1], dot_lq(0)),
    );

    // Triangle (N=3)
    let tri_m = ["m1sq", "m2sq", "m3sq"];
    let tri_s = ["s1", "s2", "s3"];
    emit("tri/scalar", &fam(&tri_m, &tri_s, vec![1, 1, 1], one()));
    emit("tri/dot211", &fam(&tri_m, &tri_s, vec![2, 1, 1], one()));
    emit("tri/dot121", &fam(&tri_m, &tri_s, vec![1, 2, 1], one()));
    emit("tri/dot221", &fam(&tri_m, &tri_s, vec![2, 2, 1], one()));
    emit("tri/pinch", &fam(&tri_m, &tri_s, vec![1, 0, 1], one()));
    emit("tri/num_q1", &fam(&tri_m, &tri_s, vec![1, 1, 1], dot_lq(0)));
    emit(
        "tri/num_q1q2",
        &fam(&tri_m, &tri_s, vec![1, 1, 1], &dot_lq(0) * &dot_lq(1)),
    );
    emit("tri/num_ll", &fam(&tri_m, &tri_s, vec![1, 1, 1], dot_ll()));
    emit(
        "tri/num_ll_q1",
        &fam(&tri_m, &tri_s, vec![1, 1, 1], &dot_ll() * &dot_lq(0)),
    );
    emit(
        "tri/dot211_num_q1",
        &fam(&tri_m, &tri_s, vec![2, 1, 1], dot_lq(0)),
    );

    // Box (N=4)
    let bx_m = ["m1sq", "m2sq", "m3sq", "m4sq"];
    let bx_s = ["p1", "p2", "p3", "p4", "sinv", "tinv"];
    emit("box/scalar", &fam(&bx_m, &bx_s, vec![1, 1, 1, 1], one()));
    emit("box/dot2111", &fam(&bx_m, &bx_s, vec![2, 1, 1, 1], one()));
    emit("box/pinch", &fam(&bx_m, &bx_s, vec![1, 0, 1, 1], one()));
    emit(
        "box/num_q1",
        &fam(&bx_m, &bx_s, vec![1, 1, 1, 1], dot_lq(0)),
    );
    emit(
        "box/num_q1q2q3",
        &fam(
            &bx_m,
            &bx_s,
            vec![1, 1, 1, 1],
            &(&dot_lq(0) * &dot_lq(1)) * &dot_lq(2),
        ),
    );
    emit("box/num_ll", &fam(&bx_m, &bx_s, vec![1, 1, 1, 1], dot_ll()));
    emit(
        "box/dot2111_num_ll",
        &fam(&bx_m, &bx_s, vec![2, 1, 1, 1], dot_ll()),
    );

    // Pentagon / hexagon scalar (N>4 path, numeric invariants) -- locks the untouched high-N path.
    let pent_m = ["m1sq", "m2sq", "m3sq", "m4sq", "m5sq"];
    let pent_inv: Vec<Atom> = [3, 5, 7, 9, 4, 6, 8, 5, 7, 6]
        .iter()
        .map(|&x| Atom::num(x))
        .collect();
    emit(
        "pent/scalar",
        &IntegralFamily {
            propagators: pent_m
                .iter()
                .map(|m| Propagator {
                    momentum: Atom::Zero,
                    mass_sq: v(m),
                })
                .collect(),
            isps: vec![],
            kinematics: Kinematics {
                invariants: pent_inv,
            },
            targets: vec![Integral {
                propagator_exponents: vec![1; 5],
                isp_exponents: vec![],
            }],
            numerator: one(),
        },
    );
}
