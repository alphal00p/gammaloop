//! Deep-dive: does the on-shell-massless-leg fix (off-shell delta regularization) work as the number
//! of massless legs grows, and for massless INTERNAL lines? Each config runs in its OWN process
//! (arg = config index) so a panic in one doesn't poison Symbolica state for the others.
//! Gitignored, local-only.
//!   for i in $(seq 0 N); do cargo run --release --example massless_legs -p oneloop -- $i; done
//!
//! Prints, for the chosen config, either "PANIC" or the reduction as
//!   RESULT <name> | TERM coeff=( .. ) <MASTER args...>   (parseable, masters have numeric args)

use oneloop::masters::MasterIntegral;
use oneloop::symbols::S;
use oneloop::{Integral, IntegralFamily, Kinematics, Propagator, reduce};
use symbolica::atom::{Atom, AtomCore};
use symbolica::function;

fn frac(a: i64, b: i64) -> Atom {
    Atom::num(a) / Atom::num(b)
}
fn kq(j: usize) -> Atom {
    let q = symbolica::symbol!(format!("oneloop::q{}", j + 1));
    function!(S.dot, Atom::var(S.k), Atom::var(q))
}
fn kk() -> Atom {
    function!(S.dot, Atom::var(S.k), Atom::var(S.k))
}

fn fam(masses: Vec<Atom>, invs: Vec<Atom>, num: Atom) -> IntegralFamily {
    IntegralFamily {
        propagators: masses
            .into_iter()
            .map(|m| Propagator {
                momentum: Atom::Zero,
                mass_sq: m,
            })
            .collect(),
        isps: vec![],
        kinematics: Kinematics { invariants: invs },
        targets: vec![Integral {
            propagator_exponents: vec![1; 3],
            isp_exponents: vec![],
        }],
        numerator: num,
    }
}

fn master_str(m: &MasterIntegral) -> String {
    let s = |x: &Atom| x.expand().to_string();
    match m {
        MasterIntegral::Tadpole { m_sq } => format!("A0 {}", s(m_sq)),
        MasterIntegral::Bubble { p_sq, m1_sq, m2_sq } => {
            format!("B0 {} {} {}", s(p_sq), s(m1_sq), s(m2_sq))
        }
        MasterIntegral::Triangle {
            p1_sq,
            p2_sq,
            p12_sq,
            m1_sq,
            m2_sq,
            m3_sq,
        } => {
            format!(
                "C0 {} {} {} {} {} {}",
                s(p1_sq),
                s(p2_sq),
                s(p12_sq),
                s(m1_sq),
                s(m2_sq),
                s(m3_sq)
            )
        }
        MasterIntegral::Box { .. } => "BOX".into(),
    }
}

fn main() {
    gammalooprs::initialisation::initialise().expect("activate Symbolica license");
    let idx: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);

    // (name, masses, invariants [s01,s02,s12], numerator)
    let mm = || vec![Atom::num(1), Atom::num(1), Atom::num(1)]; // massive internal m^2=1
    let m0 = || vec![Atom::Zero, Atom::Zero, Atom::Zero]; // massless internal
    let s = |x| Atom::num(x);
    let configs: Vec<(&str, Vec<Atom>, Vec<Atom>, Atom)> = vec![
        // ---- 2 massless legs, MASSIVE internal (gg>h). Baseline (validated). ----
        (
            "2ml_massive_scalar",
            mm(),
            vec![Atom::Zero, Atom::Zero, frac(2, 5)],
            Atom::num(1),
        ),
        (
            "2ml_massive_rank1",
            mm(),
            vec![Atom::Zero, Atom::Zero, frac(2, 5)],
            kq(0),
        ),
        (
            "2ml_massive_rank2",
            mm(),
            vec![Atom::Zero, Atom::Zero, frac(2, 5)],
            &kq(0) * &kq(0),
        ),
        // ---- 3 massless legs, MASSIVE internal: C0(0,0,0;m,m,m), collinear point ----
        (
            "3ml_massive_scalar",
            mm(),
            vec![Atom::Zero, Atom::Zero, Atom::Zero],
            Atom::num(1),
        ),
        (
            "3ml_massive_rank1",
            mm(),
            vec![Atom::Zero, Atom::Zero, Atom::Zero],
            kq(0),
        ),
        (
            "3ml_massive_rank2",
            mm(),
            vec![Atom::Zero, Atom::Zero, Atom::Zero],
            &kq(0) * &kq(0),
        ),
        (
            "3ml_massive_kk",
            mm(),
            vec![Atom::Zero, Atom::Zero, Atom::Zero],
            kk(),
        ),
        // ---- 2 massless legs, MASSLESS internal (e+e->ddx): one off-shell leg s=1 ----
        (
            "2ml_massless_scalar",
            m0(),
            vec![Atom::Zero, Atom::Zero, s(1)],
            Atom::num(1),
        ),
        (
            "2ml_massless_rank1q1",
            m0(),
            vec![Atom::Zero, Atom::Zero, s(1)],
            kq(0),
        ),
        (
            "2ml_massless_rank1q2",
            m0(),
            vec![Atom::Zero, Atom::Zero, s(1)],
            kq(1),
        ),
        (
            "2ml_massless_rank2",
            m0(),
            vec![Atom::Zero, Atom::Zero, s(1)],
            &kq(0) * &kq(0),
        ),
        // ---- ALL massless (external + internal): fully scaleless triangle ----
        (
            "allml_scalar",
            m0(),
            vec![Atom::Zero, Atom::Zero, Atom::Zero],
            Atom::num(1),
        ),
        (
            "allml_rank2",
            m0(),
            vec![Atom::Zero, Atom::Zero, Atom::Zero],
            &kq(0) * &kq(0),
        ),
        // ---- PATH-INDEPENDENCE check for 3-massless-leg rank-2 (should -> 0): reduce at small
        //      OFF-SHELL invariants (direct path, NO regularization), symmetric vs asymmetric.
        //      If 0 is genuine (not a single-delta artifact) the value shrinks toward 0 both ways.
        (
            "3ml_r2_sym_1e-1",
            mm(),
            vec![frac(1, 10), frac(1, 10), frac(1, 10)],
            &kq(0) * &kq(0),
        ),
        (
            "3ml_r2_sym_1e-2",
            mm(),
            vec![frac(1, 100), frac(1, 100), frac(1, 100)],
            &kq(0) * &kq(0),
        ),
        (
            "3ml_r2_asym_1e-2",
            mm(),
            vec![frac(3, 100), frac(2, 100), frac(1, 100)],
            &kq(0) * &kq(0),
        ),
        (
            "3ml_r2_asym_1e-3",
            mm(),
            vec![frac(3, 1000), frac(2, 1000), frac(1, 1000)],
            &kq(0) * &kq(0),
        ),
    ];

    if idx >= configs.len() {
        println!("MAXIDX {}", configs.len());
        return;
    }
    let (name, masses, invs, num) = &configs[idx];
    let f = fam(masses.clone(), invs.clone(), num.clone());
    let r = reduce(&f); // fix regularizes zero invariants automatically
    if r.terms.is_empty() {
        println!("RESULT {name} | SCALELESS_OR_ZERO (0 terms)");
    }
    for (c, m) in &r.terms {
        println!(
            "RESULT {name} | TERM coeff=( {} ) {}",
            c.expand(),
            master_str(m)
        );
    }
}
