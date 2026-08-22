//! Reduce every dot(k,.) monomial (rank <= 6, delta-regularized on-shell massive triangle) for wloop_assemble.py (see docs/10).

use oneloop::masters::MasterIntegral;
use oneloop::symbols::S;
use oneloop::{Integral, IntegralFamily, Kinematics, Propagator, reduce};
use symbolica::atom::{Atom, AtomCore};
use symbolica::function;

fn family(s: &Atom, mwsq: &Atom, delta: &Atom, numerator: Atom) -> IntegralFamily {
    IntegralFamily {
        propagators: (0..3)
            .map(|_| Propagator { momentum: Atom::Zero, mass_sq: mwsq.clone() })
            .collect(),
        isps: vec![],
        // off-shell regularization q1^2=q2^2=delta; the 1/delta poles cancel in the delta->0 assembly
        kinematics: Kinematics { invariants: vec![delta.clone(), s.clone(), delta.clone()] },
        targets: vec![Integral { propagator_exponents: vec![1, 1, 1], isp_exponents: vec![] }],
        numerator,
    }
}

fn master_line(m: &MasterIntegral) -> String {
    let s = |x: &Atom| x.expand().to_string();
    match m {
        MasterIntegral::Tadpole { m_sq } => format!("A0 {}", s(m_sq)),
        MasterIntegral::Bubble { p_sq, m1_sq, m2_sq } => format!("B0 {} {} {}", s(p_sq), s(m1_sq), s(m2_sq)),
        MasterIntegral::Triangle { p1_sq, p2_sq, p12_sq, m1_sq, m2_sq, m3_sq } =>
            format!("C0 {} {} {} {} {} {}", s(p1_sq), s(p2_sq), s(p12_sq), s(m1_sq), s(m2_sq), s(m3_sq)),
        MasterIntegral::Box { .. } => unreachable!("triangle cannot emit a box"),
    }
}

fn main() {
    gammalooprs::initialisation::initialise().expect("activate Symbolica license");
    let a: Vec<i64> = std::env::args().skip(1).filter_map(|x| x.parse().ok()).collect();
    let s = Atom::num(a[0]) / Atom::num(a[1]);
    let mwsq = Atom::num(a[2]) / Atom::num(a[3]);
    let delta = Atom::num(a[4]) / Atom::num(a[5]);
    println!("S {} {}  MWSQ {} {}  DELTA {} {}", a[0], a[1], a[2], a[3], a[4], a[5]);

    let k = Atom::var(S.k);
    let q1 = symbolica::symbol!("oneloop::q1");
    let q2 = symbolica::symbol!("oneloop::q2");
    let ll = function!(S.dot, k.clone(), k.clone());
    let lq1 = function!(S.dot, k.clone(), Atom::var(q1));
    let lq2 = function!(S.dot, k.clone(), Atom::var(q2));

    for aa in 0..=3 {
        for bb in 0..=(6 - 2 * aa) {
            for cc in 0..=(6 - 2 * aa - bb) {
                let mut num = Atom::num(1);
                for _ in 0..aa { num = &num * &ll; }
                for _ in 0..bb { num = &num * &lq1; }
                for _ in 0..cc { num = &num * &lq2; }
                let r = reduce(&family(&s, &mwsq, &delta, num));
                println!("MONO {aa} {bb} {cc}");
                for (c, m) in &r.terms {
                    println!("TERM coeff=( {} ) {}", c.expand(), master_line(m));
                }
                println!("ENDMONO");
            }
        }
    }
}
