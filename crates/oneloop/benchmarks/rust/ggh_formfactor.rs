//! Emit the gg->h massive-top-triangle reductions for the form-factor benchmark.
//!
//! End-to-end validation of the reducer on a *physical* loop amplitude: gluon-fusion
//! Higgs production through a top loop.  The transverse-projected numerator -- the
//! g_{mu nu} contraction of the top-loop Dirac trace, which extracts the (gauge-
//! invariant) form factor -- is, in d dimensions,
//!
//!   g_{mu nu} N^{mu nu} = 4m[(4-d) k^2 + (6-2d) (k.q1) + 2 (k.q2) + (2-d) s/2] + 4 m^3 d
//!
//! with the on-shell gluons q1^2 = q2^2 = 0 and q1.q2 = s/2 (verified to 14 digits
//! against explicit Dirac matrices in ggh_formfactor.py).  Here we reduce each
//! dot(k,.) monomial {1, k.q1, k.q2, k^2} for the on-shell massive triangle and print
//!     coeff(d)  *  master(numeric args)
//! for ggh_formfactor.py, which evaluates the masters with OneLOop, assembles
//! g_{mu nu} M(eps), takes the finite part, and forms
//!     A_{1/2} = g_{mu nu} M / [(d-1) q1.q2]
//! to compare against the analytic 2/tau^2 [tau + (tau-1) f(tau)] (= 1.376 at m_t=173).
//!
//!   cargo run --release --example ggh_formfactor -p oneloop -- [s] [mtsq]
//! Defaults: s = 125^2 = 15625 (m_H^2), mtsq = 173^2 = 29929 (m_t^2).

use oneloop::masters::MasterIntegral;
use oneloop::symbols::S;
use oneloop::{Integral, IntegralFamily, Kinematics, Propagator, reduce};
use symbolica::atom::{Atom, AtomCore};
use symbolica::function;

// gg->h top triangle: 3 top propagators (mass^2 = mtsq); external legs g(q1), g(q2),
// H(q1+q2).  Propagator momenta k, k+q1, k+q1+q2, so the lex pairwise invariants are
//   s_01 = q1^2 = 0,  s_02 = (q1+q2)^2 = s,  s_12 = q2^2 = 0.
fn ggh_family(s: &Atom, mtsq: &Atom, numerator: Atom) -> IntegralFamily {
    IntegralFamily {
        propagators: (0..3)
            .map(|_| Propagator {
                momentum: Atom::Zero,
                mass_sq: mtsq.clone(),
            })
            .collect(),
        isps: vec![],
        kinematics: Kinematics {
            invariants: vec![Atom::Zero, s.clone(), Atom::Zero],
        },
        targets: vec![Integral {
            propagator_exponents: vec![1, 1, 1],
            isp_exponents: vec![],
        }],
        numerator,
    }
}

fn master_line(m: &MasterIntegral) -> String {
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
        } => format!(
            "C0 {} {} {} {} {} {}",
            s(p1_sq),
            s(p2_sq),
            s(p12_sq),
            s(m1_sq),
            s(m2_sq),
            s(m3_sq)
        ),
        MasterIntegral::Box { .. } => unreachable!("triangle reduction cannot emit a box"),
    }
}

fn emit(label: &str, s: &Atom, mtsq: &Atom, numerator: Atom) {
    let fam = ggh_family(s, mtsq, numerator);
    let r = reduce(&fam);
    println!("NUM {label}");
    for (c, m) in &r.terms {
        // coeff is a rational function of d only (kinematics are numeric); print for sympy.
        println!("TERM coeff=( {} ) {}", c.expand(), master_line(m));
    }
    println!("ENDNUM");
}

fn main() {
    gammalooprs::initialisation::initialise().expect("activate Symbolica license");
    let args: Vec<String> = std::env::args().collect();
    let s_val: i64 = args.get(1).and_then(|x| x.parse().ok()).unwrap_or(15625);
    let mtsq_val: i64 = args.get(2).and_then(|x| x.parse().ok()).unwrap_or(29929);
    let s = Atom::num(s_val);
    let mtsq = Atom::num(mtsq_val);
    println!("S {s_val}");
    println!("MTSQ {mtsq_val}");

    let k = Atom::var(S.k);
    let q1 = symbolica::symbol!("oneloop::q1");
    let q2 = symbolica::symbol!("oneloop::q2");
    let kq1 = function!(S.dot, k.clone(), Atom::var(q1));
    let kq2 = function!(S.dot, k.clone(), Atom::var(q2));
    emit("one", &s, &mtsq, Atom::num(1));
    emit("lq1", &s, &mtsq, kq1.clone());
    emit("lq2", &s, &mtsq, kq2.clone());
    emit("ll", &s, &mtsq, function!(S.dot, k.clone(), k.clone()));
    emit("lq1q2", &s, &mtsq, &kq1 * &kq2);
}
