//! Emit the gg->h massive-top-triangle dot(k,.) monomial reductions for the form-factor benchmark (see docs/09).

use oneloop::masters::MasterIntegral;
use oneloop::symbols::S;
use oneloop::{Integral, IntegralFamily, Kinematics, Propagator, reduce};
use symbolica::atom::{Atom, AtomCore};
use symbolica::function;

// gg->h top triangle: 3 top propagators (mass^2=mtsq), invariants s_01=q1^2=0, s_02=(q1+q2)^2=s, s_12=q2^2=0.
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
