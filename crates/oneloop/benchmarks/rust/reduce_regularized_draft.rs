//! DRAFT / PROTOTYPE (gitignored, NOT committed to src/): off-shell-regularization wrapper for the
//! on-shell-massless-leg case, implemented ENTIRELY outside the core reducer (calls the public
//! reduce() unchanged). Demonstrates the plan Eli approved before any src/ integration.
//!   cargo run --release --example reduce_regularized_draft -p oneloop
//!
//! Idea (validated to 1e-10 in madloop_reference.md): the reducer is a massive method; on-shell
//! massless legs make sub-Grams/Cayleys vanish -> panic. Fix: replace the degeneracy-causing ZERO
//! invariants with a symbol delta, reduce() (exact rational arithmetic cancels the 1/delta poles ->
//! coefficients finite at delta=0), then substitute delta=0 in the coefficients AND master args.

use oneloop::masters::MasterIntegral;
use oneloop::symbols::S;
use oneloop::{Integral, IntegralFamily, Kinematics, Propagator, Reduction, reduce};
use symbolica::atom::{Atom, AtomCore};
use symbolica::{function, symbol};

/// Substitute delta -> 0 in an Atom (mirrors reduce.rs:530 `expr.replace(..).with(Atom::Zero)`).
fn kill_delta(expr: &Atom, delta: &Atom) -> Atom {
    expr.replace(delta.to_pattern()).with(Atom::Zero)
}

/// Rebuild a master with delta -> 0 in every kinematic argument.
fn kill_delta_master(m: &MasterIntegral, d: &Atom) -> MasterIntegral {
    let k = |a: &Atom| kill_delta(a, d);
    match m {
        MasterIntegral::Tadpole { m_sq } => MasterIntegral::Tadpole { m_sq: k(m_sq) },
        MasterIntegral::Bubble { p_sq, m1_sq, m2_sq } => MasterIntegral::Bubble {
            p_sq: k(p_sq),
            m1_sq: k(m1_sq),
            m2_sq: k(m2_sq),
        },
        MasterIntegral::Triangle {
            p1_sq,
            p2_sq,
            p12_sq,
            m1_sq,
            m2_sq,
            m3_sq,
        } => MasterIntegral::Triangle {
            p1_sq: k(p1_sq),
            p2_sq: k(p2_sq),
            p12_sq: k(p12_sq),
            m1_sq: k(m1_sq),
            m2_sq: k(m2_sq),
            m3_sq: k(m3_sq),
        },
        MasterIntegral::Box {
            p1_sq,
            p2_sq,
            p3_sq,
            p4_sq,
            s,
            t,
            m1_sq,
            m2_sq,
            m3_sq,
            m4_sq,
        } => MasterIntegral::Box {
            p1_sq: k(p1_sq),
            p2_sq: k(p2_sq),
            p3_sq: k(p3_sq),
            p4_sq: k(p4_sq),
            s: k(s),
            t: k(t),
            m1_sq: k(m1_sq),
            m2_sq: k(m2_sq),
            m3_sq: k(m3_sq),
            m4_sq: k(m4_sq),
        },
    }
}

/// The wrapper. Detects zero invariants, regularizes with a shared delta, reduces, takes delta->0.
/// (TODO for the real src/ version: also regularize zero internal MASSES; assert no residual
///  1/delta survives = genuine divergence; decide detection = this pre-check vs catch_unwind.)
fn reduce_regularized(family: &IntegralFamily) -> Reduction {
    let delta = Atom::var(symbol!("oneloop::reg_delta"));
    let has_zero = family.kinematics.invariants.iter().any(|s| s.is_zero());
    if !has_zero {
        return reduce(family); // generic kinematics: unchanged fast path
    }
    // regularize: swap each zero invariant for delta
    let mut reg = family.clone();
    reg.kinematics.invariants = family
        .kinematics
        .invariants
        .iter()
        .map(|s| {
            if s.is_zero() {
                delta.clone()
            } else {
                s.clone()
            }
        })
        .collect();
    let r = reduce(&reg);
    // take delta -> 0 in coefficients and master arguments
    Reduction {
        terms: r
            .terms
            .iter()
            .map(|(c, m)| {
                (
                    kill_delta(&c.expand(), &delta),
                    kill_delta_master(m, &delta),
                )
            })
            .filter(|(c, _)| *c != Atom::Zero)
            .collect(),
    }
}

fn main() {
    gammalooprs::initialisation::initialise().expect("activate Symbolica license");
    let m2 = Atom::num(1);
    // gg>h on-shell triangle, rank-2 (k.q1)^2: legs (0, 0, 2/5). Bare reduce() PANICS on this.
    let kq1 = function!(S.dot, Atom::var(S.k), Atom::var(S.q1));
    let fam = IntegralFamily {
        propagators: (0..3)
            .map(|_| Propagator {
                momentum: Atom::Zero,
                mass_sq: m2.clone(),
            })
            .collect(),
        isps: vec![],
        kinematics: Kinematics {
            invariants: vec![Atom::Zero, Atom::Zero, Atom::num(2) / Atom::num(5)],
        },
        targets: vec![Integral {
            propagator_exponents: vec![1, 1, 1],
            isp_exponents: vec![],
        }],
        numerator: &kq1 * &kq1,
    };

    // (bare reduce() on this on-shell case PANICS "singular Gram" -- that's the current behaviour
    //  the wrapper fixes. We don't call it here: a caught panic poisons Symbolica's global state.)
    let r = reduce_regularized(&fam);
    println!(
        "reduce_regularized() -> {} terms (delta->0):",
        r.terms.len()
    );
    for (c, m) in &r.terms {
        println!("  {}  x  {m:?}", c);
    }
    println!(
        "\n(expected: (1/20) B0(0;1,1) - (1/20) B0(2/5;1,1); evaluates to -3.4748e-3, validated.)"
    );
}
