//! Blast-radius probe for the off-shell-regularization fix: does reduce() accept a SYMBOLIC
//! invariant delta (so the on-shell leg -> delta, reduce, coefficients rational in delta, limit
//! downstream)?  If yes, the fix needs NO reducer change. Gitignored, local-only.
//!   cargo run --release --example symbolic_delta -p oneloop

use oneloop::symbols::S;
use oneloop::{Integral, IntegralFamily, Kinematics, Propagator, reduce};
use symbolica::atom::Atom;
use symbolica::{function, symbol};

fn main() {
    gammalooprs::initialisation::initialise().expect("activate Symbolica license");
    let delta = Atom::var(symbol!("oneloop::delta"));
    let m2 = Atom::num(1); // massive internal
    // gg>h triangle, the two on-shell legs kept SYMBOLIC = delta, third leg = 2/5.
    let fam = IntegralFamily {
        propagators: vec![
            Propagator {
                momentum: Atom::Zero,
                mass_sq: m2.clone(),
            },
            Propagator {
                momentum: Atom::Zero,
                mass_sq: m2.clone(),
            },
            Propagator {
                momentum: Atom::Zero,
                mass_sq: m2.clone(),
            },
        ],
        isps: vec![],
        kinematics: Kinematics {
            invariants: vec![delta.clone(), delta.clone(), Atom::num(2) / Atom::num(5)],
        },
        targets: vec![Integral {
            propagator_exponents: vec![1, 1, 1],
            isp_exponents: vec![],
        }],
        numerator: {
            let kq1 = function!(S.dot, Atom::var(S.k), Atom::var(S.q1));
            &kq1 * &kq1
        },
    };
    let res = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| reduce(&fam)));
    match res {
        Err(_) => println!(
            "PANIC — reduce() does NOT accept symbolic invariants (fix needs reducer change)"
        ),
        Ok(r) => {
            println!(
                "OK — reduce() ACCEPTS symbolic delta. {} terms:",
                r.terms.len()
            );
            for (c, m) in r.terms.iter().take(6) {
                let cs = c.to_string();
                let has_delta = cs.contains("delta");
                let short = if cs.len() > 70 {
                    format!("{}…", &cs[..70])
                } else {
                    cs
                };
                println!(
                    "  coeff{}: {short}   -> {m:?}",
                    if has_delta { " [has delta]" } else { "" }
                );
            }
        }
    }
}
