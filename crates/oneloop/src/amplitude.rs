use symbolica::atom::Atom;

use crate::family::IntegralFamily;
use crate::masters::{MasterBasis, OneLoopMasters};
use crate::reduce::reduce;

pub fn amplitude(family: &IntegralFamily) -> Atom {
    let basis = OneLoopMasters;
    reduce(family)
        .terms
        .iter()
        .fold(Atom::Zero, |acc, (coeff, master)| {
            acc + coeff * basis.symbol(master)
        })
}

#[cfg(test)]
mod tests {
    use super::amplitude;
    use crate::family::{Integral, IntegralFamily, Kinematics, Propagator};
    use crate::symbols::S;
    use symbolica::atom::Atom;
    use symbolica::function;

    fn massless_bubble_family() -> IntegralFamily {
        IntegralFamily {
            propagators: vec![
                Propagator {
                    momentum: Atom::Zero,
                    mass_sq: Atom::Zero,
                },
                Propagator {
                    momentum: Atom::Zero,
                    mass_sq: Atom::Zero,
                },
            ],
            isps: vec![],
            kinematics: Kinematics {
                invariants: vec![Atom::var(S.psq)],
                masses_sq: vec![],
            },
            targets: vec![Integral {
                propagator_exponents: vec![1, 1],
                isp_exponents: vec![],
            }],
        }
    }

    #[test]
    fn massless_bubble_amplitude_is_symbolic_b0() {
        crate::ensure_symbolica_license();
        let amp = amplitude(&massless_bubble_family());
        let want = function!(S.b0, Atom::var(S.psq), Atom::Zero, Atom::Zero);
        assert_eq!(amp, want);
    }
}
