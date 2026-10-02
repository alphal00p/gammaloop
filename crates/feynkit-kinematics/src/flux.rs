//! Initial-state normalization shared by symbolic and numerical calculations.

use std::ops::{Add, Mul, Sub};

/// The denominator multiplying an invariant squared matrix element.
///
/// Decays use `2 E` in the frame where the rate is measured. Scattering uses
/// `4 sqrt((p1.p2)^2 - m1^2 m2^2)`. Spin/color averages, final-state symmetry
/// factors and unit conversions belong to the caller.
#[derive(Clone, Debug)]
pub enum InitialStateFlux<T> {
    Decay {
        energy: T,
    },
    Scattering {
        momentum_dot: T,
        mass_squared_product: T,
    },
}

impl<T> InitialStateFlux<T>
where
    T: Clone + Add<Output = T> + Mul<Output = T> + Sub<Output = T>,
{
    /// Evaluate using the scalar's square-root operation.
    ///
    /// Supplying this operation keeps symbolic branch assumptions with
    /// Symbolica and numerical precision with the caller's scalar type.
    /// Inputs must describe physical incoming states; no threshold clipping
    /// or absolute-value replacement is performed.
    pub fn denominator(self, sqrt: impl FnOnce(T) -> T) -> T {
        match self {
            Self::Decay { energy } => energy.clone() + energy,
            Self::Scattering {
                momentum_dot,
                mass_squared_product,
            } => {
                // The Kallen form removes the rational 1/4 inside the root
                // for massless symbolic kinematics, leaving sqrt(s^2).
                let twice_dot = momentum_dot.clone() + momentum_dot;
                let twice_mass_product = mass_squared_product.clone() + mass_squared_product;
                let root = sqrt(
                    twice_dot.clone() * twice_dot
                        - (twice_mass_product.clone() + twice_mass_product),
                );
                root.clone() + root
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::InitialStateFlux;
    use crate::FourMomentum;
    use symbolica::{
        atom::{Atom, AtomCore},
        parse,
    };

    #[test]
    fn massive_massless_and_decay_flux_agree_across_scalar_types() {
        let first = FourMomentum::from_args(5.0, 0.0, 0.0, 4.0);
        let second = FourMomentum::from_args(13.0, 0.0, 0.0, -12.0);
        let flux = InitialStateFlux::Scattering {
            momentum_dot: first.dot(&second),
            mass_squared_product: first.mass_squared() * second.mass_squared(),
        }
        .denominator(f64::sqrt);
        assert_eq!(flux, 448.0);
        assert_eq!(first.flux(Some(&second)), flux);
        assert_eq!(first.flux(None), 10.0);
        let symbolic = InitialStateFlux::Scattering {
            momentum_dot: Atom::num(113),
            mass_squared_product: Atom::num(225),
        }
        .denominator(|x| x.sqrt());
        assert_eq!(symbolic, Atom::num(448));
        assert_eq!(
            InitialStateFlux::Scattering {
                momentum_dot: 50.0,
                mass_squared_product: 0.0
            }
            .denominator(f64::sqrt),
            200.0
        );
        assert_eq!(
            InitialStateFlux::Decay {
                energy: parse!("E")
            }
            .denominator(|x| x.sqrt()),
            parse!("2*E")
        );
        assert_eq!(
            InitialStateFlux::Decay { energy: 13.0 }.denominator(f64::sqrt),
            26.0
        );
    }
}
