//! Symbolic lowering and generalized residues shared by FeynKit and GammaLoop.

use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    function,
    id::Replacement,
    symbol,
};

use crate::{CffError, CffExpression, EdgeId, EnergySurface, HSurface, Surface, SurfaceCache};

impl EnergySurface {
    /// Lower a surface using the caller's on-shell and external energy coordinates.
    pub fn to_atom_with(
        &self,
        energy: impl Fn(EdgeId) -> Atom,
        external_energy: impl Fn(EdgeId) -> Atom,
    ) -> Atom {
        self.energies
            .iter()
            .fold(Atom::Zero, |sum, edge| sum + energy(*edge))
            + self
                .external_shift
                .iter()
                .fold(Atom::Zero, |sum, (edge, sign)| {
                    sum + Atom::num(*sign) * external_energy(*edge)
                })
    }
}

impl HSurface {
    /// Lower positive/negative energies and the external shift without changing their signs.
    pub fn to_atom_with(
        &self,
        energy: impl Fn(EdgeId) -> Atom,
        external_energy: impl Fn(EdgeId) -> Atom,
    ) -> Atom {
        self.positive_energies
            .iter()
            .fold(Atom::Zero, |sum, edge| sum + energy(*edge))
            - self
                .negative_energies
                .iter()
                .fold(Atom::Zero, |sum, edge| sum + energy(*edge))
            + self
                .external_shift
                .iter()
                .fold(Atom::Zero, |sum, (edge, sign)| {
                    sum + Atom::num(*sign) * external_energy(*edge)
                })
    }
}

impl Surface {
    pub fn to_atom_with(
        &self,
        energy: impl Fn(EdgeId) -> Atom,
        external_energy: impl Fn(EdgeId) -> Atom,
    ) -> Atom {
        match self {
            Self::Energy(surface) => surface.to_atom_with(energy, external_energy),
            Self::H(surface) => surface.to_atom_with(energy, external_energy),
            Self::Unit => Atom::num(1),
            Self::Infinite => crate::SurfaceId::Infinite.into(),
        }
    }
}

impl SurfaceCache {
    /// Replace the canonical eta/H placeholders by their energy expressions.
    pub fn replacements_with(
        &self,
        energy: impl Fn(EdgeId) -> Atom,
        external_energy: impl Fn(EdgeId) -> Atom,
    ) -> Vec<Replacement> {
        self.energy_surfaces()
            .iter()
            .enumerate()
            .map(|(i, surface)| {
                Replacement::new(
                    Atom::from(crate::EnergySurfaceId(i)),
                    surface.to_atom_with(&energy, &external_energy),
                )
            })
            .chain(self.h_surfaces().iter().enumerate().map(|(i, surface)| {
                Replacement::new(
                    Atom::from(crate::HSurfaceId(i)),
                    surface.to_atom_with(&energy, &external_energy),
                )
            }))
            .collect()
    }
}

impl CffExpression {
    /// Bare denominator expression in canonical eta/H coordinates, without energy or measure factors.
    pub fn to_atom(&self) -> Atom {
        self.orientations
            .iter()
            .fold(Atom::Zero, |sum, orientation| {
                sum + orientation.expression.to_atom_inverse()
            })
    }

    /// GammaLoop's loop-energy integration phase and remaining spatial measure.
    pub fn measure_normalization(loop_number: usize) -> Atom {
        (-Atom::i()).pow(loop_number as i64)
            / (Atom::var(Symbol::PI) * 2).pow(3 * loop_number as i64)
    }

    /// Each uncontracted internal propagator contributes -1/(2 E).
    pub fn inverse_energy_product(energies: impl IntoIterator<Item = Atom>) -> Atom {
        energies
            .into_iter()
            .fold(Atom::num(1), |product, energy| product / (-2 * energy))
    }
}

/// A regular coefficient divided by a positive integer power of a simple surface zero.
///
/// The supplied root is assumed to be a simple zero of `surface`. Derivatives of
/// an implicit surface at that root are retained, exactly as in GammaLoop's
/// higher-order LU residues.
pub struct SurfacePole {
    pub surface: Atom,
    pub order: usize,
}

impl SurfacePole {
    /// Laurent coefficient of `coefficient / surface^order` about `variable=root`.
    /// Only negative Laurent powers down to the pole order are supported.
    pub fn laurent_coefficient(
        &self,
        coefficient: &Atom,
        variable: Symbol,
        root: &Atom,
        laurent_power: i64,
    ) -> Result<Atom, CffError> {
        if self.order == 0 || laurent_power >= 0 || laurent_power.unsigned_abs() > self.order as u64
        {
            return Err(CffError::Invariant(
                "a Laurent power must lie between minus the positive pole order and -1".into(),
            ));
        }
        if root.contains_symbol(variable) {
            return Err(CffError::Invariant(
                "a residue root must not depend on its integration variable".into(),
            ));
        }
        let order = i64::try_from(self.order)
            .map_err(|_| CffError::Invariant("pole order exceeds i64".into()))?;
        let root_value = self.surface.replace(variable).with(root.to_pattern());
        if matches!(root_value.as_view(), AtomView::Num(_)) && root_value != Atom::Zero {
            return Err(CffError::Invariant(
                "the supplied residue root is not a zero of the surface".into(),
            ));
        }
        if self
            .surface
            .derivative(variable)
            .replace(variable)
            .with(root.to_pattern())
            .expand()
            == Atom::Zero
        {
            return Err(CffError::Invariant(
                "the supplied residue root must be simple".into(),
            ));
        }
        let expansion = self
            .surface
            .series(variable, root.clone(), (order, 1))
            .map_err(|error| CffError::Invariant(format!("cannot expand pole surface: {error}")))?
            .to_atom()
            - root_value;
        let delta = symbol!("feynkit_cff::residue_delta");
        let mut regular =
            coefficient * expansion.pow(-order) * (Atom::var(variable) - root).pow(order);
        let derivative_order = (order + laurent_power) as usize;
        for _ in 0..derivative_order {
            regular = regular.derivative(variable);
        }
        // Taking the constant term after differentiation evaluates removable
        // singularities without an undefined direct substitution at the root.
        regular = regular.replace(Atom::var(variable) - root).with(delta);
        let mut result = regular
            .series(delta, Atom::Zero, (0, 1))
            .map_err(|error| CffError::Invariant(format!("cannot evaluate pole residue: {error}")))?
            .to_atom();
        for factor in 2..=derivative_order {
            result /= Atom::num(factor as i64);
        }
        Ok(result.replace(variable).with(root.to_pattern()))
    }

    pub fn residue(
        &self,
        coefficient: &Atom,
        variable: Symbol,
        root: &Atom,
    ) -> Result<Atom, CffError> {
        self.laurent_coefficient(coefficient, variable, root, -1)
    }
}

/// An oriented cut of `(energy^2 - on_shell_energy^2 + i prescription 0)^(-power)`.
///
/// `on_shell_energy` is the positive mass-shell energy, including the mass.
/// The generalized distribution Delta(n,x) acts as f^(n-1)(0)/(n-1)!,
/// not as an ordinary derivative of a Dirac delta. See arXiv:2203.11038,
/// equations (2.24) and (2.29).
#[derive(Clone, Debug)]
pub struct CutPropagator {
    pub energy: Atom,
    pub on_shell_energy: Atom,
    pub power: usize,
    pub orientation: i8,
    pub prescription: i8,
    /// Overall replacement factor; -2 pi i gives reciprocal-propagator cuts.
    pub normalization: Atom,
}

impl CutPropagator {
    pub fn validate(&self) -> Result<(), CffError> {
        if self.power == 0 || self.power > i64::MAX as usize {
            return Err(CffError::Invariant(
                "cut propagator power must be a positive integer representable as i64".into(),
            ));
        }
        if self.on_shell_energy == Atom::Zero {
            return Err(CffError::Invariant(
                "cut on-shell energy must be nonzero so the energy roots are distinct".into(),
            ));
        }
        if ![-1, 1].contains(&self.orientation) || ![-1, 1].contains(&self.prescription) {
            return Err(CffError::Invariant(
                "cut orientation and prescription must be +1 or -1".into(),
            ));
        }
        Ok(())
    }

    /// Covariant generalized mass-shell distribution with explicit energy orientation.
    pub fn to_atom(&self) -> Result<Atom, CffError> {
        self.validate()?;
        Ok(&self.normalization
            * Atom::num(self.prescription as i64)
            * function!(
                symbol!("CutDelta"),
                self.power as i64,
                self.energy.pow(2) - self.on_shell_energy.pow(2),
                &self.energy,
                self.orientation as i64
            ))
    }

    /// Energy-space form; the uncut factor must remain inside the derivative action.
    pub fn to_energy_atom(&self) -> Result<Atom, CffError> {
        self.validate()?;
        let root = Atom::num(self.orientation as i64) * &self.on_shell_energy;
        Ok(&self.normalization
            * Atom::num((self.prescription * self.orientation) as i64)
            * function!(symbol!("Delta"), self.power as i64, &self.energy - &root)
            / (&self.energy + root).pow(self.power as i64))
    }

    /// Act on the complete remaining coefficient before setting the energy on shell.
    /// The integration variable must equal the propagator's energy; explicitly
    /// route to that coordinate before applying sequential cut distributions.
    pub fn apply(&self, coefficient: &Atom, variable: Symbol) -> Result<Atom, CffError> {
        self.validate()?;
        if self.energy != Atom::var(variable) || self.on_shell_energy.contains_symbol(variable) {
            return Err(CffError::Invariant("cut action requires an independent energy variable and an on-shell energy independent of it".into()));
        }
        let root = Atom::num(self.orientation as i64) * &self.on_shell_energy;
        let mut regular = coefficient / (&self.energy + &root).pow(self.power as i64);
        for order in 1..self.power {
            regular = regular.derivative(variable) / Atom::num(order as i64);
        }
        Ok(&self.normalization
            * Atom::num((self.prescription * self.orientation) as i64)
            * regular.replace(variable).with(root.to_pattern()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::parse;

    #[test]
    fn residue_rejects_invalid_orders_and_non_simple_roots() {
        let x = symbol!("cut_test::invalid_x");
        let pole = SurfacePole {
            surface: Atom::var(x).pow(2),
            order: 1,
        };
        assert!(pole.residue(&Atom::num(1), x, &Atom::Zero).is_err());
        assert!(pole.residue(&Atom::num(1), x, &Atom::num(1)).is_err());
        assert!(
            pole.laurent_coefficient(&Atom::num(1), x, &Atom::Zero, 0)
                .is_err()
        );
    }

    #[test]
    fn energy_lowering_and_normalization_preserve_runtime_signs() {
        let h = HSurface {
            positive_energies: vec![EdgeId::new(0)],
            negative_energies: vec![EdgeId::new(1)],
            external_shift: vec![(EdgeId::new(2), -1)].into(),
            vertex_set: crate::VertexSet::default(),
        };
        let e = |edge: EdgeId| function!(symbol!("cut_test::E"), edge.index() as i64);
        let p = |edge: EdgeId| function!(symbol!("cut_test::P"), edge.index() as i64);
        assert_eq!(
            h.to_atom_with(e, p),
            e(EdgeId::new(0)) - e(EdgeId::new(1)) - p(EdgeId::new(2))
        );
        assert_eq!(
            CffExpression::inverse_energy_product([Atom::num(2), Atom::num(3)]),
            Atom::num(1) / 24
        );
        assert_eq!(CffExpression::inverse_energy_product([]), Atom::num(1));
        assert_eq!(CffExpression::measure_normalization(0), Atom::num(1));
        assert_eq!(
            CffExpression::measure_normalization(1),
            -Atom::i() / (Atom::num(2) * Atom::var(Symbol::PI)).pow(3)
        );
    }

    #[test]
    fn raised_cut_action_equals_energy_residue_for_both_orientations() {
        let x = symbol!("cut_test::x");
        let energy = parse!("sqrt(cut_test::k^2+cut_test::m^2)");
        let coefficient = (Atom::var(x) + 3).pow(4);
        for power in 1..=3 {
            for orientation in [-1, 1] {
                let cut = CutPropagator {
                    energy: Atom::var(x),
                    on_shell_energy: energy.clone(),
                    power,
                    orientation,
                    prescription: 1,
                    normalization: -2 * Atom::var(Symbol::PI) * Atom::i(),
                };
                let pole = SurfacePole {
                    surface: Atom::var(x).pow(2) - energy.pow(2),
                    order: power,
                };
                let root = Atom::num(orientation as i64) * &energy;
                let residue = pole.residue(&coefficient, x, &root).unwrap();
                let expected = &cut.normalization * Atom::num(orientation as i64) * residue;
                assert_eq!(
                    (cut.apply(&coefficient, x).unwrap() - expected)
                        .together()
                        .expand(),
                    Atom::Zero
                );
            }
        }
    }

    #[test]
    fn triple_cut_differentiates_the_uncut_factor_and_coefficient() {
        let x = symbol!("cut_test::q0");
        let cut = CutPropagator {
            energy: Atom::var(x),
            on_shell_energy: Atom::num(2),
            power: 3,
            orientation: 1,
            prescription: 1,
            normalization: Atom::num(1),
        };
        assert_eq!(
            cut.apply(&Atom::var(x).pow(2), x).unwrap(),
            Atom::num(-1) / 128
        );
        let mut opposite = cut.clone();
        opposite.prescription = -1;
        assert_eq!(
            opposite.apply(&Atom::var(x).pow(2), x).unwrap(),
            Atom::num(1) / 128
        );
    }

    #[test]
    fn implicit_surface_laurent_coefficients_keep_surface_derivatives() {
        let x = symbol!("cut_test::t");
        let root = Atom::var(symbol!("cut_test::root"));
        let surface = function!(symbol!("cut_test::eta"), x);
        let coefficient = function!(symbol!("cut_test::f"), x);
        let pole = SurfacePole {
            surface: surface.clone(),
            order: 2,
        };
        let derivative = surface.derivative(x).replace(x).with(root.to_pattern());
        let expected = coefficient.derivative(x).replace(x).with(root.to_pattern())
            / derivative.pow(2)
            - coefficient.replace(x).with(root.to_pattern())
                * surface
                    .derivative(x)
                    .derivative(x)
                    .replace(x)
                    .with(root.to_pattern())
                / derivative.pow(3);
        assert_eq!(
            (pole.residue(&coefficient, x, &root).unwrap() - expected)
                .together()
                .expand(),
            Atom::Zero
        );
    }
}
