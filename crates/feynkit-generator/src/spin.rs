//! External-state completeness relations shared by FeynKit and GammaLoop.

use feynkit_graph::symbols;
use feynkit_model::{Model, Particle};
use idenso::{dirac::AGS, representations::Bispinor};
use spenso::network::library::symbolic::ETS;
use spenso::structure::{
    abstract_index::AbstractIndex,
    representation::{Minkowski, RepName},
    slot::DummyAind,
};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
    id::Replacement,
};
use thiserror::Error;

/// Invalid or unsupported external-state polarization sum.
#[derive(Debug, Error)]
pub enum SpinSumError {
    #[error("polarization sums for UFO spin {0} are not implemented")]
    UnsupportedSpin(i64),
    #[error("expected UFO spin {expected}, got {actual}")]
    WrongSpin { expected: i64, actual: i64 },
    #[error("the axial reference must have nonzero momentum contraction")]
    OrthogonalReference,
    #[error("a momentum must be an unindexed symbol or function call")]
    InvalidMomentum,
    #[error("a spin vector requires a massive Dirac fermion")]
    SpinVectorRequiresMassiveFermion,
    #[error("a selected spin state cannot also be spin averaged")]
    PolarizedAverage,
}

/// An axial-gauge reference contracted with the two open Lorentz slots.
///
/// The reference need not be null. All contractions use the same Lorentz
/// metric as the momentum and open slots.
pub struct AxialReference {
    pub left: Atom,
    pub right: Atom,
    pub momentum_dot_reference: Atom,
    pub norm_squared: Atom,
}

/// Completeness relations for scalars, Dirac fermions, and vector particles.
///
/// Inputs and outputs are ordinary Symbolica expressions in Spenso notation.
/// Kinematics, on-shell assumptions, and the choice of axial reference stay
/// with the caller. A covariant massless vector sum assumes a gauge-invariant
/// amplitude; a physical sum requires an explicit reference.
pub struct SpinSum {
    spin: i64,
    mass: Atom,
    massless: bool,
    antiparticle: bool,
    average: bool,
    covariant: bool,
}

impl SpinSum {
    pub fn new(particle: &Particle, model: &Model) -> Result<Self, SpinSumError> {
        if !(1..=3).contains(&particle.spin) {
            return Err(SpinSumError::UnsupportedSpin(particle.spin));
        }
        let massless = particle.is_massless(model);
        Ok(Self {
            spin: particle.spin,
            mass: particle
                .symbolic_mass(model)
                .replace(symbolica::symbol!("UFO::ZERO"))
                .with(Atom::Zero),
            massless,
            antiparticle: particle.is_antiparticle(),
            average: false,
            covariant: false,
        })
    }

    /// Divide by the number of physical spin states in four dimensions.
    pub fn averaged(mut self, average: bool) -> Self {
        self.average = average;
        self
    }

    /// Use `-g` for vector sums, including massive vectors in Feynman gauge.
    pub fn covariant(mut self, covariant: bool) -> Self {
        self.covariant = covariant;
        self
    }

    /// Build a four-dimensional completeness tensor from unindexed momenta.
    ///
    /// The two indices are bare Symbolica indices, not Spenso slots. Fermions
    /// use `bis(4, index)` and vectors use `mink(4, index)`. A momentum may be
    /// a symbol (`p`) or a labeled function (`Q(1)`). For massless vectors,
    /// omitting the reference selects the covariant sum. A massive Dirac spin
    /// vector selects one physical state, with `p.s = 0` and `s.s = -1` imposed
    /// by the caller. It uses `(slash(p) +/- m)(1 + gamma5 slash(s))/2` for
    /// particles and antiparticles alike and cannot be combined with averaging.
    pub fn expression(
        &self,
        momentum: &Atom,
        indices: [Atom; 2],
        reference: Option<&Atom>,
        spin_vector: Option<&Atom>,
    ) -> Result<Atom, SpinSumError> {
        use symbolica::atom::AtomView;
        if !matches!(momentum.as_view(), AtomView::Var(_) | AtomView::Fun(_))
            || [reference, spin_vector]
                .into_iter()
                .flatten()
                .any(|vector| !matches!(vector.as_view(), AtomView::Var(_) | AtomView::Fun(_)))
        {
            return Err(SpinSumError::InvalidMomentum);
        }
        if spin_vector.is_some() {
            if self.spin != 2 || self.massless {
                return Err(SpinSumError::SpinVectorRequiresMassiveFermion);
            }
            if self.average {
                return Err(SpinSumError::PolarizedAverage);
            }
        }
        let lorentz = Minkowski {}.new_rep(4);
        match self.spin {
            1 => Ok(Atom::one()),
            2 => {
                let spinor = Bispinor {}.new_rep(4);
                let [left, right] = indices.map(|index| spinor.pattern(index));
                let end = if spin_vector.is_some() {
                    spinor.pattern(AbstractIndex::new_dummy().to_atom())
                } else {
                    right.clone()
                };
                let slash = function!(
                    AGS.gamma,
                    &left,
                    &end,
                    lorentz.vector(momentum.as_view(), [])
                );
                let completeness = self.fermion(slash, &left, &end)?;
                Ok(if let Some(spin) = spin_vector {
                    let middle = spinor.pattern(AbstractIndex::new_dummy().to_atom());
                    let projection = function!(ETS.metric, &end, &right)
                        + function!(AGS.gamma5, &end, &middle)
                            * function!(
                                AGS.gamma,
                                &middle,
                                &right,
                                lorentz.vector(spin.as_view(), [])
                            );
                    completeness * projection / Atom::num(2)
                } else {
                    completeness
                })
            }
            3 => {
                let components = indices
                    .clone()
                    .map(|index| lorentz.vector(momentum.as_view(), [index]));
                let reference = reference.map(|reference| AxialReference {
                    left: lorentz.vector(reference.as_view(), [indices[0].clone()]),
                    right: lorentz.vector(reference.as_view(), [indices[1].clone()]),
                    momentum_dot_reference: lorentz.inner_product(momentum, reference),
                    norm_squared: lorentz.inner_product(reference, reference),
                });
                let [left, right] = indices.map(|index| lorentz.pattern(index));
                self.vector(
                    [&components[0], &components[1]],
                    &left,
                    &right,
                    reference.as_ref(),
                )
            }
            spin => Err(SpinSumError::UnsupportedSpin(spin)),
        }
    }

    pub fn averaging_factor(&self) -> Atom {
        let states = match (self.average, self.spin, self.massless) {
            (false, _, _) | (_, 1, _) => 1,
            (_, 2, _) | (_, 3, true) => 2,
            _ => 3,
        };
        Atom::one() / Atom::num(states)
    }

    /// Sum this particle's paired wavefunctions on one labeled edge.
    ///
    /// Unpaired wavefunctions and other edge labels are left unchanged. The
    /// momentum and reference conventions are those of [`Self::expression`].
    pub fn apply(
        &self,
        expression: &Atom,
        momentum: &Atom,
        edge: usize,
        reference: Option<&Atom>,
        spin_vector: Option<&Atom>,
    ) -> Result<Atom, SpinSumError> {
        let indices = [
            symbolica::symbol!("feynkit::spin_left_").into(),
            symbolica::symbol!("feynkit::spin_right_").into(),
        ];
        let tensor = self.expression(momentum, indices.clone(), reference, spin_vector)?;
        let slots = indices.map(|index| match self.spin {
            2 => Bispinor {}.new_rep(4).pattern(index),
            _ => Minkowski {}.new_rep(4).pattern(index),
        });
        Ok(match self.replacement(edge, slots, tensor) {
            Some(replacement) => expression.replace_multiple(&[replacement]),
            None => expression.clone(),
        })
    }

    /// Replace a paired external wavefunction with its completeness tensor.
    ///
    /// `indices` are full Spenso slots (possibly patterns) appearing in `tensor`.
    /// Keeping the tensor supplied by the caller lets symbolic calculations
    /// and GammaLoop's runtime gauge reference share the same wavefunction rule.
    /// Scalars have no wavefunctions and require no replacement.
    pub fn replacement(
        &self,
        edge: usize,
        indices: [Atom; 2],
        tensor: Atom,
    ) -> Option<Replacement> {
        let (ket, bra) = match (self.spin, self.antiparticle) {
            (1, _) => return None,
            (2, false) => (symbols::u(), symbols::ubar()),
            (2, true) => (symbols::v(), symbols::vbar()),
            (3, _) => (symbols::epsilon(), symbols::epsilonbar()),
            _ => unreachable!("SpinSum::new validates the spin"),
        };
        let pair = function!(ket, edge, &indices[0]) * function!(bra, edge, &indices[1]);
        Some(Replacement::new(pair.to_pattern(), tensor))
    }

    /// Return `(slash(p) +/- m identity)` with optional spin averaging.
    ///
    /// `slashed_momentum` must have the supplied open spinor slots. Passing
    /// an already contracted slash lets callers use their own dummy indices.
    pub fn fermion(
        &self,
        slashed_momentum: Atom,
        left: &Atom,
        right: &Atom,
    ) -> Result<Atom, SpinSumError> {
        if self.spin != 2 {
            return Err(SpinSumError::WrongSpin {
                expected: 2,
                actual: self.spin,
            });
        }
        let sign = if self.antiparticle { -1 } else { 1 };
        Ok(
            (slashed_momentum + Atom::num(sign) * &self.mass * function!(ETS.metric, left, right))
                * self.averaging_factor(),
        )
    }

    /// Return the covariant massless, massive Proca, or axial projector.
    ///
    /// `momentum` contains the momentum components in `left` and `right`.
    /// A supplied axial reference is used only for a massless particle.
    pub fn vector(
        &self,
        momentum: [&Atom; 2],
        left: &Atom,
        right: &Atom,
        reference: Option<&AxialReference>,
    ) -> Result<Atom, SpinSumError> {
        if self.spin != 3 {
            return Err(SpinSumError::WrongSpin {
                expected: 3,
                actual: self.spin,
            });
        }
        let mut projector = -function!(ETS.metric, left, right);
        if !self.covariant {
            if !self.massless {
                projector += momentum[0] * momentum[1] / self.mass.clone().pow(2);
            } else if let Some(reference) = reference {
                if reference.momentum_dot_reference.is_zero() {
                    return Err(SpinSumError::OrthogonalReference);
                }
                projector += (momentum[0] * &reference.right + &reference.left * momentum[1])
                    / &reference.momentum_dot_reference
                    - &reference.norm_squared * momentum[0] * momentum[1]
                        / reference.momentum_dot_reference.clone().pow(2);
            }
        }
        Ok(projector * self.averaging_factor())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use idenso::{dirac::GammaSimplifier, shorthands::schoonschip::Schoonschip};
    use symbolica::{parse, symbol};

    #[test]
    fn wavefunction_pairs_preserve_species_edges_and_unpaired_states() {
        let model =
            Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap();
        let p = parse!("p");
        let indices = [parse!("i"), parse!("j")];
        for (pdg, ket, bra) in [
            (11, symbols::u(), symbols::ubar()),
            (-11, symbols::v(), symbols::vbar()),
            (22, symbols::epsilon(), symbols::epsilonbar()),
            (23, symbols::epsilon(), symbols::epsilonbar()),
        ] {
            let particle = model
                .particle_by_id(model.particle_id_by_pdg(pdg).unwrap())
                .unwrap();
            let sum = SpinSum::new(particle, &model).unwrap().averaged(true);
            let slots = indices.clone().map(|index| match particle.spin {
                2 => Bispinor {}.new_rep(4).pattern(index),
                _ => Minkowski {}.new_rep(4).pattern(index),
            });
            let pair = function!(ket, 7, &slots[0]) * function!(bra, 7, &slots[1]);
            let spectator = function!(ket, 8, &slots[0]) * function!(bra, 8, &slots[1]);
            let unpaired = function!(ket, 7, &slots[0]);
            let input = &pair * &spectator + &unpaired;
            let tensor = sum.expression(&p, indices.clone(), None, None).unwrap();
            let output = sum.apply(&input, &p, 7, None, None).unwrap();
            assert!(
                (output.clone() - tensor * spectator - unpaired)
                    .expand()
                    .is_zero()
            );
            assert_eq!(sum.apply(&output, &p, 7, None, None).unwrap(), output);
            let wrong =
                function!(symbols::v(), 7, &slots[0]) * function!(symbols::ubar(), 7, &slots[1]);
            assert_eq!(sum.apply(&wrong, &p, 7, None, None).unwrap(), wrong);
        }
        let scalar = model
            .particle_by_id(model.particle_id_by_pdg(25).unwrap())
            .unwrap();
        let input = parse!("a+b");
        assert_eq!(
            SpinSum::new(scalar, &model)
                .unwrap()
                .apply(&input, &p, 7, None, None)
                .unwrap(),
            input
        );
    }

    #[test]
    fn fermion_and_antifermion_have_opposite_mass_terms() {
        let model =
            Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap();
        let p = parse!("p");
        let indices = [parse!("i"), parse!("j")];
        let particle = model
            .particle_by_id(model.particle_id_by_pdg(13).unwrap())
            .unwrap();
        let antiparticle = model.particle_by_id(particle.antiparticle).unwrap();
        let u = SpinSum::new(particle, &model)
            .unwrap()
            .expression(&p, indices.clone(), None, None)
            .unwrap();
        let v = SpinSum::new(antiparticle, &model)
            .unwrap()
            .expression(&p, indices.clone(), None, None)
            .unwrap();
        let identity = Bispinor {}.new_rep(4).id(&indices[0], &indices[1]);
        assert!(
            (u.clone() - v - Atom::num(2) * particle.symbolic_mass(&model) * identity)
                .expand()
                .is_zero()
        );
        let averaged = SpinSum::new(particle, &model)
            .unwrap()
            .averaged(true)
            .expression(&p, indices, None, None)
            .unwrap();
        assert!((Atom::num(2) * averaged - u).expand().is_zero());
    }

    #[test]
    fn vector_projectors_are_transverse_and_count_physical_states() {
        let model =
            Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap();
        let p = symbol!("spin_test::p").to_atom();
        let n = symbol!("spin_test::n").to_atom();
        let mu = parse!("mu");
        let nu = parse!("nu");
        let rep = Minkowski {}.new_rep(4);
        let p2 = rep.inner_product(&p, &p);
        for (pdg, states, reference) in [(22, 2, Some(&n)), (23, 3, None)] {
            let particle = model
                .particle_by_id(model.particle_id_by_pdg(pdg).unwrap())
                .unwrap();
            let projector = SpinSum::new(particle, &model)
                .unwrap()
                .expression(&p, [mu.clone(), nu.clone()], reference, None)
                .unwrap();
            let mass_squared = if pdg == 22 {
                Atom::Zero
            } else {
                particle.symbolic_mass(&model).pow(2)
            };
            let longitudinal = (projector.clone() * rep.vector(p.as_view(), [mu.clone()]))
                .expand()
                .to_dots()
                .replace(p2.to_pattern())
                .with(mass_squared.to_pattern())
                .together();
            assert!(longitudinal.is_zero(), "{longitudinal}");
            if let Some(reference) = reference {
                let axial = (projector.clone() * rep.vector(reference.as_view(), [mu.clone()]))
                    .expand()
                    .to_dots()
                    .together();
                assert!(axial.is_zero(), "{axial}");
            }
            let trace = (projector * rep.id(&mu, &nu))
                .expand()
                .to_dots()
                .replace(p2.to_pattern())
                .with(mass_squared.to_pattern())
                .together();
            assert_eq!(trace, Atom::num(-states));
        }
    }

    #[test]
    fn spinor_trace_and_scalar_sum() {
        let model =
            Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap();
        for (pdg, expected) in [(13, 4), (-13, -4)] {
            let particle = model
                .particle_by_id(model.particle_id_by_pdg(pdg).unwrap())
                .unwrap();
            let projector = SpinSum::new(particle, &model)
                .unwrap()
                .expression(&parse!("p"), [parse!("i"), parse!("i")], None, None)
                .unwrap();
            assert_eq!(
                projector.simplify_gamma().expand(),
                expected * particle.symbolic_mass(&model)
            );
        }
        let scalar = model
            .particle_by_id(model.particle_id_by_pdg(25).unwrap())
            .unwrap();
        assert_eq!(
            SpinSum::new(scalar, &model)
                .unwrap()
                .averaged(true)
                .expression(&parse!("p"), [parse!("i"), parse!("j")], None, None)
                .unwrap(),
            Atom::one()
        );
    }
}
