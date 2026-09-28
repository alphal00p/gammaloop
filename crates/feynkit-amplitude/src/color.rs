//! External color completeness using the same representations as generated rules.

use feynkit_model::Particle;
use idenso::representations::{ColorAdjoint, ColorFundamental, ColorSextet};
use spenso::{
    network::library::symbolic::ETS,
    structure::representation::{LibraryRep, RepName, Representation},
};
use symbolica::atom::Atom;
use thiserror::Error;

/// Non-singlet UFO color spaces shared by generation, completeness and runtime indices.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColorRepresentation {
    Fundamental,
    AntiFundamental,
    Sextet,
    AntiSextet,
    Adjoint,
}

impl ColorRepresentation {
    /// Resolve a supported signed UFO color code.
    /// Singlets and unsupported codes return `None`; [`ColorSum::new`] validates
    /// the distinction when constructing a completeness tensor.
    pub fn from_ufo(color: i64) -> Option<Self> {
        match color {
            3 => Some(ColorRepresentation::Fundamental),
            -3 => Some(ColorRepresentation::AntiFundamental),
            6 => Some(ColorRepresentation::Sextet),
            -6 => Some(ColorRepresentation::AntiSextet),
            8 => Some(ColorRepresentation::Adjoint),
            _ => None,
        }
    }

    /// Exchange a color space with its conjugate space.
    pub fn dual(self) -> Self {
        match self {
            Self::Fundamental => Self::AntiFundamental,
            Self::AntiFundamental => Self::Fundamental,
            Self::Sextet => Self::AntiSextet,
            Self::AntiSextet => Self::Sextet,
            Self::Adjoint => Self::Adjoint,
        }
    }

    /// Return the typed Spenso representation, with its UFO dimension.
    pub fn representation(self) -> Representation<LibraryRep> {
        match self {
            Self::Fundamental => ColorFundamental {}.new_rep(3).cast(),
            Self::AntiFundamental => ColorFundamental {}.dual().new_rep(3).cast(),
            Self::Sextet => ColorSextet {}.new_rep(6).cast(),
            Self::AntiSextet => ColorSextet {}.dual().new_rep(6).cast(),
            Self::Adjoint => ColorAdjoint {}.new_rep(8).cast(),
        }
    }
}

/// An unsupported UFO color representation.
#[derive(Debug, Error)]
#[error("color sums for UFO representation {0} are not implemented")]
pub struct ColorSumError(pub i64);

/// Color completeness for a particle and its conjugate color space.
///
/// This is a Spenso identity, not a second color-algebra implementation.
/// Singlets contribute one. The signed UFO representation selects the dual
/// orientation for antiparticles; averaging divides by the representation size.
pub struct ColorSum {
    representation: Option<ColorRepresentation>,
    dimension: u64,
    average: bool,
}

impl ColorSum {
    pub fn new(particle: &Particle) -> Result<Self, ColorSumError> {
        let representation = ColorRepresentation::from_ufo(particle.color);
        if representation.is_none() && particle.color != 1 {
            return Err(ColorSumError(particle.color));
        }
        Ok(Self {
            representation,
            dimension: particle.color.unsigned_abs(),
            average: false,
        })
    }

    /// Divide by the dimension of the particle's color space.
    pub fn averaged(mut self, average: bool) -> Self {
        self.average = average;
        self
    }

    /// Return the completeness tensor with bare particle and conjugate indices.
    pub fn expression(&self, indices: [Atom; 2]) -> Atom {
        let Some(representation) = self.representation else {
            return Atom::one();
        };
        let [left, right] = indices;
        let identity = ETS.metric(
            representation.representation().to_symbolic([left]),
            representation.dual().representation().to_symbolic([right]),
        );
        if self.average {
            identity / Atom::num(self.dimension)
        } else {
            identity
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use feynkit_model::Model;
    use idenso::tensor::{ContractionSettings, SymbolicTensor};
    use symbolica::parse;

    #[test]
    fn color_completeness_has_unit_averaged_trace_in_every_supported_representation() {
        let model =
            Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap();
        let mut particle = model
            .particle_by_id(model.particle_id_by_pdg(5).unwrap())
            .unwrap()
            .clone();
        for color in [1, 3, -3, 6, -6, 8] {
            particle.color = color;
            let sum = ColorSum::new(&particle).unwrap();
            let trace = SymbolicTensor::infer(sum.expression([parse!("i"), parse!("i")]))
                .unwrap()
                .contract(ContractionSettings::default().without_rank_one_tensors())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression();
            assert_eq!(trace, Atom::num(color.unsigned_abs()));
            let average =
                SymbolicTensor::infer(sum.averaged(true).expression([parse!("i"), parse!("i")]))
                    .unwrap()
                    .contract(ContractionSettings::default().without_rank_one_tensors())
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression();
            assert_eq!(average, Atom::one());
        }
        particle.color = 10;
        assert!(matches!(ColorSum::new(&particle), Err(ColorSumError(10))));
    }

    #[test]
    fn particle_and_antiparticle_exchange_dual_ports() {
        let model =
            Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap();
        let particle = model
            .particle_by_id(model.particle_id_by_pdg(5).unwrap())
            .unwrap();
        let antiparticle = model.particle_by_id(particle.antiparticle).unwrap();
        let quark = ColorSum::new(particle)
            .unwrap()
            .expression([parse!("i"), parse!("j")]);
        let antiquark = ColorSum::new(antiparticle)
            .unwrap()
            .expression([parse!("j"), parse!("i")]);
        assert_eq!(quark, antiquark);
    }
}
