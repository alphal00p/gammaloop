use std::fmt::Display;

use crate::momentum::{FourMomentum, Polarization, Rotatable, Rotation, ThreeMomentum};
use crate::signature::ExternalSignature;
use crate::utils::{FloatLike, F};
use crate::{Externals, Polarizations, Settings};
use bincode::{Decode, Encode};
use derive_more::{From, Into};
use serde::{Deserialize, Serialize};
use spenso::complex::Complex;
use std::ops::Index;
use symbolica::domains::float::NumericalFloatLike;
use typed_index_collections::TiVec;
use uuid::Uuid;

#[derive(From, Into, Copy, Clone, Hash, Eq, Ord, PartialEq, PartialOrd, Encode, Decode, Debug)]
pub struct LoopIndex(pub usize);
#[derive(
    From,
    Into,
    Copy,
    Clone,
    Hash,
    Eq,
    Ord,
    PartialEq,
    PartialOrd,
    Encode,
    Decode,
    Debug,
    Serialize,
    Deserialize,
)]
pub struct ExternalIndex(pub usize);

#[derive(From, Into, Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct LoopMomenta<T>(Vec<ThreeMomentum<T>>);

pub type Subspace<'a> = Option<&'a [LoopIndex]>; // None means full space

impl<T> LoopMomenta<T> {
    pub fn iter(&self) -> std::slice::Iter<'_, ThreeMomentum<T>> {
        self.0.iter()
    }

    pub fn first(&self) -> Option<&ThreeMomentum<T>> {
        self.0.first()
    }
}

impl<T> IntoIterator for LoopMomenta<T> {
    type Item = ThreeMomentum<T>;
    type IntoIter = std::vec::IntoIter<ThreeMomentum<T>>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<T> FromIterator<ThreeMomentum<T>> for LoopMomenta<T> {
    fn from_iter<I: IntoIterator<Item = ThreeMomentum<T>>>(iter: I) -> Self {
        LoopMomenta(iter.into_iter().collect())
    }
}

impl<T> Index<LoopIndex> for LoopMomenta<T> {
    type Output = ThreeMomentum<T>;

    fn index(&self, index: LoopIndex) -> &Self::Output {
        &self.0[index.0]
    }
}

impl<T: FloatLike> LoopMomenta<F<T>> {
    pub fn hyper_radius(&self, subspace: Subspace) -> F<T> {
        let zero = self.0[0].px.zero();
        match subspace {
            None => self.iter().fold(zero, |acc, x| acc + x.norm_squared()),
            Some(subspace) => subspace
                .iter()
                .fold(zero, |acc, &i| acc + self[i].norm_squared()),
        }
    }

    pub fn rescale(&self, factor: &F<T>, subspace: Subspace) -> Self {
        match subspace {
            None => LoopMomenta::from_iter(self.iter().map(|k| k * factor)),
            Some(subspace) => LoopMomenta::from_iter(subspace.iter().map(|&i| &self[i] * factor)),
        }
    }
}

pub type ExternalThreeMomenta<T> = TiVec<ExternalIndex, ThreeMomentum<T>>;
pub type ExternalFourMomenta<T> = TiVec<ExternalIndex, FourMomentum<T>>;
pub type PolarizationVectors<T> = TiVec<ExternalIndex, Polarization<T>>; // should be the same length as #externals

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BareMomentumSample<T: FloatLike> {
    pub loop_moms: LoopMomenta<F<T>>,
    pub external_moms: ExternalFourMomenta<F<T>>,
    pub polarizations: PolarizationVectors<Complex<F<T>>>,
    pub jacobian: F<T>,
}

pub struct MomentumSample<T: FloatLike> {
    pub sample: BareMomentumSample<T>,
    pub rotated_sample: Option<BareMomentumSample<T>>,
    pub uuid: Uuid,
}

impl<T: FloatLike> Display for MomentumSample<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Sample")?;
        write!(f, "\n\tloop momenta: ")?;
        for (index, loop_mom) in self.sample.loop_moms.iter().enumerate() {
            write!(f, "\n\t\tloop momentum {}: {}", index, loop_mom)?;
        }
        write!(f, "\n\texternal momenta: ")?;
        for (index, external_mom) in self.sample.external_moms.iter().enumerate() {
            write!(f, "\n\t\texternal momentum {}: {}", index, external_mom)?;
        }
        write!(f, "\n\tpolarizations: ")?;
        for (index, polarization) in self.sample.polarizations.iter().enumerate() {
            write!(f, "\n\t\tpolarization {}: {}", index, polarization)?;
        }
        write!(f, "\n\tjacobian: {:+e}", self.sample.jacobian)
    }
}

impl<T: FloatLike> BareMomentumSample<T> {
    pub fn new(
        loop_moms: LoopMomenta<F<T>>,
        external_moms: &Externals,
        jacobian: F<T>,
        polarizations: &Polarizations,
        external_signature: &ExternalSignature,
    ) -> Self {
        let polarizations = match polarizations {
            Polarizations::None => PolarizationVectors::new(),
            Polarizations::Constant { polarizations } => polarizations
                .iter()
                .map(|p| p.map(|c| c.map(|f| F::from_ff64(f))))
                .collect(),
        };

        let external_moms = external_moms.new_get_dependent_externals(external_signature);
        Self {
            polarizations,
            loop_moms,
            external_moms,
            jacobian,
        }
    }

    pub fn one(&self) -> F<T> {
        if let Some(f) = self.loop_moms.first() {
            f.px.one()
        } else if let Some(f) = self.external_moms.first() {
            return f.spatial.px.one();
        } else {
            panic!("No momenta in sample")
        }
    }

    pub fn zero(&self) -> F<T> {
        if let Some(f) = self.loop_moms.first() {
            f.px.zero()
        } else if let Some(f) = self.external_moms.first() {
            return f.spatial.px.zero();
        } else {
            panic!("No momenta in sample")
        }
    }

    /// Cast the sample to a different precision
    #[inline]
    fn cast_sample<T2: FloatLike>(&self) -> BareMomentumSample<T2>
    where
        F<T2>: From<F<T>>,
    {
        BareMomentumSample {
            loop_moms: self.loop_moms.iter().map(ThreeMomentum::cast).collect(),
            external_moms: self.external_moms.iter().map(FourMomentum::cast).collect(),
            polarizations: self
                .polarizations
                .iter()
                .map(Polarization::complex_cast)
                .collect(),
            jacobian: self.jacobian.clone().into(),
        }
    }

    pub fn higher_precision(&self) -> BareMomentumSample<T::Higher>
    where
        T::Higher: FloatLike,
    {
        BareMomentumSample {
            loop_moms: self.loop_moms.iter().map(ThreeMomentum::higher).collect(),
            external_moms: self
                .external_moms
                .iter()
                .map(FourMomentum::higher)
                .collect(),
            polarizations: self
                .polarizations
                .iter()
                .map(Polarization::higher)
                .collect(),
            jacobian: self.jacobian.higher(),
        }
    }

    pub fn lower_precision(&self) -> BareMomentumSample<T::Lower>
    where
        T::Lower: FloatLike,
    {
        BareMomentumSample {
            loop_moms: self.loop_moms.iter().map(ThreeMomentum::lower).collect(),
            external_moms: self.external_moms.iter().map(FourMomentum::lower).collect(),
            polarizations: self.polarizations.iter().map(Polarization::lower).collect(),
            jacobian: self.jacobian.lower(),
        }
    }

    #[inline]
    /// Rotation for stability checks
    pub fn get_rotated_sample_cached(
        &self,
        rotation: &Rotation,
        rotated_externals: ExternalFourMomenta<F<T>>,
        rotated_polarizations: PolarizationVectors<Complex<F<T>>>,
    ) -> Self {
        Self {
            loop_moms: self.loop_moms.iter().map(|l| l.rotate(rotation)).collect(),
            external_moms: rotated_externals,
            polarizations: rotated_polarizations,
            jacobian: self.jacobian.clone(),
        }
    }

    #[inline]
    pub fn get_rotated_sample(&self, rotation: &Rotation) -> Self {
        Self {
            loop_moms: self.loop_moms.iter().map(|l| l.rotate(rotation)).collect(),
            external_moms: self
                .external_moms
                .iter()
                .map(|l| l.rotate(rotation))
                .collect(),
            polarizations: self
                .polarizations
                .iter()
                .map(|l| l.rotate(rotation))
                .collect(),
            jacobian: self.jacobian.clone(),
        }
    }
}

impl<T: FloatLike> MomentumSample<T> {
    pub fn possibly_rotated_sample(&self) -> &BareMomentumSample<T> {
        if let Some(rot) = self.rotated_sample.as_ref() {
            rot
        } else {
            &self.sample
        }
    }

    pub fn numerator_sample(&self, settings: &Settings) -> (&BareMomentumSample<T>, Option<Uuid>) {
        if settings.stability.rotate_numerator {
            (self.possibly_rotated_sample(), self.uuid())
        } else {
            (&self.sample, self.uuid())
        }
    }

    pub fn uuid(&self) -> Option<Uuid> {
        if self.rotated_sample.is_some() {
            None
        } else {
            Some(self.uuid)
        }
    }

    pub fn loop_moms(&self) -> &LoopMomenta<F<T>> {
        if let Some(rotated_sample) = &self.rotated_sample {
            &rotated_sample.loop_moms
        } else {
            &self.sample.loop_moms
        }
    }

    #[allow(clippy::type_complexity)]
    pub fn loop_mom_pair(&self) -> (&LoopMomenta<F<T>>, Option<&LoopMomenta<F<T>>>) {
        (
            &self.sample.loop_moms,
            self.rotated_sample.as_ref().map(|s| &s.loop_moms),
        )
    }

    pub fn external_moms(&self) -> &ExternalFourMomenta<F<T>> {
        if let Some(rotated_sample) = &self.rotated_sample {
            &rotated_sample.external_moms
        } else {
            &self.sample.external_moms
        }
    }

    #[allow(clippy::type_complexity)]
    pub fn external_mom_pair(
        &self,
    ) -> (
        &ExternalFourMomenta<F<T>>,
        Option<&ExternalFourMomenta<F<T>>>,
    ) {
        (
            self.sample.external_moms.as_ref(),
            self.rotated_sample
                .as_ref()
                .map(|s| s.external_moms.as_ref()),
        )
    }

    pub fn polarizations(&self) -> &PolarizationVectors<Complex<F<T>>> {
        if let Some(rotated_sample) = &self.rotated_sample {
            &rotated_sample.polarizations
        } else {
            &self.sample.polarizations
        }
    }

    #[allow(clippy::type_complexity)]
    pub fn polarizations_pair(
        &self,
    ) -> (
        &PolarizationVectors<Complex<F<T>>>,
        Option<&PolarizationVectors<Complex<F<T>>>>,
    ) {
        (
            self.sample.polarizations.as_ref(),
            self.rotated_sample
                .as_ref()
                .map(|s| s.polarizations.as_ref()),
        )
    }

    pub fn jacobian(&self) -> F<T> {
        if let Some(rotated_sample) = &self.rotated_sample {
            rotated_sample.jacobian.clone()
        } else {
            self.sample.jacobian.clone()
        }
    }

    pub fn new(
        loop_moms: LoopMomenta<F<T>>,
        external_moms: &Externals,
        jacobian: F<T>,
        polarizations: &Polarizations,
        external_signature: &ExternalSignature,
    ) -> Self {
        Self {
            sample: BareMomentumSample::new(
                loop_moms,
                external_moms,
                jacobian,
                polarizations,
                external_signature,
            ),
            rotated_sample: None,
            uuid: Uuid::new_v4(),
        }
    }

    pub fn one(&self) -> F<T> {
        self.sample.one()
    }

    pub fn zero(&self) -> F<T> {
        self.sample.zero()
    }

    /// Cast the sample to a different precision
    #[inline]
    fn cast_sample<T2: FloatLike>(&self) -> MomentumSample<T2>
    where
        F<T2>: From<F<T>>,
    {
        MomentumSample {
            sample: self.sample.cast_sample(),
            rotated_sample: self.rotated_sample.as_ref().map(|s| s.cast_sample()),
            uuid: self.uuid,
        }
    }

    pub fn higher_precision(&self) -> MomentumSample<T::Higher>
    where
        T::Higher: FloatLike,
    {
        MomentumSample {
            sample: self.sample.higher_precision(),
            rotated_sample: self.rotated_sample.as_ref().map(|s| s.higher_precision()),
            uuid: self.uuid,
        }
    }

    pub fn lower_precision(&self) -> MomentumSample<T::Lower>
    where
        T::Lower: FloatLike,
    {
        MomentumSample {
            sample: self.sample.lower_precision(),
            rotated_sample: self.rotated_sample.as_ref().map(|s| s.lower_precision()),
            uuid: self.uuid,
        }
    }

    #[inline]
    /// Rotation for stability checks
    pub fn get_rotated_sample_cached(
        &self,
        rotation: &Rotation,
        rotated_externals: ExternalFourMomenta<F<T>>,
        rotated_polarizations: PolarizationVectors<Complex<F<T>>>,
    ) -> Self {
        Self {
            sample: self.sample.clone(),
            rotated_sample: Some(self.sample.get_rotated_sample_cached(
                rotation,
                rotated_externals,
                rotated_polarizations,
            )),
            uuid: self.uuid,
        }
    }

    #[inline]
    pub fn get_rotated_sample(&self, rotation: &Rotation) -> Self {
        Self {
            sample: self.sample.clone(),
            rotated_sample: Some(self.sample.get_rotated_sample(rotation)),
            uuid: self.uuid,
        }
    }
}
