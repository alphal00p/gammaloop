use std::fmt::Display;

use crate::momentum::{FourMomentum, Polarization, Rotatable, Rotation, ThreeMomentum};
use crate::utils::{FloatLike, Length, F};
use crate::{
    define_index, define_indexed_vec, DependentMomentaConstructor, Externals, Polarizations,
    Settings,
};
use bincode_trait_derive::{Decode, Encode};
use derive_more::{From, Into};
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use std::ops::{Add, Index, IndexMut, Sub};
use symbolica::domains::float::NumericalFloatLike;
use tabled::settings::Style;
use typed_index_collections::TiVec;
use uuid::Uuid;

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
pub struct LoopIndex(pub usize);

// #[derive(
//     From,
//     Into,
//     Copy,
//     Clone,
//     Hash,
//     Eq,
//     Ord,
//     PartialEq,
//     PartialOrd,
//     Encode,
//     Decode,
//     Debug,
//     Serialize,
//     Deserialize,
// )]
// pub struct ExternalIndex(pub usize);

define_index!(
    pub struct ExternalIndex;
);

// define_indexed_vec!(
//     EdgeIndex;

//     pub struct ExternalThr;
// );

impl Display for LoopIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Display for ExternalIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[derive(From, Into, Serialize, Deserialize, Clone, PartialEq, Debug, Encode, Decode)]
pub struct LoopMomenta<T>(pub Vec<ThreeMomentum<T>>);

impl<T: Display> Display for LoopMomenta<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut table_builder = tabled::builder::Builder::new();
        for (i, mom) in self.0.iter().enumerate() {
            let mut row = Vec::new();
            row.push(i.to_string());
            for m in mom {
                row.push(m.to_string());
            }
            table_builder.push_record(row);
        }

        write!(f, "{}", table_builder.build().with(Style::modern_rounded()))
    }
}

impl<T> Length for LoopMomenta<T> {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }
}

pub type Subspace<'a> = Option<&'a [LoopIndex]>; // None means full space

impl<T> LoopMomenta<T> {
    pub(crate) fn iter(&self) -> std::slice::Iter<'_, ThreeMomentum<T>> {
        self.0.iter()
    }

    pub(crate) fn first(&self) -> Option<&ThreeMomentum<T>> {
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

impl<T> IndexMut<LoopIndex> for LoopMomenta<T> {
    fn index_mut(&mut self, index: LoopIndex) -> &mut Self::Output {
        &mut self.0[index.0]
    }
}

impl<T: FloatLike> LoopMomenta<F<T>> {
    pub(crate) fn hyper_radius_squared(&self, subspace: Subspace) -> F<T> {
        let zero = self.0[0].px.zero();
        match subspace {
            None => self.iter().fold(zero, |acc, x| acc + x.norm_squared()),
            Some(subspace) => subspace
                .iter()
                .fold(zero, |acc, &i| acc + self[i].norm_squared()),
        }
    }

    pub(crate) fn rescale(&self, factor: &F<T>, subspace: Subspace) -> Self {
        match subspace {
            None => LoopMomenta::from_iter(self.iter().map(|k| k * factor)),
            // this branch is wrong
            Some(subspace) => LoopMomenta::from_iter(subspace.iter().map(|&i| &self[i] * factor)),
        }
    }

    pub(crate) fn rotate(&self, rotation: &Rotation) -> Self {
        LoopMomenta::from_iter(self.iter().map(|k| k.rotate(rotation)))
    }
}

impl LoopMomenta<F<f64>> {
    pub(crate) fn cast<T: FloatLike>(&self) -> LoopMomenta<F<T>> {
        LoopMomenta::from_iter(self.iter().map(|m| m.map(&|x| F::from_ff64(x))))
    }
}

impl<T: FloatLike> Sub<&LoopMomenta<F<T>>> for &LoopMomenta<F<T>> {
    type Output = LoopMomenta<F<T>>;

    fn sub(self, rhs: &LoopMomenta<F<T>>) -> Self::Output {
        LoopMomenta::from_iter(self.iter().zip(rhs.iter()).map(|(l, r)| l - r))
    }
}

impl<T: FloatLike> Add<&LoopMomenta<F<T>>> for &LoopMomenta<F<T>> {
    type Output = LoopMomenta<F<T>>;

    fn add(self, rhs: &LoopMomenta<F<T>>) -> Self::Output {
        LoopMomenta::from_iter(self.iter().zip(rhs.iter()).map(|(l, r)| l + r))
    }
}

// define_indexed_vec!(
//     MyIdx;
//     pub struct ExternalMomentumtwo<ThreeMomentum<T>>;
// );
pub type ExternalThreeMomenta<T> = TiVec<ExternalIndex, ThreeMomentum<T>>;
// define_indexed_vec!()
pub type ExternalFourMomenta<T> = TiVec<ExternalIndex, FourMomentum<T>>;
// impl<T: Display> Display for ExternalFourMomenta<T> {
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         let mut table_builder = tabled::builder::Builder::new();
//         for (i, mom) in self.iter_enumerated() {
//             let mut row = Vec::new();
//             row.push(i.to_string());
//             for m in mom {
//                 row.push(m.to_string());
//             }
//             table_builder.push_record(row);
//         }

//         write!(f, "{}", table_builder.build().with(Style::modern_rounded()))
//     }
// }
pub type PolarizationVectors<T> = TiVec<ExternalIndex, Polarization<T>>; // should be the same length as #externals

#[derive(Debug, Clone)]
pub struct BareMomentumSample<T: FloatLike> {
    pub loop_moms: LoopMomenta<F<T>>,
    pub external_moms: ExternalFourMomenta<F<T>>,
    // pub polarizations: PolarizationVectors<Complex<F<T>>>,
    pub jacobian: F<T>,
    pub orientation: Option<usize>,
}

#[derive(Debug, Clone)]
pub struct MomentumSample<T: FloatLike> {
    pub sample: BareMomentumSample<T>,
    pub rotated_sample: Option<BareMomentumSample<T>>,
    pub uuid: Uuid,
}

impl<T: FloatLike> Display for MomentumSample<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut table = tabled::builder::Builder::new();
        if self.rotated_sample.is_some() {
            table.push_record(["Rotated Sample"]);
        } else {
            table.push_record(["Sample"]);
        }
        table.push_record(["Loop Momenta", "p_x", "p_y", "p_z"]);
        // table

        for (index, loop_mom) in self.loop_moms().0.iter().enumerate() {
            table.push_record([
                index.to_string(),
                loop_mom.px.to_string(),
                loop_mom.py.to_string(),
                loop_mom.pz.to_string(),
            ]);
        }

        table.push_record(["External Momenta", "E", "p_x", "p_y", "p_z"]);
        for (index, external_mom) in self.external_moms().iter_enumerated() {
            table.push_record([
                index.to_string(),
                external_mom.temporal.to_string(),
                external_mom.spatial.px.to_string(),
                external_mom.spatial.py.to_string(),
                external_mom.spatial.pz.to_string(),
            ]);
        }

        table.push_record(["Jacobian".into(), format!("{:+e}", self.sample.jacobian)]);
        table.build().with(Style::rounded()).fmt(f)
    }
}

impl<T: FloatLike> BareMomentumSample<T> {
    pub(crate) fn new(
        loop_moms: LoopMomenta<F<T>>,
        external_moms: &Externals,
        jacobian: F<T>,
        dependent_momenta_constructor: DependentMomentaConstructor,
        orientation: Option<usize>,
    ) -> Self {
        let external_moms = external_moms.get_dependent_externals(dependent_momenta_constructor);

        Self {
            loop_moms,
            external_moms,
            jacobian,
            orientation,
        }
    }

    pub(crate) fn one(&self) -> F<T> {
        if let Some(f) = self.loop_moms.first() {
            f.px.one()
        } else if let Some(f) = self.external_moms.first() {
            return f.spatial.px.one();
        } else {
            panic!("No momenta in sample")
        }
    }

    pub(crate) fn zero(&self) -> F<T> {
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

            jacobian: self.jacobian.clone().into(),
            orientation: self.orientation,
        }
    }

    pub(crate) fn higher_precision(&self) -> BareMomentumSample<T::Higher>
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
            jacobian: self.jacobian.higher(),
            orientation: self.orientation,
        }
    }

    pub(crate) fn lower_precision(&self) -> BareMomentumSample<T::Lower>
    where
        T::Lower: FloatLike,
    {
        BareMomentumSample {
            loop_moms: self.loop_moms.iter().map(ThreeMomentum::lower).collect(),
            external_moms: self.external_moms.iter().map(FourMomentum::lower).collect(),
            jacobian: self.jacobian.lower(),
            orientation: self.orientation,
        }
    }

    #[inline]
    /// Rotation for stability checks
    pub(crate) fn get_rotated_sample_cached(
        &self,
        rotation: &Rotation,
        rotated_externals: ExternalFourMomenta<F<T>>,
    ) -> Self {
        Self {
            loop_moms: self.loop_moms.iter().map(|l| l.rotate(rotation)).collect(),
            external_moms: rotated_externals,
            jacobian: self.jacobian.clone(),
            orientation: self.orientation,
        }
    }

    #[inline]
    pub(crate) fn get_rotated_sample(&self, rotation: &Rotation) -> Self {
        Self {
            loop_moms: self.loop_moms.iter().map(|l| l.rotate(rotation)).collect(),
            external_moms: self
                .external_moms
                .iter()
                .map(|l| l.rotate(rotation))
                .collect(),
            jacobian: self.jacobian.clone(),
            orientation: self.orientation,
        }
    }

    #[inline]
    pub(crate) fn rescaled_loop_momenta(&self, factor: &F<T>, subspace: Subspace) -> Self {
        Self {
            loop_moms: self.loop_moms.rescale(factor, subspace),
            external_moms: self.external_moms.clone(),
            jacobian: self.jacobian.clone(),
            orientation: self.orientation,
        }
    }
}

impl<T: FloatLike> MomentumSample<T> {
    pub(crate) fn possibly_rotated_sample(&self) -> &BareMomentumSample<T> {
        if let Some(rot) = self.rotated_sample.as_ref() {
            rot
        } else {
            &self.sample
        }
    }

    pub(crate) fn numerator_sample(
        &self,
        settings: &Settings,
    ) -> (&BareMomentumSample<T>, Option<Uuid>) {
        if settings.stability.rotate_numerator {
            (self.possibly_rotated_sample(), self.uuid())
        } else {
            (&self.sample, self.uuid())
        }
    }

    pub(crate) fn uuid(&self) -> Option<Uuid> {
        if self.rotated_sample.is_some() {
            None
        } else {
            Some(self.uuid)
        }
    }

    pub(crate) fn loop_moms(&self) -> &LoopMomenta<F<T>> {
        if let Some(rotated_sample) = &self.rotated_sample {
            &rotated_sample.loop_moms
        } else {
            &self.sample.loop_moms
        }
    }

    #[allow(clippy::type_complexity)]
    pub(crate) fn loop_mom_pair(&self) -> (&LoopMomenta<F<T>>, Option<&LoopMomenta<F<T>>>) {
        (
            &self.sample.loop_moms,
            self.rotated_sample.as_ref().map(|s| &s.loop_moms),
        )
    }

    pub(crate) fn external_moms(&self) -> &ExternalFourMomenta<F<T>> {
        if let Some(rotated_sample) = &self.rotated_sample {
            &rotated_sample.external_moms
        } else {
            &self.sample.external_moms
        }
    }

    #[allow(clippy::type_complexity)]
    pub(crate) fn external_mom_pair(
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

    pub(crate) fn jacobian(&self) -> F<T> {
        if let Some(rotated_sample) = &self.rotated_sample {
            rotated_sample.jacobian.clone()
        } else {
            self.sample.jacobian.clone()
        }
    }

    pub(crate) fn new(
        loop_moms: LoopMomenta<F<T>>,
        external_moms: &Externals,
        jacobian: F<T>,
        dependent_momenta_constructor: DependentMomentaConstructor,
        orientation: Option<usize>,
    ) -> Self {
        Self {
            sample: BareMomentumSample::new(
                loop_moms,
                external_moms,
                jacobian,
                dependent_momenta_constructor,
                orientation,
            ),
            rotated_sample: None,
            uuid: Uuid::new_v4(),
        }
    }

    pub(crate) fn one(&self) -> F<T> {
        self.sample.one()
    }

    pub(crate) fn zero(&self) -> F<T> {
        self.sample.zero()
    }

    /// Cast the sample to a different precision
    #[inline]
    pub(crate) fn cast_sample<T2: FloatLike>(&self) -> MomentumSample<T2>
    where
        F<T2>: From<F<T>>,
    {
        MomentumSample {
            sample: self.sample.cast_sample(),
            rotated_sample: self.rotated_sample.as_ref().map(|s| s.cast_sample()),
            uuid: self.uuid,
        }
    }

    pub(crate) fn higher_precision(&self) -> MomentumSample<T::Higher>
    where
        T::Higher: FloatLike,
    {
        MomentumSample {
            sample: self.sample.higher_precision(),
            rotated_sample: self.rotated_sample.as_ref().map(|s| s.higher_precision()),
            uuid: self.uuid,
        }
    }

    pub(crate) fn lower_precision(&self) -> MomentumSample<T::Lower>
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
    pub(crate) fn get_rotated_sample_cached(
        &self,
        rotation: &Rotation,
        rotated_externals: ExternalFourMomenta<F<T>>,
    ) -> Self {
        Self {
            sample: self.sample.clone(),
            rotated_sample: Some(
                self.sample
                    .get_rotated_sample_cached(rotation, rotated_externals),
            ),
            uuid: self.uuid,
        }
    }

    #[inline]
    pub(crate) fn get_rotated_sample(&self, rotation: &Rotation) -> Self {
        Self {
            sample: self.sample.clone(),
            rotated_sample: Some(self.sample.get_rotated_sample(rotation)),
            uuid: self.uuid,
        }
    }

    #[inline]
    pub(crate) fn rescaled_loop_momenta(&self, factor: &F<T>, subspace: Subspace) -> Self {
        Self {
            sample: self.sample.rescaled_loop_momenta(factor, subspace),
            rotated_sample: self
                .rotated_sample
                .as_ref()
                .map(|s| s.rescaled_loop_momenta(factor, subspace)),
            uuid: self.uuid,
        }
    }
}
