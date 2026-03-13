use crate::momentum::FourMomentum;
use crate::utils::{F, FloatLike, into_complex_ff64};
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;
use spenso::algebra::complex::Complex;
use std::ops::{Deref, DerefMut};

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct CutInfo {
    pub particle_pdgs: (SmallVec<[isize; 2]>, SmallVec<[isize; 4]>),
    pub cut_id: usize,
    pub graph_id: usize,
}

pub type Event = GenericEvent<f64>;
pub type EventGroup = GenericEventGroup<f64>;
pub type EventGroupList = GenericEventGroupList<f64>;

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct GenericEventGroup<T: FloatLike>(pub Vec<GenericEvent<T>>);

impl<T: FloatLike> Deref for GenericEventGroup<T> {
    type Target = Vec<GenericEvent<T>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: FloatLike> DerefMut for GenericEventGroup<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<T: FloatLike> GenericEventGroup<T> {
    pub fn to_f64(&self) -> EventGroup {
        GenericEventGroup(self.0.iter().map(GenericEvent::to_f64).collect())
    }

    pub fn from_f64(event_group: &EventGroup) -> Self {
        GenericEventGroup(event_group.0.iter().map(GenericEvent::from_f64).collect())
    }
}

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct GenericEventGroupList<T: FloatLike>(pub Vec<GenericEventGroup<T>>);

impl<T: FloatLike> Deref for GenericEventGroupList<T> {
    type Target = Vec<GenericEventGroup<T>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: FloatLike> DerefMut for GenericEventGroupList<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<T: FloatLike> GenericEventGroupList<T> {
    pub fn to_f64(&self) -> EventGroupList {
        GenericEventGroupList(self.0.iter().map(GenericEventGroup::to_f64).collect())
    }

    pub fn from_f64(event_group_list: &EventGroupList) -> Self {
        GenericEventGroupList(
            event_group_list
                .0
                .iter()
                .map(GenericEventGroup::from_f64)
                .collect(),
        )
    }
}

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct GenericReweightInfo<T: FloatLike> {
    pub weights: SmallVec<[Complex<F<T>>; 1]>,
}

impl<T: FloatLike> GenericReweightInfo<T> {
    pub fn to_f64(&self) -> GenericReweightInfo<f64> {
        GenericReweightInfo {
            weights: self.weights.iter().map(into_complex_ff64).collect(),
        }
    }

    pub fn from_f64(reweight_info: &GenericReweightInfo<f64>) -> Self {
        GenericReweightInfo {
            weights: reweight_info
                .weights
                .iter()
                .map(|weight| Complex::new(F::from_ff64(weight.re), F::from_ff64(weight.im)))
                .collect(),
        }
    }
}

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct GenericEvent<T: FloatLike> {
    #[allow(clippy::type_complexity)]
    pub kinematic_configuration: (
        SmallVec<[FourMomentum<F<T>>; 2]>,
        SmallVec<[FourMomentum<F<T>>; 4]>,
    ),
    pub cut_info: CutInfo,
    // Stores the integrand contribution associated with this event.
    pub integrand_weight: Complex<F<T>>,
    // Stores the observable-specific contribution associated with this event.
    pub observable_weight: Complex<F<T>>,
    // Contains additional partial weights for future use to do post-processing reweighting,
    // e.g. for scale variation or couplings adjustments
    pub reweight_info: GenericReweightInfo<T>,
}

impl<T: FloatLike> GenericEvent<T> {
    pub fn to_f64(&self) -> Event {
        GenericEvent {
            kinematic_configuration: (
                self.kinematic_configuration
                    .0
                    .iter()
                    .map(FourMomentum::to_f64)
                    .collect(),
                self.kinematic_configuration
                    .1
                    .iter()
                    .map(FourMomentum::to_f64)
                    .collect(),
            ),
            cut_info: self.cut_info.clone(),
            integrand_weight: into_complex_ff64(&self.integrand_weight),
            observable_weight: into_complex_ff64(&self.observable_weight),
            reweight_info: self.reweight_info.to_f64(),
        }
    }

    pub fn from_f64(event: &Event) -> Self {
        GenericEvent {
            kinematic_configuration: (
                event
                    .kinematic_configuration
                    .0
                    .iter()
                    .map(FourMomentum::from_ff64)
                    .collect(),
                event
                    .kinematic_configuration
                    .1
                    .iter()
                    .map(FourMomentum::from_ff64)
                    .collect(),
            ),
            cut_info: event.cut_info.clone(),
            integrand_weight: Complex::new(
                F::from_ff64(event.integrand_weight.re),
                F::from_ff64(event.integrand_weight.im),
            ),
            observable_weight: Complex::new(
                F::from_ff64(event.observable_weight.re),
                F::from_ff64(event.observable_weight.im),
            ),
            reweight_info: GenericReweightInfo::from_f64(&event.reweight_info),
        }
    }
}
