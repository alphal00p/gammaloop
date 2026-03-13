use crate::momentum::FourMomentum;
use crate::utils::{ArbPrec, F, FloatLike, f128};
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;
use spenso::algebra::complex::Complex;

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct CutInfo {
    pub particle_pdgs: (SmallVec<[isize; 2]>, SmallVec<[isize; 4]>),
    pub cut_id: usize,
    pub graph_id: usize,
}

pub enum Events {
    F64(Vec<Event<f64>>),
    F128(Vec<Event<f128>>),
    ArbPrec(Vec<Event<ArbPrec>>),
}

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct Event<T: FloatLike> {
    #[allow(clippy::type_complexity)]
    pub kinematic_configuration: (
        SmallVec<[FourMomentum<F<T>>; 2]>,
        SmallVec<[FourMomentum<F<T>>; 4]>,
    ),
    pub cut_info: CutInfo,
    // Stores the integrand weight associated with this event
    pub weight: Complex<F<T>>,
    // Contains additional partial weights for future use to do post-processing reweighting,
    // e.g. for scale variation or couplings adjustments
    pub reweight_info: SmallVec<[Complex<F<T>>; 1]>,
}
