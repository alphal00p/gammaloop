use crate::momentum::FourMomentum;
use crate::utils::{F, FloatLike};
use smallvec::{SmallVec, smallvec};
use symbolica::domains::float::Real;

#[derive(Debug, Clone)]
pub struct Jet<T: FloatLike> {
    pub momentum: FourMomentum<F<T>>,
    pub constituents: SmallVec<[usize; 8]>,
    pt2: F<T>,
    rapidity: F<T>,
    phi: F<T>,
}

impl<T: FloatLike> Jet<T> {
    pub(crate) fn from_parts(
        momentum: FourMomentum<F<T>>,
        constituents: SmallVec<[usize; 8]>,
        pt2: F<T>,
        rapidity: F<T>,
        phi: F<T>,
    ) -> Self {
        Self {
            momentum,
            constituents,
            pt2,
            rapidity,
            phi,
        }
    }

    pub fn pt(&self) -> F<T> {
        self.pt2.clone().sqrt()
    }

    pub fn pt2(&self) -> F<T> {
        self.pt2.clone()
    }

    pub fn rapidity(&self) -> F<T> {
        self.rapidity.clone()
    }

    pub fn phi(&self) -> F<T> {
        self.phi.clone()
    }

    pub fn constituent_indices(&self) -> &[usize] {
        &self.constituents
    }

    pub fn to_f64(&self) -> Jet<f64> {
        Jet::from_parts(
            self.momentum.to_f64(),
            self.constituents.clone(),
            self.pt2.into_ff64(),
            self.rapidity.into_ff64(),
            self.phi.into_ff64(),
        )
    }

    pub fn from_f64(jet: &Jet<f64>) -> Self {
        Jet::from_parts(
            FourMomentum::from_ff64(&jet.momentum),
            jet.constituents.clone(),
            F::from_ff64(jet.pt2),
            F::from_ff64(jet.rapidity),
            F::from_ff64(jet.phi),
        )
    }
}

#[derive(Debug, Clone, Default)]
pub struct ClusteringResult<T: FloatLike> {
    pub jets: SmallVec<[Jet<T>; 8]>,
}

impl<T: FloatLike> ClusteringResult<T> {
    pub fn len(&self) -> usize {
        self.jets.len()
    }

    pub fn is_empty(&self) -> bool {
        self.jets.is_empty()
    }

    pub fn leading_jet(&self) -> Option<&Jet<T>> {
        self.jets.first()
    }

    pub fn ordered_pt(&self) -> SmallVec<[F<T>; 8]> {
        self.jets.iter().map(Jet::pt).collect()
    }

    pub fn to_f64(&self) -> ClusteringResult<f64> {
        ClusteringResult {
            jets: self.jets.iter().map(Jet::to_f64).collect(),
        }
    }

    pub fn from_f64(result: &ClusteringResult<f64>) -> Self {
        ClusteringResult {
            jets: result.jets.iter().map(Jet::from_f64).collect(),
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct PseudoJet<T: FloatLike> {
    pub(crate) momentum: FourMomentum<F<T>>,
    pub(crate) constituents: SmallVec<[usize; 8]>,
    pt2: F<T>,
    rapidity: F<T>,
    phi: F<T>,
}

impl<T: FloatLike> PseudoJet<T> {
    pub(crate) fn from_momentum(index: usize, momentum: FourMomentum<F<T>>) -> Self {
        Self::new(momentum, smallvec![index])
    }

    pub(crate) fn new(momentum: FourMomentum<F<T>>, constituents: SmallVec<[usize; 8]>) -> Self {
        let zero = momentum.temporal.value.zero();
        let mut jet = Self {
            momentum,
            constituents,
            pt2: zero.clone(),
            rapidity: zero.clone(),
            phi: zero,
        };
        jet.recompute_kinematics();
        jet
    }

    pub(crate) fn merged(lhs: &Self, rhs: &Self) -> Self {
        let mut constituents = lhs.constituents.clone();
        constituents.extend_from_slice(&rhs.constituents);
        constituents.sort_unstable();
        Self::new(lhs.momentum.clone() + rhs.momentum.clone(), constituents)
    }

    pub(crate) fn jet(self) -> Jet<T> {
        Jet::from_parts(
            self.momentum,
            self.constituents,
            self.pt2,
            self.rapidity,
            self.phi,
        )
    }

    pub(crate) fn pt(&self) -> F<T> {
        self.pt2.clone().sqrt()
    }

    pub(crate) fn pt2(&self) -> F<T> {
        self.pt2.clone()
    }

    pub(crate) fn rapidity(&self) -> F<T> {
        self.rapidity.clone()
    }

    pub(crate) fn phi(&self) -> F<T> {
        self.phi.clone()
    }

    fn recompute_kinematics(&mut self) {
        self.pt2 = transverse_momentum_squared(&self.momentum);
        self.phi = phi(&self.momentum);
        self.rapidity = rapidity(&self.momentum, &self.pt2);
    }
}

fn transverse_momentum_squared<T: FloatLike>(momentum: &FourMomentum<F<T>>) -> F<T> {
    momentum.spatial.px.square() + momentum.spatial.py.square()
}

fn phi<T: FloatLike>(momentum: &FourMomentum<F<T>>) -> F<T> {
    let zero = momentum.spatial.px.zero();
    let two_pi = momentum.spatial.px.TAU();
    let mut phi = momentum.spatial.py.atan2(&momentum.spatial.px);
    if phi < zero {
        phi += two_pi.clone();
    }
    if phi >= two_pi {
        phi -= two_pi;
    }
    phi
}

fn rapidity<T: FloatLike>(momentum: &FourMomentum<F<T>>, pt2: &F<T>) -> F<T> {
    momentum.rapidity_with_pt2(pt2)
}
