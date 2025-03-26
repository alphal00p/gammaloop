use std::{
    borrow::Borrow,
    fmt::{Display, LowerExp},
    ops::{Add, AddAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign},
};

use bincode::{Decode, Encode};
use eyre::Context;
use momtrop::vector::Vector;
use serde::{Deserialize, Serialize};
use serde_repr::{Deserialize_repr, Serialize_repr};
use spenso::{
    complex::RealOrComplexTensor,
    contraction::{Contract, RefZero},
    data::{
        DataIterator, DataTensor, DenseTensor, HasTensorData, SetTensorData, SparseTensor,
        StorageTensor,
    },
    iterators::IteratableTensor,
    parametric::{
        atomcore::TensorAtomOps, EvalTensor, FlatCoefficent, MixedTensor, ParamOrConcrete,
    },
    shadowing::Shadowable,
    structure::{
        abstract_index::AbstractIndex,
        representation::{BaseRepName, Bispinor, Euclidean, Minkowski, PhysReps, RepName},
        slot::{DualSlotTo, Slot},
        CastStructure, IndexLess, NamedStructure, TensorStructure, ToSymbolic, VecStructure,
    },
    symbolica_utils::NoArgs,
    upgrading_arithmetic::FallibleAdd,
};
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    coefficient::Coefficient,
    domains::{
        float::{NumericalFloatLike, Real, RealNumberLike, SingleFloat},
        integer::IntegerRing,
        rational::{Rational, RationalField},
    },
    evaluate::{ExpressionEvaluator, FunctionMap},
    poly::{polynomial::MultivariatePolynomial, Exponent},
};

use spenso::complex::Complex;

use crate::{
    utils::{ApproxEq, FloatLike, RefDefault, F},
    RotationSetting,
};

#[derive(Debug, PartialEq, Eq, Clone, Copy, Serialize, Deserialize)]
pub struct Energy<T> {
    pub value: T,
}

impl<T> Energy<T> {
    pub fn map_ref<U>(&self, f: &impl Fn(&T) -> U) -> Energy<U> {
        Energy {
            value: f(&self.value),
        }
    }

    pub fn map<U>(self, f: &impl Fn(T) -> U) -> Energy<U> {
        Energy {
            value: f(self.value),
        }
    }
}

impl<T: FloatLike> ApproxEq<Energy<F<T>>, F<T>> for Energy<F<T>> {
    fn approx_eq(&self, other: &Energy<F<T>>, threshold: &F<T>) -> bool {
        self.value.approx_eq(&other.value, threshold)
    }
}

impl<T: FloatLike> Energy<F<T>> {
    pub fn higher(&self) -> Energy<F<T::Higher>>
    where
        T::Higher: FloatLike,
    {
        Energy {
            value: self.value.higher(),
        }
    }

    pub fn lower(&self) -> Energy<F<T::Lower>>
    where
        T::Lower: FloatLike,
    {
        Energy {
            value: self.value.lower(),
        }
    }

    pub fn from_ff64(energy: Energy<F<f64>>) -> Self {
        Energy {
            value: F::from_ff64(energy.value),
        }
    }
}

impl<T: FloatLike> From<Energy<T>> for Energy<F<T>> {
    fn from(value: Energy<T>) -> Self {
        Energy {
            value: F(value.value),
        }
    }
}

impl<T: Real> Energy<T> {
    pub fn zero(&self) -> Self {
        Energy {
            value: self.value.zero(),
        }
    }
}

impl<T> Add<Energy<T>> for Energy<T>
where
    T: Add<T, Output = T>,
{
    type Output = Energy<T>;
    fn add(self, rhs: Energy<T>) -> Self::Output {
        Energy {
            value: self.value + rhs.value,
        }
    }
}

impl<T> Add<&Energy<T>> for Energy<T>
where
    T: for<'a> Add<&'a T, Output = T>,
{
    type Output = Energy<T>;
    fn add(self, rhs: &Energy<T>) -> Self::Output {
        Energy {
            value: self.value + &rhs.value,
        }
    }
}

impl<'b, T> Add<&Energy<T>> for &'b Energy<T>
where
    &'b T: for<'a> Add<&'a T, Output = T>,
{
    type Output = Energy<T>;
    fn add(self, rhs: &Energy<T>) -> Self::Output {
        Energy {
            value: &self.value + &rhs.value,
        }
    }
}

impl<T> Add<Energy<T>> for &Energy<T>
where
    T: for<'a> Add<&'a T, Output = T>,
{
    type Output = Energy<T>;
    fn add(self, rhs: Energy<T>) -> Self::Output {
        rhs + self
    }
}

impl<T> AddAssign<Energy<T>> for Energy<T>
where
    T: AddAssign<T>,
{
    fn add_assign(&mut self, rhs: Energy<T>) {
        self.value += rhs.value;
    }
}

impl<'a, T> AddAssign<&'a Energy<T>> for Energy<T>
where
    T: AddAssign<&'a T>,
{
    fn add_assign(&mut self, rhs: &'a Energy<T>) {
        self.value += &rhs.value;
    }
}

impl<T> Sub<Energy<T>> for Energy<T>
where
    T: Sub<T, Output = T>,
{
    type Output = Energy<T>;
    fn sub(self, rhs: Energy<T>) -> Self::Output {
        Energy {
            value: self.value - rhs.value,
        }
    }
}

impl<T> SubAssign<Energy<T>> for Energy<T>
where
    T: SubAssign<T>,
{
    fn sub_assign(&mut self, rhs: Energy<T>) {
        self.value -= rhs.value;
    }
}

impl<'a, T> SubAssign<&'a Energy<T>> for Energy<T>
where
    T: SubAssign<&'a T>,
{
    fn sub_assign(&mut self, rhs: &'a Energy<T>) {
        self.value -= &rhs.value;
    }
}
impl<T> Mul<Energy<T>> for Energy<T>
where
    T: Mul<T, Output = T>,
{
    type Output = T;
    fn mul(self, rhs: Energy<T>) -> Self::Output {
        self.value * rhs.value
    }
}

impl<T> MulAssign<Energy<T>> for Energy<T>
where
    T: MulAssign<T>,
{
    fn mul_assign(&mut self, rhs: Energy<T>) {
        self.value *= rhs.value;
    }
}

impl<T> Neg for Energy<T>
where
    T: Neg<Output = T>,
{
    type Output = Energy<T>;
    fn neg(self) -> Self::Output {
        Energy { value: -self.value }
    }
}

impl<T> Energy<T> {
    pub fn new(value: T) -> Self {
        Energy { value }
    }

    // pub fn from_three_momentum(three_momentum: &ThreeMomentum<T>) -> Self
    // where
    //     T: std::ops::Mul<Output = T> + std::ops::Add<Output = T> + Copy,
    // {
    //     let px2 = three_momentum.px * three_momentum.px;
    //     let py2 = three_momentum.py * three_momentum.py;
    //     let pz2 = three_momentum.pz * three_momentum.pz;
    //     let p2 = px2 + py2 + pz2;
    //     let value = (p2 + T::default()).sqrt();
    //     Energy { value }
    // }
}

// impl<U: Borrow<T>, T> Borrow<Energy<T>> for Energy<U> {
//     fn borrow(&self) -> &Energy<T> {

//     }
// }

impl<U, T: RefZero<U>> RefZero<Energy<U>> for Energy<T>
where
    Energy<T>: Borrow<Energy<U>>,
{
    fn ref_zero(&self) -> Energy<U> {
        Energy {
            value: self.value.ref_zero(),
        }
    }
}

impl Energy<Atom> {
    pub fn new_parametric(id: usize) -> Self {
        let value = Atom::parse(&format!("E_{}", id)).unwrap();
        Energy { value }
    }
}

impl<T: Display> Display for Energy<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "E: {}", self.value)
    }
}

impl<T: LowerExp> LowerExp for Energy<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "E: {:e}", self.value)
    }
}

impl<T: Default> Default for Energy<T> {
    fn default() -> Self {
        Energy {
            value: T::default(),
        }
    }
}

impl<T> From<T> for Energy<T> {
    fn from(value: T) -> Self {
        Energy { value }
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy, Serialize, Deserialize)]
pub struct ThreeMomentum<T> {
    pub px: T,
    pub py: T,
    pub pz: T,
}

impl<T> ThreeMomentum<T> {
    pub fn map_ref<U>(&self, f: &impl Fn(&T) -> U) -> ThreeMomentum<U> {
        ThreeMomentum {
            px: f(&self.px),
            py: f(&self.py),
            pz: f(&self.pz),
        }
    }

    pub fn map<U>(self, f: &impl Fn(T) -> U) -> ThreeMomentum<U> {
        ThreeMomentum {
            px: f(self.px),
            py: f(self.py),
            pz: f(self.pz),
        }
    }
}

impl<T: FloatLike> ApproxEq<ThreeMomentum<F<T>>, F<T>> for ThreeMomentum<F<T>> {
    fn approx_eq(&self, other: &ThreeMomentum<F<T>>, threshold: &F<T>) -> bool {
        F::approx_eq_iterator(
            [&self.px, &self.py, &self.pz],
            [&other.px, &other.py, &other.pz],
            threshold,
        )
    }
}

pub struct ThreeRotation<T> {
    pub map: fn(ThreeMomentum<T>) -> ThreeMomentum<T>,
    pub inv_map: fn(ThreeMomentum<T>) -> ThreeMomentum<T>,
}

impl<T: Neg<Output = T>> ThreeRotation<T> {
    pub fn half_pi_x() -> Self {
        let map = |mut momentum: ThreeMomentum<T>| -> ThreeMomentum<T> {
            std::mem::swap(&mut momentum.py, &mut momentum.pz);
            momentum.pz = -momentum.pz;
            momentum
        };

        let inv_map = |mut momentum: ThreeMomentum<T>| -> ThreeMomentum<T> {
            momentum.pz = -momentum.pz;
            std::mem::swap(&mut momentum.py, &mut momentum.pz);
            momentum
        };

        ThreeRotation { map, inv_map }
    }

    pub fn half_pi_y() -> Self {
        let map = |mut momentum: ThreeMomentum<T>| -> ThreeMomentum<T> {
            std::mem::swap(&mut momentum.px, &mut momentum.pz);
            momentum.px = -momentum.px;
            momentum
        };

        let inv_map = |mut momentum: ThreeMomentum<T>| -> ThreeMomentum<T> {
            momentum.px = -momentum.px;
            std::mem::swap(&mut momentum.px, &mut momentum.pz);
            momentum
        };

        ThreeRotation { map, inv_map }
    }

    pub fn half_pi_z() -> Self {
        let map = |mut momentum: ThreeMomentum<T>| -> ThreeMomentum<T> {
            std::mem::swap(&mut momentum.px, &mut momentum.py);
            momentum.py = -momentum.py;
            momentum
        };

        let inv_map = |mut momentum: ThreeMomentum<T>| -> ThreeMomentum<T> {
            momentum.py = -momentum.py;
            std::mem::swap(&mut momentum.px, &mut momentum.py);
            momentum
        };

        ThreeRotation { map, inv_map }
    }
}

impl<T: FloatLike> ThreeMomentum<F<T>> {
    pub fn higher(&self) -> ThreeMomentum<F<T::Higher>>
    where
        T::Higher: FloatLike,
    {
        ThreeMomentum {
            px: self.px.higher(),
            py: self.py.higher(),
            pz: self.pz.higher(),
        }
    }

    pub fn lower(&self) -> ThreeMomentum<F<T::Lower>>
    where
        T::Lower: FloatLike,
    {
        ThreeMomentum {
            px: self.px.lower(),
            py: self.py.lower(),
            pz: self.pz.lower(),
        }
    }

    pub fn from_ff64(three_mom: ThreeMomentum<F<f64>>) -> Self {
        ThreeMomentum {
            px: F::from_ff64(three_mom.px),
            py: F::from_ff64(three_mom.py),
            pz: F::from_ff64(three_mom.pz),
        }
    }
}

impl<T: FloatLike> From<ThreeMomentum<T>> for ThreeMomentum<F<T>> {
    fn from(value: ThreeMomentum<T>) -> Self {
        ThreeMomentum {
            px: value.px.into(),
            py: value.py.into(),
            pz: value.pz.into(),
        }
    }
}

impl<T> IntoIterator for ThreeMomentum<T> {
    type Item = T;
    type IntoIter = std::array::IntoIter<T, 3>;

    fn into_iter(self) -> Self::IntoIter {
        let [px, py, pz] = [self.px, self.py, self.pz];
        [px, py, pz].into_iter()
    }
}

impl<'a, T> IntoIterator for &'a ThreeMomentum<T> {
    type Item = &'a T;
    type IntoIter = std::array::IntoIter<&'a T, 3>;

    fn into_iter(self) -> Self::IntoIter {
        let [px, py, pz] = [&self.px, &self.py, &self.pz];
        [px, py, pz].into_iter()
    }
}

impl<T: Real> RefDefault for ThreeMomentum<T> {
    fn default(&self) -> Self {
        let zero = self.px.zero();
        ThreeMomentum {
            px: zero.clone(),
            py: zero.clone(),
            pz: zero.clone(),
        }
    }
}

impl<T: Real> ThreeMomentum<T> {
    pub fn zero(&self) -> Self {
        let zero = self.px.zero();
        ThreeMomentum {
            px: zero.clone(),
            py: zero.clone(),
            pz: zero.clone(),
        }
    }
}

impl<T: RefZero> RefZero for ThreeMomentum<T> {
    fn ref_zero(&self) -> Self {
        ThreeMomentum {
            px: self.px.ref_zero(),
            py: self.py.ref_zero(),
            pz: self.pz.ref_zero(),
        }
    }
}

impl<T> ThreeMomentum<T> {
    pub fn new(px: T, py: T, pz: T) -> Self {
        ThreeMomentum { px, py, pz }
    }

    pub fn into_dense(self, index: AbstractIndex) -> DenseTensor<T, VecStructure>
    where
        T: Clone,
    {
        let structure =
            VecStructure::from_iter(vec![PhysReps::new_slot(Euclidean {}.into(), 3, index)]);
        DenseTensor::from_data(vec![self.px, self.py, self.pz], structure).unwrap()
    }

    pub fn into_dense_param(self, index: AbstractIndex) -> DenseTensor<T, VecStructure>
    where
        T: Clone,
    {
        let structure =
            VecStructure::from_iter(vec![PhysReps::new_slot(Euclidean {}.into(), 3, index)]);
        DenseTensor::from_data(vec![self.px, self.py, self.pz], structure).unwrap()
    }
}

impl<T: FloatLike> ThreeMomentum<F<T>> {
    pub fn to_f64(&self) -> ThreeMomentum<F<f64>> {
        ThreeMomentum {
            px: F(self.px.to_f64()),
            py: F(self.py.to_f64()),
            pz: F(self.pz.to_f64()),
        }
    }
    /// Compute the phi-angle separation with p2.
    pub fn getdelphi(&self, p2: &ThreeMomentum<F<T>>) -> F<T> {
        let pt1 = self.pt();
        let pt2 = p2.pt();
        if pt1.is_zero() {
            return pt1.max_value();
        }
        if pt2.is_zero() {
            return pt2.max_value();
        }

        let mut tmp = self.px.clone() * &p2.px + self.py.clone() * &p2.py;

        tmp /= pt1 * pt2;
        if tmp.norm() > tmp.one() + tmp.epsilon() {
            panic!("Cosine larger than 1. in phase-space cuts.")
        }
        if tmp.norm() > tmp.one() {
            (tmp.clone() / tmp.norm()).acos()
        } else {
            tmp.acos()
        }
    }

    /// Compute the deltaR separation with momentum p2.
    #[inline]
    pub fn delta_r(&self, p2: &ThreeMomentum<F<T>>) -> F<T>
    where
        T: Real,
    {
        let delta_eta = self.pseudo_rap() - p2.pseudo_rap();
        let delta_phi = self.getdelphi(p2);
        (delta_eta.square() + delta_phi.square()).sqrt()
    }

    pub fn rotate_mut(&mut self, alpha: &F<T>, beta: &F<T>, gamma: &F<T>) {
        let sin_alpha = alpha.sin();
        let cos_alpha = alpha.cos();
        let sin_beta = beta.sin();
        let cos_beta = beta.cos();
        let sin_gamma = gamma.sin();
        let cos_gamma = gamma.cos();

        let px = self.px.clone();
        let py = self.py.clone();
        let pz = self.pz.clone();

        self.px = cos_gamma.clone() * &cos_beta * &px
            + (-(cos_alpha.clone()) * &sin_gamma + sin_alpha.clone() * &sin_beta * &cos_gamma)
                * &py
            + (sin_alpha.clone() * &sin_gamma + cos_alpha.clone() * &sin_beta * &cos_gamma) * &pz;

        self.py = sin_gamma.clone() * &cos_beta * &px
            + (cos_alpha.clone() * &cos_gamma + sin_alpha.clone() * &sin_beta * &sin_gamma) * &py
            + (-sin_alpha.clone() * &cos_gamma + cos_alpha.clone() * &sin_beta * &sin_gamma) * &pz;

        self.pz =
            -sin_beta * &px + cos_beta.clone() * &sin_alpha * &py + cos_alpha * &cos_beta * &pz;
    }

    /// Compute transverse momentum.
    #[inline]
    pub fn pt(&self) -> F<T> {
        (self.px.square() + self.py.square()).sqrt()
    }

    /// Compute pseudorapidity.
    #[inline]
    pub fn pseudo_rap(&self) -> F<T> {
        let pt = self.pt();
        if pt.less_than_epsilon() && self.pz.norm().less_than_epsilon() {
            if self.pz.positive() {
                return pt.max_value();
            } else {
                return pt.min_value();
            }
        }
        let th = pt.atan2(&self.pz);
        let two = pt.from_i64(2);
        -(th / two).tan().ln()
    }
}

impl<T: Neg<Output = T> + Clone> ThreeMomentum<T> {
    pub fn perform_pi2_rotation_x_mut(&mut self) {
        self.pz = -self.pz.clone();
        std::mem::swap(&mut self.pz, &mut self.py);
    }

    pub fn perform_pi2_rotation_x(&self) -> Self {
        Self {
            px: self.px.clone(),
            py: -self.pz.clone(),
            pz: self.py.clone(),
        }
    }

    pub fn perform_pi2_rotation_y_mut(&mut self) {
        self.px = -self.px.clone();
        std::mem::swap(&mut self.px, &mut self.pz);
    }

    pub fn perform_pi2_rotation_y(&self) -> Self {
        Self {
            px: self.pz.clone(),
            py: self.py.clone(),
            pz: -self.px.clone(),
        }
    }

    pub fn perform_pi2_rotation_z_mut(&mut self) {
        self.py = -self.py.clone();
        std::mem::swap(&mut self.px, &mut self.py);
    }

    pub fn perform_pi2_rotation_z(&self) -> Self {
        Self {
            px: -self.py.clone(),
            py: self.px.clone(),
            pz: self.pz.clone(),
        }
    }
}

impl<T> Add<ThreeMomentum<T>> for ThreeMomentum<T>
where
    T: Add<T, Output = T>,
{
    type Output = ThreeMomentum<T>;
    fn add(self, rhs: ThreeMomentum<T>) -> Self::Output {
        ThreeMomentum {
            px: self.px + rhs.px,
            py: self.py + rhs.py,
            pz: self.pz + rhs.pz,
        }
    }
}

impl<T> Add<&ThreeMomentum<T>> for ThreeMomentum<T>
where
    T: for<'a> Add<&'a T, Output = T>,
{
    type Output = ThreeMomentum<T>;
    fn add(self, rhs: &ThreeMomentum<T>) -> Self::Output {
        ThreeMomentum {
            px: self.px + &rhs.px,
            py: self.py + &rhs.py,
            pz: self.pz + &rhs.pz,
        }
    }
}

impl<'b, T> Add<&ThreeMomentum<T>> for &'b ThreeMomentum<T>
where
    &'b T: for<'a> Add<&'a T, Output = T>,
{
    type Output = ThreeMomentum<T>;
    fn add(self, rhs: &ThreeMomentum<T>) -> Self::Output {
        ThreeMomentum {
            px: &self.px + &rhs.px,
            py: &self.py + &rhs.py,
            pz: &self.pz + &rhs.pz,
        }
    }
}

impl<T> Add<ThreeMomentum<T>> for &ThreeMomentum<T>
where
    T: for<'a> Add<&'a T, Output = T>,
{
    type Output = ThreeMomentum<T>;
    fn add(self, rhs: ThreeMomentum<T>) -> Self::Output {
        rhs + self
    }
}

impl<T> AddAssign<ThreeMomentum<T>> for ThreeMomentum<T>
where
    T: AddAssign<T>,
{
    fn add_assign(&mut self, rhs: ThreeMomentum<T>) {
        self.px += rhs.px;
        self.py += rhs.py;
        self.pz += rhs.pz;
    }
}

impl<'a, T> AddAssign<&'a ThreeMomentum<T>> for ThreeMomentum<T>
where
    T: AddAssign<&'a T>,
{
    fn add_assign(&mut self, rhs: &'a ThreeMomentum<T>) {
        self.px += &rhs.px;
        self.py += &rhs.py;
        self.pz += &rhs.pz;
    }
}

impl<T> Sub<ThreeMomentum<T>> for ThreeMomentum<T>
where
    T: Sub<T, Output = T>,
{
    type Output = ThreeMomentum<T>;
    fn sub(self, rhs: ThreeMomentum<T>) -> Self::Output {
        ThreeMomentum {
            px: self.px - rhs.px,
            py: self.py - rhs.py,
            pz: self.pz - rhs.pz,
        }
    }
}

impl<T> SubAssign<ThreeMomentum<T>> for ThreeMomentum<T>
where
    T: SubAssign<T>,
{
    fn sub_assign(&mut self, rhs: ThreeMomentum<T>) {
        self.px -= rhs.px;
        self.py -= rhs.py;
        self.pz -= rhs.pz;
    }
}

impl<'a, T> SubAssign<&'a ThreeMomentum<T>> for ThreeMomentum<T>
where
    T: SubAssign<&'a T>,
{
    fn sub_assign(&mut self, rhs: &'a ThreeMomentum<T>) {
        self.px -= &rhs.px;
        self.py -= &rhs.py;
        self.pz -= &rhs.pz;
    }
}

impl<T> Mul<ThreeMomentum<T>> for ThreeMomentum<T>
where
    T: Mul<T, Output = T> + Add<T, Output = T>,
{
    type Output = T;
    fn mul(self, rhs: ThreeMomentum<T>) -> Self::Output {
        self.px * rhs.px + self.py * rhs.py + self.pz * rhs.pz
    }
}

impl<T> Mul<&ThreeMomentum<T>> for ThreeMomentum<T>
where
    T: for<'a> Mul<&'a T, Output = T> + Add<T, Output = T>,
{
    type Output = T;
    fn mul(self, rhs: &ThreeMomentum<T>) -> Self::Output {
        self.px * &rhs.px + self.py * &rhs.py + self.pz * &rhs.pz
    }
}

impl<T> Mul<ThreeMomentum<T>> for &ThreeMomentum<T>
where
    T: for<'b> Mul<&'b T, Output = T> + Add<T, Output = T>,
{
    type Output = T;
    fn mul(self, rhs: ThreeMomentum<T>) -> Self::Output {
        rhs * self
    }
}

impl<T> Mul<T> for ThreeMomentum<T>
where
    T: Mul<T, Output = T> + Clone,
{
    type Output = ThreeMomentum<T>;
    fn mul(self, rhs: T) -> Self::Output {
        ThreeMomentum {
            px: self.px * rhs.clone(),
            py: self.py * rhs.clone(),
            pz: self.pz * rhs,
        }
    }
}

impl<T> Mul<&T> for ThreeMomentum<T>
where
    T: for<'a> Mul<&'a T, Output = T> + Clone,
{
    type Output = ThreeMomentum<T>;
    fn mul(self, rhs: &T) -> Self::Output {
        ThreeMomentum {
            px: self.px * rhs,
            py: self.py * rhs,
            pz: self.pz * rhs,
        }
    }
}

impl<T> Mul<T> for &ThreeMomentum<T>
where
    T: Mul<T, Output = T> + Clone,
{
    type Output = ThreeMomentum<T>;
    fn mul(self, rhs: T) -> Self::Output {
        ThreeMomentum {
            px: self.px.clone() * rhs.clone(),
            py: self.py.clone() * rhs.clone(),
            pz: self.pz.clone() * rhs,
        }
    }
}

impl<T> Mul<&T> for &ThreeMomentum<T>
where
    T: for<'b> Mul<&'b T, Output = T> + Clone,
{
    type Output = ThreeMomentum<T>;
    fn mul(self, rhs: &T) -> Self::Output {
        ThreeMomentum {
            px: self.px.clone() * rhs,
            py: self.py.clone() * rhs,
            pz: self.pz.clone() * rhs,
        }
    }
}

impl<T> Neg for ThreeMomentum<T>
where
    T: Neg<Output = T>,
{
    type Output = ThreeMomentum<T>;
    fn neg(self) -> Self::Output {
        ThreeMomentum {
            px: -self.px,
            py: -self.py,
            pz: -self.pz,
        }
    }
}

impl<T> Neg for &ThreeMomentum<T>
where
    T: Neg<Output = T> + Clone,
{
    type Output = ThreeMomentum<T>;
    fn neg(self) -> Self::Output {
        ThreeMomentum {
            px: -self.px.clone(),
            py: -self.py.clone(),
            pz: -self.pz.clone(),
        }
    }
}

impl<T: Display> Display for ThreeMomentum<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "px: {}, py: {}, pz: {}", self.px, self.py, self.pz)
    }
}

impl<T: LowerExp> LowerExp for ThreeMomentum<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "px: {:e}, py: {:e}, pz: {:e}", self.px, self.py, self.pz)
    }
}

impl<T: Default> Default for ThreeMomentum<T> {
    fn default() -> Self {
        ThreeMomentum {
            px: T::default(),
            py: T::default(),
            pz: T::default(),
        }
    }
}

impl<T> ThreeMomentum<T> {
    pub fn norm_squared(&self) -> T
    where
        T: for<'a> Mul<&'a T, Output = T> + Add<T, Output = T> + Clone,
    {
        self.px.clone() * &self.px + self.py.clone() * &self.py + self.pz.clone() * &self.pz
    }

    pub fn on_shell_energy(&self, mass: Option<T>) -> Energy<T>
    where
        T: Mul<T, Output = T> + Add<T, Output = T> + Clone + std::ops::Add<Output = T> + Real,
    {
        let energy_squared = self.on_shell_energy_squared(mass);
        Energy {
            value: energy_squared.value.sqrt(),
        }
    }

    pub fn on_shell_energy_squared(&self, mass: Option<T>) -> Energy<T>
    where
        T: for<'a> Mul<&'a T, Output = T>
            + Add<T, Output = T>
            + Clone
            + std::ops::Add<Output = T>
            + Display,
    {
        let p2 = self.norm_squared();
        if let Some(mass) = mass {
            // println!("mass: {}", mass);
            Energy {
                value: p2 + mass.clone() * &mass,
            }
        } else {
            Energy { value: p2 }
        }
    }

    pub fn norm(&self) -> T
    where
        T: for<'a> Mul<&'a T, Output = T> + Add<T> + Real,
    {
        self.norm_squared().sqrt()
    }

    pub fn into_four_momentum_parametric(self, id: usize) -> FourMomentum<T, Atom> {
        let energy = Energy::new_parametric(id);
        FourMomentum {
            temporal: energy,
            spatial: self,
        }
    }

    pub fn into_on_shell_four_momentum(self, mass: Option<T>) -> FourMomentum<T, T>
    where
        T: Mul<T> + Add<T> + std::ops::Add<Output = T> + Real,
    {
        FourMomentum::new_on_shell(self, mass)
    }

    pub fn cast<U>(&self) -> ThreeMomentum<U>
    where
        T: Clone + Into<U>,
    {
        ThreeMomentum {
            px: (self.px.clone()).into(),
            py: (self.py.clone()).into(),
            pz: (self.pz.clone()).into(),
        }
    }

    pub fn into_f64(&self) -> ThreeMomentum<f64>
    where
        T: FloatLike,
    {
        ThreeMomentum {
            px: self.px.to_f64(),
            py: self.py.to_f64(),
            pz: self.pz.to_f64(),
        }
    }
}

impl<T> From<[T; 3]> for ThreeMomentum<T> {
    fn from(data: [T; 3]) -> Self {
        let [px, py, pz] = data;
        ThreeMomentum { px, py, pz }
    }
}

impl<T> From<ThreeMomentum<T>> for [T; 3] {
    fn from(data: ThreeMomentum<T>) -> Self {
        [data.px, data.py, data.pz]
    }
}

impl<T> From<(T, T, T)> for ThreeMomentum<T> {
    fn from(data: (T, T, T)) -> Self {
        let (px, py, pz) = data;
        ThreeMomentum { px, py, pz }
    }
}

impl<T> From<ThreeMomentum<T>> for (T, T, T) {
    fn from(data: ThreeMomentum<T>) -> Self {
        (data.px, data.py, data.pz)
    }
}

#[derive(Default, Debug, PartialEq, Eq, Clone, Copy, Serialize, Deserialize)]
pub struct FourMomentum<T, U = T> {
    pub temporal: Energy<U>,
    pub spatial: ThreeMomentum<T>,
}

impl<T> FourMomentum<T> {
    pub fn map_ref<U>(&self, f: &impl Fn(&T) -> U) -> FourMomentum<U> {
        FourMomentum {
            temporal: self.temporal.map_ref(f),
            spatial: self.spatial.map_ref(f),
        }
    }

    pub fn map<U>(self, f: &impl Fn(T) -> U) -> FourMomentum<U> {
        FourMomentum {
            temporal: self.temporal.map(f),
            spatial: self.spatial.map(f),
        }
    }
}

impl<T: FloatLike> ApproxEq<FourMomentum<F<T>>, F<T>> for FourMomentum<F<T>> {
    fn approx_eq(&self, other: &FourMomentum<F<T>>, threshold: &F<T>) -> bool {
        self.temporal.approx_eq(&other.temporal, threshold)
            && self.spatial.approx_eq(&other.spatial, threshold)
    }
}

impl<T: FloatLike> ApproxEq<Polarization<Complex<F<T>>>, F<T>> for FourMomentum<F<T>> {
    fn approx_eq(&self, other: &Polarization<Complex<F<T>>>, tolerance: &F<T>) -> bool {
        if other.tensor.size().unwrap() != 4 {
            false
        } else {
            self.into_iter()
                .zip(other.tensor.iter_flat())
                .all(|(a, (_, b))| a.approx_eq(b, tolerance))
        }
    }

    fn approx_eq_res(
        &self,
        other: &Polarization<Complex<F<T>>>,
        tolerance: &F<T>,
    ) -> color_eyre::Result<()> {
        if other.tensor.size().unwrap() != 4 {
            Err(eyre::eyre!("Polarization tensor has wrong size."))
        } else {
            self.into_iter()
                .zip(other.tensor.iter_flat())
                .try_for_each(|(a, (i, b))| {
                    a.approx_eq_res(b, tolerance).wrap_err(format!(
                        "Polarization tensor element {} does not match. ",
                        i
                    ))
                })
        }
    }
}

impl<T: FloatLike> FourMomentum<F<T>> {
    pub fn from_ff64(four_momentum: &FourMomentum<F<f64>>) -> Self {
        let temporal = Energy::from_ff64(four_momentum.temporal);
        let spatial = ThreeMomentum::from_ff64(four_momentum.spatial);
        FourMomentum { temporal, spatial }
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy, Serialize, Deserialize)]
#[serde(untagged)]
pub enum ExternalMomenta<T> {
    Dependent(Dep),
    Independent([T; 4]),
}

impl<T: FloatLike> Rotatable for ExternalMomenta<F<T>> {
    fn rotate(&self, rotation: &Rotation) -> Self {
        match self {
            ExternalMomenta::Dependent(_) => ExternalMomenta::Dependent(Dep::Dep),
            ExternalMomenta::Independent(data) => {
                let fm = FourMomentum::from(data.clone());
                fm.rotate(rotation).into()
            }
        }
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy, Serialize, Deserialize)]
pub enum Dep {
    #[serde(rename = "dependent")]
    Dep,
}

#[derive(Error, Debug)]
pub enum ExternalMomentaError {
    #[error("Dependent momenta cannot be converted to FourMomentum")]
    Dependent,
}

impl<T> TryFrom<ExternalMomenta<T>> for FourMomentum<T> {
    type Error = ExternalMomentaError;
    fn try_from(value: ExternalMomenta<T>) -> Result<Self, Self::Error> {
        match value {
            ExternalMomenta::Dependent(_) => Err(ExternalMomentaError::Dependent),
            ExternalMomenta::Independent(data) => Ok(FourMomentum::from(data)),
        }
    }
}

impl<T> From<FourMomentum<T>> for ExternalMomenta<T> {
    fn from(value: FourMomentum<T>) -> Self {
        ExternalMomenta::Independent(value.into())
    }
}

impl<T> From<[T; 4]> for ExternalMomenta<T> {
    fn from(data: [T; 4]) -> Self {
        ExternalMomenta::Independent(data)
    }
}

impl<T: RefZero, U: RefZero> RefZero for FourMomentum<T, U> {
    fn ref_zero(&self) -> Self {
        FourMomentum {
            temporal: self.temporal.ref_zero(),
            spatial: self.spatial.ref_zero(),
        }
    }
}

impl<T: RefZero, U: RefZero> RefZero<FourMomentum<T, U>> for &FourMomentum<T, U> {
    fn ref_zero(&self) -> FourMomentum<T, U> {
        FourMomentum {
            temporal: self.temporal.ref_zero(),
            spatial: self.spatial.ref_zero(),
        }
    }
}

#[derive(Debug, PartialEq, Clone, Copy, Serialize, Deserialize)]
pub enum PolType {
    U,
    V,
    UBar,
    VBar,
    Scalar,
    Epsilon,
    EpsilonBar,
}

impl PolType {
    pub fn bar(self) -> Self {
        match self {
            Self::Epsilon => Self::EpsilonBar,
            Self::Scalar => Self::Scalar,
            Self::U => Self::UBar,
            Self::UBar => Self::U,
            Self::EpsilonBar => Self::Epsilon,
            Self::VBar => Self::V,
            Self::V => Self::VBar,
        }
    }
}

impl Display for PolType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::U => write!(f, "u"),
            Self::V => write!(f, "v"),
            Self::VBar => write!(f, "vbar"),
            Self::UBar => write!(f, "ubar"),
            Self::Scalar => write!(f, ""),
            Self::Epsilon => write!(f, "ϵ"),
            Self::EpsilonBar => write!(f, "ϵbar"),
        }
    }
}

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct Polarization<T> {
    pub tensor: DenseTensor<T, IndexLess<PhysReps>>,
    pol_type: PolType,
}

impl<T: FloatLike> ApproxEq<Polarization<Complex<F<T>>>, F<T>> for Polarization<Complex<F<T>>> {
    fn approx_eq(&self, other: &Polarization<Complex<F<T>>>, tolerance: &F<T>) -> bool {
        self.pol_type == other.pol_type
            && self
                .tensor
                .flat_iter()
                .zip(other.tensor.flat_iter())
                .all(|((_, a), (_, b))| a.approx_eq(b, tolerance))
    }

    fn approx_eq_res(
        &self,
        other: &Polarization<Complex<F<T>>>,
        tolerance: &F<T>,
    ) -> color_eyre::Result<()> {
        if self.pol_type != other.pol_type {
            Err(eyre::eyre!(
                "Polarization types do not match. self: {}, other: {}",
                self,
                other
            ))
        } else {
            self.tensor
                .flat_iter()
                .zip(other.tensor.flat_iter())
                .try_for_each(|((i, a), (_, b))| {
                    a.approx_eq_res(b, tolerance).wrap_err(format!(
                        "Polarization tensor element {} does not match. self: {}, other: {}",
                        i, self, other
                    ))
                })
        }
    }
}

impl<T: for<'a> std::ops::AddAssign<&'a T>> AddAssign<Polarization<T>> for Polarization<T> {
    fn add_assign(&mut self, rhs: Polarization<T>) {
        self.tensor += rhs.tensor;
    }
}

impl<T: for<'a> SubAssign<&'a T>> SubAssign<Polarization<T>> for Polarization<T> {
    fn sub_assign(&mut self, rhs: Polarization<T>) {
        self.tensor -= rhs.tensor;
    }
}

impl<T: FloatLike> Display for Polarization<F<T>> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Pol {}: {}", self.pol_type, self.tensor)
    }
}

impl<T: FloatLike> Display for Polarization<Complex<F<T>>> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Pol {}: {}", self.pol_type, self.tensor)
    }
}

impl<T: FloatLike> LowerExp for Polarization<Complex<F<T>>> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Pol {}: {:+e}", self.pol_type, self.tensor)
    }
}

impl<T> Polarization<T> {
    pub fn map<U>(&self, f: impl Fn(&T) -> U) -> Polarization<U> {
        Polarization {
            tensor: self.tensor.map_data_ref(f),
            pol_type: self.pol_type,
        }
    }
}

impl<T: Clone> Polarization<T> {
    pub fn is_scalar(&self) -> bool {
        self.pol_type == PolType::Scalar
    }
    pub fn scalar(value: T) -> Self {
        let structure = IndexLess::new(vec![]);
        Polarization {
            tensor: DenseTensor {
                data: vec![value],
                structure,
            },
            pol_type: PolType::Scalar,
        }
    }

    pub fn lorentz(value: [T; 4]) -> Self {
        let structure = IndexLess::new(vec![Minkowski::rep(4).cast()]);
        Polarization {
            tensor: DenseTensor {
                data: value.to_vec(),
                structure,
            },
            pol_type: PolType::Epsilon,
        }
    }

    pub fn bispinor_u(value: [T; 4]) -> Self {
        let structure = IndexLess::new(vec![Bispinor::rep(4).cast()]);

        Polarization {
            tensor: DenseTensor {
                data: value.to_vec(),
                structure,
            },
            pol_type: PolType::U,
        }
    }

    pub fn bispinor_v(value: [T; 4]) -> Self {
        let structure = IndexLess::new(vec![Bispinor::rep(4).cast()]);

        Polarization {
            tensor: DenseTensor {
                data: value.to_vec(),
                structure,
            },
            pol_type: PolType::V,
        }
    }

    pub fn shadow(&self) -> DenseTensor<Atom, IndexLess<PhysReps>> {
        self.tensor
            .structure
            .clone()
            .to_dense_labeled(|_, i| FlatCoefficent::<NoArgs> {
                index: i,
                name: Some(Symbol::new(self.pol_type.to_string())),
                args: None,
            })
            .unwrap()
    }
}

impl<T> Index<usize> for Polarization<T> {
    type Output = T;
    fn index(&self, index: usize) -> &Self::Output {
        &self.tensor[index.into()]
    }
}

impl<T> RefZero for Polarization<T>
where
    T: RefZero + Clone,
{
    fn ref_zero(&self) -> Self {
        Polarization {
            tensor: self.tensor.ref_zero(),
            pol_type: self.pol_type,
        }
    }
}

impl<'a, T, U> Mul<&'a U> for Polarization<T>
where
    T: Clone + MulAssign<&'a U>,
{
    type Output = Polarization<T>;
    fn mul(mut self, rhs: &'a U) -> Self::Output {
        for i in self.tensor.data.iter_mut() {
            *i *= rhs;
        }
        self
    }
}

impl<'a, T> MulAssign<&'a T> for Polarization<T>
where
    T: Clone + MulAssign<&'a T>,
{
    fn mul_assign(&mut self, rhs: &'a T) {
        for i in self.tensor.data.iter_mut() {
            *i *= rhs;
        }
    }
}

impl<T, U, Out> FallibleAdd<Polarization<U>> for Polarization<T>
where
    T: FallibleAdd<U, Output = Out>,
{
    type Output = Polarization<Out>;
    fn add_fallible(&self, rhs: &Polarization<U>) -> Option<Self::Output> {
        Some(Polarization {
            tensor: self.tensor.add_fallible(&rhs.tensor)?,
            pol_type: self.pol_type,
        })
    }
}

impl<T> Polarization<T> {
    pub fn cast<U>(&self) -> Polarization<U>
    where
        T: Clone,
        U: Clone + From<T>,
    {
        Polarization {
            tensor: self.tensor.cast(),
            pol_type: self.pol_type,
        }
    }
}

impl<T> Polarization<Complex<T>> {
    pub fn complex_cast<U>(&self) -> Polarization<Complex<U>>
    where
        T: Clone,
        U: Clone + From<T>,
    {
        Polarization {
            tensor: self
                .tensor
                .map_data_ref(|d| d.map_ref(|r| r.clone().into())),
            pol_type: self.pol_type,
        }
    }
}

impl<T: FloatLike> Polarization<Complex<F<T>>> {
    pub fn higher(&self) -> Polarization<Complex<F<T::Higher>>>
    where
        T::Higher: FloatLike,
    {
        Polarization {
            tensor: self.tensor.map_data_ref(|t| t.map_ref(|r| r.higher())),
            pol_type: self.pol_type,
        }
    }

    pub fn lower(&self) -> Polarization<Complex<F<T::Lower>>>
    where
        T::Lower: FloatLike,
    {
        Polarization {
            tensor: self.tensor.map_data_ref(|t| t.map_ref(|r| r.lower())),
            pol_type: self.pol_type,
        }
    }
}

impl<T: FloatLike, U: FloatLike> From<FourMomentum<T, U>> for FourMomentum<F<T>, F<U>> {
    fn from(value: FourMomentum<T, U>) -> Self {
        FourMomentum {
            temporal: value.temporal.into(),
            spatial: value.spatial.into(),
        }
    }
}

impl<T: FloatLike, U: FloatLike> FourMomentum<F<T>, F<U>> {
    pub fn higher(&self) -> FourMomentum<F<T::Higher>, F<U::Higher>>
    where
        T::Higher: FloatLike,
        U::Higher: FloatLike,
    {
        FourMomentum {
            temporal: self.temporal.higher(),
            spatial: self.spatial.higher(),
        }
    }

    pub fn lower(&self) -> FourMomentum<F<T::Lower>, F<U::Lower>>
    where
        T::Lower: FloatLike,
        U::Lower: FloatLike,
    {
        FourMomentum {
            temporal: self.temporal.lower(),
            spatial: self.spatial.lower(),
        }
    }
}

impl<T, U> FourMomentum<T, U> {
    pub fn new(energy: Energy<U>, three_momentum: ThreeMomentum<T>) -> Self {
        FourMomentum {
            temporal: energy,
            spatial: three_momentum,
        }
    }
}

impl<T: Real> RefDefault for FourMomentum<T, T> {
    fn default(&self) -> Self {
        let zero = self.temporal.value.zero();
        FourMomentum {
            temporal: Energy::new(zero.clone()),
            spatial: ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone()),
        }
    }
}

impl<T> Mul<T> for &FourMomentum<T, T>
where
    T: Mul<T, Output = T> + Clone,
{
    type Output = FourMomentum<T, T>;
    fn mul(self, rhs: T) -> Self::Output {
        FourMomentum {
            temporal: Energy {
                value: self.temporal.value.clone() * rhs.clone(),
            },
            spatial: self.spatial.clone() * rhs,
        }
    }
}

impl<T> Mul<&T> for &FourMomentum<T, T>
where
    T: for<'b> Mul<&'b T, Output = T> + Clone,
{
    type Output = FourMomentum<T, T>;
    fn mul(self, rhs: &T) -> Self::Output {
        FourMomentum {
            temporal: Energy {
                value: self.temporal.value.clone() * rhs,
            },
            spatial: self.spatial.clone() * rhs,
        }
    }
}

impl<T> FourMomentum<T, T> {
    pub fn zero(&self) -> Self
    where
        T: Real,
    {
        let zero = self.temporal.value.zero();
        FourMomentum {
            temporal: Energy::new(zero.clone()),
            spatial: ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone()),
        }
    }

    pub fn square(&self) -> T
    where
        T: for<'a> Mul<&'a T, Output = T> + Add<T, Output = T> + Clone + Sub<T, Output = T>,
    {
        let temporal = self.temporal.value.clone();
        let spatial = self.spatial.norm_squared();
        temporal * &self.temporal.value - spatial
    }

    pub fn norm(&self) -> T
    where
        T: Real,
    {
        self.square().sqrt()
    }

    pub fn from_args(energy: T, px: T, py: T, pz: T) -> Self {
        let energy = Energy::new(energy);
        let three_momentum = ThreeMomentum { px, py, pz };
        FourMomentum {
            temporal: energy,
            spatial: three_momentum,
        }
    }
    pub fn new_on_shell(three_momentum: ThreeMomentum<T>, mass: Option<T>) -> Self
    where
        T: Mul<T> + Add<T> + std::ops::Add<Output = T> + Real,
    {
        let energy = three_momentum.on_shell_energy(mass);
        // println!("{}", energy);
        FourMomentum {
            temporal: energy,
            spatial: three_momentum,
        }
    }

    pub fn into_dense(self, index: AbstractIndex) -> DenseTensor<T, VecStructure>
    where
        T: Clone,
    {
        let structure =
            VecStructure::from_iter([PhysReps::new_slot(Minkowski {}.into(), 4, index)]);
        DenseTensor::from_data(
            vec![
                self.temporal.value,
                self.spatial.px,
                self.spatial.py,
                self.spatial.pz,
            ],
            structure,
        )
        .unwrap()
    }

    pub fn into_dense_named(
        self,
        index: AbstractIndex,
        name: Symbol,
        num: usize,
    ) -> DenseTensor<T, NamedStructure<Symbol, usize>>
    where
        T: Clone,
    {
        let structure =
            VecStructure::from_iter([PhysReps::new_slot(Minkowski {}.into(), 4, index)])
                .to_named(name, Some(num));
        DenseTensor::from_data(
            vec![
                self.temporal.value,
                self.spatial.px,
                self.spatial.py,
                self.spatial.pz,
            ],
            structure,
        )
        .unwrap()
    }

    pub fn cast<U>(&self) -> FourMomentum<U, U>
    where
        T: Clone + Into<U>,
    {
        FourMomentum {
            temporal: Energy::new((self.temporal.value.clone()).into()),
            spatial: ThreeMomentum::cast(&self.spatial),
        }
    }

    pub fn boost(&self, boost_vector: &FourMomentum<T>) -> FourMomentum<T>
    where
        T: Real + SingleFloat + PartialOrd,
    {
        let b2 = boost_vector.spatial.norm_squared();
        let one = b2.one();
        let zero = one.zero();
        let gamma = (one.clone() - &b2).sqrt().inv();

        let bp = self.spatial.clone() * &boost_vector.spatial;
        let gamma2 = if b2 > zero {
            (gamma.clone() - &one) / b2
        } else {
            zero
        };
        let factor = gamma2 * &bp + gamma.clone() * &self.temporal.value;

        FourMomentum::from_args(
            (bp + &self.temporal.value) * &gamma,
            self.spatial.px.mul_add(&factor, &boost_vector.spatial.px),
            self.spatial.py.mul_add(&factor, &boost_vector.spatial.py),
            self.spatial.pz.mul_add(&factor, &boost_vector.spatial.pz),
        )
    }
}

impl<T: FloatLike> FourMomentum<F<T>, F<T>> {
    /// Compute the phi-angle separation with p2.
    pub fn getdelphi(&self, p2: &FourMomentum<F<T>>) -> F<T>
    where
        T: Real,
    {
        self.spatial.getdelphi(&p2.spatial)
    }

    /// Compute the deltaR separation with momentum p2.
    #[inline]
    pub fn delta_r(&self, p2: &FourMomentum<F<T>>) -> F<T>
    where
        T: Real,
    {
        self.spatial.delta_r(&p2.spatial)
    }

    pub fn pt(&self) -> F<T>
    where
        T: Real,
    {
        self.spatial.pt()
    }

    pub fn to_f64(&self) -> FourMomentum<F<f64>, F<f64>> {
        FourMomentum {
            temporal: Energy {
                value: F(self.temporal.value.to_f64()),
            },
            spatial: self.spatial.to_f64().cast(),
        }
    }

    pub fn pol_one(&self) -> [F<T>; 4]
    where
        T: FloatLike,
    {
        // definition from helas_ref A.2

        // debug!("pol_one in: {}", self);

        let pt = self.pt();
        let p = self.spatial.norm();

        let (e1, e2, e3) = if pt.is_zero() {
            (pt.one(), pt.zero(), pt.zero())
        } else {
            (
                &self.spatial.px * &self.spatial.pz / (&pt * &p),
                &self.spatial.py * &self.spatial.pz / (&pt * &p),
                -(&pt / &p),
            )
        };

        // debug!(
        //     " (pt.zero(), e1, e2, e3) {} {} {} {}",
        //     pt.zero(),
        //     e1,
        //     e2,
        //     e3
        // );
        [pt.zero(), e1, e2, e3]

        // debug!("pol :{pol}");
    }

    pub fn pol_two(&self) -> [F<T>; 4]
    where
        T: FloatLike,
    {
        // definition from helas_ref A.2
        let pt = self.pt();
        let (e1, e2, e3) = if pt.is_zero() {
            if self.spatial.pz.positive() {
                (pt.zero(), pt.one(), pt.zero())
            } else {
                (pt.zero(), -pt.one(), pt.zero())
            }
        } else {
            (-(&self.spatial.py / &pt), &self.spatial.px / &pt, pt.zero())
        };
        [pt.zero(), e1, e2, e3]
    }

    pub fn pol_three(&self) -> Polarization<F<T>>
    where
        T: FloatLike,
    {
        // definition from helas_ref A.2
        let m = self.norm();
        let p = self.spatial.norm();
        let emp = &self.temporal.value / (&m * &p);
        let e0 = p.square() / &self.temporal.value;
        let e1 = &self.spatial.px / &emp;
        let e2 = &self.spatial.py / &emp;
        let e3 = &self.spatial.pz / &emp;

        Polarization::lorentz([e0, e1, e2, e3])
    }

    pub fn pol(&self, lambda: Helicity) -> Polarization<Complex<F<T>>> {
        if lambda.is_zero() {
            self.pol_three().cast()
        } else {
            let one = self.temporal.value.one();
            let sqrt_2_inv: F<T> = (&one + &one).sqrt().inv();

            let [eone0, eone1, eone2, eone3] = self.pol_one();

            let [etwo0, etwo1, etwo2, etwo3] = self.pol_two();

            Polarization::lorentz([
                Complex {
                    re: -lambda * eone0 * &sqrt_2_inv,
                    im: -etwo0 * &sqrt_2_inv,
                },
                Complex {
                    re: -lambda * eone1 * &sqrt_2_inv,
                    im: -etwo1 * &sqrt_2_inv,
                },
                Complex {
                    re: -lambda * eone2 * &sqrt_2_inv,
                    im: -etwo2 * &sqrt_2_inv,
                },
                Complex {
                    re: -lambda * eone3 * &sqrt_2_inv,
                    im: -etwo3 * &sqrt_2_inv,
                },
            ])
        }
    }

    pub fn omega(&self, lambda: Sign) -> Complex<F<T>> {
        match lambda {
            Sign::Positive => (&self.temporal.value + self.spatial.norm()).complex_sqrt(),
            Sign::Negative => (&self.temporal.value - self.spatial.norm()).complex_sqrt(),
        }
    }

    pub fn u(&self, lambda: Sign) -> Polarization<Complex<F<T>>> {
        let xi = self.xi(lambda);
        Polarization::bispinor_u([
            self.omega(-lambda) * &xi[0],
            self.omega(-lambda) * &xi[1],
            self.omega(lambda) * &xi[0],
            self.omega(lambda) * &xi[1],
        ])
    }

    pub fn v(&self, lambda: Sign) -> Polarization<Complex<F<T>>> {
        let xi = self.xi(-lambda);
        Polarization::bispinor_v([
            (-lambda) * self.omega(lambda) * &xi[0],
            (-lambda) * self.omega(lambda) * &xi[1],
            lambda * self.omega(-lambda) * &xi[0],
            lambda * self.omega(-lambda) * &xi[1],
        ])
    }

    pub fn xi(&self, lambda: Sign) -> [Complex<F<T>>; 2] {
        if self.spatial.pz == -self.spatial.norm() {
            let zero: Complex<F<T>> = self.temporal.value.zero().into();
            let one = zero.one();
            match lambda {
                Sign::Positive => [zero, one],
                Sign::Negative => [-one, zero],
            }
        } else {
            let prefactor: F<T> = ((F::from_f64(2.)
                * self.spatial.norm()
                * (self.spatial.norm() + &self.spatial.pz))
                .sqrt())
            .inv();
            let mut xi: [Complex<F<T>>; 2] = [
                Complex::new_re(&prefactor * (self.spatial.norm() + &self.spatial.pz)),
                Complex::new(self.spatial.px.clone(), self.spatial.py.clone()) * &prefactor,
            ]; //plus

            if matches!(lambda, Sign::Negative) {
                xi.swap(0, 1);
                xi[0].re = -xi[0].re.clone();
            }
            xi
        }
    }
}

impl<T: FloatLike> Polarization<Complex<F<T>>> {
    pub fn bar(&self) -> Self {
        let mut tensor = self.tensor.map_data_ref(Complex::conj);

        if matches!(
            self.pol_type,
            PolType::U | PolType::V | PolType::UBar | PolType::VBar
        ) {
            tensor.data.swap(0, 2);
            tensor.data.swap(1, 3);
        }
        Polarization {
            tensor,
            pol_type: self.pol_type.bar(),
        }
    }
}

impl<T> IntoIterator for Polarization<T> {
    type IntoIter = std::vec::IntoIter<Self::Item>;
    type Item = T;

    fn into_iter(self) -> Self::IntoIter {
        self.tensor.data.into_iter()
    }
}

impl<'a, T> IntoIterator for &'a Polarization<T> {
    type IntoIter = std::slice::Iter<'a, T>;
    type Item = &'a T;

    fn into_iter(self) -> Self::IntoIter {
        self.tensor.data.iter()
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy, Serialize, Deserialize, Hash)]
#[repr(i8)]
pub enum Sign {
    Positive = 1,
    Negative = -1,
}

pub trait Pow<E> {
    fn pow(&self, exponent: E) -> Self;
}

macro_rules! impl_pow_for_sign {
    ($type:ty, $exponent:ty,$pos:expr) => {
        impl Pow<$exponent> for $type {
            fn pow(&self, exponent: $exponent) -> Self {
                match exponent {
                    0 => $pos,
                    a => {
                        if a % 2 == 0 {
                            $pos
                        } else {
                            *self
                        }
                    }
                }
            }
        }
    };
}

impl_pow_for_sign!(Sign, i32, Sign::Positive);
impl_pow_for_sign!(Sign, u32, Sign::Positive);
impl_pow_for_sign!(Sign, i64, Sign::Positive);
impl_pow_for_sign!(Sign, u64, Sign::Positive);
impl_pow_for_sign!(Sign, isize, Sign::Positive);
impl_pow_for_sign!(Sign, usize, Sign::Positive);
impl_pow_for_sign!(Sign, i128, Sign::Positive);
impl_pow_for_sign!(Sign, u128, Sign::Positive);

use thiserror::Error;

#[derive(Error, Debug)]
pub enum SignError {
    #[error("Invalid value for Sign")]
    InvalidValue,
    #[error("Zero is not a valid value for Sign")]
    ZeroValue,
}

impl TryFrom<SignOrZero> for Sign {
    type Error = SignError;
    fn try_from(value: SignOrZero) -> Result<Self, Self::Error> {
        match value {
            SignOrZero::Zero => Err(SignError::ZeroValue),
            SignOrZero::Plus => Ok(Sign::Positive),
            SignOrZero::Minus => Ok(Sign::Negative),
        }
    }
}
// impl Serialize for Sign {
//     fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
//     where
//         S: Serializer,
//     {
//         let value = match *self {
//             Sign::Positive => 1,
//             Sign::Negative => -1,
//         };
//         serializer.serialize_i32(value)
//     }
// }

// impl<'de> Deserialize<'de> for Sign {
//     fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
//     where
//         D: Deserializer<'de>,
//     {
//         let value = i32::deserialize(deserializer)?;
//         match value {
//             1 => Ok(Sign::Positive),
//             -1 => Ok(Sign::Negative),
//             _ => Err(de::Error::custom(
//                 "Expected 1 for Positive or -1 for Negative",
//             )),
//         }
//     }
// }

impl<T: Neg<Output = T>> Mul<T> for Sign {
    type Output = T;
    fn mul(self, rhs: T) -> Self::Output {
        match self {
            Sign::Positive => rhs,
            Sign::Negative => -rhs,
        }
    }
}

impl Neg for Sign {
    type Output = Self;
    fn neg(self) -> Self::Output {
        match self {
            Sign::Positive => Sign::Negative,
            Sign::Negative => Sign::Positive,
        }
    }
}

#[derive(
    Debug,
    PartialEq,
    Eq,
    Clone,
    Copy,
    Serialize_repr,
    Deserialize_repr,
    Encode,
    Decode,
    PartialOrd,
    Ord,
    Hash,
)]
#[repr(i8)]
pub enum SignOrZero {
    Zero = 0,
    Plus = 1,
    Minus = -1,
}

impl_pow_for_sign!(SignOrZero, i32, SignOrZero::Plus);
impl_pow_for_sign!(SignOrZero, u32, SignOrZero::Plus);
impl_pow_for_sign!(SignOrZero, i64, SignOrZero::Plus);
impl_pow_for_sign!(SignOrZero, u64, SignOrZero::Plus);
impl_pow_for_sign!(SignOrZero, isize, SignOrZero::Plus);
impl_pow_for_sign!(SignOrZero, usize, SignOrZero::Plus);
impl_pow_for_sign!(SignOrZero, i128, SignOrZero::Plus);
impl_pow_for_sign!(SignOrZero, u128, SignOrZero::Plus);

impl From<linnet::half_edge::involution::Orientation> for SignOrZero {
    fn from(value: linnet::half_edge::involution::Orientation) -> Self {
        match value {
            linnet::half_edge::involution::Orientation::Default => SignOrZero::Plus,
            linnet::half_edge::involution::Orientation::Reversed => SignOrZero::Minus,
            linnet::half_edge::involution::Orientation::Undirected => SignOrZero::Zero,
        }
    }
}

impl From<linnet::half_edge::involution::Flow> for SignOrZero {
    fn from(value: linnet::half_edge::involution::Flow) -> Self {
        match value {
            linnet::half_edge::involution::Flow::Source => SignOrZero::Plus,
            linnet::half_edge::involution::Flow::Sink => SignOrZero::Minus,
        }
    }
}

impl TryFrom<i8> for SignOrZero {
    type Error = SignError;
    fn try_from(value: i8) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(SignOrZero::Zero),
            1 => Ok(SignOrZero::Plus),
            -1 => Ok(SignOrZero::Minus),
            _ => Err(SignError::InvalidValue),
        }
    }
}

impl Display for SignOrZero {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SignOrZero::Zero => write!(f, "."),
            SignOrZero::Plus => write!(f, "+"),
            SignOrZero::Minus => write!(f, "-"),
        }
    }
}

pub type Helicity = SignOrZero;

#[derive(
    Debug, PartialEq, Eq, Clone, Serialize, Deserialize, Encode, Decode, PartialOrd, Ord, Hash,
)]
pub struct Signature(pub Vec<SignOrZero>);

impl Display for Signature {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for sign in &self.0 {
            write!(f, "{}", sign)?;
        }
        Ok(())
    }
}

impl FromIterator<SignOrZero> for Signature {
    fn from_iter<I: IntoIterator<Item = SignOrZero>>(iter: I) -> Self {
        Signature(iter.into_iter().collect())
    }
}

impl FromIterator<i8> for Signature {
    fn from_iter<I: IntoIterator<Item = i8>>(iter: I) -> Self {
        Signature(
            iter.into_iter()
                .map(|x| match x {
                    0 => SignOrZero::Zero,
                    1 => SignOrZero::Plus,
                    -1 => SignOrZero::Minus,
                    _ => panic!("Invalid value for Signature"),
                })
                .collect(),
        )
    }
}

impl FromIterator<isize> for Signature {
    fn from_iter<I: IntoIterator<Item = isize>>(iter: I) -> Self {
        Signature(
            iter.into_iter()
                .map(|x| match x {
                    0 => SignOrZero::Zero,
                    1 => SignOrZero::Plus,
                    -1 => SignOrZero::Minus,
                    _ => panic!("Invalid value for Signature"),
                })
                .collect(),
        )
    }
}

impl From<Vec<i8>> for Signature {
    fn from(value: Vec<i8>) -> Self {
        Signature::from_iter(value)
    }
}

impl IntoIterator for Signature {
    type Item = SignOrZero;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a> IntoIterator for &'a Signature {
    type Item = SignOrZero;
    type IntoIter = std::iter::Copied<std::slice::Iter<'a, Self::Item>>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter().copied()
    }
}

impl Index<usize> for Signature {
    type Output = SignOrZero;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl Signature {
    pub fn validate_basis<T>(&self, basis: &[T]) -> bool {
        self.len() == basis.len()
    }

    pub fn panic_validate_basis<T>(&self, basis: &[T]) {
        if !self.validate_basis(basis) {
            panic!(
                "Invalid basis for Signature, expected length {}, got length {}",
                self.len(),
                basis.len()
            );
        }
    }

    pub fn to_momtrop_format(&self) -> Vec<isize> {
        self.0
            .iter()
            .map(|x| match x {
                SignOrZero::Zero => 0,
                SignOrZero::Plus => 1,
                SignOrZero::Minus => -1,
            })
            .collect()
    }
    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn iter(&self) -> std::slice::Iter<SignOrZero> {
        self.0.iter()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    pub fn apply<T>(&self, basis: &[T]) -> T
    where
        T: RefZero + Clone + Neg<Output = T> + AddAssign<T>,
    {
        // self.panic_validate_basis(basis);
        let mut result = basis[0].ref_zero();
        for (&sign, t) in self.0.iter().zip(basis.iter().cloned()) {
            result += sign * t;
        }
        result
    }

    pub fn apply_iter<I, T>(&self, basis: I) -> Option<T>
    where
        I: IntoIterator,
        I::Item: RefZero<T>,
        T: Clone + SubAssign<I::Item> + AddAssign<I::Item>,
    {
        let mut basis_iter = basis.into_iter();
        let mut signature_iter = self.into_iter();

        while let (Some(sign), Some(item)) = (signature_iter.next(), basis_iter.next()) {
            if sign.is_sign() {
                // Initialize the result based on the first non-zero sign
                let mut result = item.ref_zero();
                match sign {
                    SignOrZero::Zero => {
                        panic!("unreachable");
                        // return None;
                    }
                    SignOrZero::Plus => {
                        result += item;
                    }
                    SignOrZero::Minus => {
                        result -= item;
                    }
                }

                // Continue processing the rest of the iterator
                while let (Some(sign), Some(item)) = (signature_iter.next(), basis_iter.next()) {
                    match sign {
                        SignOrZero::Zero => {}
                        SignOrZero::Plus => {
                            result += item;
                        }
                        SignOrZero::Minus => {
                            result -= item;
                        }
                    }
                }

                return Some(result);
            }
        }

        // Return None if no non-zero sign was found
        None
    }

    pub fn label_with(&self, label: &str) -> String {
        let mut result = String::new();
        let mut first = true;
        for (i, sign) in self.0.iter().enumerate() {
            if !first {
                result.push_str(&sign.to_string());
            } else {
                first = false;
            }
            if sign.is_sign() {
                result.push_str(&format!("{}_{}", label, i));
            }
        }
        result
    }
}

#[test]
fn test_signature() {
    let sig = Signature(vec![SignOrZero::Plus, SignOrZero::Minus]);
    let basis: Vec<i32> = vec![1, 2];
    assert_eq!(sig.apply(&basis), 1 - 2);
    assert_eq!(sig.apply_iter(basis.iter()), Some(-1));

    let basis: [FourMomentum<i32>; 4] = [
        FourMomentum::from_args(1, 1, 0, 0),
        FourMomentum::from_args(1, 0, 1, 0),
        FourMomentum::from_args(1, 0, 0, 1),
        FourMomentum::from_args(1, 1, 1, 1),
    ];

    let sig = Signature(vec![
        SignOrZero::Plus,
        SignOrZero::Minus,
        SignOrZero::Zero,
        SignOrZero::Plus,
    ]);

    assert_eq!(sig.apply(&basis), FourMomentum::from_args(1, 2, 0, 1));
    let sig = Signature(vec![
        SignOrZero::Zero,
        SignOrZero::Zero,
        SignOrZero::Zero,
        SignOrZero::Zero,
    ]);
    assert_eq!(sig.apply_iter(basis.iter()), None);
    let sig = Signature(vec![
        SignOrZero::Zero,
        SignOrZero::Zero,
        SignOrZero::Zero,
        SignOrZero::Minus,
    ]);
    assert_eq!(sig.apply_iter(basis.iter()), Some(-basis[3]));
}

// impl Serialize for Helicity {
//     fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
//     where
//         S: Serializer,
//     {
//         match *self {
//             Helicity::Sign(ref sign) => sign.serialize(serializer),
//             Helicity::Zero => serializer.serialize_i32(0),
//         }
//     }
// }

// impl<'de> Deserialize<'de> for Helicity {
//     fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
//     where
//         D: Deserializer<'de>,
//     {
//         let value = i32::deserialize(deserializer)?;
//         match value {
//             1 => Ok(Helicity::Sign(Sign::Positive)),
//             -1 => Ok(Helicity::Sign(Sign::Negative)),
//             0 => Ok(Helicity::Zero),
//             _ => Err(de::Error::custom(
//                 "Expected 1 for Positive, -1 for Negative, or 0 for Zero",
//             )),
//         }
//     }
// }
#[allow(non_upper_case_globals)]
impl Helicity {
    pub fn is_zero(&self) -> bool {
        matches!(self, Helicity::Zero)
    }

    pub fn is_sign(&self) -> bool {
        matches!(self, Helicity::Plus | Helicity::Minus)
    }

    pub fn is_positive(&self) -> bool {
        matches!(self, Helicity::Plus)
    }

    pub fn is_negative(&self) -> bool {
        matches!(self, Helicity::Minus)
    }
}

impl<T: Neg<Output = T> + RefZero> Mul<T> for Helicity {
    type Output = T;
    fn mul(self, rhs: T) -> Self::Output {
        match self {
            Helicity::Plus => rhs,
            Helicity::Minus => -rhs,
            Helicity::Zero => rhs.ref_zero(),
        }
    }
}

impl Neg for Helicity {
    type Output = Self;
    fn neg(self) -> Self::Output {
        match self {
            Helicity::Plus => Helicity::Minus,
            Helicity::Minus => Helicity::Plus,
            Helicity::Zero => Helicity::Zero,
        }
    }
}

impl<T> IntoIterator for FourMomentum<T> {
    type Item = T;
    type IntoIter = std::array::IntoIter<T, 4>;

    fn into_iter(self) -> Self::IntoIter {
        let [e, px, py, pz] = [
            self.temporal.value,
            self.spatial.px,
            self.spatial.py,
            self.spatial.pz,
        ];
        [e, px, py, pz].into_iter()
    }
}

impl<'a, T> IntoIterator for &'a FourMomentum<T> {
    type Item = &'a T;
    type IntoIter = std::array::IntoIter<&'a T, 4>;

    fn into_iter(self) -> Self::IntoIter {
        let [e, px, py, pz] = [
            &self.temporal.value,
            &self.spatial.px,
            &self.spatial.py,
            &self.spatial.pz,
        ];
        [e, px, py, pz].into_iter()
    }
}

impl<T> From<[T; 4]> for FourMomentum<T, T> {
    fn from(data: [T; 4]) -> Self {
        let [t, px, py, pz] = data;
        FourMomentum {
            temporal: Energy::new(t),
            spatial: ThreeMomentum { px, py, pz },
        }
    }
}

impl<T> From<FourMomentum<T, T>> for [T; 4] {
    fn from(data: FourMomentum<T, T>) -> Self {
        [
            data.temporal.value,
            data.spatial.px,
            data.spatial.py,
            data.spatial.pz,
        ]
    }
}

impl<T> From<(T, T, T, T)> for FourMomentum<T, T> {
    fn from(data: (T, T, T, T)) -> Self {
        let (t, px, py, pz) = data;
        FourMomentum {
            temporal: Energy::new(t),
            spatial: ThreeMomentum { px, py, pz },
        }
    }
}

impl<T> From<FourMomentum<T, T>> for (T, T, T, T) {
    fn from(data: FourMomentum<T, T>) -> Self {
        (
            data.temporal.value,
            data.spatial.px,
            data.spatial.py,
            data.spatial.pz,
        )
    }
}

impl<T> FourMomentum<T, Atom> {
    pub fn into_dense_param(
        self,
        index: AbstractIndex,
    ) -> DenseTensor<MultivariatePolynomial<RationalField, T>, VecStructure>
    where
        T: Clone + Into<Coefficient> + Exponent,
    {
        let structure =
            VecStructure::from_iter([PhysReps::new_slot(Minkowski {}.into(), 4, index)]);
        let energy = self
            .temporal
            .value
            .to_polynomial(&RationalField::new(IntegerRing {}), None);

        let px: MultivariatePolynomial<RationalField, _> =
            Atom::new_num(self.spatial.px).to_polynomial(&RationalField::new(IntegerRing {}), None);
        let py =
            Atom::new_num(self.spatial.py).to_polynomial(&RationalField::new(IntegerRing {}), None);
        let pz =
            Atom::new_num(self.spatial.pz).to_polynomial(&RationalField::new(IntegerRing {}), None);

        DenseTensor::from_data(vec![energy, px, py, pz], structure).unwrap()
    }
}

impl<T, U> Add<FourMomentum<T, U>> for FourMomentum<T, U>
where
    T: Add<T, Output = T>,
    U: Add<U, Output = U>,
{
    type Output = FourMomentum<T, U>;
    fn add(self, rhs: FourMomentum<T, U>) -> Self::Output {
        FourMomentum {
            temporal: self.temporal + rhs.temporal,
            spatial: self.spatial + rhs.spatial,
        }
    }
}

impl<T, U> AddAssign<FourMomentum<T, U>> for FourMomentum<T, U>
where
    T: AddAssign<T>,
    U: AddAssign<U>,
{
    fn add_assign(&mut self, rhs: FourMomentum<T, U>) {
        self.temporal += rhs.temporal;
        self.spatial += rhs.spatial;
    }
}

impl<'a, T, U> AddAssign<&'a FourMomentum<T, U>> for FourMomentum<T, U>
where
    T: AddAssign<&'a T>,
    U: AddAssign<&'a U>,
{
    fn add_assign(&mut self, rhs: &'a FourMomentum<T, U>) {
        self.temporal += &rhs.temporal;
        self.spatial += &rhs.spatial;
    }
}

impl<T, U> Sub<FourMomentum<T, U>> for FourMomentum<T, U>
where
    T: Sub<T, Output = T>,
    U: Sub<U, Output = U>,
{
    type Output = FourMomentum<T, U>;
    fn sub(self, rhs: FourMomentum<T, U>) -> Self::Output {
        FourMomentum {
            temporal: self.temporal - rhs.temporal,
            spatial: self.spatial - rhs.spatial,
        }
    }
}

impl<T, U> Neg for FourMomentum<T, U>
where
    T: Neg<Output = T>,
    U: Neg<Output = U>,
{
    type Output = FourMomentum<T, U>;
    fn neg(self) -> Self::Output {
        FourMomentum {
            temporal: -self.temporal,
            spatial: -self.spatial,
        }
    }
}

impl<T, U> Neg for &FourMomentum<T, U>
where
    T: Neg<Output = T> + Clone,
    U: Neg<Output = U> + Clone,
{
    type Output = FourMomentum<T, U>;
    fn neg(self) -> Self::Output {
        FourMomentum {
            temporal: -self.temporal.clone(),
            spatial: -self.spatial.clone(),
        }
    }
}

impl<T, U> SubAssign<FourMomentum<T, U>> for FourMomentum<T, U>
where
    T: SubAssign<T>,
    U: SubAssign<U>,
{
    fn sub_assign(&mut self, rhs: FourMomentum<T, U>) {
        self.temporal -= rhs.temporal;
        self.spatial -= rhs.spatial;
    }
}

impl<'a, T, U> SubAssign<&'a FourMomentum<T, U>> for FourMomentum<T, U>
where
    T: SubAssign<&'a T>,
    U: SubAssign<&'a U>,
{
    fn sub_assign(&mut self, rhs: &'a FourMomentum<T, U>) {
        self.temporal -= &rhs.temporal;
        self.spatial -= &rhs.spatial;
    }
}

impl<T: Display, U: Display> Display for FourMomentum<T, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}, {}", self.temporal, self.spatial)
    }
}

impl<T: LowerExp, U: LowerExp> LowerExp for FourMomentum<T, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:e}, {:e}", self.temporal, self.spatial)
    }
}

impl From<Vector<F<f64>, 3>> for ThreeMomentum<F<f64>> {
    fn from(value: Vector<F<f64>, 3>) -> Self {
        ThreeMomentum::new(value[0], value[1], value[2])
    }
}

pub enum LorentzTransformation<T> {
    Boost(ThreeMomentum<T>),
    Rotation(Rotation),
    General {
        boost: ThreeMomentum<T>,
        rotation: Rotation,
    },
}

pub trait LorentzTransformable<T>: Rotatable {
    fn lorentz_transform(&self, transformation: &LorentzTransformation<T>) -> Self;
}

// pub struct IrriducibleLorentzRep {
//     m: usize,
//     n: usize,
// }

// impl IrriducibleLorentzRep {
//     pub fn scalar() -> Self {
//         Self { m: 0, n: 0 }
//     }

//     pub fn vector() -> Self {
//         Self { m: 1, n: 1 }
//     }

//     pub fn tensor() -> Self {
//         Self { m: 2, n: 2 }
//     }

//     pub fn left_weyl() -> Self {
//         Self { m: 1, n: 0 }
//     }

//     pub fn right_weyl() -> Self {
//         Self { m: 0, n: 1 }
//     }
// }

pub enum LorentzRep {
    // Irriducible(IrriducibleLorentzRep),
    // DirectSum(Vec<LorentzRep>),
    Scalar,
    Vector,
    Bispinor,
}

#[derive(Clone, Copy)]
pub enum RotationMethod {
    EulerAngles(f64, f64, f64),
    Pi2X,
    Pi2Y,
    Pi2Z,
    Identity,
}

impl From<RotationMethod> for Rotation {
    fn from(method: RotationMethod) -> Self {
        Rotation::new(method)
    }
}

#[derive(Clone)]
pub struct Rotation {
    pub method: RotationMethod,
    pub lorentz_rotation: EvalTensor<ExpressionEvaluator<Complex<F<f64>>>, VecStructure>,
    pub bispinor_rotation: EvalTensor<ExpressionEvaluator<Complex<F<f64>>>, VecStructure>,
    // phantom: PhantomData<T>,
}

impl Rotation {
    pub fn is_identity(&self) -> bool {
        matches!(self.method, RotationMethod::Identity)
    }
    pub fn setting(&self) -> RotationSetting {
        match self.method {
            RotationMethod::EulerAngles(alpha, beta, gamma) => {
                RotationSetting::EulerAngles { alpha, beta, gamma }
            }
            RotationMethod::Pi2X => RotationSetting::Pi2X,
            RotationMethod::Pi2Y => RotationSetting::Pi2Y,
            RotationMethod::Pi2Z => RotationSetting::Pi2Z,
            RotationMethod::Identity => RotationSetting::None,
        }
    }
    pub fn new(method: RotationMethod) -> Self {
        let mu = Minkowski::slot(4, 1);

        let al = Minkowski::slot(4, 3);
        let mud = mu.dual();

        let shadow: NamedStructure<String, ()> =
            VecStructure::<PhysReps>::from_iter([mu.into()]).to_named("eps".to_string(), None);
        let shadow_t: MixedTensor<_, VecStructure> =
            ParamOrConcrete::param(shadow.to_shell().expanded_shadow().unwrap().into())
                .cast_structure();

        let rotation: MixedTensor<f64, _> = method
            .lorentz_tensor(mud, al)
            .try_into_dense()
            .unwrap()
            .into();

        let fn_map = FunctionMap::new();

        let lorentz_eval: EvalTensor<ExpressionEvaluator<Rational>, VecStructure> = shadow_t
            .contract(&rotation)
            .unwrap()
            .try_into_parametric()
            .unwrap()
            .to_evaluation_tree(
                &fn_map,
                &shadow_t.try_into_parametric().unwrap().tensor.data(),
            )
            .unwrap()
            .linearize(Some(1));

        let i = Bispinor::slot(4, 1);

        let j = Bispinor::slot(4, 3);

        let shadow: NamedStructure<String, ()> =
            VecStructure::from_iter([i.cast::<PhysReps>()]).to_named("u".to_string(), None);
        let shadow_t: MixedTensor<_, VecStructure> =
            ParamOrConcrete::param(shadow.to_shell().expanded_shadow().unwrap().into())
                .cast_structure();

        let rotation: MixedTensor<f64, _> =
            ParamOrConcrete::Concrete(RealOrComplexTensor::Complex(method.bispinor_tensor(i, j)));

        let res = shadow_t
            .contract(&rotation)
            .unwrap()
            .try_into_parametric()
            .unwrap();

        let fn_map = FunctionMap::new();
        let mut params = shadow_t.try_into_parametric().unwrap().tensor.data();
        params.push(Atom::new_var(Atom::I));

        let spinor_eval: EvalTensor<ExpressionEvaluator<Rational>, VecStructure> = res
            .to_evaluation_tree(&fn_map, &params)
            .unwrap()
            .linearize(Some(1));

        Self {
            method,
            lorentz_rotation: lorentz_eval.map_coeff(&|f| Complex::new_re(F::from_f64(f.into()))),
            bispinor_rotation: spinor_eval.map_coeff(&|f| Complex::new_re(F::from_f64(f.into()))),
        }
    }
}

impl RotationMethod {
    pub fn generator(&self, i: AbstractIndex, j: AbstractIndex) -> DataTensor<f64, VecStructure> {
        let structure = VecStructure::from_iter([
            PhysReps::new_slot(Minkowski {}.into(), 4, i),
            PhysReps::new_slot(Minkowski {}.into(), 4, j),
        ]);
        let zero = 0.;
        match self {
            RotationMethod::Identity => {
                let omega = SparseTensor::empty(structure);
                omega.into()
            }
            RotationMethod::Pi2X => {
                let mut omega = SparseTensor::empty(structure);
                omega.set(&[2, 3], -zero.PIHALF()).unwrap();
                omega.set(&[3, 2], zero.PIHALF()).unwrap();
                omega.into()
            }
            RotationMethod::Pi2Y => {
                let mut omega = SparseTensor::empty(structure);
                omega.set(&[1, 3], zero.PIHALF()).unwrap();
                omega.set(&[3, 1], -zero.PIHALF()).unwrap();
                omega.into()
            }
            RotationMethod::Pi2Z => {
                let mut omega = SparseTensor::empty(structure);
                omega.set(&[1, 2], -zero.PIHALF()).unwrap();
                omega.set(&[2, 1], zero.PIHALF()).unwrap();
                omega.into()
            }
            RotationMethod::EulerAngles(alpha, beta, gamma) => DenseTensor::from_data(
                vec![
                    // row 0
                    zero, zero, zero, zero, // row 1
                    zero, zero, -gamma, *beta, // row 2
                    zero, *gamma, zero, -alpha, // row 3
                    zero, -beta, *alpha, zero,
                ],
                structure,
            )
            .unwrap()
            .into(),
        }
    }

    pub fn lorentz_tensor(
        &self,
        i: Slot<Minkowski>,
        j: Slot<Minkowski>,
    ) -> DataTensor<f64, VecStructure> {
        let structure = VecStructure::from_iter([i.cast::<PhysReps>(), j.cast()]);

        match self {
            RotationMethod::Identity => {
                let rot = DenseTensor::from_data(
                    vec![
                        1., 0., 0., 0., // row 1
                        0., -1., 0., 0., // row 2
                        0., 0., -1., 0., // row 3
                        0., 0., 0., -1.,
                    ],
                    structure,
                )
                .unwrap();
                rot.into()
            }
            RotationMethod::Pi2X => {
                let rot = DenseTensor::from_data(
                    vec![
                        1., 0., 0., 0., // row 1
                        0., -1., 0., 0., // row 2
                        0., 0., 0., -1., // row 3
                        0., 0., 1., 0.,
                    ],
                    structure,
                )
                .unwrap();
                rot.into()
            }
            RotationMethod::Pi2Y => {
                let rot = DenseTensor::from_data(
                    vec![
                        1., 0., 0., 0., // row 1
                        0., 0., 0., 1., // row 2
                        0., 0., -1., 0., // row 3
                        0., -1., 0., 0.,
                    ],
                    structure,
                )
                .unwrap();
                rot.into()
            }
            RotationMethod::Pi2Z => {
                let rot = DenseTensor::from_data(
                    vec![
                        1., 0., 0., 0., // row 1
                        0., 0., -1., 0., // row 2
                        0., 1., 0., 0., // row 3
                        0., 0., 0., -1.,
                    ],
                    structure,
                )
                .unwrap();
                rot.into()
            }
            RotationMethod::EulerAngles(alpha, beta, gamma) => DenseTensor::from_data(
                vec![
                    // row 0
                    1.,
                    0.,
                    0.,
                    0.,
                    // row 1
                    0.,
                    -gamma.cos() * beta.cos(),
                    alpha.sin() * beta.sin() * gamma.cos() - alpha.cos() * gamma.sin(),
                    alpha.sin() * gamma.sin() + alpha.cos() * beta.sin() * gamma.cos(),
                    // row 2
                    0.,
                    gamma.sin() * beta.cos(),
                    -alpha.cos() * gamma.cos() - alpha.sin() * beta.sin() * gamma.sin(),
                    -alpha.sin() * gamma.cos() + alpha.cos() * beta.sin() * gamma.sin(),
                    // row 3
                    0.,
                    -beta.sin(),
                    alpha.sin() * beta.cos(),
                    -alpha.cos() * beta.cos(),
                ],
                structure,
            )
            .unwrap()
            .into(),
            // Rotation::EulerAngles(alpha, beta, gamma) => DenseTensor::from_data(
            //     vec![
            //         // row 0
            //         zero.one(),
            //         zero.clone(),
            //         zero.clone(),
            //         zero.clone(),
            //         // row 1
            //         zero.clone(),
            //         alpha.cos() * beta.cos(),
            //         alpha.sin() * beta.sin() * gamma.cos() - alpha.cos() * gamma.sin(),
            //         alpha.sin() * gamma.sin() + alpha.cos() * beta.sin() * gamma.cos(),
            //         // row 2
            //         zero.clone(),
            //         gamma.sin() * beta.cos(),
            //         alpha.cos() * gamma.cos() + alpha.sin() * beta.sin() * gamma.sin(),
            //         -alpha.sin() * gamma.cos() + alpha.cos() * beta.sin() * gamma.sin(),
            //         // row 3
            //         zero.clone(),
            //         -beta.sin(),
            //         alpha.sin() * beta.cos(),
            //         alpha.cos() * beta.cos(),
            //     ],
            //     structure,
            // )
            // .unwrap()
            // .into(),
        }
    }

    pub fn bispinor_tensor(
        &self,
        i: Slot<Bispinor>,
        j: Slot<Bispinor>,
    ) -> DataTensor<Complex<f64>, VecStructure> {
        let structure = VecStructure::from_iter([i.cast::<PhysReps>(), j.cast()]);
        let zero = 0.; // F::new_zero();
        let zeroc = Complex::new_re(zero);

        match self {
            RotationMethod::Identity => {
                let mut rot = SparseTensor::empty(structure);
                rot.set(&[0, 0], zeroc.one()).unwrap();
                rot.set(&[1, 1], zeroc.one()).unwrap();
                rot.set(&[2, 2], zeroc.one()).unwrap();
                rot.set(&[3, 3], zeroc.one()).unwrap();
                rot.into()
            }
            RotationMethod::Pi2X => {
                let mut rot = SparseTensor::empty(structure);
                let sqrt2_half = zero.SQRT_2_HALF();
                let sqrt2_halfim = Complex::new_im(sqrt2_half);
                let sqrt2_halfre = Complex::new_re(sqrt2_half);

                rot.set(&[0, 0], sqrt2_halfim).unwrap();
                rot.set(&[0, 1], sqrt2_halfim).unwrap();
                rot.set(&[1, 0], sqrt2_halfim).unwrap();
                rot.set(&[1, 1], sqrt2_halfre).unwrap();
                rot.set(&[2, 2], sqrt2_halfim).unwrap();
                rot.set(&[2, 3], sqrt2_halfim).unwrap();
                rot.set(&[3, 2], sqrt2_halfim).unwrap();
                rot.set(&[3, 3], sqrt2_halfre).unwrap();
                rot.into()
            }
            RotationMethod::Pi2Y => {
                let mut rot = SparseTensor::empty(structure);
                let sqrt2_half = zero.SQRT_2_HALF();
                let sqrt2_halfim = Complex::new_im(sqrt2_half);
                let sqrt2_halfre = Complex::new_re(sqrt2_half);
                let nsqrt2_half = -sqrt2_halfre;

                rot.set(&[0, 0], sqrt2_halfim).unwrap();
                rot.set(&[0, 1], nsqrt2_half).unwrap();
                rot.set(&[1, 0], nsqrt2_half).unwrap();
                rot.set(&[1, 1], sqrt2_halfre).unwrap();
                rot.set(&[2, 2], sqrt2_halfim).unwrap();
                rot.set(&[2, 3], nsqrt2_half).unwrap();
                rot.set(&[3, 2], nsqrt2_half).unwrap();
                rot.set(&[3, 3], sqrt2_halfre).unwrap();
                rot.into()
            }
            RotationMethod::Pi2Z => {
                let mut rot = SparseTensor::empty(structure);
                let sqrt2_half = zero.SQRT_2_HALF();
                let sqrt2_halfc = Complex::new(sqrt2_half, sqrt2_half);
                let sqrt2_halfcc = sqrt2_halfc.conj();

                rot.set(&[0, 0], sqrt2_halfc).unwrap();
                rot.set(&[1, 1], sqrt2_halfcc).unwrap();
                rot.set(&[2, 2], sqrt2_halfc).unwrap();
                rot.set(&[3, 3], sqrt2_halfcc).unwrap();
                rot.into()
            }
            RotationMethod::EulerAngles(alpha, beta, gamma) => {
                let norm = alpha.square() + beta.square() + gamma.square();

                if alpha.is_zero() && beta.is_zero() {
                    let normhalf = &norm / 2.; //F::from_f64(2.);
                    let cos_phihalf = normhalf.cos();
                    let sin_phihalf = normhalf.sin();

                    let e = Complex::new(cos_phihalf, sin_phihalf);
                    let econj = e.conj();
                    return DenseTensor::from_data(
                        vec![
                            // row 0
                            e, zeroc, zeroc, zeroc, // row 1
                            zeroc, econj, zeroc, zeroc, // row 2
                            zeroc, zeroc, e, zeroc, // row 3
                            zeroc, zeroc, zeroc, econj,
                        ],
                        structure,
                    )
                    .unwrap()
                    .into();
                }

                let complex_phi = Complex::new(*alpha, *beta);

                let normhalf = &norm / 2.; //F::from_f64(2.);
                let cos_phihalf = normhalf.cos();
                let sin_phihalf = normhalf.sin();

                let a_00 = Complex::new(gamma / norm * cos_phihalf, sin_phihalf);
                let a_01 =
                    Complex::new_im((norm - gamma.square() / norm) * sin_phihalf) / complex_phi;
                let a_10 = complex_phi * Complex::new_im(sin_phihalf / norm);
                let a_11 = Complex::new(cos_phihalf, -gamma / norm * sin_phihalf);

                DenseTensor::from_data(
                    vec![
                        // row 0
                        a_00, a_01, zeroc, zeroc, // row 1
                        a_10, a_11, zeroc, zeroc, // row 2
                        zeroc, zeroc, a_00, a_01, // row 3
                        zeroc, zeroc, a_10, a_11,
                    ],
                    structure,
                )
                .unwrap()
                .into()
            }
        }
    }
}

pub trait Rotatable {
    fn rotate(&self, rotation: &Rotation) -> Self;
}

impl<T: FloatLike> Rotatable for ThreeMomentum<F<T>> {
    fn rotate(&self, rotation: &Rotation) -> Self {
        match rotation.method {
            RotationMethod::Identity => self.clone(),
            RotationMethod::Pi2X => ThreeMomentum::perform_pi2_rotation_x(self),
            RotationMethod::Pi2Y => ThreeMomentum::perform_pi2_rotation_y(self),
            RotationMethod::Pi2Z => ThreeMomentum::perform_pi2_rotation_z(self),
            RotationMethod::EulerAngles(alpha, beta, gamma) => {
                let mut result = self.clone();
                result.rotate_mut(&F::from_f64(alpha), &F::from_f64(beta), &F::from_f64(gamma));
                result
            }
        }
    }
}

impl<T: FloatLike> Rotatable for FourMomentum<F<T>> {
    fn rotate(&self, rotation: &Rotation) -> Self {
        Self {
            temporal: self.temporal.clone(),
            spatial: self.spatial.rotate(rotation),
        }
    }
}

impl<T: FloatLike> Rotatable for Polarization<Complex<F<T>>> {
    fn rotate(&self, rotation: &Rotation) -> Self {
        let rotated = match self.pol_type {
            PolType::Epsilon | PolType::EpsilonBar => rotation
                .lorentz_rotation
                .clone()
                .map_coeff(&|t| t.map(|f| F::from_ff64(f)))
                .evaluate(&self.tensor.data)
                .try_into_dense()
                .unwrap()
                .cast_structure(),
            PolType::U | PolType::UBar | PolType::V | PolType::VBar | PolType::Scalar => {
                self.tensor.clone()
            }
        };

        Polarization {
            tensor: rotated,
            pol_type: self.pol_type,
        }
    }
}

#[cfg(test)]
mod tests {

    use core::f64;

    use eyre::Context;
    use spenso::{
        arithmetic::ScalarMul,
        contraction::Contract,
        iterators::IteratableTensor,
        structure::{slot::DualSlotTo, TensorStructure},
        ufo,
        upgrading_arithmetic::{FallibleAdd, FallibleSub},
    };

    use crate::utils::F;

    use super::*;
    use crate::utils::ApproxEq;

    #[test]
    fn polarization() {
        let mom = FourMomentum::from_args(F(1.), F(1.), F(0.), F(0.));

        let pol: Polarization<F<f64>> = Polarization::lorentz(mom.pol_one());
        let pol2 = pol.clone();

        print!("{}", pol.add_fallible(&pol2).unwrap());

        println!("pol_one: {:?}", pol);

        let structure = pol.tensor.structure.clone();

        println!("{}", structure.flat_index([2]).unwrap());

        pol.tensor
            .iter_flat()
            .for_each(|(i, d)| println!("{}{}", i, d));
    }

    #[test]
    fn serialization() {
        let helicity = Helicity::Plus;
        let minus = Helicity::Minus;
        let hels = vec![helicity, minus, Helicity::Zero, helicity, minus];
        let serialized = serde_json::to_string(&hels).unwrap();
        println!("{}", serialized);
        let deserialized: Vec<Helicity> = serde_json::from_str("[1,-1,0,1]").unwrap();

        println!("{:?}", deserialized);
    }

    #[test]
    fn eps() {
        let mom = FourMomentum::from_args(
            F(156.2565490076973),
            F(-108.59017233120495),
            F(-100.2097685717689),
            F(50.81619686346704),
        );

        let pol = mom.pol(SignOrZero::Plus);
        println!("{}", pol);

        let mom = FourMomentum::from_args(
            F(441.7831721921727),
            F(237.7632083614734),
            F(368.76330416225863),
            F(-51.523329523343534),
        );

        let pol = mom.pol(SignOrZero::Plus);
        println!("{}", pol);

        let mom = FourMomentum::from_args(F(485.0355), F(0.), F(0.), F(485.0355));

        let pol = mom.pol(SignOrZero::Plus);

        println!("{}", pol);

        println!("pt{}", mom.pt());

        let mom = FourMomentum::from_args(F(485.0355), F(-107.6044), F(-431.5174), F(193.5805));

        let pol = mom.pol(SignOrZero::Plus).bar();

        // Helicity=           1           1           1           1
        //  W(*,1)=               (5.4707403520912967,0.0000000000000000)               (0.0000000000000000,0.0000000000000000)               (31.622776601683793,0.0000000000000000)               (0.0000000000000000,0.0000000000000000)
        //  W(*,2)=               (0.0000000000000000,0.0000000000000000)             (-0.70710678118654757,0.0000000000000000)              (0.0000000000000000,0.70710678118654757)               (0.0000000000000000,0.0000000000000000)
        //  W(*,3)=               (17.333409072341226,0.0000000000000000)              (6.3994499859138338,-25.663202641304650)               (2.9986797695150327,0.0000000000000000)              (1.1071048475630934,-4.4397340569457056)
        //  W(*,4)=               (0.0000000000000000,0.0000000000000000)        (6.82818793416064135E-002,0.68609707126238750)            (0.27382536157480858,-0.17108713804717862)              (0.64834964048112464,0.0000000000000000)

        println!("{}", pol);
        println!("pt{}", mom.pt());
    }

    #[test]
    fn spinors_degenerate() {
        let mom = FourMomentum::from_args(F(2.), F(0.), F(0.), F(-2.));

        assert_eq!(
            [Complex::new_zero(), Complex::new_re(F(1.))],
            mom.xi(Sign::Positive)
        );

        assert_eq!(
            [Complex::new_re(F(-1.)), Complex::new_re(F(0.))],
            mom.xi(Sign::Negative)
        );

        let u_p = mom.u(Sign::Positive);
        let u_p_target: Polarization<Complex<F<_>>> =
            Polarization::bispinor_u([F(0.), F(0.), F(0.), F(2.)]).cast();

        u_p.approx_eq_res(&u_p_target, &F(0.001))
            .wrap_err("u+((2,0,0,-2)) does not match target: (0,0,0,2)")
            .unwrap();

        let mut u_p_bar_target: Polarization<Complex<F<_>>> =
            Polarization::bispinor_u([F(0.), F(2.), F(0.), F(0.)]).cast();

        u_p_bar_target.pol_type = PolType::UBar;

        u_p.bar()
            .approx_eq_res(&u_p_bar_target, &F(0.001))
            .wrap_err("u+bar((2,0,0,-2)) does not match target: (0,2,0,0)")
            .unwrap();

        let u_m = mom.u(Sign::Negative);
        let u_m_target: Polarization<Complex<F<_>>> =
            Polarization::bispinor_u([F(-2.), F(0.), F(0.), F(0.)]).cast();

        u_m.approx_eq_res(&u_m_target, &F(0.001))
            .wrap_err("u-((2,0,0,-2)) does not match target: (-2,0,0,0)")
            .unwrap();

        let mut u_m_bar_target: Polarization<Complex<F<_>>> =
            Polarization::bispinor_u([F(0.), F(0.), F(-2.), F(0.)]).cast();

        u_m_bar_target.pol_type = PolType::UBar;

        u_m.bar()
            .approx_eq_res(&u_m_bar_target, &F(0.001))
            .wrap_err("u-bar((2,0,0,-2)) does not match target: (0,0,-2,0)")
            .unwrap();

        let v_p = mom.v(Sign::Positive);
        let v_p_target: Polarization<Complex<F<_>>> =
            Polarization::bispinor_v([F(2.), F(0.), F(0.), F(0.)]).cast();

        v_p.approx_eq_res(&v_p_target, &F(0.001))
            .wrap_err("v+((2,0,0,-2)) does not match target: (2,0,0,0)")
            .unwrap();

        let mut v_p_bar_target: Polarization<Complex<F<_>>> =
            Polarization::bispinor_v([F(0.), F(0.), F(2.), F(0.)]).cast();

        v_p_bar_target.pol_type = PolType::VBar;

        v_p.bar()
            .approx_eq_res(&v_p_bar_target, &F(0.001))
            .wrap_err("v+bar((2,0,0,-2)) does not match target: (0,0,2,0)")
            .unwrap();

        let v_m = mom.v(Sign::Negative);
        let v_m_target: Polarization<Complex<F<_>>> =
            Polarization::bispinor_v([F(0.), F(0.), F(0.), F(-2.)]).cast();

        v_m.approx_eq_res(&v_m_target, &F(0.001))
            .wrap_err("v-((2,0,0,-2)) does not match target: (0,0,0,-2)")
            .unwrap();

        let mut v_m_bar_target: Polarization<Complex<F<_>>> =
            Polarization::bispinor_v([F(0.), F(-2.), F(0.), F(0.)]).cast();

        v_m_bar_target.pol_type = PolType::VBar;

        v_m.bar()
            .approx_eq_res(&v_m_bar_target, &F(0.001))
            .wrap_err("v-bar((2,0,0,-2)) does not match target: (0,-2,0,0)")
            .unwrap();
    }

    #[test]
    fn spinors() {
        let mom = FourMomentum::from_args(F(4.), F(0.), F(4.), F(0.));

        assert!(
            Complex::approx_eq_slice(
                &[
                    Complex::new_re(F(f64::consts::FRAC_1_SQRT_2)),
                    Complex::new_im(F(f64::consts::FRAC_1_SQRT_2)),
                ],
                &mom.xi(Sign::Positive),
                &F(0.001),
            ),
            "xi+((4,0,4,0)) does not match target: (1/sqrt(2),i/sqrt(2))"
        );

        assert!(
            Complex::approx_eq_slice(
                &[
                    Complex::new_im(F(f64::consts::FRAC_1_SQRT_2)),
                    Complex::new_re(F(f64::consts::FRAC_1_SQRT_2)),
                ],
                &mom.xi(Sign::Negative),
                &F(0.001),
            ),
            "xi-((4,0,4,0)) does not match target: (i/sqrt(2),1/sqrt(2))"
        );

        let zero: Complex<_> = F(0.).into();
        let u_p = mom.u(Sign::Positive);
        let u_p_target: Polarization<Complex<F<_>>> =
            Polarization::bispinor_u([zero, zero, Complex::new_re(F(2.)), Complex::new_im(F(2.))])
                .cast();

        u_p.approx_eq_res(&u_p_target, &F(0.001))
            .wrap_err("u+((4,0,4,0)) does not match target: (0,0,2,i2)")
            .unwrap();

        let mut u_p_bar_target: Polarization<Complex<F<_>>> =
            Polarization::bispinor_u([Complex::new_re(F(2.)), Complex::new_im(-F(2.)), zero, zero])
                .cast();

        u_p_bar_target.pol_type = PolType::UBar;

        u_p.bar()
            .approx_eq_res(&u_p_bar_target, &F(0.001))
            .wrap_err("u+bar((4,0,4,0)) does not match target: (2,-2i,0,0)")
            .unwrap();

        let u_m = mom.u(Sign::Negative);
        let u_m_target: Polarization<Complex<F<_>>> =
            Polarization::bispinor_u([Complex::new_im(F(2.)), Complex::new_re(F(2.)), zero, zero])
                .cast();

        u_m.approx_eq_res(&u_m_target, &F(0.001))
            .wrap_err("u-((4,0,4,0)) does not match target: (2i,2,0,0)")
            .unwrap();

        let mut u_m_bar_target: Polarization<Complex<F<_>>> =
            Polarization::bispinor_u([zero, zero, Complex::new_im(-F(2.)), Complex::new_re(F(2.))])
                .cast();

        u_m_bar_target.pol_type = PolType::UBar;

        u_m.bar()
            .approx_eq_res(&u_m_bar_target, &F(0.001))
            .wrap_err("u-bar((4,0,4,0)) does not match target: (0,0,-2i,2)")
            .unwrap();

        let v_p = mom.v(Sign::Positive);
        let v_p_target: Polarization<Complex<F<_>>> = Polarization::bispinor_v([
            Complex::new_im(-F(2.)),
            Complex::new_re(-F(2.)),
            zero,
            zero,
        ])
        .cast();

        v_p.approx_eq_res(&v_p_target, &F(0.001))
            .wrap_err("v+((4,0,4,0)) does not match target: (-2i,-2,0,0)")
            .unwrap();

        let mut v_p_bar_target: Polarization<Complex<F<_>>> =
            Polarization::bispinor_v([zero, zero, Complex::new_im(F(2.)), Complex::new_re(F(-2.))])
                .cast();

        v_p_bar_target.pol_type = PolType::VBar;

        v_p.bar()
            .approx_eq_res(&v_p_bar_target, &F(0.001))
            .wrap_err("v+bar((4,0,4,0)) does not match target: (0,0,2i,-2)")
            .unwrap();

        let v_m = mom.v(Sign::Negative);
        let v_m_target: Polarization<Complex<F<_>>> = Polarization::bispinor_v([
            zero,
            zero,
            Complex::new_re(-F(2.)),
            Complex::new_im(-F(2.)),
        ])
        .cast();

        v_m.approx_eq_res(&v_m_target, &F(0.001))
            .wrap_err("v-((4,0,4,0)) does not match target: (0,0,-2,-2i)")
            .unwrap();

        let mut v_m_bar_target: Polarization<Complex<F<_>>> =
            Polarization::bispinor_v([Complex::new_re(-F(2.)), Complex::new_im(F(2.)), zero, zero])
                .cast();

        v_m_bar_target.pol_type = PolType::VBar;

        v_m.bar()
            .approx_eq_res(&v_m_bar_target, &F(0.001))
            .wrap_err("v-bar((4,0,4,0)) does not match target: (-2,2i,0,0)")
            .unwrap();
    }

    #[test]
    fn vectors_degenerate() {
        let mom = FourMomentum::from_args(F(2.), F(0.), F(0.), F(-2.));

        assert!(
            F::approx_eq_slice(&mom.pol_one(), &[0., 1., 0., 0.].map(F), &F(0.001)),
            "pol_one((2,0,0,-2)) does not match target: (0,1,0,0)"
        );

        assert!(
            F::approx_eq_slice(&mom.pol_two(), &[0., 0., -1., 0.].map(F), &F(0.001)),
            "pol_two((2,0,0,-2)) does not match target: (0,0,-1,0)"
        );

        let e_p = mom.pol(SignOrZero::Plus);

        let zero = F(0.).into();
        let e_p_target = Polarization::lorentz([
            zero,
            Complex::new_re(F(-f64::consts::FRAC_1_SQRT_2)),
            Complex::new_im(F(f64::consts::FRAC_1_SQRT_2)),
            zero,
        ]);

        let mut e_p_bar_target = Polarization::lorentz([
            zero,
            Complex::new_re(F(-f64::consts::FRAC_1_SQRT_2)),
            Complex::new_im(F(-f64::consts::FRAC_1_SQRT_2)),
            zero,
        ]);

        e_p_bar_target.pol_type = PolType::EpsilonBar;

        e_p.approx_eq_res(&e_p_target, &F(0.0001))
            .wrap_err("e+((2,0,0,-2)) does not match target: (0,-1/sqrt(2),i/sqrt(2),0)")
            .unwrap();

        e_p.bar()
            .approx_eq_res(&e_p_bar_target, &F(0.0001))
            .wrap_err("e+bar((2,0,0,-2)) does not match target: (0,-1/sqrt(2),-i/sqrt(2),0)")
            .unwrap();

        let e_m = mom.pol(SignOrZero::Minus);

        let e_m_target = Polarization::lorentz([
            zero,
            Complex::new_re(F(f64::consts::FRAC_1_SQRT_2)),
            Complex::new_im(F(f64::consts::FRAC_1_SQRT_2)),
            zero,
        ]);

        let mut e_m_bar_target = Polarization::lorentz([
            zero,
            Complex::new_re(F(f64::consts::FRAC_1_SQRT_2)),
            Complex::new_im(F(-f64::consts::FRAC_1_SQRT_2)),
            zero,
        ]);

        e_m_bar_target.pol_type = PolType::EpsilonBar;

        e_m.approx_eq_res(&e_m_target, &F(0.0001))
            .wrap_err("e-((2,0,0,-2)) does not match target: (0,1/sqrt(2),i/sqrt(2),0)")
            .unwrap();

        e_m.bar()
            .approx_eq_res(&e_m_bar_target, &F(0.0001))
            .wrap_err("e-bar((2,0,0,-2)) does not match target: (0,1/sqrt(2),-i/sqrt(2),0)")
            .unwrap();
    }

    #[test]
    fn vectors() {
        let mom = FourMomentum::from_args(F(2.), F(0.), F(0.), F(-2.));

        assert!(
            F::approx_eq_slice(&mom.pol_one(), &[0., 1., 0., 0.].map(F), &F(0.001)),
            "pol_one((2,0,0,-2)) does not match target: (0,1,0,0)"
        );

        assert!(
            F::approx_eq_slice(&mom.pol_two(), &[0., 0., -1., 0.].map(F), &F(0.001)),
            "pol_two((2,0,0,-2)) does not match target: (0,0,-1,0)"
        );

        let e_p = mom.pol(SignOrZero::Plus);

        let zero = F(0.).into();
        let e_p_target = Polarization::lorentz([
            zero,
            Complex::new_re(F(-f64::consts::FRAC_1_SQRT_2)),
            Complex::new_im(F(f64::consts::FRAC_1_SQRT_2)),
            zero,
        ]);

        let mut e_p_bar_target = Polarization::lorentz([
            zero,
            Complex::new_re(F(-f64::consts::FRAC_1_SQRT_2)),
            Complex::new_im(F(-f64::consts::FRAC_1_SQRT_2)),
            zero,
        ]);

        e_p_bar_target.pol_type = PolType::EpsilonBar;

        e_p.approx_eq_res(&e_p_target, &F(0.0001))
            .wrap_err("e+((2,0,0,-2)) does not match target: (0,-1/sqrt(2),i/sqrt(2),0)")
            .unwrap();

        e_p.bar()
            .approx_eq_res(&e_p_bar_target, &F(0.0001))
            .wrap_err("e+bar((2,0,0,-2)) does not match target: (0,-1/sqrt(2),-i/sqrt(2),0)")
            .unwrap();

        let e_m = mom.pol(SignOrZero::Minus);

        let e_m_target = Polarization::lorentz([
            zero,
            Complex::new_re(F(f64::consts::FRAC_1_SQRT_2)),
            Complex::new_im(F(f64::consts::FRAC_1_SQRT_2)),
            zero,
        ]);

        let mut e_m_bar_target = Polarization::lorentz([
            zero,
            Complex::new_re(F(f64::consts::FRAC_1_SQRT_2)),
            Complex::new_im(F(-f64::consts::FRAC_1_SQRT_2)),
            zero,
        ]);

        e_m_bar_target.pol_type = PolType::EpsilonBar;

        e_m.approx_eq_res(&e_m_target, &F(0.0001))
            .wrap_err("e-((2,0,0,-2)) does not match target: (0,1/sqrt(2),i/sqrt(2),0)")
            .unwrap();

        e_m.bar()
            .approx_eq_res(&e_m_bar_target, &F(0.0001))
            .wrap_err("e-bar((2,0,0,-2)) does not match target: (0,1/sqrt(2),-i/sqrt(2),0)")
            .unwrap();
    }

    #[test]
    fn rotations() {
        let mom = FourMomentum::from_args(F(2.), F(3.), F(1.), F(2.));

        let momc: [F<f64>; 4] = mom.into();
        let e = Polarization::lorentz(momc.map(Complex::new_re));
        let u = mom.u(Sign::Negative);

        let rot: Rotation = RotationMethod::EulerAngles(std::f64::consts::PI / 2., 0., 0.).into();
        let rotx = RotationMethod::Pi2X.into();

        let rotye: Rotation = RotationMethod::EulerAngles(0., std::f64::consts::PI / 2., 0.).into();
        let roty = RotationMethod::Pi2Y.into();

        let rotze: Rotation = RotationMethod::EulerAngles(0., 0., std::f64::consts::PI / 2.).into();
        let rotz = RotationMethod::Pi2Z.into();

        let rotid = RotationMethod::Identity.into();
        // println!(
        //     "{}",
        //     rot.bispinor_tensor(
        //         Bispinor::new_slot_selfless(4, 2),
        //         Bispinor::new_slot_selfless(4, 1),
        //     )
        // );

        let mom_rot = mom.rotate(&rot);
        let mom_id = mom.rotate(&rotid);
        let other_mom_rot = mom.rotate(&rotx);
        let mom_rotye = mom.rotate(&rotye);
        let mom_roty = mom.rotate(&roty);
        let mom_rotze = mom.rotate(&rotze);
        let mom_rotz = mom.rotate(&rotz);

        println!("id:\t{}", mom_id);
        mom_id
            .approx_eq_res(&mom, &F(0.001))
            .wrap_err("mom is not identical to mom transformed with identity")
            .unwrap();
        println!("orig:\t{}", mom);

        println!("rotxe:\t{}", mom_rot);
        println!("pi2x:\t{}", other_mom_rot);
        mom_rot.approx_eq_res(&other_mom_rot, &F(0.001)).unwrap();

        println!("rotye:\t{}", mom_rotye);
        println!("pi2y:\t{}", mom_roty);

        mom_rotye.approx_eq_res(&mom_roty, &F(0.001)).unwrap();
        println!("rotze:\t{}", mom_rotze);
        println!("pi2z:\t{}", mom_rotz);
        mom_rotze.approx_eq_res(&mom_rotz, &F(0.001)).unwrap();

        let e_rot = e.rotate(&rot);
        println!("rotxe:{}", e_rot);
        mom_rot.approx_eq_res(&e_rot, &F(0.001)).unwrap();

        let e_id = e.rotate(&rotid);
        mom_id.approx_eq_res(&e_id, &F(0.001)).unwrap();
        let e_rotye = e.rotate(&rotye);

        println!("rotye:{}", e_rotye);
        mom_rotye.approx_eq_res(&e_rotye, &F(0.001)).unwrap();
        let e_roty = e.rotate(&roty);
        mom_roty.approx_eq_res(&e_roty, &F(0.001)).unwrap();
        let e_rotze = e.rotate(&rotze);
        mom_rotze.approx_eq_res(&e_rotze, &F(0.001)).unwrap();
        let e_rotz = e.rotate(&rotz);
        mom_rotz.approx_eq_res(&e_rotz, &F(0.001)).unwrap();
        let u_rot = u.rotate(&rot);
        let u_id = u.rotate(&rotid);

        let other_e_rot = e.rotate(&rotx);

        let other_u_rot = u.rotate(&rotx);

        println!("id:{}", e_id);
        println!("orig:{}", e);
        println!("euler:{}", e_rot);
        println!("Pi2X:{}", other_e_rot);
        println!("rotye:{}", e_rotye);
        println!("pi2y:{}", e_roty);
        println!("rotze:{}", e_rotze);
        println!("pi2z:{}", e_rotz);

        println!("id:{}", u_id);
        println!("orig:{}", u);
        println!("euler:{}", u_rot);
        println!("Pi2X {}", other_u_rot);
    }

    #[test]
    fn omega() {
        let mu = PhysReps::new_slot(Minkowski {}.into(), 4, 0);

        let nu = PhysReps::new_slot(Minkowski {}.into(), 4, 1);

        let mud = mu.dual();

        let nud = nu.dual();

        let i = PhysReps::new_slot(Bispinor {}.into(), 4, 2);

        let j = PhysReps::new_slot(Bispinor {}.into(), 4, 3);

        let k = PhysReps::new_slot(Bispinor {}.into(), 4, 4);

        let zero = Atom::new_num(0);
        let theta_x = Atom::parse("theta_x").unwrap();
        let theta_y = Atom::parse("theta_y").unwrap();
        let theta_z = Atom::parse("theta_z").unwrap();

        let omega = DenseTensor::from_data(
            vec![
                zero.clone(),
                zero.clone(),
                zero.clone(),
                zero.clone(),
                //
                zero.clone(),
                zero.clone(),
                theta_z.clone(),
                -theta_y.clone(),
                //
                zero.clone(),
                -theta_z.clone(),
                zero.clone(),
                theta_x.clone(),
                //
                zero.clone(),
                theta_y.clone(),
                -theta_x.clone(),
                zero.clone(),
            ],
            VecStructure::from_iter([mu, nu]),
        )
        .unwrap();

        println!("{}", omega);

        let gammamujk: SparseTensor<Complex<f64>> =
            ufo::gamma_data_weyl(VecStructure::from_iter([mud, j, k]));

        let gammamukj: SparseTensor<Complex<f64>> =
            ufo::gamma_data_weyl(VecStructure::from_iter([mud, k, i]));

        let gammanuki: SparseTensor<Complex<f64>> =
            ufo::gamma_data_weyl(VecStructure::from_iter([nud, k, i]));

        let gammanujk: SparseTensor<Complex<f64>> =
            ufo::gamma_data_weyl(VecStructure::from_iter([nud, j, k]));

        let sigma = (gammamujk
            .contract(&gammanuki)
            .unwrap()
            .sub_fallible(&gammanujk.contract(&gammamukj).unwrap())
            .unwrap())
        .scalar_mul(&Complex::new_im(0.25))
        .unwrap();

        println!("{}", sigma);

        println!("{}", omega.contract(&sigma).unwrap());
    }
}
