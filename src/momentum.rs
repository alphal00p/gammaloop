use std::{
    borrow::Borrow,
    fmt::{Display, LowerExp},
    ops::{Add, AddAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign},
};

use bincode::{Decode, Encode};
use momtrop::vector::Vector;
use serde::{Deserialize, Serialize};
use serde_repr::{Deserialize_repr, Serialize_repr};
use spenso::{
    arithmetic::ScalarMul,
    contraction::RefZero,
    data::DenseTensor,
    parametric::FlatCoefficent,
    structure::{
        AbstractIndex, BaseRepName, Bispinor, Euclidean, IndexLess, Lorentz, NamedStructure,
        NoArgs, PhysReps, RepName, ToSymbolic, VecStructure,
    },
    upgrading_arithmetic::{FallibleAdd, FallibleMul},
};
use symbolica::{
    atom::{Atom, Symbol},
    coefficient::Coefficient,
    domains::{
        float::{NumericalFloatLike, Real, RealNumberLike, SingleFloat},
        integer::IntegerRing,
        rational::RationalField,
    },
    poly::{polynomial::MultivariatePolynomial, Exponent},
    state::State,
};

use spenso::complex::Complex;

use crate::utils::{FloatLike, RefDefault, F};

#[derive(Debug, PartialEq, Eq, Clone, Copy, Serialize, Deserialize)]
pub struct Energy<T> {
    pub value: T,
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

impl<'b, T> Add<Energy<T>> for &'b Energy<T>
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

    pub fn rotate_mut(&mut self, alpha: F<T>, beta: F<T>, gamma: F<T>) {
        let sin_alpha = alpha.sin();
        let cos_alpha = alpha.cos();
        let sin_beta = beta.sin();
        let cos_beta = beta.cos();
        let sin_gamma = gamma.sin();
        let cos_gamma = gamma.cos();

        let px = self.px.clone();
        let py = self.py.clone();
        let pz = self.pz.clone();

        self.px = cos_alpha.clone() * &cos_beta * &px
            + (-(cos_alpha.clone()) * &sin_gamma + sin_alpha.clone() * &sin_beta * &cos_gamma)
                * &py
            + (sin_alpha.clone() * &sin_gamma + cos_alpha.clone() * &sin_beta * &cos_gamma) * &pz;

        self.py = sin_gamma.clone() * &cos_beta * &px
            + (cos_alpha.clone() * &cos_gamma + sin_alpha.clone() * &sin_beta * &sin_gamma) * &py
            + (-sin_alpha.clone() * &cos_gamma + cos_alpha.clone() * &sin_beta * &sin_gamma) * &pz;

        self.pz =
            -sin_beta * &px + cos_beta.clone() * &sin_alpha * &py + cos_alpha * &cos_beta * &pz;
    }

    pub fn rotate(&self, alpha: F<T>, beta: F<T>, gamma: F<T>) -> Self {
        let mut result = self.clone();
        result.rotate_mut(alpha, beta, gamma);
        result
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

impl<'b, T> Add<ThreeMomentum<T>> for &'b ThreeMomentum<T>
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

impl<'a, T> Mul<ThreeMomentum<T>> for &'a ThreeMomentum<T>
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

impl<'a, T> Mul<T> for &'a ThreeMomentum<T>
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

impl<'a, T> Mul<&T> for &'a ThreeMomentum<T>
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

impl<'a, T> Neg for &'a ThreeMomentum<T>
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
        T: for<'a> Mul<&'a T, Output = T> + Add<T, Output = T> + Clone + std::ops::Add<Output = T>,
    {
        let p2 = self.norm_squared();
        if let Some(mass) = mass {
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

#[derive(Debug, PartialEq, Clone, Copy)]
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

#[derive(Debug, PartialEq, Clone)]
pub struct Polarization<T> {
    pub tensor: DenseTensor<T, IndexLess<PhysReps>>,
    pol_type: PolType,
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

impl<T> Polarization<T> {
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
        let structure = IndexLess::new(vec![Lorentz::new_dimed_rep_selfless(4).cast()]);
        let [v1, v2, v3, v4] = value;
        Polarization {
            tensor: DenseTensor {
                data: vec![v1, v2, v3, v4],
                structure,
            },
            pol_type: PolType::Epsilon,
        }
    }

    pub fn bispinor_u(value: [T; 4]) -> Self {
        let structure = IndexLess::new(vec![Bispinor::new_dimed_rep_selfless(4).cast()]);
        let [v1, v2, v3, v4] = value;
        Polarization {
            tensor: DenseTensor {
                data: vec![v1, v2, v3, v4],
                structure,
            },
            pol_type: PolType::U,
        }
    }

    pub fn bispinor_v(value: [T; 4]) -> Self {
        let structure = IndexLess::new(vec![Bispinor::new_dimed_rep_selfless(4).cast()]);
        let [v1, v2, v3, v4] = value;
        Polarization {
            tensor: DenseTensor {
                data: vec![v1, v2, v3, v4],
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
                name: Some(State::get_symbol(self.pol_type.to_string())),
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

impl<T, U, Out> Mul<&U> for Polarization<T>
where
    T: FallibleMul<U, Output = Out>,
{
    type Output = Polarization<Out>;
    fn mul(self, rhs: &U) -> Self::Output {
        Polarization {
            tensor: self.tensor.scalar_mul(rhs).unwrap(),
            pol_type: self.pol_type,
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

impl<T: Neg<Output = T>> Neg for Polarization<T> {
    type Output = Polarization<T>;
    fn neg(self) -> Self::Output {
        Polarization {
            tensor: -self.tensor,
            pol_type: self.pol_type,
        }
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

impl<'a, T> Mul<T> for &'a FourMomentum<T, T>
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

impl<'a, T> Mul<&T> for &'a FourMomentum<T, T>
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
        FourMomentum {
            temporal: energy,
            spatial: three_momentum,
        }
    }

    pub fn into_dense(self, index: AbstractIndex) -> DenseTensor<T, VecStructure>
    where
        T: Clone,
    {
        let structure = VecStructure::from_iter([PhysReps::new_slot(Lorentz {}.into(), 4, index)]);
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
        let structure = VecStructure::from_iter([PhysReps::new_slot(Lorentz {}.into(), 4, index)])
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
            boost_vector.spatial.px.mul_add(&factor, &self.spatial.px),
            boost_vector.spatial.py.mul_add(&factor, &self.spatial.py),
            boost_vector.spatial.pz.mul_add(&factor, &self.spatial.pz),
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

    pub fn pol_one(&self) -> Polarization<F<T>>
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
        Polarization::lorentz([pt.zero(), e1, e2, e3])

        // debug!("pol :{pol}");
    }

    pub fn pol_two(&self) -> Polarization<F<T>>
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
        Polarization::lorentz([pt.zero(), e1, e2, e3])
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
            let sqrt_2_inv: Complex<F<T>> = (&one + &one).sqrt().inv().into();
            let i = one.i();

            let eone: Polarization<Complex<F<T>>> = -lambda * self.pol_one().cast();

            let etwo: Polarization<Complex<F<T>>> = self.pol_two().cast();

            (eone.add_fallible(&(-(etwo * &i)))).unwrap() * &sqrt_2_inv
        }
    }

    pub fn omega(&self, lambda: Sign) -> Complex<F<T>> {
        match lambda {
            Sign::Positive => (&self.temporal.value + self.spatial.norm()).complex_sqrt(),
            Sign::Negative => (&self.temporal.value - self.spatial.norm()).complex_sqrt(),
        }
    }

    pub fn u(&self, lambda: Sign) -> Polarization<Complex<F<T>>> {
        let zero: Complex<F<T>> = self.temporal.value.zero().into();
        Polarization::bispinor_u(match lambda {
            Sign::Positive => [
                zero.clone(),
                self.omega(Sign::Negative),
                zero.clone(),
                self.omega(Sign::Positive),
            ],
            Sign::Negative => [
                -self.omega(Sign::Positive),
                zero.clone(),
                self.omega(Sign::Negative),
                zero.clone(),
            ],
        })
    }

    pub fn v(&self, lambda: Sign) -> Polarization<Complex<F<T>>> {
        let zero: Complex<F<T>> = self.temporal.value.zero().into();
        Polarization::bispinor_v(match lambda {
            Sign::Negative => [
                zero.clone(),
                self.omega(Sign::Negative),
                zero.clone(),
                -self.omega(Sign::Positive),
            ],
            Sign::Positive => [
                self.omega(Sign::Positive),
                zero.clone(),
                -self.omega(Sign::Negative),
                zero.clone(),
            ],
        })
    }
}

impl<T: FloatLike> Polarization<Complex<F<T>>> {
    pub fn bar(&self) -> Self {
        Polarization {
            tensor: self.tensor.map_data_ref(Complex::conj),
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

#[derive(Debug, PartialEq, Eq, Clone, Copy, Serialize, Deserialize)]
#[repr(i8)]
pub enum Sign {
    Positive = 1,
    Negative = -1,
}

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

impl Display for SignOrZero {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SignOrZero::Zero => write!(f, ""),
            SignOrZero::Plus => write!(f, "+"),
            SignOrZero::Minus => write!(f, "-"),
        }
    }
}

pub type Helicity = SignOrZero;

#[derive(
    Debug, PartialEq, Eq, Clone, Serialize, Deserialize, Encode, Decode, PartialOrd, Ord, Hash,
)]
pub struct Signature(Vec<SignOrZero>);

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
            self.spatial.px,
            self.spatial.py,
            self.spatial.pz,
            self.temporal.value,
        ];
        [e, px, py, pz].into_iter()
    }
}

impl<'a, T> IntoIterator for &'a FourMomentum<T> {
    type Item = &'a T;
    type IntoIter = std::array::IntoIter<&'a T, 4>;

    fn into_iter(self) -> Self::IntoIter {
        let [e, px, py, pz] = [
            &self.spatial.px,
            &self.spatial.py,
            &self.spatial.pz,
            &self.temporal.value,
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
        let structure = VecStructure::from_iter([PhysReps::new_slot(Lorentz {}.into(), 4, index)]);
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

impl<'a, T, U> Neg for &'a FourMomentum<T, U>
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

impl From<Vector<f64, 3>> for ThreeMomentum<F<f64>> {
    fn from(value: Vector<f64, 3>) -> Self {
        ThreeMomentum::new(F(value[0]), F(value[1]), F(value[2]))
    }
}
#[cfg(test)]
mod tests {
    use spenso::{
        iterators::IteratableTensor, structure::TensorStructure, upgrading_arithmetic::FallibleAdd,
    };

    use crate::utils::F;

    use super::*;

    #[test]
    fn polarization() {
        let mom = FourMomentum::from_args(F(1.), F(1.), F(0.), F(0.));

        let pol = mom.pol_one();
        let pol2 = pol.clone();

        print!("{}", pol.add_fallible(&pol2).unwrap());

        println!("pol_one: {:?}", pol);

        let structure = pol.tensor.structure.clone();

        println!("{}", structure.flat_index(&[2]).unwrap());

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
}
