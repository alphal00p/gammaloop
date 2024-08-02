use std::{
    fmt::{Display, LowerExp},
    ops::{Add, AddAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign},
};

use log::debug;

use serde::{Deserialize, Serialize};
use spenso::{
    arithmetic::ScalarMul,
    contraction::RefZero,
    data::DenseTensor,
    structure::{
        AbstractIndex, BaseRepName, Bispinor, Euclidean, IndexLess, Lorentz, NamedStructure,
        PhysReps, RepName, VecStructure,
    },
    upgrading_arithmetic::{FallibleAdd, FallibleMul},
};
use symbolica::{
    atom::{Atom, Symbol},
    coefficient::Coefficient,
    domains::{
        float::{NumericalFloatComparison, NumericalFloatLike, Real},
        integer::IntegerRing,
        rational::RationalField,
    },
    poly::{polynomial::MultivariatePolynomial, Exponent},
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

impl<T: RefZero> RefZero for Energy<T> {
    fn ref_zero(&self) -> Self {
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

#[derive(Debug, PartialEq, Clone)]
pub struct Polarization<T> {
    tensor: DenseTensor<T, IndexLess<PhysReps>>,
}

impl<T: FloatLike> Display for Polarization<F<T>> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Pol: {}", self.tensor)
    }
}

impl<T: FloatLike> Display for Polarization<Complex<F<T>>> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Pol: {}", self.tensor)
    }
}

impl<T> Polarization<T> {
    pub fn scalar(value: T) -> Self {
        let structure = IndexLess::new(vec![]);
        Polarization {
            tensor: DenseTensor {
                data: vec![value],
                structure,
            },
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
        }
    }

    pub fn bispinor(value: [T; 4]) -> Self {
        let structure = IndexLess::new(vec![Bispinor::new_dimed_rep_selfless(4).cast()]);
        let [v1, v2, v3, v4] = value;
        Polarization {
            tensor: DenseTensor {
                data: vec![v1, v2, v3, v4],
                structure,
            },
        }
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
        })
    }
}

impl<T: Neg<Output = T>> Neg for Polarization<T> {
    type Output = Polarization<T>;
    fn neg(self) -> Self::Output {
        Polarization {
            tensor: -self.tensor,
        }
    }
}

impl<T> Polarization<T> {
    fn cast<U>(&self) -> Polarization<U>
    where
        T: Clone,
        U: Clone + From<T>,
    {
        Polarization {
            tensor: self.tensor.cast(),
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
        T: Real + NumericalFloatComparison,
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

        debug!("pol_one in: {}", self);
        let pt = self.pt();
        let p = self.spatial.norm();

        let (e1, e2, e3) = if pt.is_zero() {
            (pt.one(), pt.zero(), pt.zero())
        } else {
            (
                &self.spatial.px * &self.spatial.pz / (&pt * &p),
                &self.spatial.pz * &self.spatial.pz / (&pt * &p),
                -(&pt / &p),
            )
        };

        debug!(
            " (pt.zero(), e1, e2, e3) {} {} {} {}",
            pt.zero(),
            e1,
            e2,
            e3
        );

        let pol = Polarization::lorentz([pt.zero(), e1, e2, e3]);

        debug!("pol :{pol}");
        pol
    }

    pub fn pol_two(&self) -> Polarization<F<T>>
    where
        T: FloatLike,
    {
        // definition from helas_ref A.2
        let pt = self.pt();
        let (e1, e2, e3) = if pt.is_zero() {
            (pt.zero(), pt.one(), pt.zero())
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

    pub fn pol(&self, lambda: SignOrZero) -> Polarization<Complex<F<T>>> {
        if lambda.is_zero() {
            self.pol_three().cast()
        } else {
            let one = self.temporal.value.one();
            let sqrt_2_inv: Complex<F<T>> = (&one + &one).sqrt().inv().into();
            let i = one.i();

            i.add_fallible(&i).unwrap();

            let eone: Polarization<Complex<F<T>>> = lambda * self.pol_one().cast();

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
        Polarization::bispinor(match lambda {
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
        Polarization::bispinor(match lambda {
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
        }
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy, Serialize, Deserialize)]
pub enum Sign {
    Positive,
    Negative,
}

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

#[derive(Debug, PartialEq, Eq, Clone, Copy, Serialize, Deserialize)]
pub enum SignOrZero {
    Sign(Sign),
    Zero,
}

impl SignOrZero {
    fn is_zero(&self) -> bool {
        matches!(self, SignOrZero::Zero)
    }

    #[allow(dead_code)]
    fn is_sign(&self) -> bool {
        matches!(self, SignOrZero::Sign(_))
    }
}

impl<T: Neg<Output = T> + RefZero> Mul<T> for SignOrZero {
    type Output = T;
    fn mul(self, rhs: T) -> Self::Output {
        match self {
            SignOrZero::Sign(s) => s * rhs,
            SignOrZero::Zero => rhs.ref_zero(),
        }
    }
}

impl Neg for SignOrZero {
    type Output = Self;
    fn neg(self) -> Self::Output {
        match self {
            SignOrZero::Sign(s) => SignOrZero::Sign(-s),
            SignOrZero::Zero => SignOrZero::Zero,
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
#[cfg(test)]
mod tests {
    use spenso::{
        iterators::IteratableTensor, structure::TensorStructure, upgrading_arithmetic::FallibleAdd,
    };

    use crate::utils::F;

    use super::FourMomentum;

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
}
