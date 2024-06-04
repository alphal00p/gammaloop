use std::{
    fmt::{Display, LowerExp},
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use num::{Float, Zero};

use spenso::{AbstractIndex, DenseTensor, Representation, VecStructure};
use symbolica::{
    atom::Atom,
    coefficient::Coefficient,
    domains::rational::RationalField,
    poly::{polynomial::MultivariatePolynomial, Exponent},
};

use crate::utils::FloatLike;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct Energy<T> {
    pub value: T,
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

impl<T> Zero for Energy<T>
where
    T: Zero,
{
    fn zero() -> Self {
        Energy { value: T::zero() }
    }

    fn is_zero(&self) -> bool {
        self.value.is_zero()
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

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct ThreeMomentum<T> {
    pub px: T,
    pub py: T,
    pub pz: T,
}

impl<T: Clone> From<[T; 3]> for ThreeMomentum<T> {
    fn from(data: [T; 3]) -> Self {
        ThreeMomentum {
            px: data[0].clone(),
            py: data[1].clone(),
            pz: data[2].clone(),
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
            VecStructure::new(vec![(index, Representation::Euclidean(3.into())).into()]);
        DenseTensor::from_data(&[self.px, self.py, self.pz], structure).unwrap()
    }

    pub fn into_dense_param(self, index: AbstractIndex) -> DenseTensor<T, VecStructure>
    where
        T: Clone,
    {
        let structure =
            VecStructure::new(vec![(index, Representation::Euclidean(3.into())).into()]);
        DenseTensor::from_data(&[self.px, self.py, self.pz], structure).unwrap()
    }
}

impl<T: FloatLike> ThreeMomentum<T> {
    pub fn rotate_mut(&mut self, alpha: T, beta: T, gamma: T) {
        let sin_alpha = alpha.sin();
        let cos_alpha = alpha.cos();
        let sin_beta = beta.sin();
        let cos_beta = beta.cos();
        let sin_gamma = gamma.sin();
        let cos_gamma = gamma.cos();

        let px = self.px.clone();
        let py = self.py.clone();
        let pz = self.pz.clone();

        self.px = cos_alpha * cos_beta * px
            + (-cos_alpha * sin_gamma + sin_alpha * sin_beta * cos_gamma) * py
            + (sin_alpha * sin_gamma + cos_alpha * sin_beta * cos_gamma) * pz;

        self.py = sin_gamma * cos_beta * px
            + (cos_alpha * cos_gamma + sin_alpha * sin_beta * sin_gamma) * py
            + (-sin_alpha * cos_gamma + cos_alpha * sin_beta * sin_gamma) * pz;

        self.pz = -sin_beta * px + cos_beta * sin_alpha * py + cos_alpha * cos_beta * pz;
    }

    pub fn rotate(&self, alpha: T, beta: T, gamma: T) -> Self {
        let mut result = self.clone();
        result.rotate_mut(alpha, beta, gamma);
        result
    }
}

impl<T: Neg<Output = T> + Clone> ThreeMomentum<T> {
    pub fn perform_pi2_rotation_x_mut(&mut self) {
        let px = self.px.clone();
        self.px = -px;
    }

    pub fn perform_pi2_rotation_x(&self) -> Self {
        Self {
            px: -self.px.clone(),
            py: self.py.clone(),
            pz: self.pz.clone(),
        }
    }

    pub fn perform_pi2_rotation_y_mut(&mut self) {
        let py = self.py.clone();
        self.py = -py;
    }

    pub fn perform_pi2_rotation_y(&self) -> Self {
        Self {
            px: self.px.clone(),
            py: -self.py.clone(),
            pz: self.pz.clone(),
        }
    }

    pub fn perform_pi2_rotation_z_mut(&mut self) {
        let pz = self.pz.clone();
        self.pz = -pz;
    }

    pub fn perform_pi2_rotation_z(&self) -> Self {
        Self {
            px: self.px.clone(),
            py: self.py.clone(),
            pz: -self.pz.clone(),
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
    pub fn norm_squared(self) -> T
    where
        T: Mul<T, Output = T> + Add<T, Output = T> + Clone,
    {
        self.px.clone() * self.px + self.py.clone() * self.py + self.pz.clone() * self.pz
    }

    pub fn on_shell_energy(self, mass: Option<T>) -> Energy<T>
    where
        T: Mul<T, Output = T> + Add<T, Output = T> + Clone + std::ops::Add<Output = T> + Float,
    {
        let energy_squared = self.on_shell_energy_squared(mass);
        Energy {
            value: energy_squared.value.sqrt(),
        }
    }

    pub fn on_shell_energy_squared(self, mass: Option<T>) -> Energy<T>
    where
        T: Mul<T, Output = T> + Add<T, Output = T> + Clone + std::ops::Add<Output = T> + Float,
    {
        let p2 = self.norm_squared();
        if let Some(mass) = mass {
            Energy {
                value: p2 + mass * mass,
            }
        } else {
            Energy { value: p2 }
        }
    }

    pub fn norm(self) -> T
    where
        T: Mul<T> + Add<T> + Float,
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
        T: Mul<T> + Add<T> + std::ops::Add<Output = T> + Float,
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
}

#[derive(Default, Debug, PartialEq, Eq, Clone, Copy)]
pub struct FourMomentum<T, U = T> {
    pub temporal: Energy<U>,
    pub spatial: ThreeMomentum<T>,
}

impl<T, U> FourMomentum<T, U> {
    pub fn new(energy: Energy<U>, three_momentum: ThreeMomentum<T>) -> Self {
        FourMomentum {
            temporal: energy,
            spatial: three_momentum,
        }
    }
}
impl<T> FourMomentum<T, T> {
    pub fn square(self) -> T
    where
        T: Mul<T, Output = T> + Add<T, Output = T> + Clone + Sub<T, Output = T>,
    {
        let temporal = self.temporal.value.clone();
        let spatial = self.spatial.norm_squared();
        temporal * self.temporal.value - spatial
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
        T: Mul<T> + Add<T> + std::ops::Add<Output = T> + Float,
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
        let structure = VecStructure::new(vec![(index, Representation::Lorentz(4.into())).into()]);
        DenseTensor::from_data(
            &[
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
}

impl<T> FourMomentum<T, Atom> {
    pub fn into_dense_param(
        self,
        index: AbstractIndex,
    ) -> DenseTensor<MultivariatePolynomial<RationalField, T>, VecStructure>
    where
        T: Clone + Into<Coefficient> + Exponent,
    {
        let structure = VecStructure::new(vec![(index, Representation::Lorentz(4.into())).into()]);
        let energy = self
            .temporal
            .value
            .to_polynomial(&RationalField::new(), None);

        let px: MultivariatePolynomial<RationalField, _> =
            Atom::new_num(self.spatial.px).to_polynomial(&RationalField::new(), None);
        let py = Atom::new_num(self.spatial.py).to_polynomial(&RationalField::new(), None);
        let pz = Atom::new_num(self.spatial.pz).to_polynomial(&RationalField::new(), None);

        DenseTensor::from_data(&[energy, px, py, pz], structure).unwrap()
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
