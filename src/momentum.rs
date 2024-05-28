use std::{
    fmt::{Display, LowerExp},
    ops::{Add, AddAssign, Mul, Sub, SubAssign},
};

use num::Float;

use rand::seq::index;
use spenso::{AbstractIndex, DenseTensor, Representation, VecStructure};
use symbolica::{
    atom::Atom,
    coefficient::Coefficient,
    domains::rational::RationalField,
    poly::{polynomial::MultivariatePolynomial, Exponent},
};

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct Energy<T> {
    pub value: T,
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

impl<T> ThreeMomentum<T> {
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

    pub fn on_shell_energy(self, mass: T) -> Energy<T>
    where
        T: Mul<T, Output = T> + Add<T, Output = T> + Clone + std::ops::Add<Output = T> + Float,
    {
        let p2 = self.norm_squared();
        let value = (p2 + mass * mass).sqrt();
        Energy { value }
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
            energy,
            three_momentum: self,
        }
    }

    pub fn into_on_shell_four_momentum(self, mass: T) -> FourMomentum<T, T>
    where
        T: Mul<T> + Add<T> + std::ops::Add<Output = T> + Float,
    {
        FourMomentum::new_on_shell(self, mass)
    }
}

#[derive(Default, Debug, PartialEq, Eq, Clone, Copy)]
pub struct FourMomentum<T, U> {
    pub energy: Energy<U>,
    pub three_momentum: ThreeMomentum<T>,
}

impl<T, U> FourMomentum<T, U> {
    pub fn new(energy: Energy<U>, three_momentum: ThreeMomentum<T>) -> Self {
        FourMomentum {
            energy,
            three_momentum,
        }
    }
}
impl<T> FourMomentum<T, T> {
    pub fn new_on_shell(three_momentum: ThreeMomentum<T>, mass: T) -> Self
    where
        T: Mul<T> + Add<T> + std::ops::Add<Output = T> + Float,
    {
        let energy = three_momentum.on_shell_energy(mass);
        FourMomentum {
            energy,
            three_momentum,
        }
    }

    pub fn into_dense(self, index: AbstractIndex) -> DenseTensor<T, VecStructure>
    where
        T: Clone,
    {
        let structure = VecStructure::new(vec![(index, Representation::Lorentz(4.into())).into()]);
        DenseTensor::from_data(
            &[
                self.energy.value,
                self.three_momentum.px,
                self.three_momentum.py,
                self.three_momentum.pz,
            ],
            structure,
        )
        .unwrap()
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
        let energy = self.energy.value.to_polynomial(&RationalField::new(), None);

        let px: MultivariatePolynomial<RationalField, _> =
            Atom::new_num(self.three_momentum.px).to_polynomial(&RationalField::new(), None);
        let py = Atom::new_num(self.three_momentum.py).to_polynomial(&RationalField::new(), None);
        let pz = Atom::new_num(self.three_momentum.pz).to_polynomial(&RationalField::new(), None);

        DenseTensor::from_data(&[energy, px, py, pz], structure).unwrap()
    }
}

impl<T: Display, U: Display> Display for FourMomentum<T, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}, {}", self.energy, self.three_momentum)
    }
}

impl<T: LowerExp, U: LowerExp> LowerExp for FourMomentum<T, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:e}, {:e}", self.energy, self.three_momentum)
    }
}
