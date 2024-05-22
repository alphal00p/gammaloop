use std::{
    fmt::{Display, LowerExp},
    ops::{Add, AddAssign, Mul, Sub, SubAssign},
};

use num::Float;

use crate::tensor::{AbstractIndex, DenseTensor, Representation, VecStructure};

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

impl<T: Clone> From<(ThreeMomentum<T>, AbstractIndex)> for DenseTensor<T, VecStructure> {
    fn from(value: (ThreeMomentum<T>, AbstractIndex)) -> Self {
        let (three_momentum, index) = value;
        let structure =
            VecStructure::new(vec![(index, Representation::Euclidean(3.into())).into()]);
        DenseTensor::from_data(
            &[three_momentum.px, three_momentum.py, three_momentum.pz],
            structure,
        )
        .unwrap()
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
}

#[derive(Default, Debug, PartialEq, Eq, Clone, Copy)]
pub struct FourMomentum<T> {
    pub energy: Energy<T>,
    pub three_momentum: ThreeMomentum<T>,
}

impl<T> FourMomentum<T> {
    pub fn new(energy: Energy<T>, three_momentum: ThreeMomentum<T>) -> Self {
        FourMomentum {
            energy,
            three_momentum,
        }
    }

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
}

impl<T: Display> Display for FourMomentum<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}, {}", self.energy, self.three_momentum)
    }
}

impl<T: LowerExp> LowerExp for FourMomentum<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:e}, {:e}", self.energy, self.three_momentum)
    }
}

impl<T: Clone> From<(FourMomentum<T>, AbstractIndex)> for DenseTensor<T, VecStructure> {
    fn from(value: (FourMomentum<T>, AbstractIndex)) -> Self {
        let (four_mom, index) = value;
        let structure = VecStructure::new(vec![(index, Representation::Lorentz(4.into())).into()]);
        DenseTensor::from_data(
            &[
                four_mom.energy.value,
                four_mom.three_momentum.px,
                four_mom.three_momentum.py,
                four_mom.three_momentum.pz,
            ],
            structure,
        )
        .unwrap()
    }
}
