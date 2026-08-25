use std::{
    fmt::{Display, LowerExp},
    ops::{Add, AddAssign, Mul, Neg, Sub, SubAssign},
};

use serde::{Deserialize, Serialize};

use crate::KinematicScalar;

const LARGE_RAPIDITY_SENTINEL: usize = 100_000;

/// The temporal component of a four-momentum.
#[repr(transparent)]
#[derive(Debug, Default, PartialEq, Eq, Clone, Copy, Hash, Serialize, Deserialize)]
pub struct Energy<T> {
    /// The energy value.
    pub value: T,
}

impl<T> Energy<T> {
    /// Wrap an energy value.
    pub const fn new(value: T) -> Self {
        Self { value }
    }

    /// Unwrap the energy value.
    pub fn into_inner(self) -> T {
        self.value
    }

    /// Transform the energy's scalar type.
    pub fn map<U>(self, map: impl FnOnce(T) -> U) -> Energy<U> {
        Energy::new(map(self.value))
    }

    /// Transform a borrowed energy value.
    pub fn map_ref<U>(&self, map: impl FnOnce(&T) -> U) -> Energy<U> {
        Energy::new(map(&self.value))
    }
}

impl<T> From<T> for Energy<T> {
    fn from(value: T) -> Self {
        Self::new(value)
    }
}

impl<T: Display> Display for Energy<T> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "E: {}", self.value)
    }
}

impl<T: LowerExp> LowerExp for Energy<T> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "E: {:e}", self.value)
    }
}

impl<T> Add for Energy<T>
where
    T: Add<Output = T>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.value + rhs.value)
    }
}

impl<T> AddAssign for Energy<T>
where
    T: AddAssign,
{
    fn add_assign(&mut self, rhs: Self) {
        self.value += rhs.value;
    }
}

impl<T> Sub for Energy<T>
where
    T: Sub<Output = T>,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.value - rhs.value)
    }
}

impl<T> SubAssign for Energy<T>
where
    T: SubAssign,
{
    fn sub_assign(&mut self, rhs: Self) {
        self.value -= rhs.value;
    }
}

impl<T> Neg for Energy<T>
where
    T: Neg<Output = T>,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.value)
    }
}

/// A Cartesian spatial momentum `(px, py, pz)`.
#[derive(Debug, Default, PartialEq, Eq, Clone, Copy, Hash, Serialize, Deserialize)]
pub struct ThreeMomentum<T> {
    /// Momentum along the x-axis.
    pub px: T,
    /// Momentum along the y-axis.
    pub py: T,
    /// Momentum along the z-axis.
    pub pz: T,
}

impl<T> ThreeMomentum<T> {
    /// Construct a Cartesian spatial momentum.
    pub const fn new(px: T, py: T, pz: T) -> Self {
        Self { px, py, pz }
    }

    /// Transform the momentum's scalar type.
    pub fn map<U>(self, mut map: impl FnMut(T) -> U) -> ThreeMomentum<U> {
        ThreeMomentum::new(map(self.px), map(self.py), map(self.pz))
    }

    /// Transform borrowed components into another scalar type.
    pub fn map_ref<U>(&self, mut map: impl FnMut(&T) -> U) -> ThreeMomentum<U> {
        ThreeMomentum::new(map(&self.px), map(&self.py), map(&self.pz))
    }

    /// Borrow the components in Cartesian order.
    pub fn components(&self) -> [&T; 3] {
        [&self.px, &self.py, &self.pz]
    }

    /// Compute the Euclidean norm squared without restricting the scalar type to floats.
    pub fn norm_squared(&self) -> T
    where
        T: Clone + Mul<Output = T> + Add<Output = T>,
    {
        self.px.clone() * self.px.clone()
            + self.py.clone() * self.py.clone()
            + self.pz.clone() * self.pz.clone()
    }

    /// Compute the Euclidean dot product.
    pub fn dot(&self, rhs: &Self) -> T
    where
        T: Clone + Mul<Output = T> + Add<Output = T>,
    {
        self.px.clone() * rhs.px.clone()
            + self.py.clone() * rhs.py.clone()
            + self.pz.clone() * rhs.pz.clone()
    }

    /// Compute the Cartesian cross product.
    pub fn cross(&self, rhs: &Self) -> Self
    where
        T: Clone + Mul<Output = T> + Sub<Output = T>,
    {
        Self {
            px: self.py.clone() * rhs.pz.clone() - self.pz.clone() * rhs.py.clone(),
            py: self.pz.clone() * rhs.px.clone() - self.px.clone() * rhs.pz.clone(),
            pz: self.px.clone() * rhs.py.clone() - self.py.clone() * rhs.px.clone(),
        }
    }
}

impl<T: KinematicScalar> ThreeMomentum<T> {
    /// Construct the zero vector using the active scalar's precision.
    pub fn zero(representative: &T) -> Self {
        let zero = representative.zero();
        Self::new(zero.clone(), zero.clone(), zero)
    }

    /// Compute the Euclidean norm.
    pub fn norm(&self) -> T {
        self.norm_squared().sqrt()
    }

    /// Compute transverse momentum squared.
    pub fn pt_squared(&self) -> T {
        self.px.clone() * self.px.clone() + self.py.clone() * self.py.clone()
    }

    /// Compute transverse momentum.
    pub fn pt(&self) -> T {
        self.pt_squared().sqrt()
    }

    /// Return the azimuthal angle in `[0, 2*pi)`.
    pub fn phi(&self) -> T {
        let zero = self.px.zero();
        let two_pi = self.px.pi() * self.px.from_usize(2);
        let mut phi = self.py.atan2(&self.px);
        if phi < zero {
            phi += two_pi.clone();
        }
        if phi >= two_pi {
            phi -= self.px.pi() * self.px.from_usize(2);
        }
        phi
    }

    /// Compute pseudorapidity, using a finite sentinel on the beam axis.
    pub fn pseudorapidity(&self) -> T {
        let pt = self.pt();
        let zero = pt.zero();
        if pt == zero {
            let abs_pz = self.pz.norm();
            if abs_pz == zero {
                return zero;
            }
            let sentinel = abs_pz.from_usize(LARGE_RAPIDITY_SENTINEL) + abs_pz;
            return if self.pz > self.pz.zero() {
                sentinel
            } else {
                -sentinel
            };
        }

        let theta = pt.atan2(&self.pz);
        let two = theta.from_usize(2);
        -(theta / two).tan().log()
    }

    /// Compute the wrapped absolute azimuthal separation in `[0, pi]`.
    pub fn delta_phi(&self, rhs: &Self) -> T {
        wrapped_delta_phi(self.phi(), rhs.phi())
    }

    /// Compute separation in pseudorapidity and azimuth.
    pub fn delta_r(&self, rhs: &Self) -> T {
        let delta_eta = self.pseudorapidity() - rhs.pseudorapidity();
        let delta_phi = self.delta_phi(rhs);
        (delta_eta.clone() * delta_eta + delta_phi.clone() * delta_phi).sqrt()
    }

    /// Compute the positive on-shell energy for this momentum and an optional mass.
    pub fn on_shell_energy(&self, mass: Option<&T>) -> Energy<T> {
        let mut energy_squared = self.norm_squared();
        if let Some(mass) = mass {
            energy_squared += mass.clone() * mass.clone();
        }
        Energy::new(energy_squared.sqrt())
    }

    /// Promote this spatial momentum to a positive-energy on-shell four-momentum.
    pub fn on_shell(self, mass: Option<&T>) -> FourMomentum<T> {
        let temporal = self.on_shell_energy(mass);
        FourMomentum::new(temporal, self)
    }
}

impl<T> From<[T; 3]> for ThreeMomentum<T> {
    fn from([px, py, pz]: [T; 3]) -> Self {
        Self::new(px, py, pz)
    }
}

impl<T> From<ThreeMomentum<T>> for [T; 3] {
    fn from(momentum: ThreeMomentum<T>) -> Self {
        [momentum.px, momentum.py, momentum.pz]
    }
}

impl<T> From<(T, T, T)> for ThreeMomentum<T> {
    fn from((px, py, pz): (T, T, T)) -> Self {
        Self::new(px, py, pz)
    }
}

impl<T> From<ThreeMomentum<T>> for (T, T, T) {
    fn from(momentum: ThreeMomentum<T>) -> Self {
        (momentum.px, momentum.py, momentum.pz)
    }
}

impl<T> IntoIterator for ThreeMomentum<T> {
    type Item = T;
    type IntoIter = std::array::IntoIter<T, 3>;

    fn into_iter(self) -> Self::IntoIter {
        [self.px, self.py, self.pz].into_iter()
    }
}

impl<'a, T> IntoIterator for &'a ThreeMomentum<T> {
    type Item = &'a T;
    type IntoIter = std::array::IntoIter<&'a T, 3>;

    fn into_iter(self) -> Self::IntoIter {
        [&self.px, &self.py, &self.pz].into_iter()
    }
}

impl<T> Add for ThreeMomentum<T>
where
    T: Add<Output = T>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.px + rhs.px, self.py + rhs.py, self.pz + rhs.pz)
    }
}

impl<T> Add<&Self> for ThreeMomentum<T>
where
    T: for<'a> Add<&'a T, Output = T>,
{
    type Output = Self;

    fn add(self, rhs: &Self) -> Self::Output {
        Self::new(self.px + &rhs.px, self.py + &rhs.py, self.pz + &rhs.pz)
    }
}

impl<T> AddAssign for ThreeMomentum<T>
where
    T: AddAssign,
{
    fn add_assign(&mut self, rhs: Self) {
        self.px += rhs.px;
        self.py += rhs.py;
        self.pz += rhs.pz;
    }
}

impl<T> Sub for ThreeMomentum<T>
where
    T: Sub<Output = T>,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.px - rhs.px, self.py - rhs.py, self.pz - rhs.pz)
    }
}

impl<T> SubAssign for ThreeMomentum<T>
where
    T: SubAssign,
{
    fn sub_assign(&mut self, rhs: Self) {
        self.px -= rhs.px;
        self.py -= rhs.py;
        self.pz -= rhs.pz;
    }
}

impl<T> Neg for ThreeMomentum<T>
where
    T: Neg<Output = T>,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.px, -self.py, -self.pz)
    }
}

impl<T> Neg for &ThreeMomentum<T>
where
    T: Clone + Neg<Output = T>,
{
    type Output = ThreeMomentum<T>;

    fn neg(self) -> Self::Output {
        ThreeMomentum::new(-self.px.clone(), -self.py.clone(), -self.pz.clone())
    }
}

impl<T> Mul<T> for ThreeMomentum<T>
where
    T: Clone + Mul<Output = T>,
{
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        Self::new(self.px * rhs.clone(), self.py * rhs.clone(), self.pz * rhs)
    }
}

impl<T: Display> Display for ThreeMomentum<T> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "px: {}, py: {}, pz: {}",
            self.px, self.py, self.pz
        )
    }
}

impl<T: LowerExp> LowerExp for ThreeMomentum<T> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "px: {:e}, py: {:e}, pz: {:e}",
            self.px, self.py, self.pz
        )
    }
}

/// A contravariant four-momentum `(E, px, py, pz)` with mostly-minus metric.
#[derive(Debug, Default, PartialEq, Eq, Clone, Copy, Hash, Serialize, Deserialize)]
pub struct FourMomentum<T> {
    /// Temporal component.
    pub temporal: Energy<T>,
    /// Cartesian spatial component.
    pub spatial: ThreeMomentum<T>,
}

impl<T> FourMomentum<T> {
    /// Construct a four-momentum from temporal and spatial components.
    pub const fn new(temporal: Energy<T>, spatial: ThreeMomentum<T>) -> Self {
        Self { temporal, spatial }
    }

    /// Construct a four-momentum from `(E, px, py, pz)`.
    pub const fn from_args(energy: T, px: T, py: T, pz: T) -> Self {
        Self::new(Energy::new(energy), ThreeMomentum::new(px, py, pz))
    }

    /// Transform the momentum's scalar type.
    pub fn map<U>(self, mut map: impl FnMut(T) -> U) -> FourMomentum<U> {
        FourMomentum::from_args(
            map(self.temporal.value),
            map(self.spatial.px),
            map(self.spatial.py),
            map(self.spatial.pz),
        )
    }

    /// Transform borrowed components into another scalar type.
    pub fn map_ref<U>(&self, mut map: impl FnMut(&T) -> U) -> FourMomentum<U> {
        FourMomentum::from_args(
            map(&self.temporal.value),
            map(&self.spatial.px),
            map(&self.spatial.py),
            map(&self.spatial.pz),
        )
    }

    /// Borrow the components in `(E, px, py, pz)` order.
    pub fn components(&self) -> [&T; 4] {
        [
            &self.temporal.value,
            &self.spatial.px,
            &self.spatial.py,
            &self.spatial.pz,
        ]
    }

    /// Compute the Minkowski scalar product with mostly-minus metric.
    pub fn dot(&self, rhs: &Self) -> T
    where
        T: Clone + Mul<Output = T> + Add<Output = T> + Sub<Output = T>,
    {
        self.temporal.value.clone() * rhs.temporal.value.clone() - self.spatial.dot(&rhs.spatial)
    }

    /// Compute invariant mass squared with mostly-minus metric.
    pub fn mass_squared(&self) -> T
    where
        T: Clone + Mul<Output = T> + Add<Output = T> + Sub<Output = T>,
    {
        self.dot(self)
    }
}

impl<T: KinematicScalar> FourMomentum<T> {
    /// Construct the zero four-vector using the active scalar's precision.
    pub fn zero(representative: &T) -> Self {
        let zero = representative.zero();
        Self::from_args(zero.clone(), zero.clone(), zero.clone(), zero)
    }

    /// Return the positive square root of invariant mass squared.
    pub fn mass(&self) -> T {
        self.mass_squared().sqrt()
    }

    /// Compute transverse momentum squared.
    pub fn pt_squared(&self) -> T {
        self.spatial.pt_squared()
    }

    /// Compute transverse momentum.
    pub fn pt(&self) -> T {
        self.spatial.pt()
    }

    /// Return the azimuthal angle in `[0, 2*pi)`.
    pub fn phi(&self) -> T {
        self.spatial.phi()
    }

    /// Compute pseudorapidity from the spatial momentum.
    pub fn pseudorapidity(&self) -> T {
        self.spatial.pseudorapidity()
    }

    /// Compute rapidity, using a finite sentinel for a massless beam-axis momentum.
    pub fn rapidity(&self) -> T {
        self.rapidity_with_pt_squared(&self.pt_squared())
    }

    pub(crate) fn rapidity_with_pt_squared(&self, pt_squared: &T) -> T {
        let zero = self.temporal.value.zero();
        let pz = self.spatial.pz.clone();
        let abs_pz = pz.norm();
        let energy = self.temporal.value.clone();

        if *pt_squared == zero.clone() && energy == abs_pz {
            let sentinel = abs_pz.from_usize(LARGE_RAPIDITY_SENTINEL) + abs_pz;
            return if pz >= zero { sentinel } else { -sentinel };
        }

        let spatial_norm_squared = self.spatial.norm_squared();
        let mass_squared = energy.clone() * energy.clone() - spatial_norm_squared;
        let effective_mass_squared = if mass_squared < zero {
            self.temporal.value.zero()
        } else {
            mass_squared
        };
        let energy_plus_abs_pz = energy + abs_pz;
        let two = energy_plus_abs_pz.from_usize(2);
        let half = energy_plus_abs_pz.one() / two;
        let denominator = energy_plus_abs_pz.clone() * energy_plus_abs_pz;
        let mut rapidity =
            ((pt_squared.clone() + effective_mass_squared) / denominator).log() * half;
        if pz > self.spatial.pz.zero() {
            rapidity = -rapidity;
        }
        rapidity
    }

    /// Compute the wrapped absolute azimuthal separation in `[0, pi]`.
    pub fn delta_phi(&self, rhs: &Self) -> T {
        self.spatial.delta_phi(&rhs.spatial)
    }

    /// Compute separation in rapidity and azimuth, as used by hadron-collider jets.
    pub fn delta_r(&self, rhs: &Self) -> T {
        let delta_rapidity = self.rapidity() - rhs.rapidity();
        let delta_phi = self.delta_phi(rhs);
        (delta_rapidity.clone() * delta_rapidity + delta_phi.clone() * delta_phi).sqrt()
    }
}

impl<T> From<[T; 4]> for FourMomentum<T> {
    fn from([energy, px, py, pz]: [T; 4]) -> Self {
        Self::from_args(energy, px, py, pz)
    }
}

impl<T> From<FourMomentum<T>> for [T; 4] {
    fn from(momentum: FourMomentum<T>) -> Self {
        [
            momentum.temporal.value,
            momentum.spatial.px,
            momentum.spatial.py,
            momentum.spatial.pz,
        ]
    }
}

impl<T> From<(T, T, T, T)> for FourMomentum<T> {
    fn from((energy, px, py, pz): (T, T, T, T)) -> Self {
        Self::from_args(energy, px, py, pz)
    }
}

impl<T> From<FourMomentum<T>> for (T, T, T, T) {
    fn from(momentum: FourMomentum<T>) -> Self {
        (
            momentum.temporal.value,
            momentum.spatial.px,
            momentum.spatial.py,
            momentum.spatial.pz,
        )
    }
}

impl<T> IntoIterator for FourMomentum<T> {
    type Item = T;
    type IntoIter = std::array::IntoIter<T, 4>;

    fn into_iter(self) -> Self::IntoIter {
        <[T; 4]>::from(self).into_iter()
    }
}

impl<'a, T> IntoIterator for &'a FourMomentum<T> {
    type Item = &'a T;
    type IntoIter = std::array::IntoIter<&'a T, 4>;

    fn into_iter(self) -> Self::IntoIter {
        self.components().into_iter()
    }
}

impl<T> Add for FourMomentum<T>
where
    T: Add<Output = T>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.temporal + rhs.temporal, self.spatial + rhs.spatial)
    }
}

impl<T> Add<&Self> for FourMomentum<T>
where
    T: for<'a> Add<&'a T, Output = T>,
{
    type Output = Self;

    fn add(self, rhs: &Self) -> Self::Output {
        Self::from_args(
            self.temporal.value + &rhs.temporal.value,
            self.spatial.px + &rhs.spatial.px,
            self.spatial.py + &rhs.spatial.py,
            self.spatial.pz + &rhs.spatial.pz,
        )
    }
}

impl<T> AddAssign for FourMomentum<T>
where
    T: AddAssign,
{
    fn add_assign(&mut self, rhs: Self) {
        self.temporal += rhs.temporal;
        self.spatial += rhs.spatial;
    }
}

impl<T> Sub for FourMomentum<T>
where
    T: Sub<Output = T>,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.temporal - rhs.temporal, self.spatial - rhs.spatial)
    }
}

impl<T> SubAssign for FourMomentum<T>
where
    T: SubAssign,
{
    fn sub_assign(&mut self, rhs: Self) {
        self.temporal -= rhs.temporal;
        self.spatial -= rhs.spatial;
    }
}

impl<T> Neg for FourMomentum<T>
where
    T: Neg<Output = T>,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.temporal, -self.spatial)
    }
}

impl<T> Neg for &FourMomentum<T>
where
    T: Clone + Neg<Output = T>,
{
    type Output = FourMomentum<T>;

    fn neg(self) -> Self::Output {
        FourMomentum::from_args(
            -self.temporal.value.clone(),
            -self.spatial.px.clone(),
            -self.spatial.py.clone(),
            -self.spatial.pz.clone(),
        )
    }
}

impl<T> Mul<T> for FourMomentum<T>
where
    T: Clone + Mul<Output = T>,
{
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        Self::from_args(
            self.temporal.value * rhs.clone(),
            self.spatial.px * rhs.clone(),
            self.spatial.py * rhs.clone(),
            self.spatial.pz * rhs,
        )
    }
}

impl<T: Display> Display for FourMomentum<T> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{}, {}", self.temporal, self.spatial)
    }
}

impl<T: LowerExp> LowerExp for FourMomentum<T> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{:e}, {:e}", self.temporal, self.spatial)
    }
}

pub(crate) fn wrapped_delta_phi<T: KinematicScalar>(lhs: T, rhs: T) -> T {
    let pi = lhs.pi();
    let two_pi = pi.clone() * pi.from_usize(2);
    let delta = (lhs - rhs).norm();
    if delta > pi { two_pi - delta } else { delta }
}

#[cfg(test)]
mod tests {
    use numerica::domains::float::{DoubleFloat, FloatLike, Real, RealLike};

    use super::*;

    fn assert_close(lhs: f64, rhs: f64) {
        assert!((lhs - rhs).abs() < 1.0e-12, "{lhs} != {rhs}");
    }

    #[test]
    fn on_shell_construction_preserves_mass() {
        let spatial = ThreeMomentum::new(3.0, 4.0, 0.0);
        let momentum = spatial.on_shell(Some(&12.0));

        assert_close(momentum.temporal.value, 13.0);
        assert_close(momentum.mass_squared(), 144.0);
    }

    #[test]
    fn azimuthal_separation_wraps_at_two_pi() {
        let epsilon = 0.01_f64;
        let left = ThreeMomentum::new(epsilon.cos(), epsilon.sin(), 0.0);
        let angle = 0.0_f64.pi() * 2.0 - epsilon;
        let right = ThreeMomentum::new(angle.cos(), angle.sin(), 0.0);

        assert_close(left.delta_phi(&right), 2.0 * epsilon);
    }

    #[test]
    fn vector_arithmetic_supports_symbolica_high_precision_scalars() {
        let zero = DoubleFloat::default();
        let one = zero.one();
        let momentum = FourMomentum::from_args(
            one.from_usize(5),
            one.from_usize(3),
            one.from_usize(4),
            zero,
        );

        assert_close(momentum.mass_squared().to_f64(), 0.0);
    }
}
