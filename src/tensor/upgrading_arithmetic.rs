use duplicate::duplicate;
use duplicate::duplicate_item;
use num::Complex;

pub trait SmallestUpgrade<T> {
    type LCM;
    fn upgrade(self) -> Self::LCM;
}

// impl<T, U> SmallestUpgrade<Up<T>> for Up<U>
// where
//     T: From<U>,
// {
//     type LCM = T;
//     fn upgrade(self) -> Self::LCM {
//         T::from(self.up)
//     }
// } We can't do this because of possible future impls, means that any specialization is forbidden

// impl<T> SmallestUpgrade<T> for T {
//     type LCM = T;
//     fn upgrade(self) -> Self::LCM {
//         self
//     }
// } We don't want this, so that we can specialize binary operations

// impl<T, U> SmallestUpgrade<T> for U
// where
//     T: SmallestUpgrade<U, LCM = U>,
// {
//     type LCM = U;
//     fn upgrade(self) -> Self::LCM {
//         self
//     }
// } This should work but doesn't

// duplicate! {
//     [
//   U T;
//   [ f64 ]    [Complex<f64>];
//   [ &f64]    [Complex<f64>];
// ]
// impl SmallestUpgrade<U> for T {
//     type LCM = T;
//     fn upgrade(self) -> Self::LCM {
//         self
//     }
// }
// impl<'a> SmallestUpgrade<U> for &'a T {
//     type LCM = &'a T;
//     fn upgrade(self) -> Self::LCM {
//         self
//     }
// }
// }

// duplicate! {
//     [
//   U T;
//   [ f64 ]    [Complex<f64>];
//   [ f64]    [&Complex<f64>];
// ]
// impl SmallestUpgrade<T> for U {
//     type LCM = Complex<f64>;
//     fn upgrade(self) -> Self::LCM {
//         println!("U -> T");
//         Complex::new(self, 0.0)
//     }
// }
// impl<'a> SmallestUpgrade<T> for &'a U {
//     type LCM = Complex<f64>;
//     fn upgrade(self) -> Self::LCM {
//         println!("U -> T");
//         Complex::new(*self, 0.0)
//     }
// }
// }

// pub trait UpgradingMul<T> {
//     type Output;
//     fn up_mul(self, rhs: T) -> Self::Output;
// }

// impl<T, U, Up1, Up2, Out> UpgradingMul<U> for T
// where
//     T: SmallestUpgrade<U, LCM = Up1>,
//     U: SmallestUpgrade<T, LCM = Up2>,
//     Up1: std::ops::Mul<Up2, Output = Out>,
// {
//     type Output = Out;
//     fn up_mul(self, rhs: U) -> Self::Output {
//         self.upgrade() * rhs.upgrade()
//     }
// }

// impl<T, U, Out> UpgradingMul<U> for T
// where
//     T: std::ops::Mul<U, Output = Out>,
// {
//     type Output = Out;
//     fn up_mul(self, rhs: U) -> Self::Output {
//         self * rhs
//     }
// }

// pub trait UpgradingAdd<T> {
//     type Output;
//     fn up_add(self, rhs: T) -> Self::Output;
// }

// impl<T, U, Up1, Up2, Out> UpgradingAdd<U> for T
// where
//     T: SmallestUpgrade<U, LCM = Up1>,
//     U: SmallestUpgrade<T, LCM = Up2>,
//     Up1: std::ops::Add<Up2, Output = Out>,
// {
//     type Output = Out;
//     fn up_add(self, rhs: U) -> Self::Output {
//         self.upgrade() + rhs.upgrade()
//     }
// }
