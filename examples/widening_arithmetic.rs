// use std::ops::Mul;

// use num::Complex;

// pub trait SmallestUpgrade<T> {
//     type LCM;
//     fn upgrade(self) -> Self::LCM;
// }

// pub struct Up<T> {
//     up: T,
// }

// impl<T> Up<T> {
//     pub fn new(up: T) -> Self {
//         Up { up }
//     }
// }

// impl<T, U, Out> Mul<Up<U>> for Up<T>
// where
//     T: Mul<U, Output = Out>,
// {
//     type Output = Up<Out>;
//     fn mul(self, rhs: Up<U>) -> Self::Output {
//         Up {
//             up: self.up * rhs.up,
//         }
//     }
// }

// impl<T, U> SmallestUpgrade<Up<T>> for Up<U>
// where
//     T: From<U>,
// {
//     type LCM = Up<T>;
//     fn upgrade(self) -> Self::LCM {
//         T::from(self.up)
//     }
// }

// // impl SmallestUpgrade<Up<f64>> for Up<Complex<f64>> {
// //     type LCM = Up<Complex<f64>>;
// //     fn upgrade(self) -> Self::LCM {
// //         self
// //     }
// // }

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

// fn main() {
//     let a = Up::new(1.0);
//     let b = Up::new(2.0);
//     let c = a * b;

//     let ca = Up::new(Complex::new(1.0, 0.0));
//     let f = UpgradingMul::up_mul(b, ca);
// }

fn main() {}
