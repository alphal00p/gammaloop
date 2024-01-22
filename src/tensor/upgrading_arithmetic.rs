use num::Complex;

// pub struct Up<T> {
//     up: T,
// }
pub trait SmallestUpgrade<T> {
    type LCM;
    fn upgrade(self) -> Self::LCM;
}

// impl SmallestUpgrade<A> for B {
//     type LCM = B;
// }

// impl<T, U> SmallestUpgrade<Up<T>> for Up<U>
// where
//     T: From<U>,
// {
//     type LCM = T;
//     fn upgrade(self) -> Self::LCM {
//         T::from(self.up)
//     }
// }

impl<T> SmallestUpgrade<T> for T {
    type LCM = T;
    fn upgrade(self) -> Self::LCM {
        self
    }
}

// impl<T, U> SmallestUpgrade<T> for U
// where
//     T: SmallestUpgrade<U, LCM = U>,
// {
//     type LCM = U;
//     fn upgrade(self) -> Self::LCM {
//         self
//     }
// }
impl SmallestUpgrade<f64> for Complex<f64> {
    type LCM = Complex<f64>;
    fn upgrade(self) -> Self::LCM {
        self
    }
}

impl SmallestUpgrade<&f64> for Complex<f64> {
    type LCM = Complex<f64>;
    fn upgrade(self) -> Self::LCM {
        self
    }
}

impl<'a> SmallestUpgrade<&f64> for &'a Complex<f64> {
    type LCM = &'a Complex<f64>;
    fn upgrade(self) -> Self::LCM {
        self
    }
}

impl<'a> SmallestUpgrade<f64> for &'a Complex<f64> {
    type LCM = &'a Complex<f64>;
    fn upgrade(self) -> Self::LCM {
        self
    }
}

impl SmallestUpgrade<Complex<f64>> for f64 {
    type LCM = Complex<f64>;
    fn upgrade(self) -> Self::LCM {
        Complex::new(self, 0.0)
    }
}

impl<'a> SmallestUpgrade<&Complex<f64>> for &'a f64 {
    type LCM = Complex<f64>;
    fn upgrade(self) -> Self::LCM {
        Complex::new(*self, 0.0)
    }
}

impl<'a> SmallestUpgrade<Complex<f64>> for &'a f64 {
    type LCM = Complex<f64>;
    fn upgrade(self) -> Self::LCM {
        Complex::new(*self, 0.0)
    }
}

impl SmallestUpgrade<&Complex<f64>> for f64 {
    type LCM = Complex<f64>;
    fn upgrade(self) -> Self::LCM {
        Complex::new(self, 0.0)
    }
}

pub trait UpgradingMul<T> {
    type Output;
    fn up_mul(self, rhs: T) -> Self::Output;
}

impl<T, U, Up1, Up2, Out> UpgradingMul<U> for T
where
    T: SmallestUpgrade<U, LCM = Up1>,
    U: SmallestUpgrade<T, LCM = Up2>,
    Up1: std::ops::Mul<Up2, Output = Out>,
{
    type Output = Out;
    fn up_mul(self, rhs: U) -> Self::Output {
        self.upgrade() * rhs.upgrade()
    }
}

pub trait UpgradingAdd<T> {
    type Output;
    fn up_add(self, rhs: T) -> Self::Output;
}

impl<T, U, Up1, Up2, Out> UpgradingAdd<U> for T
where
    T: SmallestUpgrade<U, LCM = Up1>,
    U: SmallestUpgrade<T, LCM = Up2>,
    Up1: std::ops::Add<Up2, Output = Out>,
{
    type Output = Out;
    fn up_add(self, rhs: U) -> Self::Output {
        self.upgrade() + rhs.upgrade()
    }
}

#[test]
fn test_upgrading_mul() {
    let a = 1.0;
    let b = Complex::new(1.0, 1.0);
    let c = UpgradingMul::up_mul(&a, b);
    let d = b.up_mul(a);
    assert_eq!(c, d);
    assert_eq!(c, Complex::new(1.0, 1.0));
}
