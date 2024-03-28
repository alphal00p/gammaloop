use std::borrow::Cow;
use std::ops::Mul;

use duplicate::duplicate;

use symbolica::domains::float::Complex;

use symbolica::domains::float::Real;
use symbolica::representations::Atom;

use symbolica::state::State;

#[macro_export]
macro_rules! forward_ref_bino {
    (impl $imp:ident, $method:ident for $t:ty, $u:ty,$out:ty) => {
        impl<'a> $imp<$u> for &'a $t {
            type Output = $out;

            #[inline]
            fn $method(self, other: $u) -> Option<Self::Output> {
                $imp::$method(self, &other)
            }
        }

        impl $imp<&$u> for $t {
            type Output = $out;

            #[inline]
            fn $method(self, other: &$u) -> Option<Self::Output> {
                $imp::$method(&self, other)
            }
        }

        impl $imp<$u> for $t {
            type Output = $out;

            #[inline]
            fn $method(self, other: $u) -> Option<Self::Output> {
                $imp::$method(&self, &other)
            }
        }
    };
}
pub trait SmallestUpgrade<T> {
    type LCM;
    fn upgrade(self) -> Self::LCM;
}

pub trait TryFromUpgrade<T> {
    fn try_from_upgrade(value: T) -> Option<Self>
    where
        Self: Sized;
}

pub trait TrySmallestUpgrade<T> {
    type LCM: Clone;

    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>;
}

impl<T, U> TryFromUpgrade<T> for U
where
    T: TrySmallestUpgrade<U, LCM = U>,
    U: Clone,
{
    fn try_from_upgrade(value: T) -> Option<Self> {
        let cow = value.try_upgrade()?;
        Some(cow.into_owned())
    }
}

pub trait TryIntoUpgrade<T> {
    fn try_into_upgrade(self) -> Option<T>;
}

impl<T, U> TryIntoUpgrade<U> for T
where
    U: TryFromUpgrade<T>,
{
    fn try_into_upgrade(self) -> Option<U> {
        U::try_from_upgrade(self)
    }
}

// duplicate! {
//     [ num;
//     [f64] ;
//     [i32] ;]

// impl TrySmallestUpgrade<num> for Complex<num> {
//     type LCM = Complex<num>;

//     fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
//         where
//             Self::LCM: Clone {
//         Some(Cow::Borrowed(self))
//     }
// }

// impl TrySmallestUpgrade<Complex<num>> for num {
//     type LCM = Complex<num>;

//     fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
//         where
//             Self::LCM: Clone {
//         Some(Cow::Owned(Complex::from(*self)))
//     }
// }
// }

// impl<T, U> TrySmallestUpgrade<U> for T
// where
//     T: Borrow<T>,
//     U: Borrow<T>,
// {
//     type LCM = T;

//     fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
//     where
//         Self::LCM: Clone,
//     {
//         Some(Cow::Borrowed(self))
//     }
// } can't do this because of future impls GRR.

duplicate! {
    [smaller larger;
    [i16] [i16];
    [i32] [i32];
    [f64] [f64];
    [Complex<f64>] [Complex<f64>];
    [f64] [Complex<f64>];
    [f64] [Atom];
    [Atom] [Atom];
    [Complex<f64>] [Atom];
    [f32] [f64];
    [i32] [f64];]

impl TrySmallestUpgrade<smaller> for larger {
    type LCM = larger;
    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
        where
            Self::LCM: Clone {
        Some(Cow::Borrowed(self))
    }
}

impl<'a> TrySmallestUpgrade<&'a smaller> for larger {
    type LCM = larger;
    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
        where
            Self::LCM: Clone {
        Some(Cow::Borrowed(self))
    }
}

impl<'a,'b> TrySmallestUpgrade<&'a smaller> for &'b larger {
    type LCM = larger;
    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
        where
            Self::LCM: Clone {
        Some(Cow::Borrowed(*self))
    }
}

impl<'b> TrySmallestUpgrade<smaller> for &'b larger {
    type LCM = larger;
    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
        where
            Self::LCM: Clone {
        Some(Cow::Borrowed(*self))
    }
}
}

duplicate! {
    [smaller larger;
    [f32] [f64];
    [i32] [f64];]

impl TrySmallestUpgrade<larger> for smaller {
    type LCM = larger;


    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
        where
            Self::LCM: Clone {
        Some(Cow::Owned(larger::from(*self)))
    }
}

}

impl<T> TrySmallestUpgrade<Complex<T>> for T
where
    T: Real,
{
    type LCM = Complex<T>;

    fn try_upgrade(&self) -> Option<Cow<Self::LCM>> {
        let new = Complex::new(*self, T::zero());
        Some(Cow::Owned(new))
    }
}

impl TrySmallestUpgrade<Atom> for f64 {
    type LCM = Atom;
    fn try_upgrade(&self) -> Option<Cow<Self::LCM>> {
        let rugrat = rug::Rational::from_f64(*self)?;
        let natrat = symbolica::domains::rational::Rational::from_large(rugrat);
        let symrat = Atom::new_num(symbolica::coefficient::Coefficient::from(natrat));

        Some(Cow::Owned(symrat))
    }
}

impl TrySmallestUpgrade<Atom> for Complex<f64> {
    type LCM = Atom;
    fn try_upgrade(&self) -> Option<Cow<Self::LCM>> {
        let real: Cow<'_, Atom> = <f64 as TrySmallestUpgrade<Atom>>::try_upgrade(&self.re)?;
        let imag: Cow<'_, Atom> = <f64 as TrySmallestUpgrade<Atom>>::try_upgrade(&self.im)?;
        let i = Atom::new_var(State::I);
        let symrat = (i * imag.as_ref()) + real.as_ref();

        Some(Cow::Owned(symrat))
    }
}

duplicate! {
    [smaller larger;
    [f64] [Atom];
    [Complex<f64>] [Atom];
    [f64][Complex<f64>];
    [f32] [f64];
    [i32] [f64];]
impl<'a> TrySmallestUpgrade<&'a larger> for smaller {
    type LCM = larger;
    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
        where
            Self::LCM: Clone {
        <smaller as TrySmallestUpgrade<larger>>::try_upgrade(self)
    }
}

impl<'a,'b> TrySmallestUpgrade<&'a larger> for &'b smaller {
    type LCM = larger;
    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
        where
            Self::LCM: Clone {
       <smaller as TrySmallestUpgrade<larger>>::try_upgrade(*self)
    }}

impl<'b> TrySmallestUpgrade<larger> for &'b smaller {
    type LCM = larger;
    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
        where
            Self::LCM: Clone {
        <smaller as TrySmallestUpgrade<larger>>::try_upgrade(*self)
    }
}

}

pub trait FallibleMul<T> {
    type Output;
    fn mul_fallible(self, rhs: T) -> Option<Self::Output>;
}
impl<T, U> FallibleMul<T> for U
where
    U: TrySmallestUpgrade<T>,
    T: TrySmallestUpgrade<U, LCM = <U as TrySmallestUpgrade<T>>::LCM>,
    U::LCM: Clone,
    for<'a, 'b> &'a U::LCM: std::ops::Mul<&'b U::LCM, Output = U::LCM>,
{
    type Output = U::LCM;

    fn mul_fallible(self, rhs: T) -> Option<Self::Output> {
        let lhs = self.try_upgrade()?;
        let rhs = rhs.try_upgrade()?;
        Some(lhs.as_ref().mul(rhs.as_ref()))
    }
}

pub trait FallibleAdd<T> {
    type Output;
    fn add_fallible(self, rhs: T) -> Option<Self::Output>;
}

impl<T, U, Out> FallibleAdd<T> for U
where
    U: TrySmallestUpgrade<T, LCM = Out>,
    T: TrySmallestUpgrade<U, LCM = Out>,
    Out: Clone,
    for<'a, 'b> &'a Out: std::ops::Add<&'b Out, Output = Out>,
{
    type Output = Out;

    fn add_fallible(self, rhs: T) -> Option<Self::Output> {
        let lhs = self.try_upgrade()?;
        let rhs = rhs.try_upgrade()?;
        Some(lhs.as_ref() + rhs.as_ref())
    }
}

pub trait FallibleSub<T> {
    type Output;
    fn sub_fallible(self, rhs: T) -> Option<Self::Output>;
}

impl<T, U> FallibleSub<T> for U
where
    U: TrySmallestUpgrade<T>,
    T: TrySmallestUpgrade<U, LCM = <U as TrySmallestUpgrade<T>>::LCM>,
    U::LCM: Clone,
    for<'a, 'b> &'a U::LCM: std::ops::Sub<&'b U::LCM, Output = U::LCM>,
{
    type Output = U::LCM;

    fn sub_fallible(self, rhs: T) -> Option<Self::Output> {
        let lhs = self.try_upgrade()?;
        let rhs = rhs.try_upgrade()?;
        Some(lhs.as_ref() - rhs.as_ref())
    }
}

pub trait FallibleAddAssign<T> {
    fn add_assign_fallible(&mut self, rhs: T);
}

impl<T, U> FallibleAddAssign<T> for U
where
    U: TrySmallestUpgrade<T, LCM = U>,
    T: TrySmallestUpgrade<U, LCM = U>,
    U::LCM: Clone,
    for<'a, 'b> &'a U::LCM: std::ops::Add<&'b U::LCM, Output = U::LCM>,
{
    fn add_assign_fallible(&mut self, rhs: T) {
        let lhs = self.try_upgrade().unwrap();
        let rhs = rhs.try_upgrade().unwrap();
        let out = lhs.as_ref() + rhs.as_ref();
        *self = out;
    }
}

pub trait FallibleSubAssign<T> {
    fn sub_assign_fallible(&mut self, rhs: T);
}

impl<T, U> FallibleSubAssign<T> for U
where
    U: TrySmallestUpgrade<T, LCM = U>,
    T: TrySmallestUpgrade<U, LCM = U>,
    U::LCM: Clone,
    for<'a, 'b> &'a U::LCM: std::ops::Sub<&'b U::LCM, Output = U::LCM>,
{
    fn sub_assign_fallible(&mut self, rhs: T) {
        let lhs = self.try_upgrade().unwrap();
        let rhs = rhs.try_upgrade().unwrap();
        let out = lhs.as_ref() - rhs.as_ref();
        *self = out;
    }
}

// impl<T, U> SmallestUpgrade<U> for T
// where
//     U: From<T>,
// {
//     type LCM = U;
//     fn upgrade(self) -> Self::LCM {
//         U::from(self)
//     }
// }

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
