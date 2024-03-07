use std::borrow::Cow;
use std::ops::Mul;

use duplicate::duplicate;

use num::traits::Zero;
use num::Complex;

use symbolica::representations::Atom;

use symbolica::state::State;
use symbolica::state::Workspace;

/// For a given binary operation on &T and &U, we can implement the same operation on any combination of &T,T and &U,U.
#[macro_export]
macro_rules! forward_ref_binop {
    (impl $imp:ident, $method:ident for $t:ty, $u:ty,$out:ty) => {
        impl<'a> $imp<$u> for &'a $t {
            type Output = $out;

            #[inline]
            fn $method(self, other: $u, ws: &Workspace, state: &State) -> Option<Self::Output> {
                $imp::$method(self, &other, ws, state)
            }
        }

        impl $imp<&$u> for $t {
            type Output = $out;

            #[inline]
            fn $method(self, other: &$u, ws: &Workspace, state: &State) -> Option<Self::Output> {
                $imp::$method(&self, other, ws, state)
            }
        }

        impl $imp<$u> for $t {
            type Output = $out;

            #[inline]
            fn $method(self, other: $u, ws: &Workspace, state: &State) -> Option<Self::Output> {
                $imp::$method(&self, &other, ws, state)
            }
        }
    };
}

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

pub trait TrySmallestUpgrade<T> {
    type LCM;

    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
    where
        Self::LCM: Clone;
}

duplicate! {
    [ num;
    [f64] ;
    [i32] ;]

impl TrySmallestUpgrade<num> for Complex<num> {
    type LCM = Complex<num>;


    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
        where
            Self::LCM: Clone {
        Some(Cow::Borrowed(self))
    }
}

impl TrySmallestUpgrade<Complex<num>> for num {
    type LCM = Complex<num>;


    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
        where
            Self::LCM: Clone {
        Some(Cow::Owned(Complex::from(*self)))
    }
}
}

duplicate! {
    [smaller larger;
    [f64] [f64];
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
        Some(Cow::Borrowed(self))
    }
}


impl<'b> TrySmallestUpgrade<smaller> for &'b larger {
    type LCM = larger;


    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
        where
            Self::LCM: Clone {
        Some(Cow::Borrowed(self))
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

impl TrySmallestUpgrade<Atom> for f64 {
    type LCM = Atom;
    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
    where
        Self::LCM: Clone,
    {
        let rugrat = rug::Rational::from_f64(*self)?;
        let natrat = symbolica::domains::rational::Rational::from_large(rugrat);
        let symrat = Atom::new_num(symbolica::coefficient::Coefficient::from(natrat));

        Some(Cow::Owned(symrat))
    }
}

impl TrySmallestUpgrade<Atom> for Complex<f64> {
    type LCM = Atom;
    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
    where
        Self::LCM: Clone,
    {
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

impl<T, U> FallibleAdd<T> for U
where
    U: TrySmallestUpgrade<T>,
    T: TrySmallestUpgrade<U, LCM = <U as TrySmallestUpgrade<T>>::LCM>,
    U::LCM: Clone,
    for<'a, 'b> &'a U::LCM: std::ops::Add<&'b U::LCM, Output = U::LCM>,
{
    type Output = U::LCM;

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

#[test]
fn test_fallible_mul() {
    let a: i32 = 4;
    let b: f64 = 4.;
    let mut c = a.mul_fallible(b).unwrap();
    c.add_assign_fallible(&a);
    let d: Option<f64> = b.mul_fallible(&a);
    let a: &i32 = &a;
    let e: Option<f64> = a.mul_fallible(&b);
    assert_eq!(c, 20.);
    assert_eq!(d, Some(16.));
    assert_eq!(e, Some(16.));

    let atom = &Atom::parse("a", &mut State::get_global_state().write().unwrap()).unwrap();

    let f = atom.mul_fallible(4.);

    // print!("{}", f.unwrap());

    let g = Complex::new(0.1, 3.);

    let mut h = atom.sub_fallible(g).unwrap();

    h.add_assign_fallible(atom);
    let f = atom.mul_fallible(atom);

    print!("{}", h);
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

pub trait SmallestUpgradeSymbolic<T> {
    type LCMS;
    fn upgrade_sym(self, ws: &Workspace, state: &State) -> Option<Self::LCMS>;
}

/// If possble turn into an Atom.
pub trait SymbolicInto {
    fn into_sym(self) -> Option<Atom>;
}

/// First turns into a rational, then into an Atom.
impl SymbolicInto for f64 {
    fn into_sym(self) -> Option<Atom> {
        let rugrat = rug::Rational::from_f64(self)?;
        let natrat = symbolica::domains::rational::Rational::from_large(rugrat);
        let symrat = Atom::new_num(symbolica::coefficient::Coefficient::from(natrat));

        Some(symrat)
    }
}

#[allow(clippy::used_underscore_binding)]
impl SymbolicInto for &f64 {
    fn into_sym(self) -> Option<Atom> {
        SymbolicInto::into_sym(*self)
    }
}

/// Uses `f64::into_sym` to turn real and imaginary parts into Atoms, then adds them together, multiplying the imaginary part by i.
impl SymbolicInto for Complex<f64> {
    fn into_sym(self) -> Option<Atom> {
        let real = self.re.into_sym()?;
        let imag = self.im.into_sym()?;
        let i = Atom::new_var(State::I);
        let symrat = (i * &imag) + &real;

        Some(symrat)
    }
}

impl SymbolicInto for &Complex<f64> {
    fn into_sym(self) -> Option<Atom> {
        SymbolicInto::into_sym(*self)
    }
}

/// Smart multiplication, that automatically upgrades types when necessary.
///
/// Normal multiplication is used for types that support it.
pub trait SymbolicMul<T> {
    type Output;
    fn mul_sym(self, rhs: T, ws: &Workspace, state: &State) -> Option<Self::Output>;
}

// impl<T, Out> SymbolicMul<T> for T
// where
//     T: std::ops::Mul<T, Output = Out>,
// {
//     type Output = Out;
//     fn mul_sym(self, rhs: T, _ws: &Workspace, _state: &State) -> Option<Self::Output> {
//         Some(self * rhs)
//     }
// } can't do this because of future impls GRR.

pub trait SymbolicAdd<T> {
    type Output;
    fn add_sym(self, rhs: T, ws: &Workspace, state: &State) -> Option<Self::Output>;
}

// impl<T, U, Out> SymbolicMul<U> for T
// where
//     T: std::ops::Mul<U, Output = Out>,
// {
//     type Output = Out;
//     fn mul_sym(self, rhs: U, _ws: &Workspace, _state: &State) -> Self::Output {
//         self * rhs
//     }
// } can't do this either because of future impls GRR.

pub trait SymbolicSub<T> {
    type Output;
    fn sub_sym(self, other: T, ws: &Workspace, state: &State) -> Option<Self::Output>;
}

pub trait SymbolicAddAssign<T> {
    fn add_assign_sym(&mut self, rhs: T, ws: &Workspace, state: &State);
}

impl SymbolicAddAssign<Atom> for Atom {
    fn add_assign_sym(&mut self, rhs: Atom, ws: &Workspace, state: &State) {
        SymbolicAddAssign::<&Atom>::add_assign_sym(self, &rhs, ws, state);
    }
}

impl SymbolicAddAssign<&Atom> for Atom {
    fn add_assign_sym(&mut self, rhs: &Atom, _ws: &Workspace, _state: &State) {
        let out = self.clone() + rhs;
        *self = out;
    }
}

pub trait SymbolicSubAssign<T> {
    fn sub_assign_sym(&mut self, rhs: T, ws: &Workspace, state: &State);
}

impl SymbolicSubAssign<&Atom> for Atom {
    fn sub_assign_sym(&mut self, rhs: &Atom, _ws: &Workspace, _state: &State) {
        let out = self.clone() - rhs;
        *self = out;
    }
}
impl SymbolicSubAssign<Atom> for Atom {
    fn sub_assign_sym(&mut self, rhs: Atom, ws: &Workspace, state: &State) {
        SymbolicSubAssign::<&Atom>::sub_assign_sym(self, &rhs, ws, state);
    }
}

pub trait SymbolicNeg {
    #[must_use]
    fn neg_sym(self, ws: &Workspace, state: &State) -> Self;
}

impl SymbolicNeg for Atom {
    fn neg_sym(self, _ws: &Workspace, _state: &State) -> Self {
        -self
    }
}

pub trait SymbolicZero {
    fn sym_zero(state: &State, ws: &Workspace) -> Self;
}

impl SymbolicZero for f64 {
    fn sym_zero(_state: &State, _ws: &Workspace) -> Self {
        0.
    }
}

impl SymbolicZero for Complex<f64> {
    fn sym_zero(_state: &State, _ws: &Workspace) -> Self {
        Complex::new(0., 0.)
    }
}

impl SymbolicZero for Atom {
    fn sym_zero(_state: &State, _ws: &Workspace) -> Self {
        Atom::new_num(0)
    }
}

duplicate! {[
  U                 T               Out;
  [ f64]            [ f64]          [ f64];
  [ Complex<f64>]   [Complex<f64>]  [Complex<f64>];
  [ f64]            [Complex<f64>]  [Complex<f64>];
  [ Complex<f64>]   [f64]           [Complex<f64>];
]

impl<'a> SymbolicMul<&U> for &'a T{
    type Output = Out;
    fn mul_sym<'b>(self, rhs: &U, _ws: &Workspace, _state: &State) -> Option<Self::Output> {
        Some(self*rhs)
    }
}

forward_ref_binop!(impl SymbolicMul, mul_sym for T, U, Out);


impl<'a> SymbolicAdd<&U> for &'a T{
    type Output = Out;
    fn add_sym<'b>(self, rhs: &U, _ws: &Workspace, _state: &State) -> Option<Self::Output> {
        Some(self+rhs)
    }
}

forward_ref_binop!(impl SymbolicAdd, add_sym for T, U, Out);



impl<'a> SymbolicSub<&U> for &'a T{
    type Output = Out;
    fn sub_sym<'b>(self, rhs: &U, _ws: &'b Workspace, _state: &'b State) -> Option<Self::Output> {
        Some(self-rhs)
    }
}

forward_ref_binop!(impl SymbolicSub, sub_sym for T, U, Out);


}

impl<'a> SymbolicMul<&Atom> for &'a Atom {
    type Output = Atom;
    fn mul_sym(self, rhs: &Atom, _ws: &Workspace, _state: &State) -> Option<Self::Output> {
        Some(self * rhs)
    }
}

forward_ref_binop!(impl SymbolicMul, mul_sym for Atom, Atom, Atom);

impl<'a> SymbolicAdd<&Atom> for &'a Atom {
    type Output = Atom;
    fn add_sym(self, rhs: &Atom, _ws: &Workspace, _state: &State) -> Option<Self::Output> {
        Some(self + rhs)
    }
}

forward_ref_binop!(impl SymbolicAdd, add_sym for Atom, Atom, Atom);

impl<'a> SymbolicSub<&Atom> for &'a Atom {
    type Output = Atom;
    fn sub_sym(self, rhs: &Atom, _ws: &Workspace, _state: &State) -> Option<Self::Output> {
        Some(self - rhs)
    }
}

forward_ref_binop!(impl SymbolicSub, sub_sym for Atom, Atom, Atom);

duplicate! {
    [ num;
    [f64] ;
    [Complex<f64>] ;]

    impl SymbolicNeg for num{
        fn neg_sym(self, _ws: &Workspace, _state: &State) -> Self {
            -self
        }
    }

    impl SymbolicAddAssign<num> for num{
        fn add_assign_sym(&mut self, rhs: num, _ws: &Workspace, _state: &State) {
            *self+=rhs;
        }
    }

    impl SymbolicAddAssign<&num> for num{
        fn add_assign_sym(&mut self, rhs: &num, _ws: &Workspace, _state: &State) {
            *self+=rhs;
        }
    }

    impl SymbolicSubAssign<num> for num{
        fn sub_assign_sym(&mut self, rhs: num, _ws: &Workspace, _state: &State) {
            *self-=rhs;
        }
    }

    impl SymbolicSubAssign<&num> for num{
        fn sub_assign_sym(&mut self, rhs: &num, _ws: &Workspace, _state: &State) {
            *self-=rhs;
        }
    }

    impl<'a> SymbolicMul<&Atom> for &'a num{
        type Output = Atom;
        fn mul_sym(self, rhs: &Atom, ws: &Workspace, state: &State) -> Option<Self::Output> {
            self.into_sym()?.mul_sym(rhs, ws, state)
        }
    }

    forward_ref_binop! {impl SymbolicMul, mul_sym for num, Atom, Atom}

    impl<'a> SymbolicMul<&num> for &'a Atom{
        type Output = Atom;
        fn mul_sym(self, rhs: &num, ws: &Workspace, state: &State) -> Option<Self::Output> {
            SymbolicMul::<&Atom>::mul_sym(rhs,self, ws, state)
        }
    }

    forward_ref_binop! {impl SymbolicMul, mul_sym for Atom, num, Atom}

    impl<'a> SymbolicAdd<&Atom> for &'a num{
        type Output = Atom;
        fn add_sym(self, rhs: &Atom, ws: &Workspace, state: &State) -> Option<Self::Output> {
            self.into_sym()?.add_sym(rhs, ws, state)
        }
    }

    forward_ref_binop! {impl SymbolicAdd, add_sym for num, Atom, Atom}

    impl<'a> SymbolicAdd<&num> for &'a Atom{
        type Output = Atom;
        fn add_sym(self, rhs: &num, ws: &Workspace, state: &State) -> Option<Self::Output> {
            SymbolicAdd::<&Atom>::add_sym(rhs,self, ws, state)
        }
    }

    forward_ref_binop! {impl SymbolicAdd, add_sym for Atom, num, Atom}

    impl<'a> SymbolicSub<&Atom> for &'a num{
        type Output = Atom;
        fn sub_sym(self, rhs: &Atom, ws: &Workspace, state: &State) -> Option<Self::Output> {
            self.into_sym()?.sub_sym(rhs, ws, state)
        }
    }

    forward_ref_binop! {impl SymbolicSub, sub_sym for num, Atom, Atom}

    impl<'a> SymbolicSub<&num> for &'a Atom{
        type Output = Atom;
        fn sub_sym(self, rhs: &num, ws: &Workspace, state: &State) -> Option<Self::Output> {
            SymbolicSub::<&Atom>::sub_sym(rhs,self, ws, state)
        }
    }

    forward_ref_binop! {impl SymbolicSub, sub_sym for Atom, num, Atom}

}

#[test]
fn symbolic() {
    let mut state = State::get_global_state().write().unwrap();
    let ws: Workspace = Workspace::new();

    let a = Atom::parse("c", &mut state).unwrap();
    let b = Atom::parse("d", &mut state).unwrap();
    let c = a.add_sym(&b, &ws, &state).unwrap();
    let _d = c.add_sym(4., &ws, &state).unwrap();

    let _e = (4.).add_sym(4., &ws, &state);
    let _f = Complex::new(2., 3.).add_sym(4., &ws, &state);
    let _g = Complex::new(2., 3.).add_sym(Complex::new(2., 3.), &ws, &state);
    let _h = (3.).add_sym(Complex::new(1., 2.), &ws, &state);
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
