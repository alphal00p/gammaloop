use std::mem::transmute;

use ref_ops::{RefAdd, RefDiv, RefMul, RefSub};
use symbolica::{
    coefficient::Coefficient,
    domains::{
        float::{Complex as SymComplex, Float, FloatLike, Real, SingleFloat},
        rational::Rational,
    },
    evaluate::{
        CompileOptions, CompiledComplexEvaluator, CompiledNumber, EvaluatorLoader, ExportNumber,
        ExportSettings, ExpressionEvaluator,
    },
};

use crate::algebra::{
    algebraic_traits::{RefOne, RefZero},
    complex::RealOrComplexRef,
};
use rand::Rng;

use super::{Complex, RealOrComplex};

macro_rules! from_via_rational {
    ($type:ty) => {
        impl From<Complex<$type>> for Coefficient {
            fn from(value: Complex<$type>) -> Self {
                Coefficient::Complex(SymComplex::new(
                    Rational::from(value.re),
                    Rational::from(value.im),
                ))
            }
        }

        impl From<&Complex<$type>> for Coefficient {
            fn from(value: &Complex<$type>) -> Self {
                Coefficient::Complex(SymComplex::new(
                    Rational::from(value.re),
                    Rational::from(value.im),
                ))
            }
        }

        impl From<RealOrComplex<$type>> for Coefficient {
            fn from(value: RealOrComplex<$type>) -> Self {
                match value {
                    RealOrComplex::Real(r) => Coefficient::Complex(Rational::from(r).into()),
                    RealOrComplex::Complex(c) => Coefficient::Complex(SymComplex::new(
                        Rational::from(c.re),
                        Rational::from(c.im),
                    )),
                }
            }
        }

        impl From<&RealOrComplex<$type>> for Coefficient {
            fn from(value: &RealOrComplex<$type>) -> Self {
                match value {
                    RealOrComplex::Real(r) => Coefficient::Complex(Rational::from(*r).into()),
                    RealOrComplex::Complex(c) => Coefficient::Complex(SymComplex::new(
                        Rational::from(c.re),
                        Rational::from(c.im),
                    )),
                }
            }
        }

        impl From<RealOrComplexRef<'_, $type>> for Coefficient {
            fn from(value: RealOrComplexRef<'_, $type>) -> Self {
                match value {
                    RealOrComplexRef::Real(r) => Coefficient::Complex(Rational::from(*r).into()),
                    RealOrComplexRef::Complex(c) => Coefficient::Complex(SymComplex::new(
                        Rational::from(c.re),
                        Rational::from(c.im),
                    )),
                }
            }
        }
    };
}

from_via_rational!(i32);
from_via_rational!(i64);
from_via_rational!(i128);
from_via_rational!(isize);
from_via_rational!(u32);
from_via_rational!(u64);
from_via_rational!(u128);
from_via_rational!(usize);

pub trait ToFloat {
    fn to_float(&self) -> Float;
}

impl ToFloat for f32 {
    fn to_float(&self) -> Float {
        Float::with_val(53, *self as f64)
    }
}
impl ToFloat for f64 {
    fn to_float(&self) -> Float {
        Float::with_val(53, *self)
    }
}

impl<F: ToFloat> From<Complex<F>> for Coefficient {
    fn from(value: Complex<F>) -> Self {
        Coefficient::Float(SymComplex {
            re: value.re.to_float(),
            im: value.im.to_float(),
        })
    }
}

impl<F: ToFloat> From<&Complex<F>> for Coefficient {
    fn from(value: &Complex<F>) -> Self {
        Coefficient::Float(SymComplex {
            re: value.re.to_float(),
            im: value.im.to_float(),
        })
    }
}

impl<F: ToFloat> From<RealOrComplex<F>> for Coefficient {
    fn from(value: RealOrComplex<F>) -> Self {
        match value {
            RealOrComplex::Real(r) => {
                Coefficient::Float(SymComplex::new(r.to_float(), r.to_float().zero()))
            }
            RealOrComplex::Complex(c) => {
                Coefficient::Float(SymComplex::new(c.re.to_float(), c.im.to_float()))
            }
        }
    }
}

impl<F: ToFloat> From<&RealOrComplex<F>> for Coefficient {
    fn from(value: &RealOrComplex<F>) -> Self {
        match value {
            RealOrComplex::Real(r) => {
                Coefficient::Float(SymComplex::new(r.to_float(), r.to_float().zero()))
            }
            RealOrComplex::Complex(c) => {
                Coefficient::Float(SymComplex::new(c.re.to_float(), c.im.to_float()))
            }
        }
    }
}

impl<F: ToFloat> From<RealOrComplexRef<'_, F>> for Coefficient {
    fn from(value: RealOrComplexRef<'_, F>) -> Self {
        match value {
            RealOrComplexRef::Real(r) => {
                Coefficient::Float(SymComplex::new(r.to_float(), r.to_float().zero()))
            }
            RealOrComplexRef::Complex(c) => {
                Coefficient::Float(SymComplex::new(c.re.to_float(), c.im.to_float()))
            }
        }
    }
}

impl<T: ExportNumber + SingleFloat> ExportNumber for Complex<T> {
    fn export(&self) -> String {
        if self.im.is_zero() {
            self.re.export()
        } else {
            format!("{}, {}", self.re.export(), self.im.export())
        }
    }

    fn is_real(&self) -> bool {
        self.im.is_zero()
    }
}
impl<T: SingleFloat> SingleFloat for Complex<T>
where
    T: for<'a> RefMul<&'a T, Output = T>
        + for<'a> RefAdd<&'a T, Output = T>
        + for<'a> RefSub<&'a T, Output = T>
        + for<'a> RefDiv<&'a T, Output = T>,
{
    fn is_finite(&self) -> bool {
        self.re.is_finite() && self.im.is_finite()
    }

    fn is_one(&self) -> bool {
        self.re.is_one() && self.im.is_zero()
    }

    fn is_zero(&self) -> bool {
        self.re.is_zero() && self.im.is_zero()
    }

    fn from_rational(&self, rat: &Rational) -> Self {
        Complex {
            re: self.re.from_rational(rat),
            im: self.im.zero(),
        }
    }
}

impl<T: FloatLike> FloatLike for Complex<T>
where
    T: for<'a> RefMul<&'a T, Output = T>
        + for<'a> RefAdd<&'a T, Output = T>
        + for<'a> RefSub<&'a T, Output = T>
        + for<'a> RefDiv<&'a T, Output = T>,
{
    #[inline(always)]
    fn is_fully_zero(&self) -> bool {
        self.re.is_fully_zero() && self.im.is_fully_zero()
    }
    #[inline]
    fn mul_add(&self, a: &Self, b: &Self) -> Self {
        self.clone() * a + b.clone()
    }

    #[inline]
    fn neg(&self) -> Self {
        Complex {
            re: -self.re.clone(),
            im: -self.im.clone(),
        }
    }

    #[inline]
    fn zero(&self) -> Self {
        Complex {
            re: self.re.zero(),
            im: self.im.zero(),
        }
    }

    fn new_zero() -> Self {
        Complex {
            re: T::new_zero(),
            im: T::new_zero(),
        }
    }

    fn one(&self) -> Self {
        Complex {
            re: self.re.one(),
            im: self.im.zero(),
        }
    }

    fn pow(&self, e: u64) -> Self {
        // TODO: use binary exponentiation
        let mut r = self.one();
        for _ in 0..e {
            r *= self;
        }
        r
    }

    fn inv(&self) -> Self {
        self.inv()
    }

    fn from_usize(&self, a: usize) -> Self {
        Complex {
            re: self.re.from_usize(a),
            im: self.im.zero(),
        }
    }

    fn from_i64(&self, a: i64) -> Self {
        Complex {
            re: self.re.from_i64(a),
            im: self.im.zero(),
        }
    }

    #[inline(always)]
    fn get_precision(&self) -> u32 {
        self.re.get_precision().min(self.im.get_precision())
    }

    #[inline(always)]
    fn get_epsilon(&self) -> f64 {
        (2.0f64).powi(-(self.get_precision() as i32))
    }

    #[inline(always)]
    fn fixed_precision(&self) -> bool {
        self.re.fixed_precision() || self.im.fixed_precision()
    }

    fn sample_unit<R: Rng + ?Sized>(&self, rng: &mut R) -> Self {
        Complex {
            re: self.re.sample_unit(rng),
            im: self.im.zero(),
        }
    }
}

#[cfg(feature = "shadowing")]
impl<T: Real> Real for Complex<T>
where
    T: for<'a> RefMul<&'a T, Output = T>
        + for<'a> RefAdd<&'a T, Output = T>
        + for<'a> RefSub<&'a T, Output = T>
        + for<'a> RefDiv<&'a T, Output = T>
        + RefOne
        + RefZero,
{
    #[inline(always)]
    fn conj(&self) -> Self {
        Complex::new(self.re.clone(), -self.im.clone())
    }
    fn i(&self) -> Option<Self> {
        Some(Complex {
            re: self.re.zero(),
            im: self.re.one(),
        })
    }
    fn e(&self) -> Self {
        Complex {
            re: self.re.e(),
            im: self.im.zero(),
        }
    }

    fn pi(&self) -> Self {
        Complex {
            re: self.re.pi(),
            im: self.im.zero(),
        }
    }

    fn phi(&self) -> Self {
        Complex {
            re: self.re.phi(),
            im: self.im.zero(),
        }
    }

    fn euler(&self) -> Self {
        Complex {
            re: self.re.euler(),
            im: self.im.zero(),
        }
    }

    #[inline]
    fn norm(&self) -> Self {
        Complex::new(self.norm_squared().sqrt(), self.im.zero())
    }

    #[inline]
    fn sqrt(&self) -> Self {
        let (r, phi) = (self.norm_squared().sqrt(), self.im.atan2(&self.re));
        let phi = phi.ref_div(&self.re.from_usize(2));
        let r = r.sqrt();
        Complex::new(r.ref_mul(&phi.cos()), r.ref_mul(&phi.sin()))
    }

    #[inline]
    fn log(&self) -> Self {
        Complex::new(self.norm_squared().sqrt().log(), self.im.atan2(&self.re))
    }

    #[inline]
    fn exp(&self) -> Self {
        let r = self.re.exp();
        Complex::new(r.clone() * self.im.cos(), r * self.im.sin())
    }

    #[inline]
    fn sin(&self) -> Self {
        Complex::new(
            self.re.sin() * self.im.cosh(),
            self.re.cos() * self.im.sinh(),
        )
    }

    #[inline]
    fn cos(&self) -> Self {
        Complex::new(
            self.re.cos() * self.im.cosh(),
            -self.re.sin() * self.im.sinh(),
        )
    }

    #[inline]
    fn tan(&self) -> Self {
        let (r, i) = (self.re.clone() + &self.re, self.im.clone() + &self.im);
        let m = r.cos() + i.cosh();
        Self::new(r.sin() / &m, i.sinh() / m)
    }

    #[inline]
    fn asin(&self) -> Self {
        let i = self.ref_i();
        -i.clone() * ((self.one() - self.clone() * self).sqrt() + i * self).log()
    }

    #[inline]
    fn acos(&self) -> Self {
        let i = self.ref_i();
        -i.clone() * (i * (self.one() - self.clone() * self).sqrt() + self).log()
    }

    #[inline]
    fn atan2(&self, x: &Self) -> Self {
        // TODO: pick proper branch
        let r = self.clone() / x;
        let i = self.ref_i();
        let one = self.one();
        let two = one.clone() + &one;
        // TODO: add edge cases
        ((&one + &i * &r).log() - (&one - &i * r).log()) / (two * i)
    }

    #[inline]
    fn sinh(&self) -> Self {
        Complex::new(
            self.re.sinh() * self.im.cos(),
            self.re.cosh() * self.im.sin(),
        )
    }

    #[inline]
    fn cosh(&self) -> Self {
        Complex::new(
            self.re.cosh() * self.im.cos(),
            self.re.sinh() * self.im.sin(),
        )
    }

    #[inline]
    fn tanh(&self) -> Self {
        let (two_re, two_im) = (self.re.clone() + &self.re, self.im.clone() + &self.im);
        let m = two_re.cosh() + two_im.cos();
        Self::new(two_re.sinh() / &m, two_im.sin() / m)
    }

    #[inline]
    fn asinh(&self) -> Self {
        let one = self.one();
        (self.clone() + (one + self.clone() * self).sqrt()).log()
    }

    #[inline]
    fn acosh(&self) -> Self {
        let one = self.one();
        let two = one.clone() + &one;
        &two * (((self.clone() + &one) / &two).sqrt() + ((self.clone() - one) / &two).sqrt()).log()
    }

    #[inline]
    fn atanh(&self) -> Self {
        let one = self.one();
        let two = one.clone() + &one;
        // TODO: add edge cases
        ((&one + self).log() - (one - self).log()) / two
    }

    #[inline]
    fn powf(&self, e: &Self) -> Self {
        if e.re == self.re.zero() && e.im == self.im.zero() {
            self.one()
        } else if e.im == self.im.zero() {
            let (r, phi) = (self.norm_squared().sqrt(), self.im.atan2(&self.re));
            let r = r.powf(&e.re);
            let phi = phi * &e.re;
            Complex::new(r.ref_mul(&phi.cos()), r.ref_mul(&phi.sin()))
        } else {
            (e * self.log()).exp()
        }
    }
}

#[derive(Clone, Debug, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
pub struct CompiledComplexEvaluatorSpenso(CompiledComplexEvaluator);

impl CompiledComplexEvaluatorSpenso {
    pub fn evaluate(&mut self, args: &[Complex<f64>], out: &mut [Complex<f64>]) {
        unsafe {
            self.0.evaluate(
                transmute::<&[Complex<f64>], &[SymComplex<f64>]>(args),
                transmute::<&mut [Complex<f64>], &mut [SymComplex<f64>]>(out),
            );
        }
    }
}

impl EvaluatorLoader<Complex<f64>> for CompiledComplexEvaluatorSpenso {
    fn load_with_settings(
        file: impl AsRef<std::path::Path>,
        function_name: &str,
        settings: (),
    ) -> Result<Self, String> {
        Ok(CompiledComplexEvaluatorSpenso(
            CompiledComplexEvaluator::load_with_settings(file, function_name, settings)?,
        ))
    }
}

impl CompiledNumber for Complex<f64> {
    type Evaluator = CompiledComplexEvaluatorSpenso;
    type Settings = ();
    const SUFFIX: &'static str = "complexf64";
    fn get_default_compile_options() -> CompileOptions {
        CompileOptions::default()
    }

    fn export_cpp<N: ExportNumber + SingleFloat>(
        eval: &ExpressionEvaluator<N>,
        function_name: &str,
        settings: ExportSettings,
    ) -> Result<String, String> {
        SymComplex::<f64>::export_cpp(eval, function_name, settings)
    }
}
