use crate::{ParameterizationMapping, ParameterizationMode, Settings, MAX_LOOP};
use colored::Colorize;
use hyperdual::Hyperdual;
use itertools::{izip, Itertools};
use lorentz_vector::{Field, LorentzVector, RealNumberLike};
use num::Complex;
use num::ToPrimitive;
use num_traits::{Float, FloatConst, FromPrimitive, Num, NumAssign, NumCast, Signed};
use num_traits::{Inv, One, Zero};
use smartstring::{LazyCompact, SmartString};
use statrs::function::gamma::{gamma, gamma_lr, gamma_ur};
use std::cmp::{Ord, Ordering};
use std::ops::Neg;

#[allow(unused_imports)]
use log::{debug, info};
use symbolica::printer::{AtomPrinter, PrintOptions};
use symbolica::representations::Atom;

use git_version::git_version;
pub const GIT_VERSION: &str = git_version!();
pub const VERSION: &str = "0.0.1";

#[allow(unused)]
const MAX_DIMENSION: usize = MAX_LOOP * 3;

pub const PINCH_TEST_THRESHOLD: f64 = 1e-10;

pub trait FloatConvertFrom<U> {
    fn convert_from(x: &U) -> Self;
}

impl FloatConvertFrom<f64> for f64 {
    fn convert_from(x: &f64) -> f64 {
        *x
    }
}

impl FloatConvertFrom<f128::f128> for f64 {
    fn convert_from(x: &f128::f128) -> f64 {
        (*x).to_f64().unwrap()
    }
}

impl FloatConvertFrom<f128::f128> for f128::f128 {
    fn convert_from(x: &f128::f128) -> f128::f128 {
        *x
    }
}

impl FloatConvertFrom<f64> for f128::f128 {
    fn convert_from(x: &f64) -> f128::f128 {
        f128::f128::from_f64(*x).unwrap()
    }
}

pub trait FloatLike:
    From<f64>
    + FloatConvertFrom<f64>
    + FloatConvertFrom<f128::f128>
    + Num
    + FromPrimitive
    + Float
    + Field
    + RealNumberLike
    + Signed
    + FloatConst
    + std::fmt::LowerExp
    + 'static
    + Signum
{
}

impl FloatLike for f64 {}
impl FloatLike for f128::f128 {}

/// An iterator which iterates two other iterators simultaneously
#[derive(Clone, Debug)]
#[must_use = "iterator adaptors are lazy and do nothing unless consumed"]
pub struct ZipEq<I, J> {
    a: I,
    b: J,
}

/// An iterator which iterates two other iterators simultaneously and checks
/// if the sizes are equal in debug mode.
#[allow(unused)]
pub fn zip_eq<I, J>(i: I, j: J) -> ZipEq<I::IntoIter, J::IntoIter>
where
    I: IntoIterator,
    J: IntoIterator,
{
    ZipEq {
        a: i.into_iter(),
        b: j.into_iter(),
    }
}

impl<I, J> Iterator for ZipEq<I, J>
where
    I: Iterator,
    J: Iterator,
{
    type Item = (I::Item, J::Item);

    fn next(&mut self) -> Option<Self::Item> {
        match (self.a.next(), self.b.next()) {
            (None, None) => None,
            (Some(a), Some(b)) => Some((a, b)),
            (None, Some(_)) => {
                #[cfg(debug_assertions)]
                panic!("Unequal length of iterators; first iterator finished first");
                #[cfg(not(debug_assertions))]
                None
            }
            (Some(_), None) => {
                #[cfg(debug_assertions)]
                panic!("Unequal length of iterators; second iterator finished first");
                #[cfg(not(debug_assertions))]
                None
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let sa = self.a.size_hint();
        let sb = self.b.size_hint();
        (sa.0.min(sb.0), sa.1.zip(sb.1).map(|(ua, ub)| ua.min(ub)))
    }
}

impl<I, J> ExactSizeIterator for ZipEq<I, J>
where
    I: ExactSizeIterator,
    J: ExactSizeIterator,
{
}

pub fn parse_python_expression(
    expression: &str,
    sb_state: &mut symbolica::state::State,
    sb_workspace: &symbolica::state::Workspace,
) -> Atom {
    let processed_string = String::from(expression)
        .replace("**", "^")
        .replace("cmath.sqrt", "sqrt")
        .replace("cmath.pi", "pi")
        .replace("math.sqrt", "sqrt")
        .replace("math.pi", "pi");
    Atom::parse(processed_string.as_str(), sb_state, sb_workspace)
        .map_err(|e| {
            format!(
                "Failed to parse expression : '{}'\nError: {}",
                processed_string,
                e.to_string()
            )
        })
        .unwrap()
}

pub fn to_str_expression(
    expression: &Atom,
    sb_state: &symbolica::state::State,
) -> SmartString<LazyCompact> {
    format!(
        "{}",
        AtomPrinter::new_with_options(
            expression.as_view(),
            PrintOptions {
                terms_on_new_line: false,
                color_top_level_sum: false,
                color_builtin_functions: false,
                print_finite_field: false,
                explicit_rational_polynomial: false,
                number_thousands_separator: None,
                multiplication_operator: '*',
                square_brackets_for_function: false,
                num_exp_as_superscript: false,
                latex: false
            },
            sb_state
        )
    )
    .into()
}

/// Format a mean ± sdev as mean(sdev) with the correct number of digits.
/// Based on the Python package gvar.
pub fn format_uncertainty(mean: f64, sdev: f64) -> String {
    fn ndec(x: f64, offset: usize) -> i32 {
        let mut ans = (offset as f64 - x.log10()) as i32;
        if ans > 0 && x * 10.0.powi(ans) >= [0.5, 9.5, 99.5][offset] {
            ans -= 1;
        }
        if ans < 0 {
            0
        } else {
            ans
        }
    }
    let v = mean;
    let dv = sdev.abs();

    // special cases
    if v.is_nan() || dv.is_nan() {
        format!("{:e} ± {:e}", v, dv)
    } else if dv.is_infinite() {
        format!("{:e} ± inf", v)
    } else if v == 0. && (dv >= 1e5 || dv < 1e-4) {
        if dv == 0. {
            "0(0)".to_owned()
        } else {
            let e = format!("{:.1e}", dv);
            let mut ans = e.split('e');
            let e1 = ans.next().unwrap();
            let e2 = ans.next().unwrap();
            "0.0(".to_owned() + e1 + ")e" + e2
        }
    } else if v == 0. {
        if dv >= 9.95 {
            format!("0({:.0})", dv)
        } else if dv >= 0.995 {
            format!("0.0({:.1})", dv)
        } else {
            let ndecimal = ndec(dv, 2);
            format!(
                "{:.*}({:.0})",
                ndecimal as usize,
                v,
                dv * 10.0.powi(ndecimal)
            )
        }
    } else if dv == 0. {
        let e = format!("{:e}", v);
        let mut ans = e.split('e');
        let e1 = ans.next().unwrap();
        let e2 = ans.next().unwrap();
        if e2 != "0" {
            e1.to_owned() + "(0)e" + e2
        } else {
            e1.to_owned() + "(0)"
        }
    } else if dv > 1e4 * v.abs() {
        format!("{:.1e} ± {:.2e}", v, dv)
    } else if v.abs() >= 1e6 || v.abs() < 1e-5 {
        // exponential notation for large |self.mean|
        let exponent = v.abs().log10().floor();
        let fac = 10.0.powf(exponent);
        let mantissa = format_uncertainty(v / fac, dv / fac);
        let e = format!("{:.0e}", fac);
        let mut ee = e.split('e');
        mantissa + "e" + ee.nth(1).unwrap()
    }
    // normal cases
    else if dv >= 9.95 {
        if v.abs() >= 9.5 {
            format!("{:.0}({:.0})", v, dv)
        } else {
            let ndecimal = ndec(v.abs(), 1);
            format!("{:.*}({:.*})", ndecimal as usize, v, ndecimal as usize, dv)
        }
    } else if dv >= 0.995 {
        if v.abs() >= 0.95 {
            format!("{:.1}({:.1})", v, dv)
        } else {
            let ndecimal = ndec(v.abs(), 1);
            format!("{:.*}({:.*})", ndecimal as usize, v, ndecimal as usize, dv)
        }
    } else {
        let ndecimal = ndec(v.abs(), 1).max(ndec(dv, 2));
        format!(
            "{:.*}({:.0})",
            ndecimal as usize,
            v,
            dv * 10.0.powi(ndecimal)
        )
    }
}

/// Compare two slices, selecting on length first
#[allow(unused)]
pub fn compare_slice<T: Ord>(slice1: &[T], slice2: &[T]) -> Ordering {
    match slice1.len().cmp(&slice2.len()) {
        Ordering::Equal => (),
        non_eq => return non_eq,
    }

    let l = slice1.len();
    // Slice to the loop iteration range to enable bound check
    // elimination in the compiler
    let lhs = &slice1[..l];
    let rhs = &slice2[..l];

    for i in 0..l {
        match lhs[i].cmp(&rhs[i]) {
            Ordering::Equal => (),
            non_eq => return non_eq,
        }
    }

    Ordering::Equal
}

pub trait Signum {
    fn multiply_sign(&self, sign: i8) -> Self;
}

impl Signum for f128::f128 {
    #[inline]
    fn multiply_sign(&self, sign: i8) -> f128::f128 {
        match sign {
            1 => self.clone(),
            0 => f128::f128::zero(),
            -1 => self.neg(),
            _ => unreachable!("Sign should be -1,0,1"),
        }
    }
}

impl Signum for f64 {
    #[inline]
    fn multiply_sign(&self, sign: i8) -> f64 {
        match sign {
            1 => self.clone(),
            0 => f64::zero(),
            -1 => self.neg(),
            _ => unreachable!("Sign should be -1,0,1"),
        }
    }
}

impl<T: Num + Neg<Output = T> + Copy> Signum for Complex<T> {
    #[inline]
    fn multiply_sign(&self, sign: i8) -> Complex<T> {
        match sign {
            1 => *self,
            0 => Complex::zero(),
            -1 => Complex::new(-self.re, -self.im),
            _ => unreachable!("Sign should be -1,0,1"),
        }
    }
}

impl<T: FloatLike, const U: usize> Signum for Hyperdual<T, U> {
    #[inline]
    fn multiply_sign(&self, sign: i8) -> Hyperdual<T, U> {
        match sign {
            1 => *self,
            0 => Hyperdual::zero(),
            -1 => -*self,
            _ => unreachable!("Sign should be -1,0,1"),
        }
    }
}

impl<T: Field> Signum for LorentzVector<T> {
    #[inline]
    fn multiply_sign(&self, sign: i8) -> LorentzVector<T> {
        match sign {
            1 => *self,
            0 => LorentzVector::default(),
            -1 => -self,
            _ => unreachable!("Sign should be -1,0,1"),
        }
    }
}

#[allow(unused)]
#[inline]
/// Invert with better precision
pub fn finv<T: Float>(c: Complex<T>) -> Complex<T> {
    let norm = c.norm();
    c.conj() / norm / norm
}

#[allow(unused)]
#[inline]
pub fn powi<T: Float + NumAssign>(c: Complex<T>, n: usize) -> Complex<T> {
    let mut c1 = Complex::<T>::one();
    for _ in 0..n {
        c1 *= c;
    }
    c1
}

#[allow(unused)]
#[inline]
pub fn powf<T: 'static + Float + NumAssign + std::fmt::Debug, const N: usize>(
    h: Hyperdual<T, N>,
    n: T,
) -> Hyperdual<T, N> {
    let r = Float::powf(h.real(), n - T::one());
    let rr = n * r;
    h.map_dual(r * h.real(), |x| *x * rr)
}

#[allow(unused)]
pub fn evaluate_signature<T: Field>(
    signature: &[i8],
    momenta: &[LorentzVector<T>],
) -> LorentzVector<T> {
    let mut momentum = LorentzVector::default();
    for (&sign, mom) in zip_eq(signature, momenta) {
        match sign {
            0 => {}
            1 => momentum += mom,
            -1 => momentum -= mom,
            _ => {
                #[cfg(debug_assertions)]
                panic!("Sign should be -1,0,1")
            }
        }
    }

    momentum
}

#[allow(unused)]
#[inline]
pub fn pinch_dampening_function<T: FloatLike>(
    dampening_arg: T,
    delta_t: T,
    powers: (i32, i32),
    multiplier: f64,
) -> T {
    // Make sure the function is even in t-tstar
    assert!(powers.1 % 2 == 0);
    let a = dampening_arg.powi(powers.0);
    return a / (a + Into::<T>::into(multiplier as f64) * delta_t.powi(powers.1));
}

pub fn h<T: FloatLike>(
    t: T,
    tstar: Option<T>,
    sigma: Option<T>,
    h_function_settings: &crate::HFunctionSettings,
) -> T {
    let sig = if let Some(s) = sigma {
        s
    } else {
        Into::<T>::into(h_function_settings.sigma as f64)
    };
    let power = h_function_settings.power;
    match h_function_settings.function {
        crate::HFunction::Exponential => {
            (-(t * t) / (sig * sig)).exp() * Into::<T>::into(2 as f64)
                / (<T as FloatConst>::PI().sqrt() * sig)
        }
        crate::HFunction::PolyExponential => {
            // Result of \int_0^{\infty} dt (t/sigma)^{-p} exp(2-t^2/sigma^2-sigma^2/t^2)
            let normalisation = match power {
                None | Some(0) => <T as FloatConst>::PI().sqrt() * sig / Into::<T>::into(2 as f64),
                Some(1) => Into::<T>::into(0.84156821507077141791912486734584) * sig,
                Some(3) => Into::<T>::into(1.0334768470686885731753571058796) * sig,
                Some(4) => Into::<T>::into(1.3293403881791370204736256125059) * sig,
                Some(6) => Into::<T>::into(2.8802375077214635443595221604294) * sig,
                Some(7) => Into::<T>::into(4.7835669713476085553643210523305) * sig,
                Some(9) => Into::<T>::into(16.225745976182285657187445130217) * sig,
                Some(10) => Into::<T>::into(32.735007058911249129163030707957) * sig,
                Some(12) => Into::<T>::into(155.83746592258341696260606919938) * sig,
                Some(13) => Into::<T>::into(364.65850035656604157775795299621) * sig,
                Some(15) => Into::<T>::into(2257.6375530154730006506618195504) * sig,
                Some(16) => Into::<T>::into(5939.8044185378636927153327426791) * sig,
                _ => panic!(
                    "Value {} of power in poly exponential h function not supported",
                    power.unwrap()
                ),
            };
            let prefactor = match power {
                None | Some(0) => normalisation.inv(),
                Some(p) => (t / sig).powi(-(p as i32)) / normalisation,
            };
            prefactor
                * (Into::<T>::into(2 as f64) - (t * t) / (sig * sig) - (sig * sig) / (t * t)).exp()
        }
        crate::HFunction::PolyLeftRightExponential => {
            // Result of \int_0^{\infty} dt (t/sigma)^{-p} exp( -((t^2/sigma^2 +1)/ (t/sigma) -2) )
            let normalisation = match power {
                None | Some(0) => Into::<T>::into(2.0669536941373771463507142117592) * sig,
                Some(1) => Into::<T>::into(1.6831364301415428358382497346917) * sig,
                Some(3) => Into::<T>::into(3.7500901242789199821889639464509) * sig,
                Some(4) => Into::<T>::into(9.567133942695217110728642104661) * sig,
                Some(6) => Into::<T>::into(139.3731017521535023682282031464) * sig,
                Some(7) => Into::<T>::into(729.31700071313208315551590599241) * sig,
                Some(9) => Into::<T>::into(32336.242742929754092264781379699) * sig,
                Some(10) => Into::<T>::into(263205.21704946897873941957467669) * sig,
                Some(12) => Into::<T>::into(2.4275037178930974606209829109376e7) * sig,
                Some(13) => Into::<T>::into(2.694265921644288712310551611566e8) * sig,
                Some(15) => Into::<T>::into(9.0407420577601240420566809613266e12) * sig,
                Some(16) => Into::<T>::into(1.4525174802464911684668046368611e14) * sig,
                _ => panic!(
                    "Value {} of power in poly exponential h function not supported",
                    power.unwrap()
                ),
            };
            let prefactor = match power {
                None | Some(0) => normalisation.inv(),
                Some(p) => (t / sig).powi(-(p as i32)) / normalisation,
            };
            prefactor
                * (Into::<T>::into(2 as f64) - ((t * t) / (sig * sig) + T::one()) / (t / sig)).exp()
        }
        crate::HFunction::ExponentialCT => {
            let delta_t_sq = (t - tstar.unwrap()) * (t - tstar.unwrap());
            let tstar_sq = tstar.unwrap() * tstar.unwrap();
            // info!("dampener: {}", dampener);
            // info!("delta_t_sq: {}", delta_t_sq);
            // info!("tstar_sq: {}", tstar_sq);
            // info!(
            //     "Exp arg: {}",
            //     -sig.inv() * (delta_t_sq / tstar_sq + sig * sig * (dampener * dampener))
            // );
            // info!(
            //     "result: {}",
            //     (-sig.inv() * (delta_t_sq / tstar_sq + sig * sig * (dampener * dampener))).exp()
            // );
            if h_function_settings.enabled_dampening {
                let dampener = delta_t_sq / (delta_t_sq - tstar_sq);
                (-sig.inv() * (delta_t_sq / tstar_sq + sig * sig * (dampener * dampener))).exp()
            } else {
                (-sig.inv() * (delta_t_sq / tstar_sq)).exp()
            }
        }
    }
}

/// Calculate the determinant of any complex-valued input matrix using LU-decomposition.
/// Original C-code by W. Gong and D.E. Soper.
#[allow(unused)]
pub fn determinant<T: Float + RealNumberLike>(
    bb: &Vec<Complex<T>>,
    dimension: usize,
) -> Complex<T> {
    // Define matrix related variables.
    let mut determinant = Complex::new(T::one(), T::zero());
    let mut indx = [0; MAX_DIMENSION];
    let mut d = 1; // initialize parity parameter

    // Inintialize the matrix to be decomposed with the transferred matrix b.
    let mut aa = bb.clone();

    // Define parameters used in decomposition.
    let mut imax = 0;
    let mut flag = 1;
    let mut dumc;
    let mut sum;

    let mut aamax;
    let mut dumr;
    let mut vv = [T::zero(); MAX_DIMENSION];

    // Get the implicit scaling information.
    for i in 0..dimension {
        aamax = T::zero();
        for j in 0..dimension {
            let r = aa[i * dimension + j].norm_sqr();
            if r > aamax {
                aamax = r;
            }
        }
        // Set a flag to check if the determinant is zero.
        if aamax.is_zero() {
            flag = 0;
        }
        // Save the scaling.
        vv[i] = aamax.inv();
    }
    if flag == 1 {
        for j in 0..dimension {
            for i in 0..j {
                sum = aa[i * dimension + j];
                for k in 0..i {
                    sum = sum - aa[i * dimension + k] * aa[k * dimension + j];
                }
                aa[i * dimension + j] = sum;
            }
            //Initialize for the search for largest pivot element.
            aamax = T::zero();
            for i in j..dimension {
                sum = aa[i * dimension + j];
                for k in 0..j {
                    sum = sum - aa[i * dimension + k] * aa[k * dimension + j];
                }
                aa[i * dimension + j] = sum;
                // Figure of merit for the pivot.
                dumr = vv[i] * sum.norm_sqr();
                // Is it better than the best so far?
                if dumr >= aamax {
                    imax = i;
                    aamax = dumr;
                }
            }
            // See if we need to interchange rows.
            if j != imax {
                for k in 0..dimension {
                    dumc = aa[imax * dimension + k];
                    aa[imax * dimension + k] = aa[j * dimension + k];
                    aa[j * dimension + k] = dumc;
                }
                // Change the parity of d.
                d = -d;
                // Interchange the scale factor.
                vv[imax] = vv[j];
            }
            indx[j] = imax;
            if j + 1 != dimension {
                dumc = aa[j * dimension + j].inv();
                for i in j + 1..dimension {
                    aa[i * dimension + j] = aa[i * dimension + j] * dumc;
                }
            }
        }
    }
    // Calculate the determinant using the decomposed matrix.
    if flag == 0 {
        determinant = Complex::default();
    } else {
        // Multiply the diagonal elements.
        for diagonal in 0..dimension {
            determinant = determinant * aa[diagonal * dimension + diagonal];
        }
        determinant = determinant * <T as NumCast>::from(d).unwrap();
    }
    determinant
}

#[allow(unused)]
pub fn next_combination_with_replacement(state: &mut [usize], max_entry: usize) -> bool {
    for i in (0..state.len()).rev() {
        if state[i] < max_entry {
            state[i] += 1;
            for j in i + 1..state.len() {
                state[j] = state[i]
            }
            return true;
        }
    }
    false
}

#[allow(unused)]
pub fn compute_momentum<T: FloatLike>(
    signature: &(Vec<isize>, Vec<isize>),
    loop_moms: &Vec<LorentzVector<T>>,
    external_moms: &Vec<LorentzVector<T>>,
) -> LorentzVector<T> {
    let mut res = LorentzVector::default();
    for (i_l, sign) in signature.0.iter().enumerate() {
        res += loop_moms[i_l] * Into::<T>::into(*sign as f64);
    }
    for (i_l, sign) in signature.1.iter().enumerate() {
        res += external_moms[i_l] * Into::<T>::into(*sign as f64);
    }
    res
}

// Bilinear form for E-surface defined as sqrt[(k+p1)^2+m1sq] + sqrt[(k+p2)^2+m2sq] + e_shift
// The Bilinear system then reads 4 k.a.k + 4 k.n + C = 0
#[allow(unused)]
pub fn one_loop_e_surface_bilinear_form<T: FloatLike>(
    p1: &[T; 3],
    p2: &[T; 3],
    m1_sq: T,
    m2_sq: T,
    e_shift: T,
) -> ([[T; 3]; 3], [T; 3], T) {
    let two = Into::<T>::into(2 as f64);
    let e_shift_sq = e_shift * e_shift;
    let p1_sq = p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2];
    let p2_sq = p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2];

    let mut a = [[T::zero(); 3]; 3];
    a[0][0] = (p1[0] - p2[0] - e_shift) * (p2[0] - p1[0] - e_shift);
    a[0][1] = (p1[0] - p2[0]) * (p2[1] - p1[1]);
    a[1][0] = a[0][1];
    a[0][2] = (p1[0] - p2[0]) * (p2[2] - p1[2]);
    a[2][0] = a[0][2];
    a[1][1] = (p1[1] - p2[1] - e_shift) * (p2[1] - p1[1] - e_shift);
    a[1][2] = (p1[1] - p2[1]) * (p2[2] - p1[2]);
    a[2][1] = a[1][2];
    a[2][2] = (p1[2] - p2[2] - e_shift) * (p2[2] - p1[2] - e_shift);

    let mut b = [T::zero(); 3];
    b[0] = (p2[0] - p1[0]) * (m1_sq - m2_sq + p1_sq - p2_sq) + e_shift_sq * (p1[0] + p2[0]);
    b[1] = (p2[1] - p1[1]) * (m1_sq - m2_sq + p1_sq - p2_sq) + e_shift_sq * (p1[1] + p2[1]);
    b[2] = (p2[2] - p1[2]) * (m1_sq - m2_sq + p1_sq - p2_sq) + e_shift_sq * (p1[2] + p2[2]);

    let c = -e_shift_sq * e_shift_sq + two * e_shift_sq * (m1_sq + m2_sq + p1_sq + p2_sq)
        - (m1_sq - m2_sq + p1_sq - p2_sq) * (m1_sq - m2_sq + p1_sq - p2_sq);

    (a, b, c)
}

#[allow(unused)]
pub fn one_loop_e_surface_exists<T: FloatLike>(
    p1: &[T; 3],
    p2: &[T; 3],
    m1_sq: T,
    m2_sq: T,
    e_shift: T,
) -> (bool, bool) {
    let p_norm_sq = (p1[0] - p2[0]) * (p1[0] - p2[0])
        + (p1[1] - p2[1]) * (p1[1] - p2[1])
        + (p1[2] - p2[2]) * (p1[2] - p2[2]);

    // /!\ In alphaLoop this should be done without numerical check but purely symbolically, or at least
    // one must make sure there is no weird transition from non-pinched to pinched for the 2->2 massless E-surface sandwich,
    // i.e. such cases should be *both* existing and pinched!
    if e_shift > Into::<T>::into(PINCH_TEST_THRESHOLD) {
        return (false, false);
    }
    let test = (e_shift * e_shift - p_norm_sq)
        - (m1_sq.sqrt() + m2_sq.sqrt()) * (m1_sq.sqrt() + m2_sq.sqrt());
    if test.abs() < Into::<T>::into(PINCH_TEST_THRESHOLD) {
        return (false, true);
    } else {
        if test < T::zero() {
            return (false, false);
        } else {
            return (true, false);
        }
    }
}

#[allow(unused)]
pub fn one_loop_eval_e_surf<T: FloatLike>(
    k: &[T; 3],
    p1: &[T; 3],
    p2: &[T; 3],
    m1_sq: T,
    m2_sq: T,
    e_shift: T,
) -> T {
    ((k[0] + p1[0]) * (k[0] + p1[0])
        + (k[1] + p1[1]) * (k[1] + p1[1])
        + (k[2] + p1[2]) * (k[2] + p1[2])
        + m1_sq)
        .sqrt()
        + ((k[0] + p2[0]) * (k[0] + p2[0])
            + (k[1] + p2[1]) * (k[1] + p2[1])
            + (k[2] + p2[2]) * (k[2] + p2[2])
            + m2_sq)
            .sqrt()
        + e_shift
}

#[allow(unused)]
pub fn one_loop_eval_e_surf_k_derivative<T: FloatLike>(
    k: &[T; 3],
    p1: &[T; 3],
    p2: &[T; 3],
    m1_sq: T,
    m2_sq: T,
) -> [T; 3] {
    let e1 = ((k[0] + p1[0]) * (k[0] + p1[0])
        + (k[1] + p1[1]) * (k[1] + p1[1])
        + (k[2] + p1[2]) * (k[2] + p1[2])
        + m1_sq)
        .sqrt();
    let e2 = ((k[0] + p2[0]) * (k[0] + p2[0])
        + (k[1] + p2[1]) * (k[1] + p2[1])
        + (k[2] + p2[2]) * (k[2] + p2[2])
        + m2_sq)
        .sqrt();
    [
        (k[0] + p1[0]) / e1 + (k[0] + p2[0]) / e2,
        (k[1] + p1[1]) / e1 + (k[1] + p2[1]) / e2,
        (k[2] + p1[2]) / e1 + (k[2] + p2[2]) / e2,
    ]
}

#[allow(unused)]
pub fn one_loop_get_e_surf_t_scaling<T: FloatLike>(
    k: &[T; 3],
    p1: &[T; 3],
    p2: &[T; 3],
    m1_sq: T,
    m2_sq: T,
    e_shift: T,
) -> [T; 2] {
    let (a, b, c_coef) = one_loop_e_surface_bilinear_form(p1, p2, m1_sq, m2_sq, e_shift);
    let mut a_coef = T::zero();
    for i in 0..=2 {
        for j in 0..=2 {
            a_coef += k[i] * a[i][j] * k[j];
        }
    }
    a_coef *= Into::<T>::into(4 as f64);
    let mut b_coef = T::zero();
    for i in 0..=2 {
        b_coef += k[i] * b[i];
    }
    b_coef *= Into::<T>::into(4 as f64);
    let discr = b_coef * b_coef - Into::<T>::into(4 as f64) * a_coef * c_coef;
    if discr < T::zero() {
        [T::zero(), T::zero()]
    } else {
        [
            (-b_coef + discr.sqrt()) / (Into::<T>::into(2 as f64) * a_coef),
            (-b_coef - discr.sqrt()) / (Into::<T>::into(2 as f64) * a_coef),
        ]
    }
}

pub fn box_muller<T: FloatLike>(x1: T, x2: T) -> (T, T) {
    let r = (-Into::<T>::into(2.) * x1.ln()).sqrt();
    let theta = Into::<T>::into(2.) * T::PI() * x2;
    (r * theta.cos(), r * theta.sin())
}

pub fn compute_surface_and_volume<T: FloatLike>(n_dim: usize, radius: T) -> (T, T) {
    let mut surface = Into::<T>::into(2.0);
    let mut volume = T::one();
    for i in 1..n_dim + 1 {
        (surface, volume) = (
            Into::<T>::into(2.0) * <T as FloatConst>::PI() * volume,
            surface / Into::<T>::into(i as f64),
        );
    }
    (
        surface * radius.powi(n_dim as i32),
        volume * radius.powi(n_dim as i32),
    )
}

pub fn get_n_dim_for_n_loop_momenta(
    settings: &Settings,
    n_loop_momenta: usize,
    force_radius: bool,
) -> usize {
    match settings.parameterization.mode {
        ParameterizationMode::HyperSphericalFlat => {
            // Because we use Box-Muller, we need to have an even number of angular dimensions
            let mut n_dim = 3 * n_loop_momenta;
            if n_dim % 2 == 1 {
                n_dim += 1;
            }
            // Then if the radius is not forced, then we need to add one more dimension
            if !force_radius {
                n_dim += 1;
            }
            return n_dim;
        }
        ParameterizationMode::HyperSpherical
        | ParameterizationMode::Cartesian
        | ParameterizationMode::Spherical => {
            if force_radius {
                return 3 * n_loop_momenta - 1;
            } else {
                return 3 * n_loop_momenta;
            }
        }
    }
}

pub fn global_parameterize<T: FloatLike>(
    x: &[T],
    e_cm_squared: T,
    settings: &Settings,
    force_radius: bool,
) -> (Vec<[T; 3]>, T) {
    match settings.parameterization.mode {
        ParameterizationMode::HyperSpherical | ParameterizationMode::HyperSphericalFlat => {
            let e_cm = e_cm_squared.sqrt() * Into::<T>::into(settings.parameterization.shifts[0].0);
            let mut jac = T::one();
            // rescale the input to the desired range
            let mut x_r = Vec::with_capacity(x.len());
            if !force_radius {
                x_r.push(x[0]);
            } else {
                let lo = Into::<T>::into(settings.parameterization.input_rescaling[0][0].0);
                let hi = Into::<T>::into(settings.parameterization.input_rescaling[0][0].1);
                x_r.push(lo + x[0] * (hi - lo));
                jac *= Into::<T>::into(hi - lo);
            }
            let lo = Into::<T>::into(settings.parameterization.input_rescaling[0][1].0);
            let hi = Into::<T>::into(settings.parameterization.input_rescaling[0][1].1);
            x_r.push(lo + x[1] * (hi - lo));
            jac *= Into::<T>::into(hi - lo);
            for xi in &x[2..] {
                let lo = Into::<T>::into(settings.parameterization.input_rescaling[0][2].0);
                let hi = Into::<T>::into(settings.parameterization.input_rescaling[0][2].1);
                x_r.push(lo + *xi * (hi - lo));
                jac *= Into::<T>::into(hi - lo);
            }

            let radius: T = if force_radius {
                x[0]
            } else {
                match settings.parameterization.mapping {
                    ParameterizationMapping::Log => {
                        // r = e_cm * ln(1 + b*x/(1-x))
                        let b = Into::<T>::into(settings.parameterization.b);
                        let radius = e_cm * (T::one() + b * x_r[0] / (T::one() - x_r[0])).ln();
                        jac *=
                            e_cm * b / (T::one() - x_r[0]) / (T::one() + x_r[0] * (b - T::one()));
                        radius
                    }
                    ParameterizationMapping::Linear => {
                        // r = e_cm * b * x/(1-x)
                        let b = Into::<T>::into(settings.parameterization.b);
                        let radius = e_cm * b * x_r[0] / (T::one() - x_r[0]);
                        jac *= <T as num_traits::Float>::powi(e_cm * b + radius, 2) / e_cm / b;
                        radius
                    }
                }
            };
            match settings.parameterization.mode {
                ParameterizationMode::HyperSpherical => {
                    let phi = Into::<T>::into(2.) * <T as FloatConst>::PI() * x_r[1];
                    jac *= Into::<T>::into(2.) * <T as FloatConst>::PI();

                    let mut cos_thetas = Vec::with_capacity(x.len() - 2);
                    let mut sin_thetas = Vec::with_capacity(x.len() - 2);

                    for (i, xi) in x_r[2..].iter().enumerate() {
                        let cos_theta = -T::one() + Into::<T>::into(2.) * *xi;
                        jac *= Into::<T>::into(2.);
                        let sin_theta = (T::one() - cos_theta * cos_theta).sqrt();
                        if i > 0 {
                            jac *= sin_theta.powi(i as i32);
                        }
                        cos_thetas.push(cos_theta);
                        sin_thetas.push(sin_theta);
                    }

                    let mut concatenated_vecs = Vec::with_capacity(x.len() / 3);
                    let mut base = radius;
                    for (_i, (cos_theta, sin_theta)) in
                        cos_thetas.iter().zip(sin_thetas.iter()).enumerate()
                    {
                        concatenated_vecs.push(base * cos_theta);
                        base *= *sin_theta;
                    }
                    concatenated_vecs.push(base * phi.cos());
                    concatenated_vecs.push(base * phi.sin());

                    jac *= radius.powi((x.len() - 1) as i32); // hyperspherical coords

                    (
                        concatenated_vecs
                            .chunks(3)
                            .map(|v| [v[0], v[1], v[2]])
                            .collect(),
                        jac,
                    )
                }
                ParameterizationMode::HyperSphericalFlat => {
                    // As we will use Box Muller we expect an even number of random variables
                    assert!(x_r[1..].len() % 2 == 0);
                    let mut normal_distributed_xs = vec![];
                    for x_pair in x_r[1..].chunks(2) {
                        let (z1, z2) = box_muller(x_pair[0], x_pair[1]);
                        normal_distributed_xs.push(z1);
                        normal_distributed_xs.push(z2);
                    }
                    // ignore the last variable generated if we had to pad to get the 3*n_loop_momenta
                    if normal_distributed_xs.len() % 3 != 0 {
                        normal_distributed_xs.pop();
                    }
                    let curr_norm = normal_distributed_xs[..]
                        .iter()
                        .map(|&x| x * x)
                        .sum::<T>()
                        .sqrt();
                    let surface =
                        compute_surface_and_volume(normal_distributed_xs.len() - 1, radius).0;
                    jac *= surface;
                    let rescaling_factor = radius / curr_norm;
                    (
                        normal_distributed_xs
                            .chunks(3)
                            .map(|v| {
                                [
                                    v[0] * rescaling_factor,
                                    v[1] * rescaling_factor,
                                    v[2] * rescaling_factor,
                                ]
                            })
                            .collect(),
                        jac,
                    )
                }
                _ => unreachable!(),
            }
        }
        ParameterizationMode::Cartesian | ParameterizationMode::Spherical => {
            if force_radius {
                panic!("Cannot force radius for non-hyperspherical parameterization.");
            }
            let mut jac = T::one();
            let mut vecs = Vec::with_capacity(x.len() / 3);
            for (i, xi) in x.chunks(3).enumerate() {
                let (vec_i, jac_i) = parameterize3d(xi, e_cm_squared, i, settings);
                vecs.push(vec_i);
                jac *= jac_i;
            }
            (vecs, jac)
        }
    }
}

#[allow(unused)]
pub fn global_inv_parameterize<T: FloatLike>(
    moms: &Vec<LorentzVector<T>>,
    e_cm_squared: T,
    settings: &Settings,
    force_radius: bool,
) -> (Vec<T>, T) {
    match settings.parameterization.mode {
        ParameterizationMode::HyperSpherical => {
            let e_cm = e_cm_squared.sqrt() * Into::<T>::into(settings.parameterization.shifts[0].0);
            let mut inv_jac = T::one();
            let mut xs = Vec::with_capacity(moms.len() * 3);

            let cartesian_xs = moms
                .iter()
                .map(|lv| [lv.x, lv.y, lv.z])
                .flatten()
                .collect::<Vec<_>>();

            let mut k_r_sq = cartesian_xs.iter().map(|xi| *xi * xi).sum::<T>();
            // cover the degenerate case
            if k_r_sq.is_zero() {
                return (vec![T::zero(); cartesian_xs.len()], T::zero());
            }
            let k_r = k_r_sq.sqrt();
            if force_radius {
                xs.push(k_r);
            } else {
                match settings.parameterization.mapping {
                    ParameterizationMapping::Log => {
                        let b = Into::<T>::into(settings.parameterization.b);
                        let x1 = T::one() - b / (-T::one() + b + (k_r / e_cm).exp());
                        inv_jac /= e_cm * b / (T::one() - x1) / (T::one() + x1 * (b - T::one()));
                        xs.push(x1);
                    }
                    ParameterizationMapping::Linear => {
                        let b = Into::<T>::into(settings.parameterization.b);
                        inv_jac /= <T as num_traits::Float>::powi(e_cm * b + k_r, 2) / e_cm / b;
                        xs.push(k_r / (e_cm * b + k_r));
                    }
                }
            };

            let y = cartesian_xs[cartesian_xs.len() - 2];
            let x = cartesian_xs[cartesian_xs.len() - 1];
            let xphi = if x < T::zero() {
                T::one() + Into::<T>::into(0.5) * T::FRAC_1_PI() * T::atan2(x, y)
            } else {
                Into::<T>::into(0.5) * T::FRAC_1_PI() * T::atan2(x, y)
            };
            xs.push(xphi);
            inv_jac /= Into::<T>::into(2.) * <T as FloatConst>::PI();

            for (i, &x) in cartesian_xs[..cartesian_xs.len() - 2].iter().enumerate() {
                xs.push(Into::<T>::into(0.5) * (T::one() + x / k_r_sq.sqrt()));
                inv_jac /= Into::<T>::into(2.);
                if i > 0 {
                    inv_jac /= (T::one() - (x * x / k_r_sq)).sqrt().powi(i as i32);
                }
                k_r_sq -= x * x;
            }

            inv_jac /= k_r.powi((cartesian_xs.len() - 1) as i32);

            let lo = Into::<T>::into(settings.parameterization.input_rescaling[0][0].0);
            let hi = Into::<T>::into(settings.parameterization.input_rescaling[0][0].1);
            xs[0] = (xs[0] - Into::<T>::into(lo)) / Into::<T>::into(hi - lo);
            inv_jac /= Into::<T>::into(hi - lo);

            let lo = Into::<T>::into(settings.parameterization.input_rescaling[0][1].0);
            let hi = Into::<T>::into(settings.parameterization.input_rescaling[0][1].1);
            xs[1] = (xs[1] - Into::<T>::into(lo)) / Into::<T>::into(hi - lo);
            inv_jac /= Into::<T>::into(hi - lo);

            let lo = Into::<T>::into(settings.parameterization.input_rescaling[0][2].0);
            let hi = Into::<T>::into(settings.parameterization.input_rescaling[0][2].1);
            for x in &mut xs[2..] {
                *x -= Into::<T>::into(lo) / Into::<T>::into(hi - lo);
                inv_jac /= Into::<T>::into(hi - lo);
            }
            (xs, inv_jac)
        }
        ParameterizationMode::HyperSphericalFlat => {
            panic!("Inverse of flat hyperspherical sampling is not available since it is not bijective.");
        }
        ParameterizationMode::Cartesian | ParameterizationMode::Spherical => {
            if force_radius {
                panic!("Cannot force radius for non-hyperspherical parameterization.");
            }
            let mut inv_jac = T::one();
            let mut xs = Vec::with_capacity(moms.len() * 3);
            for (i, mom) in moms.iter().enumerate() {
                let (xs_i, inv_jac_i) = inv_parametrize3d(mom, e_cm_squared, i, settings);
                xs.extend(xs_i);
                inv_jac *= inv_jac_i;
            }
            (xs, inv_jac)
        }
    }
}

/// Map a vector in the unit hypercube to the infinite hypercube.
/// Also compute the Jacobian.
pub fn parameterize3d<T: FloatLike>(
    x: &[T],
    e_cm_squared: T,
    loop_index: usize,
    settings: &Settings,
) -> ([T; 3], T) {
    let e_cm =
        e_cm_squared.sqrt() * Into::<T>::into(settings.parameterization.shifts[loop_index].0);
    let mut l_space = [T::zero(); 3];
    let mut jac = T::one();

    // rescale the input to the desired range
    let mut x_r = [T::zero(); 3];
    for (xd, xi, &(lo, hi)) in izip!(
        &mut x_r,
        x,
        &settings.parameterization.input_rescaling[loop_index]
    ) {
        let lo = Into::<T>::into(lo);
        let hi = Into::<T>::into(hi);
        *xd = lo + *xi * (hi - lo);
        jac *= Into::<T>::into(hi - lo);
    }

    match settings.parameterization.mode {
        ParameterizationMode::Cartesian => match settings.parameterization.mapping {
            ParameterizationMapping::Log => {
                for i in 0..3 {
                    let x = x_r[i];
                    l_space[i] = e_cm * (x / (T::one() - x)).ln();
                    jac *= e_cm / (x - x * x);
                }
            }
            ParameterizationMapping::Linear => {
                for i in 0..3 {
                    let x = x_r[i];
                    l_space[i] = e_cm * (T::one() / (T::one() - x) - T::one() / x);
                    jac *=
                        e_cm * (T::one() / (x * x) + T::one() / ((T::one() - x) * (T::one() - x)));
                }
            }
        },
        ParameterizationMode::Spherical => {
            let radius = match settings.parameterization.mapping {
                ParameterizationMapping::Log => {
                    // r = e_cm * ln(1 + b*x/(1-x))
                    let x = x_r[0];
                    let b = Into::<T>::into(settings.parameterization.b);
                    let radius = e_cm * (T::one() + b * x / (T::one() - x)).ln();
                    jac *= e_cm * b / (T::one() - x) / (T::one() + x * (b - T::one()));

                    radius
                }
                ParameterizationMapping::Linear => {
                    // r = e_cm * b * x/(1-x)
                    let b = Into::<T>::into(settings.parameterization.b);
                    let radius = e_cm * b * x_r[0] / (T::one() - x_r[0]);
                    jac *= <T as num_traits::Float>::powi(e_cm * b + radius, 2) / e_cm / b;
                    radius
                }
            };
            let phi = Into::<T>::into(2.) * <T as FloatConst>::PI() * x_r[1];
            jac *= Into::<T>::into(2.) * <T as FloatConst>::PI();

            let cos_theta = -T::one() + Into::<T>::into(2.) * x_r[2]; // out of range
            jac *= Into::<T>::into(2.);
            let sin_theta = (T::one() - cos_theta * cos_theta).sqrt();

            l_space[0] = radius * sin_theta * phi.cos();
            l_space[1] = radius * sin_theta * phi.sin();
            l_space[2] = radius * cos_theta;

            jac *= radius * radius; // spherical coord
        }
        _ => {
            panic!(
                "Inappropriate parameterization mapping specified for parameterize: {:?}.",
                settings.parameterization.mode.clone()
            );
        }
    }

    // add a shift such that k=l is harder to be picked up by integrators such as cuhre
    l_space[0] += e_cm * Into::<T>::into(settings.parameterization.shifts[loop_index].1);
    l_space[1] += e_cm * Into::<T>::into(settings.parameterization.shifts[loop_index].2);
    l_space[2] += e_cm * Into::<T>::into(settings.parameterization.shifts[loop_index].3);

    (l_space, jac)
}

pub fn inv_parametrize3d<T: FloatLike>(
    mom: &LorentzVector<T>,
    e_cm_squared: T,
    loop_index: usize,
    settings: &Settings,
) -> ([T; 3], T) {
    if settings.parameterization.mode != ParameterizationMode::Spherical {
        panic!("Inverse mapping is only implemented for spherical coordinates");
    }

    let mut jac = T::one();
    let e_cm =
        e_cm_squared.sqrt() * Into::<T>::into(settings.parameterization.shifts[loop_index].0);

    let x: T = mom.x - e_cm * Into::<T>::into(settings.parameterization.shifts[loop_index].1);
    let y: T = mom.y - e_cm * Into::<T>::into(settings.parameterization.shifts[loop_index].2);
    let z: T = mom.z - e_cm * Into::<T>::into(settings.parameterization.shifts[loop_index].3);

    let k_r_sq = x * x + y * y + z * z;
    let k_r = k_r_sq.sqrt();

    let x2 = if y < T::zero() {
        T::one() + Into::<T>::into(0.5) * T::FRAC_1_PI() * T::atan2(y, x)
    } else {
        Into::<T>::into(0.5) * T::FRAC_1_PI() * T::atan2(y, x)
    };

    // cover the degenerate case
    if k_r_sq.is_zero() {
        return ([T::zero(), x2, T::zero()], T::zero());
    }

    let x1 = match settings.parameterization.mapping {
        ParameterizationMapping::Log => {
            let b = Into::<T>::into(settings.parameterization.b);
            let x1 = T::one() - b / (-T::one() + b + (k_r / e_cm).exp());
            jac /= e_cm * b / (T::one() - x1) / (T::one() + x1 * (b - T::one()));
            x1
        }
        ParameterizationMapping::Linear => {
            let b = Into::<T>::into(settings.parameterization.b);
            jac /= <T as num_traits::Float>::powi(e_cm * b + k_r, 2) / e_cm / b;
            k_r / (e_cm * b + k_r)
        }
    };

    let x3 = Into::<T>::into(0.5) * (T::one() + z / k_r);

    jac /= Into::<T>::into(2.) * <T as FloatConst>::PI();
    jac /= Into::<T>::into(2.);
    jac /= k_r * k_r;

    let mut x = [x1, x2, x3];
    for (xi, &(lo, hi)) in x
        .iter_mut()
        .zip_eq(&settings.parameterization.input_rescaling[loop_index])
    {
        *xi = (*xi - Into::<T>::into(lo)) / Into::<T>::into(hi - lo);
        jac /= Into::<T>::into(hi - lo);
    }

    (x, jac)
}

pub const MINUTE: usize = 60;
pub const HOUR: usize = 3_600;
pub const DAY: usize = 86_400;
pub const WEEK: usize = 604_800;
pub fn format_wdhms(seconds: usize) -> String {
    let mut compound_duration = vec![];
    if seconds == 0 {
        compound_duration.push(format!("{}", "0s"));
        return compound_duration.join(" ");
    }

    let mut sec = seconds % WEEK;
    // weeks
    let ws = seconds / WEEK;
    if ws != 0 {
        compound_duration.push(format!("{ws}w"));
    }

    // days
    let ds = sec / DAY;
    sec %= DAY;
    if ds != 0 {
        compound_duration.push(format!("{ds}d"));
    }

    // hours
    let hs = sec / HOUR;
    sec %= HOUR;
    if hs != 0 {
        compound_duration.push(format!("{hs}h"));
    }

    // minutes
    let ms = sec / MINUTE;
    sec %= MINUTE;
    if ms != 0 {
        compound_duration.push(format!("{ms}m"));
    }

    // seconds
    if sec != 0 {
        compound_duration.push(format!("{sec}s"));
    }

    compound_duration.join(" ")
}

#[allow(unused)]
pub fn inverse_gamma_lr(a: f64, p: f64, n_iter: usize) -> f64 {
    // this algorithm is taken from https://dl.acm.org/doi/pdf/10.1145/22721.23109

    // get an estimate for x0 to start newton iterations.
    let q = 1.0 - p;

    if 1.0 - 1.0e-8 <= a && a <= 1.0 + 1.0e-8 {
        return -q.ln();
    }

    let gamma_a = gamma(a);
    let b = q * gamma_a;
    let c = 0.57721566490153286060651209008240;

    let mut x0 = 0.5;
    if a < 1.0 {
        if b > 0.6 || (b >= 0.45 && a >= 0.3) {
            let u: f64;
            if b * q > 10e-8 {
                u = (p * gamma(a + 1.0)).powf(a.recip());
            } else {
                u = (-q / a - c).exp();
            }
            x0 = u / (1.0 - u / (a + 1.0));
        } else if a < 0.3 && 0.35 <= b && b <= 0.6 {
            let t = (-c - b).exp();
            let u = t * t.exp();
            x0 = t * u.exp();
        } else if (0.15 <= b && b <= 0.35) || (0.15 <= b && b < 0.45 && a >= 0.3) {
            let y = -b.ln();
            let u = y - (1.0 - a) * y.ln();
            x0 = y - (1.0 - a) * y.ln() - (1.0 + (1.0 - a) / (1.0 + u)).ln();
        } else if 0.01 < b && b < 0.15 {
            let y = -b.ln();
            let u = y - (1.0 - a) * y.ln();
            x0 = y
                - (1.0 - a) * u.ln()
                - ((u * u + 2.0 * (3.0 - a) * u + (2.0 - a) * (3.0 - a))
                    / (u * u + (5.0 - a) * u + 2.0))
                    .ln();
        } else if b <= 0.01 {
            let y = -b.ln();
            let c1 = (a - 1.0) * y.ln();
            let c2 = (a - 1.0) * (1.0 + c1);
            let c3 = (a - 1.0) * (-0.5 * c1 * c1 + (a - 2.0) * c1 + (3.0 * a - 5.0) * 0.5);
            let c4 = (a - 1.0)
                * (1.0 / 3.0 * c1 * c1 * c1 - (3.0 * a - 5.0) * 0.5 * c1 * c1
                    + (a * a - 6.0 * a + 7.0) * c1
                    + (11.0 * a * a - 46.0 * a + 47.0) / 6.0);
            let c5 = (a - 1.0)
                * (-0.25 * c1 * c1 * c1 * c1
                    + (11.0 * a - 7.0) / 6.0 * c1 * c1 * c1
                    + (-3.0 * a * a - 13.0) * c1 * c1
                    + (2.0 * a * a * a - 25.0 * a * a + 72.0 * a - 61.0) * 0.5 * c1
                    + (25.0 * a * a * a - 195.0 * a * a + 477.0 * a - 379.0) / 12.0);
            x0 = y + c1 + c2 / (y) + c3 / (y * y) + c4 / (y * y * y) + c5 / (y * y * y * y);

            if b <= 1.0e-28 {
                return x0;
            }
        }
    } else {
        let pref;
        let tau;
        if p < 0.5 {
            pref = -1.0;
            tau = p;
        } else {
            pref = 1.0;
            tau = q;
        }
        let t = (-2.0 * tau.ln()).sqrt();

        let a_0 = 3.31125922108741;
        let a_1 = 11.6616720288968;
        let a_2 = 4.28342155967104;
        let a_3 = 0.213623493715853;

        let b_1 = 6.61053765625462;
        let b_2 = 6.40691597760039;
        let b_3 = 1.27364489782223;
        let b_4 = 0.3611708101884203e-1;

        let t2 = t * t;
        let t3 = t2 * t;
        let t4 = t3 * t;

        let numerator = a_0 + a_1 * t + a_2 * t2 + a_3 * t3;
        let denominator = 1.0 + b_1 * t + b_2 * t2 + b_3 * t3 + b_4 * t4;

        let s = pref * (t - numerator / denominator);
        let s2 = s * s;
        let s3 = s * s2;
        let s4 = s * s3;
        let s5 = s * s4;

        let a_sqrt = a.sqrt();

        let w = a + s * a_sqrt + (s2 - 1.0) / 3.0 + (s3 - 7.0 * s) / (36.0 * a_sqrt)
            - (3.0 * s4 + 7.0 * s2 - 16.0) / (810.0 * a)
            + (9.0 * s5 + 256.0 * s3 - 433.0 * s) / (38880.0 * a * a_sqrt);

        if a >= 500.0 && (1.0 - w / a).abs() < 1.0e-6 {
            return w;
        } else {
            if p > 0.5 {
                if w < 3.0 * a {
                    x0 = w;
                } else {
                    let d = 2f64.max(a * (a - 1.0));
                    if b > 10f64.powf(-d) {
                        let u = -b.ln() + (a - 1.0) * w.ln() - (1.0 + (1.0 - a) / (1.0 + w)).ln();
                        x0 = -b.ln() + (a - 1.0) * u.ln() - (1.0 + (1.0 - a) / (1.0 + u)).ln();
                    } else {
                        let y = -b.ln();
                        let c1 = (a - 1.0) * y.ln();
                        let c2 = (a - 1.0) * (1.0 + c1);
                        let c3 =
                            (a - 1.0) * (-0.5 * c1 * c1 + (a - 2.0) * c1 + (3.0 * a - 5.0) * 0.5);
                        let c4 = (a - 1.0)
                            * (1.0 / 3.0 * c1 * c1 * c1 - (3.0 * a - 5.0) * 0.5 * c1 * c1
                                + (a * a - 6.0 * a + 7.0) * c1
                                + (11.0 * a * a - 46.0 * a + 47.0) / 6.0);
                        let c5 = (a - 1.0)
                            * (-0.25 * c1 * c1 * c1 * c1
                                + (11.0 * a - 7.0) / 6.0 * c1 * c1 * c1
                                + (-3.0 * a * a - 13.0) * c1 * c1
                                + (2.0 * a * a * a - 25.0 * a * a + 72.0 * a - 61.0) * 0.5 * c1
                                + (25.0 * a * a * a - 195.0 * a * a + 477.0 * a - 379.0) / 12.0);
                        x0 = y
                            + c1
                            + c2 / (y)
                            + c3 / (y * y)
                            + c4 / (y * y * y)
                            + c5 / (y * y * y * y);
                    }
                }
            } else {
                // this part is heavily simplified from the paper, if any issues occur this estimate
                // will need more refinement.
                let v = (p * gamma(a + 1.0)).ln();
                x0 = ((v + w) / a).exp();
            }
        }
    }

    // start iteration
    let mut x_n = x0;
    for _ in 0..n_iter {
        let r = x_n.powf(a - 1.0) * (-x_n).exp() / gamma_a;
        let t_n;
        if p <= 0.5 {
            t_n = (gamma_lr(a, x_n) - p) / r;
        } else {
            t_n = -(gamma_ur(a, x_n) - q) / r;
        }
        let w_n = (a - 1.0 - x_n) / 2.0;

        let h_n;
        if t_n.abs() <= 0.1 && (w_n * t_n).abs() <= 0.1 {
            h_n = t_n + w_n * t_n * t_n;
        } else {
            h_n = t_n;
        }

        x_n -= h_n;
    }

    return x_n;
}

#[allow(unused)]
pub fn inv_3x3_sig_matrix(mat: [[isize; 3]; 3]) -> [[isize; 3]; 3] {
    let denom = -mat[0][2] * mat[1][1] * mat[2][0]
        + mat[0][1] * mat[1][2] * mat[2][0]
        + mat[0][2] * mat[1][0] * mat[2][1]
        - mat[0][0] * mat[1][2] * mat[2][1]
        - mat[0][1] * mat[1][0] * mat[2][2]
        + mat[0][0] * mat[1][1] * mat[2][2];
    if denom != 1 && denom != -1 {
        panic!("Non invertible signature matrix.");
    }
    let mut inv_mat = [[0; 3]; 3];
    inv_mat[0][0] = (-mat[1][2] * mat[2][1] + mat[1][1] * mat[2][2]) * denom;
    inv_mat[0][1] = (mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2]) * denom;
    inv_mat[0][2] = (-mat[0][2] * mat[1][1] + mat[0][1] * mat[1][2]) * denom;
    inv_mat[1][0] = (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]) * denom;
    inv_mat[1][1] = (-mat[0][2] * mat[2][0] + mat[0][0] * mat[2][2]) * denom;
    inv_mat[1][2] = (mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2]) * denom;
    inv_mat[2][0] = (-mat[1][1] * mat[2][0] + mat[1][0] * mat[2][1]) * denom;
    inv_mat[2][1] = (mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1]) * denom;
    inv_mat[2][2] = (-mat[0][1] * mat[1][0] + mat[0][0] * mat[1][1]) * denom;

    inv_mat
}

pub fn print_banner() {
    info!(
        "\n{}{}\n",
        format!(
            "{}",
            r#"                                        _                       
                                       | |                      
   __ _  __ _ _ __ ___  _ __ ___   __ _| |     ___   ___  _ __  
  / _` |/ _` | '_ ` _ \| '_ ` _ \ / _` | |    / _ \ / _ \| '_ \ 
 | (_| | (_| | | | | | | | | | | | (_| | |___| (_) | (_) | |_) |
  \__, |\__,_|_| |_| |_|_| |_| |_|\__,_|______\___/ \___/| .__/ 
   __/ |                                                 | |    
"#
        )
        .bold()
        .blue(),
        format!(
            r#"  |___/    {}                    |_|    "#,
            format!("{:-26}", GIT_VERSION).green(),
        )
        .bold()
        .blue(),
    );
}
