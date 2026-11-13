//! Exact Laurent polynomials in group invariants with complex rational
//! coefficients: the values of closed colour networks.

use std::ops::{Add, AddAssign, Mul, Neg};

use symbolica::domains::rational::Rational;

/// Variables of the invariant polynomials. The named higher invariants of
/// [`super::names`] follow [`FIRST_NAMED`].
pub(crate) const NA: u16 = 0;
pub(crate) const NR: u16 = 1;
pub(crate) const I2R: u16 = 2;
pub(crate) const CR: u16 = 3;
pub(crate) const CA: u16 = 4;
pub(crate) const FIRST_NAMED: u16 = 5;

const SLOTS: usize = 8;
const EMPTY: u16 = u16::MAX;

thread_local! {
    /// Set when a monomial overflowed (more than [`SLOTS`] variables or an
    /// exponent outside i8): the polynomials computed since are garbage, and
    /// the reduction that sees the flag fails instead of panicking.
    static OVERFLOW: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Whether a monomial overflowed since the last call; clears the flag.
pub(crate) fn take_overflow() -> bool {
    OVERFLOW.with(|o| o.replace(false))
}

fn overflowed() {
    OVERFLOW.with(|o| o.set(true));
}

/// A monomial: up to eight (variable, exponent) pairs, sorted by variable.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub(crate) struct Mono([(u16, i8); SLOTS]);

impl Mono {
    pub(crate) const ONE: Mono = Mono([(EMPTY, 0); SLOTS]);

    pub(crate) fn var(var: u16, exp: i8) -> Mono {
        let mut m = Self::ONE;
        if exp != 0 {
            m.0[0] = (var, exp);
        }
        m
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = (u16, i8)> + '_ {
        self.0.iter().copied().take_while(|(v, _)| *v != EMPTY)
    }

    pub(crate) fn exponent(&self, var: u16) -> i8 {
        self.iter().find(|(v, _)| *v == var).map_or(0, |(_, e)| e)
    }

    pub(crate) fn mul(&self, other: &Mono) -> Mono {
        let mut out = Self::ONE;
        let mut n = 0;
        let (mut a, mut b) = (self.iter().peekable(), other.iter().peekable());
        loop {
            let next = match (a.peek().copied(), b.peek().copied()) {
                (None, None) => break,
                (Some(x), None) => {
                    a.next();
                    x
                }
                (None, Some(y)) => {
                    b.next();
                    y
                }
                (Some(x), Some(y)) if x.0 < y.0 => {
                    a.next();
                    x
                }
                (Some(x), Some(y)) if y.0 < x.0 => {
                    b.next();
                    y
                }
                (Some(x), Some(y)) => {
                    a.next();
                    b.next();
                    (x.0, x.1.checked_add(y.1).unwrap_or_else(|| {
                        overflowed();
                        0
                    }))
                }
            };
            if next.1 != 0 {
                if n == SLOTS {
                    overflowed();
                    break;
                }
                out.0[n] = next;
                n += 1;
            }
        }
        out
    }

    pub(crate) fn pow(&self, e: i8) -> Mono {
        let mut out = *self;
        for slot in out.0.iter_mut().take_while(|(v, _)| *v != EMPTY) {
            slot.1 = slot.1.checked_mul(e).unwrap_or_else(|| {
                overflowed();
                1
            });
        }
        if e == 0 { Self::ONE } else { out }
    }
}

/// A complex rational number.
#[derive(Clone, PartialEq, Eq, Debug)]
pub(crate) struct C {
    pub re: Rational,
    pub im: Rational,
}

impl Default for C {
    fn default() -> C {
        C::one()
    }
}

impl C {
    pub(crate) fn zero() -> C {
        C { re: Rational::zero(), im: Rational::zero() }
    }
    pub(crate) fn one() -> C {
        C { re: Rational::one(), im: Rational::zero() }
    }
    pub(crate) fn i() -> C {
        C { re: Rational::zero(), im: Rational::one() }
    }
    pub(crate) fn rat(n: i64, d: i64) -> C {
        C { re: Rational::new(n, d), im: Rational::zero() }
    }
    pub(crate) fn imag(n: i64, d: i64) -> C {
        C { re: Rational::zero(), im: Rational::new(n, d) }
    }
    pub(crate) fn is_zero(&self) -> bool {
        self.re.is_zero() && self.im.is_zero()
    }
    pub(crate) fn mul(&self, o: &C) -> C {
        if self.im.is_zero() && o.im.is_zero() {
            return C { re: &self.re * &o.re, im: Rational::zero() };
        }
        C {
            re: &(&self.re * &o.re) - &(&self.im * &o.im),
            im: &(&self.re * &o.im) + &(&self.im * &o.re),
        }
    }
    pub(crate) fn add_assign(&mut self, o: &C) {
        self.re += &o.re;
        if !o.im.is_zero() {
            self.im += &o.im;
        }
    }
    pub(crate) fn neg(&self) -> C {
        C { re: -self.re.clone(), im: -self.im.clone() }
    }
    pub(crate) fn inv(&self) -> C {
        if self.im.is_zero() {
            return C { re: self.re.inv(), im: Rational::zero() };
        }
        let d = &(&self.re * &self.re) + &(&self.im * &self.im);
        C { re: &self.re / &d, im: -(&self.im / &d) }
    }
}

/// A sparse polynomial: sorted monomials with nonzero coefficients.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub(crate) struct Poly(pub(crate) Vec<(Mono, C)>);

impl Poly {
    pub(crate) fn zero() -> Poly {
        Poly(Vec::new())
    }
    pub(crate) fn one() -> Poly {
        Poly(vec![(Mono::ONE, C::one())])
    }
    pub(crate) fn constant(c: C) -> Poly {
        if c.is_zero() { Poly::zero() } else { Poly(vec![(Mono::ONE, c)]) }
    }
    pub(crate) fn rat(n: i64, d: i64) -> Poly {
        Poly::constant(C::rat(n, d))
    }
    pub(crate) fn var(v: u16) -> Poly {
        Poly(vec![(Mono::var(v, 1), C::one())])
    }
    pub(crate) fn mono(m: Mono, c: C) -> Poly {
        if c.is_zero() { Poly::zero() } else { Poly(vec![(m, c)]) }
    }
    pub(crate) fn is_zero(&self) -> bool {
        self.0.is_empty()
    }
    /// The number if this is a constant.
    pub(crate) fn as_constant(&self) -> Option<C> {
        match self.0.as_slice() {
            [] => Some(C::zero()),
            [(m, c)] if *m == Mono::ONE => Some(c.clone()),
            _ => None,
        }
    }
    pub(crate) fn scale(&self, c: &C) -> Poly {
        if c.is_zero() {
            return Poly::zero();
        }
        Poly(self.0.iter().map(|(m, x)| (*m, x.mul(c))).collect())
    }
    pub(crate) fn mul_mono(&self, m: &Mono) -> Poly {
        let mut v: Vec<_> = self.0.iter().map(|(n, x)| (n.mul(m), x.clone())).collect();
        v.sort_by(|a, b| a.0.cmp(&b.0));
        Poly(v)
    }
    pub(crate) fn add_scaled(&mut self, other: &Poly, c: &C) {
        if c.is_zero() || other.is_zero() {
            return;
        }
        let mut out = Vec::with_capacity(self.0.len() + other.0.len());
        let mut a = std::mem::take(&mut self.0).into_iter().peekable();
        let mut b = other.0.iter().peekable();
        loop {
            match (a.peek(), b.peek()) {
                (None, None) => break,
                (Some(_), None) => out.push(a.next().unwrap()),
                (None, Some(_)) => {
                    let (m, x) = b.next().unwrap();
                    out.push((*m, x.mul(c)));
                }
                (Some((ma, _)), Some((mb, _))) => {
                    if ma < mb {
                        out.push(a.next().unwrap());
                    } else if mb < ma {
                        let (m, x) = b.next().unwrap();
                        out.push((*m, x.mul(c)));
                    } else {
                        let (m, mut x) = a.next().unwrap();
                        let (_, y) = b.next().unwrap();
                        x.add_assign(&y.mul(c));
                        if !x.is_zero() {
                            out.push((m, x));
                        }
                    }
                }
            }
        }
        self.0 = out;
    }
    pub(crate) fn mul(&self, other: &Poly) -> Poly {
        if self.is_zero() || other.is_zero() {
            return Poly::zero();
        }
        if let Some(c) = other.as_constant() {
            return self.scale(&c);
        }
        if let Some(c) = self.as_constant() {
            return other.scale(&c);
        }
        let mut terms: Vec<(Mono, C)> = Vec::with_capacity(self.0.len() * other.0.len());
        for (ma, xa) in &self.0 {
            for (mb, xb) in &other.0 {
                terms.push((ma.mul(mb), xa.mul(xb)));
            }
        }
        Self::collect(terms)
    }
    pub(crate) fn collect(mut terms: Vec<(Mono, C)>) -> Poly {
        terms.sort_by(|a, b| a.0.cmp(&b.0));
        let mut out: Vec<(Mono, C)> = Vec::with_capacity(terms.len());
        for (m, x) in terms {
            match out.last_mut() {
                Some((lm, lx)) if *lm == m => lx.add_assign(&x),
                _ => out.push((m, x)),
            }
        }
        out.retain(|(_, x)| !x.is_zero());
        Poly(out)
    }
    pub(crate) fn pow(&self, n: u32) -> Poly {
        let mut out = Poly::one();
        for _ in 0..n {
            out = out.mul(self);
        }
        out
    }
    /// Substitute a variable by a polynomial (nonnegative exponents only
    /// for non-monomial replacements).
    pub(crate) fn substitute(&self, var: u16, value: &Poly) -> Poly {
        let mut out = Poly::zero();
        for (m, x) in &self.0 {
            let e = m.exponent(var);
            if e == 0 {
                out.add_scaled(&Poly(vec![(*m, x.clone())]), &C::one());
                continue;
            }
            let rest = m.mul(&Mono::var(var, -e));
            let power = if e > 0 {
                value.pow(e as u32)
            } else {
                // Only monomial values can be inverted.
                let [(vm, vx)] = value.0.as_slice() else {
                    panic!("cannot invert a polynomial colour invariant")
                };
                Poly(vec![(vm.pow(-1), vx.inv())]).pow((-e) as u32)
            };
            out.add_scaled(&power.mul_mono(&rest), x);
        }
        out
    }
    /// Write N_A T_R as N_R C_R (the same number): closures of R lines carry
    /// the dimension and Casimir of R, as the local rules write them.
    pub(crate) fn with_fundamental_closures(&self) -> Poly {
        let terms = self
            .0
            .iter()
            .map(|(m, x)| {
                let (na, i2r) = (m.exponent(NA), m.exponent(I2R));
                let k = na.min(i2r);
                if k <= 0 {
                    return (*m, x.clone());
                }
                let shift = Mono::var(NA, -k).mul(&Mono::var(I2R, -k)).mul(&Mono::var(NR, k)).mul(&Mono::var(CR, k));
                (m.mul(&shift), x.clone())
            })
            .collect();
        Poly::collect(terms)
    }
}

impl Add<&Poly> for Poly {
    type Output = Poly;
    fn add(mut self, rhs: &Poly) -> Poly {
        self.add_scaled(rhs, &C::one());
        self
    }
}

impl AddAssign<&Poly> for Poly {
    fn add_assign(&mut self, rhs: &Poly) {
        self.add_scaled(rhs, &C::one());
    }
}

impl Mul<&Poly> for &Poly {
    type Output = Poly;
    fn mul(self, rhs: &Poly) -> Poly {
        Poly::mul(self, rhs)
    }
}

impl Neg for Poly {
    type Output = Poly;
    fn neg(self) -> Poly {
        Poly(self.0.into_iter().map(|(m, x)| (m, x.neg())).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polynomial_arithmetic_collects_and_cancels() {
        let a = Poly::var(CA) + &Poly::rat(1, 2);
        let b = Poly::var(CA) + &Poly::rat(-1, 2);
        let p = a.mul(&b);
        assert_eq!(p, Poly::var(CA).mul(&Poly::var(CA)) + &Poly::rat(-1, 4));
        let q = p.substitute(CA, &Poly::rat(1, 2));
        assert!(q.is_zero());
        let i = Poly::constant(C::i());
        assert_eq!(i.mul(&i), Poly::rat(-1, 1));
        let inv = Poly::mono(Mono::var(NA, -1), C::one()).mul(&Poly::var(NA));
        assert_eq!(inv, Poly::one());
    }

    /// A monomial outgrowing its slots or its exponents raises the overflow
    /// flag instead of panicking (a panic would abort a parallel caller).
    #[test]
    fn monomial_overflow_sets_a_flag() {
        take_overflow();
        let mut m = Mono::ONE;
        for v in 0..SLOTS as u16 {
            m = m.mul(&Mono::var(v, 1));
        }
        assert!(!take_overflow());
        let _ = m.mul(&Mono::var(SLOTS as u16, 1));
        assert!(take_overflow());
        let _ = Mono::var(CA, 100).mul(&Mono::var(CA, 100));
        assert!(take_overflow());
        let _ = Mono::var(CA, 100).pow(2);
        assert!(take_overflow());
        assert!(!take_overflow());
    }
}
