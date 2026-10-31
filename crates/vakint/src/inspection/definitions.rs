//! Defining integrals of registered constants.
//!
//! A definition is an exact Symbolica atom built from inert `vakint` heads that
//! Symbolica never evaluates:
//!
//! - `Integral(f, Interval(x1, a1, b1), …, Interval(xn, an, bn))` is the
//!   iterated definite integral with the first interval outermost, as in
//!   Mathematica's `Integrate`. Inner bounds may depend on outer variables.
//! - `LoopIntegral(D, f, k(1), …, k(L))` is `∏ ∫ d^D k_i / π^(D/2) f` over
//!   Euclidean momenta, with `f` rational in `EuclideanDot` and analytic
//!   continuation in `D` where needed.
//! - `SeriesCoefficient(f, x, x0, n)` is the coefficient of `(x-x0)^n`, not the
//!   n-th derivative.
//!
//! Every other normalization (signs, `exp(L γ ep)`, projections) is explicit
//! outside these heads. Bound variables are `t`, `x(i)`, `k(i)`, `ep`
//! (`d = 4-2 ep`) and `delta` (`D = 2-2 delta`), all in the `vakint` namespace.

use std::sync::LazyLock;

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::coefficient::{Coefficient, CoefficientView};
use symbolica::domains::rational::Rational;
use symbolica::printer::PrintOptions;
use symbolica::transcendental::{
    bessel_i, bessel_k, euler_gamma, gamma, imaginary_part, polylog, real_part, zeta,
};

use super::InspectionError;
use crate::symbols::S;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{VakintError, get_integer_from_atom};

struct Heads {
    integral: Symbol,
    interval: Symbol,
    loop_integral: Symbol,
    series_coefficient: Symbol,
    t: Symbol,
    delta: Symbol,
}

static HEADS: LazyLock<Heads> = LazyLock::new(|| {
    crate::Vakint::initialize_vakint_symbols();
    Heads {
        integral: vk_symbol!("Integral"),
        interval: vk_symbol!("Interval"),
        loop_integral: vk_symbol!("LoopIntegral"),
        series_coefficient: vk_symbol!("SeriesCoefficient"),
        t: vk_symbol!("t"),
        delta: vk_symbol!("delta"),
    }
});

/// Inert heads that only definitions may use; a constant must not be named after them.
pub(crate) fn is_definition_head(symbol: Symbol) -> bool {
    let heads = &*HEADS;
    [
        heads.integral,
        heads.interval,
        heads.loop_integral,
        heads.series_coefficient,
        S.euclidean_dot,
    ]
    .contains(&symbol)
}

fn undefined(symbol: AtomView, reason: impl Into<String>) -> VakintError {
    InspectionError::UndefinedConstant {
        symbol: symbol.to_canonical_string(),
        reason: reason.into(),
    }
    .into()
}

fn call(symbol: Symbol, arguments: Vec<Atom>) -> Atom {
    symbol.call(arguments.as_slice())
}

fn t() -> Atom {
    Atom::var(HEADS.t)
}

fn x(index: usize) -> Atom {
    call(S.x, vec![Atom::num(index as i64)])
}

fn k(index: usize) -> Atom {
    call(S.k, vec![Atom::num(index as i64)])
}

fn epsilon() -> Atom {
    Atom::var(S.form_epsilon)
}

fn exp(argument: Atom) -> Atom {
    call(Symbol::EXP, vec![argument])
}

fn log(argument: Atom) -> Atom {
    call(Symbol::LOG, vec![argument])
}

fn factorial(n: u32) -> Atom {
    (2..=n).fold(Atom::num(1), |product, factor| {
        product * Atom::num(i64::from(factor))
    })
}

fn sqrt3() -> Atom {
    Atom::num(3).pow(Atom::num((1, 2)))
}

fn integral(integrand: Atom, intervals: Vec<(Atom, Atom, Atom)>) -> Atom {
    let mut arguments = vec![integrand];
    arguments.extend(
        intervals
            .into_iter()
            .map(|(variable, lower, upper)| call(HEADS.interval, vec![variable, lower, upper])),
    );
    call(HEADS.integral, arguments)
}

/// Iterated integral over 1 > x(1) > x(2) > … > x(n) > 0.
fn simplex_integral(integrand: Atom, depth: usize) -> Atom {
    integral(
        integrand,
        (1..=depth)
            .map(|index| {
                let upper = if index == 1 {
                    Atom::num(1)
                } else {
                    x(index - 1)
                };
                (x(index), Atom::Zero, upper)
            })
            .collect(),
    )
}

fn loop_integral(dimension: Atom, integrand: Atom, loops: usize) -> Atom {
    let mut arguments = vec![dimension, integrand];
    arguments.extend((1..=loops).map(k));
    call(HEADS.loop_integral, arguments)
}

fn series_coefficient(function: Atom, variable: Atom, order: i64) -> Atom {
    call(
        HEADS.series_coefficient,
        vec![function, variable, Atom::Zero, Atom::num(order)],
    )
}

fn projection(value: Atom, imaginary: bool) -> Atom {
    call(
        if imaginary {
            imaginary_part()
        } else {
            real_part()
        },
        vec![value],
    )
}

/// `π = ∫_0^1 4/(1+t²) dt`.
fn pi() -> Atom {
    integral(
        Atom::num(4) / (Atom::num(1) + t().pow(Atom::num(2))),
        vec![(t(), Atom::Zero, Atom::num(1))],
    )
}

/// `γ = -∫_0^∞ exp(-t) log(t) dt`.
fn euler_mascheroni() -> Atom {
    -integral(
        exp(-t()) * log(t()),
        vec![(t(), Atom::Zero, Atom::num(Coefficient::positive_infinity()))],
    )
}

/// `Li_n(z) = z/(n-1)! ∫_0^1 (-log t)^(n-1)/(1-z t) dt`, also `ζ(n) = Li_n(1)`.
fn polylogarithm(order: u32, point: Atom) -> Atom {
    let integrand =
        (-log(t())).pow(Atom::num(i64::from(order) - 1)) / (Atom::num(1) - point.clone() * t());
    point / factorial(order - 1) * integral(integrand, vec![(t(), Atom::Zero, Atom::num(1))])
}

/// `Cl_2(π/3) = -∫_0^(π/3) log(2 sin(t/2)) dt`.
fn clausen() -> Atom {
    -integral(
        log(Atom::num(2) * call(Symbol::SIN, vec![t() / Atom::num(2)])),
        vec![(t(), Atom::Zero, Atom::var(Symbol::PI) / Atom::num(3))],
    )
}

/// `Cl_4(π/3) = Im Li_4(e^(iπ/3)) = √3/12 ∫_0^1 (-log t)³/(1-t+t²) dt`.
fn clausen4() -> Atom {
    sqrt3() / Atom::num(12)
        * integral(
            (-log(t())).pow(Atom::num(3)) / (Atom::num(1) - t() + t().pow(Atom::num(2))),
            vec![(t(), Atom::Zero, Atom::num(1))],
        )
}

/// `log(q) = ∫_1^q dt/t` for a positive rational `q`.
fn logarithm(point: Atom) -> Atom {
    integral(Atom::num(1) / t(), vec![(t(), Atom::num(1), point)])
}

/// `H(a1,…,an; z)` at `z = (1 + sign i √3)/2` on the straight path `t = z s`.
///
/// Each letter contributes `z f_a(z s_i)` with `f_0 = 1/t`, `f_1 = 1/(1-t)`
/// and `f_-1 = 1/(1+t)`. Trailing zeros give the innermost principal
/// logarithm `log(z s_k)^m/m!`; an all-zero word is `log(z)^n/n!`.
pub(crate) fn hpl(word: &[i8], sign: i8, imaginary: bool) -> Atom {
    let z = (Atom::num(1) + Atom::num(i64::from(sign)) * Atom::i() * sqrt3()) / Atom::num(2);
    let depth = word
        .iter()
        .rposition(|letter| *letter != 0)
        .map_or(0, |p| p + 1);
    let trailing = (word.len() - depth) as u32;
    let value = if depth == 0 {
        log(z).pow(Atom::num(i64::from(trailing))) / factorial(trailing)
    } else {
        let mut integrand = Atom::num(1);
        for (index, letter) in word[..depth].iter().enumerate() {
            let s = x(index + 1);
            integrand *= match letter {
                0 => Atom::num(1) / s,
                1 => z.clone() / (Atom::num(1) - z.clone() * s),
                _ => z.clone() / (Atom::num(1) + z.clone() * s),
            };
        }
        if trailing > 0 {
            integrand = integrand * log(z.clone() * x(depth)).pow(Atom::num(i64::from(trailing)))
                / factorial(trailing);
        }
        simplex_integral(integrand, depth)
    };
    projection(value, imaginary)
}

/// `G(a1,…,an; 1) = ∫_0^1 dt1/(t1-a1) ∫_0^t1 dt2/(t2-a2) …`, without regularization.
fn gpl(letters: &[Atom], imaginary: bool) -> Result<Atom, String> {
    let first_is_one = letters
        .first()
        .is_some_and(|letter| letter.expand() == Atom::num(1));
    let last_is_zero = letters
        .last()
        .is_none_or(|letter| letter.expand().is_zero());
    if first_is_one || last_is_zero {
        return Err(
            "the iterated integral needs regularization (first letter one or last letter zero)"
                .into(),
        );
    }
    let integrand = letters
        .iter()
        .enumerate()
        .fold(Atom::num(1), |product, (index, letter)| {
            product / (x(index + 1) - letter)
        });
    Ok(projection(
        simplex_integral(integrand, letters.len()),
        imaginary,
    ))
}

/// A unit-mass Euclidean propagator, `q_E² + 1` or `q_E²` for a massless line.
fn propagator(momentum: &[i64], massive: bool) -> Atom {
    let q = momentum
        .iter()
        .enumerate()
        .filter(|(_, coefficient)| **coefficient != 0)
        .fold(Atom::Zero, |sum, (index, coefficient)| {
            sum + Atom::num(*coefficient) * k(index + 1)
        });
    let square = call(S.euclidean_dot, vec![q.clone(), q]);
    if massive {
        square + Atom::num(1)
    } else {
        square
    }
}

/// Positive Euclidean rows `D_E` of an ordered affine vacuum family.
/// A stored Minkowski row `D_M = c0 + Σ c_j q_j` becomes `D_E = Σ c_j q_j^E - c0`.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct EuclideanFamily {
    loops: usize,
    rows: Vec<Atom>,
}

impl EuclideanFamily {
    fn from_momenta(loops: usize, rows: &[(&[i64], bool)]) -> Self {
        Self {
            loops,
            rows: rows
                .iter()
                .map(|(momentum, massive)| propagator(momentum, *massive))
                .collect(),
        }
    }

    /// Read the external-data family schema (`coordinates`, `denominators`).
    pub(crate) fn from_json(family: &serde_json::Value) -> Result<Self, String> {
        if family
            .get("coordinate_index_base")
            .is_some_and(|base| *base != "0" && *base != 0)
        {
            return Err("family coordinates must use index base zero".into());
        }
        if let Some(convention) = family.get("denominator_convention")
            && *convention != "constant + sum_j coefficients[j] * coordinates[j]"
        {
            return Err(format!("unsupported denominator convention {convention}"));
        }
        if family
            .get("external_momenta")
            .and_then(|value| value.as_array())
            .is_some_and(|momenta| !momenta.is_empty())
        {
            return Err("constant definitions require a vacuum family".into());
        }
        let coordinates = family["coordinates"]
            .as_array()
            .ok_or("family has no coordinates")?
            .iter()
            .map(|coordinate| {
                let index = |key: &str| {
                    coordinate[key]
                        .as_u64()
                        .map(|value| value as usize)
                        .ok_or_else(|| format!("family coordinate has no {key} index"))
                };
                if coordinate["kind"] != "loop_loop" {
                    return Err(format!("unsupported family coordinate {coordinate}"));
                }
                Ok((index("left")?, index("right")?))
            })
            .collect::<Result<Vec<_>, String>>()?;
        let loops = coordinates
            .iter()
            .map(|(left, right)| left.max(right) + 1)
            .max()
            .ok_or("family has no coordinates")?;
        let products = coordinates
            .iter()
            .map(|(left, right)| call(S.euclidean_dot, vec![k(left + 1), k(right + 1)]))
            .collect::<Vec<_>>();
        let rational = |value: &serde_json::Value| -> Result<Atom, String> {
            let text = match value {
                serde_json::Value::String(text) => text.clone(),
                serde_json::Value::Number(number) => number.to_string(),
                _ => return Err(format!("invalid family coefficient {value}")),
            };
            let atom = vk_parse!(&text)?;
            match atom.as_view() {
                AtomView::Num(number)
                    if matches!(
                        number.get_coeff_view(),
                        CoefficientView::Natural(..) | CoefficientView::Large(..)
                    ) && number.get_coeff_view().is_real() =>
                {
                    Ok(atom)
                }
                _ => Err(format!("family coefficient {text} is not a real rational")),
            }
        };
        let rows = family["denominators"]
            .as_array()
            .ok_or("family has no denominators")?
            .iter()
            .map(|row| {
                let coefficients = row["coefficients"]
                    .as_array()
                    .ok_or("family row has no coefficients")?;
                if coefficients.len() != products.len() {
                    return Err("family row length differs from its coordinates".to_string());
                }
                let mut value = -rational(&row["constant"])?;
                for (coefficient, product) in coefficients.iter().zip(&products) {
                    value += rational(coefficient)? * product;
                }
                Ok(value)
            })
            .collect::<Result<Vec<_>, String>>()?;
        Ok(Self { loops, rows })
    }

    /// `∏ D_E^(-a)` and the Wick sign exponent `Σ a`.
    fn integrand(&self, powers: &[i64]) -> Result<(Atom, i64), String> {
        if powers.len() != self.rows.len() {
            return Err(format!(
                "{} propagator powers for a family with {} rows",
                powers.len(),
                self.rows.len()
            ));
        }
        let integrand = self
            .rows
            .iter()
            .zip(powers)
            .filter(|(_, power)| **power != 0)
            .fold(Atom::num(1), |product, (row, power)| {
                product * row.pow(Atom::num(-power))
            });
        Ok((integrand, powers.iter().sum()))
    }

    /// The positive Euclidean integral in `D = 2`, without any extra factor.
    pub(crate) fn d2_period(&self, powers: &[i64]) -> Result<Atom, String> {
        Ok(loop_integral(
            Atom::num(2),
            self.integrand(powers)?.0,
            self.loops,
        ))
    }

    /// `[δ^n] exp(L γ δ) × (positive Euclidean integral in D = 2-2δ)`.
    pub(crate) fn d2_jet(&self, powers: &[i64], order: i64) -> Result<Atom, String> {
        let delta = Atom::var(HEADS.delta);
        let dimension = Atom::num(2) - Atom::num(2) * &delta;
        let prefactor = exp(Atom::num(self.loops as i64) * Atom::var(euler_gamma()) * &delta);
        Ok(series_coefficient(
            prefactor * loop_integral(dimension, self.integrand(powers)?.0, self.loops),
            delta,
            order,
        ))
    }

    /// `[ep^n]` of a public FMFT/MATAD master in the Minkowski measure
    /// `(i π^(d/2) e^(-γ ep))^(-L)`: `(-1)^(Σ a) e^(L γ ep) LoopIntegral(4-2ep, ∏ D_E^(-a))`.
    pub(crate) fn public_master(&self, powers: &[i64], order: i64) -> Result<Atom, String> {
        let (integrand, sign) = self.integrand(powers)?;
        Ok(series_coefficient(
            Atom::num(if sign.rem_euclid(2) == 0 { 1 } else { -1 })
                * euclidean_master(integrand, self.loops),
            epsilon(),
            order,
        ))
    }
}

/// `e^(L γ ep) LoopIntegral(4-2ep, f)`: the Euclidean MATAD measure with unit masses.
fn euclidean_master(integrand: Atom, loops: usize) -> Atom {
    let dimension = Atom::num(4) - Atom::num(2) * epsilon();
    exp(Atom::num(loops as i64) * Atom::var(euler_gamma()) * epsilon())
        * loop_integral(dimension, integrand, loops)
}

/// Rows of the public four-loop families as loop-momentum combinations of k1…k4.
/// They transcribe `4L_h` (eb225dae…) and `4L_x` (693719b1…) of the external
/// data manifest, in family order; every row is massive.
const FAMILY_4L_H: [[i64; 4]; 10] = [
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [1, 0, -1, 0],
    [0, 1, -1, 0],
    [1, 0, -1, -1],
    [0, 1, -1, -1],
    [0, 0, 1, 1],
    [1, -1, 0, 0],
];
const FAMILY_4L_X: [[i64; 4]; 10] = [
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [1, 0, -1, 0],
    [0, 1, -1, 0],
    [1, 0, -1, -1],
    [0, 1, -1, -1],
    [1, 1, -1, -1],
    [0, 0, 1, 1],
];

/// FMFT master names with their representative family and powers in the public manifest.
pub(crate) const FMFT_MASTERS: [(&str, bool, [i64; 10]); 19] = [
    ("PR0", false, [1, 1, 1, 1, 1, 1, 1, 1, 1, 0]),
    ("PR1", true, [0, 0, 0, 0, 0, 1, 1, 1, 1, 0]),
    ("PR2", true, [0, 0, 0, 0, 1, 1, 1, 1, 1, 0]),
    ("PR3", true, [0, 0, 0, 1, 0, 1, 1, 1, 1, 0]),
    ("PR4", true, [0, 1, 0, 0, 1, 1, 1, 0, 1, 0]),
    ("PR4d", true, [0, 3, 0, 0, 1, 1, 1, 0, 1, 0]),
    ("PR5", true, [0, 0, 1, 0, 1, 1, 1, 1, 1, 0]),
    ("PR6", true, [0, 0, 0, 1, 1, 1, 1, 1, 1, 0]),
    ("PR7", true, [0, 1, 0, 0, 1, 1, 1, 1, 1, 0]),
    ("PR8", true, [0, 1, 0, 1, 1, 0, 1, 1, 1, 0]),
    ("PR9", true, [0, 1, 1, 0, 1, 1, 1, 1, 1, 0]),
    ("PR9d", true, [0, 2, 1, 0, 1, 1, 1, 1, 1, 0]),
    ("PR10", true, [0, 1, 0, 1, 1, 1, 1, 1, 1, 0]),
    ("PR11", true, [1, 1, 0, 1, 1, 1, 1, 1, 1, 0]),
    ("PR11d", true, [2, 1, 1, 0, 1, 1, 1, 1, 1, 0]),
    ("PR12", true, [1, 1, 1, 1, 1, 1, 1, 1, 1, 0]),
    ("PR13", true, [0, 1, 1, 1, 0, 1, 1, 1, 1, 0]),
    ("PR14", true, [0, 0, 1, 1, 1, 1, 1, 1, 1, 0]),
    ("PR15", true, [0, 1, 1, 1, 1, 1, 1, 1, 1, 0]),
];

/// Finite FMFT table literals that are complete Laurent coefficients of one master.
/// `PR3ep2` (not part of the PR3 expansion) and `PR13ep1` (truncated by an
/// `Oep` remainder) are deliberately absent.
pub(crate) const FMFT_TABLE: [(&str, &str, i64); 33] = [
    ("PR0ep0", "PR0", 0),
    ("PR0ep1", "PR0", 1),
    ("PR1ep1", "PR1", 1),
    ("PR1ep2", "PR1", 2),
    ("PR2ep1", "PR2", 1),
    ("PR2ep2", "PR2", 2),
    ("PR3ep1", "PR3", 1),
    ("PR4fin", "PR4", 0),
    ("PR4ep1", "PR4", 1),
    ("PR4dfin", "PR4d", 0),
    ("PR4dep1", "PR4d", 1),
    ("PR5ep0", "PR5", 0),
    ("PR5ep1", "PR5", 1),
    ("PR6ep0", "PR6", 0),
    ("PR6ep1", "PR6", 1),
    ("PR7ep1", "PR7", 1),
    ("PR8ep0", "PR8", 0),
    ("PR8ep1", "PR8", 1),
    ("PR9ep0", "PR9", 0),
    ("PR9ep1", "PR9", 1),
    ("PR9dep0", "PR9d", 0),
    ("PR10ep0", "PR10", 0),
    ("PR10ep1", "PR10", 1),
    ("PR11ep0", "PR11", 0),
    ("PR11ep1", "PR11", 1),
    ("PR11dep0", "PR11d", 0),
    ("PR12ep0", "PR12", 0),
    ("PR12ep1", "PR12", 1),
    ("PR13ep0", "PR13", 0),
    ("PR14ep0", "PR14", 0),
    ("PR14ep1", "PR14", 1),
    ("PR15ep0", "PR15", 0),
    ("PR15ep1", "PR15", 1),
];

/// `[ep^order]` of an FMFT four-loop master, independent of loaded data.
pub(crate) fn fmft_master(name: &str, order: i64) -> Option<Atom> {
    let (_, family_h, powers) = FMFT_MASTERS.iter().find(|(master, _, _)| *master == name)?;
    let rows = if *family_h {
        &FAMILY_4L_H
    } else {
        &FAMILY_4L_X
    };
    let rows = rows
        .iter()
        .map(|momentum| (momentum.as_slice(), true))
        .collect::<Vec<_>>();
    EuclideanFamily::from_momenta(4, &rows)
        .public_master(powers, order)
        .ok()
}

type MatadLines = &'static [([i64; 3], bool)];

/// MATAD master topologies in MATAD's own line labels (matad-ng.hh), as
/// momentum combinations of k1…k3 and a massive flag. MATAD masters are
/// Euclidean with `e^(γ ep)` per loop and no sign.
const MATAD_MASTERS: [(&str, usize, MatadLines); 10] = [
    (
        "miT111",
        2,
        &[([1, 0, 0], true), ([0, 1, 0], true), ([1, 1, 0], true)],
    ),
    // Tetrahedron, k = (p1, p3, p6): p2 = p1+p3, p4 = p1+p6, p5 = p6-p3.
    (
        "miD6",
        3,
        &[
            ([1, 0, 0], true),
            ([1, 1, 0], true),
            ([0, 1, 0], true),
            ([1, 0, 1], true),
            ([0, -1, 1], true),
            ([0, 0, 1], true),
        ],
    ),
    (
        "miD5",
        3,
        &[
            ([1, 0, 0], true),
            ([1, 1, 0], true),
            ([0, 1, 0], true),
            ([1, 0, 1], true),
            ([0, -1, 1], true),
            ([0, 0, 1], false),
        ],
    ),
    (
        "miD4",
        3,
        &[
            ([1, 0, 0], true),
            ([1, 1, 0], true),
            ([0, 1, 0], false),
            ([1, 0, 1], true),
            ([0, -1, 1], false),
            ([0, 0, 1], true),
        ],
    ),
    (
        "miDM",
        3,
        &[
            ([1, 0, 0], true),
            ([1, 1, 0], true),
            ([0, 1, 0], true),
            ([1, 0, 1], false),
            ([0, -1, 1], false),
            ([0, 0, 1], false),
        ],
    ),
    // DN labels, k = (p1, p2, p5): p3 = p5-p2, p4 = p1+p5, p6 = p1+p5-p2.
    (
        "miDN",
        3,
        &[
            ([1, 0, 0], true),
            ([0, 1, 0], true),
            ([0, -1, 1], false),
            ([1, 0, 1], false),
            ([0, 0, 1], false),
            ([1, -1, 1], false),
        ],
    ),
    // BN1 labels, k = (p1, p3, p4): p2 = p4-p3-p1, p5 = p4-p1, p6 = p1+p3.
    (
        "miD3",
        3,
        &[
            ([1, 0, 0], false),
            ([-1, -1, 1], false),
            ([0, 1, 0], true),
            ([0, 0, 1], true),
            ([-1, 0, 1], false),
            ([1, 1, 0], true),
        ],
    ),
    (
        "miBN",
        3,
        &[
            ([0, 1, 0], true),
            ([0, 0, 1], true),
            ([-1, 0, 1], true),
            ([1, 1, 0], true),
        ],
    ),
    (
        "miBN1",
        3,
        &[
            ([0, 1, 0], true),
            ([0, 0, 1], true),
            ([-1, 0, 1], false),
            ([1, 1, 0], true),
        ],
    ),
    // E3 labels, k = (p1, p3, p4): p2 = p1+p3, p5 = p4-p2.
    (
        "miE3",
        3,
        &[
            ([1, 0, 0], true),
            ([1, 1, 0], true),
            ([0, 1, 0], true),
            ([0, 0, 1], false),
            ([-1, -1, 1], false),
        ],
    ),
];

/// The Euclidean integrand and loop count of a MATAD master topology.
fn matad_topology(name: &str) -> Option<(Atom, usize)> {
    let (_, loops, lines) = MATAD_MASTERS
        .iter()
        .find(|(master, _, _)| *master == name)?;
    let integrand = lines
        .iter()
        .fold(Atom::num(1), |product, (momentum, massive)| {
            product / propagator(&momentum[..*loops], *massive)
        });
    Some((integrand, *loops))
}

/// `[ep^order]` of a MATAD master integral, `e^(L γ ep) LoopIntegral(4-2ep, …)`.
pub(crate) fn matad_master(name: &str, order: i64) -> Option<Atom> {
    let (integrand, loops) = matad_topology(name)?;
    Some(series_coefficient(
        euclidean_master(integrand, loops),
        epsilon(),
        order,
    ))
}

/// `Xbar(order)` of `matad.rs`: `[ep^order]` of the MATAD prefactor times the master.
///
/// `matad.rs` writes `miD5`, `miD6 = (Xbar(0) + ep Xbar(1) + … + 2ζ₃/ep)/((1-2ep)(1-ep))`,
/// `miBN = e^(3γep) Γ(1+ep)³ (4/ep³ - 44/(3ep²) + Σ_(k≥1) ep^k BNbar(k))/((1-ep)(1-2ep)(1-3ep)(2-3ep))`
/// and `miT111 = e^(2γep) Γ(1+ep)² (3√3/(2ep²) + Σ_(k≥0) ep^k T111bar(k))/(√3 (1-2ep)(ep-1))`.
/// With `mi = e^(L γ ep) LoopIntegral(4-2ep, …)` the explicit `e^(L γ ep)` of BN and
/// T111 cancels, and the pole terms do not contribute at the admitted orders.
pub(crate) fn matad_bar(family: &str, order: i64) -> Option<Atom> {
    let (master, prefactor, minimum_order) = matad_bar_prefactor(family)?;
    if order < minimum_order {
        return None;
    }
    let (integrand, loops) = matad_topology(master)?;
    let ep = epsilon();
    let dimension = Atom::num(4) - Atom::num(2) * &ep;
    Some(series_coefficient(
        prefactor * loop_integral(dimension, integrand, loops),
        ep,
        order,
    ))
}

/// The MATAD master, the `matad.rs` prefactor and the lowest bar order of a bar family.
fn matad_bar_prefactor(family: &str) -> Option<(&'static str, Atom, i64)> {
    let ep = epsilon();
    // a - b ep
    let linear = |a: i64, b: i64| Atom::num(a) - Atom::num(b) * &ep;
    let inverse_gamma = |power: i64| call(gamma(), vec![Atom::num(1) + &ep]).pow(Atom::num(-power));
    Some(match family {
        "D5bar" | "D6bar" => (
            if family == "D5bar" { "miD5" } else { "miD6" },
            linear(1, 2) * linear(1, 1) * exp(Atom::num(3) * Atom::var(euler_gamma()) * &ep),
            0,
        ),
        "BNbar" => (
            "miBN",
            linear(1, 1) * linear(1, 2) * linear(1, 3) * linear(2, 3) * inverse_gamma(3),
            1,
        ),
        "T111bar" => (
            "miT111",
            sqrt3() * linear(1, 2) * (ep.clone() - Atom::num(1)) * inverse_gamma(2),
            0,
        ),
        _ => return None,
    })
}

fn integer_arguments(arguments: &[AtomView]) -> Option<Vec<i64>> {
    arguments
        .iter()
        .map(|argument| get_integer_from_atom(*argument))
        .collect()
}

/// The definition of a built-in constant, or `None` for an unregistered identity.
///
/// Recognition and argument validation are those of
/// [`super::constants::builtin_definition`], so the two can never disagree.
pub(crate) fn builtin(source: AtomView<'_>) -> Result<Option<Atom>, VakintError> {
    if super::constants::builtin_definition(source)?.is_none() {
        return Ok(None);
    }
    let unsupported = || undefined(source, "no defining integral is implemented");
    let definition = match source {
        AtomView::Var(variable) => {
            let symbol = variable.get_symbol();
            let name = symbol.get_stripped_name();
            if symbol == Symbol::PI {
                pi()
            } else if symbol == Symbol::E {
                exp(Atom::num(1))
            } else if symbol == euler_gamma() || symbol == vk_symbol!("EulerGamma") {
                euler_mascheroni()
            } else if symbol == vk_symbol!("sqrt3") {
                sqrt3()
            } else if symbol == vk_symbol!("cl2") {
                clausen()
            } else if symbol == vk_symbol!("S2") {
                // S2 = 4/(9 √3) Cl_2(π/3), with the radical rationalized.
                Atom::num((4, 27)) * sqrt3() * clausen()
            } else if symbol == vk_symbol!("D6") {
                matad_master("miD6", 0).ok_or_else(unsupported)?
            } else if symbol == vk_symbol!("T1ep") {
                matad_master("miT111", 1).ok_or_else(unsupported)?
            } else if symbol == vk_symbol!("T1ep2") {
                matad_master("miT111", 2).ok_or_else(unsupported)?
            } else if symbol == vk_symbol!("cl4") {
                clausen4()
            } else if let Some((_, family, order)) = super::constants::MATAD_BARS
                .iter()
                .find(|(bar, _, _)| symbol == vk_symbol!(bar))
            {
                matad_bar(family, *order).ok_or_else(unsupported)?
            } else if let Some(order) = name
                .strip_prefix('z')
                .and_then(|order| order.parse::<u32>().ok())
            {
                polylogarithm(order, Atom::num(1))
            } else {
                return Err(unsupported());
            }
        }
        AtomView::Fun(function) => {
            let symbol = function.get_symbol();
            let arguments = function.iter().collect::<Vec<_>>();
            if symbol == vk_symbol!("HR") || symbol == vk_symbol!("HI") {
                let (kernel, imaginary) =
                    crate::hpl::projection_kernel(source)?.ok_or_else(unsupported)?;
                let AtomView::Fun(kernel_function) = kernel.as_view() else {
                    return Err(unsupported());
                };
                let tags = kernel_function.iter().collect::<Vec<_>>();
                let [AtomView::Fun(indices), sign] = tags.as_slice() else {
                    return Err(unsupported());
                };
                let word = integer_arguments(&indices.iter().collect::<Vec<_>>())
                    .ok_or_else(unsupported)?
                    .into_iter()
                    .map(|letter| letter as i8)
                    .collect::<Vec<_>>();
                let sign = get_integer_from_atom(*sign).ok_or_else(unsupported)? as i8;
                hpl(&word, sign, imaginary)
            } else if symbol == vk_symbol!("GPLReal") || symbol == vk_symbol!("GPLImag") {
                let (_, letters) = arguments.split_last().ok_or_else(unsupported)?;
                let letters = letters
                    .iter()
                    .map(|letter| letter.to_owned())
                    .collect::<Vec<_>>();
                gpl(&letters, symbol == vk_symbol!("GPLImag"))
                    .map_err(|reason| undefined(source, reason))?
            } else if symbol == vk_symbol!("PolylogReal") || symbol == vk_symbol!("PolylogImag") {
                let [order, point] = arguments.as_slice() else {
                    return Err(unsupported());
                };
                let order = get_integer_from_atom(*order).ok_or_else(unsupported)? as u32;
                projection(
                    polylogarithm(order, point.to_owned()),
                    symbol == vk_symbol!("PolylogImag"),
                )
            } else if symbol == vk_symbol!("Zeta") || symbol == zeta() {
                let order = arguments
                    .first()
                    .and_then(|order| get_integer_from_atom(*order))
                    .ok_or_else(unsupported)?;
                polylogarithm(order as u32, Atom::num(1))
            } else if symbol == vk_symbol!("Sqrt") || symbol == Symbol::SQRT {
                arguments
                    .first()
                    .ok_or_else(unsupported)?
                    .to_owned()
                    .pow(Atom::num((1, 2)))
            } else if symbol == polylog() {
                let [order, point] = arguments.as_slice() else {
                    return Err(unsupported());
                };
                let order = get_integer_from_atom(*order).ok_or_else(unsupported)? as u32;
                polylogarithm(order, point.to_owned())
            } else if symbol == Symbol::LOG {
                let [point] = arguments.as_slice() else {
                    return Err(unsupported());
                };
                logarithm(point.to_owned())
            } else {
                return Err(unsupported());
            }
        }
        _ => return Err(unsupported()),
    };
    validate(definition.as_view()).map_err(|detail| {
        VakintError::from(InspectionError::InvalidDefinition {
            symbol: source.to_canonical_string(),
            detail,
        })
    })?;
    Ok(Some(definition))
}

type NumericalTable = ahash::HashMap<
    Atom,
    (
        Atom,
        symbolica::id::Condition<symbolica::id::PatternRestriction>,
    ),
>;

/// Whether the leaf at `order` is the only content of that Laurent coefficient.
fn complete_literal(table: &NumericalTable, original: &Atom, order: i64) -> bool {
    let Some((target, _)) = table.get(original) else {
        return false;
    };
    target
        .coefficient_list::<i32>(&[Atom::var(S.form_epsilon)])
        .into_iter()
        .find(|(monomial, _)| {
            crate::precision::legacy::coefficient_power(monomial.as_view(), S.form_epsilon)
                .is_ok_and(|power| power == order)
        })
        .is_some_and(|(_, coefficient)| {
            matches!(coefficient.as_view(), AtomView::Num(number)
                if matches!(number.get_coeff_view(), CoefficientView::Float(..)))
        })
}

/// Definition of a historical FORM table leaf `NumericalSource(…)`.
///
/// FMFT and MATAD leaves denote a complete Laurent coefficient of a named
/// master; other backends are runtime inputs or estimates without a defining
/// integral and return [`InspectionError::UndefinedConstant`].
pub(crate) fn legacy_source(source: AtomView<'_>) -> Result<Option<Atom>, VakintError> {
    let Some(identity) = crate::precision::legacy::source_identity(source)? else {
        return Ok(None);
    };
    let name = match identity.original.as_view() {
        AtomView::Var(variable) => variable.get_symbol().get_stripped_name().to_owned(),
        _ => String::new(),
    };
    let definition = match identity.backend.as_str() {
        "FMFT" if ["D6", "T1ep", "T1ep2"].contains(&name.as_str()) && identity.order == 0 => {
            builtin(identity.original.as_view())?
        }
        "FMFT" => FMFT_TABLE
            .iter()
            .find(|(entry, _, _)| *entry == name)
            .filter(|_| {
                identity.ordinal == 0
                    && identity.order == 0
                    && complete_literal(
                        &crate::fmft_numerics::MASTERS_NUMERIC_SUBSTITUTIONS,
                        &identity.original,
                        identity.order,
                    )
            })
            .and_then(|(_, master, order)| fmft_master(master, *order)),
        "MATAD" => (identity.ordinal == 0
            && complete_literal(
                &crate::matad_numerics::DIRECT_SUBSTITUTIONS,
                &identity.original,
                identity.order,
            ))
        .then(|| matad_master(&name, identity.order))
        .flatten(),
        "AlphaLoop" => {
            return Err(undefined(
                source,
                "AlphaLoop table literals carry no encoded master integral or normalization",
            ));
        }
        "Input" => {
            return Err(undefined(
                source,
                "a runtime numerical input literal is not a mathematical constant",
            ));
        }
        "PySecDec" => {
            return Err(undefined(
                source,
                "a PySecDec numerical estimate has no registered defining integral",
            ));
        }
        _ => None,
    };
    let Some(definition) = definition else {
        return Err(undefined(
            source,
            format!(
                "{} table literal {} of {} at epsilon^{} is not a complete coefficient of a named master",
                identity.backend, identity.ordinal, identity.original, identity.order
            ),
        ));
    };
    Ok(Some(if identity.component == 1 {
        projection(definition, true)
    } else {
        definition
    }))
}

fn json_integer(value: &serde_json::Value) -> Option<i64> {
    value
        .as_i64()
        .or_else(|| value.as_str().and_then(|text| text.trim().parse().ok()))
}

fn json_powers(value: &serde_json::Value) -> Option<Vec<i64>> {
    value.as_array()?.iter().map(json_integer).collect()
}

/// Definitions of the external constants recorded before `symbolic_definition`
/// existed, reconstructed from the record's own definition metadata.
///
/// Returns `None` for any other symbol. Metadata that contradicts the symbol is
/// an error, never a silently different definition.
pub(crate) fn external_fallback(
    symbol: &Atom,
    provenance: &serde_json::Value,
    families: &[serde_json::Value],
) -> Result<Option<Atom>, String> {
    let metadata = &provenance["definition_metadata"];
    if let AtomView::Fun(function) = symbol.as_view() {
        let head = function.get_symbol();
        let logarithmic = head == vk_symbol!("K5LogMoment");
        if head != vk_symbol!("K5Moment") && !logarithmic {
            return Ok(None);
        }
        let arguments = function.iter().collect::<Vec<_>>();
        let power = match integer_arguments(&arguments).as_deref() {
            Some([power]) if *power >= 0 => *power,
            _ => return Err(format!("{symbol} needs one nonnegative integer power")),
        };
        if !metadata.is_null() {
            let expected = |key: &str, value: serde_json::Value| {
                metadata.get(key).is_none_or(|actual| {
                    *actual == value || (value.is_i64() && json_integer(actual) == value.as_i64())
                })
            };
            if !(expected("power", power.into())
                && expected("Bessel_power", 5.into())
                && expected("modified_Bessel_order", 0.into())
                && expected("logarithmic", logarithmic.into())
                && (!logarithmic
                    || (expected("logarithmic_moment_factor", 3.into())
                        && expected("logarithmic_factor", "log(t)+EulerGamma-log(2)".into()))))
            {
                return Err(format!("{symbol} metadata disagrees with its symbol"));
            }
        }
        let mut integrand =
            t().pow(Atom::num(power)) * call(bessel_k(), vec![Atom::Zero, t()]).pow(Atom::num(5));
        if logarithmic {
            integrand *= log(t()) + Atom::var(euler_gamma()) - log(Atom::num(2));
        }
        let value = integral(
            integrand,
            vec![(t(), Atom::Zero, Atom::num(Coefficient::positive_infinity()))],
        );
        return Ok(Some(if logarithmic {
            Atom::num(3) * value
        } else {
            value
        }));
    }
    let AtomView::Var(variable) = symbol.as_view() else {
        return Ok(None);
    };
    let jet = match variable.get_symbol() {
        name if [
            vk_symbol!("PR7D2"),
            vk_symbol!("PR10D2"),
            vk_symbol!("PR15D2"),
        ]
        .contains(&name) =>
        {
            None
        }
        name if [vk_symbol!("PR7D2Eps1"), vk_symbol!("PR10D2Eps1")].contains(&name) => Some(1),
        _ => return Ok(None),
    };
    // Jets keep their portable integral beside the Taylor metadata.
    let specification = metadata
        .get("portable_integral_definition")
        .unwrap_or(metadata);
    if let Some(order) = jet
        && (metadata
            .get("Taylor_order")
            .is_some_and(|value| json_integer(value) != Some(order))
            || metadata
                .get("Taylor_coefficient_not_factorial_scaled_derivative")
                .is_some_and(|value| *value != true))
    {
        return Err(format!(
            "{symbol} Taylor metadata disagrees with its symbol"
        ));
    }
    if specification.get("dimension").is_some_and(|dimension| {
        if jet.is_some() {
            *dimension != "2-2*delta"
        } else {
            *dimension != 2 && *dimension != "2"
        }
    }) {
        return Err(format!("{symbol} metadata has an unexpected dimension"));
    }
    let manifest_family =
        |sha: &serde_json::Value| families.iter().find(|family| family["sha256"] == *sha);
    let (family, powers) = if let Some(ordered) = specification.get("ordered_family") {
        let family = EuclideanFamily::from_json(ordered)?;
        if let Some(listed) = ordered.get("sha256").and_then(&manifest_family)
            && EuclideanFamily::from_json(listed)? != family
        {
            return Err(format!(
                "{symbol} embedded family disagrees with the manifest family of the same digest"
            ));
        }
        let powers = json_powers(&specification["signed_index_tuple"])
            .ok_or_else(|| format!("{symbol} metadata has no signed index tuple"))?;
        (family, powers)
    } else {
        let family = manifest_family(&specification["family_sha256"])
            .ok_or_else(|| format!("{symbol} metadata names no manifest family"))?;
        let powers = json_powers(&specification["powers"])
            .ok_or_else(|| format!("{symbol} metadata has no propagator powers"))?;
        (EuclideanFamily::from_json(family)?, powers)
    };
    if powers.iter().any(|power| *power < 0) {
        return Err(format!("{symbol} metadata has numerator powers"));
    }
    match jet {
        Some(order) => family.d2_jet(&powers, order),
        None => family.d2_period(&powers),
    }
    .map(Some)
}

fn deep_expand(atom: AtomView) -> Atom {
    let rebuilt = match atom {
        AtomView::Fun(function) => function.get_symbol().call(
            function
                .iter()
                .map(deep_expand)
                .collect::<Vec<_>>()
                .as_slice(),
        ),
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            deep_expand(base).pow(deep_expand(exponent))
        }
        AtomView::Mul(product) => product
            .iter()
            .fold(Atom::num(1), |result, factor| result * deep_expand(factor)),
        AtomView::Add(sum) => sum
            .iter()
            .fold(Atom::Zero, |result, term| result + deep_expand(term)),
        _ => atom.to_owned(),
    };
    rebuilt.expand()
}

/// Propagators `(v·k)² + c` of a loop integral: integer momentum vector, constant, power.
struct Propagators {
    dimension: Atom,
    loops: usize,
    rows: Vec<(Vec<Rational>, Atom, i64)>,
}

fn propagators(loop_integral: AtomView) -> Option<Propagators> {
    let AtomView::Fun(function) = loop_integral else {
        return None;
    };
    let arguments = function.iter().collect::<Vec<_>>();
    let [dimension, integrand, momenta @ ..] = arguments.as_slice() else {
        return None;
    };
    let loops = momenta.len();
    let index = |momentum: AtomView| momenta.iter().position(|candidate| *candidate == momentum);
    let factors = match integrand {
        AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
        _ => vec![*integrand],
    };
    let mut rows = Vec::new();
    for factor in factors {
        let (base, power) = match factor {
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                (base, -get_integer_from_atom(exponent)?)
            }
            _ => (factor, -1),
        };
        let terms = match base {
            AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
            _ => vec![base],
        };
        let mut quadratic = vec![vec![Rational::zero(); loops]; loops];
        let mut constant = Atom::Zero;
        for term in terms {
            let (coefficient, dot) = match term {
                AtomView::Num(_) => {
                    constant += term.to_owned();
                    continue;
                }
                AtomView::Mul(product) => {
                    let parts = product.iter().collect::<Vec<_>>();
                    let [coefficient @ AtomView::Num(_), dot] = parts.as_slice() else {
                        return None;
                    };
                    (Rational::try_from(*coefficient).ok()?, *dot)
                }
                _ => (Rational::one(), term),
            };
            let AtomView::Fun(dot) = dot else {
                return None;
            };
            let pair = dot.iter().collect::<Vec<_>>();
            if dot.get_symbol() != S.euclidean_dot || pair.len() != 2 {
                return None;
            }
            let (i, j) = (index(pair[0])?, index(pair[1])?);
            if i == j {
                quadratic[i][i] += &coefficient;
            } else {
                let half = &coefficient / &Rational::from(2);
                quadratic[i][j] += &half;
                quadratic[j][i] += &half;
            }
        }
        // A rank-one form v vᵀ with an integer vector v.
        let pivot = (0..loops).find(|i| !quadratic[*i][*i].is_zero())?;
        let root = (1..=16)
            .map(Rational::from)
            .find(|root| root * root == quadratic[pivot][pivot])?;
        let vector = (0..loops)
            .map(|j| &quadratic[pivot][j] / &root)
            .collect::<Vec<_>>();
        if (0..loops).any(|i| (0..loops).any(|j| quadratic[i][j] != &vector[i] * &vector[j])) {
            return None;
        }
        rows.push((vector, constant, power));
    }
    Some(Propagators {
        dimension: dimension.to_owned(),
        loops,
        rows,
    })
}

/// Solve `x B = w` for the coefficients of `w` in the basis rows `B`, if square and invertible.
fn inverse(matrix: &[Vec<Rational>]) -> Option<Vec<Vec<Rational>>> {
    let n = matrix.len();
    let mut work = matrix
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut row = row.clone();
            row.extend((0..n).map(|j| {
                if i == j {
                    Rational::one()
                } else {
                    Rational::zero()
                }
            }));
            row
        })
        .collect::<Vec<_>>();
    for column in 0..n {
        let pivot = (column..n).find(|row| !work[*row][column].is_zero())?;
        work.swap(pivot, column);
        let scale = work[column][column].clone();
        for entry in work[column].iter_mut() {
            *entry = &*entry / &scale;
        }
        let pivot_row = work[column].clone();
        for (row, values) in work.iter_mut().enumerate() {
            if row != column && !values[column].is_zero() {
                let factor = values[column].clone();
                for (entry, above) in values.iter_mut().zip(&pivot_row) {
                    *entry = &*entry - &(&factor * above);
                }
            }
        }
    }
    Some(work.into_iter().map(|row| row[n..].to_vec()).collect())
}

/// Row reduction over the rationals; returns the rank and the determinant.
fn eliminate(matrix: &[Vec<Rational>]) -> (usize, Rational) {
    let mut work = matrix.to_vec();
    let mut rank = 0;
    let mut determinant = Rational::one();
    let columns = work.first().map_or(0, Vec::len);
    for column in 0..columns {
        let Some(pivot) = (rank..work.len()).find(|row| !work[*row][column].is_zero()) else {
            determinant = Rational::zero();
            continue;
        };
        if pivot != rank {
            work.swap(pivot, rank);
            determinant = -determinant;
        }
        let pivot_row = work[rank].clone();
        determinant = &determinant * &pivot_row[column];
        for row in work.iter_mut().skip(rank + 1) {
            let factor = &row[column] / &pivot_row[column];
            for (entry, above) in row.iter_mut().zip(&pivot_row) {
                *entry = &*entry - &(&factor * above);
            }
        }
        rank += 1;
    }
    (rank, determinant)
}

fn determinant(matrix: &[Vec<Rational>]) -> Rational {
    eliminate(matrix).1
}

fn multiply(left: &[Vec<Rational>], right: &[Vec<Rational>]) -> Vec<Vec<Rational>> {
    left.iter()
        .map(|row| {
            (0..right[0].len())
                .map(|j| {
                    row.iter()
                        .zip(right)
                        .fold(Rational::zero(), |sum, (a, column)| sum + &(a * &column[j]))
                })
                .collect()
        })
        .collect()
}

/// Whether two loop integrals agree under a change of loop momenta `k → T k` with
/// `|det T| = 1`, which leaves the Euclidean measure invariant, up to propagator order.
fn equivalent_loop_integrals(left: &Propagators, right: &Propagators) -> bool {
    if left.dimension != right.dimension
        || left.loops != right.loops
        || left.rows.len() != right.rows.len()
    {
        return false;
    }
    // Greedily choose a basis of left rows.
    let mut basis = Vec::new();
    for (index, row) in left.rows.iter().enumerate() {
        let mut candidate = basis
            .iter()
            .map(|&chosen: &usize| left.rows[chosen].0.clone())
            .collect::<Vec<_>>();
        candidate.push(row.0.clone());
        let rank = eliminate(&candidate).0;
        if rank > basis.len() {
            basis.push(index);
        }
        if basis.len() == left.loops {
            break;
        }
    }
    if basis.len() != left.loops {
        return false;
    }
    let Some(basis_inverse) = inverse(
        &basis
            .iter()
            .map(|&row| left.rows[row].0.clone())
            .collect::<Vec<_>>(),
    ) else {
        return false;
    };
    let compatible =
        |l: usize, r: usize| left.rows[l].1 == right.rows[r].1 && left.rows[l].2 == right.rows[r].2;
    fn search(
        left: &Propagators,
        right: &Propagators,
        basis: &[usize],
        basis_inverse: &[Vec<Rational>],
        images: &mut Vec<(usize, Vec<Rational>)>,
        compatible: &dyn Fn(usize, usize) -> bool,
    ) -> bool {
        if images.len() == basis.len() {
            // Row vectors map as v ↦ v M with M = B⁻¹ W, i.e. k ↦ M k on loop momenta.
            let images_matrix = images
                .iter()
                .map(|(_, row)| row.clone())
                .collect::<Vec<_>>();
            let map = multiply(basis_inverse, &images_matrix);
            let unimodular = determinant(&map);
            if unimodular != 1 && unimodular != -1 {
                return false;
            }
            let mut used = images.iter().map(|(r, _)| *r).collect::<Vec<_>>();
            for (index, row) in left.rows.iter().enumerate() {
                if basis.contains(&index) {
                    continue;
                }
                let image = multiply(std::slice::from_ref(&row.0), &map).remove(0);
                let negated = image.iter().map(|value| -value.clone()).collect::<Vec<_>>();
                let Some(target) = (0..right.rows.len()).find(|r| {
                    !used.contains(r)
                        && compatible(index, *r)
                        && (right.rows[*r].0 == image || right.rows[*r].0 == negated)
                }) else {
                    return false;
                };
                used.push(target);
            }
            return true;
        }
        let source = basis[images.len()];
        for target in 0..right.rows.len() {
            if images.iter().any(|(used, _)| *used == target) || !compatible(source, target) {
                continue;
            }
            for sign in [1, -1] {
                let image = right.rows[target]
                    .0
                    .iter()
                    .map(|value| value * &Rational::from(sign))
                    .collect::<Vec<_>>();
                images.push((target, image));
                if search(left, right, basis, basis_inverse, images, compatible) {
                    return true;
                }
                images.pop();
            }
        }
        false
    }
    search(
        left,
        right,
        &basis,
        &basis_inverse,
        &mut Vec::new(),
        &compatible,
    )
}

fn unique_loop_integral(definition: &Atom) -> Option<Atom> {
    let mut found = Vec::new();
    definition.visitor(&mut |part| {
        if let AtomView::Fun(function) = part
            && function.get_symbol() == HEADS.loop_integral
        {
            found.push(part.to_owned());
            return false;
        }
        true
    });
    (found.len() == 1).then(|| found.remove(0))
}

/// Whether two definitions denote the same integral: the same atom after expanding
/// every argument, or the same expression around one loop integral whose propagators
/// agree under a linear change of loop momenta with unit Jacobian.
pub(crate) fn agree(left: &Atom, right: &Atom) -> bool {
    if left == right || deep_expand(left.as_view()) == deep_expand(right.as_view()) {
        return true;
    }
    let (Some(left_loop), Some(right_loop)) =
        (unique_loop_integral(left), unique_loop_integral(right))
    else {
        return false;
    };
    let (Some(left_rows), Some(right_rows)) = (
        propagators(left_loop.as_view()),
        propagators(right_loop.as_view()),
    ) else {
        return false;
    };
    if !equivalent_loop_integrals(&left_rows, &right_rows) {
        return false;
    }
    let placeholder = Atom::var(vk_symbol!("RelabeledLoopIntegral"));
    let left = left
        .replace(left_loop.to_pattern())
        .with(placeholder.to_pattern());
    let right = right
        .replace(right_loop.to_pattern())
        .with(placeholder.to_pattern());
    deep_expand(left.as_view()) == deep_expand(right.as_view())
}

/// Check that every variable is bound and only definition heads occur.
///
/// A definition never references another named constant: its only free symbols
/// are `𝜋`, `𝑒` and `γ`, and the other heads are the inert definition heads and
/// native elementary or special functions. `zeta` and `polylog` are admitted
/// only as the complete definition of that standard function at exact numbers;
/// anywhere else they must be written as their integrals.
pub(crate) fn validate(definition: AtomView) -> Result<(), String> {
    if let AtomView::Fun(function) = definition
        && [zeta(), polylog()].contains(&function.get_symbol())
    {
        return function
            .iter()
            .try_for_each(|argument| check(argument, &mut Vec::new()));
    }
    let mut scope = Vec::new();
    check(definition, &mut scope)
}

fn named_constant_hint(symbol: Symbol) -> &'static str {
    if symbol.get_namespace() == crate::NAMESPACE || [zeta(), polylog()].contains(&symbol) {
        " (a definition must not reference another named constant; write its defining integral)"
    } else {
        ""
    }
}

fn bound_variable(view: AtomView) -> Option<Atom> {
    let heads = &*HEADS;
    match view {
        AtomView::Var(variable)
            if [heads.t, heads.delta, S.form_epsilon].contains(&variable.get_symbol()) =>
        {
            Some(view.to_owned())
        }
        AtomView::Fun(function)
            if [S.x, S.k].contains(&function.get_symbol()) && function.get_nargs() == 1 =>
        {
            let index = get_integer_from_atom(function.iter().next()?)?;
            (index >= 1).then(|| view.to_owned())
        }
        _ => None,
    }
}

fn bind(variable: AtomView, scope: &[Atom]) -> Result<Atom, String> {
    let variable = bound_variable(variable)
        .ok_or_else(|| format!("{variable} is not an admissible bound variable"))?;
    if scope.contains(&variable) {
        return Err(format!("{variable} is bound twice"));
    }
    Ok(variable)
}

fn check(view: AtomView, scope: &mut Vec<Atom>) -> Result<(), String> {
    let heads = &*HEADS;
    match view {
        AtomView::Num(number) => match number.get_coeff_view() {
            CoefficientView::Float(..) => Err("definitions must be exact".into()),
            CoefficientView::Indeterminate => Err("indeterminate number".into()),
            _ => Ok(()),
        },
        AtomView::Var(variable) => {
            let symbol = variable.get_symbol();
            if [Symbol::PI, Symbol::E, euler_gamma()].contains(&symbol)
                || scope.contains(&view.to_owned())
            {
                Ok(())
            } else {
                Err(format!("free symbol {view}{}", named_constant_hint(symbol)))
            }
        }
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            check(base, scope)?;
            check(exponent, scope)
        }
        AtomView::Mul(product) => product.iter().try_for_each(|part| check(part, scope)),
        AtomView::Add(sum) => sum.iter().try_for_each(|part| check(part, scope)),
        AtomView::Fun(function) => {
            let symbol = function.get_symbol();
            let arguments = function.iter().collect::<Vec<_>>();
            let depth = scope.len();
            let result = if symbol == heads.integral {
                let Some((integrand, intervals)) = arguments.split_first() else {
                    return Err("empty Integral".into());
                };
                if intervals.is_empty() {
                    return Err("Integral without an Interval".into());
                }
                for interval in intervals {
                    let AtomView::Fun(interval) = interval else {
                        return Err(format!("{interval} is not an Interval"));
                    };
                    let parts = interval.iter().collect::<Vec<_>>();
                    let [variable, lower, upper] = parts.as_slice() else {
                        return Err("Interval needs a variable and two bounds".into());
                    };
                    if interval.get_symbol() != heads.interval {
                        return Err(format!("{} is not an Interval", interval.as_view()));
                    }
                    if matches!(variable, AtomView::Fun(f) if f.get_symbol() == S.k) {
                        return Err("loop momenta cannot be integration variables".into());
                    }
                    let variable = bind(*variable, scope)?;
                    check(*lower, scope)?;
                    check(*upper, scope)?;
                    scope.push(variable);
                }
                check(*integrand, scope)
            } else if symbol == heads.loop_integral {
                let [dimension, integrand, momenta @ ..] = arguments.as_slice() else {
                    return Err("LoopIntegral needs a dimension and an integrand".into());
                };
                if momenta.is_empty() {
                    return Err("LoopIntegral without loop momenta".into());
                }
                check(*dimension, scope)?;
                for momentum in momenta {
                    if !matches!(momentum, AtomView::Fun(f) if f.get_symbol() == S.k) {
                        return Err(format!("{momentum} is not a loop momentum vakint::k(i)"));
                    }
                    let momentum = bind(*momentum, scope)?;
                    scope.push(momentum);
                }
                check(*integrand, scope)
            } else if symbol == heads.series_coefficient {
                let [function, variable, center, order] = arguments.as_slice() else {
                    return Err("SeriesCoefficient needs four arguments".into());
                };
                if get_integer_from_atom(*order).is_none() {
                    return Err("SeriesCoefficient order must be an integer".into());
                }
                if !matches!(variable, AtomView::Var(_)) {
                    return Err("SeriesCoefficient variable must be a symbol".into());
                }
                let variable = bind(*variable, scope)?;
                check(*center, scope)?;
                scope.push(variable);
                check(*function, scope)
            } else if symbol == S.euclidean_dot {
                if arguments.len() != 2 {
                    return Err("EuclideanDot needs two momenta".into());
                }
                for momentum in &arguments {
                    if !matches!(momentum, AtomView::Fun(f) if f.get_symbol() == S.k)
                        || !scope.contains(&momentum.to_owned())
                    {
                        return Err(format!("{momentum} is not a bound loop momentum"));
                    }
                }
                Ok(())
            } else if symbol == S.x && bound_variable(view).is_some() {
                if scope.contains(&view.to_owned()) {
                    Ok(())
                } else {
                    Err(format!("free variable {view}"))
                }
            } else if [
                real_part(),
                imaginary_part(),
                Symbol::EXP,
                Symbol::LOG,
                Symbol::SIN,
                gamma(),
                bessel_k(),
                bessel_i(),
            ]
            .contains(&symbol)
            {
                arguments
                    .iter()
                    .try_for_each(|argument| check(*argument, scope))
            } else {
                Err(format!(
                    "head {} is not admitted in a definition{}",
                    symbol,
                    named_constant_hint(symbol)
                ))
            };
            scope.truncate(depth);
            result
        }
    }
}

/// Plain-language normalization statements, one per class of definition.
const PLAIN_INTEGRAL: &str = "Plain real number: the convergent integral shown, first interval outermost; no loop measure and no further factor.";
const ALGEBRAIC: &str = "Exact algebraic or elementary number; no integral is needed.";
const CLAUSEN: &str = "Clausen value at π/3, a plain real number: Cl2(π/3) = Im Li2(e^(iπ/3)) = -∫_0^(π/3) log(2 sin(t/2)) dt, Cl4(π/3) = Im Li4(e^(iπ/3)); the MATAD constant S2 = 4/(9√3) Cl2(π/3).";
const HPL: &str = "Real (HR) or imaginary (HI) part of the harmonic polylogarithm H(w; z) at z = (1 ± i√3)/2 = e^(±iπ/3): iterated integral on the straight path t = z s with 1 > s1 > … > sk > 0 (first interval outermost), letters 1/t, 1/(1-t), 1/(1+t) for 0, 1, -1, trailing zeros as log(z s)^m/m! on the principal branch; a plain number.";
const GPL: &str = "Real or imaginary part of the Goncharov polylogarithm G(a1,…,an; 1) = ∫_(1>t1>…>tn>0) ∏ dt_i/(t_i - a_i), first interval outermost, letters among 0, ±1 and the sixth roots of unity; it converges without shuffle regularization and is a plain number.";
const POLYLOG_PROJECTION: &str = "Real or imaginary part of Li_n(z) = z/(n-1)! ∫_0^1 (-log t)^(n-1)/(1 - z t) dt at the complex point shown; a plain number.";
const BESSEL_MOMENT: &str = "One-dimensional Bessel moment ∫_0^∞ t^k K_0(t)^5 dt (position-space form of the equal-mass four-loop banana in D = 2); K5LogMoment carries the factor 3 and log t + γ - log 2 as shown. A plain number.";
const BESSEL_TAYLOR: &str = "Coefficient of ε^j of the Taylor expansion (not the j-th derivative) of the ε-deformed Bessel moment (e^γ/2)^(3ε) ∫_0^∞ t^(k+3ε) K_ε(t)^5 dt; a plain number.";
/// Labels of the loop-integral normalization statements, read off the definition.
const PUBLIC_MASTER: &str = "Public Vakint/FMFT master normalization";
const MATAD_MASTER: &str = "MATAD normalization";
const EVEN_MASTER: &str = "Public Vakint/FMFT and MATAD normalization (identical for even Σa)";
const MATAD_BAR: &str = "MATAD bar";
/// The external-data acceptance tests select the D = 2 constants by these labels.
pub(crate) const D2_PERIOD: &str = "D = 2 period";
pub(crate) const D2_JET: &str = "D = 2-2δ jet";
const LOOP_INTEGRAL: &str = "Loop integral with explicit factors";

/// A loop-integral definition read off its structure:
/// `[SeriesCoefficient(] sign e^exponent other LoopIntegral(D, f, k(1), …, k(L)) [+ added, x, x0, n)]`.
struct LoopDefinition {
    /// Variable, point and order of the outer `SeriesCoefficient`.
    expansion: Option<[Atom; 3]>,
    dimension: Atom,
    loops: usize,
    /// Squared mass (constant term) and power `a` of every denominator.
    denominators: Vec<(Atom, i64)>,
    /// Whether the integrand has factors other than denominators.
    numerator: bool,
    negative: bool,
    exponent: Option<Atom>,
    other: Atom,
    /// Terms added to the loop-integral term inside the expansion.
    added: Option<Atom>,
}

/// Factors of a product with their integer exponents, distributing `(a b)^n`.
fn factor_powers(view: AtomView, power: i64, factors: &mut Vec<(Atom, i64)>) {
    match view {
        AtomView::Mul(product) => product
            .iter()
            .for_each(|factor| factor_powers(factor, power, factors)),
        AtomView::Pow(pow) => {
            let (base, exponent) = pow.get_base_exp();
            match get_integer_from_atom(exponent) {
                Some(exponent) if !matches!(base, AtomView::Num(_)) => {
                    factor_powers(base, power * exponent, factors)
                }
                _ => factors.push((view.to_owned(), power)),
            }
        }
        _ => factors.push((view.to_owned(), power)),
    }
}

/// The constant term of a propagator denominator: its squared mass.
fn squared_mass(denominator: AtomView) -> Atom {
    match denominator {
        AtomView::Add(sum) => sum
            .iter()
            .filter(|term| !term.contains_symbol(S.euclidean_dot))
            .fold(Atom::Zero, |total, term| total + term.to_owned()),
        view if view.contains_symbol(S.euclidean_dot) => Atom::Zero,
        view => view.to_owned(),
    }
}

fn parse_loop_definition(definition: AtomView) -> Option<LoopDefinition> {
    let (body, expansion) = match definition {
        AtomView::Fun(function)
            if function.get_symbol() == HEADS.series_coefficient && function.get_nargs() == 4 =>
        {
            let arguments = function.iter().collect::<Vec<_>>();
            (
                arguments[0],
                Some([
                    arguments[1].to_owned(),
                    arguments[2].to_owned(),
                    arguments[3].to_owned(),
                ]),
            )
        }
        view => (view, None),
    };
    let terms = match body {
        AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
        view => vec![view],
    };
    let (integral_terms, added_terms): (Vec<_>, Vec<_>) = terms
        .into_iter()
        .partition(|term| term.contains_symbol(HEADS.loop_integral));
    let [term] = integral_terms[..] else {
        return None;
    };
    let added = (!added_terms.is_empty()).then(|| {
        added_terms
            .iter()
            .fold(Atom::Zero, |total, term| total + term.to_owned())
    });
    let mut factors = Vec::new();
    factor_powers(term, 1, &mut factors);
    let (mut integral, mut negative, mut exponent, mut other) =
        (None, false, None::<Atom>, Atom::num(1));
    for (factor, power) in factors {
        let mut add_exponent = |summand: Atom| {
            exponent = Some(exponent.take().unwrap_or(Atom::Zero) + summand * Atom::num(power));
        };
        match factor.as_view() {
            AtomView::Fun(function) if function.get_symbol() == HEADS.loop_integral => {
                if power != 1 || integral.is_some() {
                    return None;
                }
                integral = Some(function.iter().map(|a| a.to_owned()).collect::<Vec<_>>());
            }
            view if view.contains_symbol(HEADS.loop_integral) => return None,
            AtomView::Num(_) if factor == Atom::num(-1) => negative ^= power % 2 != 0,
            AtomView::Var(variable) if variable.get_symbol() == Symbol::E => {
                add_exponent(Atom::num(1))
            }
            AtomView::Pow(pow) if matches!(pow.get_base_exp().0, AtomView::Var(base) if base.get_symbol() == Symbol::E) => {
                add_exponent(pow.get_base_exp().1.to_owned())
            }
            AtomView::Fun(function)
                if function.get_symbol() == Symbol::EXP && function.get_nargs() == 1 =>
            {
                add_exponent(function.iter().next().unwrap().to_owned())
            }
            _ => other *= factor.pow(Atom::num(power)),
        }
    }
    let arguments = integral?;
    if arguments.len() < 3 {
        return None;
    }
    let mut integrand = Vec::new();
    factor_powers(arguments[1].as_view(), 1, &mut integrand);
    let mut denominators = Vec::new();
    let mut numerator = false;
    for (factor, power) in integrand {
        if power < 0 && factor.contains_symbol(S.euclidean_dot) {
            denominators.push((squared_mass(factor.as_view()), -power));
        } else if factor != Atom::num(1) {
            numerator = true;
        }
    }
    Some(LoopDefinition {
        expansion,
        dimension: arguments[0].clone(),
        loops: arguments.len() - 2,
        denominators,
        numerator,
        negative,
        exponent,
        other,
        added,
    })
}

/// The power `m` if `term = c x^m` with `c` free of `x`.
fn monomial_power(term: AtomView, variable: &Atom) -> Option<i64> {
    let mut factors = Vec::new();
    factor_powers(term, 1, &mut factors);
    factors.into_iter().try_fold(0, |total, (factor, power)| {
        if factor == *variable {
            Some(total + power)
        } else if factor.contains(variable.as_view()) {
            None
        } else {
            Some(total)
        }
    })
}

fn equal(left: &Atom, right: &Atom) -> bool {
    (left.clone() - right.clone()).expand().is_zero()
}

/// The normalization statement of a loop-integral definition, derived from its
/// structure: the coefficient taken, the measure and dimension, the masses and
/// `Σa`, and every factor in front of the `LoopIntegral`. The label names the
/// convention that these factors realise; it never depends on the symbol name.
fn loop_statement(definition: &Atom) -> Option<String> {
    let parsed = parse_loop_definition(definition.as_view())?;
    let show = |atom: &Atom| atom.printer(PrintOptions::file_no_namespace()).to_string();
    let loops = parsed.loops;
    let propagators = parsed.denominators.len();
    let power_sum: i64 = parsed.denominators.iter().map(|(_, power)| power).sum();
    let odd = power_sum.rem_euclid(2) == 1;
    let variable = parsed.expansion.as_ref().map(|[variable, _, _]| variable);
    let per_loop = |variable: &Atom| {
        parsed.exponent.as_ref().is_some_and(|exponent| {
            equal(
                exponent,
                &(Atom::num(loops as i64) * Atom::var(euler_gamma()) * variable),
            )
        })
    };
    let plain = parsed.other == Atom::num(1) && parsed.added.is_none();
    let prefactor = parsed
        .exponent
        .as_ref()
        .map_or(Atom::num(1), |exponent| exp(exponent.clone()))
        * &parsed.other
        * Atom::num(if parsed.negative { -1 } else { 1 });
    let matad = ["D5bar", "BNbar", "T111bar"]
        .into_iter()
        .find_map(|family| {
            let (master, candidate, _) = matad_bar_prefactor(family)?;
            (matad_topology(master)?.1 == loops && equal(&prefactor, &candidate)).then_some(family)
        });
    let dimension_is = |expected: Atom| equal(&parsed.dimension, &expected);
    let label = match (matad, variable) {
        (Some(family), _) => format!(
            "{MATAD_BAR} (matad.rs {} prefactor)",
            match family {
                "D5bar" => "D5/D6",
                "BNbar" => "BN",
                _ => "T111",
            }
        ),
        (None, None)
            if dimension_is(Atom::num(2))
                && plain
                && parsed.exponent.is_none()
                && !parsed.negative =>
        {
            D2_PERIOD.into()
        }
        (None, Some(variable))
            if plain
                && per_loop(variable)
                && variable == &Atom::var(HEADS.delta)
                && dimension_is(Atom::num(2) - Atom::num(2) * variable)
                && !parsed.negative =>
        {
            D2_JET.into()
        }
        (None, Some(variable))
            if plain
                && per_loop(variable)
                && dimension_is(Atom::num(4) - Atom::num(2) * variable) =>
        {
            match (odd, parsed.negative) {
                (true, true) => PUBLIC_MASTER.into(),
                (true, false) => MATAD_MASTER.into(),
                (false, false) => EVEN_MASTER.into(),
                (false, true) => LOOP_INTEGRAL.into(),
            }
        }
        _ => LOOP_INTEGRAL.into(),
    };

    let mut text = format!("{label}: ");
    match &parsed.expansion {
        Some([variable, point, order]) => text.push_str(&format!(
            "coefficient of {}^{} of the series expansion in {} around {} (a series coefficient, not a derivative) of ",
            show(variable),
            show(order),
            show(variable),
            show(point)
        )),
        None => text.push_str("exactly "),
    }
    let mut factors = Vec::new();
    if parsed.negative {
        factors.push(if odd {
            "the sign -1 = (-1)^Σa".to_string()
        } else {
            "the sign -1".to_string()
        });
    }
    if let Some(exponent) = &parsed.exponent {
        factors.push(match variable {
            Some(variable) if per_loop(variable) => format!(
                "e^({}), i.e. e^(γ {}) per loop",
                show(exponent),
                show(variable)
            ),
            _ => format!("e^({})", show(exponent)),
        });
    }
    if parsed.other != Atom::num(1) {
        factors.push(format!("the factor {}", show(&parsed.other)));
    }
    for factor in factors {
        text.push_str(&factor);
        text.push_str(" times ");
    }
    let unit = parsed
        .denominators
        .iter()
        .filter(|(mass, _)| *mass == Atom::num(1))
        .count();
    let massless = parsed
        .denominators
        .iter()
        .filter(|(mass, _)| mass.is_zero())
        .count();
    let masses = if unit == propagators {
        "all with unit mass".to_string()
    } else if unit + massless == propagators {
        format!("{unit} with unit mass and {massless} massless")
    } else {
        "with the squared masses written in them".to_string()
    };
    text.push_str(&format!(
        "the {loops}-loop Euclidean integral with measure d^Dk/π^(D/2) per loop, D = {}, over {propagators} propagators, {masses}, whose powers sum to Σa = {power_sum}{}",
        show(&parsed.dimension),
        if parsed.numerator {
            ", with the numerator shown"
        } else {
            ""
        }
    ));
    if let Some(added) = &parsed.added {
        text.push_str(&format!(
            "; the term {} is added inside the expansion",
            show(added)
        ));
        let contributes = match &parsed.expansion {
            Some([variable, point, order]) if point.is_zero() => {
                let terms = match added.as_view() {
                    AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
                    view => vec![view],
                };
                terms.into_iter().any(|term| {
                    monomial_power(term, variable)
                        .is_none_or(|power| Some(power) == get_integer_from_atom(order.as_view()))
                })
            }
            _ => true,
        };
        if !contributes {
            text.push_str(" (it does not contribute to this coefficient)");
        }
    }
    text.push_str(". ");
    if !parsed.negative && odd {
        text.push_str("No sign (-1)^Σa. ");
    }
    if let Some(variable) = variable
        && parsed.exponent.is_none()
    {
        text.push_str(&format!("No e^(γ {}) factor", show(variable)));
        text.push_str(if matches!(matad, Some("BNbar" | "T111bar")) {
            ": the explicit e^(-Lγε) of matad.rs has cancelled the e^(γε) per loop of the MATAD master. "
        } else {
            ". "
        });
    }
    if label == PUBLIC_MASTER {
        text.push_str(
            "This is the Euclidean form of the Minkowski measure (iπ^(d/2) e^(-γε))^(-1) d^dk per loop. ",
        );
    } else if label == MATAD_MASTER {
        text.push_str(
            "The public Vakint/FMFT master with the same propagators carries in addition the sign (-1)^Σa. ",
        );
    } else if label == D2_PERIOD {
        text.push_str(
            "No regulator: the integral is finite as written (positive, hat convention). ",
        );
    } else if label == D2_JET {
        text.push_str("Positive (hat) sign convention. ");
    }
    text.push_str("No other factor.");
    Some(text)
}

fn named(symbol: AtomView, names: &[&str]) -> bool {
    let head = match symbol {
        AtomView::Var(variable) => variable.get_symbol(),
        AtomView::Fun(function) => function.get_symbol(),
        _ => return false,
    };
    names.iter().any(|name| head == vk_symbol!(name))
}

/// A plain-language statement of the normalization conventions of a definition.
///
/// For a loop-integral definition the statement is derived from the definition
/// itself (see [`loop_statement`]), so it cannot contradict it. The other
/// classes (elementary integrals, HPL/GPL projections, Bessel moments) carry a
/// fixed statement. Without a definition it records why.
pub(crate) fn normalisation(
    symbol: AtomView,
    definition: Option<&Atom>,
    unavailable: Option<&str>,
) -> String {
    let Some(definition) = definition else {
        return format!(
            "No defining integral is registered ({}); the numerical value is the recorded source value.",
            unavailable.unwrap_or("no definition")
        );
    };
    if definition.contains_symbol(HEADS.loop_integral) {
        return loop_statement(definition).unwrap_or_else(|| {
            format!(
                "{LOOP_INTEGRAL}: the expression shown, with the Euclidean measure d^Dk/π^(D/2) per loop inside each LoopIntegral and every other factor written out."
            )
        });
    }
    if definition.contains_symbol(bessel_k()) {
        return if definition.contains_symbol(HEADS.series_coefficient) {
            BESSEL_TAYLOR.into()
        } else {
            BESSEL_MOMENT.into()
        };
    }
    if named(symbol, &["HR", "HI"]) {
        return HPL.into();
    }
    if named(symbol, &["GPLReal", "GPLImag"]) {
        return GPL.into();
    }
    if named(symbol, &["PolylogReal", "PolylogImag"]) {
        return POLYLOG_PROJECTION.into();
    }
    if named(symbol, &["cl2", "cl4", "S2"]) {
        return CLAUSEN.into();
    }
    if !definition.contains_symbol(HEADS.integral) {
        return ALGEBRAIC.into();
    }
    PLAIN_INTEGRAL.into()
}

#[cfg(test)]
#[path = "definitions_tests.rs"]
mod tests;
