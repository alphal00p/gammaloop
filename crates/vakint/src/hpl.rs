//! Historical transcendental notation shared by FORM and RustRed evaluation.
//!
//! HR/HI denote the real/imaginary projections of HPLs at exp(±i*pi/3).
//! Symbolica evaluates words over {-1, 0, 1} through weight six at the requested
//! precision. The exact tagged constant retains its word and root sign; its
//! metadata must not be visited as ordinary numerical function arguments.
//! Real HR/HI projections remain intact during symbolic expansion, then each
//! projection is evaluated as a whole so products remain exactly real.
//! Unsupported words or arguments return errors instead of using finite tables.

use std::collections::BTreeMap;

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::coefficient::CoefficientView;
use symbolica::domains::algebraic::AlgebraicExtension;
use symbolica::domains::float::Complex;
#[cfg(test)]
use symbolica::domains::float::Float;
use symbolica::domains::rational::{Q, Rational};
use symbolica::id::Replacement;
use symbolica::transcendental::{
    GplLetter, endpoint_one_gpl, hpl_at_sixth_root, imaginary_part, polylog, real_part,
    sixth_root_hpl, zeta,
};

#[cfg(test)]
use crate::VakintSettings;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{VakintError, get_integer_from_atom};

fn failure(detail: impl Into<String>) -> VakintError {
    VakintError::EvaluationError(detail.into())
}

/// Real period projections are scalar functions, even at arity two.
pub(crate) fn is_projection_symbol(symbol: Symbol) -> bool {
    [
        "HR",
        "HI",
        "GPLReal",
        "GPLImag",
        "PolylogReal",
        "PolylogImag",
    ]
    .iter()
    .any(|name| symbol == vk_symbol!(name))
}

fn projection(word: &[i8], sign: i8, imaginary: bool) -> Result<Atom, VakintError> {
    let positive = sixth_root_hpl(word, sign).map_err(failure)?;
    let negative = sixth_root_hpl(word, -sign).map_err(failure)?;
    Ok(if imaginary {
        (positive - negative) / (Atom::num(2) * vk_parse!("𝑖").unwrap())
    } else {
        (positive + negative) / Atom::num(2)
    })
}

/// Namespaced user functions with the same short names are unrelated.
fn native_constant(part: AtomView) -> Result<Option<Atom>, VakintError> {
    if let AtomView::Var(variable) = part {
        let symbol = variable.get_symbol();
        for n in 2..=6 {
            if symbol == vk_symbol!(format!("z{n}")) {
                return Ok(Some(zeta().call(Atom::num(n))));
            }
        }
        if symbol == vk_symbol!("sqrt3") {
            return Ok(Some(Atom::num(3).pow(Atom::num((1, 2)))));
        }
        if symbol == vk_symbol!("cl2") || symbol == vk_symbol!("S2") {
            // H(1,0; exp(i*pi/3)) has imaginary part -Cl_2(pi/3).
            let clausen = -projection(&[1, 0], 1, true)?;
            return Ok(Some(if symbol == vk_symbol!("S2") {
                // Rationalize the positive radical so exact cancellations do
                // not depend on identifying sqrt(1/3) with 1/sqrt(3).
                Atom::num((4, 27)) * clausen * Atom::num(3).pow(Atom::num((1, 2)))
            } else {
                clausen
            }));
        }
        return Ok(None);
    }
    let AtomView::Fun(fun) = part else {
        return Ok(None);
    };
    let symbol = fun.get_symbol();
    if symbol == hpl_at_sixth_root() {
        // Symbolica's conservative is_constant predicate does not descend
        // through these symbolic tags. Its exact constant hook evaluates them.
        return Ok(Some(part.to_owned()));
    }
    let args = fun.iter().collect::<Vec<_>>();
    if symbol == polylog() {
        if let [order, point] = args.as_slice()
            && get_integer_from_atom(*order).is_some_and(|order| order >= 1)
            && let Some((positive, conjugate)) = polylog_pair(*point)
        {
            let real = vk_symbol!("PolylogReal").call((order.to_owned(), &positive));
            let imaginary = vk_symbol!("PolylogImag").call((order.to_owned(), positive));
            return Ok(Some(
                real + if conjugate {
                    -Atom::i() * imaginary
                } else {
                    Atom::i() * imaginary
                },
            ));
        }
        return Ok(None);
    }
    if symbol == vk_symbol!("Zeta") || symbol == zeta() {
        if let [argument] = args.as_slice()
            && get_integer_from_atom(*argument).is_some()
        {
            return Ok(Some(zeta().call(argument.to_owned())));
        }
        return Ok(None);
    }
    if symbol == vk_symbol!("Sqrt") || symbol == Symbol::SQRT {
        if let [argument] = args.as_slice()
            && get_integer_from_atom(*argument).is_some()
        {
            return Ok(Some(argument.to_owned().pow(Atom::num((1, 2)))));
        }
        return Ok(None);
    }
    Ok(projection_kernel(part)?.map(|(kernel, _)| kernel))
}

/// Return the native complex period and the real projection to take from it.
/// Tracked evaluation projects the ball itself, preserving an exactly zero
/// imaginary component without subtracting independently rounded conjugates.
pub(crate) fn projection_kernel(part: AtomView) -> Result<Option<(Atom, bool)>, VakintError> {
    let AtomView::Fun(fun) = part else {
        return Ok(None);
    };
    let symbol = fun.get_symbol();
    let args = fun.iter().collect::<Vec<_>>();
    if symbol == vk_symbol!("PolylogReal") || symbol == vk_symbol!("PolylogImag") {
        let [order, point] = args.as_slice() else {
            return Err(failure(
                "polylogarithm projection requires an order and argument",
            ));
        };
        if get_integer_from_atom(*order).is_none_or(|order| order < 1) {
            return Err(failure(
                "polylogarithm projection order must be a positive integer",
            ));
        }
        let Some((positive, conjugate)) = polylog_pair(*point) else {
            return Err(failure("unsupported polylogarithm projection argument"));
        };
        if conjugate {
            return Err(failure(
                "polylogarithm projection uses the positive-imaginary representative",
            ));
        }
        let kernel = if positive == vk_parse!("(1+3^(1/2)*𝑖)/2").unwrap() {
            // Li_n(z) = H(0,...,0,1;z). The existing sixth-root
            // HPL kernel provides bounded high-precision evaluation on the
            // unit circle without a general polylogarithm continuation.
            let order = get_integer_from_atom(*order).unwrap();
            if order > 6 {
                return Err(failure(
                    "sixth-root polylogarithm projections support orders one through six",
                ));
            }
            let mut word = vec![0; order as usize - 1];
            word.push(1);
            sixth_root_hpl(&word, 1).map_err(failure)?
        } else if get_integer_from_atom(*order) == Some(3) {
            // (1-conj(omega))/(1+conj(omega)) = i/sqrt(3).
            // Repeated integration of dLi_n(t(z))/dz =
            // -(f_1+f_-1) Li_(n-1)(t(z)), with Li_n(t(0))=zeta(n),
            // gives a short expression entirely in the bounded HPL kernel.
            let h = |word: &[i8]| sixth_root_hpl(word, -1).map_err(failure);
            let mut value = zeta().call(3) - zeta().call(2) * (h(&[1])? + h(&[-1])?);
            let log_two = Symbol::LOG.call(2);
            for first in [-1, 1] {
                for second in [-1, 1] {
                    value -= &log_two * h(&[first, second])?;
                    value -= h(&[first, second, 0])?;
                    value += h(&[first, second, -1])?;
                }
            }
            value
        } else {
            polylog().call((order.to_owned(), positive))
        };
        return Ok(Some((kernel, symbol == vk_symbol!("PolylogImag"))));
    }
    if symbol == vk_symbol!("GPLReal") || symbol == vk_symbol!("GPLImag") {
        let (endpoint, letters) = args
            .split_last()
            .ok_or_else(|| failure("GPL projection has no endpoint"))?;
        if *endpoint != Atom::num(1).as_view() {
            return Err(failure("GPL projection endpoint must be one"));
        }
        let alphabet = [
            (Atom::num(0), GplLetter::Zero),
            (Atom::num(1), GplLetter::One),
            (Atom::num(-1), GplLetter::MinusOne),
            (vk_parse!("(1+3^(1/2)*𝑖)/2").unwrap(), GplLetter::R1),
            (vk_parse!("(-1+3^(1/2)*𝑖)/2").unwrap(), GplLetter::R2),
            (vk_parse!("(-1-3^(1/2)*𝑖)/2").unwrap(), GplLetter::R4),
            (vk_parse!("(1-3^(1/2)*𝑖)/2").unwrap(), GplLetter::R5),
        ]
        .map(|(letter, tag)| (letter.expand(), tag));
        let letters = letters
            .iter()
            .map(|letter| {
                // Registered constants can spell this exact endpoint with
                // the minus sign outside the sum. Compare expanded algebra,
                // retaining exact recognition rather than a numerical test.
                let letter = letter.expand();
                alphabet
                    .iter()
                    .find_map(|(value, tag)| (value == &letter).then_some(*tag))
                    .ok_or_else(|| failure(format!("unsupported GPL letter {letter}")))
            })
            .collect::<Result<Vec<_>, _>>()?;
        return Ok(Some((
            endpoint_one_gpl(&letters).map_err(failure)?,
            symbol == vk_symbol!("GPLImag"),
        )));
    }
    let imaginary = symbol == vk_symbol!("HI");
    if !imaginary && symbol != vk_symbol!("HR") {
        return Ok(None);
    }
    let (point, indices) = args
        .split_last()
        .ok_or_else(|| failure("HPL expression has no argument"))?;
    let positive = vk_parse!("(1+3^(1/2)*𝑖)/2").unwrap();
    let negative = vk_parse!("(1-3^(1/2)*𝑖)/2").unwrap();
    let sign = if *point == Atom::var(vk_symbol!("om1")).as_view() || *point == positive.as_view() {
        1
    } else if *point == negative.as_view() {
        -1
    } else {
        return Err(failure(format!(
            "unsupported HPL argument {point}; expected exp(±i*pi/3)"
        )));
    };
    let word = indices
        .iter()
        .map(|&index| match get_integer_from_atom(index) {
            Some(value @ -1..=1) => Ok(value as i8),
            _ => Err(failure(format!(
                "invalid HPL index {index}; expected -1, 0, or 1"
            ))),
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Some((
        sixth_root_hpl(&word, sign).map_err(failure)?,
        imaginary,
    )))
}

// These are the nonreal conjugate pairs occurring in the authenticated master
// formulas. Integer polylogarithms obey conjugation here, away from the real
// branch cut. Keep each real component opaque before rounded arithmetic.
fn polylog_pair(point: AtomView) -> Option<(Atom, bool)> {
    for spelling in ["(1+3^(1/2)*𝑖)/2", "(1/3)^(1/2)*𝑖", "3^(1/2)*𝑖/3"] {
        let positive = vk_parse!(spelling).unwrap();
        let negative = if spelling.starts_with("(1+") {
            vk_parse!("(1-3^(1/2)*𝑖)/2").unwrap()
        } else {
            -&positive
        };
        let canonical = if spelling.starts_with("(1+") {
            positive.clone()
        } else {
            vk_parse!("3^(1/2)*𝑖/3").unwrap()
        };
        if point == positive.as_view() {
            return Some((canonical, false));
        }
        if point == negative.as_view() {
            return Some((canonical, true));
        }
    }
    None
}

fn replace_constants(
    expression: AtomView,
    mut materialize: impl FnMut(AtomView, Atom) -> Result<Atom, VakintError>,
) -> Result<Atom, VakintError> {
    let mut constants = BTreeMap::new();
    let mut error = None;
    expression.visitor(&mut |part| {
        if error.is_some() {
            return false;
        }
        match native_constant(part) {
            Ok(Some(native)) => {
                constants.insert(part.to_owned(), native);
                false
            }
            Ok(None) => true,
            Err(failure) => {
                error = Some(failure);
                false
            }
        }
    });
    if let Some(error) = error {
        return Err(error);
    }
    let replacements = constants
        .into_iter()
        .map(|(source, native)| {
            let value = materialize(source.as_view(), native)?;
            Ok(Replacement::new(source, value))
        })
        .collect::<Result<Vec<_>, VakintError>>()?;
    Ok(expression.replace_multiple(&replacements))
}

/// Expand exact constant denominator products before numerical evaluation.
/// Symbolica's factorized `together` intentionally preserves factors such as
/// `(sqrt(3)-i)*(sqrt(3)+i)`. Their expanded polynomial is exactly four;
/// evaluating the reciprocals separately would introduce an uncertain
/// imaginary zero into an otherwise real master coefficient.
fn normalize_constant_denominators(expression: Atom) -> Result<Atom, VakintError> {
    let mut has_float = false;
    expression.visitor(&mut |part| {
        if let AtomView::Num(number) = part {
            has_float |= matches!(number.get_coeff_view(), CoefficientView::Float(_, _));
        }
        !has_float
    });
    if has_float {
        return Ok(expression);
    }
    let mut needed = false;
    expression.visitor(&mut |part| {
        // Rational conversion treats functions, including private source
        // identities and exact evaluator tags, as opaque indeterminates.
        if matches!(part, AtomView::Fun(_)) {
            return false;
        }
        if let AtomView::Pow(power) = part
            && get_integer_from_atom(power.get_exp()).is_some_and(|exponent| exponent < 0)
            && power.get_base().is_constant()
        {
            power.get_base().visitor(&mut |base| {
                if matches!(base, AtomView::Fun(_)) {
                    return false;
                }
                if let AtomView::Num(number) = base {
                    needed |= !number.get_coeff_view().is_real();
                }
                true
            });
        }
        true
    });
    if !needed {
        return Ok(expression);
    }
    let complex = AlgebraicExtension::complex(Q);
    let polynomial = expression
        .try_to_rational_polynomial::<_, _, u32>(&complex, &complex, None)
        .map_err(|error| failure(format!("exact constant denominator: {error}")))?;
    Ok(polynomial
        .to_expression_with_coeff_map(|_, coefficient, out| {
            *out = Atom::num(Complex::new(
                coefficient.poly().get_constant(),
                coefficient
                    .poly()
                    .coefficient(&[1])
                    .unwrap_or(Rational::zero()),
            ));
        })
        .expand())
}

/// Normalize elementary constants while keeping real HPL projections opaque.
/// Expanding their conjugate sums before Laurent algebra would distribute a
/// real polynomial into complex terms that only cancel after rounded arithmetic.
pub(crate) fn normalize_constants(expression: AtomView) -> Result<Atom, VakintError> {
    // Use one exact basis for the positive radical appearing in the sunset
    // formulas. Symbolica deliberately does not infer radical identities in
    // generic complex algebra; here both rationals have their positive roots.
    let inverse_root = vk_parse!("(1/3)^(1/2)").unwrap();
    let rationalized = vk_parse!("3^(1/2)/3").unwrap();
    let expression = expression
        .replace(inverse_root.to_pattern())
        .with(rationalized.clone());
    let expression = replace_constants(expression.as_view(), |source, native| {
        if let AtomView::Fun(function) = source {
            let symbol = function.get_symbol();
            if ["GPLReal", "GPLImag", "PolylogReal", "PolylogImag"]
                .iter()
                .any(|name| symbol == vk_symbol!(name))
            {
                return Ok(source.to_owned());
            }
            if symbol == vk_symbol!("HR") || symbol == vk_symbol!("HI") {
                // native_constant has already validated the word and endpoint.
                let mut arguments = function
                    .iter()
                    .map(|part| part.to_owned())
                    .collect::<Vec<_>>();
                let point = arguments.pop().unwrap();
                let conjugate = point == vk_parse!("(1-3^(1/2)*𝑖)/2").unwrap();
                arguments.push(Atom::var(vk_symbol!("om1")));
                let projection = symbol.call(&arguments);
                return Ok(if conjugate && symbol == vk_symbol!("HI") {
                    -projection
                } else {
                    projection
                });
            }
        }
        if let AtomView::Var(variable) = source {
            let symbol = variable.get_symbol();
            if symbol == vk_symbol!("cl2") || symbol == vk_symbol!("S2") {
                let clausen = -vk_parse!("HI(1,0,om1)").unwrap();
                return Ok(if symbol == vk_symbol!("S2") {
                    Atom::num((4, 27)) * clausen * Atom::num(3).pow(Atom::num((1, 2)))
                } else {
                    clausen
                });
            }
        }
        Ok(native)
    })?;
    // Replacing legacy Sqrt(3) inside a denominator can create this radical.
    normalize_constant_denominators(
        expression
            .replace(inverse_root.to_pattern())
            .with(rationalized),
    )
}

/// Export pure constants using Symbolica's own arbitrary-precision callbacks.
/// Finite snapshot-bound coefficient identities intentionally require their
/// owning Vakint registry rather than installing mutable global callbacks.
pub(crate) fn to_symbolica_expression(expression: AtomView) -> Result<Atom, VakintError> {
    let normalized = normalize_constants(expression)?;
    replace_constants(normalized.as_view(), |source, native| {
        if let Some((kernel, imaginary)) = projection_kernel(source)? {
            Ok(if imaginary {
                imaginary_part()
            } else {
                real_part()
            }
            .call(kernel))
        } else {
            Ok(native)
        }
    })
}

/// Replace individual constants without evaluating the surrounding master
/// formula, tensor spectators, mass dependence, or Laurent structure.
#[cfg(test)]
pub(crate) fn substitute_numerically(
    expression: AtomView,
    settings: &VakintSettings,
) -> Result<Atom, VakintError> {
    let parameters = ahash::HashMap::<Atom, Complex<Float>>::default();
    replace_constants(expression, |source, native| {
        // The symbolic normalizer decomposes supported conjugate pairs into
        // real projections. This legacy one-constant evaluator can evaluate
        // the original native polylogarithm directly.
        let native = if matches!(source, AtomView::Fun(function) if function.get_symbol() == polylog())
        {
            source.to_owned()
        } else {
            native
        };
        let value = native
            .evaluate_with_prec(&parameters, settings.get_binary_precision())
            .map_err(|error| {
                failure(format!("numerical evaluation of {source} failed: {error}"))
            })?;
        if !value.re.is_finite() || !value.im.is_finite() {
            return Err(failure(format!(
                "numerical evaluation of {source} produced a non-finite value"
            )));
        }
        if let Some((_, imaginary)) = projection_kernel(source)? {
            Ok(Atom::num(if imaginary { value.im } else { value.re }))
        } else {
            Ok(Atom::num(value))
        }
    })
}

#[cfg(test)]
mod tests;
