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
use symbolica::domains::float::{Complex, Float};
use symbolica::id::Replacement;
use symbolica::transcendental::{hpl_at_sixth_root, sixth_root_hpl, zeta};

use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{VakintError, VakintSettings, get_integer_from_atom};

fn failure(detail: impl Into<String>) -> VakintError {
    VakintError::EvaluationError(detail.into())
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
    Ok(Some(projection(&word, sign, imaginary)?))
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

/// Normalize elementary constants while keeping real HPL projections opaque.
/// Expanding their conjugate sums before Laurent algebra would distribute a
/// real polynomial into complex terms that only cancel after rounded arithmetic.
pub(crate) fn normalize_constants(expression: AtomView) -> Result<Atom, VakintError> {
    replace_constants(expression, |source, native| {
        if let AtomView::Fun(function) = source {
            let symbol = function.get_symbol();
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
    })
}

/// Replace individual constants without evaluating the surrounding master
/// formula, tensor spectators, mass dependence, or Laurent structure.
pub(crate) fn substitute_numerically(
    expression: AtomView,
    settings: &VakintSettings,
) -> Result<Atom, VakintError> {
    let parameters = ahash::HashMap::<Atom, Complex<Float>>::default();
    replace_constants(expression, |source, native| {
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
        Ok(Atom::num(value))
    })
}

#[cfg(test)]
mod tests;
