//! Exact definitions of registered constants, independent of a data snapshot.
//!
//! Finite four-loop primitive periods belong to the external master registry.
//! Real HPL/GPL projections remain opaque until the shared numerical evaluator
//! materializes them, so surrounding algebra preserves their reality.

use crate::symbols::S;

use std::collections::BTreeSet;

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::coefficient::CoefficientView;
use symbolica::transcendental::{euler_gamma, polylog, zeta};

#[cfg(test)]
use crate::utils::vakint_macros::vk_symbol;
use crate::{VakintError, get_integer_from_atom, hpl, matad};

pub(crate) const NAMED_CONSTANTS: [&str; 12] = [
    "z2",
    "z3",
    "z4",
    "z5",
    "z6",
    "sqrt3",
    "S2",
    "cl2",
    "D6",
    "T1ep",
    "T1ep2",
    "EulerGamma",
];

/// Higher MATAD expansion coefficients exposed as named constants: the public
/// name, the `matad.rs` coefficient family and its order. Their exact values
/// are the tabulated HPL combinations of `matad::MASTERS_SUBSTITUTION`.
/// `BNbar4` is a convenience consistent with `matad.rs`; it is not
/// independent of `Li_6(1/2)`, `s6` and the classical constants.
pub(crate) const MATAD_BARS: [(&str, &str, i64); 6] = [
    ("D5bar1", "D5bar", 1),
    ("D5bar2", "D5bar", 2),
    ("D6bar1", "D6bar", 1),
    ("D6bar2", "D6bar", 2),
    ("T111bar4", "T111bar", 4),
    ("BNbar4", "BNbar", 4),
];

/// `Cl_4(π/3) = Im Li_4(e^(iπ/3))`.
fn clausen4() -> Atom {
    S.polylog_imag.call((Atom::num(4), positive_sixth_root()))
}

/// The alternating Euler sum `s6 = Σ_(m≥k≥1) (-1)^(m+k)/(m^5 k) = ζ(-5,-1) + ζ(6)`
/// (diagonal included), evaluated as `π^6/945 + G(0,0,0,0,-1,1; 1)`.
fn alternating_euler_sum6() -> Atom {
    Atom::var(Symbol::PI).pow(Atom::num(6)) / Atom::num(945)
        + S.gpl_real
            .call(&[0_i64, 0, 0, 0, -1, 1, 1].map(Atom::num)[..])
}

fn positive_sixth_root() -> Atom {
    (Atom::num(1) + Atom::num(3).pow((1, 2)) * Atom::i()) / Atom::num(2)
}

/// A positive rational other than one: its logarithm is a real constant.
fn positive_rational(argument: AtomView) -> bool {
    let AtomView::Num(number) = argument else {
        return false;
    };
    let value = number.get_coeff_view();
    matches!(
        value,
        CoefficientView::Natural(..) | CoefficientView::Large(..)
    ) && value.is_real()
        && value > CoefficientView::Natural(0, 1, 0, 1)
        && value != CoefficientView::Natural(1, 1, 0, 1)
}

fn invalid(source: AtomView, detail: &str) -> VakintError {
    VakintError::EvaluationError(format!("constant {source}: {detail}"))
}

fn matad_coefficient(name: Symbol, power: i64) -> Result<Atom, VakintError> {
    let key = name.call(power);
    matad::MASTERS_SUBSTITUTION
        .get(&key)
        .cloned()
        .ok_or_else(|| invalid(key.as_view(), "the exact MATAD definition is missing"))
}

fn sunset_coefficient(second_order: bool) -> Result<Atom, VakintError> {
    let b0 = matad_coefficient(S.matad_bars[4].1, 0)?;
    let b1 = matad_coefficient(S.matad_bars[4].1, 1)?;
    let z2 = zeta().call(2);
    let z3 = zeta().call(3);
    let inverse_sqrt3 = Atom::num(3).pow((1, 2)) / Atom::num(3);

    // MATAD has miT111 = A(ep)/( (1-2ep)(ep-1) ) *
    // (3/(2ep²) + (b0 + ep*b1 + ep²*b2 + ...)/sqrt(3)), where
    // A = exp(2γ ep)Γ(1+ep)² = 1 + ζ₂ep² - 2ζ₃ep³/3
    //     + (ζ₄+ζ₂²)ep⁴/2 + ... . T1ep/T1ep2 are coefficients of
    // miT111 itself; the external sunset master is -miT111.
    let coefficient = if second_order {
        let b2 = matad_coefficient(S.matad_bars[4].1, 2)?;
        -(Atom::num((93, 2)) + Atom::num((21, 2)) * &z2 - Atom::num(3) * z3
            + Atom::num((3, 4)) * (zeta().call(4) + z2.pow(2))
            + ((Atom::num(7) + z2) * b0 + Atom::num(3) * b1 + b2) * inverse_sqrt3)
    } else {
        -(Atom::num((45, 2)) + Atom::num((9, 2)) * z2 - z3
            + (Atom::num(3) * b0 + b1) * inverse_sqrt3)
    };
    hpl::normalize_constants(coefficient.as_view())
}

/// Return an exact definition, or `None` for an unregistered identity.
///
/// A recognized evaluator terminal can return itself. Callers must stop
/// rewriting at that fixed point. Qualified symbol identity is checked before
/// normalization; a user namespace with the same short name is unrelated.
pub(crate) fn builtin_definition(source: AtomView<'_>) -> Result<Option<Atom>, VakintError> {
    match source {
        AtomView::Var(variable) => {
            let symbol = variable.get_symbol();
            if symbol == Symbol::PI || symbol == Symbol::E || symbol == euler_gamma() {
                return Ok(Some(source.to_owned()));
            }
            if symbol == S.euler_gamma {
                return Ok(Some(Atom::var(euler_gamma())));
            }
            if symbol == S.d6 {
                // The pole 2ζ₃/ep multiplies 1/((1-2ep)(1-ep)) = 1+3ep+... .
                let value =
                    matad_coefficient(S.matad_bars[2].1, 0)? + Atom::num(6) * zeta().call(3);
                return hpl::normalize_constants(value.as_view()).map(Some);
            }
            if symbol == S.t1ep || symbol == S.t1ep2 {
                return sunset_coefficient(symbol == S.t1ep2).map(Some);
            }
            if symbol == S.cl4 {
                return Ok(Some(clausen4()));
            }
            if symbol == S.s6 {
                return hpl::normalize_constants(alternating_euler_sum6().as_view()).map(Some);
            }
            if let Some((_, family, order)) =
                S.matad_bars.iter().find(|(name, _, _)| symbol == *name)
            {
                let value = matad_coefficient(*family, *order)?;
                return hpl::normalize_constants(value.as_view()).map(Some);
            }
            if S.named_constants.contains(&symbol) {
                return hpl::normalize_constants(source).map(Some);
            }
        }
        AtomView::Fun(function) => {
            let symbol = function.get_symbol();
            if hpl::is_projection_symbol(symbol) {
                return hpl::normalize_constants(source).map(Some);
            }
            if symbol == S.zeta || symbol == zeta() {
                let arguments = function.iter().collect::<Vec<_>>();
                if let [argument] = arguments.as_slice()
                    && get_integer_from_atom(*argument).is_some_and(|order| order >= 2)
                {
                    return hpl::normalize_constants(source).map(Some);
                }
                return Err(invalid(
                    source,
                    "expected zeta at an integer order of at least two",
                ));
            }
            if symbol == S.mzv {
                let arguments = function.iter().collect::<Vec<_>>();
                // Symbolica's GPL evaluator stops at weight six; the
                // evaluable terminal sums the MZV directly at any weight.
                return integer_indices(&arguments)
                    .filter(|indices| mzv_word(indices).is_some())
                    .map(|indices| Some(crate::p_integral::mzv::mzv_value(&indices)))
                    .ok_or_else(|| {
                        invalid(
                            source,
                            "expected MZV(n_1,...,n_k) with integers n_1 >= 2, n_i >= 1",
                        )
                    });
            }
            if symbol == S.sqrt || symbol == Symbol::SQRT {
                let arguments = function.iter().collect::<Vec<_>>();
                if let [argument] = arguments.as_slice()
                    && get_integer_from_atom(*argument).is_some_and(|value| value >= 0)
                {
                    return hpl::normalize_constants(source).map(Some);
                }
                return Err(invalid(
                    source,
                    "expected the square root of a nonnegative integer",
                ));
            }
            if symbol == Symbol::LOG {
                // log(q) of a positive rational is a real constant (log 2, log 3).
                // Other logarithms are ordinary expression parts.
                let arguments = function.iter().collect::<Vec<_>>();
                return Ok(match arguments.as_slice() {
                    [argument] if positive_rational(*argument) => Some(source.to_owned()),
                    _ => None,
                });
            }
            if symbol == polylog() {
                let arguments = function.iter().collect::<Vec<_>>();
                // A public constant must be real. Ordinary polylog calls can
                // instead contain parameters, complex arguments or tracked
                // finite leaves: leave those calls to the general evaluator,
                // whose expression visitor must descend into their arguments.
                if let [order, AtomView::Num(argument)] = arguments.as_slice()
                    && let Some(order) = get_integer_from_atom(*order)
                    && order >= 1
                {
                    let value = argument.get_coeff_view();
                    let one = CoefficientView::Natural(1, 1, 0, 1);
                    if matches!(
                        value,
                        CoefficientView::Natural(..) | CoefficientView::Large(..)
                    ) && value.is_real()
                        && (value < one || (order >= 2 && value == one))
                    {
                        return Ok(Some(source.to_owned()));
                    }
                }
                return Ok(None);
            }
        }
        _ => {}
    }
    Ok(None)
}

fn integer_indices(arguments: &[AtomView]) -> Option<Vec<i64>> {
    arguments
        .iter()
        .map(|argument| get_integer_from_atom(*argument))
        .collect()
}

/// The multiple zeta value `MZV(n_1, ..., n_k) = Σ_(m_1 > ... > m_k >= 1) Π m_i^(-n_i)`
/// (the datamine convention used by Forcer) as a GPL at one:
/// `(-1)^k G(0^(n_1-1), 1, ..., 0^(n_k-1), 1; 1)`. Returns the sign and the
/// letters, or `None` for a divergent or malformed index list.
pub(crate) fn mzv_word(indices: &[i64]) -> Option<(i64, Vec<i64>)> {
    if indices.is_empty() || indices[0] < 2 || indices.iter().any(|&n| n < 1) {
        return None;
    }
    let letters = indices
        .iter()
        .flat_map(|&n| std::iter::repeat_n(0, n as usize - 1).chain([1]))
        .collect();
    Some((
        if indices.len().is_multiple_of(2) {
            1
        } else {
            -1
        },
        letters,
    ))
}

/// Discover named constants and the concrete projections in the built-in data.
/// The active registry can extend this set with calls occurring in loaded data.
pub(crate) fn builtin_symbols() -> Vec<Atom> {
    let mut symbols = S
        .named_constants
        .iter()
        .map(|name| Atom::var(*name))
        .collect::<BTreeSet<_>>();
    symbols.extend([
        Atom::var(Symbol::PI),
        Atom::var(Symbol::E),
        Atom::var(euler_gamma()),
        Atom::var(S.cl4),
        Atom::var(S.s6),
    ]);
    symbols.extend(S.matad_bars.iter().map(|(name, _, _)| Atom::var(*name)));
    for expression in matad::MASTERS_SUBSTITUTION.values() {
        expression.visitor(&mut |part| {
            if let AtomView::Fun(function) = part
                && (function.get_symbol() == S.hr || function.get_symbol() == S.hi)
            {
                symbols.insert(part.to_owned());
                return false;
            }
            true
        });
    }
    let r4 = -(Atom::one() + Atom::i() * Atom::num(3).pow((1, 2))) / Atom::num(2);
    symbols.insert(
        S.gpl_real
            .call(&[Atom::Zero, Atom::Zero, r4.clone(), Atom::one(), Atom::one()][..]),
    );
    symbols.insert(
        S.gpl_imag
            .call(&[Atom::Zero, Atom::one(), Atom::one(), r4, Atom::one()][..]),
    );
    symbols.into_iter().collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::vakint_macros::vk_parse;
    use crate::{Vakint, VakintSettings};
    use symbolica::domains::float::{Complex, Float, Real};
    use symbolica::domains::rational::Rational;
    use symbolica::id::Replacement;

    #[test]
    fn lower_loop_aliases_agree_with_the_complete_matad_prefactors() {
        Vakint::initialize_vakint_symbols();
        let replacements = matad::MASTERS_SUBSTITUTION
            .iter()
            .map(|(source, target)| Replacement::new(source.clone(), target.clone()))
            .collect::<Vec<_>>();
        for (name, master, power) in [
            ("D6", "miD6", 0),
            ("T1ep", "miT111", 1),
            ("T1ep2", "miT111", 2),
        ] {
            let source = &matad::MASTERS_EXPANSION[&Atom::var(vk_symbol!(master))].0;
            let expanded = Vakint::normalize_gamma_functions(source.as_view())
                .replace(vk_parse!("Oep(x_,y_)").unwrap().to_pattern())
                .with(Atom::Zero)
                .series(
                    vk_symbol!("ep"),
                    Atom::Zero.as_view(),
                    Rational::from(power + 1),
                )
                .unwrap()
                .to_atom()
                .replace_multiple(&replacements);
            let expected = expanded
                .coefficient_list::<i32>(&[Atom::var(vk_symbol!("ep"))])
                .into_iter()
                .find(|(monomial, _)| *monomial == Atom::var(vk_symbol!("ep")).pow(power))
                .unwrap()
                .1;
            let expected = hpl::normalize_constants(expected.as_view()).unwrap();
            let actual = builtin_definition(Atom::var(vk_symbol!(name)).as_view())
                .unwrap()
                .unwrap();
            let residual = ((actual - expected) * Atom::num(3).pow((1, 2)))
                .expand()
                .together();
            assert!(residual.is_zero(), "{name}: exact residual {residual}");
        }
    }

    #[test]
    fn recovered_constants_match_independent_numerical_values() {
        Vakint::initialize_vakint_symbols();
        let settings = VakintSettings {
            run_time_decimal_precision: 80,
            ..Default::default()
        };
        for (name, reference) in [
            (
                "D6",
                "-10.03527847976878917191470068515890023865033349600275134054452580066391412413269541906349",
            ),
            (
                "T1ep",
                "-24.20892802120359267872133821957094892580149308596048255918527447258483053931479496514833",
            ),
            (
                "T1ep2",
                "-38.71759974491583887231661364177794370941749496295656315949521883198933793534638589756375",
            ),
        ] {
            let expression = builtin_definition(Atom::var(vk_symbol!(name)).as_view())
                .unwrap()
                .unwrap();
            expression.visitor(&mut |part| {
                if let AtomView::Num(number) = part {
                    assert!(
                        !matches!(number.get_coeff_view(), CoefficientView::Float(..)),
                        "{name}"
                    );
                }
                true
            });
            let actual = hpl::substitute_numerically(expression.as_view(), &settings)
                .unwrap()
                .evaluate_with_prec(&ahash::HashMap::<Atom, Complex<Float>>::default(), 256)
                .unwrap();
            assert!(actual.im.as_raw().is_zero(), "{name}");
            assert!(
                (actual.re - Float::parse(reference, Some(256)).unwrap()).norm()
                    < Float::parse("1e-65", Some(256)).unwrap(),
                "{name}"
            );
        }
    }

    #[test]
    fn foreign_names_and_finite_primitives_are_not_builtin_definitions() {
        Vakint::initialize_vakint_symbols();
        for name in [
            "user_space::HR(1,0,vakint::om1)",
            "user_space::D6",
            "user_space::polylog(2,1/2)",
            "PR4fin",
            "PR4dfin",
            "PR11dep0",
        ] {
            assert!(
                builtin_definition(vk_parse!(name).unwrap().as_view())
                    .unwrap()
                    .is_none(),
                "{name}"
            );
        }
        for name in [
            "HR(2,om1)",
            "HI(1,0,1/2)",
            "Zeta(1)",
            "Zeta(2,3)",
            "Sqrt(-1)",
        ] {
            assert!(
                builtin_definition(vk_parse!(name).unwrap().as_view()).is_err(),
                "{name}"
            );
        }
    }

    #[test]
    fn advertised_endpoint_periods_can_be_inspected() {
        let engine = Vakint::new().unwrap();
        let mut checked = 0;
        for symbol in engine.constant_symbols().unwrap() {
            if !matches!(symbol.as_atom().as_view(), AtomView::Fun(function)
                if function.get_symbol() == vk_symbol!("GPLReal")
                    || function.get_symbol() == vk_symbol!("GPLImag"))
            {
                continue;
            }
            let constant = engine.get_constant(&symbol, 40).unwrap();
            assert!(constant.numerical_value().as_raw().is_finite(), "{symbol}");
            assert_eq!(
                constant.available_precision(),
                crate::AvailablePrecision::ArbitraryPrecision
            );
            assert!(constant.precision().imaginary().exact_zero());
            checked += 1;
        }
        assert!(
            checked >= 2,
            "both advertised endpoint projections must be covered"
        );
    }

    #[test]
    fn discovery_contains_real_projections_and_exact_aliases() {
        Vakint::initialize_vakint_symbols();
        let symbols = builtin_symbols();
        assert_eq!(symbols.iter().collect::<BTreeSet<_>>().len(), symbols.len());
        for name in [
            "D6",
            "T1ep",
            "T1ep2",
            "S2",
            "HR(1,0,1,0,om1)",
            "HI(1,0,om1)",
        ] {
            assert!(symbols.contains(&vk_parse!(name).unwrap()), "{name}");
        }
    }

    #[test]
    fn ordinary_polylog_calls_do_not_hide_argument_dependencies() {
        Vakint::initialize_vakint_symbols();
        for expression in [
            "symbolica::polylog(4,x)",
            "symbolica::polylog(4,RustRedFiniteLiteral(tag,source,0))",
            "symbolica::polylog(4,0.5`30)",
            "symbolica::polylog(2,𝑖)",
            "symbolica::polylog(2,2)",
        ] {
            assert!(
                builtin_definition(vk_parse!(expression).unwrap().as_view())
                    .unwrap()
                    .is_none(),
                "{expression}"
            );
        }
        assert!(
            builtin_definition(vk_parse!("symbolica::polylog(4,1/2)").unwrap().as_view())
                .unwrap()
                .is_some()
        );
    }
}
