//! Exact comparison of reduced p-integrals in the tag basis.
//!
//! Forcer writes its reductions with the G-scheme constants `PG(x1,x2)` and
//! the masters `PMaster(name)`; RustRed writes Gamma tags `PGamma(m)` and the
//! same masters. Every `PG(x1,x2)` is a rational function of `ε` times a
//! product of `Γ(1+cε)` with integer `c` (`finalize::pg_value`), so the common
//! exact basis is the monomials `Π Γ(1+mε)^e × Π PMaster(name)^k` with
//! coefficients rational in `ε` (and in whatever the user's coefficients
//! hold). Forcer's own relations between `PG`s (`Masters.prc`, the `Masters-1`
//! section) are identities in this basis.
//!
//! The basis assumes that the `Γ(1+mε)` are algebraically independent over
//! `Q(ε)`. A difference that is not zero here while the numerical values
//! agree is a Γ identity, never a false agreement (M2 plan, section 10,
//! risk 22).

use std::collections::{BTreeMap, BTreeSet};

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::transcendental::gamma;

use super::finalize::{gamma_one_plus, pg_value};
use super::laurent::unsupported_power;
use super::{masters, p_gamma_symbol, p_master_symbol, pg_symbol};
use crate::symbols::S;
use crate::{VakintError, get_integer_from_atom};

/// `Π Γ(1+mε)^e`: `(m, e)` with `m ≠ 0`, `e ≠ 0`, sorted by `m`.
pub type GammaMonomial = Vec<(i64, i64)>;

/// `Π PMaster(name)^k`: `(name, k)` with `k > 0`, sorted by name.
pub type MasterMonomial = Vec<(String, i64)>;

/// A reduced value in the tag basis: the coefficient of every monomial,
/// normalized by `together()`, zero coefficients dropped.
pub type TagBasis = BTreeMap<(GammaMonomial, MasterMonomial), Atom>;

fn failure(detail: impl std::fmt::Display) -> VakintError {
    VakintError::EvaluationError(format!("p-integral tag basis: {detail}"))
}

/// `(n, c)` of a Γ argument `n + c·ε` with integers `n` and `c`.
fn integer_linear_argument(argument: AtomView, epsilon: Symbol) -> Option<(i64, i64)> {
    let argument = argument.to_owned().expand();
    let e = Atom::var(epsilon);
    let c = get_integer_from_atom(argument.coefficient(e.as_view()).expand().as_view())?;
    let n = get_integer_from_atom(
        argument
            .replace(e.to_pattern())
            .with(Atom::Zero.to_pattern())
            .expand()
            .as_view(),
    )?;
    (Atom::num(n) + Atom::num(c) * &e - &argument)
        .expand()
        .is_zero()
        .then_some((n, c))
}

/// Every `Γ(n + cε)` written as a rational function of `ε` times
/// `Γ(1+cε)` (`Γ(n)` itself for `c = 0`); any other Γ is an error.
fn normalize_gammas(expression: &Atom, epsilon: Symbol) -> Result<Atom, VakintError> {
    let mut unsupported = None;
    expression.visitor(&mut |part| {
        if let AtomView::Fun(function) = part
            && function.get_symbol() == gamma()
        {
            let argument = function.iter().next();
            if function.get_nargs() != 1
                || !integer_linear_argument(argument.unwrap(), epsilon)
                    .is_some_and(|(n, c)| c != 0 || n > 0)
            {
                unsupported = Some(part.to_owned());
            }
            return false;
        }
        unsupported.is_none()
    });
    if let Some(part) = unsupported {
        return Err(failure(format!("{part} is not a finite Γ(n + cε)")));
    }
    Ok(expression
        .replace(gamma().call(Atom::var(S.x_)).to_pattern())
        .with_map(move |matched| {
            let argument = matched.get(S.x_).unwrap().to_atom();
            let (n, c) =
                integer_linear_argument(argument.as_view(), epsilon).expect("validated Γ argument");
            let e = Atom::var(epsilon);
            let z = Atom::num(c) * &e;
            let mut value = if c == 0 {
                Atom::num(1)
            } else {
                gamma_one_plus(c, &e)
            };
            for m in 1..n {
                value *= Atom::num(m) + &z;
            }
            for m in n..=0 {
                value /= Atom::num(m) + &z;
            }
            value
        }))
}

/// `m` of `Γ(1+mε)`, if `atom` is one with `m ≠ 0`.
fn gamma_tag(atom: AtomView, epsilon: Symbol) -> Option<i64> {
    let AtomView::Fun(function) = atom else {
        return None;
    };
    if function.get_symbol() != gamma() || function.get_nargs() != 1 {
        return None;
    }
    match integer_linear_argument(function.iter().next().unwrap(), epsilon)? {
        (1, m) if m != 0 => Some(m),
        _ => None,
    }
}

/// Every `PG(x1,x2)` has two integer arguments, every `PGamma(m)` one nonzero
/// integer argument and every `PMaster(name)` a tabulated name in vakint's
/// namespace (`masters::p_master_names`); anything else is an error.
pub(crate) fn check_symbols(expression: &Atom) -> Result<(), VakintError> {
    let known = masters::p_master_names()?;
    let mut malformed = None;
    expression.visitor(&mut |part| {
        if let AtomView::Fun(function) = part {
            let integers = function
                .iter()
                .map(get_integer_from_atom)
                .collect::<Option<Vec<_>>>();
            let valid = if function.get_symbol() == pg_symbol() {
                integers.is_some_and(|x| x.len() == 2)
            } else if function.get_symbol() == p_gamma_symbol() {
                integers.is_some_and(|m| m.len() == 1 && m[0] != 0)
            } else if function.get_symbol() == p_master_symbol() {
                masters::p_master_name(part).is_some_and(|name| known.contains(&name.as_str()))
            } else {
                return malformed.is_none();
            };
            if !valid {
                malformed = Some(part.to_owned());
            }
            return false;
        }
        malformed.is_none()
    });
    match malformed {
        Some(part) => Err(VakintError::EvaluationError(format!(
            "p-integral: malformed {part}"
        ))),
        None => Ok(()),
    }
}

/// `reduced` in the tag basis: a Laurent polynomial in the `PG` and `PGamma`
/// and a polynomial in the `PMaster`, with coefficients rational in `epsilon`
/// and in whatever other symbols they hold. Everything else is an error:
/// unknown masters, malformed `PG`/`PGamma`/`PMaster` ([`check_symbols`]), Γ
/// functions other than a finite `Γ(n + cε)`, powers that are not `i16`
/// integers, negative master powers, and coefficients that still contain a Γ
/// or a master (`1/(1+PGamma(1))`, `foo(PMaster(no))`).
pub fn tag_basis(reduced: &Atom, epsilon: Symbol) -> Result<TagBasis, VakintError> {
    check_symbols(reduced)?;
    // Before the rational factors of a huge `PG` power get expanded.
    let tags = [pg_symbol(), p_gamma_symbol(), p_master_symbol()];
    if let Some(power) = unsupported_power(reduced, &tags) {
        return Err(failure(format!("{power} is not an i16 integer power")));
    }
    let substituted = reduced
        .replace(
            pg_symbol()
                .call((Atom::var(S.x_), Atom::var(S.y_)))
                .to_pattern(),
        )
        .with_map(move |matched| {
            let x = get_integer_from_atom(matched.get(S.x_).unwrap().to_atom().as_view());
            let y = get_integer_from_atom(matched.get(S.y_).unwrap().to_atom().as_view());
            pg_value(x.expect("validated PG"), y.expect("validated PG"), epsilon)
        })
        .replace(p_gamma_symbol().call(Atom::var(S.x_)).to_pattern())
        .with_map(move |matched| {
            let m = get_integer_from_atom(matched.get(S.x_).unwrap().to_atom().as_view());
            gamma_one_plus(m.expect("validated PGamma"), &Atom::var(epsilon))
        });
    let normalized = normalize_gammas(
        &crate::Vakint::normalize_gamma_functions(substituted.as_view()),
        epsilon,
    )?;
    if !normalized.is_finite() {
        return Err(failure(format!("{reduced} is not finite")));
    }
    let mut variables = BTreeSet::new();
    let mut unknown = None;
    normalized.visitor(&mut |part| {
        if let AtomView::Fun(function) = part {
            if function.get_symbol() == gamma() {
                if gamma_tag(part, epsilon).is_some() {
                    variables.insert(part.to_owned());
                } else {
                    unknown = Some(format!("unexpected {part}"));
                }
                return false;
            }
            if function.get_symbol() == p_master_symbol() {
                variables.insert(part.to_owned());
                return false;
            }
            if [pg_symbol(), p_gamma_symbol()].contains(&function.get_symbol()) {
                unknown = Some(format!("malformed {part}"));
                return false;
            }
        }
        true
    });
    if let Some(detail) = unknown {
        return Err(failure(detail));
    }
    let variables = variables.into_iter().collect::<Vec<_>>();
    let expanded = normalized.expand();
    if let Some(power) = unsupported_power(&expanded, &[gamma(), p_master_symbol()]) {
        return Err(failure(format!("{power} is not an i16 integer power")));
    }
    let mut basis = TagBasis::new();
    for (key, coefficient) in expanded.coefficient_list::<i16>(&variables) {
        if [gamma(), p_master_symbol()]
            .iter()
            .any(|symbol| coefficient.contains_symbol(*symbol))
        {
            return Err(failure(format!(
                "the coefficient {coefficient} still contains a Γ or a master"
            )));
        }
        let factors = match key.as_view() {
            AtomView::Mul(product) => product.iter().map(|f| f.to_owned()).collect(),
            _ => vec![key.clone()],
        };
        let mut gammas = GammaMonomial::new();
        let mut masters = MasterMonomial::new();
        for factor in factors {
            let (base, exponent) = match factor.as_view() {
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    (base.to_owned(), get_integer_from_atom(exponent))
                }
                AtomView::Num(_) if factor == Atom::num(1) => continue,
                _ => (factor.clone(), Some(1)),
            };
            let exponent = exponent.ok_or_else(|| failure(format!("monomial factor {factor}")))?;
            if let Some(m) = gamma_tag(base.as_view(), epsilon) {
                gammas.push((m, exponent));
            } else if let Some(name) = masters::p_master_name(base.as_view()) {
                if exponent < 0 {
                    return Err(failure(format!("negative master power {factor}")));
                }
                masters.push((name, exponent));
            } else {
                return Err(failure(format!("monomial factor {factor}")));
            }
        }
        gammas.sort();
        masters.sort();
        let entry = basis.entry((gammas, masters)).or_default();
        *entry = (&*entry + &coefficient).together();
    }
    basis.retain(|_, coefficient| !coefficient.is_zero());
    Ok(basis)
}

/// The difference `a - b` in the tag basis: empty iff `a` and `b` agree
/// exactly.
pub fn tag_basis_difference(a: &Atom, b: &Atom, epsilon: Symbol) -> Result<TagBasis, VakintError> {
    tag_basis(&(a - b), epsilon)
}

#[cfg(test)]
mod tests {
    use symbolica::atom::{Atom, AtomCore};
    use symbolica::function;

    use super::*;
    use crate::Vakint;
    use crate::utils::vakint_macros::{vk_parse, vk_symbol};

    fn epsilon() -> Symbol {
        Vakint::initialize_vakint_symbols();
        vk_symbol!("ε")
    }

    fn pg(x1: i64, x2: i64) -> Atom {
        pg_symbol().call((x1, x2))
    }

    fn master(name: &str) -> Atom {
        function!(p_master_symbol(), Atom::var(vk_symbol!(name)))
    }

    fn tag(m: i64) -> Atom {
        p_gamma_symbol().call((m,))
    }

    #[test]
    fn forcer_g_scheme_relations_are_identities_in_the_tag_basis() {
        let e = epsilon();
        let ep = Atom::var(e);
        let one = Atom::num(1);
        // Masters.prc:472-473: G(2,2) = G(3,1) G(2,0) InvG(1,0) (1-4ε)/(1-3ε).
        let relation = pg(2, 2)
            - pg(3, 1) * pg(2, 0) / pg(1, 0) * (&one - Atom::num(4) * &ep)
                / (&one - Atom::num(3) * &ep);
        assert_eq!(tag_basis(&relation, e).unwrap(), TagBasis::new());
        // Masters.prc:474-475: G(3,1) G(0,0) = G(4,0) G(3,0) (1-5ε)/(1-2ε).
        let relation = pg(3, 1) * pg(0, 0)
            - pg(4, 0) * pg(3, 0) * (&one - Atom::num(5) * &ep) / (&one - Atom::num(2) * &ep);
        assert_eq!(tag_basis(&relation, e).unwrap(), TagBasis::new());
        // Masters.prc: GschemeConstants(0,0) = rat(1,ep).
        assert_eq!(
            tag_basis(&(pg(0, 0) - &one / &ep), e).unwrap(),
            TagBasis::new()
        );
        // Without the rational factor the relation fails.
        let wrong = pg(2, 2) - pg(3, 1) * pg(2, 0) / pg(1, 0);
        assert!(!tag_basis(&wrong, e).unwrap().is_empty());
    }

    #[test]
    fn tags_and_masters_are_the_basis_monomials() {
        let e = epsilon();
        let ep = Atom::var(e);
        let mass = vk_parse!("mass").unwrap();
        let reduced = Atom::num(3) / &ep * tag(2).pow(Atom::num(-1)) * master("no")
            + tag(-1).pow(Atom::num(2)) * tag(2) * master("t105").pow(Atom::num(2)) * &mass
            + Atom::num(5) * master("no") / tag(2)
            + Atom::num(7);
        let basis = tag_basis(&reduced, e).unwrap();
        let expected = TagBasis::from([
            ((vec![], vec![]), Atom::num(7)),
            (
                (vec![(2, -1)], vec![("no".to_string(), 1)]),
                (Atom::num(3) / &ep + Atom::num(5)).together(),
            ),
            (
                (vec![(-1, 2), (2, 1)], vec![("t105".to_string(), 2)]),
                mass.clone(),
            ),
        ]);
        assert_eq!(basis.len(), expected.len(), "{basis:?}");
        for (key, value) in &expected {
            let found = basis
                .get(key)
                .unwrap_or_else(|| panic!("{key:?} in {basis:?}"));
            assert!(
                (found - value).expand().is_zero(),
                "{key:?}: {found} != {value}"
            );
        }
        // A raw Γ(n+mε) is the same monomial as its PGamma.
        let raw = vk_parse!("Gamma(3+2*ε)").unwrap();
        let tagged =
            tag(2) * (Atom::num(1) + Atom::num(2) * &ep) * (Atom::num(2) + Atom::num(2) * &ep);
        assert!(tag_basis_difference(&raw, &tagged, e).unwrap().is_empty());
    }

    #[test]
    fn malformed_and_unknown_symbols_are_refused() {
        let e = epsilon();
        for bad in [
            master("nothere"),
            function!(
                p_master_symbol(),
                Atom::var(symbolica::symbol!("other::no"))
            ),
            function!(p_master_symbol(), Atom::var(vk_symbol!("no")), 2),
            p_gamma_symbol().call((0,)),
            pg_symbol().call((1,)),
            vk_parse!("Gamma(1/2+ε)").unwrap(),
        ] {
            assert!(tag_basis(&bad, e).is_err(), "{bad}");
        }
    }

    #[test]
    fn what_the_basis_cannot_hold_is_refused() {
        let e = epsilon();
        let ep = Atom::var(e);
        let one = Atom::num(1);
        let half = Atom::num(1) / Atom::num(2);
        let foo = vk_symbol!("foo");
        for bad in [
            // Γ poles; Symbolica evaluates them to complex infinity.
            vk_parse!("Gamma(0)").unwrap(),
            vk_parse!("Gamma(0)+1").unwrap(),
            &ep * vk_parse!("Gamma(-1)").unwrap(),
            // Not a monomial in the tags and masters.
            &one / (&one + tag(1)),
            &one / (&one + master("no")),
            function!(foo, tag(1)),
            function!(foo, master("no")),
            function!(foo, pg(1, 0)),
            tag(1).pow(&half),
            master("no").pow(&half),
            tag(1).pow(ep.clone()),
            // Negative master powers.
            master("no").pow(Atom::num(-1)),
            tag(1) * master("no").pow(Atom::num(-2)),
            // Beyond the i16 exponents of the tags.
            tag(1).pow(Atom::num(32768)),
            tag(1).pow(Atom::num(-32769)),
            tag(1).pow(Atom::num(40000)),
            (tag(1).pow(Atom::num(30000)) + &one).pow(Atom::num(2)),
            pg(1, 0).pow(Atom::num(40000)),
        ] {
            let error = tag_basis(&bad, e).expect_err(&bad.to_string());
            assert!(
                matches!(error, VakintError::EvaluationError(_)),
                "{bad}: {error}"
            );
        }
        // Γ(n + cε) with its ε-poles is fine, and so are the ends of the range.
        assert!(tag_basis(&vk_parse!("Gamma(-1+ε)").unwrap(), e).is_ok());
        for power in [i64::from(i16::MAX), i64::from(i16::MIN)] {
            let basis = tag_basis(&tag(1).pow(Atom::num(power)), e).unwrap();
            assert!(basis.contains_key(&(vec![(1, power)], vec![])), "{basis:?}");
        }
    }
}
