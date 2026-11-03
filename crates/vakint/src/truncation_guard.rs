//! Fail-closed truncation of substituted master integrals (the `Oep` guard).
//!
//! Every backend substitutes a master integral by the Laurent coefficients it
//! knows plus an explicit unknown tail `ε^k·Oep(k, master)`, where `ε^k` is the
//! first order that is not substituted. Reduction coefficients can carry
//! (spurious) `1/ε` poles even when the dimension stays parametric in the IBPs,
//! so a coefficient `c/ε^p` multiplying such a master moves the unknown tail
//! down to `ε^(k-p)`. After the result is expanded and truncated at the
//! requested order, a surviving `Oep` therefore means that the result depends
//! on master coefficients that are not available. Such a result is rejected
//! with [`MissingMasterOrders`] instead of being returned with a silently
//! missing contribution. A tail beyond the requested order is truncated
//! exactly, so the guard never changes a correct result.
//!
//! The same marker has always been used by the MATAD and FMFT master tables
//! (`ep^3*Oep(3,D6)`, `ep^2*Oep(2,I_PR6)`, ...); this module makes the final
//! check one shared, backend-independent function.

use std::collections::BTreeMap;
use std::fmt;

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::printer::{AtomPrinter, PrintOptions};

use crate::get_integer_from_atom;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};

/// One master whose unknown Laurent tail survived truncation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MissingMasterOrder {
    /// The master integral, in the naming of the backend that substituted it.
    pub master: String,
    /// The first Laurent order `ε^k` of the master that was not available to
    /// the substitution (its `Oep(k, master)` marker).
    pub first_missing_order: i64,
    /// The highest Laurent order of the master that the requested result
    /// needs, when the requested order is known.
    pub last_needed_order: Option<i64>,
    /// The lowest order of the result that the missing coefficients affect.
    pub first_affected_result_order: i64,
}

/// The result needs master Laurent coefficients beyond those available.
///
/// This is raised by every evaluation backend instead of returning a result
/// that silently treats unknown master coefficients as zero.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MissingMasterOrders {
    /// The backend (and data snapshot, when relevant) that substituted the masters.
    pub backend: String,
    /// The highest order `ε^n` requested from this evaluation, when known.
    pub requested_order: Option<i64>,
    /// The masters with missing Laurent coefficients, sorted by name and order.
    pub masters: Vec<MissingMasterOrder>,
}

impl fmt::Display for MissingMasterOrders {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.requested_order {
            Some(order) => write!(f, "the result through 𝒪(ε^{order})")?,
            None => write!(f, "the requested result")?,
        }
        write!(
            f,
            " needs master Laurent coefficients beyond the available expansion depth \
             ({}; an unknown 𝒪(ε^k) master tail survived truncation, e.g. because a \
             reduction coefficient with an ε pole multiplies a master known only to a \
             lower order):",
            self.backend
        )?;
        for missing in &self.masters {
            write!(f, "\n  - master {}: missing ", missing.master)?;
            match missing.last_needed_order {
                Some(last) if last > missing.first_missing_order => write!(
                    f,
                    "Laurent coefficients ε^{}..ε^{last}",
                    missing.first_missing_order
                )?,
                _ => write!(f, "Laurent coefficient ε^{}", missing.first_missing_order)?,
            }
            write!(
                f,
                " (affects the result from ε^{})",
                missing.first_affected_result_order
            )?;
        }
        Ok(())
    }
}

impl std::error::Error for MissingMasterOrders {}

/// The function symbol of an unknown master tail, `Oep(k, master)`.
pub(crate) fn unknown_tail_symbol() -> Symbol {
    vk_symbol!("Oep")
}

/// The marker `Oep(first_missing_order, master)` without its epsilon power.
pub(crate) fn unknown_tail_marker(first_missing_order: i64, master: &Atom) -> Atom {
    unknown_tail_symbol().call((first_missing_order, master))
}

/// The unknown tail `ε^k·Oep(k, master)` of a master known through `ε^(k-1)`.
pub(crate) fn unknown_tail(epsilon: Symbol, first_missing_order: i64, master: &Atom) -> Atom {
    Atom::var(epsilon).pow(Atom::num(first_missing_order))
        * unknown_tail_marker(first_missing_order, master)
}

/// Whether any unknown master tail is present.
pub(crate) fn has_unknown_tail(expression: AtomView) -> bool {
    expression.contains_symbol(unknown_tail_symbol())
}

/// The canonical name of a marker's master without Vakint's namespace.
pub(crate) fn plain_name(master: &Atom) -> String {
    AtomPrinter::new_with_options(master.as_view(), PrintOptions::file_no_namespace()).to_string()
}

/// alphaLoop's `integrateduv.frm` masters end in an `alarm…` symbol times
/// `ε^k` (in the normalization of its `Masters()` substitution).
const ALPHALOOP_ALARMS: [(&str, &str, i64); 8] = [
    ("alarmt1m3", "uvid(1,1)", 4),
    ("alarmt2m21", "uvid(2,1)", 4),
    ("alarmt2m22", "uvid(2,2)", 2),
    ("alarmt3m111111", "uvid(3,1,1,1,1,1,1,1)", 1),
    ("alarmt3m011111", "uvid(3,1,0,1,1,1,1,1)", 1),
    ("alarmt3m011101", "uvid(3,1,0,1,1,1,0,1)", 2),
    ("alarmt3m001111", "uvid(3,1,0,0,1,1,1,1)", 1),
    ("alarmt3m001011", "uvid(3,1,0,0,1,0,1,1)", 2),
];

/// Replace alphaLoop's `alarm…` tail symbols by `Oep(k, alphaLoop master)`.
///
/// FORM has already multiplied each alarm symbol by its `ε^k`, so only the
/// marker is substituted. An unexpected `alarm…` symbol still becomes a
/// marker and fails closed.
pub(crate) fn alphaloop_alarms_to_unknown_tails(expression: &Atom) -> Atom {
    let mut alarms = Vec::new();
    expression.visitor(&mut |part| {
        if let AtomView::Var(variable) = part
            && variable
                .get_symbol()
                .get_stripped_name()
                .starts_with("alarm")
        {
            alarms.push(variable.get_symbol());
        }
        true
    });
    alarms.sort();
    alarms.dedup();
    let replacements: Vec<_> = alarms
        .into_iter()
        .map(|alarm| {
            let name = alarm.get_stripped_name();
            let (master, order) = ALPHALOOP_ALARMS
                .iter()
                .find(|(known, _, _)| *known == name)
                .map(|(_, master, order)| {
                    let master = vk_parse!(master).expect("alphaLoop master names parse");
                    (master, *order)
                })
                .unwrap_or_else(|| (Atom::var(alarm), 0));
            symbolica::id::Replacement::new(
                Atom::var(alarm).to_pattern(),
                unknown_tail_marker(order, &master),
            )
        })
        .collect();
    if replacements.is_empty() {
        expression.clone()
    } else {
        expression.replace_multiple(&replacements)
    }
}

/// Reject `expression` if an unknown master tail survived truncation.
///
/// `expression` must already be truncated at `requested_order` in `epsilon`.
/// The success path is a single symbol search; the diagnostic analysis only
/// runs when a tail survived. `describe` names a master in the backend's
/// terms (the default is [`plain_name`]).
pub(crate) fn ensure_known_orders(
    expression: &Atom,
    epsilon: Symbol,
    requested_order: Option<i64>,
    backend: &str,
    describe: impl Fn(&Atom) -> String,
) -> Result<(), MissingMasterOrders> {
    if !has_unknown_tail(expression.as_view()) {
        return Ok(());
    }
    Err(missing_orders(
        expression,
        epsilon,
        requested_order,
        backend,
        describe,
    ))
}

/// Analyse the surviving tails of an expression that contains one.
fn missing_orders(
    expression: &Atom,
    epsilon: Symbol,
    requested_order: Option<i64>,
    backend: &str,
    describe: impl Fn(&Atom) -> String,
) -> MissingMasterOrders {
    // (master, first missing order) -> lowest affected result order
    let mut tails = BTreeMap::<(String, i64), i64>::new();
    let expanded = expression.expand();
    let terms: Vec<AtomView> = match expanded.as_view() {
        AtomView::Add(sum) => sum.iter().collect(),
        term => vec![term],
    };
    for term in terms {
        let power = epsilon_power(term, epsilon);
        term.visitor(&mut |part| {
            if let AtomView::Fun(function) = part
                && function.get_symbol() == unknown_tail_symbol()
            {
                let arguments: Vec<_> = function.iter().collect();
                let (order, master) = match arguments.as_slice() {
                    [order, master] => (get_integer_from_atom(*order), master.to_owned()),
                    _ => (None, part.to_owned()),
                };
                // A malformed marker still fails closed; it is reported at
                // the order where it appears.
                let order = order.unwrap_or(power);
                tails
                    .entry((describe(&master), order))
                    .and_modify(|lowest| *lowest = (*lowest).min(power))
                    .or_insert(power);
                return false;
            }
            true
        });
    }
    let masters = tails
        .into_iter()
        .map(
            |((master, first_missing_order), first_affected_result_order)| MissingMasterOrder {
                master,
                first_missing_order,
                // A coefficient multiplying ε^k·Oep(k,M) at result order j
                // needs master orders k..k+(n-j) for a result through ε^n.
                last_needed_order: requested_order.map(|requested| {
                    first_missing_order
                        + requested.saturating_sub(first_affected_result_order).max(0)
                }),
                first_affected_result_order,
            },
        )
        .collect();
    MissingMasterOrders {
        backend: backend.to_owned(),
        requested_order,
        masters,
    }
}

/// The integer power of `epsilon` in an expanded product term.
fn epsilon_power(term: AtomView, epsilon: Symbol) -> i64 {
    match term {
        AtomView::Var(variable) if variable.get_symbol() == epsilon => 1,
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            match base {
                AtomView::Var(variable) if variable.get_symbol() == epsilon => {
                    get_integer_from_atom(exponent).unwrap_or(0)
                }
                _ => 0,
            }
        }
        AtomView::Mul(product) => product
            .iter()
            .map(|factor| epsilon_power(factor, epsilon))
            .sum(),
        _ => 0,
    }
}

/// The missing-order report of an evaluation that must fail closed.
#[cfg(test)]
pub(crate) fn expect_missing_orders(error: crate::VakintError) -> MissingMasterOrders {
    match error {
        crate::VakintError::MissingMasterOrders(missing) => missing,
        other => panic!("expected missing master Laurent orders, got: {other}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ep() -> Symbol {
        vk_symbol!("ep")
    }

    /// `coefficient·(known + ε^k·Oep(k,M))`, truncated through ε^requested.
    fn truncated(coefficient: &str, known: &str, first_missing: i64, requested: i64) -> Atom {
        let master = vk_parse!("M").unwrap();
        let value = vk_parse!(known).unwrap() + unknown_tail(ep(), first_missing, &master);
        (vk_parse!(coefficient).unwrap() * value)
            .series(
                ep(),
                Atom::Zero.as_view(),
                symbolica::domains::rational::Rational::from(requested),
            )
            .unwrap()
            .to_atom()
    }

    #[test]
    fn a_pole_exposing_an_unknown_order_is_rejected_with_master_and_orders() {
        // M = 1/ep + 2 + O(ep): the 1/ep^2 coefficient needs M through ep^2.
        let result = truncated("1/ep^2", "1/ep+2", 1, 0);
        let error = ensure_known_orders(&result, ep(), Some(0), "test", plain_name).unwrap_err();
        assert_eq!(
            error.masters,
            [MissingMasterOrder {
                master: "M".into(),
                first_missing_order: 1,
                last_needed_order: Some(2),
                first_affected_result_order: -1,
            }]
        );
        let message = error.to_string();
        for part in ["𝒪(ε^0)", "master M", "ε^1..ε^2", "from ε^-1", "test"] {
            assert!(message.contains(part), "{part}: {message}");
        }
    }

    #[test]
    fn a_tail_beyond_the_requested_order_is_truncated_exactly() {
        // Enough depth: the same pole with M known through ep^2.
        let result = truncated("1/ep^2", "1/ep+2+3*ep+4*ep^2", 3, 0);
        assert!(ensure_known_orders(&result, ep(), Some(0), "test", plain_name).is_ok());
        let expected = vk_parse!("ep^-3+2*ep^-2+3*ep^-1+4").unwrap();
        assert!((result - expected).expand().is_zero());
    }
}
