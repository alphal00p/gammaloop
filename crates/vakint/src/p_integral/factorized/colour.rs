//! The colour split: every local factor (edge or vertex numerator in source
//! notation) is written as `Σ_rows colour_row × residual_row`, the colour
//! rows of all locals are multiplied out into branches, and the colour
//! network of each branch is closed on its external legs and reduced to one
//! coefficient.
//!
//! This is the same decomposition as the beta crate's `color.rs`
//! (`split_local` → Cartesian product of rows → `ColorBranch { coefficient,
//! diagram }`, coefficient `δ_ab C_ab / N_A`), done on SOURCE atoms (the
//! uncooked DOT notation) instead of cooked tensors, so that rows are keyed by
//! portable canonical strings and the result does not depend on process-local
//! cooking hashes. It is written to be reused by the reference chain (lane
//! refpath): one implementation of the split for the helper.
//!
//! Differences to the beta crate (stated, not hidden): no branch merging
//! (`merge_channels`): one value per non-zero channel; rows are normalised by
//! taking each colour product with unit numeric coefficient (the sign stays in
//! the residual), not by the beta's residual-key sign rule.

use std::collections::BTreeMap;

use idenso::{
    CookMode, CookSettings,
    color::ColorSimplifySettings,
    tensor::{AlgebraContraction, AlgebraSettings, SymbolicTensor},
};
use spenso::network::library::symbolic::ETS;
use symbolica::atom::{Atom, AtomCore, AtomView};

use super::{FactorizedError, Result, as_slot, contains_colour_slot, is_colour_rep, names};

/// How the colour coefficient of a branch is reduced.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ColourMode {
    /// Group invariants with symbolic dimensions (idenso's closed-network
    /// reducer, color.h identities): `coad(8, ..)` / `cof(3, ..)` are read as
    /// `coad(dA, ..)` / `cof(Nc, ..)` (`spenso::dA`, `spenso::Nc`), the result
    /// is a polynomial in `cas(2, ..)`, `idx(2, ..)`, `dA`, `Nc` and named
    /// higher invariants, divided by `dA`. Symbolic SU(N) / general group.
    #[default]
    Invariants,
    /// SU(3) numbers, as the beta crate's production `contract_color`
    /// (idenso `simplify_algebra` with colour, `with_cof_dimension_invariants`,
    /// dot contraction), divided by `N_A = 8`.
    Su3,
}

/// One colour row of a local factor.
#[derive(Clone, Debug, PartialEq)]
pub struct ColourRow {
    /// Product of the colour tensors (unit numeric coefficient); `1` when the
    /// local has no colour.
    pub colour: Atom,
    /// The colour-free rest: `local = Σ_rows colour × residual`.
    pub residual: Atom,
}

/// Is `view` one colour tensor: a function with a colour slot argument?
fn is_colour_tensor(view: AtomView<'_>) -> bool {
    let AtomView::Fun(f) = view else {
        return false;
    };
    let symbol = f.get_symbol();
    let name = symbol.get_name();
    if name == names::STRUCTURE_F || name == names::GENERATOR_T {
        return true;
    }
    f.iter()
        .filter_map(as_slot)
        .any(|slot| slot.rep.starts_with("dind(") || is_colour_rep(&slot.rep))
}

fn product(factors: impl IntoIterator<Item = Atom>) -> Atom {
    factors.into_iter().fold(Atom::num(1), |acc, f| acc * f)
}

/// Split one local factor into colour rows, ordered by the canonical string of
/// the colour part. A term whose colour sits inside a nested sum or power is
/// expanded once; if colour and the rest still do not separate, an error.
pub fn split_local_colour(local: &Atom) -> Result<Vec<ColourRow>> {
    let mut rows = BTreeMap::<String, (Atom, Vec<Atom>)>::new();
    let mut add_term =
        |term: AtomView<'_>, expanded: bool| -> std::result::Result<bool, FactorizedError> {
            let factors = match term {
                AtomView::Mul(m) => m.iter().collect::<Vec<_>>(),
                other => vec![other],
            };
            let mut colour = Vec::new();
            let mut rest = Vec::new();
            for factor in factors {
                if is_colour_tensor(factor) {
                    colour.push(factor.to_owned());
                } else if contains_colour_slot(factor) {
                    if expanded {
                        return Err(FactorizedError::Colour {
                            detail: format!(
                                "colour does not factor out of {}",
                                term.to_owned().to_canonical_string()
                            ),
                        });
                    }
                    return Ok(false);
                } else {
                    rest.push(factor.to_owned());
                }
            }
            let colour = product(colour);
            let key = colour.to_canonical_string();
            rows.entry(key)
                .or_insert_with(|| (colour, Vec::new()))
                .1
                .push(product(rest));
            Ok(true)
        };
    let terms = match local.as_view() {
        AtomView::Add(a) => a.iter().map(|t| t.to_owned()).collect::<Vec<_>>(),
        _ => vec![local.clone()],
    };
    for term in terms {
        if !add_term(term.as_view(), false)? {
            let expanded = term.expand();
            let parts = match expanded.as_view() {
                AtomView::Add(a) => a.iter().map(|t| t.to_owned()).collect::<Vec<_>>(),
                _ => vec![expanded.clone()],
            };
            for part in parts {
                add_term(part.as_view(), true)?;
            }
        }
    }
    Ok(rows
        .into_values()
        .map(|(colour, residuals)| ColourRow {
            colour,
            residual: Atom::add_many(residuals.iter()),
        })
        .collect())
}

/// The colour slots of a product, with multiplicities (by canonical string).
fn colour_slots(colour: &Atom) -> BTreeMap<String, (usize, Atom)> {
    let mut slots = BTreeMap::<String, (usize, Atom)>::new();
    let mut add = |view: AtomView<'_>| {
        if let Some(slot) = as_slot(view)
            && (is_colour_rep(&slot.rep) || slot.rep.starts_with("dind("))
        {
            // dind(cof(N, i)) and cof(N, i) are the two ends of one index.
            let base = match view {
                AtomView::Fun(f) if f.get_symbol().get_name() == names::DIND => {
                    f.iter().next().unwrap().to_owned()
                }
                _ => view.to_owned(),
            };
            let entry = slots.entry(base.to_canonical_string()).or_insert((0, base));
            entry.0 += 1;
        }
    };
    let factors = match colour.as_view() {
        AtomView::Mul(m) => m.iter().collect::<Vec<_>>(),
        other => vec![other],
    };
    for factor in factors {
        let (base, times) = match factor {
            AtomView::Pow(p) => {
                let (base, exponent) = p.get_base_exp();
                (base, i64::try_from(exponent).unwrap_or(1).max(1) as usize)
            }
            other => (other, 1),
        };
        if let AtomView::Fun(f) = base {
            for _ in 0..times {
                for arg in f.iter() {
                    add(arg);
                }
            }
        }
    }
    slots
}

/// `δ_ab C_ab / N_A` of a branch's colour product `C` with two open adjoint
/// slots (the external legs); `C` itself when it has no open slot. Other open
/// colour structures are not supported in v1.
pub fn colour_coefficient(colour: &Atom, mode: ColourMode) -> Result<Atom> {
    let err = |detail: String| FactorizedError::Colour { detail };
    if colour.is_zero() {
        return Ok(Atom::Zero);
    }
    let open = colour_slots(colour)
        .into_values()
        .filter(|(count, _)| *count == 1)
        .map(|(_, slot)| slot)
        .collect::<Vec<_>>();
    if open.is_empty() && !contains_colour_slot(colour.as_view()) {
        return Ok(colour.clone());
    }
    let (closed, dimension) = match open.as_slice() {
        [a, b] => {
            let (sa, sb) = (as_slot(a.as_view()).unwrap(), as_slot(b.as_view()).unwrap());
            if sa.rep != names::COAD || sb.rep != names::COAD || sa.dim != sb.dim {
                return Err(err(format!(
                    "open colour slots {} and {} are not two adjoint legs",
                    a.to_canonical_string(),
                    b.to_canonical_string()
                )));
            }
            (colour * ETS.metric(a.clone(), b.clone()), sa.dim)
        }
        [] => (colour.clone(), Atom::num(1)),
        _ => {
            return Err(err(format!(
                "{} open colour slots; v1 closes exactly two adjoint legs",
                open.len()
            )));
        }
    };
    let coefficient = match mode {
        ColourMode::Invariants => reduce_invariants(&closed, &dimension)?,
        ColourMode::Su3 => reduce_su3(&closed)? / &dimension,
    };
    let coefficient = coefficient.expand();
    if contains_colour_slot(coefficient.as_view()) {
        return Err(err(format!(
            "the reduction left an open colour factor in {}",
            coefficient.to_canonical_string()
        )));
    }
    Ok(coefficient)
}

/// Symbolic dimensions `dA`, `Nc` for the numeric `coad(8)` / `cof(3)`, then
/// idenso's closed-network reducer; divided by the (symbolic) `N_A` when the
/// network was closed on two adjoint legs.
fn reduce_invariants(closed: &Atom, dimension: &Atom) -> Result<Atom> {
    let adjoint = Atom::var(symbolica::symbol!("spenso::dA"));
    let fundamental = Atom::var(idenso::color::CS.nc);
    let symbolic = map_colour_dimensions(closed.as_view(), &adjoint, &fundamental);
    let value = idenso::color::network::try_reduce_closed_color_network(symbolic.as_view())
        .map_err(|detail| FactorizedError::Colour {
            detail: format!("closed-network reducer failed: {detail}"),
        })?
        .ok_or_else(|| FactorizedError::Colour {
            detail: format!(
                "not a closed colour network the reducer supports: {}",
                symbolic.to_canonical_string()
            ),
        })?;
    if dimension.is_one() {
        Ok(value)
    } else {
        Ok(value / adjoint)
    }
}

fn map_colour_dimensions(view: AtomView<'_>, adjoint: &Atom, fundamental: &Atom) -> Atom {
    view.replace_map(|part, _ctx, out| {
        let AtomView::Fun(f) = part else {
            return;
        };
        let symbol = f.get_symbol();
        let name = symbol.get_name();
        let dim = match name {
            names::COAD => adjoint,
            names::COF => fundamental,
            _ => return,
        };
        let args = f.iter().collect::<Vec<_>>();
        if args.len() != 2 {
            return;
        }
        let mut builder = symbolica::atom::FunctionBuilder::new(f.get_symbol());
        builder = builder.add_arg(dim).add_arg(args[1]);
        **out = builder.finish();
    })
}

/// The beta crate's SU(3) production reduction of a closed colour network.
fn reduce_su3(closed: &Atom) -> Result<Atom> {
    let err = |detail: String| FactorizedError::Colour { detail };
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let cooked = cooking
        .try_cook(closed.as_view())
        .map_err(|e| err(format!("cooking: {e}")))?;
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
        contract: AlgebraContraction::Dots,
        ..Default::default()
    };
    let mut reduced = SymbolicTensor::infer_indexed(cooked)
        .and_then(|t| t.undo_trace())
        .and_then(|t| t.undo_chain())
        .and_then(|t| t.simplify_algebra(&settings))
        .map_err(|e| err(format!("colour algebra: {e}")))?;
    for _ in 0..3 {
        let value = cooking.uncook(reduced.expression().as_view());
        if !contains_colour_slot(value.as_view()) {
            break;
        }
        reduced = SymbolicTensor::infer_indexed(reduced.expression().clone())
            .and_then(|t| t.simplify_algebra(&settings))
            .map_err(|e| err(format!("colour algebra: {e}")))?;
    }
    Ok(cooking.uncook(reduced.expression().as_view()))
}
