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
//! Differences to the beta crate (stated, not hidden): branch merging
//! ([`merge_colour_channels`], the beta's `merge_channels`) only as an option
//! of `from_feynkit` (default: one value per non-zero channel); rows are normalised by
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
use symbolica::domains::rational::Rational;

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

/// One term of a merged colour decomposition ([`merge_colour_channels`]):
/// per multi-row local its weights over the local's rows (exact rationals,
/// the first non-zero weight is 1) and the colour coefficient of the term.
#[derive(Clone, Debug, PartialEq)]
pub struct MergedColourTerm {
    pub weights: Vec<Vec<Rational>>,
    pub coefficient: Atom,
}

/// An exact rational number atom (`None`: a symbol, a complex or a float).
pub(crate) fn rational_number(view: AtomView<'_>) -> Option<Rational> {
    use symbolica::coefficient::CoefficientView;
    let AtomView::Num(number) = view else {
        return None;
    };
    match number.get_coeff_view() {
        CoefficientView::Natural(n, d, 0, _) => Some(Rational::from((n, d))),
        CoefficientView::Large(re, im) if im.is_zero() => Some(re.to_rat()),
        _ => None,
    }
}

/// An exact complex rational number atom `(re, im)` (`None`: a float).
fn complex_number(view: AtomView<'_>) -> Option<(Rational, Rational)> {
    use symbolica::coefficient::CoefficientView;
    let AtomView::Num(number) = view else {
        return None;
    };
    match number.get_coeff_view() {
        CoefficientView::Natural(n, d, ni, di) => {
            Some((Rational::from((n, d)), Rational::from((ni, di))))
        }
        CoefficientView::Large(re, im) => Some((re.to_rat(), im.to_rat())),
        _ => None,
    }
}

/// `coefficient = Σ_m c_m M_m`: the colour monomials `M_m` (products of
/// invariants, `1` in SU(3) mode) with exact rational `c_m`, keyed by the
/// canonical string of `M_m`. The imaginary unit is treated as one more
/// factor of a monomial (`(a + b 𝑖) M = a M + b (𝑖 M)`): the weights stay
/// real rationals (a quark loop with three gluons gives imaginary channel
/// coefficients, `Tr(T^a T^b T^c) ∝ d^{abc} + 𝑖 f^{abc}`). `None` when a
/// term's number is not exact.
fn colour_monomials(coefficient: &Atom) -> Option<BTreeMap<String, (Atom, Rational)>> {
    let expanded = coefficient.expand();
    let mut out = BTreeMap::<String, (Atom, Rational)>::new();
    let mut add = |term: AtomView<'_>| -> Option<()> {
        let factors = match term {
            AtomView::Mul(m) => m.iter().collect::<Vec<_>>(),
            other => vec![other],
        };
        let (mut re, mut im) = (Rational::one(), Rational::zero());
        let mut monomial = Atom::num(1);
        for factor in factors {
            match factor {
                AtomView::Num(_) => {
                    let (a, b) = complex_number(factor)?;
                    // (re + im i)(a + b i)
                    let next_re = &(&re * &a) - &(&im * &b);
                    let next_im = &(&re * &b) + &(&im * &a);
                    (re, im) = (next_re, next_im);
                }
                other => monomial *= other.to_owned(),
            }
        }
        for (value, unit) in [(re, Atom::num(1)), (im, Atom::i())] {
            if value.is_zero() {
                continue;
            }
            let monomial = &monomial * unit;
            let entry = out
                .entry(monomial.to_canonical_string())
                .or_insert_with(|| (monomial, Rational::zero()));
            entry.1 += &value;
        }
        Some(())
    };
    match expanded.as_view() {
        AtomView::Add(sum) => {
            for term in sum.iter() {
                add(term)?;
            }
        }
        view if view.is_zero() => {}
        view => add(view)?,
    }
    out.retain(|_, (_, value)| !value.is_zero());
    Some(out)
}

/// Exact rank factorisation of the channel tensor `c(k_1, ..., k_m)` over Q
/// (the beta crate's `rank_factors`, Gauss-Jordan on the unfolding of the
/// first axis, recursively): terms `(prod_v w_v(k_v), coefficient)`.
fn rank_factors(
    dims: &[usize],
    entries: &BTreeMap<Vec<usize>, Rational>,
) -> Vec<(Vec<Vec<Rational>>, Rational)> {
    let Some((&rows, rest)) = dims.split_first() else {
        return entries
            .get(&Vec::new())
            .filter(|value| !value.is_zero())
            .map(|value| vec![(Vec::new(), value.clone())])
            .unwrap_or_default();
    };
    let columns = entries
        .keys()
        .map(|key| key[1..].to_vec())
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let original = (0..rows)
        .map(|row| {
            columns
                .iter()
                .map(|column| {
                    let mut key = vec![row];
                    key.extend_from_slice(column);
                    entries.get(&key).cloned().unwrap_or_else(Rational::zero)
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut reduced = original.clone();
    let mut pivots = Vec::new();
    for column in 0..columns.len() {
        let rank = pivots.len();
        if rank == rows {
            break;
        }
        let Some(pivot) = (rank..rows).find(|row| !reduced[*row][column].is_zero()) else {
            continue;
        };
        reduced.swap(rank, pivot);
        let scale = reduced[rank][column].clone();
        for value in reduced[rank].iter_mut() {
            *value = &*value / &scale;
        }
        for row in 0..rows {
            if row == rank || reduced[row][column].is_zero() {
                continue;
            }
            let factor = reduced[row][column].clone();
            for position in 0..columns.len() {
                let eliminated = &factor * &reduced[rank][position];
                reduced[row][position] = &reduced[row][position] - &eliminated;
            }
        }
        pivots.push(column);
    }
    let mut merged = Vec::new();
    for (index, pivot) in pivots.iter().enumerate() {
        let weights = (0..rows)
            .map(|row| original[row][*pivot].clone())
            .collect::<Vec<_>>();
        let slice = columns
            .iter()
            .zip(&reduced[index])
            .filter(|(_, value)| !value.is_zero())
            .map(|(column, value)| (column.clone(), value.clone()))
            .collect::<BTreeMap<_, _>>();
        for (mut rest_weights, coefficient) in rank_factors(rest, &slice) {
            rest_weights.insert(0, weights.clone());
            merged.push((rest_weights, coefficient));
        }
    }
    merged
}

/// MERGED COLOUR TERMS (helper VH lane share, experiment A1; the beta crate's
/// `merge_channels` / `merge_symbolic_channels`): the non-zero colour channels
/// of a graph, `channels[i] = (row of every multi-row local, coefficient)`
/// with `dims[v]` rows at multi-row local `v`, written as a sum of product
/// terms `coefficient_t × Π_v (Σ_k w_{t,v,k} row_{v,k})`.
///
/// Exact and simple with SYMBOLIC colour: every coefficient is split into its
/// colour monomials `Σ_m c_m M_m` with rational `c_m` (in SU(3) mode one
/// monomial, `1`), the monomials become one more axis of the tensor
/// `c(k_1, .., k_m, m)` of rationals, and that tensor is rank-factorised over
/// Q. The local weights are therefore RATIONAL numbers (each local factor
/// stays a colour-free tensor with rational coefficients, as CONTRACT-v1
/// requires) and the coefficient of a term is `Σ_m w_m M_m`, a polynomial in
/// the invariants. Not chosen: a rank factorisation over the field of
/// rational functions of the invariants (fewer terms possible, but weights
/// would carry colour symbols into the factors); merging only proportional
/// channels (a special case of this one). The identity
/// `Σ_channels c(k) Π_v row_{v,k_v} = Σ_t coefficient_t Π_v Σ_k w_{t,v,k} row_{v,k}`
/// holds exactly (test `merged_colour_terms_reconstruct_every_channel`).
///
/// `None` when a coefficient is not a rational combination of monomials.
/// Each weight vector is normalised to its first non-zero entry.
pub fn merge_colour_channels(
    dims: &[usize],
    channels: &[(Vec<usize>, Atom)],
) -> Option<Vec<MergedColourTerm>> {
    let mut monomials = BTreeMap::<String, Atom>::new();
    let mut split = Vec::with_capacity(channels.len());
    for (choices, coefficient) in channels {
        if choices.len() != dims.len() || choices.iter().zip(dims).any(|(k, n)| k >= n) {
            return None;
        }
        let parts = colour_monomials(coefficient)?;
        for (key, (monomial, _)) in &parts {
            monomials
                .entry(key.clone())
                .or_insert_with(|| monomial.clone());
        }
        split.push((choices, parts));
    }
    let position = monomials
        .keys()
        .enumerate()
        .map(|(index, key)| (key.clone(), index))
        .collect::<BTreeMap<_, _>>();
    let mut entries = BTreeMap::new();
    for (choices, parts) in &split {
        for (key, (_, value)) in parts {
            let mut index = choices.to_vec();
            index.push(position[key]);
            entries.insert(index, value.clone());
        }
    }
    let mut axes = dims.to_vec();
    axes.push(monomials.len());
    let monomials = monomials.into_values().collect::<Vec<_>>();
    let mut merged = Vec::new();
    for (mut weights, mut value) in rank_factors(&axes, &entries) {
        let colour = weights.pop()?;
        for vector in weights.iter_mut() {
            let first = vector.iter().find(|weight| !weight.is_zero())?.clone();
            for weight in vector.iter_mut() {
                *weight = &*weight / &first;
            }
            value = &value * &first;
        }
        let coefficient = Atom::add_many(
            colour
                .iter()
                .zip(&monomials)
                .filter(|(weight, _)| !weight.is_zero())
                .map(|(weight, monomial)| Atom::num(&value * weight) * monomial),
        );
        merged.push(MergedColourTerm {
            weights,
            coefficient: coefficient.expand(),
        });
    }
    Some(merged)
}
