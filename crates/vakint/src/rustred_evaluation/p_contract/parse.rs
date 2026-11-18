//! Factor atoms (CONTRACT-v1 source notation) into kernel tensors.
//!
//! A Lorentz factor is read recursively without expanding it: a sum of
//! momenta on one index (`-Q(2,μ) + Q(3,μ)`, the vertex rules) becomes ONE
//! momentum leg with the combined canonical vector; products contract their
//! repeated indices at once (slot forms for momentum pairs). Momenta
//! `Q(e, μ)` become the canonical vector of edge `e` (the match's linear map
//! applied to the factors, idea (a)).
//!
//! A fermion factor (one with spinor slots `bis(4, ..)`) is split into its
//! gammas `(row, column, Lorentz slot)`, spinor identities and the rest; the
//! gammas of one closed loop are chained into a trace word
//! ([`super::trace::Letter`]): a gamma whose Lorentz index is a dummy of its
//! own factor contracted with one momentum `Q(e, x)` is the slash of edge `e`,
//! every other gamma is an open index.

use ahash::AHashMap as HashMap;
use symbolica::atom::{Atom, AtomView, Symbol};
use symbolica::coefficient::CoefficientView;
use symbolica::domains::rational::Rational;

use super::kernel::{DPoly, End, KernelCounters, Kinematics, Mono, Slots, Tensor, contract, leg};
use super::trace::Letter;
use crate::p_integral::factorized::{as_slot, names};

/// Reads factor atoms against one graph's canonical momenta.
pub(crate) struct Parser {
    pub(crate) dimension: Symbol,
    /// Edge id → interned canonical vector (`None`: the zero vector).
    pub(crate) edge_vectors: HashMap<i64, Option<u16>>,
    /// Interned Lorentz slots (`mink(D, ..)` atoms).
    pub(crate) indices: HashMap<Atom, u32>,
}

/// An exact rational number atom.
pub(crate) fn rational(view: AtomView<'_>) -> Option<Rational> {
    let AtomView::Num(number) = view else {
        return None;
    };
    // Machine-size exact rationals only (FeynGen rules and trace results);
    // anything else is refused by the caller.
    match number.get_coeff_view() {
        CoefficientView::Natural(n, d, 0, _) => Some(Rational::from((n, d))),
        _ => None,
    }
}

impl Parser {
    pub(crate) fn index(&mut self, slot: &Atom) -> u32 {
        let next = self.indices.len() as u32;
        *self.indices.entry(slot.clone()).or_insert(next)
    }

    /// The edge of a bare momentum `Q(e, mink(D))` (idenso's dot notation:
    /// the representation without an index).
    fn bare_momentum<'b>(&self, view: AtomView<'b>) -> Option<AtomView<'b>> {
        let AtomView::Fun(function) = view else {
            return None;
        };
        if function.get_symbol().get_name() != names::MOMENTUM || function.get_nargs() != 2 {
            return None;
        }
        let mut args = function.iter();
        let edge = args.next()?;
        let AtomView::Fun(rep) = args.next()? else {
            return None;
        };
        (rep.get_symbol().get_name() == names::MINK && rep.get_nargs() == 1).then_some(edge)
    }

    fn momentum_vector(&self, edge: AtomView<'_>) -> Result<Option<u16>, String> {
        let id = i64::try_from(edge).map_err(|_| format!("a momentum of edge {edge}"))?;
        self.edge_vectors
            .get(&id)
            .copied()
            .ok_or_else(|| format!("Q({id}, ..) of an unknown edge"))
    }

    /// The tensor of a Lorentz factor (no spinor slot).
    pub(crate) fn parse(
        &mut self,
        view: AtomView<'_>,
        kinematics: &mut Kinematics,
        slots: &Slots,
        counters: &mut KernelCounters,
    ) -> Result<Tensor, String> {
        match view {
            AtomView::Num(_) => {
                let value = rational(view)
                    .ok_or_else(|| format!("a non-rational number {view} in a factor"))?;
                Ok(Tensor::scalar(DPoly::constant(value)))
            }
            AtomView::Var(variable) => {
                if variable.get_symbol() == self.dimension {
                    Ok(Tensor::scalar(DPoly::d_power(1)))
                } else {
                    Err(format!("an unexpected scalar {view} in a factor"))
                }
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                let name = symbol.get_name();
                let args = function.iter().collect::<Vec<_>>();
                match (name, args.as_slice()) {
                    (names::METRIC, [a, b])
                        if self.bare_momentum(*a).is_some() || self.bare_momentum(*b).is_some() =>
                    {
                        // idenso's dot notation: `g(Q(e, mink(D)), Q(f, mink(D)))` is
                        // `Q_e · Q_f`, `g(Q(e, mink(D)), mink(D, μ))` is `Q_e^μ`.
                        let mut ends = Vec::with_capacity(2);
                        for argument in [*a, *b] {
                            match self.bare_momentum(argument) {
                                Some(edge) => match self.momentum_vector(edge)? {
                                    None => return Ok(Tensor::default()),
                                    Some(vector) => ends.push(End::Mom(vector)),
                                },
                                None => {
                                    let Some(slot) =
                                        as_slot(argument).filter(|s| s.rep == names::MINK)
                                    else {
                                        return Err(format!(
                                            "a metric argument {argument} in {view}"
                                        ));
                                    };
                                    ends.push(End::Idx(self.index(&slot.atom)));
                                }
                            }
                        }
                        match (ends[0], ends[1]) {
                            (End::Mom(u), End::Mom(v)) => {
                                let form = kinematics.pair_form(u, v);
                                let mut tensor = Tensor::one();
                                let poly = tensor.terms.remove(&Vec::new()).unwrap_or_default();
                                let poly =
                                    super::kernel::poly_mul_form(&poly, &form, slots, counters)?;
                                tensor.add(Vec::new(), poly);
                                Ok(tensor)
                            }
                            (x, y) => Ok(Tensor::structure(vec![leg(x, y)])),
                        }
                    }
                    (names::METRIC, [a, b]) => {
                        let (Some(sa), Some(sb)) = (as_slot(*a), as_slot(*b)) else {
                            return Err(format!("a metric without two slots: {view}"));
                        };
                        if sa.rep != names::MINK || sb.rep != names::MINK {
                            return Err(format!(
                                "a non-Lorentz metric in a Lorentz factor: {view}"
                            ));
                        }
                        let (ia, ib) = (self.index(&sa.atom), self.index(&sb.atom));
                        if ia == ib {
                            return Ok(Tensor::scalar(DPoly::d_power(1)));
                        }
                        Ok(Tensor::structure(vec![leg(End::Idx(ia), End::Idx(ib))]))
                    }
                    (names::MOMENTUM, [edge, slot]) => {
                        let Some(slot) = as_slot(*slot).filter(|s| s.rep == names::MINK) else {
                            return Err(format!("a momentum without a Lorentz slot: {view}"));
                        };
                        match self.momentum_vector(*edge)? {
                            None => Ok(Tensor::default()),
                            Some(vector) => {
                                let index = self.index(&slot.atom);
                                Ok(Tensor::structure(vec![leg(
                                    End::Mom(vector),
                                    End::Idx(index),
                                )]))
                            }
                        }
                    }
                    _ => Err(format!("an unsupported tensor {view} in a Lorentz factor")),
                }
            }
            AtomView::Add(sum) => {
                let mut total = Tensor::default();
                for term in sum.iter() {
                    total.add_tensor(self.parse(term, kinematics, slots, counters)?);
                }
                merge_momentum_legs(total, kinematics)
            }
            AtomView::Mul(product) => {
                let mut total = Tensor::one();
                for factor in product.iter() {
                    let tensor = self.parse(factor, kinematics, slots, counters)?;
                    total = contract(&total, &tensor, kinematics, slots, counters)?;
                    if total.terms.is_empty() {
                        break;
                    }
                }
                Ok(total)
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let exponent = i64::try_from(exponent)
                    .ok()
                    .filter(|e| *e >= 0)
                    .ok_or_else(|| format!("a non-polynomial power {view} in a factor"))?;
                let base = self.parse(base, kinematics, slots, counters)?;
                let mut total = Tensor::one();
                for _ in 0..exponent {
                    total = contract(&total, &base, kinematics, slots, counters)?;
                }
                Ok(total)
            }
        }
    }
}

/// A sum whose every term is one momentum leg on the same index with a
/// rational coefficient becomes one leg with the combined vector (the vertex
/// rules' `(-Q(2,μ) + Q(3,μ))`); any other sum is kept as it is.
fn merge_momentum_legs(tensor: Tensor, kinematics: &mut Kinematics) -> Result<Tensor, String> {
    if tensor.terms.len() < 2 {
        return Ok(tensor);
    }
    let mut index = None;
    let mut combined = vec![Rational::zero(); kinematics.loops + 1];
    for (structure, poly) in &tensor.terms {
        let [(End::Mom(vector), End::Idx(i))] = structure.as_slice() else {
            return Ok(tensor);
        };
        if index.is_some_and(|seen| seen != *i) {
            return Ok(tensor);
        }
        index = Some(*i);
        let [(mono, value)] = poly.iter().collect::<Vec<_>>()[..] else {
            return Ok(tensor);
        };
        if *mono != Mono::default() || value.0.len() != 1 {
            return Ok(tensor);
        }
        for (target, entry) in combined
            .iter_mut()
            .zip(&kinematics.vectors[*vector as usize])
        {
            *target += &(entry * &value.0[0]);
        }
    }
    let index = index.expect("at least two terms");
    Ok(match kinematics.intern(combined)? {
        None => Tensor::default(),
        Some(vector) => Tensor::structure(vec![leg(End::Mom(vector), End::Idx(index))]),
    })
}

/// A gamma matrix of a fermion factor: `(γ^slot)_{row, column}`.
#[derive(Clone, Debug)]
pub(crate) struct Gamma {
    pub(crate) row: Atom,
    pub(crate) column: Atom,
    pub(crate) lorentz: Atom,
}

/// A fermion factor split into its spinor content and its Lorentz rest.
#[derive(Clone, Debug)]
pub(crate) struct FermionFactor {
    pub(crate) gammas: Vec<Gamma>,
    /// Spinor identities `g(bis, bis)`.
    pub(crate) identities: Vec<(Atom, Atom)>,
    /// Momenta `Q(e, slot)` of the factor (edge id, Lorentz slot).
    pub(crate) momenta: Vec<(i64, Atom)>,
    /// The remaining Lorentz / scalar factors (numbers, metrics, ...).
    pub(crate) rest: Vec<Atom>,
}

/// Split one fermion factor (a single product: v1 QCD rules) into its parts.
pub(crate) fn split_fermion_factor(view: AtomView<'_>) -> Result<FermionFactor, String> {
    let factors = match view {
        AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
        AtomView::Add(_) => {
            return Err(format!(
                "a fermion factor that is a sum (not in v1 QCD rules): {view}"
            ));
        }
        other => vec![other],
    };
    let mut out = FermionFactor {
        gammas: Vec::new(),
        identities: Vec::new(),
        momenta: Vec::new(),
        rest: Vec::new(),
    };
    for factor in factors {
        let AtomView::Fun(function) = factor else {
            out.rest.push(factor.to_owned());
            continue;
        };
        let symbol = function.get_symbol();
        let args = function.iter().collect::<Vec<_>>();
        match (symbol.get_name(), args.as_slice()) {
            (names::GAMMA, [row, column, lorentz]) => out.gammas.push(Gamma {
                row: row.to_owned(),
                column: column.to_owned(),
                lorentz: lorentz.to_owned(),
            }),
            (names::METRIC, [a, b])
                if as_slot(*a).is_some_and(|s| s.rep == names::BIS)
                    && as_slot(*b).is_some_and(|s| s.rep == names::BIS) =>
            {
                out.identities.push((a.to_owned(), b.to_owned()));
            }
            (names::MOMENTUM, [edge, slot]) => {
                let id = i64::try_from(*edge).map_err(|_| format!("a momentum {factor}"))?;
                out.momenta.push((id, slot.to_owned()));
            }
            _ => out.rest.push(factor.to_owned()),
        }
    }
    Ok(out)
}

/// The trace word of one closed fermion loop from its factors: the gammas
/// chained row → column (matrices multiply against the fermion flow; the
/// trace is cyclic), slashes recognised, and the Lorentz rest (unconsumed
/// momenta and other factors) returned for the ordinary contraction.
pub(crate) fn trace_word(factors: Vec<FermionFactor>) -> Result<(Vec<Letter>, Vec<Atom>), String> {
    // Spinor identities: map each slot to a representative.
    let mut representative = HashMap::<Atom, Atom>::default();
    fn find(map: &HashMap<Atom, Atom>, slot: &Atom) -> Atom {
        let mut current = slot.clone();
        while let Some(next) = map.get(&current) {
            if *next == current {
                break;
            }
            current = next.clone();
        }
        current
    }
    for factor in &factors {
        for (a, b) in &factor.identities {
            let (ra, rb) = (find(&representative, a), find(&representative, b));
            if ra != rb {
                representative.insert(ra, rb);
            }
        }
    }
    let mut gammas = Vec::new();
    let mut rest = Vec::new();
    let symbol_q = symbolica::symbol!(names::MOMENTUM);
    for factor in factors {
        let mut momenta = factor.momenta;
        for gamma in factor.gammas {
            // A slash: the gamma's Lorentz slot is carried by exactly one
            // momentum of the same factor.
            let carriers = momenta
                .iter()
                .enumerate()
                .filter(|(_, (_, slot))| *slot == gamma.lorentz)
                .map(|(position, _)| position)
                .collect::<Vec<_>>();
            let letter = if let [position] = carriers[..] {
                let (edge, _) = momenta.remove(position);
                Letter::Slash(edge)
            } else {
                Letter::Index(gamma.lorentz.clone())
            };
            gammas.push((
                find(&representative, &gamma.row),
                find(&representative, &gamma.column),
                letter,
            ));
        }
        for (edge, slot) in momenta {
            rest.push(symbolica::function!(symbol_q, edge, &slot));
        }
        rest.extend(factor.rest);
    }
    if gammas.is_empty() {
        return Err("a fermion loop without gamma matrices".into());
    }
    // Chain: the gamma whose row is the previous column.
    let mut by_row = HashMap::<Atom, usize>::default();
    for (position, (row, _, _)) in gammas.iter().enumerate() {
        if by_row.insert(row.clone(), position).is_some() {
            return Err("two gammas share a row spinor index".into());
        }
    }
    let mut word = Vec::with_capacity(gammas.len());
    let mut used = vec![false; gammas.len()];
    let mut current = 0usize;
    loop {
        if used[current] {
            break;
        }
        used[current] = true;
        word.push(gammas[current].2.clone());
        let column = &gammas[current].1;
        current = *by_row
            .get(column)
            .ok_or("an open spinor chain (not a closed loop)")?;
    }
    if current != 0 || used.iter().any(|u| !u) {
        return Err("the gammas of one fermion loop do not form one cycle".into());
    }
    Ok((word, rest))
}
