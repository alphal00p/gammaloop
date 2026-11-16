//! Direct lowering of a massive closure onto the RustRed p library.
//!
//! The ordinary path ([`crate::p_integral::MassiveClosure::evaluate_batch`])
//! receives `numerator · topo(...)` as one expression: the vacuum numerator
//! with the cut routing substituted and fully expanded. It then splits that
//! expression by topology, normalizes and matches the topology, rewrites the
//! numerator into the canonical momenta of the matched family (another
//! substitution and expansion), splits it into homogeneous loop-dot parts
//! ([`crate::p_integral::parts::split_parts`]), lowers every part onto the
//! slots of the family's presentation and folds the lowered channels by
//! spectator monomial. Every one of those steps builds and re-parses Symbolica
//! expressions term by term.
//!
//! Everything structural is a property of the topology alone: the routing,
//! the match and the slot forms are linear maps of the loop momenta. Here they
//! are computed once per closure and applied to the vacuum numerator by exact
//! polynomial arithmetic:
//!
//! 1. the cut graph is normalized and matched once (numerator 1); the cut
//!    routing composed with the matched momentum images gives every vacuum
//!    loop momentum as a vector over the canonical `k(1..L)` and the head;
//! 2. the numerator is read once, without expanding it: every loop-dependent
//!    `dot` becomes a linear form in the canonical scalar products `k(i)·k(j)`,
//!    `k(i)·Q` plus a multiple of `dot(head, head)`, sums of them are added
//!    and products multiplied as sparse polynomials with exact rational
//!    coefficients, and the loop-free factors are kept as interned spectator
//!    monomials;
//! 3. the canonical terms are classified exactly as the ordinary split
//!    classifies them (the spectator images of the ordinary numerator
//!    rewrites, [`super::p_integral::spectator_terms`]) and merged as the
//!    ordinary canonical numerator merges them;
//! 4. every distinct canonical loop monomial is lowered once onto the slots
//!    and accumulated into the folded channels of
//!    [`super::p_integral::reduce_requests_folded`].
//!
//! The folded channels, keys and rational inputs are those of the ordinary
//! path, so the reductions and values are the same. Whenever an input falls
//! outside what this reproduces exactly (a complex number with real and
//! imaginary parts, an inexact or large number, a loop momentum outside a
//! scalar product, a merged canonical coefficient that the ordinary path
//! would not read as a rational number, ...) the caller falls back to the
//! ordinary path (`None`).

use std::collections::{BTreeMap, BTreeSet};

use ahash::AHashMap as HashMap;
use std::rc::Rc;

use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::coefficient::CoefficientView;
use symbolica::domains::integer::Integer;
use symbolica::domains::rational::Rational;
use symbolica::function;

use super::p_data::{PLibraryHandle, PTop};
use super::p_integral::{Channel, FoldedInputs, spectator_terms};
use crate::p_integral::{P_FAMILIES, PropagatorRequest, external_head, parts};
use crate::symbols::S;
use crate::utils::vakint_macros::vk_symbol;
use crate::{
    EvaluationMethod, ReplacementRules, Topology, Vakint, VakintError, VakintSettings, VakintTerm,
};

/// Bits per canonical scalar product exponent in a packed monomial.
const BITS: u32 = 6;
/// The power of `dot(head, head)` sits in the top byte.
const J_SHIFT: u32 = 120;
/// The phase bit of a spectator id: the term carries the imaginary unit.
const IMAGINARY: u32 = 1 << 31;

/// Statistics of one direct lowering (census).
#[derive(Clone, Debug, Default)]
pub(crate) struct DirectStats {
    pub numerator_terms: usize,
    pub canonical_terms: usize,
    pub merged_terms: usize,
    pub distinct_monomials: usize,
    pub spectators: usize,
    pub lowered: usize,
    pub folded: usize,
    pub read_seconds: f64,
    pub classify_seconds: f64,
    pub monomial_seconds: f64,
    pub scatter_seconds: f64,
    pub finish_seconds: f64,
    pub match_seconds: f64,
    pub reduce_seconds: f64,
    pub finalize_seconds: f64,
}

/// Why the direct lowering declined a closure (the ordinary path runs).
#[derive(Clone, Debug)]
pub(crate) struct Declined(pub &'static str);

type Key = (u32, u128);
type Poly = HashMap<Key, Rational>;

/// The interned loop-free spectator monomials and their products.
#[derive(Default)]
struct Spectators {
    atoms: Vec<Atom>,
    ids: HashMap<Atom, u32>,
    products: HashMap<(u32, u32), u32>,
}

impl Spectators {
    fn new() -> Self {
        let mut spectators = Self::default();
        spectators.intern(Atom::num(1));
        spectators
    }

    fn intern(&mut self, atom: Atom) -> u32 {
        if let Some(id) = self.ids.get(&atom) {
            return *id;
        }
        let id = self.atoms.len() as u32;
        self.atoms.push(atom.clone());
        self.ids.insert(atom, id);
        id
    }

    /// The product of two keyed spectators with their phases, and the sign
    /// `i·i = -1`.
    fn multiply(&mut self, left: u32, right: u32) -> (u32, bool) {
        let phase = (left & IMAGINARY) ^ (right & IMAGINARY);
        let negate = left & IMAGINARY != 0 && right & IMAGINARY != 0;
        let (a, b) = (left & !IMAGINARY, right & !IMAGINARY);
        let id = if a == 0 {
            b
        } else if b == 0 {
            a
        } else {
            let pair = (a.min(b), a.max(b));
            if let Some(id) = self.products.get(&pair) {
                *id
            } else {
                let product = &self.atoms[a as usize] * &self.atoms[b as usize];
                let id = self.intern(product);
                self.products.insert(pair, id);
                id
            }
        };
        (id | phase, negate)
    }
}

/// A rational, real or purely imaginary, number factor.
fn number(view: AtomView) -> Result<(Rational, bool), Declined> {
    let AtomView::Num(number) = view else {
        return Err(Declined("not a number"));
    };
    match number.get_coeff_view() {
        CoefficientView::Natural(nr, dr, 0, _) => Ok((Rational::from((nr, dr)), false)),
        CoefficientView::Natural(0, _, ni, di) => Ok((Rational::from((ni, di)), true)),
        _ => Err(Declined("a complex, large or inexact number")),
    }
}

fn fits_natural(value: &Rational) -> bool {
    matches!(value.numerator_ref(), Integer::Single(_))
        && matches!(value.denominator_ref(), Integer::Single(_))
}

/// vk5l-physrun: polynomials over a coefficient ring.
type PolyC<C> = HashMap<Key, C>;

thread_local! {
    /// Set when an [`i128`] coefficient operation overflowed: the read is
    /// repeated with exact rationals.
    static OVERFLOW: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// vk5l-physrun: the coefficients of the numerator reader: exact rationals,
/// or machine integers (`i128`, checked: an overflow or a non-integer input
/// makes the reader decline and the exact rational reader runs instead).
trait Coef: Clone + std::fmt::Debug {
    fn zero() -> Self;
    fn one() -> Self;
    fn from_i64(value: i64) -> Self;
    fn from_rational(value: &Rational) -> Result<Self, Declined>;
    fn to_rational(&self) -> Rational;
    fn is_zero(&self) -> bool;
    fn add(&mut self, other: &Self);
    fn mul(&self, other: &Self) -> Self;
    fn neg(self) -> Self;
}

impl Coef for Rational {
    fn zero() -> Self {
        Rational::zero()
    }
    fn one() -> Self {
        Rational::one()
    }
    fn from_i64(value: i64) -> Self {
        Rational::from(value)
    }
    fn from_rational(value: &Rational) -> Result<Self, Declined> {
        Ok(value.clone())
    }
    fn to_rational(&self) -> Rational {
        self.clone()
    }
    fn is_zero(&self) -> bool {
        Rational::is_zero(self)
    }
    fn add(&mut self, other: &Self) {
        *self += other;
    }
    fn mul(&self, other: &Self) -> Self {
        self * other
    }
    fn neg(self) -> Self {
        -self
    }
}

impl Coef for i128 {
    fn zero() -> Self {
        0
    }
    fn one() -> Self {
        1
    }
    fn from_i64(value: i64) -> Self {
        value as i128
    }
    fn from_rational(value: &Rational) -> Result<Self, Declined> {
        match (value.numerator_ref(), value.denominator_ref()) {
            (Integer::Single(n), Integer::Single(1)) => Ok(*n as i128),
            _ => Err(Declined("a non-integer coefficient for the machine-integer reader")),
        }
    }
    fn to_rational(&self) -> Rational {
        Rational::from(Integer::from(*self))
    }
    fn is_zero(&self) -> bool {
        *self == 0
    }
    fn add(&mut self, other: &Self) {
        match self.checked_add(*other) {
            Some(sum) => *self = sum,
            None => OVERFLOW.with(|flag| flag.set(true)),
        }
    }
    fn mul(&self, other: &Self) -> Self {
        self.checked_mul(*other).unwrap_or_else(|| {
            OVERFLOW.with(|flag| flag.set(true));
            0
        })
    }
    fn neg(self) -> Self {
        self.checked_neg().unwrap_or_else(|| {
            OVERFLOW.with(|flag| flag.set(true));
            0
        })
    }
}

fn add_into<C: Coef>(target: &mut HashMap<Key, C>, key: Key, value: C) {
    use std::collections::hash_map::Entry;
    match target.entry(key) {
        Entry::Occupied(mut entry) => {
            entry.get_mut().add(&value);
            if entry.get().is_zero() {
                entry.remove();
            }
        }
        Entry::Vacant(entry) => {
            if !value.is_zero() {
                entry.insert(value);
            }
        }
    }
}

/// The degree of a packed monomial in canonical scalar products.
fn degree(packed: u128, count: usize) -> u32 {
    (0..count)
        .map(|index| ((packed >> (BITS * index as u32)) & ((1 << BITS) - 1)) as u32)
        .sum()
}

/// vk5l-physrun: `VAKINT_DIRECT_READER=prefix` multiplies prefix products
/// into each group's coefficient (top-down) instead of the Horner evaluation.
fn horner_reader() -> bool {
    static ON: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| {
        !matches!(std::env::var("VAKINT_DIRECT_READER").as_deref(), Ok("prefix"))
    });
    *ON
}

/// vk5l-physrun: `VAKINT_DIRECT_PREFIX_TRIE=off` restores one product per
/// group in [`Reader::read`] (comparison runs).
fn prefix_trie_enabled() -> bool {
    static ON: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| {
        !matches!(std::env::var("VAKINT_DIRECT_PREFIX_TRIE").as_deref(), Ok("off" | "0"))
    });
    *ON
}

/// The numerator reader of one closure.
struct Reader<'a, C: Coef> {
    /// Vacuum (or cut) loop momentum `k(i)` -> its vector over the canonical
    /// `k(1..L)` and the head.
    momenta: &'a HashMap<Atom, Vec<C>>,
    loops: usize,
    sigma: i64,
    /// The number of canonical scalar products `k(i)·k(j)` (i <= j), `k(i)·Q`.
    count: usize,
    spectators: &'a mut Spectators,
    cache: HashMap<Atom, Rc<PolyC<C>>>,
    terms: usize,
}

impl<C: Coef> Reader<'_, C> {
    fn pair_index(&self, i: usize, j: usize) -> usize {
        let (i, j) = (i.min(j), i.max(j));
        // Upper triangle by rows: (0,0..L-1), (1,1..L-1), ...
        i * self.loops - i * i.saturating_sub(1) / 2 + (j - i)
    }

    fn unit() -> PolyC<C> {
        PolyC::<C>::from([((0, 0), C::one())])
    }

    fn multiply(&mut self, left: &PolyC<C>, right: &PolyC<C>) -> Result<PolyC<C>, Declined> {
        let mut product = PolyC::<C>::with_capacity(left.len().max(right.len()));
        self.multiply_into(&mut product, left, right)?;
        Ok(product)
    }

    /// vk5l-physrun: `target += left · right`, without a product map of its own.
    fn multiply_into(&mut self, product: &mut PolyC<C>, left: &PolyC<C>, right: &PolyC<C>) -> Result<(), Declined> {
        for ((ls, lp), lv) in left {
            for ((rs, rp), rv) in right {
                let packed = lp + rp;
                if degree(packed, self.count) > 60 || (packed >> J_SHIFT) > 250 {
                    return Err(Declined("a numerator degree beyond the packed range"));
                }
                let (spectator, negate) = self.spectators.multiply(*ls, *rs);
                let mut value = lv.mul(rv);
                if negate {
                    value = value.neg();
                }
                add_into(product, (spectator, packed), value);
            }
        }
        Ok(())
    }

    fn power(&mut self, base: &PolyC<C>, exponent: i64) -> Result<PolyC<C>, Declined> {
        let mut result = Self::unit();
        for _ in 0..exponent {
            result = self.multiply(&result, base)?;
        }
        Ok(result)
    }

    /// A linear momentum (a loop momentum `k(i)` or a sum of them) as its
    /// vector.
    fn vector(&self, momentum: AtomView) -> Result<Vec<C>, Declined> {
        if let Some(vector) = self.momenta.get(&momentum.to_owned()) {
            return Ok(vector.clone());
        }
        let mut total = vec![C::zero(); self.loops + 1];
        let mut add_term = |term: AtomView| -> Result<(), Declined> {
            let (scale, head) = match term {
                AtomView::Mul(product) => {
                    let factors = product.iter().collect::<Vec<_>>();
                    match factors.as_slice() {
                        [n @ AtomView::Num(_), head] | [head, n @ AtomView::Num(_)] => {
                            let (value, imaginary) = number(*n)?;
                            if imaginary {
                                return Err(Declined("an imaginary momentum coefficient"));
                            }
                            (C::from_rational(&value)?, *head)
                        }
                        _ => return Err(Declined("a nonlinear momentum")),
                    }
                }
                other => (C::one(), other),
            };
            let vector = self
                .momenta
                .get(&head.to_owned())
                .ok_or(Declined("an unknown momentum in a scalar product"))?;
            for (slot, value) in total.iter_mut().zip(vector) {
                slot.add(&scale.mul(value));
            }
            Ok(())
        };
        match momentum {
            AtomView::Add(sum) => sum.iter().try_for_each(&mut add_term)?,
            other => add_term(other)?,
        }
        Ok(total)
    }

    /// A loop-dependent `dot(a, b)` as a linear polynomial.
    fn dot(&mut self, a: AtomView, b: AtomView) -> Result<PolyC<C>, Declined> {
        let (va, vb) = (self.vector(a)?, self.vector(b)?);
        let l = self.loops;
        let mut poly = PolyC::<C>::new();
        for i in 0..l {
            for j in 0..l {
                let value = va[i].mul(&vb[j]);
                if !value.is_zero() {
                    let index = self.pair_index(i, j);
                    add_into(&mut poly, (0, 1u128 << (BITS * index as u32)), value);
                }
            }
        }
        let external = l * (l + 1) / 2;
        let sigma = C::from_i64(self.sigma);
        for i in 0..l {
            let mut value = va[i].mul(&vb[l]);
            value.add(&va[l].mul(&vb[i]));
            let value = value.mul(&sigma);
            if !value.is_zero() {
                add_into(&mut poly, (0, 1u128 << (BITS * (external + i) as u32)), value);
            }
        }
        let value = va[l].mul(&vb[l]);
        if !value.is_zero() {
            add_into(&mut poly, (0, 1u128 << J_SHIFT), value);
        }
        Ok(poly)
    }

    /// Whether a factor depends on a loop momentum.
    fn has_loop(view: AtomView) -> bool {
        view.contains_symbol(S.k)
    }

    fn cached(&mut self, view: AtomView) -> Result<Rc<PolyC<C>>, Declined> {
        let cacheable = !matches!(view, AtomView::Add(_)) || view.get_byte_size() <= 4096;
        if cacheable && let Some(poly) = self.cache.get(&view.to_owned()) {
            return Ok(Rc::clone(poly));
        }
        let poly = Rc::new(self.read(view)?);
        if cacheable {
            self.cache.insert(view.to_owned(), Rc::clone(&poly));
        }
        Ok(poly)
    }

    /// The polynomial of an expression.
    fn read(&mut self, view: AtomView) -> Result<PolyC<C>, Declined> {
        match view {
            AtomView::Add(sum) => {
                // Summands grouped by their loop-dependent factors: the
                // loop-free coefficients of a group are added first and the
                // group's loop product is multiplied once.
                let mut order = Vec::<Atom>::new();
                let mut groups = HashMap::<Atom, PolyC<C>>::new();
                for term in sum.iter() {
                    self.terms += 1;
                    let factors = match term {
                        AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
                        other => vec![other],
                    };
                    let (loops, rest): (Vec<_>, Vec<_>) =
                        factors.into_iter().partition(|factor| Self::has_loop(*factor));
                    let key = match loops.as_slice() {
                        [] => Atom::num(1),
                        [single] => single.to_owned(),
                        many => Atom::mul_many(many.iter().copied()),
                    };
                    let coefficient = self.term_factors(rest)?;
                    match groups.get_mut(&key) {
                        Some(group) => {
                            for (k, value) in coefficient {
                                add_into(group, k, value);
                            }
                        }
                        None => {
                            order.push(key.clone());
                            groups.insert(key, coefficient);
                        }
                    }
                }
                if !prefix_trie_enabled() {
                    let mut total = PolyC::<C>::new();
                    for key in order {
                        let coefficient = groups.remove(&key).unwrap_or_default();
                        if coefficient.is_empty() {
                            continue;
                        }
                        let product = if key == Atom::num(1) {
                            coefficient
                        } else {
                            let loops = self.term(key.as_view())?;
                            self.multiply(&coefficient, &loops)?
                        };
                        for (k, value) in product {
                            add_into(&mut total, k, value);
                        }
                    }
                    return Ok(total);
                }
                // vk5l-physrun: the loop products of the groups share factors.
                // Each group's loop key is cut into its factors, the factors
                // are numbered by how many groups use them (most used first),
                // every group lists its factor numbers in increasing order, and
                // the groups are visited in the order of those lists: a stack
                // keeps the product of the current prefix, so a prefix shared
                // by consecutive groups is multiplied once. Exact sums: the
                // total is the same polynomial.
                let mut factor_ids = HashMap::<Atom, u32>::new();
                let mut factor_atoms = Vec::<Atom>::new();
                let mut uses = Vec::<usize>::new();
                let mut lists = Vec::<(Vec<u32>, PolyC<C>)>::with_capacity(order.len());
                let mut total = PolyC::<C>::new();
                for key in order {
                    let coefficient = groups.remove(&key).unwrap_or_default();
                    if coefficient.is_empty() {
                        continue;
                    }
                    if key == Atom::num(1) {
                        for (k, value) in coefficient {
                            add_into(&mut total, k, value);
                        }
                        continue;
                    }
                    let factors = match key.as_view() {
                        AtomView::Mul(product) => product.iter().map(|f| f.to_owned()).collect::<Vec<_>>(),
                        other => vec![other.to_owned()],
                    };
                    let mut list = Vec::with_capacity(factors.len());
                    for factor in factors {
                        let next = factor_atoms.len() as u32;
                        let id = *factor_ids.entry(factor.clone()).or_insert_with(|| {
                            factor_atoms.push(factor);
                            uses.push(0);
                            next
                        });
                        uses[id as usize] += 1;
                        list.push(id);
                    }
                    lists.push((list, coefficient));
                }
                // Rank by use (most used first), ties by first appearance.
                let mut ranked = (0..factor_atoms.len() as u32).collect::<Vec<_>>();
                ranked.sort_by_key(|id| (std::cmp::Reverse(uses[*id as usize]), *id));
                let mut rank = vec![0u32; factor_atoms.len()];
                for (position, id) in ranked.iter().enumerate() {
                    rank[*id as usize] = position as u32;
                }
                for (list, _) in &mut lists {
                    for id in list.iter_mut() {
                        *id = rank[*id as usize];
                    }
                    list.sort_unstable();
                }
                lists.sort_by(|left, right| left.0.cmp(&right.0));
                if horner_reader() {
                    // vk5l-physrun: Horner over the trie of sorted factor lists,
                    // bottom-up: value(node) = Σ coefficients of the groups ending
                    // at the node + Σ_children factor · value(child); every factor
                    // multiplies a merged subtree sum once. Exact: the same total.
                    let factors = ranked
                        .iter()
                        .map(|id| self.cached(factor_atoms[*id as usize].as_view()))
                        .collect::<Result<Vec<_>, Declined>>()?;
                    let value = self.horner(&lists, 0, &factors)?;
                    for (k, v) in value {
                        add_into(&mut total, k, v);
                    }
                    return Ok(total);
                }
                let mut stack: Vec<(u32, Rc<PolyC<C>>)> = Vec::new();
                let unit = Rc::new(Self::unit());
                for (list, coefficient) in lists {
                    let common = stack
                        .iter()
                        .zip(&list)
                        .take_while(|((id, _), other)| id == *other)
                        .count();
                    stack.truncate(common);
                    for &id in &list[common..] {
                        let factor = self.cached(factor_atoms[ranked[id as usize] as usize].as_view())?;
                        let previous = stack.last().map(|(_, poly)| Rc::clone(poly)).unwrap_or_else(|| Rc::clone(&unit));
                        let product = self.multiply(&previous, &factor)?;
                        stack.push((id, Rc::new(product)));
                    }
                    let loops = stack.last().map(|(_, poly)| Rc::clone(poly)).unwrap_or_else(|| Rc::clone(&unit));
                    self.multiply_into(&mut total, &coefficient, &loops)?;
                }
                Ok(total)
            }
            AtomView::Mul(_) => self.term(view),
            _ if !Self::has_loop(view) => self.term(view),
            AtomView::Fun(function) if function.get_symbol() == S.dot => self.dot_factor(view),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let exponent = crate::get_integer_from_atom(exponent)
                    .filter(|exponent| *exponent > 0)
                    .ok_or(Declined("a non-polynomial loop factor"))?;
                let base = self.cached(base)?;
                self.power(&base, exponent)
            }
            _ => Err(Declined("a loop momentum outside a scalar product")),
        }
    }


    /// vk5l-physrun: the Horner value of the sorted factor lists `lists`
    /// (all sharing their first `depth` factors): Σ of the coefficients of the
    /// lists that end at `depth`, plus Σ over the next factor `f` of
    /// `factors[f] · horner(lists continuing with f, depth + 1)`.
    #[allow(clippy::type_complexity)]
    fn horner(
        &mut self,
        lists: &[(Vec<u32>, PolyC<C>)],
        depth: usize,
        factors: &[Rc<PolyC<C>>],
    ) -> Result<PolyC<C>, Declined> {
        let mut value = PolyC::<C>::new();
        let mut start = 0;
        while start < lists.len() && lists[start].0.len() == depth {
            for (k, v) in &lists[start].1 {
                add_into(&mut value, *k, v.clone());
            }
            start += 1;
        }
        while start < lists.len() {
            let factor = lists[start].0[depth];
            let mut end = start + 1;
            while end < lists.len() && lists[end].0[depth] == factor {
                end += 1;
            }
            let below = self.horner(&lists[start..end], depth + 1, factors)?;
            self.multiply_into(&mut value, &factors[factor as usize], &below)?;
            start = end;
        }
        Ok(value)
    }

    /// One summand: numbers, spectators and loop-dependent factors.
    fn term(&mut self, term: AtomView) -> Result<PolyC<C>, Declined> {
        let factors = match term {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            other => vec![other],
        };
        self.term_factors(factors)
    }

    /// The product of `factors` (one summand, or part of it).
    fn term_factors(&mut self, factors: Vec<AtomView<'_>>) -> Result<PolyC<C>, Declined> {
        let mut value = C::one();
        let mut imaginary = false;
        let mut spectators = Vec::new();
        let mut loops = Vec::new();
        for factor in factors {
            match factor {
                AtomView::Num(_) => {
                    let (number, phase) = number(factor)?;
                    value = value.mul(&C::from_rational(&number)?);
                    if phase {
                        if imaginary {
                            value = value.neg();
                        }
                        imaginary = !imaginary;
                    }
                }
                _ if !Self::has_loop(factor) => spectators.push(factor),
                _ => loops.push(self.cached(factor)?),
            }
        }
        let spectator = match spectators.as_slice() {
            [] => 0,
            [single] => self.spectators.intern(single.to_owned()),
            many => self.spectators.intern(Atom::mul_many(many.iter().copied())),
        };
        let spectator = if imaginary {
            spectator | IMAGINARY
        } else {
            spectator
        };
        // Smallest factors first.
        loops.sort_by_key(|poly| poly.len());
        let mut product = PolyC::<C>::from([((spectator, 0), value)]);
        for poly in loops {
            product = self.multiply(&product, &poly)?;
        }
        Ok(product)
    }

    fn dot_factor(&mut self, view: AtomView) -> Result<PolyC<C>, Declined> {
        let AtomView::Fun(function) = view else {
            unreachable!()
        };
        let arguments = function.iter().collect::<Vec<_>>();
        let [a, b] = arguments.as_slice() else {
            return Err(Declined("a dot with other than two arguments"));
        };
        self.dot(*a, *b)
    }
}

/// vk5l-physrun: the canonical polynomial of a closure numerator (times the
/// cut graph's factor) read with coefficients `C`, as exact rationals.
#[allow(clippy::too_many_arguments)]
fn read_numerator<C: Coef>(
    vacuum_momenta: &HashMap<Atom, Vec<Rational>>,
    cut_momenta: &HashMap<Atom, Vec<Rational>>,
    loops: usize,
    sigma: i64,
    count: usize,
    spectators: &mut Spectators,
    numerator: AtomView,
    factor: &Atom,
    stats: &mut DirectStats,
) -> Result<Poly, Declined> {
    let convert = |momenta: &HashMap<Atom, Vec<Rational>>| -> Result<HashMap<Atom, Vec<C>>, Declined> {
        momenta
            .iter()
            .map(|(momentum, vector)| {
                Ok((
                    momentum.clone(),
                    vector.iter().map(C::from_rational).collect::<Result<Vec<_>, _>>()?,
                ))
            })
            .collect()
    };
    let vacuum = convert(vacuum_momenta)?;
    let cut = convert(cut_momenta)?;
    let mut reader = Reader::<C> {
        momenta: &vacuum,
        loops,
        sigma,
        count,
        spectators,
        cache: HashMap::new(),
        terms: 0,
    };
    let numerator_poly = reader.read(numerator)?;
    stats.numerator_terms += reader.terms;
    let polynomial = if *factor == Atom::num(1) {
        numerator_poly
    } else {
        let mut factor_reader = Reader::<C> {
            momenta: &cut,
            loops,
            sigma,
            count,
            spectators: reader.spectators,
            cache: HashMap::new(),
            terms: 0,
        };
        let factor_poly = factor_reader.read(factor.as_view())?;
        factor_reader.multiply(&numerator_poly, &factor_poly)?
    };
    Ok(polynomial
        .into_iter()
        .map(|(key, value)| (key, value.to_rational()))
        .collect())
}

/// The integer (rational) coefficients of a linear momentum over the
/// canonical `k(1..L)` and the user head.
fn canonical_vector(
    momentum: &Atom,
    loops: usize,
    head: &Atom,
) -> Result<Vec<Rational>, Declined> {
    let expanded = momentum.expand();
    let mut vector = vec![Rational::zero(); loops + 1];
    let mut add = |term: AtomView| -> Result<(), Declined> {
        let (scale, vector_atom) = match term {
            AtomView::Mul(product) => {
                let factors = product.iter().collect::<Vec<_>>();
                match factors.as_slice() {
                    [n @ AtomView::Num(_), v] | [v, n @ AtomView::Num(_)] => {
                        let (value, imaginary) = number(*n)?;
                        if imaginary {
                            return Err(Declined("an imaginary routing coefficient"));
                        }
                        (value, *v)
                    }
                    _ => return Err(Declined("a nonlinear routing")),
                }
            }
            other => (Rational::one(), other),
        };
        let slot = if vector_atom == head.as_view() {
            loops
        } else {
            let AtomView::Fun(function) = vector_atom else {
                return Err(Declined("a routing term is not a momentum"));
            };
            if function.get_symbol() != S.k || function.get_nargs() != 1 {
                return Err(Declined("a routing term is not a canonical loop momentum"));
            }
            let index = crate::get_integer_from_atom(function.iter().next().unwrap())
                .ok_or(Declined("a symbolic loop momentum index"))?;
            if index < 1 || index as usize > loops {
                return Err(Declined("a loop momentum outside the canonical range"));
            }
            index as usize - 1
        };
        vector[slot] += &scale;
        Ok(())
    };
    match expanded.as_view() {
        AtomView::Add(sum) => sum.iter().try_for_each(&mut add)?,
        v if v.is_zero() => {}
        other => add(other)?,
    }
    Ok(vector)
}

/// A matched cut graph: everything the lowering needs that does not depend
/// on the numerator.
pub(crate) struct MatchedCut {
    pub request_family: &'static crate::p_integral::PFamilyData,
    pub powers: Vec<i64>,
    pub external: Atom,
    pub head: Atom,
    pub sigma: i64,
    /// The normalization factor of the cut graph (tree factors, negative
    /// powers), in the cut momenta.
    pub factor: Atom,
    pub rules: ReplacementRules,
    pub scaleless: bool,
}

/// Normalize and match the cut graph `topology` once, as the ordinary path
/// does for every closure term. `Ok(None)` when the normalized graph is
/// scaleless (an exact zero).
pub(crate) fn match_cut(
    vakint: &Vakint,
    settings: &VakintSettings,
    topology: &Atom,
    closure_head: &Atom,
) -> Result<Result<Option<MatchedCut>, Declined>, VakintError> {
    let mut term = VakintTerm {
        integral: topology.clone(),
        numerator: Atom::num(1),
        vectors: Vec::new(),
    };
    if !term.normalize_propagator()? {
        return Ok(Ok(None));
    }
    vakint.validate_settings(settings)?;
    let external = term
        .propagator_external()?
        .ok_or_else(|| crate::p_integral::invalid("a massive closure requires a massless p-integral"))?;
    if &external != closure_head {
        return Err(crate::p_integral::invalid(format!(
            "the cut graph uses {external}, but the closing line uses {closure_head}"
        )));
    }
    let factor = term.numerator.clone();
    term.numerator = Atom::num(1);
    let matched = term.prepare_integral_specs(vakint, settings, None)?;
    if !matches!(matched.canonical_topology, Topology::Propagator(_)) {
        return Ok(Err(Declined("not a p-integral topology")));
    }
    match vakint.evaluation_method(settings, &matched.canonical_topology) {
        Some(EvaluationMethod::RustRed(_)) if crate::p_integral::massive::rustred_max_rank().is_none() => {}
        _ => return Ok(Err(Declined("not reduced by RustRed alone"))),
    }
    let one = Atom::num(1);
    let request = matched
        .propagator_request(one.as_view())
        .ok_or_else(|| crate::p_integral::invalid("a p-integral match has a propagator request"))?;
    let (head, sigma) = external_head(&request.external)?;
    let scaleless = parts::split(
        request.family,
        request.powers,
        Atom::num(1).as_view(),
        &request.external,
    )?
    .is_none();
    Ok(Ok(Some(MatchedCut {
        request_family: request.family,
        powers: request.powers.to_vec(),
        external: request.external.clone(),
        head,
        sigma,
        factor,
        rules: matched,
        scaleless,
    })))
}

/// vk5l-physrun (shared path, rule 17): a plain p-integral request (any loop
/// count) through the direct reader and the Horner slot substitution: its
/// canonical numerator is already in the canonical momenta and the user head,
/// so the cut is the request itself (no replacement rules, identity routing).
/// `Err(Declined)`: the ordinary split and lowering run.
pub(crate) fn lower_request(
    settings: &VakintSettings,
    handle: &PLibraryHandle,
    request: &PropagatorRequest,
    stats: &mut DirectStats,
) -> Result<Result<BTreeMap<Channel, FoldedInputs>, Declined>, VakintError> {
    let (head, sigma) = external_head(&request.external)?;
    let scaleless = parts::split(
        request.family,
        request.powers,
        Atom::num(1).as_view(),
        &request.external,
    )?
    .is_none();
    let cut = MatchedCut {
        request_family: request.family,
        powers: request.powers.to_vec(),
        external: request.external.clone(),
        head: head.clone(),
        sigma,
        factor: Atom::num(1),
        rules: ReplacementRules::default(),
        scaleless,
    };
    let mut routing = (1..=request.family.loops as i64)
        .map(|index| (function!(S.k, index), function!(S.k, index)))
        .collect::<Vec<_>>();
    routing.push((head.clone(), head));
    lower_closure(settings, handle, &cut, &routing, request.canonical_numerator, stats)
}

/// The folded channels of one closure, lowered directly; `Err(Declined)` for
/// the ordinary path.
#[allow(clippy::too_many_arguments)]
pub(crate) fn lower_closure(
    settings: &VakintSettings,
    handle: &PLibraryHandle,
    cut: &MatchedCut,
    routing: &[(Atom, Atom)],
    numerator: AtomView,
    stats: &mut DirectStats,
) -> Result<Result<BTreeMap<Channel, FoldedInputs>, Declined>, VakintError> {
    if cut.scaleless {
        return Ok(Ok(BTreeMap::new()));
    }
    // vk5l-physrun: a negative line power is a numerator factor of the
    // ordinary split (`negative_powers` in parts::split_with), which this
    // lowering does not form: the ordinary path takes the request.
    if cut.powers.iter().any(|power| *power < 0) {
        return Ok(Err(Declined("a negative line power")));
    }
    if !crate::utils::has_only_exact_numbers(numerator) {
        return Ok(Err(Declined("inexact numbers")));
    }
    let family = P_FAMILIES
        .iter()
        .position(|family| family.name == cut.request_family.name)
        .ok_or_else(|| crate::p_integral::invalid(format!("unknown family {}", cut.request_family.name)))?;
    let top: PTop = handle.top(family)?;
    let presentation = handle.library.presentation(top.presentation);
    let loops = cut.request_family.loops;
    if usize::from(presentation.loops) != loops
        || top.arity as u32 * BITS > 128
        || (loops * (loops + 1) / 2 + loops) as u32 * BITS > J_SHIFT
    {
        return Ok(Err(Declined("a presentation beyond the packed range")));
    }
    // The cut momenta in canonical coordinates: each cut loop momentum through
    // the matched images.
    let canonical = |momentum: &Atom| -> Result<Result<Vec<Rational>, Declined>, VakintError> {
        let mut probe = VakintTerm {
            integral: Atom::num(1),
            numerator: momentum.clone(),
            vectors: Vec::new(),
        };
        probe.apply_numerator_replacement_rules(&cut.rules, settings)?;
        Ok(canonical_vector(&probe.numerator, loops, &cut.head))
    };
    let mut cut_momenta = HashMap::new();
    for index in 1..=loops as i64 {
        let momentum = function!(S.k, index);
        match canonical(&momentum)? {
            Ok(vector) => {
                cut_momenta.insert(momentum, vector);
            }
            Err(declined) => return Ok(Err(declined)),
        }
    }
    let mut head_vector = vec![Rational::zero(); loops + 1];
    head_vector[loops] = Rational::one();
    cut_momenta.insert(cut.head.clone(), head_vector);
    let mut vacuum_momenta = HashMap::new();
    for (vacuum, image) in routing {
        let reader_momenta = &cut_momenta;
        let expanded = image.expand();
        let mut total = vec![Rational::zero(); loops + 1];
        let mut add = |term: AtomView| -> Result<(), Declined> {
            let (scale, head) = match term {
                AtomView::Mul(product) => {
                    let factors = product.iter().collect::<Vec<_>>();
                    match factors.as_slice() {
                        [n @ AtomView::Num(_), v] | [v, n @ AtomView::Num(_)] => {
                            let (value, imaginary) = number(*n)?;
                            if imaginary {
                                return Err(Declined("an imaginary routing coefficient"));
                            }
                            (value, *v)
                        }
                        _ => return Err(Declined("a nonlinear routing")),
                    }
                }
                other => (Rational::one(), other),
            };
            let vector = reader_momenta
                .get(&head.to_owned())
                .ok_or(Declined("a routing image outside the cut momenta"))?;
            for (slot, value) in total.iter_mut().zip(vector) {
                *slot += &(&scale * value);
            }
            Ok(())
        };
        let outcome = match expanded.as_view() {
            AtomView::Add(sum) => sum.iter().try_for_each(&mut add),
            other => add(other),
        };
        if let Err(declined) = outcome {
            return Ok(Err(declined));
        }
        vacuum_momenta.insert(vacuum.clone(), total);
    }
    let count = loops * (loops + 1) / 2 + loops;
    let phase = std::time::Instant::now();
    let mut spectators = Spectators::new();
    // vk5l-physrun: machine-integer coefficients first (no gcd per operation,
    // small entries); exact rationals when an input is not an integer or an
    // operation overflows. Same polynomial either way.
    let polynomial = match read_numerator::<i128>(
        &vacuum_momenta,
        &cut_momenta,
        loops,
        cut.sigma,
        count,
        &mut spectators,
        numerator,
        &cut.factor,
        stats,
    ) {
        Ok(polynomial) if !OVERFLOW.with(|flag| flag.replace(false)) => polynomial,
        _ => {
            OVERFLOW.with(|flag| flag.set(false));
            match read_numerator::<Rational>(
                &vacuum_momenta,
                &cut_momenta,
                loops,
                cut.sigma,
                count,
                &mut spectators,
                numerator,
                &cut.factor,
                stats,
            ) {
                Ok(polynomial) => polynomial,
                Err(declined) => return Ok(Err(declined)),
            }
        }
    };
    stats.canonical_terms += polynomial.len();
    stats.read_seconds += phase.elapsed().as_secs_f64();
    let phase = std::time::Instant::now();

    // Classify the canonical terms as the ordinary split does: the spectator
    // image of the ordinary numerator rewrites, split into
    // (monomial, ε-basis, value) terms, times `dot(head, head)^j`, merged per
    // canonical loop monomial and symbolic part.
    let epsilon = vk_symbol!(&settings.epsilon_symbol);
    let head_square = function!(S.dot, &cut.head, &cut.head);
    let mut images = HashMap::<u32, Rc<Vec<(Atom, Atom, Rational)>>>::new();
    let mut monomials = Spectators::new();
    let mut bases = Spectators::new();
    let mut powered = HashMap::<(u32, u32), u32>::new();
    // (loop monomial, monomial id, basis id) -> (real, imaginary) value.
    let mut merged = HashMap::<(u128, u32, u32), (Rational, Rational)>::new();
    let loop_mask = (1u128 << J_SHIFT) - 1;
    for ((spectator, packed), value) in polynomial {
        let id = spectator & !IMAGINARY;
        let imaginary = spectator & IMAGINARY != 0;
        let image = match images.get(&id) {
            Some(image) => Rc::clone(image),
            None => {
                let atom = &spectators.atoms[id as usize];
                let transformed = match spectator_image(settings, &cut.rules, &cut.head, atom)? {
                    Ok(transformed) => transformed,
                    Err(declined) => return Ok(Err(declined)),
                };
                let terms = Rc::new(spectator_terms(&transformed, epsilon));
                images.insert(id, Rc::clone(&terms));
                terms
            }
        };
        let j = (packed >> J_SHIFT) as u32;
        let loop_monomial = packed & loop_mask;
        for (monomial, basis, coefficient) in image.iter() {
            let monomial_id = monomials.intern(monomial.clone());
            let monomial_id = *powered.entry((monomial_id, j)).or_insert_with(|| {
                if j == 0 {
                    monomial_id
                } else {
                    let atom = &monomials.atoms[monomial_id as usize]
                        * head_square.pow(Atom::num(j as i64));
                    monomials.intern(atom)
                }
            });
            let basis_id = bases.intern(basis.clone());
            let entry = merged
                .entry((loop_monomial, monomial_id, basis_id))
                .or_insert_with(|| (Rational::zero(), Rational::zero()));
            let contribution = &value * coefficient;
            if imaginary {
                entry.1 += &contribution;
            } else {
                entry.0 += &contribution;
            }
        }
    }
    stats.spectators += images.len();
    stats.merged_terms += merged.len();
    stats.classify_seconds += phase.elapsed().as_secs_f64();
    let phase = std::time::Instant::now();

    // Lower every distinct loop monomial once and accumulate the folded
    // channels.
    let total_power: i64 = cut.powers.iter().filter(|power| **power > 0).sum();
    let base = (0..top.arity)
        .map(|slot| {
            if slot < top.line_count {
                cut.powers[slot].max(0)
            } else {
                0
            }
        })
        .collect::<Vec<_>>();
    let mut forms = Vec::with_capacity(count);
    for index in 0..count {
        let (a, b) = canonical_rows(index, loops);
        match presentation.express_scalar_product(&a, &b) {
            Ok(form) => forms.push(form),
            Err(_) => return Ok(Err(Declined("a scalar product the presentation cannot express"))),
        }
    }
    // Gather the canonical loop polynomial of every bucket
    // (monomial, phase, scale power, basis), then substitute the slot forms
    // into each bucket polynomial at once (multivariate Horner): sums are
    // merged after every multiplication, so the product of forms of each
    // monomial is never formed on its own. The lowered keys and exact values
    // are those of lowering every monomial and summing.
    let mut key_ids = HashMap::<Vec<i64>, usize>::new();
    let mut keys = Vec::<Vec<i64>>::new();
    let imaginary_unit = Atom::i();
    let mut bucket_ids = HashMap::<(u32, bool, i64, u32), usize>::new();
    let mut bucket_terms = Vec::<Vec<(u128, Rational)>>::new();
    let mut distinct = ahash::AHashSet::<u128>::new();
    for ((loop_monomial, monomial_id, basis_id), (real, imaginary)) in merged {
        let (value, phase) = match (real.is_zero(), imaginary.is_zero()) {
            (true, true) => continue,
            (false, true) => (real, false),
            (true, false) => (imaginary, true),
            (false, false) => {
                return Ok(Err(Declined("a canonical coefficient with real and imaginary parts")));
            }
        };
        if !fits_natural(&value) {
            return Ok(Err(Declined("a canonical coefficient beyond a machine rational")));
        }
        let r = degree(loop_monomial, count) as i64;
        let scale_power = 2 * loops as i64 - total_power + r;
        distinct.insert(loop_monomial);
        let next = bucket_terms.len();
        let bucket = *bucket_ids
            .entry((monomial_id, phase, scale_power, basis_id))
            .or_insert_with(|| {
                bucket_terms.push(Vec::new());
                next
            });
        bucket_terms[bucket].push((loop_monomial, value));
    }
    stats.distinct_monomials += distinct.len();
    let mut buckets = Vec::<Vec<Rational>>::with_capacity(bucket_terms.len());
    for mut terms in bucket_terms {
        let substituted = substitute(&mut terms, count, &forms);
        let mut bucket = Vec::new();
        for (monomial, value) in substituted {
            if value.is_zero() {
                continue;
            }
            let key = (0..top.arity)
                .map(|slot| {
                    base[slot] - ((monomial >> (BITS * slot as u32)) & ((1 << BITS) - 1)) as i64
                })
                .collect::<Vec<_>>();
            let next = keys.len();
            let id = *key_ids.entry(key.clone()).or_insert_with(|| {
                keys.push(key);
                next
            });
            if bucket.len() <= id {
                bucket.resize(id + 1, Rational::zero());
            }
            bucket[id] += &value;
            stats.lowered += 1;
        }
        buckets.push(bucket);
    }
    stats.scatter_seconds += phase.elapsed().as_secs_f64();
    let phase = std::time::Instant::now();
    let integral = super::p_integral::integral_name(&PropagatorRequest {
        family: cut.request_family,
        powers: &cut.powers,
        canonical_numerator: Atom::num(1).as_view(),
        external: cut.external.clone(),
    });
    let mut folded = BTreeMap::<Channel, FoldedInputs>::new();
    for ((monomial_id, phase, scale_power, basis_id), bucket) in bucket_ids {
        let inputs = std::mem::take(&mut buckets[bucket])
            .into_iter()
            .enumerate()
            .filter(|(_, value)| !value.is_zero())
            .map(|(id, value)| ((top.presentation, keys[id].clone()), value))
            .collect::<BTreeMap<_, _>>();
        if inputs.is_empty() {
            continue;
        }
        let monomial = &monomials.atoms[monomial_id as usize];
        let coefficient = if phase {
            monomial * &imaginary_unit
        } else {
            monomial.clone()
        };
        let target = folded
            .entry(Channel {
                loops,
                head: cut.head.clone(),
                scale_power,
                coefficient,
            })
            .or_default();
        target.integrals = BTreeSet::from([integral.clone()]);
        let basis = bases.atoms[basis_id as usize].clone();
        let previous = target.by_basis.insert(basis, inputs);
        if previous.is_some() {
            return Ok(Err(Declined("two spectator classes share a channel and basis")));
        }
    }
    stats.folded += folded.len();
    stats.finish_seconds += phase.elapsed().as_secs_f64();
    Ok(Ok(folded))
}

/// The momentum rows `(a, b)` over `(k_1..k_L, q)` of canonical scalar
/// product `index` (the packing order of [`Reader`]).
fn canonical_rows(index: usize, loops: usize) -> (Vec<i64>, Vec<i64>) {
    let unit = |position: usize| {
        let mut row = vec![0; loops + 1];
        row[position] = 1;
        row
    };
    let mut position = 0;
    for i in 0..loops {
        for j in i..loops {
            if position == index {
                return (unit(i), unit(j));
            }
            position += 1;
        }
    }
    (unit(index - position), unit(loops))
}

/// `Σ c·Π form_d^n_d` of a polynomial in the canonical scalar products
/// `0..vars` (packed monomials), as a polynomial in the slots (packed slot
/// exponents): Horner in the highest variable, recursively.
fn substitute(
    terms: &mut [(u128, Rational)],
    vars: usize,
    forms: &[rustred::propagator::presentation::AffineSlotForm],
) -> HashMap<u128, Rational> {
    let mut result = HashMap::<u128, Rational>::new();
    if terms.is_empty() {
        return result;
    }
    // The highest variable present.
    let exponent = |packed: u128, var: usize| ((packed >> (BITS * var as u32)) & ((1 << BITS) - 1)) as u32;
    let Some(var) = (0..vars).rev().find(|var| terms.iter().any(|(m, _)| exponent(*m, *var) > 0)) else {
        let mut total = Rational::zero();
        for (_, value) in terms.iter() {
            total += value;
        }
        if !total.is_zero() {
            result.insert(0, total);
        }
        return result;
    };
    terms.sort_unstable_by_key(|(m, _)| std::cmp::Reverse(exponent(*m, var)));
    let top = exponent(terms[0].0, var);
    let mask = !(((1u128 << BITS) - 1) << (BITS * var as u32));
    let mut start = 0;
    for power in (0..=top).rev() {
        if !result.is_empty() {
            result = multiply_form(&result, &forms[var]);
        }
        let mut end = start;
        while end < terms.len() && exponent(terms[end].0, var) == power {
            end += 1;
        }
        let slice = &mut terms[start..end];
        for (monomial, _) in slice.iter_mut() {
            *monomial &= mask;
        }
        let inner = substitute(slice, var, forms);
        for (monomial, value) in inner {
            match result.get_mut(&monomial) {
                Some(sum) => *sum += &value,
                None => {
                    result.insert(monomial, value);
                }
            }
        }
        start = end;
    }
    result.retain(|_, value| !value.is_zero());
    result
}

/// `polynomial · form` on packed slot exponents.
fn multiply_form(
    polynomial: &HashMap<u128, Rational>,
    form: &rustred::propagator::presentation::AffineSlotForm,
) -> HashMap<u128, Rational> {
    let mut next = HashMap::<u128, Rational>::with_capacity(polynomial.len() * 2);
    for (monomial, value) in polynomial {
        if !form.constant.is_zero() {
            *next.entry(*monomial).or_insert_with(Rational::zero) += &(value * &form.constant);
        }
        for (slot, coefficient) in &form.terms {
            let shifted = monomial + (1u128 << (BITS * u32::from(*slot)));
            *next.entry(shifted).or_insert_with(Rational::zero) += &(value * coefficient);
        }
    }
    next
}

/// A loop-free spectator monomial after the numerator rewrites of the
/// ordinary path (scalar admission, the matched substitutions, dot notation
/// and expansion).
fn spectator_image(
    settings: &VakintSettings,
    rules: &ReplacementRules,
    head: &Atom,
    spectator: &Atom,
) -> Result<Result<Atom, Declined>, VakintError> {
    // The matched substitutions only rewrite input loop momenta `k(j)`, which a
    // spectator does not contain, so the image depends on the spectator and
    // the head alone and is shared by every closure of the thread.
    thread_local! {
        static IMAGES: std::cell::RefCell<HashMap<(Atom, Atom), Atom>> =
            std::cell::RefCell::new(HashMap::new());
    }
    if spectator.contains_symbol(S.k) {
        return Ok(Err(Declined("a spectator with a loop momentum")));
    }
    let key = (spectator.clone(), head.clone());
    if let Some(image) = IMAGES.with(|images| images.borrow().get(&key).cloned()) {
        return Ok(Ok(image));
    }
    let mut term = VakintTerm {
        integral: Atom::num(1),
        numerator: crate::p_integral::scalar_numerator(spectator.as_view(), head)?,
        vectors: Vec::new(),
    };
    term.apply_numerator_replacement_rules(rules, settings)?;
    let image = Vakint::convert_to_dot_notation(term.numerator.as_view()).expand();
    if image.contains_symbol(S.k) {
        return Ok(Err(Declined("a spectator acquired a loop momentum")));
    }
    IMAGES.with(|images| {
        let mut images = images.borrow_mut();
        if images.len() > 1 << 16 {
            images.clear();
        }
        images.insert(key, image.clone());
    });
    Ok(Ok(image))
}
