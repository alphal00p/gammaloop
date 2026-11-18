//! The kernel of the contracting lowering: exact polynomials in `D`, packed
//! slot monomials with their zero-sector test, the slot forms of momentum
//! pairs, and open tensors (sums of index structures with slot-polynomial
//! coefficients) with their contraction.
//!
//! Representation (helper VH, numerator assignment, lane engine):
//!
//! * [`DPoly`]: a dense polynomial in the dimension `D` with exact rational
//!   coefficients (index = power of `D`).
//! * [`Mono`]: one monomial of the running polynomial: the slot exponents
//!   packed six bits per slot in a `u128` (the packing of
//!   [`super::super::p_lowering`]), `r` = the number of loop-dependent
//!   canonical scalar products it was built from (fixes the scale power of
//!   its channel) and `j` = the power of `head·head` (kept symbolic, as the
//!   reference reader keeps it).
//! * [`Poly`]: `Mono → DPoly`, sparse.
//! * [`Tensor`]: `Structure → Poly`, where a structure is the sorted list of
//!   the open legs of a term: metrics `g(i, j)` between open Lorentz indices
//!   and momentum legs `v^i` (a canonical momentum vector carrying an open
//!   index). Scalar products are never kept: the moment two momenta are
//!   contracted their pair is replaced by its slot form ([`Kinematics::pair_form`]).
//!
//! Contraction rules ([`resolve`]): an index that occurs twice joins the two
//! legs that carry it: `g(x,i) g(i,y) → g(x,y)`, `g(x,i) v^i → v^x`,
//! `v^i w^i → (v·w)` (a slot form), `g(i,i) → D`.
//!
//! Early zero-sector pruning ([`Slots::is_zero_sector`]): numerator factors
//! only raise slot exponents, so a monomial whose cancelled lines already make
//! its sector zero (RustRed's census, lane prune's oracle, closed under
//! cancelling lines) has only zero-sector descendants: it is dropped the
//! moment it is formed.

use std::rc::Rc;
use std::sync::Arc;

use ahash::AHashMap as HashMap;
use symbolica::domains::rational::Rational;

use super::super::p_zero::ZeroSectors;

/// Bits per slot exponent of a packed slot monomial.
pub(crate) const BITS: u32 = 6;
/// Exponents stay at most `2^(BITS-1) - 1`: sums of two never carry, and the
/// field-wise cancellation test has a free top bit.
pub(crate) const MAX_EXPONENT: u32 = (1 << (BITS - 1)) - 1;

/// A polynomial in `D` with exact rational coefficients, dense, no trailing
/// zeros (empty = 0).
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct DPoly(pub(crate) Vec<Rational>);

impl DPoly {
    pub(crate) fn zero() -> Self {
        Self(Vec::new())
    }

    pub(crate) fn constant(value: Rational) -> Self {
        if value.is_zero() {
            Self::zero()
        } else {
            Self(vec![value])
        }
    }

    pub(crate) fn one() -> Self {
        Self(vec![Rational::one()])
    }

    /// `D^power`.
    pub(crate) fn d_power(power: usize) -> Self {
        let mut coefficients = vec![Rational::zero(); power + 1];
        coefficients[power] = Rational::one();
        Self(coefficients)
    }

    pub(crate) fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    fn trim(&mut self) {
        while self.0.last().is_some_and(Rational::is_zero) {
            self.0.pop();
        }
    }

    pub(crate) fn add_assign(&mut self, other: &DPoly) {
        if self.0.len() < other.0.len() {
            self.0.resize(other.0.len(), Rational::zero());
        }
        for (target, value) in self.0.iter_mut().zip(&other.0) {
            *target += value;
        }
        self.trim();
    }

    pub(crate) fn mul(&self, other: &DPoly) -> DPoly {
        if self.is_zero() || other.is_zero() {
            return DPoly::zero();
        }
        if other.0.len() == 1 {
            return self.scale(&other.0[0]);
        }
        if self.0.len() == 1 {
            return other.scale(&self.0[0]);
        }
        let mut out = vec![Rational::zero(); self.0.len() + other.0.len() - 1];
        for (i, a) in self.0.iter().enumerate() {
            if a.is_zero() {
                continue;
            }
            for (j, b) in other.0.iter().enumerate() {
                out[i + j] += &(a * b);
            }
        }
        let mut product = DPoly(out);
        product.trim();
        product
    }

    pub(crate) fn scale(&self, factor: &Rational) -> DPoly {
        if factor.is_zero() {
            return DPoly::zero();
        }
        DPoly(self.0.iter().map(|value| value * factor).collect())
    }

    /// `self · D^power`.
    pub(crate) fn shift(&self, power: u32) -> DPoly {
        if power == 0 || self.is_zero() {
            return self.clone();
        }
        let mut out = vec![Rational::zero(); power as usize];
        out.extend(self.0.iter().cloned());
        DPoly(out)
    }
}

/// One monomial of the running polynomial (module documentation).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) struct Mono {
    pub(crate) slots: u128,
    pub(crate) r: u16,
    pub(crate) j: u16,
}

pub(crate) type Poly = HashMap<Mono, DPoly>;

/// The scalar polynomial `c` (constant monomial).
pub(crate) fn poly_constant(value: DPoly) -> Poly {
    let mut poly = Poly::default();
    if !value.is_zero() {
        poly.insert(Mono::default(), value);
    }
    poly
}

/// `target += poly`.
pub(crate) fn poly_add_into(target: &mut Poly, poly: Poly) {
    if target.is_empty() {
        *target = poly;
        return;
    }
    for (mono, value) in poly {
        add_term(target, mono, value);
    }
}

fn add_term(target: &mut Poly, mono: Mono, value: DPoly) {
    use std::collections::hash_map::Entry;
    match target.entry(mono) {
        Entry::Occupied(mut entry) => {
            entry.get_mut().add_assign(&value);
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

/// Counters of the kernel (statistics of one lowering).
#[derive(Clone, Debug, Default)]
pub(crate) struct KernelCounters {
    /// Candidate monomials formed by a product (before the zero test).
    pub(crate) formed: u64,
    /// Candidate monomials dropped because their sector is already zero.
    pub(crate) dropped: u64,
}

/// The slot layout of the matched presentation and its zero sectors.
pub(crate) struct Slots {
    pub(crate) arity: usize,
    pub(crate) lines: usize,
    /// Line powers of the request (0 for a pinched line), 0 for ISPs.
    pub(crate) base: Vec<i64>,
    /// Bit `BITS - 1` of every field.
    high: u128,
    /// Per line field `2^(BITS-1) - base`, 0 for ISP fields: `exp + cancel`
    /// sets the field's top bit iff `exp >= base` (the line is cancelled).
    cancel: u128,
    zero: Option<Arc<ZeroSectors>>,
    pub(crate) prune: bool,
}

impl Slots {
    pub(crate) fn new(
        arity: usize,
        lines: usize,
        base: Vec<i64>,
        zero: Option<Arc<ZeroSectors>>,
        prune: bool,
    ) -> Result<Self, String> {
        if arity as u32 * BITS > 128 {
            return Err(format!(
                "a presentation of {arity} slots beyond the packed range"
            ));
        }
        let mut high = 0u128;
        let mut cancel = 0u128;
        for (slot, power) in base.iter().enumerate() {
            let shift = BITS * slot as u32;
            high |= 1u128 << (shift + BITS - 1);
            if slot < lines {
                if !(0..=i64::from(MAX_EXPONENT) + 1).contains(power) {
                    return Err(format!("a line power {power} beyond the packed range"));
                }
                cancel |= ((1u128 << (BITS - 1)) - *power as u128) << shift;
            }
        }
        if let Some(zero) = &zero
            && zero.lines != lines
        {
            return Err("the zero-sector table has another line count".into());
        }
        Ok(Self {
            arity,
            lines,
            base,
            high,
            cancel,
            zero,
            prune,
        })
    }

    /// `a + b` of two packed monomials, `None` when a field would exceed
    /// [`MAX_EXPONENT`].
    #[inline]
    pub(crate) fn add(&self, a: u128, b: u128) -> Option<u128> {
        let sum = a + b;
        (sum & self.high == 0).then_some(sum)
    }

    /// The support of a packed monomial: bit `s` for every line slot whose
    /// remaining power `base[s] - exp[s]` is positive.
    #[inline]
    pub(crate) fn support(&self, slots: u128) -> u64 {
        let cancelled = (slots + self.cancel) & self.high;
        let mut support = 0u64;
        for line in 0..self.lines {
            if (cancelled >> (BITS * line as u32 + BITS - 1)) & 1 == 0 {
                support |= 1 << line;
            }
        }
        support
    }

    /// Whether the sector of a packed monomial (and of all its descendants)
    /// is zero: only when pruning is on and the oracle is known.
    #[inline]
    pub(crate) fn is_zero_sector(&self, slots: u128) -> bool {
        match &self.zero {
            Some(zero) if self.prune => zero.is_zero(self.support(slots)),
            _ => false,
        }
    }

    /// The exponent of `slot` in a packed monomial.
    #[inline]
    pub(crate) fn exponent(slots: u128, slot: usize) -> i64 {
        ((slots >> (BITS * slot as u32)) & ((1 << BITS) - 1)) as i64
    }

    /// The integral key of a packed monomial: `base - exponent` per slot.
    pub(crate) fn key(&self, slots: u128) -> Vec<i64> {
        (0..self.arity)
            .map(|slot| self.base[slot] - Self::exponent(slots, slot))
            .collect()
    }
}

/// One slot form `c + Σ r_s D_s` of a canonical scalar product (RustRed's
/// `AffineSlotForm`).
#[derive(Clone, Debug)]
pub(crate) struct SlotForm {
    pub(crate) constant: Rational,
    pub(crate) terms: Vec<(u16, Rational)>,
}

/// The canonical momenta of the matched presentation, the slot forms of the
/// canonical scalar products, and the interned momentum vectors with the
/// memo of their pair forms.
pub(crate) struct Kinematics {
    pub(crate) loops: usize,
    /// `external = sigma · head`: the sign of `k(i)·Q` entries.
    pub(crate) sigma: i64,
    /// Forms of `k_i·k_j` (i <= j, upper triangle by rows) then `k_i·Q`.
    pub(crate) forms: Vec<SlotForm>,
    /// Interned canonical vectors over `(k_1..k_L, head)`.
    pub(crate) vectors: Vec<Vec<Rational>>,
    ids: HashMap<Vec<Rational>, u16>,
    pairs: HashMap<(u16, u16), Rc<Vec<(Mono, Rational)>>>,
}

impl Kinematics {
    pub(crate) fn new(loops: usize, sigma: i64, forms: Vec<SlotForm>) -> Self {
        Self {
            loops,
            sigma,
            forms,
            vectors: Vec::new(),
            ids: HashMap::default(),
            pairs: HashMap::default(),
        }
    }

    /// The id of a vector (`None` for the zero vector).
    pub(crate) fn intern(&mut self, vector: Vec<Rational>) -> Result<Option<u16>, String> {
        if vector.iter().all(Rational::is_zero) {
            return Ok(None);
        }
        if let Some(id) = self.ids.get(&vector) {
            return Ok(Some(*id));
        }
        let id = u16::try_from(self.vectors.len()).map_err(|_| "too many momentum vectors")?;
        self.vectors.push(vector.clone());
        self.ids.insert(vector, id);
        Ok(Some(id))
    }

    fn pair_index(&self, i: usize, j: usize) -> usize {
        let (i, j) = (i.min(j), i.max(j));
        i * self.loops - i * i.saturating_sub(1) / 2 + (j - i)
    }

    /// The slot form of `v_a · v_b`: the loop-dependent part as slot forms of
    /// the canonical products (`r = 1`) and the `head·head` part as `j = 1`
    /// (`r = 0`), exactly the split of the reference reader
    /// (`p_lowering::Reader::dot` followed by the slot substitution).
    pub(crate) fn pair_form(&mut self, a: u16, b: u16) -> Rc<Vec<(Mono, Rational)>> {
        let key = (a.min(b), a.max(b));
        if let Some(form) = self.pairs.get(&key) {
            return Rc::clone(form);
        }
        let (va, vb) = (&self.vectors[a as usize], &self.vectors[b as usize]);
        let l = self.loops;
        let mut terms = HashMap::<Mono, Rational>::default();
        let add_form = |terms: &mut HashMap<Mono, Rational>, form: &SlotForm, weight: &Rational| {
            if weight.is_zero() {
                return;
            }
            let mut push = |mono: Mono, value: Rational| {
                let entry = terms.entry(mono).or_insert_with(Rational::zero);
                *entry += &value;
            };
            if !form.constant.is_zero() {
                push(
                    Mono {
                        slots: 0,
                        r: 1,
                        j: 0,
                    },
                    &form.constant * weight,
                );
            }
            for (slot, coefficient) in &form.terms {
                push(
                    Mono {
                        slots: 1u128 << (BITS * u32::from(*slot)),
                        r: 1,
                        j: 0,
                    },
                    coefficient * weight,
                );
            }
        };
        for i in 0..l {
            for j in 0..l {
                let weight = &va[i] * &vb[j];
                let index = self.pair_index(i, j);
                add_form(&mut terms, &self.forms[index], &weight);
            }
        }
        let external = l * (l + 1) / 2;
        let sigma = Rational::from(self.sigma);
        for i in 0..l {
            let weight = &(&(&va[i] * &vb[l]) + &(&va[l] * &vb[i])) * &sigma;
            add_form(&mut terms, &self.forms[external + i], &weight);
        }
        let head = &va[l] * &vb[l];
        if !head.is_zero() {
            let entry = terms
                .entry(Mono {
                    slots: 0,
                    r: 0,
                    j: 1,
                })
                .or_insert_with(Rational::zero);
            *entry += &head;
        }
        let mut form = terms
            .into_iter()
            .filter(|(_, value)| !value.is_zero())
            .collect::<Vec<_>>();
        form.sort_by(|a, b| a.0.cmp(&b.0));
        let form = Rc::new(form);
        self.pairs.insert(key, Rc::clone(&form));
        form
    }
}

/// `a · b` of two polynomials, dropping zero-sector monomials.
pub(crate) fn poly_mul(
    a: &Poly,
    b: &Poly,
    slots: &Slots,
    counters: &mut KernelCounters,
) -> Result<Poly, String> {
    // A constant side is a scale.
    for (constant_side, other) in [(a, b), (b, a)] {
        if constant_side.len() == 1
            && let Some(value) = constant_side.get(&Mono::default())
        {
            let mut out = Poly::with_capacity(other.len());
            for (mono, coefficient) in other {
                let product = coefficient.mul(value);
                if !product.is_zero() {
                    out.insert(*mono, product);
                }
            }
            return Ok(out);
        }
    }
    let mut out = Poly::with_capacity(a.len().max(b.len()) * 2);
    for (ma, ca) in a {
        for (mb, cb) in b {
            counters.formed += 1;
            let packed = slots
                .add(ma.slots, mb.slots)
                .ok_or("a slot exponent beyond the packed range")?;
            if slots.is_zero_sector(packed) {
                counters.dropped += 1;
                continue;
            }
            let mono = Mono {
                slots: packed,
                r: ma.r + mb.r,
                j: ma.j + mb.j,
            };
            add_term(&mut out, mono, ca.mul(cb));
        }
    }
    Ok(out)
}

/// `a · form` of a polynomial and a pair form, dropping zero-sector
/// monomials.
pub(crate) fn poly_mul_form(
    a: &Poly,
    form: &[(Mono, Rational)],
    slots: &Slots,
    counters: &mut KernelCounters,
) -> Result<Poly, String> {
    let mut out = Poly::with_capacity(a.len() * form.len().max(1));
    for (ma, ca) in a {
        for (mf, cf) in form {
            counters.formed += 1;
            let packed = slots
                .add(ma.slots, mf.slots)
                .ok_or("a slot exponent beyond the packed range")?;
            if mf.slots != 0 && slots.is_zero_sector(packed) {
                counters.dropped += 1;
                continue;
            }
            let mono = Mono {
                slots: packed,
                r: ma.r + mf.r,
                j: ma.j + mf.j,
            };
            add_term(&mut out, mono, ca.scale(cf));
        }
    }
    Ok(out)
}

/// One end of a leg: a momentum vector (interned id) or an open Lorentz
/// index (interned id).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) enum End {
    Mom(u16),
    Idx(u32),
}

/// A leg: a metric `(Idx, Idx)` or a momentum leg `(Mom, Idx)`, normalised
/// `(min, max)`.
pub(crate) type Leg = (End, End);

pub(crate) fn leg(a: End, b: End) -> Leg {
    if a <= b { (a, b) } else { (b, a) }
}

/// The open legs of a term, sorted.
pub(crate) type Structure = Vec<Leg>;

/// An open tensor: structures with slot-polynomial coefficients.
#[derive(Clone, Debug, Default)]
pub(crate) struct Tensor {
    pub(crate) terms: HashMap<Structure, Poly>,
}

impl Tensor {
    /// The scalar 1.
    pub(crate) fn one() -> Self {
        Self::scalar(DPoly::one())
    }

    pub(crate) fn scalar(value: DPoly) -> Self {
        let mut terms = HashMap::default();
        if !value.is_zero() {
            terms.insert(Vec::new(), poly_constant(value));
        }
        Self { terms }
    }

    /// One structure with coefficient 1.
    pub(crate) fn structure(structure: Structure) -> Self {
        let mut terms = HashMap::default();
        terms.insert(structure, poly_constant(DPoly::one()));
        Self { terms }
    }

    pub(crate) fn add(&mut self, structure: Structure, poly: Poly) {
        if poly.is_empty() {
            return;
        }
        match self.terms.get_mut(&structure) {
            Some(sum) => {
                poly_add_into(sum, poly);
                if sum.is_empty() {
                    self.terms.remove(&structure);
                }
            }
            None => {
                self.terms.insert(structure, poly);
            }
        }
    }

    pub(crate) fn add_tensor(&mut self, other: Tensor) {
        for (structure, poly) in other.terms {
            self.add(structure, poly);
        }
    }

    /// Number of (structure, monomial) terms.
    pub(crate) fn size(&self) -> usize {
        self.terms.values().map(|poly| poly.len()).sum()
    }

    /// The open indices of the tensor (union over its structures).
    pub(crate) fn open_indices(&self) -> Vec<u32> {
        let mut open = Vec::new();
        for structure in self.terms.keys() {
            for (a, b) in structure {
                for end in [a, b] {
                    if let End::Idx(index) = end {
                        open.push(*index);
                    }
                }
            }
        }
        open.sort_unstable();
        open.dedup();
        open
    }
}

/// A structure after its repeated indices are contracted.
pub(crate) struct Resolved {
    pub(crate) structure: Structure,
    /// Momentum pairs that met (scalar products to multiply in).
    pub(crate) products: Vec<(u16, u16)>,
    /// Closed metric loops (`g(i,i) = D`).
    pub(crate) d_power: u32,
}

/// Contract every index that occurs twice in `legs` (module documentation).
pub(crate) fn resolve(mut legs: Vec<Leg>) -> Result<Resolved, String> {
    let mut products = Vec::new();
    let mut d_power = 0;
    loop {
        // Find an index with two occurrences.
        let mut found = None;
        'search: for p in 0..legs.len() {
            for side in 0..2 {
                let end = if side == 0 { legs[p].0 } else { legs[p].1 };
                let End::Idx(index) = end else { continue };
                // The same leg's other side.
                if side == 0 && legs[p].1 == end {
                    found = Some((p, 0, p, 1));
                    break 'search;
                }
                for q in p + 1..legs.len() {
                    for other_side in 0..2 {
                        let other = if other_side == 0 {
                            legs[q].0
                        } else {
                            legs[q].1
                        };
                        if other == End::Idx(index) {
                            found = Some((p, side, q, other_side));
                            break 'search;
                        }
                    }
                }
            }
        }
        let Some((p, sp, q, sq)) = found else { break };
        if p == q {
            legs.swap_remove(p);
            d_power += 1;
            continue;
        }
        let x = if sp == 0 { legs[p].1 } else { legs[p].0 };
        let y = if sq == 0 { legs[q].1 } else { legs[q].0 };
        // Remove the higher position first.
        legs.swap_remove(q.max(p));
        legs.swap_remove(q.min(p));
        match (x, y) {
            (End::Mom(a), End::Mom(b)) => products.push((a.min(b), a.max(b))),
            _ => legs.push(leg(x, y)),
        }
    }
    // Every index must now occur once.
    let mut seen = legs
        .iter()
        .flat_map(|(a, b)| [*a, *b])
        .filter(|end| matches!(end, End::Idx(_)))
        .collect::<Vec<_>>();
    seen.sort_unstable();
    if seen.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err("an index occurs more than twice".into());
    }
    legs.sort_unstable();
    Ok(Resolved {
        structure: legs,
        products,
        d_power,
    })
}

/// `a · b` of two tensors: every pair of structures is joined and its
/// repeated indices contracted; the coefficient is `poly_a · poly_b · Π pair
/// forms · D^loops`, with early zero-sector pruning.
pub(crate) fn contract(
    a: &Tensor,
    b: &Tensor,
    kinematics: &mut Kinematics,
    slots: &Slots,
    counters: &mut KernelCounters,
) -> Result<Tensor, String> {
    let mut out = Tensor::default();
    for (sa, pa) in &a.terms {
        for (sb, pb) in &b.terms {
            let mut legs = Vec::with_capacity(sa.len() + sb.len());
            legs.extend_from_slice(sa);
            legs.extend_from_slice(sb);
            let resolved = resolve(legs)?;
            let mut poly = poly_mul(pa, pb, slots, counters)?;
            for (u, v) in &resolved.products {
                if poly.is_empty() {
                    break;
                }
                let form = kinematics.pair_form(*u, *v);
                poly = poly_mul_form(&poly, &form, slots, counters)?;
            }
            if resolved.d_power > 0 {
                for value in poly.values_mut() {
                    *value = value.shift(resolved.d_power);
                }
            }
            out.add(resolved.structure, poly);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(n: i64) -> Rational {
        Rational::from(n)
    }

    #[test]
    fn d_polynomials_multiply_and_cancel() {
        let a = DPoly(vec![r(-1), r(1)]); // D - 1
        let b = DPoly(vec![r(1), r(1)]); // D + 1
        assert_eq!(a.mul(&b), DPoly(vec![r(-1), r(0), r(1)]));
        let mut c = a.clone();
        c.add_assign(&DPoly(vec![r(1), r(-1)]));
        assert!(c.is_zero());
        assert_eq!(DPoly::one().shift(2), DPoly::d_power(2));
    }

    #[test]
    fn index_chains_resolve_to_metrics_products_and_traces() {
        // g(1,2) g(2,3) v^3 w^1 -> v.w
        let legs = vec![
            leg(End::Idx(1), End::Idx(2)),
            leg(End::Idx(2), End::Idx(3)),
            leg(End::Mom(0), End::Idx(3)),
            leg(End::Mom(1), End::Idx(1)),
        ];
        let resolved = resolve(legs).unwrap();
        assert!(resolved.structure.is_empty());
        assert_eq!(resolved.products, vec![(0, 1)]);
        // g(1,2) g(2,1) -> D
        let resolved = resolve(vec![
            leg(End::Idx(1), End::Idx(2)),
            leg(End::Idx(2), End::Idx(1)),
        ])
        .unwrap();
        assert!(resolved.structure.is_empty() && resolved.products.is_empty());
        assert_eq!(resolved.d_power, 1);
        // g(1,2) v^2 with 1 open -> v^1
        let resolved = resolve(vec![
            leg(End::Idx(1), End::Idx(2)),
            leg(End::Mom(4), End::Idx(2)),
        ])
        .unwrap();
        assert_eq!(resolved.structure, vec![leg(End::Mom(4), End::Idx(1))]);
    }

    #[test]
    fn the_cancellation_test_is_field_wise() {
        // Three lines with powers 1, 2, 0 and one ISP.
        let slots = Slots::new(4, 3, vec![1, 2, 0, 0], None, true).unwrap();
        let pack = |e: [u128; 4]| e[0] | e[1] << BITS | e[2] << (2 * BITS) | e[3] << (3 * BITS);
        assert_eq!(slots.support(pack([0, 0, 0, 0])), 0b011);
        assert_eq!(slots.support(pack([1, 0, 0, 5])), 0b010);
        assert_eq!(slots.support(pack([0, 1, 0, 0])), 0b011);
        assert_eq!(slots.support(pack([0, 2, 0, 0])), 0b001);
        assert_eq!(slots.support(pack([31, 31, 31, 31])), 0);
        assert_eq!(slots.key(pack([1, 1, 0, 2])), vec![0, 1, 0, -2]);
        assert!(slots.add(pack([31, 0, 0, 0]), pack([1, 0, 0, 0])).is_none());
        assert!(slots.add(pack([30, 0, 0, 0]), pack([1, 0, 0, 0])).is_some());
    }
}
