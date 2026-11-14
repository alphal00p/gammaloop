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

use std::collections::{BTreeMap, BTreeSet, HashMap};
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

fn add_into(target: &mut Poly, key: Key, value: Rational) {
    use std::collections::hash_map::Entry;
    match target.entry(key) {
        Entry::Occupied(mut entry) => {
            *entry.get_mut() += &value;
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

/// The numerator reader of one closure.
struct Reader<'a> {
    /// Vacuum (or cut) loop momentum `k(i)` -> its vector over the canonical
    /// `k(1..L)` and the head.
    momenta: &'a HashMap<Atom, Vec<Rational>>,
    loops: usize,
    sigma: i64,
    /// The number of canonical scalar products `k(i)·k(j)` (i <= j), `k(i)·Q`.
    count: usize,
    spectators: &'a mut Spectators,
    cache: HashMap<Atom, Rc<Poly>>,
    terms: usize,
}

impl Reader<'_> {
    fn pair_index(&self, i: usize, j: usize) -> usize {
        let (i, j) = (i.min(j), i.max(j));
        // Upper triangle by rows: (0,0..L-1), (1,1..L-1), ...
        i * self.loops - i * i.saturating_sub(1) / 2 + (j - i)
    }

    fn unit() -> Poly {
        Poly::from([((0, 0), Rational::one())])
    }

    fn multiply(&mut self, left: &Poly, right: &Poly) -> Result<Poly, Declined> {
        let mut product = Poly::with_capacity(left.len().max(right.len()));
        for ((ls, lp), lv) in left {
            for ((rs, rp), rv) in right {
                let packed = lp + rp;
                if degree(packed, self.count) > 60 || (packed >> J_SHIFT) > 250 {
                    return Err(Declined("a numerator degree beyond the packed range"));
                }
                let (spectator, negate) = self.spectators.multiply(*ls, *rs);
                let mut value = lv * rv;
                if negate {
                    value = -value;
                }
                add_into(&mut product, (spectator, packed), value);
            }
        }
        Ok(product)
    }

    fn power(&mut self, base: &Poly, exponent: i64) -> Result<Poly, Declined> {
        let mut result = Self::unit();
        for _ in 0..exponent {
            result = self.multiply(&result, base)?;
        }
        Ok(result)
    }

    /// A linear momentum (a loop momentum `k(i)` or a sum of them) as its
    /// vector.
    fn vector(&self, momentum: AtomView) -> Result<Vec<Rational>, Declined> {
        if let Some(vector) = self.momenta.get(&momentum.to_owned()) {
            return Ok(vector.clone());
        }
        let mut total = vec![Rational::zero(); self.loops + 1];
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
                            (value, *head)
                        }
                        _ => return Err(Declined("a nonlinear momentum")),
                    }
                }
                other => (Rational::one(), other),
            };
            let vector = self
                .momenta
                .get(&head.to_owned())
                .ok_or(Declined("an unknown momentum in a scalar product"))?;
            for (slot, value) in total.iter_mut().zip(vector) {
                *slot += &(&scale * value);
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
    fn dot(&mut self, a: AtomView, b: AtomView) -> Result<Poly, Declined> {
        let (va, vb) = (self.vector(a)?, self.vector(b)?);
        let l = self.loops;
        let mut poly = Poly::new();
        for i in 0..l {
            for j in 0..l {
                let value = &va[i] * &vb[j];
                if !value.is_zero() {
                    let index = self.pair_index(i, j);
                    add_into(&mut poly, (0, 1u128 << (BITS * index as u32)), value);
                }
            }
        }
        let external = l * (l + 1) / 2;
        let sigma = Rational::from(self.sigma);
        for i in 0..l {
            let value = &(&(&va[i] * &vb[l]) + &(&va[l] * &vb[i])) * &sigma;
            if !value.is_zero() {
                add_into(&mut poly, (0, 1u128 << (BITS * (external + i) as u32)), value);
            }
        }
        let value = &va[l] * &vb[l];
        if !value.is_zero() {
            add_into(&mut poly, (0, 1u128 << J_SHIFT), value);
        }
        Ok(poly)
    }

    /// Whether a factor depends on a loop momentum.
    fn has_loop(view: AtomView) -> bool {
        view.contains_symbol(S.k)
    }

    fn cached(&mut self, view: AtomView) -> Result<Rc<Poly>, Declined> {
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
    fn read(&mut self, view: AtomView) -> Result<Poly, Declined> {
        match view {
            AtomView::Add(sum) => {
                let mut total = Poly::new();
                for term in sum.iter() {
                    self.terms += 1;
                    let poly = self.cached_term(term)?;
                    for (key, value) in poly.iter() {
                        add_into(&mut total, *key, value.clone());
                    }
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

    fn cached_term(&mut self, term: AtomView) -> Result<Rc<Poly>, Declined> {
        match term {
            AtomView::Mul(_) => Ok(Rc::new(self.term(term)?)),
            other => self.cached(other),
        }
    }

    /// One summand: numbers, spectators and loop-dependent factors.
    fn term(&mut self, term: AtomView) -> Result<Poly, Declined> {
        let factors = match term {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            other => vec![other],
        };
        let mut value = Rational::one();
        let mut imaginary = false;
        let mut spectators = Vec::new();
        let mut loops = Vec::new();
        for factor in factors {
            match factor {
                AtomView::Num(_) => {
                    let (number, phase) = number(factor)?;
                    value *= &number;
                    if phase {
                        if imaginary {
                            value = -value;
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
        let mut product = Poly::from([((spectator, 0), value)]);
        for poly in loops {
            product = self.multiply(&product, &poly)?;
        }
        Ok(product)
    }

    fn dot_factor(&mut self, view: AtomView) -> Result<Poly, Declined> {
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

/// The folded channels of one closure, lowered directly; `Err(Declined)` for
/// the ordinary path.
#[allow(clippy::too_many_arguments)]
pub(crate) fn lower_closure(
    _vakint: &Vakint,
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
    let mut spectators = Spectators::new();
    let polynomial = {
        let mut reader = Reader {
            momenta: &vacuum_momenta,
            loops,
            sigma: cut.sigma,
            count,
            spectators: &mut spectators,
            cache: HashMap::new(),
            terms: 0,
        };
        let numerator_poly = match reader.read(numerator) {
            Ok(poly) => poly,
            Err(declined) => return Ok(Err(declined)),
        };
        stats.numerator_terms += reader.terms;
        if cut.factor == Atom::num(1) {
            numerator_poly
        } else {
            let mut factor_reader = Reader {
                momenta: &cut_momenta,
                loops,
                sigma: cut.sigma,
                count,
                spectators: reader.spectators,
                cache: HashMap::new(),
                terms: 0,
            };
            let factor = match factor_reader.read(cut.factor.as_view()) {
                Ok(poly) => poly,
                Err(declined) => return Ok(Err(declined)),
            };
            match factor_reader.multiply(&numerator_poly, &factor) {
                Ok(poly) => poly,
                Err(declined) => return Ok(Err(declined)),
            }
        }
    };
    stats.canonical_terms += polynomial.len();

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
    // Lowered monomials as (key id, value); keys numbered once per closure,
    // buckets dense over key ids.
    let mut lowered_cache = HashMap::<u128, Rc<Vec<(usize, Rational)>>>::new();
    let mut key_ids = HashMap::<Vec<i64>, usize>::new();
    let mut keys = Vec::<Vec<i64>>::new();
    let imaginary_unit = Atom::i();
    // (monomial id, phase, scale power, basis id) -> bucket.
    let mut bucket_ids = HashMap::<(u32, bool, i64, u32), usize>::new();
    let mut buckets = Vec::<Vec<Rational>>::new();
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
        let lowered = match lowered_cache.get(&loop_monomial) {
            Some(lowered) => Rc::clone(lowered),
            None => {
                let lowered = lower_monomial(loop_monomial, count, &forms, &base, top.arity)
                    .into_iter()
                    .map(|(key, value)| {
                        let next = keys.len();
                        let id = *key_ids.entry(key.clone()).or_insert_with(|| {
                            keys.push(key);
                            next
                        });
                        (id, value)
                    })
                    .collect::<Vec<_>>();
                let lowered = Rc::new(lowered);
                lowered_cache.insert(loop_monomial, Rc::clone(&lowered));
                lowered
            }
        };
        let next = buckets.len();
        let bucket = *bucket_ids
            .entry((monomial_id, phase, scale_power, basis_id))
            .or_insert_with(|| {
                buckets.push(Vec::new());
                next
            });
        let bucket = &mut buckets[bucket];
        if bucket.len() < keys.len() {
            bucket.resize(keys.len(), Rational::zero());
        }
        let unit = value == Rational::one();
        for (id, coefficient) in lowered.iter() {
            if unit {
                bucket[*id] += coefficient;
            } else {
                bucket[*id] += &(&value * coefficient);
            }
        }
        stats.lowered += lowered.len();
    }
    stats.distinct_monomials += lowered_cache.len();
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

/// `Π form_d^n_d` of a packed loop monomial as (slot key, value): the key is
/// the line powers lowered by the slot exponents (ISPs at negative powers),
/// as [`super::p_integral`]'s `lower_part`.
fn lower_monomial(
    packed: u128,
    count: usize,
    forms: &[rustred::propagator::presentation::AffineSlotForm],
    base: &[i64],
    arity: usize,
) -> Vec<(Vec<i64>, Rational)> {
    let mut product = HashMap::<u128, Rational>::from([(0u128, Rational::one())]);
    for index in 0..count {
        let exponent = (packed >> (BITS * index as u32)) & ((1 << BITS) - 1);
        let form = &forms[index];
        for _ in 0..exponent {
            let mut next = HashMap::<u128, Rational>::with_capacity(product.len() * 2);
            for (monomial, value) in &product {
                if !form.constant.is_zero() {
                    *next.entry(*monomial).or_insert_with(Rational::zero) += &(value * &form.constant);
                }
                for (slot, coefficient) in &form.terms {
                    let shifted = monomial + (1u128 << (BITS * u32::from(*slot)));
                    *next.entry(shifted).or_insert_with(Rational::zero) += &(value * coefficient);
                }
            }
            next.retain(|_, value| !value.is_zero());
            product = next;
        }
    }
    let mut lowered = product
        .into_iter()
        .map(|(monomial, value)| {
            let key = (0..arity)
                .map(|slot| {
                    base[slot] - ((monomial >> (BITS * slot as u32)) & ((1 << BITS) - 1)) as i64
                })
                .collect::<Vec<_>>();
            (key, value)
        })
        .collect::<Vec<_>>();
    lowered.sort_by(|left, right| left.0.cmp(&right.0));
    lowered
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
