//! The numerator split shared by the p-integral backends.
//!
//! A backend request (`docs/p_integrals.md`, section 5) is one canonical
//! family, its line powers and a scalar numerator in the canonical `k(j)` and
//! the user head `p(i)`. [`split_parts`] turns it into homogeneous parts:
//! negative powers join the numerator, every loop-momentum monomial is written
//! in the canonical scalar products `k(i)·k(j)` and `k(i)·Q` (the user head is
//! `σ·Q`), and the monomials are grouped by degree `r` and by the factor of
//! their coefficient that is not a rational number. Each [`PPart`] retains
//! its own scale power, Wick sign and non-rational coefficient (`1+ε`,
//! `mass²`). Forcer sets `Q·Q = 1`, so parts of different degree must remain
//! distinguishable: independent scalar tags let them share one expression
//! and reduction call. The tags are extracted before each part's scale,
//! sign and coefficient are restored.
//!
//! `None` means the sector is scaleless (an exact zero), decided here once for
//! every backend: the scaleless rule of `docs/p_integrals.md` (section 1,
//! item 6) on the graph with the pinched lines contracted, or remaining lines
//! that do not span the loop momenta.
//!
//! The Wick sign `(-1)^(Σa+r)` of a part is applied on the Forcer path only
//! (M2 plan, conventions C3/C4): Forcer reduces Euclidean integrals, while the
//! RustRed reduction already returns Minkowski-signed values.

use std::collections::hash_map::Entry;
use std::collections::{BTreeMap, HashMap};
use std::rc::Rc;
use std::sync::OnceLock;

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::coefficient::{Coefficient, CoefficientView};
use symbolica::function;

use super::{PFamilyData, PropagatorRequest, external_head};
use crate::symbols::S;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{Vakint, VakintError, get_integer_from_atom, get_prop_with_id};

/// One canonical line: endpoints and the integer coefficients of
/// `k(1..=L)` and of the canonical `p(1)` in its momentum.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct CanonicalLine {
    pub from: i64,
    pub to: i64,
    pub loops: Vec<i64>,
    pub external: i64,
}

impl CanonicalLine {
    /// The momentum with the canonical `p(1)` written as `external`.
    pub(crate) fn momentum(&self, external: &Atom) -> Atom {
        let mut momentum = external * Atom::num(self.external);
        for (index, coefficient) in self.loops.iter().enumerate() {
            momentum += function!(S.k, Atom::num(index as i64 + 1)) * Atom::num(*coefficient);
        }
        momentum
    }
}

/// The canonical lines of a registry family, in line order.
pub(crate) fn canonical_lines(family: &PFamilyData) -> Result<Vec<CanonicalLine>, VakintError> {
    // Cache only the immutable registry definitions. A caller-provided family
    // with the same name but different lines must still be validated below.
    static LINES: [OnceLock<Vec<CanonicalLine>>; super::P_FAMILIES.len()] =
        [const { OnceLock::new() }; super::P_FAMILIES.len()];
    let Some(index) = super::P_FAMILIES.iter().position(|registered| {
        registered.name == family.name
            && registered.loops == family.loops
            && registered.lines == family.lines
    }) else {
        return canonical_lines_uncached(family);
    };
    if let Some(lines) = LINES[index].get() {
        return Ok(lines.clone());
    }
    let lines = canonical_lines_uncached(family)?;
    Ok(LINES[index].get_or_init(|| lines).clone())
}

fn canonical_lines_uncached(family: &PFamilyData) -> Result<Vec<CanonicalLine>, VakintError> {
    let invalid = |detail: String| {
        VakintError::InvalidGenericExpression(format!("{}: {detail}", family.name))
    };
    Vakint::initialize_vakint_symbols();
    let expression = vk_parse!(family.lines).map_err(|error| invalid(error.to_string()))?;
    let heads = (1..=family.loops)
        .map(|i| function!(S.k, Atom::num(i as i64)))
        .chain([function!(S.p, Atom::num(1))])
        .collect::<Vec<_>>();
    (1..=family.line_count())
        .map(|id| {
            let prop = get_prop_with_id(expression.as_view(), id)
                .ok_or_else(|| invalid(format!("missing line {id}")))?;
            let node = |name: &str| {
                get_integer_from_atom(prop[&vk_symbol!(name)].as_view())
                    .ok_or_else(|| invalid(format!("line {id} has a non-integer vertex")))
            };
            let momentum = &prop[&vk_symbol!("q_")];
            let coefficients = heads
                .iter()
                .map(|head| {
                    let mut value = momentum.clone();
                    for other in &heads {
                        let image = if other == head {
                            Atom::num(1)
                        } else {
                            Atom::Zero
                        };
                        value = value.replace(other.to_pattern()).with(image.to_pattern());
                    }
                    get_integer_from_atom(value.expand().as_view())
                        .ok_or_else(|| invalid(format!("line {id} is not an integer combination")))
                })
                .collect::<Result<Vec<_>, _>>()?;
            Ok(CanonicalLine {
                from: node("nl_")?,
                to: node("nr_")?,
                external: coefficients[family.loops],
                loops: coefficients[..family.loops].to_vec(),
            })
        })
        .collect()
}

/// Exact rational numbers for the small linear systems of a routing and the
/// rational coefficients of a numerator.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Fraction(pub i128, pub i128);

impl Fraction {
    pub(crate) fn new(numerator: i128, denominator: i128) -> Self {
        fn gcd(a: i128, b: i128) -> i128 {
            if b == 0 { a.abs() } else { gcd(b, a % b) }
        }
        let divisor = gcd(numerator, denominator).max(1) * denominator.signum();
        Fraction(numerator / divisor, denominator / divisor)
    }
    pub(crate) fn is_zero(self) -> bool {
        self.0 == 0
    }
    pub(crate) fn sub(self, other: Self) -> Self {
        Self::new(self.0 * other.1 - other.0 * self.1, self.1 * other.1)
    }
    pub(crate) fn mul(self, other: Self) -> Self {
        Self::new(self.0 * other.0, self.1 * other.1)
    }
    pub(crate) fn div(self, other: Self) -> Self {
        Self::new(self.0 * other.1, self.1 * other.0)
    }
}

/// The loop momenta in terms of the remaining (positive-power) lines and the
/// external momentum: row `j` holds the coefficients of `p1..pN` and of `Q`
/// in `k(j+1)`. `None` when the remaining lines do not span the loop
/// momenta: some loop carries no propagator and the integral is scaleless.
pub(crate) fn loop_momenta_in_lines(
    lines: &[&CanonicalLine],
    loops: usize,
) -> Option<Vec<Vec<Fraction>>> {
    // Gauss-Jordan elimination on [A | I] with A the loop coefficients of
    // the remaining lines (rows), keeping the row operations to read off the
    // left inverse X with X A = 1; then k = X (p - b Q).
    let rows = lines.len();
    let mut matrix = lines
        .iter()
        .enumerate()
        .map(|(row, line)| {
            let mut entries = line
                .loops
                .iter()
                .map(|&value| Fraction::new(value as i128, 1))
                .collect::<Vec<_>>();
            entries.extend((0..rows).map(|column| Fraction::new((column == row) as i128, 1)));
            entries
        })
        .collect::<Vec<_>>();
    let mut pivot_rows = Vec::new();
    for column in 0..loops {
        let pivot = (pivot_rows.len()..rows).find(|&row| !matrix[row][column].is_zero())?;
        let target = pivot_rows.len();
        matrix.swap(target, pivot);
        let pivot_value = matrix[target][column];
        for entry in matrix[target].iter_mut() {
            *entry = entry.div(pivot_value);
        }
        for row in 0..rows {
            let factor = matrix[row][column];
            if row != target && !factor.is_zero() {
                let pivot_row = matrix[target].clone();
                for (entry, pivot_entry) in matrix[row].iter_mut().zip(pivot_row) {
                    *entry = entry.sub(factor.mul(pivot_entry));
                }
            }
        }
        pivot_rows.push(target);
    }
    Some(
        (0..loops)
            .map(|j| {
                let inverse = &matrix[j][loops..];
                let mut row = inverse.to_vec();
                let shift = inverse
                    .iter()
                    .zip(lines)
                    .fold(Fraction::new(0, 1), |sum, (x, line)| {
                        sum.sub(x.mul(Fraction::new(-(line.external as i128), 1)))
                    });
                // k = X p - (X b) Q, with `shift` = X b.
                row.push(Fraction::new(0, 1).sub(shift));
                row
            })
            .collect(),
    )
}

/// The contracted graph of a sector: union-find representatives of the
/// vertices once all pinched lines are contracted. It decides scalelessness
/// only; Forcer itself always gets the uncontracted family graph.
pub(crate) struct ContractedGraph {
    representative: BTreeMap<i64, i64>,
}

impl ContractedGraph {
    pub(crate) fn new(lines: &[CanonicalLine], pinched: &[bool]) -> Self {
        let mut representative = BTreeMap::new();
        for line in lines {
            representative.insert(line.from, line.from);
            representative.insert(line.to, line.to);
        }
        let mut graph = ContractedGraph { representative };
        for (line, &pinched) in lines.iter().zip(pinched) {
            if pinched {
                let (a, b) = (graph.find(line.from), graph.find(line.to));
                graph.representative.insert(a.max(b), a.min(b));
            }
        }
        graph
    }

    pub(crate) fn find(&self, mut vertex: i64) -> i64 {
        while self.representative[&vertex] != vertex {
            vertex = self.representative[&vertex];
        }
        vertex
    }

    /// The scaleless rule of `docs/p_integrals.md` (section 1, item 6) on the
    /// remaining lines: the external vertices coincide, a self-loop remains,
    /// or removing one vertex isolates a part without external vertex.
    pub(crate) fn is_scaleless(&self, remaining: &[&CanonicalLine]) -> bool {
        let (q_in, q_out) = (self.find(1), self.find(2));
        if q_in == q_out {
            return true;
        }
        let edges = remaining
            .iter()
            .map(|line| (self.find(line.from), self.find(line.to)))
            .collect::<Vec<_>>();
        if edges.iter().any(|(a, b)| a == b) {
            return true;
        }
        let mut vertices = edges.iter().flat_map(|&(a, b)| [a, b]).collect::<Vec<_>>();
        vertices.sort();
        vertices.dedup();
        vertices.iter().any(|&removed| {
            // Components of the graph without `removed`: each must reach an
            // external vertex (unless it is `removed`'s own external vertex).
            let mut component = BTreeMap::new();
            for &start in vertices.iter().filter(|&&v| v != removed) {
                if component.contains_key(&start) {
                    continue;
                }
                let mut stack = vec![start];
                component.insert(start, start);
                while let Some(vertex) = stack.pop() {
                    for &(a, b) in &edges {
                        for (x, y) in [(a, b), (b, a)] {
                            if x == vertex && y != removed && !component.contains_key(&y) {
                                component.insert(y, start);
                                stack.push(y);
                            }
                        }
                    }
                }
            }
            let mut roots = component.values().copied().collect::<Vec<_>>();
            roots.sort();
            roots.dedup();
            roots.iter().any(|root| {
                component
                    .iter()
                    .filter(|(_, r)| *r == root)
                    .all(|(&vertex, _)| vertex != q_in && vertex != q_out)
            })
        })
    }
}

/// One canonical scalar product of a numerator monomial.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum PDot {
    /// `k(i)·k(j)`, 1-based, in the order written.
    Loops(usize, usize),
    /// `k(i)·Q`, 1-based, with `Q` the canonical `p(1)`.
    External(usize),
}

/// A monomial in canonical scalar products: `(dot, exponent > 0)` factors in
/// the order of the numerator; empty for `1`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct PMonomial(pub Vec<(PDot, i64)>);

impl PMonomial {
    /// The degree in scalar products: the power of `Q²` the monomial carries.
    pub(crate) fn degree(&self) -> i64 {
        self.0.iter().map(|(_, exponent)| exponent).sum()
    }
}

/// One homogeneous part of a request.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct PPart {
    /// `w = 2L - Σa + r`: the value carries `(Q²)^w`.
    pub scale_power: i64,
    /// `(-1)^(Σa+r)`, applied on the Forcer path only (module documentation).
    pub wick_sign: i64,
    /// The factor of the coefficients that is not a rational number (`1` if
    /// there is none), left outside the reduction.
    pub coefficient: Atom,
    /// `Σ c·monomial` with rational `c` (the sign of `σ` included), every
    /// monomial of degree `r`, in the order of the numerator.
    pub polynomial: Vec<(Fraction, PMonomial)>,
}

/// A request split into parts, with the sector data Forcer renders.
#[derive(Clone, Debug)]
pub(crate) struct PSplit {
    /// The canonical lines of the family.
    pub lines: Vec<CanonicalLine>,
    /// `(line index, power)` of the lines of positive power, in line order.
    pub remaining: Vec<(usize, i64)>,
    /// [`loop_momenta_in_lines`] of the remaining lines.
    pub loop_momenta: Vec<Vec<Fraction>>,
    pub parts: Vec<PPart>,
}

pub(crate) fn tensor_numerator(detail: impl std::fmt::Display) -> VakintError {
    VakintError::InvalidNumerator(format!(
        "p-integral: tensor p-integrals are not supported yet: {detail}"
    ))
}

/// A real rational coefficient as an exact fraction, if it is one.
pub(crate) fn rational(value: AtomView) -> Option<Fraction> {
    let AtomView::Num(number) = value else {
        return None;
    };
    match number.get_coeff_view() {
        CoefficientView::Natural(numerator, denominator, 0, _) => {
            Some(Fraction::new(numerator as i128, denominator as i128))
        }
        _ => None,
    }
}

/// Split a numerator coefficient into an exact rational and the rest.
pub(crate) fn split_rational(coefficient: &Atom) -> (Fraction, Atom) {
    if let Some(value) = rational(coefficient.as_view()) {
        return (value, Atom::num(1));
    }
    if let AtomView::Mul(product) = coefficient.as_view() {
        let mut value = Fraction::new(1, 1);
        let mut rest = Atom::num(1);
        for factor in product.iter() {
            match rational(factor) {
                Some(number) => value = value.mul(number),
                None => rest *= factor.to_owned(),
            }
        }
        return (value, rest);
    }
    (Fraction::new(1, 1), coefficient.clone())
}

/// A loop index `j` of `k(j)` with `1 <= j <= loops`.
fn loop_index(momentum: AtomView, loops: usize) -> Option<usize> {
    let AtomView::Fun(function) = momentum else {
        return None;
    };
    if function.get_symbol() != S.k || function.get_nargs() != 1 {
        return None;
    }
    let index = get_integer_from_atom(function.iter().next().unwrap())?;
    (1..=loops as i64)
        .contains(&index)
        .then_some(index as usize)
}

/// The canonical scalar product of one loop-dependent `dot(...)` and its
/// sign: the user's head `p(i)` is `σ·Q`.
pub(crate) fn scalar_product(
    dot: AtomView,
    head: &Atom,
    sigma: i64,
    loops: usize,
) -> Result<(i64, PDot), VakintError> {
    let AtomView::Fun(function) = dot else {
        return Err(VakintError::InvalidNumerator(format!(
            "p-integral: unexpected numerator factor {dot}"
        )));
    };
    let arguments = function.iter().collect::<Vec<_>>();
    let [a, b] = arguments.as_slice() else {
        return Err(tensor_numerator(dot));
    };
    let (loop_argument, other) = match (loop_index(*a, loops), loop_index(*b, loops)) {
        (Some(i), Some(j)) => return Ok((1, PDot::Loops(i, j))),
        (Some(i), None) => (i, *b),
        (None, Some(j)) => (j, *a),
        (None, None) => {
            return Err(VakintError::InvalidNumerator(format!(
                "p-integral: unknown loop momentum in {dot}"
            )));
        }
    };
    if other == head.as_view() {
        Ok((sigma, PDot::External(loop_argument)))
    } else {
        Err(tensor_numerator(format!(
            "{dot} contracts a loop momentum with a vector other than the external momentum {head}"
        )))
    }
}

/// The sign and canonical form of a monomial in loop-dependent dot products.
fn canonical_monomial(
    monomial: AtomView,
    head: &Atom,
    sigma: i64,
    loops: usize,
) -> Result<(i64, PMonomial), VakintError> {
    let factors = match monomial {
        AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
        _ => vec![monomial],
    };
    let mut sign = 1;
    let mut dots = Vec::new();
    for factor in factors {
        let (base, exponent) = match factor {
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                (base, get_integer_from_atom(exponent))
            }
            AtomView::Num(_) if rational(factor) == Some(Fraction::new(1, 1)) => continue,
            _ => (factor, Some(1)),
        };
        let exponent = exponent.filter(|exponent| *exponent > 0).ok_or_else(|| {
            VakintError::InvalidNumerator(format!(
                "p-integral: the numerator must be polynomial in the loop momenta, found {factor}"
            ))
        })?;
        let (factor_sign, dot) = scalar_product(base, head, sigma, loops)?;
        if exponent % 2 == 1 {
            sign *= factor_sign;
        }
        dots.push((dot, exponent));
    }
    Ok((sign, PMonomial(dots)))
}

/// The directory of the opt-in `VAKINT_DUMP_SPLIT_INPUTS` diagnostic, read
/// once per process.
fn split_dump_directory() -> Option<&'static std::path::Path> {
    static DIRECTORY: std::sync::LazyLock<Option<std::path::PathBuf>> =
        std::sync::LazyLock::new(|| std::env::var_os("VAKINT_DUMP_SPLIT_INPUTS").map(Into::into));
    DIRECTORY.as_deref()
}

/// Write one [`split`] input and its time to the dump directory, as
/// `split-<pid>-<n>.txt` with one `key value` line per field; the numerator
/// and the external momentum are canonical strings (they parse back).
fn dump_split_input(
    directory: &std::path::Path,
    family: &PFamilyData,
    powers: &[i64],
    canonical_numerator: AtomView,
    external: &Atom,
    elapsed: std::time::Duration,
) {
    static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let index = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let powers = powers
        .iter()
        .map(i64::to_string)
        .collect::<Vec<_>>()
        .join(",");
    let text = format!(
        "seconds {:.6}\nfamily {}\npowers {powers}\nexternal {}\nnumerator {}\n",
        elapsed.as_secs_f64(),
        family.name,
        external.to_canonical_string(),
        canonical_numerator.to_canonical_string(),
    );
    let path = directory.join(format!("split-{}-{index:07}.txt", std::process::id()));
    let _ = std::fs::write(path, text);
}

/// The split of `topo(NAME(external, powers...))` with a scalar numerator in
/// the canonical `k(j)` and the user's head; `None` when the sector is
/// scaleless (an exact zero).
///
/// With `VAKINT_DUMP_SPLIT_INPUTS=<dir>` every input is also written to
/// `<dir>` (a diagnostic corpus for the split tests; off by default).
pub(crate) fn split(
    family: &PFamilyData,
    powers: &[i64],
    canonical_numerator: AtomView,
    external: &Atom,
) -> Result<Option<PSplit>, VakintError> {
    let Some(directory) = split_dump_directory() else {
        return split_with(
            family,
            powers,
            canonical_numerator,
            external,
            numerator_parts,
        );
    };
    let started = std::time::Instant::now();
    let result = split_with(
        family,
        powers,
        canonical_numerator,
        external,
        numerator_parts,
    );
    dump_split_input(
        directory,
        family,
        powers,
        canonical_numerator,
        external,
        started.elapsed(),
    );
    result
}

/// [`split`] with the parts of [`expanded_parts`] only: the reference of the
/// grouped split.
#[cfg(test)]
pub(crate) fn split_reference(
    family: &PFamilyData,
    powers: &[i64],
    canonical_numerator: AtomView,
    external: &Atom,
) -> Result<Option<PSplit>, VakintError> {
    split_with(
        family,
        powers,
        canonical_numerator,
        external,
        expanded_parts,
    )
}

/// How the parts of a numerator are computed from it, its negative-power
/// factors and the labels of its sector.
type NumeratorParts = fn(AtomView, &[Atom], &PartCollector) -> Result<Vec<PPart>, VakintError>;

/// [`split`] with the parts of `parts`.
fn split_with(
    family: &PFamilyData,
    powers: &[i64],
    canonical_numerator: AtomView,
    external: &Atom,
    parts: NumeratorParts,
) -> Result<Option<PSplit>, VakintError> {
    let lines = canonical_lines(family)?;
    if powers.len() != lines.len() {
        return Err(VakintError::InvalidGenericExpression(format!(
            "{}: {} powers for {} lines",
            family.name,
            powers.len(),
            lines.len()
        )));
    }
    let (head, sigma) = external_head(external)?;
    let loops = family.loops;

    // A negative power is a numerator factor in the user's momenta.
    let negative_powers = lines
        .iter()
        .zip(powers)
        .filter(|(_, power)| **power < 0)
        .map(|(line, &power)| {
            let momentum = line.momentum(external);
            function!(S.dot, &momentum, &momentum).pow(Atom::num(-power))
        })
        .collect::<Vec<_>>();
    let pinched = powers.iter().map(|&power| power <= 0).collect::<Vec<_>>();
    let graph = ContractedGraph::new(&lines, &pinched);
    let remaining = powers
        .iter()
        .enumerate()
        .filter(|(_, power)| **power > 0)
        .map(|(index, power)| (index, *power))
        .collect::<Vec<_>>();
    let remaining_lines = remaining
        .iter()
        .map(|(index, _)| &lines[*index])
        .collect::<Vec<_>>();
    if graph.is_scaleless(&remaining_lines) {
        return Ok(None);
    }
    let Some(loop_momenta) = loop_momenta_in_lines(&remaining_lines, loops) else {
        return Ok(None);
    };
    let collector = PartCollector {
        head: &head,
        sigma,
        loops,
        total_power: remaining.iter().map(|(_, power)| power).sum(),
    };
    let parts = parts(canonical_numerator, &negative_powers, &collector)?;
    Ok(Some(PSplit {
        lines,
        remaining,
        loop_momenta,
        parts,
    }))
}

/// How the parts of one sector are labelled: the user's head `σ·Q`, the loop
/// count and the total power `Σa` of the remaining lines.
struct PartCollector<'a> {
    head: &'a Atom,
    sigma: i64,
    loops: usize,
    total_power: i64,
}

impl PartCollector<'_> {
    /// The loop-dependent scalar products of `expressions`, sorted, each
    /// checked to be a canonical scalar product.
    fn loop_dots<'b>(
        &self,
        expressions: impl IntoIterator<Item = AtomView<'b>>,
    ) -> Result<Vec<Atom>, VakintError> {
        let mut loop_dots = Vec::new();
        let mut error = None;
        for expression in expressions {
            expression.visitor(&mut |part| {
                if error.is_some() {
                    return false;
                }
                if let AtomView::Fun(function) = part
                    && function.get_symbol() == S.dot
                    && part.contains_symbol(S.k)
                {
                    match scalar_product(part, self.head, self.sigma, self.loops) {
                        Ok(_) => loop_dots.push(part.to_owned()),
                        Err(e) => error = Some(e),
                    }
                    return false;
                }
                true
            });
            if let Some(error) = error {
                return Err(error);
            }
        }
        loop_dots.sort();
        loop_dots.dedup();
        Ok(loop_dots)
    }

    /// Homogeneous parts of `Σ coefficient·monomial`: loop-momentum
    /// monomials, grouped by degree and by the non-rational factor of their
    /// coefficient.
    fn parts(
        &self,
        monomials: impl IntoIterator<Item = (Atom, Atom)>,
    ) -> Result<Vec<PPart>, VakintError> {
        let mut groups = Vec::<((i64, Atom), Vec<(Fraction, PMonomial)>)>::new();
        let mut group_indices = HashMap::<(i64, Atom), usize>::new();
        for (monomial, coefficient) in monomials {
            if coefficient.contains_symbol(S.k) {
                return Err(tensor_numerator(format!(
                    "the loop momenta must appear in scalar products only, found {coefficient}"
                )));
            }
            let (monomial_sign, canonical) =
                canonical_monomial(monomial.as_view(), self.head, self.sigma, self.loops)?;
            let (value, rest) = split_rational(&coefficient);
            let term = (
                Fraction::new(monomial_sign as i128, 1).mul(value),
                canonical,
            );
            let key = (term.1.degree(), rest);
            // Keep first-seen part order and every polynomial's existing order;
            // only the lookup uses hashing, never the emitted sequence.
            if let Some(&index) = group_indices.get(&key) {
                groups[index].1.push(term);
            } else {
                group_indices.insert(key.clone(), groups.len());
                groups.push((key, vec![term]));
            }
        }
        Ok(groups
            .into_iter()
            .map(|((degree, coefficient), polynomial)| {
                let exponent = self.total_power + degree;
                PPart {
                    scale_power: 2 * self.loops as i64 - self.total_power + degree,
                    wick_sign: if exponent % 2 == 0 { 1 } else { -1 },
                    coefficient,
                    polynomial,
                }
            })
            .collect())
    }
}

/// The parts of [`grouped_parts`] where it applies, else those of
/// [`expanded_parts`]: the same parts either way.
fn numerator_parts(
    canonical_numerator: AtomView,
    negative_powers: &[Atom],
    collector: &PartCollector,
) -> Result<Vec<PPart>, VakintError> {
    match grouped_parts(canonical_numerator, negative_powers, collector) {
        Some(parts) => Ok(parts),
        None => expanded_parts(canonical_numerator, negative_powers, collector),
    }
}

/// The original parts: expand the whole numerator with its negative-power
/// factors, then collect it by loop-momentum monomial.
fn expanded_parts(
    canonical_numerator: AtomView,
    negative_powers: &[Atom],
    collector: &PartCollector,
) -> Result<Vec<PPart>, VakintError> {
    let mut numerator = canonical_numerator.to_owned();
    for factor in negative_powers {
        numerator *= factor.clone();
    }
    let numerator = Vakint::convert_to_dot_notation(numerator.as_view()).expand();
    let loop_dots = collector.loop_dots([numerator.as_view()])?;
    collector.parts(numerator.coefficient_list::<i8>(&loop_dots))
}

/// The parts of [`expanded_parts`], without expanding the whole numerator:
/// the numerator is converted to a polynomial in the loop dots
/// ([`LoopPolynomials`]), with each distinct loop factor and product of loop
/// factors converted once and the loop-free factors of a sum's terms summed
/// and expanded once per loop factor.
///
/// Every coefficient is the canonical (normalized, merged) sum of the same
/// products of terms that the expansion in [`expanded_parts`] forms, and the
/// monomials are listed in the lexicographic exponent order of
/// `coefficient_list`, so the parts are those of [`expanded_parts`] when
/// - all numbers are exact, every power is a real rational power of a
///   symbol or function or an integer power of a sum or product (no power
///   of a number), and there is no `abs`: then normalized products are
///   associative, multiply powers as the expansion does and never form a
///   sum factor ([`admits_loop_polynomials`]);
/// - every coefficient term expands to itself: the second expansion in
///   [`expanded_parts`] changes nothing;
/// - no dot-notation pass applies to a coefficient;
/// - every loop dot is a canonical scalar product, the loop momenta appear in
///   loop dots to positive integer powers only, and the parts can be
///   collected (the errors of [`expanded_parts`]).
///
/// `None` (use [`expanded_parts`]) when one of these does not hold.
fn grouped_parts(
    canonical_numerator: AtomView,
    negative_powers: &[Atom],
    collector: &PartCollector,
) -> Option<Vec<PPart>> {
    grouped_parts_with(
        canonical_numerator,
        negative_powers,
        collector,
        PENDING_PRODUCT_BYTES,
    )
}

/// [`grouped_parts`], merging the pending products of a polynomial operation
/// once they exceed `pending_bound` bytes ([`ProductSums`]; `usize::MAX`
/// merges every monomial once, at the end).
fn grouped_parts_with(
    canonical_numerator: AtomView,
    negative_powers: &[Atom],
    collector: &PartCollector,
    pending_bound: usize,
) -> Option<Vec<PPart>> {
    if !admits_loop_polynomials(canonical_numerator) {
        return None;
    }
    let loop_dots = collector
        .loop_dots(
            std::iter::once(canonical_numerator).chain(negative_powers.iter().map(Atom::as_view)),
        )
        .ok()?;
    let atom_polynomial = || -> Option<LoopPolynomial> {
        let mut polynomials = LoopPolynomials {
            loop_dots: &loop_dots,
            factors: HashMap::new(),
            products: HashMap::new(),
            pending_bound,
        };
        let mut polynomial = polynomials.sum(sum_terms(canonical_numerator))?;
        for factor in negative_powers {
            let factor = polynomials.factor(factor.as_view())?;
            polynomial = polynomials.mul(&polynomial, &factor)?;
        }
        Some(polynomial)
    };
    let numeric = if *SPLIT_NUMERIC_LEGACY {
        None
    } else {
        numeric_loop_polynomial(&loop_dots, canonical_numerator, negative_powers)
    };
    let polynomial = match numeric {
        Some(polynomial) => {
            if *SPLIT_NUMERIC_CHECK {
                let reference = atom_polynomial();
                eprintln!(
                    "SPLIT_CHECK stage=split_numeric bytes={} equal={}",
                    canonical_numerator.get_byte_size(),
                    reference.as_ref() == Some(&polynomial)
                );
            }
            polynomial
        }
        None => atom_polynomial()?,
    };
    let mut monomials = polynomial.into_iter().collect::<Vec<_>>();
    monomials.sort_by(|(left, _), (right, _)| left.cmp(right));
    let mut fixed_factors = Vec::new();
    let mut coefficients = Vec::with_capacity(monomials.len());
    for (exponents, coefficient) in monomials {
        if !sum_terms(coefficient.as_view())
            .into_iter()
            .all(|term| expands_to_itself(term, &mut fixed_factors))
            || Vakint::dot_notation_passes_may_apply(coefficient.as_view())
        {
            return None;
        }
        // The monomial as `coefficient_list` writes it.
        let mut monomial = Atom::num(1);
        for (&exponent, dot) in exponents.iter().zip(&loop_dots) {
            if exponent != 0 {
                monomial = monomial * dot.pow(Atom::num(i64::from(exponent)));
            }
        }
        coefficients.push((monomial, coefficient));
    }
    collector.parts(coefficients).ok()
}

/// Whether all numbers of `expression` are exact and every power has a
/// real rational exponent on a symbol or function, or an integer one on a
/// sum or product. Normalized products of its terms then merge powers of
/// equal bases by adding the exponents, `x^a*x^b = x^(a+b)` (`x^0 = 1`,
/// `x^1 = x`), as the powers `t^n` of terms in the expansion multiply them,
/// `(x^a)^n = x^(n*a)`; this is associative and never forms a sum factor.
///
/// The two simplifications of a power that a product does not make are not
/// admitted: a power of a number (or of a power), since a product merges
/// `2^(1/3)*2^(1/3) = 2^(2/3)` but a power renormalizes the number,
/// `(2^(1/3))^2 = 4^(1/3)`; and `abs`, whose even powers of a real argument
/// a power rewrites, `abs(x)^2 = x^2`, but a product keeps,
/// `abs(x)*abs(x) = abs(x)^2`.
fn admits_loop_polynomials(expression: AtomView) -> bool {
    if !crate::utils::has_only_exact_numbers(expression) {
        return false;
    }
    let mut admitted = true;
    expression.visitor(&mut |part| {
        admitted &= match part {
            AtomView::Pow(power) => match power.get_base_exp() {
                (AtomView::Var(_) | AtomView::Fun(_), AtomView::Num(number)) => matches!(
                    number.get_coeff_view(),
                    CoefficientView::Natural(_, _, 0, _)
                ),
                (AtomView::Add(_) | AtomView::Mul(_), AtomView::Num(number)) => matches!(
                    number.get_coeff_view(),
                    CoefficientView::Natural(_, 1, 0, _)
                ),
                _ => false,
            },
            AtomView::Fun(function) => function.get_symbol() != Symbol::ABS,
            _ => true,
        };
        admitted
    });
    admitted
}

/// The exponents of a monomial over the sorted loop dots.
type Exponents = Vec<u8>;

/// A polynomial in the loop dots: nonzero coefficients that are expanded
/// loop-free sums.
type LoopPolynomial = HashMap<Exponents, Atom>;

/// The loop-dot polynomials of the subexpressions of one numerator, each
/// distinct loop factor (by its bytes) and product of loop factors converted
/// once.
struct LoopPolynomials<'a> {
    loop_dots: &'a [Atom],
    factors: HashMap<&'a [u8], Rc<LoopPolynomial>>,
    products: HashMap<Vec<&'a [u8]>, Rc<LoopPolynomial>>,
    /// The bound of [`ProductSums`].
    pending_bound: usize,
}

/// The products of one polynomial operation that [`ProductSums`] holds
/// unmerged, in bytes, before it merges them (at least the size of what it
/// holds already, so every merge at least doubles what it has merged).
const PENDING_PRODUCT_BYTES: usize = 64 << 20;

/// The products of a polynomial operation collected by monomial: each
/// coefficient is the one canonical sum of its products ([`Atom::add_many`];
/// a single product is kept as it is).
///
/// Once the pending products exceed the bound, every monomial with a merged
/// sum or more than one pending product is merged into its running sum, so
/// the unmerged products of a large operation are never all held at once.
/// Merging in batches gives the one merge of all products: the products are
/// products of terms of expanded sums with exact numbers
/// ([`admits_loop_polynomials`]), so no merge turns a term into a sum (that
/// needs a number times one sum, [`crate::utils::ExactSum`]) and the exact
/// coefficients of the terms with the same non-number part add up to the
/// same number in any grouping. A product that is a number times one sum
/// stops the batches: everything after it is merged once, at the end.
struct ProductSums {
    /// The merged sum and the pending products of each monomial.
    monomials: HashMap<Exponents, (Atom, Vec<Atom>)>,
    pending_bound: usize,
    pending_bytes: usize,
    held_bytes: usize,
    batches: bool,
}

impl ProductSums {
    fn new(pending_bound: usize) -> Self {
        Self {
            monomials: HashMap::new(),
            pending_bound,
            pending_bytes: 0,
            held_bytes: 0,
            batches: pending_bound != usize::MAX,
        }
    }

    fn push(&mut self, exponents: &Exponents, product: Atom) {
        if self.batches && crate::utils::is_scaled_sum(product.as_view()) {
            self.batches = false;
        }
        self.pending_bytes += product.as_view().get_byte_size();
        match self.monomials.get_mut(exponents) {
            Some((_, pending)) => pending.push(product),
            None => {
                self.monomials
                    .insert(exponents.clone(), (Atom::Zero, vec![product]));
            }
        }
        if self.batches && self.pending_bytes > self.pending_bound.max(self.held_bytes) {
            self.merge();
        }
    }

    /// Merge the pending products of every monomial that has a merged sum or
    /// more than one of them.
    fn merge(&mut self) {
        let mut held_bytes = 0;
        for (merged, pending) in self.monomials.values_mut() {
            if !merged.is_zero() || pending.len() > 1 {
                let sum = std::mem::take(merged);
                *merged = Atom::add_many(std::iter::once(sum).chain(pending.drain(..)));
            }
            held_bytes += merged.as_view().get_byte_size()
                + pending
                    .iter()
                    .map(|product| product.as_view().get_byte_size())
                    .sum::<usize>();
        }
        self.held_bytes = held_bytes;
        self.pending_bytes = 0;
    }

    fn finish(self) -> LoopPolynomial {
        let mut polynomial = HashMap::with_capacity(self.monomials.len());
        for (exponents, (merged, pending)) in self.monomials {
            let coefficient = if !merged.is_zero() {
                if pending.is_empty() {
                    merged
                } else {
                    Atom::add_many(std::iter::once(merged).chain(pending))
                }
            } else {
                match <[Atom; 1]>::try_from(pending) {
                    Ok([product]) => product,
                    Err(products) => Atom::add_many(products.iter()),
                }
            };
            if !coefficient.is_zero() {
                polynomial.insert(exponents, coefficient);
            }
        }
        polynomial
    }
}

impl<'a> LoopPolynomials<'a> {
    /// `dot^exponent` for a loop dot `dot`.
    fn monomial(&self, dot: AtomView, exponent: u8) -> Option<LoopPolynomial> {
        let index = self.loop_dots.binary_search(&dot.to_owned()).ok()?;
        let mut exponents = vec![0; self.loop_dots.len()];
        exponents[index] = exponent;
        Some(HashMap::from([(exponents, Atom::num(1))]))
    }

    fn mul(&self, left: &LoopPolynomial, right: &LoopPolynomial) -> Option<LoopPolynomial> {
        let mut products = ProductSums::new(self.pending_bound);
        for (left_exponents, left_coefficient) in left {
            let left_terms = sum_terms(left_coefficient.as_view());
            for (right_exponents, right_coefficient) in right {
                // `coefficient_list` stores the exponents as `i8`.
                let exponents = left_exponents
                    .iter()
                    .zip(right_exponents)
                    .map(|(left, right)| left.checked_add(*right).filter(|sum| *sum <= 127))
                    .collect::<Option<Exponents>>()?;
                let right_terms = sum_terms(right_coefficient.as_view());
                for &left in &left_terms {
                    for &right in &right_terms {
                        products.push(&exponents, left * right);
                    }
                }
            }
        }
        Some(products.finish())
    }

    /// The polynomial of a factor that holds a loop momentum: a loop dot, a
    /// sum, or a positive integer power of either.
    fn factor(&mut self, factor: AtomView<'a>) -> Option<Rc<LoopPolynomial>> {
        if let Some(polynomial) = self.factors.get(factor.get_data()) {
            return Some(polynomial.clone());
        }
        let polynomial = match factor {
            AtomView::Fun(function) if function.get_symbol() == S.dot => {
                self.monomial(factor, 1)?
            }
            AtomView::Add(_) | AtomView::Mul(_) => self.sum(sum_terms(factor))?,
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let exponent = match exponent {
                    AtomView::Num(number) => match number.get_coeff_view() {
                        CoefficientView::Natural(exponent, 1, 0, _) => u8::try_from(exponent)
                            .ok()
                            .filter(|exponent| (1..=127).contains(exponent))?,
                        _ => return None,
                    },
                    _ => return None,
                };
                match base {
                    AtomView::Fun(function) if function.get_symbol() == S.dot => {
                        self.monomial(base, exponent)?
                    }
                    AtomView::Add(_) => {
                        let base = self.factor(base)?;
                        let mut power = (*base).clone();
                        for _ in 1..exponent {
                            power = self.mul(&power, &base)?;
                        }
                        power
                    }
                    _ => return None,
                }
            }
            _ => return None,
        };
        let polynomial = Rc::new(polynomial);
        self.factors.insert(factor.get_data(), polynomial.clone());
        Some(polynomial)
    }

    /// The polynomial of a product of loop factors.
    fn product(&mut self, factors: &[AtomView<'a>]) -> Option<Rc<LoopPolynomial>> {
        let key = factors.iter().map(AtomView::get_data).collect::<Vec<_>>();
        if let Some(polynomial) = self.products.get(&key) {
            return Some(polynomial.clone());
        }
        let polynomial = match factors {
            [] => Rc::new(HashMap::from([(
                vec![0; self.loop_dots.len()],
                Atom::num(1),
            )])),
            [factor] => self.factor(*factor)?,
            [rest @ .., last] => {
                let rest = self.product(rest)?;
                let last = self.factor(*last)?;
                Rc::new(self.mul(&rest, &last)?)
            }
        };
        self.products.insert(key, polynomial.clone());
        Some(polynomial)
    }

    /// The polynomial of a sum: its terms grouped by their loop factors, the
    /// loop-free factors of each group summed and expanded once, and every
    /// coefficient one canonical sum of the products of the terms of those
    /// sums with the coefficient terms of the group's loop factors.
    fn sum(&mut self, terms: Vec<AtomView<'a>>) -> Option<LoopPolynomial> {
        let mut groups = Vec::<(Vec<AtomView<'a>>, Vec<Atom>)>::new();
        let mut group_indices = HashMap::<Vec<&'a [u8]>, usize>::new();
        for term in terms {
            let (loop_factors, spectators): (Vec<_>, Vec<_>) = product_factors(term)
                .into_iter()
                .partition(|factor| factor.contains_symbol(S.k));
            let spectator = if spectators.is_empty() {
                Atom::num(1)
            } else {
                Atom::mul_many(spectators)
            };
            let key = loop_factors.iter().map(AtomView::get_data).collect();
            match group_indices.entry(key) {
                Entry::Occupied(entry) => groups[*entry.get()].1.push(spectator),
                Entry::Vacant(entry) => {
                    entry.insert(groups.len());
                    groups.push((loop_factors, vec![spectator]));
                }
            }
        }
        let mut products = ProductSums::new(self.pending_bound);
        for (loop_factors, spectators) in groups {
            let spectators = Atom::add_many(spectators.iter()).expand();
            if spectators.is_zero() {
                continue;
            }
            let spectator_terms = sum_terms(spectators.as_view());
            for (exponents, coefficient) in self.product(&loop_factors)?.iter() {
                let coefficient_terms = sum_terms(coefficient.as_view());
                for &spectator in &spectator_terms {
                    for &term in &coefficient_terms {
                        products.push(exponents, spectator * term);
                    }
                }
            }
        }
        Some(products.finish())
    }
}

/// Whether the numeric loop polynomials ([`NumericLoopPolynomials`]) are
/// disabled (`VAKINT_ALGEBRA_LEGACY=split_numeric`) or also checked against
/// the atom loop polynomials (`VAKINT_ALGEBRA_CHECK=split_numeric`, reported
/// as `SPLIT_CHECK equal=...` on stderr).
fn split_numeric_switch(variable: &str) -> bool {
    std::env::var(variable)
        .map(|value| {
            value
                .split(',')
                .any(|entry| entry.trim() == "split_numeric" || entry.trim() == "all")
        })
        .unwrap_or(false)
}

static SPLIT_NUMERIC_LEGACY: std::sync::LazyLock<bool> =
    std::sync::LazyLock::new(|| split_numeric_switch("VAKINT_ALGEBRA_LEGACY"));
static SPLIT_NUMERIC_CHECK: std::sync::LazyLock<bool> =
    std::sync::LazyLock::new(|| split_numeric_switch("VAKINT_ALGEBRA_CHECK"));

/// A small exact rational `n/d` (`d > 0`, coprime): the numbers of the loop
/// polynomials of [`NumericLoopPolynomials`]. Arithmetic is checked; an
/// overflow makes the numeric split decline (`None`).
#[derive(Clone, Copy, Debug, PartialEq)]
struct SmallRational(i64, i64);

impl SmallRational {
    const ONE: Self = Self(1, 1);

    fn new(numerator: i128, denominator: i128) -> Option<Self> {
        fn gcd(mut a: i128, mut b: i128) -> i128 {
            while b != 0 {
                (a, b) = (b, a % b);
            }
            a.abs()
        }
        if denominator == 0 {
            return None;
        }
        if numerator == 0 {
            return Some(Self(0, 1));
        }
        let divisor = gcd(numerator, denominator);
        let (mut numerator, mut denominator) = (numerator / divisor, denominator / divisor);
        if denominator < 0 {
            (numerator, denominator) = (-numerator, -denominator);
        }
        Some(Self(
            i64::try_from(numerator).ok()?,
            i64::try_from(denominator).ok()?,
        ))
    }

    fn mul(self, other: Self) -> Option<Self> {
        Self::new(
            i128::from(self.0) * i128::from(other.0),
            i128::from(self.1) * i128::from(other.1),
        )
    }

    fn add(self, other: Self) -> Option<Self> {
        Self::new(
            i128::from(self.0) * i128::from(other.1) + i128::from(other.0) * i128::from(self.1),
            i128::from(self.1) * i128::from(other.1),
        )
    }

    fn is_zero(self) -> bool {
        self.0 == 0
    }

    /// A real rational number that fits, else `None`.
    fn from_view(number: CoefficientView) -> Option<Self> {
        match number {
            CoefficientView::Natural(numerator, denominator, 0, _) => {
                Self::new(i128::from(numerator), i128::from(denominator))
            }
            _ => None,
        }
    }

    fn from_coefficient(number: &Coefficient) -> Option<Self> {
        let atom = Atom::num(number.clone());
        match atom.as_view() {
            AtomView::Num(number) => Self::from_view(number.get_coeff_view()),
            _ => None,
        }
    }

    fn coefficient(self) -> Coefficient {
        Coefficient::from(symbolica::domains::rational::Rational::from((
            self.0, self.1,
        )))
    }
}

/// A loop-free expanded sum of a loop polynomial: small rationals per
/// distinct non-number monomial (an id of [`Monomials`]), sorted by id.
type SmallCoefficient = Vec<(u32, SmallRational)>;

/// A polynomial in the loop dots with [`SmallCoefficient`]s.
type SmallPolynomial = HashMap<Exponents, SmallCoefficient>;

/// The coefficients of the whole numerator: exact (complex) numbers per
/// monomial id, sorted by id.
type NumericCoefficient = Vec<(u32, Coefficient)>;

/// A polynomial in the loop dots with [`NumericCoefficient`]s.
type NumericPolynomial = HashMap<Exponents, NumericCoefficient>;

fn add_small(
    coefficient: &mut SmallCoefficient,
    monomial: u32,
    value: SmallRational,
) -> Option<()> {
    match coefficient.binary_search_by_key(&monomial, |(id, _)| *id) {
        Ok(index) => coefficient[index].1 = coefficient[index].1.add(value)?,
        Err(index) => coefficient.insert(index, (monomial, value)),
    }
    Some(())
}

fn add_numeric(coefficient: &mut NumericCoefficient, monomial: u32, value: Coefficient) {
    match coefficient.binary_search_by_key(&monomial, |(id, _)| *id) {
        Ok(index) => {
            let sum = std::mem::take(&mut coefficient[index].1) + value;
            coefficient[index].1 = sum;
        }
        Err(index) => coefficient.insert(index, (monomial, value)),
    }
}

fn clean_small(mut polynomial: SmallPolynomial) -> SmallPolynomial {
    polynomial.retain(|_, coefficient| {
        coefficient.retain(|(_, value)| !value.is_zero());
        coefficient.shrink_to_fit();
        !coefficient.is_empty()
    });
    polynomial.shrink_to_fit();
    polynomial
}

fn clean_numeric(mut polynomial: NumericPolynomial) -> NumericPolynomial {
    polynomial.retain(|_, coefficient| {
        coefficient.retain(|(_, value)| !value.is_zero());
        !coefficient.is_empty()
    });
    polynomial
}

/// The distinct loop-free monomials (normalized products without a number
/// factor; id 0 is `1`) and their products.
struct Monomials {
    atoms: Vec<Atom>,
    ids: HashMap<Atom, u32>,
    products: HashMap<(u32, u32), Option<(Coefficient, u32)>>,
    small_products: HashMap<(u32, u32), Option<(SmallRational, u32)>>,
}

impl Monomials {
    fn new() -> Self {
        Self {
            atoms: vec![Atom::num(1)],
            ids: HashMap::from([(Atom::num(1), 0)]),
            products: HashMap::new(),
            small_products: HashMap::new(),
        }
    }

    fn intern(&mut self, monomial: Atom) -> u32 {
        if let Some(&id) = self.ids.get(&monomial) {
            return id;
        }
        let id = self.atoms.len() as u32;
        self.atoms.push(monomial.clone());
        self.ids.insert(monomial, id);
        id
    }

    /// A normalized term as its exact number and the id of its monomial;
    /// `None` if the term holds an inexact number or its non-number part is a
    /// sum (a number times a sum, which a sum would not keep as a term).
    fn term(&mut self, term: AtomView) -> Option<(Coefficient, u32)> {
        let exact = |number: CoefficientView| -> Option<Coefficient> {
            match number {
                CoefficientView::Natural(..) | CoefficientView::Large(..) => {
                    Some(number.to_owned())
                }
                _ => None,
            }
        };
        match term {
            AtomView::Num(number) => Some((exact(number.get_coeff_view())?, 0)),
            AtomView::Add(_) => None,
            AtomView::Mul(product) => {
                let mut value = Coefficient::one();
                let mut others = Vec::new();
                for factor in product.iter() {
                    match factor {
                        AtomView::Num(number) => {
                            value = value * exact(number.get_coeff_view())?;
                        }
                        AtomView::Add(_) => return None,
                        other => others.push(other),
                    }
                }
                let monomial = match others[..] {
                    [] => Atom::num(1),
                    [single] => single.to_owned(),
                    _ => Atom::mul_many(others.iter().copied()),
                };
                Some((value, self.intern(monomial)))
            }
            other => Some((Coefficient::one(), self.intern(other.to_owned()))),
        }
    }

    /// [`Self::term`] with a small real rational number.
    fn small_term(&mut self, term: AtomView) -> Option<(SmallRational, u32)> {
        let (value, monomial) = self.term(term)?;
        Some((SmallRational::from_coefficient(&value)?, monomial))
    }

    /// The normalized product of two monomials, as a number and a monomial.
    fn product(&mut self, left: u32, right: u32) -> Option<(Coefficient, u32)> {
        if left == 0 {
            return Some((Coefficient::one(), right));
        }
        if right == 0 {
            return Some((Coefficient::one(), left));
        }
        let key = (left.min(right), left.max(right));
        if let Some(product) = self.products.get(&key) {
            return product.clone();
        }
        let atom = &self.atoms[key.0 as usize] * &self.atoms[key.1 as usize];
        let product = self.term(atom.as_view());
        self.products.insert(key, product.clone());
        product
    }

    /// [`Self::product`] with a small real rational number.
    fn small_product(&mut self, left: u32, right: u32) -> Option<(SmallRational, u32)> {
        if left == 0 {
            return Some((SmallRational::ONE, right));
        }
        if right == 0 {
            return Some((SmallRational::ONE, left));
        }
        let key = (left.min(right), left.max(right));
        if let Some(product) = self.small_products.get(&key) {
            return *product;
        }
        let product = self.product(left, right).and_then(|(value, monomial)| {
            Some((SmallRational::from_coefficient(&value)?, monomial))
        });
        self.small_products.insert(key, product);
        product
    }

    /// `Σ value·monomial` as the canonical sum of its terms.
    fn to_atom(&self, coefficient: &NumericCoefficient) -> Atom {
        let terms = coefficient
            .iter()
            .filter(|(_, value)| !value.is_zero())
            .map(|(monomial, value)| &Atom::num(value.clone()) * &self.atoms[*monomial as usize])
            .collect::<Vec<_>>();
        match <[Atom; 1]>::try_from(terms) {
            Ok([term]) => term,
            Err(terms) => Atom::add_many(terms.iter()),
        }
    }
}

/// [`LoopPolynomials`] with every coefficient held as exact numbers per
/// distinct loop-free monomial: the products of terms that it pushes one by
/// one as atoms ([`ProductSums`]) become number products, with each distinct
/// product of two monomials normalized once ([`Monomials::product`]). Loop
/// factors and their products carry small real rationals
/// ([`SmallPolynomial`]); only the numerator's own sum, whose spectators may
/// be complex, carries [`Coefficient`]s.
///
/// The coefficients are those of [`LoopPolynomials`]: a normalized sum is
/// the canonical sum of its terms, a normalized product of two terms is the
/// product of their numbers times the normalized product of their
/// non-number parts (exact numbers, [`admits_loop_polynomials`]), and like
/// monomials add their numbers. `None` (use [`LoopPolynomials`]) where a
/// term is not a number times a monomial, a loop factor has a number that is
/// not a small real rational, or an operation overflows.
struct NumericLoopPolynomials<'a> {
    loop_dots: &'a [Atom],
    monomials: Monomials,
    factors: HashMap<&'a [u8], Rc<SmallPolynomial>>,
    products: HashMap<Vec<&'a [u8]>, Rc<SmallPolynomial>>,
}

/// The terms of a sum grouped by their loop factors (in first-seen order),
/// with the product of the other factors of each term.
fn loop_groups<'a>(terms: Vec<AtomView<'a>>) -> Vec<(Vec<AtomView<'a>>, Vec<Atom>)> {
    let mut groups = Vec::<(Vec<AtomView<'a>>, Vec<Atom>)>::new();
    let mut group_indices = HashMap::<Vec<&'a [u8]>, usize>::new();
    for term in terms {
        let (loop_factors, spectators): (Vec<_>, Vec<_>) = product_factors(term)
            .into_iter()
            .partition(|factor| factor.contains_symbol(S.k));
        let spectator = if spectators.is_empty() {
            Atom::num(1)
        } else {
            Atom::mul_many(spectators)
        };
        let key = loop_factors.iter().map(AtomView::get_data).collect();
        match group_indices.entry(key) {
            Entry::Occupied(entry) => groups[*entry.get()].1.push(spectator),
            Entry::Vacant(entry) => {
                entry.insert(groups.len());
                groups.push((loop_factors, vec![spectator]));
            }
        }
    }
    groups
}

impl<'a> NumericLoopPolynomials<'a> {
    fn one(&self) -> SmallPolynomial {
        HashMap::from([(vec![0; self.loop_dots.len()], vec![(0, SmallRational::ONE)])])
    }

    fn monomial(&self, dot: AtomView, exponent: u8) -> Option<SmallPolynomial> {
        let index = self.loop_dots.binary_search(&dot.to_owned()).ok()?;
        let mut exponents = vec![0; self.loop_dots.len()];
        exponents[index] = exponent;
        Some(HashMap::from([(exponents, vec![(0, SmallRational::ONE)])]))
    }

    fn add_exponents(left: &Exponents, right: &Exponents) -> Option<Exponents> {
        // `coefficient_list` stores the exponents as `i8`.
        left.iter()
            .zip(right)
            .map(|(left, right)| left.checked_add(*right).filter(|sum| *sum <= 127))
            .collect()
    }

    fn mul(&mut self, left: &SmallPolynomial, right: &SmallPolynomial) -> Option<SmallPolynomial> {
        let mut product = SmallPolynomial::new();
        for (left_exponents, left_coefficient) in left {
            for (right_exponents, right_coefficient) in right {
                let exponents = Self::add_exponents(left_exponents, right_exponents)?;
                let target = product.entry(exponents).or_default();
                for &(left_monomial, left_value) in left_coefficient {
                    for &(right_monomial, right_value) in right_coefficient {
                        let (scale, monomial) = self
                            .monomials
                            .small_product(left_monomial, right_monomial)?;
                        add_small(target, monomial, scale.mul(left_value)?.mul(right_value)?)?;
                    }
                }
            }
        }
        Some(clean_small(product))
    }

    /// The whole numerator's polynomial times a loop factor's.
    fn mul_numeric(
        &mut self,
        left: &NumericPolynomial,
        right: &SmallPolynomial,
    ) -> Option<NumericPolynomial> {
        let mut product = NumericPolynomial::new();
        for (left_exponents, left_coefficient) in left {
            for (right_exponents, right_coefficient) in right {
                let exponents = Self::add_exponents(left_exponents, right_exponents)?;
                let target = product.entry(exponents).or_default();
                for (left_monomial, left_value) in left_coefficient {
                    for &(right_monomial, right_value) in right_coefficient {
                        let (scale, monomial) =
                            self.monomials.product(*left_monomial, right_monomial)?;
                        add_numeric(
                            target,
                            monomial,
                            scale * left_value.clone() * right_value.coefficient(),
                        );
                    }
                }
            }
        }
        Some(clean_numeric(product))
    }

    fn factor(&mut self, factor: AtomView<'a>) -> Option<Rc<SmallPolynomial>> {
        if let Some(polynomial) = self.factors.get(factor.get_data()) {
            return Some(polynomial.clone());
        }
        let polynomial = match factor {
            AtomView::Fun(function) if function.get_symbol() == S.dot => {
                self.monomial(factor, 1)?
            }
            AtomView::Add(_) | AtomView::Mul(_) => self.small_sum(sum_terms(factor))?,
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let exponent = match exponent {
                    AtomView::Num(number) => match number.get_coeff_view() {
                        CoefficientView::Natural(exponent, 1, 0, _) => u8::try_from(exponent)
                            .ok()
                            .filter(|exponent| (1..=127).contains(exponent))?,
                        _ => return None,
                    },
                    _ => return None,
                };
                match base {
                    AtomView::Fun(function) if function.get_symbol() == S.dot => {
                        self.monomial(base, exponent)?
                    }
                    AtomView::Add(_) => {
                        let base = self.factor(base)?;
                        let mut power = (*base).clone();
                        for _ in 1..exponent {
                            power = self.mul(&power, &base)?;
                        }
                        power
                    }
                    _ => return None,
                }
            }
            _ => return None,
        };
        let polynomial = Rc::new(polynomial);
        self.factors.insert(factor.get_data(), polynomial.clone());
        Some(polynomial)
    }

    fn product(&mut self, factors: &[AtomView<'a>]) -> Option<Rc<SmallPolynomial>> {
        let key = factors.iter().map(AtomView::get_data).collect::<Vec<_>>();
        if let Some(polynomial) = self.products.get(&key) {
            return Some(polynomial.clone());
        }
        let polynomial = match factors {
            [] => Rc::new(self.one()),
            [factor] => self.factor(*factor)?,
            [rest @ .., last] => {
                let rest = self.product(rest)?;
                let last = self.factor(*last)?;
                Rc::new(self.mul(&rest, &last)?)
            }
        };
        self.products.insert(key, polynomial.clone());
        Some(polynomial)
    }

    /// [`LoopPolynomials::sum`] of a sum inside a loop factor, with small
    /// real rational spectators.
    fn small_sum(&mut self, terms: Vec<AtomView<'a>>) -> Option<SmallPolynomial> {
        let mut result = SmallPolynomial::new();
        for (loop_factors, spectators) in loop_groups(terms) {
            let spectators = Atom::add_many(spectators.iter()).expand();
            if spectators.is_zero() {
                continue;
            }
            let spectators = sum_terms(spectators.as_view())
                .into_iter()
                .map(|term| self.monomials.small_term(term))
                .collect::<Option<Vec<_>>>()?;
            let polynomial = self.product(&loop_factors)?;
            for &(weight, spectator_monomial) in &spectators {
                for (exponents, coefficient) in polynomial.iter() {
                    let target = match result.get_mut(exponents) {
                        Some(target) => target,
                        None => result.entry(exponents.clone()).or_default(),
                    };
                    for &(term_monomial, value) in coefficient {
                        let (scale, product) = self
                            .monomials
                            .small_product(spectator_monomial, term_monomial)?;
                        add_small(target, product, scale.mul(weight)?.mul(value)?)?;
                    }
                }
            }
        }
        Some(clean_small(result))
    }

    /// [`LoopPolynomials::sum`] of the numerator: the spectators of each
    /// group of terms with the same loop factors summed and expanded once,
    /// and their terms multiplied with the terms of the group's loop
    /// polynomial as numbers and cached monomial products, straight into the
    /// result.
    fn sum(&mut self, terms: Vec<AtomView<'a>>) -> Option<NumericPolynomial> {
        let mut result = NumericPolynomial::new();
        for (loop_factors, spectators) in loop_groups(terms) {
            let spectators = Atom::add_many(spectators.iter()).expand();
            if spectators.is_zero() {
                continue;
            }
            let spectators = sum_terms(spectators.as_view())
                .into_iter()
                .map(|term| self.monomials.term(term))
                .collect::<Option<Vec<_>>>()?;
            let polynomial = self.product(&loop_factors)?;
            for (weight, spectator_monomial) in &spectators {
                for (exponents, coefficient) in polynomial.iter() {
                    let target = match result.get_mut(exponents) {
                        Some(target) => target,
                        None => result.entry(exponents.clone()).or_default(),
                    };
                    for &(term_monomial, value) in coefficient {
                        let (scale, product) =
                            self.monomials.product(*spectator_monomial, term_monomial)?;
                        add_numeric(
                            target,
                            product,
                            scale * weight.clone() * value.coefficient(),
                        );
                    }
                }
            }
        }
        Some(clean_numeric(result))
    }

    /// The atom polynomial of [`LoopPolynomials`].
    fn to_atoms(&self, polynomial: NumericPolynomial) -> LoopPolynomial {
        polynomial
            .into_iter()
            .map(|(exponents, coefficient)| (exponents, self.monomials.to_atom(&coefficient)))
            .filter(|(_, coefficient)| !coefficient.is_zero())
            .collect()
    }
}

/// The loop polynomial of a numerator with its negative-power factors, by
/// [`NumericLoopPolynomials`]; `None` where those do not apply.
fn numeric_loop_polynomial<'a>(
    loop_dots: &'a [Atom],
    canonical_numerator: AtomView<'a>,
    negative_powers: &'a [Atom],
) -> Option<LoopPolynomial> {
    let mut polynomials = NumericLoopPolynomials {
        loop_dots,
        monomials: Monomials::new(),
        factors: HashMap::new(),
        products: HashMap::new(),
    };
    let mut polynomial = polynomials.sum(sum_terms(canonical_numerator))?;
    for factor in negative_powers {
        let factor = polynomials.factor(factor.as_view())?;
        polynomial = polynomials.mul_numeric(&polynomial, &factor)?;
    }
    // The cached loop polynomials are no longer needed for the conversion.
    polynomials.factors.clear();
    polynomials.products.clear();
    Some(polynomials.to_atoms(polynomial))
}

/// The terms of a sum, or the expression itself.
fn sum_terms(expression: AtomView) -> Vec<AtomView> {
    match expression {
        AtomView::Add(sum) => sum.iter().collect(),
        term => vec![term],
    }
}

/// The factors of a product, or the expression itself.
fn product_factors(expression: AtomView) -> Vec<AtomView> {
    match expression {
        AtomView::Mul(product) => product.iter().collect(),
        factor => vec![factor],
    }
}

/// Whether expanding the term `term` (no sum, no sum factor) returns it
/// unchanged: every factor does. A number, symbol, function or power of a
/// symbol or function to a number is kept as is; other factors, such as a
/// power of a sum, are expanded once and remembered in `fixed` (a short
/// list: a numerator has a few distinct denominators; `false` beyond).
fn expands_to_itself(term: AtomView, fixed: &mut Vec<Atom>) -> bool {
    const MAX_FIXED_FACTORS: usize = 64;
    let mut factor_expands_to_itself = |factor: AtomView| match factor {
        AtomView::Num(_) | AtomView::Var(_) | AtomView::Fun(_) => true,
        AtomView::Add(_) | AtomView::Mul(_) => false,
        AtomView::Pow(power)
            if matches!(
                power.get_base_exp(),
                (AtomView::Var(_) | AtomView::Fun(_), AtomView::Num(_))
            ) =>
        {
            true
        }
        AtomView::Pow(_) => {
            if fixed.iter().any(|known| known.as_view() == factor) {
                return true;
            }
            if fixed.len() == MAX_FIXED_FACTORS {
                return false;
            }
            let expanded = factor.expand();
            if expanded.as_view() != factor {
                return false;
            }
            fixed.push(expanded);
            true
        }
    };
    match term {
        AtomView::Add(_) => false,
        AtomView::Mul(product) => product.iter().all(factor_expands_to_itself),
        factor => factor_expands_to_itself(factor),
    }
}

/// The homogeneous parts of a backend request, `None` when it is scaleless.
/// Forcer renders [`split`] directly; the RustRed p backend lowers these
/// parts onto the family slots.
#[allow(dead_code)]
pub(crate) fn split_parts(request: &PropagatorRequest) -> Result<Option<Vec<PPart>>, VakintError> {
    Ok(split(
        request.family,
        request.powers,
        request.canonical_numerator,
        &request.external,
    )?
    .map(|split| split.parts))
}

#[cfg(test)]
mod tests {
    use symbolica::atom::Atom;

    use super::*;
    use crate::p_integral::P_FAMILIES;

    fn parse(text: &str) -> Atom {
        Vakint::initialize_vakint_symbols();
        vk_parse!(text).unwrap()
    }

    fn parts(family: &str, powers: &[i64], numerator: &str, external: &str) -> Option<Vec<PPart>> {
        let family = P_FAMILIES.iter().find(|f| f.name == family).unwrap();
        let numerator = parse(numerator);
        let request = PropagatorRequest {
            family,
            powers,
            canonical_numerator: numerator.as_view(),
            external: parse(external),
        };
        split_parts(&request).unwrap()
    }

    #[test]
    fn canonical_line_cache_preserves_every_registry_family() {
        for family in &P_FAMILIES {
            let expected = canonical_lines_uncached(family).unwrap();
            assert_eq!(
                canonical_lines(family).unwrap(),
                expected,
                "{}",
                family.name
            );
            let mut local_copy = canonical_lines(family).unwrap();
            local_copy[0].from += 100;
            assert_eq!(
                canonical_lines(family).unwrap(),
                expected,
                "{}",
                family.name
            );
        }
    }

    #[test]
    fn canonical_line_cache_does_not_alias_custom_families_or_errors() {
        let registered = &P_FAMILIES[0];
        let original = canonical_lines(registered).unwrap();
        let custom = PFamilyData {
            name: registered.name,
            loops: 1,
            forcer_topology: None,
            vacuum_parent: "custom cache test",
            lines: "topo(prop(1,edge(7,9),2*k(1),0,pow(1))*prop(2,edge(9,7),k(1)+p(1),0,pow(2)))",
            isps: &[],
        };
        let expected = vec![
            CanonicalLine {
                from: 7,
                to: 9,
                loops: vec![2],
                external: 0,
            },
            CanonicalLine {
                from: 9,
                to: 7,
                loops: vec![1],
                external: 1,
            },
        ];
        assert_eq!(canonical_lines(&custom).unwrap(), expected);
        assert_eq!(canonical_lines(&custom).unwrap(), expected);
        for (loops, lines) in [
            (2, registered.lines),
            (1, "topo(prop(1,edge(1,2),k(1),0,pow(1)))"),
            (
                1,
                "topo(prop(1,edge(1,2),k(1)/2,0,pow(1))*prop(2,edge(2,1),k(1)-p(1),0,pow(2)))",
            ),
        ] {
            let invalid = PFamilyData {
                loops,
                lines,
                ..custom
            };
            let cached = canonical_lines(&invalid).unwrap_err();
            let uncached = canonical_lines_uncached(&invalid).unwrap_err();
            assert_eq!(
                std::mem::discriminant(&cached),
                std::mem::discriminant(&uncached)
            );
            assert_eq!(cached.to_string(), uncached.to_string());
        }
        assert_eq!(canonical_lines(registered).unwrap(), original);
    }

    // Preserve the former linear grouping as a test oracle for both part
    // ordering and polynomial term ordering. Routing and degree conventions
    // are unchanged; input is the already expanded scalar numerator.
    fn linear_parts(
        numerator: &Atom,
        external: &Atom,
        loops: usize,
        total_power: i64,
    ) -> Vec<PPart> {
        let (head, sigma) = external_head(external).unwrap();
        let mut loop_dots = Vec::new();
        numerator.visitor(&mut |part| {
            if let AtomView::Fun(function) = part
                && function.get_symbol() == S.dot
                && part.contains_symbol(S.k)
            {
                scalar_product(part, &head, sigma, loops).unwrap();
                loop_dots.push(part.to_owned());
                return false;
            }
            true
        });
        loop_dots.sort();
        loop_dots.dedup();
        let mut groups = Vec::<((i64, Atom), Vec<(Fraction, PMonomial)>)>::new();
        for (monomial, coefficient) in numerator.coefficient_list::<i8>(&loop_dots) {
            let (sign, canonical) =
                canonical_monomial(monomial.as_view(), &head, sigma, loops).unwrap();
            let (value, rest) = split_rational(&coefficient);
            let term = (Fraction::new(sign as i128, 1).mul(value), canonical);
            let key = (term.1.degree(), rest);
            match groups.iter_mut().find(|(group, _)| *group == key) {
                Some((_, terms)) => terms.push(term),
                None => groups.push((key, vec![term])),
            }
        }
        groups
            .into_iter()
            .map(|((degree, coefficient), polynomial)| PPart {
                scale_power: 2 * loops as i64 - total_power + degree,
                wick_sign: if (total_power + degree) % 2 == 0 {
                    1
                } else {
                    -1
                },
                coefficient,
                polynomial,
            })
            .collect()
    }

    #[test]
    fn indexed_part_lookup_preserves_linear_order_with_spectators_and_negative_powers() {
        let family = &P_FAMILIES[1];
        let external = parse("-p(3)");
        let first = parse("dot(k(1),p(3))");
        let second = parse("dot(k(2),p(3))");
        let epsilon = parse("ε");
        let mut numerator = parse("3/7*dot(k(1),k(2)) - 2/5");
        for index in 1..=32 {
            let power = index;
            let spectator = function!(S.x, index) * (Atom::one() + &epsilon * Atom::num(index));
            numerator +=
                spectator * (first.pow(power) / Atom::num(3) + second.pow(power) / Atom::num(3));
        }
        for (powers, extra_numerator) in [
            ([1, 1, 1, 1, 2], Atom::one()),
            ([-1, 1, 1, 1, 2], parse("dot(k(1),k(1))")),
        ] {
            let actual = split(family, &powers, numerator.as_view(), &external)
                .unwrap()
                .unwrap();
            let total_power = powers.iter().filter(|power| **power > 0).sum();
            let expanded = (&numerator * extra_numerator).expand();
            let expected = linear_parts(&expanded, &external, family.loops, total_power);
            assert!(
                expected.len() > 16,
                "fixture must exercise many distinct groups"
            );
            assert!(expected.iter().any(|part| part.polynomial.len() > 1));
            assert_eq!(actual.parts, expected);
        }
    }

    #[test]
    fn parts_are_homogeneous_in_canonical_scalar_products() {
        // The user head p(3) is -Q: odd powers of k·p(3) flip the sign.
        let split = parts(
            "P2L",
            &[1, 1, 1, 1, 2],
            "dot(k(1),p(3))^2*dot(k(2),k(2)) + 3*dot(k(1),p(3)) + (1+ε)*dot(k(2),p(3)) - 1/2",
            "-p(3)",
        )
        .unwrap();
        let q = |i| PDot::External(i);
        let one = |n, d| Fraction::new(n, d);
        // Degree 3: w = 4 - 6 + 3, Wick sign (-1)^(6+3).
        let cubic = split.iter().find(|part| part.scale_power == 1).unwrap();
        assert_eq!(cubic.wick_sign, -1);
        assert_eq!(cubic.coefficient, Atom::num(1));
        assert_eq!(
            cubic.polynomial,
            [(
                one(1, 1),
                PMonomial(vec![(q(1), 2), (PDot::Loops(2, 2), 1)])
            )]
        );
        // Degree 1, rational and (1+ε) coefficients in separate parts.
        let linear = split
            .iter()
            .filter(|part| part.scale_power == -1)
            .collect::<Vec<_>>();
        assert_eq!(linear.len(), 2);
        assert!(linear.iter().all(|part| part.wick_sign == -1));
        let rational = linear
            .iter()
            .find(|part| part.coefficient == Atom::num(1))
            .unwrap();
        assert_eq!(
            rational.polynomial,
            [(one(-3, 1), PMonomial(vec![(q(1), 1)]))]
        );
        let shifted = linear
            .iter()
            .find(|part| part.coefficient != Atom::num(1))
            .unwrap();
        assert_eq!(shifted.coefficient, parse("1+ε"));
        assert_eq!(
            shifted.polynomial,
            [(one(-1, 1), PMonomial(vec![(q(2), 1)]))]
        );
        // Degree 0: w = 4 - 6, Wick sign +1.
        let constant = split.iter().find(|part| part.scale_power == -2).unwrap();
        assert_eq!(constant.wick_sign, 1);
        assert_eq!(constant.polynomial, [(one(-1, 2), PMonomial::default())]);
        assert_eq!(split.len(), 4);
    }

    #[test]
    fn negative_powers_join_the_numerator_and_scaleless_sectors_are_none() {
        // A negative power on line 1 (momentum k(1)) is the numerator k(1)²:
        // one scalar product.
        let split = parts("P1L", &[-1, 1], "1", "p(1)");
        assert_eq!(split, None, "a single remaining line is scaleless");
        let split = parts("P2L", &[-1, 1, 1, 1, 1], "1", "p(1)").unwrap();
        assert_eq!(split.len(), 1);
        // w = 2L - Σa + r = 4 - 4 + 1.
        assert_eq!(split[0].scale_power, 1);
        assert_eq!(
            split[0].polynomial,
            [(Fraction::new(1, 1), PMonomial(vec![(PDot::Loops(1, 1), 1)]))]
        );
        // A pinched loop (the two lines of k(2) gone) is scaleless.
        assert_eq!(
            parts("P2L", &[1, 0, 0, 1, 1], "dot(k(1),p(1))", "p(1)"),
            None
        );
    }

    type SplitSummary = Result<
        Option<(
            Vec<CanonicalLine>,
            Vec<(usize, i64)>,
            Vec<Vec<Fraction>>,
            Vec<PPart>,
        )>,
        String,
    >;

    fn summary(split: Result<Option<PSplit>, VakintError>) -> SplitSummary {
        split
            .map(|split| {
                split.map(|split| {
                    (
                        split.lines,
                        split.remaining,
                        split.loop_momenta,
                        split.parts,
                    )
                })
            })
            .map_err(|error| error.to_string())
    }

    /// [`split`] with the grouped parts only, an error where they do not
    /// apply. The grouped parts are the same whether their products are
    /// merged in the production batches, after every doubling of the pending
    /// products (bound 0) or once at the end (`usize::MAX`, the merge before
    /// the batches).
    fn split_grouped_only(
        family: &PFamilyData,
        powers: &[i64],
        numerator: AtomView,
        external: &Atom,
    ) -> SplitSummary {
        let with_bound = |bound: usize| {
            let parts: NumeratorParts = match bound {
                0 => |numerator, negative_powers, collector| {
                    grouped_parts_with(numerator, negative_powers, collector, 0)
                        .ok_or_else(not_grouped)
                },
                usize::MAX => |numerator, negative_powers, collector| {
                    grouped_parts_with(numerator, negative_powers, collector, usize::MAX)
                        .ok_or_else(not_grouped)
                },
                _ => |numerator, negative_powers, collector| {
                    grouped_parts(numerator, negative_powers, collector).ok_or_else(not_grouped)
                },
            };
            summary(split_with(family, powers, numerator, external, parts))
        };
        let grouped = with_bound(PENDING_PRODUCT_BYTES);
        for bound in [0, usize::MAX] {
            assert_eq!(with_bound(bound), grouped, "bound {bound}: {numerator}");
        }
        grouped
    }

    fn not_grouped() -> VakintError {
        VakintError::InvalidNumerator("the grouped split does not apply".into())
    }

    /// The parts of [`split`] before the grouped split and the skip of the
    /// dot-notation passes: expand the numerator with its negative-power
    /// factors, always run the passes, expand again and collect with
    /// `coefficient_list`.
    fn original_parts(
        canonical_numerator: AtomView,
        negative_powers: &[Atom],
        collector: &PartCollector,
    ) -> Result<Vec<PPart>, VakintError> {
        let mut numerator = canonical_numerator.to_owned();
        for factor in negative_powers {
            numerator *= factor.clone();
        }
        let numerator = Vakint::dot_notation_passes(numerator.expand()).expand();
        collected_parts(numerator.as_view(), &[], collector)
    }

    /// The last stage of [`original_parts`]: the parts of a numerator that is
    /// multiplied by its negative-power factors, expanded and in dot notation
    /// already (`negative_powers` is ignored).
    fn collected_parts(
        numerator: AtomView,
        _negative_powers: &[Atom],
        collector: &PartCollector,
    ) -> Result<Vec<PPart>, VakintError> {
        let loop_dots = collector.loop_dots([numerator])?;
        collector.parts(numerator.coefficient_list::<i8>(&loop_dots))
    }

    /// The split before the grouped split and the skip of the dot-notation
    /// passes ([`original_parts`]).
    fn split_original(
        family: &PFamilyData,
        powers: &[i64],
        numerator: AtomView,
        external: &Atom,
    ) -> SplitSummary {
        summary(split_with(
            family,
            powers,
            numerator,
            external,
            original_parts,
        ))
    }

    /// [`split`] equals [`split_reference`] and the original split
    /// ([`split_original`]), errors included; with `grouped` the grouped
    /// parts apply and are the reference parts, else they defer to the
    /// reference.
    fn assert_split_matches_reference(
        family: &PFamilyData,
        powers: &[i64],
        numerator: AtomView,
        external: &Atom,
        grouped: bool,
    ) -> SplitSummary {
        let reference = summary(split_reference(family, powers, numerator, external));
        assert_eq!(
            split_original(family, powers, numerator, external),
            reference,
            "{} {powers:?} {numerator}",
            family.name
        );
        assert_eq!(
            summary(split(family, powers, numerator, external)),
            reference,
            "{} {powers:?} {numerator}",
            family.name
        );
        let grouped_only = split_grouped_only(family, powers, numerator, external);
        if grouped {
            assert_eq!(grouped_only, reference, "{} {powers:?}", family.name);
        } else {
            assert!(
                grouped_only.is_err(),
                "{} {powers:?} {numerator}",
                family.name
            );
        }
        reference
    }

    fn family(name: &str) -> &'static PFamilyData {
        P_FAMILIES
            .iter()
            .find(|family| family.name == name)
            .unwrap()
    }

    #[test]
    fn grouped_split_equals_the_expanded_split() {
        let _real = symbolica::symbol!("vakint::parts_real_spectator"; Real);
        let cases: &[(&str, &[i64], &str, &str, bool)] = &[
            ("P1L", &[1, 1], "0", "p(1)", true),
            ("P1L", &[1, 1], "1", "p(1)", true),
            (
                "P2L",
                &[1, 1, 1, 1, 2],
                "dot(k(1),p(3))^2*dot(k(2),k(2))+3*dot(k(1),p(3))+(1+ε)*dot(k(2),p(3))-1/2",
                "-p(3)",
                true,
            ),
            // Repeated loop factors with different spectators, and a negative power.
            (
                "P2L",
                &[-2, 1, 1, 1, 1],
                "x*(dot(k(1),p(1))+dot(k(2),p(1)))^2+y*(dot(k(1),p(1))+dot(k(2),p(1)))^2\
                 +(x+y)*dot(k(1),k(2))+z",
                "p(1)",
                true,
            ),
            // Spectator sums and powers of sums.
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "(1+ε)^2*(dot(k(1),p(1))+dot(p(1),p(1)))*dot(k(2),p(1))\
                 +(nf+ε)*(dot(k(1),p(1))+dot(p(1),p(1)))*dot(k(2),p(1))\
                 +(nf+ε)*(2-ε)^-1*dot(k(1),k(1))",
                "p(1)",
                true,
            ),
            // The loop dot k(1)·p(1) cancels: the reference never sees it.
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "(dot(k(1),p(1))+y)*x-x*dot(k(1),p(1))+dot(k(2),k(2))",
                "p(1)",
                true,
            ),
            // One monomial from several loop factors; a product and a power
            // at the top.
            (
                "P2L",
                &[1, 1, 1, 1, 2],
                "x*(dot(k(1),p(1))+dot(k(2),p(1)))^2+y*dot(k(1),p(1))^2+z",
                "-p(1)",
                true,
            ),
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "x*(dot(k(1),p(1))+dot(p(1),p(1)))^3",
                "-p(1)",
                true,
            ),
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "(dot(k(1),p(1))+ε)^2",
                "p(1)",
                true,
            ),
            // Nested sums of products with loop-free factors at every level,
            // a repeated inner loop factor and a merged rational power.
            (
                "P2L",
                &[-1, 1, 1, 1, 1],
                "((dot(k(1),p(1))+dot(k(2),p(1)))*ε^2*x+(dot(k(1),p(1))+dot(k(2),p(1)))*(1+ε)*y\
                 +z*dot(k(1),k(1)))*(dot(k(2),k(2))+(2-ε)^-1*dot(p(1),p(1)))\
                 *(dot(k(1),k(2))+nf)+(dot(k(1),p(1))+dot(k(2),p(1)))^2*(ε-1)",
                "-p(1)",
                true,
            ),
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "x^(1/2)*dot(k(1),p(1))*(x^(1/2)*dot(k(2),p(1))+1)",
                "p(1)",
                true,
            ),
            // Rational powers of symbols and functions multiply as in the
            // expansion's powers of terms, (x^(1/3))^3 = x.
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "(x^(1/3)*dot(k(1),p(1))+1)^3",
                "p(1)",
                true,
            ),
            (
                "P2L",
                &[1, 1, 1, -1, 1],
                "x^(2/3)*(x^(1/3)*dot(k(1),p(1))+x^(-1/3)*dot(k(2),p(1))+f(y)^(-1/3))^2\
                 *(f(y)^(1/3)*dot(k(1),k(2))+1)",
                "-p(1)",
                true,
            ),
            // Indexed vectors that no dot-notation pass matches.
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "p(0,mu)*p(0,nu)*dot(k(1),p(1))*ε+g(mu,nu)*dot(p(0),p(0))*dot(k(1),k(1))",
                "p(1)",
                true,
            ),
            // Deferred: a product squares a denominator, which the second
            // expansion rewrites; a product squares an indexed vector, which a
            // pass contracts; inexact numbers; a non-integer power of a sum.
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "(2-ε)^-1*(dot(k(1),p(1))*(2-ε)^-1+1)*dot(k(2),k(2))",
                "p(1)",
                false,
            ),
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "p(0,mu)*(p(0,mu)*dot(k(1),p(1))+dot(k(2),p(1)))",
                "p(1)",
                false,
            ),
            ("P2L", &[1, 1, 1, 1, 1], "1.5*dot(k(1),p(1))", "p(1)", false),
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "(1+ε)^(1/2)*dot(k(1),p(1))",
                "p(1)",
                false,
            ),
            // Deferred: powers of numbers, which a product merges as
            // 2^(1/3)*2^(1/3) = 2^(2/3) but the expansion's power of a term
            // as (2^(1/3))^2 = 4^(1/3); also square roots, and complex
            // exponents.
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "(2^(1/3)*dot(k(1),p(1))+1)^2",
                "p(1)",
                false,
            ),
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "(3^(1/2)*dot(k(1),k(2))+dot(k(2),p(1))*2^(1/3))^3",
                "p(1)",
                false,
            ),
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "(2^(1/2)*dot(k(1),p(1))+2^(1/3))^3",
                "p(1)",
                false,
            ),
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "3^(1/2)*dot(k(1),p(1))*(3^(1/2)*dot(k(2),p(1))+1)",
                "p(1)",
                false,
            ),
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "x^(2𝑖)*dot(k(1),p(1))*(x^(2𝑖)*dot(k(2),p(1))+1)",
                "p(1)",
                false,
            ),
            // Deferred: `abs`, whose even powers of a real argument a power
            // of a term rewrites, (abs(r))^2 = r^2, but a product keeps,
            // abs(r)*abs(r) = abs(r)^2.
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "(abs(parts_real_spectator)*dot(k(1),p(1))+abs(x))^2",
                "p(1)",
                false,
            ),
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "abs(parts_real_spectator)*dot(k(1),p(1))\
                 *(abs(parts_real_spectator)*dot(k(2),p(1))+1)",
                "p(1)",
                false,
            ),
            // Deferred errors: a loop momentum outside a scalar product and a
            // loop momentum contracted with another external vector.
            (
                "P2L",
                &[1, 1, 1, 1, 1],
                "f(dot(k(1),p(1)))*dot(k(2),p(1))",
                "p(1)",
                false,
            ),
            ("P2L", &[1, 1, 1, 1, 1], "dot(k(1),p(2))", "p(1)", false),
        ];
        for &(name, powers, numerator, external, grouped) in cases {
            let numerator = parse(numerator);
            let _ = assert_split_matches_reference(
                family(name),
                powers,
                numerator.as_view(),
                &parse(external),
                grouped,
            );
        }
    }

    /// A [`split`] input as `dump_split_input` writes it.
    fn split_record(text: &str) -> (&'static PFamilyData, Vec<i64>, Atom, Atom) {
        let field = |key: &str| {
            text.lines()
                .find_map(|line| line.strip_prefix(key)?.strip_prefix(' '))
                .unwrap_or_else(|| panic!("a split record without {key}"))
        };
        let powers = field("powers")
            .split(',')
            .map(|power| power.parse().unwrap())
            .collect();
        (
            family(field("family")),
            powers,
            parse(field("external")),
            parse(field("numerator")),
        )
    }

    /// Real five-loop inputs: requests of the nf^3 self-energy diagrams
    /// 58440, 58844 and 59566 on the RustRed path, written by
    /// `VAKINT_DUMP_SPLIT_INPUTS` (records separated by `---` lines).
    #[test]
    fn grouped_split_reproduces_five_loop_requests() {
        let records = include_str!("../../tests/fixtures/p_split_requests.txt")
            .split("\n---\n")
            .collect::<Vec<_>>();
        assert!(records.len() >= 10);
        for record in records {
            let (family, powers, external, numerator) = split_record(record);
            let split = assert_split_matches_reference(
                family,
                &powers,
                numerator.as_view(),
                &external,
                true,
            );
            assert!(split.unwrap().is_some_and(|split| !split.3.is_empty()));
        }
    }

    /// Whether `expression` holds a power of a number, such as `3^(1/2)`,
    /// which the grouped parts defer to the expanded parts.
    fn has_number_power(expression: AtomView) -> bool {
        let mut found = false;
        expression.visitor(&mut |part| {
            found |= matches!(part, AtomView::Pow(power)
                if matches!(power.get_base_exp().0, AtomView::Num(_)));
            !found
        });
        found
    }

    /// Opt-in diagnostic: one dumped request named by `VAKINT_SPLIT_RECORD`,
    /// split in the mode `VAKINT_SPLIT_RECORD_MODE` (`split` by default,
    /// `grouped`, `grouped_unbounded` (one merge per monomial), `reference`, `expand` or `expand_merging` of the numerator
    /// with its negative powers), with the time and the process high-water
    /// mark printed (one mode per process: the mark only grows).
    #[test]
    fn split_one_dumped_request() {
        let Some(path) = std::env::var_os("VAKINT_SPLIT_RECORD") else {
            return;
        };
        let high_water_mark = || {
            std::fs::read_to_string("/proc/self/status")
                .ok()
                .and_then(|status| {
                    status
                        .lines()
                        .find_map(|line| Some(line.strip_prefix("VmHWM:")?.trim().to_string()))
                })
                .unwrap_or_default()
        };
        let (family, powers, external, numerator) =
            split_record(&std::fs::read_to_string(path).unwrap());
        let mode = std::env::var("VAKINT_SPLIT_RECORD_MODE").unwrap_or_else(|_| "split".into());
        eprintln!(
            "{} {powers:?}: numerator {} bytes, high-water mark {}",
            family.name,
            numerator.as_view().get_byte_size(),
            high_water_mark()
        );
        let started = std::time::Instant::now();
        let outcome = match mode.as_str() {
            "split" => format!(
                "{} parts",
                summary(split(family, &powers, numerator.as_view(), &external))
                    .map_or(0, |split| split.map_or(0, |split| split.3.len()))
            ),
            "grouped" | "grouped_unbounded" => {
                let parts: NumeratorParts = if mode == "grouped" {
                    |numerator, negative_powers, collector| {
                        grouped_parts(numerator, negative_powers, collector).ok_or_else(not_grouped)
                    }
                } else {
                    |numerator, negative_powers, collector| {
                        grouped_parts_with(numerator, negative_powers, collector, usize::MAX)
                            .ok_or_else(not_grouped)
                    }
                };
                format!(
                    "applies {}",
                    split_with(family, &powers, numerator.as_view(), &external, parts).is_ok()
                )
            }
            "reference" => format!(
                "{} parts",
                summary(split_reference(
                    family,
                    &powers,
                    numerator.as_view(),
                    &external
                ))
                .map_or(0, |split| split.map_or(0, |split| split.3.len()))
            ),
            "expand" | "expand_merging" => {
                let mut product = numerator.clone();
                for (line, &power) in canonical_lines(family).unwrap().iter().zip(&powers) {
                    if power < 0 {
                        let momentum = line.momentum(&external);
                        product *= function!(S.dot, &momentum, &momentum).pow(Atom::num(-power));
                    }
                }
                let expanded = if mode == "expand" {
                    product.expand()
                } else {
                    crate::utils::expand_merging(product.as_view())
                };
                format!("{} bytes", expanded.as_view().get_byte_size())
            }
            other => panic!("unknown mode {other}"),
        };
        eprintln!(
            "{mode}: {outcome} in {:.1} s, high-water mark {}",
            started.elapsed().as_secs_f64(),
            high_water_mark()
        );
    }

    /// The same on a whole `VAKINT_DUMP_SPLIT_INPUTS` directory named by
    /// `VAKINT_SPLIT_CORPUS` (opt-in: real corpora are far too large to
    /// vendor), against the original split ([`original_parts`], run stage by
    /// stage). The grouped parts must apply to every request without a power
    /// of a number. Also checks that the dot-notation passes are no-ops where
    /// they are skipped, counts where the second expansion changes anything,
    /// and prints timings.
    #[test]
    fn grouped_split_reproduces_a_dumped_corpus() {
        let Some(directory) = std::env::var_os("VAKINT_SPLIT_CORPUS") else {
            eprintln!("VAKINT_SPLIT_CORPUS is not set: no corpus to compare");
            return;
        };
        let mut paths = std::fs::read_dir(directory)
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .collect::<Vec<_>>();
        paths.sort();
        let seconds = |started: std::time::Instant| started.elapsed().as_secs_f64();
        let mut totals = [0.0; 5];
        let (mut skipped, mut second_expansion_changed, mut number_powers) = (0, 0, 0);
        let mut declined = Vec::new();
        for path in &paths {
            let (family, powers, external, numerator) =
                split_record(&std::fs::read_to_string(path).unwrap());
            // The original split, stage by stage.
            let mut product = numerator.clone();
            for (line, &power) in canonical_lines(family).unwrap().iter().zip(&powers) {
                if power < 0 {
                    let momentum = line.momentum(&external);
                    product *= function!(S.dot, &momentum, &momentum).pow(Atom::num(-power));
                }
            }
            let started = std::time::Instant::now();
            let expanded = product.expand();
            let expand_seconds = seconds(started);
            assert_eq!(
                crate::utils::expand_merging(product.as_view()),
                expanded,
                "{}",
                path.display()
            );
            let started = std::time::Instant::now();
            let passed = Vakint::dot_notation_passes(expanded.clone());
            let passes_seconds = seconds(started);
            if !Vakint::dot_notation_passes_may_apply(expanded.as_view()) {
                skipped += 1;
                assert_eq!(passed, expanded, "{}", path.display());
            }
            let started = std::time::Instant::now();
            let reexpanded = passed.expand();
            let second_expansion_seconds = seconds(started);
            second_expansion_changed += usize::from(reexpanded != passed);
            let started = std::time::Instant::now();
            let original = summary(split_with(
                family,
                &powers,
                reexpanded.as_view(),
                &external,
                collected_parts,
            ));
            let original_seconds =
                expand_seconds + passes_seconds + second_expansion_seconds + seconds(started);
            let started = std::time::Instant::now();
            let actual = summary(split(family, &powers, numerator.as_view(), &external));
            let split_seconds = seconds(started);
            assert_eq!(actual, original, "{}", path.display());
            let grouped = split_grouped_only(family, &powers, numerator.as_view(), &external);
            if grouped != original {
                assert!(grouped.is_err(), "{}", path.display());
                if has_number_power(numerator.as_view()) {
                    number_powers += 1;
                } else {
                    declined.push(path.file_name().unwrap().to_string_lossy().into_owned());
                }
            }
            for (total, value) in totals.iter_mut().zip([
                original_seconds,
                split_seconds,
                expand_seconds,
                passes_seconds,
                second_expansion_seconds,
            ]) {
                *total += value;
            }
            eprintln!(
                "{}: original {original_seconds:.4} s, split {split_seconds:.4} s \
                 (expansion {expand_seconds:.4} s, passes {passes_seconds:.4} s, \
                 second expansion {second_expansion_seconds:.4} s)",
                path.file_name().unwrap().to_string_lossy()
            );
        }
        eprintln!(
            "{} requests: original {:.3} s, split {:.3} s; in the original: expansion {:.3} s, \
             passes {:.3} s, second expansion {:.3} s; passes skipped on {skipped}, \
             second expansion changed {second_expansion_changed}, deferred for a power of \
             a number {number_powers}",
            paths.len(),
            totals[0],
            totals[1],
            totals[2],
            totals[3],
            totals[4]
        );
        assert!(
            declined.is_empty(),
            "the grouped split does not apply to {declined:?}"
        );
    }

    /// Random nested two-loop numerators: sums of products of loop dots,
    /// their powers and sums with loop-free factors (complex rationals,
    /// rational powers of numbers, symbols and functions, denominators,
    /// functions of sums, indexed vectors), with negative line powers and
    /// `±p(1)`. [`split`] equals the original split, whether or not the
    /// grouped parts apply.
    #[test]
    fn split_equals_the_original_split_on_random_numerators() {
        const DOTS: &[&str] = &[
            "dot(k(1),k(1))",
            "dot(k(1),k(2))",
            "dot(k(2),k(2))",
            "dot(k(1),p(1))",
            "dot(k(2),p(1))",
        ];
        const SPECTATORS: &[&str] = &[
            "x",
            "y",
            "nf",
            "ε",
            "ε^2",
            "ε^-1",
            "(2-ε)^-1",
            "(1-2*ε)^-1",
            "(3+ε+ε^2)^-1",
            "dot(p(1),p(1))",
            "dot(p(1),p(1))^-1",
            "1/3",
            "-2",
            "5/7",
            "2𝑖",
            "(1/2+3𝑖/4)",
            "x^(1/2)",
            "x^(-1/2)",
            "x^(1/3)",
            "x^(-2/3)",
            "f(y)^(1/3)",
            "2^(1/2)",
            "3^(1/2)",
            "2^(1/3)",
            "4^(1/3)",
            "f(x+y)",
            "g(mu,nu)",
            "p(0,mu)*p(0,nu)",
            "(1+ε)",
            "(x+y)",
            "(1+ε)^2",
            "dot(p(0),p(0))",
            "𝜋^-2",
        ];
        struct Random(u64);
        impl Random {
            fn below(&mut self, n: usize) -> usize {
                self.0 ^= self.0 << 13;
                self.0 ^= self.0 >> 7;
                self.0 ^= self.0 << 17;
                (self.0 % n as u64) as usize
            }
            fn pick(&mut self, items: &[&'static str]) -> &'static str {
                items[self.below(items.len())]
            }
            /// A sum of one to three products of one to three factors,
            /// nested to `depth` further levels.
            fn numerator(&mut self, depth: usize) -> String {
                let mut sum = String::new();
                for term in 0..1 + self.below(3) {
                    if term > 0 {
                        sum.push(if self.below(3) == 0 { '-' } else { '+' });
                    }
                    for factor in 0..1 + self.below(3) {
                        if factor > 0 {
                            sum.push('*');
                        }
                        let factor = match self.below(if depth == 0 { 6 } else { 10 }) {
                            0 | 1 => self.pick(SPECTATORS).to_string(),
                            2 | 3 => self.pick(DOTS).to_string(),
                            4 => format!("{}^{}", self.pick(DOTS), 2 + self.below(2)),
                            5 => format!("({}+{})", self.pick(DOTS), self.pick(SPECTATORS)),
                            8 if depth == 1 => format!("({})^2", self.numerator(0)),
                            8 => format!("({}+{})^2", self.pick(DOTS), self.pick(SPECTATORS)),
                            _ => format!("({})", self.numerator(depth - 1)),
                        };
                        sum.push_str(&factor);
                    }
                }
                sum
            }
        }
        let powers: &[&[i64]] = &[
            &[1, 1, 1, 1, 1],
            &[-1, 1, 1, 1, 1],
            &[1, 1, 1, 1, -2],
            &[1, 2, 1, 1, 1],
            &[1, 1, -1, 1, 1],
        ];
        let externals = [parse("p(1)"), parse("-p(1)")];
        let mut random = Random(0x9e37_79b9_7f4a_7c15);
        let (cases, mut grouped) = (600, 0);
        for _ in 0..cases {
            let numerator = parse(&random.numerator(2));
            let powers = powers[random.below(powers.len())];
            let external = &externals[random.below(externals.len())];
            let family = family("P2L");
            let original = split_original(family, powers, numerator.as_view(), external);
            assert_eq!(
                summary(split(family, powers, numerator.as_view(), external)),
                original,
                "{powers:?} {external} {numerator}"
            );
            let grouped_only = split_grouped_only(family, powers, numerator.as_view(), external);
            if grouped_only.is_ok() {
                assert_eq!(grouped_only, original, "{powers:?} {external} {numerator}");
                grouped += 1;
            }
        }
        // Many cases exercise the grouped parts (about half: powers of
        // numbers, indexed vectors that a pass contracts and squared
        // denominators defer).
        assert!(3 * grouped > cases, "{grouped} of {cases} grouped");
    }
}
