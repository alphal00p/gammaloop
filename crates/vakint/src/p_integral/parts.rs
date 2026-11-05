//! The numerator split shared by the p-integral backends.
//!
//! A backend request (`docs/p_integrals.md`, section 5) is one canonical
//! family, its line powers and a scalar numerator in the canonical `k(j)` and
//! the user head `p(i)`. [`split_parts`] turns it into homogeneous parts:
//! negative powers join the numerator, every loop-momentum monomial is written
//! in the canonical scalar products `k(i)·k(j)` and `k(i)·Q` (the user head is
//! `σ·Q`), and the monomials are grouped by degree `r` and by the factor of
//! their coefficient that is not a rational number. Each [`PPart`] is then
//! reduced on its own: Forcer sets `Q·Q = 1`, so parts of different degree
//! cannot share one computation, and a non-rational factor (`1+ε`, `mass²`)
//! stays outside the reduction.
//!
//! `None` means the sector is scaleless (an exact zero), decided here once for
//! every backend: the scaleless rule of `docs/p_integrals.md` (section 1,
//! item 6) on the graph with the pinched lines contracted, or remaining lines
//! that do not span the loop momenta.
//!
//! The Wick sign `(-1)^(Σa+r)` of a part is applied on the Forcer path only
//! (M2 plan, conventions C3/C4): Forcer reduces Euclidean integrals, while the
//! RustRed reduction already returns Minkowski-signed values.

use std::collections::BTreeMap;

use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::coefficient::CoefficientView;
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

/// The split of `topo(NAME(external, powers...))` with a scalar numerator in
/// the canonical `k(j)` and the user's head; `None` when the sector is
/// scaleless (an exact zero).
pub(crate) fn split(
    family: &PFamilyData,
    powers: &[i64],
    canonical_numerator: AtomView,
    external: &Atom,
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
    let mut numerator = canonical_numerator.to_owned();
    for (line, &power) in lines.iter().zip(powers) {
        if power < 0 {
            let momentum = line.momentum(external);
            numerator *= function!(S.dot, &momentum, &momentum).pow(Atom::num(-power));
        }
    }
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
    let total_power: i64 = remaining.iter().map(|(_, power)| power).sum();

    // Homogeneous parts: loop-momentum monomials, grouped by degree and by
    // the non-rational factor of their coefficient.
    let numerator = Vakint::convert_to_dot_notation(numerator.as_view()).expand();
    let mut loop_dots = Vec::new();
    let mut error = None;
    numerator.visitor(&mut |part| {
        if error.is_some() {
            return false;
        }
        if let AtomView::Fun(function) = part
            && function.get_symbol() == S.dot
            && part.contains_symbol(S.k)
        {
            match scalar_product(part, &head, sigma, loops) {
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
    loop_dots.sort();
    loop_dots.dedup();
    let mut groups = Vec::<((i64, Atom), Vec<(Fraction, PMonomial)>)>::new();
    for (monomial, coefficient) in numerator.coefficient_list::<i8>(&loop_dots) {
        if coefficient.contains_symbol(S.k) {
            return Err(tensor_numerator(format!(
                "the loop momenta must appear in scalar products only, found {coefficient}"
            )));
        }
        let (monomial_sign, canonical) =
            canonical_monomial(monomial.as_view(), &head, sigma, loops)?;
        let (value, rest) = split_rational(&coefficient);
        let term = (Fraction::new(monomial_sign as i128, 1).mul(value), canonical);
        let key = (term.1.degree(), rest);
        match groups.iter_mut().find(|(group, _)| *group == key) {
            Some((_, terms)) => terms.push(term),
            None => groups.push((key, vec![term])),
        }
    }
    let parts = groups
        .into_iter()
        .map(|((degree, coefficient), polynomial)| {
            let exponent = total_power + degree;
            PPart {
                scale_power: 2 * loops as i64 - total_power + degree,
                wick_sign: if exponent % 2 == 0 { 1 } else { -1 },
                coefficient,
                polynomial,
            }
        })
        .collect();
    Ok(Some(PSplit {
        lines,
        remaining,
        loop_momenta,
        parts,
    }))
}

/// The homogeneous parts of a backend request, `None` when it is scaleless.
/// Forcer renders [`split`] directly; the RustRed p backend lowers these
/// parts onto the family slots.
#[allow(dead_code)]
pub(crate) fn split_parts(
    request: &PropagatorRequest,
) -> Result<Option<Vec<PPart>>, VakintError> {
    Ok(split(
        request.family,
        request.powers,
        request.canonical_numerator,
        &request.external,
    )?
    .map(|split| split.parts))
}
