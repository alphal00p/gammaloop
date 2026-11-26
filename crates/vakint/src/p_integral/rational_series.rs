//! Bounded exact rational expansion with a regulator-free coefficient ring.
//!
//! Only denominators D(epsilon)*C(spectators) are admitted. The recurrence
//! therefore multiplies spectator polynomials by rational numbers; it never
//! constructs powers of unexpanded AtomField series. All other inputs retain
//! the original general series implementation.
use std::sync::Arc;

use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::CoefficientView,
    domains::{
        integer::{IntegerRing, Z},
        rational::{Q, Rational, RationalField},
        rational_polynomial::RationalPolynomial,
    },
    poly::{PolyVariable, polynomial::MultivariatePolynomial},
};

use super::Laurent;

type Polynomial = MultivariatePolynomial<RationalField, u32>;
const MAX_NODES: usize = 16_384;
const MAX_SUPPORT: usize = 32_768;
const MAX_VARIABLES: usize = 64;
const MAX_DEGREE: usize = 256;
const MAX_DEPTH: usize = 128;
const MAX_SERIES: usize = 64;
const MAX_WORK: usize = 1_048_576;
const MAX_SPECTATOR_KEYS: usize = 64;

#[derive(Clone)]
struct Shape {
    numerator: usize,
    denominator: usize,
    degree: usize,
    denominator_key: Option<Atom>,
}
impl Shape {
    fn leaf() -> Self {
        Self {
            numerator: 1,
            denominator: 1,
            degree: 1,
            denominator_key: None,
        }
    }
    fn one() -> Self {
        Self {
            numerator: 1,
            denominator: 1,
            degree: 0,
            denominator_key: None,
        }
    }
    fn bounded(self) -> Option<Self> {
        (self.numerator <= MAX_SUPPORT
            && self.denominator <= MAX_SUPPORT
            && self.degree <= MAX_DEGREE)
            .then_some(self)
    }
    fn denominator_product(left: Option<Atom>, right: Option<Atom>) -> Option<Atom> {
        match (left, right) {
            (None, key) | (key, None) => key,
            (Some(a), Some(b)) => Some(a * b),
        }
    }
    fn mul(self, rhs: Self) -> Option<Self> {
        Self {
            numerator: self.numerator.checked_mul(rhs.numerator)?,
            denominator: self.denominator.checked_mul(rhs.denominator)?,
            degree: self.degree.checked_add(rhs.degree)?,
            denominator_key: Self::denominator_product(self.denominator_key, rhs.denominator_key),
        }
        .bounded()
    }
    fn add(self, rhs: Self) -> Option<Self> {
        // Expanded coefficients repeat one denominator for every epsilon power.
        // Reuse that exact structural denominator instead of estimating q^N.
        if self.denominator_key == rhs.denominator_key {
            return Self {
                numerator: self.numerator.checked_add(rhs.numerator)?,
                denominator: self.denominator.max(rhs.denominator),
                degree: self.degree.max(rhs.degree),
                denominator_key: self.denominator_key,
            }
            .bounded();
        }
        Self {
            numerator: self
                .numerator
                .checked_mul(rhs.denominator)?
                .checked_add(rhs.numerator.checked_mul(self.denominator)?)?,
            denominator: self.denominator.checked_mul(rhs.denominator)?,
            degree: self.degree.checked_add(rhs.degree)?,
            denominator_key: Self::denominator_product(self.denominator_key, rhs.denominator_key),
        }
        .bounded()
    }
}

fn scan(view: AtomView<'_>, epsilon: Symbol, nodes: &mut usize, depth: usize) -> Option<Shape> {
    *nodes = nodes.checked_add(1)?;
    if *nodes > MAX_NODES || depth > MAX_DEPTH {
        return None;
    }
    match view {
        AtomView::Num(n) => matches!(
            n.get_coeff_view(),
            CoefficientView::Natural(..) | CoefficientView::Large(..)
        )
        .then_some(Shape::leaf()),
        // An epsilon-free symbol is a coefficient, even when it has an
        // evaluator. The regulator itself remains the first polynomial variable.
        AtomView::Var(_) => Some(Shape::leaf()),
        AtomView::Fun(function) => {
            if view.contains_symbol(epsilon) {
                return None;
            }
            // Keep the whole epsilon-free call as a coefficient. Its analytic
            // hooks are irrelevant to expansion in epsilon; inspect arguments
            // only for exactness and the existing bounded syntax policy.
            for argument in function {
                scan(argument, epsilon, nodes, depth + 1)?;
            }
            Some(Shape::leaf())
        }
        AtomView::Add(sum) => {
            let mut shape = None;
            for term in sum {
                let next = scan(term, epsilon, nodes, depth + 1)?;
                shape = Some(match shape {
                    None => next,
                    Some(s) => Shape::add(s, next)?,
                });
            }
            shape
        }
        AtomView::Mul(product) => {
            let mut shape = Shape::one();
            for term in product {
                shape = shape.mul(scan(term, epsilon, nodes, depth + 1)?)?;
            }
            Some(shape)
        }
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            let mut shape = scan(base, epsilon, nodes, depth + 1)?;
            // Reject floats even in an epsilon-free exponent.
            scan(exponent, epsilon, nodes, depth + 1)?;
            // The native converter extracts integer content even from formal
            // fractional/symbolic exponents. Keep them out of this bounded lane.
            let power = crate::get_integer_from_atom(exponent)?;
            if power.unsigned_abs() > 64 {
                return None;
            }
            if power < 0 {
                std::mem::swap(&mut shape.numerator, &mut shape.denominator);
                shape.denominator_key = Some(base.to_owned());
            }
            let mut result = Shape::one();
            for _ in 0..power.unsigned_abs() {
                result = result.mul(shape.clone())?;
            }
            Some(result)
        }
    }
}

// A specialized preflight for Laurent monomials in regulator-free spectator
// atoms, with rational Q[epsilon] coefficients. Repeated/scaled univariate
// denominators share their exact LCM instead of multiplying support estimates.
#[derive(Clone)]
struct UnivariateShape {
    spectator_support: usize,
    // Exact Laurent-monomial support while small. Beyond the key cap retain
    // the existing conservative support count, never relax MAX_SUPPORT.
    spectator_keys: Option<Arc<std::collections::HashSet<Atom>>>,
    numerator_degree: usize,
    spectator_degree: usize,
    denominator: Polynomial,
    only_epsilon: bool,
    height: usize,
}

impl UnivariateShape {
    fn singleton_key(key: Atom) -> Option<Arc<std::collections::HashSet<Atom>>> {
        Some(Arc::new(std::collections::HashSet::from([key])))
    }

    fn combine_keys(
        &self,
        rhs: &Self,
        multiply: bool,
    ) -> Option<Arc<std::collections::HashSet<Atom>>> {
        let (Some(left), Some(right)) = (&self.spectator_keys, &rhs.spectator_keys) else {
            return None;
        };
        let mut keys = std::collections::HashSet::new();
        if multiply {
            // Do not allocate or enumerate a large Cartesian product merely
            // to refine the bound. The conservative count remains available.
            if left.len().checked_mul(right.len())? > MAX_SPECTATOR_KEYS {
                return None;
            }
            for a in left.iter() {
                for b in right.iter() {
                    keys.insert(a * b);
                    if keys.len() > MAX_SPECTATOR_KEYS {
                        return None;
                    }
                }
            }
        } else {
            for key in left.iter().chain(right.iter()) {
                keys.insert(key.clone());
                if keys.len() > MAX_SPECTATOR_KEYS {
                    return None;
                }
            }
        }
        Some(Arc::new(keys))
    }

    fn bounded(self) -> Option<Self> {
        (self.numerator_degree <= MAX_DEGREE
            && self.spectator_degree <= MAX_DEGREE
            && self.denominator.degree(0) <= MAX_DEGREE as u32
            && self
                .spectator_support
                .checked_mul(self.numerator_degree.checked_add(1)?)?
                <= MAX_SUPPORT)
            .then_some(self)
    }

    fn mul(self, rhs: Self) -> Option<Self> {
        // Bound the univariate product before allocating it.
        if self
            .denominator
            .degree(0)
            .checked_add(rhs.denominator.degree(0))?
            > MAX_DEGREE as u32
        {
            return None;
        }
        let conservative_support = self.spectator_support.checked_mul(rhs.spectator_support)?;
        let numerator_degree = self.numerator_degree.checked_add(rhs.numerator_degree)?;
        let spectator_degree = self.spectator_degree.checked_add(rhs.spectator_degree)?;
        if numerator_degree > MAX_DEGREE || spectator_degree > MAX_DEGREE {
            return None;
        }
        let spectator_keys = self.combine_keys(&rhs, true);
        let spectator_support = spectator_keys
            .as_ref()
            .map_or(conservative_support, |keys| keys.len());
        if spectator_support.checked_mul(numerator_degree.checked_add(1)?)? > MAX_SUPPORT {
            return None;
        }
        Self {
            spectator_support,
            spectator_keys,
            numerator_degree,
            spectator_degree,
            denominator: &self.denominator * &rhs.denominator,
            only_epsilon: self.only_epsilon && rhs.only_epsilon,
            height: self.height.max(rhs.height),
        }
        .bounded()
    }

    fn add(self, rhs: Self) -> Option<Self> {
        let left_degree = self.denominator.degree(0) as usize;
        let right_degree = rhs.denominator.degree(0) as usize;
        let denominator = if self.denominator == rhs.denominator {
            self.denominator.clone()
        } else {
            let gcd = self.denominator.gcd(&rhs.denominator);
            let cofactor = self.denominator.try_div(&gcd)?;
            if cofactor.degree(0).checked_add(rhs.denominator.degree(0))? > MAX_DEGREE as u32 {
                return None;
            }
            &cofactor * &rhs.denominator
        };
        let degree = denominator.degree(0) as usize;
        let only_epsilon = self.only_epsilon && rhs.only_epsilon;
        let conservative_support = if only_epsilon {
            1
        } else {
            self.spectator_support.checked_add(rhs.spectator_support)?
        };
        let spectator_keys = self.combine_keys(&rhs, false);
        let spectator_support = spectator_keys
            .as_ref()
            .map_or(conservative_support, |keys| keys.len());
        Self {
            spectator_support,
            spectator_keys,
            numerator_degree: self
                .numerator_degree
                .checked_add(degree.checked_sub(left_degree)?)?
                .max(
                    rhs.numerator_degree
                        .checked_add(degree.checked_sub(right_degree)?)?,
                ),
            spectator_degree: self.spectator_degree.max(rhs.spectator_degree),
            denominator,
            only_epsilon,
            height: self.height.max(rhs.height),
        }
        .bounded()
    }
}

struct UnivariateScan {
    epsilon: Symbol,
    one: Polynomial,
    // Repeated denominator polynomials and tensor leaves are scanned once.
    // The original unique-node budget and a separate total-visit budget bound
    // both cached state and traversal work; neither depends on cancellation.
    memo: std::collections::HashMap<Atom, UnivariateShape>,
    visits: usize,
}

impl UnivariateScan {
    fn new(epsilon: Symbol) -> Self {
        Self {
            epsilon,
            one: Polynomial::new(&Q, None, Arc::new(vec![PolyVariable::Symbol(epsilon)])).one(),
            memo: std::collections::HashMap::new(),
            visits: 0,
        }
    }

    fn leaf(&self, numerator_degree: usize, only_epsilon: bool) -> UnivariateShape {
        UnivariateShape {
            spectator_support: 1,
            spectator_keys: if only_epsilon {
                UnivariateShape::singleton_key(Atom::one())
            } else {
                None
            },
            numerator_degree,
            spectator_degree: usize::from(!only_epsilon),
            denominator: self.one.clone(),
            only_epsilon,
            height: 0,
        }
    }

    fn scan(&mut self, view: AtomView<'_>, depth: usize) -> Option<UnivariateShape> {
        self.visits = self.visits.checked_add(1)?;
        if depth > MAX_DEPTH || self.visits > MAX_WORK {
            return None;
        }
        let atom = view.to_owned();
        if let Some(shape) = self.memo.get(&atom) {
            return (depth.checked_add(shape.height)? <= MAX_DEPTH + 1).then(|| shape.clone());
        }
        if self.memo.len() >= MAX_NODES {
            return None;
        }
        let mut shape = match view {
            AtomView::Num(number) => {
                if !matches!(
                    number.get_coeff_view(),
                    CoefficientView::Natural(..) | CoefficientView::Large(..)
                ) {
                    return None;
                }
                self.leaf(0, true)
            }
            AtomView::Var(variable) => {
                let is_epsilon = variable.get_symbol() == self.epsilon;
                let mut shape = self.leaf(usize::from(is_epsilon), is_epsilon);
                if !is_epsilon {
                    shape.spectator_keys = UnivariateShape::singleton_key(atom.clone());
                }
                shape
            }
            AtomView::Fun(function) => {
                if view.contains_symbol(self.epsilon) {
                    return None;
                }
                // Keep the exactness checks for opaque function arguments.
                let mut result = self.leaf(0, false);
                result.spectator_keys = UnivariateShape::singleton_key(atom.clone());
                for argument in function {
                    result.height = result.height.max(self.scan(argument, depth + 1)?.height);
                }
                result
            }
            AtomView::Add(sum) => {
                let mut terms = sum.iter();
                let mut result = self.scan(terms.next()?, depth + 1)?;
                for term in terms {
                    result = result.add(self.scan(term, depth + 1)?)?;
                }
                result
            }
            AtomView::Mul(product) => {
                let mut result = self.leaf(0, true);
                for factor in product {
                    result = result.mul(self.scan(factor, depth + 1)?)?;
                }
                result
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let mut shape = self.scan(base, depth + 1)?;
                let height = shape.height.max(self.scan(exponent, depth + 1)?.height);
                let power = crate::get_integer_from_atom(exponent)?;
                if power.unsigned_abs() > 64 {
                    return None;
                }
                if power < 0 {
                    if shape.only_epsilon {
                        // This conversion is only univariate and has already
                        // passed the support and degree bounds recursively.
                        let rational: RationalPolynomial<IntegerRing, u32> = base
                            .try_to_rational_polynomial(&Q, &Z, self.one.variables().clone())
                            .ok()?;
                        if rational.numerator.variables().len() != 1 || rational.numerator.is_zero()
                        {
                            return None;
                        }
                        let denominator = rational.numerator.map_coeff(|c| c.into(), Q);
                        let leading = denominator.lcoeff();
                        shape = UnivariateShape {
                            spectator_support: 1,
                            spectator_keys: UnivariateShape::singleton_key(Atom::one()),
                            numerator_degree: rational.denominator.degree(0) as usize,
                            spectator_degree: 0,
                            denominator: denominator.div_coeff(&leading),
                            only_epsilon: true,
                            height,
                        }
                        .bounded()?;
                    } else if !base.contains_symbol(self.epsilon) && shape.spectator_support == 1 {
                        // Inverting an epsilon-free monomial still contributes
                        // one Laurent monomial. A spectator sum must fall back.
                        let spectator_degree = shape.spectator_degree;
                        let spectator_keys = shape.spectator_keys.as_ref().map(|keys| {
                            Arc::new(keys.iter().map(|key| key.pow(Atom::num(-1))).collect())
                        });
                        shape = self.leaf(0, false);
                        shape.spectator_degree = spectator_degree;
                        shape.spectator_keys = spectator_keys;
                    } else {
                        return None;
                    }
                }
                // Keep the original power-expansion budget. Tracking repeated
                // monomials in sums must not turn a large polynomial power into
                // an admitted input merely because its expanded keys overlap.
                if shape
                    .spectator_support
                    .checked_pow(power.unsigned_abs() as u32)?
                    > MAX_SUPPORT
                {
                    return None;
                }
                let mut result = self.leaf(0, true);
                result.height = height;
                for _ in 0..power.unsigned_abs() {
                    result = result.mul(shape.clone())?;
                }
                result
            }
        }
        .bounded()?;
        shape.height = shape.height.checked_add(1)?;
        if depth.checked_add(shape.height)? > MAX_DEPTH + 1 {
            return None;
        }
        // Recursive calls may have filled the cache since the entry check.
        if self.memo.len() >= MAX_NODES {
            return None;
        }
        self.memo.insert(atom, shape.clone());
        Some(shape)
    }
}

/// The exact series of a rational function of ε with spectator coefficients:
/// one [`RationalPart`], or — when the whole input exceeds the part's support
/// bound (spectator monomials × numerator degree, [`MAX_SUPPORT`]) — one part
/// per spectator monomial ([`Self::grouped`]), the series being the sum of the
/// parts' series.
pub(super) struct RationalSeries {
    parts: Vec<RationalPart>,
}

/// Groups of a split input at most (one per distinct spectator monomial).
const MAX_GROUPS: usize = 65_536;

impl RationalSeries {
    #[cfg(test)]
    pub(super) fn validate_spectators(expression: &Atom, epsilon: Symbol) -> Option<()> {
        RationalPart::validate_spectators(expression, epsilon)
    }

    pub(super) fn new(expression: &Atom, epsilon: Symbol) -> Option<Self> {
        match RationalPart::new(expression, epsilon) {
            Some(part) => Some(Self { parts: vec![part] }),
            None => Self::grouped(expression, epsilon),
        }
    }

    /// The input as a sum over its spectator monomials `s`: `Σ_s s·R_s(ε)`,
    /// every `R_s` a univariate rational function of ε (its own denominators),
    /// each prepared as one [`RationalPart`] with every guard of the single
    /// path (a group the part refuses refuses the whole). Admitted only when
    /// the whole input passes the spectator-coordinate and semantic-zero
    /// checks ([`RationalPart::validate_spectators`]). The valuation is the
    /// least of the parts' valuations: a lower bound of the sum's (spectator
    /// monomials may still be related, `ζ(2)³` and `ζ(6)`), so a caller that
    /// expands its other factors that much deeper stays exact.
    fn grouped(expression: &Atom, epsilon: Symbol) -> Option<Self> {
        RationalPart::validate_spectators(expression, epsilon)?;
        // an overall numeric (complex) phase stays a factor of every group
        let factors = match expression.as_view() {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            view => vec![view],
        };
        let (phase, rest): (Vec<_>, Vec<_>) =
            factors.into_iter().partition(|factor| matches!(factor, AtomView::Num(_)));
        let phase = Atom::mul_many(phase);
        let expanded = Atom::mul_many(rest).expand();
        let summands = match expanded.as_view() {
            AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
            view => vec![view],
        };
        let mut index = std::collections::HashMap::<Atom, usize>::new();
        let mut groups: Vec<(Atom, Vec<Atom>)> = Vec::new();
        for summand in summands {
            let factors = match summand {
                AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
                view => vec![view],
            };
            let (key, coefficient): (Vec<_>, Vec<_>) = factors
                .into_iter()
                .partition(|f| !matches!(f, AtomView::Num(_)) && !f.contains_symbol(epsilon));
            let key = Atom::mul_many(key);
            let coefficient = Atom::mul_many(coefficient);
            match index.get(&key) {
                Some(&i) => groups[i].1.push(coefficient),
                None => {
                    if groups.len() >= MAX_GROUPS {
                        return None;
                    }
                    index.insert(key.clone(), groups.len());
                    groups.push((key, vec![coefficient]));
                }
            }
        }
        if groups.len() < 2 {
            return None;
        }
        let parts = groups
            .into_iter()
            .map(|(key, coefficients)| {
                RationalPart::new(&(&phase * &key * Atom::add_many(coefficients.iter())), epsilon)
            })
            .collect::<Option<Vec<_>>>()?;
        Some(Self { parts })
    }

    pub(super) fn valuation(&self) -> Option<i64> {
        self.parts.iter().filter_map(|part| part.valuation()).min()
    }

    pub(super) fn series(&self, order: i64) -> Option<Laurent> {
        let Some(valuation) = self.valuation() else {
            return Some(Laurent {
                valuation: 0,
                coefficients: Vec::new(),
            });
        };
        if let [part] = self.parts.as_slice() {
            return part.series(order);
        }
        let length = usize::try_from(order.checked_sub(valuation)?.checked_add(1)?.max(0)).ok()?;
        let mut coefficients = vec![Atom::Zero; length];
        for part in &self.parts {
            let series = part.series(order)?;
            for (i, coefficient) in series.coefficients.iter().enumerate() {
                let at = series.valuation + i as i64 - valuation;
                if let Ok(at) = usize::try_from(at)
                    && at < length
                {
                    coefficients[at] = &coefficients[at] + coefficient;
                }
            }
        }
        Some(Laurent {
            valuation,
            coefficients: coefficients.into_iter().map(|c| c.expand()).collect(),
        })
    }
}

/// One exact rational series (module documentation).
pub(super) struct RationalPart {
    numerator: Vec<Polynomial>,
    denominator: Vec<Rational>,
    prefactor: Atom,
    valuation: Option<i64>,
    numerator_start: usize,
}

impl RationalPart {
    /// Check original spectator coordinates and inverse bases before rational
    /// conversion can enlarge the variable map or cancel denominator factors.
    fn validate_spectators(expression: &Atom, epsilon: Symbol) -> Option<()> {
        let mut variables = std::collections::HashSet::from([Atom::var(epsilon)]);
        let mut denominators = std::collections::HashSet::new();
        let mut valid = true;
        let mut visits = 0usize;
        expression.visitor(&mut |part| {
            visits += 1;
            if !valid || visits > MAX_WORK {
                valid = false;
                return false;
            }
            match part {
                AtomView::Var(_) | AtomView::Fun(_) => {
                    // An epsilon-free call is one coordinate. Its metadata
                    // arguments were checked by the syntax scan, not variables
                    // of the polynomial coefficient ring.
                    variables.insert(part.to_owned());
                    valid = variables.len() <= MAX_VARIABLES;
                    return false;
                }
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    if crate::get_integer_from_atom(exponent).is_some_and(|power| power < 0)
                        && denominators.insert(base.to_owned())
                    {
                        let mut semantic_coefficients = false;
                        base.visitor(&mut |part| {
                            semantic_coefficients |= match part {
                                AtomView::Fun(_) => true,
                                AtomView::Var(variable) if variable.get_symbol() != epsilon => {
                                    let symbol = variable.get_symbol();
                                    symbol.get_normalization_function().is_some()
                                        || symbol.get_evaluation_info().is_some()
                                        || symbol.get_series_function().is_some()
                                        || symbol.get_derivative_function().is_some()
                                }
                                _ => false,
                            };
                            !semantic_coefficients
                        });
                        // A formally common factor may be identically zero
                        // after restoring function semantics. Reject it before
                        // polynomial GCD cancellation can hide the division.
                        // Pure Q[epsilon] denominators do not need this proof.
                        if semantic_coefficients && base.check_zero().is_true() {
                            valid = false;
                            return false;
                        }
                    }
                }
                _ => {}
            }
            valid
        });
        valid.then_some(())
    }

    pub(super) fn new(expression: &Atom, epsilon: Symbol) -> Option<Self> {
        if scan(expression.as_view(), epsilon, &mut 0, 0).is_none() {
            UnivariateScan::new(epsilon).scan(expression.as_view(), 0)?;
        }
        Self::validate_spectators(expression, epsilon)?;
        // A global Gaussian rational phase belongs to the epsilon-free factor.
        // Mixed complex sums currently fall back rather than changing domains.
        let factors = match expression.as_view() {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            view => vec![view],
        };
        let (phase, rest): (Vec<_>, Vec<_>) = factors.into_iter().partition(|factor| {
            matches!(factor, AtomView::Num(_)) && factor.has_complex_coefficients()
        });
        let mut phase = Atom::mul_many(phase);
        let mut expression = Atom::mul_many(rest);
        if expression.as_view().has_complex_coefficients() {
            // Reuse the master finalizer's exact common-phase policy for a
            // phase distributed across an expanded sum. Do not invent a
            // different meaning for opaque complex labels or finite literals.
            let common =
                crate::rustred_evaluation::masters::common_exact_numeric_phase(&expression);
            if common.is_one() {
                return None;
            }
            expression = (&expression / &common).expand();
            phase *= common;
        }
        let initial = Arc::new(vec![PolyVariable::Symbol(epsilon)]);
        let rational: RationalPolynomial<IntegerRing, u32> = expression
            .try_to_rational_polynomial(&Q, &Z, initial)
            .ok()?;
        let variables = rational.numerator.variables();
        // Conversion can introduce opaque variables; none may hide epsilon.
        if variables.len() > MAX_VARIABLES
            || variables.first() != Some(&PolyVariable::Symbol(epsilon))
            || variables.as_ref() != rational.denominator.variables().as_ref()
            || variables
                .iter()
                .skip(1)
                .any(|v| Atom::from(v.clone()).contains_symbol(epsilon))
            || rational.numerator.nterms() > MAX_SUPPORT
            || rational.denominator.nterms() > MAX_SUPPORT
            || rational.numerator.degree(0) > MAX_DEGREE as u32
            || rational.denominator.degree(0) > MAX_DEGREE as u32
        {
            return None;
        }
        // Native conversion clears denominators into integer polynomials.
        // Embed their coefficients in Q for the exact division recurrence.
        let mut numerator = rational
            .numerator
            .map_coeff(|coefficient| coefficient.into(), Q)
            .to_univariate(0)
            .coefficients()
            .to_vec();
        let denominator = rational
            .denominator
            .map_coeff(|coefficient| coefficient.into(), Q)
            .to_univariate(0);
        let denominator_start = denominator
            .coefficients()
            .iter()
            .position(|c| !c.is_zero())?;
        let common = &denominator.coefficients()[denominator_start];
        let common_atom = common.to_expression();
        if common_atom.check_zero().is_true() {
            return None;
        }
        // This polynomial proportionality is a sufficient identity after
        // restoring the spectator atoms, whether or not those atoms satisfy
        // further relations. The recurrence only multiplies them by rationals;
        // restored zero proofs below determine the actual leading epsilon order.
        let mut numeric_denominator = Vec::new();
        for coefficient in &denominator.coefficients()[denominator_start..] {
            let ratio = coefficient.lcoeff() / common.lcoeff();
            if *coefficient != common.clone().mul_coeff(ratio.clone()) {
                return None;
            }
            numeric_denominator.push(ratio);
        }
        // As in Symbolica's normalize_trailing_coefficient, only a proof of
        // zero can shift the valuation. Keep opaque function semantics intact.
        let mut numerator_start = 0;
        while numerator_start < numerator.len() {
            let coefficient = &numerator[numerator_start];
            if !coefficient.is_zero() && !coefficient.to_expression().check_zero().is_true() {
                break;
            }
            numerator[numerator_start] = coefficient.zero();
            numerator_start += 1;
        }
        let valuation = (numerator_start < numerator.len())
            .then(|| numerator_start as i64 - denominator_start as i64);
        Some(Self {
            numerator,
            denominator: numeric_denominator,
            prefactor: phase / common_atom,
            valuation,
            numerator_start,
        })
    }

    pub(super) fn valuation(&self) -> Option<i64> {
        self.valuation
    }

    pub(super) fn series(&self, order: i64) -> Option<Laurent> {
        let Some(valuation) = self.valuation else {
            return Some(Laurent {
                valuation: 0,
                coefficients: Vec::new(),
            });
        };
        let length = usize::try_from(order.checked_sub(valuation)?.checked_add(1)?.max(0)).ok()?;
        let support = self.numerator.iter().map(|c| c.nterms()).sum::<usize>();
        if length > MAX_SERIES || length.checked_mul(length)?.checked_mul(support)? > MAX_WORK {
            return None;
        }
        let prototype = &self.numerator[self.numerator_start];
        let mut coefficients: Vec<Polynomial> = Vec::with_capacity(length);
        for n in 0..length {
            let mut value = self
                .numerator
                .get(self.numerator_start + n)
                .cloned()
                .unwrap_or_else(|| prototype.zero());
            // q0 == 1 by construction. The other qk are exact rationals.
            for k in 1..=n.min(self.denominator.len() - 1) {
                value = &value
                    - &coefficients[n - k]
                        .clone()
                        .mul_coeff(self.denominator[k].clone());
            }
            coefficients.push(value);
        }
        Some(Laurent {
            valuation,
            coefficients: coefficients
                .into_iter()
                .map(|coefficient| (coefficient.to_expression() * &self.prefactor).expand())
                .collect(),
        })
    }
}

#[cfg(test)]
mod shared_denominator_tests {
    use super::*;
    use crate::utils::vakint_macros::{vk_parse, vk_symbol};

    #[test]
    fn spectator_coordinate_limit_precedes_common_denominator_conversion() {
        crate::Vakint::initialize_vakint_symbols();
        let epsilon = vk_symbol!("ε");
        let head = vk_symbol!("coordinate_limit_spectator");
        let input = Atom::add_many((0..MAX_VARIABLES).map(|i| head.call(i as i64).pow(-1)));
        // The shared-denominator shape is small despite the distinct inverse
        // coordinates. Reject before clearing their common denominator.
        assert!(
            UnivariateScan::new(epsilon)
                .scan(input.as_view(), 0)
                .is_some()
        );
        assert!(RationalSeries::validate_spectators(&input, epsilon).is_none());
        assert!(RationalSeries::new(&input, epsilon).is_none());
    }

    #[test]
    fn opaque_function_arguments_are_not_extra_spectator_coordinates() {
        crate::Vakint::initialize_vakint_symbols();
        let epsilon = vk_symbol!("ε");
        let arguments = (0..MAX_VARIABLES + 5)
            .map(|i| Atom::var(vk_symbol!(&format!("coordinate_metadata_{i}"))))
            .collect::<Vec<_>>();
        let spectator = vk_symbol!("coordinate_metadata").call_args(arguments);
        let e = Atom::var(epsilon);
        let input = &spectator / (Atom::one() - &e);
        let actual = RationalSeries::new(&input, epsilon)
            .unwrap()
            .series(2)
            .unwrap()
            .to_atom(epsilon);
        let expected = &spectator * (Atom::one() + &e + e.pow(2));
        assert!((actual - expected).expand().is_zero());
    }

    #[test]
    fn scaled_and_overlapping_epsilon_denominators_use_their_lcm() {
        crate::Vakint::initialize_vakint_symbols();
        let epsilon = vk_symbol!("ε");
        let expression = vk_parse!("a/(2-6*ε+4*ε^2)+b/(3-6*ε)+c/(5-15*ε+10*ε^2)").unwrap();
        let shape = UnivariateScan::new(epsilon)
            .scan(expression.as_view(), 0)
            .unwrap();
        assert_eq!(shape.denominator.degree(0), 2);
        assert_eq!(shape.spectator_support, 3);
        assert_eq!(shape.numerator_degree, 1);
        let prepared = RationalSeries::new(&expression, epsilon).unwrap();
        let actual = prepared.series(4).unwrap().to_atom(epsilon);
        let reference = expression
            .series(epsilon, Atom::Zero.as_view(), Rational::from(4))
            .unwrap()
            .to_atom();
        assert!((actual - reference).expand().together().is_zero());
    }

    #[test]
    fn repeated_scaled_kernels_no_longer_explode_the_support_estimate() {
        crate::Vakint::initialize_vakint_symbols();
        let epsilon = vk_symbol!("ε");
        let expression = Atom::add_many((1..=24).map(|n| {
            vk_parse!(&format!(
                "f({n})*ε^{}/({n}-{}*ε+{}*ε^2)",
                n % 5,
                3 * n,
                2 * n
            ))
            .unwrap()
        }));
        assert!(scan(expression.as_view(), epsilon, &mut 0, 0).is_none());
        let shape = UnivariateScan::new(epsilon)
            .scan(expression.as_view(), 0)
            .unwrap();
        assert_eq!(shape.denominator.degree(0), 2);
        assert_eq!(shape.spectator_support, 24);
        let prepared = RationalSeries::new(&expression, epsilon).unwrap();
        let actual = prepared.series(3).unwrap().to_atom(epsilon);
        let reference = expression
            .series(epsilon, Atom::Zero.as_view(), Rational::from(3))
            .unwrap()
            .to_atom();
        assert!((actual - reference).expand().together().is_zero());
    }

    #[test]
    fn bounded_univariate_scan_keeps_mixed_denominators_and_large_work_out() {
        crate::Vakint::initialize_vakint_symbols();
        let epsilon = vk_symbol!("ε");
        for expression in [
            "1/(a+ε)",
            "1/(a+b)",
            "sin(ε)",
            "1/(1-ε)^65",
            "(1+ε^64)^64",
            "(a+b)^32",
            "((a*b)^64*(c*d)^64)^64",
            "((a*b)^(-64)*(c*d)^(-64))^(-64)",
            "ε^(1/2)",
        ] {
            let input = vk_parse!(expression).unwrap();
            assert!(
                UnivariateScan::new(epsilon)
                    .scan(input.as_view(), 0)
                    .is_none(),
                "{expression}"
            );
        }
    }
}

#[cfg(test)]
mod spectator_support_tests {
    use super::*;
    use crate::utils::vakint_macros::{vk_parse, vk_symbol};

    #[test]
    fn repeated_spectator_monomials_are_counted_once_across_scaled_kernels() {
        crate::Vakint::initialize_vakint_symbols();
        let epsilon = vk_symbol!("ε");
        let expression = Atom::add_many((1..=600).map(|n| {
            let spectator = if n % 2 == 0 { "a" } else { "b" };
            vk_parse!(&format!(
                "{spectator}*ε^{}/({n}-{}*ε+{}*ε^2)",
                n % 64,
                3 * n,
                2 * n
            ))
            .unwrap()
        }));
        assert!(scan(expression.as_view(), epsilon, &mut 0, 0).is_none());
        let shape = UnivariateScan::new(epsilon)
            .scan(expression.as_view(), 0)
            .unwrap();
        let AtomView::Add(sum) = expression.as_view() else {
            panic!("expected an expanded sum");
        };
        assert!(
            sum.iter().count() * (shape.numerator_degree + 1) > MAX_SUPPORT,
            "fixture must exceed the old repeated-term support bound"
        );
        assert_eq!(shape.spectator_support, 2);
        assert_eq!(shape.spectator_keys.as_ref().unwrap().len(), 2);
        assert_eq!(shape.denominator.degree(0), 2);
        let prepared = RationalSeries::new(&expression, epsilon).unwrap();
        let actual = prepared.series(1).unwrap().to_atom(epsilon);
        let reference = expression
            .series(epsilon, Atom::Zero.as_view(), Rational::from(1))
            .unwrap()
            .to_atom();
        assert!((actual - reference).expand().together().is_zero());
    }

    #[test]
    fn spectator_support_keeps_laurent_powers_and_cross_products() {
        crate::Vakint::initialize_vakint_symbols();
        let epsilon = vk_symbol!("ε");
        let expression = vk_parse!("(a+b)*(a^(-1)+b^(-1))/(1-ε)").unwrap();
        let shape = UnivariateScan::new(epsilon)
            .scan(expression.as_view(), 0)
            .unwrap();
        assert_eq!(shape.spectator_support, 3);
        let expected: std::collections::HashSet<_> = ["1", "a/b", "b/a"]
            .map(|s| vk_parse!(s).unwrap())
            .into_iter()
            .collect();
        assert_eq!(shape.spectator_keys.as_deref(), Some(&expected));
        let actual = RationalSeries::new(&expression, epsilon)
            .unwrap()
            .series(2)
            .unwrap()
            .to_atom(epsilon);
        let reference = expression
            .series(epsilon, Atom::Zero.as_view(), Rational::from(2))
            .unwrap()
            .to_atom();
        assert!((actual - reference).expand().together().is_zero());
    }

    #[test]
    fn spectator_key_cap_retains_conservative_support_and_degree_guards() {
        crate::Vakint::initialize_vakint_symbols();
        let epsilon = vk_symbol!("ε");
        let expression =
            Atom::add_many((0..=64).map(|n| vk_parse!(&format!("a^{n}*b^{}", 64 - n)).unwrap()));
        let shape = UnivariateScan::new(epsilon)
            .scan(expression.as_view(), 0)
            .unwrap();
        assert!(shape.spectator_keys.is_none());
        assert_eq!(shape.spectator_support, 65);
        for expression in [
            "((a*b)^64*(c*d)^64)^64",
            "((a*b)^(-64)*(c*d)^(-64))^(-64)",
            "(a+b+c)^64",
        ] {
            assert!(
                UnivariateScan::new(epsilon)
                    .scan(vk_parse!(expression).unwrap().as_view(), 0)
                    .is_none(),
                "{expression}"
            );
        }
    }
}

#[cfg(test)]
mod decline_probe {
    use super::*;

    /// Diagnostics (ignored): the stages of [`RationalSeries::new`] on a
    /// declined input dumped by `VAKINT_SERIES_DECLINE_LOG` (path in
    /// `VAKINT_RATIONAL_PROBE`).
    #[test]
    #[ignore]
    fn probe_a_declined_input() {
        crate::Vakint::initialize_vakint_symbols();
        let path = std::env::var("VAKINT_RATIONAL_PROBE").unwrap();
        let text = std::fs::read_to_string(path).unwrap();
        let line = text.lines().find(|l| l.starts_with("expression=")).unwrap();
        let expression = symbolica::try_parse!(&line["expression=".len()..]).unwrap();
        let epsilon = symbolica::symbol!("gammalooprs::ε");
        let mut nodes = 0;
        eprintln!("PROBE scan: {}", scan(expression.as_view(), epsilon, &mut nodes, 0).is_some());
        eprintln!("PROBE scan nodes {nodes}");
        let mut univariate = UnivariateScan::new(epsilon);
        let shape = univariate.scan(expression.as_view(), 0);
        eprintln!(
            "PROBE univariate: {} memo {} visits {}",
            shape.is_some(),
            univariate.memo.len(),
            univariate.visits
        );
        if let AtomView::Mul(product) = expression.as_view() {
            for factor in product.iter() {
                let mut u = UnivariateScan::new(epsilon);
                let s = u.scan(factor, 0);
                eprintln!(
                    "PROBE factor bytes {} ok {} memo {} {}",
                    factor.get_byte_size(),
                    s.is_some(),
                    u.memo.len(),
                    s.map(|s| format!(
                        "support {} keys {:?} numdeg {} specdeg {} dendeg {}",
                        s.spectator_support,
                        s.spectator_keys.as_ref().map(|k| k.len()),
                        s.numerator_degree,
                        s.spectator_degree,
                        s.denominator.degree(0)
                    ))
                    .unwrap_or_default()
                );
                if let AtomView::Add(sum) = factor {
                    eprintln!("PROBE  sum terms {}", sum.get_nargs());
                    let mut u = UnivariateScan::new(epsilon);
                    let mut acc: Option<UnivariateShape> = None;
                    for (i, term) in sum.iter().enumerate() {
                        let Some(t) = u.scan(term, 1) else {
                            eprintln!("PROBE  term {i} declined alone (memo {})", u.memo.len());
                            break;
                        };
                        acc = match acc {
                            None => Some(t),
                            Some(a) => match a.clone().add(t.clone()) {
                                Some(x) => Some(x),
                                None => {
                                    eprintln!(
                                        "PROBE  add declined at term {i}: acc support {} numdeg {} dendeg {}; term support {} numdeg {} dendeg {}",
                                        a.spectator_support, a.numerator_degree, a.denominator.degree(0),
                                        t.spectator_support, t.numerator_degree, t.denominator.degree(0)
                                    );
                                    break;
                                }
                            },
                        };
                    }
                }
            }
        }
        eprintln!("PROBE validate: {:?}", RationalSeries::validate_spectators(&expression, epsilon).is_some());
    }
}
