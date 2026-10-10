//! Bounded, exact finite Laurent jets for the offline soft-overlap audit.
//! No native Atom::series, statistical zero test, or numerator expansion.

use super::soft_energy_certificate::{EnergyZeroOutcome, SoftEnergyAlgebra};
use crate::utils::GS;
use color_eyre::eyre::{Result, ensure, eyre};
use std::{collections::BTreeMap, sync::Arc};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    domains::{atom::AtomField, rational::Rational},
    poly::series::Series,
};

type Coefficients = BTreeMap<i64, Atom>;
const MAX_ORDER: i64 = 32;
const MAX_POWER: u64 = 16;
const MAX_BYTES: usize = 64 * 1024 * 1024;
const MAX_CACHE_BYTES: usize = 256 * 1024 * 1024;

pub(crate) struct ExactSoftJet {
    lambda: Symbol,
    endpoint: i64,
    lower: BTreeMap<Atom, i64>,
    leading: BTreeMap<Atom, (i64, Atom)>,
    coefficients: BTreeMap<(Atom, i64), Coefficients>,
    algebra: SoftEnergyAlgebra,
    cache_bytes: usize,
    cache_limit: usize,
}

impl ExactSoftJet {
    pub(crate) fn new(lambda: Symbol, endpoint: i64) -> Self {
        Self {
            lambda,
            endpoint,
            lower: BTreeMap::new(),
            leading: BTreeMap::new(),
            coefficients: BTreeMap::new(),
            algebra: SoftEnergyAlgebra::default(),
            cache_bytes: 0,
            cache_limit: MAX_CACHE_BYTES,
        }
    }

    pub(crate) fn coefficients(&mut self, atom: &Atom) -> Result<Coefficients> {
        self.through(atom, self.endpoint)
    }

    /// Adapt the exact finite map to the retained-numerator precision owner.
    /// Even an empty map carries only O(lambda^(depth+1)), never exact zero.
    pub(crate) fn series(
        atom: &Atom,
        lambda: Symbol,
        depth: Rational,
    ) -> Result<Series<AtomField>> {
        ensure!(
            depth.is_integer(),
            "unproven: noninteger exact soft-jet endpoint"
        );
        let endpoint = depth
            .numerator()
            .to_i64()
            .ok_or_else(|| eyre!("soft-jet endpoint overflow"))?;
        let coefficients = Self::new(lambda, endpoint).coefficients(atom)?;
        let first = coefficients.keys().next().copied().unwrap_or(endpoint + 1);
        let field = AtomField {
            statistical_zero_test: false,
            ..Default::default()
        };
        let template = Series::new(
            &field,
            None,
            Arc::new(lambda.into()),
            Atom::Zero,
            Rational::from((endpoint - first + 1).max(1)),
        );
        let mut terms = coefficients.into_iter();
        let mut series = if let Some((power, coefficient)) = terms.next() {
            template.monomial(coefficient, power.into())
        } else {
            template.monomial(Atom::one(), (endpoint + 1).into())
        };
        for (power, coefficient) in terms {
            series = &series + &template.monomial(coefficient, power.into());
        }
        series.truncate_absolute_order((endpoint + 1).into());
        ensure!(
            series.absolute_order() == endpoint + 1,
            "exact soft-jet adapter lost absolute precision"
        );
        Ok(series)
    }

    fn lower_bound(&mut self, atom: &Atom) -> Result<i64> {
        if let Some(lower) = self.lower.get(atom) {
            return Ok(*lower);
        }
        let lower = if !atom.contains_symbol(self.lambda) {
            0
        } else {
            match atom.as_view() {
                AtomView::Var(variable) if variable.get_symbol() == self.lambda => 1,
                AtomView::Add(sum) => {
                    let mut lower = i64::MAX;
                    for term in sum {
                        lower = lower.min(self.lower_bound(&term.to_owned())?);
                    }
                    lower
                }
                AtomView::Mul(product) => {
                    let mut lower = 0;
                    for factor in product {
                        lower += self.lower_bound(&factor.to_owned())?;
                    }
                    lower
                }
                AtomView::Pow(power) => {
                    let exponent = Self::power(power.get_exp())?;
                    let base = power.get_base().to_owned();
                    exponent
                        * if exponent < 0 {
                            self.leading(&base)?.0
                        } else {
                            self.lower_bound(&base)?
                        }
                }
                AtomView::Fun(call)
                    if call.get_symbol() == GS.on_shell_energy && call.get_nargs() == 2 =>
                {
                    let valuation = self.leading(&call.get(1).to_owned())?.0;
                    ensure!(
                        valuation % 2 == 0,
                        "unproven: half-integer energy valuation"
                    );
                    valuation / 2
                }
                _ => return Err(eyre!("unproven: unknown Taylor-dependent atom {atom}")),
            }
        };
        ensure!(
            lower.abs() <= MAX_ORDER,
            "unproven: Laurent lower-bound budget"
        );
        let bytes = atom
            .as_view()
            .get_byte_size()
            .saturating_add(size_of::<i64>());
        if self.reserve_cache(bytes) {
            self.lower.insert(atom.clone(), lower);
        }
        Ok(lower)
    }

    fn leading(&mut self, atom: &Atom) -> Result<(i64, Atom)> {
        if let Some(leading) = self.leading.get(atom) {
            return Ok(leading.clone());
        }
        let lower = self.lower_bound(atom)?;
        for extra in [0, 1, 2, 4, 8, 16, MAX_ORDER] {
            for (power, coefficient) in self.through(atom, lower + extra)? {
                if matches!(
                    self.algebra.scalar_zero(&coefficient)?,
                    EnergyZeroOutcome::Zero(_)
                ) {
                    continue;
                }
                ensure!(
                    self.leading_nonzero(&coefficient)?,
                    "unproven: no exact nonzero witness for leading denominator coefficient {coefficient}"
                );
                let bytes = atom
                    .as_view()
                    .get_byte_size()
                    .saturating_add(coefficient.as_view().get_byte_size())
                    .saturating_add(size_of::<i64>());
                if self.reserve_cache(bytes) {
                    self.leading
                        .insert(atom.clone(), (power, coefficient.clone()));
                }
                return Ok((power, coefficient));
            }
        }
        Err(eyre!(
            "unproven: leading coefficient not found within exact finite probe"
        ))
    }

    fn through(&mut self, atom: &Atom, endpoint: i64) -> Result<Coefficients> {
        ensure!(
            endpoint.abs() <= 2 * MAX_ORDER,
            "unproven: finite-jet endpoint budget"
        );
        let cache_key = (atom.clone(), endpoint);
        if let Some(coefficients) = self.coefficients.get(&cache_key) {
            return Ok(coefficients.clone());
        }
        let lower = self.lower_bound(atom)?;
        let mut result = BTreeMap::new();
        if endpoint >= lower {
            if !atom.contains_symbol(self.lambda) {
                if !atom.is_zero() {
                    result.insert(0, atom.clone());
                }
            } else {
                match atom.as_view() {
                    AtomView::Var(variable) if variable.get_symbol() == self.lambda => {
                        result.insert(1, Atom::one());
                    }
                    AtomView::Add(sum) => {
                        for term in sum {
                            for (power, coefficient) in self.through(&term.to_owned(), endpoint)? {
                                *result.entry(power).or_insert(Atom::Zero) += coefficient;
                            }
                        }
                    }
                    AtomView::Mul(product) => {
                        let factors = product
                            .iter()
                            .map(|factor| factor.to_owned())
                            .collect::<Vec<_>>();
                        result = self.product(&factors, endpoint)?;
                    }
                    AtomView::Pow(power) => {
                        let exponent = Self::power(power.get_exp())?;
                        let base = power.get_base().to_owned();
                        if exponent >= 0 {
                            result = self.product(&vec![base; exponent as usize], endpoint)?;
                        } else {
                            let (valuation, _) = self.leading(&base)?;
                            let count = -exponent;
                            let inverse =
                                self.inverse(&base, endpoint + (count - 1) * valuation)?;
                            result.insert(0, Atom::one());
                            for completed in 1..=count {
                                result = Self::convolve(
                                    &result,
                                    &inverse,
                                    endpoint + (count - completed) * valuation,
                                )?;
                            }
                        }
                    }
                    AtomView::Fun(call)
                        if call.get_symbol() == GS.on_shell_energy && call.get_nargs() == 2 =>
                    {
                        result = self.energy(&call.get(1).to_owned(), endpoint)?;
                    }
                    _ => return Err(eyre!("unproven: unsupported exact soft jet {atom}")),
                }
            }
        }
        result.retain(|_, coefficient| !coefficient.is_zero());
        Self::bounded(&result)?;
        let bytes = result.values().fold(
            atom.as_view()
                .get_byte_size()
                .saturating_add(size_of::<i64>()),
            |bytes, coefficient| {
                bytes
                    .saturating_add(coefficient.as_view().get_byte_size())
                    .saturating_add(size_of::<i64>())
            },
        );
        if self.reserve_cache(bytes) {
            self.coefficients.insert(cache_key, result.clone());
        }
        Ok(result)
    }

    /// Memoization is an optimization, not a certificate-size limit. Account
    /// keys as well as values across all three tables. Eviction is safe during
    /// recursion because in-flight coefficients and bounds are owned locally.
    fn reserve_cache(&mut self, bytes: usize) -> bool {
        if bytes > self.cache_limit {
            return false;
        }
        if self.cache_bytes.saturating_add(bytes) > self.cache_limit {
            self.lower.clear();
            self.leading.clear();
            self.coefficients.clear();
            self.cache_bytes = 0;
        }
        self.cache_bytes += bytes;
        true
    }

    fn product(&mut self, factors: &[Atom], endpoint: i64) -> Result<Coefficients> {
        let lower = factors
            .iter()
            .map(|factor| self.lower_bound(factor))
            .collect::<Result<Vec<_>>>()?;
        let total = lower.iter().sum::<i64>();
        let mut remaining = total;
        let mut result = BTreeMap::from([(0, Atom::one())]);
        for (factor, bound) in factors.iter().zip(lower) {
            let coefficients = self.through(factor, endpoint - total + bound)?;
            remaining -= bound;
            result = Self::convolve(&result, &coefficients, endpoint - remaining)?;
        }
        Ok(result)
    }

    fn inverse(&mut self, atom: &Atom, endpoint: i64) -> Result<Coefficients> {
        let (valuation, first) = self.leading(atom)?;
        let count = endpoint + valuation;
        if count < 0 {
            return Ok(BTreeMap::new());
        }
        ensure!(count <= MAX_ORDER, "unproven: inverse-jet depth budget");
        let base = self.through(atom, valuation + count)?;
        let mut inverse = vec![Atom::one() / &first];
        for degree in 1..=count {
            let mut sum = Atom::Zero;
            for offset in 1..=degree {
                if let Some(coefficient) = base.get(&(valuation + offset)) {
                    sum += coefficient * &inverse[(degree - offset) as usize];
                }
            }
            inverse.push(-sum / &first);
        }
        Ok(inverse
            .into_iter()
            .enumerate()
            .map(|(index, coefficient)| (-valuation + index as i64, coefficient))
            .collect())
    }

    fn energy(&mut self, invariant: &Atom, endpoint: i64) -> Result<Coefficients> {
        let (valuation, first) = self.leading(invariant)?;
        ensure!(
            valuation % 2 == 0,
            "unproven: half-integer energy valuation"
        );
        let energy_order = valuation / 2;
        let count = endpoint - energy_order;
        if count < 0 {
            return Ok(BTreeMap::new());
        }
        ensure!(count <= MAX_ORDER, "unproven: energy-jet depth budget");
        let base = self.through(invariant, valuation + count)?;
        let root = GS.on_shell_energy.call_args([Atom::Zero, first]);
        let root = self.algebra.normalize_energies(&root, self.lambda)?;
        let mut coefficients = vec![root.clone()];
        for degree in 1..=count {
            let mut convolution = Atom::Zero;
            for offset in 1..degree {
                convolution +=
                    &coefficients[offset as usize] * &coefficients[(degree - offset) as usize];
            }
            coefficients.push(
                (base
                    .get(&(valuation + degree))
                    .cloned()
                    .unwrap_or(Atom::Zero)
                    - convolution)
                    / (&root * 2),
            );
        }
        Ok(coefficients
            .into_iter()
            .enumerate()
            .map(|(index, coefficient)| (energy_order + index as i64, coefficient))
            .collect())
    }

    fn convolve(left: &Coefficients, right: &Coefficients, endpoint: i64) -> Result<Coefficients> {
        let mut result = BTreeMap::new();
        for (lp, lc) in left {
            for (rp, rc) in right {
                if lp + rp <= endpoint {
                    *result.entry(lp + rp).or_insert(Atom::Zero) += lc * rc;
                }
            }
        }
        Self::bounded(&result)?;
        Ok(result)
    }

    fn bounded(coefficients: &Coefficients) -> Result<()> {
        ensure!(
            coefficients.len() <= 2 * MAX_ORDER as usize + 1
                && coefficients
                    .values()
                    .map(|coefficient| coefficient.as_view().get_byte_size())
                    .sum::<usize>()
                    <= MAX_BYTES,
            "unproven: finite-jet coefficient budget"
        );
        Ok(())
    }

    fn leading_nonzero(&self, coefficient: &Atom) -> Result<bool> {
        if self.algebra.scalar_nonzero(coefficient)? {
            return Ok(true);
        }
        // A norm can vanish because another algebraic branch has a zero.
        // The independent positive-energy exact point can still certify this
        // coefficient nonzero. Unsupported points provide no evidence.
        for point in crate::uv::overlap_control::ExactControlPoint::points() {
            if point.nonzero_scalar_witness(coefficient).unwrap_or(false) {
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn power(atom: AtomView<'_>) -> Result<i64> {
        let exponent = i64::try_from(atom).map_err(|_| eyre!("unproven: noninteger jet power"))?;
        ensure!(
            exponent.unsigned_abs() <= MAX_POWER,
            "unproven: jet power budget"
        );
        Ok(exponent)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::{function, parse_lit, symbol};

    #[test]
    fn exact_inverse_retains_a_tiny_leading_coefficient() -> Result<()> {
        let t = symbol!("exact_soft_jet_test::t");
        let tiny = parse_lit!(1 / 10 ^ 30);
        let mut jets = ExactSoftJet::new(t, 1);
        let result = jets.coefficients(&(Atom::one() / (&tiny + Atom::var(t))))?;
        assert_eq!(result[&0], parse_lit!(10 ^ 30));
        assert_eq!(result[&1], -parse_lit!(10 ^ 60));
        Ok(())
    }

    #[test]
    fn exact_jets_lift_product_precision_and_do_not_drop_delayed_terms() -> Result<()> {
        let t = symbol!("exact_soft_jet_test::t");
        let x = Atom::var(t);
        let expression = (x.pow(9) + x.pow(11)) / x.pow(12);
        let coefficients = ExactSoftJet::new(t, -2).coefficients(&expression)?;
        assert_eq!(coefficients, BTreeMap::from([(-3, Atom::one())]));
        let inverse =
            ExactSoftJet::new(t, -7).coefficients(&(Atom::one() / (x.pow(9) + x.pow(11))))?;
        assert_eq!(
            inverse,
            BTreeMap::from([(-9, Atom::one()), (-7, -Atom::one())])
        );
        let algebraic_zero = parse_lit!((a + b) ^ 2 - a ^ 2 - 2 * a * b - b ^ 2);
        let delayed = Atom::one() / (algebraic_zero + x.pow(9));
        let delayed = ExactSoftJet::new(t, -9).coefficients(&delayed)?;
        assert_eq!(delayed, BTreeMap::from([(-9, Atom::one())]));
        let zero = ExactSoftJet::series(&x.pow(9), t, Rational::from(-2))?;
        assert!(zero.is_zero());
        assert_eq!(zero.absolute_order(), Rational::from(-1));
        Ok(())
    }

    #[test]
    fn exact_jet_cache_eviction_preserves_recomputed_coefficients() -> Result<()> {
        let lambda = symbol!("exact_soft_jet_cache_test::lambda");
        let t = Atom::var(lambda);
        let first = Atom::one() / (parse_lit!(a) + &t);
        let second = Atom::one() / (parse_lit!(b) + &t);
        let mut reference = ExactSoftJet::new(lambda, 2);
        let expected_first = reference.coefficients(&first)?;
        assert!(!reference.leading.is_empty());
        let first_cache_bytes = reference.cache_bytes;
        let expected_second = ExactSoftJet::new(lambda, 2).coefficients(&second)?;

        let mut bounded = ExactSoftJet::new(lambda, 2);
        bounded.cache_limit = first_cache_bytes;
        assert_eq!(bounded.coefficients(&first)?, expected_first);
        assert!(bounded.coefficients.contains_key(&(first.clone(), 2)));
        assert_eq!(bounded.coefficients(&second)?, expected_second);
        assert!(!bounded.coefficients.contains_key(&(first.clone(), 2)));
        assert!(bounded.cache_bytes <= bounded.cache_limit);
        assert_eq!(bounded.coefficients(&first)?, expected_first);
        assert!(bounded.cache_bytes <= bounded.cache_limit);
        Ok(())
    }

    #[test]
    fn exact_jet_entry_larger_than_cache_is_computed_without_caching() -> Result<()> {
        let lambda = symbol!("exact_soft_jet_cache_test::lambda");
        let t = Atom::var(lambda);
        let expression = Atom::one() / (parse_lit!(a) + t);
        let expected = ExactSoftJet::new(lambda, 2).coefficients(&expression)?;
        let mut uncached = ExactSoftJet::new(lambda, 2);
        uncached.cache_limit = 1;
        for _ in 0..2 {
            assert_eq!(uncached.coefficients(&expression)?, expected);
            assert!(uncached.lower.is_empty());
            assert!(uncached.leading.is_empty());
            assert!(uncached.coefficients.is_empty());
            assert_eq!(uncached.cache_bytes, 0);
        }
        Ok(())
    }

    #[test]
    fn exact_energy_jet_has_the_positive_soft_branch() -> Result<()> {
        crate::initialisation::test_initialise()?;
        let t = symbol!("exact_soft_jet_test::t");
        let x = Atom::var(t);
        let energy = function!(GS.on_shell_energy, 0, x.pow(2) * (1 + &x));
        let coefficients = ExactSoftJet::new(t, 3).coefficients(&energy)?;
        let e = function!(GS.on_shell_energy, 0, 1);
        assert_eq!(coefficients[&1], e.clone());
        assert!(matches!(
            SoftEnergyAlgebra::default().scalar_zero(&(coefficients[&2].clone() - &e / 2))?,
            EnergyZeroOutcome::Zero(_)
        ));
        assert!(matches!(
            SoftEnergyAlgebra::default().scalar_zero(&(coefficients[&3].clone() + &e / 8))?,
            EnergyZeroOutcome::Zero(_)
        ));
        Ok(())
    }
}
