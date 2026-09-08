//! Vacuum projection through universal tensor kernels and opaque coefficients.
//!
//! All source coefficients remain Symbolica ASTs. The only distributive
//! operations are convolution of loop degrees and contraction of tensor slots.
//! The FORM callback receives a coefficient-one universal tensor monomial.

use std::collections::{BTreeMap, BTreeSet};

use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol};
use symbolica::coefficient::CoefficientView;

use crate::utils::vakint_macros::vk_symbol;
use crate::{
    PySecDecOptions, Vakint, VakintError, VakintSettings, lorentz::LorentzTensor, symbols::S,
};

type Degrees = BTreeMap<i64, u32>;
type TensorTerms = BTreeMap<Degrees, Atom>;

/// A coefficient ring over scalar-product monomials. This same owning
/// representation is used for projected vacuum outputs and for scalar
/// nonvacuum inputs; the latter may have dot(k,p) in their monomial keys.
#[derive(Clone, Default)]
pub(crate) struct ScalarTerms(pub(crate) BTreeMap<Atom, Atom>);

impl ScalarTerms {
    /// Reuse the existing AtomField series for rational epsilon dependence.
    /// Generic epsilon-dependent functions need coefficient protection before
    /// Symbolica's function-series branch, which expands derivatives.
    pub(crate) fn epsilon_series(
        input: AtomView<'_>,
        settings: &VakintSettings,
        depth: symbolica::poly::series::SeriesDepth,
    ) -> Result<symbolica::poly::series::Series<symbolica::domains::atom::AtomField>, VakintError>
    {
        fn check(view: AtomView<'_>, epsilon: Symbol) -> Result<(), VakintError> {
            if !view.get_all_symbols(true).contains(&epsilon) {
                return Ok(());
            }
            match view {
                AtomView::Fun(_) => Err(VakintError::InvalidNumerator(
                    "generic epsilon-dependent functions require protected coefficient series"
                        .into(),
                )),
                AtomView::Add(sum) => sum.iter().try_for_each(|term| check(term, epsilon)),
                AtomView::Mul(product) => {
                    product.iter().try_for_each(|factor| check(factor, epsilon))
                }
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    check(base, epsilon)?;
                    check(exponent, epsilon)
                }
                _ => Ok(()),
            }
        }
        let epsilon = vk_symbol!(&settings.epsilon_symbol);
        check(input, epsilon)?;
        input
            .series(epsilon, 0, depth)
            .map_err(|error| VakintError::SymbolicaError(error.to_string()))
    }

    fn insert(&mut self, key: Atom, value: Atom) {
        if value.is_zero() {
            return;
        }
        let coefficient = self.0.entry(key.clone()).or_insert_with(Atom::zero);
        *coefficient += value;
        // Extracting loop monomials can leave a zero external coefficient in
        // a factored form, also inside a nonzero coefficient. Certify additive
        // subtrees before their parents, replacing only proven zeros. Keep all
        // collection on scratch candidates so surviving factors stay intact.
        // Collection can cycle between equivalent forms; never expand them.
        *coefficient = coefficient.replace_map_bottom_up(|term, _, out| {
            if !matches!(term, AtomView::Add(_)) {
                return;
            }
            let mut candidate = term.to_owned();
            let mut seen = BTreeSet::new();
            while !candidate.is_zero() && seen.insert(candidate.clone()) {
                // Numerical-content GCDs are safe only in the exact real
                // rational domain. Functions stay opaque, as in collect_num;
                // other powers can create complex or infinite coefficients.
                let mut rational = true;
                candidate.visitor(&mut |part| {
                    match part {
                        AtomView::Fun(_) => return false,
                        AtomView::Num(number) => {
                            let number = number.get_coeff_view();
                            rational &= matches!(
                                number,
                                CoefficientView::Natural(..) | CoefficientView::Large(..)
                            ) && number.is_real();
                        }
                        AtomView::Pow(power) => {
                            rational &= i64::try_from(power.get_base_exp().1)
                                .is_ok_and(|exponent| exponent > 0);
                        }
                        _ => {}
                    }
                    rational
                });
                candidate = if rational {
                    candidate.collect_num()
                } else {
                    candidate.collect_by_coefficient()
                };
            }
            if candidate.is_zero() {
                **out = Atom::zero();
            }
        });
        if coefficient.is_zero() {
            self.0.remove(&key);
        }
    }

    pub(crate) fn parse(view: AtomView<'_>) -> Result<Self, VakintError> {
        // Select only complete loop-dot atoms. Symbolica's existing AtomField
        // collector retains every external scalar/tensor AST as a coefficient;
        // it never converts the whole numerator into a polynomial over Q(i).
        fn inventory(view: AtomView<'_>, dots: &mut BTreeSet<Atom>) -> Result<i32, VakintError> {
            if !view.get_all_symbols(true).contains(&S.k) {
                return Ok(0);
            }
            let overflow = || {
                VakintError::InvalidNumerator("total loop scalar-product degree exceeds i32".into())
            };
            match view {
                AtomView::Fun(dot) if dot.get_symbol() == S.dot && dot.get_nargs() == 2 => {
                    for vector in dot {
                        let AtomView::Fun(vector) = vector else {
                            return Err(VakintError::InvalidNumerator(
                                "expected elementary vector".into(),
                            ));
                        };
                        if ![S.k, S.p].contains(&vector.get_symbol())
                            || vector.get_nargs() != 1
                            || i64::try_from(vector.get(0)).is_err()
                        {
                            return Err(VakintError::InvalidNumerator(
                                "expected elementary k(i) or p(i)".into(),
                            ));
                        }
                    }
                    dots.insert(view.to_owned());
                    Ok(1)
                }
                AtomView::Add(sum) => sum
                    .iter()
                    .try_fold(0i32, |degree, term| Ok(degree.max(inventory(term, dots)?))),
                AtomView::Mul(product) => product.iter().try_fold(0i32, |degree, factor| {
                    degree
                        .checked_add(inventory(factor, dots)?)
                        .ok_or_else(overflow)
                }),
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    let Ok(exponent @ 0..=2147483646) = i64::try_from(exponent) else {
                        return Err(VakintError::InvalidNumerator(
                            "loop dependence requires a nonnegative integer power below i32::MAX"
                                .into(),
                        ));
                    };
                    inventory(base, dots)?
                        .checked_mul(exponent as i32)
                        .ok_or_else(overflow)
                }
                _ => Err(VakintError::InvalidNumerator(format!(
                    "expected polynomial loop scalar products after contraction: {view}"
                ))),
            }
        }
        let mut dots = BTreeSet::new();
        inventory(view, &mut dots)?;
        let mut result = Self::default();
        if dots.is_empty() {
            result.insert(Atom::one(), view.to_owned());
            return Ok(result);
        }
        let dots = dots.into_iter().collect::<Vec<_>>();
        let mut used = view.get_all_symbols(true);
        let mut serial = 0;
        let mut aliases = BTreeMap::<Atom, Atom>::new();
        // AtomField retains products of sums, but its polynomial collector
        // can rewrite an unselected c^(a+b) as c^a*c^b. Hide maximal external
        // composite coefficients too, so their complete AST survives intact.
        let protected = view.replace_map(|term, _, out| {
            if term.get_all_symbols(true).contains(&S.k) {
                if matches!(term, AtomView::Fun(dot) if dot.get_symbol() == S.dot) {
                    // A selected complete dot is one variable: do not visit
                    // its momentum arguments or replace their integer IDs.
                    **out = term.to_owned();
                }
            } else if !matches!(term, AtomView::Num(_) | AtomView::Var(_)) {
                let alias = aliases.entry(term.to_owned()).or_insert_with(|| {
                    loop {
                        let symbol = vk_symbol!(format!("scalar_coefficient_{serial}"));
                        serial += 1;
                        if used.insert(symbol) {
                            break Atom::var(symbol);
                        }
                    }
                });
                **out = alias.clone();
            }
        });
        for (monomial, mut coefficient) in protected.coefficient_list::<i32>(&dots) {
            for (original, alias) in &aliases {
                coefficient = coefficient
                    .replace(alias.to_pattern())
                    .with(original.to_pattern());
            }
            if coefficient.get_all_symbols(true).contains(&S.k) {
                return Err(VakintError::InvalidNumerator(
                    "loop dependence remains outside a polynomial scalar-product monomial".into(),
                ));
            }
            result.insert(monomial, coefficient);
        }
        Ok(result)
    }

    pub(crate) fn expression(&self) -> Atom {
        self.0
            .iter()
            .fold(Atom::zero(), |sum, (monomial, coefficient)| {
                sum + coefficient * monomial.clone()
            })
    }

    /// Batch one synthetic scalar kernel, restoring each coefficient only
    /// after the backend returns. Callers must reserve topology symbols too.
    /// For general epsilon-dependent functions use protected series aliases;
    /// the direct series calls here cover rational epsilon dependence.
    #[allow(clippy::type_complexity)]
    pub(crate) fn backend_kernel(
        &self,
        settings: &VakintSettings,
        reserved: &[Atom],
    ) -> Result<(Atom, Vec<(Atom, Atom)>, i64), VakintError> {
        use symbolica::poly::series::SeriesDepth;
        let mut used = self.expression().get_all_symbols(true);
        for expression in reserved {
            used.extend(expression.get_all_symbols(true));
        }
        let mut kernel = Atom::zero();
        let mut aliases = Vec::new();
        let mut extra_orders = 0;
        let mut serial = 0;
        for (monomial, coefficient) in &self.0 {
            if coefficient.is_zero() {
                continue;
            }
            let series =
                Self::epsilon_series(coefficient.as_view(), settings, SeriesDepth::relative(1))?;
            let valuation = series.get_trailing_exponent();
            if !valuation.is_integer() {
                return Err(VakintError::InvalidNumerator(
                    "the numerator requires integer Laurent powers of epsilon".into(),
                ));
            }
            let valuation = valuation.numerator().to_i64().ok_or_else(|| {
                VakintError::InvalidNumerator("epsilon valuation exceeds i64".into())
            })?;
            let pole_order = valuation.checked_neg().ok_or_else(|| {
                VakintError::InvalidNumerator("epsilon pole order exceeds i64".into())
            })?;
            extra_orders = extra_orders.max(pole_order);
            let symbol = loop {
                let symbol = vk_symbol!(format!("vacuum_coefficient_{serial}"));
                serial += 1;
                if used.insert(symbol) {
                    break symbol;
                }
            };
            let alias = Atom::var(symbol);
            kernel += &alias * monomial.clone();
            aliases.push((alias, coefficient.clone()));
        }
        Ok((kernel, aliases, extra_orders))
    }

    /// Every coefficient is evaluated at the supplied point before one
    /// numerical integral is constructed; ERROR therefore stays correlated.
    /// The owner supplies the inclusive coefficient order: the requested
    /// answer's order plus a valid pole bound for its denominator class.
    pub(crate) fn numerical_kernel(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        options: &PySecDecOptions,
        coefficient_order: i64,
        reserved: &[Atom],
    ) -> Result<(Atom, PySecDecOptions, i64), VakintError> {
        use symbolica::{
            domains::float::{Complex, RealLike},
            poly::series::SeriesDepth,
        };

        let _ = vakint.scalar_numerator(settings, self.expression().as_view())?;
        let params = vakint.params_from_complex_f64(settings, &options.numerical_parameters);
        let externals = options
            .numerical_external_momenta
            .iter()
            .map(|(name, value)| {
                let id = name
                    .strip_prefix('p')
                    .and_then(|id| id.parse::<usize>().ok())
                    .ok_or_else(|| {
                        VakintError::EvaluationError(format!(
                            "external momentum name must be p followed by an integer: {name}"
                        ))
                    })?;
                Ok((id, *value))
            })
            .collect::<Result<std::collections::HashMap<_, _>, VakintError>>()?;
        let externals = vakint.externals_from_f64(settings, &externals);
        let mut values = Vec::new();
        let mut common_power = 0;
        for (monomial, coefficient) in &self.0 {
            let series = Self::epsilon_series(
                coefficient.as_view(),
                settings,
                SeriesDepth::absolute(coefficient_order),
            )?;
            for (power, coefficient) in series.terms() {
                if !power.is_integer() {
                    return Err(VakintError::InvalidNumerator(
                        "the numerator requires integer Laurent powers of epsilon".into(),
                    ));
                }
                let power = power.numerator().to_i64().ok_or_else(|| {
                    VakintError::InvalidNumerator("epsilon power exceeds i64".into())
                })?;
                // Evaluate each epsilon-independent coefficient. Passing the
                // full Laurent polynomial would unnecessarily narrow powers
                // through the numeric helper's i8 coefficient collector.
                let evaluated = Vakint::full_numerical_evaluation_without_error(
                    settings,
                    coefficient.as_view(),
                    &Default::default(),
                    &params,
                    Some(&externals),
                )?;
                let mut value = Complex::new(0.0, 0.0);
                for (coefficient_power, term) in evaluated.0 {
                    if coefficient_power != 0 {
                        return Err(VakintError::InvalidNumerator(
                            "contract Lorentz traces before constructing scalar coefficient sectors".into(),
                        ));
                    }
                    value.re += term.re.to_f64();
                    value.im += term.im.to_f64();
                }
                if !value.re.is_finite() || !value.im.is_finite() {
                    return Err(VakintError::EvaluationError(
                        "a numerical numerator coefficient is not finite".into(),
                    ));
                }
                if value.re == 0. && value.im == 0. {
                    continue;
                }
                common_power = common_power.min(power);
                values.push((monomial.clone(), power, value));
            }
        }
        let mut used = self.expression().get_all_symbols(true);
        for expression in reserved {
            used.extend(expression.get_all_symbols(true));
        }
        let mut options = options.clone();
        let mut kernel = Atom::zero();
        let epsilon = Atom::var(vk_symbol!(&settings.epsilon_symbol));
        let mut serial = 0;
        for (monomial, power, value) in values {
            let (alias, name) = loop {
                let symbol = vk_symbol!(format!("vacuum_numerical_coefficient_{serial}"));
                serial += 1;
                let alias = Atom::var(symbol);
                let name = crate::utils::undress_vakint_symbols(&alias.to_canonical_string());
                if used.insert(symbol) && !options.numerical_parameters.contains_key(&name) {
                    break (alias, name);
                }
            };
            options.numerical_parameters.insert(name, value);
            let shifted_power = power.checked_sub(common_power).ok_or_else(|| {
                VakintError::InvalidNumerator("epsilon power shift exceeds i64".into())
            })?;
            kernel += alias * epsilon.clone().pow(Atom::num(shifted_power)) * monomial.clone();
        }
        Ok((kernel, options, common_power))
    }
}

/// One instance per VakintExpression operation; its slot family and cache
/// never escape the owning operation. Include every expression's numerator
/// when choosing the family, so a cache hit cannot collide with user indices.
pub(crate) struct VacuumProjection {
    slot: Symbol,
    dimension: Atom,
    kernels: BTreeMap<Degrees, ScalarTerms>,
}

impl VacuumProjection {
    pub(crate) fn new(settings: &VakintSettings, inputs: &[Atom]) -> Self {
        let mut serial = 0;
        let slot = loop {
            let candidate = vk_symbol!(format!("vacuum_projection_slot_{serial}"));
            if !candidate.is_symmetric()
                && !candidate.is_antisymmetric()
                && !candidate.is_cyclesymmetric()
                && !candidate.is_linear()
                && candidate.get_normalization_function().is_none()
                && inputs
                    .iter()
                    .all(|input| !input.get_all_symbols(true).contains(&candidate))
            {
                break candidate;
            }
            serial += 1;
        };
        Self {
            slot,
            dimension: Atom::num(4)
                - Atom::num(2) * Atom::var(vk_symbol!(&settings.epsilon_symbol)),
            kernels: BTreeMap::new(),
        }
    }

    fn slot(&self, loop_id: i64, ordinal: u32) -> Atom {
        FunctionBuilder::new(S.tensor_index)
            .add_arg(Atom::var(self.slot))
            .add_arg(Atom::num(loop_id))
            .add_arg(Atom::num(i64::from(ordinal)))
            .finish()
    }

    fn metric(left: Atom, right: Atom) -> Atom {
        FunctionBuilder::new(S.g)
            .add_arg(left)
            .add_arg(right)
            .finish()
    }

    fn loop_id(vector: AtomView<'_>) -> Option<i64> {
        let AtomView::Fun(vector) = vector else {
            return None;
        };
        (vector.get_symbol() == S.k && vector.get_nargs() == 1)
            .then(|| i64::try_from(vector.get(0)).ok())
            .flatten()
    }

    fn angular_dependence(view: AtomView<'_>) -> bool {
        match view {
            AtomView::Fun(dot)
                if dot.get_symbol() == S.dot
                    && dot.get_nargs() == 2
                    && dot.iter().all(|vector| Self::loop_id(vector).is_some()) =>
            {
                false
            }
            AtomView::Fun(function) => {
                function.get_symbol() == S.k || function.iter().any(Self::angular_dependence)
            }
            AtomView::Add(sum) => sum.iter().any(Self::angular_dependence),
            AtomView::Mul(product) => product.iter().any(Self::angular_dependence),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                Self::angular_dependence(base) || Self::angular_dependence(exponent)
            }
            _ => false,
        }
    }

    fn product(&self, left: &TensorTerms, right: &TensorTerms) -> Result<TensorTerms, VakintError> {
        let mut result = TensorTerms::new();
        for (a, coefficient_a) in left {
            for (b, coefficient_b) in right {
                let mut degrees = a.clone();
                let mut shifted = coefficient_b.clone();
                for (loop_id, rank) in b {
                    let offset = *a.get(loop_id).unwrap_or(&0);
                    // Descending substitution is simultaneous for a positive
                    // shift: no newly-created target can be a later source.
                    for ordinal in (0..*rank).rev() {
                        shifted = shifted
                            .replace(self.slot(*loop_id, ordinal).to_pattern())
                            .with(
                                self.slot(
                                    *loop_id,
                                    ordinal.checked_add(offset).ok_or_else(|| {
                                        VakintError::InvalidNumerator(
                                            "tensor degree overflows u32".into(),
                                        )
                                    })?,
                                )
                                .to_pattern(),
                            );
                    }
                    degrees.insert(
                        *loop_id,
                        offset.checked_add(*rank).ok_or_else(|| {
                            VakintError::InvalidNumerator("tensor degree overflows u32".into())
                        })?,
                    );
                }
                *result.entry(degrees).or_insert_with(Atom::zero) += coefficient_a * shifted;
            }
        }
        result.retain(|_, coefficient| !coefficient.is_zero());
        Ok(result)
    }

    fn tensor_terms(&self, view: AtomView<'_>) -> Result<TensorTerms, VakintError> {
        if !Self::angular_dependence(view) {
            return Ok(BTreeMap::from([(Degrees::new(), view.to_owned())]));
        }
        match view {
            AtomView::Add(sum) => {
                let mut result = TensorTerms::new();
                for term in sum {
                    for (degree, coefficient) in self.tensor_terms(term)? {
                        *result.entry(degree).or_insert_with(Atom::zero) += coefficient;
                    }
                }
                result.retain(|_, coefficient| !coefficient.is_zero());
                Ok(result)
            }
            AtomView::Mul(product) => {
                let mut result = BTreeMap::from([(Degrees::new(), Atom::one())]);
                for factor in product {
                    result = self.product(&result, &self.tensor_terms(factor)?)?;
                }
                Ok(result)
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let exponent = u32::try_from(i64::try_from(exponent).map_err(|_| {
                    VakintError::InvalidNumerator("nonpolynomial loop tensor".into())
                })?)
                .map_err(|_| VakintError::InvalidNumerator("nonpolynomial loop tensor".into()))?;
                let mut base = self.tensor_terms(base)?;
                let mut exponent = exponent;
                let mut result = BTreeMap::from([(Degrees::new(), Atom::one())]);
                while exponent != 0 {
                    if exponent % 2 == 1 {
                        result = self.product(&result, &base)?;
                    }
                    exponent /= 2;
                    if exponent != 0 {
                        base = self.product(&base, &base)?;
                    }
                }
                Ok(result)
            }
            AtomView::Fun(vector) if vector.get_symbol() == S.k && vector.get_nargs() == 2 => {
                let loop_id = i64::try_from(vector.get(0)).map_err(|_| {
                    VakintError::InvalidNumerator("loop momentum id must be an integer".into())
                })?;
                Ok(BTreeMap::from([(
                    BTreeMap::from([(loop_id, 1)]),
                    Self::metric(vector.get(1).to_owned(), self.slot(loop_id, 0)),
                )]))
            }
            AtomView::Fun(dot) if dot.get_symbol() == S.dot && dot.get_nargs() == 2 => {
                let (loop_id, external) = if let Some(loop_id) = Self::loop_id(dot.get(0)) {
                    (loop_id, dot.get(1))
                } else if let Some(loop_id) = Self::loop_id(dot.get(1)) {
                    (loop_id, dot.get(0))
                } else {
                    return Err(VakintError::InvalidNumerator(
                        "unsupported loop scalar product".into(),
                    ));
                };
                let AtomView::Fun(external) = external else {
                    return Err(VakintError::InvalidNumerator(
                        "expected elementary external vector".into(),
                    ));
                };
                if external.get_symbol() != S.p
                    || external.get_nargs() != 1
                    || i64::try_from(external.get(0)).is_err()
                {
                    return Err(VakintError::InvalidNumerator(
                        "expected elementary external vector".into(),
                    ));
                }
                Ok(BTreeMap::from([(
                    BTreeMap::from([(loop_id, 1)]),
                    FunctionBuilder::new(S.p)
                        .add_arg(external.get(0))
                        .add_arg(self.slot(loop_id, 0))
                        .finish(),
                )]))
            }
            _ => Err(VakintError::InvalidNumerator(format!(
                "unsupported polynomial loop-tensor factor: {view}"
            ))),
        }
    }

    /// The owner must certify vacuum denominators before calling this
    /// angular projector. Nonvacuum scalar inputs use ScalarTerms directly.
    pub(crate) fn project(
        &mut self,
        input: AtomView<'_>,
        mut universal_projector: impl FnMut(Atom) -> Result<Atom, VakintError>,
    ) -> Result<ScalarTerms, VakintError> {
        // This is a contraction, never a call to FORM on the input. It also
        // normalizes closed tensor powers and checks literal dummy incidence.
        let input = LorentzTensor::parse(input, &self.dimension)?.expression;
        let mut result = ScalarTerms::default();
        for (degrees, coefficient) in self.tensor_terms(input.as_view())? {
            let odd = degrees
                .values()
                .fold(false, |odd, rank| odd ^ (rank % 2 == 1));
            if odd {
                continue;
            }
            if !self.kernels.contains_key(&degrees) {
                let mut monomial = Atom::one();
                for (loop_id, rank) in &degrees {
                    for ordinal in 0..*rank {
                        monomial *= FunctionBuilder::new(S.k)
                            .add_arg(Atom::num(*loop_id))
                            .add_arg(self.slot(*loop_id, ordinal))
                            .finish();
                    }
                }
                let kernel = if degrees.is_empty() {
                    Atom::one()
                } else {
                    universal_projector(monomial)?
                };
                self.kernels
                    .insert(degrees.clone(), ScalarTerms::parse(kernel.as_view())?);
            }
            for (kernel_monomial, metric_coefficient) in &self.kernels[&degrees].0 {
                let coefficient = LorentzTensor::parse(
                    (&coefficient * metric_coefficient).as_view(),
                    &self.dimension,
                )?
                .expression;
                // Scalar invariants retained above now enter only their
                // monomial key. This also deduplicates before backend calls.
                for (coefficient_monomial, coefficient) in
                    ScalarTerms::parse(coefficient.as_view())?.0
                {
                    result.insert(kernel_monomial * coefficient_monomial, coefficient);
                }
            }
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::vakint_macros::vk_parse;

    fn atom(source: &str) -> Atom {
        let _ = &*S;
        vk_parse!(source).unwrap()
    }

    fn rank_two_kernel(left: Atom, right: Atom, loops: &str, dimension: &Atom) -> Atom {
        VacuumProjection::metric(left, right) * atom(loops) / dimension
    }

    #[test]
    fn projected_zero_external_sector_preserves_the_nonzero_spectator() {
        let settings = VakintSettings::default();
        let spectator = atom("(a+2*b)*(c+3*d)");
        let input = &spectator
            * atom(
                "-2*(-dot(p(3),p(3))-dot(p(3),p(4))+dot(p(3),k(1)))
                 -2*(-dot(p(3),p(4))-dot(p(4),p(4))+dot(p(4),k(1)))
                 -(dot(p(3),p(3))+2*dot(p(3),p(4))-2*dot(p(3),k(1))
                   +dot(p(4),p(4))-2*dot(p(4),k(1))+dot(k(1),k(1)))
                 -dot(p(3),p(3))-2*dot(p(3),p(4))-dot(p(4),p(4))",
            );
        let result = VacuumProjection::new(&settings, std::slice::from_ref(&input))
            .project(input.as_view(), |_| {
                panic!("odd angular terms vanish and scalar loop products need no kernel")
            })
            .unwrap();
        assert_eq!(result.0.len(), 1);
        assert_eq!(result.expression(), -spectator * atom("dot(k(1),k(1))"));
    }

    #[test]
    fn scalar_zero_certificate_handles_cycles_and_accumulated_cancellation() {
        let complex: Atom = Atom::num(2) + Atom::i() * 2;
        for spectator in [
            atom("(a+2*b)*(c+3*d)"),
            &complex * atom("a") + atom("3*b"),
            atom("(1.25`80*a+2.5`80*b)*(c+d)"),
            atom("(2*(a+b))^(1/2)+3*c"),
            atom("(2*(a+b))^-1+3*c"),
        ] {
            assert_eq!(
                ScalarTerms::parse(spectator.as_view())
                    .unwrap()
                    .expression(),
                spectator
            );
        }
        let mut terms = ScalarTerms::default();
        terms.insert(Atom::one(), atom("2*(a+b)"));
        terms.insert(Atom::one(), atom("-2*a-2*b"));
        assert!(terms.0.is_empty());
        terms.insert(Atom::one(), &complex * atom("a+b"));
        terms.insert(Atom::one(), -&complex * atom("a") - &complex * atom("b"));
        assert!(terms.0.is_empty());
    }

    #[test]
    fn scalar_zero_subtrees_disappear_inside_complex_coefficients() {
        let zero = atom("2*(-a-b)+2*(-b-c)+2*a+4*b+2*c");
        let spectator = atom("(u+2*v)*(w+3*x)");
        let complex = Atom::i() * atom("eps*(2+log(mu)^2)+log(mu)/pi^2");
        let input = &spectator + &complex * &zero;
        assert_eq!(
            ScalarTerms::parse(input.as_view()).unwrap().expression(),
            spectator
        );

        let settings = VakintSettings::default();
        let positive = atom(
            "2*(-dot(p(3),p(3))-dot(p(3),p(4))+dot(p(3),k(1)))
             +2*(-dot(p(3),p(4))-dot(p(4),p(4))+dot(p(4),k(1)))
             +2*dot(p(3),p(3))+4*dot(p(3),p(4))-2*dot(p(3),k(1))
             +2*dot(p(4),p(4))-2*dot(p(4),k(1))+dot(k(1),k(1))",
        );
        let negative = atom(
            "-2*(-dot(p(3),p(3))-dot(p(3),p(4))+dot(p(3),k(1)))
             -2*(-dot(p(3),p(4))-dot(p(4),p(4))+dot(p(4),k(1)))
             -(dot(p(3),p(3))+2*dot(p(3),p(4))-2*dot(p(3),k(1))
               +dot(p(4),p(4))-2*dot(p(4),k(1))+dot(k(1),k(1)))
             -dot(p(3),p(3))-2*dot(p(3),p(4))-dot(p(4),p(4))",
        );
        let input = &spectator
            * (atom("eps") * (&complex * &positive + &positive) + atom("log(mu)") * &positive)
            * negative;
        let result = VacuumProjection::new(&settings, std::slice::from_ref(&input))
            .project(input.as_view(), |_| {
                panic!("canceled angular terms require no tensor kernel")
            })
            .unwrap();
        assert_eq!(result.0.len(), 1);
        assert_eq!(
            result.expression(),
            -spectator
                * (atom("eps") * (complex + Atom::one()) + atom("log(mu)"))
                * atom("dot(k(1),k(1))^2")
        );
    }

    #[test]
    fn same_rank_sums_merge_before_products_and_spectators_stay_factorized() {
        let settings = VakintSettings::default();
        let spectator = atom("(a+b)*(c+d)");
        let input =
            &spectator * atom("(dot(k(1),p(1))+dot(k(1),p(2)))*(dot(k(1),p(3))+dot(k(1),p(4)))");
        let mut projection = VacuumProjection::new(&settings, std::slice::from_ref(&input));
        let kernel = rank_two_kernel(
            projection.slot(1, 0),
            projection.slot(1, 1),
            "dot(k(1),k(1))",
            &projection.dimension,
        );
        let expected_monomial = FunctionBuilder::new(S.k)
            .add_arg(Atom::num(1))
            .add_arg(projection.slot(1, 0))
            .finish()
            * FunctionBuilder::new(S.k)
                .add_arg(Atom::num(1))
                .add_arg(projection.slot(1, 1))
                .finish();
        let mut calls = 0;
        let result = projection
            .project(input.as_view(), |monomial| {
                calls += 1;
                assert_eq!(monomial, expected_monomial);
                Ok(kernel.clone())
            })
            .unwrap();
        assert_eq!(calls, 1);
        assert_eq!(result.0.len(), 1);
        assert_eq!(
            result.0.values().next().unwrap(),
            &(&spectator * atom("dot(p(1),p(3))+dot(p(1),p(4))+dot(p(2),p(3))+dot(p(2),p(4))",)
                / &projection.dimension),
        );
        // The comparison above is exact structural Atom equality; there is
        // no expanded source numerator or baseline projection in this test.
    }

    #[test]
    fn two_loop_rank_two_has_dot12_over_dimension_and_keeps_open_indices() {
        let settings = VakintSettings::default();
        let input = atom("(a+b)*(c+d)*k(1,mu)*k(2,nu)");
        let mut projection = VacuumProjection::new(&settings, std::slice::from_ref(&input));
        let kernel = rank_two_kernel(
            projection.slot(1, 0),
            projection.slot(2, 0),
            "dot(k(1),k(2))",
            &projection.dimension,
        );
        let result = projection
            .project(input.as_view(), |_| Ok(kernel.clone()))
            .unwrap();
        assert_eq!(
            result.expression(),
            atom("(a+b)*(c+d)*g(mu,nu)*dot(k(1),k(2))") / &projection.dimension
        );
    }

    #[test]
    fn declared_spin_tensor_slots_survive_the_universal_kernel_boundary() {
        let settings = VakintSettings::default();
        let input = atom(
            "(a+b)*(c+d)*(tensor(gamma(mu),mu)*k(1,mu)+tensor(gamma(nu),nu)*k(1,nu))*k(1,rho)",
        );
        let mut projection = VacuumProjection::new(&settings, std::slice::from_ref(&input));
        let kernel = rank_two_kernel(
            projection.slot(1, 0),
            projection.slot(1, 1),
            "dot(k(1),k(1))",
            &projection.dimension,
        );
        let result = projection
            .project(input.as_view(), |monomial| {
                assert!(!monomial.get_all_symbols(true).contains(&S.tensor));
                Ok(kernel.clone())
            })
            .unwrap();
        assert_eq!(
            result.expression(),
            atom("2*(a+b)*(c+d)*tensor(gamma(rho),rho)*dot(k(1),k(1))") / &projection.dimension
        );
    }

    #[test]
    fn closed_spin_tensor_projection_preserves_coefficients_up_to_dummy_renaming() {
        let settings = VakintSettings::default();
        let spectator = atom("(a+b)*(c+d)");
        let input = &spectator * atom("tensor(gamma(mu),mu)*k(1,mu)*tensor(gamma(nu),nu)*k(1,nu)");
        let mut projection = VacuumProjection::new(&settings, std::slice::from_ref(&input));
        let kernel = rank_two_kernel(
            projection.slot(1, 0),
            projection.slot(1, 1),
            "dot(k(1),k(1))",
            &projection.dimension,
        );
        let result = projection
            .project(input.as_view(), |_| Ok(kernel.clone()))
            .unwrap();
        let coefficient = &result.0[&atom("dot(k(1),k(1))")];
        assert!(
            LorentzTensor::parse(coefficient.as_view(), &projection.dimension)
                .unwrap()
                .is_scalar()
        );
        // A closed contraction may retain either original dummy or a projector
        // slot. Compare the complete tensor coefficient after alpha-renaming.
        let tensor = atom("tensor(gamma(index_),index_)").to_pattern();
        let indices = coefficient
            .pattern_match(&tensor, None, None)
            .map(|matched| matched[&vk_symbol!("index_")].clone())
            .collect::<BTreeSet<_>>();
        assert_eq!(indices.len(), 1);
        let restored = coefficient
            .replace(indices.iter().next().unwrap().to_pattern())
            .with(atom("mink(4,fresh_dummy)"));
        assert_eq!(
            restored,
            &spectator * atom("tensor(gamma(mink(4,fresh_dummy)),mink(4,fresh_dummy))^2")
                / &projection.dimension
        );

        // Force the marker boundary separately: a metric with one fresh slot
        // must rename the opaque leaf and its declared slot consistently.
        let marker = projection.slot(1, 0);
        let input = &spectator
            * atom("tensor(gamma(mu),mu)")
            * VacuumProjection::metric(atom("mu"), marker.clone());
        let marked = LorentzTensor::parse(input.as_view(), &projection.dimension)
            .unwrap()
            .expression;
        assert!(marked.get_all_symbols(true).contains(&S.tensor_index));
        let restored = marked
            .replace(marker.to_pattern())
            .with(atom("mink(4,fresh_dummy)"));
        assert_eq!(
            restored,
            spectator * atom("tensor(gamma(mink(4,fresh_dummy)),mink(4,fresh_dummy))")
        );
    }

    #[test]
    fn odd_total_rank_vanishes_without_calling_form() {
        let settings = VakintSettings::default();
        let input = atom("(a+b)*(c+d)*dot(k(1),p(1))");
        let mut projection = VacuumProjection::new(&settings, std::slice::from_ref(&input));
        let result = projection
            .project(input.as_view(), |_| panic!("odd kernel needs no FORM"))
            .unwrap();
        assert!(result.expression().is_zero());
    }

    #[test]
    fn invariant_powers_do_not_increase_projector_rank() {
        let settings = VakintSettings::default();
        let input = atom("dot(k(1),k(1))^20*dot(k(1),p(1))*dot(k(1),p(2))");
        let mut projection = VacuumProjection::new(&settings, std::slice::from_ref(&input));
        let kernel = rank_two_kernel(
            projection.slot(1, 0),
            projection.slot(1, 1),
            "dot(k(1),k(1))",
            &projection.dimension,
        );
        let result = projection
            .project(input.as_view(), |monomial| {
                let AtomView::Mul(product) = monomial.as_view() else {
                    panic!("rank two")
                };
                assert_eq!(product.get_nargs(), 2);
                Ok(kernel.clone())
            })
            .unwrap();
        assert_eq!(
            result.expression(),
            atom("dot(k(1),k(1))^21*dot(p(1),p(2))") / &projection.dimension
        );
    }

    #[test]
    fn scalar_sectors_keep_external_products_and_deduplicate_loop_monomials() {
        let input = atom("(a+b)*(c+d)*(dot(k(1),k(2))+dot(k(1),k(1)))");
        let terms = ScalarTerms::parse(input.as_view()).unwrap();
        assert_eq!(terms.0.len(), 2);
        assert!(
            terms
                .0
                .values()
                .all(|coefficient| *coefficient == atom("(a+b)*(c+d)"))
        );
        assert!(ScalarTerms::parse(atom("dot(k(1),p(1))^-1").as_view()).is_err());
        let external = atom("c^(a+b)*f((a+b)*(c+d))");
        let input = &external * atom("dot(k(1),p(1))+dot(k(1),k(2))");
        let terms = ScalarTerms::parse(input.as_view()).unwrap();
        assert!(terms.0.values().all(|coefficient| *coefficient == external));
    }

    #[test]
    fn projection_is_idempotent_for_internal_scalar_products() {
        let settings = VakintSettings::default();
        let input = atom("(a+b)*(c+d)*dot(k(1),k(2))*dot(p(1),p(2))");
        let mut projection = VacuumProjection::new(&settings, std::slice::from_ref(&input));
        let result = projection
            .project(input.as_view(), |_| panic!("already projected"))
            .unwrap();
        assert_eq!(result.expression(), input);
    }

    #[test]
    fn rational_epsilon_series_uses_existing_atom_coefficient_ring() {
        use symbolica::poly::series::SeriesDepth;
        let epsilon = vk_symbol!("eps");
        let source = atom("(a+b)*(c+d)/(eps*(4-2*eps))");
        let leading = source.series(epsilon, 0, SeriesDepth::relative(1)).unwrap();
        assert_eq!(leading.get_trailing_exponent(), (-1, 1));
        let result = source.series(epsilon, 0, SeriesDepth::absolute(0)).unwrap();
        // Compare the series coefficients individually: scalar ASTs remain
        // products of sums, while epsilon alone is series-expanded.
        let terms = result.terms().collect::<Vec<_>>();
        assert_eq!(terms.len(), 2);
        assert_eq!(*terms[0].1, atom("(a+b)*(c+d)/4"));
        assert_eq!(*terms[1].1, atom("(a+b)*(c+d)/8"));
    }

    #[test]
    fn rank_four_kernel_has_the_fixed_minkowski_0011_component() {
        let settings = VakintSettings::default();
        let input = atom("k(1,mu)*k(1,nu)*k(1,rho)*k(1,sigma)");
        let mut projection = VacuumProjection::new(&settings, std::slice::from_ref(&input));
        projection.dimension = Atom::num(4);
        let slots = (0..4)
            .map(|ordinal| projection.slot(1, ordinal))
            .collect::<Vec<_>>();
        let metric =
            |a: usize, b: usize| VacuumProjection::metric(slots[a].clone(), slots[b].clone());
        let kernel = atom("dot(k(1),k(1))^2")
            * (metric(0, 1) * metric(2, 3)
                + metric(0, 2) * metric(1, 3)
                + metric(0, 3) * metric(1, 2))
            / Atom::num(24);
        let result = projection
            .project(input.as_view(), |_| Ok(kernel.clone()))
            .unwrap();
        let mut coefficient = result.0[&atom("dot(k(1),k(1))^2")].clone();
        // Fixed Minkowski component mu=nu=0, rho=sigma=1: only g00*g11 survives.
        for (metric, value) in [
            ("g(mu,nu)", 1),
            ("g(rho,sigma)", -1),
            ("g(mu,rho)", 0),
            ("g(nu,sigma)", 0),
            ("g(mu,sigma)", 0),
            ("g(nu,rho)", 0),
        ] {
            coefficient = coefficient
                .replace(atom(metric).to_pattern())
                .with(Atom::num(value));
        }
        assert_eq!(coefficient, atom("-1/24"));
    }

    #[test]
    fn mixed_rank_two_two_has_both_fixed_minkowski_components() {
        let settings = VakintSettings::default();
        let input = atom("k(1,mu)*k(1,nu)*k(2,rho)*k(2,sigma)");
        let mut projection = VacuumProjection::new(&settings, std::slice::from_ref(&input));
        projection.dimension = Atom::num(4);
        let slots = [
            projection.slot(1, 0),
            projection.slot(1, 1),
            projection.slot(2, 0),
            projection.slot(2, 1),
        ];
        let metric =
            |a: usize, b: usize| VacuumProjection::metric(slots[a].clone(), slots[b].clone());
        // The two metric traces give A=(5*q11*q22-2*q12^2)/72 and
        // B=(4*q12^2-q11*q22)/72 in four dimensions.
        let kernel = atom("(5*dot(k(1),k(1))*dot(k(2),k(2))-2*dot(k(1),k(2))^2)/72")
            * metric(0, 1)
            * metric(2, 3)
            + atom("(4*dot(k(1),k(2))^2-dot(k(1),k(1))*dot(k(2),k(2)))/72")
                * (metric(0, 2) * metric(1, 3) + metric(0, 3) * metric(1, 2));
        let result = projection
            .project(input.as_view(), |_| Ok(kernel.clone()))
            .unwrap();
        assert_eq!(result.0.len(), 2);
        let metrics = [
            "g(mu,nu)",
            "g(rho,sigma)",
            "g(mu,rho)",
            "g(nu,sigma)",
            "g(mu,sigma)",
            "g(nu,rho)",
        ];
        for (values, expected) in [
            ([1, -1, 0, 0, 0, 0], ["-5/72", "1/36"]),
            ([0, 0, 1, -1, 0, 0], ["1/72", "-1/18"]),
        ] {
            for (monomial, expected) in [
                ("dot(k(1),k(1))*dot(k(2),k(2))", expected[0]),
                ("dot(k(1),k(2))^2", expected[1]),
            ] {
                let mut coefficient = result.0[&atom(monomial)].clone();
                for (metric, value) in metrics.iter().zip(values) {
                    coefficient = coefficient
                        .replace(atom(metric).to_pattern())
                        .with(Atom::num(value));
                }
                assert_eq!(coefficient, atom(expected));
            }
        }
    }

    #[test]
    fn validation_precedes_odd_projection_and_bounds_collector_degrees() {
        let settings = VakintSettings::default();
        let input = atom("1+k(1,mu)");
        let mut projection = VacuumProjection::new(&settings, std::slice::from_ref(&input));
        assert!(
            projection
                .project(input.as_view(), |_| panic!("invalid tensor"))
                .is_err()
        );
        for input in [
            "dot(k(1),p(1))^-1",
            "dot(k(1),p(1))^(1/2)",
            "f(dot(k(1),p(1)))",
            "(dot(k(1),k(1))+dot(k(1),k(2)))^2147483646*dot(k(1),k(1))^2",
        ] {
            assert!(
                ScalarTerms::parse(atom(input).as_view()).is_err(),
                "{input}"
            );
        }
    }

    #[test]
    fn analytic_aliases_preserve_coefficients_and_account_for_epsilon_poles() {
        let settings = VakintSettings {
            epsilon_symbol: "eps".into(),
            ..VakintSettings::default()
        };
        let input = atom("(a+b)*(c+d)*dot(k(1),k(1))/eps^2+(a+b)*dot(k(1),k(2))/eps");
        let terms = ScalarTerms::parse(input.as_view()).unwrap();
        let (kernel, aliases, extra_orders) = terms.backend_kernel(&settings, &[]).unwrap();
        assert_eq!(extra_orders, 2);
        assert_eq!(aliases.len(), 2);
        for name in ["a", "b", "c", "d", "eps"] {
            assert!(!kernel.get_all_symbols(true).contains(&vk_symbol!(name)));
        }
        let restored = aliases.iter().fold(kernel, |kernel, (alias, coefficient)| {
            kernel
                .replace(alias.to_pattern())
                .with(coefficient.to_pattern())
        });
        assert_eq!(restored, input);
        assert!(
            ScalarTerms::epsilon_series(
                atom("f(eps*(a+b)*(c+d))").as_view(),
                &settings,
                symbolica::poly::series::SeriesDepth::absolute(1),
            )
            .is_err()
        );
    }
}
