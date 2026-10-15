use std::{ops::Deref, sync::Arc};

use color_eyre::eyre::{bail, ensure};
use spenso::{
    network::parsing::{
        AtomStructureExt, ParseSettings, SchoonschipExpansionMode, ShorthandParsing,
        StrictTensorFilter,
    },
    structure::representation::{Minkowski, RepName},
};

use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomType, AtomView, Symbol},
    coefficient::{Coefficient, CoefficientView},
    domains::atom::AtomField,
    function,
    poly::{polynomial::MultivariatePolynomial, series::SeriesDepth},
    symbol,
};

use crate::{
    cff::expression::OrientationID,
    utils::{GS, TENSORLIB, W_},
};

use super::ParsingNet;
pub type ParsingNetError = spenso::network::TensorNetworkError<
    spenso::structure::IndexlessNamedStructure<
        Symbol,
        Vec<Atom>,
        spenso::structure::representation::LibraryRep,
        super::aind::Aind,
    >,
    Symbol,
>;

pub trait NumeratorAtomExt {
    /// Collect a bounded Laurent polynomial in literal, independent keys using
    /// exact zero tests. Complete scalar sums are admitted as indivisible keys;
    /// their bodies and key-independent subtrees remain opaque coefficients.
    /// Return None for unsupported dependence or an exceeded algebra budget.
    fn coefficient_list_exact(&self, keys: &[Atom]) -> Option<Vec<(Atom, Atom)>>;

    /// Recover Horner forms and common factors without distributing graph
    /// numerators or extracting denominators across branch guards.
    fn collect_compact_factors(&self) -> Atom;

    /// Cancel removable scalar denominators without distributing numerator
    /// factors. Collect a bounded polynomial of opaque numerator calls, without
    /// distributing products or powers of their sums. Guards and function
    /// arguments retain their boundaries.
    /// Rational normalization requires exact rational coefficients, has a
    /// preflight algebra budget and must reduce size.
    fn cancel_scalar_poles(&self, numerator_keys: &[Atom]) -> Atom;

    /// Truncate through an absolute integer order, preserving an existing root
    /// product's variable-independent factors outside the coefficient sum.
    /// Other roots retain native series behavior, including additive zeros.
    /// Regroup the result in the caller's owned, Taylor-independent numerator
    /// coefficient keys. These must be symbols or functions occurring linearly,
    /// outside other functions and inverses. Their definitions and all dependent
    /// energy arguments remain the caller's responsibility; this method never
    /// makes a Taylor-dependent numerator opaque.
    fn series_preserving_factors(
        &self,
        variable: Symbol,
        expansion_point: AtomView<'_>,
        depth: i64,
        numerator_family_keys: &[Atom],
    ) -> color_eyre::Result<Atom>;

    fn map_mink_dim<'a>(&self, dim: impl Into<AtomOrView<'a>>) -> Atom;

    fn unwrap_function(&self, symbol: Symbol) -> Atom;

    #[allow(clippy::result_large_err)]
    fn parse_into_net(&self) -> Result<ParsingNet, ParsingNetError>;
}

impl NumeratorAtomExt for Atom {
    fn coefficient_list_exact(&self, keys: &[Atom]) -> Option<Vec<(Atom, Atom)>> {
        self.as_view().coefficient_list_exact(keys)
    }

    fn cancel_scalar_poles(&self, numerator_keys: &[Atom]) -> Atom {
        self.as_view().cancel_scalar_poles(numerator_keys)
    }

    fn collect_compact_factors(&self) -> Atom {
        self.as_view().collect_compact_factors()
    }

    fn series_preserving_factors(
        &self,
        variable: Symbol,
        expansion_point: AtomView<'_>,
        depth: i64,
        numerator_family_keys: &[Atom],
    ) -> color_eyre::Result<Atom> {
        self.as_view().series_preserving_factors(
            variable,
            expansion_point,
            depth,
            numerator_family_keys,
        )
    }

    fn map_mink_dim<'a>(&self, dim: impl Into<AtomOrView<'a>>) -> Atom {
        self.as_view().map_mink_dim(dim)
    }

    fn unwrap_function(&self, symbol: Symbol) -> Atom {
        self.as_view().unwrap_function(symbol)
    }

    fn parse_into_net(&self) -> Result<ParsingNet, ParsingNetError> {
        self.as_view().parse_into_net()
    }
}

/// Upper bounds before rational polynomial conversion, including the expanded
/// denominator of a sum. Saturating arithmetic makes compact high powers cheap
/// to reject; a post-conversion byte limit would be too late.
struct ScalarRationalSize {
    numerator_terms: usize,
    denominator_terms: usize,
    numerator_degree: usize,
    denominator_degree: usize,
}

impl ScalarRationalSize {
    fn estimate(atom: AtomView<'_>) -> Option<Self> {
        let leaf = Self {
            numerator_terms: 1,
            denominator_terms: 1,
            numerator_degree: usize::from(!matches!(atom, AtomView::Num(_))),
            denominator_degree: 0,
        };
        let size = match atom {
            AtomView::Fun(fun) => {
                if fun.get_symbol() != GS.on_shell_energy
                    && atom.is_tensorial(StrictTensorFilter::ContainsReps)
                {
                    return None;
                }
                leaf
            }
            AtomView::Var(_) if atom.is_tensorial(StrictTensorFilter::ContainsReps) => return None,
            AtomView::Num(number) => {
                // Integrated child coefficients can contain floating zeta
                // values. Rational cancellation cannot represent that domain;
                // leave it untouched, including inside sums and inverses.
                if !matches!(
                    number.get_coeff_view(),
                    CoefficientView::Natural(..) | CoefficientView::Large(..)
                ) {
                    return None;
                }
                leaf
            }
            AtomView::Var(_) => leaf,
            AtomView::Pow(power) => {
                let exponent = i64::try_from(power.get_exp()).ok()?;
                if exponent.unsigned_abs() > 64 {
                    return None;
                }
                let mut base = Self::estimate(power.get_base())?;
                if exponent < 0 {
                    std::mem::swap(&mut base.numerator_terms, &mut base.denominator_terms);
                    std::mem::swap(&mut base.numerator_degree, &mut base.denominator_degree);
                }
                let exponent = exponent.unsigned_abs() as u32;
                Self {
                    numerator_terms: base.numerator_terms.saturating_pow(exponent),
                    denominator_terms: base.denominator_terms.saturating_pow(exponent),
                    numerator_degree: base.numerator_degree.saturating_mul(exponent as usize),
                    denominator_degree: base.denominator_degree.saturating_mul(exponent as usize),
                }
            }
            AtomView::Add(_) | AtomView::Mul(_) => {
                let (add, mut arguments) = match atom {
                    AtomView::Add(sum) => (true, sum.iter()),
                    AtomView::Mul(product) => (false, product.iter()),
                    _ => unreachable!(),
                };
                let first = Self::estimate(arguments.next()?)?;
                arguments.try_fold(first, |left, right| {
                    let right = Self::estimate(right)?;
                    Self {
                        numerator_terms: if add {
                            left.numerator_terms
                                .saturating_mul(right.denominator_terms)
                                .saturating_add(
                                    right.numerator_terms.saturating_mul(left.denominator_terms),
                                )
                        } else {
                            left.numerator_terms.saturating_mul(right.numerator_terms)
                        },
                        denominator_terms: left
                            .denominator_terms
                            .saturating_mul(right.denominator_terms),
                        numerator_degree: if add {
                            left.numerator_degree
                                .saturating_add(right.denominator_degree)
                                .max(
                                    right
                                        .numerator_degree
                                        .saturating_add(left.denominator_degree),
                                )
                        } else {
                            left.numerator_degree.saturating_add(right.numerator_degree)
                        },
                        denominator_degree: left
                            .denominator_degree
                            .saturating_add(right.denominator_degree),
                    }
                    .bounded()
                })?
            }
        };
        size.bounded()
    }

    fn bounded(self) -> Option<Self> {
        (self.numerator_terms.saturating_add(self.denominator_terms) <= 1024
            && self
                .numerator_degree
                .saturating_add(self.denominator_degree)
                <= 64)
            .then_some(self)
    }
}

impl NumeratorAtomExt for AtomView<'_> {
    fn coefficient_list_exact(&self, keys: &[Atom]) -> Option<Vec<(Atom, Atom)>> {
        type Polynomial = MultivariatePolynomial<AtomField, i32>;

        fn collect(atom: AtomView<'_>, keys: &[Atom], template: &Polynomial) -> Option<Polynomial> {
            if let Some(index) = keys.iter().position(|key| key.as_view() == atom) {
                let mut exponents = vec![0; keys.len()];
                exponents[index] = 1;
                return Some(template.monomial(Atom::num(1), exponents));
            }
            if !keys.iter().any(|key| atom.contains(key)) {
                return Some(template.constant(atom.to_owned()));
            }
            match atom {
                AtomView::Add(sum) => {
                    let mut result = template.zero();
                    for term in sum {
                        let term = collect(term, keys, template)?;
                        if result.nterms().checked_add(term.nterms())? > 128 {
                            return None;
                        }
                        result = &result + &term;
                    }
                    Some(result)
                }
                AtomView::Mul(product) => {
                    let mut result = template.one();
                    for factor in product {
                        let factor = collect(factor, keys, template)?;
                        if result.nterms().checked_mul(factor.nterms())? > 128 {
                            return None;
                        }
                        for left in result.exponents_iter() {
                            for right in factor.exponents_iter() {
                                let degree = left.iter().zip(right).try_fold(
                                    0u32,
                                    |degree, (left, right)| {
                                        degree.checked_add(left.checked_add(*right)?.unsigned_abs())
                                    },
                                )?;
                                if degree > 64 {
                                    return None;
                                }
                            }
                        }
                        result = &result * &factor;
                    }
                    Some(result)
                }
                AtomView::Pow(power) => {
                    let index = keys
                        .iter()
                        .position(|key| key.as_view() == power.get_base())?;
                    let exponent = i64::try_from(power.get_exp()).ok()?;
                    if exponent.unsigned_abs() > 64 {
                        return None;
                    }
                    let mut exponents = vec![0; keys.len()];
                    exponents[index] = i32::try_from(exponent).ok()?;
                    Some(template.monomial(Atom::num(1), exponents))
                }
                _ => None,
            }
        }

        if self.get_byte_size() > 1024 * 1024
            || keys.len() > 128
            || keys.iter().enumerate().any(|(index, key)| {
                let admitted = match key.as_view() {
                    AtomView::Var(_) | AtomView::Fun(_) => true,
                    AtomView::Add(_) => {
                        // A literal sum can be a scalar coefficient owner, but
                        // cannot hide free tensor ports or a lazy branch guard.
                        // Inspect its native interface without simplifying or
                        // rewriting the key used by the arithmetic below.
                        let guarded = [
                            Symbol::IF,
                            OrientationID::symbol(),
                            GS.theta,
                            GS.orientation_delta,
                        ]
                        .into_iter()
                        .any(|head| key.contains_symbol(head));
                        !guarded
                            && idenso::CookSettings::indices()
                                .with_mode(idenso::CookMode::ReversibleEncoding)
                                .try_cook(key.as_view())
                                .ok()
                                .and_then(|cooked| {
                                    idenso::tensor::SymbolicTensor::infer(cooked).ok()
                                })
                                .is_some_and(|tensor| tensor.is_scalar())
                    }
                    _ => false,
                };
                !admitted
                    || keys[..index]
                        .iter()
                        .any(|other| key.contains(other) || other.contains(key))
            })
        {
            return None;
        }
        if keys.is_empty() {
            return Some(if self.is_zero() {
                vec![]
            } else {
                vec![(Atom::num(1), self.to_owned())]
            });
        }
        // Native coefficient_list uses AtomField's statistical zero test. That
        // can discard small nonzero coefficients or classify the same stored
        // coefficient inconsistently. Regrouping must use literal exact zeros.
        let field = AtomField {
            statistical_zero_test: false,
            ..AtomField::new()
        };
        // These columns are internal arithmetic handles, not symbolic aliases.
        // A composite key is not a native polynomial indeterminate; reconstruct
        // every output monomial from the original immutable keys below.
        let variables = (0..keys.len())
            .map(symbolica::poly::PolyVariable::Temporary)
            .collect::<Vec<_>>();
        let template = Polynomial::new(&field, None, Arc::new(variables));
        let polynomial = collect(*self, keys, &template)?;
        Some(
            polynomial
                .into_iter()
                .map(|term| {
                    let key = keys
                        .iter()
                        .zip(term.exponents)
                        .fold(Atom::num(1), |product, (key, exponent)| {
                            product * key.pow(*exponent)
                        });
                    (key, term.coefficient.clone())
                })
                .collect(),
        )
    }

    fn cancel_scalar_poles(&self, numerator_keys: &[Atom]) -> Atom {
        // These are optimization budgets, never Taylor truncation cutoffs.
        // Bound family collection as well as the later rational conversion.
        if self.get_byte_size() > 256 * 1024 || numerator_keys.len() > 128 {
            return self.to_owned();
        }
        // Even a fixed residue key can contain inner lazy guards. Never pull a
        // denominator out of its guarded evaluation context.
        if [
            OrientationID::symbol(),
            GS.theta,
            GS.orientation_delta,
            Symbol::IF,
        ]
        .into_iter()
        .any(|symbol| self.contains_symbol(symbol))
        {
            return self.to_owned();
        }
        let started = std::time::Instant::now();
        let contains_family =
            |part: AtomView<'_>| numerator_keys.iter().any(|key| part.contains(key));
        let literal_family =
            |part: AtomView<'_>| numerator_keys.iter().any(|key| key.as_view() == part);
        let monomial_factor = |part: AtomView<'_>| {
            literal_family(part)
                || matches!(part, AtomView::Pow(power)
                if literal_family(power.get_base())
                    && i64::try_from(power.get_exp()).is_ok_and(|n| (0..=64).contains(&n)))
        };
        let mut eligible = numerator_keys
            .iter()
            .all(|key| matches!(key.as_view(), AtomView::Var(_) | AtomView::Fun(_)));
        let mut monomials = 1usize;
        self.visitor(&mut |part| {
            if !eligible || literal_family(part) || !contains_family(part) {
                return false;
            }
            match part {
                AtomView::Add(sum) => {
                    monomials = monomials.saturating_add(sum.get_nargs().saturating_sub(1));
                    eligible = monomials <= 128;
                }
                AtomView::Mul(product) => {
                    let factors = product
                        .iter()
                        .filter(|factor| contains_family(*factor))
                        .collect::<Vec<_>>();
                    // One sum with scalar spectators is a linear jet. Multiple
                    // numerator factors must already be literal monomials.
                    eligible = factors.len() <= 1 || factors.into_iter().all(monomial_factor);
                }
                AtomView::Pow(_) => eligible = monomial_factor(part),
                _ => eligible = false,
            }
            eligible
        });
        if !eligible {
            return self.to_owned();
        }
        let Some(groups) = self.coefficient_list_exact(numerator_keys) else {
            return self.to_owned();
        };
        let mut attempted = 0usize;
        let mut reduced = 0usize;
        let result = Atom::add_many(groups.into_iter().map(|(key, scalar)| {
            let mut has_inverse = false;
            scalar.visitor(&mut |part| match part {
                AtomView::Fun(_) => false,
                AtomView::Pow(power) => {
                    has_inverse |= i64::try_from(power.get_exp()).is_ok_and(|n| n < 0);
                    true
                }
                _ => true,
            });
            let coefficient = if has_inverse
                && scalar.as_view().get_byte_size() <= 16 * 1024
                && ScalarRationalSize::estimate(scalar.as_view()).is_some()
            {
                attempted += 1;
                // Functions are indeterminates to rational normalization: in
                // particular, E(owner, P)^2 is never replaced by P here.
                // In this exact domain, together already cancels common
                // factors; cancel would repeat rational conversion and GCDs.
                let cancelled = scalar.together().collect_compact_factors();
                if cancelled.is_zero()
                    || cancelled.as_view().get_byte_size() < scalar.as_view().get_byte_size()
                {
                    reduced += 1;
                    cancelled
                } else {
                    scalar
                }
            } else {
                scalar
            };
            key * coefficient
        }));
        let result = if result.as_view().get_byte_size() < self.get_byte_size() {
            result
        } else {
            self.to_owned()
        };
        crate::debug_tags!(#generation, #profile, #uv, #numerator, #summary;
            stage = "scalar_pole_cancellation", attempted, reduced,
            input_bytes = self.get_byte_size(), output_bytes = result.as_view().get_byte_size(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Bounded scalar coefficient cancellation"
        );
        result
    }

    fn collect_compact_factors(&self) -> Atom {
        // A complex float serializes the precision of its zero component too.
        // Equal finite coefficients can otherwise become unequal literal keys
        // inside a function. Canonicalize only that exact zero's metadata;
        // neither the nonzero component nor its precision is rounded.
        let mut atom = self.replace_map(|view, _, output| {
            let AtomView::Num(number) = view else {
                return;
            };
            let CoefficientView::Float(real, imaginary) = number.get_coeff_view() else {
                return;
            };
            let real_zero = real.is_zero();
            let imaginary_zero = imaginary.is_zero();
            if real_zero && imaginary_zero {
                **output = Atom::Zero;
            } else if real_zero || imaginary_zero {
                let mut real = real.to_float();
                let mut imaginary = imaginary.to_float();
                if real.prec() != imaginary.prec() {
                    if real_zero {
                        real = symbolica::domains::float::Float::new(imaginary.prec());
                    } else {
                        imaginary = symbolica::domains::float::Float::new(real.prec());
                    }
                    **output = Atom::num(symbolica::domains::float::Complex::new(real, imaginary));
                }
            }
        });
        // Energy powers retain their owner through enclosing Taylor operations.
        // Replacing E(owner, P)^2 by P here would discard the mass-rearrangement
        // boundary needed by an outer U operation.
        // This call is one collection scope. A residue can still contain inner
        // selectors: only share scalar inverses when the entire scope is free
        // of guards, so a finite inactive branch cannot become 0 * infinity.
        let guarded = [
            OrientationID::symbol(),
            GS.theta,
            GS.orientation_delta,
            Symbol::IF,
        ]
        .into_iter()
        .any(|head| atom.contains_symbol(head));
        // Extract only direct numeric content. In particular, implicit unit
        // coefficients participate in the GCD, and term bodies stay literal.
        // Ordering by those bodies makes opposite sums choose the same owner.
        let numeric_content = |value: AtomView<'_>| {
            let AtomView::Add(sum) = value else {
                return value.to_owned();
            };
            let mut terms = Vec::with_capacity(sum.get_nargs());
            for term in sum.iter() {
                let number = match term {
                    AtomView::Num(number) => number.get_coeff_view().to_owned(),
                    AtomView::Mul(product) => match product.get_coefficient() {
                        Some(AtomView::Num(number)) => number.get_coeff_view().to_owned(),
                        _ => Coefficient::one(),
                    },
                    _ => Coefficient::one(),
                };
                if !matches!(number, Coefficient::Complex(_)) {
                    return value.to_owned();
                }
                terms.push((term / Atom::num(&number).as_view(), number));
            }
            terms.sort_by(|left, right| left.0.cmp(&right.0));
            let mut content = terms[0].1.clone();
            for (_, number) in &terms[1..] {
                content = content.gcd(number);
            }
            let first = Atom::num(&terms[0].1) / Atom::num(&content);
            if let AtomView::Num(first) = first.as_view()
                && first.get_coeff_view().to_owned().is_negative()
            {
                content = -content;
            }
            if content.is_one() {
                return value.to_owned();
            }
            let content = Atom::num(content);
            &content * Atom::add_many(sum.iter().map(|term| term / content.as_view()))
        };
        let split_number = |value: &Atom| {
            let number = match value.as_view() {
                AtomView::Mul(product) => product
                    .get_coefficient()
                    .map(|number| number.to_owned())
                    .unwrap_or_else(Atom::one),
                _ => Atom::one(),
            };
            let body = value / &number;
            (number, body)
        };
        loop {
            // Collect only complete factors after Taylor and residue mapping.
            // Opaque powers keep distinct inverse denominators, their owners,
            // and numerator powers intact; functions keep their arguments intact.
            // Wrap original functions too so an input using the temporary head
            // is restored verbatim by the single outer unwrapping pass.
            let opaque = symbol!("gammalooprs::uv::opaque_factor");
            // A positive scalar power and its linear base share one factor owner.
            // The complete summed base stays opaque, so no polynomial expansion of
            // a graph numerator is needed to use A*A^n = A^(n+1).
            let mut powered_sums = Vec::new();
            let share_positive_powers =
                !guarded && !atom.contains_symbol(symbol!("gammalooprs::uv::numerator_family"));
            if share_positive_powers {
                atom.visitor(&mut |part| match part {
                    AtomView::Fun(_) => false,
                    AtomView::Pow(power) => {
                        if i64::try_from(power.get_exp()).is_ok_and(|exponent| exponent > 0)
                            && !part.is_tensorial(StrictTensorFilter::ContainsReps)
                        {
                            let (_, base) = split_number(&numeric_content(power.get_base()));
                            if matches!(base.as_view(), AtomView::Add(_))
                                && !powered_sums.contains(&base)
                            {
                                powered_sums.push(base);
                            }
                        }
                        false
                    }
                    _ => true,
                });
            }
            let mut occurrence = 0usize;
            let mut protected = atom.replace_map(|view, context, out| {
                // Sharing a summed vector across contractions can distribute a
                // vanishing contracted factor across separately evaluated terms.
                // Keep these tensor sums local to their original occurrences.
                let tensor_sum = matches!(view, AtomView::Add(_))
                    && context.parent_type == Some(AtomType::Mul)
                    && view.is_tensorial(StrictTensorFilter::ContainsReps);
                if share_positive_powers
                    && !tensor_sum
                    && !view.is_tensorial(StrictTensorFilter::ContainsReps)
                {
                    if let AtomView::Pow(power) = view
                        && i64::try_from(power.get_exp()).is_ok_and(|exponent| exponent > 0)
                    {
                        let (number, base) = split_number(&numeric_content(power.get_base()));
                        if powered_sums.contains(&base) {
                            **out = (number * function!(opaque, base)).pow(power.get_exp());
                            return;
                        }
                        // Traversal still wraps function bases. Atomic scalar
                        // powers retain their native relationship to linear terms.
                        if matches!(power.get_base(), AtomView::Var(_) | AtomView::Fun(_)) {
                            return;
                        }
                    } else if matches!(view, AtomView::Add(_)) {
                        let (number, base) = split_number(&numeric_content(view));
                        if powered_sums.contains(&base) {
                            **out = number * function!(opaque, base);
                            return;
                        }
                    }
                }
                if tensor_sum || matches!(view, AtomView::Pow(_) | AtomView::Fun(_)) {
                    let mut branch_local = tensor_sum;
                    let protect_inverse =
                        guarded || view.is_tensorial(StrictTensorFilter::ContainsReps);
                    view.visitor(&mut |part| {
                        branch_local |= match part {
                            AtomView::Pow(power) => {
                                protect_inverse
                                    && !i64::try_from(power.get_base_exp().1)
                                        .is_ok_and(|exponent| exponent >= 0)
                            }
                            AtomView::Fun(fun) => [
                                OrientationID::symbol(),
                                GS.theta,
                                GS.orientation_delta,
                                Symbol::IF,
                                symbol!("gammalooprs::uv::numerator_family"),
                            ]
                            .contains(&fun.get_symbol()),
                            _ => false,
                        };
                        !branch_local
                    });
                    // Guarded inverses remain occurrence-local, including inverses
                    // nested in functions or powers. Numerator families and tensor
                    // sums must still form independent contraction networks even
                    // when their scalar denominators can share a factor.
                    // Normalize a scalar integer power before it becomes opaque:
                    // (-2*A)^2 and 4*A^2 have the same complete power owner. Only
                    // its literal numeric prefactor leaves that owner. Inverses
                    // remain opaque to subsequent common-factor collection.
                    let normalized = if let AtomView::Pow(power) = view
                        && !branch_local
                        && !view.is_tensorial(StrictTensorFilter::ContainsReps)
                        && let Ok(exponent) = i64::try_from(power.get_exp())
                    {
                        numeric_content(power.get_base()).pow(exponent)
                    } else {
                        view.to_owned()
                    };
                    let (number, body) = split_number(&normalized);
                    **out = number
                        * if branch_local {
                            occurrence += 1;
                            function!(opaque, body, occurrence)
                        } else {
                            function!(opaque, body)
                        };
                }
            });
            // Original functions and powers are opaque before numeric distribution.
            // This pass can cancel k - (k-p), but cannot enter a denominator or a
            // fractional-power body. The same boundary restores canonical signs
            // after Horner collection introduces new scalar sums.
            let normalize = |value: Atom| {
                if guarded {
                    return value;
                }
                value
                    .expand_num()
                    .replace_map_bottom_up(|view, context, out| {
                        if context.function_level == 0 && matches!(view, AtomView::Add(_)) {
                            **out = numeric_content(view);
                        }
                    })
            };
            protected = normalize(protected);
            loop {
                let collected =
                    normalize(protected.collect_horner::<Symbol>(None).collect_factors());
                // Strictly decreasing size bounds the iteration and prevents two
                // equivalent factor orders from alternating indefinitely.
                if collected == protected
                    || collected.as_view().get_byte_size() >= protected.as_view().get_byte_size()
                {
                    break;
                }
                protected = collected;
            }
            let result = protected.replace_map(|view, _, out| {
                if let AtomView::Fun(fun) = view
                    && fun.get_symbol() == opaque
                {
                    out.set_from_view(
                        &fun.iter()
                            .next()
                            .expect("opaque factor has an expression argument"),
                    );
                }
            });
            // Restoring a factored base can expose a new linear cancellation.
            // Re-enter the same protected scope only after a strict size decrease;
            // no alternative variable orders or grouping searches are attempted.
            if result.as_view().get_byte_size() >= atom.as_view().get_byte_size() {
                return result;
            }
            atom = result;
        }
    }

    fn series_preserving_factors(
        &self,
        variable: Symbol,
        expansion_point: AtomView<'_>,
        depth: i64,
        numerator_family_keys: &[Atom],
    ) -> color_eyre::Result<Atom> {
        for (index, key) in numerator_family_keys.iter().enumerate() {
            ensure!(
                matches!(key.as_view(), AtomView::Var(_) | AtomView::Fun(_))
                    && !key.contains_symbol(variable),
                "numerator-family key must be a Taylor-independent symbol or function: {key}"
            );
            ensure!(
                numerator_family_keys[..index]
                    .iter()
                    .all(|other| !key.contains(other) && !other.contains(key)),
                "numerator-family keys must be distinct and cannot contain one another: {key}"
            );
        }
        let contains_family =
            |atom: AtomView<'_>| numerator_family_keys.iter().any(|key| atom.contains(key));
        let mut independent = Atom::num(1);
        let mut dependent = self.to_owned();
        if let AtomView::Mul(product) = self {
            dependent = Atom::num(1);
            for factor in product.iter() {
                if factor.contains_symbol(variable) || contains_family(factor) {
                    dependent *= factor;
                } else {
                    independent *= factor;
                }
            }
        }
        // Keep native whole-product precision lifting and cross-arm zero
        // detection before regrouping the surviving coefficient families.
        let started = std::time::Instant::now();
        let native_series =
            dependent.series(variable, expansion_point, SeriesDepth::absolute(depth))?;
        crate::debug_tags!(#generation, #profile, #uv, #numerator, #summary;
            stage = "factor_preserving_native_series_done", depth,
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            dependent_bytes = dependent.as_view().get_byte_size(),
            coefficient_count = native_series.terms().count(),
            coefficient_bytes = native_series.terms().map(|(_, c)| c.as_view().get_byte_size()).sum::<usize>(),
            max_coefficient_bytes = native_series.terms().map(|(_, c)| c.as_view().get_byte_size()).max().unwrap_or(0),
            "Native series before atom conversion and factor collection"
        );
        let series = native_series.to_atom();
        if numerator_family_keys.is_empty() {
            return Ok((independent * series).collect_compact_factors());
        }
        // Reject nonlinear or hidden occurrences before polynomial collection
        // can distribute powers or products of coefficient families.
        let mut invalid = None;
        series.visitor(&mut |atom| {
            if invalid.is_some()
                || numerator_family_keys
                    .iter()
                    .any(|key| key.as_view() == atom)
            {
                return false;
            }
            let valid = match atom {
                AtomView::Mul(product) => {
                    product
                        .iter()
                        .filter(|factor| contains_family(*factor))
                        .count()
                        <= 1
                }
                AtomView::Fun(_) | AtomView::Pow(_) => !contains_family(atom),
                _ => true,
            };
            if !valid {
                invalid = Some(atom.to_owned());
            }
            valid
        });
        if let Some(invalid) = invalid {
            bail!(
                "numerator-family keys must occur linearly outside functions and inverses: {invalid}"
            );
        }
        let mut grouped = Atom::Zero;
        for (key, coefficient) in series.coefficient_list::<u32>(numerator_family_keys) {
            ensure!(
                key.is_one() || numerator_family_keys.contains(&key),
                "series must be linear in its numerator-family keys, found {key}"
            );
            ensure!(
                !contains_family(coefficient.as_view()),
                "numerator-family key remains hidden in series coefficient {coefficient}"
            );
            grouped += key * coefficient;
        }
        Ok((independent * grouped).collect_compact_factors())
    }

    fn map_mink_dim<'a>(&self, dim: impl Into<AtomOrView<'a>>) -> Atom {
        self.replace(Minkowski {}.to_symbolic([W_.d_, W_.a___]))
            .with(Minkowski {}.to_symbolic([dim.into().into_owned(), Atom::var(W_.a___)]))
    }

    fn unwrap_function(&self, symbol: Symbol) -> Atom {
        self.replace(function!(symbol, W_.a___)).with(W_.a___)
    }

    fn parse_into_net(&self) -> Result<ParsingNet, ParsingNetError> {
        ParsingNet::try_from_view(
            *self,
            TENSORLIB.read().unwrap().deref(),
            &ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode {
                        inner_products: false,
                        expand_schoonship: true,
                        expand_inside_chains: true,
                    },
                    trace: true,
                    chain: true,
                },
                ..Default::default()
            },
        )
    }

    // fn parse_into_only_lib_net<T: TensorLibraryData + Clone + Default>(
    //     &self,
    //     one: T,
    //     zero: T,
    // ) -> Result<ParsingNet, ParsingNetError> {
    //     let mut lib = hep_lib(one, zero);

    //     ParsingNet::try_from_view(*self, &lib);
    // }
}

#[cfg(test)]
mod tests {
    use idenso::IndexTooling;
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::{
        network::tags::SPENSO_TAG,
        structure::{
            representation::{Minkowski, RepName},
            slot::{DummyAind, IsAbstractSlot, Slot},
        },
    };
    use symbolica::{
        atom::{Atom, AtomCore},
        function, parse_lit, symbol,
    };

    use crate::{initialisation::test_initialise, numerator::aind::Aind, utils::GS};

    use super::NumeratorAtomExt;

    #[test]
    fn compact_factors_canonicalize_only_literal_zero_precision_metadata() {
        use symbolica::{
            atom::AtomView,
            coefficient::CoefficientView,
            domains::float::{Complex, Float},
        };

        test_initialise().unwrap();
        let opaque = symbol!("compact_zero_precision::opaque");
        let spectator = parse_lit!((x + y) ^ 11 * (z + w) ^ (-3));
        let tiny = Atom::num((1, 2)).pow(120);
        for imaginary_nonzero in [false, true] {
            let nonzero = Float::with_val(128, 5);
            let expected = if imaginary_nonzero {
                Atom::num(Complex::new(Float::new(128), nonzero.clone()))
            } else {
                Atom::num(Complex::new(nonzero.clone(), Float::new(128)))
            };
            for zero_precision in [53, 96, 256] {
                let coefficient = if imaginary_nonzero {
                    Atom::num(Complex::new(Float::new(zero_precision), nonzero.clone()))
                } else {
                    Atom::num(Complex::new(nonzero.clone(), Float::new(zero_precision)))
                };
                assert_ne!(coefficient, expected);
                assert_eq!(
                    coefficient.to_canonical_string(),
                    expected.to_canonical_string()
                );
                let normalized = coefficient.collect_compact_factors();
                assert_eq!(normalized, expected);
                assert_eq!(normalized.collect_compact_factors(), normalized);
                let (AtomView::Num(before), AtomView::Num(after)) =
                    (coefficient.as_view(), normalized.as_view())
                else {
                    panic!("literal coefficient")
                };
                let (
                    CoefficientView::Float(before_re, before_im),
                    CoefficientView::Float(after_re, after_im),
                ) = (before.get_coeff_view(), after.get_coeff_view())
                else {
                    panic!("floating coefficient")
                };
                assert_eq!(
                    if imaginary_nonzero {
                        before_im.0
                    } else {
                        before_re.0
                    },
                    if imaginary_nonzero {
                        after_im.0
                    } else {
                        after_re.0
                    },
                    "all nonzero value and precision bytes must survive unchanged"
                );
                let source =
                    &spectator * (function!(opaque, &coefficient) - function!(opaque, &expected));
                assert!(!source.is_zero());
                assert!(source.collect_compact_factors().is_zero());
                assert_eq!(
                    (&spectator * function!(opaque, &coefficient)).collect_compact_factors(),
                    &spectator * function!(opaque, &expected),
                    "powers and inverse spectators retain their factorization"
                );
            }
            let direction = if imaginary_nonzero {
                Atom::i()
            } else {
                Atom::one()
            };
            let distinct = &expected + &tiny * direction;
            assert_ne!(distinct, expected, "the finite perturbation is represented");
            let normalized_distinct = distinct.collect_compact_factors();
            assert_ne!(normalized_distinct, expected);
            let (AtomView::Num(before), AtomView::Num(after)) =
                (distinct.as_view(), normalized_distinct.as_view())
            else {
                panic!("literal perturbed coefficient")
            };
            let (
                CoefficientView::Float(before_re, before_im),
                CoefficientView::Float(after_re, after_im),
            ) = (before.get_coeff_view(), after.get_coeff_view())
            else {
                panic!("floating perturbed coefficient")
            };
            assert_eq!(
                if imaginary_nonzero {
                    before_im.0
                } else {
                    before_re.0
                },
                if imaginary_nonzero {
                    after_im.0
                } else {
                    after_re.0
                },
                "the small finite perturbation remains bit-exact"
            );
            assert!(
                !(&spectator * (function!(opaque, &distinct) - function!(opaque, &expected)))
                    .collect_compact_factors()
                    .is_zero(),
                "normalization must not erase a nonzero difference"
            );
        }
    }

    #[test]
    fn exact_scalar_sum_keys_regroup_nested_coefficients_without_opening_binomials() {
        test_initialise().unwrap();
        let [a, b, c, d, e, f, g] = [
            parse_lit!(p1 - p2),
            parse_lit!(p1 + p2),
            parse_lit!(k0 + k3),
            parse_lit!(k1 + k2),
            parse_lit!(p0 + p3),
            parse_lit!(k1 - k2),
            parse_lit!(-k0 + k3),
        ];
        let source: Atom = 4 * &a * (&b * &c - &d * &e)
            + 4 * &b * (&f * &e - &a * &c)
            + 4 * &e * (&g * &e + &a * &d)
            - 4 * &e * (&g * &e + &b * &f);
        assert!(
            !source.is_zero(),
            "the input retains its nested factorization"
        );
        let spectator = parse_lit!((x + y) ^ 11) * idenso::color_cas!(2, idenso::cof!(3));
        for value in [source.clone(), &spectator * &source] {
            let before = value.clone();
            assert!(
                value
                    .coefficient_list_exact(std::slice::from_ref(&e))
                    .unwrap()
                    .is_empty()
            );
            assert_eq!(value, before);
        }
        // A nearby nonzero expression must retain the complete binomial and
        // spectator. No polynomial identity of their bodies is requested.
        let retained = &source + &spectator * &a * &e;
        assert_eq!(
            retained
                .coefficient_list_exact(std::slice::from_ref(&e))
                .unwrap(),
            vec![(e, spectator * a)]
        );
    }

    #[test]
    fn exact_scalar_sum_keys_preserve_signed_powers_and_closed_dummy_scopes() {
        test_initialise().unwrap();
        let slot = spenso::mink!(4, Atom::from(Aind::new_dummy()));
        let key = Atom::one() + GS.emr_mom(EdgeIndex(1), &slot) * GS.emr_mom(EdgeIndex(2), &slot);
        let spectator = parse_lit!((x + y) ^ 11)
            * function!(
                symbol!("literal_scalar_factor::opaque"),
                parse_lit!((z + w) ^ 9)
            );
        let source: Atom = &spectator * (key.pow(-2) + 3 * key.pow(3));
        let before = source.clone();
        let rows = source
            .coefficient_list_exact(std::slice::from_ref(&key))
            .unwrap();
        let coefficients = rows
            .iter()
            .cloned()
            .collect::<std::collections::BTreeMap<_, _>>();
        assert_eq!(coefficients.len(), 2);
        assert_eq!(coefficients.get(&key.pow(-2)), Some(&spectator));
        assert_eq!(coefficients.get(&key.pow(3)), Some(&(3 * &spectator)));
        let reconstructed = Atom::add_many(
            rows.iter()
                .map(|(factor, coefficient)| factor * coefficient),
        );
        assert!(
            (&reconstructed - &source)
                .coefficient_list_exact(std::slice::from_ref(&key))
                .unwrap()
                .is_empty(),
            "the literal-key rows exactly reconstruct the original coefficient"
        );
        assert_eq!(source, before);
        for exponent in [-2, -1, 2, 3, 11] {
            let value = key.pow(exponent);
            assert_eq!(
                value
                    .coefficient_list_exact(std::slice::from_ref(&key))
                    .unwrap(),
                vec![(value, Atom::one())]
            );
        }
    }

    #[test]
    fn exact_scalar_sum_keys_reject_open_ports_guards_and_non_laurent_dependence() {
        test_initialise().unwrap();
        let slot = spenso::mink!(4, Atom::from(Aind::new_dummy()));
        let open = GS.emr_mom(EdgeIndex(1), &slot) + GS.emr_mom(EdgeIndex(2), &slot);
        assert!(
            open.coefficient_list_exact(std::slice::from_ref(&open))
                .is_none()
        );
        let key = parse_lit!(literal_scalar_factor::a + literal_scalar_factor::b);
        let extra = parse_lit!(literal_scalar_factor::c);
        for source in [
            (&key * &extra + 1).pow(-1),
            (&key * &extra + 1).pow(2),
            key.pow(Atom::num((1, 2))),
            key.pow(Atom::num((-1, 2))),
            function!(symbol!("literal_scalar_factor::hidden"), &key),
        ] {
            assert!(
                source
                    .coefficient_list_exact(std::slice::from_ref(&key))
                    .is_none(),
                "unsupported selected dependence must remain unproven: {source}"
            );
        }
        for guard in [
            crate::cff::expression::OrientationID(3).atom(),
            function!(GS.theta, &extra),
            function!(GS.orientation_delta, &extra, 1),
            function!(symbolica::atom::Symbol::IF, &extra, 1, 0),
        ] {
            let guarded = &key + guard;
            assert!(
                guarded
                    .coefficient_list_exact(std::slice::from_ref(&guarded))
                    .is_none()
            );
        }
        assert!(
            key.coefficient_list_exact(&[key.clone(), parse_lit!(literal_scalar_factor::a)])
                .is_none(),
            "overlapping keys cannot become independent polynomial columns"
        );
    }

    #[test]
    fn compact_factors_normalize_numeric_content_before_integer_power_ownership() {
        test_initialise().unwrap();
        let base = parse_lit!(p0 * k0 - p1 * k1 - p2 * k2 - p3 * k3);
        let spectator = parse_lit!((a + b) ^ 11 * (c + d));
        for coefficient in [Atom::one(), Atom::i() * Atom::num((4, 3))] {
            let difference = &coefficient
                * (-Atom::num(2) * base.pow(2)
                    + Atom::num((1, 2)) * (-Atom::num(2) * &base).expand_num().pow(2));
            assert!(
                (&spectator * difference)
                    .collect_compact_factors()
                    .is_zero()
            );
        }
        let opposite = (-&base).expand_num();
        assert!(
            (base.pow(2) * opposite + base.pow(3))
                .collect_compact_factors()
                .is_zero()
        );
        let retained = (&spectator * base.pow(3)).collect_compact_factors();
        assert!(retained.contains(&parse_lit!((a + b) ^ 11)));
        assert!(retained.contains(&parse_lit!(c + d)));
    }

    #[test]
    fn compact_factors_keep_signed_integer_and_fractional_power_domains_distinct() {
        test_initialise().unwrap();
        let base = parse_lit!(a - b);
        let opposite = parse_lit!(-a + b);
        for exponent in [-2, -1, 2, 3, 11] {
            let left = base.pow(exponent);
            let right = Atom::num(if exponent % 2 == 0 { 1 } else { -1 }) * opposite.pow(exponent);
            let normalized = left.collect_compact_factors();
            assert_eq!(
                normalized,
                right.collect_compact_factors(),
                "power {exponent}"
            );
            assert_eq!(normalized.collect_compact_factors(), normalized);
        }
        let fractional = parse_lit!((-a + b) ^ (1 / 2));
        assert_eq!(fractional.collect_compact_factors(), fractional);
        let separate_denominators = parse_lit!(x / (a + b) + y / (c + d));
        assert_eq!(
            separate_denominators.collect_compact_factors(),
            separate_denominators
        );
    }

    #[test]
    fn compact_factors_preserve_atomic_mass_power_relations() {
        test_initialise().unwrap();
        let difference =
            parse_lit!(m * (v - m ^ 2) - m * (m ^ 2 - k2) - m * (v - m ^ 2 + k2) + m ^ 3);
        let spectator = parse_lit!((a + b) ^ 11 * opaque(c + 2 * (d + e)) / (x + y));
        assert!(
            (&spectator * &difference)
                .collect_compact_factors()
                .is_zero()
        );
        assert_eq!(
            (&spectator * (&difference + 1)).collect_compact_factors(),
            spectator
        );
    }

    #[test]
    fn compact_factors_share_scalar_inverses_only_without_guards() {
        test_initialise().unwrap();
        let a = parse_lit!(a);
        let b = parse_lit!(b);
        // The complete inverse stays opaque, including its summed denominator.
        let denominator = parse_lit!((D1 + D2) ^ -1);
        let input = &a * &denominator + &b * &denominator;
        let compact = input.collect_compact_factors();
        assert_eq!(compact, (&a + &b) * &denominator);
        assert_eq!(compact.collect_compact_factors(), compact);
        for guard in [
            crate::cff::expression::OrientationID(3).atom(),
            function!(GS.theta, &a),
            function!(GS.orientation_delta, &a, &b),
            function!(symbolica::atom::Symbol::IF, &a, &b, 0),
        ] {
            let guarded = &guard * &a * &denominator + &b * &denominator;
            assert_eq!(guarded.collect_compact_factors(), guarded);
        }
    }

    #[test]
    fn compact_factors_share_denominators_without_merging_numerator_networks() {
        test_initialise().unwrap();
        let numerator = function!(symbol!("gammalooprs::uv::numerator_family"), 1);
        let denominator = parse_lit!((D1 + D2) ^ -1);
        let first = &numerator * parse_lit!(a);
        let second = &numerator * parse_lit!(b);
        let input = &first * &denominator + &second * &denominator;
        assert_eq!(
            input.collect_compact_factors(),
            (first + second) * &denominator
        );
        // An inverse containing a tensor contraction is not a scalar factor.
        let tensor = spenso::tensor!(compact_test, spenso::mink!(4, mu));
        let inverse = (&tensor + spenso::tensor!(other_compact_test, spenso::mink!(4, mu))).pow(-1);
        let tensor_input = &inverse * parse_lit!(a) + &inverse * parse_lit!(b);
        assert_eq!(tensor_input.collect_compact_factors(), tensor_input);
    }

    #[test]
    fn series_preserves_independent_tensor_factors() {
        test_initialise().unwrap();
        let t = symbol!("series_t");
        let slot: Slot<Minkowski, Aind> = Minkowski {}.new_rep(4).slot(Aind::new_dummy());
        let tensor =
            GS.emr_vec(EdgeIndex(0), slot.to_atom()) * GS.emr_vec(EdgeIndex(1), slot.to_atom());
        for spectator in [parse_lit!(g), tensor * parse_lit!(u + v)] {
            let source = &spectator * parse_lit!((a + b * series_t) * (c + d * series_t));
            let result = source
                .series_preserving_factors(t, Atom::Zero.as_view(), 1, &[])
                .unwrap();
            assert_eq!(
                result,
                &spectator * parse_lit!(a * c + (a * d + b * c) * series_t)
            );
            assert_eq!(
                result
                    .series_preserving_factors(t, Atom::Zero.as_view(), 1, &[])
                    .unwrap(),
                result
            );
        }
    }

    #[test]
    fn series_preserves_laurent_order_with_mapped_numerator() {
        let t = symbol!("series_t");
        let source = parse_lit!(g * (a + b * series_t) * (c + d * series_t) / series_t ^ 2);
        for (depth, expected) in [
            (
                -1,
                parse_lit!(g * (a * c / series_t ^ 2 + (a * d + b * c) / series_t)),
            ),
            (
                0,
                parse_lit!(g * (a * c / series_t ^ 2 + (a * d + b * c) / series_t + b * d)),
            ),
        ] {
            assert_eq!(
                source
                    .series_preserving_factors(t, Atom::Zero.as_view(), depth, &[])
                    .unwrap(),
                expected
            );
        }
    }

    #[test]
    fn series_preserves_native_additive_zero() {
        let t = symbol!("series_t");
        for source in [
            Atom::Zero,
            parse_lit!(g * (a + b * series_t) / series_t - g * a / series_t - g * b),
        ] {
            let result = source
                .series_preserving_factors(t, Atom::Zero.as_view(), 0, &[])
                .unwrap();
            assert!(result.is_zero());
            assert_eq!(result, source.series(t, Atom::Zero, 0).unwrap().to_atom());
            assert!(result.replace(t).with(Atom::num(1)).is_zero());
        }
    }

    #[test]
    fn series_preserves_nonzero_expansion_point() {
        let t = symbol!("series_t");
        let source = parse_lit!(g * (a + b * series_t) * (c + d * series_t));
        let result = source
            .series_preserving_factors(t, Atom::num(2).as_view(), 1, &[])
            .unwrap();
        assert_eq!(
            result,
            parse_lit!(
                g * ((a + 2 * b) * (c + 2 * d)
                    + (series_t - 2) * (b * (c + 2 * d) + d * (a + 2 * b)))
            )
        );
    }

    #[test]
    fn series_compacts_spectators_without_merging_additive_denominators() {
        let t = symbol!("series_t");
        let source =
            parse_lit!(s1 * g * (a + b * series_t) / D1 + s2 * g * (c + d * series_t) / D2);
        for (depth, expected) in [
            (0, parse_lit!(g * (a * s1 / D1 + c * s2 / D2))),
            (
                1,
                parse_lit!(
                    g * (a * s1 / D1 + c * s2 / D2 + series_t * (b * s1 / D1 + d * s2 / D2))
                ),
            ),
        ] {
            let result = source
                .series_preserving_factors(t, Atom::Zero.as_view(), depth, &[])
                .unwrap();
            // This is a small scalar oracle, with no graph numerator. Different
            // Horner layouts are valid; retain the separate inverse factors.
            assert!((&result - expected).expand().is_zero());
            for denominator in [parse_lit!(D1 ^ -1), parse_lit!(D2 ^ -1)] {
                assert!(result.contains(&denominator));
            }
            let native = source.series(t, Atom::Zero, depth).unwrap().to_atom();
            assert!(result.as_view().get_byte_size() <= native.as_view().get_byte_size());
        }
    }

    #[test]
    fn series_regroups_owned_numerator_coefficient_families() {
        let t = symbol!("series_t");
        let family = symbol!("series_owned_coefficient"; Scalar);
        let keys = [0, 1, 2].map(|order| function!(family, order));
        let [n0, n1, n2] = &keys;
        let spectator = parse_lit!(g * (u + v));
        let source = &spectator
            * (n0 + n1 * Atom::var(t) + n2 * Atom::var(t).pow(2))
            * parse_lit!((a + b * series_t) / series_t ^ 2);
        for (depth, expected) in [
            (
                -1,
                n0 * parse_lit!(a / series_t ^ 2 + b / series_t) + n1 * parse_lit!(a / series_t),
            ),
            (
                0,
                n0 * parse_lit!(a / series_t ^ 2 + b / series_t)
                    + n1 * parse_lit!(a / series_t + b)
                    + n2 * parse_lit!(a),
            ),
        ] {
            let result = source
                .series_preserving_factors(t, Atom::Zero.as_view(), depth, &keys)
                .unwrap();
            assert_eq!(result, &spectator * expected);
            assert_eq!(
                result
                    .series_preserving_factors(t, Atom::Zero.as_view(), depth, &keys)
                    .unwrap(),
                result
            );
            let native = source.series(t, Atom::Zero, depth).unwrap().to_atom();
            assert!((result - native).expand().is_zero());
        }
    }

    #[test]
    fn series_family_regrouping_preserves_zero_and_nonzero_center() {
        let t = symbol!("series_t");
        let family = symbol!("series_owned_coefficient"; Scalar);
        let keys = [0, 1, 2].map(|order| function!(family, order));
        let [n0, n1, n2] = &keys;
        let source = n0 * parse_lit!((a + b * series_t) / series_t)
            - n0 * parse_lit!(a / series_t)
            - n0 * parse_lit!(b);
        for zero in [Atom::Zero, source] {
            let result = zero
                .series_preserving_factors(t, Atom::Zero.as_view(), 0, &keys)
                .unwrap();
            assert!(result.is_zero());
            assert!(result.replace(t).with(Atom::num(1)).is_zero());
        }

        let shift = parse_lit!(series_t - 2);
        let source = (n0 + n1 * &shift + n2 * shift.pow(2))
            * (parse_lit!(a) + parse_lit!(b) * &shift)
            / shift.pow(2);
        let result = source
            .series_preserving_factors(t, Atom::num(2).as_view(), 0, &keys)
            .unwrap();
        assert_eq!(
            result,
            n0 * (parse_lit!(a) / shift.pow(2) + parse_lit!(b) / &shift)
                + n1 * (parse_lit!(a) / &shift + parse_lit!(b))
                + n2 * parse_lit!(a)
        );
        assert!(
            (result - source.series(t, Atom::num(2), 0).unwrap().to_atom())
                .expand()
                .is_zero()
        );
    }

    #[test]
    fn series_family_specialization_preserves_ose_derivatives() {
        let t = symbol!("series_t");
        let family = symbol!("series_ose_coefficient"; Scalar);
        let keys = [0, 1, 2].map(|order| {
            function!(
                family,
                order,
                parse_lit!(sigma),
                parse_lit!(z),
                parse_lit!(tau),
                parse_lit!(w)
            )
        });
        // Scalar diagnostic: each energy retains its own sign, integer node
        // and fixed shift. Nonconstant OSEs must contribute their derivatives.
        let q = parse_lit!(sigma * (1 + series_t) ^ (1 / 2) + z * M + x);
        let r = parse_lit!(tau * (4 + 2 * series_t) ^ (1 / 2) + w * M + y);
        let q0 = parse_lit!(sigma + z * M + x);
        let r0 = parse_lit!(2 * tau + w * M + y);
        let coefficients = [
            &q0 * &r0,
            &q0 * parse_lit!(tau / 2) + parse_lit!(sigma / 2) * &r0,
            -&q0 * parse_lit!(tau / 16) + parse_lit!(sigma * tau / 4) - parse_lit!(sigma / 8) * &r0,
        ];
        let truncated_numerator =
            &keys[0] + &keys[1] * Atom::var(t) + &keys[2] * Atom::var(t).pow(2);
        let denominator = parse_lit!(series_t ^ 2 * (1 - series_t));
        let retained = (truncated_numerator / &denominator)
            .series_preserving_factors(t, Atom::Zero.as_view(), 0, &keys)
            .unwrap();
        assert_eq!(
            retained,
            &keys[0] * parse_lit!(1 / series_t ^ 2 + 1 / series_t + 1)
                + &keys[1] * parse_lit!(1 / series_t + 1)
                + &keys[2]
        );
        let rebuilt = keys
            .iter()
            .zip(&coefficients)
            .fold(retained, |atom, (key, body)| {
                atom.replace(key.to_pattern()).with(body.to_pattern())
            });
        let source = q * r / denominator;
        let native = source.series(t, Atom::Zero, 0).unwrap().to_atom();
        assert!((&rebuilt - native).expand().is_zero());
        // The first row makes the first energy vanish at the expansion point.
        // Other rows keep the two repeated-occurrence nodes independent.
        for values in [
            [1, -1, 1, 2, 1, 0, 0],
            [-1, 2, 1, -3, 2, 3, -1],
            [1, 0, -1, 2, -2, 1, 4],
        ] {
            let variables = ["sigma", "z", "tau", "w", "M", "x", "y"].map(|name| symbol!(name));
            let (specialized, rebuilt) = variables.into_iter().zip(values).fold(
                (source.clone(), rebuilt.clone()),
                |(source, result), (variable, value)| {
                    (
                        source.replace(variable).with(Atom::num(value)),
                        result.replace(variable).with(Atom::num(value)),
                    )
                },
            );
            assert!(
                (rebuilt - specialized.series(t, Atom::Zero, 0).unwrap().to_atom())
                    .expand()
                    .is_zero()
            );
        }
    }

    #[test]
    fn series_rejects_invalid_numerator_families() {
        let t = symbol!("series_t");
        let key = parse_lit!(n0);
        for keys in [
            vec![parse_lit!(n0 + n1)],
            vec![parse_lit!(n0(series_t))],
            vec![key.clone(), key.clone()],
            vec![key.clone(), parse_lit!(outer(n0))],
        ] {
            assert!(
                key.series_preserving_factors(t, Atom::Zero.as_view(), 0, &keys)
                    .is_err()
            );
        }
        for source in [
            parse_lit!(n0 ^ 2),
            parse_lit!((n0 + 1) ^ 1000),
            parse_lit!(1 / n0),
            parse_lit!(outer(n0)),
        ] {
            assert!(
                source
                    .series_preserving_factors(
                        t,
                        Atom::Zero.as_view(),
                        0,
                        std::slice::from_ref(&key),
                    )
                    .is_err()
            );
        }
    }

    #[test]
    fn series_with_families_preserves_native_errors() {
        let t = symbol!("series_t");
        let source = parse_lit!(n0 * singular_argument(1 / series_t));
        let native = source.series(t, Atom::Zero, 0).unwrap_err();
        let error = source
            .series_preserving_factors(t, Atom::Zero.as_view(), 0, &[parse_lit!(n0)])
            .unwrap_err();
        assert_eq!(
            error.downcast_ref::<symbolica::poly::series::SeriesError>(),
            Some(&native)
        );
    }

    #[test]
    fn dummy_parsing() {
        test_initialise().unwrap();

        let e_mass = parse_lit!(M_e);

        let m2 = &e_mass * &e_mass;

        let mink: Slot<Minkowski, Aind> = Minkowski {}.new_rep(4).slot(Aind::new_dummy());

        let e = EdgeIndex(0);

        let sqrt = symbol!("sqrt_scalar", tag = SPENSO_TAG.broadcast);

        let a = function!(
            sqrt,
            (GS.emr_vec(e, mink.to_atom()) * GS.emr_vec(e, mink.to_atom()) + m2).pow(Atom::num(2))
        );

        let net = a.parse_into_net().unwrap();

        println!("{}", net.dot_pretty())
    }

    #[test]
    fn canonizations() {
        test_initialise().unwrap();

        let a = parse_lit!(
            ((-2 * spenso::projp(
                spenso::bis(4, gammalooprs::edge(0)),
                spenso::bis(4, gammalooprs::hedge(2))
            ) + spenso::projm(
                spenso::bis(4, gammalooprs::edge(0)),
                spenso::bis(4, gammalooprs::hedge(2))
            )) * -1𝑖
                / 6
                * UFO::sw
                ^ 2 + -1𝑖 / 2 * UFO::cw
                ^ 2 * spenso::projm(
                    spenso::bis(4, gammalooprs::edge(0)),
                    spenso::bis(4, gammalooprs::hedge(2))
                ))
                * ((-2
                    * spenso::projp(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(8))
                    )
                    + spenso::projm(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(8))
                    ))
                    * -1𝑖
                    / 6
                    * UFO::sw
                    ^ 2 + -1𝑖 / 2 * UFO::cw
                    ^ 2 * spenso::projm(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(8))
                    ))
                * (-1 * UFO::MZ
                    ^ 2 * spenso::g(
                        spenso::mink(4, gammalooprs::hedge(4)),
                        spenso::mink(4, gammalooprs::hedge(5))
                    ) + gammalooprs::K(1, spenso::mink(4, gammalooprs::hedge(4)))
                        * gammalooprs::K(1, spenso::mink(4, gammalooprs::hedge(5))))
                * (-1 * gammalooprs::P(0, spenso::mink(4, gammalooprs::edge(5, 1)))
                    + gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(5, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(5, 1))))
                * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(1, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(1, 1))))
                * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(4, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(4, 1))))
                * 1
                / 3
                * UFO::MZ
                ^ (-2) * UFO::cw
                ^ (-2) * UFO::ee
                ^ 4 * UFO::sw
                ^ (-2)
                    * gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(2, 1)))
                    * gammalooprs::e(0, spenso::mink(4, gammalooprs::hedge(1)))
                    * gammalooprs::ebar(0, spenso::mink(4, gammalooprs::hedge(0)))
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(10)),
                        spenso::bis(4, gammalooprs::hedge(11)),
                        spenso::mink(4, gammalooprs::edge(5, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(11)),
                        spenso::bis(4, gammalooprs::hedge(7)),
                        spenso::mink(4, gammalooprs::hedge(1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(2)),
                        spenso::bis(4, gammalooprs::hedge(3)),
                        spenso::mink(4, gammalooprs::edge(2, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(3)),
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::mink(4, gammalooprs::hedge(5))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(6)),
                        spenso::bis(4, gammalooprs::edge(0)),
                        spenso::mink(4, gammalooprs::hedge(4))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(7)),
                        spenso::bis(4, gammalooprs::hedge(6)),
                        spenso::mink(4, gammalooprs::edge(4, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(8)),
                        spenso::bis(4, gammalooprs::hedge(9)),
                        spenso::mink(4, gammalooprs::edge(1, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(9)),
                        spenso::bis(4, gammalooprs::hedge(10)),
                        spenso::mink(4, gammalooprs::hedge(0))
                    )
        );
        println!("a:{}", a);
        let b = parse_lit!(
            ((-2 * spenso::projp(
                spenso::bis(4, gammalooprs::edge(0)),
                spenso::bis(4, gammalooprs::hedge(6))
            ) + spenso::projm(
                spenso::bis(4, gammalooprs::edge(0)),
                spenso::bis(4, gammalooprs::hedge(6))
            )) * -1𝑖
                / 6
                * UFO::sw
                ^ 2 + -1𝑖 / 2 * UFO::cw
                ^ 2 * spenso::projm(
                    spenso::bis(4, gammalooprs::edge(0)),
                    spenso::bis(4, gammalooprs::hedge(6))
                ))
                * ((-2
                    * spenso::projp(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(3))
                    )
                    + spenso::projm(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(3))
                    ))
                    * -1𝑖
                    / 6
                    * UFO::sw
                    ^ 2 + -1𝑖 / 2 * UFO::cw
                    ^ 2 * spenso::projm(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(3))
                    ))
                * (-1 * UFO::MZ
                    ^ 2 * spenso::g(
                        spenso::mink(4, gammalooprs::hedge(4)),
                        spenso::mink(4, gammalooprs::hedge(5))
                    ) + gammalooprs::K(1, spenso::mink(4, gammalooprs::hedge(4)))
                        * gammalooprs::K(1, spenso::mink(4, gammalooprs::hedge(5))))
                * (-1 * gammalooprs::P(0, spenso::mink(4, gammalooprs::edge(5, 1)))
                    + gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(5, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(5, 1))))
                * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(1, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(1, 1))))
                * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(4, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(4, 1))))
                * 1
                / 3
                * UFO::MZ
                ^ (-2) * UFO::cw
                ^ (-2) * UFO::ee
                ^ 4 * UFO::sw
                ^ (-2)
                    * gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(2, 1)))
                    * gammalooprs::e(0, spenso::mink(4, gammalooprs::hedge(1)))
                    * gammalooprs::ebar(0, spenso::mink(4, gammalooprs::hedge(0)))
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(10)),
                        spenso::bis(4, gammalooprs::hedge(9)),
                        spenso::mink(4, gammalooprs::hedge(0))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(11)),
                        spenso::bis(4, gammalooprs::hedge(10)),
                        spenso::mink(4, gammalooprs::edge(5, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(2)),
                        spenso::bis(4, gammalooprs::edge(0)),
                        spenso::mink(4, gammalooprs::hedge(4))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(3)),
                        spenso::bis(4, gammalooprs::hedge(2)),
                        spenso::mink(4, gammalooprs::edge(2, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(6)),
                        spenso::bis(4, gammalooprs::hedge(7)),
                        spenso::mink(4, gammalooprs::edge(4, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(7)),
                        spenso::bis(4, gammalooprs::hedge(11)),
                        spenso::mink(4, gammalooprs::hedge(1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(8)),
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::mink(4, gammalooprs::hedge(5))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(9)),
                        spenso::bis(4, gammalooprs::hedge(8)),
                        spenso::mink(4, gammalooprs::edge(1, 1))
                    )
        );

        println!("b:{}", b);

        println!("ratio:{}", &a / &b);

        let ac = a
            .canonize(Aind::Dummy)
            .expect("test expression should canonicalize");
        let bc = b
            .canonize(Aind::Dummy)
            .expect("test expression should canonicalize");
        println!("ac:{}", ac);
        println!("bc:{}", bc);
        println!("ratio canonized:{}", ac / bc);
    }

    // #[test]
    // fn test_can() {
    //     let a = parse_lit!(T(a, b, c) * T(c, d, e) * T(d, b, f)(K(e) + P(e)) * (K(f) + P(f)));
    //     let b = parse_lit!(T(a, d) * T(d, c));

    //     let indices = vec![
    //         (parse_lit!(a), 1),
    //         (parse_lit!(a), 1),
    //         (parse_lit!(a), 1),
    //         (parse_lit!(b), 1),
    //         (parse_lit!(b), 1),
    //         (parse_lit!(c), 1),
    //         (parse_lit!(d), 1),
    //         (parse_lit!(e), 1),
    //         (parse_lit!(f), 1),
    //     ];

    //     let ac = a.canonize_tensors(&indices);
    //     println!("{}", ac.unwrap());
    //     let bc = b.canonize_tensors(&indices);
    //     println!("{}", bc.unwrap());
    // }
}

#[cfg(test)]
#[path = "scalar_pole_tests.rs"]
mod scalar_pole_tests;
