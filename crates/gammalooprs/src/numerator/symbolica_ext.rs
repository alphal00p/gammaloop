use std::{ops::Deref, sync::Arc};

use color_eyre::eyre::{bail, ensure};
use idenso::{
    color::{CS, ColorSimplifier},
    representations::{ColorAdjoint, ColorFundamental},
};
use spenso::{
    network::parsing::{
        AtomStructureExt, ParseSettings, SchoonschipExpansionMode, ShorthandParsing,
        StrictTensorFilter,
    },
    structure::representation::{Minkowski, RepName},
};

use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomType, AtomView, Symbol},
    coefficient::CoefficientView,
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
    /// exact zero tests. Key-independent subtrees remain opaque coefficients.
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

    fn to_param_color(&self) -> Atom;
    // fn wrap_color(&self, symbol: Symbol) -> Atom;
    fn kill_color(&self) -> Atom;

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

    fn to_param_color(&self) -> Atom {
        self.as_view().to_param_color()
    }
    fn kill_color(&self) -> Atom {
        self.wrap_color(GS.killing_func)
    }

    fn map_mink_dim<'a>(&self, dim: impl Into<AtomOrView<'a>>) -> Atom {
        self.as_view().map_mink_dim(dim)
    }
    // fn wrap_color(&self, symbol: Symbol) -> Atom {
    //     self.as_view().wrap_color(symbol)
    // }

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
                if fun.get_symbol() != GS.energy_surface
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
                !matches!(key.as_view(), AtomView::Var(_) | AtomView::Fun(_))
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
        let variables = keys
            .iter()
            .cloned()
            .map(TryInto::try_into)
            .collect::<Result<Vec<_>, _>>()
            .ok()?;
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
                let cancelled = scalar.together().cancel().collect_compact_factors();
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
        // Energy powers retain their owner through enclosing Taylor operations.
        // Replacing E(owner, P)^2 by P here would discard the mass-rearrangement
        // boundary needed by an outer U operation.
        let atom = self;
        // Collect only complete factors after Taylor and residue mapping.
        // Opaque powers keep distinct inverse denominators, their owners,
        // and numerator powers intact; functions keep their arguments intact.
        // Wrap original functions too so an input using the temporary head
        // is restored verbatim by the single outer unwrapping pass.
        let opaque = symbol!("gammalooprs::uv::opaque_factor");
        let mut occurrence = 0usize;
        let mut protected = atom.replace_map(|view, context, out| {
            // Sharing a summed vector across contractions can distribute a
            // vanishing contracted factor across separately evaluated terms.
            // Keep these tensor sums local to their original occurrences.
            let tensor_sum = matches!(view, AtomView::Add(_))
                && context.parent_type == Some(AtomType::Mul)
                && view.is_tensorial(StrictTensorFilter::ContainsReps);
            if tensor_sum || matches!(view, AtomView::Pow(_) | AtomView::Fun(_)) {
                let mut branch_local = tensor_sum;
                view.visitor(&mut |part| {
                    branch_local |= match part {
                        AtomView::Pow(power) => !i64::try_from(power.get_base_exp().1)
                            .is_ok_and(|exponent| exponent >= 0),
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
                // Identical inverses must remain inside their branch guards,
                // including inverses nested in a function or numerator power.
                // Selectors stay with their contributions so independent
                // scalar contractions do not become one combined network.
                **out = if branch_local {
                    occurrence += 1;
                    function!(opaque, view, occurrence)
                } else {
                    function!(opaque, view)
                };
            }
        });
        loop {
            let collected = protected.collect_horner::<Symbol>(None).collect_factors();
            // Strictly decreasing size bounds the iteration and prevents two
            // equivalent factor orders from alternating indefinitely.
            if collected == protected
                || collected.as_view().get_byte_size() >= protected.as_view().get_byte_size()
            {
                break;
            }
            protected = collected;
        }
        protected.replace_map(|view, _, out| {
            if let AtomView::Fun(fun) = view
                && fun.get_symbol() == opaque
            {
                out.set_from_view(
                    &fun.iter()
                        .next()
                        .expect("opaque factor has an expression argument"),
                );
            }
        })
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

    fn kill_color(&self) -> Atom {
        self.wrap_color(GS.killing_func)
    }

    fn to_param_color(&self) -> Atom {
        let adj = ColorAdjoint {};
        let fund = ColorFundamental {};
        self.replace(adj.to_symbolic([W_.d_, W_.a_]))
            .with(adj.to_symbolic([CS.nc * CS.nc - 1, Atom::var(W_.a_)]))
            .replace(fund.to_symbolic([W_.d_, W_.a_]))
            .with(fund.to_symbolic([CS.nc, W_.a_]))
    }
    fn map_mink_dim<'a>(&self, dim: impl Into<AtomOrView<'a>>) -> Atom {
        self.replace(Minkowski {}.to_symbolic([W_.d_, W_.a___]))
            .with(Minkowski {}.to_symbolic([dim.into().into_owned(), Atom::var(W_.a___)]))
    }
    // fn wrap_color(&self, symbol: Symbol) -> Atom {
    //     self.expand_color()
    //         .into_iter()
    //         .fold(Atom::Zero, |a, (c, s)| a + function!(symbol, c) * s)
    // }

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

    use crate::{
        dot,
        graph::{FeynmanGraph, Graph, parse::IntoGraph},
        initialisation::test_initialise,
        numerator::aind::Aind,
        utils::GS,
        uv::UltravioletGraph,
    };

    use super::NumeratorAtomExt;

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
    fn canonize_color() {
        test_initialise().unwrap();
        let gls: Vec<Graph> = dot!(
            digraph{
            num = "1";

            ext0 [style=invis];
            2:0-> ext0 [id=0 dir=none is_cut=0  particle="a"];
            ext1 [style=invis];
            ext1-> 3:1 [id=1 dir=none is_cut=0  particle="a"];
            0:2-> 1:3 [id=2   particle="d"];
            0:4-> 1:5 [id=3 dir=none   particle="g"];
            3:6-> 0:7 [id=4   particle="d"];
            1:8-> 2:9 [id=5   particle="d"];
            2:10-> 3:11 [id=6   particle="d"];
        }

        digraph GL8{
            num = 1;
        0[int_id=V_74];
        1[int_id=V_74];
        2[int_id=V_71];
        3[int_id=V_71];
        ext0 [style=invis];
        2:0-> ext0 [id=0 dir=none is_cut=0  particle=a];
        ext1 [style=invis];
        ext1-> 3:1 [id=1 dir=none is_cut=0  particle=a];
        0:2-> 1:3 [id=2   particle=d];
        0:4-> 1:5 [id=3 dir=none   particle=g];
        0:6-> 3:7 [id=4 dir=back   particle="d~"];
        1:8-> 2:9 [id=5   particle=d];
        2:10-> 3:11 [id=6   particle=d];
        }

        )
        .unwrap();

        for g in gls {
            let mut numerator = g.numerator(&g.no_dummy(), &g.empty_subgraph());

            // TODO Check if we include overall factor in main
            numerator.state.expr *= &g.global_prefactor.num * &g.global_prefactor.projector; // * &gl5.overall_factor;
            // numerator.state.expr = numerator.state.expr.replace_multiple(&cpl_reps);

            let numerator_color_simplified = numerator
                .clone()
                .color_simplify()
                .get_single_atom()
                .unwrap()
                .canonize(Aind::Dummy)
                .expect("test expression should canonicalize");

            println!("numerator_color_simplified:{numerator_color_simplified}");
            println!("numerator:{}", numerator.state.expr);
        }
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
