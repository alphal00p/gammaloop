//! Combine scalar spectators with their integer common-mass scaling before IBP.

use super::numerator::MatchedScalarInput;
use rustred::algebra::{
    Coefficient, CoefficientContext, ExactAlgebraLimits, SYMBOLICA_COEFFICIENT_EXPONENT_LIMIT,
};
use std::collections::BTreeSet;
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::coefficient::{Coefficient as SymbolicaCoefficient, CoefficientView};
use symbolica::domains::float::Complex;
use symbolica::prelude::{Q, Z};

/// Formal coefficientwise conjugation fixes native parameters. Only rational
/// expressions with exact Gaussian coefficients and real symbolic denominators
/// are admitted; no analytic reality assumptions are made about unknown atoms.
fn formal_coefficient_reality(part: AtomView<'_>, parameters: &BTreeSet<Symbol>) -> Option<bool> {
    match part {
        AtomView::Num(number) => match number.get_coeff_view() {
            value @ (CoefficientView::Natural(..) | CoefficientView::Large(..)) => {
                Some(value.is_real())
            }
            _ => None,
        },
        AtomView::Var(variable) => parameters.contains(&variable.get_symbol()).then_some(true),
        AtomView::Add(sum) => sum.iter().try_fold(true, |real, term| {
            Some(formal_coefficient_reality(term, parameters)? && real)
        }),
        AtomView::Mul(product) => product.iter().try_fold(true, |real, factor| {
            Some(formal_coefficient_reality(factor, parameters)? && real)
        }),
        AtomView::Pow(power) => {
            let AtomView::Num(exponent) = power.get_exp() else {
                return None;
            };
            let exponent_view = exponent.get_coeff_view();
            if !matches!(
                exponent_view,
                CoefficientView::Natural(..) | CoefficientView::Large(..)
            ) || !exponent_view.is_integer()
            {
                return None;
            }
            let exponent = i64::try_from(power.get_exp()).ok()?;
            let real = formal_coefficient_reality(power.get_base(), parameters)?;
            (exponent >= 0 || real).then_some(real)
        }
        _ => None,
    }
}

/// Conservative numerator/denominator degree envelope, including intermediate
/// products. Failure merely keeps the factor opaque: this is the existing
/// Native exponent representation boundary, never a user work quota.
fn representable_degree(part: AtomView<'_>) -> Option<u64> {
    let maximum = u64::from(SYMBOLICA_COEFFICIENT_EXPONENT_LIMIT);
    let degree = match part {
        AtomView::Num(_) => 0,
        AtomView::Var(_) => 1,
        AtomView::Add(sum) => sum.iter().try_fold(0u64, |degree, term| {
            degree.checked_add(representable_degree(term)?)
        })?,
        AtomView::Mul(product) => product.iter().try_fold(0u64, |degree, factor| {
            degree.checked_add(representable_degree(factor)?)
        })?,
        AtomView::Pow(power) => {
            let exponent = i64::try_from(power.get_exp()).ok()?.unsigned_abs();
            if exponent > maximum {
                return None;
            }
            representable_degree(power.get_base())?.checked_mul(exponent)?
        }
        _ => return None,
    };
    (degree <= maximum).then_some(degree)
}

impl MatchedScalarInput {
    /// Exact numeric common-mass factors belong in the native weight, rather
    /// than fragmenting spectator identities. Native coefficients are real
    /// Q(d), so Q(i) factors use the two canonical carriers one and i. Formal
    /// polynomial numerators in native parameters may also be split over real
    /// Q(d) denominators. Complex symbolic denominators, unknown symbols,
    /// functions and precision-tracking floats remain opaque spectators.
    pub(super) fn channel_coefficients(
        spectator: Atom,
        coefficient: &Coefficient,
        context: &CoefficientContext,
        limits: ExactAlgebraLimits,
    ) -> Result<Vec<(Atom, Coefficient)>, super::RustRedEvaluationError> {
        let factors = match spectator.as_view() {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            view => vec![view],
        };
        let mut parameters = None;
        let mut scalar = Atom::one();
        let mut formal = Atom::one();
        let mut formal_degree = 0u64;
        for factor in factors {
            let eligible = match factor {
                AtomView::Num(number) => matches!(
                    number.get_coeff_view(),
                    CoefficientView::Natural(..) | CoefficientView::Large(..)
                ),
                AtomView::Fun(_) => false,
                _ => {
                    let parameters = parameters.get_or_insert_with(|| {
                        context
                            .parameter_names()
                            .iter()
                            .flat_map(|name| {
                                context
                                    .parameter(name)
                                    .unwrap()
                                    .to_expression()
                                    .get_all_symbols(true)
                            })
                            .collect::<BTreeSet<_>>()
                    });
                    formal_coefficient_reality(factor, parameters).is_some()
                }
            };
            let degree = if eligible {
                representable_degree(factor)
            } else {
                None
            };
            let combined_degree = degree.and_then(|degree| formal_degree.checked_add(degree));
            if let Some(degree) = combined_degree
                && degree <= u64::from(SYMBOLICA_COEFFICIENT_EXPONENT_LIMIT)
            {
                formal *= factor.to_owned();
                formal_degree = degree;
            } else {
                scalar *= factor.to_owned();
            }
        }
        let imaginary_unit = Atom::num(SymbolicaCoefficient::Complex(Complex::new(
            0.into(),
            1.into(),
        )));
        let parts = if let AtomView::Num(number) = formal.as_view() {
            let SymbolicaCoefficient::Complex(number) = number.get_coeff_view().to_owned() else {
                unreachable!("Natural and Large views contain exact complex rationals");
            };
            [Atom::num(number.re), Atom::num(number.im)]
        } else {
            // Project expanded numerators before rational collection. Asking
            // together() to collect a complex expression first could factor
            // even a real denominator d²+1 into (d-i)(d+i), which Q(d) rejects.
            let formal = formal.expand();
            let conjugate = formal.replace_map(|part, _, output| {
                if let AtomView::Num(number) = part
                    && !number.get_coeff_view().is_real()
                {
                    **output = Atom::num(number.get_coeff_view().to_owned().conjugate());
                }
            });
            [
                ((&formal + &conjugate) / Atom::num(2)).expand().together(),
                ((formal - conjugate) / (Atom::num(2) * &imaginary_unit))
                    .expand()
                    .together(),
            ]
        };
        let mut channels = Vec::with_capacity(2);
        for (part, carrier) in parts.into_iter().zip([Atom::one(), imaginary_unit]) {
            if part.is_zero() {
                continue;
            }
            if representable_degree(part.as_view()).is_none() {
                return Ok(vec![(spectator, coefficient.clone())]);
            }
            let rational = part
                .try_to_rational_polynomial(&Q, &Z, Some(context.one().get_variables().clone()))
                .map_err(super::data::failure)?;
            let weight = context
                .try_mul(coefficient, &rational, limits)
                .map_err(super::RustRedEvaluationError::from_native)?;
            channels.push((&scalar * carrier, weight));
        }
        Ok(channels)
    }

    pub(super) fn mass_channel(
        spectator: &Atom,
        mass_squared: &Atom,
        exponent: i64,
        anchor: i64,
        diagnostic_max_degree: Option<usize>,
    ) -> (Atom, i64) {
        let Some(shift) = exponent.checked_sub(anchor) else {
            return (spectator.clone(), exponent);
        };
        if diagnostic_max_degree.is_some_and(|max_degree| {
            !usize::try_from(shift.unsigned_abs()).is_ok_and(|shift| shift <= max_degree)
        }) {
            return (spectator.clone(), exponent);
        }
        let mut exact = true;
        for expression in [spectator, mass_squared] {
            expression.visitor(&mut |part| {
                if let AtomView::Num(number) = part
                    && matches!(number.get_coeff_view(), CoefficientView::Float(..))
                {
                    exact = false;
                }
                exact
            });
        }
        if !exact {
            // Preserve raw floating inputs and their original precision. Staged
            // source identities remain symbolic and need no reevaluation here.
            return (spectator.clone(), exponent);
        }
        // Use one existing key's exponent as origin: raising a numeric mass to
        // -sum(a) could create an enormous rational for large dotted inputs.
        // An optional diagnostic can bound this difference; ordinary grouping
        // relies on exact integer scaling and the shared RAM/time budget.
        // Atom normalization combines powers without expanding compound masses;
        // (m²)^(Ld/2+anchor) stays outside the unchanged d-only native coefficients.
        (spectator * mass_squared.clone().pow(shift), anchor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;
    use symbolica::{domains::rational::Rational, parse};

    #[test]
    fn nonunit_numeric_mass_shifts_do_not_fragment_ordered_channels() {
        let context = CoefficientContext::try_new(["d"]).unwrap();
        let (spectator, exponent) =
            MatchedScalarInput::mass_channel(&parse!("x"), &parse!("7/3"), -6, -5, None);
        let channels = MatchedScalarInput::channel_coefficients(
            spectator,
            &context.one(),
            &context,
            ExactAlgebraLimits::default(),
        )
        .unwrap();
        assert_eq!(channels.len(), 1);
        let (spectator, coefficient) = &channels[0];
        assert_eq!(spectator, &parse!("x"));
        assert_eq!(exponent, -5);
        assert_eq!(coefficient.to_expression(), parse!("3/7"));
    }

    #[test]
    fn exact_complex_prefactors_use_two_real_native_weights() {
        let context = CoefficientContext::try_new(["d"]).unwrap();
        let weight = context.parameter("d").unwrap();
        for spectator in [
            parse!("-8*𝑖*x/3"),
            parse!("(2/3-8*𝑖/3)*x"),
            parse!("(123456789012345678901234567890/17+98765432109876543210987654321*𝑖/19)*x"),
            Atom::Zero,
        ] {
            let channels = MatchedScalarInput::channel_coefficients(
                spectator.clone(),
                &weight,
                &context,
                ExactAlgebraLimits::default(),
            )
            .unwrap();
            assert!(channels.len() <= 2);
            let restored = channels
                .iter()
                .map(|(carrier, coefficient)| carrier * coefficient.to_expression())
                .sum::<Atom>();
            assert!(
                (restored - spectator * weight.to_expression())
                    .expand()
                    .together()
                    .is_zero()
            );
            for (carrier, _) in channels {
                assert!(carrier == parse!("x") || carrier == parse!("𝑖*x"));
            }
        }
    }

    #[test]
    fn different_complex_prefactors_fuse_and_cancel_in_the_same_carriers() {
        let context = CoefficientContext::try_new(["d"]).unwrap();
        let mut combined = BTreeMap::<Atom, Atom>::new();
        // Two independently emitted contributions to the same integral:
        // (a+i*b) I + (-a+i*c) I = i*(b+c) I.
        for spectator in [parse!("(2/3+3*𝑖/5)*x"), parse!("(-2/3+7*𝑖/5)*x")] {
            for (carrier, coefficient) in MatchedScalarInput::channel_coefficients(
                spectator,
                &context.one(),
                &context,
                ExactAlgebraLimits::default(),
            )
            .unwrap()
            {
                *combined.entry(carrier).or_insert(Atom::Zero) += coefficient.to_expression();
            }
        }
        combined.retain(|_, coefficient| !coefficient.is_zero());
        assert_eq!(combined, BTreeMap::from([(parse!("𝑖*x"), Atom::num(2))]));
    }

    #[test]
    fn complex_parameter_numerators_share_carriers_with_separated_forms() {
        let context = CoefficientContext::try_new(["d"]).unwrap();
        let d = context.parameter("d").unwrap().to_expression();
        let imaginary = parse!("𝑖");
        let split = |spectator| {
            MatchedScalarInput::channel_coefficients(
                spectator,
                &context.one(),
                &context,
                ExactAlgebraLimits::default(),
            )
            .unwrap()
        };
        for denominator in [&d + Atom::one(), d.clone().pow(2) + Atom::one()] {
            let combined = split((&d + &imaginary) / &denominator);
            let mut separated = BTreeMap::<Atom, Atom>::new();
            for part in [&d / &denominator, &imaginary / &denominator] {
                for (carrier, coefficient) in split(part) {
                    *separated.entry(carrier).or_insert(Atom::Zero) += coefficient.to_expression();
                }
            }
            assert_eq!(combined.len(), 2);
            assert_eq!(separated.len(), 2);
            for (carrier, coefficient) in combined {
                assert!(carrier == Atom::one() || carrier == imaginary);
                assert!(
                    (coefficient.to_expression() - &separated[&carrier])
                        .together()
                        .expand()
                        .is_zero()
                );
            }
        }
    }

    #[test]
    fn complex_symbolic_denominators_and_unknown_factors_remain_opaque() {
        let context = CoefficientContext::try_new(["d"]).unwrap();
        let d = context.parameter("d").unwrap().to_expression();
        for spectator in [
            (&d + parse!("𝑖")).pow(-1),
            (&d + parse!("𝑖")).pow(-2),
            parse!("x+𝑖"),
            symbolica::symbol!("opaque").call(d.clone()),
            d.pow(parse!("1/2")),
        ] {
            let channels = MatchedScalarInput::channel_coefficients(
                spectator.clone(),
                &context.one(),
                &context,
                ExactAlgebraLimits::default(),
            )
            .unwrap();
            assert_eq!(channels.len(), 1);
            assert_eq!(channels[0].0, spectator);
            assert_eq!(channels[0].1, context.one());
        }
    }

    #[test]
    fn powers_outside_native_exponent_representation_remain_opaque() {
        let context = CoefficientContext::try_new(["d"]).unwrap();
        let d = context.parameter("d").unwrap().to_expression();
        let beyond = u64::from(SYMBOLICA_COEFFICIENT_EXPONENT_LIMIT) + 1;
        for spectator in [
            d.clone().pow(Atom::num(beyond)),
            d.clone().pow(Atom::num(u64::from(u32::MAX) + 1)),
            (d.pow(Atom::num(beyond / 2)) + Atom::one()).pow(2),
        ] {
            let channels = MatchedScalarInput::channel_coefficients(
                spectator.clone(),
                &context.one(),
                &context,
                ExactAlgebraLimits::default(),
            )
            .unwrap();
            assert_eq!(channels.len(), 1);
            assert_eq!(channels[0].0, spectator);
            assert_eq!(channels[0].1, context.one());
        }
    }

    #[test]
    fn complex_mass_shifts_preserve_the_exact_full_coefficient() {
        let context = CoefficientContext::try_new(["d"]).unwrap();
        let mass = parse!("(2+𝑖)/3");
        let (spectator, exponent) =
            MatchedScalarInput::mass_channel(&parse!("x"), &mass, -6, -5, None);
        let restored = MatchedScalarInput::channel_coefficients(
            spectator,
            &context.one(),
            &context,
            ExactAlgebraLimits::default(),
        )
        .unwrap()
        .into_iter()
        .map(|(carrier, coefficient)| carrier * coefficient.to_expression())
        .sum::<Atom>();
        assert_eq!(exponent, -5);
        assert!(
            (restored * mass.clone().pow(exponent) - parse!("x") * mass.pow(-6))
                .together()
                .expand()
                .is_zero()
        );
    }

    #[test]
    fn complex_float_spectators_keep_their_original_precision() {
        use symbolica::domains::float::Float;
        let context = CoefficientContext::try_new(["d"]).unwrap();
        let spectator = Atom::num(Complex::new(
            Float::parse("1.234567890123456789", Some(73)).unwrap(),
            Float::parse("2.345678901234567891", Some(71)).unwrap(),
        )) * parse!("x");
        let channels = MatchedScalarInput::channel_coefficients(
            spectator.clone(),
            &context.one(),
            &context,
            ExactAlgebraLimits::default(),
        )
        .unwrap();
        assert_eq!(channels.len(), 1);
        assert_eq!(channels[0].0, spectator);
        assert_eq!(channels[0].1, context.one());
    }

    fn channel_sum(mass: &Atom, rows: &[(i64, Atom, Atom)]) -> (Atom, usize) {
        let mut channels = BTreeMap::<(Atom, i64), Atom>::new();
        let anchor = -rows[0].0;
        for (degree, spectator, value) in rows {
            let channel =
                MatchedScalarInput::mass_channel(spectator, mass, -*degree, anchor, Some(512));
            *channels.entry(channel).or_insert(Atom::Zero) += value;
        }
        let n = channels.len();
        let result = channels
            .into_iter()
            .map(|((spectator, exponent), value)| {
                spectator * mass.clone().pow(parse!("10-5*ep") + exponent) * value
            })
            .sum();
        (result, n)
    }

    #[test]
    fn homogeneous_channels_preserve_the_full_five_loop_laurent_scaling() {
        let mass = parse!("mu");
        let rows = [
            (10, Atom::one(), parse!("2/ep^2+3/ep+5")),
            (11, mass.clone(), parse!("-3*(-1/ep^2+7)")),
            (12, mass.clone().pow(2), parse!("22/ep")),
        ];
        let (result, channels) = channel_sum(&mass, &rows);
        assert_eq!(channels, 1);
        let ep = symbolica::symbol!("ep");
        let expanded = result
            .series(ep, Atom::Zero.as_view(), Rational::from(0))
            .unwrap()
            .to_atom();
        let expected = parse!("5/ep^2+(25-25*log(mu))/ep-16-125*log(mu)+125/2*log(mu)^2");
        assert!((expanded - expected).expand().together().is_zero());
        // The diagnostic has engineering degree nine instead of ten, hence
        // one additional mass squared. No rank inference enters grouping.
        let diagnostic = rows.map(|(d, s, f)| (d - 1, s, f));
        let (other, channels) = channel_sum(&mass, &diagnostic);
        assert_eq!(channels, 1);
        assert!((other - &mass * result).expand().together().is_zero());
    }

    #[test]
    fn compound_and_nonunit_masses_agree_with_independent_per_integral_scaling() {
        for mass in [parse!("7/3"), parse!("x^2+y^2"), parse!("4*x^2")] {
            let rows = [
                (10, Atom::one(), parse!("f1")),
                (11, mass.clone(), parse!("-3*f2")),
                (12, mass.clone().pow(2), parse!("2*f3")),
            ];
            let before = rows
                .iter()
                .map(|(d, s, f)| s * mass.clone().pow(parse!("10-5*ep") - *d) * f)
                .sum::<Atom>();
            let (after, _) = channel_sum(&mass, &rows);
            // An exact integer dimension specialization also catches incorrect
            // scaling of nonunit numeric and compound masses without requiring
            // branch-sensitive symbolic power simplifications from the CAS.
            for epsilon in [parse!("0"), parse!("1/5"), parse!("-2/5")] {
                let residual = (before.clone() - &after)
                    .replace(parse!("ep"))
                    .with(epsilon);
                assert!(residual.expand().together().is_zero());
            }
        }
    }

    #[test]
    fn canceled_unknowns_share_a_channel_with_exact_rational_spectators() {
        let mass = parse!("mu");
        let rows = [
            (11, &mass / parse!("d-5"), parse!("unknown")),
            (10, parse!("1/(d-5)"), parse!("-unknown")),
        ];
        let (result, channels) = channel_sum(&mass, &rows);
        assert_eq!(channels, 1);
        assert!(result.is_zero());
    }

    #[test]
    fn compound_powers_stay_compact_and_staged_sources_are_not_evaluated() {
        let mass = parse!("x^2+y^2");
        let channel =
            MatchedScalarInput::mass_channel(&mass.clone().pow(500), &mass, -510, 0, Some(512));
        assert_eq!(channel, (mass.pow(-10), 0));
        let source = parse!("source_float(17)");
        let mass = parse!("mu");
        assert_eq!(
            MatchedScalarInput::mass_channel(
                &(&source * mass.clone().pow(2)),
                &mass,
                -12,
                -10,
                Some(2)
            ),
            (source.clone(), -10)
        );
        // No huge rational is materialized for a large absolute index; unsafe
        // differences and shifts beyond an explicit diagnostic grant stay put.
        let numeric_mass = parse!("7/3");
        for (exponent, anchor, limit) in [
            (i64::MIN, i64::MIN, 6),
            (i64::MIN, i64::MAX, 6),
            (-1_000_000_000_000, 0, 6),
            (-12, -10, 1),
        ] {
            assert_eq!(
                MatchedScalarInput::mass_channel(
                    &source,
                    &numeric_mass,
                    exponent,
                    anchor,
                    Some(limit)
                ),
                (source.clone(), exponent)
            );
        }
    }

    #[test]
    fn ordinary_symbolic_mass_channels_ignore_dormant_degree_thresholds() {
        let source = parse!("source_float(17)");
        let mass = parse!("mu");
        let exponent = -1_000_000_000_000_i64;
        let expected = (&source * mass.clone().pow(exponent), 0);
        assert_eq!(
            MatchedScalarInput::mass_channel(&source, &mass, exponent, 0, None),
            expected
        );
        assert_eq!(
            MatchedScalarInput::mass_channel(&source, &mass, exponent, 0, Some(6)),
            (source.clone(), exponent)
        );
        assert_eq!(
            MatchedScalarInput::mass_channel(&source, &mass, i64::MIN, i64::MAX, None),
            (source, i64::MIN)
        );
    }

    #[test]
    fn raw_float_components_keep_the_original_mass_factorization_and_precision() {
        use symbolica::domains::float::Float;
        let value = Atom::num(Float::parse("1.234567890123456789", Some(73)).unwrap());
        for (spectator, mass) in [
            (value.clone(), parse!("mu")),
            (parse!("x"), &value + parse!("mu")),
            (
                symbolica::symbol!("wrapped").call(value.clone()),
                parse!("mu"),
            ),
        ] {
            let original = (spectator.clone(), -12);
            for diagnostic_max_degree in [None, Some(6)] {
                assert_eq!(
                    MatchedScalarInput::mass_channel(
                        &spectator,
                        &mass,
                        -12,
                        -10,
                        diagnostic_max_degree
                    ),
                    original
                );
            }
        }
    }
}
