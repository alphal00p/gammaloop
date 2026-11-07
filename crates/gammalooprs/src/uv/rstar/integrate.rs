//! Integration of the remaining, Taylor-expanded single-mass vacuum graphs.
//!
//! The existing vacuum tensor projector acts before the marked line is cut.
//! Its external tensor coefficients accompany the scalar numerator through the
//! exact p-integral reduction and the final Laurent expansion.

use std::collections::BTreeSet;

use eyre::{Result, WrapErr, ensure, eyre};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function,
    id::Replacement,
    poly::series::SeriesDepth,
};
use vakint::{
    EvaluationMethod, RustRedMasterMode, TensorReductionMethod, Vakint, VakintExpression,
    VakintSettings, p_integral::MassiveClosure, symbols::S, vakint_parse,
};

use crate::uv::approx::integrated::Integrated;

/// Integrate explicit vacuum topologies with exactly one massive denominator.
///
/// The mass and power must already be those of the selected surviving line.
/// Massless vacuum terms vanish; local terms (`integral = 1`) are retained.
/// Proper UV counterterms are supplied by the ordinary massive R operation.
/// The marked momentum needs a unit coefficient in the supplied loop basis;
/// the remaining massless graph must be covered by Vakint's p-integral backend.
/// The requested Laurent order applies after multiplying insertion coefficients:
/// their poles increase the depth needed from the remaining integral.
pub fn integrate_single_mass(
    vakint: &Vakint,
    settings: &VakintSettings,
    input: &VakintExpression,
) -> Result<Atom> {
    let mut tensor_settings = settings.clone();
    tensor_settings.tensor_reduction_method = TensorReductionMethod::FeynKit;
    tensor_settings.use_dot_product_notation = true;
    let mut result = Atom::Zero;
    for source in &input.0 {
        if source.numerator.is_zero() {
            continue;
        }
        if source.integral == Atom::one() {
            ensure!(
                !source.numerator.contains_symbol(S.k),
                "a local R* insertion cannot contain an unintegrated loop momentum"
            );
            result += &source.numerator;
            continue;
        }
        let lines = VacuumLine::parse(&source.integral)?;
        let massive = lines
            .iter()
            .filter(|line| !line.mass_squared.is_zero() && line.power > 0)
            .collect::<Vec<_>>();
        if massive.is_empty() {
            continue;
        }
        ensure!(
            massive.len() == 1,
            "an R* remainder needs exactly one massive propagator, found {}",
            massive.len()
        );
        ensure!(
            lines
                .iter()
                .all(|line| { line.mass_squared.is_zero() || std::ptr::eq(line, massive[0]) }),
            "an R* remainder contains another massive propagator with nonpositive power"
        );
        let cut = VacuumCut::new(&lines, massive[0], &source.numerator)?;
        // A massive self-loop times a massless vacuum component is scaleless.
        if cut.scaleless {
            continue;
        }
        let mut term = source.clone();
        term.tensor_reduce(vakint, &tensor_settings)
            .wrap_err("while applying the existing vacuum tensor projector to an R* remainder")?;
        // Lorentz representation metadata is constant during the following
        // Laurent expansion; only scalar evanescent coefficients may expand.
        term.numerator = Integrated::restore_lorentz_dimensions(&term.numerator);
        if settings.verify_numerator_identification {
            // Strict matching admits only authenticated scalar symbols. Keep
            // arbitrary tensor spectators outside that check, as before.
            result += cut.evaluate_split(vakint, settings, &term.numerator)?;
            continue;
        }
        // The p-integral owner separates homogeneous loop monomials from their
        // tensor spectators. Keeping their coefficients here lets its Laurent
        // expansion include insertion poles before truncating the full product.
        let numerator = term.numerator.replace_multiple(&cut.routing).expand();
        let input = numerator * &cut.integrand;
        result += cut
            .closure
            .evaluate(vakint, settings, input.as_view())
            .wrap_err("while integrating an R* remainder through a massless p-integral")?;
    }
    Ok(result.expand())
}

struct VacuumLine {
    id: Atom,
    ends: Atom,
    momentum: Atom,
    mass_squared: Atom,
    power: i64,
}

impl VacuumLine {
    fn parse(integral: &Atom) -> Result<Vec<Self>> {
        let AtomView::Fun(topology) = integral.as_view() else {
            return Err(eyre!("R* integration requires an explicit topo(prop(...))"));
        };
        ensure!(
            topology.get_symbol() == S.topo && topology.get_nargs() == 1,
            "R* integration requires an explicit topo(prop(...))"
        );
        let argument = topology.iter().next().unwrap();
        let factors = match argument {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            _ => vec![argument],
        };
        factors
            .into_iter()
            .map(|factor| {
                let AtomView::Fun(prop) = factor else {
                    return Err(eyre!("invalid R* propagator {factor}"));
                };
                ensure!(
                    prop.get_symbol() == S.prop && prop.get_nargs() == 5,
                    "invalid R* propagator {factor}"
                );
                let arguments = prop.iter().collect::<Vec<_>>();
                ensure!(
                    !arguments[2].contains_symbol(S.p),
                    "the R* remainder must be Taylor-expanded to a vacuum before integration"
                );
                Ok(Self {
                    id: arguments[0].to_owned(),
                    ends: arguments[1].to_owned(),
                    momentum: arguments[2].to_owned(),
                    mass_squared: arguments[3].to_owned(),
                    power: arguments[4]
                        .try_into()
                        .map_err(|_| eyre!("an R* propagator must have integer power"))?,
                })
            })
            .collect()
    }
}

struct VacuumCut {
    closure: MassiveClosure,
    routing: Vec<Replacement>,
    integrand: Atom,
    scaleless: bool,
    loops: usize,
}

impl VacuumCut {
    fn evaluate_split(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        numerator: &Atom,
    ) -> Result<Atom> {
        let mut invariants = BTreeSet::new();
        numerator.visitor(&mut |part| {
            if let AtomView::Fun(function) = part
                && function.get_symbol() == S.dot
                && part.contains_symbol(S.k)
            {
                invariants.insert(part.to_owned());
                return false;
            }
            true
        });
        let mut result = Atom::Zero;
        for (scalar, coefficient) in
            numerator.coefficient_list::<i16>(&invariants.into_iter().collect::<Vec<_>>())
        {
            ensure!(
                !coefficient.contains_symbol(S.k),
                "the vacuum tensor projector left an uncontracted loop momentum: {coefficient}"
            );
            let scalar = scalar.replace_multiple(&self.routing).expand() * &self.integrand;
            result += self
                .evaluate(vakint, settings, &scalar, &coefficient)
                .wrap_err("while integrating an R* remainder through a massless p-integral")?;
        }
        Ok(result)
    }

    fn evaluate(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        scalar: &Atom,
        coefficient: &Atom,
    ) -> Result<Atom> {
        let retained = settings
            .evaluation_order
            .0
            .iter()
            .find_map(|method| match method {
                EvaluationMethod::Forcer(options) => Some(options.master_mode),
                _ => None,
            })
            == Some(RustRedMasterMode::Retain);
        if retained {
            return Ok(coefficient * self.closure.evaluate(vakint, settings, scalar.as_view())?);
        }
        let epsilon = vakint_parse!(&settings.epsilon_symbol)
            .map_err(|error| eyre!(error))?
            .get_symbol()
            .ok_or_else(|| eyre!("Vakint epsilon must be a symbol"))?;
        let series = coefficient
            .series(epsilon, 0, SeriesDepth::relative(1))
            .wrap_err("could not determine the Laurent order of an R* insertion coefficient")?;
        if series.is_zero() {
            return Ok(Atom::Zero);
        }
        let exponent = series.get_trailing_exponent();
        ensure!(
            exponent.is_integer(),
            "an R* insertion coefficient needs integer epsilon powers"
        );
        let valuation = exponent
            .numerator()
            .to_i64()
            .ok_or_else(|| eyre!("the R* insertion's epsilon power is too large"))?;
        // A pole in the already integrated insertion promotes higher orders of
        // the remaining loop integral. Keep tensor spectators outside Vakint,
        // but supply precisely the extra depth their Laurent coefficient needs.
        let mut closure_settings = settings.clone();
        closure_settings.number_of_terms_in_epsilon_expansion = settings
            .number_of_terms_in_epsilon_expansion
            .checked_sub(valuation.min(0))
            .ok_or_else(|| eyre!("the R* insertion requires an overflowing Laurent depth"))?;
        let integrated = self
            .closure
            .evaluate(vakint, &closure_settings, scalar.as_view())?;
        let target = settings.number_of_terms_in_epsilon_expansion - self.loops as i64 - 1;
        Ok((coefficient * integrated)
            .series(epsilon, 0, target)
            .wrap_err("could not truncate the complete R* insertion and integral product")?
            .to_atom())
    }

    fn new(lines: &[VacuumLine], massive: &VacuumLine, numerator: &Atom) -> Result<Self> {
        let mut heads = BTreeSet::new();
        for line in lines {
            line.momentum.visitor(&mut |part| {
                if let AtomView::Fun(function) = part
                    && function.get_symbol() == S.k
                    && function.get_nargs() == 1
                {
                    heads.insert(part.to_owned());
                    return false;
                }
                true
            });
        }
        let heads = heads.into_iter().collect::<Vec<_>>();
        for line in lines {
            for (head, coefficient) in line.momentum.coefficient_list::<i16>(&heads) {
                ensure!(
                    heads.contains(&head) && i64::try_from(coefficient.as_view()).is_ok(),
                    "R* vacuum routing requires integer linear loop momenta: {}",
                    line.momentum
                );
            }
        }
        let (pivot, sign) = heads
            .iter()
            .find_map(|head| {
                let coefficient = massive.momentum.coefficient(head.as_view());
                let sign = i64::try_from(coefficient.as_view()).ok()?;
                (sign == 1 || sign == -1).then_some((head, sign))
            })
            .ok_or_else(|| {
                eyre!(
                    "the marked R* momentum needs a unit loop coefficient for a measure-preserving cut: {}",
                    massive.momentum
                )
            })?;
        // The closing momentum is integrated, whereas every pre-existing p(i)
        // remains an external spectator, including indexed and compact forms.
        let mut external_ids = BTreeSet::new();
        numerator.visitor(&mut |part| {
            if let AtomView::Fun(function) = part
                && function.get_symbol() == S.p
                && let Some(index) = function.iter().next()
            {
                external_ids.insert(index.to_owned());
                return false;
            }
            true
        });
        let hard_index = (1_i64..)
            .find(|index| !external_ids.contains(&Atom::num(*index)))
            .expect("a finite numerator leaves unused momentum identifiers");
        let hard = function!(S.p, hard_index);
        let mut routing = heads
            .iter()
            .filter(|head| *head != pivot)
            .enumerate()
            .map(|(index, head)| Replacement::new(head.to_pattern(), function!(S.k, index + 1)))
            .collect::<Vec<_>>();
        let rest = (&massive.momentum - Atom::num(sign) * pivot)
            .expand()
            .replace_multiple(&routing);
        routing.push(Replacement::new(
            pivot.to_pattern(),
            Atom::num(sign) * (&hard - rest),
        ));
        let mut propagators = Atom::one();
        let mut radial = Atom::one();
        let mut has_loop = false;
        let mut has_hard = false;
        for line in lines.iter().filter(|line| !std::ptr::eq(*line, massive)) {
            let momentum = line.momentum.replace_multiple(&routing).expand();
            ensure!(
                !momentum.is_zero(),
                "R* vacuum routing produced a zero-momentum propagator"
            );
            has_loop |= momentum.contains_symbol(S.k);
            has_hard |= momentum.contains_symbol(S.p);
            propagators *= function!(S.prop, &line.id, &line.ends, &momentum, 0, line.power);
            // Only used if the cut is a tree, where every momentum is ±p.
            // Multi-loop cuts leave tree contraction to Vakint's frontend.
            if !momentum.contains_symbol(S.k) {
                let coefficient = momentum.coefficient(hard.as_view());
                let square = coefficient.pow(2) * function!(S.dot, &hard, &hard);
                radial *= square.pow(-line.power);
            }
        }
        Ok(Self {
            closure: MassiveClosure {
                momentum: hard,
                mass_squared: massive.mass_squared.clone(),
                power: massive.power,
            },
            routing,
            integrand: if has_loop {
                function!(S.topo, propagators)
            } else {
                radial
            },
            scaleless: has_loop && !has_hard,
            loops: heads.len(),
        })
    }
}

#[cfg(test)]
mod tests {
    use std::sync::LazyLock;

    use super::*;
    use symbolica::{coefficient::CoefficientView, prelude::RealLike};
    use vakint::{EvaluationOrder, ForcerOptions, LoopNormalizationFactor, VakintTerm};

    static VAKINT: LazyLock<Vakint> = LazyLock::new(|| Vakint::new().unwrap());

    fn parse(input: &str) -> Atom {
        vakint_parse!(input).unwrap()
    }

    fn settings() -> VakintSettings {
        VakintSettings {
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
            number_of_terms_in_epsilon_expansion: 4,
            n_cores: 1,
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::Forcer(ForcerOptions {
                master_mode: RustRedMasterMode::Symbolic,
            })]),
            ..Default::default()
        }
    }

    fn integrate(input: &str) -> Atom {
        integrate_single_mass(
            &VAKINT,
            &settings(),
            &VakintExpression::try_from(parse(input)).unwrap(),
        )
        .unwrap()
    }

    fn assert_equal(actual: Atom, expected: Atom) {
        assert!(
            (&actual - &expected).together().is_zero(),
            "actual: {actual}\nexpected: {expected}"
        );
    }

    fn assert_numerically_equal(actual: Atom, expected: Atom) {
        let difference = (actual - expected).expand();
        let mut spectators = BTreeSet::new();
        difference.visitor(&mut |part| {
            if matches!(part, AtomView::Var(_) | AtomView::Fun(_)) {
                spectators.insert(part.to_owned());
                return false;
            }
            true
        });
        for (_, coefficient) in
            difference.coefficient_list::<i16>(&spectators.into_iter().collect::<Vec<_>>())
        {
            if coefficient.is_zero() {
                continue;
            }
            let AtomView::Num(number) = coefficient.as_view() else {
                panic!("non-numeric residual coefficient: {coefficient}")
            };
            let CoefficientView::Float(real, imaginary) = number.get_coeff_view() else {
                panic!("nonzero exact residual coefficient: {coefficient}")
            };
            // Default precision is thirty decimal digits. Summing before or
            // after numerical materialization can differ in the last digits.
            assert!(real.to_float().to_f64().abs() < 1e-24, "{difference}");
            assert!(imaginary.to_float().to_f64().abs() < 1e-24, "{difference}");
        }
    }

    // The former integration path, retained only as an independent ordering
    // oracle: integrate each loop monomial before multiplying its spectator,
    // explicitly promoting the integral depth for spectator poles.
    fn split_reference(input: &VakintExpression, settings: &VakintSettings) -> Result<Atom> {
        let mut tensor_settings = settings.clone();
        tensor_settings.tensor_reduction_method = TensorReductionMethod::FeynKit;
        tensor_settings.use_dot_product_notation = true;
        let mut result = Atom::Zero;
        for source in &input.0 {
            if source.numerator.is_zero() {
                continue;
            }
            if source.integral == Atom::one() {
                result += &source.numerator;
                continue;
            }
            let lines = VacuumLine::parse(&source.integral)?;
            let Some(massive) = lines
                .iter()
                .find(|line| !line.mass_squared.is_zero() && line.power > 0)
            else {
                continue;
            };
            let cut = VacuumCut::new(&lines, massive, &source.numerator)?;
            if cut.scaleless {
                continue;
            }
            let mut term = source.clone();
            term.tensor_reduce(&VAKINT, &tensor_settings)?;
            term.numerator = Integrated::restore_lorentz_dimensions(&term.numerator);
            result += cut.evaluate_split(&VAKINT, settings, &term.numerator)?;
        }
        Ok(result.expand())
    }

    #[test]
    fn complete_numerator_matches_split_integration_with_tensor_spectators_and_poles() {
        let topology = parse(
            "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-k(2),0,1)*prop(3,edge(2,1),k(2),1,2))",
        );
        let polynomial = parse(
            "1/ε+dot(k(1),k(1))/(3-2*ε)+dot(k(1),k(2))/ε^2+dot(p(1),p(1))*dot(k(2),k(2))/(ε*(1-ε))",
        );
        let spectators = [
            Atom::one(),
            parse("g(71,72)"),
            parse("spenso::gamma(spenso::bis(4,a),spenso::bis(4,b),spenso::mink(4,mu))"),
            function!(idenso::color::CS.cas, 2, parse("spenso::coad(8)")).pow(2),
        ];
        for mode in [
            RustRedMasterMode::Symbolic,
            RustRedMasterMode::Retain,
            RustRedMasterMode::Numerical,
        ] {
            let settings = VakintSettings {
                number_of_terms_in_epsilon_expansion: 3,
                evaluation_order: EvaluationOrder(vec![EvaluationMethod::Forcer(ForcerOptions {
                    master_mode: mode,
                })]),
                ..settings()
            };
            for spectator in &spectators {
                let input =
                    VakintExpression::try_from(spectator * &polynomial * &topology).unwrap();
                let actual = integrate_single_mass(&VAKINT, &settings, &input).unwrap();
                assert!(!actual.contains_symbol(symbolica::atom::Symbol::DERIVATIVE));
                let expected = split_reference(&input, &settings).unwrap();
                if mode == RustRedMasterMode::Numerical {
                    assert_numerically_equal(actual, expected);
                } else {
                    assert_equal(actual, expected);
                }
            }
        }
        let mut cancelled = VakintExpression::try_from(topology).unwrap();
        cancelled.0[0].numerator = &polynomial - &polynomial;
        assert!(
            integrate_single_mass(&VAKINT, &settings(), &cancelled)
                .unwrap()
                .is_zero()
        );
    }

    #[test]
    fn cut_momentum_avoids_indexed_and_compact_external_identifiers() {
        let topology = parse("topo(prop(1,edge(1,1),k(1),1,2))");
        let lines = VacuumLine::parse(&topology).unwrap();
        let numerator = parse(
            "p(1,71)*p(2,spenso::mink(4,72))*spenso::dot(p(3,spenso::mink(4)),p(3,spenso::mink(4)))",
        );
        let cut = VacuumCut::new(&lines, &lines[0], &numerator).unwrap();
        assert_eq!(cut.closure.momentum, parse("p(4)"));
    }

    #[test]
    fn strict_numerator_verification_keeps_scalar_spectators_outside_matching() {
        let input = VakintExpression::try_from(parse(
            "c*(1/ε+dot(k(1),k(2))/(3-2*ε))*topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-k(2),0,1)*prop(3,edge(2,1),k(2),1,2))",
        ))
        .unwrap();
        let expected = integrate_single_mass(&VAKINT, &settings(), &input).unwrap();
        let strict = VakintSettings {
            verify_numerator_identification: true,
            ..settings()
        };
        assert_equal(
            integrate_single_mass(&VAKINT, &strict, &input).unwrap(),
            expected,
        );
    }

    #[test]
    fn one_loop_cut_keeps_the_tree_denominator_and_marked_power() {
        let actual = integrate("topo(prop(1,edge(1,2),k(7),mass2,2)*prop(2,edge(2,1),k(7),0,1))");
        let expected = MassiveClosure {
            momentum: parse("p(1)"),
            mass_squared: parse("mass2"),
            power: 2,
        }
        .evaluate(&VAKINT, &settings(), parse("1/dot(p(1),p(1))").as_view())
        .unwrap();
        assert_equal(actual, expected);
    }

    #[test]
    fn tensor_prepass_retains_external_coefficients_even_with_the_cut_momentum_id() {
        let actual = integrate("spectator*dot(k(3),p(1))^2*topo(prop(1,edge(1,1),k(3),1,3))");
        let expected = MassiveClosure {
            momentum: parse("p(1)"),
            mass_squared: Atom::one(),
            power: 3,
        }
        .evaluate(&VAKINT, &settings(), parse("dot(p(1),p(1))").as_view())
        .unwrap()
            * parse("spectator*dot(p(1),p(1))/(4-2*ε)");
        assert_equal(
            actual,
            expected
                .series(parse("ε").get_symbol().unwrap(), 0, 2)
                .unwrap()
                .to_atom(),
        );
    }

    #[test]
    fn insertion_poles_promote_the_required_higher_orders_of_the_closure() {
        for (terms, numerator, expected) in [
            (2, "1/ε", "1/ε^2+𝜋^2/12"),
            (1, "1/ε^2", "1/ε^3+𝜋^2/(12*ε)"),
            (1, "1/(ε^2*(4-2*ε))", "1/(4*ε^3)+1/(8*ε^2)+(1/16+𝜋^2/48)/ε"),
            (
                2,
                "dot(p(1),p(1))/(ε*(1-ε))",
                "dot(p(1),p(1))*(1/ε^2+1/ε+1+𝜋^2/12)",
            ),
        ] {
            let settings = VakintSettings {
                number_of_terms_in_epsilon_expansion: terms,
                ..settings()
            };
            let input = VakintExpression::try_from(
                parse(numerator) * parse("topo(prop(1,edge(1,1),k(1),1,2))"),
            )
            .unwrap();
            let actual = integrate_single_mass(&VAKINT, &settings, &input)
                .unwrap()
                .replace(parse("mursq").to_pattern())
                .with(Atom::one());
            assert_equal(actual, parse(expected));
        }
    }

    #[test]
    fn massless_remainders_vanish_but_local_insertions_survive() {
        let mut input =
            VakintExpression::try_from(parse("local_mass*topo(prop(1,edge(1,1),k(1),0,2))"))
                .unwrap();
        input.0.push(VakintTerm {
            integral: Atom::one(),
            numerator: parse("local_mass/ε"),
            vectors: Vec::new(),
        });
        assert_equal(
            integrate_single_mass(&VAKINT, &settings(), &input).unwrap(),
            parse("local_mass/ε"),
        );
    }

    #[test]
    fn a_marked_linear_combination_routes_with_unit_jacobian() {
        let canonical = parse(
            "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-k(2),0,1)*prop(3,edge(2,1),k(2),1,2))",
        );
        let rerouted = canonical.replace_multiple(&[
            Replacement::new(parse("k(1)").to_pattern(), parse("k(5)")),
            Replacement::new(parse("k(2)").to_pattern(), parse("k(5)-k(8)")),
        ]);
        let evaluate = |input| {
            integrate_single_mass(
                &VAKINT,
                &settings(),
                &VakintExpression::try_from(input).unwrap(),
            )
            .unwrap()
        };
        let actual = evaluate(canonical);
        assert_equal(evaluate(rerouted), actual.clone());
        let epsilon = parse("ε").get_symbol().unwrap();
        let series = actual.series(epsilon, 0, 0).unwrap();
        let (_, leading_pole) = series.terms().find(|(power, _)| *power == -2).unwrap();
        assert_equal(leading_pole.clone(), parse("1/2"));
    }

    #[test]
    fn a_pure_hard_bridge_is_contracted_by_the_p_integral_frontend() {
        let actual = integrate(
            "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-k(2),0,1)*prop(3,edge(2,3),k(2),1,2)*prop(4,edge(3,1),k(2),0,1))",
        );
        let expected = MassiveClosure {
            momentum: parse("p(1)"),
            mass_squared: Atom::one(),
            power: 2,
        }
        .evaluate(
            &VAKINT,
            &settings(),
            parse(
                "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,1))/dot(p(1),p(1))",
            )
            .as_view(),
        )
        .unwrap();
        assert_equal(actual, expected);
    }

    #[test]
    fn two_massive_propagators_are_rejected() {
        let input = VakintExpression::try_from(parse(
            "topo(prop(1,edge(1,2),k(1),1,1)*prop(2,edge(2,1),k(1),1,1))",
        ))
        .unwrap();
        assert!(integrate_single_mass(&VAKINT, &settings(), &input).is_err());
    }
}
