//! Checked numerical coefficients with symbolic spectators retained.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::CoefficientView,
    domains::float::{Complex, Float},
    id::Replacement,
};

use crate::{
    Vakint, VakintError, VakintSettings, get_integer_from_atom, hpl,
    inspection::{PreparedExpression, constants},
    rustred_evaluation::masters::{MasterCatalog, combine_exact},
    utils::vakint_macros::vk_symbol,
};

use super::{attach_origins, evaluate_coefficient, legacy};

fn failure(message: impl Into<String>) -> VakintError {
    VakintError::NumericalPrecision(message.into())
}

pub(crate) fn snapshot_reference(part: AtomView) -> bool {
    matches!(part, AtomView::Fun(function) if ["RustRedMaster", "RustRedCoefficient", "RustRedFiniteLiteral", "RustRedConstant"]
        .iter().any(|name| function.get_symbol() == vk_symbol!(name)))
}

fn private_reference(part: AtomView) -> bool {
    snapshot_reference(part) || legacy::is_source(part)
}

/// Classification is structural: no trial numerical evaluations or retries.
struct NumericalParts<'a> {
    parameters: &'a HashMap<Atom, Complex<Float>>,
    catalog: Option<&'a MasterCatalog>,
    numeric: BTreeMap<Atom, bool>,
}

impl NumericalParts<'_> {
    fn is_numeric(&mut self, part: AtomView) -> Result<bool, VakintError> {
        let owned = part.to_owned();
        if legacy::is_source(part) {
            return Ok(true);
        }
        if self.parameters.contains_key(&owned) {
            return Ok(true);
        }
        if let Some(value) = self.numeric.get(&owned) {
            return Ok(*value);
        }
        let value = match part {
            AtomView::Num(_) => true,
            AtomView::Var(variable) => {
                if constants::builtin_definition(part)?.is_some()
                    || variable.get_symbol().get_evaluation_info().is_some()
                {
                    true
                } else if let Some(catalog) = self.catalog {
                    // This resolves snapshot-bound named periods without
                    // evaluating their numerical payloads. Unknown variables
                    // remain unchanged and are therefore spectators.
                    catalog.prepare_evaluation(owned.clone())?.expression != owned
                } else {
                    false
                }
            }
            AtomView::Fun(function) => {
                if function.get_symbol() == vk_symbol!("RustRedCoefficient")
                    || function.get_symbol() == vk_symbol!("RustRedConstant")
                {
                    // resolve_numerical validated its snapshot and epsilon
                    // order before classification. Arguments are identity tags.
                    true
                } else if self
                    .catalog
                    .is_some_and(|catalog| catalog.constants.contains_key(&owned))
                    || hpl::projection_kernel(part)?.is_some()
                {
                    true
                } else {
                    let symbol = function.get_symbol();
                    let info = symbol.get_evaluation_info();
                    let tags = info.map_or(0, |info| info.get_tag_count());
                    let mut constant_arguments = true;
                    for argument in function.iter().skip(tags) {
                        constant_arguments &= self.is_numeric(argument)?;
                    }
                    constant_arguments
                        && (symbol.is_builtin()
                            || info.is_some()
                            || constants::builtin_definition(part)?.is_some())
                }
            }
            AtomView::Pow(power) => {
                self.is_numeric(power.get_base())? && self.is_numeric(power.get_exp())?
            }
            AtomView::Add(sum) => {
                let mut numeric = true;
                for term in sum.iter() {
                    numeric &= self.is_numeric(term)?;
                }
                numeric
            }
            AtomView::Mul(product) => {
                let mut numeric = true;
                for factor in product.iter() {
                    numeric &= self.is_numeric(factor)?;
                }
                numeric
            }
        };
        self.numeric.insert(owned, value);
        Ok(value)
    }

    fn spectators(
        &mut self,
        part: AtomView,
        spectators: &mut BTreeSet<Atom>,
    ) -> Result<(), VakintError> {
        if self.is_numeric(part)? {
            return Ok(());
        }
        match part {
            AtomView::Add(sum) => {
                for term in sum.iter() {
                    self.spectators(term, spectators)?;
                }
            }
            AtomView::Mul(product) => {
                for factor in product.iter() {
                    self.spectators(factor, spectators)?;
                }
            }
            AtomView::Pow(power)
                if get_integer_from_atom(power.get_exp())
                    .is_some_and(|exponent| i32::try_from(exponent).is_ok())
                    && matches!(power.get_base(), AtomView::Var(_) | AtomView::Fun(_)) =>
            {
                self.spectators(power.get_base(), spectators)?;
            }
            _ => {
                let mut references = BTreeSet::new();
                part.visitor(&mut |child| {
                    if private_reference(child) {
                        references.insert(child.to_owned());
                        false
                    } else {
                        true
                    }
                });
                if !references.is_empty() {
                    return Err(failure(format!(
                        "cannot numerically materialize master coefficients inside the unevaluated spectator {part}; unresolved references: {references:?}"
                    )));
                }
                spectators.insert(part.to_owned());
            }
        }
        Ok(())
    }
}

fn prepare_coefficient(
    expression: Atom,
    catalog: Option<&MasterCatalog>,
    parameters: &HashMap<Atom, Complex<Float>>,
) -> Result<PreparedExpression, VakintError> {
    if let Some(catalog) = catalog {
        return catalog.prepare_evaluation(expression);
    }
    // FORM callers without a RustRed snapshot still use the shared exact
    // definitions. Finite declared parameters remain numerical evaluator inputs.
    let mut expression = combine_exact(Vakint::normalize_gamma_functions(expression.as_view()));
    for _ in 0..128 {
        let mut replacements = BTreeMap::new();
        let mut error = None;
        expression.visitor(&mut |part| {
            if error.is_some()
                || parameters.contains_key(&part.to_owned())
                || legacy::is_source(part)
            {
                return false;
            }
            match constants::builtin_definition(part) {
                Ok(Some(definition)) => {
                    if definition.as_view() != part {
                        replacements.insert(part.to_owned(), definition);
                    }
                    false
                }
                Ok(None) => true,
                Err(value) => {
                    error = Some(value);
                    false
                }
            }
        });
        if let Some(error) = error {
            return Err(error);
        }
        if replacements.is_empty() {
            return Ok(PreparedExpression {
                expression,
                sources: BTreeMap::new(),
                origins: Vec::new(),
            });
        }
        expression = combine_exact(
            expression.replace_multiple(
                replacements
                    .into_iter()
                    .map(|(source, target)| Replacement::new(source, target))
                    .collect::<Vec<_>>(),
            ),
        );
    }
    Err(failure("constant definitions did not reach a fixed point"))
}

fn epsilon_power(monomial: AtomView, epsilon: Symbol) -> Option<i64> {
    match monomial {
        AtomView::Var(variable) if variable.get_symbol() == epsilon => Some(1),
        AtomView::Pow(power) if power.get_base() == Atom::var(epsilon).as_view() => {
            get_integer_from_atom(power.get_exp())
        }
        AtomView::Mul(product) => product.iter().try_fold(0_i64, |sum, factor| {
            sum.checked_add(epsilon_power(factor, epsilon)?)
        }),
        _ if !monomial.contains_symbol(epsilon) => Some(0),
        _ => None,
    }
}

/// Evaluate complete constant coefficients, retaining unbound symbolic factors.
pub(crate) fn materialize(
    expression: AtomView,
    settings: &VakintSettings,
    parameters: &HashMap<Atom, Complex<Float>>,
    catalog: Option<&MasterCatalog>,
) -> Result<Atom, VakintError> {
    for identity in parameters.keys() {
        crate::inspection::validate_parameter(identity.as_view(), catalog)?;
    }
    if settings.run_time_decimal_precision == 0 {
        return Err(failure("requested decimal precision must be positive"));
    }
    // Raw approximate literals must keep independent source identities before
    // expansion can accidentally cancel their equal midpoints.
    let expression = legacy::stage_table(
        Atom::var(vk_symbol!("partial_numerical_input")).as_view(),
        expression,
        vk_symbol!(&settings.epsilon_symbol),
        "Input",
    )?;
    let expression = combine_exact(Vakint::normalize_gamma_functions(expression.as_view()));
    let expression = if let Some(catalog) = catalog {
        catalog.resolve_numerical(expression, settings)?
    } else {
        let mut references = BTreeSet::new();
        expression.visitor(&mut |part| {
            if snapshot_reference(part) {
                references.insert(part.to_owned());
                false
            } else {
                !legacy::is_source(part)
            }
        });
        if !references.is_empty() {
            return Err(failure(format!(
                "cannot materialize master references without their data snapshot: {references:?}"
            )));
        }
        expression
    };
    if expression.is_zero() {
        return Ok(Atom::Zero);
    }
    let mut classification = NumericalParts {
        parameters,
        catalog,
        numeric: BTreeMap::new(),
    };
    let mut spectators = BTreeSet::new();
    classification.spectators(expression.as_view(), &mut spectators)?;

    // Function indeterminates also represent composite negative/fractional
    // powers safely. Every temporary call is absent from the input, and all
    // maps are request-local; only one global function symbol is registered.
    let mut occupied = BTreeSet::new();
    expression.visitor(&mut |part| {
        occupied.insert(part.to_owned());
        !private_reference(part)
    });
    let mut replacements = Vec::new();
    let mut restore = Vec::new();
    let mut variables = Vec::new();
    let mut index = 0_i64;
    for spectator in spectators {
        let temporary = loop {
            let candidate = vk_symbol!("RustRedMaterializeSpectator").call(index);
            index = index
                .checked_add(1)
                .ok_or_else(|| failure("too many symbolic spectators"))?;
            if !occupied.contains(&candidate) {
                break candidate;
            }
        };
        replacements.push(Replacement::new(spectator.clone(), temporary.clone()));
        restore.push(Replacement::new(temporary.clone(), spectator));
        variables.push(temporary);
    }
    let grouped = expression.replace_multiple(&replacements);
    let mut result = Atom::Zero;
    for (monomial, coefficient) in grouped.coefficient_list::<i32>(&variables) {
        let monomial = monomial.replace_multiple(&restore);
        let prepared = prepare_coefficient(coefficient, catalog, parameters)?;
        // A resolved rational is already an exact numerical coefficient.
        // Retain it across partial evaluation boundaries (for example nested
        // forest integrations); only final numerical evaluation needs a Float.
        // Staged Float literals and finite sources cannot enter this branch.
        if let AtomView::Num(number) = prepared.expression.as_view()
            && matches!(
                number.get_coeff_view(),
                CoefficientView::Natural(..) | CoefficientView::Large(..)
            )
        {
            result += prepared.expression * monomial;
            continue;
        }
        let (value, mut warnings) = evaluate_coefficient(
            prepared.expression.as_view(),
            &prepared.sources,
            parameters,
            settings,
        )?;
        attach_origins(
            &mut warnings,
            &prepared.origins,
            epsilon_power(monomial.as_view(), vk_symbol!(&settings.epsilon_symbol)),
            catalog.map_or("built-in definitions and declared parameters", |catalog| {
                catalog.version.as_str()
            }),
        );
        for warning in warnings {
            log::warn!("{warning}; symbolic factor {monomial}");
        }
        let coefficient = value
            .to_atom()
            .map_err(|error| failure(format!("coefficient multiplying {monomial}: {error}")))?;
        result += coefficient * monomial;
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::vakint_macros::vk_parse;

    fn settings() -> VakintSettings {
        Vakint::initialize_vakint_symbols();
        VakintSettings {
            run_time_decimal_precision: 40,
            ..Default::default()
        }
    }

    #[test]
    fn tensor_spectators_and_laurent_orders_survive_materialization() {
        let settings = settings();
        let source = vk_parse!("(z3+Zeta(3))*ep^-2*p(1,mu)").unwrap();
        let result = materialize(source.as_view(), &settings, &HashMap::new(), None).unwrap();
        let terms = result
            .coefficient_list::<i32>(&[vk_parse!("ep").unwrap(), vk_parse!("p(1,mu)").unwrap()]);
        assert_eq!(terms.len(), 1);
        assert_eq!(terms[0].0, vk_parse!("ep^-2*p(1,mu)").unwrap());
        assert!(matches!(terms[0].1.as_view(), AtomView::Num(_)));
        assert!(!result.contains_symbol(vk_symbol!("z3")));
        assert!(!result.contains_symbol(vk_symbol!("Zeta")));
    }

    #[test]
    fn complete_coefficients_cancel_before_approximation() {
        let settings = settings();
        let source = vk_parse!("(10^80*Zeta(3)-10^80*z3+2)*x").unwrap();
        let result = materialize(source.as_view(), &settings, &HashMap::new(), None).unwrap();
        let terms = result.coefficient_list::<i32>(&[vk_parse!("x").unwrap()]);
        assert_eq!(terms.len(), 1);
        let value = terms[0]
            .1
            .evaluate_with_prec(&ahash::HashMap::<Atom, Complex<Float>>::default(), 128)
            .unwrap();
        assert_eq!(value.re, Float::with_val(128, 2));
        assert!(value.im.as_raw().is_zero());
    }

    #[test]
    fn public_partial_evaluation_distinguishes_exact_and_approximate_polynomial_zero() {
        let settings = settings();
        let evaluate = |expression: &str| {
            Vakint::partial_numerical_evaluation(
                &settings,
                vk_parse!(expression).unwrap().as_view(),
                &Default::default(),
                &Default::default(),
                None,
            )
        };
        assert!(
            evaluate("((z3+1)^2-Zeta(3)^2-2*Zeta(3)-1)*x")
                .unwrap()
                .is_zero()
        );
        // Equal central values do not assert an exact relation among these
        // three independently supplied finite estimates.
        let error = evaluate("((1.25`30+z3)^2-z3^2-1.5625`30-2.5`30*z3)*x").unwrap_err();
        assert!(matches!(error, VakintError::NumericalPrecision(_)));
        assert!(error.to_string().contains("uncertain"), "{error}");
    }

    #[test]
    fn composite_denominators_fractional_powers_and_temporary_name_collisions_are_preserved() {
        let settings = settings();
        let source = vk_parse!("z3/(1+x)+z3*y^(1/2)+z3*RustRedMaterializeSpectator(0)").unwrap();
        let result = materialize(source.as_view(), &settings, &HashMap::new(), None).unwrap();
        assert!(!result.contains_symbol(vk_symbol!("z3")));
        for factor in ["(1+x)^-1", "y^(1/2)", "RustRedMaterializeSpectator(0)"] {
            assert!(
                result
                    .pattern_match(&vk_parse!(factor).unwrap().to_pattern(), None, None)
                    .next()
                    .is_some(),
                "{factor}: {result}"
            );
        }
    }

    #[test]
    fn uncertain_zero_is_not_silently_exported() {
        let settings = settings();
        let source = vk_parse!("(a-b)*x").unwrap();
        let parameters = ["a", "b"]
            .map(|name| {
                (
                    Atom::var(vk_symbol!(name)),
                    Complex::new(Float::with_val(128, 1), Float::new(128)),
                )
            })
            .into_iter()
            .collect();
        assert!(matches!(
            materialize(source.as_view(), &settings, &parameters, None),
            Err(VakintError::NumericalPrecision(_))
        ));
        let cancelled = vk_parse!("(a-a)*x").unwrap();
        assert!(
            materialize(cancelled.as_view(), &settings, &parameters, None)
                .unwrap()
                .is_zero()
        );
    }

    #[test]
    fn unbound_master_references_are_not_treated_as_spectators() {
        let settings = settings();
        let source = vk_parse!("x*RustRedMaster(family,1)").unwrap();
        assert!(materialize(source.as_view(), &settings, &HashMap::new(), None).is_err());
    }
}
