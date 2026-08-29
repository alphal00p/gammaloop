use ::rustred::algebra::CoefficientContext;
use ::rustred::family::ScalarProductCoordinate;
use ::rustred::tensor::{
    TensorGuardOrigin, TensorHeads, TensorLane, TensorLimits, TensorMomenta, TensorProjection,
    TensorService,
};
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::coefficient::CoefficientView;
use symbolica::function;
use symbolica::id::Pattern;

use super::{RustRedOptions, TensorReductionError};
use crate::symbols::S;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{
    Topology, Vakint, VakintExpression, VakintSettings, VakintTerm, get_integer_from_atom,
    get_prop_with_id,
};

mod momentum_labels;
mod vacuum_family;

use momentum_labels::MomentumLabels;
use vacuum_family::{PhysicalVacuumPropagator, VacuumFamilyCache, VacuumFamilySpec};

pub(super) fn reduce(
    expression: &mut VakintExpression,
    vakint: &Vakint,
    settings: &VakintSettings,
    _options: RustRedOptions,
) -> Result<(), TensorReductionError> {
    let heads = TensorHeads::try_new(S.k, S.p, S.g, S.dot).map_err(|source| {
        TensorReductionError::RustRedTensor {
            term: 0,
            source: source.into(),
        }
    })?;

    // Keep the operation transactional: topology canonicalization and native
    // projection become visible only if every input term is admitted.
    let mut projected = expression.clone();
    let mut family_cache = VacuumFamilyCache::default();
    let momentum_labels = MomentumLabels::new();
    for (term_index, term) in projected.0.iter_mut().enumerate() {
        project_term(
            term,
            term_index,
            vakint,
            settings,
            heads,
            momentum_labels,
            &mut family_cache,
        )?;
    }
    *expression = projected;
    Ok(())
}

fn project_term(
    term: &mut VakintTerm,
    term_index: usize,
    vakint: &Vakint,
    settings: &VakintSettings,
    heads: TensorHeads,
    momentum_labels: MomentumLabels,
    family_cache: &mut VacuumFamilyCache,
) -> Result<(), TensorReductionError> {
    // Reject caller-spelled private heads before Vakint's canonicalizer can
    // interpret the surrounding expression, then repeat the check on the
    // normalized boundary value immediately before adapter tagging.
    momentum_labels.reject_reserved_input(term_index, term.numerator.as_view())?;
    let replacement_rules = vakint
        .topologies
        .match_topologies_to_user_input(term.integral.as_view(), settings.allow_unknown_integrals)?
        .ok_or_else(|| TensorReductionError::RustRedUnrecognizedIntegral {
            term: term_index,
            integral: term.integral.to_string(),
        })?;

    if matches!(replacement_rules.canonical_topology, Topology::Unknown(_)) {
        return Err(TensorReductionError::RustRedUnknownTopology { term: term_index });
    }
    let matched = replacement_rules.canonical_topology.get_integral();
    let loop_count = matched.n_loops;
    let propagator_count = matched.n_props;
    validate_explicit_input_routing(term, term_index, loop_count, &replacement_rules)?;

    // This is the existing Vakint matcher/routing engine. In particular,
    // numerator substitutions are simultaneous; the adapter never invents a
    // second routing or dispatches on a topology name.
    term.canonicalize(settings, &replacement_rules, false)?;
    let canonical_loop_ids = canonical_loop_ids(term_index, loop_count)?;
    let loop_momenta = canonical_loop_ids
        .iter()
        .map(|label| function!(S.k, label.clone()))
        .collect::<Vec<_>>();
    let mut physical_propagators = Vec::with_capacity(propagator_count);
    for propagator_id in 1..=propagator_count {
        // Contracted registry topologies retain the parent topology's
        // `n_props`; the canonical expression is authoritative for which
        // physical rows survive the pinch.
        let Some(propagator) = get_prop_with_id(term.integral.as_view(), propagator_id) else {
            continue;
        };
        let momentum = propagator
            .get(&vk_symbol!("q_"))
            .expect("Vakint's canonical propagator matcher always captures momentum");
        let momentum_row =
            canonical_momentum_row(term_index, propagator_id, momentum, &loop_momenta)?;
        let mass = propagator
            .get(&vk_symbol!("mUVsq_"))
            .expect("Vakint's canonical propagator matcher always captures mass");
        let power = propagator
            .get(&vk_symbol!("pow_"))
            .expect("Vakint's canonical propagator matcher always captures power");
        let power = get_integer_from_atom(power.as_view()).ok_or_else(|| {
            TensorReductionError::RustRedUnsupportedPower {
                term: term_index,
                propagator: propagator_id,
                power: power.to_string(),
            }
        })?;
        physical_propagators.push(PhysicalVacuumPropagator::new(
            propagator_id,
            momentum_row,
            mass.clone(),
            power,
        ));
    }

    let bridge = family_cache.bind(
        term_index,
        VacuumFamilySpec::new(loop_count, physical_propagators),
    )?;
    let retained_denominator_count = bridge
        .physical_denominator_ids()
        .len()
        .checked_add(bridge.appended_coordinate_ordinals().len())
        .ok_or_else(|| TensorReductionError::RustRedUnsupportedOutput {
            term: term_index,
            detail: "retained vacuum denominator count overflowed".to_owned(),
        })?;
    if retained_denominator_count != bridge.presentation().family().denominator_count() {
        return Err(TensorReductionError::RustRedUnsupportedOutput {
            term: term_index,
            detail: "vacuum-family role witness does not cover every denominator".to_owned(),
        });
    }
    // Normalize indexed contractions before crossing the native boundary.
    // The private labels retain Vakint's loop/external type even when both use
    // the same numeric ID, and collect external vectors that occur only in a
    // pre-existing scalar `dot`.
    let dot_numerator = Vakint::convert_to_dot_notation(term.numerator.as_view());
    momentum_labels.reject_reserved_input(term_index, dot_numerator.as_view())?;
    let (rustred_numerator, external_labels) = momentum_labels.encode(dot_numerator.as_view());
    let loop_labels = momentum_labels.loop_labels(term_index, loop_count)?;
    let service = TensorService::try_new(
        bridge.presentation(),
        TensorLane::Auto,
        heads,
        TensorMomenta::new(loop_labels.clone(), external_labels),
        TensorLimits::default(),
    )
    .map_err(|source| TensorReductionError::RustRedTensor {
        term: term_index,
        source,
    })?;

    let dimension = parse_dimension(settings)?;
    let projection = service
        .project(&rustred_numerator, bridge.base_integral())
        .map_err(|source| TensorReductionError::RustRedTensor {
            term: term_index,
            source,
        })?;
    let projected_numerator = projection_to_vakint(
        term_index,
        &projection,
        bridge.presentation().family().coefficient_context(),
        &dimension,
        bridge.common_mass_squared(),
        &loop_labels,
    )?;
    let projected_numerator = momentum_labels.decode(term_index, projected_numerator.as_view())?;
    term.numerator = if settings.use_dot_product_notation {
        projected_numerator
    } else {
        Vakint::convert_from_dot_notation(projected_numerator.as_view())
    };
    term.vectors = VakintTerm::identify_vectors_in_numerator(term.numerator.as_view())?;
    Ok(())
}

fn validate_explicit_input_routing(
    term: &VakintTerm,
    term_index: usize,
    loop_count: usize,
    replacement_rules: &crate::ReplacementRules,
) -> Result<(), TensorReductionError> {
    // Short topology notation carries the registry's canonical routing by
    // construction. A graph-shaped input, however, must not let the broad
    // topology matcher erase an inadmissible denominator momentum before the
    // RustRed family proof is built.
    if !contains_function_head(term.integral.as_view(), S.prop) {
        return Ok(());
    }
    if loop_count != 1 {
        return Err(TensorReductionError::RustRedExplicitRoutingUnsupported {
            term: term_index,
            loop_count,
        });
    }
    let input_propagator = replacement_rules
        .edge_ids_canonical_to_input_map
        .get(&1)
        .copied()
        .ok_or(TensorReductionError::RustRedMissingPropagator {
            term: term_index,
            propagator: 1,
        })?;
    let propagator = get_prop_with_id(term.integral.as_view(), input_propagator).ok_or(
        TensorReductionError::RustRedMissingPropagator {
            term: term_index,
            propagator: input_propagator,
        },
    )?;
    let momentum = propagator
        .get(&vk_symbol!("q_"))
        .expect("Vakint's explicit propagator matcher always captures momentum");
    if !is_bare_loop_momentum(momentum.as_view()) {
        return Err(TensorReductionError::RustRedUnsupportedMomentum {
            term: term_index,
            propagator: input_propagator,
            momentum: momentum.to_string(),
        });
    }
    Ok(())
}

fn canonical_loop_ids(
    term_index: usize,
    loop_count: usize,
) -> Result<Vec<Atom>, TensorReductionError> {
    (1..=loop_count)
        .map(|label| {
            i64::try_from(label).map(Atom::num).map_err(|_| {
                TensorReductionError::RustRedUnsupportedOutput {
                    term: term_index,
                    detail: format!("canonical loop label {label} does not fit in i64"),
                }
            })
        })
        .collect()
}

fn canonical_momentum_row(
    term_index: usize,
    propagator_id: usize,
    momentum: &Atom,
    loop_momenta: &[Atom],
) -> Result<Vec<i64>, TensorReductionError> {
    let mut row = vec![0; loop_momenta.len()];
    for (monomial, coefficient) in momentum.coefficient_list::<i16>(loop_momenta) {
        if monomial == Atom::num(1) {
            if coefficient.is_zero() {
                continue;
            }
            return Err(TensorReductionError::RustRedUnsupportedMomentum {
                term: term_index,
                propagator: propagator_id,
                momentum: momentum.to_string(),
            });
        }
        let Some(loop_index) = loop_momenta
            .iter()
            .position(|candidate| candidate == &monomial)
        else {
            return Err(TensorReductionError::RustRedUnsupportedMomentum {
                term: term_index,
                propagator: propagator_id,
                momentum: momentum.to_string(),
            });
        };
        row[loop_index] = get_integer_from_atom(coefficient.as_view()).ok_or_else(|| {
            TensorReductionError::RustRedUnsupportedMomentum {
                term: term_index,
                propagator: propagator_id,
                momentum: momentum.to_string(),
            }
        })?;
    }
    if row.iter().all(|coefficient| *coefficient == 0) {
        return Err(TensorReductionError::RustRedUnsupportedMomentum {
            term: term_index,
            propagator: propagator_id,
            momentum: momentum.to_string(),
        });
    }
    Ok(row)
}

fn contains_function_head(atom: AtomView<'_>, head: symbolica::atom::Symbol) -> bool {
    match atom {
        AtomView::Fun(function) => {
            function.get_symbol() == head
                || function
                    .iter()
                    .any(|argument| contains_function_head(argument, head))
        }
        AtomView::Add(addition) => addition
            .iter()
            .any(|argument| contains_function_head(argument, head)),
        AtomView::Mul(product) => product
            .iter()
            .any(|argument| contains_function_head(argument, head)),
        AtomView::Pow(power) => {
            contains_function_head(power.get_base(), head)
                || contains_function_head(power.get_exp(), head)
        }
        AtomView::Num(_) | AtomView::Var(_) => false,
    }
}

fn is_bare_loop_momentum(momentum: AtomView<'_>) -> bool {
    let AtomView::Fun(function) = momentum else {
        return false;
    };
    function.get_symbol() == S.k
        && function.get_nargs() == 1
        && function
            .iter()
            .next()
            .and_then(get_integer_from_atom)
            .is_some()
}

fn is_exact_real_scalar(atom: AtomView<'_>) -> bool {
    match atom {
        AtomView::Var(_) => true,
        AtomView::Num(number) => match number.get_coeff_view() {
            CoefficientView::Natural(_, _, imaginary, _) => imaginary == 0,
            CoefficientView::Large(_, imaginary) => imaginary.is_zero(),
            _ => false,
        },
        _ => false,
    }
}

fn parse_dimension(settings: &VakintSettings) -> Result<Atom, TensorReductionError> {
    let epsilon = vk_parse!(settings.epsilon_symbol.as_str()).map_err(|error| {
        TensorReductionError::RustRedInvalidDimension {
            symbol: settings.epsilon_symbol.clone(),
            detail: error.to_string(),
        }
    })?;
    if !matches!(epsilon.as_view(), AtomView::Var(_)) {
        return Err(TensorReductionError::RustRedInvalidDimension {
            symbol: settings.epsilon_symbol.clone(),
            detail: "expected one Symbolica symbol, not a value or expression".to_owned(),
        });
    }
    // Retain Vakint/FORM's established canonical spelling while representing
    // the exact same dimension d = 4 - 2 epsilon.
    Ok(Atom::num(-1) * (Atom::num(2) * epsilon - Atom::num(4)))
}

fn projection_to_vakint(
    term_index: usize,
    projection: &TensorProjection,
    coefficients: &CoefficientContext,
    dimension: &Atom,
    mass: &Atom,
    loop_labels: &[Atom],
) -> Result<Atom, TensorReductionError> {
    let dimension_parameter = coefficients
        .parameter("d")
        .expect("the adapter registers the dimension parameter")
        .to_expression();
    let mass_parameter = coefficients
        .parameter("m2")
        .expect("the adapter registers the mass parameter")
        .to_expression();
    let substitutions: [(Pattern, Pattern); 2] = [
        (dimension_parameter.to_pattern(), dimension.to_pattern()),
        (mass_parameter.to_pattern(), mass.to_pattern()),
    ];
    for guard in projection.guards() {
        if guard.origin() != TensorGuardOrigin::RankTwoProjectorDimension {
            return Err(TensorReductionError::RustRedUnsupportedOutput {
                term: term_index,
                detail: format!("unmapped tensor guard {:?}", guard.origin()),
            });
        }
        let specialized = substitutions.iter().fold(
            guard.polynomial().to_expression(),
            |atom, (source, target)| atom.replace(source.clone()).with(target.clone()),
        );
        if specialized.is_zero() {
            return Err(TensorReductionError::RustRedSingularDimension {
                term: term_index,
                dimension: dimension.to_string(),
            });
        }
    }

    let mut output_terms = Vec::new();
    output_terms
        .try_reserve_exact(projection.terms().len())
        .map_err(|_| TensorReductionError::RustRedUnsupportedOutput {
            term: term_index,
            detail: "could not reserve projected output terms".to_owned(),
        })?;
    for projected in projection.terms() {
        let coefficient = substitutions.iter().fold(
            projected.coefficient().to_expression(),
            |atom, (source, target)| atom.replace(source.clone()).with(target.clone()),
        );
        let mut factors = Vec::new();
        let factor_count = 3usize
            .checked_add(projected.scalar_products().len())
            .ok_or_else(|| TensorReductionError::RustRedUnsupportedOutput {
                term: term_index,
                detail: "projected factor count overflowed".to_owned(),
            })?;
        factors.try_reserve_exact(factor_count).map_err(|_| {
            TensorReductionError::RustRedUnsupportedOutput {
                term: term_index,
                detail: format!("could not reserve {factor_count} projected factors"),
            }
        })?;
        factors.push(coefficient);
        factors.push(projected.scalar_spectator().clone());
        factors.push(projected.outside_tensor().clone());
        for scalar_product in projected.scalar_products() {
            factors.push(scalar_product_to_vakint(
                term_index,
                *scalar_product,
                loop_labels,
            )?);
        }
        output_terms.push(Atom::mul_many(factors));
    }
    Ok(Atom::add_many(output_terms))
}

fn scalar_product_to_vakint(
    term_index: usize,
    coordinate: ScalarProductCoordinate,
    loop_labels: &[Atom],
) -> Result<Atom, TensorReductionError> {
    match coordinate {
        ScalarProductCoordinate::LoopLoop { left, right } => {
            let left = loop_labels.get(left).ok_or_else(|| {
                TensorReductionError::RustRedUnsupportedOutput {
                    term: term_index,
                    detail: format!("loop scalar-product index {left} is out of range"),
                }
            })?;
            let right = loop_labels.get(right).ok_or_else(|| {
                TensorReductionError::RustRedUnsupportedOutput {
                    term: term_index,
                    detail: format!("loop scalar-product index {right} is out of range"),
                }
            })?;
            Ok(S.dot(left, right))
        }
        ScalarProductCoordinate::LoopExternal { .. } => {
            Err(TensorReductionError::RustRedUnsupportedOutput {
                term: term_index,
                detail: "the single-scale vacuum bridge produced a loop-external scalar product"
                    .to_owned(),
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reserved_input_label_failure_does_not_publish_prior_term_projection() {
        let vakint = Vakint::new().unwrap();
        let settings = VakintSettings {
            form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
            ..VakintSettings::default()
        };
        let valid = VakintExpression::try_from(
            vk_parse!("k(1,mu)*k(1,nu)*topo(I1L(muvsq,1))")
                .unwrap()
                .as_view(),
        )
        .unwrap()
        .0
        .pop()
        .unwrap();
        let attack = VakintExpression::try_from(
            vk_parse!("k(rustred_loop_label(1),mu)*topo(I1L(muvsq,1))")
                .unwrap()
                .as_view(),
        )
        .unwrap()
        .0
        .pop()
        .unwrap();
        let mut expression = VakintExpression(vec![valid, attack]);
        let original = expression.clone();

        let error = reduce(&mut expression, &vakint, &settings, RustRedOptions::new()).unwrap_err();

        assert!(matches!(
            error,
            TensorReductionError::RustRedReservedMomentumLabel { term: 1, .. }
        ));
        assert_eq!(expression, original);
    }
}
