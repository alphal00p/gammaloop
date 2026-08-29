use ::rustred::algebra::{CoefficientContext, ExactAlgebraLimits};
use ::rustred::family::presentation::{
    CommonMassScale, DenominatorRole, FamilyConventions, FamilyPresentation, MetricConvention,
    MomentumCombination, MomentumRouting, PhysicalPropagator, PropagatorConvention,
};
use ::rustred::family::{
    AffineDenominator, IntegralFamily, IntegralFamilyError, IntegralKey, ScalarProductCoordinate,
};
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

pub(super) fn reduce(
    expression: &mut VakintExpression,
    vakint: &Vakint,
    settings: &VakintSettings,
    _options: RustRedOptions,
) -> Result<(), TensorReductionError> {
    // The abstract one-loop bridge is independent of the matched mass value
    // and propagator power, so authenticate it once per expression rather
    // than once per term.
    let (presentation, coefficient_context) = one_loop_family(0)?;
    let heads = TensorHeads::try_new(S.k, S.p, S.g, S.dot).map_err(|source| {
        TensorReductionError::RustRedTensor {
            term: 0,
            source: source.into(),
        }
    })?;
    let loop_labels = vec![Atom::num(1)];
    let service = TensorService::try_new(
        &presentation,
        TensorLane::SingleScaleVacuum,
        heads,
        TensorMomenta::new(loop_labels.clone(), Vec::new()),
        TensorLimits::default(),
    )
    .map_err(|source| TensorReductionError::RustRedTensor { term: 0, source })?;

    // Keep the operation transactional: topology canonicalization and native
    // projection become visible only if every input term is admitted.
    let mut projected = expression.clone();
    for (term_index, term) in projected.0.iter_mut().enumerate() {
        project_term(
            term,
            term_index,
            vakint,
            settings,
            &service,
            &coefficient_context,
            &loop_labels,
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
    service: &TensorService<'_>,
    coefficient_context: &CoefficientContext,
    loop_labels: &[Atom],
) -> Result<(), TensorReductionError> {
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
    if matched.n_loops != 1 || matched.n_props != 1 {
        return Err(TensorReductionError::RustRedUnsupportedFamily {
            term: term_index,
            loop_count: matched.n_loops,
            propagator_count: matched.n_props,
        });
    }
    validate_explicit_input_routing(term, term_index, &replacement_rules)?;

    // This is the existing Vakint matcher/routing engine. In particular,
    // numerator substitutions are simultaneous; the adapter never invents a
    // second routing or dispatches on a topology name.
    term.canonicalize(settings, &replacement_rules, false)?;
    let propagator = get_prop_with_id(term.integral.as_view(), 1).ok_or(
        TensorReductionError::RustRedMissingPropagator {
            term: term_index,
            propagator: 1,
        },
    )?;
    let momentum = propagator
        .get(&vk_symbol!("q_"))
        .expect("Vakint's canonical propagator matcher always captures momentum");
    let expected_momentum = function!(S.k, Atom::num(1));
    if momentum != &expected_momentum {
        return Err(TensorReductionError::RustRedUnsupportedMomentum {
            term: term_index,
            propagator: 1,
            momentum: momentum.to_string(),
        });
    }
    let mass = propagator
        .get(&vk_symbol!("mUVsq_"))
        .expect("Vakint's canonical propagator matcher always captures mass");
    if mass.is_zero() || !is_exact_real_scalar(mass.as_view()) {
        return Err(TensorReductionError::RustRedUnsupportedMass {
            term: term_index,
            propagator: 1,
            mass: mass.to_string(),
        });
    }
    let power = propagator
        .get(&vk_symbol!("pow_"))
        .expect("Vakint's canonical propagator matcher always captures power");
    let power = get_integer_from_atom(power.as_view()).ok_or_else(|| {
        TensorReductionError::RustRedUnsupportedPower {
            term: term_index,
            propagator: 1,
            power: power.to_string(),
        }
    })?;

    let dimension = parse_dimension(settings)?;
    let base_integral =
        IntegralKey::try_new([power]).map_err(|source| TensorReductionError::RustRedTensor {
            term: term_index,
            source: source.into(),
        })?;
    let projection = service
        .project(&term.numerator, &base_integral)
        .map_err(|source| TensorReductionError::RustRedTensor {
            term: term_index,
            source,
        })?;
    let projected_numerator = projection_to_vakint(
        term_index,
        &projection,
        coefficient_context,
        &dimension,
        mass,
        loop_labels,
    )?;
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
    replacement_rules: &crate::ReplacementRules,
) -> Result<(), TensorReductionError> {
    // Short topology notation carries the registry's canonical routing by
    // construction. A graph-shaped input, however, must not let the broad
    // topology matcher erase an inadmissible denominator momentum before the
    // RustRed family proof is built.
    if !contains_function_head(term.integral.as_view(), S.prop) {
        return Ok(());
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

fn one_loop_family(
    term_index: usize,
) -> Result<(FamilyPresentation, CoefficientContext), TensorReductionError> {
    let coefficients = CoefficientContext::try_new(["d", "m2"])
        .map_err(|source| TensorReductionError::RustRedCoefficientContext { source })?;
    let dimension = coefficients
        .parameter("d")
        .expect("the adapter registers the dimension parameter");
    let mass = coefficients
        .parameter("m2")
        .expect("the adapter registers the mass parameter");
    let negative_mass = coefficients
        .try_neg(&mass, ExactAlgebraLimits::default())
        .map_err(|source| TensorReductionError::RustRedFamily {
            term: term_index,
            source: IntegralFamilyError::ExactAlgebra(source),
        })?;
    let family = IntegralFamily::new(
        "vakint-native-one-loop-vacuum",
        vec!["k1".to_owned()],
        Vec::new(),
        coefficients.clone(),
        dimension,
        vec![AffineDenominator::new(
            negative_mass,
            vec![coefficients.one()],
        )],
        Vec::new(),
        vec![coefficients.zero()],
    )
    .map_err(|source| TensorReductionError::RustRedFamily {
        term: term_index,
        source,
    })?;
    let presentation = FamilyPresentation::try_new(
        family,
        vec![DenominatorRole::Physical(PhysicalPropagator::new(
            "propagator-1".to_owned(),
            MomentumCombination::new(vec![coefficients.one()], Vec::new()),
            mass.clone(),
        ))],
        MomentumRouting::new(
            vec!["k1".to_owned()],
            Vec::new(),
            vec![vec![coefficients.one()]],
            vec![Vec::new()],
            Vec::new(),
        ),
        FamilyConventions::new(
            MetricConvention::MinkowskiMostlyMinus,
            PropagatorConvention::MOMENTUM_SQUARED_MINUS_MASS_SQUARED,
        ),
        Some(CommonMassScale::new(mass)),
    )
    .map_err(|source| TensorReductionError::RustRedPresentation {
        term: term_index,
        source,
    })?;
    Ok((presentation, coefficients))
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
            Ok(S.dot(function!(S.k, left.clone()), function!(S.k, right.clone())))
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
