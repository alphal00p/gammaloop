//! Apply an externally hosted package to a native five-loop DOT fixture.
//! The package's exact entry union decides whether this example is admitted.

#[path = "five_loop_cancellation_fixtures/mod.rs"]
mod cancellation_fixtures;
#[path = "saved_application/dot.rs"]
mod saved_application_dot;

use std::collections::{BTreeMap, BTreeSet};
use std::{path::PathBuf, time::Instant};

use feynkit_graph::{FeynmanDiagram, IntegralFamily};
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder};
use symbolica::id::Replacement;
use symbolica::prelude::{Integer, Rational};
use vakint::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions, Vakint,
    VakintExpression, VakintSettings, diagram_integral::DiagramIntegralOptions, vakint_parse,
    vakint_symbol,
};

struct ApplicationLogger;
impl log::Log for ApplicationLogger {
    fn enabled(&self, metadata: &log::Metadata<'_>) -> bool {
        matches!(
            metadata.target(),
            "vakint::saved_application" | "vakint::five_loop_match"
        )
    }
    fn log(&self, record: &log::Record<'_>) {
        if self.enabled(record.metadata()) {
            eprintln!("{}", record.args());
        }
    }
    fn flush(&self) {}
}
static LOGGER: ApplicationLogger = ApplicationLogger;

fn rerouted_fixture(input: &Atom) -> Result<Atom, Box<dyn std::error::Error>> {
    let routing = [
        ("k(1,idx___)", "k(2,idx___)+k(3,idx___)"),
        ("k(2,idx___)", "-k(1,idx___)"),
        ("k(3,idx___)", "k(3,idx___)"),
        ("k(4,idx___)", "k(5,idx___)"),
        ("k(5,idx___)", "k(4,idx___)"),
    ]
    .iter()
    .map(|(source, target)| {
        Ok(Replacement::new(
            vakint_parse!(source)?.to_pattern(),
            vakint_parse!(target)?.to_pattern(),
        ))
    })
    .collect::<Result<Vec<_>, Box<dyn std::error::Error>>>()?;
    Ok(input
        .replace_multiple(&routing)
        .replace(vakint_parse!("prop(id_,edge(a_,b_),q_,mass_,power_)")?.to_pattern())
        .with(vakint_parse!("prop(100-id_,edge(100-b_,100-a_),-q_,mass_,power_)")?.to_pattern()))
}

fn expanded_propagator_momenta(input: &Atom) -> Atom {
    input.replace_map(|term, _context, out| {
        if let AtomView::Fun(function) = term
            && function.get_symbol() == vakint_symbol!("prop")
        {
            let arguments = function.iter().collect::<Vec<_>>();
            **out = FunctionBuilder::new(function.get_symbol())
                .add_arg(arguments[0])
                .add_arg(arguments[1])
                .add_arg(arguments[2].expand())
                .add_arg(arguments[3])
                .add_arg(arguments[4])
                .finish();
        }
    })
}

fn verify_routing_inverse(
    original: &Atom,
    transformed: &Atom,
) -> Result<(), Box<dyn std::error::Error>> {
    let inverse = [
        ("k(1,idx___)", "-k(2,idx___)"),
        ("k(2,idx___)", "k(1,idx___)-k(3,idx___)"),
        ("k(3,idx___)", "k(3,idx___)"),
        ("k(4,idx___)", "k(5,idx___)"),
        ("k(5,idx___)", "k(4,idx___)"),
    ]
    .iter()
    .map(|(source, target)| {
        Ok(Replacement::new(
            vakint_parse!(source)?.to_pattern(),
            vakint_parse!(target)?.to_pattern(),
        ))
    })
    .collect::<Result<Vec<_>, Box<dyn std::error::Error>>>()?;
    let restored = transformed
        .replace(vakint_parse!("prop(id_,edge(a_,b_),q_,mass_,power_)")?.to_pattern())
        .with(vakint_parse!("prop(100-id_,edge(100-b_,100-a_),-q_,mass_,power_)")?.to_pattern())
        .replace_multiple(&inverse);
    if !(expanded_propagator_momenta(&restored) - expanded_propagator_momenta(original))
        .expand()
        .together()
        .is_zero()
    {
        return Err("the alternate-basis fixture failed exact inverse reconstruction".into());
    }
    println!("inverse_fixture_reconstruction=passed");
    Ok(())
}

fn admitted_denominator_numerator(
    vakint: &Vakint,
    settings: &VakintSettings,
    degree: usize,
) -> Result<Atom, Box<dyn std::error::Error>> {
    let topology = vakint_parse!(format!(
        "topo(I5L_32745(1,0,0,0,0,1,0,0,1,1,0,1,{}))",
        5 + degree,
    ))?;
    let canonical = vakint.to_canonical(settings, topology.as_view(), false)?;
    // Short-form loop names denote Vakint's canonical LMB, not the external
    // package's original k1..k5. Resolve the parent's physical q1 line through
    // its retained propagator ID before constructing the cancelling numerator.
    let square = physical_square(&canonical, 1)?;
    Ok((square - Atom::one()).pow(Atom::num(degree)) * canonical)
}

fn physical_square(canonical: &Atom, id: usize) -> Result<Atom, Box<dyn std::error::Error>> {
    let pattern = vakint_parse!(format!("prop({id},edge(a_,b_),q_,mass_,power_)"))?;
    let matched = canonical
        .pattern_match(&pattern.to_pattern(), None, None)
        .next()
        .ok_or("canonical parent is missing the requested physical line")?;
    let momentum = matched[&vakint_symbol!("q_")].clone();
    Ok(vakint_symbol!("dot").call((&momentum, &momentum)))
}

fn homogeneous_tadpole_control(
    vakint: &Vakint,
    settings: &VakintSettings,
) -> Result<Atom, Box<dyn std::error::Error>> {
    let parent = vakint_parse!("topo(I5L_32745(1,0,0,0,0,1,0,0,1,1,0,1,10))")?;
    let parent = vakint.to_canonical(settings, parent.as_view(), false)?;
    // Both physical rows are absent from the positive support. The common
    // masses cancel in D1-D2, leaving a homogeneous momentum polynomial.
    let numerator = (physical_square(&parent, 1)? - physical_square(&parent, 2)?).pow(5);
    Ok(numerator * parent)
}

type GuardFactorCache = BTreeMap<Atom, BTreeSet<Atom>>;

/// Test-only domain comparison. Native factorization supplies every factor;
/// this only interprets products/powers as nonzero obligations. Returned
/// application guards are never replaced or discarded.
fn guard_domain(
    values: &[Atom],
    cache: &mut GuardFactorCache,
) -> Result<BTreeSet<Atom>, Box<dyn std::error::Error>> {
    fn collect(
        value: AtomView,
        out: &mut BTreeSet<Atom>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        if value.is_zero() {
            return Err("a retained nonzero condition vanished".into());
        }
        match value {
            AtomView::Num(_) => {
                Rational::try_from(value)
                    .map_err(|_| "only proved nonzero rational units may be omitted")?;
            }
            AtomView::Mul(product) => {
                for factor in product.iter() {
                    collect(factor, out)?;
                }
            }
            AtomView::Pow(power) if Integer::try_from(power.get_exp()).is_ok() => {
                if Integer::try_from(power.get_exp())? != Integer::from(0) {
                    collect(power.get_base(), out)?;
                }
            }
            _ => {
                out.insert(value.to_owned());
            }
        }
        Ok(())
    }
    let mut result = BTreeSet::new();
    for condition in values {
        let condition = condition.expand().together();
        if !cache.contains_key(&condition) {
            let mut factors = BTreeSet::new();
            collect(condition.factor().as_view(), &mut factors)?;
            cache.insert(condition.clone(), factors);
        }
        result.extend(cache[&condition].iter().cloned());
    }
    Ok(result)
}

fn require_same_value(
    left: &vakint::RustRedSymbolicEvaluation,
    right: &vakint::RustRedSymbolicEvaluation,
    cache: &mut GuardFactorCache,
    guard_discrepancies: &mut Vec<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    if !(&left.value - &right.value).expand().together().is_zero()
        || left.families != right.families
    {
        return Err("equivalent input changed the exact qualified master expression".into());
    }
    let left_guards = guard_domain(&left.nonzero_conditions, cache)?;
    let right_guards = guard_domain(&right.nonzero_conditions, cache)?;
    if left_guards != right_guards {
        let bounded = |values: Vec<&Atom>| {
            values
                .iter()
                .take(8)
                .map(|value| {
                    value
                        .to_canonical_string()
                        .chars()
                        .take(512)
                        .collect::<String>()
                })
                .collect::<Vec<_>>()
        };
        let left_only = left_guards.difference(&right_guards).collect::<Vec<_>>();
        let right_only = right_guards.difference(&left_guards).collect::<Vec<_>>();
        let diagnostic = format!(
            "equivalent value has distinct native-factored guard domains: left_only_count={} right_only_count={} left_only_prefix={:?} right_only_prefix={:?}",
            left_only.len(),
            right_only.len(),
            bounded(left_only),
            bounded(right_only),
        );
        eprintln!("guard_domain_discrepancy={diagnostic}");
        guard_discrepancies.push(diagnostic);
        return Ok(());
    }
    println!(
        "native_factored_guard_domain_parity=passed factors={}",
        left_guards.len()
    );
    Ok(())
}

fn check_independent_tadpoles(
    result: &vakint::RustRedSymbolicEvaluation,
) -> Result<(), Box<dyn std::error::Error>> {
    use rustred::{
        family::IntegralKey,
        foundry::artifact::ClosedArtifact,
        reduction::{Reducer, ReductionLimits},
    };
    let artifact =
        ClosedArtifact::decode_durable(include_bytes!("../artifacts/1L/ibps/rules.rrbin"))?;
    let label = result
        .families
        .iter()
        .find_map(|(label, fingerprint)| {
            (fingerprint == artifact.family_fingerprint()).then_some(label)
        })
        .ok_or("five-tadpole result did not use the independently shipped one-loop family")?;
    let dimension = artifact.family().dimension().to_expression();
    let physical_dimension = vakint_parse!("4-2*eps")?;
    let mut expected = Atom::one();
    for power in [1, 1, 1, 1, 5] {
        let (reduced, conditions, _) = Reducer::reduce_artifact_unit_mass_guarded(
            &artifact,
            &IntegralKey::try_new([power])?,
            ReductionLimits::default(),
        )?;
        let mut factor = Atom::Zero;
        for (integral, coefficient) in reduced.terms() {
            let master = vakint_parse!(format!(
                "RustRedMaster({},{})",
                label.to_canonical_string(),
                integral.powers()[0]
            ))?;
            factor += coefficient
                .to_expression()
                .replace(dimension.to_pattern())
                .with(physical_dimension.to_pattern())
                * master;
        }
        for condition in conditions {
            let condition = condition
                .to_expression()
                .replace(dimension.to_pattern())
                .with(physical_dimension.to_pattern());
            if !result.nonzero_conditions.contains(&condition) {
                return Err(format!("independent factor condition was lost: {condition}").into());
            }
        }
        expected *= factor;
    }
    if !(&result.value - expected).expand().together().is_zero() {
        return Err(
            "five-loop result differs from five independent guarded one-loop reductions".into(),
        );
    }
    println!("independent_five_tadpole_parity=passed");
    Ok(())
}

/// Native five-tadpole chain, with four zero-flow bridges authenticated from
/// its stored routing. Only one bridge is restored in the second input.
fn massive_bridge_fixtures() -> Result<Vec<(String, Atom, Atom, Atom)>, Box<dyn std::error::Error>>
{
    let diagram = FeynmanDiagram::from_dot(
        Model::phi_3_4(),
        "digraph five_tadpoles {\
         a -> a [particle=\"phi\"]; a -> b [particle=\"phi\"];\
         b -> b [particle=\"phi\"]; b -> c [particle=\"phi\"];\
         c -> c [particle=\"phi\"]; c -> d [particle=\"phi\"];\
         d -> d [particle=\"phi\"]; d -> e [particle=\"phi\"];\
         e -> e [particle=\"phi\"]; }",
    )?;
    assert_eq!(diagram.loop_count(), 5);
    let original = diagram.propagator_family(&Kinematics::new())?;
    let mut edges = diagram.edges().collect::<Vec<_>>();
    edges.sort_by_key(|(id, _, _)| *id);
    let mut powers = vec![Integer::from(0); original.denominators().len()];
    let mut bridges = Vec::new();
    let mut tadpoles = Vec::new();
    for (slot, (id, _, _)) in edges.iter().enumerate() {
        if diagram.loop_momentum_basis().edge_signatures[id]
            .loops
            .integer_coefficients()
            .iter()
            .all(|value| *value == 0)
        {
            bridges.push(slot);
        } else {
            tadpoles.push(slot);
            powers[slot] = Integer::from(1);
        }
    }
    assert_eq!(bridges.len(), 4);
    assert_eq!(tadpoles.len(), 5);
    powers[tadpoles[4]] = Integer::from(5);
    let symbols = original.denominators()[bridges[0]].get_all_symbols(true);
    assert_eq!(symbols.len(), 1);
    let mass = Atom::var(*symbols.iter().next().unwrap());
    let mut fixtures = Vec::new();
    for mass_value in [Some(1), Some(2), None] {
        let substitutions = mass_value
            .map(|value| (mass.clone(), Atom::num(value)))
            .into_iter()
            .collect::<Vec<_>>();
        let replacements = substitutions
            .iter()
            .map(|(source, target)| Replacement::new(source.clone(), target.clone()))
            .collect::<Vec<_>>();
        let family = IntegralFamily::new(
            original.loop_momenta().to_vec(),
            vec![],
            original
                .denominators()
                .iter()
                .map(|denominator| denominator.replace_multiple(&replacements))
                .collect(),
            original.kinematics(),
        )?;
        let convert = |bridge_power: i64| -> Result<Atom, Box<dyn std::error::Error>> {
            let mut powers = powers.clone();
            powers[bridges[0]] = Integer::from(bridge_power);
            Ok(Atom::from(VakintExpression::from_diagram(
                &diagram,
                &family,
                &Atom::one(),
                &DiagramIntegralOptions {
                    powers: Some(powers),
                    parameter_substitutions: substitutions.clone(),
                    ..Default::default()
                },
            )?))
        };
        let baseline = convert(0)?;
        let bridge = convert(1)?;
        let denominator = family.denominators()[bridges[0]].clone();
        if !(&bridge - &baseline * denominator.clone().pow(-1))
            .expand()
            .together()
            .is_zero()
        {
            return Err("native massive bridge extraction changed the five-loop input".into());
        }
        fixtures.push((
            format!("mass_{mass_value:?}"),
            baseline,
            bridge,
            denominator,
        ));
    }
    Ok(fixtures)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    log::set_logger(&LOGGER)?;
    log::set_max_level(log::LevelFilter::Debug);
    let args = std::env::args_os()
        .skip(1)
        .map(PathBuf::from)
        .collect::<Vec<_>>();
    let scalar_dot = saved_application_dot::ScalarDotInput::from_arguments(&args)?;
    if args.is_empty() || (scalar_dot.is_none() && !(1..=2).contains(&args.len())) {
        return Err("usage: five_loop_saved_application APPLICATION_DIRECTORY [VACUUM_DOT|--admitted-tadpoles|--admitted-rank11|--admitted-suite|--dot PATH --powers CSV]".into());
    }
    let mut settings = VakintSettings {
        epsilon_symbol: "eps".into(),
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                substitute_masters: false,
            },
        )]),
        form_exe_path: "/five-loop-example-must-not-use-form".into(),
        integral_normalization_factor: LoopNormalizationFactor::Custom("1".into()),
        ..Default::default()
    };
    let catalog_started = Instant::now();
    let vakint = Vakint::new()?;
    println!(
        "vakint_catalogue_construction_seconds={:.6}",
        catalog_started.elapsed().as_secs_f64()
    );
    println!(
        "five_loop_catalogue_before_load={:?}",
        vakint.topologies.five_loop_catalogue_statistics()
    );
    if args.len() == 1 && args[0].as_os_str() == "--catalogue-only" {
        return Ok(());
    }
    if args.len() == 1 && args[0].as_os_str() == "--native-fixtures" {
        let fixtures = cancellation_fixtures::fixtures(&vakint)?;
        cancellation_fixtures::verify_fixture_shapes(&fixtures)?;
        println!(
            "native_positive_line_cancellation_fixtures={}",
            fixtures.len()
        );
        println!(
            "native_massive_bridge_fixtures={}",
            massive_bridge_fixtures()?.len()
        );
        return Ok(());
    }
    if args.len() == 1 && args[0].as_os_str() == "--routing-fixtures" {
        return Err("routing fixtures require APPLICATION_DIRECTORY --routing-fixtures; --native-fixtures needs no package".into());
    }
    let started = Instant::now();
    let vakint = vakint.with_rustred_application(&mut settings, &args[0])?;
    println!(
        "load_and_prepare_seconds={:.6}",
        started.elapsed().as_secs_f64()
    );
    println!(
        "five_loop_catalogue_after_load={:?}",
        vakint.topologies.five_loop_catalogue_statistics()
    );
    println!(
        "products={:?}",
        vakint.rustred_application_product_statistics()?
    );
    println!(
        "product_charts={:?}",
        vakint.rustred_application_product_chart_statistics()?
    );
    println!(
        "owner_preparation={:?}",
        vakint.rustred_application_owner_preparation_statistics()?
    );
    if args.len() == 1 {
        return Ok(());
    }
    let fixtures_started = Instant::now();
    let cancellation_controls = if args.get(1).is_some_and(|arg| {
        arg.as_os_str() == "--admitted-suite" || arg.as_os_str() == "--routing-fixtures"
    }) {
        let fixtures = cancellation_fixtures::fixtures(&vakint)?;
        cancellation_fixtures::verify_fixture_shapes(&fixtures)?;
        println!(
            "native_positive_line_cancellation_fixtures={}",
            fixtures.len()
        );
        fixtures
    } else {
        Vec::new()
    };
    println!(
        "native_fixture_preparation_seconds={:.6}",
        fixtures_started.elapsed().as_secs_f64()
    );
    if args[1].as_os_str() == "--routing-fixtures" {
        let bridges = massive_bridge_fixtures()?;
        println!("native_massive_bridge_fixtures={}", bridges.len());
        let mut scalar_settings = settings.clone();
        scalar_settings.use_dot_product_notation = true;
        for degree in [5, 11] {
            let original = admitted_denominator_numerator(&vakint, &settings, degree)?;
            let transformed = rerouted_fixture(&original)?;
            let indexed = Vakint::convert_from_dot_notation(transformed.as_view());
            println!("scope_degree={degree} stage=first_canonicalization");
            let direct = vakint.to_canonical(&settings, transformed.as_view(), false)?;
            let indexed = vakint.to_canonical(&settings, indexed.as_view(), false)?;
            println!("scope_degree={degree} stage=direct_dot_tensor_prepass");
            let direct = vakint.tensor_reduce(&scalar_settings, direct.as_view())?;
            println!("scope_degree={degree} stage=raw_indexed_tensor_prepass");
            let indexed = vakint.tensor_reduce(&scalar_settings, indexed.as_view())?;
            if !(&direct - indexed).expand().together().is_zero() {
                return Err(format!(
                    "initial canonicalization broke indexed scalar scope at degree{degree}"
                )
                .into());
            }
            let expanded = original.expand();
            println!("scope_degree={degree} stage=expanded_internal_scalar_tensor_prepass");
            let first = vakint.tensor_reduce(&scalar_settings, expanded.as_view())?;
            println!("scope_degree={degree} stage=internal_scalar_output_restore");
            let restored = vakint.tensor_reduce(&scalar_settings, first.as_view())?;
            if !(restored - expanded).expand().together().is_zero() {
                return Err(format!(
                    "internal tensor output lost expanded scalar scope at degree{degree}"
                )
                .into());
            }
            println!("full_indexed_and_expanded_scalar_scope_degree={degree}:passed");
        }
        let input = admitted_denominator_numerator(&vakint, &settings, 5)?;
        let canonical = vakint.to_canonical(&settings, input.as_view(), false)?;
        let transformed = rerouted_fixture(&canonical)?;
        verify_routing_inverse(&canonical, &transformed)?;
        let rematched = vakint.to_canonical(&settings, transformed.as_view(), false)?;
        let tensor_reduced = vakint.tensor_reduce(&settings, rematched.as_view())?;
        let scalar = vakint.tensor_reduce(&scalar_settings, tensor_reduced.as_view())?;
        if !(&scalar - &rematched).expand().together().is_zero() {
            return Err(
                "the native closed-scalar bridge changed the routed five-loop numerator".into(),
            );
        }
        println!("native_closed_scalar_scope_restoration=passed");
        for (stage, expression) in [
            ("constructed", input),
            ("canonical", canonical),
            ("transformed", transformed),
            ("rematched", rematched),
            ("tensor_reduced", tensor_reduced),
            ("scalar_restored", scalar),
        ] {
            println!(
                "stage={stage} expression={}",
                expression.to_canonical_string()
            );
        }
        return Ok(());
    }
    let input = if let Some(input) = &scalar_dot {
        let prepared = input.prepare()?;
        println!(
            "scalar_dot_loops={} physical_denominators={} powers={:?} mass_squared=1 numerator=1",
            prepared.diagram.loop_count(),
            prepared.family.denominators().len(),
            input.powers
        );
        prepared.expression
    } else if args[1].as_os_str() == "--admitted-tadpoles"
        || args[1].as_os_str() == "--admitted-suite"
    {
        // Generation54 point00162: 15-slot powers
        // [0,0,0,0,1,0,0,1,1,0,0,1,0,0,5]. The explicit parent physical
        // slots place q12 and q15 at positions 11 and 12 of this short form.
        vakint_parse!("topo(I5L_32745(1,0,0,0,0,1,0,0,1,1,0,1,5))")?
    } else if args[1].as_os_str() == "--admitted-rank11" {
        // Generation54 point00000. This is an inverse-denominator numerator,
        // not a claim that all homogeneous momentum polynomials are admitted.
        admitted_denominator_numerator(&vakint, &settings, 11)?
    } else {
        let diagram =
            FeynmanDiagram::from_dot(Model::phi_3_4(), &std::fs::read_to_string(&args[1])?)?;
        if diagram.loop_count() != 5 {
            return Err("the DOT fixture must have five loops".into());
        }
        let family = diagram.propagator_family(&Kinematics::new())?;
        if family.denominators().len() != 12 {
            return Err("this degree-ten example requires a maximal twelve-line fixture".into());
        }
        let mut powers = vec![1.into(); family.denominators().len()];
        powers[0] = 2.into();
        powers[1] = 2.into();
        // The common masses cancel: momentum degree ten, scalar signed difference 9.
        // Each lowered/routed endpoint must still pass the original package union.
        let numerator = (&family.denominators()[0] - &family.denominators()[1]).pow(Atom::num(5));
        let input = VakintExpression::from_diagram(
            &diagram,
            &family,
            &numerator,
            &DiagramIntegralOptions {
                powers: Some(powers),
                ..Default::default()
            },
        )?;
        Atom::from(input)
    };
    let mut inputs = vec![input];
    let mut guard_factors = GuardFactorCache::new();
    let mut guard_discrepancies = Vec::new();
    if args[1].as_os_str() == "--admitted-suite" {
        inputs.push(admitted_denominator_numerator(&vakint, &settings, 5)?);
        inputs.push(admitted_denominator_numerator(&vakint, &settings, 11)?);
    }
    for (case, input) in inputs.into_iter().enumerate() {
        println!("case={case}");
        let started = Instant::now();
        let first = vakint.evaluate_rustred_symbolic(&settings, input.as_view())?;
        println!(
            "first_application_seconds={:.6}",
            started.elapsed().as_secs_f64()
        );
        let started = Instant::now();
        let warm = vakint.evaluate_rustred_symbolic(&settings, input.as_view())?;
        println!(
            "warm_application_seconds={:.6}",
            started.elapsed().as_secs_f64()
        );
        if first.value != warm.value
            || first.nonzero_conditions != warm.nonzero_conditions
            || first.families != warm.families
        {
            return Err("warm application changed the exact result".into());
        }
        println!("phase=expanded_input");
        let expanded = vakint.evaluate_rustred_symbolic(&settings, input.expand().as_view())?;
        require_same_value(
            &first,
            &expanded,
            &mut guard_factors,
            &mut guard_discrepancies,
        )?;
        let standard = vakint.evaluate(&settings, input.expand().as_view())?;
        if !(&first.value - standard).expand().together().is_zero() {
            return Err("standard evaluation and detailed scalar handoff disagree".into());
        }
        println!("standard_and_detailed_expanded_value_parity=passed");
        let long_form = vakint.to_canonical(&settings, input.as_view(), false)?;
        let transformed = rerouted_fixture(&long_form)?;
        println!("canonical_fixture={}", long_form.to_canonical_string());
        println!("rerouted_fixture={}", transformed.to_canonical_string());
        verify_routing_inverse(&long_form, &transformed)?;
        println!("phase=rerouted_input");
        let rerouted = vakint.evaluate_rustred_symbolic(&settings, transformed.as_view())?;
        require_same_value(
            &first,
            &rerouted,
            &mut guard_factors,
            &mut guard_discrepancies,
        )?;
        println!("expanded_numerator_and_simultaneous_routing_value_parity=passed");
        let indexed_transformed = Vakint::convert_from_dot_notation(transformed.as_view());
        let indexed = vakint.evaluate_rustred_symbolic(&settings, indexed_transformed.as_view())?;
        require_same_value(
            &first,
            &indexed,
            &mut guard_factors,
            &mut guard_discrepancies,
        )?;
        println!("indexed_simultaneous_routing_value_parity=passed");
        println!("cache={:?}", vakint.rustred_application_cache_statistics()?);
        println!(
            "cache_work={:?}",
            vakint.rustred_application_cache_work_statistics()?
        );
        println!(
            "owner_preparation={:?}",
            vakint.rustred_application_owner_preparation_statistics()?
        );
        println!(
            "factor_cache={:?}",
            vakint.rustred_application_factor_cache_statistics()?
        );
        println!(
            "witness_cache={:?}",
            vakint.rustred_application_witness_cache_statistics()?
        );
        println!("nonzero_condition_count={}", first.nonzero_conditions.len());
        let output = first.value.to_canonical_string();
        println!("output_bytes={}", output.len());
        println!("{output}");
        if case == 0
            && (args[1].as_os_str() == "--admitted-tadpoles"
                || args[1].as_os_str() == "--admitted-suite")
        {
            check_independent_tadpoles(&first)?;
            let spectator = vakint_parse!("dot(p(1),p(2))*g(71,72)")?;
            let decorated_input = &input * &spectator;
            let decorated =
                vakint.evaluate_rustred_symbolic(&settings, decorated_input.as_view())?;
            let mut dot_settings = settings.clone();
            dot_settings.use_dot_product_notation = true;
            let mut decorated_dots =
                vakint.evaluate_rustred_symbolic(&dot_settings, decorated_input.as_view())?;
            decorated_dots.value =
                Vakint::convert_from_dot_notation(decorated_dots.value.as_view());
            require_same_value(
                &decorated,
                &decorated_dots,
                &mut guard_factors,
                &mut guard_discrepancies,
            )?;
            let expected = Vakint::convert_from_dot_notation((&first.value * spectator).as_view());
            if !(&decorated.value - expected).expand().together().is_zero() {
                return Err(
                    "internal compact handoff changed external scalar/tensor spectators".into(),
                );
            }
            println!("external_spectator_presentation_parity=passed");
            let massive = vakint_parse!("topo(I5L_32745(4,0,0,0,0,1,0,0,1,1,0,1,5))")?;
            let massive = vakint.evaluate_rustred_symbolic(&settings, massive.as_view())?;
            let expected = &first.value * vakint_parse!("4^(1-5*eps)")?;
            if !(&massive.value - expected).expand().together().is_zero() {
                return Err("nonunit common mass changed the five-loop scaling law".into());
            }
            println!("nonunit_mass_scaling=passed");
        }
        // Settings remain the sole enable path, including on a clone with a warm cache.
        let clone = vakint.clone();
        let mut disabled = settings.clone();
        disabled.rustred_artifact_path = None;
        if clone
            .evaluate_rustred_symbolic(&disabled, input.as_view())
            .is_ok()
        {
            return Err("removing the settings path did not disable a cached application".into());
        }
        let mut tightened = settings.clone();
        tightened
            .rustred_application_options
            .binary_limits
            .max_program_bytes = 0;
        if clone
            .evaluate_rustred_symbolic(&tightened, input.as_view())
            .is_ok()
        {
            return Err("a warm cache bypassed changed application options".into());
        }
        vakint.clear_rustred_application_cache()?;
        if vakint
            .rustred_application_witness_cache_statistics()?
            .retained_entries
            != 0
        {
            return Err("clearing the application left frontend witnesses cached".into());
        }
        if vakint
            .rustred_application_factor_cache_statistics()?
            .cached_factors
            != 0
        {
            return Err("clearing the application left lower-factor entries cached".into());
        }
        println!("cached_path_and_options_guards=passed");
    }
    if args[1].as_os_str() == "--admitted-suite" {
        for fixture in cancellation_controls {
            let cancelled =
                vakint.evaluate_rustred_symbolic(&settings, fixture.cancelled.as_view())?;
            let pinched = vakint.evaluate_rustred_symbolic(&settings, fixture.pinched.as_view())?;
            require_same_value(
                &cancelled,
                &pinched,
                &mut guard_factors,
                &mut guard_discrepancies,
            )?;
            let expanded = vakint
                .evaluate_rustred_symbolic(&settings, fixture.cancelled.expand().as_view())?;
            require_same_value(
                &cancelled,
                &expanded,
                &mut guard_factors,
                &mut guard_discrepancies,
            )?;
            println!(
                "positive_line_cancellation_value_parity={}:passed expected_key={:?} mass_squared={} positive_lines={:?} independent_dot_bytes={}",
                fixture.label,
                fixture.expected_key.powers(),
                fixture.mass_squared,
                fixture.positive_cancellations,
                fixture.pinched_dot.len(),
            );
        }
        for (label, baseline, bridge, denominator) in massive_bridge_fixtures()? {
            let baseline = vakint.evaluate_rustred_symbolic(&settings, baseline.as_view())?;
            let actual = vakint.evaluate_rustred_symbolic(&settings, bridge.as_view())?;
            let expected = &baseline.value * denominator.clone().pow(-1);
            if !(&actual.value - expected).expand().together().is_zero()
                || actual.families != baseline.families
            {
                return Err(format!("five-loop massive bridge parity failed: {label}").into());
            }
            if !denominator.get_all_symbols(true).is_empty()
                && !actual.nonzero_conditions.iter().any(|condition| {
                    (condition - &denominator).expand().together().is_zero()
                        || (condition + &denominator).expand().together().is_zero()
                })
            {
                return Err("the symbolic massive bridge lost its original mass pole guard".into());
            }
            println!("massive_bridge_parity={label}:passed");
        }
        // This particular scalar twelve-line corner (D=12) is outside the
        // pinned generation-54 entry union; this is not a global D cutoff.
        // neither a loaded owner nor cached terminal values grant admission.
        let unsupported = vakint_parse!("topo(I5L_32745(1,1,1,1,1,1,1,1,1,1,1,1,1))")?;
        match vakint.evaluate_rustred_symbolic(&settings, unsupported.as_view()) {
            Err(vakint::VakintError::RustRedEvaluation(
                vakint::RustRedEvaluationError::SavedApplicationEntryNotAdmitted { .. },
            )) => {}
            Err(error) => {
                return Err(
                    format!("D=12 control failed for an unexpected reason: {error}").into(),
                );
            }
            Ok(_) => {
                return Err("the scalar D=12 corner escaped the generation-54 entry union".into());
            }
        }
        println!("unsupported_entry_after_cache_reuse=refused");
        println!("phase=homogeneous_degree10_tadpole_control");
        let homogeneous = homogeneous_tadpole_control(&vakint, &settings)?;
        let first = vakint.evaluate_rustred_symbolic(&settings, homogeneous.as_view())?;
        let warm = vakint.evaluate_rustred_symbolic(&settings, homogeneous.as_view())?;
        require_same_value(&first, &warm, &mut guard_factors, &mut guard_discrepancies)?;
        let expanded =
            vakint.evaluate_rustred_symbolic(&settings, homogeneous.expand().as_view())?;
        require_same_value(
            &first,
            &expanded,
            &mut guard_factors,
            &mut guard_discrepancies,
        )?;
        println!("homogeneous_degree10_tadpole_value_parity=passed");
    }
    if !guard_discrepancies.is_empty() {
        return Err(format!(
            "{} native-factored guard-domain discrepancies require provenance audit",
            guard_discrepancies.len()
        )
        .into());
    }
    println!("guard_domain_parity_all_comparisons=passed");
    Ok(())
}
