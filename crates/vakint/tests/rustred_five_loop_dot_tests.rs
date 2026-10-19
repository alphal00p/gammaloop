#![cfg(feature = "feynkit-ingress")]

#[allow(dead_code)] // Shared CLI helpers also support file-based input outside this suite.
#[path = "../examples/saved_application/dot.rs"]
mod saved_application_dot;

use std::{collections::BTreeSet, path::PathBuf, time::Instant};

use feynkit_graph::{DiagramEdge, DiagramVertex, FeynmanDiagram};
use feynkit_model::Model;
use rustred::{family::IntegralKey, solver::FiniteRootAdmission};
use rustred_app::saved_application_domains;
use saved_application_dot::{Result, ScalarDotInput};
use symbolica::atom::{Atom, AtomCore};
use vakint::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions, Vakint,
    VakintExpression, VakintSettings, vakint_parse,
};

fn arguments(values: &[&str]) -> Vec<PathBuf> {
    values.iter().map(PathBuf::from).collect()
}

fn banana_dot() -> Result<String> {
    let model = Model::phi_3_4();
    let particle = model.particle_id("phi")?;
    let mut builder = FeynmanDiagram::builder(model, "six_line_scalar_banana");
    let mut vertex = |name: &str| {
        builder.add_vertex(DiagramVertex {
            name: name.into(),
            interaction: None,
            numerator: Atom::one(),
        })
    };
    let left = vertex("left");
    let right = vertex("right");
    for _ in 0..6 {
        builder.add_edge(left, right, DiagramEdge::new(particle, false))?;
    }
    Ok(builder.build()?.to_dot()?)
}

#[test]
fn scalar_dot_arguments_preserve_named_modes_and_accept_both_option_orders() -> Result<()> {
    for mode in [
        "--admitted-tadpoles",
        "--admitted-rank11",
        "--admitted-suite",
        "old.dot",
    ] {
        assert!(ScalarDotInput::from_arguments(&arguments(&["package", mode]))?.is_none());
    }
    for values in [
        ["package", "--dot", "input.dot", "--powers", "1, 0,-2"],
        ["package", "--powers", "1, 0,-2", "--dot", "input.dot"],
    ] {
        assert_eq!(
            ScalarDotInput::from_arguments(&arguments(&values))?,
            Some(ScalarDotInput {
                path: PathBuf::from("input.dot"),
                powers: vec![1, 0, -2],
            })
        );
    }
    Ok(())
}

#[test]
fn scalar_dot_arguments_reject_missing_duplicate_and_invalid_values() {
    for values in [
        vec!["package", "--dot"],
        vec!["package", "--dot", "input.dot"],
        vec!["package", "--dot", "input.dot", "--powers", ""],
        vec!["package", "--dot", "input.dot", "--powers", "1,true"],
        vec!["package", "--dot", "input.dot", "--powers", "1,,2"],
        vec!["package", "--dot", "input.dot", "--unknown", "1"],
        vec![
            "package",
            "--dot",
            "input.dot",
            "--dot",
            "other.dot",
            "--powers",
            "1",
        ],
        vec![
            "package",
            "--dot",
            "input.dot",
            "--powers",
            "1",
            "--powers",
            "2",
        ],
    ] {
        assert!(
            ScalarDotInput::from_arguments(&arguments(&values)).is_err(),
            "{values:?}"
        );
    }
}

#[test]
fn scalar_dot_roundtrip_preserves_native_routing_and_refuses_wrong_power_count() -> Result<()> {
    Vakint::initialize_vakint_symbols();
    let dot = banana_dot()?;
    let prepared = ScalarDotInput::prepare_text(&dot, &[1, 1, 1, 1, 1, 4])?;
    assert_eq!(prepared.diagram.loop_count(), 5);
    assert_eq!(prepared.diagram.edges().count(), 6);
    assert_eq!(prepared.family.loop_momenta().len(), 5);
    assert!(prepared.family.external_momenta().is_empty());
    assert_eq!(prepared.family.denominators().len(), 6);
    prepared
        .diagram
        .loop_momentum_basis()
        .validate(&prepared.diagram)?;
    assert_eq!(prepared.diagram.loop_momentum_basis().loop_edges.len(), 5);
    let parsed = VakintExpression::try_from(prepared.expression.as_view())?;
    assert_eq!(parsed.0.len(), 1);
    assert_eq!(parsed.0[0].numerator, Atom::one());
    assert!(ScalarDotInput::prepare_text(&dot, &[1, 1, 1, 1, 4]).is_err());
    assert!(ScalarDotInput::prepare_text("not a DOT graph", &[1]).is_err());
    Ok(())
}

fn guards(conditions: &[Atom]) -> BTreeSet<Atom> {
    conditions
        .iter()
        .map(|condition| condition.expand().together())
        .collect()
}

#[test]
#[ignore = "requires VAKINT_SAVED_APPLICATION with the admitted scalar banana; symbolic reduction only"]
fn external_r0_connected_dot_matches_known_admitted_key_and_warm_repeat() -> Result<()> {
    let package = std::env::var_os("VAKINT_SAVED_APPLICATION")
        .ok_or("set VAKINT_SAVED_APPLICATION to an external completed runtime package")?;
    let package = PathBuf::from(package);
    let known_key = IntegralKey::try_new([0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 4])?;
    let domains = saved_application_domains::<15>(&package)?;
    let max_regions = domains.len();
    let admission = FiniteRootAdmission::try_new(
        domains
            .into_iter()
            .filter(|domain| !domain.is_empty)
            .map(|domain| domain.region),
        max_regions,
    )?;
    admission.validate_entry(&known_key)?;

    let mut settings = VakintSettings {
        epsilon_symbol: "eps".into(),
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: vakint::RustRedMasterMode::Retain,
            },
        )]),
        form_exe_path: "/native-dot-test-must-not-use-form".into(),
        integral_normalization_factor: LoopNormalizationFactor::Custom("1".into()),
        ..Default::default()
    };
    let load_started = Instant::now();
    let engine = Vakint::new()?.with_rustred_application(&mut settings, &package)?;
    println!(
        "native_dot_package_load_seconds={:.6}",
        load_started.elapsed().as_secs_f64()
    );
    let preparation_started = Instant::now();
    let dot = banana_dot()?;
    let prepared = ScalarDotInput::prepare_text(&dot, &[1, 1, 1, 1, 1, 4])?;
    prepared
        .diagram
        .loop_momentum_basis()
        .validate(&prepared.diagram)?;
    assert_eq!(prepared.diagram.loop_count(), 5);
    let expression = VakintExpression::try_from(prepared.expression.as_view())?;
    let witness = engine
        .topologies
        .match_topologies_to_user_input(expression.0[0].integral.as_view(), false)?
        .ok_or("the native DOT has no verified simultaneous momentum-routing witness")?;
    println!("native_dot_matched_routing={witness}");
    println!(
        "native_dot_input_preparation_seconds={:.6}",
        preparation_started.elapsed().as_secs_f64()
    );

    // Independent common-family input for the admitted key. The six signed
    // momenta sum to zero and span the five loop variables with unit Jacobian.
    // This is fixture data, never topology-dependent evaluator dispatch.
    let reference = vakint_parse!(
        "topo(\
        prop(1,edge(1,2),k(5),1,1)*\
        prop(2,edge(1,2),k(1)-k(3),1,1)*\
        prop(3,edge(1,2),k(2)-k(3),1,1)*\
        prop(4,edge(1,2),k(3)-k(5),1,1)*\
        prop(5,edge(1,2),-k(1)-k(2)+k(4),1,1)*\
        prop(6,edge(1,2),k(3)-k(4),1,4))"
    )?;
    let first_started = Instant::now();
    let first = engine.evaluate_rustred_symbolic(&settings, prepared.expression.as_view())?;
    println!(
        "native_dot_first_seconds={:.6}",
        first_started.elapsed().as_secs_f64()
    );
    let warm_started = Instant::now();
    let warm = engine.evaluate_rustred_symbolic(&settings, prepared.expression.as_view())?;
    println!(
        "native_dot_warm_seconds={:.6}",
        warm_started.elapsed().as_secs_f64()
    );
    assert_eq!(first.value, warm.value);
    assert_eq!(first.nonzero_conditions, warm.nonzero_conditions);
    assert_eq!(first.families, warm.families);
    let expected = engine.evaluate_rustred_symbolic(&settings, reference.as_view())?;
    assert!(
        (&first.value - &expected.value)
            .expand()
            .together()
            .is_zero()
    );
    assert_eq!(first.families, expected.families);
    assert_eq!(
        guards(&first.nonzero_conditions),
        guards(&expected.nonzero_conditions)
    );
    assert!(
        !first.value.is_zero(),
        "this massive scalar control must be nonzero"
    );
    assert!(!first.families.is_empty());
    println!("native_dot_known_key_exact_parity=passed");
    println!(
        "native_dot_cache={:?}",
        engine.rustred_application_cache_statistics()?
    );
    Ok(())
}
