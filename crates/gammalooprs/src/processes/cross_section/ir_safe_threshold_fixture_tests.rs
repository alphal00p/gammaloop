use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::Cursor,
    path::PathBuf,
    time::{SystemTime, UNIX_EPOCH},
};

use crate::{
    GammaLoopContextContainer,
    feyngen::GenerationType,
    graph::{
        Graph, autogen::Autogen, feynman_graph::FeynmanGraph, parse::from_dot::IntoGraph,
        threshold_counterterms::ThresholdCountertermSpec,
    },
    initialisation::test_initialise,
    integrands::process::{
        ProcessIntegrand, ProcessIntegrandImpl, cross_section::CrossSectionIntegrand,
        evaluate_profile_momentum_point,
    },
    model::Model,
    momentum::ThreeMomentum,
    observables::ThresholdCountertermComponentOccurrence,
    processes::{
        CrossSection, DotExportSettings, ProcessDefinition,
        threshold_counterterms::{ThresholdCountertermComponentKind, ThresholdCountertermSide},
    },
    settings::{
        GlobalSettings, RuntimeSettings,
        global::{GenerationSettings, OrientationPattern},
    },
    utils::{
        F,
        fitting::{constant_dropped_fit_points, log_log_slope_constant_dropped},
        load_generic_model,
    },
};
use linnet::half_edge::involution::{EdgeIndex, Orientation};
use spenso::algebra::complex::Complex;
use symbolica::{numerical_integration::Sample, state::State};

const GL297_DOT: &str =
    include_str!("../../../../../tests/resources/graphs/ir_safe_thresholds/GL297.dot");
const GL638_DOT: &str =
    include_str!("../../../../../tests/resources/graphs/ir_safe_thresholds/GL638.dot");
const EMPTY_DIRECTIVES: &str = "schema_version = 1\n";
const GL297_FORCE_CUTS: &[&[&str]] = &[
    &["e2", "e11", "e14"],
    &["e2", "e9", "e12", "e14"],
    &["e2", "e7", "e13", "e14"],
    &["e2", "e6", "e11", "e13"],
    &["e2", "e6", "e9", "e12", "e13"],
    &["e2", "e6", "e7"],
    &["e2", "e4", "e11", "e12"],
    &["e2", "e4", "e9"],
    &["e2", "e4", "e7", "e12", "e13"],
];
const GL638_FORCE_CUTS: &[&[&str]] = &[
    &["e2", "e6", "e12", "e13"],
    &["e2", "e6", "e10"],
    &["e2", "e6", "e7", "e13", "e14"],
    &["e2", "e4", "e12"],
    &["e2", "e4", "e10", "e13"],
    &["e2", "e4", "e7", "e14"],
];

struct TemporaryDirectory(PathBuf);

impl TemporaryDirectory {
    fn new(prefix: &str) -> Self {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "gammalooprs-{prefix}-{}-{unique}",
            std::process::id()
        ));
        fs::create_dir_all(&path).unwrap();
        Self(path)
    }
}

impl Drop for TemporaryDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

const GL297_DIRECTIVES: &str = r#"
schema_version = 1

[[cuts]]
edges = [2, 11, 14]

  [[cuts.thresholds]]
  edges = [5, 7]

    [[cuts.thresholds.counterterms]]
    name = "forced_1l"
    subspace = [5]

  [[cuts.thresholds]]
  edges = [3, 9]

    [[cuts.thresholds.counterterms]]
    name = "forced_1l"
    subspace = [3]
"#;

const GL297_CURE_DIRECTIVES: &str = r#"
schema_version = 1

[[cuts]]
edges = [2, 4, 9]

  [[cuts.thresholds]]
  edges = [5, 11, 13]

    [[cuts.thresholds.counterterms]]
    name = "forced_1l"
    subspace = [5]

  [[cuts.thresholds]]
  edges = [5, 7]

    [[cuts.thresholds.counterterms]]
    name = "forced_1l"
    subspace = [5]

[[cuts]]
edges = [2, 6, 7]

  [[cuts.thresholds]]
  edges = [3, 11, 12]

    [[cuts.thresholds.counterterms]]
    name = "forced_1l"
    subspace = [3]

  [[cuts.thresholds]]
  edges = [3, 9]

    [[cuts.thresholds.counterterms]]
    name = "forced_1l"
    subspace = [3]
"#;

const GL638_TARGET_DIRECTIVES: &str = r#"
schema_version = 1

[[cuts]]
edges = [2, 4, 12]

  # Keep this fixture focused on the duplicated (7,8) geometry. An all-disabled explicit list is
  # authoritative and therefore suppresses the unrelated legacy-default left threshold.
  [[cuts.thresholds]]
  edges = [8, 12, 14]

    [[cuts.thresholds.counterterms]]
    name = "fixture_disabled"
    disable = true

  # The verified orientation also exposes a right threshold equal to another physical-cut
  # geometry. It is unrelated to the two historical right partners exercised by this fixture.
  [[cuts.thresholds]]
  edges = [2, 6, 10]

    [[cuts.thresholds.counterterms]]
    name = "fixture_disabled"
    disable = true

  [[cuts.thresholds]]
  edges = [7, 8]

    [[cuts.thresholds.counterterms]]
    name = "intrinsic_1l"
    subspace = [7]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, eset(7, 8))^2 / (eta(effective, eset(7, 8))^2 + eta(effective, eset(8, 12, 14))^2)"

    [[cuts.thresholds.counterterms]]
    name = "embedded_2l"
    subspace = [3, 7]
    # Revalidated by enumerating every generated parent: omission is genuinely ambiguous.
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, eset(8, 12, 14))^2 / (eta(effective, eset(7, 8))^2 + eta(effective, eset(8, 12, 14))^2)"
"#;

const GL638_FULL_DIRECTIVES: &str = r#"
schema_version = 1

[[cuts]]
edges = [2, 6, 12, 13]

  [[cuts.thresholds]]
  edges = [7, 8]

    [[cuts.thresholds.counterterms]]
    name = "intrinsic_1l"
    subspace = [7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, eset(7, 8))^2 / (eta(effective, eset(7, 8))^2 + eta(effective, eset(8, 12, 14))^2)"

    [[cuts.thresholds.counterterms]]
    name = "embedded_2l"
    subspace = [3, 7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, eset(8, 12, 14))^2 / (eta(effective, eset(7, 8))^2 + eta(effective, eset(8, 12, 14))^2)"

[[cuts]]
edges = [2, 6, 10]

  [[cuts.thresholds]]
  edges = [8, 12, 14]

    [[cuts.thresholds.counterterms]]
    name = "shared_1l"
    subspace = [7]
    parent_lmb = [3, 4, 7, 10]

  [[cuts.thresholds]]
  edges = [7, 8]

    [[cuts.thresholds.counterterms]]
    name = "intrinsic_1l"
    subspace = [7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, eset(7, 8))^2 / (eta(effective, eset(7, 8))^2 + eta(effective, eset(8, 12, 14))^2)"

    [[cuts.thresholds.counterterms]]
    name = "embedded_2l"
    subspace = [3, 7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, eset(8, 12, 14))^2 / (eta(effective, eset(7, 8))^2 + eta(effective, eset(8, 12, 14))^2)"

[[cuts]]
edges = [2, 4, 12]

  [[cuts.thresholds]]
  edges = [7, 8]

    [[cuts.thresholds.counterterms]]
    name = "intrinsic_1l"
    subspace = [7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, eset(7, 8))^2 / (eta(effective, eset(7, 8))^2 + eta(effective, eset(8, 12, 14))^2)"

    [[cuts.thresholds.counterterms]]
    name = "embedded_2l"
    subspace = [3, 7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, eset(8, 12, 14))^2 / (eta(effective, eset(7, 8))^2 + eta(effective, eset(8, 12, 14))^2)"

[[cuts]]
edges = [2, 4, 10, 13]

  [[cuts.thresholds]]
  edges = [7, 8]

    [[cuts.thresholds.counterterms]]
    name = "intrinsic_1l"
    subspace = [7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, eset(7, 8))^2 / (eta(effective, eset(7, 8))^2 + eta(effective, eset(2, 4, 10, 13))^2)"

    [[cuts.thresholds.counterterms]]
    name = "embedded_2l"
    subspace = [3, 7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, eset(2, 4, 10, 13))^2 / (eta(effective, eset(7, 8))^2 + eta(effective, eset(2, 4, 10, 13))^2)"
"#;

fn runtime_settings() -> RuntimeSettings {
    toml::from_str(
        r#"
[kinematics]
e_cm = 1000.0

[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
  [500.0, 0.0, 0.0, 500.0],
  [500.0, 0.0, 0.0, -500.0],
]
helicities = ["summed_averaged", "summed_averaged"]
"#,
    )
    .unwrap()
}

fn generation_settings(
    force_cut: &[&str],
    orientation_pattern: &str,
    generate_integrated_uv: bool,
) -> GenerationSettings {
    let mut settings = GenerationSettings {
        force_cuts: vec![force_cut.iter().map(|edge| (*edge).to_string()).collect()],
        orientation_pattern: OrientationPattern::from_user_pattern(orientation_pattern).unwrap(),
        ..Default::default()
    };
    settings.uv.generate_integrated = generate_integrated_uv;
    settings
}

fn fixture_pool() -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .stack_size(256 * 1024 * 1024)
        .build()
        .unwrap()
}

fn selected_signature(cross_section: &CrossSection) -> String {
    let orientations = &cross_section.supergraphs[0]
        .derived_data
        .global_cff_expression
        .as_ref()
        .expect("preprocessing must store the selected CFF")
        .orientations;
    assert_eq!(
        orientations.len(),
        1,
        "fixture must generate exactly one orientation"
    );
    format!(
        "({})",
        orientations
            .first()
            .expect("validated one selected orientation")
            .data
            .orientation
            .iter()
            .map(|(_, orientation)| match orientation {
                Orientation::Default => "+",
                Orientation::Reversed => "-",
                Orientation::Undirected => "0",
            })
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn preprocess_fixture(
    dot: &str,
    directives: &str,
    force_cut: &[&str],
    orientation_pattern: &str,
    generate_integrated_uv: bool,
) -> CrossSection {
    let settings = generation_settings(force_cut, orientation_pattern, generate_integrated_uv);
    preprocess_fixture_with_settings(dot, directives, &settings)
}

fn preprocess_fixture_with_settings(
    dot: &str,
    directives: &str,
    settings: &GenerationSettings,
) -> CrossSection {
    test_initialise().unwrap();
    let model = load_generic_model("sm");
    let mut graph: Graph = dot.into_graph(&model).unwrap();
    graph.threshold_counterterms =
        Autogen::explicit(toml::from_str::<ThresholdCountertermSpec>(directives).unwrap());
    let process_definition =
        ProcessDefinition::from_graph_list(&[graph.clone()], GenerationType::CrossSection, &model)
            .unwrap();
    let mut cross_section =
        CrossSection::from_graph_list("NNLO_fixture".to_string(), vec![graph], &model).unwrap();
    let runtime = runtime_settings();
    let pool = fixture_pool();
    cross_section
        .preprocess(
            &model,
            &process_definition,
            settings,
            (&runtime).into(),
            &pool,
        )
        .unwrap();
    cross_section
}

fn build_fixture_integrand_with_settings(
    dot: &str,
    directives: &str,
    settings: &GenerationSettings,
) -> (CrossSection, Model) {
    let mut cross_section = preprocess_fixture_with_settings(dot, directives, settings);
    let model = load_generic_model("sm");
    let runtime = runtime_settings();
    let global = GlobalSettings {
        generation: settings.clone(),
        ..Default::default()
    };
    cross_section
        .build_integrand(
            &model,
            "NNLO_fixture",
            &global,
            (&runtime).into(),
            &fixture_pool(),
        )
        .unwrap();
    (cross_section, model)
}

fn gl297_soft_loop_momenta(
    direction: &str,
    lambda: f64,
    _top_mass: f64,
) -> Vec<ThreeMomentum<F<f64>>> {
    let momentum = |px, py, pz| ThreeMomentum::new(F(px), F(py), F(pz));
    let displacement = momentum(211.0 * lambda, -157.0 * lambda, 193.0 * lambda);
    let k1 = momentum(-91.0, 64.0, 83.0);
    let k2 = momentum(37.0, 126.0, -109.0);
    match direction {
        // e13 = -K1 + K3 because the incoming spatial momenta sum to zero.
        "e13" => vec![momentum(110.0, -73.0, 57.0), k1, k2, k1 + displacement],
        // e12 = K0 + K1 - K2 - K3 in the same centre-of-mass frame.
        "e12" => {
            let k3 = momentum(74.0, -52.0, 137.0);
            vec![k2 + k3 - k1 + displacement, k1, k2, k3]
        }
        _ => unreachable!("the GL297 fixture has only the two verified soft-gluon directions"),
    }
}

fn gl297_correlated_soft_threshold_loop_momenta(
    direction: &str,
    lambda: f64,
    top_mass: f64,
) -> Vec<ThreeMomentum<F<f64>>> {
    let momentum = |px, py, pz| ThreeMomentum::new(F(px), F(py), F(pz));
    let threshold_magnitude = ((500.0_f64).powi(2) - top_mass.powi(2)).sqrt();
    let threshold_momentum = momentum(
        threshold_magnitude * 3.0 / 13.0,
        -threshold_magnitude * 4.0 / 13.0,
        threshold_magnitude * 12.0 / 13.0,
    );
    let hard_displacement = momentum(47.0 * lambda, 29.0 * lambda, -31.0 * lambda);
    let soft_displacement = momentum(211.0 * lambda, -157.0 * lambda, 193.0 * lambda);
    let k1_generic = momentum(-91.0, 64.0, 83.0);
    let k2 = momentum(37.0, 126.0, -109.0);

    match direction {
        // q13=-K1+K3=lambda*d while both hard top momenta approach K*.
        "e13" => {
            let k1 = threshold_momentum + hard_displacement;
            let k3 = k1 + soft_displacement;
            vec![momentum(110.0, -73.0, 57.0), k1, k2, k3]
        }
        // q12=K0+K1-K2-K3=lambda*d and p3=K3+lambda*d, with K3 -> K*.
        "e12" => {
            let k3 = threshold_momentum + hard_displacement;
            let k0 = k2 + k3 - k1_generic + soft_displacement;
            vec![k0, k1_generic, k2, k3]
        }
        _ => unreachable!("the GL297 fixture has only the two verified soft-gluon directions"),
    }
}

fn gl297_full_generation_settings() -> GenerationSettings {
    let mut settings =
        generation_settings(GL297_FORCE_CUTS[0], "(+,+,-,+,+,+,+,-,0,-,0,-,-,+,+)", true);
    settings.force_cuts = GL297_FORCE_CUTS
        .iter()
        .map(|cut| cut.iter().map(|edge| (*edge).to_string()).collect())
        .collect();
    settings
}

type Gl297Profile = (f64, f64, Vec<f64>);
type Gl297Profiles = BTreeMap<&'static str, Gl297Profile>;
type Gl297ApproachFits = BTreeMap<&'static str, [(f64, f64); 3]>;
type Gl297MomentumPath = fn(&str, f64, f64) -> Vec<ThreeMomentum<F<f64>>>;

fn gl297_profiles(
    integrand: &mut CrossSectionIntegrand,
    model: &Model,
    evaluate_thresholds: bool,
    top_mass: f64,
    momentum_path: Gl297MomentumPath,
) -> Gl297Profiles {
    integrand.settings.subtraction.disable_threshold_subtraction = !evaluate_thresholds;
    let lambdas = constant_dropped_fit_points(&F(1.0e-2), &F(1.0e-6), 7).unwrap();

    ["e12", "e13"]
        .into_iter()
        .map(|direction| {
            let magnitudes = lambdas
                .iter()
                .map(|lambda| {
                    let result = evaluate_profile_momentum_point(
                        integrand,
                        model,
                        0,
                        Some(0),
                        momentum_path(direction, lambda.0, top_mass),
                        true,
                    )
                    .unwrap();
                    assert!(
                        !result.evaluation_metadata.is_nan,
                        "{direction} at lambda={} produced a non-finite result",
                        lambda.0,
                    );
                    F(result.integrand_result.norm_squared().sqrt().0)
                })
                .collect::<Vec<_>>();
            let fit = log_log_slope_constant_dropped(&lambdas, &magnitudes).unwrap();
            (
                direction,
                (
                    fit.slope.0,
                    fit.r_squared.0,
                    magnitudes.into_iter().map(|value| value.0).collect(),
                ),
            )
        })
        .collect()
}

fn gl297_correlated_approach_fits(top_mass: f64) -> Gl297ApproachFits {
    let lambdas = constant_dropped_fit_points(&F(1.0e-2), &F(1.0e-6), 7).unwrap();
    ["e12", "e13"]
        .into_iter()
        .map(|direction| {
            let mut magnitudes = [Vec::new(), Vec::new(), Vec::new()];
            for lambda in &lambdas {
                let loop_momenta =
                    gl297_correlated_soft_threshold_loop_momenta(direction, lambda.0, top_mass);
                let energy = |momentum: &ThreeMomentum<F<f64>>| {
                    (momentum.norm_squared() + F(top_mass.powi(2))).sqrt().0
                };
                let (soft_momentum, eta_first, eta_second) = match direction {
                    "e13" => {
                        let soft_momentum = loop_momenta[3] - loop_momenta[1];
                        (
                            soft_momentum,
                            2.0 * energy(&loop_momenta[1]) - 1000.0,
                            energy(&loop_momenta[1])
                                + energy(&loop_momenta[3])
                                + soft_momentum.norm_squared().sqrt().0
                                - 1000.0,
                        )
                    }
                    "e12" => {
                        let soft_momentum =
                            loop_momenta[0] + loop_momenta[1] - loop_momenta[2] - loop_momenta[3];
                        let p3 = loop_momenta[0] + loop_momenta[1] - loop_momenta[2];
                        (
                            soft_momentum,
                            2.0 * energy(&p3) - 1000.0,
                            energy(&p3)
                                + energy(&loop_momenta[3])
                                + soft_momentum.norm_squared().sqrt().0
                                - 1000.0,
                        )
                    }
                    _ => unreachable!(),
                };
                assert_ne!(eta_first, 0.0);
                assert_ne!(eta_second, 0.0);
                magnitudes[0].push(soft_momentum.norm_squared().sqrt());
                magnitudes[1].push(F(eta_first.abs()));
                magnitudes[2].push(F(eta_second.abs()));
            }
            let fits = magnitudes.map(|values| {
                let fit = log_log_slope_constant_dropped(&lambdas, &values).unwrap();
                (fit.slope.0, fit.r_squared.0)
            });
            (direction, fits)
        })
        .collect()
}

fn gl638_correlated_soft_threshold_loop_momenta(
    lambda: f64,
    top_mass: f64,
) -> Vec<ThreeMomentum<F<f64>>> {
    let momentum = |px, py, pz| ThreeMomentum::new(F(px), F(py), F(pz));
    let threshold_magnitude = ((500.0_f64).powi(2) - top_mass.powi(2)).sqrt();
    let threshold_momentum = momentum(
        threshold_magnitude * 3.0 / 13.0,
        -threshold_magnitude * 4.0 / 13.0,
        threshold_magnitude * 12.0 / 13.0,
    );
    let k0_displacement = momentum(47.0 * lambda, 29.0 * lambda, -31.0 * lambda);
    let k3_displacement = momentum(258.0 * lambda, -128.0 * lambda, 162.0 * lambda);
    vec![
        threshold_momentum + k0_displacement,
        momentum(-59.0, 101.0, -73.0),
        momentum(127.0, 41.0, -109.0),
        threshold_momentum + k3_displacement,
    ]
}

type Gl638ApproachFits = BTreeMap<&'static str, BTreeMap<&'static str, (f64, f64)>>;

fn gl638_correlated_approach_fits(
    integrand: &mut CrossSectionIntegrand,
    model: &Model,
    top_mass: f64,
) -> Gl638ApproachFits {
    let lambdas = constant_dropped_fit_points(&F(1.0e-2), &F(1.0e-6), 7).unwrap();
    let energy =
        |momentum: &ThreeMomentum<F<f64>>| (momentum.norm_squared() + F(top_mass.powi(2))).sqrt();

    [("negative", -1.0), ("positive", 1.0)]
        .into_iter()
        .map(|(branch, sign)| {
            let mut magnitudes = [
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            ];
            for lambda in &lambdas {
                let loop_momenta = gl638_correlated_soft_threshold_loop_momenta(
                    sign * lambda.0,
                    top_mass,
                );
                let q13 = loop_momenta[0] - loop_momenta[3];
                magnitudes[0].push(q13.norm_squared().sqrt());
                magnitudes[1]
                    .push((F(2.0) * energy(&loop_momenta[3]) - F(1000.0)).abs());
                magnitudes[2].push(
                    (energy(&loop_momenta[0])
                        + energy(&loop_momenta[3])
                        + q13.norm_squared().sqrt()
                        - F(1000.0))
                    .abs(),
                );

                let result = evaluate_profile_momentum_point(
                    integrand,
                    model,
                    0,
                    Some(0),
                    loop_momenta,
                    true,
                )
                .unwrap();
                assert!(
                    !result.evaluation_metadata.is_nan
                        && result.integrand_result.re.0.is_finite()
                        && result.integrand_result.im.0.is_finite(),
                    "GL638 {branch} branch at |lambda|={} produced a non-finite result",
                    lambda.0,
                );

                let mut original = Complex::new(F(0.0), F(0.0));
                let mut counterterms = Complex::new(F(0.0), F(0.0));
                let mut event_sum = Complex::new(F(0.0), F(0.0));
                let mut event_scale = 0.0_f64;
                let mut event_count = 0;
                for event_group in result.event_groups.iter() {
                    for event in event_group.iter() {
                        event_count += 1;
                        let decomposition = event
                            .additional_weights
                            .threshold_counterterms
                            .as_ref()
                            .expect("the cured GL638 profile must record its CT decomposition");
                        assert_eq!(event.weight, decomposition.total());
                        original += &decomposition.original;
                        event_sum += &event.weight;
                        event_scale += event.weight.norm_squared().sqrt().0;
                        for component in &decomposition.components {
                            assert!(
                                component
                                    .multiplier_values
                                    .iter()
                                    .all(|value| value.0.is_finite())
                                    && component.effective_multiplier.0.is_finite()
                                    && component.weighted.re.0.is_finite()
                                    && component.weighted.im.0.is_finite(),
                                "GL638 {branch} branch at |lambda|={} recorded a non-finite component",
                                lambda.0,
                            );
                            counterterms += &component.weighted;
                        }
                    }
                }
                assert!(event_count > 0, "the GL638 profile must generate events");
                let closure = (event_sum - result.integrand_result)
                    .norm_squared()
                    .sqrt()
                    .0;
                let scale = event_scale.max(f64::MIN_POSITIVE);
                assert!(
                    closure / scale < 1.0e-12,
                    "GL638 {branch} branch event sum does not close at |lambda|={}: relative residual {}",
                    lambda.0,
                    closure / scale,
                );
                magnitudes[3].push(original.norm_squared().sqrt());
                magnitudes[4].push(counterterms.norm_squared().sqrt());
                magnitudes[5].push(result.integrand_result.norm_squared().sqrt());
            }

            let fits = [
                "q13",
                "eta_5_10",
                "eta_5_12_13",
                "original",
                "counterterms",
                "total",
            ]
            .into_iter()
            .zip(magnitudes)
            .enumerate()
            .map(|(quantity_index, (quantity, values))| {
                let (fit_lambdas, fit_values) = if quantity_index < 3 {
                    (&lambdas[..], &values[..])
                } else {
                    (&lambdas[3..], &values[3..])
                };
                let fit = log_log_slope_constant_dropped(fit_lambdas, fit_values).unwrap();
                (quantity, (fit.slope.0, fit.r_squared.0))
            })
            .collect();
            (branch, fits)
        })
        .collect()
}

#[test]
fn gl297_selected_orientation_resolves_forced_one_loop_subspaces_with_full_uv() {
    std::thread::Builder::new()
        .name("gl297-ir-safe-threshold-fixture".to_string())
        .stack_size(64 * 1024 * 1024)
        .spawn(|| {
            let cross_section = preprocess_fixture(
                GL297_DOT,
                GL297_DIRECTIVES,
                &["e2", "e11", "e14"],
                "(+,+,-,+,+,+,+,-,0,-,0,-,-,+,+)",
                true,
            );
            assert_eq!(
                selected_signature(&cross_section),
                "(+,+,-,+,+,+,+,-,0,-,0,-,-,+,+)"
            );

            let graph = &cross_section.supergraphs[0];
            assert_eq!(graph.cuts.len(), 1, "fixture must retain only the target cut");
            let resolved = graph
                .derived_data
                .resolved_threshold_counterterms
                .as_ref()
                .expect("threshold directives must resolve");
            let all_lmbs = graph.derived_data.lmbs.as_ref().unwrap();
            let forced = resolved
                .variants
                .iter()
                .filter(|variant| variant.name == "forced_1l")
                .collect::<Vec<_>>();
            assert_eq!(forced.len(), 2);
            assert_eq!(
                forced
                    .iter()
                    .map(|variant| {
                        assert_eq!(variant.subspace_loop_count, 1);
                        let mut edges = variant
                            .subspace
                            .iter_basis_edges(all_lmbs)
                            .map(usize::from)
                            .collect::<Vec<_>>();
                        edges.sort_unstable();
                        edges
                    })
                    .collect::<BTreeSet<_>>(),
                BTreeSet::from([vec![3], vec![5]])
            );
            assert!(
                forced
                    .iter()
                    .any(|variant| variant.requested_subspace.as_deref()
                        == Some(&[EdgeIndex::from(3)])),
                "the non-generation-basis edge e3 must resolve through a genuine one-loop parent embedding",
            );

            // Revalidated against the current GL297 DOT: both requests use the same parent, but
            // topology assigns them to opposite cut sides and clips each fundamental cycle to a
            // distinct physical one-loop subgraph.
            let all_edges = graph
                .graph
                .iter_edges()
                .map(|(_, edge, _)| edge)
                .collect::<Vec<_>>();
            let geometries = forced
                .iter()
                .map(|variant| {
                    let requested = variant.requested_subspace.as_ref().unwrap()[0].0;
                    let mut active = variant
                        .subspace
                        .contains(&all_edges, &graph.graph)
                        .map(usize::from)
                        .collect::<Vec<_>>();
                    active.sort_unstable();
                    (
                        requested,
                        (
                            variant.side,
                            variant
                                .resolved_parent_lmb
                                .iter()
                                .map(|edge| edge.0)
                                .collect::<Vec<_>>(),
                            active,
                        ),
                    )
                })
                .collect::<BTreeMap<_, _>>();
            assert_eq!(
                geometries,
                BTreeMap::from([
                    (
                        3,
                        (
                            ThresholdCountertermSide::Right,
                            vec![3, 5, 11, 14],
                            vec![3, 4, 9, 12],
                        ),
                    ),
                    (
                        5,
                        (
                            ThresholdCountertermSide::Left,
                            vec![3, 5, 11, 14],
                            vec![5, 6, 7, 13],
                        ),
                    ),
                ])
            );
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn cross_section_dot_export_materializes_defaults_and_reimports_as_legacy() {
    std::thread::Builder::new()
        .name("cross-section-threshold-dot-export".to_string())
        .stack_size(64 * 1024 * 1024)
        .spawn(|| {
            test_initialise().unwrap();
            let model = load_generic_model("sm");
            let graph: Graph = GL297_DOT.into_graph(&model).unwrap();
            let mut cross_section = CrossSection::from_graph_list(
                "NNLO_fixture".to_string(),
                vec![graph.clone()],
                &model,
            )
            .unwrap();
            let wrapper = &cross_section.supergraphs[0];
            assert!(wrapper.graph.threshold_counterterms.autogenerated);

            let mut ordinary_before = String::new();
            wrapper
                .write_dot_fmt(&mut ordinary_before, &DotExportSettings::default())
                .unwrap();
            assert!(!ordinary_before.contains("threshold_counterterms"));
            let mut normalized_before = String::new();
            wrapper
                .write_dot_fmt(
                    &mut normalized_before,
                    &DotExportSettings {
                        include_autogenerated_fields: true,
                        ..DotExportSettings::default()
                    },
                )
                .unwrap();
            let normalized_before_graph: Graph =
                normalized_before.as_str().into_graph(&model).unwrap();
            assert!(!normalized_before_graph.threshold_counterterms.autogenerated);
            assert!(
                normalized_before_graph
                    .threshold_counterterms
                    .cuts
                    .is_empty()
            );

            let mut settings = generation_settings(
                &["e2", "e11", "e14"],
                "(+,+,-,+,+,+,+,-,0,-,0,-,-,+,+)",
                false,
            );
            // This regression exercises DOT provenance and the homogeneous/generalized runtime
            // branch only. UV subtraction and integrated threshold pieces are covered by the
            // dedicated cure/runtime fixtures and would merely duplicate expensive evaluators.
            settings.uv.softct = false;
            settings.uv.subtract_uv = false;
            settings.threshold_subtraction.disable_integrated_ct = true;
            let process_definition =
                ProcessDefinition::from_graph_list(&[graph], GenerationType::CrossSection, &model)
                    .unwrap();
            cross_section
                .preprocess(
                    &model,
                    &process_definition,
                    &settings,
                    (&runtime_settings()).into(),
                    &fixture_pool(),
                )
                .unwrap();
            let wrapper = &cross_section.supergraphs[0];
            let resolved = wrapper
                .derived_data
                .resolved_threshold_counterterms
                .as_ref()
                .unwrap();
            assert!(resolved.legacy_equivalent);
            let expected = resolved.materialized_spec(
                &wrapper.graph.threshold_counterterms,
                wrapper.derived_data.lmbs.as_ref().unwrap(),
            );
            assert!(!expected.cuts.is_empty());

            let mut ordinary_after = String::new();
            wrapper
                .write_dot_fmt(&mut ordinary_after, &DotExportSettings::default())
                .unwrap();
            assert!(!ordinary_after.contains("threshold_counterterms"));
            let mut normalized_after = String::new();
            wrapper
                .write_dot_fmt(
                    &mut normalized_after,
                    &DotExportSettings {
                        include_autogenerated_fields: true,
                        ..DotExportSettings::default()
                    },
                )
                .unwrap();
            let normalized_graph: Graph = normalized_after.as_str().into_graph(&model).unwrap();
            assert_eq!(*normalized_graph.threshold_counterterms, expected);
            assert!(wrapper.graph.threshold_counterterms.autogenerated);
            assert!(wrapper.graph.threshold_counterterms.cuts.is_empty());

            let normalized_process_definition = ProcessDefinition::from_graph_list(
                std::slice::from_ref(&normalized_graph),
                GenerationType::CrossSection,
                &model,
            )
            .unwrap();
            let mut round_tripped = CrossSection::from_graph_list(
                "NNLO_fixture".to_string(),
                vec![normalized_graph],
                &model,
            )
            .unwrap();
            round_tripped
                .preprocess(
                    &model,
                    &normalized_process_definition,
                    &settings,
                    (&runtime_settings()).into(),
                    &fixture_pool(),
                )
                .unwrap();
            assert!(
                round_tripped.supergraphs[0]
                    .derived_data
                    .resolved_threshold_counterterms
                    .as_ref()
                    .unwrap()
                    .legacy_equivalent,
                "materialized defaults must remain on the homogeneous legacy fast path",
            );

            let global = GlobalSettings {
                generation: settings,
                ..Default::default()
            };
            cross_section
                .build_integrand(
                    &model,
                    "NNLO_fixture",
                    &global,
                    (&runtime_settings()).into(),
                    &fixture_pool(),
                )
                .unwrap();
            let ProcessIntegrand::CrossSection(autogenerated_integrand) =
                cross_section.integrand.as_ref().unwrap()
            else {
                unreachable!("GL297 is a cross-section fixture")
            };
            let autogenerated_counterterm =
                &autogenerated_integrand.data.graph_terms[0].counterterm;
            assert!(
                autogenerated_counterterm.metadata_registry.is_none(),
                "absent/autogenerated legacy defaults must not allocate a metadata registry",
            );
            assert!(
                autogenerated_counterterm.variant_subspaces.is_none(),
                "absent/autogenerated legacy defaults must retain the homogeneous subspace path",
            );

            round_tripped
                .build_integrand(
                    &model,
                    "NNLO_fixture",
                    &global,
                    (&runtime_settings()).into(),
                    &fixture_pool(),
                )
                .unwrap();
            let ProcessIntegrand::CrossSection(materialized_integrand) =
                round_tripped.integrand.as_ref().unwrap()
            else {
                unreachable!("GL297 is a cross-section fixture")
            };
            let materialized_counterterm =
                &materialized_integrand.data.graph_terms[0].counterterm;
            assert!(
                materialized_counterterm.metadata_registry.is_some(),
                "explicit materialized defaults retain provenance metadata after re-import",
            );
            assert!(
                materialized_counterterm.variant_subspaces.is_none(),
                "semantic legacy equivalence must keep materialized defaults on the homogeneous subspace path",
            );
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn gl297_full_cut_forced_one_loop_subspaces_restore_soft_scaling() {
    std::thread::Builder::new()
        .name("gl297-ir-safe-threshold-cure".to_string())
        .stack_size(64 * 1024 * 1024)
        .spawn(|| {
            // A numerical LU regression must keep all nine process-valid cuts: its IR finiteness
            // is a dual cancellation across cuts. The orientation is selected only to avoid
            // generating other acyclic flows.
            let settings = gl297_full_generation_settings();
            let (mut legacy_cross_section, legacy_model) =
                build_fixture_integrand_with_settings(GL297_DOT, EMPTY_DIRECTIVES, &settings);
            assert_eq!(legacy_cross_section.supergraphs[0].cuts.len(), 9);
            assert_eq!(
                selected_signature(&legacy_cross_section),
                "(+,+,-,+,+,+,+,-,0,-,0,-,-,+,+)",
            );
            let legacy_graph = &legacy_cross_section.supergraphs[0];
            let legacy_resolved = legacy_graph
                .derived_data
                .resolved_threshold_counterterms
                .as_ref()
                .unwrap();
            let legacy_lmbs = legacy_graph.derived_data.lmbs.as_ref().unwrap();
            let top_mass = legacy_model
                .get_parameter("MT")
                .value
                .as_ref()
                .unwrap()
                .re
                .0;
            for (direction, fits) in gl297_correlated_approach_fits(top_mass) {
                for (slope, r_squared) in fits {
                    assert!(
                        (slope - 1.0).abs() < 0.002 && r_squared > 0.9999,
                        "{direction} correlated soft/threshold kinematics do not scale linearly: {fits:?}",
                    );
                }
            }
            for (cut_edges, threshold_edges, side, basis) in [
                (&[2, 4, 9][..], &[5, 11, 13][..], ThresholdCountertermSide::Left, &[5, 11][..]),
                (&[2, 4, 9][..], &[5, 7][..], ThresholdCountertermSide::Left, &[5, 11][..]),
                (&[2, 6, 7][..], &[3, 11, 12][..], ThresholdCountertermSide::Right, &[3, 11][..]),
                (&[2, 6, 7][..], &[3, 9][..], ThresholdCountertermSide::Right, &[3, 11][..]),
            ] {
                let variant = legacy_resolved
                    .variants
                    .iter()
                    .find(|variant| {
                        variant.associations.iter().any(|association| {
                            association
                                .cut_edges
                                .iter()
                                .map(|edge| edge.0)
                                .eq(cut_edges.iter().copied())
                                && association
                                    .threshold_edges
                                    .iter()
                                    .map(|edge| edge.0)
                                    .eq(threshold_edges.iter().copied())
                        })
                    })
                    .expect("each problematic legacy threshold association must resolve");
                assert_eq!(variant.side, side);
                assert_eq!(variant.subspace_loop_count, 2);
                assert_eq!(
                    variant
                        .subspace
                        .iter_basis_edges(legacy_lmbs)
                        .map(usize::from)
                        .collect::<Vec<_>>(),
                    basis,
                );
            }
            let ProcessIntegrand::CrossSection(legacy_integrand) =
                legacy_cross_section.integrand.as_mut().unwrap()
            else {
                unreachable!("GL297 is a cross-section fixture")
            };
            legacy_integrand.warm_up(&legacy_model).unwrap();
            let bare = gl297_profiles(
                legacy_integrand,
                &legacy_model,
                false,
                top_mass,
                gl297_soft_loop_momenta,
            );
            let correlated_bare = gl297_profiles(
                legacy_integrand,
                &legacy_model,
                false,
                top_mass,
                gl297_correlated_soft_threshold_loop_momenta,
            );
            let default = gl297_profiles(
                legacy_integrand,
                &legacy_model,
                true,
                top_mass,
                gl297_soft_loop_momenta,
            );
            let correlated_default = gl297_profiles(
                legacy_integrand,
                &legacy_model,
                true,
                top_mass,
                gl297_correlated_soft_threshold_loop_momenta,
            );

            let (mut forced_cross_section, forced_model) =
                build_fixture_integrand_with_settings(GL297_DOT, GL297_CURE_DIRECTIVES, &settings);
            assert_eq!(forced_cross_section.supergraphs[0].cuts.len(), 9);
            assert_eq!(
                selected_signature(&forced_cross_section),
                "(+,+,-,+,+,+,+,-,0,-,0,-,-,+,+)",
            );
            let forced_graph = &forced_cross_section.supergraphs[0];
            let forced_resolved = forced_graph
                .derived_data
                .resolved_threshold_counterterms
                .as_ref()
                .unwrap();
            let forced_lmbs = forced_graph.derived_data.lmbs.as_ref().unwrap();
            let forced_variants = forced_resolved
                .variants
                .iter()
                .filter(|variant| variant.name == "forced_1l")
                .collect::<Vec<_>>();
            assert_eq!(forced_variants.len(), 4);
            assert!(forced_variants.iter().all(|variant| {
                variant.subspace_loop_count == 1
                    && variant
                        .subspace
                        .iter_basis_edges(forced_lmbs)
                        .map(usize::from)
                        .collect::<Vec<_>>()
                        == match variant.side {
                            ThresholdCountertermSide::Left => vec![5],
                            ThresholdCountertermSide::Right => vec![3],
                            ThresholdCountertermSide::Amplitude => return false,
                        }
            }));
            let ProcessIntegrand::CrossSection(forced_integrand) =
                forced_cross_section.integrand.as_mut().unwrap()
            else {
                unreachable!("GL297 is a cross-section fixture")
            };
            forced_integrand.warm_up(&forced_model).unwrap();
            let forced = gl297_profiles(
                forced_integrand,
                &forced_model,
                true,
                top_mass,
                gl297_soft_loop_momenta,
            );
            let correlated_forced = gl297_profiles(
                forced_integrand,
                &forced_model,
                true,
                top_mass,
                gl297_correlated_soft_threshold_loop_momenta,
            );

            for direction in ["e12", "e13"] {
                let (bare_slope, bare_r_squared, _) = &bare[direction];
                let (default_slope, default_r_squared, _) = &default[direction];
                let (forced_slope, forced_r_squared, _) = &forced[direction];
                assert!(
                    [bare_r_squared, default_r_squared, forced_r_squared]
                        .into_iter()
                        .all(|r_squared| *r_squared > 0.9999),
                    "{direction} slope fits are not asymptotic enough: bare={:?}, default={:?}, forced={:?}",
                    bare[direction],
                    default[direction],
                    forced[direction],
                );
                assert!((bare_slope + 1.0).abs() < 0.02);
                assert!((default_slope + 3.0).abs() < 0.02);
                assert!((forced_slope + 1.0).abs() < 0.02);
                assert!(*default_slope < *bare_slope - 1.8);
                assert!((forced_slope - bare_slope).abs() < 0.02);
            }
            // This trajectory approaches the relevant soft direction and both associated
            // threshold surfaces at the same linear rate. It independently confirms that the
            // legacy maximal-subspace counterterms introduce two extra inverse powers, while
            // the forced one-loop variants restore the bare/IR-safe scaling.
            for direction in ["e12", "e13"] {
                let (bare_slope, bare_r_squared, _) = &correlated_bare[direction];
                let (default_slope, default_r_squared, _) = &correlated_default[direction];
                let (forced_slope, forced_r_squared, _) = &correlated_forced[direction];
                assert!(
                    [bare_r_squared, default_r_squared, forced_r_squared]
                        .into_iter()
                        .all(|r_squared| *r_squared > 0.9995),
                    "{direction} correlated fits are not asymptotic enough: bare={:?}, default={:?}, forced={:?}",
                    correlated_bare[direction],
                    correlated_default[direction],
                    correlated_forced[direction],
                );
                assert!((bare_slope + 1.0).abs() < 0.03);
                assert!((default_slope + 3.0).abs() < 0.02);
                assert!((forced_slope + 1.0).abs() < 0.02);
                assert!(*default_slope < *bare_slope - 1.8);
                assert!((forced_slope - bare_slope).abs() < 0.03);
            }
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn gl638_cartesian_structure_and_full_cut_runtime_roundtrip() {
    std::thread::Builder::new()
        .name("gl638-ir-safe-threshold-fixture".to_string())
        .stack_size(64 * 1024 * 1024)
        .spawn(|| {
            // Integrated UV is disabled only for this GL638 smoke fixture because direct-import
            // integrated-UV reconstruction currently fails in Vakint. Threshold local and
            // integrated pieces themselves remain enabled and are both generated.
            let cross_section = preprocess_fixture(
                GL638_DOT,
                GL638_TARGET_DIRECTIVES,
                &["e2", "e4", "e12"],
                "(+,+,+,+,-,-,-,+,-,0,+,0,+,-,+)",
                false,
            );
            assert_eq!(
                selected_signature(&cross_section),
                "(+,+,+,+,-,-,-,+,-,0,+,0,+,-,+)"
            );

            let graph = &cross_section.supergraphs[0];
            assert_eq!(
                graph.cuts.len(),
                1,
                "fixture must retain only the target cut"
            );
            assert_eq!(graph.derived_data.cut_group_data.cut_groups.len(), 1);
            let resolved = graph
                .derived_data
                .resolved_threshold_counterterms
                .as_ref()
                .expect("threshold directives must resolve");
            let group = &resolved.cross_section_cut_groups[crate::processes::CutGroupId::from(0)];
            assert_eq!(group.left.len(), 2);
            assert_eq!(
                group.right.len(),
                2,
                "resolved right variants: {:?}",
                group
                    .right
                    .iter()
                    .map(|variant_id| {
                        let variant = &resolved.variants[*variant_id];
                        variant
                            .associations
                            .iter()
                            .map(|association| {
                                association
                                    .threshold_edges
                                    .iter()
                                    .map(|edge| edge.0)
                                    .collect::<Vec<_>>()
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>(),
            );

            let left = group
                .left
                .iter()
                .map(|variant_id| &resolved.variants[*variant_id])
                .collect::<Vec<_>>();
            assert_eq!(
                left.iter()
                    .map(|variant| variant.name.as_str())
                    .collect::<BTreeSet<_>>(),
                BTreeSet::from(["embedded_2l", "intrinsic_1l"]),
            );
            assert_eq!(
                left.iter()
                    .map(|variant| variant.subspace_loop_count)
                    .collect::<BTreeSet<_>>(),
                BTreeSet::from([1, 2]),
            );
            assert!(left.iter().all(|variant| variant.multiplier.is_some()));
            assert!(
                left.iter().all(|variant| {
                    variant
                        .resolved_parent_lmb
                        .iter()
                        .map(|edge| edge.0)
                        .collect::<Vec<_>>()
                        == vec![3, 4, 7, 10]
                }),
                "the global two-loop request must rebase the complete cut group into the verified generation parent",
            );

            let right = group
                .right
                .iter()
                .map(|variant_id| &resolved.variants[*variant_id])
                .collect::<Vec<_>>();
            assert!(right.iter().all(|variant| variant.name == "default"));
            assert!(right.iter().all(|variant| {
                variant
                    .resolved_parent_lmb
                    .iter()
                    .map(|edge| edge.0)
                    .collect::<Vec<_>>()
                    == vec![3, 4, 7, 10]
            }));
            assert_eq!(
                right
                    .iter()
                    .flat_map(|variant| &variant.associations)
                    .map(|association| {
                        association
                            .threshold_edges
                            .iter()
                            .copied()
                            .map(usize::from)
                            .collect::<Vec<_>>()
                    })
                    .collect::<BTreeSet<_>>(),
                BTreeSet::from([vec![5, 10], vec![5, 12, 13]])
            );

            let generated =
                &graph.derived_data.threshold_counterterms[crate::processes::CutGroupId::from(0)];
            assert_eq!(generated.left_thresholds.len(), 2);
            assert_eq!(generated.right_thresholds.len(), 2);
            assert_eq!(generated.iterated.iter().count(), 4);

            // The original-side identity is a rescaling-map statement, not an extra integrand.
            // There is exactly one O_L*O_R container, four one-sided variant containers, and the
            // internal 2x2 Cartesian product. L/I splitting remains internal to each container.
            let original_terms = 1;
            let single_terms = generated.left_thresholds.len() + generated.right_thresholds.len();
            let pair_terms = generated.iterated.iter().count();
            assert_eq!((original_terms, single_terms, pair_terms), (1, 4, 4));
            assert_eq!(original_terms + single_terms + pair_terms, 9);

            // Numerical LU evaluation retains all six process-valid cuts. The target-only state
            // above is a structural Cartesian-product test and is not an IR-complete integrand.
            let mut full_generation = generation_settings(
                GL638_FORCE_CUTS[0],
                "(+,+,+,+,-,-,-,+,-,0,+,0,+,-,+)",
                false,
            );
            full_generation.force_cuts = GL638_FORCE_CUTS
                .iter()
                .map(|cut| cut.iter().map(|edge| (*edge).to_string()).collect())
                .collect();
            let mut cross_section = preprocess_fixture_with_settings(
                GL638_DOT,
                GL638_FULL_DIRECTIVES,
                &full_generation,
            );
            assert_eq!(cross_section.supergraphs[0].cuts.len(), 6);
            assert_eq!(
                selected_signature(&cross_section),
                "(+,+,+,+,-,-,-,+,-,0,+,0,+,-,+)",
            );

            let graph = &cross_section.supergraphs[0];
            let resolved = graph
                .derived_data
                .resolved_threshold_counterterms
                .as_ref()
                .unwrap();
            let expected_variant_cuts = BTreeSet::from([
                vec![2, 4, 10, 13],
                vec![2, 4, 12],
                vec![2, 6, 10],
                vec![2, 6, 12, 13],
            ]);
            for name in ["intrinsic_1l", "embedded_2l"] {
                assert_eq!(
                    resolved
                        .variants
                        .iter()
                        .filter(|variant| variant.name == name)
                        .flat_map(|variant| &variant.associations)
                        .map(|association| {
                            association
                                .cut_edges
                                .iter()
                                .map(|edge| edge.0)
                                .collect::<Vec<_>>()
                        })
                        .collect::<BTreeSet<_>>(),
                    expected_variant_cuts,
                );
            }
            let shared_variant = resolved
                .variants
                .iter()
                .find(|variant| variant.name == "shared_1l")
                .expect("the shared (8,12,14) surface must be forced into its one-loop subspace");
            assert_eq!(shared_variant.subspace_loop_count, 1);
            assert_eq!(
                shared_variant
                    .subspace
                    .iter_basis_edges(graph.derived_data.lmbs.as_ref().unwrap())
                    .map(usize::from)
                    .collect::<Vec<_>>(),
                vec![7],
            );
            assert!(shared_variant.associations.iter().any(|association| {
                association
                    .cut_edges
                    .iter()
                    .map(|edge| edge.0)
                    .eq([2, 6, 10])
                    && association
                        .threshold_edges
                        .iter()
                        .map(|edge| edge.0)
                        .eq([8, 12, 14])
            }));
            let target_left_variant_ids = resolved
                .variants
                .iter_enumerated()
                .filter(|(_, variant)| {
                    variant.associations.iter().any(|association| {
                        association
                            .cut_edges
                            .iter()
                            .map(|edge| edge.0)
                            .eq([2, 4, 12])
                            && association
                                .threshold_edges
                                .iter()
                                .map(|edge| edge.0)
                                .eq([7, 8])
                    })
                })
                .map(|(variant_id, _)| variant_id)
                .collect::<Vec<_>>();
            assert_eq!(target_left_variant_ids.len(), 2);
            let target_cut_group_id = resolved.variants[target_left_variant_ids[0]]
                .cut_group_id
                .unwrap();
            assert!(target_left_variant_ids.iter().all(|variant_id| {
                resolved.variants[*variant_id].cut_group_id == Some(target_cut_group_id)
            }));
            let target_group = &resolved.cross_section_cut_groups[target_cut_group_id];
            let target_right_variant_ids = target_group
                .right
                .iter()
                .copied()
                .filter(|variant_id| {
                    resolved.variants[*variant_id]
                        .associations
                        .iter()
                        .any(|association| {
                            matches!(
                                association
                                    .threshold_edges
                                    .iter()
                                    .map(|edge| edge.0)
                                    .collect::<Vec<_>>()
                                    .as_slice(),
                                [5, 10] | [5, 12, 13]
                            )
                        })
                })
                .collect::<Vec<_>>();
            assert_eq!(target_right_variant_ids.len(), 2);
            let target_generated =
                &graph.derived_data.threshold_counterterms[target_cut_group_id];
            assert!(target_left_variant_ids
                .iter()
                .all(|variant_id| target_generated.left_variant_ids.contains(variant_id)));
            assert!(target_right_variant_ids
                .iter()
                .all(|variant_id| target_generated.right_variant_ids.contains(variant_id)));

            // Filtering out unrelated legacy variants leaves one original, four one-sided terms,
            // and the intended internal 2x2 Cartesian product. The full integrand still retains
            // every other threshold needed for cross-cut cancellation.
            assert_eq!(1 + 4 + target_left_variant_ids.len() * target_right_variant_ids.len(), 9);

            let model = load_generic_model("sm");
            let mut runtime = runtime_settings();
            runtime.general.generate_events = true;
            runtime.general.store_additional_weights_in_event = true;
            let global = GlobalSettings {
                generation: full_generation,
                ..Default::default()
            };
            let pool = fixture_pool();
            cross_section
                .build_integrand(
                    &model,
                    "NNLO_fixture",
                    &global,
                    (&runtime).into(),
                    &pool,
                )
                .unwrap();
            let n_dim = cross_section.supergraphs[0].graph.get_loop_number() * 3;
            let samples = (0..2)
                .map(|sample_index| {
                    Sample::Continuous(
                        F(1.0),
                        (0..n_dim)
                            .map(|axis| {
                                F(0.12 + ((axis * 7 + sample_index * 5) % 17) as f64 * 0.043)
                            })
                            .collect(),
                    )
                })
                .collect::<Vec<_>>();
            let integrand = cross_section
                .integrand
                .as_mut()
                .expect("fixture integrand must be built");
            let metadata_registry = match &*integrand {
                ProcessIntegrand::CrossSection(cross_section_integrand) => cross_section_integrand
                    .data
                    .graph_terms[0]
                    .threshold_counterterm_metadata()
                    .expect("explicit GL638 variants must allocate component metadata")
                    .clone(),
                _ => unreachable!("GL638 is a cross-section fixture"),
            };
            let target_left_variant_ids = target_left_variant_ids
                .iter()
                .map(|variant_id| variant_id.0)
                .collect::<BTreeSet<_>>();
            let target_right_variant_ids = target_right_variant_ids
                .iter()
                .map(|variant_id| variant_id.0)
                .collect::<BTreeSet<_>>();
            let target_component_ids = metadata_registry
                .components
                .iter()
                .filter(|component| component.cut_group_id == Some(target_cut_group_id.0))
                .filter(|component| match component.variant_ids.as_slice() {
                    [variant_id] => {
                        target_left_variant_ids.contains(variant_id)
                            || target_right_variant_ids.contains(variant_id)
                    }
                    [left_variant_id, right_variant_id] => {
                        target_left_variant_ids.contains(left_variant_id)
                            && target_right_variant_ids.contains(right_variant_id)
                    }
                    _ => false,
                })
                .map(|component| component.component_id)
                .collect::<BTreeSet<_>>();
            assert_eq!(target_component_ids.len(), 24);
            assert!(metadata_registry.components.len() > target_component_ids.len());
            assert_eq!(
                target_component_ids
                    .iter()
                    .map(|component_id| metadata_registry.components[*component_id].kind)
                    .collect::<BTreeSet<_>>(),
                BTreeSet::from([
                    ThresholdCountertermComponentKind::Local,
                    ThresholdCountertermComponentKind::Integrated,
                    ThresholdCountertermComponentKind::LocalLocal,
                    ThresholdCountertermComponentKind::LocalIntegrated,
                    ThresholdCountertermComponentKind::IntegratedLocal,
                    ThresholdCountertermComponentKind::IntegratedIntegrated,
                ]),
            );
            integrand.warm_up(&model).unwrap();
            let top_mass = model.get_parameter("MT").value.as_ref().unwrap().re.0;
            let ProcessIntegrand::CrossSection(cross_section_integrand) = &mut *integrand else {
                unreachable!("GL638 is a cross-section fixture")
            };
            let correlated_fits = gl638_correlated_approach_fits(
                cross_section_integrand,
                &model,
                top_mass,
            );
            for branch in ["negative", "positive"] {
                let fits = &correlated_fits[branch];
                for quantity in ["q13", "eta_5_10", "eta_5_12_13"] {
                    let (slope, r_squared) = fits[quantity];
                    assert!(
                        (slope - 1.0).abs() < 0.002 && r_squared > 0.9999,
                        "GL638 {branch} {quantity} does not vanish linearly: {fits:?}",
                    );
                }
                for quantity in ["original", "counterterms", "total"] {
                    let (slope, r_squared) = fits[quantity];
                    assert!(
                        (slope + 1.0).abs() < 0.003 && r_squared > 0.99999,
                        "GL638 {branch} cured {quantity} has the wrong correlated IR scaling: {fits:?}",
                    );
                }
            }
            let before_save = integrand
                .evaluate_samples_raw(
                    &model,
                    &samples,
                    0,
                    false,
                    false,
                    Complex::new(F(0.0), F(0.0)),
                )
                .unwrap()
                .samples;
            assert!(before_save.iter().all(|result| {
                !result.evaluation_metadata.is_nan
                    && result.integrand_result.re.0.is_finite()
                    && result.integrand_result.im.0.is_finite()
            }));
            let mut observed_component_ids = BTreeSet::new();
            let mut observed_target_component_kinds = BTreeSet::new();
            let mut event_count = 0;
            for result in &before_save {
                for event_group in result.event_groups.iter() {
                    for event in event_group.iter() {
                        event_count += 1;
                        let decomposition = event
                            .additional_weights
                            .threshold_counterterms
                            .as_ref()
                            .expect("explicit threshold variants must record event decomposition");
                        assert_eq!(event.weight, decomposition.total());
                        for component in &decomposition.components {
                            let metadata = metadata_registry
                                .components
                                .get(component.component_id)
                                .expect("event component ID must resolve in the graph registry");
                            assert_eq!(metadata.component_id, component.component_id);
                            assert_eq!(
                                metadata.variant_ids.len(),
                                component.multiplier_values.len(),
                            );
                            assert!(
                                component
                                    .multiplier_values
                                    .iter()
                                    .all(|value| value.0.is_finite())
                                    && component.effective_multiplier.0.is_finite(),
                                "every recorded multiplier value must be finite",
                            );
                            assert_eq!(
                                component.effective_multiplier,
                                component
                                    .multiplier_values
                                    .iter()
                                    .fold(F(1.0), |product, value| product * *value),
                            );
                            assert!(
                                component.weighted.re.0.is_finite()
                                    && component.weighted.im.0.is_finite(),
                                "every recorded weighted component must be finite",
                            );
                            let ThresholdCountertermComponentOccurrence::LocalUnitarity {
                                overlap_groups,
                                left_threshold_order,
                                right_threshold_order,
                                lu_cut_order,
                            } = &component.occurrence
                            else {
                                panic!("a cross-section fixture must record LU occurrences")
                            };
                            assert!(lu_cut_order.is_some());
                            let expected_variant_count = metadata.kind.variant_count();
                            assert_eq!(overlap_groups.len(), expected_variant_count);
                            assert_eq!(
                                usize::from(left_threshold_order.is_some())
                                    + usize::from(right_threshold_order.is_some()),
                                expected_variant_count,
                            );
                            if component.evaluation_skipped {
                                assert!(component.bare.is_none());
                                assert_eq!(component.effective_multiplier, F(0.0));
                                assert_eq!(component.weighted, Complex::new(F(0.0), F(0.0)));
                            } else {
                                let bare = component.bare.as_ref().unwrap();
                                assert!(bare.re.0.is_finite() && bare.im.0.is_finite());
                            }
                            observed_component_ids.insert(component.component_id);
                            if target_component_ids.contains(&component.component_id) {
                                observed_target_component_kinds.insert(metadata.kind);
                            }
                        }
                    }
                }
            }
            assert!(event_count > 0, "the selected samples must produce events");
            assert!(
                target_component_ids.is_subset(&observed_component_ids),
                "the selected samples must exercise every intended target-cut Cartesian component: missing {:?}",
                target_component_ids
                    .difference(&observed_component_ids)
                    .collect::<Vec<_>>(),
            );
            assert_eq!(
                observed_target_component_kinds,
                BTreeSet::from([
                    ThresholdCountertermComponentKind::Local,
                    ThresholdCountertermComponentKind::Integrated,
                    ThresholdCountertermComponentKind::LocalLocal,
                    ThresholdCountertermComponentKind::LocalIntegrated,
                    ThresholdCountertermComponentKind::IntegratedLocal,
                    ThresholdCountertermComponentKind::IntegratedIntegrated,
                ]),
            );

            let save_root = TemporaryDirectory::new("gl638-ir-safe-threshold");
            cross_section.save(&save_root.0, true).unwrap();

            let mut state_bytes = Vec::new();
            State::export(&mut state_bytes).unwrap();
            let state_map = State::import(&mut Cursor::new(state_bytes), None).unwrap();
            let context = GammaLoopContextContainer {
                model: &model,
                state_map: &state_map,
            };
            let mut loaded =
                CrossSection::load(save_root.0.join("NNLO_fixture"), context).unwrap();
            let loaded_integrand = loaded
                .integrand
                .as_mut()
                .expect("saved fixture must reload its integrand");
            loaded_integrand.warm_up(&model).unwrap();
            let after_load = loaded_integrand
                .evaluate_samples_raw(
                    &model,
                    &samples,
                    0,
                    false,
                    false,
                    Complex::new(F(0.0), F(0.0)),
                )
                .unwrap()
                .samples;
            assert_eq!(before_save.len(), after_load.len());
            for (before, after) in before_save.iter().zip(&after_load) {
                assert_eq!(before.integrand_result, after.integrand_result);
                assert!(!after.evaluation_metadata.is_nan);
            }
        })
        .unwrap()
        .join()
        .unwrap();
}
