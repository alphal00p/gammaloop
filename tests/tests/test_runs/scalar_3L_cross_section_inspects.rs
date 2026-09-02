use super::utils::*;
use super::*;
use gammaloop_api::commands::evaluate_samples::{EvaluateSamplesPrecise, evaluate_sample_precise};
use gammalooprs::{
    graph::Autogen,
    integrands::evaluation::PreciseEvaluationResultOutput,
    processes::ProcessCollection,
    utils::{ArbPrec, FloatLike, symbolica_ext::DOD},
};
use ndarray::Array2;
use spenso::algebra::algebraic_traits::IsZero;
use std::time::Instant;
use symbolica::{domains::float::Real, parse};

// The scalar LU acceptance matrix always exercises the complete subtraction:
// local and integrated UV counterterms together with threshold counterterms.
// Feature-isolation diagnostics belong in separately named tests so they cannot
// silently weaken this production-route comparison.
// The current three CFF routes compare orientation-local and explicit-sum 3D UV
// with the projected local-4D construction. Proper LTD will append a fourth
// route to this same local-comparison matrix.

const DEFAULT_FINAL_STATES: &str = "{scalar_0 scalar_0, scalar_0 scalar_0 scalar_1}";
const SAMPLE_POINT: [f64; 9] = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];

#[derive(Clone, Copy)]
struct Scalar3LGraphCase {
    graph: &'static str,
    final_states: &'static str,
    numerator: NumeratorChoice,
}

impl Scalar3LGraphCase {
    const fn default(graph: &'static str) -> Self {
        Self {
            graph,
            final_states: DEFAULT_FINAL_STATES,
            numerator: NumeratorChoice::DefaultHigherEnergyProbe,
        }
    }

    const fn with_numerator(self, numerator: NumeratorChoice) -> Self {
        Self { numerator, ..self }
    }

    const fn with_squared_edge_numerator(self, edge: usize) -> Self {
        Self {
            numerator: NumeratorChoice::SquaredEdge {
                edge,
                dummy_index: 1,
            },
            ..self
        }
    }

    const fn with_quartic_edge_numerator(self, edge: usize) -> Self {
        Self {
            numerator: NumeratorChoice::QuarticEdge {
                edge,
                first_dummy_index: 1,
                second_dummy_index: 2,
            },
            ..self
        }
    }

    const fn with_final_states(self, final_states: &'static str) -> Self {
        Self {
            final_states,
            ..self
        }
    }
}

#[derive(Clone, Copy)]
enum NumeratorChoice {
    DefaultHigherEnergyProbe,
    SingleComponent {
        edge: usize,
        component: usize,
    },
    Dot {
        left_edge: usize,
        right_edge: usize,
        dummy_index: usize,
    },
    ComponentProduct {
        left_edge: usize,
        right_edge: usize,
        component: usize,
    },
    SquaredEdge {
        edge: usize,
        dummy_index: usize,
    },
    QuarticEdge {
        edge: usize,
        first_dummy_index: usize,
        second_dummy_index: usize,
    },
    QuarticEnergy {
        edge: usize,
    },
    DotPower {
        left_edge: usize,
        right_edge: usize,
        first_dummy_index: usize,
        power: usize,
    },
}

impl NumeratorChoice {
    fn owned_edge_factors(self) -> Vec<(usize, symbolica::atom::Atom)> {
        let component =
            |edge: usize, index: &str| parse!(&format!("gammalooprs::Q({edge},spenso::{index})"));
        let mink_component =
            |edge: usize, dummy_index: usize| component(edge, &format!("mink(4,{dummy_index})"));
        match self {
            Self::DefaultHigherEnergyProbe => vec![
                (1, mink_component(1, 1)),
                (2, mink_component(2, 1)),
                (5, mink_component(5, 2)),
                (6, mink_component(6, 2)),
            ],
            Self::SingleComponent {
                edge,
                component: component_index,
            } => vec![(edge, component(edge, &format!("cind({component_index})")))],
            Self::Dot {
                left_edge,
                right_edge,
                dummy_index,
            } => {
                vec![
                    (left_edge, mink_component(left_edge, dummy_index)),
                    (right_edge, mink_component(right_edge, dummy_index)),
                ]
            }
            Self::ComponentProduct {
                left_edge,
                right_edge,
                component: component_index,
            } => vec![
                (
                    left_edge,
                    component(left_edge, &format!("cind({component_index})")),
                ),
                (
                    right_edge,
                    component(right_edge, &format!("cind({component_index})")),
                ),
            ],
            Self::SquaredEdge { edge, dummy_index } => vec![
                (edge, mink_component(edge, dummy_index)),
                (edge, mink_component(edge, dummy_index)),
            ],
            Self::QuarticEdge {
                edge,
                first_dummy_index,
                second_dummy_index,
            } => [first_dummy_index, second_dummy_index]
                .into_iter()
                .flat_map(|dummy_index| {
                    [
                        (edge, mink_component(edge, dummy_index)),
                        (edge, mink_component(edge, dummy_index)),
                    ]
                })
                .collect(),
            Self::QuarticEnergy { edge } => {
                (0..4).map(|_| (edge, component(edge, "cind(0)"))).collect()
            }
            Self::DotPower {
                left_edge,
                right_edge,
                first_dummy_index,
                power,
            } => (0..power)
                .flat_map(|offset| {
                    let dummy_index = first_dummy_index + offset;
                    [
                        (left_edge, mink_component(left_edge, dummy_index)),
                        (right_edge, mink_component(right_edge, dummy_index)),
                    ]
                })
                .collect(),
        }
    }

    fn label(self) -> String {
        match self {
            Self::DefaultHigherEnergyProbe | Self::Dot { .. } => {
                "quadratic_pair_numerator".to_string()
            }
            Self::SingleComponent { edge, component } => {
                format!("q{edge}_component_{component}_numerator")
            }
            Self::ComponentProduct { component, .. } => {
                format!("component_{component}_pair_numerator")
            }
            Self::SquaredEdge { edge, .. } => format!("q{edge}_squared_numerator"),
            Self::QuarticEdge { edge, .. } => format!("q{edge}_quartic_numerator"),
            Self::QuarticEnergy { edge } => format!("q{edge}_energy_quartic_numerator"),
            Self::DotPower { power, .. } => format!("dot_power_{power}_numerator"),
        }
    }
}

#[derive(Clone, Copy)]
enum ScalarLocalUvRoute {
    OrientationLocal3d,
    OrientationLocal3dParametric,
    Explicit3d,
    Projected4d,
}

fn setup_scalar_3l_cross_section_cli(
    test_name: &str,
    local_uv_route: ScalarLocalUvRoute,
    graph_commands: &[&str],
    owned_numerators: &[(&str, &str, NumeratorChoice)],
    integrand_commands: &[&str],
    lmb_multichanneling: bool,
) -> Result<gammaloop_integration_tests::CLIState> {
    let diagnostic_bool = |name: &str, default: bool| {
        std::env::var(name)
            .ok()
            .map(|value| value == "true")
            .unwrap_or(default)
    };
    let subtract_uv = diagnostic_bool("GL_SCALAR_DIAGNOSTIC_SUBTRACT_UV", true);
    let generate_integrated = diagnostic_bool("GL_SCALAR_DIAGNOSTIC_INTEGRATED_UV", true);
    let enable_thresholds = diagnostic_bool("GL_SCALAR_DIAGNOSTIC_THRESHOLDS", true);
    let use_summed_evaluator = diagnostic_bool("GL_SCALAR_DIAGNOSTIC_SUMMED_EVALUATOR", true);
    let (explicit_orientation_sum_only, local_uv_from_expanded_4d, summed_evaluator) =
        match local_uv_route {
            ScalarLocalUvRoute::OrientationLocal3d => (false, false, use_summed_evaluator),
            ScalarLocalUvRoute::OrientationLocal3dParametric => (false, false, false),
            ScalarLocalUvRoute::Explicit3d => (true, false, false),
            ScalarLocalUvRoute::Projected4d => (true, true, false),
        };
    let evaluator_method = if summed_evaluator {
        "Summed"
    } else {
        "SingleParametric"
    };
    let mut cli = get_test_cli(
        None,
        get_tests_workspace_path().join(test_name),
        Some(test_name.to_string()),
        true,
    )?;

    run_commands(
        &mut cli,
        &[
            "import model scalars-default.json",
            "remove processes",
            &format!(
                "set global kv global.generation.explicit_orientation_sum_only={explicit_orientation_sum_only} global.generation.evaluator.compile=false global.generation.evaluator.summed={summed_evaluator} global.generation.uv.subtract_uv={subtract_uv} global.generation.uv.generate_integrated={generate_integrated} global.generation.uv.local_uv_cts_from_expanded_4d_integrands={local_uv_from_expanded_4d} global.generation.threshold_subtraction.enable_thresholds={enable_thresholds} global.generation.threshold_subtraction.check_esurface_at_generation=false"
            ),
            &format!(
                r#"set default-runtime string '
[general]
evaluator_method = "{evaluator_method}"
generate_events = true
store_additional_weights_in_event = true
mu_r = 3.0
m_uv = 20.0

[subtraction]
disable_threshold_subtraction = {}

[sampling]
graphs = "summed"
orientations = "summed"
lmb_multichanneling = {lmb_multichanneling}
lmb_channel_weight = "ose"
lmb_channels = "summed"
coordinate_system = "spherical"
mapping = "linear"
b = 1.0

[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
    [1.0, 0.0, 0.0, 0.0]
]
helicities = [0]
'"#,
                !enable_thresholds
            ),
        ],
    )?;
    run_commands(&mut cli, graph_commands)?;
    // Keep every momentum-dependent probe on the graph element that owns it.
    // Splitting a dot product across its two edges preserves the shared Lorentz
    // contraction while allowing each UV subgraph to retain only its own factor.
    for (process_name, integrand_name, numerator) in owned_numerators {
        let (process_id, canonical_integrand_name) = cli.state.find_integrand_ref(
            Some(&ProcessRef::Unqualified((*process_name).to_string())),
            Some(&(*integrand_name).to_string()),
        )?;
        let process = &mut cli.state.process_list.processes[process_id];
        let ProcessCollection::CrossSections(cross_sections) = &mut process.collection else {
            return Err(eyre::eyre!(
                "scalar probe process '{process_name}' is not a cross section"
            ));
        };
        let cross_section = cross_sections
            .get_mut(&canonical_integrand_name)
            .expect("resolved scalar probe integrand must exist");
        let [supergraph] = cross_section.supergraphs.as_mut_slice() else {
            return Err(eyre::eyre!(
                "scalar probe process '{process_name}' must select exactly one supergraph, found {}",
                cross_section.supergraphs.len()
            ));
        };
        for (requested_edge_id, factor) in numerator.owned_edge_factors() {
            let edge_id = supergraph
                .graph
                .underlying
                .iter_edges()
                .find_map(|(_, edge_id, _)| (edge_id.0 == requested_edge_id).then_some(edge_id))
                .ok_or_else(|| {
                    eyre::eyre!(
                        "scalar probe process '{process_name}' has no edge {requested_edge_id}"
                    )
                })?;
            let edge = &mut supergraph.graph.underlying[edge_id];
            let edge_numerator = &edge.num.value * factor;
            let edge_dod = edge_numerator.edge_dod(edge_id) - 2;
            edge.num = Autogen::explicit(edge_numerator);
            edge.dod = Autogen::explicit(edge_dod);
        }
        if std::env::var_os("GL_SCALAR_DIAGNOSTIC_PRINT_GRAPH").is_some() {
            eprintln!(
                "SCALAR GRAPH {process_name}: {}",
                supergraph.graph.debug_dot()
            );
        }
    }
    run_commands(&mut cli, &["set model mass_scalar_1=1.0"])?;
    run_commands(&mut cli, integrand_commands)?;
    run_commands(&mut cli, &["set model mass_scalar_1=0.1"])?;

    Ok(cli)
}

fn generate_graph_command(
    process: &str,
    integrand: &str,
    graph: &str,
    final_states: &str,
) -> String {
    format!(
        "generate xs scalar_1 > {final_states} | scalar_0 scalar_1 [{{{{3}}}}] --allowed-vertex-interactions V_3_SCALAR_001 V_3_SCALAR_000 -p {process} -i {integrand} -o --select-graphs {graph}"
    )
}

fn run_scalar_3l_cross_section_case(case: Scalar3LGraphCase) -> Result<()> {
    run_scalar_3l_cross_section_case_impl("scalar_3l_all", case, true, true)
}

fn run_scalar_3l_cross_section_numerator_only_case(
    test_scope: &str,
    case: Scalar3LGraphCase,
) -> Result<()> {
    run_scalar_3l_cross_section_case_impl(test_scope, case, false, false)
}

fn run_scalar_3l_cross_section_case_impl(
    test_scope: &str,
    case: Scalar3LGraphCase,
    include_no_numerator: bool,
    exercise_orientation_local_3d: bool,
) -> Result<()> {
    let diagnostic_only_orientation_local_numerator =
        std::env::var_os("GL_SCALAR_DIAGNOSTIC_ONLY_ORIENTATION_LOCAL_NUMERATOR").is_some();
    let include_no_numerator = include_no_numerator && !diagnostic_only_orientation_local_numerator;
    let numerator_label = case.numerator.label();
    let graph_label = case.graph.to_ascii_lowercase();
    let no_numerator_process = format!("{test_scope}_{graph_label}_no_num");
    let numerator_process = format!("{test_scope}_{graph_label}_{numerator_label}");

    let mut graph_commands = Vec::new();
    let mut integrand_commands = Vec::new();
    let mut evaluations = Vec::new();

    if include_no_numerator {
        graph_commands.push(generate_graph_command(
            &no_numerator_process,
            "no_numerator",
            case.graph,
            case.final_states,
        ));
        integrand_commands.push(format!(
            "generate existing -p {no_numerator_process} -i no_numerator"
        ));
        evaluations.push((
            no_numerator_process.clone(),
            "no_numerator".to_string(),
            "no_numerator".to_string(),
        ));
    }

    graph_commands.push(generate_graph_command(
        &numerator_process,
        "numerator",
        case.graph,
        case.final_states,
    ));
    integrand_commands.push(format!(
        "generate existing -p {numerator_process} -i numerator"
    ));
    evaluations.push((
        numerator_process.clone(),
        "numerator".to_string(),
        numerator_label,
    ));

    let graph_command_refs = graph_commands
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>();
    let integrand_command_refs = integrand_commands
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>();
    let arb_points = Array2::from_shape_vec((1, SAMPLE_POINT.len()), SAMPLE_POINT.to_vec())?;
    let arb_discrete_dims = Array2::from_shape_vec((1, 0), Vec::<usize>::new())?;
    let evaluate_arb = |route: &str,
                        cli: &mut gammaloop_integration_tests::CLIState,
                        process: &str,
                        integrand: &str|
     -> Result<(
        Complex<F<ArbPrec>>,
        Vec<(usize, usize, Complex<F<ArbPrec>>)>,
    )> {
        let (process_id, resolved_integrand_name) = cli.state.find_integrand_ref(
            Some(&ProcessRef::Unqualified(process.to_string())),
            Some(&integrand.to_string()),
        )?;
        let result = evaluate_sample_precise(
            &mut cli.state,
            &EvaluateSamplesPrecise {
                process_id: Some(process_id),
                integrand_name: Some(resolved_integrand_name),
                use_arb_prec: true,
                minimal_output: false,
                return_generated_events: Some(true),
                momentum_space: false,
                points: arb_points.view(),
                integrator_weights: None,
                discrete_dims: Some(arb_discrete_dims.view()),
                graph_names: None,
                orientations: None,
            },
        )?;
        match result.sample.evaluation {
            PreciseEvaluationResultOutput::Arb(result) => {
                for (group_index, group) in result.event_groups.iter().enumerate() {
                    for event in group.iter() {
                        eprintln!(
                            "ARB EVENT route={route} group={group_index} cut={} orientation={:?} weight={:e}",
                            event.cut_info.cut_id, event.cut_info.orientation_id, event.weight,
                        );
                    }
                }
                let events = result
                    .event_groups
                    .iter()
                    .enumerate()
                    .flat_map(|(group_index, group)| {
                        group.iter().map(move |event| {
                            (group_index, event.cut_info.cut_id, event.weight.clone())
                        })
                    })
                    .collect();
                eprintln!("ARB DIAG route={route} total={:e}", result.integrand_result);
                Ok((result.integrand_result, events))
            }
            evaluation => Err(eyre::eyre!(
                "{process} precise scalar local-UV comparison requested Arb but got {:?}",
                evaluation.precision()
            )),
        }
    };

    let total_started = Instant::now();
    let owned_numerators = [(numerator_process.as_str(), "numerator", case.numerator)];
    let profile_orientation_local_3d = (include_no_numerator
        || diagnostic_only_orientation_local_numerator)
        && exercise_orientation_local_3d;
    // Keep every route in a fresh graph and evaluator state. Regenerating the
    // explicit route in the orientation-local state would retain route-local caches.
    let localized_3d_results = if exercise_orientation_local_3d {
        let generation_started = Instant::now();
        let mut localized_3d = setup_scalar_3l_cross_section_cli(
            &format!(
                "{test_scope}_{}_cff_orientation_local_3d",
                case.graph.to_ascii_lowercase()
            ),
            if profile_orientation_local_3d {
                ScalarLocalUvRoute::OrientationLocal3dParametric
            } else {
                ScalarLocalUvRoute::OrientationLocal3d
            },
            &graph_command_refs,
            &owned_numerators,
            &integrand_command_refs,
            true,
        )?;
        println!(
            "scalar {} local-UV route orientation-local local-3D: setup and generation {:?}",
            case.graph,
            generation_started.elapsed()
        );
        if diagnostic_only_orientation_local_numerator {
            clean_test(&localized_3d.cli_settings.state.folder);
            return Ok(());
        }
        let mut results = Vec::with_capacity(evaluations.len());
        for (process, integrand, label) in &evaluations {
            let evaluation_started = Instant::now();
            let result = evaluate_xspace_process_with_events(
                &mut localized_3d,
                process,
                integrand,
                &SAMPLE_POINT,
                &[],
            )?;
            println!(
                "scalar {} {label} local-UV route orientation-local local-3D: evaluation {:?}",
                case.graph,
                evaluation_started.elapsed()
            );
            let arb_started = Instant::now();
            let (arb_result, _) = evaluate_arb(
                "orientation-local-3d",
                &mut localized_3d,
                process,
                integrand,
            )?;
            println!(
                "scalar {} {label} local-UV route orientation-local local-3D: Arb evaluation {:?}",
                case.graph,
                arb_started.elapsed()
            );
            results.push((result, arb_result));
        }
        if profile_orientation_local_3d {
            let generation = &localized_3d.cli_settings.global.generation;
            assert!(!generation.explicit_orientation_sum_only);
            assert!(!generation.uv.local_uv_cts_from_expanded_4d_integrands);
            assert!(generation.uv.generate_integrated);
            assert!(generation.threshold_subtraction.enable_thresholds);

            for (process, integrand, label) in &evaluations {
                let profile_started = Instant::now();
                let analysis = Profile::UltraViolet(UltraVioletProfile {
                    process: Some(ProcessRef::Unqualified(process.clone())),
                    integrand_name: Some(integrand.clone()),
                    graph: Some(case.graph.to_string()),
                    min_scale_exponent: 4.0,
                    max_scale_exponent: 8.0,
                    n_points: 6,
                    per_orientation: true,
                    seed: Some(9320),
                    ..Default::default()
                })
                .run(&mut localized_3d.state, &localized_3d.cli_settings)?
                .unwrap_uv();
                let subsets = analysis
                    .graphs
                    .iter()
                    .flat_map(|graph| &graph.lmbs)
                    .flat_map(|lmb| &lmb.subsets)
                    .collect::<Vec<_>>();
                assert!(
                    !subsets.is_empty(),
                    "scalar {} {label} only-divergent UV profile selected no limits",
                    case.graph,
                );
                let fitted_subsets = subsets
                    .iter()
                    .filter(|subset| {
                        subset
                            .per_orientation_inspect_entries
                            .iter()
                            .flatten()
                            .any(|entry| entry.analysis.is_some())
                    })
                    .count();
                if case.graph == "GL00" && integrand == "no_numerator" {
                    assert!(
                        fitted_subsets >= 2,
                        "GL00 must nontrivially profile both a component UV limit and its disconnected union",
                    );
                }
                let orientation_failures = analysis
                    .pass_fail(-0.9)
                    .failures
                    .into_iter()
                    .filter(|failure| failure.orientation_label.is_some())
                    .collect::<Vec<_>>();
                assert!(
                    orientation_failures.is_empty(),
                    "scalar {} {label} has non-integrable active orientations: {orientation_failures:#?}",
                    case.graph,
                );
                println!(
                    "scalar {} {label} orientation-local UV profile: {} selected subsets, {} fitted subsets, {:?}",
                    case.graph,
                    subsets.len(),
                    fitted_subsets,
                    profile_started.elapsed(),
                );
                for table in analysis
                    .per_orientation_tables_per_graph(-0.9)
                    .into_iter()
                    .flatten()
                {
                    println!("{table}");
                }
            }
        }
        clean_test(&localized_3d.cli_settings.state.folder);
        Some(results)
    } else {
        None
    };

    let generation_started = Instant::now();
    let mut cff_3d = setup_scalar_3l_cross_section_cli(
        &format!(
            "{test_scope}_{}_cff_local_3d",
            case.graph.to_ascii_lowercase()
        ),
        ScalarLocalUvRoute::Explicit3d,
        &graph_command_refs,
        &owned_numerators,
        &integrand_command_refs,
        exercise_orientation_local_3d,
    )?;
    if exercise_orientation_local_3d {
        println!(
            "scalar {} local-UV route explicit-sum local-3D: setup and generation {:?}",
            case.graph,
            generation_started.elapsed()
        );
    }

    let generation_started = Instant::now();
    let mut cff_4d = setup_scalar_3l_cross_section_cli(
        &format!(
            "{test_scope}_{}_cff_local_4d",
            case.graph.to_ascii_lowercase()
        ),
        ScalarLocalUvRoute::Projected4d,
        &graph_command_refs,
        &owned_numerators,
        &integrand_command_refs,
        exercise_orientation_local_3d,
    )?;
    if exercise_orientation_local_3d {
        println!(
            "scalar {} local-UV route projected local-4D: setup and generation {:?}",
            case.graph,
            generation_started.elapsed()
        );
    }
    for (evaluation_index, (process, integrand, label)) in evaluations.into_iter().enumerate() {
        let evaluation_started = Instant::now();
        let cff_3d_result = evaluate_xspace_process_with_events(
            &mut cff_3d,
            &process,
            &integrand,
            &SAMPLE_POINT,
            &[],
        )?;
        if exercise_orientation_local_3d {
            println!(
                "scalar {} local-UV route explicit-sum local-3D: evaluation {:?}",
                case.graph,
                evaluation_started.elapsed()
            );
        }
        let evaluation_started = Instant::now();
        let cff_4d_result = evaluate_xspace_process_with_events(
            &mut cff_4d,
            &process,
            &integrand,
            &SAMPLE_POINT,
            &[],
        )?;
        if exercise_orientation_local_3d {
            println!(
                "scalar {} local-UV route projected local-4D: evaluation {:?}",
                case.graph,
                evaluation_started.elapsed()
            );
        }
        assert_evaluation_outputs_match(
            &cff_4d_result.sample.evaluation,
            &cff_3d_result.sample.evaluation,
            &format!(
                "scalar 3L cross-section {} {label} rich inspect parity: CFF local-4D vs CFF local-3D",
                case.graph
            ),
        );
        if let Some(localized_3d_results) = &localized_3d_results {
            let (localized_3d_result, localized_3d_arb) = &localized_3d_results[evaluation_index];
            let route_totals = [
                (
                    "orientation-local local-3D",
                    complex_ff64(&localized_3d_result.sample.evaluation.integrand_result),
                ),
                (
                    "explicit-sum local-3D",
                    complex_ff64(&cff_3d_result.sample.evaluation.integrand_result),
                ),
                (
                    "projected local-4D",
                    complex_ff64(&cff_4d_result.sample.evaluation.integrand_result),
                ),
            ];
            for (route, total) in route_totals {
                assert!(
                    total.re.is_finite() && total.im.is_finite() && total.re.hypot(total.im) > 0.0,
                    "scalar {} {label} {route} total must be finite and nonzero, got {total:e}",
                    case.graph
                );
            }
            let arb_started = Instant::now();
            let (explicit_3d_arb, explicit_3d_events) =
                evaluate_arb("explicit-3d", &mut cff_3d, &process, &integrand)?;
            let (projected_4d_arb, projected_4d_events) =
                evaluate_arb("projected-4d", &mut cff_4d, &process, &integrand)?;
            if integrand == "numerator" {
                for (
                    (explicit_group, explicit_cut, explicit_weight),
                    (projected_group, projected_cut, projected_weight),
                ) in explicit_3d_events.iter().zip(&projected_4d_events)
                {
                    assert_eq!(
                        (explicit_group, explicit_cut),
                        (projected_group, projected_cut)
                    );
                    let delta = projected_weight.clone() - explicit_weight.clone();
                    if !delta.norm().re.is_zero() {
                        eprintln!(
                            "ARB DIAG event group={explicit_group} cut={explicit_cut} explicit={explicit_weight:e} projected={projected_weight:e} delta={delta:e}"
                        );
                    }
                }
            }
            let precision_tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
            let f64_input_tolerance = F::<ArbPrec>::from_f64(1.0e-14);
            for (route, actual) in [
                ("orientation-local local-3D", localized_3d_arb),
                ("projected local-4D", &projected_4d_arb),
            ] {
                let distance = (actual.clone() - explicit_3d_arb.clone()).norm().re;
                let actual_norm = actual.norm().re;
                let reference_norm = explicit_3d_arb.norm().re;
                let scale = if actual_norm > reference_norm {
                    actual_norm
                } else {
                    reference_norm
                };
                let relative_distance = if scale.is_zero() {
                    distance.clone()
                } else {
                    distance.clone() / scale.clone()
                };
                // The sample is supplied as f64, so a cancellation to zero can
                // preserve route-dependent f64 roundoff after both evaluators
                // are upcast to Arb. Keep reporting the precision-scaled Arb
                // comparison, but accept that known input floor at 1e-14.
                let f64_scale = if scale > F::<ArbPrec>::from_f64(1.0) {
                    scale
                } else {
                    F::<ArbPrec>::from_f64(1.0)
                };
                let f64_scaled_distance = distance / f64_scale;
                assert!(
                    relative_distance <= precision_tolerance.clone()
                        || f64_scaled_distance <= f64_input_tolerance.clone(),
                    "scalar {} {label} {route} differs from precise explicit-sum local-3D: actual={actual:e}, reference={explicit_3d_arb:e}, precision-scaled relative delta={relative_distance:e}, precision tolerance={precision_tolerance:e}, f64-scaled delta={f64_scaled_distance:e}, f64-input tolerance={f64_input_tolerance:e}",
                    case.graph,
                );
            }
            println!(
                "scalar {} local-UV three-route Arb comparison: {:?}",
                case.graph,
                arb_started.elapsed()
            );
        }
    }

    if exercise_orientation_local_3d {
        println!(
            "scalar {} local-UV three-route acceptance: total {:?}",
            case.graph,
            total_started.elapsed()
        );
    }

    clean_test(&cff_3d.cli_settings.state.folder);
    clean_test(&cff_4d.cli_settings.state.folder);
    Ok(())
}

fn run_scalar_3l_quadratic_energy_numerator_case(
    case: Scalar3LGraphCase,
    edge: usize,
) -> Result<()> {
    let numerator_case = case.with_squared_edge_numerator(edge);
    run_scalar_3l_cross_section_numerator_only_case(
        &format!("scalar_3l_quadratic_energy_q{edge}_squared"),
        numerator_case,
    )
}

fn run_scalar_3l_gl16_quadratic_energy_three_route_case() -> Result<()> {
    let case = Scalar3LGraphCase::default("GL16").with_squared_edge_numerator(7);
    run_scalar_3l_cross_section_case_impl(
        "scalar_3l_quadratic_energy_q7_squared",
        case,
        false,
        true,
    )
}

fn run_scalar_3l_gl16_quartic_energy_three_route_case() -> Result<()> {
    let case = Scalar3LGraphCase::default("GL16").with_quartic_edge_numerator(7);
    run_scalar_3l_cross_section_case_impl("scalar_3l_quartic_energy_q7", case, false, true)
}

fn run_scalar_3l_quartic_energy_numerator_case(case: Scalar3LGraphCase, edge: usize) -> Result<()> {
    let numerator_case = case.with_quartic_edge_numerator(edge);
    run_scalar_3l_cross_section_numerator_only_case(
        &format!("scalar_3l_quartic_energy_q{edge}"),
        numerator_case,
    )
}

#[test]
#[serial]
fn scalar_3l_gl00_higher_power_cff_is_invariant_under_nonzero_sampling_scale() -> Result<()> {
    let case = Scalar3LGraphCase::default("GL00");
    let process = "scalar_3l_gl00_higher_power_sampling_scale";
    let integrand = "numerator";
    let numerator = NumeratorChoice::DotPower {
        left_edge: 1,
        right_edge: 2,
        first_dummy_index: 1,
        power: 3,
    };
    let sampling_mode =
        "set global kv global.generation.uniform_numerator_sampling_scale=beyond_quadratic"
            .to_string();
    let graph_command = generate_graph_command(process, integrand, case.graph, case.final_states);
    let integrand_command = format!("generate existing -p {process} -i {integrand}");
    let mut cli = setup_scalar_3l_cross_section_cli(
        "scalar_3l_gl00_higher_power_sampling_scale",
        ScalarLocalUvRoute::Explicit3d,
        &[&sampling_mode, &graph_command],
        &[(process, integrand, numerator)],
        &[&integrand_command],
        false,
    )?;

    cli.run_command(&format!(
        "set process -p {process} -i {integrand} kv general.numerator_sampling_scale=0.75"
    ))?;
    let first = complex_ff64(
        &evaluate_xspace_process_with_events(&mut cli, process, integrand, &SAMPLE_POINT, &[])?
            .sample
            .evaluation
            .integrand_result,
    );
    cli.run_command(&format!(
        "set process -p {process} -i {integrand} kv general.numerator_sampling_scale=2.25"
    ))?;
    let second = complex_ff64(
        &evaluate_xspace_process_with_events(&mut cli, process, integrand, &SAMPLE_POINT, &[])?
            .sample
            .evaluation
            .integrand_result,
    );
    let first_norm = (first.re * first.re + first.im * first.im).sqrt();
    let second_norm = (second.re * second.re + second.im * second.im).sqrt();
    let scale = first_norm.max(second_norm);
    assert!(scale > 1.0e-18, "sampling-scale acceptance point is zero");
    assert!(
        complex_distance(first, second) <= 1.0e-8 * scale,
        "GammaLoop higher-power CFF changed with auxiliary M: {first} vs {second}"
    );

    cli.run_command(&format!(
        "set process -p {process} -i {integrand} kv general.numerator_sampling_scale=0.0"
    ))?;
    let error =
        match evaluate_xspace_process_with_events(&mut cli, process, integrand, &SAMPLE_POINT, &[])
        {
            Ok(_) => panic!("an evaluator that uses M accepted a zero sampling scale"),
            Err(error) => error,
        };
    assert!(format!("{error:#}").contains("sampling scale M"));

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

macro_rules! scalar_3l_graph_test {
    ($name:ident, $graph:literal) => {
        #[test]
        #[serial]
        fn $name() -> Result<()> {
            run_scalar_3l_cross_section_case(Scalar3LGraphCase::default($graph))
        }
    };
    ($name:ident, $case:expr) => {
        #[test]
        #[serial]
        fn $name() -> Result<()> {
            run_scalar_3l_cross_section_case($case)
        }
    };
}

macro_rules! quadratic_energy_graph_test {
    ($name:ident, q1_only: $q1_case:expr) => {
        mod $name {
            use super::*;

            #[test]
            #[serial]
            fn q1_squared() -> Result<()> {
                run_scalar_3l_quadratic_energy_numerator_case($q1_case, 1)
            }
        }
    };
    ($name:ident, q7_only: $q7_case:expr) => {
        mod $name {
            use super::*;

            #[test]
            #[serial]
            fn q7_squared() -> Result<()> {
                run_scalar_3l_quadratic_energy_numerator_case($q7_case, 7)
            }
        }
    };
    ($name:ident, $graph:literal) => {
        mod $name {
            use super::*;

            #[test]
            #[serial]
            fn q1_squared() -> Result<()> {
                run_scalar_3l_quadratic_energy_numerator_case(Scalar3LGraphCase::default($graph), 1)
            }

            #[test]
            #[serial]
            fn q7_squared() -> Result<()> {
                run_scalar_3l_quadratic_energy_numerator_case(Scalar3LGraphCase::default($graph), 7)
            }
        }
    };
    ($name:ident, $case:expr) => {
        mod $name {
            use super::*;

            #[test]
            #[serial]
            fn q1_squared() -> Result<()> {
                run_scalar_3l_quadratic_energy_numerator_case($case, 1)
            }

            #[test]
            #[serial]
            fn q7_squared() -> Result<()> {
                run_scalar_3l_quadratic_energy_numerator_case($case, 7)
            }
        }
    };
    ($name:ident, q1: $q1_case:expr, q7: $q7_case:expr) => {
        mod $name {
            use super::*;

            #[test]
            #[serial]
            fn q1_squared() -> Result<()> {
                run_scalar_3l_quadratic_energy_numerator_case($q1_case, 1)
            }

            #[test]
            #[serial]
            fn q7_squared() -> Result<()> {
                run_scalar_3l_quadratic_energy_numerator_case($q7_case, 7)
            }
        }
    };
}

macro_rules! quartic_energy_graph_test {
    ($name:ident, $graph:literal) => {
        mod $name {
            use super::*;

            #[test]
            #[serial]
            fn q1_quartic() -> Result<()> {
                run_scalar_3l_quartic_energy_numerator_case(Scalar3LGraphCase::default($graph), 1)
            }
        }
    };
    ($name:ident, $case:expr) => {
        mod $name {
            use super::*;

            #[test]
            #[serial]
            fn q1_quartic() -> Result<()> {
                run_scalar_3l_quartic_energy_numerator_case($case, 1)
            }
        }
    };
}

mod default_scalar_3l_cross_section_inspects {
    use super::*;

    scalar_3l_graph_test!(scalar_3l_cross_section_gl00_inspects_match, "GL00");
    scalar_3l_graph_test!(scalar_3l_cross_section_gl02_inspects_match, "GL02");
    scalar_3l_graph_test!(scalar_3l_cross_section_gl04_inspects_match, "GL04");
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl08_inspects_match,
        Scalar3LGraphCase::default("GL08")
    );
    scalar_3l_graph_test!(scalar_3l_cross_section_gl09_inspects_match, "GL09");
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl24_inspects_match,
        Scalar3LGraphCase::default("GL24").with_numerator(NumeratorChoice::Dot {
            left_edge: 1,
            right_edge: 2,
            dummy_index: 1,
        })
    );

    #[test]
    #[serial]
    fn scalar_3l_cross_section_gl00_first_owned_dot_inspects_match() -> Result<()> {
        run_scalar_3l_cross_section_case_impl(
            "scalar_3l_gl00_first_owned_dot",
            Scalar3LGraphCase::default("GL00").with_numerator(NumeratorChoice::Dot {
                left_edge: 1,
                right_edge: 2,
                dummy_index: 1,
            }),
            false,
            true,
        )
    }

    #[test]
    #[serial]
    fn scalar_3l_cross_section_gl00_second_owned_dot_inspects_match() -> Result<()> {
        run_scalar_3l_cross_section_case_impl(
            "scalar_3l_gl00_second_owned_dot",
            Scalar3LGraphCase::default("GL00").with_numerator(NumeratorChoice::Dot {
                left_edge: 5,
                right_edge: 6,
                dummy_index: 1,
            }),
            false,
            true,
        )
    }

    #[test]
    #[serial]
    fn scalar_3l_cross_section_gl00_second_owned_temporal_product_inspects_match() -> Result<()> {
        run_scalar_3l_cross_section_case_impl(
            "scalar_3l_gl00_second_owned_temporal_product",
            Scalar3LGraphCase::default("GL00").with_numerator(NumeratorChoice::ComponentProduct {
                left_edge: 5,
                right_edge: 6,
                component: 0,
            }),
            false,
            true,
        )
    }

    #[test]
    #[serial]
    fn scalar_3l_cross_section_gl00_q5_temporal_component_inspects_match() -> Result<()> {
        run_scalar_3l_cross_section_case_impl(
            "scalar_3l_gl00_q5_temporal_component",
            Scalar3LGraphCase::default("GL00").with_numerator(NumeratorChoice::SingleComponent {
                edge: 5,
                component: 0,
            }),
            false,
            true,
        )
    }

    #[test]
    #[serial]
    fn scalar_3l_cross_section_gl00_q6_temporal_component_inspects_match() -> Result<()> {
        run_scalar_3l_cross_section_case_impl(
            "scalar_3l_gl00_q6_temporal_component",
            Scalar3LGraphCase::default("GL00").with_numerator(NumeratorChoice::SingleComponent {
                edge: 6,
                component: 0,
            }),
            false,
            true,
        )
    }

    #[test]
    #[serial]
    fn scalar_3l_cross_section_gl02_dod_zero_triangle_inspects_match() -> Result<()> {
        run_scalar_3l_cross_section_case_impl(
            "scalar_3l_gl02_dod_zero_triangle",
            Scalar3LGraphCase::default("GL02"),
            false,
            true,
        )
    }

    mod quadratic_energy_numerators {
        use super::*;

        mod scalar_3l_cross_section_gl00_quadratic_energy_inspects_match {
            use super::*;

            #[test]
            #[serial]
            fn q7_squared() -> Result<()> {
                run_scalar_3l_cross_section_case_impl(
                    "scalar_3l_quadratic_energy_q7_squared",
                    Scalar3LGraphCase::default("GL00").with_squared_edge_numerator(7),
                    false,
                    true,
                )
            }
        }
        mod scalar_3l_cross_section_gl16_quadratic_energy_inspects_match {
            use super::*;

            #[test]
            #[serial]
            fn q7_squared() -> Result<()> {
                run_scalar_3l_gl16_quadratic_energy_three_route_case()
            }
        }
    }

    mod quartic_energy_numerators {
        use super::*;

        mod scalar_3l_cross_section_gl02_quartic_energy_inspects_match {
            use super::*;

            #[test]
            #[serial]
            fn q1_quartic() -> Result<()> {
                run_scalar_3l_cross_section_case_impl(
                    "scalar_3l_quartic_energy_q1",
                    Scalar3LGraphCase::default("GL02")
                        .with_numerator(NumeratorChoice::QuarticEnergy { edge: 1 }),
                    false,
                    true,
                )
            }
        }
    }
}

mod slow {
    use super::*;

    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl01_inspects_match,
        Scalar3LGraphCase::default("GL01")
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl03_inspects_match,
        Scalar3LGraphCase::default("GL03").with_numerator(NumeratorChoice::Dot {
            left_edge: 7,
            right_edge: 8,
            dummy_index: 1,
        })
    );
    scalar_3l_graph_test!(scalar_3l_cross_section_gl05_inspects_match, "GL05");
    scalar_3l_graph_test!(scalar_3l_cross_section_gl06_inspects_match, "GL06");
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl07_inspects_match,
        Scalar3LGraphCase::default("GL07")
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl10_inspects_match,
        Scalar3LGraphCase::default("GL10")
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl11_inspects_match,
        Scalar3LGraphCase::default("GL11")
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl12_inspects_match,
        Scalar3LGraphCase::default("GL12")
    );
    scalar_3l_graph_test!(scalar_3l_cross_section_gl13_inspects_match, "GL13");
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl14_inspects_match,
        Scalar3LGraphCase::default("GL14")
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl15_inspects_match,
        Scalar3LGraphCase::default("GL15")
    );
    scalar_3l_graph_test!(scalar_3l_cross_section_gl16_inspects_match, "GL16");
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl17_inspects_match,
        Scalar3LGraphCase::default("GL17")
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl18_inspects_match,
        Scalar3LGraphCase::default("GL18")
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl19_inspects_match,
        Scalar3LGraphCase::default("GL19")
    );
    scalar_3l_graph_test!(scalar_3l_cross_section_gl20_inspects_match, "GL20");
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl21_inspects_match,
        Scalar3LGraphCase::default("GL21")
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl22_inspects_match,
        Scalar3LGraphCase::default("GL22").with_numerator(NumeratorChoice::Dot {
            left_edge: 1,
            right_edge: 2,
            dummy_index: 1,
        })
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl23_inspects_match,
        Scalar3LGraphCase::default("GL23").with_numerator(NumeratorChoice::Dot {
            left_edge: 1,
            right_edge: 2,
            dummy_index: 1,
        })
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl25_inspects_match,
        Scalar3LGraphCase::default("GL25").with_numerator(NumeratorChoice::Dot {
            left_edge: 1,
            right_edge: 2,
            dummy_index: 1,
        })
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl26_inspects_match,
        Scalar3LGraphCase::default("GL26")
    );
    scalar_3l_graph_test!(scalar_3l_cross_section_gl27_inspects_match, "GL27");
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl28_inspects_match,
        Scalar3LGraphCase::default("GL28").with_numerator(NumeratorChoice::Dot {
            left_edge: 1,
            right_edge: 2,
            dummy_index: 1,
        })
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl29_inspects_match,
        Scalar3LGraphCase::default("GL29")
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl30_inspects_match,
        Scalar3LGraphCase::default("GL30").with_numerator(NumeratorChoice::Dot {
            left_edge: 1,
            right_edge: 2,
            dummy_index: 1,
        })
    );
    scalar_3l_graph_test!(scalar_3l_cross_section_gl31_inspects_match, "GL31");
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl32_inspects_match,
        Scalar3LGraphCase::default("GL32").with_numerator(NumeratorChoice::Dot {
            left_edge: 1,
            right_edge: 2,
            dummy_index: 1,
        })
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl33_inspects_match,
        Scalar3LGraphCase::default("GL33")
    );
    scalar_3l_graph_test!(scalar_3l_cross_section_gl34_inspects_match, "GL34");
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl35_inspects_match,
        Scalar3LGraphCase::default("GL35")
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl36_inspects_match,
        Scalar3LGraphCase::default("GL36")
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl37_inspects_match,
        Scalar3LGraphCase::default("GL37")
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl38_inspects_match,
        Scalar3LGraphCase::default("GL38")
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl39_inspects_match,
        Scalar3LGraphCase::default("GL39").with_numerator(NumeratorChoice::Dot {
            left_edge: 1,
            right_edge: 2,
            dummy_index: 1,
        })
    );
    scalar_3l_graph_test!(scalar_3l_cross_section_gl40_inspects_match, "GL40");
    scalar_3l_graph_test!(scalar_3l_cross_section_gl41_inspects_match, "GL41");
    scalar_3l_graph_test!(scalar_3l_cross_section_gl42_inspects_match, "GL42");
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl43_inspects_match,
        Scalar3LGraphCase::default("GL43").with_numerator(NumeratorChoice::Dot {
            left_edge: 1,
            right_edge: 2,
            dummy_index: 1,
        })
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl44_inspects_match,
        Scalar3LGraphCase::default("GL44")
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl45_inspects_match,
        Scalar3LGraphCase::default("GL45")
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl46_inspects_match,
        Scalar3LGraphCase::default("GL46")
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl47_inspects_match,
        Scalar3LGraphCase::default("GL47").with_numerator(NumeratorChoice::Dot {
            left_edge: 1,
            right_edge: 2,
            dummy_index: 1,
        })
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl48_inspects_match,
        Scalar3LGraphCase::default("GL48")
            .with_final_states("scalar_0 scalar_0")
            .with_numerator(NumeratorChoice::Dot {
                left_edge: 7,
                right_edge: 8,
                dummy_index: 1,
            })
    );

    mod quadratic_energy_numerators {
        use super::*;

        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl00_quadratic_energy_inspects_match,
            q1_only: Scalar3LGraphCase::default("GL00")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl01_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL01")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl02_quadratic_energy_inspects_match,
            q1_only: Scalar3LGraphCase::default("GL02")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl03_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL03")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl04_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL04")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl05_quadratic_energy_inspects_match,
            q1_only: Scalar3LGraphCase::default("GL05")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl06_quadratic_energy_inspects_match,
            "GL06"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl07_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL07")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl08_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL08")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl09_quadratic_energy_inspects_match,
            "GL09"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl10_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL10")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl11_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL11")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl12_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL12")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl13_quadratic_energy_inspects_match,
            "GL13"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl14_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL14")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl15_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL15")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl16_quadratic_energy_inspects_match,
            q7_only: Scalar3LGraphCase::default("GL16")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl17_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL17")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl18_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL18")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl19_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL19")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl20_quadratic_energy_inspects_match,
            "GL20"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl21_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL21")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl22_quadratic_energy_inspects_match,
            "GL22"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl23_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL23")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl24_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL24")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl25_quadratic_energy_inspects_match,
            "GL25"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl26_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL26")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl27_quadratic_energy_inspects_match,
            q1_only: Scalar3LGraphCase::default("GL27")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl28_quadratic_energy_inspects_match,
            "GL28"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl29_quadratic_energy_inspects_match,
            q1_only: Scalar3LGraphCase::default("GL29")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl30_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL30")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl31_quadratic_energy_inspects_match,
            q1_only: Scalar3LGraphCase::default("GL31")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl32_quadratic_energy_inspects_match,
            "GL32"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl33_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL33")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl34_quadratic_energy_inspects_match,
            q1_only: Scalar3LGraphCase::default("GL34")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl35_quadratic_energy_inspects_match,
            q1_only: Scalar3LGraphCase::default("GL35")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl36_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL36")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl37_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL37")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl38_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL38")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl39_quadratic_energy_inspects_match,
            "GL39"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl40_quadratic_energy_inspects_match,
            "GL40"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl41_quadratic_energy_inspects_match,
            "GL41"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl42_quadratic_energy_inspects_match,
            q1_only: Scalar3LGraphCase::default("GL42")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl43_quadratic_energy_inspects_match,
            "GL43"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl44_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL44")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl45_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL45")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl46_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL46")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl47_quadratic_energy_inspects_match,
            "GL47"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl48_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL48").with_final_states("scalar_0 scalar_0")
        );
    }

    mod quartic_energy_numerators {
        use super::*;

        // This all-routes rank-four stress remains valuable, but its first
        // orientation-local generation alone can exceed ten minutes. The base
        // suite covers the same retained-cut component bridge with a fast
        // exact-residue oracle and keeps GL16 rank two in all three routes.
        mod scalar_3l_cross_section_gl16_quartic_energy_inspects_match {
            use super::*;

            #[test]
            #[serial]
            fn q7_quartic() -> Result<()> {
                run_scalar_3l_gl16_quartic_energy_three_route_case()
            }
        }

        mod scalar_3l_cross_section_gl24_quartic_energy_inspects_match {
            use super::*;

            #[test]
            #[serial]
            fn q1_quartic() -> Result<()> {
                run_scalar_3l_cross_section_case_impl(
                    "scalar_3l_quartic_energy_q1",
                    Scalar3LGraphCase::default("GL24").with_quartic_edge_numerator(1),
                    false,
                    true,
                )
            }
        }

        quartic_energy_graph_test!(
            scalar_3l_cross_section_gl00_quartic_energy_inspects_match,
            Scalar3LGraphCase::default("GL00")
        );
        quartic_energy_graph_test!(
            scalar_3l_cross_section_gl06_quartic_energy_inspects_match,
            "GL06"
        );
        quartic_energy_graph_test!(
            scalar_3l_cross_section_gl10_quartic_energy_inspects_match,
            Scalar3LGraphCase::default("GL10")
        );
        quartic_energy_graph_test!(
            scalar_3l_cross_section_gl20_quartic_energy_inspects_match,
            "GL20"
        );
        quartic_energy_graph_test!(
            scalar_3l_cross_section_gl30_quartic_energy_inspects_match,
            Scalar3LGraphCase::default("GL30")
        );
        quartic_energy_graph_test!(
            scalar_3l_cross_section_gl41_quartic_energy_inspects_match,
            "GL41"
        );
        quartic_energy_graph_test!(
            scalar_3l_cross_section_gl42_quartic_energy_inspects_match,
            "GL42"
        );
    }
}
