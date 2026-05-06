use super::utils::*;
use super::*;

// These flags keep the exhaustive 3L scalar cross-section coverage alive while
// raised threshold counterterms and disconnected UV unions remain known broader
// limitations. Turning either flag to false promotes the corresponding fallback
// to a hard test failure.
const ALLOW_DISABLING_UV_SUBTRACTION: bool = true;
const ALLOW_DISABLING_THRESHOLD_SUBTRACTION: bool = true;
// Toggle LTD into the exhaustive slow sweep. With this true, each case compares
// CFF local UV from 3D expansions, CFF local UV from expanded 4D integrands, and
// LTD local UV from expanded 4D integrands before snapshotting the agreed CFF
// result.
const CONSIDER_COMPARISON_WITH_LTD: bool = true;

const DEFAULT_FINAL_STATES: &str = "{scalar_0 scalar_0, scalar_0 scalar_0 scalar_1}";
const SAMPLE_POINT: [f64; 9] = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];

#[derive(Clone, Copy)]
struct Scalar3LGraphCase {
    graph: &'static str,
    final_states: &'static str,
    numerator: NumeratorChoice,
    disable_uv_subtraction: bool,
    disable_threshold_subtraction: bool,
}

impl Scalar3LGraphCase {
    const fn default(graph: &'static str) -> Self {
        Self {
            graph,
            final_states: DEFAULT_FINAL_STATES,
            numerator: NumeratorChoice::DefaultHigherEnergyProbe,
            disable_uv_subtraction: false,
            disable_threshold_subtraction: false,
        }
    }

    const fn without_uv_subtraction(self) -> Self {
        Self {
            disable_uv_subtraction: true,
            ..self
        }
    }

    const fn without_threshold_subtraction(self) -> Self {
        Self {
            disable_threshold_subtraction: true,
            ..self
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

    fn subtract_uv(self) -> bool {
        if self.disable_uv_subtraction {
            assert!(
                ALLOW_DISABLING_UV_SUBTRACTION,
                "{} requests UV-subtraction fallback, but ALLOW_DISABLING_UV_SUBTRACTION is false",
                self.graph
            );
            false
        } else {
            true
        }
    }

    fn enable_thresholds(self) -> bool {
        if self.disable_threshold_subtraction {
            assert!(
                ALLOW_DISABLING_THRESHOLD_SUBTRACTION,
                "{} requests threshold-subtraction fallback, but ALLOW_DISABLING_THRESHOLD_SUBTRACTION is false",
                self.graph
            );
            false
        } else {
            true
        }
    }
}

#[derive(Clone, Copy)]
enum NumeratorChoice {
    DefaultHigherEnergyProbe,
    Dot {
        left_edge: usize,
        right_edge: usize,
        dummy_index: usize,
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
}

impl NumeratorChoice {
    fn expression(self) -> String {
        match self {
            Self::DefaultHigherEnergyProbe => [mink_dot(1, 2, 1), mink_dot(5, 6, 2)].join("*"),
            Self::Dot {
                left_edge,
                right_edge,
                dummy_index,
            } => mink_dot(left_edge, right_edge, dummy_index),
            Self::SquaredEdge { edge, dummy_index } => mink_dot(edge, edge, dummy_index),
            Self::QuarticEdge {
                edge,
                first_dummy_index,
                second_dummy_index,
            } => [
                mink_dot(edge, edge, first_dummy_index),
                mink_dot(edge, edge, second_dummy_index),
            ]
            .join("*"),
        }
    }

    fn label(self) -> String {
        match self {
            Self::DefaultHigherEnergyProbe | Self::Dot { .. } => {
                "quadratic_pair_numerator".to_string()
            }
            Self::SquaredEdge { edge, .. } => format!("q{edge}_squared_numerator"),
            Self::QuarticEdge { edge, .. } => format!("q{edge}_quartic_numerator"),
        }
    }
}

fn setup_scalar_3l_cross_section_cli(
    test_name: &str,
    representation: &str,
    local_uv_from_expanded_4d: bool,
    graph_commands: &[&str],
    integrand_commands: &[&str],
    subtract_uv: bool,
    enable_thresholds: bool,
) -> Result<gammaloop_integration_tests::CLIState> {
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
                "set global kv global.3d_representation={representation} global.generation.explicit_orientation_sum_only=true global.generation.evaluator.compile=false global.generation.uv.subtract_uv={subtract_uv} global.generation.uv.generate_integrated=false global.generation.uv.local_uv_cts_from_expanded_4d_integrands={local_uv_from_expanded_4d} global.generation.threshold_subtraction.enable_thresholds={enable_thresholds} global.generation.threshold_subtraction.check_esurface_at_generation=false"
            ),
            &format!(
                r#"set default-runtime string '
[general]
evaluator_method = "Summed"
generate_events = true
store_additional_weights_in_event = true
mu_r = 3.0
m_uv = 20.0

[subtraction]
disable_threshold_subtraction = {}

[sampling]
graphs = "summed"
orientations = "summed"
lmb_multichanneling = false
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
    numerator: Option<&str>,
) -> String {
    let numerator = numerator
        .map(|expr| format!(" --global-prefactor-num '{expr}'"))
        .unwrap_or_default();
    format!(
        "generate xs scalar_1 > {final_states} | scalar_0 scalar_1 [{{{{3}}}}] --allowed-vertex-interactions V_3_SCALAR_001 V_3_SCALAR_000 -p {process} -i {integrand} -o --select-graphs {graph}{numerator}"
    )
}

fn assert_weight_benchmark(
    case: Scalar3LGraphCase,
    integrand_label: &str,
    result: &gammalooprs::integrands::evaluation::SingleSampleEvaluationResult,
) {
    let weight = complex_ff64(&result.sample.evaluation.integrand_result);
    let subtraction_label = format!(
        "uv_{}_threshold_{}",
        if case.subtract_uv() { "on" } else { "off" },
        if case.enable_thresholds() {
            "on"
        } else {
            "off"
        }
    );
    insta::assert_snapshot!(
        format!(
            "scalar_3l_cross_section_{}_{}_{}",
            case.graph.to_ascii_lowercase(),
            integrand_label,
            subtraction_label
        ),
        format!("re = {:.17e}\nim = {:.17e}", weight.re, weight.im)
    );
}

fn run_scalar_3l_cross_section_case(case: Scalar3LGraphCase) -> Result<()> {
    run_scalar_3l_cross_section_case_impl("scalar_3l_all", case, true)
}

fn run_scalar_3l_cross_section_numerator_only_case(
    test_scope: &str,
    case: Scalar3LGraphCase,
) -> Result<()> {
    run_scalar_3l_cross_section_case_impl(test_scope, case, false)
}

fn run_scalar_3l_cross_section_case_impl(
    test_scope: &str,
    case: Scalar3LGraphCase,
    include_no_numerator: bool,
) -> Result<()> {
    let numerator = case.numerator.expression();
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
            None,
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
        Some(&numerator),
    ));
    integrand_commands.push(format!(
        "generate existing -p {numerator_process} -i numerator"
    ));
    evaluations.push((numerator_process, "numerator".to_string(), numerator_label));

    let graph_command_refs = graph_commands
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>();
    let integrand_command_refs = integrand_commands
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>();

    let subtract_uv = case.subtract_uv();
    let enable_thresholds = case.enable_thresholds();
    let mut cff_3d = setup_scalar_3l_cross_section_cli(
        &format!(
            "{test_scope}_{}_cff_local_3d",
            case.graph.to_ascii_lowercase()
        ),
        "CFF",
        false,
        &graph_command_refs,
        &integrand_command_refs,
        subtract_uv,
        enable_thresholds,
    )?;
    let mut cff_4d = setup_scalar_3l_cross_section_cli(
        &format!(
            "{test_scope}_{}_cff_local_4d",
            case.graph.to_ascii_lowercase()
        ),
        "CFF",
        true,
        &graph_command_refs,
        &integrand_command_refs,
        subtract_uv,
        enable_thresholds,
    )?;
    let mut ltd_4d = if CONSIDER_COMPARISON_WITH_LTD {
        Some(setup_scalar_3l_cross_section_cli(
            &format!(
                "{test_scope}_{}_ltd_local_4d",
                case.graph.to_ascii_lowercase()
            ),
            "LTD",
            true,
            &graph_command_refs,
            &integrand_command_refs,
            subtract_uv,
            enable_thresholds,
        )?)
    } else {
        None
    };

    for (process, integrand, label) in evaluations {
        let cff_3d_result = evaluate_xspace_process_with_events(
            &mut cff_3d,
            &process,
            &integrand,
            &SAMPLE_POINT,
            &[],
        )?;
        let cff_4d_result = evaluate_xspace_process_with_events(
            &mut cff_4d,
            &process,
            &integrand,
            &SAMPLE_POINT,
            &[],
        )?;
        assert_evaluation_outputs_match(
            &cff_4d_result.sample.evaluation,
            &cff_3d_result.sample.evaluation,
            &format!(
                "scalar 3L cross-section {} {label} rich inspect parity: CFF local-4D vs CFF local-3D",
                case.graph
            ),
        );
        if let Some(ltd_4d) = ltd_4d.as_mut() {
            let ltd_4d_result = evaluate_xspace_process_with_events(
                ltd_4d,
                &process,
                &integrand,
                &SAMPLE_POINT,
                &[],
            )?;
            assert_evaluation_outputs_match(
                &ltd_4d_result.sample.evaluation,
                &cff_3d_result.sample.evaluation,
                &format!(
                    "scalar 3L cross-section {} {label} rich inspect parity: LTD local-4D vs CFF local-3D",
                    case.graph
                ),
            );
        }
        assert_weight_benchmark(case, &label, &cff_3d_result);
    }

    clean_test(&cff_3d.cli_settings.state.folder);
    clean_test(&cff_4d.cli_settings.state.folder);
    if let Some(ltd_4d) = ltd_4d {
        clean_test(&ltd_4d.cli_settings.state.folder);
    }
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

fn run_scalar_3l_quartic_energy_numerator_case(case: Scalar3LGraphCase, edge: usize) -> Result<()> {
    let numerator_case = case.with_quartic_edge_numerator(edge);
    run_scalar_3l_cross_section_numerator_only_case(
        &format!("scalar_3l_quartic_energy_q{edge}"),
        numerator_case,
    )
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

    scalar_3l_graph_test!(scalar_3l_cross_section_gl02_inspects_match, "GL02");

    mod quadratic_energy_numerators {
        use super::*;

        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl00_quadratic_energy_inspects_match,
            q7_only: Scalar3LGraphCase::default("GL00").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl16_quadratic_energy_inspects_match,
            q1_only: Scalar3LGraphCase::default("GL16")
        );
    }

    mod quartic_energy_numerators {
        use super::*;

        quartic_energy_graph_test!(
            scalar_3l_cross_section_gl02_quartic_energy_inspects_match,
            "GL02"
        );
        quartic_energy_graph_test!(
            scalar_3l_cross_section_gl16_quartic_energy_inspects_match,
            "GL16"
        );
        quartic_energy_graph_test!(
            scalar_3l_cross_section_gl24_quartic_energy_inspects_match,
            Scalar3LGraphCase::default("GL24").without_threshold_subtraction()
        );
    }
}

mod slow {
    use super::*;

    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl00_inspects_match,
        Scalar3LGraphCase::default("GL00").without_uv_subtraction()
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl01_inspects_match,
        Scalar3LGraphCase::default("GL01")
            .without_uv_subtraction()
            .without_threshold_subtraction()
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl03_inspects_match,
        Scalar3LGraphCase::default("GL03")
            .without_uv_subtraction()
            .without_threshold_subtraction()
            .with_numerator(NumeratorChoice::Dot {
                left_edge: 7,
                right_edge: 8,
                dummy_index: 1,
            })
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl04_inspects_match,
        Scalar3LGraphCase::default("GL04")
            .without_uv_subtraction()
            .without_threshold_subtraction()
    );
    scalar_3l_graph_test!(scalar_3l_cross_section_gl05_inspects_match, "GL05");
    scalar_3l_graph_test!(scalar_3l_cross_section_gl06_inspects_match, "GL06");
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl07_inspects_match,
        Scalar3LGraphCase::default("GL07").without_uv_subtraction()
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl08_inspects_match,
        Scalar3LGraphCase::default("GL08").without_uv_subtraction()
    );
    scalar_3l_graph_test!(scalar_3l_cross_section_gl09_inspects_match, "GL09");
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl10_inspects_match,
        Scalar3LGraphCase::default("GL10").without_uv_subtraction()
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl11_inspects_match,
        Scalar3LGraphCase::default("GL11").without_uv_subtraction()
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl12_inspects_match,
        Scalar3LGraphCase::default("GL12").without_uv_subtraction()
    );
    scalar_3l_graph_test!(scalar_3l_cross_section_gl13_inspects_match, "GL13");
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl14_inspects_match,
        Scalar3LGraphCase::default("GL14").without_uv_subtraction()
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl15_inspects_match,
        Scalar3LGraphCase::default("GL15").without_uv_subtraction()
    );
    scalar_3l_graph_test!(scalar_3l_cross_section_gl16_inspects_match, "GL16");
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl17_inspects_match,
        Scalar3LGraphCase::default("GL17").without_uv_subtraction()
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl18_inspects_match,
        Scalar3LGraphCase::default("GL18").without_uv_subtraction()
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl19_inspects_match,
        Scalar3LGraphCase::default("GL19").without_uv_subtraction()
    );
    scalar_3l_graph_test!(scalar_3l_cross_section_gl20_inspects_match, "GL20");
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl21_inspects_match,
        Scalar3LGraphCase::default("GL21").without_uv_subtraction()
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
        Scalar3LGraphCase::default("GL23")
            .without_uv_subtraction()
            .with_numerator(NumeratorChoice::Dot {
                left_edge: 1,
                right_edge: 2,
                dummy_index: 1,
            })
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl24_inspects_match,
        Scalar3LGraphCase::default("GL24")
            .without_threshold_subtraction()
            .with_numerator(NumeratorChoice::Dot {
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
        Scalar3LGraphCase::default("GL26").without_uv_subtraction()
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
    scalar_3l_graph_test!(scalar_3l_cross_section_gl29_inspects_match, "GL29");
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl30_inspects_match,
        Scalar3LGraphCase::default("GL30")
            .without_uv_subtraction()
            .without_threshold_subtraction()
            .with_numerator(NumeratorChoice::Dot {
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
        Scalar3LGraphCase::default("GL33").without_uv_subtraction()
    );
    scalar_3l_graph_test!(scalar_3l_cross_section_gl34_inspects_match, "GL34");
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl35_inspects_match,
        Scalar3LGraphCase::default("GL35").without_threshold_subtraction()
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl36_inspects_match,
        Scalar3LGraphCase::default("GL36").without_uv_subtraction()
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl37_inspects_match,
        Scalar3LGraphCase::default("GL37").without_uv_subtraction()
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl38_inspects_match,
        Scalar3LGraphCase::default("GL38").without_uv_subtraction()
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
        Scalar3LGraphCase::default("GL44").without_uv_subtraction()
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl45_inspects_match,
        Scalar3LGraphCase::default("GL45").without_uv_subtraction()
    );
    scalar_3l_graph_test!(
        scalar_3l_cross_section_gl46_inspects_match,
        Scalar3LGraphCase::default("GL46").without_uv_subtraction()
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
            q1_only: Scalar3LGraphCase::default("GL00").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl01_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL01")
                .without_uv_subtraction()
                .without_threshold_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl02_quadratic_energy_inspects_match,
            q1_only: Scalar3LGraphCase::default("GL02")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl03_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL03")
                .without_uv_subtraction()
                .without_threshold_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl04_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL04").without_uv_subtraction()
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
            Scalar3LGraphCase::default("GL07").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl08_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL08").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl09_quadratic_energy_inspects_match,
            "GL09"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl10_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL10").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl11_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL11").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl12_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL12").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl13_quadratic_energy_inspects_match,
            "GL13"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl14_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL14").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl15_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL15").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl16_quadratic_energy_inspects_match,
            q7_only: Scalar3LGraphCase::default("GL16")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl17_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL17").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl18_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL18").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl19_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL19").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl20_quadratic_energy_inspects_match,
            "GL20"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl21_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL21").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl22_quadratic_energy_inspects_match,
            "GL22"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl23_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL23").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl24_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL24").without_threshold_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl25_quadratic_energy_inspects_match,
            "GL25"
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl26_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL26").without_uv_subtraction()
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
                .without_uv_subtraction()
                .without_threshold_subtraction()
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
            Scalar3LGraphCase::default("GL33").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl34_quadratic_energy_inspects_match,
            q1_only: Scalar3LGraphCase::default("GL34")
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl35_quadratic_energy_inspects_match,
            q1_only: Scalar3LGraphCase::default("GL35").without_threshold_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl36_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL36").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl37_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL37").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl38_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL38").without_uv_subtraction()
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
            Scalar3LGraphCase::default("GL44").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl45_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL45").without_uv_subtraction()
        );
        quadratic_energy_graph_test!(
            scalar_3l_cross_section_gl46_quadratic_energy_inspects_match,
            Scalar3LGraphCase::default("GL46").without_uv_subtraction()
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

        quartic_energy_graph_test!(
            scalar_3l_cross_section_gl00_quartic_energy_inspects_match,
            Scalar3LGraphCase::default("GL00").without_uv_subtraction()
        );
        quartic_energy_graph_test!(
            scalar_3l_cross_section_gl06_quartic_energy_inspects_match,
            "GL06"
        );
        quartic_energy_graph_test!(
            scalar_3l_cross_section_gl10_quartic_energy_inspects_match,
            Scalar3LGraphCase::default("GL10").without_uv_subtraction()
        );
        quartic_energy_graph_test!(
            scalar_3l_cross_section_gl20_quartic_energy_inspects_match,
            "GL20"
        );
        quartic_energy_graph_test!(
            scalar_3l_cross_section_gl30_quartic_energy_inspects_match,
            Scalar3LGraphCase::default("GL30")
                .without_uv_subtraction()
                .without_threshold_subtraction()
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
