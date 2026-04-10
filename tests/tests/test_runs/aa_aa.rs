use super::utils::*;
use super::*;
use gammaloop_api::CLISettings;
use gammaloop_api::commands::duplicate::DuplicateIntegrand;
use gammaloop_api::commands::evaluate_samples::{
    EvaluateSamples, EvaluateSamplesPrecise, evaluate_sample_precise,
};
use gammaloop_api::session::CliSessionState;
use gammaloop_api::state::{RunHistory, State};
use gammalooprs::{
    initialisation::initialise,
    integrands::evaluation::PreciseEvaluationResultOutput,
    momentum::{Dep, ExternalMomenta, Helicity, SignOrZero},
    settings::{RuntimeSettings, runtime::kinematic::Externals},
    utils::{ArbPrec, F, FloatLike, PrecisionUpgradable, f128},
};
use ndarray::Array2;
use std::collections::BTreeMap;
use symbolica::domains::float::Float as SymbolicaFloat;
use symbolica::domains::float::Real;

const AA_AA_PROCESS: &str = "aa_aa_all_helicities";
const AA_AA_ASSEMBLY_PROCESS: &str = "aa_aa_all_helicities_assembly";
const AA_AA_PRIMARY_INTEGRAND: &str = "1L_ppmm";
const AA_AA_GRAPH_COUNT: usize = 3;
const HISTOGRAM_TARGET_N_SIGMA: f64 = 5.0;

const AA_AA_HELICITIES_ALL: [(&str, &str); 5] = [
    ("1L_ppmm", "[+1,+1,-1,-1]"),
    ("1L_mpmm", "[-1,+1,-1,-1]"),
    ("1L_mmmm", "[-1,-1,-1,-1]"),
    ("1L_pmmp", "[+1,-1,-1,+1]"),
    ("1L_pmpm", "[+1,-1,+1,-1]"),
];

const AA_AA_HELICITIES_INTEGRATED: [(&str, &str); 2] =
    [("1L_ppmm", "[+1,+1,-1,-1]"), ("1L_pmpm", "[+1,-1,+1,-1]")];

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum InspectTargetKind {
    F64Symjit,
    F64Assembly,
    Quad,
    Arb,
}

impl InspectTargetKind {
    fn from_target_tag(tag: &str) -> Result<Self> {
        match tag {
            "F64_SYMJIT_TARGET" => Ok(Self::F64Symjit),
            "F64_ASSEMBLY_TARGET" => Ok(Self::F64Assembly),
            "QUAD_TARGET" => Ok(Self::Quad),
            "ARB_TARGET" => Ok(Self::Arb),
            _ => Err(eyre::eyre!("Unsupported inspect target tag '{tag}'")),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct InspectTargetKey {
    kind: InspectTargetKind,
    kinematics: String,
    integrand: String,
    graph_id: usize,
}

impl InspectTargetKey {
    fn new(
        kind: InspectTargetKind,
        kinematics: impl Into<String>,
        integrand: impl Into<String>,
        graph_id: usize,
    ) -> Self {
        Self {
            kind,
            kinematics: kinematics.into(),
            integrand: integrand.into(),
            graph_id,
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct IntegratedHistogramTarget {
    re_avg: f64,
    re_err: f64,
    im_avg: f64,
    im_err: f64,
}

fn parse_aa_aa_helicities(helicities: &str) -> Vec<Helicity> {
    match helicities {
        "[+1,+1,-1,-1]" => vec![
            Helicity::Signed(SignOrZero::Plus),
            Helicity::Signed(SignOrZero::Plus),
            Helicity::Signed(SignOrZero::Minus),
            Helicity::Signed(SignOrZero::Minus),
        ],
        "[-1,+1,-1,-1]" => vec![
            Helicity::Signed(SignOrZero::Minus),
            Helicity::Signed(SignOrZero::Plus),
            Helicity::Signed(SignOrZero::Minus),
            Helicity::Signed(SignOrZero::Minus),
        ],
        "[-1,-1,-1,-1]" => vec![
            Helicity::Signed(SignOrZero::Minus),
            Helicity::Signed(SignOrZero::Minus),
            Helicity::Signed(SignOrZero::Minus),
            Helicity::Signed(SignOrZero::Minus),
        ],
        "[+1,-1,-1,+1]" => vec![
            Helicity::Signed(SignOrZero::Plus),
            Helicity::Signed(SignOrZero::Minus),
            Helicity::Signed(SignOrZero::Minus),
            Helicity::Signed(SignOrZero::Plus),
        ],
        "[+1,-1,+1,-1]" => vec![
            Helicity::Signed(SignOrZero::Plus),
            Helicity::Signed(SignOrZero::Minus),
            Helicity::Signed(SignOrZero::Plus),
            Helicity::Signed(SignOrZero::Minus),
        ],
        other => panic!("unsupported aa_aa helicity pattern {other}"),
    }
}

fn benchmark_resource_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("resources")
        .join("benchmarks")
        .join(name)
}

fn example_aa_aa_state_folder() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("integration-tests crate must live under the workspace root")
        .join("examples/cli/aa_aa/1L/state")
}

fn load_inspect_targets() -> Result<BTreeMap<InspectTargetKey, (String, String)>> {
    let mut targets = BTreeMap::new();
    for line in
        std::fs::read_to_string(benchmark_resource_path("aa_aa_local_inspect_targets.txt"))?.lines()
    {
        if line.is_empty() || line.ends_with("_START") || line.ends_with("_END") {
            continue;
        }
        let parts = line.split('|').collect_vec();
        if parts.len() != 6 {
            return Err(eyre::eyre!("Malformed inspect target line '{line}'"));
        }
        let graph_id = parts[3]
            .parse()
            .map_err(|err| eyre::eyre!("Invalid graph id in inspect target '{line}': {err}"))?;
        targets.insert(
            InspectTargetKey::new(
                InspectTargetKind::from_target_tag(parts[0])?,
                parts[1],
                parts[2],
                graph_id,
            ),
            (parts[4].to_string(), parts[5].to_string()),
        );
    }
    Ok(targets)
}

fn load_integrated_targets() -> Result<BTreeMap<(String, String, usize), IntegratedHistogramTarget>>
{
    let mut targets = BTreeMap::new();
    for line in
        std::fs::read_to_string(benchmark_resource_path("aa_aa_integrated_targets.txt"))?.lines()
    {
        if line.is_empty() || line.ends_with("_START") || line.ends_with("_END") {
            continue;
        }
        let parts = line.split('|').collect_vec();
        if parts.len() != 8 || parts[0] != "INTEGRATED_TARGET" {
            return Err(eyre::eyre!("Malformed integrated target line '{line}'"));
        }
        let graph_id = parts[3]
            .parse()
            .map_err(|err| eyre::eyre!("Invalid graph id in integrated target '{line}': {err}"))?;
        targets.insert(
            (parts[1].to_string(), parts[2].to_string(), graph_id),
            IntegratedHistogramTarget {
                re_avg: parts[4].parse().map_err(|err| {
                    eyre::eyre!("Invalid real average in integrated target '{line}': {err}")
                })?,
                re_err: parts[5].parse().map_err(|err| {
                    eyre::eyre!("Invalid real error in integrated target '{line}': {err}")
                })?,
                im_avg: parts[6].parse().map_err(|err| {
                    eyre::eyre!("Invalid imaginary average in integrated target '{line}': {err}")
                })?,
                im_err: parts[7].parse().map_err(|err| {
                    eyre::eyre!("Invalid imaginary error in integrated target '{line}': {err}")
                })?,
            },
        );
    }
    Ok(targets)
}

fn inspect_target_entry<'a>(
    targets: &'a BTreeMap<InspectTargetKey, (String, String)>,
    kind: InspectTargetKind,
    kinematics: &str,
    integrand: &str,
    graph_id: usize,
) -> &'a (String, String) {
    targets
        .get(&InspectTargetKey::new(kind, kinematics, integrand, graph_id))
        .unwrap_or_else(|| {
            panic!(
                "missing inspect target for kind={kind:?}, kinematics={kinematics}, integrand={integrand}, graph_id={graph_id}"
            )
        })
}

fn f64_inspect_target(
    targets: &BTreeMap<InspectTargetKey, (String, String)>,
    kind: InspectTargetKind,
    kinematics: &str,
    integrand: &str,
    graph_id: usize,
) -> Complex<f64> {
    let (re, im) = inspect_target_entry(targets, kind, kinematics, integrand, graph_id);
    Complex::new(
        re.parse().expect("f64 real inspect target must parse"),
        im.parse().expect("f64 imaginary inspect target must parse"),
    )
}

fn precise_inspect_target<T>(
    targets: &BTreeMap<InspectTargetKey, (String, String)>,
    kind: InspectTargetKind,
    kinematics: &str,
    integrand: &str,
    graph_id: usize,
) -> Result<Complex<F<T>>>
where
    T: gammalooprs::utils::FloatLike + From<SymbolicaFloat>,
{
    let (re, im) = inspect_target_entry(targets, kind, kinematics, integrand, graph_id);
    decimal_complex(re, im)
}

fn exact_f64_as_quad(value: Complex<f64>) -> Complex<F<f128>> {
    Complex::new(
        F(f128::from(SymbolicaFloat::from(value.re))),
        F(f128::from(SymbolicaFloat::from(value.im))),
    )
}

fn exact_quad_as_arb(value: &Complex<F<f128>>) -> Complex<F<ArbPrec>> {
    Complex::new(F(value.re.0.higher()), F(value.im.0.higher()))
}

fn quad_vs_f64_tolerance() -> F<f128> {
    F(f128::from(SymbolicaFloat::from(f64::EPSILON.sqrt())))
}

fn arb_vs_quad_tolerance() -> F<ArbPrec> {
    F(f128::default().epsilon().sqrt().higher())
}

fn quad_integrand_result(result: PreciseEvaluationResultOutput) -> Result<Complex<F<f128>>> {
    match result {
        PreciseEvaluationResultOutput::Quad(result) => Ok(result.integrand_result),
        _ => Err(eyre::eyre!("Expected a Quad precise evaluation result")),
    }
}

fn arb_integrand_result(result: PreciseEvaluationResultOutput) -> Result<Complex<F<ArbPrec>>> {
    match result {
        PreciseEvaluationResultOutput::Arb(result) => Ok(result.integrand_result),
        _ => Err(eyre::eyre!("Expected an ArbPrec precise evaluation result")),
    }
}

fn assert_histogram_estimate_compatible(
    actual_value: f64,
    actual_error: f64,
    target_value: f64,
    target_error: f64,
    context: &str,
) {
    let sigma = (actual_error * actual_error + target_error * target_error)
        .sqrt()
        .max(1.0e-30);
    let delta = (actual_value - target_value).abs();
    assert!(
        delta <= HISTOGRAM_TARGET_N_SIGMA * sigma,
        "{context}: actual={actual_value:e} ± {actual_error:e}, target={target_value:e} ± {target_error:e}, |Δ|={delta:e}, allowed={}σ={:e}",
        HISTOGRAM_TARGET_N_SIGMA,
        HISTOGRAM_TARGET_N_SIGMA * sigma,
    );
}

fn setup_aa_aa_cli(test_name: &str) -> Result<gammaloop_integration_tests::CLIState> {
    let mut cli = get_test_cli(
        Some("aa_aa_1L.toml".into()),
        get_tests_workspace_path().join(test_name),
        Some(test_name.to_string()),
        true,
    )?;
    cli.run_command("import model sm-default.json")?;
    cli.run_command("run set_model_parameters")?;
    Ok(cli)
}

fn ensure_aa_aa_helicity_family(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    helicities: &[(&str, &str)],
) -> Result<()> {
    let (source_integrand, source_helicities) = helicities
        .first()
        .copied()
        .expect("helicity family must not be empty");
    let process_id = cli
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified(process.to_string())))?;
    let source_settings = cli
        .state
        .process_list
        .get_integrand_mut(process_id, source_integrand)?
        .get_mut_settings();
    let Externals::Constant {
        helicities: source_existing,
        ..
    } = &mut source_settings.kinematics.externals;
    *source_existing = parse_aa_aa_helicities(source_helicities);
    for (integrand, helicities) in helicities.iter().skip(1) {
        let integrand_exists = cli
            .state
            .find_integrand_ref(
                Some(&ProcessRef::Unqualified(process.to_string())),
                Some(&integrand.to_string()),
            )
            .is_ok();
        if !integrand_exists {
            DuplicateIntegrand {
                process: Some(ProcessRef::Unqualified(process.to_string())),
                integrand_name: Some(source_integrand.to_string()),
                output_process_name: process.to_string(),
                output_integrand_name: integrand.to_string(),
            }
            .run(&mut cli.state)?;
        }
        let settings = cli
            .state
            .process_list
            .get_integrand_mut(process_id, integrand)?
            .get_mut_settings();
        let Externals::Constant {
            helicities: existing,
            ..
        } = &mut settings.kinematics.externals;
        *existing = parse_aa_aa_helicities(helicities);
    }
    for (integrand, _) in helicities {
        let settings = cli
            .state
            .process_list
            .get_integrand_mut(process_id, integrand)?
            .get_mut_settings();
        settings.general.m_uv = 20.0;
        settings.general.mu_r = 20.0;
    }
    Ok(())
}

fn load_example_aa_aa_cli() -> Result<gammaloop_integration_tests::CLIState> {
    initialise()?;
    let state_folder = example_aa_aa_state_folder();

    let mut state = State::load(state_folder.clone(), None, None)?;
    state.activate_loaded_integrand_backends(true)?;

    let mut cli_settings: CLISettings = toml::from_str(&std::fs::read_to_string(
        state_folder.join("global_settings.toml"),
    )?)?;
    cli_settings.state.folder = state_folder.clone();
    cli_settings.session.read_only_state = true;
    cli_settings.global.n_cores.integrate = 1;
    cli_settings.sync_settings()?;

    let default_runtime_settings: RuntimeSettings = toml::from_str(&std::fs::read_to_string(
        state_folder.join("default_runtime_settings.toml"),
    )?)?;
    let run_history = RunHistory::load(state_folder.join("run.toml"))?;
    let mut cli = gammaloop_integration_tests::CLIState {
        state,
        cli_settings,
        default_runtime_settings,
        run_history,
        session_state: CliSessionState::default(),
    };
    ensure_aa_aa_helicity_family(&mut cli, AA_AA_PROCESS, &AA_AA_HELICITIES_ALL)?;
    Ok(cli)
}

fn set_aa_aa_kinematics(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    label: &str,
) -> Result<()> {
    let momenta = match label {
        "a" => vec![
            ExternalMomenta::Independent([F(86.5), F(0.0), F(0.0), F(86.5)]),
            ExternalMomenta::Independent([F(86.5), F(0.0), F(0.0), F(-86.5)]),
            ExternalMomenta::Independent([F(86.5), F(0.0), F(-79.27855952273603), F(-34.6)]),
            ExternalMomenta::Dependent(Dep::Dep),
        ],
        "b" => vec![
            ExternalMomenta::Independent([F(103.8), F(0.0), F(0.0), F(103.8)]),
            ExternalMomenta::Independent([F(103.8), F(0.0), F(0.0), F(-103.8)]),
            ExternalMomenta::Independent([F(103.8), F(0.0), F(-95.13427142728324), F(-41.52)]),
            ExternalMomenta::Dependent(Dep::Dep),
        ],
        "f" => vec![
            ExternalMomenta::Independent([
                F(1751.8893288527534),
                F(0.0),
                F(0.0),
                F(1751.8893288527534),
            ]),
            ExternalMomenta::Independent([
                F(1751.8893288527534),
                F(0.0),
                F(0.0),
                F(-1751.8893288527534),
            ]),
            ExternalMomenta::Independent([
                F(1751.8893288527534),
                F(0.0),
                F(-1605.6330917306252),
                F(-700.755731541101),
            ]),
            ExternalMomenta::Dependent(Dep::Dep),
        ],
        other => {
            return Err(eyre::eyre!("Unsupported aa_aa kinematics label '{other}'"));
        }
    };
    let process_id = cli
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified(process.to_string())))?;
    for (integrand, _) in AA_AA_HELICITIES_ALL {
        if cli
            .state
            .find_integrand_ref(
                Some(&ProcessRef::Unqualified(process.to_string())),
                Some(&integrand.to_string()),
            )
            .is_ok()
        {
            let settings = cli
                .state
                .process_list
                .get_integrand_mut(process_id, integrand)?
                .get_mut_settings();
            let Externals::Constant {
                momenta: existing, ..
            } = &mut settings.kinematics.externals;
            *existing = momenta.clone();
        }
    }
    Ok(())
}

fn set_compilation_mode(
    cli: &mut gammaloop_integration_tests::CLIState,
    compilation_mode: &str,
) -> Result<()> {
    cli.run_command(&format!(
        "set global string '\n[global.generation.compile]\ncompilation_mode = \"{compilation_mode}\"\noptimization_level = \"O0\"\nfast_math = false\nunsafe_math = false\ncustom = []\n'"
    ))
}

fn generate_aa_aa_helicity_family(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    compilation_mode: &str,
    helicities: &[(&str, &str)],
) -> Result<()> {
    let first_integrand = helicities
        .first()
        .map(|(integrand, _)| *integrand)
        .expect("helicity family must not be empty");
    set_compilation_mode(cli, compilation_mode)?;
    cli.run_command(&format!(
        "import graphs ./examples/cli/aa_aa/1L/aa_aa_1L.dot -p {process} -i {first_integrand} -o"
    ))?;
    cli.run_command(&format!("generate existing -p {process}"))?;
    ensure_aa_aa_helicity_family(cli, process, helicities)?;
    Ok(())
}

fn set_single_precision_level(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    integrand: &str,
    precision: &str,
) -> Result<()> {
    let thresholds = if precision == "Double" {
        "required_precision_for_re = 1.0e-10\nrequired_precision_for_im = 1.0e-10\nescalate_for_large_weight_threshold = 0.9"
    } else {
        "required_precision_for_re = 1.0e-30\nrequired_precision_for_im = 1.0e-30\nescalate_for_large_weight_threshold = -1.0"
    };
    cli.run_command(&format!(
        "set process -p {process} -i {integrand} string '\n[stability]\ncheck_on_norm = true\nescalate_if_exact_zero = false\nloop_momenta_norm_escalation_factor = -1.0\n\n[[stability.rotation_axis]]\ntype = \"z\"\n\n[[stability.levels]]\nprecision = \"{precision}\"\n{thresholds}\n'"
    ))
}

fn evaluate_graph_momentum_sample_f64(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    integrand: &str,
    graph_id: usize,
    point: &[f64],
) -> Result<Complex<f64>> {
    let process_id = cli
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified(process.to_string())))?;
    let graph_name = cli
        .state
        .process_list
        .get_integrand(process_id, integrand)?
        .require_generated()?
        .graph_name_by_id(graph_id)
        .expect("graph id must exist")
        .to_string();
    let points = Array2::from_shape_vec((1, point.len()), point.to_vec())?;
    let result = EvaluateSamples {
        process_id: Some(process_id),
        integrand_name: Some(integrand.to_string()),
        use_arb_prec: false,
        minimal_output: true,
        return_generated_events: None,
        momentum_space: true,
        points: points.view(),
        integrator_weights: None,
        discrete_dims: None,
        graph_names: Some(vec![Some(graph_name)]),
        orientations: None,
    }
    .run(&mut cli.state)?;
    Ok(result.samples[0]
        .evaluation
        .integrand_result
        .map(|entry| entry.0))
}

fn default_momentum_space_point_for(
    cli: &gammaloop_integration_tests::CLIState,
    process: &str,
    integrand: &str,
) -> Result<Vec<f64>> {
    let process_id = cli
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified(process.to_string())))?;
    let integrand = cli
        .state
        .process_list
        .get_integrand(process_id, integrand)?
        .require_generated()?;
    let n_dim = integrand.get_n_dim();
    let seed = [0.11, -0.07, 0.19, -0.13, 0.05, 0.29];
    Ok((0..n_dim).map(|index| seed[index % seed.len()]).collect())
}

fn evaluate_graph_momentum_sample_precise(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    integrand: &str,
    graph_id: usize,
    point: &[f64],
    use_arb_prec: bool,
) -> Result<PreciseEvaluationResultOutput> {
    let process_id = cli
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified(process.to_string())))?;
    let graph_name = cli
        .state
        .process_list
        .get_integrand(process_id, integrand)?
        .require_generated()?
        .graph_name_by_id(graph_id)
        .expect("graph id must exist")
        .to_string();
    let points = Array2::from_shape_vec((1, point.len()), point.to_vec())?;
    Ok(evaluate_sample_precise(
        &mut cli.state,
        &EvaluateSamplesPrecise {
            process_id: Some(process_id),
            integrand_name: Some(integrand.to_string()),
            use_arb_prec,
            minimal_output: true,
            return_generated_events: None,
            momentum_space: true,
            points: points.view(),
            integrator_weights: None,
            discrete_dims: None,
            graph_names: Some(vec![Some(graph_name)]),
            orientations: None,
        },
    )?
    .sample
    .evaluation)
}

fn add_graph_id_observable(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
) -> Result<()> {
    cli.run_command(&format!(
        "set process -p {process} string '\n[quantities.graph_id]\ntype = \"graph_id\"\n\n[observables.graph_id_hist_real]\nquantity = \"graph_id\"\nkind = \"discrete\"\nphase = \"real\"\ndomain = {{ type = \"explicit_range\", min = 0, max = 2 }}\n\n[observables.graph_id_hist_imag]\nquantity = \"graph_id\"\nkind = \"discrete\"\nphase = \"imag\"\ndomain = {{ type = \"explicit_range\", min = 0, max = 2 }}\n'"
    ))
}

#[test]
#[serial]
#[ignore = "manual aa_aa inspect regression; Quad is currently known-buggy"]
fn aa_aa_local_inspect_precisions_and_backends() -> Result<()> {
    let targets = load_inspect_targets()?;
    let mut symjit_cli = setup_aa_aa_cli("aa_aa_local_inspect_precisions_and_backends_symjit")?;
    generate_aa_aa_helicity_family(
        &mut symjit_cli,
        AA_AA_PROCESS,
        "symjit",
        &AA_AA_HELICITIES_ALL,
    )?;
    let mut assembly_cli = setup_aa_aa_cli("aa_aa_local_inspect_precisions_and_backends")?;
    generate_aa_aa_helicity_family(
        &mut assembly_cli,
        AA_AA_ASSEMBLY_PROCESS,
        "assembly",
        &AA_AA_HELICITIES_ALL,
    )?;

    set_aa_aa_kinematics(&mut symjit_cli, AA_AA_PROCESS, "a")?;
    let point_a =
        default_momentum_space_point_for(&symjit_cli, AA_AA_PROCESS, AA_AA_PRIMARY_INTEGRAND)?;
    set_aa_aa_kinematics(&mut symjit_cli, AA_AA_PROCESS, "b")?;
    let point_b =
        default_momentum_space_point_for(&symjit_cli, AA_AA_PROCESS, AA_AA_PRIMARY_INTEGRAND)?;

    let quad_compat_tolerance = quad_vs_f64_tolerance();
    let arb_compat_tolerance = arb_vs_quad_tolerance();

    for (kinematics, point) in [("a", &point_a), ("b", &point_b)] {
        set_aa_aa_kinematics(&mut symjit_cli, AA_AA_PROCESS, kinematics)?;
        set_aa_aa_kinematics(&mut assembly_cli, AA_AA_ASSEMBLY_PROCESS, kinematics)?;
        for (integrand, _) in AA_AA_HELICITIES_ALL {
            for graph_id in 0..AA_AA_GRAPH_COUNT {
                let context =
                    format!("kinematics={kinematics}, integrand={integrand}, graph_id={graph_id}");

                set_single_precision_level(&mut symjit_cli, AA_AA_PROCESS, integrand, "Double")?;
                let symjit_result = evaluate_graph_momentum_sample_f64(
                    &mut symjit_cli,
                    AA_AA_PROCESS,
                    integrand,
                    graph_id,
                    point,
                )?;
                let symjit_target = f64_inspect_target(
                    &targets,
                    InspectTargetKind::F64Symjit,
                    kinematics,
                    integrand,
                    graph_id,
                );
                assert_complex_approx_eq(
                    symjit_result,
                    symjit_target,
                    &format!("{context}: symjit f64 benchmark"),
                );

                set_single_precision_level(
                    &mut assembly_cli,
                    AA_AA_ASSEMBLY_PROCESS,
                    integrand,
                    "Double",
                )?;
                let assembly_result = evaluate_graph_momentum_sample_f64(
                    &mut assembly_cli,
                    AA_AA_ASSEMBLY_PROCESS,
                    integrand,
                    graph_id,
                    point,
                )?;
                let assembly_target = f64_inspect_target(
                    &targets,
                    InspectTargetKind::F64Assembly,
                    kinematics,
                    integrand,
                    graph_id,
                );
                assert_complex_approx_eq(
                    assembly_result,
                    assembly_target,
                    &format!("{context}: assembly f64 benchmark"),
                );
                assert_complex_approx_eq(
                    assembly_result,
                    symjit_result,
                    &format!("{context}: assembly vs symjit f64"),
                );

                set_single_precision_level(&mut symjit_cli, AA_AA_PROCESS, integrand, "Quad")?;
                let quad_result = quad_integrand_result(evaluate_graph_momentum_sample_precise(
                    &mut symjit_cli,
                    AA_AA_PROCESS,
                    integrand,
                    graph_id,
                    point,
                    false,
                )?)?;
                let quad_target = precise_inspect_target::<f128>(
                    &targets,
                    InspectTargetKind::Quad,
                    kinematics,
                    integrand,
                    graph_id,
                )?;
                assert_complex_approx_eq_precise(
                    &quad_result,
                    &quad_target,
                    &quad_compat_tolerance,
                    &format!("{context}: quad benchmark"),
                );
                assert_complex_approx_eq_precise(
                    &quad_result,
                    &exact_f64_as_quad(symjit_result),
                    &quad_compat_tolerance,
                    &format!("{context}: quad vs f64 compatibility"),
                );

                set_single_precision_level(&mut symjit_cli, AA_AA_PROCESS, integrand, "Arb")?;
                let arb_result = arb_integrand_result(evaluate_graph_momentum_sample_precise(
                    &mut symjit_cli,
                    AA_AA_PROCESS,
                    integrand,
                    graph_id,
                    point,
                    true,
                )?)?;
                let arb_target = precise_inspect_target::<ArbPrec>(
                    &targets,
                    InspectTargetKind::Arb,
                    kinematics,
                    integrand,
                    graph_id,
                )?;
                assert_complex_approx_eq_precise(
                    &arb_result,
                    &arb_target,
                    &arb_compat_tolerance,
                    &format!("{context}: arb benchmark"),
                );
                assert_complex_approx_eq_precise(
                    &arb_result,
                    &exact_quad_as_arb(&quad_result),
                    &arb_compat_tolerance,
                    &format!("{context}: arb vs quad compatibility"),
                );
            }
        }
    }

    clean_test(&assembly_cli.cli_settings.state.folder);
    clean_test(&symjit_cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
#[ignore = "manual aa_aa integrated regression"]
fn aa_aa_integrated_graph_histogram_bins() -> Result<()> {
    let targets = load_integrated_targets()?;
    let mut cli = load_example_aa_aa_cli()?;
    let workspace_root = get_tests_workspace_path().join("aa_aa_integrated_graph_histogram_bins");

    clean_test(&workspace_root);
    add_graph_id_observable(&mut cli, AA_AA_PROCESS)?;
    cli.run_command(
        "set process -p aa_aa_all_helicities kv integrator.target_relative_accuracy=0.0 integrator.n_start=100000 integrator.n_increase=100000 integrator.n_max=1000000 integrator.seed=1337",
    )?;

    for (kinematics, integrand) in [
        ("a", AA_AA_HELICITIES_INTEGRATED[0].0),
        ("a", AA_AA_HELICITIES_INTEGRATED[1].0),
        ("f", AA_AA_HELICITIES_INTEGRATED[0].0),
        ("f", AA_AA_HELICITIES_INTEGRATED[1].0),
    ] {
        set_aa_aa_kinematics(&mut cli, AA_AA_PROCESS, kinematics)?;
        let result = Integrate {
            process: vec![ProcessRef::Unqualified(AA_AA_PROCESS.to_string())],
            integrand_name: vec![integrand.to_string()],
            workspace_path: Some(
                workspace_root.join(format!("{kinematics}_{integrand}_workspace")),
            ),
            n_cores: Some(1),
            restart: true,
            ..Default::default()
        }
        .run(&mut cli.state, &cli.cli_settings)?;

        let slot_key = format!("{AA_AA_PROCESS}@{integrand}");
        let bundle = result
            .slot_observables(&slot_key)
            .unwrap_or_else(|| panic!("missing observables bundle for slot '{slot_key}'"));
        let real_hist = bundle
            .histograms
            .get("graph_id_hist_real")
            .expect("missing graph_id_hist_real histogram");
        let imag_hist = bundle
            .histograms
            .get("graph_id_hist_imag")
            .expect("missing graph_id_hist_imag histogram");

        for graph_id in 0..AA_AA_GRAPH_COUNT {
            let target = targets
                .get(&(kinematics.to_string(), integrand.to_string(), graph_id))
                .unwrap_or_else(|| {
                    panic!(
                        "missing integrated target for kinematics={kinematics}, integrand={integrand}, graph_id={graph_id}"
                    )
                });
            let (actual_re_avg, actual_re_err) =
                discrete_histogram_bin_average_and_error(real_hist, graph_id as isize)?;
            let (actual_im_avg, actual_im_err) =
                discrete_histogram_bin_average_and_error(imag_hist, graph_id as isize)?;
            assert_histogram_estimate_compatible(
                actual_re_avg,
                actual_re_err,
                target.re_avg,
                target.re_err,
                &format!(
                    "integrated real bin kinematics={kinematics}, integrand={integrand}, graph_id={graph_id}"
                ),
            );
            assert_histogram_estimate_compatible(
                actual_im_avg,
                actual_im_err,
                target.im_avg,
                target.im_err,
                &format!(
                    "integrated imag bin kinematics={kinematics}, integrand={integrand}, graph_id={graph_id}"
                ),
            );
        }
    }

    clean_test(&workspace_root);
    Ok(())
}
