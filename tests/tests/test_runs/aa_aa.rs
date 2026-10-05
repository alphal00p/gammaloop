use super::utils::*;
use super::*;
use gammaloop_api::commands::duplicate::DuplicateIntegrand;
use gammalooprs::{
    momentum::{Dep, ExternalMomenta, Helicity, SignOrZero},
    settings::runtime::kinematic::Externals,
    utils::F,
};
use std::collections::BTreeMap;

const AA_AA_PROCESS: &str = "aa_aa_all_helicities";
const AA_AA_ASSEMBLY_PROCESS: &str = "aa_aa_all_helicities_assembly";

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
    gammaloop_integration_tests::workspace_root()
        .join("tests/resources")
        .join("benchmarks")
        .join(name)
}

// The f rows use independent MSbar Feynman-parameter box integrals at mu_r=20,
// including -GC_2^4*S'(4)/(32*pi^2), where S(D)=D(D+2)<N_rank4>.
// Their errors cover contour/order convergence and rounding. Graph 2 is
// spacelike and has exactly zero absorptive part. All targets store M=i*A.
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

fn add_graph_id_observable(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
) -> Result<()> {
    cli.run_command(&format!(
        "set process -p {process} string '\n[quantities.graph_id]\ntype = \"graph_id\"\n\n[observables.graph_id_hist_real]\nquantity = \"graph_id\"\nkind = \"discrete\"\nphase = \"real\"\ndomain = {{ type = \"explicit_range\", min = 0, max = 2 }}\n\n[observables.graph_id_hist_imag]\nquantity = \"graph_id\"\nkind = \"discrete\"\nphase = \"imag\"\ndomain = {{ type = \"explicit_range\", min = 0, max = 2 }}\n'"
    ))
}

mod important {
    use super::*;

    #[test]
    #[serial]
    fn aa_aa_local_inspect_backend_consistency() -> Result<()> {
        let point = vec![0.11, -0.07, 0.19];
        let cases = [
            ("a", "1L_ppmm", 0),
            ("a", "1L_mpmm", 1),
            ("a", "1L_mmmm", 2),
            ("a", "1L_pmmp", 0),
            ("a", "1L_pmpm", 1),
            ("b", "1L_ppmm", 2),
            ("b", "1L_mpmm", 0),
            ("b", "1L_mmmm", 1),
            ("b", "1L_pmmp", 2),
            ("b", "1L_pmpm", 0),
        ];

        let mut symjit_cli = setup_aa_aa_cli("aa_aa_local_inspect_backend_consistency_symjit")?;
        generate_aa_aa_helicity_family(
            &mut symjit_cli,
            AA_AA_PROCESS,
            "symjit",
            &AA_AA_HELICITIES_ALL,
        )?;
        let mut assembly_cli = setup_aa_aa_cli("aa_aa_local_inspect_backend_consistency_assembly")?;
        generate_aa_aa_helicity_family(
            &mut assembly_cli,
            AA_AA_ASSEMBLY_PROCESS,
            "assembly",
            &AA_AA_HELICITIES_ALL,
        )?;

        for (integrand, _) in AA_AA_HELICITIES_ALL {
            set_single_precision_level(&mut symjit_cli, AA_AA_PROCESS, integrand, "Double")?;
            set_single_precision_level(
                &mut assembly_cli,
                AA_AA_ASSEMBLY_PROCESS,
                integrand,
                "Double",
            )?;
        }

        for (kinematics, integrand, graph_id) in cases {
            let context =
                format!("kinematics={kinematics}, integrand={integrand}, graph_id={graph_id}");

            set_aa_aa_kinematics(&mut symjit_cli, AA_AA_PROCESS, kinematics)?;
            let (_, symjit_result) = Inspect {
                process: Some(ProcessRef::Unqualified(AA_AA_PROCESS.to_string())),
                integrand_name: Some(integrand.to_string()),
                point: point.clone(),
                momentum_space: true,
                graph_id: Some(graph_id),
                ..Default::default()
            }
            .run(&mut symjit_cli.state)?;
            set_aa_aa_kinematics(&mut assembly_cli, AA_AA_ASSEMBLY_PROCESS, kinematics)?;
            let (_, assembly_result) = Inspect {
                process: Some(ProcessRef::Unqualified(AA_AA_ASSEMBLY_PROCESS.to_string())),
                integrand_name: Some(integrand.to_string()),
                point: point.clone(),
                momentum_space: true,
                graph_id: Some(graph_id),
                ..Default::default()
            }
            .run(&mut assembly_cli.state)?;
            assert_complex_approx_eq(
                assembly_result,
                symjit_result,
                format!("{context}: assembly vs symjit inspect"),
            );
        }

        clean_test(&assembly_cli.cli_settings.state.folder);
        clean_test(&symjit_cli.cli_settings.state.folder);
        Ok(())
    }
}

#[test]
#[serial]
fn aa_aa_integrated_graph_histogram_bins() -> Result<()> {
    let targets = load_integrated_targets()?;
    let mut cli = setup_aa_aa_cli("aa_aa_integrated_graph_histogram_bins_generated")?;
    // The historical saved state was generated with SM-default parameters;
    // its targets use aEWM1=132.507, rather than the inspection card's 137.036.
    cli.run_command("set model aEWM1=132.507")?;
    {
        use symbolica::atom::{Atom, AtomCore, Symbol};

        let ee = cli.state.model.get_parameter("ee");
        let alpha = cli.state.model.get_parameter("aEW");
        let inverse_alpha = cli.state.model.get_parameter("aEWM1");
        let ee_symbol: Atom = ee.name.into();
        let alpha_symbol: Atom = alpha.name.into();
        let inverse_alpha_symbol: Atom = inverse_alpha.name.into();
        // Four photon-top vertices give exactly GC_2^4 =
        // 256*pi^2/(81*aEWM1^2), with no residual coupling phase.
        let coupling = cli
            .state
            .model
            .get_coupling("GC_2")
            .expression
            .replace(ee_symbol.to_pattern())
            .with(ee.expression.as_ref().unwrap().to_pattern())
            .replace(alpha_symbol.to_pattern())
            .with(alpha.expression.as_ref().unwrap().to_pattern());
        assert_eq!(
            coupling.pow(4) * inverse_alpha_symbol.pow(2),
            Atom::num((256, 81)) * Atom::var(Symbol::PI).pow(2),
        );
    }
    generate_aa_aa_helicity_family(
        &mut cli,
        AA_AA_PROCESS,
        "symjit",
        &AA_AA_HELICITIES_INTEGRATED,
    )?;
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
            // GammaLoop returns A=-i*M; the historical targets store M.
            // Compare i*A, rotating the component errors with their means.
            assert_histogram_estimate_compatible(
                -actual_im_avg,
                actual_im_err,
                target.re_avg,
                target.re_err,
                &format!(
                    "integrated real bin kinematics={kinematics}, integrand={integrand}, graph_id={graph_id}"
                ),
            );
            assert_histogram_estimate_compatible(
                actual_re_avg,
                actual_re_err,
                target.im_avg,
                target.im_err,
                &format!(
                    "integrated imag bin kinematics={kinematics}, integrand={integrand}, graph_id={graph_id}"
                ),
            );
        }
    }

    clean_test(&workspace_root);
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
