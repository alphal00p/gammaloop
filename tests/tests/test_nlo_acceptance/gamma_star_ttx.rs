use std::{
    collections::{BTreeMap, BTreeSet},
    f64::consts::PI,
};

use color_eyre::{Result, eyre::eyre};
use gammaloop_api::{
    commands::{
        Inspect, Integrate,
        evaluate_samples::{EvaluateSamplesPrecise, evaluate_samples_precise},
        integrate::RendererOption,
    },
    state::ProcessRef,
};
use gammaloop_integration_tests::{clean_test, get_example_cli, get_tests_workspace_path};
use gammalooprs::{
    graph::edge::PossibleParticle,
    integrands::{
        HasIntegrand, evaluation::PreciseEvaluationResultOutput, process::ProcessIntegrand,
    },
    settings::global::VectorPolarizationSumGauge,
    utils::{ArbPrec, F, FloatLike},
    uv::{
        CTIdentifier, CTRenormalizationRule, UltravioletGraph,
        settings::{ApproximationType, FinalIntegrandDimension},
    },
};
use linnet::half_edge::subgraph::SubSetOps;
use ndarray::Array2;
use spenso::algebra::{algebraic_traits::IsZero, complex::Complex};
use symbolica::domains::float::{Real, SingleFloat};
use which::which;

#[test]
fn gamma_star_ttx_msbar_nlo_matches_the_published_absolute_cross_sections() -> Result<()> {
    gamma_star_ttx_nlo_matches_reference(false, false)
}

#[test]
fn gamma_star_ttx_os_nlo_matches_the_analytic_inclusive_cross_section_and_auxiliary_scales()
-> Result<()> {
    gamma_star_ttx_nlo_matches_reference(true, false)
}

#[test]
fn gamma_star_ttx_os_sampling_probes_compare_route_maps() -> Result<()> {
    gamma_star_ttx_nlo_matches_reference(true, true)
}

// Pole-mass R-ratio, Chetyrkin, Kuehn and Steinhauser, hep-ph/9606230,
// equations (33)-(35): https://arxiv.org/pdf/hep-ph/9606230 .
// The current normalization follows from equation (7.1) of 2203.11038v1:
// https://arxiv.org/pdf/2203.11038v1 . The second result is the O(alpha_s)
// increment, including all real and virtual contributions, not LO + NLO.
fn pole_mass_vector_current_cross_section(e_cm: f64, mass: f64) -> (f64, f64) {
    assert!(mass > 0.0 && e_cm > 2.0 * mass);
    let dilog = |x: f64| {
        assert!((0.0..1.0).contains(&x));
        let mut sum = 0.0;
        let mut power = x;
        for n in 1..=4096 {
            sum += power / f64::from(n).powi(2);
            power *= x;
            // After n terms, sum_{k>n} x^k/k^2 is bounded above by this
            // geometric tail. This controls truncation, not f64 roundoff.
            let tail_bound = power / (f64::from(n + 1).powi(2) * (1.0 - x));
            if tail_bound <= f64::EPSILON * sum {
                return sum;
            }
        }
        panic!("dilogarithm series did not converge at x={x}");
    };
    let v = (1.0 - 4.0 * mass.powi(2) / e_cm.powi(2)).sqrt();
    let p = (1.0 - v) / (1.0 + v);
    let log_ratio = ((1.0 + v) / (1.0 - v)).ln();
    let r0 = 3.0 * v * (3.0 - v.powi(2)) / 2.0;
    let r1 = 3.0 * (3.0 - v.powi(2)) / 2.0
        * ((1.0 + v.powi(2))
            * (dilog(p.powi(2))
                + 2.0 * dilog(p)
                + log_ratio * ((1.0 + v).powi(3) / (8.0 * v.powi(2))).ln())
            + 3.0 * v * ((1.0 - v.powi(2)) / (4.0 * v)).ln()
            - v * v.ln()
            + (33.0 + 22.0 * v.powi(2) - 7.0 * v.powi(4)) / (8.0 * (3.0 - v.powi(2))) * log_ratio
            + (15.0 * v - 9.0 * v.powi(3)) / (4.0 * (3.0 - v.powi(2))));
    // alpha_QED E Q_t^2 / 2, with N_c=3 already included in R_0 and R_1.
    let normalization = e_cm * (4.0 / 9.0) / (2.0 * 132.507);
    (
        normalization * r0,
        normalization * (4.0 / 3.0) * 0.118 / PI * r1,
    )
}

#[test]
fn pole_mass_vector_current_reference_matches_published_nlo_benchmarks() {
    // Independent numerical anchors from table 4 of 2203.11038v1. Its
    // second column is 3000 GeV, as specified by the table header and text;
    // the caption's 360 GeV is a typo. Allow the published numerical precision.
    for (energy, published_nlo, rounding) in [(400.0, 0.432831, 1.0e-6), (3000.0, 0.59047, 5.0e-6)]
    {
        let (lo, nlo) = pole_mass_vector_current_cross_section(energy, 173.0);
        assert!(lo > 0.0 && (nlo - published_nlo).abs() <= rounding);
    }
    let (lo, nlo) = pole_mass_vector_current_cross_section(400.0, 173.0);
    assert!((lo - 1.387_581_507_745_258).abs() <= 1.0e-13);
    assert!((nlo - 0.432_830_086_776_986).abs() <= 1.0e-13);
}

fn gamma_star_ttx_nlo_matches_reference(on_shell: bool, sampling_probe_only: bool) -> Result<()> {
    let _serial_guard = gammaloop_integration_tests::SERIAL_TEST_LOCK.lock();
    const LO_PROCESS: &str = "gamma_star_ttx_lo";
    const MSBAR_NLO_PROCESS: &str = "gamma_star_ttx_explicit_local_3d";
    const MSBAR_ROUTES: [(&str, bool, bool); 3] = [
        ("gamma_star_ttx_orientation_local_3d", false, false),
        (MSBAR_NLO_PROCESS, true, false),
        ("gamma_star_ttx_local_4d_then_cff", true, true),
    ];
    const OS_ROUTES: [(&str, bool, bool); 3] = [
        ("gamma_star_ttx_os_orientation_local_3d", false, false),
        ("gamma_star_ttx_os_explicit_local_3d", true, false),
        ("gamma_star_ttx_os_local_4d_then_cff", true, true),
    ];
    const LO: &str = "LO";
    const NLO: &str = "NLO";
    let e_cm = if on_shell { 400.0 } else { 600.0 };
    let mu_r = if on_shell { 173.0 } else { e_cm };
    let (reference_lo, reference_nlo) = if on_shell {
        pole_mass_vector_current_cross_section(e_cm, 173.0)
    } else {
        (2.876302, 2.01520e-1)
    };
    const PUBLISHED_VERTEX: f64 = 3.48276e-1;
    const PUBLISHED_SELF_ENERGY: f64 = -1.46756e-1;
    let routes: &[(&str, bool, bool)] = if on_shell { &OS_ROUTES } else { &MSBAR_ROUTES };
    let reference_index = if on_shell { 2 } else { 1 };
    let nlo_process = routes[reference_index].0;
    let test_name = if sampling_probe_only {
        "gamma_star_ttx_os_sampling_probes_compare_route_maps"
    } else if on_shell {
        "gamma_star_ttx_os_nlo_matches_the_analytic_inclusive_cross_section_and_auxiliary_scales"
    } else {
        "gamma_star_ttx_msbar_nlo_matches_the_published_absolute_cross_sections"
    };
    let test_root = get_tests_workspace_path().join(test_name);
    clean_test(&test_root);
    let acceptance_artifacts = (on_shell && !sampling_probe_only)
        .then(|| gammaloop_integration_tests::new_test_artifact_dir("gamma_star_ttx_os_acceptance"))
        .transpose()?;

    let mut cli = get_example_cli(
        "epem_a_ddx/NLO/epem_a_ddx_NLO.toml",
        &[],
        Some(test_root.join("state")),
        None,
        true,
    )?;
    {
        let generation = &mut cli.cli_settings.global.generation;
        generation.orientation_pattern = Default::default();
        generation.explicit_orientation_sum_only = true;
        generation
            .tropical_subgraph_table
            .disable_tropical_generation = true;
        generation.threshold_subtraction.enable_thresholds = true;
        generation.threshold_subtraction.disable_integrated_ct = false;
        generation.vector_polarization_sum_gauge = VectorPolarizationSumGauge::Feynman;
        let uv = &mut generation.uv;
        uv.softct = false;
        uv.generate_integrated = true;
        uv.subtract_uv = true;
        uv.final_integrand = FinalIntegrandDimension::ThreeD;
        uv.local_uv_cts_from_expanded_4d_integrands = false;
        uv.renormalization_prescription.log_divergent = ApproximationType::MUV;
        uv.renormalization_prescription.massive_power_divergent = ApproximationType::MUV;
        uv.renormalization_prescription.massless_power_divergent = ApproximationType::MUV;
        uv.renormalization_prescription.overrides.clear();
        if on_shell {
            uv.renormalization_prescription
                .overrides
                .push(CTRenormalizationRule::new(
                    CTIdentifier::new(
                        [-6, 6].into_iter().collect(),
                        Some([6, 21].into_iter().collect()),
                    ),
                    ApproximationType::OS,
                ));
        }
        uv.vakint.normalization = "MSbar".to_string();
        uv.vakint.additional_normalization = "1".to_string();
        uv.vakint.form_exe_path = which("form")?.display().to_string();
    }
    cli.default_runtime_settings
        .subtraction
        .disable_threshold_subtraction = false;

    // Integrand preprocessing applies external spin sums, so the direct-current
    // external state must already be the runtime default before `generate
    // existing` (the reused e+e- card otherwise still supplies two helicities).
    cli.run_command(&format!(
        r#"set default-runtime string '
            [kinematics]
            e_cm = {e_cm}

            [kinematics.externals]
            type = "constant"

            [kinematics.externals.data]
            momenta = [[{e_cm}, 0.0, 0.0, 0.0]]
            helicities = ["summed"]

            [general]
            mu_r = {mu_r}
            m_uv = {mu_r}
            renormalization_localization_scale = {mu_r}
            integral_unit = "none"
            disable_flux_factor = false

            [sampling]
            graphs = "monte_carlo"
            graph_names = []
            orientations = "summed"
            sampling_multichanneling = true
            sampling_channels = "summed"
            sampling_channel_weight = "ose"
        '"#,
    ))?;

    cli.run_command("import model sm-default.json")?;
    cli.run_command("run set_model_parameters")?;
    cli.run_command(&format!(
        r#"generate xs a > t t~ | a t t~ g ghG ghG~ QCD^2==0 QED^2==2 [{{{{1}}}} QCD=0]
            --numerator-grouping group_identical_graphs_up_to_scalar_rescaling
            --symmetrize-left-right-states true
            -p {LO_PROCESS} -i {LO} --global-prefactor-num "1/2" --only-diagrams"#,
    ))?;
    cli.run_command(&format!("generate existing -p {LO_PROCESS} -i {LO}"))?;

    // The paper averages over the two incoming-current polarizations. Keep
    // that factor separate so that `summed` below applies exactly -g^{mu nu}
    // rather than the invalid on-shell photon spin sum.
    for (process, explicit_orientation_sum_only, project_local_4d_to_cff) in routes.iter().copied()
    {
        {
            let generation = &mut cli.cli_settings.global.generation;
            generation.explicit_orientation_sum_only = explicit_orientation_sum_only;
            generation.uv.local_uv_cts_from_expanded_4d_integrands = project_local_4d_to_cff;
        }
        cli.run_command(&format!(
            r#"generate xs a > t t~ | a t t~ g ghG ghG~ QCD^2==2 QED^2==2 [{{{{2}}}} QCD=1]
                --numerator-grouping group_identical_graphs_up_to_scalar_rescaling
                --symmetrize-left-right-states true
                -p {process} -i {NLO} --global-prefactor-num "1/2" --only-diagrams"#,
        ))?;
        cli.run_command(&format!("generate existing -p {process} -i {NLO}"))?;
    }
    let mut role_master_names = BTreeMap::<&str, BTreeMap<&str, Vec<String>>>::new();
    let mut reference_group_inventory = None;
    for (process, _, _) in routes.iter().copied() {
        let process_id = cli
            .state
            .resolve_process_ref(Some(&ProcessRef::Unqualified(process.to_string())))?;
        if on_shell {
            let model_parameters = cli
                .state
                .resolve_effective_model_parameter_card_for_integrand(process_id, NLO)?;
            for (parameter, expected) in [
                ("MT", 173.0),
                ("WT", 0.0),
                ("aS", 0.118),
                ("aEWM1", 132.507),
            ] {
                let (re, im) = model_parameters.data.get(parameter).ok_or_else(|| {
                    eyre!("the effective OS model parameters do not contain {parameter}")
                })?;
                assert_eq!(im.0, 0.0, "the analytic reference needs real {parameter}");
                assert_eq!(
                    re.0, expected,
                    "the analytic reference has different {parameter}"
                );
            }
        }
        let info = cli.state.get_integrand_info(
            Some(&ProcessRef::Unqualified(process.to_string())),
            Some(&NLO.to_string()),
        )?;
        let resolved = cli.state.process_list.processes[process_id].get_integrand(NLO)?;
        let ProcessIntegrand::CrossSection(integrand) = resolved.require_generated()? else {
            return Err(eyre!(
                "the direct-current acceptance must use a generated cross-section"
            ));
        };
        assert_eq!(
            integrand.get_n_dim(),
            6,
            "relative spherical sampling requires two integrated three-momenta"
        );
        let settings = &cli.cli_settings.global.generation.uv;
        let mut physical_cut_particles = BTreeSet::new();
        let mut graph_roles = Vec::new();
        for term in &integrand.data.graph_terms {
            let graph = &term.graph;
            let mut roles = BTreeSet::new();
            let grouped_cuts = term
                .cut_group_data
                .cut_groups
                .iter()
                .flat_map(|group| group.cuts.iter().copied())
                .collect::<Vec<_>>();
            assert_eq!(
                grouped_cuts.len(),
                term.cuts.len(),
                "{} lost a physical cut",
                graph.name
            );
            assert_eq!(
                grouped_cuts
                    .iter()
                    .map(|cut_id| cut_id.0)
                    .collect::<BTreeSet<_>>()
                    .len(),
                term.cuts.len(),
                "{} selected a physical cut more than once",
                graph.name
            );
            for cut in &term.cuts {
                let mut particles = graph
                    .iter_edges_of(&cut.cut.as_subgraph())
                    .map(|(_, _, edge)| match &edge.data.particle {
                        PossibleParticle::Particle(particle)
                        | PossibleParticle::MassOverriddenParticle { particle, .. } => graph
                            .model
                            .particle_by_id(*particle)
                            .expect("the physical cut particle belongs to the graph model")
                            .pdg_code
                            .abs(),
                        PossibleParticle::JustMass { .. } => {
                            panic!("physical cuts must retain their particle identities")
                        }
                    })
                    .collect::<Vec<_>>();
                particles.sort_unstable();
                physical_cut_particles.insert(particles);
            }
            // Match production's raised-cut grouping. Cutting the outer photon
            // correlator leaves at most one virtual loop, so this acceptance
            // needs no nested photon-parent/OS-child Ward policy.
            for group in &term.cut_group_data.cut_groups {
                let cut_union = group
                    .cuts
                    .iter()
                    .map(|cut_id| term.cuts[*cut_id].cut.as_subgraph())
                    .reduce(|left, right| left.union(&right))
                    .ok_or_else(|| eyre!("{} has an empty cut group", graph.name))?;
                let uncut = graph
                    .full_filter()
                    .subtract(&graph.initial_state_cut.left)
                    .subtract(&cut_union);
                for spinney in graph.classified_spinneys(
                    &uncut,
                    &graph.model,
                    settings,
                    &graph.loop_momentum_basis,
                )? {
                    if spinney.n_components() == 0 {
                        continue;
                    }
                    assert_eq!(
                        spinney.n_components(),
                        1,
                        "NLO has no disjoint virtual loops"
                    );
                    assert_eq!(
                        graph.n_loops(spinney.filter()),
                        1,
                        "NLO has no nested virtual loops"
                    );
                    let identifier = graph.ct_identifier(spinney.filter(), &graph.model);
                    if spinney.dod == 1 {
                        assert_eq!(
                            spinney.renormalization_scheme,
                            if on_shell {
                                ApproximationType::OS
                            } else {
                                ApproximationType::MUV
                            }
                        );
                        assert_eq!(
                            identifier,
                            CTIdentifier::new(
                                [-6, 6].into_iter().collect(),
                                Some([6, 21].into_iter().collect())
                            )
                        );
                        roles.insert("self-energy");
                    } else {
                        assert_eq!(spinney.dod, 0);
                        assert_eq!(spinney.renormalization_scheme, ApproximationType::MUV);
                        roles.insert("vertex");
                    }
                }
            }
            assert_eq!(
                roles.len(),
                1,
                "{} must have exactly one physical role from its virtual cuts: {roles:?}",
                graph.name
            );
            graph_roles.push(*roles.first().expect("one physical role was required"));
        }
        assert_eq!(
            physical_cut_particles,
            BTreeSet::from([vec![6, 6], vec![6, 6, 21]]),
            "the inclusive observable requires both virtual ttbar and real ttbar-g cuts"
        );

        // Native grouping may retain several masters for one physical role
        // when its conservative numerator comparison cannot prove equality.
        // Keep all of them and resolve published rows by their physical role.
        let mut masters_by_role = BTreeMap::<&str, Vec<String>>::new();
        let mut grouped_graph_ids = BTreeSet::new();
        let mut master_names = BTreeSet::new();
        let mut group_inventory = BTreeSet::new();
        for group in &info.graph_groups {
            let masters = group
                .graphs
                .iter()
                .filter(|graph| graph.is_master)
                .collect::<Vec<_>>();
            assert_eq!(
                masters.len(),
                1,
                "each native group needs exactly one master"
            );
            let master = masters[0];
            let role = graph_roles[master.graph_id];
            assert!(
                master_names.insert(master.name.clone()),
                "duplicate master {}",
                master.name
            );
            let mut members = BTreeSet::new();
            for member in &group.graphs {
                assert!(
                    grouped_graph_ids.insert(member.graph_id),
                    "a graph belongs to more than one group"
                );
                assert_eq!(
                    member.name,
                    integrand.data.graph_terms[member.graph_id].graph.name
                );
                assert_eq!(
                    graph_roles[member.graph_id], role,
                    "{} mixes physical roles",
                    master.name
                );
                assert!(
                    members.insert(member.name.clone()),
                    "duplicate graph name {}",
                    member.name
                );
            }
            masters_by_role
                .entry(role)
                .or_default()
                .push(master.name.clone());
            group_inventory.insert((role, master.name.clone(), members));
        }
        assert_eq!(
            grouped_graph_ids,
            (0..integrand.data.graph_terms.len()).collect::<BTreeSet<_>>()
        );
        assert_eq!(master_names.len(), info.graph_groups.len());
        assert_eq!(
            masters_by_role.keys().copied().collect::<BTreeSet<_>>(),
            BTreeSet::from(["self-energy", "vertex"])
        );
        for names in masters_by_role.values_mut() {
            names.sort();
        }
        if let Some(reference) = &reference_group_inventory {
            assert_eq!(
                &group_inventory, reference,
                "{process} changed native groups or physical roles between routes"
            );
        } else {
            reference_group_inventory = Some(group_inventory);
        }
        role_master_names.insert(process, masters_by_role);
    }
    let graph_names = role_master_names[nlo_process]
        .values()
        .flatten()
        .cloned()
        .collect::<Vec<_>>();
    if sampling_probe_only {
        let artifacts = gammaloop_integration_tests::new_test_artifact_dir(
            "gamma_star_ttx_os_sampling_probes",
        )?;
        // The first two rows are the largest positive and negative samples of
        // the projected-route M_uv=346 pilot at 320k evaluations (seed 1338).
        // Compare the map boundary before changing adaptive sampling settings.
        let probe_points = Array2::from_shape_vec(
            (3, 6),
            vec![
                0.6698519657474593,
                0.5672641947748677,
                0.8884692830398706,
                0.811372050762004,
                0.5632622014369728,
                0.8820961667765356,
                0.8363521033071066,
                0.8795577560012201,
                0.6799454654814769,
                0.9140268219434564,
                0.8850144423462233,
                0.6983671899069425,
                0.19,
                0.43,
                0.71,
                0.29,
                0.53,
                0.61,
            ],
        )?;
        let mut matrix = Vec::new();
        let mut mappings = Vec::new();
        let mut reference_values: BTreeMap<_, Complex<F<ArbPrec>>> = BTreeMap::new();
        let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
        for (process, _, _) in routes {
            let process_id = cli
                .state
                .resolve_process_ref(Some(&ProcessRef::Unqualified((*process).to_owned())))?;
            let graph_groups = cli.state.process_list.processes[process_id]
                .get_integrand(NLO)?
                .require_generated()?
                .graph_group_master_names()
                .into_iter()
                .map(str::to_owned)
                .collect::<Vec<_>>();
            for graph in
                std::iter::once(None).chain(graph_names.iter().map(|name| Some(name.as_str())))
            {
                let graph_sampling = if graph.is_some() {
                    "monte_carlo"
                } else {
                    "summed"
                };
                // EvaluateSamples consumes full-domain discrete group indices;
                // the integration-only graph_names filter does not select them.
                let discrete_graph = graph
                    .map(|name| {
                        graph_groups
                            .iter()
                            .position(|master| master == name)
                            .map(|index| Array2::from_elem((probe_points.nrows(), 1), index))
                            .ok_or_else(|| eyre!("sampling probe has no graph group {name}"))
                    })
                    .transpose()?;
                cli.run_command(&format!(
                    r#"set process -p {process} -i {NLO} kv
                        general.m_uv=346.0
                        general.renormalization_localization_scale=173.0
                        sampling.graphs="{graph_sampling}"
                        sampling.graph_names=[]
                        sampling.orientations="summed"
                        sampling.sampling_multichanneling=true
                        sampling.sampling_channels="summed"
                        sampling.sampling_channel_weight="map_density"
                        sampling.coordinate_system="relative_spherical""#,
                ))?;
                let result = evaluate_samples_precise(
                    &mut cli.state,
                    &EvaluateSamplesPrecise {
                        process_id: Some(process_id),
                        integrand_name: Some(NLO.to_owned()),
                        use_arb_prec: true,
                        minimal_output: false,
                        return_generated_events: Some(false),
                        momentum_space: false,
                        points: probe_points.view(),
                        integrator_weights: None,
                        discrete_dims: discrete_graph.as_ref().map(|indices| indices.view()),
                        graph_names: None,
                        orientations: None,
                    },
                )?;
                for (point, sample) in result.samples.into_iter().enumerate() {
                    let PreciseEvaluationResultOutput::Arb(value) = sample.evaluation else {
                        return Err(eyre!("sampling probes require Arb precision"));
                    };
                    let jacobian = value.parameterization_jacobian.as_ref().ok_or_else(|| {
                        eyre!("x-space sampling probe needs its remaining Jacobian")
                    })?;
                    let weighted = value.integrand_result.clone() * jacobian.clone();
                    assert!(weighted.re.is_finite() && weighted.im.is_finite());
                    let key = (graph.map(str::to_owned), point);
                    if let Some(reference) = reference_values.get(&key) {
                        let distance = (weighted.clone() - reference).norm().re;
                        let actual_norm = weighted.norm().re;
                        let reference_norm = reference.norm().re;
                        let scale = if actual_norm > reference_norm {
                            actual_norm
                        } else {
                            reference_norm
                        };
                        assert!(
                            distance <= &tolerance * scale,
                            "sampling probe {process} {graph:?} point {point} differs between routes"
                        );
                    } else {
                        reference_values.insert(key, weighted.clone());
                    }
                    matrix.push(serde_json::json!({
                        "process": process,
                        "graph": graph,
                        "point": point,
                        "integrand": format!("{:e}", value.integrand_result),
                        "parameterization_jacobian": format!("{jacobian:e}"),
                        "weighted_re": format!("{:e}", weighted.re),
                        "weighted_im": format!("{:e}", weighted.im),
                        "metadata": value.evaluation_metadata,
                    }));
                }
                if graph.is_none() {
                    let resolved =
                        cli.state.process_list.processes[process_id].get_integrand(NLO)?;
                    let ProcessIntegrand::CrossSection(integrand) = resolved.require_generated()?
                    else {
                        return Err(eyre!("sampling probes need a cross-section"));
                    };
                    for term in &integrand.data.graph_terms {
                        let bridge = term.multi_channeling_setup.sampling_bridge::<ArbPrec>()?;
                        for (channel, entry) in bridge.channels().iter().enumerate() {
                            for (point, row) in probe_points.outer_iter().enumerate() {
                                // These f64 coordinates are the exact same explicit
                                // API input boundary used by evaluate_samples_precise.
                                let xs = row
                                    .iter()
                                    .map(|x| F::<ArbPrec>::from_f64(*x).0)
                                    .collect::<Vec<_>>();
                                let mapped = bridge.forward(channel.into(), &xs)?;
                                mappings.push(serde_json::json!({
                                    "process": process,
                                    "graph": term.graph.name,
                                    "point": point,
                                    "channel": channel,
                                    "name": entry.name,
                                    "basis_id": entry.basis_id,
                                    "embedded_edges": entry.embedded_edges,
                                    "raw_coordinates": mapped.raw_coordinates.iter().map(|x| format!("{:e}", F(x.clone()))).collect::<Vec<_>>(),
                                    "jacobian": format!("{:e}", F(mapped.map.jacobian)),
                                    "partition": mapped.partition.weights.iter().map(|x| format!("{:e}", F(x.clone()))).collect::<Vec<_>>(),
                                }));
                            }
                        }
                    }
                }
            }
        }
        for point in 0..probe_points.nrows() {
            let total = &reference_values[&(None, point)];
            let graph_sum = graph_names
                .iter()
                .map(|name| reference_values[&(Some(name.clone()), point)].clone())
                .reduce(|sum, value| sum + value)
                .expect("the physical roles contain native masters");
            let distance = (graph_sum - total).norm().re;
            assert!(
                distance <= &tolerance * total.norm().re,
                "sampling probe point {point} does not sum every native graph group exactly once"
            );
        }
        std::fs::write(
            artifacts.join("route_map_matrix.json"),
            serde_json::to_vec_pretty(&serde_json::json!({
                "physical_role_masters": role_master_names,
                "points": probe_points.rows().into_iter().map(|row| row.to_vec()).collect::<Vec<_>>(),
                "weighted_values": matrix,
                "maps": mappings,
            }))?,
        )?;
        clean_test(&test_root);
        return Ok(());
    }
    {
        let generation = &mut cli.cli_settings.global.generation;
        generation.explicit_orientation_sum_only = true;
        generation.uv.local_uv_cts_from_expanded_4d_integrands = false;
    }

    cli.run_command(&format!(
        r#"set process -p {LO_PROCESS} -i {LO} kv
            sampling.graphs="summed"
            sampling.orientations="summed"
            sampling.sampling_multichanneling=true
            sampling.sampling_channels="summed"
            sampling.sampling_channel_weight="ose""#,
    ))?;
    let (_, lo_probe) = Inspect {
        process: Some(ProcessRef::Unqualified(LO_PROCESS.to_string())),
        integrand_name: Some(LO.to_string()),
        point: vec![0.19, 0.43, 0.71],
        ..Default::default()
    }
    .run(&mut cli.state)?;
    assert!(
        lo_probe.re > 0.0 && lo_probe.im.abs() <= 1.0e-12 * lo_probe.re,
        "the direct-current LO must be positive and real, got {lo_probe:e}",
    );
    cli.run_command(&format!(
        r#"set process -p {LO_PROCESS} -i {LO} kv
            integrator.integrated_phase="real"
            integrator.min_samples_for_update=5000
            integrator.n_start=5000
            integrator.n_increase=5000
            integrator.n_max=10000
            integrator.target_relative_accuracy=1.0e-12
            integrator.seed=7331
            integrator.discrete_dim_learning_rate=0.0
            integrator.continuous_dim_learning_rate=0.1"#,
    ))?;
    let lo_output = Integrate {
        process: vec![ProcessRef::Unqualified(LO_PROCESS.to_string())],
        integrand_name: vec![LO.to_string()],
        n_cores: Some(10),
        workspace_path: Some(test_root.join("lo_integration_workspace")),
        restart: true,
        renderer: RendererOption::Tabled,
        show_max_weight_info: false,
        no_stream_iterations: true,
        no_stream_updates: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;
    let lo_estimate = lo_output
        .single_slot_integral()
        .ok_or_else(|| eyre!("expected one direct-current LO integration slot"))?;
    let lo_value = lo_estimate.result.re.0;
    let lo_error = lo_estimate.error.re.0.abs();
    assert!(
        lo_estimate.result.im.0.abs()
            <= 3.0 * lo_estimate.error.im.0.abs() + 1.0e-12 * reference_lo
            && lo_estimate.error.im.0.abs() <= 0.05 * reference_lo,
        "direct-current LO must have a vanishing imaginary component: {:e} ± {:e}",
        lo_estimate.result,
        lo_estimate.error,
    );
    assert!(
        lo_value > 0.0 && lo_error / lo_value <= 0.05,
        "direct-current LO uncertainty is too large: {lo_value:e} ± {lo_error:e}",
    );
    assert!(
        (lo_value - reference_lo).abs() <= 3.0 * lo_error,
        "direct-current LO mismatch: {lo_value:e} ± {lo_error:e}, published={reference_lo:e}",
    );

    for (process, _, _) in routes {
        cli.run_command(&format!(
            r#"set process -p {process} -i {NLO} kv
                sampling.graphs="summed"
                sampling.orientations="summed"
                sampling.sampling_multichanneling=true
                sampling.sampling_channels="summed"
                sampling.sampling_channel_weight="ose""#,
        ))?;
    }
    let phase_process = nlo_process;
    let (_, probe) = Inspect {
        process: Some(ProcessRef::Unqualified(phase_process.to_string())),
        integrand_name: Some(NLO.to_string()),
        point: vec![17.0, -31.0, 43.0, -29.0, 53.0, 71.0],
        momentum_space: true,
        ..Default::default()
    }
    .run(&mut cli.state)?;
    assert!(
        probe.re.is_finite()
            && probe.im.is_finite()
            && probe.re.abs() > 0.0
            && probe.im.abs() <= 1.0e-8 * probe.re.abs(),
        "the {phase_process} NLO pointwise probe must be finite, nonzero and real, got {probe:e}",
    );

    let points = Array2::from_shape_fn((graph_names.len(), 6), |(_, coordinate)| {
        [17.0, -31.0, 43.0, -29.0, 53.0, 71.0][coordinate]
    });
    // Retaining one eighth of ArbPrec's requested precision permits substantial
    // expression-dependent loss while still testing the precision-increase path.
    let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
    // The phase-space point originates as f64, but all routes receive identical
    // values. Nonzero results still require relative agreement after Arb upcast.
    let mut route_values = Vec::new();
    for (process, _, _) in routes.iter().copied() {
        let process_id = cli
            .state
            .resolve_process_ref(Some(&ProcessRef::Unqualified(process.to_string())))?;
        let result = evaluate_samples_precise(
            &mut cli.state,
            &EvaluateSamplesPrecise {
                process_id: Some(process_id),
                integrand_name: Some(NLO.to_string()),
                use_arb_prec: true,
                minimal_output: true,
                return_generated_events: Some(false),
                momentum_space: true,
                points: points.view(),
                integrator_weights: None,
                discrete_dims: None,
                graph_names: Some(graph_names.iter().cloned().map(Some).collect()),
                orientations: Some(vec![None; graph_names.len()]),
            },
        )?;
        let values = result
            .samples
            .into_iter()
            .zip(&graph_names)
            .map(|(sample, graph_name)| match sample.evaluation {
                PreciseEvaluationResultOutput::Arb(result) => Ok(result.integrand_result),
                evaluation => Err(eyre!(
                    "{process} {graph_name} comparison requested Arb but got {:?}",
                    evaluation.precision()
                )),
            })
            .collect::<Result<Vec<_>>>()?;
        route_values.push((process, values));
    }
    let reference = &route_values[reference_index];
    for (graph_index, graph_name) in graph_names.iter().enumerate() {
        let reference_value = &reference.1[graph_index];
        if on_shell {
            assert!(
                reference_value.re.is_finite() && reference_value.im.is_finite(),
                "OS {graph_name} must evaluate finitely in Arb precision"
            );
        }
        for actual in route_values
            .iter()
            .enumerate()
            .filter_map(|(index, route)| (index != reference_index).then_some(route))
        {
            let actual_value = &actual.1[graph_index];
            let distance = (actual_value.clone() - reference_value.clone()).norm().re;
            let actual_norm = actual_value.norm().re;
            let reference_norm = reference_value.norm().re;
            let scale = if actual_norm > reference_norm {
                actual_norm
            } else {
                reference_norm
            };
            let relative_distance = if IsZero::is_zero(&scale) {
                distance.clone()
            } else {
                distance.clone() / scale.clone()
            };
            assert!(
                [actual_value, reference_value]
                    .into_iter()
                    .all(|value| { value.re.is_finite() && value.im.is_finite() })
                    && relative_distance <= tolerance,
                "{graph_name} differs between precise routes {} and {}: actual={:e}, reference={:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}",
                actual.0,
                reference.0,
                actual_value,
                reference_value,
            );
        }
    }

    // For MSbar, sum every native master belonging to each published physical
    // row through the existing runtime graph-subset selection. The
    // OS acceptance integrates every route at each auxiliary-scale choice;
    // projected 4D runs first so its established seeds remain unchanged.
    let integration_processes = if on_shell {
        routes.iter().rev().map(|route| route.0).collect::<Vec<_>>()
    } else {
        vec![nlo_process]
    };
    let mut graph_estimates = BTreeMap::new();
    let integration_cases = if on_shell {
        vec![
            ("baseline", 173.0, 173.0),
            ("m_uv", 346.0, 173.0),
            ("localization", 173.0, 346.0),
        ]
    } else {
        vec![("self-energy", e_cm, e_cm), ("vertex", e_cm, e_cm)]
    };
    let mut imaginary_sum = 0.0;
    let mut imaginary_variance = 0.0;
    let mut measurements = Vec::new();
    for (route_index, graph_process) in integration_processes.iter().copied().enumerate() {
        for (graph_index, (graph, m_uv, localization_scale)) in
            integration_cases.iter().copied().enumerate()
        {
            let graph_sampling = if on_shell { "summed" } else { "monte_carlo" };
            let graph_names = if on_shell {
                "[]".to_string()
            } else {
                serde_json::to_string(&role_master_names[graph_process][graph])?
            };
            // OS counterterms introduce soft points that differ from the
            // original massive propagator energies. Use each LMB map's
            // actual density to partition these overlapping soft regions.
            let channel_weight = if on_shell { "map_density" } else { "ose" };
            // The largest OS pilot weights occur at k1 ~ k2. Sampling the
            // common and relative momenta supplies the relative r^2 dr
            // Jacobian explicitly, including this correlated soft region.
            let coordinate_system = if on_shell {
                "relative_spherical"
            } else {
                "spherical"
            };
            // Keep the sample budget fixed while testing stronger Vegas
            // adaptation of the measured correlated soft regions. Independent
            // seeds and the 5% reference-relative bound remain unchanged.
            let n_max = if on_shell { 320_000 } else { 40_000 };
            let continuous_learning_rate = if on_shell { 0.25 } else { 0.1 };
            // Different seeds keep route/scale comparisons statistically
            // independent; the projected baseline retains seed 1337.
            let seed = 1337 + route_index * integration_cases.len() + graph_index;
            cli.run_command(&format!(
                r#"set process -p {graph_process} -i {NLO} kv
                    general.m_uv={m_uv}
                    general.renormalization_localization_scale={localization_scale}
                    sampling.graphs="{graph_sampling}"
                    sampling.graph_names={graph_names}
                    sampling.orientations="summed"
                    sampling.sampling_multichanneling=true
                    sampling.sampling_channels="summed"
                    sampling.sampling_channel_weight="{channel_weight}"
                    sampling.coordinate_system="{coordinate_system}"
                    integrator.integrated_phase="real"
                    integrator.min_samples_for_update=20000
                    integrator.n_start=20000
                    integrator.n_increase=20000
                    integrator.n_max={n_max}
                    integrator.target_relative_accuracy=1.0e-12
                    integrator.seed={seed}
                    integrator.discrete_dim_learning_rate=0.0
                    integrator.continuous_dim_learning_rate={continuous_learning_rate}"#,
            ))?;
            let started = std::time::Instant::now();
            let output = Integrate {
                process: vec![ProcessRef::Unqualified(graph_process.to_string())],
                integrand_name: vec![NLO.to_string()],
                n_cores: Some(10),
                workspace_path: Some(
                    test_root.join(format!("{graph_process}_{graph}_integration_workspace")),
                ),
                restart: true,
                renderer: RendererOption::Tabled,
                show_max_weight_info: false,
                no_stream_iterations: true,
                no_stream_updates: true,
                ..Default::default()
            }
            .run(&mut cli.state, &cli.cli_settings)?;
            let elapsed_seconds = started.elapsed().as_secs_f64();
            let estimate = output
                .single_slot_integral()
                .ok_or_else(|| eyre!("expected one {graph_process} {graph} integration slot"))?;
            assert!(
                estimate.result.im.0.abs()
                    <= 3.0 * estimate.error.im.0.abs() + 1.0e-12 * reference_nlo
                    && estimate.error.im.0.abs() <= 0.15 * reference_nlo,
                "{graph_process} {graph} must have a vanishing imaginary component: {:e} ± {:e}",
                estimate.result,
                estimate.error,
            );
            imaginary_sum += estimate.result.im.0;
            imaginary_variance += estimate.error.im.0.powi(2);
            let value = estimate.result.re.0;
            let error = estimate.error.re.0.abs();
            graph_estimates.insert((graph_process, graph), (value, error));
            measurements.push(serde_json::json!({
                "process": graph_process,
                "case": graph,
                "m_uv": m_uv,
                "renormalization_localization_scale": localization_scale,
                "seed": seed,
                "n_max": n_max,
                "sampling_channel_weight": channel_weight,
                "coordinate_system": coordinate_system,
                "continuous_dim_learning_rate": continuous_learning_rate,
                "value": value,
                "error": error,
                "imaginary_value": estimate.result.im.0,
                "imaginary_error": estimate.error.im.0.abs(),
                "elapsed_seconds": elapsed_seconds,
            }));
        }
    }
    if on_shell {
        let artifact_root = acceptance_artifacts
            .as_ref()
            .expect("the OS run has a dedicated acceptance artifact directory");
        std::fs::write(
            artifact_root.join("results.json"),
            serde_json::to_vec_pretty(&serde_json::json!({
                "schema_version": 1,
                "measurement_process_id": std::process::id(),
                "e_cm": e_cm,
                "pole_mass": 173.0,
                "mu_r": mu_r,
                "analytic_lo": reference_lo,
                "analytic_nlo": reference_nlo,
                "lo": {"value": lo_value, "error": lo_error},
                "physical_role_masters": role_master_names,
                "measurements": measurements,
            }))?,
        )?;
        // Auxiliary masses and localization kernels may change the integrand.
        // Each complete integrated sum must satisfy the analytic reference,
        // its own scale scan, and the other routes at the same scales.
        for (process, _, _) in routes.iter().copied() {
            let (baseline, baseline_error) = graph_estimates[&(process, "baseline")];
            for (scale_choice, _, _) in integration_cases.iter().copied() {
                let (value, error) = graph_estimates[&(process, scale_choice)];
                assert!(
                    error <= 0.05 * reference_nlo,
                    "OS {process} {scale_choice} uncertainty is too large: {value:e} +/- {error:e}"
                );
                assert!(
                    (value - reference_nlo).abs() <= 3.0 * error,
                    "OS {process} {scale_choice} disagrees with the inclusive pole-mass R-ratio: {value:e} +/- {error:e}, analytic={reference_nlo:e}"
                );
                if scale_choice != "baseline" {
                    let combined_error = error.hypot(baseline_error);
                    assert!(
                        (value - baseline).abs() <= 3.0 * combined_error,
                        "OS {process} {scale_choice} changes the integrated result: {value:e} versus {baseline:e}, combined error={combined_error:e}"
                    );
                }
                if process != nlo_process {
                    let (reference, reference_error) =
                        graph_estimates[&(nlo_process, scale_choice)];
                    let combined_error = error.hypot(reference_error);
                    assert!(
                        (value - reference).abs() <= 3.0 * combined_error,
                        "OS {scale_choice} differs between {process} and {nlo_process}: {value:e} versus {reference:e}, combined error={combined_error:e}"
                    );
                }
            }
        }
        clean_test(&test_root);
        return Ok(());
    }
    for (graph, target) in [
        ("self-energy", PUBLISHED_SELF_ENERGY),
        ("vertex", PUBLISHED_VERTEX),
    ] {
        let (value, error) = graph_estimates[&(nlo_process, graph)];
        assert!(
            error / target.abs() <= 0.15,
            "{graph} uncertainty is too large: {value:e} ± {error:e}, published={target:e}",
        );
        assert!(
            (value - target).abs() <= 3.0 * error,
            "{graph} mismatch: {value:e} ± {error:e}, published={target:e}",
        );
    }
    assert!(
        imaginary_sum.abs() <= 3.0 * imaginary_variance.sqrt() + 1.0e-12 * reference_nlo,
        "the inclusive NLO imaginary component must vanish: {imaginary_sum:e} ± {:e}",
        imaginary_variance.sqrt(),
    );
    let graph_sum = graph_estimates
        .values()
        .map(|(value, _)| value)
        .sum::<f64>();
    let graph_sum_error = graph_estimates
        .values()
        .map(|(_, error)| error.powi(2))
        .sum::<f64>()
        .sqrt();
    assert!(
        (graph_sum - reference_nlo).abs() <= 3.0 * graph_sum_error,
        "sum of graphwise NLO rows mismatches the paper: {graph_sum:e} ± {graph_sum_error:e}, published={reference_nlo:e}",
    );

    clean_test(&test_root);
    Ok(())
}
