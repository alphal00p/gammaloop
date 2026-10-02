use super::*;
use crate::{
    initialisation::test_initialise,
    integrands::process::{
        GaussianReferenceFunction, ProcessIntegrand, SamplingChannelBridgeAcceptanceReport,
        SamplingSupport,
    },
    processes::Amplitude,
    utils::load_generic_model,
};

#[test]
fn fermi_sampling_generated_amplitude_preserves_geometry_and_reference() -> Result<()> {
    test_initialise()?;
    let mut model = load_generic_model("sm");
    model.get_parameter_mut("MT")?.value = Some(Complex::new_re(F(0.5)));
    model.get_parameter_mut("muB")?.value = Some(Complex::new_re(F(3.0)));
    model.recompute_dependents()?;
    let graphs = Graph::from_string(
        r#"digraph fermi_sunrise {
            node [num=1]; edge [num=1];
            A -> B [id=0, particle="t", lmb_id=0];
            B -> A [id=1, particle="t", lmb_id=1];
            A -> B [id=2, particle="g"];
        }"#,
        &model,
    )?;
    let mut global: GlobalSettings = toml::from_str(
        r#"
[generation]
override_lmb_heuristics = true
[generation.uv]
subtract_uv = false
generate_integrated = false
[generation.threshold_subtraction]
enable_thresholds = false
[generation.tropical_subgraph_table]
disable_tropical_generation = true
[generation.evaluator]
compile = false
summed = false
summed_function_map = true
iterative_orientation_optimization = false
"#,
    )?;
    let settings: RuntimeSettings = toml::from_str(
        r#"
[general]
evaluator_method = "SummedFunctionMap"
[kinematics]
e_cm = 1.0
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = []
helicities = []
[sampling]
graphs = "summed"
orientations = "summed"
sampling_multichanneling = true
sampling_channels = "summed"
sampling_channel_weight = "map_density"
power = 2.0
default_channel_selection = ["lmb(0,1)"]
[sampling.channel_definitions.fermi_sunrise.pair]
around = "product(block(lmb(0),fermi(0)),block(lmb(1),fermi(1)))"
parent_lmb = [0,1]
"#,
    )?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .stack_size(256 * 1024 * 1024)
        .build()?;
    for (mode, explicit) in [
        (MediumMode::ThermodynamicEquilibrium, false),
        (MediumMode::ThermodynamicEquilibrium, true),
        (MediumMode::ZeroTemperatureEquilibrium, false),
        (MediumMode::Vacuum, false),
    ] {
        global.generation.medium.mode = mode;
        global.generation.explicit_orientation_sum_only = explicit;
        let mut amplitude = Amplitude::from_graph_list("fermi_sunrise", graphs.clone())?;
        amplitude.preprocess(&model, &global.generation, &(&settings).into(), &pool)?;
        amplitude.build_integrand(&model, "fermi_sunrise", &global, (&settings).into(), &pool)?;
        let runtime = amplitude.integrand.as_mut().unwrap();
        let mut parameterization = settings.sampling.get_parameterization_settings().unwrap();
        parameterization.sampling_channels.default_channel_selection =
            vec!["pair".to_owned(), "lmb(0,1)".to_owned()];
        let mut selected = settings.clone();
        let SamplingSettings::MultiChanneling(channels) = &mut selected.sampling else {
            unreachable!("fixture requests summed graph and channel sampling");
        };
        channels.parameterization_settings = parameterization.clone();
        *runtime.get_mut_settings() = selected.clone();
        if mode == MediumMode::Vacuum {
            let error = runtime.warm_up(&model).unwrap_err();
            assert!(
                format!("{error:#}").contains("active thermal edge"),
                "{error:#}"
            );
            continue;
        }
        runtime.warm_up(&model)?;
        let ProcessIntegrand::Amplitude(generated) = &*runtime else {
            unreachable!()
        };
        let term = &generated.data.graph_terms[0];
        for edge in [EdgeIndex(0), EdgeIndex(1)] {
            assert!(term.thermal_edges.iter().any(|edges| edges.contains(&edge)));
        }
        let bridge = term.sampling_setup().sampling_bridge::<f64>()?;
        let cube = [0.23, 0.31, 0.61, 0.72, 0.27, 0.57];
        let mapped = bridge.forward(SamplingChannelId(0), &cube)?;
        assert_eq!(mapped.map.support, SamplingSupport::Full);
        let inverse = bridge
            .inverse(SamplingChannelId(0), &mapped.raw_coordinates)?
            .unwrap();
        assert!((mapped.map.jacobian * inverse.map.inverse_jacobian - 1.0).abs() < 1.0e-10);
        assert!((mapped.partition.weights.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
        // The component oracle fixes the physical shell, independently of
        // metadata extraction and graph-to-master binding.
        let shell = SurfaceRadialMap::new(3, vec![0.0; 3], Some(0.75_f64.sqrt()), 1.0, 2.0)?;
        for (point, coordinates) in mapped
            .raw_coordinates
            .chunks_exact(3)
            .zip(cube.chunks_exact(3))
        {
            let expected = shell.forward(&coordinates.iter().copied().map(F).collect_vec())?;
            for (actual, expected) in point.iter().zip(expected.point) {
                assert!(
                    (actual - expected.0).abs() < 1.0e-11,
                    "{mode:?}, explicit={explicit}"
                );
            }
        }
        if mode == MediumMode::ThermodynamicEquilibrium && !explicit {
            let report = SamplingChannelBridgeAcceptanceReport::normalized_gaussian(
                bridge,
                8192,
                0.8,
                &[0.1, -0.2, 0.3, -0.1, 0.2, -0.3],
            )?;
            assert!((report.normalization - 1.0).abs() < 0.025, "{report:?}");
            assert!(
                (report.second_moment / report.expected_second_moment - 1.0).abs() < 0.035,
                "{report:?}"
            );
            assert_eq!(
                report.finite_sample_count,
                report.sample_count * report.channel_count
            );
        }
        let reference = GaussianReferenceFunction::centered(0.8, 2)?;
        let mut expected = 0.0;
        let mut expected_moment = 0.0;
        for channel in 0..bridge.channels().len() {
            let point = bridge.forward(SamplingChannelId(channel), &cube)?;
            let norm_squared = point.raw_coordinates.iter().map(|x| x * x).sum::<f64>();
            let gaussian = (-norm_squared / (2.0 * 0.8_f64.powi(2))).exp()
                / (2.0 * std::f64::consts::PI * 0.8_f64.powi(2)).powi(3);
            let weight = gaussian * point.map.jacobian * point.partition.weights[channel];
            expected += weight;
            expected_moment += weight * norm_squared;
        }
        let source = Sample::Continuous(F(1.0), cube.map(F).to_vec());
        let evaluated = runtime.evaluate_reference_sample_detailed(&source, &reference)?;
        let jacobian = evaluated.evaluation.parameterization_jacobian.unwrap().0;
        assert!(!evaluated.evaluation.evaluation_metadata.is_nan);
        assert!(
            (evaluated.evaluation.integrand_result.re.0 * jacobian / expected - 1.0).abs() < 1.0e-9
        );
        assert!(
            (evaluated.moments.second_moment.0 * jacobian / expected_moment - 1.0).abs() < 1.0e-9
        );
        if mode != MediumMode::ThermodynamicEquilibrium || explicit {
            continue;
        }
        // A group samples its master's proposal even if another member has no
        // surviving thermal factors on those edges. These generated clones test
        // cache ownership and reference normalization, not a new physical graph sum.
        {
            let ProcessIntegrand::Amplitude(generated) = &*runtime else {
                unreachable!()
            };
            let mut grouped = generated.clone();
            grouped.data.graph_terms[0].graph.group_id = Some(GroupId(0));
            grouped.data.graph_terms[0].graph.is_group_master = true;
            let mut member = grouped.data.graph_terms[0].clone();
            member.graph.name = "fermi_sunrise_member".to_owned();
            member.graph.is_group_master = false;
            for edges in member.thermal_edges.iter_mut() {
                edges.clear();
            }
            let error = member
                .compile_sampling_bridge::<f64>(&parameterization, &selected, &[], None)
                .unwrap_err();
            assert!(
                format!("{error:#}").contains("active thermal edge"),
                "{error:#}"
            );
            grouped.data.graph_terms.push(member);
            let mut graphs = grouped
                .data
                .graph_terms
                .iter()
                .map(|term| term.graph.clone())
                .collect_vec();
            grouped.data.graph_group_structure =
                crate::graph::parse::complete_group_parsing(&mut graphs)?;
            grouped.data.graph_to_group_id = vec![0, 0];
            grouped.warm_up_sampling()?;
            assert_eq!(grouped.data.graph_group_structure.len(), 1);
            assert!(
                grouped.data.graph_terms[1]
                    .thermal_edges
                    .iter()
                    .all(Vec::is_empty)
            );
            let master = grouped
                .get_master_graph(GroupId(0))
                .sampling_setup()
                .sampling_bridge::<f64>()?;
            for channel in 0..master.channels().len() {
                let expected = master.forward(SamplingChannelId(channel), &cube)?;
                for term in &grouped.data.graph_terms {
                    let actual = term
                        .sampling_setup()
                        .sampling_bridge::<f64>()?
                        .forward(SamplingChannelId(channel), &cube)?;
                    assert_eq!(actual.raw_coordinates, expected.raw_coordinates);
                    assert_eq!(actual.map.jacobian, expected.map.jacobian);
                    assert_eq!(actual.partition.weights, expected.partition.weights);
                }
            }
            let mut grouped = ProcessIntegrand::Amplitude(grouped);
            let evaluated = grouped.evaluate_reference_sample_detailed(&source, &reference)?;
            let jacobian = evaluated.evaluation.parameterization_jacobian.unwrap().0;
            assert!(!evaluated.evaluation.evaluation_metadata.is_nan);
            assert!(
                (evaluated.evaluation.integrand_result.re.0 * jacobian / expected - 1.0).abs()
                    < 1.0e-9
            );
            assert!(
                (evaluated.moments.second_moment.0 * jacobian / expected_moment - 1.0).abs()
                    < 1.0e-9
            );
        }
        // Changing model parameters must rebuild the same catalogue, including
        // the zero-radius onset and negative-mu orientation union.
        for mu in [0.5_f64, 0.4, 0.0, -1.0, 1.0] {
            model.get_parameter_mut("muB")?.value = Some(Complex::new_re(F(3.0 * mu)));
            model.recompute_dependents()?;
            runtime.warm_up(&model)?;
            let ProcessIntegrand::Amplitude(generated) = &*runtime else {
                unreachable!()
            };
            let bridge = generated.data.graph_terms[0]
                .sampling_setup()
                .sampling_bridge::<f64>()?;
            assert_eq!(bridge.channels().len(), 2);
            let point = bridge.forward(SamplingChannelId(0), &cube)?;
            let radius = (mu.abs() > 0.5).then(|| ((mu.abs() - 0.5) * (mu.abs() + 0.5)).sqrt());
            let oracle = SurfaceRadialMap::new(3, vec![0.0; 3], radius, 1.0, 2.0)?;
            for (point, coordinates) in point
                .raw_coordinates
                .chunks_exact(3)
                .zip(cube.chunks_exact(3))
            {
                let expected = oracle.forward(&coordinates.iter().copied().map(F).collect_vec())?;
                for (actual, expected) in point.iter().zip(expected.point) {
                    assert!((actual - expected.0).abs() < 1.0e-11, "mu={mu}");
                }
            }
            if mu.abs() > 0.5 {
                // Selecting one physical branch must not borrow the shell of
                // its opposite sign, including when runtime filtering precedes
                // a summed proposal. The model already owns the sign of mu.
                let mut term = generated.data.graph_terms[0].clone();
                for direction in [Orientation::Default, Orientation::Reversed] {
                    let id = term
                        .orientations
                        .iter_enumerated()
                        .find_map(|(id, orientation)| {
                            (orientation[EdgeIndex(0)] == direction
                                && term.thermal_edges[id].contains(&EdgeIndex(0)))
                            .then_some(id)
                        })
                        .expect("both physical directions retain a thermal edge");
                    let signed_mu = if direction == Orientation::Default {
                        mu
                    } else {
                        -mu
                    };
                    let oracle = SurfaceRadialMap::new(
                        3,
                        vec![0.0; 3],
                        if signed_mu > 0.5 { radius } else { None },
                        1.0,
                        2.0,
                    )?;
                    let expected =
                        oracle.forward(&cube[..3].iter().copied().map(F).collect_vec())?;
                    for filtered_sum in [false, true] {
                        term.orientation_filter = if filtered_sum {
                            let mut filter = SubSet::empty(term.orientations.len());
                            filter.add(id);
                            filter
                        } else {
                            SubSet::full(term.orientations.len())
                        };
                        let bridge = term.compile_sampling_bridge::<f64>(
                            &parameterization,
                            &selected,
                            &[],
                            (!filtered_sum).then_some(usize::from(id)),
                        )?;
                        let actual = bridge.forward(SamplingChannelId(0), &cube)?;
                        assert_eq!(actual.map.support, SamplingSupport::Full);
                        for (actual, expected) in
                            actual.raw_coordinates[..3].iter().zip(&expected.point)
                        {
                            assert!(
                                (actual - expected.0).abs() < 1.0e-11,
                                "mu={mu}, direction={direction:?}, filtered_sum={filtered_sum}",
                            );
                        }
                    }
                }
            }
        }
        let mut complex_model = model.clone();
        complex_model.get_parameter_mut("MT")?.value = Some(Complex::new(F(0.5), F(0.1)));
        let error = runtime.warm_up(&complex_model).unwrap_err();
        assert!(format!("{error:#}").contains("real mass"), "{error:#}");
        runtime.warm_up(&model)?;
        complex_model = model.clone();
        complex_model.get_parameter_mut("mut")?.value = Some(Complex::new(F(1.0), F(0.1)));
        let error = runtime.warm_up(&complex_model).unwrap_err();
        assert!(
            format!("{error:#}").contains("real chemical potential"),
            "{error:#}"
        );
        runtime.warm_up(&model)?;
    }
    Ok(())
}
