use color_eyre::Result;
use linnet::half_edge::involution::EdgeIndex;
use rayon::ThreadPoolBuilder;
use symbolica::evaluate::OptimizationSettings;
use symbolica::prelude::*;
use three_dimensional_reps::thermal::OCCUPATION_FORM;
use three_dimensional_reps::{MediumMode, ThermalDistributionCall, ThermalDistributionFactor};

use crate::cff::PhaseQuenchedDifference;
use crate::dot;
use crate::graph::Graph;
use crate::graph::parse::from_dot::IntoGraph;
use crate::initialisation;
use crate::integrands::process::MomentumSpaceEvaluationInput;
use crate::momentum::ThreeMomentum;
use crate::processes::Amplitude;
use crate::settings::{GlobalSettings, RuntimeSettings};
use crate::utils;
use crate::utils::{F, GS};
use crate::uv::uv_graph::UVE;

use super::FermiSurfaceSector;

struct TwoQuarkLoops {
    graph: Graph,
}

impl TwoQuarkLoops {
    fn new(mode: MediumMode) -> Result<Self> {
        initialisation::test_initialise()?;
        // The gluon bridge connects two same-flavor cycles without coupling
        // their Fermi momenta. Every declared numerator is local and constant.
        let mut graph: Graph = dot!(digraph phase_quenched_fermi_cycles {
            node [num=1]
            edge [num=1 particle="d"]
            A -> B [id=0 lmb_id=0]
            B -> A [id=1]
            C -> D [id=2 lmb_id=1]
            D -> C [id=3]
            A -> C [id=4 particle="g"]
        })?;
        graph.set_medium_mode(mode)?;
        Ok(Self { graph })
    }

    fn weight(order: i64, temperature: i64, first: i32, second: i32) -> Atom {
        let derivative = ThermalDistributionCall::atom(0, order, temperature, 0, first);
        let steps = ThermalDistributionCall::atom(2, 0, temperature, 1, second)
            + Atom::num(2) * ThermalDistributionCall::atom(3, 0, temperature, -1, -second);
        function!(GS.thermal_weight_wrapper, derivative * steps)
    }
}

#[test]
fn phase_quenched_fermi_sectors_preserve_signed_supports_and_coefficients() -> Result<()> {
    let fixture = TwoQuarkLoops::new(MediumMode::ZeroTemperatureEquilibrium)?;
    let graph = &fixture.graph;
    let projection = PhaseQuenchedDifference::new(graph)?;
    let chemical_potential = graph[EdgeIndex(0)].chemical_potential_atom().unwrap();
    assert_eq!(
        graph[EdgeIndex(2)].chemical_potential_atom().unwrap(),
        chemical_potential
    );
    let params = [
        GS.ose(EdgeIndex(0)),
        GS.ose(EdgeIndex(2)),
        GS.ose(EdgeIndex(3)),
        chemical_potential,
    ];
    // This smooth factor remains outside the complete thermal weight.
    let coefficient = GS.ose(EdgeIndex(0)).pow(2);
    for order in [1, 2, 3] {
        let original = &coefficient * TwoQuarkLoops::weight(order, 0, 1, 1);
        let projected = projection.apply(&original)?;
        let mut explicit = -&original;
        for first in [-1, 1] {
            for second in [-1, 1] {
                explicit +=
                    &coefficient * TwoQuarkLoops::weight(order, 0, first, second) / Atom::num(4);
            }
        }
        let (bulk, sectors) = FermiSurfaceSector::extract(graph, &projected)?;
        let (explicit_bulk, explicit_sectors) = FermiSurfaceSector::extract(graph, &explicit)?;
        assert!(bulk.is_zero());
        assert_eq!(bulk, explicit_bulk);
        assert_eq!(sectors.len(), 2);
        assert_eq!(explicit_sectors.len(), 2);
        for sector in &sectors {
            assert_eq!(
                sector.product.factors(),
                [ThermalDistributionFactor {
                    edge_id: EdgeIndex(0),
                    sign: 0,
                    derivative_order: order as usize,
                }]
            );
            assert_eq!(sector.product.derivative_orders, [order as usize - 1]);
            let [orientation] = sector.orientations.as_slice() else {
                panic!("each term has exactly one Fermi support");
            };
            assert!(matches!(*orientation, -1 | 1));
            let explicit_sector = explicit_sectors
                .iter()
                .find(|candidate| candidate.orientations == sector.orientations)
                .unwrap();
            assert!(
                (&sector.coefficient - &explicit_sector.coefficient)
                    .replace_multiple(&*OCCUPATION_FORM)
                    .expand()
                    .is_zero()
            );
            let mut evaluator = sector
                .coefficient
                .as_view()
                .to_evaluation_tree(&graph.param_builder.fn_map, &params)?
                .linearize(&OptimizationSettings::default())
                .map_coeff(&|value| value.re.to_f64());
            for mu in [-5.0_f64, 5.0] {
                for energy in [2.0_f64, 7.0] {
                    let input = [mu.abs(), energy, energy, mu];
                    let step = |sign: i32| {
                        let nu = f64::from(sign) * mu;
                        (if energy > nu { 1.0 } else { 0.0 })
                            - 2.0 * if energy < -nu { 1.0 } else { 0.0 }
                    };
                    let expected = mu
                        * mu
                        * ((step(1) + step(-1)) / 4.0
                            - if *orientation == 1 { step(1) } else { 0.0 });
                    let actual = evaluator.evaluate_single(&input);
                    assert!(
                        (actual - expected).abs() < 1.0e-12,
                        "order={order}, orientation={orientation}, mu={mu}, E={energy}: {actual} != {expected}"
                    );
                }
            }
        }
        assert!(sectors.iter().any(|sector| sector.orientations == [-1]));
        assert!(sectors.iter().any(|sector| sector.orientations == [1]));
    }
    Ok(())
}

#[test]
fn phase_quenched_finite_temperature_derivative_weights_match_independent_signs() -> Result<()> {
    let fixture = TwoQuarkLoops::new(MediumMode::ThermodynamicEquilibrium)?;
    let graph = &fixture.graph;
    let projection = PhaseQuenchedDifference::new(graph)?;
    let params = [
        GS.ose(EdgeIndex(0)),
        GS.ose(EdgeIndex(2)),
        GS.ose(EdgeIndex(3)),
        GS.inverse_temperature.to_atom(),
        graph[EdgeIndex(0)].chemical_potential_atom().unwrap(),
    ];
    for order in [1, 2] {
        let projected = projection.apply(&TwoQuarkLoops::weight(order, 1, 1, 1))?;
        let mut evaluator = projected
            .as_view()
            .to_evaluation_tree(&graph.param_builder.fn_map, &params)?
            .linearize(&OptimizationSettings::default())
            .map_coeff(&|value| value.re.to_f64());
        for beta in [0.3_f64, 2.0] {
            for mu in [-3.0_f64, 0.0, 3.0] {
                let weight = |first: i32, second: i32| {
                    let x = beta * (1.3 - f64::from(first) * mu) / 2.0;
                    let first_derivative = beta / (4.0 * x.cosh().powi(2));
                    let derivative = if order == 1 {
                        first_derivative
                    } else {
                        -beta * x.tanh() * first_derivative
                    };
                    let occupation = |energy: f64, sign: i32| {
                        1.0 / (1.0 + (beta * (energy - f64::from(sign) * mu)).exp())
                    };
                    derivative * (1.0 - occupation(2.1, second) - 2.0 * occupation(2.1, -second))
                };
                let expected = (weight(-1, -1) + weight(-1, 1) + weight(1, -1) + weight(1, 1))
                    / 4.0
                    - weight(1, 1);
                let actual = evaluator.evaluate_single(&[1.3, 2.1, 2.1, beta, mu]);
                assert!(
                    (actual - expected).abs() < 1.0e-13,
                    "order={order}, beta={beta}, mu={mu}: {actual} != {expected}"
                );
            }
        }
    }
    Ok(())
}

#[test]
fn phase_quenched_setting_reaches_the_production_amplitude_evaluator() -> Result<()> {
    initialisation::test_initialise()?;
    let model = utils::load_generic_model("sm");
    let mut global: GlobalSettings = toml::from_str(
        r#"
[generation]
explicit_orientation_sum_only = true
[generation.medium]
mode = "thermodynamic_equilibrium"
phase_quenched_difference = true
[generation.uv]
subtract_uv = false
generate_integrated = false
softct = false
[generation.threshold_subtraction]
enable_thresholds = false
[generation.tropical_subgraph_table]
disable_tropical_generation = true
[generation.evaluator]
compile = false
summed = false
summed_function_map = false
iterative_orientation_optimization = false
"#,
    )?;
    let settings: RuntimeSettings = toml::from_str(
        r#"
[general]
evaluator_method = "SingleParametric"
inverse_temperature = 2.0
[kinematics]
e_cm = 10.0
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = []
helicities = []
[sampling]
graphs = "summed"
orientations = "summed"
sampling_multichanneling = false
"#,
    )?;
    let pool = ThreadPoolBuilder::new()
        .num_threads(1)
        .stack_size(256 * 1024 * 1024)
        .build()?;
    let input = MomentumSpaceEvaluationInput {
        loop_momenta: vec![ThreeMomentum::new(F(2.0), F(-1.0), F(0.5))],
        integrator_weight: F(1.0),
        graph_id: Some(0),
        group_id: None,
        orientation: None,
        channel_id: None,
    };
    for enabled in [false, true] {
        global.generation.medium.phase_quenched_difference = enabled;
        // A local scalar numerator and a positive mass keep this point regular.
        // With no quark loops the complete phase-quenched difference is zero,
        // although the ordinary finite-temperature contour is nonzero.
        let graphs = Graph::from_string(
            r#"digraph phase_quenched_gluon_tadpole {
                node [num=1];
                A -> A [id=0 lmb_id=0 particle="g" mass="1" num=1];
            }"#,
            &model,
        )?;
        let mut amplitude = Amplitude::from_graph_list("phase_quenched_gluon_tadpole", graphs)?;
        amplitude.preprocess(&model, &global.generation, &(&settings).into(), &pool)?;
        amplitude.build_integrand(
            &model,
            "phase_quenched_gluon_tadpole",
            &global,
            (&settings).into(),
            &pool,
        )?;
        let runtime = amplitude.integrand.as_mut().unwrap();
        runtime.warm_up(&model)?;
        let value = runtime
            .evaluate_momentum_configuration(&model, &input, false)?
            .integrand_result;
        assert!(value.re.0.is_finite() && value.im.0.is_finite());
        if enabled {
            assert_eq!(value.re, F(0.0));
            assert_eq!(value.im, F(0.0));
        } else {
            assert!(value.re.0.abs() + value.im.0.abs() > 1.0e-8);
        }
    }
    Ok(())
}
