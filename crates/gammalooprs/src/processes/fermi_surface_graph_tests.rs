use std::f64::consts::PI;

use eyre::Result;
use linnet::half_edge::involution::EdgeIndex;
use spenso::algebra::complex::Complex;
use typed_index_collections::TiVec;

use crate::{
    graph::{FeynmanGraph, Graph, GraphGroupPosition},
    initialisation::test_initialise,
    integrands::process::{
        GraphTerm,
        amplitude::AmplitudeGraphTerm,
        evaluators::{GenericEvaluatorFloat, InputParams},
        param_builder::ParamBuilderGraph,
    },
    processes::AmplitudeGraph,
    settings::{
        GlobalSettings, RuntimeSettings,
        global::{GenerationSettings, MediumMode},
    },
    utils::{F, load_generic_model},
    uv::uv_graph::UVE,
};

#[test]
fn fermi_surface_generated_dotted_vacuum_cycles_match_mass_derivatives() -> Result<()> {
    test_initialise()?;
    let mut model = load_generic_model("sm");
    let mu = 2.0_f64;
    let mass = 1.0_f64;
    let fermi_radius = (mu * mu - mass * mass).sqrt();
    let logarithm = ((mu + fermi_radius) / mass).ln();
    // After removing the production i/(2*pi)^3 spatial prefactor, the source
    // contour's occupied tadpole is +theta(mu-E)/(2E). Raising D=q0²-q²-m²
    // differentiates with respect to m², including the moving boundary.
    let references = [
        PI * (mu * fermi_radius - mass * mass * logarithm),
        -PI * logarithm,
        PI * mu / (4.0 * mass * mass * fermi_radius),
    ];

    for (power, reference) in (1..=3).zip(references) {
        for explicit_orientation_sum_only in [false, true] {
            let edges = (0..power)
                .map(|edge| {
                    format!(
                        "v{edge} -> v{} [id={edge}, particle=\"d\", mass=1, num=1];",
                        (edge + 1) % power
                    )
                })
                .collect::<Vec<_>>()
                .join("\n");
            let source = format!(
                "digraph dotted_vacuum_{power} {{
                num=1; projector=1; overall_factor=1;
                node [num=1];
                {edges}
            }}"
            );
            let mut graph = AmplitudeGraph::new(Graph::from_string(&source, &model)?.remove(0));
            assert_eq!(graph.graph.get_loop_number(), 1);
            let mu_parameter = graph.graph[EdgeIndex(0)]
                .particle()
                .unwrap()
                .chemical_potential
                .unwrap();
            model.parameters.get_mut(&mu_parameter).unwrap().value = Some(Complex::new_re(F(mu)));
            graph.graph.param_builder.update_model_values(&model);
            let mut settings = GenerationSettings::default();
            settings.medium.mode = MediumMode::ZeroTemperatureEquilibrium;
            settings.medium.vacuum_subtraction = true;
            settings.explicit_orientation_sum_only = explicit_orientation_sum_only;
            settings.threshold_subtraction.enable_thresholds = false;
            settings.uv.subtract_uv = false;
            settings.uv.generate_integrated = false;
            graph.preprocess(&model, &settings, &(&RuntimeSettings::default()).into())?;
            if power > 1 {
                assert!(
                    graph
                        .derived_data
                        .cff_expression
                        .as_ref()
                        .unwrap()
                        .expression
                        .orientations
                        .iter()
                        .flat_map(|orientation| &orientation.variants)
                        .flat_map(|variant| &variant.thermal_weight.distributions)
                        .any(|factor| factor.derivative_order == power - 1),
                    "the native dotted cycle must generate its Fermi-surface derivative",
                );
            }
            graph
                .graph
                .param_builder
                .numerator_sampling_scale_value(Complex::new_re(F(1.0)));
            let (mut term, _) = AmplitudeGraphTerm::from_amplitude_graph(
                &graph,
                GraphGroupPosition(0),
                TiVec::new(),
                &model,
                &GlobalSettings {
                    generation: settings,
                    ..Default::default()
                },
            )?;
            let builder = &graph.graph.param_builder;
            let params = (&builder.pairs)
                .into_iter()
                .flat_map(|pair| &pair.params)
                .collect::<Vec<_>>();
            let chemical_potential = graph.graph[EdgeIndex(0)].chemical_potential_atom().unwrap();
            let mu_slot = params
                .iter()
                .position(|param| **param == chemical_potential)
                .unwrap();
            let loop_slots = graph
                .graph
                .loop_mom_params(&graph.graph.loop_momentum_basis)
                .iter()
                .map(|momentum| params.iter().position(|param| *param == momentum).unwrap())
                .collect::<Vec<_>>();
            let mut values = builder.values[0].clone();
            values[mu_slot] = Complex::new_re(F(mu));
            let orientations = term
                .original_integrand
                .production_orientation_ids()
                .iter()
                .copied()
                .zip(term.orientations.iter().cloned())
                .collect::<Vec<_>>();
            {
                let mut evaluate = <f64 as GenericEvaluatorFloat>::get_evaluator_single(
                    &mut term.original_integrand.single_parametric,
                );
                let mut evaluate_complete = |values: &mut [Complex<F<f64>>]| {
                    if explicit_orientation_sum_only {
                        return evaluate(values);
                    }
                    let mut sum = Complex::new_re(F(0.0));
                    for (residue_map, orientation) in &orientations {
                        values[builder.pairs.residue_map_id.value_range.start] =
                            Complex::new_re(F(residue_map.0 as f64));
                        InputParams::<f64>::set_orientation_values_impl(
                            values,
                            Complex::new_re(F(1.0)),
                            Complex::new_re(F(0.0)),
                            1,
                            builder.pairs.orientations.value_range.start,
                            orientation,
                        );
                        sum += evaluate(values);
                    }
                    sum
                };

                // The x=.5 boundary is a quadrature cell boundary at both resolutions.
                // The auxiliary shell density becomes constant under this radial map;
                // the remaining compactly supported bulk is smooth on each half.
                let mut previous_error = f64::INFINITY;
                for steps in [512, 1024] {
                    let mut integral = 0.0;
                    for index in 0..steps {
                        let x = (index as f64 + 0.5) / steps as f64;
                        let radius = fermi_radius * x / (1.0 - x);
                        let jacobian =
                            4.0 * PI * radius * radius * fermi_radius / (1.0 - x).powi(2);
                        values[loop_slots[0]] = Complex::new_re(F(radius));
                        values[loop_slots[1]] = Complex::new_re(F(0.0));
                        values[loop_slots[2]] = Complex::new_re(F(0.0));
                        let result = evaluate_complete(&mut values);
                        assert!(result.re.0.abs() < 1.0e-12);
                        integral += result.im.0 * (2.0 * PI).powi(3) * jacobian / steps as f64;
                    }
                    let error = (integral - reference).abs();
                    assert!(
                        error < previous_error * 0.3,
                        "D^(-{power}), explicit={explicit_orientation_sum_only}, {steps} cells: {integral} != {reference}",
                    );
                    previous_error = error;
                }
                assert!(previous_error < 2.0e-5 * reference.abs());

                // Vacuum subtraction removes the whole graph below onset, including
                // every localized derivative contribution from the dotted cases.
                values[mu_slot] = Complex::new_re(F(0.5));
                values[loop_slots[0]] = Complex::new_re(F(0.7));
                let absent = evaluate_complete(&mut values);
                assert!(absent.re.0.abs() < 1.0e-13);
                assert!(absent.im.0.abs() < 1.0e-13);
            }

            if power > 1 {
                model.parameters.get_mut(&mu_parameter).unwrap().value =
                    Some(Complex::new_re(F(mass)));
                let error = term
                    .warm_up(&RuntimeSettings::default(), &model)
                    .unwrap_err();
                assert!(error.to_string().contains("nondegenerate shell"));
            }
        }
    }
    Ok(())
}
