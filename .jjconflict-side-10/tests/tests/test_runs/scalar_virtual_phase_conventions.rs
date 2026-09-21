use super::utils::*;
use super::*;
use gammaloop_api::commands::integrate::RendererOption;
use gammalooprs::{
    graph::FeynmanGraph,
    processes::{CutId, ProcessCollection},
};
use std::f64::consts::PI;

#[test]
#[serial]
fn scalar_virtual_triangle_mirrors_match_signed_born_interference() -> Result<()> {
    let parent_mass: f64 = 1.0;
    // Feynman parameters give the UV-finite triangle integral
    // -i/(16π²) integral_{x,y>=0,x+y<=1} dx dy/(m²-M²xy).
    // The simplex integral is 2 asin²(M/(2m))/M². Three vertices -ig and
    // three propagator numerators i multiply to g³, so A_triangle=-ig³c.
    // Sewing with A_Born* = +ig gives +g⁴c below threshold. Above threshold,
    // asin(x+i0)=π/2+i acosh(x): c has a positive imaginary part on the left,
    // and its complex conjugate enters when the virtual loop is on the right.
    for virtual_mass in [2.0_f64, 0.3] {
        let above_threshold = parent_mass > 2.0 * virtual_mass;
        let x = parent_mass / (2.0 * virtual_mass);
        let (angle_re, angle_im) = if above_threshold {
            (PI / 2.0, x.acosh())
        } else {
            (x.asin(), 0.0)
        };
        let normalization = 2.0
            / (16.0 * PI * PI * parent_mass * parent_mass)
            / (8.0 * PI)
            / (2.0 * parent_mass)
            / 2.0;
        let expected_real = normalization * (angle_re * angle_re - angle_im * angle_im);
        let expected_imaginary = normalization * 2.0 * angle_re * angle_im;
        assert!(expected_real > 0.0);
        let mut mirror_results = Vec::new();
        for triangle_on_right in [false, true] {
            let side = if triangle_on_right { "right" } else { "left" };
            let region = if above_threshold { "above" } else { "below" };
            let name = format!("scalar_virtual_triangle_{region}_{side}");
            let expected = (
                expected_real,
                if triangle_on_right {
                    -expected_imaginary
                } else {
                    expected_imaginary
                },
            );
            let root = get_tests_workspace_path().join(&name);
            clean_test(&root);
            let mut cli = get_test_cli(None, root, Some(name), true)?;
            let (incoming, outgoing, cut_edges) = if triangle_on_right {
                ("d", "a", "d -> b [id=3 lmb_id=0]; d -> c [id=4];")
            } else {
                ("a", "d", "b -> d [id=3 lmb_id=0]; c -> d [id=4];")
            };
            // The loop is on a,b,c; d is the Born vertex. Final scalar_0
            // lines are the only process-selected physical cut. For m=2 the
            // heavy threshold is closed; for m=0.3 it requires a causal CT on
            // the virtual side. The labelled graph has
            // exactly one identical-final-state factor 1/2! and no extra sum
            // over the mirror: each graph is tested against its own target.
            let dot = format!(
                r#"digraph virtual_triangle {{
                overall_factor="1/2";
                edge [particle=scalar_0];
                ext0 [style=invis is_cut=0];
                a [int_id="V_3_SCALAR_112"];
                b [int_id="V_3_SCALAR_011"];
                c [int_id="V_3_SCALAR_011"];
                d [int_id="V_3_SCALAR_002"];
                ext0 -> {incoming} [particle=scalar_2];
                {outgoing} -> ext0 [particle=scalar_2];
                a -> b [id=0 particle=scalar_1 lmb_id=1];
                a -> c [id=1 particle=scalar_1];
                b -> c [id=2 particle=scalar_1];
                {cut_edges}
            }}"#,
            );
            run_commands(
                &mut cli,
                &[
                    "import model scalars-default.json",
                    &format!("set model mass_scalar_1={virtual_mass} mass_scalar_2=1.0 lam=1.0"),
                    &format!(
                        "set global kv global.generation.evaluator.compile=false global.generation.evaluator.summed=false global.generation.evaluator.iterative_orientation_optimization=false global.generation.uv.subtract_uv=false global.generation.uv.generate_integrated=false global.generation.threshold_subtraction.enable_thresholds={above_threshold} global.generation.threshold_subtraction.disable_integrated_ct=false"
                    ),
                    &format!(
                        r#"set default-runtime string '
[general]
evaluator_method = "SingleParametric"
integral_unit = "none"
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
[h_function]
function = "poly_exponential"
sigma = 1.0
[kinematics]
e_cm = 1.0
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = [[1.0,0.0,0.0,0.0]]
helicities = [0]
[integrator]
integrated_phase = "real"
n_start = 20000
n_increase = 20000
n_max = 500000
min_samples_for_update = 1000
seed = 19427
target_relative_accuracy = 0.02
'"#,
                        !above_threshold
                    ),
                    &format!(
                        "import graphs --inline-dot \"\"\"{dot}\"\"\" --process-spec 'scalar_2 > scalar_0 scalar_0' -p interference -i virtual -o"
                    ),
                    "generate existing -p interference -i virtual",
                ],
            )?;
            let ProcessCollection::CrossSections(cross_sections) =
                &cli.state.process_list.processes[0].collection
            else {
                panic!("expected a scalar cross section")
            };
            let [supergraph] = cross_sections["virtual"].supergraphs.as_slice() else {
                panic!("the interference oracle must select exactly one graph")
            };
            assert_eq!(supergraph.graph.get_loop_number(), 2);
            assert_eq!(supergraph.cuts.len(), 1, "expected only the two-scalar cut");
            assert_eq!(
                supergraph
                    .graph
                    .iter_edges_of(&supergraph.cuts[CutId(0)].cut)
                    .count(),
                2,
            );
            for point in [
                vec![0.17, -0.09, 0.13, 0.41, 0.23, -0.36],
                vec![0.24, 0.08, -0.11, -0.27, 0.35, 0.19],
            ] {
                let (_, value) = Inspect {
                    process: Some(ProcessRef::Unqualified("interference".to_string())),
                    integrand_name: Some("virtual".to_string()),
                    point,
                    momentum_space: true,
                    graph_id: Some(0),
                    ..Default::default()
                }
                .run(&mut cli.state)?;
                assert!(
                    value.re.is_finite() && value.im.is_finite(),
                    "{region}/{side}: {value:e}"
                );
                if above_threshold {
                    // The remaining propagator on the massive two-particle cut
                    // is -2ℓ·p_final<0. Its absorptive weight has a uniform sign.
                    assert!(
                        value.im * expected.1 > 0.0,
                        "{side}: incorrect causal phase {value:e}"
                    );
                } else {
                    assert!(
                        value.re > 0.0 && value.im.abs() <= 1.0e-11 * value.re,
                        "{side}: below-threshold scalar interference is not positive real: {value:e}"
                    );
                }
            }
            let output = Integrate {
                process: vec![ProcessRef::Unqualified("interference".to_string())],
                integrand_name: vec!["virtual".to_string()],
                workspace_path: Some(cli.cli_settings.state.folder.join("integration_workspace")),
                n_cores: Some(1),
                restart: true,
                renderer: RendererOption::Tabled,
                show_max_weight_info: false,
                no_stream_iterations: true,
                no_stream_updates: true,
                ..Default::default()
            }
            .run(&mut cli.state, &cli.cli_settings)?;
            let integral = single_slot_integral(&output);
            for (value, error, target) in [
                (integral.result.re.0, integral.error.re.0, expected.0),
                (integral.result.im.0, integral.error.im.0, expected.1),
            ] {
                assert!(value.is_finite() && error.is_finite());
                if target == 0.0 {
                    assert!(
                        value.abs() <= 1.0e-12 * expected_real
                            && error.abs() <= 1.0e-12 * expected_real,
                        "{side}: below-threshold imaginary component: {integral:?}"
                    );
                } else {
                    assert!(
                        value * target > 0.0 && error < 0.05 * target.abs(),
                        "{region}/{side}: {integral:?}"
                    );
                    assert!(
                        (value - target).abs() <= (5.0 * error).max(1.0e-10 * target.abs()),
                        "{region}/{side}: {value:e} ± {error:e}, signed oracle {target:e}"
                    );
                }
            }
            mirror_results.push((integral.result, integral.error));
            clean_test(&cli.cli_settings.state.folder);
        }
        let [(left, left_error), (right, right_error)] = mirror_results.as_slice() else {
            unreachable!()
        };
        assert!(
            (left.re.0 - right.re.0).abs() <= 5.0 * left_error.re.0.hypot(right_error.re.0),
            "left/right real interference differs: {mirror_results:?}"
        );
        assert!(
            (left.im.0 + right.im.0).abs()
                <= (5.0 * left_error.im.0.hypot(right_error.im.0)).max(1.0e-12 * expected_real),
            "left/right imaginary interference is not conjugated: {mirror_results:?}"
        );
    }
    Ok(())
}
