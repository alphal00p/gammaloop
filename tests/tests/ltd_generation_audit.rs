//! Process-isolated assertions about which native generator actually executes.
//!
//! The counters include all worker threads. Each parent launches exactly one
//! ignored child test, so unrelated tests can neither hide nor add CFF calls.

use std::{collections::BTreeMap, env, fs, process::Command};

use color_eyre::Result;
use gammaloop_api::state::SyncSettings;
use gammaloop_integration_tests::{
    clean_test, get_test_cli, get_tests_workspace_path, run_commands, workspace_root,
};
use gammalooprs::{processes::ProcessCollection, uv::UVOrchestrator};
use three_dimensional_reps::{
    CffGenerationContext, Generate3DExpressionOptions, MomentumSignature, ParsedGraph,
    RepresentationMode, generate_3d_expression,
    generation::audit::GenerationInvocations,
    graph_io::{ParsedGraphExternalEdge, ParsedGraphInternalEdge},
};

const CHILD_CASE: &str = "GAMMALOOP_LTD_GENERATION_AUDIT_CASE";

fn run_isolated(case: &str) -> Result<()> {
    let output = Command::new(env::current_exe()?)
        .current_dir(workspace_root())
        .args([
            "--exact",
            "generation_audit_child",
            "--ignored",
            "--nocapture",
        ])
        .env(CHILD_CASE, case)
        .output()?;
    eyre::ensure!(
        output.status.success(),
        "generation audit {case} failed: {}\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    Ok(())
}

#[test]
fn cff_invocation_audit_detects_a_positive_control() -> Result<()> {
    run_isolated("cff_control")
}

#[test]
fn ltd_disconnected_raised_sources_never_build_cff() -> Result<()> {
    run_isolated("disconnected")
}

#[test]
fn ltd_uv_and_threshold_generation_never_build_cff() -> Result<()> {
    for case in ["legacy", "hedge"] {
        run_isolated(case)?;
    }
    Ok(())
}

#[test]
fn ltd_cut_and_threshold_generation_never_build_cff() -> Result<()> {
    for case in ["cuts_legacy", "cuts_hedge"] {
        run_isolated(case)?;
    }
    Ok(())
}

#[test]
#[ignore = "launched by the process-isolated generation audits"]
fn generation_audit_child() -> Result<()> {
    let case = env::var(CHILD_CASE)?;
    assert_eq!(GenerationInvocations::snapshot().cff_builders, 0);
    if case == "cff_control" || case == "disconnected" {
        return audit_disconnected_source(case == "cff_control");
    }
    let root = get_tests_workspace_path().join(format!("ltd_generation_audit_{case}"));
    let mut cli = get_test_cli(None, &root, Some(format!("ltd_audit_{case}")), true)?;
    let cuts = case.starts_with("cuts_");
    let generation = &mut cli.cli_settings.global.generation;
    generation.three_dimensional_representations = vec![RepresentationMode::Ltd];
    generation.uv.orchestrator = if case.ends_with("legacy") {
        UVOrchestrator::LegacyDagForest
    } else {
        UVOrchestrator::HedgePoset
    };
    generation.uv.local_uv_cts_from_expanded_4d_integrands = true;
    generation.uv.subtract_uv = true;
    generation.uv.generate_integrated = true;
    generation.uv.softct = false;
    generation.evaluator.compile = false;
    generation.evaluator.summed = false;
    generation.evaluator.summed_function_map = true;
    generation.evaluator.iterative_orientation_optimization = false;
    generation.threshold_subtraction.enable_thresholds = true;
    generation.threshold_subtraction.disable_integrated_ct = false;
    generation
        .threshold_subtraction
        .check_esurface_at_generation = false;
    generation
        .tropical_subgraph_table
        .disable_tropical_generation = true;
    if cuts {
        generation.force_cuts = vec![
            vec!["e5".into(), "e7".into()],
            vec!["e6".into(), "e7".into()],
        ];
    }
    cli.cli_settings.global.n_cores.generate = 2;
    cli.cli_settings.sync_settings()?;
    run_commands(
        &mut cli,
        &[
            "import model assets/models/json/scalars/scalars_2p_3p.json",
            "set model mass_scalar_1=2.0",
            "set model mass_scalar_2=2.0001",
            if cuts {
                r#"set default-runtime string '
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = [[5.0, 0.0, 0.0, 0.0]]
helicities = [0]
'"#
            } else {
                r#"set default-runtime string '
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = [[5.0, 0.0, 0.0, 0.0], "dependent"]
helicities = [0, 0]
'"#
            },
            if cuts {
                "import graphs tests/resources/graphs/double_dotted_double_triangle.dot -p audit -i audit"
            } else {
                "import graphs tests/resources/graphs/scalar_self_energy_amp.dot -p audit -i audit"
            },
            "generate existing -p audit -i audit",
        ],
    )?;

    let process = &cli.state.process_list.processes[0];
    match &process.collection {
        ProcessCollection::Amplitudes(amplitudes) => {
            let graphs = &amplitudes["audit"].graphs;
            assert!(!cuts && !graphs.is_empty());
            assert!(
                graphs.iter().any(|graph| {
                    let source = &graph.derived_data.representations[&RepresentationMode::Ltd];
                    source.threshold_counterterms.iter().any(|counterterm| {
                        counterterm
                            .parametric
                            .iter()
                            .map(|(_, expression)| expression)
                            .any(|atom| !atom.is_zero())
                    })
                }),
                "the amplitude must actually generate threshold counterterms"
            );
            for graph in graphs {
                assert_eq!(graph.derived_data.representations.len(), 1);
                let source = &graph.derived_data.representations[&RepresentationMode::Ltd];
                assert!(
                    source
                        .all_mighty_integrand
                        .iter()
                        .map(|(_, expression)| expression)
                        .any(|atom| !atom.is_zero())
                );
            }
        }
        ProcessCollection::CrossSections(cross_sections) => {
            let graphs = &cross_sections["audit"].supergraphs;
            assert!(cuts && !graphs.is_empty());
            assert!(
                graphs.iter().any(|graph| {
                    graph
                        .derived_data
                        .threshold_counterterms
                        .iter()
                        .any(|counterterm| {
                            let sources = &counterterm.integrands[&RepresentationMode::Ltd];
                            sources
                                .left_atoms
                                .iter()
                                .chain(&sources.right_atoms)
                                .any(|source| {
                                    source
                                        .integrands
                                        .iter()
                                        .map(|(_, expression)| expression)
                                        .any(|atom| !atom.is_zero())
                                })
                        })
                }),
                "the cut fixture must actually generate nonzero threshold programs"
            );
            for graph in graphs {
                assert_eq!(graph.derived_data.expressions.len(), 1);
                assert!(
                    !graph.derived_data.cut_parametric_integrands[&RepresentationMode::Ltd]
                        .is_empty()
                );
            }
        }
    }

    let generated = GenerationInvocations::snapshot();
    assert_eq!(generated.cff_builders, 0, "{generated:#?}");
    assert!(generated.ltd.iter().any(|call| call.loops == 2));
    assert!(generated.ltd.iter().any(|call| call.repeated_poles));
    assert!(generated.ltd_residue_selections > 0);
    assert!(
        generated
            .ltd
            .iter()
            .any(|call| call.thread != std::thread::current().id()),
        "generation must exercise worker-thread recording"
    );
    if !cuts {
        assert!(
            generated.ltd.iter().any(|call| {
                call.context == CffGenerationContext::EmbeddedCffFactor && call.loops > 0
            }),
            "the UV child projection must execute: {generated:#?}"
        );
        assert!(
            generated.ltd.iter().any(|call| {
                call.context == CffGenerationContext::Standalone && call.loops == 1
            }),
            "a nontrivial outer cograph must execute: {generated:#?}"
        );
    }

    let diagnostic = root.join("ltd-residue-map.json");
    cli.run_command(&format!(
        "3drep build -p audit -i audit -g 0 --representation ltd --workspace-path {} --json-out {} --no-pretty",
        root.join("diagnostic").display(), diagnostic.display(),
    ))?;
    let value: serde_json::Value = serde_json::from_str(&fs::read_to_string(diagnostic)?)?;
    assert_eq!(value["family"], "ltd");
    let after_diagnostic = GenerationInvocations::snapshot();
    assert_eq!(after_diagnostic.cff_builders, 0, "{after_diagnostic:#?}");
    assert!(after_diagnostic.ltd.len() > generated.ltd.len());
    clean_test(&root);
    Ok(())
}

fn audit_disconnected_source(cff_control: bool) -> Result<()> {
    // Two independent loop components; the first has a genuine double pole.
    // No GammaLoop CFF oracle is constructed in the LTD child process.
    let graph = ParsedGraph {
        internal_edges: (0..4)
            .map(|edge_id| ParsedGraphInternalEdge {
                edge_id,
                tail: 2 * (edge_id / 2) + edge_id % 2,
                head: 2 * (edge_id / 2) + (edge_id + 1) % 2,
                label: format!("e{edge_id}"),
                mass_key: Some(if edge_id < 2 {
                    "m".to_owned()
                } else {
                    format!("m{edge_id}")
                }),
                signature: MomentumSignature {
                    loop_signature: if edge_id < 2 { vec![1, 0] } else { vec![0, 1] },
                    external_signature: vec![i32::from(edge_id == 3)],
                },
                had_pow: false,
            })
            .collect(),
        external_edges: vec![
            ParsedGraphExternalEdge {
                edge_id: 0,
                source: Some(2),
                destination: None,
                label: "outgoing".into(),
                external_coefficients: vec![1],
            },
            ParsedGraphExternalEdge {
                edge_id: 1,
                source: None,
                destination: Some(3),
                label: "incoming".into(),
                external_coefficients: vec![1],
            },
        ],
        initial_state_cut_edges: Vec::new(),
        loop_names: vec!["k0".into(), "k1".into()],
        external_names: vec!["p0".into()],
        node_name_to_internal: (0..4)
            .map(|node| (format!("v{node}"), node))
            .collect::<BTreeMap<_, _>>(),
    };
    let expression = generate_3d_expression(
        &graph,
        &Generate3DExpressionOptions {
            representation: if cff_control {
                RepresentationMode::Cff
            } else {
                RepresentationMode::Ltd
            },
            energy_degree_bounds: Some(Vec::new()),
            ..Default::default()
        },
    )?;
    assert_eq!(expression.energy_factor_components.len(), 2);
    assert!(!expression.expression.orientations.is_empty());
    let calls = GenerationInvocations::snapshot();
    if cff_control {
        assert!(
            calls.cff_builders > 0,
            "the audit must detect real CFF generation"
        );
        assert!(calls.ltd.is_empty());
    } else {
        assert_eq!(calls.cff_builders, 0, "{calls:#?}");
        assert_eq!(calls.ltd.len(), 2);
        assert!(calls.ltd.iter().all(|call| call.loops == 1));
        assert!(calls.ltd.iter().any(|call| call.repeated_poles));
    }
    Ok(())
}
