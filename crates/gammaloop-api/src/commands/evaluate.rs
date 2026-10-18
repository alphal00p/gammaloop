use std::{collections::BTreeMap, fs, path::PathBuf};

use clap::Args;
use gammalooprs::{
    graph::FeynmanGraph,
    processes::{AmplitudeGraph, AnalyticalEvaluationConfig},
    uv::UltravioletGraph,
};

use gammalooprs::processes::{Amplitude, ProcessCollection};

use gammalooprs::utils::vakint;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use color_eyre::Result;
use colored::Colorize;
use gammalooprs::settings::RuntimeSettings;
use idenso::tensor::SymbolicTensor;
use spenso::structure::IndexedStructure;
use symbolica::atom::{Atom, AtomCore};
use tracing::{info, warn};

use crate::{
    completion::CompletionArgExt,
    state::{ProcessRef, State},
    CLISettings,
};

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, unsendable, name = "EvaluationSettings")
)]
#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Evaluate {
    /// Process reference: `#<id>`, `name:<name>`, or `<id>/<name>`
    #[arg(
        short = 'p',
        long = "process",
        value_name = "PROCESS",
        completion_process_selector(crate::completion::SelectorKind::Amplitude)
    )]
    pub process: Option<ProcessRef>,

    /// The integrand name to evaluate
    #[arg(
        short = 'i',
        long = "integrand-name",
        value_name = "NAME",
        completion_integrand_selector(crate::completion::SelectorKind::Amplitude)
    )]
    pub graphs_group_name: Option<String>,

    /// Write the canonical symbolic evaluation result to this TOML file
    #[arg(short = 'o', long, value_hint = clap::ValueHint::FilePath)]
    pub result_path: Option<PathBuf>,

    /// Whether to evaluate numerically or not the resulting analytical expression
    #[arg(short = 'm', long = "numerical")]
    pub numerical: bool,

    /// The number of terms in the epsilon expansion to compute
    /// Defaults is automatically inferred from input graph to get to the order O(epsilon^0)
    #[arg(short = 'e', long = "n-epsilon-terms")]
    pub number_of_terms_in_epsilon_expansion: Option<usize>,
}

impl Evaluate {
    pub fn run(
        &self,
        state: &mut State,
        global_cli_settings: &CLISettings,
        default_runtime_settings: &RuntimeSettings,
    ) -> Result<Atom> {
        if let Some(path) = &self.result_path {
            global_cli_settings
                .ensure_write_target_outside_active_state(path, "write evaluation result")?;
        }

        let (process_id, integrand_name) =
            state.find_integrand_ref(self.process.as_ref(), self.graphs_group_name.as_ref())?;

        let amplitude: &Amplitude = match &state.process_list.processes[process_id].collection {
            ProcessCollection::Amplitudes(amplitudes) => amplitudes.get(&integrand_name).unwrap(),
            ProcessCollection::CrossSections(_) => {
                return Err(color_eyre::eyre::eyre!(
                    "Evaluate command does not support cross-section graphs"
                ));
            }
        };
        // Always resolve the model through the targeted integrand. `state.model` only carries the
        // shared structural/default model and would ignore any per-integrand parameter overlay.
        let model = state.resolve_model_for_integrand(process_id, &integrand_name)?;

        let mut true_settings = global_cli_settings
            .global
            .generation
            .uv
            .vakint
            .true_settings();
        let vakint = vakint()?;

        if let Some(n_terms) = self.number_of_terms_in_epsilon_expansion {
            true_settings.number_of_terms_in_epsilon_expansion = n_terms as i64;
        }

        let mut full_evaluation: Option<SymbolicTensor<IndexedStructure>> = None;

        for graph_term in amplitude.graphs.iter() {
            let g = &graph_term.graph;
            if self.number_of_terms_in_epsilon_expansion.is_none() {
                true_settings.number_of_terms_in_epsilon_expansion =
                    g.n_loops(&g.no_dummy()) as i64 + 1;
            }
            let mut complete_evaluation_for_this_graph =
                SymbolicTensor::infer_indexed(Atom::one())?;
            if g.n_externals() != 0 {
                return Err(color_eyre::eyre::eyre!(
                    "Graph named: {} has external legs. Analytical evaluation in gammaloop only supports vacuum graphs.",
                    graph_term.graph.name
                ));
            }
            let connected_components = g.connected_components(&g.full_filter());
            if connected_components.len() > 1 && !g.global_prefactor.num.expression().is_one() {
                warn!(
                    "Graph named: {} has more than one ({}) connected components and its expression contains a global numerator ({}), which will only be applied to the first connected component. Make sure this is intended.",
                    graph_term.graph.name,
                    connected_components.len(),
                    feynkit_graph::NumeratorInterface::export(&g.global_prefactor.num).0
                );
            }
            for (i_gc, gc) in connected_components.iter().enumerate() {
                info!(
                    "Evaluating connected component #{} of graph named: {}",
                    i_gc + 1,
                    graph_term.graph.name.blue()
                );
                let component = graph_term.analytical_evaluation(
                    gc,
                    AnalyticalEvaluationConfig {
                        model: &model,
                        evaluate_numerically: self.numerical,
                        vakint,
                        true_settings: &true_settings,
                        settings: &global_cli_settings.global.generation.uv.vakint,
                        run_time_settings: default_runtime_settings,
                        // Only include the overall global numerator on the first connected component.
                        include_global_numerator: i_gc == 0,
                    },
                )?;
                complete_evaluation_for_this_graph =
                    complete_evaluation_for_this_graph.multiply(&component)?;
            }
            let overall_factor = SymbolicTensor::infer_indexed(g.overall_factor.clone())?;
            let complete = complete_evaluation_for_this_graph
                .multiply(&g.global_prefactor.projector)?
                .multiply(&overall_factor)?;
            full_evaluation = Some(match full_evaluation {
                Some(sum) => sum.try_add(&complete)?,
                None => complete,
            });
        }

        // The command returns and optionally writes a source expression. Tensor
        // interfaces remain attached until this explicit presentation boundary.
        let full_evaluation = full_evaluation
            .map(|tensor| feynkit_graph::NumeratorInterface::export(&tensor).0)
            .unwrap_or(Atom::Zero);

        if let Some(p) = self.result_path.clone() {
            info!(
                "Saving evaluation result for process {} to path: {}",
                state.process_list.processes[process_id]
                    .definition
                    .folder_name
                    .green(),
                p.display()
            );
            fs::write(
                p,
                toml::to_string_pretty(&BTreeMap::from([(
                    "expression",
                    full_evaluation.to_canonical_string(),
                )]))?,
            )?;
        } else if self.numerical {
            let numerical_evaluation_result =
                AmplitudeGraph::to_numerical(full_evaluation.as_view(), &true_settings)?;
            info!(
                "Numerical evaluation of the analytical result for process {}:\n{}",
                state.process_list.processes[process_id]
                    .definition
                    .folder_name
                    .green(),
                numerical_evaluation_result
            );
        } else {
            info!(
                "Analytical result for process {}:\n{}",
                state.process_list.processes[process_id]
                    .definition
                    .folder_name
                    .green(),
                full_evaluation
            );
        }

        Ok(full_evaluation)
    }
}

#[cfg(test)]
mod tests {
    use super::Evaluate;
    use crate::{
        commands::{
            generate::Generate,
            import::{model::ImportModel, Import},
        },
        state::State,
        CLISettings,
    };
    use gammalooprs::{
        initialisation::test_initialise,
        settings::{global::MediumMode, RuntimeSettings},
    };
    use std::path::PathBuf;

    #[test]
    fn evaluate_vacuum_tensor_restores_numeric_ports_before_global_projection(
    ) -> color_eyre::Result<()> {
        use feynkit_graph::NumeratorInterface;
        use gammalooprs::{
            processes::{AnalyticalEvaluationConfig, ProcessCollection},
            utils::{vakint, GS},
        };
        use idenso::tensor::ContractSettings;
        use spenso::structure::{
            representation::{Minkowski, RepName},
            slot::IsAbstractSlot,
        };
        use symbolica::{
            atom::{Atom, AtomCore},
            parse,
        };

        test_initialise()?;
        let temp = tempfile::tempdir()?;
        let mut state = State::new_test();
        let mut cli_settings = CLISettings::default();
        cli_settings.state.folder = temp.path().join("state");
        ImportModel {
            path: "scalars-default.json".into(),
            simplify_model: true,
        }
        .run(&mut state)?;
        // Both sources have the same scalar vacuum integral. The first keeps
        // its two explicit numeric D-dimensional ports until its projector.
        for (name, numerator, projector) in [
            (
                "tensor",
                "spenso::g(spenso::mink(gammalooprs::dim,0),spenso::mink(gammalooprs::dim,1))",
                "spenso::g(spenso::mink(gammalooprs::dim,0),spenso::mink(gammalooprs::dim,1))/gammalooprs::dim",
            ),
            ("scalar", "1", "1"),
        ] {
            Import::Graphs {
                source: None,
                inline_dot: Some(format!(
                    r#"digraph {name} {{
                        num="{numerator}"; projector="{projector}";
                        overall_factor=1; overall_factor_evaluated=1;
                        node [num=1];
                        a:0 -> a:1 [particle="scalar_1" num=1 id=0 lmb_id=0,
                            lmb_rep="K(0,a___)", source="{{ufo_order:0}}", sink="{{ufo_order:1}}"];
                    }}"#
                )),
                process_spec: None,
                process: Some(crate::state::ProcessRef::Name("vacuum_tensor".into())),
                integrand_name: Some(name.into()),
                overwrite: false,
                append: false,
            }
            .run(&mut state, &cli_settings)?;
        }
        let runtime = RuntimeSettings::default();
        let model = state.resolve_model_for_integrand(0, "tensor")?;
        let ProcessCollection::Amplitudes(amplitudes) = &state.process_list.processes[0].collection
        else {
            unreachable!("vacuum fixtures are amplitudes")
        };
        let graph = &amplitudes["tensor"].graphs[0];
        let settings = &cli_settings.global.generation.uv.vakint;
        let true_settings = settings.true_settings();
        let evaluate = |graph: &gammalooprs::processes::AmplitudeGraph| {
            graph.analytical_evaluation(
                &graph.graph.full_filter(),
                AnalyticalEvaluationConfig {
                    model: &model,
                    evaluate_numerically: false,
                    vakint: vakint()?,
                    true_settings: &true_settings,
                    settings,
                    run_time_settings: &runtime,
                    include_global_numerator: true,
                },
            )
        };
        let component = evaluate(graph)?;
        let source = &graph.graph.global_prefactor.num;
        assert_eq!(component.rank(), 2);
        assert_eq!(
            component.structure().logical_slots(),
            source.structure().logical_slots()
        );
        for slot in component.structure().logical_slots() {
            assert_eq!(slot.rep(), Minkowski {}.new_rep(GS.dim).to_lib());
        }
        assert!(!component.is_zero());
        let closed = component.multiply(&graph.graph.global_prefactor.projector)?;
        assert!(
            closed.is_scalar(),
            "the declared projector must close both ports"
        );
        let mut zero_graph = graph.clone();
        zero_graph.graph.global_prefactor.num = source.with_checked_expression(Atom::Zero)?;
        let zero = evaluate(&zero_graph)?;
        assert!(zero.is_zero());
        assert_eq!(
            zero.structure().logical_slots(),
            source.structure().logical_slots()
        );

        let mut results = Vec::new();
        for name in ["tensor", "scalar"] {
            let output = Evaluate {
                process: Some(crate::state::ProcessRef::Id(0)),
                graphs_group_name: Some(name.into()),
                result_path: Some(temp.path().join(format!("{name}.toml"))),
                numerical: false,
                number_of_terms_in_epsilon_expansion: Some(2),
            }
            .run(&mut state, &cli_settings, &runtime)?;
            let saved: std::collections::BTreeMap<String, String> = toml::from_str(
                &std::fs::read_to_string(temp.path().join(format!("{name}.toml")))?,
            )?;
            assert_eq!(
                saved,
                std::collections::BTreeMap::from([(
                    "expression".to_owned(),
                    output.to_canonical_string(),
                )]),
            );
            let output =
                NumeratorInterface::admit(output, None)?.contract(ContractSettings::default())?;
            assert!(output.is_scalar());
            results.push(output.into_expression());
        }
        assert!(!results[1].is_zero());
        assert!((&results[0] - &results[1]).together().is_zero());
        // Exact contraction happens in D, with no physical-dimension shortcut.
        let source_metric = NumeratorInterface::admit(
            parse!("spenso::g(spenso::mink(gammalooprs::dim,0),spenso::mink(gammalooprs::dim,1))"),
            None,
        )?;
        assert_eq!(
            source_metric
                .multiply(&source_metric)?
                .contract(ContractSettings::default())?
                .expression(),
            &Atom::var(GS.dim)
        );
        Ok(())
    }

    #[test]
    fn evaluate_rejects_generated_thermal_graphs_after_defaults_change() -> color_eyre::Result<()> {
        test_initialise()?;
        let temp = tempfile::tempdir()?;
        for mode in [
            MediumMode::ThermodynamicEquilibrium,
            MediumMode::ZeroTemperatureEquilibrium,
        ] {
            let mut state = State::new_test();
            let mut cli_settings = CLISettings::default();
            cli_settings.state.folder = temp.path().join("state");
            let generation = &mut cli_settings.global.generation;
            generation.medium.mode = mode;
            generation.evaluator.compile = false;
            generation.evaluator.summed = false;
            generation.evaluator.iterative_orientation_optimization = false;
            generation.uv.subtract_uv = false;
            generation.uv.generate_integrated = false;
            generation.threshold_subtraction.enable_thresholds = false;
            generation
                .tropical_subgraph_table
                .disable_tropical_generation = true;
            ImportModel {
                path: "scalars-default.json".into(),
                simplify_model: true,
            }
            .run(&mut state)?;
            Import::Graphs {
                source: None,
                inline_dot: Some(
                    r#"digraph thermal_tadpole {
                        num=1; overall_factor=1; overall_factor_evaluated=1; projector=1;
                        node [num=1];
                        a:0 -> a:1 [particle="scalar_1" num=1 id=0 lmb_id=0,
                            lmb_rep="K(0,a___)", source="{ufo_order:0}", sink="{ufo_order:1}"];
                    }"#
                    .to_string(),
                ),
                process_spec: None,
                process: None,
                integrand_name: Some("tadpole".to_string()),
                overwrite: false,
                append: false,
            }
            .run(&mut state, &cli_settings)?;
            let runtime = RuntimeSettings::default();
            Generate {
                keep_sources: false,
                mode: None,
            }
            .run(
                &mut state,
                temp.path(),
                false,
                &cli_settings.global,
                &runtime,
            )?;

            // Evaluation must use the generated graph's mode, even if subsequent
            // commands select vacuum defaults for future generation.
            cli_settings.global.generation.medium.mode = MediumMode::Vacuum;
            let result_path = temp.path().join("result.toml");
            for numerical in [false, true] {
                let err = Evaluate {
                    process: None,
                    graphs_group_name: Some("tadpole".to_string()),
                    result_path: Some(result_path.clone()),
                    numerical,
                    number_of_terms_in_epsilon_expansion: None,
                }
                .run(&mut state, &cli_settings, &runtime)
                .unwrap_err();
                assert!(
                    err.to_string().contains("does not support thermal modes"),
                    "{mode:?}, numerical={numerical}: {err:?}"
                );
                assert!(!result_path.exists());
            }
        }
        Ok(())
    }

    #[test]
    fn evaluate_rejects_result_path_inside_read_only_state() {
        let mut cli_settings = CLISettings::default();
        cli_settings.state.folder = PathBuf::from("/tmp/read_only_state");
        cli_settings.session.read_only_state = true;

        let err = Evaluate {
            process: None,
            graphs_group_name: None,
            result_path: Some(PathBuf::from("/tmp/read_only_state/result.toml")),
            numerical: false,
            number_of_terms_in_epsilon_expansion: None,
        }
        .run(
            &mut State::new_test(),
            &cli_settings,
            &RuntimeSettings::default(),
        )
        .unwrap_err();

        assert!(format!("{err:?}").contains("--read-only-state"));
    }
}
