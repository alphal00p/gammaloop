use clap::Args;
use color_eyre::Result;
use colored::Colorize;
use eyre::{eyre, Context};
use gammalooprs::processes::{
    CycleSignature, GraphGroupSelectionMode, GraphGroupSelectionSpec, ParticleSignature,
    RaisedPropagatorScope, RaisedPropagatorSignature, SelectionPolarity, VertexSignature,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::str::FromStr;
use tracing::info;

use crate::{
    completion::CompletionArgExt,
    state::{
        GraphGroupSelectionContext, GraphGroupSelectionTarget, ProcessRef, SelectedGraphGroups,
        State,
    },
    CLISettings,
};

#[derive(Args, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Select {
    /// Process reference: #<id>, name:<name>, or <id>/<name>
    #[arg(
        long = "process",
        short = 'p',
        value_name = "PROCESS",
        completion_process_selector(crate::completion::SelectorKind::Any)
    )]
    pub process: Option<ProcessRef>,

    /// Integrand name inside the selected process
    #[arg(
        long = "integrand-name",
        short = 'i',
        value_name = "NAME",
        completion_integrand_selector(crate::completion::SelectorKind::Any)
    )]
    pub integrand_name: Option<String>,

    /// Apply structural filters to left/right amplitude graphs of cross-section cuts
    #[arg(long = "amplitude-graphs", default_value_t = false)]
    pub amplitude_graphs: bool,

    /// Retain graph groups with these master graph names
    #[arg(
        long = "with-graph-names",
        value_name = "MASTER_GRAPH",
        num_args = 1..,
        completion_selected_master_graph()
    )]
    pub with_graph_names: Vec<String>,

    /// Remove graph groups with these master graph names
    #[arg(
        long = "without-graph-names",
        value_name = "MASTER_GRAPH",
        num_args = 1..,
        completion_selected_master_graph()
    )]
    pub without_graph_names: Vec<String>,

    /// Retain graph groups with one of these raised-propagator signatures
    #[arg(
        long = "with-raised-propagator-signatures",
        value_name = "SIGNATURE",
        num_args = 1..,
        completion_select_raised_signature(crate::completion::SelectRaisedSignatureScope::All)
    )]
    pub with_raised_propagator_signatures: Vec<String>,

    /// Remove graph groups with one of these raised-propagator signatures
    #[arg(
        long = "without-raised-propagator-signatures",
        value_name = "SIGNATURE",
        num_args = 1..,
        completion_select_raised_signature(crate::completion::SelectRaisedSignatureScope::All)
    )]
    pub without_raised_propagator_signatures: Vec<String>,

    /// Retain graph groups with one of these massive raised-propagator signatures
    #[arg(
        long = "with-massive-raised-propagator-signatures",
        value_name = "SIGNATURE",
        num_args = 1..,
        completion_select_raised_signature(crate::completion::SelectRaisedSignatureScope::Massive)
    )]
    pub with_massive_raised_propagator_signatures: Vec<String>,

    /// Remove graph groups with one of these massive raised-propagator signatures
    #[arg(
        long = "without-massive-raised-propagator-signatures",
        value_name = "SIGNATURE",
        num_args = 1..,
        completion_select_raised_signature(crate::completion::SelectRaisedSignatureScope::Massive)
    )]
    pub without_massive_raised_propagator_signatures: Vec<String>,

    /// Retain graph groups with one of these massless raised-propagator signatures
    #[arg(
        long = "with-massless-raised-propagator-signatures",
        value_name = "SIGNATURE",
        num_args = 1..,
        completion_select_raised_signature(crate::completion::SelectRaisedSignatureScope::Massless)
    )]
    pub with_massless_raised_propagator_signatures: Vec<String>,

    /// Remove graph groups with one of these massless raised-propagator signatures
    #[arg(
        long = "without-massless-raised-propagator-signatures",
        value_name = "SIGNATURE",
        num_args = 1..,
        completion_select_raised_signature(crate::completion::SelectRaisedSignatureScope::Massless)
    )]
    pub without_massless_raised_propagator_signatures: Vec<String>,

    /// Retain graph groups where a process-valid cut touches one of these raised-propagator group signatures
    #[arg(
        long = "with-raised-cuts-signatures",
        value_name = "SIGNATURE",
        num_args = 1..,
        completion_select_raised_cut_signature(crate::completion::SelectRaisedSignatureScope::All)
    )]
    pub with_raised_cuts_signatures: Vec<String>,

    /// Remove graph groups where a process-valid cut touches one of these raised-propagator group signatures
    #[arg(
        long = "without-raised-cuts-signatures",
        value_name = "SIGNATURE",
        num_args = 1..,
        completion_select_raised_cut_signature(crate::completion::SelectRaisedSignatureScope::All)
    )]
    pub without_raised_cuts_signatures: Vec<String>,

    /// Retain graph groups where a process-valid cut touches one of these massive raised-propagator group signatures
    #[arg(
        long = "with-massive-raised-cuts-signatures",
        value_name = "SIGNATURE",
        num_args = 1..,
        completion_select_raised_cut_signature(crate::completion::SelectRaisedSignatureScope::Massive)
    )]
    pub with_massive_raised_cuts_signatures: Vec<String>,

    /// Remove graph groups where a process-valid cut touches one of these massive raised-propagator group signatures
    #[arg(
        long = "without-massive-raised-cuts-signatures",
        value_name = "SIGNATURE",
        num_args = 1..,
        completion_select_raised_cut_signature(crate::completion::SelectRaisedSignatureScope::Massive)
    )]
    pub without_massive_raised_cuts_signatures: Vec<String>,

    /// Retain graph groups where a process-valid cut touches one of these massless raised-propagator group signatures
    #[arg(
        long = "with-massless-raised-cuts-signatures",
        value_name = "SIGNATURE",
        num_args = 1..,
        completion_select_raised_cut_signature(crate::completion::SelectRaisedSignatureScope::Massless)
    )]
    pub with_massless_raised_cuts_signatures: Vec<String>,

    /// Remove graph groups where a process-valid cut touches one of these massless raised-propagator group signatures
    #[arg(
        long = "without-massless-raised-cuts-signatures",
        value_name = "SIGNATURE",
        num_args = 1..,
        completion_select_raised_cut_signature(crate::completion::SelectRaisedSignatureScope::Massless)
    )]
    pub without_massless_raised_cuts_signatures: Vec<String>,

    /// Retain graph groups matching one of these cycle signatures
    #[arg(
        long = "with-cycle-signatures",
        value_name = "SIGNATURE",
        num_args = 1..,
        completion_select_cycle_signature()
    )]
    pub with_cycle_signatures: Vec<String>,

    /// Remove graph groups matching one of these cycle signatures
    #[arg(
        long = "without-cycle-signatures",
        value_name = "SIGNATURE",
        num_args = 1..,
        completion_select_cycle_signature()
    )]
    pub without_cycle_signatures: Vec<String>,

    /// Retain graph groups containing at least one of these vertex-rule multisets
    #[arg(
        long = "with-vertices",
        value_name = "VERTEX_LIST",
        num_args = 1..,
        completion_select_vertex_signature()
    )]
    pub with_vertices: Vec<String>,

    /// Remove graph groups containing at least one of these vertex-rule multisets
    #[arg(
        long = "without-vertices",
        value_name = "VERTEX_LIST",
        num_args = 1..,
        completion_select_vertex_signature()
    )]
    pub without_vertices: Vec<String>,

    /// Retain graph groups containing all particles from at least one of these particle sets
    #[arg(
        long = "with-particles",
        value_name = "PARTICLE_SET",
        num_args = 1..,
        completion_select_particle_signature()
    )]
    pub with_particles: Vec<String>,

    /// Remove graph groups containing all particles from at least one of these particle sets
    #[arg(
        long = "without-particles",
        value_name = "PARTICLE_SET",
        num_args = 1..,
        completion_select_particle_signature()
    )]
    pub without_particles: Vec<String>,

    /// Copy the selected graph groups to this output process instead of modifying in place
    #[arg(
        long = "output_process",
        value_name = "NAME",
        completion_disable_special_value()
    )]
    pub output_process: Option<String>,

    /// Copy the selected graph groups to this output integrand instead of modifying in place
    #[arg(
        long = "output_integrand",
        value_name = "NAME",
        completion_disable_special_value()
    )]
    pub output_integrand: Option<String>,

    /// Allow copy-mode selection to replace an existing output process or integrand
    #[arg(
        long = "clear-existing-processes",
        short = 'c',
        alias = "clear",
        default_value_t = false
    )]
    pub clear_existing_processes: bool,
}

impl Select {
    pub fn run(&self, state: &mut State, global_cli_settings: &CLISettings) -> Result<()> {
        let selection = self.selection_spec(state)?;
        let target = self.selection_target();
        let selected = state.select_integrand_graph_groups(
            self.process.as_ref(),
            self.integrand_name.as_ref(),
            &selection,
            &target,
            GraphGroupSelectionContext::new(
                &global_cli_settings.global.generation,
                &global_cli_settings.state.folder,
                global_cli_settings.session.read_only_state,
            ),
        )?;

        report_selection(&selected);
        if selected.copied_to_output {
            if selected.replaced_existing_target {
                info!("{}", "Replaced existing output target.".yellow());
            }
            if selected.removed_target_artifacts {
                info!("{}", "Removed stale saved output artifacts.".yellow());
            }
            info!(
                "{} {}",
                "Created selected pregeneration integrand; generate with".yellow(),
                format!(
                    "generate existing -p #{} -i {}",
                    selected.process_id, selected.integrand_name
                )
                .green()
            );
        } else if selected.discarded_generated_integrand {
            info!(
                "{} {}",
                "Discarded generated integrand state; regenerate with".yellow(),
                format!(
                    "generate existing -p #{} -i {}",
                    selected.process_id, selected.integrand_name
                )
                .green()
            );
        }

        Ok(())
    }

    fn selection_target(&self) -> GraphGroupSelectionTarget {
        GraphGroupSelectionTarget::copy(
            self.output_process.clone(),
            self.output_integrand.clone(),
            self.clear_existing_processes,
        )
    }

    fn selection_spec(&self, state: &State) -> Result<GraphGroupSelectionSpec> {
        let mut spec = GraphGroupSelectionSpec::new()
            .with_mode(if self.amplitude_graphs {
                GraphGroupSelectionMode::CrossSectionAmplitudeGraphs
            } else {
                GraphGroupSelectionMode::MasterGraphs
            })
            .with_master_graph_names_polarity(
                SelectionPolarity::With,
                self.with_graph_names.clone(),
            )
            .with_master_graph_names_polarity(
                SelectionPolarity::Without,
                self.without_graph_names.clone(),
            );

        spec = spec
            .with_raised_propagator_signatures(
                SelectionPolarity::With,
                RaisedPropagatorScope::All,
                parse_raised_signatures(&self.with_raised_propagator_signatures)?,
            )
            .with_raised_propagator_signatures(
                SelectionPolarity::Without,
                RaisedPropagatorScope::All,
                parse_raised_signatures(&self.without_raised_propagator_signatures)?,
            )
            .with_raised_propagator_signatures(
                SelectionPolarity::With,
                RaisedPropagatorScope::Massive,
                parse_raised_signatures(&self.with_massive_raised_propagator_signatures)?,
            )
            .with_raised_propagator_signatures(
                SelectionPolarity::Without,
                RaisedPropagatorScope::Massive,
                parse_raised_signatures(&self.without_massive_raised_propagator_signatures)?,
            )
            .with_raised_propagator_signatures(
                SelectionPolarity::With,
                RaisedPropagatorScope::Massless,
                parse_raised_signatures(&self.with_massless_raised_propagator_signatures)?,
            )
            .with_raised_propagator_signatures(
                SelectionPolarity::Without,
                RaisedPropagatorScope::Massless,
                parse_raised_signatures(&self.without_massless_raised_propagator_signatures)?,
            )
            .with_raised_cut_signatures(
                SelectionPolarity::With,
                RaisedPropagatorScope::All,
                parse_raised_cut_signatures(&self.with_raised_cuts_signatures)?,
            )
            .with_raised_cut_signatures(
                SelectionPolarity::Without,
                RaisedPropagatorScope::All,
                parse_raised_cut_signatures(&self.without_raised_cuts_signatures)?,
            )
            .with_raised_cut_signatures(
                SelectionPolarity::With,
                RaisedPropagatorScope::Massive,
                parse_raised_cut_signatures(&self.with_massive_raised_cuts_signatures)?,
            )
            .with_raised_cut_signatures(
                SelectionPolarity::Without,
                RaisedPropagatorScope::Massive,
                parse_raised_cut_signatures(&self.without_massive_raised_cuts_signatures)?,
            )
            .with_raised_cut_signatures(
                SelectionPolarity::With,
                RaisedPropagatorScope::Massless,
                parse_raised_cut_signatures(&self.with_massless_raised_cuts_signatures)?,
            )
            .with_raised_cut_signatures(
                SelectionPolarity::Without,
                RaisedPropagatorScope::Massless,
                parse_raised_cut_signatures(&self.without_massless_raised_cuts_signatures)?,
            )
            .with_cycle_signatures(
                SelectionPolarity::With,
                parse_cycle_signatures(&self.with_cycle_signatures, state)?,
            )
            .with_cycle_signatures(
                SelectionPolarity::Without,
                parse_cycle_signatures(&self.without_cycle_signatures, state)?,
            )
            .with_vertex_signatures(
                SelectionPolarity::With,
                parse_vertex_signatures(&self.with_vertices, state)?,
            )
            .with_vertex_signatures(
                SelectionPolarity::Without,
                parse_vertex_signatures(&self.without_vertices, state)?,
            )
            .with_particle_signatures(
                SelectionPolarity::With,
                parse_particle_signatures(&self.with_particles, state)?,
            )
            .with_particle_signatures(
                SelectionPolarity::Without,
                parse_particle_signatures(&self.without_particles, state)?,
            );

        if spec.is_empty() {
            return Err(eyre!(
                "No graph-group selection filters were provided. Use --with-graph-names or another --with/--without selection option."
            ));
        }

        Ok(spec)
    }
}

fn parse_raised_signatures(values: &[String]) -> Result<Vec<RaisedPropagatorSignature>> {
    values
        .iter()
        .map(|value| {
            RaisedPropagatorSignature::from_str(value)
                .with_context(|| format!("While parsing raised-propagator signature '{value}'"))
        })
        .collect()
}

fn parse_raised_cut_signatures(values: &[String]) -> Result<Vec<RaisedPropagatorSignature>> {
    values
        .iter()
        .map(|value| {
            RaisedPropagatorSignature::from_str(value)
                .with_context(|| format!("While parsing raised-cut signature '{value}'"))
        })
        .collect()
}

fn parse_cycle_signatures(values: &[String], state: &State) -> Result<Vec<CycleSignature>> {
    values
        .iter()
        .map(|value| CycleSignature::parse(value, &state.model))
        .collect()
}

fn parse_vertex_signatures(values: &[String], state: &State) -> Result<Vec<VertexSignature>> {
    values
        .iter()
        .map(|value| {
            let signature = VertexSignature::parse(value)?;
            for vertex_rule_name in signature.vertex_rule_names() {
                if !state
                    .model
                    .vertex_rule_name_to_position
                    .contains_key(vertex_rule_name)
                {
                    return Err(eyre!(
                        "Unknown vertex rule '{}' in vertex selection signature '{}'.",
                        vertex_rule_name,
                        value
                    ));
                }
            }
            Ok(signature)
        })
        .collect()
}

fn parse_particle_signatures(values: &[String], state: &State) -> Result<Vec<ParticleSignature>> {
    values
        .iter()
        .map(|value| ParticleSignature::parse(value, &state.model))
        .collect()
}

fn report_selection(selected: &SelectedGraphGroups) {
    let report = &selected.report;
    if selected.copied_to_output {
        info!(
            "{} {} {} {} {} {} {} {} {} {}",
            "Selected".blue(),
            report.kept_master_graphs.len().to_string().green(),
            "graph groups from".blue(),
            selected.source_integrand_name.green(),
            "in process".blue(),
            selected.source_process_name.green(),
            "to".blue(),
            selected.integrand_name.green(),
            "in process".blue(),
            selected.process_name.green()
        );
    } else {
        info!(
            "{} {} {} {} {} {}",
            "Selected".blue(),
            report.kept_master_graphs.len().to_string().green(),
            "graph groups for integrand".blue(),
            selected.integrand_name.green(),
            "in process".blue(),
            selected.process_name.green()
        );
    }
    info!(
        "{} {}",
        "Kept master graphs:".blue(),
        format_graph_name_list(&report.kept_master_graphs).green()
    );
    info!(
        "{} {}",
        "Removed master graphs:".blue(),
        format_graph_name_list(&report.removed_master_graphs).red()
    );
    if !report.removed_graphs.is_empty()
        && report.removed_graphs.len() != report.removed_master_graphs.len()
    {
        info!(
            "{} {}",
            "Removed concrete graphs:".blue(),
            format_graph_name_list(&report.removed_graphs).red()
        );
    }
}

fn format_graph_name_list(graph_names: &[String]) -> String {
    match graph_names.len() {
        0 => "0".to_string(),
        1..=10 => format!("{} ({})", graph_names.len(), graph_names.join(", ")),
        len => format!("{len} (list elided; more than 10 graph names)"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn empty_select() -> Select {
        Select {
            process: None,
            integrand_name: None,
            amplitude_graphs: false,
            with_graph_names: Vec::new(),
            without_graph_names: Vec::new(),
            with_raised_propagator_signatures: Vec::new(),
            without_raised_propagator_signatures: Vec::new(),
            with_massive_raised_propagator_signatures: Vec::new(),
            without_massive_raised_propagator_signatures: Vec::new(),
            with_massless_raised_propagator_signatures: Vec::new(),
            without_massless_raised_propagator_signatures: Vec::new(),
            with_raised_cuts_signatures: Vec::new(),
            without_raised_cuts_signatures: Vec::new(),
            with_massive_raised_cuts_signatures: Vec::new(),
            without_massive_raised_cuts_signatures: Vec::new(),
            with_massless_raised_cuts_signatures: Vec::new(),
            without_massless_raised_cuts_signatures: Vec::new(),
            with_cycle_signatures: Vec::new(),
            without_cycle_signatures: Vec::new(),
            with_vertices: Vec::new(),
            without_vertices: Vec::new(),
            with_particles: Vec::new(),
            without_particles: Vec::new(),
            output_process: None,
            output_integrand: None,
            clear_existing_processes: false,
        }
    }

    #[test]
    fn selection_spec_rejects_empty_filters() {
        let state = State::new_test();
        let err = empty_select().selection_spec(&state).unwrap_err();
        assert!(
            err.to_string()
                .contains("No graph-group selection filters were provided"),
            "{err:?}"
        );
    }

    #[test]
    fn selection_spec_accepts_raised_cut_filters() {
        let state = State::new_test();
        let mut select = empty_select();
        select.with_raised_cuts_signatures = vec!["[2,]".to_string()];
        select.with_massive_raised_cuts_signatures = vec!["[]".to_string()];
        select.without_massless_raised_cuts_signatures = vec!["[2]".to_string()];

        let spec = select.selection_spec(&state).unwrap();
        assert!(spec.has_raised_cut_rules());
    }

    #[test]
    fn selection_spec_accepts_any_raising_keyword() {
        let state = State::new_test();
        let mut select = empty_select();
        select.without_raised_propagator_signatures = vec!["ANY_RAISING".to_string()];
        select.without_raised_cuts_signatures = vec!["ANY_RAISING".to_string()];

        let spec = select.selection_spec(&state).unwrap();
        assert!(spec.has_raised_cut_rules());
    }
}
