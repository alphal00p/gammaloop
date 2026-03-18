use std::collections::BTreeMap;

use clap::{Args, Subcommand};
use color_eyre::{eyre::eyre, Result};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use gammalooprs::{
    integrands::process::ProcessIntegrand,
    processes::{Amplitude, CrossSection, Process, ProcessCollection},
};

use crate::{
    completion::CompletionArgExt,
    state::{ProcessRef, State},
};

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Duplicate {
    Integrand(DuplicateIntegrand),
}

impl Duplicate {
    pub fn run(self, state: &mut State) -> Result<()> {
        match self {
            Duplicate::Integrand(command) => command.run(state),
        }
    }
}

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct DuplicateIntegrand {
    /// Source process reference: #<id>, name:<name>, or <id>/<name>
    #[arg(
        short = 'p',
        long = "process",
        value_name = "PROCESS",
        completion_process_selector(crate::completion::SelectorKind::Any)
    )]
    pub process: Option<ProcessRef>,

    /// Source integrand name
    #[arg(
        short = 'i',
        long = "integrand-name",
        value_name = "NAME",
        completion_integrand_selector(crate::completion::SelectorKind::Any)
    )]
    pub integrand_name: Option<String>,

    /// Destination process name
    #[arg(long = "output_process_name", value_name = "NAME")]
    pub output_process_name: String,

    /// Destination integrand name
    #[arg(long = "output_integrand_name", value_name = "NAME")]
    pub output_integrand_name: String,
}

impl DuplicateIntegrand {
    pub fn run(self, state: &mut State) -> Result<()> {
        let (source_process_id, source_integrand_name) =
            state.find_integrand_ref(self.process.as_ref(), self.integrand_name.as_ref())?;

        if self.output_process_name.trim().is_empty()
            || self.output_integrand_name.trim().is_empty()
        {
            return Err(eyre!(
                "Both --output_process_name and --output_integrand_name must be provided"
            ));
        }

        let source_process = &state.process_list.processes[source_process_id];
        let existing_destination = state
            .process_list
            .processes
            .iter()
            .position(|process| process.definition.folder_name == self.output_process_name);
        if let Some(destination_process_id) = existing_destination {
            let destination_process = &state.process_list.processes[destination_process_id];
            if destination_process
                .collection
                .get_integrand_names()
                .contains(&self.output_integrand_name.as_str())
            {
                return Err(eyre!(
                    "An integrand '{}' already exists in process '{}'",
                    self.output_integrand_name,
                    self.output_process_name
                ));
            }
        }

        match &source_process.collection {
            ProcessCollection::Amplitudes(amplitudes) => {
                let mut duplicated = amplitudes
                    .get(&source_integrand_name)
                    .cloned()
                    .ok_or_else(|| eyre!("Missing source amplitude '{}'", source_integrand_name))?;
                duplicated.name = self.output_integrand_name.clone();
                rename_process_integrand(
                    duplicated.integrand.as_mut(),
                    &self.output_integrand_name,
                );
                self.insert_amplitude_copy(state, source_process_id, duplicated)
            }
            ProcessCollection::CrossSections(cross_sections) => {
                let mut duplicated = cross_sections
                    .get(&source_integrand_name)
                    .cloned()
                    .ok_or_else(|| {
                        eyre!("Missing source cross section '{}'", source_integrand_name)
                    })?;
                duplicated.name = self.output_integrand_name.clone();
                rename_process_integrand(
                    duplicated.integrand.as_mut(),
                    &self.output_integrand_name,
                );
                self.insert_cross_section_copy(state, source_process_id, duplicated)
            }
        }
    }

    fn insert_amplitude_copy(
        &self,
        state: &mut State,
        source_process_id: usize,
        amplitude: Amplitude,
    ) -> Result<()> {
        if let Some(destination_process_id) = state
            .process_list
            .processes
            .iter()
            .position(|process| process.definition.folder_name == self.output_process_name)
        {
            let destination_process = &mut state.process_list.processes[destination_process_id];
            match &mut destination_process.collection {
                ProcessCollection::Amplitudes(collection) => {
                    collection.insert(amplitude.name.clone(), amplitude);
                    Ok(())
                }
                ProcessCollection::CrossSections(_) => Err(eyre!(
                    "Destination process '{}' exists but does not contain amplitudes",
                    self.output_process_name
                )),
            }
        } else {
            let mut definition = state.process_list.processes[source_process_id]
                .definition
                .clone();
            definition.folder_name = self.output_process_name.clone();
            definition.process_id = next_process_id(state);
            let mut collection = ProcessCollection::Amplitudes(BTreeMap::new());
            collection.add_amplitude(amplitude);
            state.process_list.add_process(Process {
                definition,
                settings_history: state.process_list.processes[source_process_id]
                    .settings_history
                    .clone(),
                collection,
            });
            Ok(())
        }
    }

    fn insert_cross_section_copy(
        &self,
        state: &mut State,
        source_process_id: usize,
        cross_section: CrossSection,
    ) -> Result<()> {
        if let Some(destination_process_id) = state
            .process_list
            .processes
            .iter()
            .position(|process| process.definition.folder_name == self.output_process_name)
        {
            let destination_process = &mut state.process_list.processes[destination_process_id];
            match &mut destination_process.collection {
                ProcessCollection::CrossSections(collection) => {
                    collection.insert(cross_section.name.clone(), cross_section);
                    Ok(())
                }
                ProcessCollection::Amplitudes(_) => Err(eyre!(
                    "Destination process '{}' exists but does not contain cross sections",
                    self.output_process_name
                )),
            }
        } else {
            let mut definition = state.process_list.processes[source_process_id]
                .definition
                .clone();
            definition.folder_name = self.output_process_name.clone();
            definition.process_id = next_process_id(state);
            let mut collection = ProcessCollection::CrossSections(BTreeMap::new());
            collection.add_cross_section(cross_section);
            state.process_list.add_process(Process {
                definition,
                settings_history: state.process_list.processes[source_process_id]
                    .settings_history
                    .clone(),
                collection,
            });
            Ok(())
        }
    }
}

fn rename_process_integrand(integrand: Option<&mut ProcessIntegrand>, new_name: &str) {
    let Some(integrand) = integrand else {
        return;
    };
    let _ = integrand.get_mut_settings();
    match integrand {
        ProcessIntegrand::Amplitude(amplitude) => {
            amplitude.data.name = new_name.to_string();
        }
        ProcessIntegrand::CrossSection(cross_section) => {
            cross_section.data.name = new_name.to_string();
        }
    }
}

fn next_process_id(state: &State) -> usize {
    state
        .process_list
        .processes
        .iter()
        .map(|process| process.definition.process_id)
        .max()
        .map(|max_id| max_id + 1)
        .unwrap_or(0)
}
