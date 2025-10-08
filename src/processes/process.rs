use ahash::HashMap;
use rayon::ThreadPool;
use std::{
    collections::BTreeMap,
    fs::{self, File},
    io::Write,
    path::{Path, PathBuf},
};
// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::{Help, Result};
use itertools::Itertools;
use log::debug;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::fmt;

use crate::{
    feyngen::NumeratorAwareGraphGroupingOption,
    gammaloop_integrand::NewIntegrand,
    numerator::GlobalPrefactor,
    settings::{runtime::LockedRuntimeSettings, GlobalSettings},
    GammaLoopContext, GammaLoopContextContainer,
};
use eyre::{eyre, Context};

use crate::{
    feyngen::{FeynGenFilters, GenerationType},
    graph::Graph,
    model::Model,
    settings::global::GenerationSettings,
};

use super::{Amplitude, CrossSection};

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ProcessDefinition {
    pub generation_type: GenerationType,
    pub initial_pdgs: Vec<i64>,
    pub final_pdgs_lists: Vec<Vec<i64>>,
    pub loop_count_range: (usize, usize),
    pub symmetrize_initial_states: bool,
    pub symmetrize_final_states: bool,
    pub symmetrize_left_right_states: bool,
    pub allow_symmetrization_of_external_fermions_in_amplitudes: bool,
    pub max_multiplicity_for_fast_cut_filter: usize,
    pub amplitude_filters: FeynGenFilters,
    pub cross_section_filters: FeynGenFilters,
    pub folder_name: String,
    pub process_id: usize,
    pub numerator_grouping: NumeratorAwareGraphGroupingOption,
    pub filter_self_loop: bool,
    pub filter_zero_flow_edges: bool,
    pub graph_prefix: String,
    pub selected_graphs: Option<Vec<String>>,
    pub vetoed_graphs: Option<Vec<String>>,
    pub loop_momentum_bases: Option<HashMap<String, Vec<String>>>,
    pub prefactor: GlobalPrefactor,
}

impl fmt::Display for ProcessDefinition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Process #{}: '{}'\nGeneration type: {}{}{}\nInitial PDGs: {:?}{}\nFinal PDGs: {}{}\nLoop count: {}\nAmplitude filters:{}{}\nCross-section filters:{}{}",
            self.process_id,
            self.folder_name,
            self.generation_type,
            if self.symmetrize_left_right_states { " (left-right symmetrized)" } else { "" },
            if self.allow_symmetrization_of_external_fermions_in_amplitudes
                && self.generation_type == GenerationType::Amplitude
                && (self.symmetrize_initial_states  || self.symmetrize_final_states || self.symmetrize_left_right_states)
                { " (allowing fermion symmetrization)" } else { "" },
            self.initial_pdgs,
            if self.symmetrize_initial_states { " (symmetrized)" } else { "" },
            if self.final_pdgs_lists.len() == 1 {
                format!("{:?}",self.final_pdgs_lists[0])
            } else {
                format!("[ {} ]", self.final_pdgs_lists.iter().map(|pdgs| format!("{:?}", pdgs)).join(" | "))
            },
            if self.symmetrize_final_states { " (symmetrized)" } else { "" },
            if self.loop_count_range.0 == self.loop_count_range.1 {
                format!("{}", self.loop_count_range.0)
            } else {
                format!("{:?}", self.loop_count_range)
            },
            if self.amplitude_filters.0.is_empty() { " None" } else {"\n"},
            if self.amplitude_filters.0.is_empty() { "".into() } else { self.amplitude_filters
                .0
                .iter()
                .map(|f| format!(" > {}", f))
                .collect::<Vec<String>>()
                .join("\n")
            },
            if self.cross_section_filters.0.is_empty() { " None" } else {"\n"},
            if self.cross_section_filters.0.is_empty() { "".into() } else { self.cross_section_filters
                .0
                .iter()
                .map(|f| format!(" > {}", f))
                .collect::<Vec<String>>()
                .join("\n")
            }
        )
    }
}

impl Default for ProcessDefinition {
    fn default() -> Self {
        Self {
            generation_type: GenerationType::Amplitude,
            initial_pdgs: vec![],
            final_pdgs_lists: vec![],
            loop_count_range: (1, 1),
            symmetrize_initial_states: false,
            symmetrize_final_states: false,
            symmetrize_left_right_states: false,
            allow_symmetrization_of_external_fermions_in_amplitudes: false,
            max_multiplicity_for_fast_cut_filter: 6,
            amplitude_filters: FeynGenFilters(vec![]),
            cross_section_filters: FeynGenFilters(vec![]),
            folder_name: "undefined_process".to_string(),
            process_id: 0,
            numerator_grouping: NumeratorAwareGraphGroupingOption::NoGrouping,
            filter_self_loop: true,
            graph_prefix: "GL".to_string(),
            selected_graphs: None,
            vetoed_graphs: None,
            loop_momentum_bases: None,
            prefactor: GlobalPrefactor::default(),
            filter_zero_flow_edges: true,
        }
    }
}

impl ProcessDefinition {
    // Best attempt at creating what process definition matches the given graphs
    pub fn from_graph_list(_graphs: &[Graph], generation_type: GenerationType) -> Result<Self> {
        // TODO: At least set correctly things like loop count range and initial/final states
        Ok(Self {
            generation_type,
            ..Self::default()
        })
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct Process {
    pub definition: ProcessDefinition,
    pub settings_history: Option<GlobalSettings>,
    pub collection: ProcessCollection,
}

impl Process {
    pub(crate) fn warm_up(&mut self, model: &Model) -> Result<()> {
        self.collection.warm_up(model)
    }
    pub fn preprocess(
        &mut self,
        model: &Model,
        settings: &GlobalSettings,
        thread_pool: &ThreadPool,
    ) -> Result<()> {
        self.collection
            .preprocess(model, &self.definition, &settings.generation, thread_pool)?;
        self.settings_history = Some(settings.clone());
        Ok(())
    }
}

impl Process {
    pub(crate) fn load_amplitude(
        path: impl AsRef<Path>,
        context: GammaLoopContextContainer,
    ) -> Result<Self> {
        let binary = fs::read(path.as_ref().join("def.bin")).context(format!(
            "Error reading def.bin in {}",
            path.as_ref().display()
        ))?;

        let settings_history = File::open(path.as_ref().join("settings_history.yaml"))
            .ok()
            .map(|a| serde_yaml::from_reader(a))
            .transpose()?;

        let (definition, _) =
            bincode::decode_from_slice_with_context(&binary, bincode::config::standard(), context)
                .context("Error decoding process definition")?;

        let mut collection = ProcessCollection::new_amplitude();

        for entry in fs::read_dir(path.as_ref())
            .context(format!("Error reading {}", path.as_ref().display()))?
        {
            let Ok(entry) = entry else {
                debug!("Error reading entry");
                continue;
            };
            if entry.file_type()?.is_file() {
                continue; //skip def.bin, amplitudes are in folders
            }
            let path = entry.path();
            debug!("loading amplitude at {}", path.display());
            let amp = Amplitude::load(path, context.clone()).context("Error loading amplitude")?;

            collection.add_amplitude(amp);
        }

        Ok(Self {
            definition,
            collection,
            settings_history,
        })
    }

    pub(crate) fn load_cross_section(
        path: impl AsRef<Path>,
        context: GammaLoopContextContainer,
    ) -> Result<Self> {
        let binary = fs::read(path.as_ref().join("def.bin"))?;
        let (definition, _) =
            bincode::decode_from_slice_with_context(&binary, bincode::config::standard(), context)?;

        let collection = ProcessCollection::new_cross_section();
        let settings_history = File::open(path.as_ref().join("settings_history.yaml"))
            .ok()
            .map(|a| serde_yaml::from_reader(a))
            .transpose()?;

        // for entry in fs::read_dir(path)? {
        //     let entry = entry?;
        //     let path = entry.path();
        //     let amp = CrossSection::load(path, context.clone())?;
        //     collection.add_amplitude(amp);
        // }

        Ok(Self {
            definition,
            collection,
            settings_history,
        })
    }

    pub fn save(&mut self, path: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        match &mut self.collection {
            ProcessCollection::Amplitudes(a) => {
                let p = path.as_ref().join("amplitudes");
                fs::create_dir_all(&p)?;
                let p = p.join(PathBuf::from(self.definition.folder_name.clone()));

                let r = fs::create_dir_all(&p).with_context(|| {
                    format!(
                        "Trying to create directory to export amplitude dot {}",
                        p.display()
                    )
                });
                if override_existing {
                    r?;
                }

                let binary = bincode::encode_to_vec(&self.definition, bincode::config::standard())?;
                fs::write(p.join("def.bin"), binary)?;

                if let Some(a) = &self.settings_history {
                    File::create(path.as_ref().join("settings_history.toml"))?
                        .write(&toml::to_string_pretty(a)?.into_bytes())?;
                }

                for (_, amp) in a {
                    amp.save(&p, override_existing)?;
                }
            }
            ProcessCollection::CrossSections(_a) => {
                todo!("Implement save for cross sections");
            }
        }

        Ok(())
    }

    pub fn compile(
        &mut self,
        path: impl AsRef<Path>,
        override_existing: bool,
        settings: &GlobalSettings,
        integrand_name: Option<String>,
        thread_pool: &ThreadPool,
    ) -> Result<()> {
        match &mut self.collection {
            ProcessCollection::Amplitudes(a) => {
                let p = path.as_ref().join("amplitudes");
                fs::create_dir_all(&p)?;
                let p = p.join(PathBuf::from(self.definition.folder_name.clone()));

                let r = fs::create_dir_all(&p).with_context(|| {
                    format!(
                        "Trying to create directory to export amplitude dot {}",
                        p.display()
                    )
                });
                if override_existing {
                    r?;
                }

                for (_, amp) in a {
                    if let Some(int_name) = integrand_name.clone() {
                        if amp.name != int_name {
                            continue;
                        }
                    }

                    amp.compile(&p, override_existing, settings, thread_pool)?;
                }
            }
            ProcessCollection::CrossSections(_a) => {
                todo!("Implement compile for cross sections");
            }
        }

        Ok(())
    }

    pub fn get_integrand(&self, integrand_name: impl AsRef<str>) -> Result<&NewIntegrand> {
        self.collection.get_integrand(integrand_name)
    }

    pub fn get_integrand_names(&self) -> Vec<&str> {
        self.collection.get_integrand_names()
    }

    pub fn get_integrand_mut(
        &mut self,
        integrand_name: impl AsRef<str>,
    ) -> Result<&mut NewIntegrand> {
        self.collection.get_integrand_mut(integrand_name)
    }

    pub(crate) fn export_dot(&self, path: impl AsRef<Path>) -> Result<()> {
        let p = path.as_ref().join("amplitudes");
        let path = p.join(PathBuf::from(self.definition.folder_name.clone()));
        fs::create_dir_all(&path)?;

        match &self.collection {
            ProcessCollection::Amplitudes(a) => {
                for (_, amp) in a {
                    let mut dot = File::create_new(path.join(&format!("{}.dot", amp.name)))
                        .with_context(|| {
                            format!(
                                "Trying to create file to export amplitude dot {}",
                                p.join(&format!("{}.dot", amp.name)).display()
                            )
                        })?;
                    amp.write_dot(&mut dot)?;
                }
            }
            ProcessCollection::CrossSections(cs) => {
                for (_, cs) in cs {
                    let mut dot = File::create_new(path.join(&format!("{}.dot", cs.name)))
                        .with_context(|| {
                            format!(
                                "Trying to create file to export cross section dot {}",
                                p.join(&format!("{}.dot", cs.name)).display()
                            )
                        })?;
                    cs.write_dot(&mut dot)?;
                }
            }
        }
        Ok(())
    }

    pub fn from_graph_list(
        process_name: String,
        integrand_name: String,
        graphs: Vec<Graph>,
        generation_type: GenerationType,
        definition: Option<ProcessDefinition>,
        sub_classes: Option<Vec<Vec<String>>>,
    ) -> Result<Self> {
        let mut proc_definition = definition.unwrap_or_default();
        proc_definition.folder_name = process_name;
        match generation_type {
            GenerationType::Amplitude => {
                let mut collection: ProcessCollection = ProcessCollection::new_amplitude();

                if let Some(_sub_classes) = sub_classes {
                    todo!("implement seperation of processes into user defined sub classes");
                } else {
                    collection.add_amplitude(Amplitude::from_graph_list(integrand_name, graphs)?);

                    // TODO: construct a better default definition from graph (i.e. at least the external IDs)
                    Ok(Self {
                        settings_history: None,
                        definition: proc_definition,
                        collection,
                    })
                }
            }
            GenerationType::CrossSection => {
                let mut collection: ProcessCollection = ProcessCollection::new_cross_section();

                if let Some(_sub_classes) = sub_classes {
                    todo!("implement seperation of processes into user defined sub classes");
                } else {
                    collection
                        .add_cross_section(CrossSection::from_graph_list(integrand_name, graphs)?);
                    // TODO: construct a better default definition from graph (i.e. at least the external IDs)
                    Ok(Self {
                        settings_history: None,
                        definition: proc_definition,
                        collection,
                    })
                }
            }
        }
    }

    pub fn generate_integrands(
        &mut self,
        model: &Model,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
        thread_pool: &ThreadPool,
    ) -> Result<()> {
        self.collection
            .generate_integrands(model, global_settings, runtime_default, thread_pool)
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub enum ProcessCollection {
    Amplitudes(BTreeMap<String, Amplitude>),
    CrossSections(BTreeMap<String, CrossSection>),
}

impl ProcessCollection {
    fn new_amplitude() -> Self {
        Self::Amplitudes(BTreeMap::new())
    }

    pub fn get_integrand_names(&self) -> Vec<&str> {
        match self {
            Self::Amplitudes(amplitudes) => amplitudes.keys().map(|a| a.as_str()).collect(),
            Self::CrossSections(cross_sections) => {
                cross_sections.keys().map(|a| a.as_str()).collect()
            }
        }
    }

    fn get_integrand(&self, name: impl AsRef<str>) -> Result<&NewIntegrand> {
        let res = match self {
            Self::Amplitudes(amplitudes) => {
                if amplitudes.contains_key(name.as_ref()) {
                    Ok(amplitudes.get(name.as_ref()).unwrap().integrand.as_ref())
                } else {
                    let names = amplitudes.keys().map(|a| a.as_str()).collect::<Vec<_>>();

                    Err(eyre!("Integrand {} does not exist", name.as_ref()))
                        .suggestion(format!("Available amplitude names: {}", names.join(", ")))
                }
            }
            Self::CrossSections(cross_sections) => {
                if cross_sections.contains_key(name.as_ref()) {
                    Ok(cross_sections
                        .get(name.as_ref())
                        .unwrap()
                        .integrand
                        .as_ref())
                } else {
                    let names = cross_sections
                        .keys()
                        .map(|a| a.as_str())
                        .collect::<Vec<_>>();

                    Err(eyre!("Integrand {} does not exist", name.as_ref())).suggestion(format!(
                        "Available cross section names: {}",
                        names.join(", ")
                    ))
                }
            }
        }?;

        match res {
            Some(integrand) => Ok(integrand),
            None => Err(eyre!(
                "Integrand {} has not yet been generated, but exists",
                name.as_ref()
            )),
        }
    }

    pub fn find_integrand(&self, name: Option<String>) -> Result<String> {
        let all_integrand_names = self.get_integrand_names();

        let integrand_name = if let Some(name) = name {
            if !all_integrand_names.contains(&name.as_str()) {
                return Err(color_eyre::eyre::eyre!(
                    "No integrand named '{}' in process, Available integrands: {:?}",
                    name,
                    all_integrand_names
                ));
            }
            name
        } else {
            if all_integrand_names.len() != 1 {
                return Err(color_eyre::eyre::eyre!(
                    "Multiple integrands in process,Please specify one of: {:?}",
                    all_integrand_names
                ));
            }
            all_integrand_names[0].to_string()
        };

        Ok(integrand_name)
    }

    fn get_integrand_mut(&mut self, name: impl AsRef<str>) -> Result<&mut NewIntegrand> {
        let res = match self {
            Self::Amplitudes(amplitudes) => {
                if amplitudes.contains_key(name.as_ref()) {
                    Ok(amplitudes
                        .get_mut(name.as_ref())
                        .unwrap()
                        .integrand
                        .as_mut())
                } else {
                    let names = amplitudes.keys().map(|a| a.as_str()).collect::<Vec<_>>();

                    Err(eyre!("Integrand {} does not exist", name.as_ref()))
                        .suggestion(format!("Available amplitude names: {}", names.join(", ")))
                }
            }
            Self::CrossSections(cross_sections) => {
                if cross_sections.contains_key(name.as_ref()) {
                    Ok(cross_sections
                        .get_mut(name.as_ref())
                        .unwrap()
                        .integrand
                        .as_mut())
                } else {
                    let names = cross_sections
                        .keys()
                        .map(|a| a.as_str())
                        .collect::<Vec<_>>();

                    Err(eyre!("Integrand {} does not exist", name.as_ref())).suggestion(format!(
                        "Available cross section names: {}",
                        names.join(", ")
                    ))
                }
            }
        }?;

        match res {
            Some(integrand) => Ok(integrand),
            None => Err(eyre!(
                "Integrand {} has not yet been generated, but exists",
                name.as_ref()
            )),
        }
    }

    fn new_cross_section() -> Self {
        Self::CrossSections(BTreeMap::new())
    }

    pub fn add_amplitude(&mut self, amplitude: Amplitude) {
        match self {
            Self::Amplitudes(amplitudes) => amplitudes.insert(amplitude.name.clone(), amplitude),
            _ => panic!("Cannot add amplitude to a cross section collection"),
        };
    }

    pub fn add_cross_section(&mut self, cross_section: CrossSection) {
        match self {
            Self::CrossSections(cross_sections) => {
                cross_sections.insert(cross_section.name.clone(), cross_section);
            }
            _ => panic!("Cannot add cross section to an amplitude collection"),
        }
    }

    fn preprocess(
        &mut self,
        model: &Model,
        process_definition: &ProcessDefinition,
        settings: &GenerationSettings,
        thread_pool: &ThreadPool,
    ) -> Result<()> {
        match self {
            Self::Amplitudes(amplitudes) => {
                for (_, amplitude) in amplitudes {
                    amplitude.preprocess(model, settings, thread_pool)?;
                }
            }
            Self::CrossSections(cross_sections) => {
                for (_, cross_section) in cross_sections {
                    cross_section.preprocess(model, process_definition)?;
                }
            }
        }
        Ok(())
    }

    fn warm_up(&mut self, model: &Model) -> Result<()> {
        match self {
            Self::Amplitudes(amplitudes) => {
                for (_, amplitude) in amplitudes {
                    amplitude.warm_up(model)?;
                }
            }
            Self::CrossSections(cross_sections) => {
                for (_, cross_section) in cross_sections {
                    cross_section.warm_up(model)?;
                }
            }
        }
        Ok(())
    }

    fn generate_integrands(
        &mut self,
        model: &Model,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
        thread_pool: &ThreadPool,
    ) -> Result<()> {
        // let mut result = HashMap::default();
        match self {
            Self::Amplitudes(amplitudes) => {
                for (_, amplitude) in amplitudes {
                    amplitude.build_integrand(
                        model,
                        global_settings,
                        runtime_default,
                        thread_pool,
                    )?;
                }
            }
            Self::CrossSections(cross_sections) => {
                for (_, cross_section) in cross_sections {
                    cross_section.build_integrand(model, runtime_default)?;
                }
            }
        }
        // result
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::{utils::test_utils::load_generic_model, GammaLoopContextContainer};

    #[test]
    fn test_proc_definition_encode() {
        let def = super::ProcessDefinition::default();
        let encoded = bincode::encode_to_vec(&def, bincode::config::standard()).unwrap();
        let model_sm = load_generic_model("sm");

        let mut state_file = std::fs::File::create("state_map.bin").unwrap();

        let _exported = symbolica::state::State::export(&mut state_file).unwrap();
        let state_map = symbolica::state::State::import(&mut state_file, None).unwrap();

        let context = GammaLoopContextContainer {
            model: &model_sm,
            state_map: &state_map,
        };

        let (decoded, _): (super::ProcessDefinition, _) =
            bincode::decode_from_slice_with_context(&encoded, bincode::config::standard(), context)
                .unwrap();
        assert_eq!(def, decoded);
    }
}
