use ahash::HashMap;
use ahash::HashSet;
use linnet::half_edge::involution::Flow;
use linnet::half_edge::involution::HedgePair;
use linnet::half_edge::involution::Orientation;
use rayon::ThreadPool;
use std::{
    collections::BTreeMap,
    fs::{self, File},
    io::Write,
    path::{Path, PathBuf},
};
use tracing::warn;
// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::{Help, Result};
use itertools::Itertools;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::fmt;
use tracing::debug;

use crate::graph::FeynmanGraph;
use crate::graph::edge::PossibleParticle;
use crate::processes::DotExportSettings;

use crate::processes::StandaloneExportSettings;
use crate::{
    GammaLoopContext, GammaLoopContextContainer,
    feyngen::NumeratorAwareGraphGroupingOption,
    integrands::process::ProcessIntegrand,
    numerator::GlobalPrefactor,
    settings::{GlobalSettings, RuntimeSettings, runtime::LockedRuntimeSettings},
};
use eyre::{Context, eyre};

use crate::{
    feyngen::{FeynGenFilters, GenerationType},
    graph::Graph,
    model::Model,
    settings::global::GenerationSettings,
};

use super::{Amplitude, CrossSection};

const SETTINGS_HISTORY_TOML: &str = "settings_history.toml";
const SETTINGS_HISTORY_YAML: &str = "settings_history.yaml";

pub struct ResolvedIntegrandRef<'a> {
    pub canonical_name: String,
    pub integrand: Option<&'a ProcessIntegrand>,
}

impl<'a> ResolvedIntegrandRef<'a> {
    pub fn get_settings(&self) -> Option<&'a RuntimeSettings> {
        self.integrand.map(ProcessIntegrand::get_settings)
    }

    pub fn require_generated(&self) -> Result<&'a ProcessIntegrand> {
        self.integrand.ok_or_else(|| {
            eyre!(
                "Integrand {} has not yet been generated, but exists",
                self.canonical_name
            )
        })
    }
}

fn create_overwriting_file(path: &Path, file_kind: &str) -> Result<File> {
    if path.exists() {
        if path.is_dir() {
            fs::remove_dir_all(path).with_context(|| {
                format!(
                    "Trying to remove existing directory before exporting {file_kind} {}",
                    path.display()
                )
            })?;
        } else {
            fs::remove_file(path).with_context(|| {
                format!(
                    "Trying to remove existing file before exporting {file_kind} {}",
                    path.display()
                )
            })?;
        }
    }

    File::create(path).with_context(|| {
        format!(
            "Trying to create file to export {file_kind} {}",
            path.display()
        )
    })
}

fn load_settings_history(path: &Path) -> Result<Option<GlobalSettings>> {
    let settings_history_toml = path.join(SETTINGS_HISTORY_TOML);
    if settings_history_toml.exists() {
        let settings_history_raw =
            fs::read_to_string(&settings_history_toml).with_context(|| {
                format!(
                    "Error reading process settings history file {}",
                    settings_history_toml.display()
                )
            })?;
        let settings_history = toml::from_str(&settings_history_raw).with_context(|| {
            format!(
                "Error parsing process settings history file {}",
                settings_history_toml.display()
            )
        })?;
        return Ok(Some(settings_history));
    }

    let settings_history_yaml = path.join(SETTINGS_HISTORY_YAML);
    if settings_history_yaml.exists() {
        warn!(
            "Using legacy process settings history file {}. Re-save state to migrate to {}.",
            settings_history_yaml.display(),
            SETTINGS_HISTORY_TOML
        );
        let settings_history = serde_yaml::from_reader(File::open(&settings_history_yaml)?)
            .with_context(|| {
                format!(
                    "Error parsing legacy process settings history file {}",
                    settings_history_yaml.display()
                )
            })?;
        return Ok(Some(settings_history));
    }

    Ok(None)
}

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
    pub loop_momentum_bases: Option<HashMap<String, Vec<usize>>>,
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
            if self.symmetrize_left_right_states {
                " (left-right symmetrized)"
            } else {
                ""
            },
            if self.allow_symmetrization_of_external_fermions_in_amplitudes
                && self.generation_type == GenerationType::Amplitude
                && (self.symmetrize_initial_states
                    || self.symmetrize_final_states
                    || self.symmetrize_left_right_states)
            {
                " (allowing fermion symmetrization)"
            } else {
                ""
            },
            self.initial_pdgs,
            if self.symmetrize_initial_states {
                " (symmetrized)"
            } else {
                ""
            },
            if self.final_pdgs_lists.len() == 1 {
                format!("{:?}", self.final_pdgs_lists[0])
            } else {
                format!(
                    "[ {} ]",
                    self.final_pdgs_lists
                        .iter()
                        .map(|pdgs| format!("{:?}", pdgs))
                        .join(" | ")
                )
            },
            if self.symmetrize_final_states {
                " (symmetrized)"
            } else {
                ""
            },
            if self.loop_count_range.0 == self.loop_count_range.1 {
                format!("{}", self.loop_count_range.0)
            } else {
                format!("{:?}", self.loop_count_range)
            },
            if self.amplitude_filters.0.is_empty() {
                " None"
            } else {
                "\n"
            },
            if self.amplitude_filters.0.is_empty() {
                "".into()
            } else {
                self.amplitude_filters
                    .0
                    .iter()
                    .map(|f| format!(" > {}", f))
                    .collect::<Vec<String>>()
                    .join("\n")
            },
            if self.cross_section_filters.0.is_empty() {
                " None"
            } else {
                "\n"
            },
            if self.cross_section_filters.0.is_empty() {
                "".into()
            } else {
                self.cross_section_filters
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
    pub fn from_graph_list(
        graphs: &[Graph],
        generation_type: GenerationType,
        model: &Model,
    ) -> Result<Self> {
        let mut initial_pdgs = HashSet::default();

        for g in graphs {
            let mut initial_pdgs_of_graph = vec![];
            match generation_type {
                GenerationType::Amplitude => {
                    for (pair, _, edge) in g.iter_edges() {
                        if matches!(
                            pair,
                            HedgePair::Unpaired {
                                hedge: _,
                                flow: Flow::Sink
                            }
                        ) {
                            if let PossibleParticle::Particle(particle) = &edge.data.particle {
                                initial_pdgs_of_graph.push(particle.0.pdg_code as i64);
                            } else {
                                debug!("Edge without particle data in initial state");
                            }
                        }
                    }
                }
                GenerationType::CrossSection => {
                    for (_, _, edge) in g.iter_edges_of(&g.initial_state_cut) {
                        if let PossibleParticle::Particle(particle) = &edge.data.particle {
                            initial_pdgs_of_graph.push(particle.0.pdg_code as i64);
                        } else {
                            debug!("Edge without particle data in initial state");
                        }
                    }
                }
            }

            initial_pdgs_of_graph.sort();
            initial_pdgs.insert(initial_pdgs_of_graph);
        }

        let initial_pdgs = if initial_pdgs.len() == 1 {
            initial_pdgs.into_iter().next().unwrap()
        } else {
            warn!("Multiple initial states found in graphs, setting initial state to empty");
            vec![]
        };

        let mut final_states = HashSet::default();
        match generation_type {
            GenerationType::Amplitude => {
                for g in graphs {
                    let mut final_pdgs_of_graph = vec![];
                    for (pair, _, edge) in g.iter_edges() {
                        if matches!(
                            pair,
                            HedgePair::Unpaired {
                                hedge: _,
                                flow: Flow::Source
                            }
                        ) {
                            if let PossibleParticle::Particle(particle) = &edge.data.particle {
                                final_pdgs_of_graph.push(particle.0.pdg_code as i64);
                            } else {
                                debug!("Edge without particle data in final state");
                            }
                        }
                    }
                    final_pdgs_of_graph.sort();
                    final_states.insert(final_pdgs_of_graph);
                }
            }
            GenerationType::CrossSection => {
                for g in graphs {
                    let (source_nodes, target_nodes) = g.get_source_and_target();
                    let st_cuts = g.all_st_cuts_for_cs(
                        source_nodes,
                        target_nodes,
                        &g.get_initial_state_tree().0,
                    );
                    for (_, cut, _) in st_cuts {
                        let mut final_pdgs_of_cut = vec![];
                        for (orientaion, edge) in cut.iter_edges(&g.underlying) {
                            if let PossibleParticle::Particle(particle) = &edge.data.particle {
                                if orientaion == Orientation::Reversed {
                                    final_pdgs_of_cut
                                        .push(particle.0.get_anti_particle(model).pdg_code as i64);
                                } else {
                                    final_pdgs_of_cut.push(particle.0.pdg_code as i64);
                                }
                            } else {
                                debug!("Edge without particle data in final state");
                            }
                        }
                        final_pdgs_of_cut.sort();
                        final_states.insert(final_pdgs_of_cut);
                    }
                }
            }
        }

        let final_pdgs_lists = final_states.into_iter().sorted().collect_vec();
        let mut min_loop_count = usize::MAX;
        let mut max_loop_count = 0usize;

        for g in graphs {
            // don't know how the looop count is really intended, for now it doesn't matter I think
            let lc = g.underlying.cyclotomatic_number(&g.full_filter());
            if lc < min_loop_count {
                min_loop_count = lc;
            }
            if lc > max_loop_count {
                max_loop_count = lc;
            }
        }

        let loop_count_range = (min_loop_count, max_loop_count);

        Ok(Self {
            generation_type,
            initial_pdgs,
            final_pdgs_lists,
            loop_count_range,
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
    pub fn warm_up(&mut self, model: &Model) -> Result<()> {
        self.collection.warm_up(model)
    }
    pub fn preprocess(
        &mut self,
        model: &Model,
        settings: &GlobalSettings,
        locked_runtime_settings: &LockedRuntimeSettings,
        thread_pool: &ThreadPool,
    ) -> Result<()> {
        self.collection.preprocess(
            model,
            &self.definition,
            &settings.generation,
            locked_runtime_settings,
            thread_pool,
        )?;
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

        let settings_history = load_settings_history(path.as_ref())?;

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
            let amp = Amplitude::load(path, context).context("Error loading amplitude")?;

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

        let mut collection = ProcessCollection::new_cross_section();
        let settings_history = load_settings_history(path.as_ref())?;

        for entry in fs::read_dir(path.as_ref())
            .context(format!("error reading {}", path.as_ref().display()))?
        {
            let Ok(entry) = entry else {
                debug!("Error reading entry");
                continue;
            };
            if entry.file_type()?.is_file() {
                continue; //skip def.bin, cross sections are in folders
            }
            let path = entry.path();
            debug!("loading cross section at {}", path.display());
            let cs = CrossSection::load(path, context).context("Error loading cross section")?;

            collection.add_cross_section(cs);
        }

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
                    File::create(p.join(SETTINGS_HISTORY_TOML))?
                        .write_all(toml::to_string_pretty(a)?.as_bytes())?;
                }

                for amp in a.values_mut() {
                    amp.save(&p, override_existing)?;
                }
            }
            ProcessCollection::CrossSections(cs) => {
                let p = path.as_ref().join("cross_sections");
                fs::create_dir_all(&p)?;
                let p = p.join(PathBuf::from(self.definition.folder_name.clone()));

                let r = fs::create_dir_all(&p).with_context(|| {
                    format!(
                        "Trying to create directory to save cross section dot {}",
                        p.display()
                    )
                });

                if override_existing {
                    r?;
                }

                let binary = bincode::encode_to_vec(&self.definition, bincode::config::standard())?;
                fs::write(p.join("def.bin"), binary)?;

                if let Some(a) = &self.settings_history {
                    File::create(p.join(SETTINGS_HISTORY_TOML))?
                        .write_all(toml::to_string_pretty(a)?.as_bytes())?;
                }

                for cs in cs.values_mut() {
                    cs.save(&p, override_existing)?;
                }
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

                for amp in a.values_mut() {
                    if let Some(int_name) = integrand_name.clone()
                        && amp.name != int_name
                    {
                        continue;
                    }

                    amp.compile(&p, override_existing, settings, thread_pool)?;
                }
            }
            ProcessCollection::CrossSections(cs) => {
                let p = path.as_ref().join("cross_sections");
                fs::create_dir_all(&p)?;
                let p = p.join(PathBuf::from(self.definition.folder_name.clone()));

                let r = fs::create_dir_all(&p).with_context(|| {
                    format!(
                        "Trying to create directory to export cross section dot {}",
                        p.display()
                    )
                });
                if override_existing {
                    r?;
                }

                for cs in cs.values_mut() {
                    if let Some(int_name) = integrand_name.clone()
                        && cs.name != int_name
                    {
                        continue;
                    }

                    cs.compile(&p, override_existing, settings, thread_pool)?;
                }
            }
        }

        Ok(())
    }

    pub fn get_integrand(
        &self,
        integrand_name: impl AsRef<str>,
    ) -> Result<ResolvedIntegrandRef<'_>> {
        self.collection.get_integrand(integrand_name)
    }

    pub fn get_integrand_names(&self) -> Vec<&str> {
        self.collection.get_integrand_names()
    }

    pub fn get_integrand_mut(
        &mut self,
        integrand_name: impl AsRef<str>,
    ) -> Result<&mut ProcessIntegrand> {
        self.collection.get_integrand_mut(integrand_name)
    }

    pub(crate) fn export_standalone(
        &self,
        path: impl AsRef<Path>,
        settings: &StandaloneExportSettings,
    ) -> Result<()> {
        match &self.collection {
            ProcessCollection::Amplitudes(a) => {
                let p = path.as_ref().join("amplitudes");
                let path = p.join(PathBuf::from(self.definition.folder_name.clone()));
                fs::create_dir_all(&path)?;
                for amp in a.values() {
                    // Create a folder for each amplitude
                    let amp_path = path.join(&amp.name);
                    fs::create_dir_all(&amp_path).with_context(|| {
                        format!(
                            "Trying to create directory for amplitude {}",
                            amp_path.display()
                        )
                    })?;

                    amp.export_standalone(&amp_path, settings)?;
                }
            }
            ProcessCollection::CrossSections(cs) => {
                let p = path.as_ref().join("cross_sections");
                let path = p.join(PathBuf::from(self.definition.folder_name.clone()));
                fs::create_dir_all(&path)?;
                for (xs_name, cs) in cs {
                    // Create a folder for each cross section
                    let cs_path = path.join(&cs.name);
                    fs::create_dir_all(&cs_path).with_context(|| {
                        format!(
                            "Trying to create directory for cross section {}",
                            cs_path.display()
                        )
                    })?;

                    let _ = (xs_name, cs);
                    warn!(
                        "Standalone evaluator export is currently implemented for amplitudes only; skipping cross section output in {}",
                        cs_path.display()
                    );
                }
            }
        }
        Ok(())
    }

    pub(crate) fn export_dot(
        &self,
        path: impl AsRef<Path>,
        settings: &DotExportSettings,
    ) -> Result<()> {
        match &self.collection {
            ProcessCollection::Amplitudes(a) => {
                let p = path.as_ref().join("amplitudes");
                let path = p.join(PathBuf::from(self.definition.folder_name.clone()));
                fs::create_dir_all(&path)?;
                for (amp_name, amp) in a {
                    // Create a folder for each amplitude
                    let amp_path = path.join(&amp.name);
                    fs::create_dir_all(&amp_path).with_context(|| {
                        format!(
                            "Trying to create directory for amplitude {}",
                            amp_path.display()
                        )
                    })?;

                    if settings.combine_diagrams {
                        // Save all graphs combined in one file
                        let output_path = amp_path.join(format!("{}_graphs.dot", amp_name.clone()));
                        let mut dot = create_overwriting_file(&output_path, "amplitude graph")?;
                        for graph in amp.graphs.iter() {
                            graph.graph.dot_serialize_io(&mut dot, settings)?;
                        }
                    } else {
                        // Save each graph in its own file
                        for graph in amp.graphs.iter() {
                            let output_path = amp_path.join(format!("{}.dot", graph.graph.name));
                            let mut dot = create_overwriting_file(&output_path, "amplitude graph")?;
                            graph.graph.dot_serialize_io(&mut dot, settings)?;
                        }
                    }
                }
            }
            ProcessCollection::CrossSections(cs) => {
                let p = path.as_ref().join("cross_sections");
                let path = p.join(PathBuf::from(self.definition.folder_name.clone()));
                fs::create_dir_all(&path)?;
                for (xs_name, cs) in cs {
                    // Create a folder for each cross section
                    let cs_path = path.join(&cs.name);
                    fs::create_dir_all(&cs_path).with_context(|| {
                        format!(
                            "Trying to create directory for cross section {}",
                            cs_path.display()
                        )
                    })?;

                    if settings.combine_diagrams {
                        // Save all graphs combined in one file
                        let output_path = cs_path.join(format!("{}_graphs.dot", xs_name.clone()));
                        let mut dot = create_overwriting_file(&output_path, "cross section graph")?;
                        for graph in cs.supergraphs.iter() {
                            graph.graph.dot_serialize_io(&mut dot, settings)?;
                        }
                    } else {
                        // Save each supergraph in its own file
                        for graph in cs.supergraphs.iter() {
                            let output_path = cs_path.join(format!("{}.dot", graph.graph.name));
                            let mut dot =
                                create_overwriting_file(&output_path, "cross section graph")?;
                            graph.graph.dot_serialize_io(&mut dot, settings)?;
                        }
                    }
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
        model: &Model,
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
                    collection.add_cross_section(CrossSection::from_graph_list(
                        integrand_name,
                        graphs,
                        model,
                    )?);
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

    fn get_integrand(&self, name: impl AsRef<str>) -> Result<ResolvedIntegrandRef<'_>> {
        let canonical_name = self.find_integrand(Some(name.as_ref().to_string()))?;
        let integrand = match self {
            Self::Amplitudes(amplitudes) => amplitudes
                .get(&canonical_name)
                .expect("resolved amplitude name must exist")
                .integrand
                .as_ref(),
            Self::CrossSections(cross_sections) => cross_sections
                .get(&canonical_name)
                .expect("resolved cross section name must exist")
                .integrand
                .as_ref(),
        };

        Ok(ResolvedIntegrandRef {
            canonical_name,
            integrand,
        })
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

    fn get_integrand_mut(&mut self, name: impl AsRef<str>) -> Result<&mut ProcessIntegrand> {
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
    pub fn remove_integrand(&mut self, integrand_name: &str) -> Result<()> {
        match self {
            Self::Amplitudes(amplitudes) => {
                amplitudes
                    .remove(integrand_name)
                    .ok_or(eyre!("No amplitude named {}", integrand_name))?;
            }
            Self::CrossSections(cross_sections) => {
                cross_sections
                    .remove(integrand_name)
                    .ok_or(eyre!("No cross section named {}", integrand_name))?;
            }
        }
        Ok(())
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
        locked_runtime_settings: &LockedRuntimeSettings,
        thread_pool: &ThreadPool,
    ) -> Result<()> {
        match self {
            Self::Amplitudes(amplitudes) => {
                for amplitude in amplitudes.values_mut() {
                    amplitude.preprocess(model, settings, locked_runtime_settings, thread_pool)?;
                }
            }
            Self::CrossSections(cross_sections) => {
                for cross_section in cross_sections.values_mut() {
                    cross_section.preprocess(
                        model,
                        process_definition,
                        settings,
                        *locked_runtime_settings,
                        thread_pool,
                    )?;
                }
            }
        }
        Ok(())
    }

    pub fn warm_up(&mut self, model: &Model) -> Result<()> {
        match self {
            Self::Amplitudes(amplitudes) => {
                for amplitude in amplitudes.values_mut() {
                    amplitude.warm_up(model)?;
                }
            }
            Self::CrossSections(cross_sections) => {
                for cross_section in cross_sections.values_mut() {
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
                for amplitude in amplitudes.values_mut() {
                    amplitude.build_integrand(
                        model,
                        global_settings,
                        runtime_default,
                        thread_pool,
                    )?;
                }
            }
            Self::CrossSections(cross_sections) => {
                for cross_section in cross_sections.values_mut() {
                    cross_section.build_integrand(
                        model,
                        global_settings,
                        runtime_default,
                        thread_pool,
                    )?;
                }
            }
        }
        // result
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::{GammaLoopContextContainer, utils::load_generic_model};

    #[test]
    fn test_proc_definition_encode() {
        let def = super::ProcessDefinition::default();
        let encoded = bincode::encode_to_vec(&def, bincode::config::standard()).unwrap();
        let model_sm = load_generic_model("sm");

        let mut state_file = std::fs::File::create("state_map.bin").unwrap();
        symbolica::state::State::export(&mut state_file).unwrap();
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
