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
    uv::{
        approx::OrientationProjection,
        export::{UVForestExportSettings, sanitize_file_component},
    },
};
use eyre::{Context, eyre};

use crate::{
    feyngen::{FeynGenFilters, GenerationType},
    graph::Graph,
    model::Model,
    settings::global::GenerationSettings,
};

use super::{
    Amplitude, CrossSection, GeneratedGraphReport, GenerationProcessKind, GenerationProgressPhase,
    NamedGraphGenerationReport, generation_progress,
};

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

fn saved_child_dirs(root: &Path, expected_binary: &str, kind: &str) -> Result<Vec<PathBuf>> {
    let mut saved_dirs = Vec::new();

    for entry in fs::read_dir(root).with_context(|| format!("Error reading {}", root.display()))? {
        let Ok(entry) = entry else {
            debug!("Error reading entry");
            continue;
        };
        if !entry.file_type()?.is_dir() {
            continue;
        }

        let path = entry.path();
        if !path.join(expected_binary).is_file() {
            debug!(
                "Skipping helper directory {} while loading {}s because '{}' is missing",
                path.display(),
                kind,
                expected_binary
            );
            continue;
        }

        saved_dirs.push(path);
    }

    saved_dirs.sort();
    Ok(saved_dirs)
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
    pub(crate) fn covariant_cut_states(&self, model: &Model) -> Result<Vec<Vec<i64>>> {
        if self.generation_type != GenerationType::CrossSection {
            return Ok(self.final_pdgs_lists.clone());
        }
        let states = model.covariant_cut_states(&self.final_pdgs_lists)?;
        let representatives = self.covariant_cut_representatives(model);
        for &member in self.final_pdgs_lists.iter().flatten() {
            if let Some(&physical) = representatives.get(&(member as isize))
                && member != physical as i64
            {
                return Err(eyre!(
                    "Final-state PDG {member} overlaps the covariant cut multiplet of requested physical PDG {physical}; separate physical-vector and unphysical diagnostic final-state requests so their observable labels remain unambiguous. For imported covariant partner graphs, supply the physical channel with --process-spec; raw graph states do not establish that intent"
                ));
            }
        }
        for filter in [&self.amplitude_filters, &self.cross_section_filters] {
            if let Some(vetoes) = filter.get_particle_vetos() {
                for (&member, &physical) in &representatives {
                    let anti = model
                        .get_particle_from_pdg(member)
                        .get_anti_particle(model)
                        .pdg_code;
                    if vetoes.contains(&(member as i64)) || vetoes.contains(&(anti as i64)) {
                        return Err(eyre!(
                            "Particle veto removes PDG {member} from the covariant cut multiplet of physical PDG {physical}; retain the complete vector/Goldstone/ghost sector for a physical cross section"
                        ));
                    }
                }
            }
        }
        Ok(states)
    }

    pub(crate) fn covariant_cut_representatives(&self, model: &Model) -> BTreeMap<isize, isize> {
        if self.generation_type != GenerationType::CrossSection {
            return BTreeMap::new();
        }
        let (_, unresolved) = self.unresolved_cut_content(model);
        model
            .covariant_cut_multiplets
            .iter()
            .filter(|(physical, _)| {
                self.final_pdgs_lists
                    .iter()
                    .any(|state| state.contains(physical))
                    || unresolved
                        .iter()
                        .any(|particle| particle.pdg_code as i64 == **physical)
            })
            .flat_map(|(&physical, members)| {
                members
                    .iter()
                    .map(move |&member| (member as isize, physical as isize))
            })
            .collect()
    }

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
                        &g.get_initial_state_tree(),
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
    ) -> Result<Vec<GeneratedGraphReport>> {
        let reports = self.collection.preprocess(
            model,
            &self.definition,
            &settings.generation,
            locked_runtime_settings,
            thread_pool,
        )?;
        self.settings_history = Some(settings.clone());
        Ok(self.attach_process_id(reports))
    }

    fn attach_process_id(
        &self,
        reports: Vec<NamedGraphGenerationReport>,
    ) -> Vec<GeneratedGraphReport> {
        reports
            .into_iter()
            .map(|report| GeneratedGraphReport {
                process_id: self.definition.process_id,
                integrand_name: report.integrand_name,
                graph_name: report.graph_name,
                stats: report.stats,
            })
            .collect()
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
        for path in saved_child_dirs(path.as_ref(), "amp.bin", "amplitude")? {
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
        for path in saved_child_dirs(path.as_ref(), "cs.bin", "cross section")? {
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
        integrand_name: Option<String>,
        thread_pool: &ThreadPool,
    ) -> Result<Vec<GeneratedGraphReport>> {
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

                let mut reports = Vec::new();
                for amp in a.values_mut() {
                    if let Some(int_name) = integrand_name.clone()
                        && amp.name != int_name
                    {
                        continue;
                    }

                    reports.extend(amp.compile(&p, override_existing, thread_pool)?);
                }
                Ok(self.attach_process_id(reports))
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

                let mut reports = Vec::new();
                for cs in cs.values_mut() {
                    if let Some(int_name) = integrand_name.clone()
                        && cs.name != int_name
                    {
                        continue;
                    }

                    reports.extend(cs.compile(&p, override_existing, thread_pool)?);
                }
                Ok(self.attach_process_id(reports))
            }
        }
    }

    pub fn activate_loaded_integrand_backends(
        &mut self,
        allow_symjit_fallback: bool,
    ) -> Result<()> {
        match &mut self.collection {
            ProcessCollection::Amplitudes(amplitudes) => {
                for (integrand_name, amplitude) in amplitudes.iter_mut() {
                    if let Some(integrand) = amplitude.integrand.as_mut()
                        && let Some(reason) =
                            integrand.activate_runtime_backends_after_load(allow_symjit_fallback)?
                    {
                        warn!(
                            "Falling back to symjit for integrand '{}' in process #{} ({}) after external compiled evaluator loading failed: {}",
                            integrand_name,
                            self.definition.process_id,
                            self.definition.folder_name,
                            reason
                        );
                    }
                }
            }
            ProcessCollection::CrossSections(cross_sections) => {
                for (integrand_name, cross_section) in cross_sections.iter_mut() {
                    if let Some(integrand) = cross_section.integrand.as_mut()
                        && let Some(reason) =
                            integrand.activate_runtime_backends_after_load(allow_symjit_fallback)?
                    {
                        warn!(
                            "Falling back to symjit for integrand '{}' in process #{} ({}) after external compiled evaluator loading failed: {}",
                            integrand_name,
                            self.definition.process_id,
                            self.definition.folder_name,
                            reason
                        );
                    }
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
                for cs in cs.values() {
                    // Create a folder for each cross section
                    let cs_path = path.join(&cs.name);
                    fs::create_dir_all(&cs_path).with_context(|| {
                        format!(
                            "Trying to create directory for cross section {}",
                            cs_path.display()
                        )
                    })?;

                    cs.export_standalone(&cs_path, settings)?;
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

    pub(crate) fn export_uv_forests(
        &self,
        path: impl AsRef<Path>,
        integrand_name: &str,
        graph_ids: &[usize],
        settings: &UVForestExportSettings,
    ) -> Result<()> {
        let generation_settings = &self
            .settings_history
            .as_ref()
            .ok_or_else(|| {
                eyre!(
                    "Cannot export UV forests for process {} without generation settings history",
                    self.definition.folder_name
                )
            })?
            .generation;
        let resolved = self.get_integrand(integrand_name)?;
        let integrand = resolved.require_generated()?;
        let integrand_path = match &self.collection {
            ProcessCollection::Amplitudes(_) => path
                .as_ref()
                .join("amplitudes")
                .join(PathBuf::from(self.definition.folder_name.clone()))
                .join(&resolved.canonical_name),
            ProcessCollection::CrossSections(_) => path
                .as_ref()
                .join("cross_sections")
                .join(PathBuf::from(self.definition.folder_name.clone()))
                .join(&resolved.canonical_name),
        };
        fs::create_dir_all(&integrand_path).with_context(|| {
            format!(
                "Trying to create directory for UV forest export {}",
                integrand_path.display()
            )
        })?;

        for &graph_id in graph_ids {
            let source = if settings.computed {
                let (graph, expression) = match &self.collection {
                    ProcessCollection::Amplitudes(amplitudes) => {
                        let source = amplitudes[&resolved.canonical_name]
                            .graphs
                            .get(graph_id)
                            .ok_or_else(|| eyre!("Missing source amplitude graph {graph_id}"))?;
                        (&source.graph, source.derived_data.cff_expression.as_ref())
                    }
                    ProcessCollection::CrossSections(cross_sections) => {
                        let source = cross_sections[&resolved.canonical_name]
                            .supergraphs
                            .get(graph_id)
                            .ok_or_else(|| {
                                eyre!("Missing source cross-section graph {graph_id}")
                            })?;
                        (
                            &source.graph,
                            source.derived_data.global_cff_expression.as_ref(),
                        )
                    }
                };
                // Generation and persistent selection preserve graph order;
                // reject stale runtime metadata instead of pairing a different
                // graph with this stored production residue map.
                if integrand.graph_name_by_id(graph_id) != Some(graph.name.as_str()) {
                    return Err(eyre!(
                        "Source/runtime graph mismatch for computed UV forest export at id {graph_id}"
                    ));
                }
                Some((
                    graph,
                    expression.ok_or_else(|| {
                        eyre!(
                            "Graph {} has no stored production CFF for computed UV forest export",
                            graph.name
                        )
                    })?,
                ))
            } else {
                None
            };
            let cff_options = source
                .map(|(graph, _)| graph.production_cff_3d_expression_options(generation_settings))
                .transpose()?;
            let orientation = source
                .zip(cff_options.as_ref())
                .map(|((_, expression), options)| {
                    OrientationProjection::exact_expression(
                        expression,
                        options,
                        &generation_settings.orientation_pattern,
                        generation_settings.explicit_orientation_sum_only,
                    )
                });
            let export = integrand.export_uv_forest_graph(
                graph_id,
                orientation,
                generation_settings,
                settings,
            )?;
            let graph_name = sanitize_file_component(&export.graph_name);
            let forest_path = integrand_path.join(format!("{graph_name}.forest.dot"));
            let mut forest_file = create_overwriting_file(&forest_path, "UV forest")?;
            forest_file.write_all(export.forest_dot.as_bytes())?;

            for term in export.node_terms {
                let node_dir = integrand_path
                    .join(format!("{graph_name}_nodes"))
                    .join(format!("forest_{:03}", term.forest_index));
                fs::create_dir_all(&node_dir).with_context(|| {
                    format!(
                        "Trying to create directory for UV forest node graph {}",
                        node_dir.display()
                    )
                })?;
                let mut dot = create_overwriting_file(
                    &node_dir.join(term.file_name()),
                    "UV forest node graph",
                )?;
                dot.write_all(term.dot.as_bytes())?;
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
    ) -> Result<Vec<GeneratedGraphReport>> {
        let reports = self.collection.generate_integrands(
            model,
            &self.definition,
            global_settings,
            runtime_default,
            thread_pool,
        )?;
        Ok(self.attach_process_id(reports))
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
    ) -> Result<Vec<NamedGraphGenerationReport>> {
        match self {
            Self::Amplitudes(amplitudes) => {
                let mut reports = Vec::new();
                for amplitude in amplitudes.values_mut() {
                    generation_progress::begin_phase(
                        GenerationProgressPhase::GraphPreprocessing,
                        GenerationProcessKind::Amplitude,
                        &process_definition.folder_name,
                        &amplitude.name,
                        amplitude.graphs.len(),
                        None,
                    );
                    reports.extend(amplitude.preprocess(
                        model,
                        settings,
                        locked_runtime_settings,
                        thread_pool,
                    )?);
                }
                Ok(reports)
            }
            Self::CrossSections(cross_sections) => {
                let mut reports = Vec::new();
                for cross_section in cross_sections.values_mut() {
                    reports.extend(cross_section.preprocess(
                        model,
                        process_definition,
                        settings,
                        *locked_runtime_settings,
                        thread_pool,
                    )?);
                }
                Ok(reports)
            }
        }
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
        process_definition: &ProcessDefinition,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
        thread_pool: &ThreadPool,
    ) -> Result<Vec<NamedGraphGenerationReport>> {
        match self {
            Self::Amplitudes(amplitudes) => {
                let mut reports = Vec::new();
                for amplitude in amplitudes.values_mut() {
                    reports.extend(amplitude.build_integrand(
                        model,
                        &process_definition.folder_name,
                        global_settings,
                        runtime_default,
                        thread_pool,
                    )?);
                }
                Ok(reports)
            }
            Self::CrossSections(cross_sections) => {
                let mut reports = Vec::new();
                for cross_section in cross_sections.values_mut() {
                    reports.extend(cross_section.build_integrand(
                        model,
                        process_definition,
                        global_settings,
                        runtime_default,
                        thread_pool,
                    )?);
                }
                Ok(reports)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::{
        fs,
        path::PathBuf,
        time::{SystemTime, UNIX_EPOCH},
    };

    use crate::{GammaLoopContextContainer, utils::load_generic_model};
    use symbolica::atom::{Atom, AtomCore};

    fn fresh_temp_dir(name: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "gammalooprs-{name}-{}-{unique}",
            std::process::id()
        ));
        fs::create_dir_all(&path).unwrap();
        path
    }

    #[test]
    fn saved_child_dirs_skip_cross_section_compile_artifact_folders() {
        let temp = fresh_temp_dir("saved-child-dirs");
        let saved_dir = temp.join("NLO");
        let compiled_dir = temp.join("cs_NLO");

        fs::create_dir_all(&saved_dir).unwrap();
        fs::write(saved_dir.join("cs.bin"), []).unwrap();
        fs::create_dir_all(compiled_dir.join("integrand").join("GL08")).unwrap();

        let dirs = super::saved_child_dirs(&temp, "cs.bin", "cross section").unwrap();

        assert_eq!(dirs, vec![saved_dir]);
        fs::remove_dir_all(temp).unwrap();
    }

    #[test]
    fn computed_uv_forest_process_export_uses_stored_sources_in_all_routes()
    -> color_eyre::Result<()> {
        use crate::{
            feyngen::GenerationType,
            graph::{Graph, GroupId},
            initialisation::test_initialise,
            processes::{
                GraphGroupSelectionSpec, Process, ProcessCollection, ProcessDefinition, ProcessList,
            },
            settings::{GlobalSettings, RuntimeSettings},
            uv::{
                UVOrchestrator, export::UVForestExportSettings, settings::FinalIntegrandDimension,
            },
        };
        use std::collections::{BTreeMap, BTreeSet};

        test_initialise()?;
        let model = load_generic_model("scalars");
        let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build()?;
        let runtime = RuntimeSettings::default();
        let directory = fresh_temp_dir("computed-uv-process-export");
        for (kind, folder, source) in [
            (
                GenerationType::Amplitude,
                "amplitudes",
                include_str!(concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/../../tests/resources/graphs/scalar_bubble.dot"
                )),
            ),
            (
                GenerationType::CrossSection,
                "cross_sections",
                include_str!(concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/../../tests/resources/graphs/mass_approach_scalar_self_energy.dot"
                )),
            ),
        ] {
            for orchestrator in [UVOrchestrator::LegacyDagForest, UVOrchestrator::HedgePoset] {
                for (mode, explicit_sum, projected) in [
                    ("direct_keyed", false, false),
                    ("direct_summed", true, false),
                    ("projected_summed", true, true),
                ] {
                    // Keep the fixture's physical graph and all default subtraction
                    // terms; only the three requested CFF/UV routes differ.
                    let mut graphs = Graph::from_string(source, &model)?;
                    let exercise_selection = kind == GenerationType::CrossSection
                        && orchestrator == UVOrchestrator::HedgePoset
                        && mode == "direct_keyed";
                    if !exercise_selection {
                        graphs.truncate(1);
                    }
                    let mut graph_name = graphs[0].name.clone();
                    let definition = ProcessDefinition::from_graph_list(&graphs, kind, &model)?;
                    let process = Process::from_graph_list(
                        "export_fixture".into(),
                        "default".into(),
                        graphs,
                        kind,
                        Some(definition),
                        None,
                        &model,
                    )?;
                    let mut processes = ProcessList {
                        processes: vec![process],
                    };
                    let mut settings = GlobalSettings::default();
                    settings.generation.uv.orchestrator = orchestrator;
                    settings.generation.explicit_orientation_sum_only = explicit_sum;
                    settings
                        .generation
                        .uv
                        .local_uv_cts_from_expanded_4d_integrands = projected;
                    assert_eq!(
                        settings.generation.uv.final_integrand,
                        FinalIntegrandDimension::ThreeD
                    );
                    assert!(settings.generation.uv.subtract_uv);
                    assert!(settings.generation.uv.generate_integrated);
                    assert!(settings.generation.uv.softct);
                    assert!(settings.generation.threshold_subtraction.enable_thresholds);
                    assert!(!settings.generation.evaluator.store_atom);
                    processes.preprocess(&model, &settings, &(&runtime).into(), &pool)?;
                    processes.generate_integrands(&model, &settings, (&runtime).into(), &pool)?;
                    // Record both physical cut inventories before selection, so
                    // removal/reindexing cannot silently redefine the expected output.
                    let expected_by_graph: Vec<Vec<BTreeSet<String>>> = match processes.processes[0]
                        .get_integrand("default")?
                        .require_generated()?
                    {
                        crate::integrands::process::ProcessIntegrand::Amplitude(_) => {
                            vec![vec![["all_none".to_string()].into_iter().collect()]]
                        }
                        crate::integrands::process::ProcessIntegrand::CrossSection(integrand) => {
                            integrand
                                .data
                                .graph_terms
                                .iter()
                                .map(|term| {
                                    term.cut_group_data
                                        .cut_groups
                                        .iter()
                                        .map(|cuts| {
                                            crate::graph::cuts::ResidueSelector {
                                                lu: Some(
                                                    cuts.lu_cut_selection(&term.graph, &term.cuts),
                                                ),
                                                left_th_cut: None,
                                                right_th_cut: None,
                                            }
                                            .generate_allowed_keys()
                                            .into_iter()
                                            .map(|key| {
                                                format!("lu_cut_{}", key.lu_cut_order.unwrap())
                                            })
                                            .collect()
                                        })
                                        .collect()
                                })
                                .collect()
                        }
                    };
                    assert!(!expected_by_graph.is_empty());
                    let case_dir = directory.join(format!("{folder}-{orchestrator}-{mode}"));
                    let source_expression =
                        |graph: &Graph,
                         source: &three_dimensional_reps::GeneratedThreeDExpression<
                            crate::cff::esurface::Esurface,
                            crate::cff::hsurface::Hsurface,
                        >|
                         -> eyre::Result<Atom> {
                            let expressions = graph
                                .cff_from_production_expression(
                                    source,
                                    &crate::graph::cuts::CutSet::empty(graph.n_hedges()),
                                    &crate::settings::global::OrientationPattern::default(),
                                )?
                                .expression_with_selectors();
                            Ok(expressions
                                .iter()
                                .fold(Atom::Zero, |sum, (_, atom)| sum + atom))
                        };
                    let mut production_expression: Option<Atom> = None;
                    let mut expected_exported_expressions: Option<BTreeMap<(usize, String), Atom>> =
                        None;
                    for phase in ["generated", "loaded", "selected"] {
                        if phase == "selected" {
                            if !exercise_selection {
                                continue;
                            }
                            // Exercise the same persistent selection owner as the CLI,
                            // after save/load, removing graph 0 and compacting graph 1.
                            let ProcessCollection::CrossSections(cross_sections) =
                                &mut processes.processes[0].collection
                            else {
                                unreachable!();
                            };
                            let cross_section = cross_sections.get_mut("default").unwrap();
                            let retained = &cross_section.supergraphs[1];
                            let removed_name = graph_name.clone();
                            graph_name = retained.graph.name.clone();
                            assert_ne!(graph_name, removed_name);
                            let expected_cuts = retained
                                .cuts
                                .iter()
                                .map(|cut| cut.cut.clone())
                                .collect::<BTreeSet<_>>();
                            let expected_numerator = retained
                                .graph
                                .production_numerator_atom_for_full_3d_expression();
                            expected_exported_expressions = None;
                            production_expression = Some(source_expression(
                                &retained.graph,
                                retained
                                    .derived_data
                                    .global_cff_expression
                                    .as_ref()
                                    .unwrap(),
                            )?);
                            let plan = cross_section.plan_graph_group_selection(
                                &GraphGroupSelectionSpec::from_master_graph_names(vec![
                                    graph_name.clone(),
                                ]),
                            )?;
                            assert_eq!(plan.retained_group_ids(), &[GroupId(1)]);
                            assert_eq!(plan.new_group_id_for_old(GroupId(1)), Some(GroupId(0)));
                            assert_eq!(plan.report().removed_graphs, vec![removed_name.clone()]);
                            cross_section.apply_graph_group_selection(&plan)?;
                            let retained = &cross_section.supergraphs[0];
                            assert!(
                                (retained
                                    .graph
                                    .production_numerator_atom_for_full_3d_expression()
                                    - expected_numerator)
                                    .collect_factors()
                                    .is_zero(),
                                "persistent selection changed the retained physical numerator"
                            );
                            assert!(
                                processes.processes[0]
                                    .get_integrand("default")?
                                    .require_generated()
                                    .is_err()
                            );
                            // Match normal regeneration: preprocess the retained source
                            // and then rebuild the runtime integrand with identical settings.
                            processes.preprocess(&model, &settings, &(&runtime).into(), &pool)?;
                            processes.generate_integrands(
                                &model,
                                &settings,
                                (&runtime).into(),
                                &pool,
                            )?;
                            assert_eq!(
                                processes.processes[0]
                                    .get_integrand("default")?
                                    .require_generated()?
                                    .graph_name_by_id(1),
                                None,
                            );
                            let crate::integrands::process::ProcessIntegrand::CrossSection(
                                integrand,
                            ) = processes.processes[0]
                                .get_integrand("default")?
                                .require_generated()?
                            else {
                                unreachable!();
                            };
                            let term = &integrand.data.graph_terms[0];
                            assert_eq!(
                                term.cuts
                                    .iter()
                                    .map(|cut| cut.cut.clone())
                                    .collect::<BTreeSet<_>>(),
                                expected_cuts,
                                "selection/regeneration changed the retained physical cuts"
                            );
                        }
                        let expected_residues =
                            &expected_by_graph[usize::from(phase == "selected")];
                        assert!(!expected_residues.is_empty());
                        let process = &processes.processes[0];
                        let (source_graph, production) = match &process.collection {
                            ProcessCollection::Amplitudes(amplitudes) => {
                                let graph = &amplitudes["default"].graphs[0];
                                (
                                    &graph.graph,
                                    graph.derived_data.cff_expression.as_ref().unwrap(),
                                )
                            }
                            ProcessCollection::CrossSections(cross_sections) => {
                                let graph = &cross_sections["default"].supergraphs[0];
                                (
                                    &graph.graph,
                                    graph.derived_data.global_cff_expression.as_ref().unwrap(),
                                )
                            }
                        };
                        assert_eq!(source_graph.name, graph_name);
                        assert_eq!(
                            process
                                .get_integrand("default")?
                                .require_generated()?
                                .graph_name_by_id(0),
                            Some(graph_name.as_str())
                        );
                        assert!(!production.expression.orientations.is_empty());
                        // Persistence retains the complete consumed expression, including
                        // surface references, energy factors, and convention prefactors;
                        // transient degree reports and internal cache/tree layout do not
                        // affect this contract.
                        // Selection/regeneration must also retain the old graph 1 source.
                        let expression = source_expression(source_graph, production)?;
                        if let Some(expected) = &production_expression {
                            assert!(
                                (&expression - expected).expand().together().is_zero(),
                                "stored production CFF changed after {phase}"
                            );
                        } else {
                            production_expression = Some(expression);
                        }
                        let export_dir = case_dir.join(phase);
                        processes.export_uv_forests(
                            &export_dir,
                            0,
                            "default",
                            &[0],
                            &UVForestExportSettings { computed: true },
                        )?;
                        let graph_dir = export_dir
                            .join("processes")
                            .join(folder)
                            .join("export_fixture/default");
                        assert!(graph_dir.join(format!("{graph_name}.forest.dot")).is_file());
                        assert_eq!(
                            fs::read_dir(&graph_dir)?
                                .map(|entry| entry
                                    .map(|entry| entry.file_name().to_string_lossy().into_owned()))
                                .collect::<std::io::Result<Vec<_>>>()?
                                .into_iter()
                                .filter(|name| name.ends_with(".forest.dot"))
                                .collect::<BTreeSet<_>>(),
                            BTreeSet::from([format!("{graph_name}.forest.dot")]),
                            "computed export retained a removed graph",
                        );
                        let mut exported_expressions = BTreeMap::<(usize, String), Atom>::new();
                        let mut forest_indices = BTreeSet::new();
                        for forest in fs::read_dir(graph_dir.join(format!("{graph_name}_nodes")))? {
                            let forest = forest?;
                            let forest_index = forest
                                .file_name()
                                .to_string_lossy()
                                .strip_prefix("forest_")
                                .unwrap()
                                .parse::<usize>()?;
                            assert!(forest_indices.insert(forest_index));
                            let expected = &expected_residues[forest_index];
                            let mut actual = BTreeSet::new();
                            for node in fs::read_dir(forest.path())? {
                                let node = node?;
                                let file = node.file_name();
                                let file = file.to_string_lossy();
                                let residue = expected
                                    .iter()
                                    .find(|residue| file.ends_with(&format!("_{residue}.dot")))
                                    .expect(
                                        "exported residue must belong to this complete physical cut",
                                    );
                                actual.insert(residue.clone());
                                let node_dot = fs::read_to_string(node.path())?;
                                assert!(node_dot.contains("forest_residue_index"));
                                for exported in Graph::from_string(&node_dot, &model)? {
                                    *exported_expressions
                                        .entry((forest_index, residue.clone()))
                                        .or_insert(Atom::Zero) += exported.global_prefactor.num;
                                }
                            }
                            assert_eq!(&actual, expected, "a physical cut lost a residue order");
                        }
                        assert_eq!(
                            forest_indices,
                            (0..expected_residues.len()).collect::<BTreeSet<_>>()
                        );
                        assert!(
                            exported_expressions
                                .values()
                                .any(|expression| !expression.is_zero()),
                            "{folder}/{orchestrator}/{mode}: the computed export has no nonzero value"
                        );
                        // Persistence preserves complete exported values for every
                        // physical cut and residue, independently of how terms are
                        // grouped into nodes or numbered in the forest traversal.
                        if let Some(expected) = &expected_exported_expressions {
                            assert_eq!(
                                exported_expressions.keys().collect::<BTreeSet<_>>(),
                                expected.keys().collect::<BTreeSet<_>>()
                            );
                            for (key, expression) in &exported_expressions {
                                assert!(
                                    (expression.collect_factors()
                                        - expected[key].collect_factors())
                                    .collect_factors()
                                    .is_zero(),
                                    "{folder}/{orchestrator}/{mode}: exported residue {key:?} changed after {phase}"
                                );
                            }
                        } else {
                            expected_exported_expressions = Some(exported_expressions);
                        }
                        if phase == "generated" {
                            let saved = case_dir.join("saved");
                            processes.processes[0].save(&saved, true)?;
                            let mut symbols = Vec::new();
                            symbolica::state::State::export(&mut symbols)?;
                            let state_map = symbolica::state::State::import(
                                &mut std::io::Cursor::new(symbols),
                                None,
                            )?;
                            let context = GammaLoopContextContainer {
                                model: &model,
                                state_map: &state_map,
                            };
                            let path = saved.join(folder).join("export_fixture");
                            processes.processes[0] = match kind {
                                GenerationType::Amplitude => {
                                    Process::load_amplitude(path, context)?
                                }
                                GenerationType::CrossSection => {
                                    Process::load_cross_section(path, context)?
                                }
                            };
                        }
                    }
                    // The public lookup must reject an ID/source mismatch rather
                    // than pairing a runtime graph with a different stored CFF.
                    let source_graph = match &mut processes.processes[0].collection {
                        ProcessCollection::Amplitudes(amplitudes) => {
                            &mut amplitudes.get_mut("default").unwrap().graphs[0].graph
                        }
                        ProcessCollection::CrossSections(cross_sections) => {
                            &mut cross_sections.get_mut("default").unwrap().supergraphs[0].graph
                        }
                    };
                    source_graph.name.push_str("_mismatch");
                    let error = processes
                        .export_uv_forests(
                            case_dir.join("mismatch"),
                            0,
                            "default",
                            &[0],
                            &UVForestExportSettings { computed: true },
                        )
                        .unwrap_err();
                    assert!(format!("{error:#}").contains("Source/runtime graph mismatch"));
                    let source = match &mut processes.processes[0].collection {
                        ProcessCollection::Amplitudes(amplitudes) => {
                            let graph = &mut amplitudes.get_mut("default").unwrap().graphs[0];
                            graph.graph.name = graph_name.clone();
                            &mut graph.derived_data.cff_expression
                        }
                        ProcessCollection::CrossSections(cross_sections) => {
                            let graph =
                                &mut cross_sections.get_mut("default").unwrap().supergraphs[0];
                            graph.graph.name = graph_name.clone();
                            &mut graph.derived_data.global_cff_expression
                        }
                    };
                    *source = None;
                    let topology_dir = case_dir.join("topology_without_source");
                    processes.export_uv_forests(
                        &topology_dir,
                        0,
                        "default",
                        &[0],
                        &UVForestExportSettings { computed: false },
                    )?;
                    let topology_graph_dir = topology_dir
                        .join("processes")
                        .join(folder)
                        .join("export_fixture/default");
                    assert!(
                        topology_graph_dir
                            .join(format!("{graph_name}.forest.dot"))
                            .is_file()
                    );
                    assert!(
                        !topology_graph_dir
                            .join(format!("{graph_name}_nodes"))
                            .exists()
                    );
                    let error = processes
                        .export_uv_forests(
                            case_dir.join("missing_source"),
                            0,
                            "default",
                            &[0],
                            &UVForestExportSettings { computed: true },
                        )
                        .unwrap_err();
                    assert!(format!("{error:#}").contains("has no stored production CFF"));
                }
            }
        }
        fs::remove_dir_all(directory)?;
        Ok(())
    }

    mod failing {
        use super::*;

        #[test]
        fn test_proc_definition_encode() {
            let def = crate::processes::ProcessDefinition::default();
            let encoded = bincode::encode_to_vec(&def, bincode::config::standard()).unwrap();
            let model_sm = load_generic_model("sm");

            let mut state_file = std::fs::File::create("state_map.bin").unwrap();
            symbolica::state::State::export(&mut state_file).unwrap();
            let state_map = symbolica::state::State::import(&mut state_file, None).unwrap();

            let context = GammaLoopContextContainer {
                model: &model_sm,
                state_map: &state_map,
            };

            let (decoded, _): (crate::processes::ProcessDefinition, _) =
                bincode::decode_from_slice_with_context(
                    &encoded,
                    bincode::config::standard(),
                    context,
                )
                .unwrap();
            assert_eq!(def, decoded);
        }
    }
}
