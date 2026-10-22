use ahash::HashSet;
use linnet::half_edge::involution::Flow;
use linnet::half_edge::involution::HedgePair;
use linnet::half_edge::involution::Orientation;
use rayon::ThreadPool;
use std::result::Result as SnapshotResult;
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File},
    io::Write,
    path::Path,
};
use tracing::warn;
// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::{Help, Result};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::fmt;
use tracing::debug;

use crate::graph::edge::PossibleParticle;
use crate::processes::DotExportSettings;

use crate::processes::StandaloneExportSettings;
use crate::{
    GammaLoopContext, GammaLoopContextContainer,
    integrands::process::ProcessIntegrand,
    settings::{GlobalSettings, RuntimeSettings, runtime::LockedRuntimeSettings},
    uv::export::UVForestExportSettings,
};
use eyre::{Context, eyre};
use feynkit_generator::{
    FilterScope, GenerationOptions, GenerationType, Process as GenerationProcess,
};

use crate::{graph::Graph, model::Model, settings::global::GenerationSettings};

use super::{
    Amplitude, CrossSection, GeneratedGraphReport, GenerationProcessKind, GenerationProgressPhase,
    IntegrandSnapshotEntry, NamedGraphGenerationReport, OutputPath, ProcessSnapshotEntry,
    ProcessSnapshotError, generation_progress,
};

const SETTINGS_HISTORY_TOML: &str = "settings_history.toml";

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
    OutputPath::new(path)?.create_file().with_context(|| {
        format!(
            "Trying to create file to export {file_kind} {}",
            path.display()
        )
    })
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, bincode::Encode, bincode::Decode)]
pub struct ProcessDefinition {
    pub generation_type: GenerationType,
    pub process: GenerationProcess,
    pub generation_options: GenerationOptions,
    pub folder_name: String,
    pub process_id: usize,
}

impl fmt::Display for ProcessDefinition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Process #{}: '{}'\nGeneration type: {}{}{}\nInitial PDGs: {:?}{}\nFinal PDGs: {}{}\nLoop count: {}\nAmplitude filters:{}{}\nCross-section filters:{}{}",
            self.process_id,
            self.folder_name,
            self.generation_type,
            if self.generation_options.symmetrizes_left_right() {
                " (left-right symmetrized)"
            } else {
                ""
            },
            if self.generation_options.symmetrizes_external_fermions()
                && self.generation_type == GenerationType::Amplitude
                && (self.generation_options.symmetrizes_initial()
                    || self.generation_options.symmetrizes_final()
                    || self.generation_options.symmetrizes_left_right())
            {
                " (allowing fermion symmetrization)"
            } else {
                ""
            },
            self.process.incoming(),
            if self.generation_options.symmetrizes_initial() {
                " (symmetrized)"
            } else {
                ""
            },
            if self.process.outgoing_alternatives().len() == 1 {
                format!("{:?}", self.process.outgoing_alternatives()[0])
            } else {
                format!(
                    "[ {} ]",
                    self.process
                        .outgoing_alternatives()
                        .iter()
                        .map(|pdgs| format!("{:?}", pdgs))
                        .join(" | ")
                )
            },
            if self.generation_options.symmetrizes_final() {
                " (symmetrized)"
            } else {
                ""
            },
            if self.generation_options.loop_count().start()
                == self.generation_options.loop_count().end()
            {
                format!("{}", self.generation_options.loop_count().start())
            } else {
                format!("{:?}", self.generation_options.loop_count())
            },
            if self
                .generation_options
                .filters(FilterScope::CutAmplitude)
                .is_empty()
            {
                " None"
            } else {
                "\n"
            },
            if self
                .generation_options
                .filters(FilterScope::CutAmplitude)
                .is_empty()
            {
                "".into()
            } else {
                self.generation_options
                    .filters(FilterScope::CutAmplitude)
                    .iter()
                    .map(|f| format!(" > {}", f))
                    .collect::<Vec<String>>()
                    .join("\n")
            },
            if self
                .generation_options
                .filters(FilterScope::Graph)
                .is_empty()
            {
                " None"
            } else {
                "\n"
            },
            if self
                .generation_options
                .filters(FilterScope::Graph)
                .is_empty()
            {
                "".into()
            } else {
                self.generation_options
                    .filters(FilterScope::Graph)
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
            process: GenerationProcess::new(Vec::<i64>::new(), Vec::<i64>::new()),
            generation_options: GenerationOptions::default()
                .graph_prefix("GL")
                .with_loop_count(1, 1)
                .expect("the default loop range is valid"),
            folder_name: "undefined_process".to_string(),
            process_id: 0,
        }
    }
}

impl ProcessDefinition {
    pub(crate) fn covariant_cut_states(&self, model: &Model) -> Result<Vec<Vec<i64>>> {
        self.process.validate_covariant_cut_filters(
            model,
            &self.generation_options,
            self.generation_type,
        )?;
        Ok(self
            .process
            .covariant_cut_states(model, self.generation_type)?)
    }

    pub(crate) fn covariant_cut_representatives(&self, model: &Model) -> BTreeMap<isize, isize> {
        self.process
            .covariant_cut_representatives(model, self.generation_type)
            .expect("process selectors were validated against the graph model")
            .into_iter()
            .map(|(member, physical)| (member as isize, physical as isize))
            .collect()
    }

    pub fn may_filter_covariant_partners(&self, model: &Model) -> bool {
        !self.covariant_cut_representatives(model).is_empty()
            && self.generation_options.may_filter_covariant_partners()
    }

    // Best attempt at creating what process definition matches the given graphs
    pub fn from_graph_list(
        graphs: &[Graph],
        generation_type: GenerationType,
        model: &Model,
    ) -> Result<Self> {
        Self::infer_graph_process(
            graphs
                .iter()
                .map(|graph| (graph, graph.finalized_cuts.as_slice())),
            generation_type,
            model,
        )
    }

    /// Infer from the runtime cut owner after finalized inventories are consumed.
    #[cfg(test)]
    pub(crate) fn from_cross_section_graphs(
        graphs: &[super::CrossSectionGraph],
        model: &Model,
    ) -> Result<Self> {
        Self::infer_graph_process(
            graphs
                .iter()
                .map(|graph| (&graph.graph, graph.cuts.raw.as_slice())),
            GenerationType::CrossSection,
            model,
        )
    }

    fn infer_graph_process<'a>(
        graphs: impl IntoIterator<Item = (&'a Graph, &'a [crate::graph::FinalizedCut])>,
        generation_type: GenerationType,
        model: &Model,
    ) -> Result<Self> {
        let graphs = graphs.into_iter().collect::<Vec<_>>();
        let mut initial_pdgs = HashSet::default();

        for &(g, _) in &graphs {
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
                                initial_pdgs_of_graph
                                    .push(model.particle_by_id(*particle)?.pdg_code);
                            } else {
                                debug!("Edge without particle data in initial state");
                            }
                        }
                    }
                }
                GenerationType::CrossSection => {
                    for (_, _, edge) in g.iter_edges_of(&g.initial_state_cut) {
                        if let PossibleParticle::Particle(particle) = &edge.data.particle {
                            initial_pdgs_of_graph.push(model.particle_by_id(*particle)?.pdg_code);
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
                for &(g, _) in &graphs {
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
                                final_pdgs_of_graph.push(model.particle_by_id(*particle)?.pdg_code);
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
                for &(g, cuts) in &graphs {
                    for finalized in cuts {
                        let mut final_pdgs_of_cut = vec![];
                        for (orientaion, edge) in finalized.cut.iter_edges(&g.underlying) {
                            if let PossibleParticle::Particle(particle) = &edge.data.particle {
                                if orientaion == Orientation::Reversed {
                                    final_pdgs_of_cut.push(
                                        model
                                            .particle_by_id(
                                                model.particle_by_id(*particle)?.antiparticle,
                                            )?
                                            .pdg_code,
                                    );
                                } else {
                                    final_pdgs_of_cut
                                        .push(model.particle_by_id(*particle)?.pdg_code);
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

        for &(g, _) in &graphs {
            // don't know how the looop count is really intended, for now it doesn't matter I think
            let lc = g.underlying.cyclotomatic_number(&g.full_filter());
            if lc < min_loop_count {
                min_loop_count = lc;
            }
            if lc > max_loop_count {
                max_loop_count = lc;
            }
        }

        let first_final = final_pdgs_lists.first().cloned().unwrap_or_default();
        let process = match generation_type {
            GenerationType::Amplitude => GenerationProcess::new(initial_pdgs, first_final),
            GenerationType::CrossSection => GenerationProcess::new(initial_pdgs, first_final)
                .with_final_state_alternatives(final_pdgs_lists)
                .map_err(|error| eyre!(error))?,
        };

        Ok(Self {
            generation_type,
            process,
            generation_options: Self::default()
                .generation_options
                .with_loop_count(min_loop_count, max_loop_count)
                .map_err(|error| eyre!(error))?,
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
    fn load_settings_history(path: &Path) -> SnapshotResult<GlobalSettings, ProcessSnapshotError> {
        let path = path.join(SETTINGS_HISTORY_TOML);
        let source = fs::read_to_string(&path).map_err(|source| ProcessSnapshotError::Read {
            path: path.clone(),
            source,
        })?;
        toml::from_str(&source)
            .map_err(|source| ProcessSnapshotError::SettingsDecode { path, source })
    }

    pub(super) fn snapshot_entry(
        &self,
    ) -> SnapshotResult<ProcessSnapshotEntry, ProcessSnapshotError> {
        let mut children = Vec::new();
        let kind = match &self.collection {
            ProcessCollection::Amplitudes(amplitudes) => {
                for (name, amplitude) in amplitudes {
                    if name != &amplitude.name {
                        return Err(ProcessSnapshotError::ChildNameMismatch {
                            process: self.definition.folder_name.clone(),
                            expected: name.clone(),
                            actual: amplitude.name.clone(),
                        });
                    }
                    if amplitude.integrand.as_ref().is_some_and(|integrand| !matches!(integrand, ProcessIntegrand::Amplitude(value) if value.data.name == *name)) {
                        return Err(ProcessSnapshotError::GeneratedChildMismatch { process: self.definition.folder_name.clone(), child: name.clone(), kind: GenerationType::Amplitude });
                    }
                    children.push(IntegrandSnapshotEntry {
                        name: name.clone(),
                        generated: amplitude.integrand.is_some(),
                    });
                }
                GenerationType::Amplitude
            }
            ProcessCollection::CrossSections(cross_sections) => {
                for (name, cross_section) in cross_sections {
                    if name != &cross_section.name {
                        return Err(ProcessSnapshotError::ChildNameMismatch {
                            process: self.definition.folder_name.clone(),
                            expected: name.clone(),
                            actual: cross_section.name.clone(),
                        });
                    }
                    if cross_section.integrand.as_ref().is_some_and(|integrand| !matches!(integrand, ProcessIntegrand::CrossSection(value) if value.data.name == *name)) {
                        return Err(ProcessSnapshotError::GeneratedChildMismatch { process: self.definition.folder_name.clone(), child: name.clone(), kind: GenerationType::CrossSection });
                    }
                    children.push(IntegrandSnapshotEntry {
                        name: name.clone(),
                        generated: cross_section.integrand.is_some(),
                    });
                }
                GenerationType::CrossSection
            }
        };
        if kind != self.definition.generation_type {
            return Err(ProcessSnapshotError::ProcessKindMismatch {
                process: self.definition.folder_name.clone(),
                expected: kind,
                actual: self.definition.generation_type,
            });
        }
        Ok(ProcessSnapshotEntry {
            kind,
            folder_name: self.definition.folder_name.clone(),
            children,
            settings_history: self.settings_history.is_some(),
        })
    }

    pub(super) fn load_snapshot(
        root: &Path,
        entry: &ProcessSnapshotEntry,
        process_id: usize,
        context: GammaLoopContextContainer<'_>,
        selected: Option<&BTreeSet<String>>,
    ) -> SnapshotResult<Self, ProcessSnapshotError> {
        let path = OutputPath::new(root)?
            .join("processes")?
            .join(entry.kind_directory())?
            .join(&entry.folder_name)?;
        let path = path.as_path();
        let definition_path = path.join("def.bin");
        let binary = fs::read(&definition_path).map_err(|source| ProcessSnapshotError::Read {
            path: definition_path.clone(),
            source,
        })?;
        let (definition, _): (ProcessDefinition, _) =
            bincode::decode_from_slice_with_context(&binary, bincode::config::standard(), context)
                .map_err(|source| ProcessSnapshotError::DefinitionDecode {
                    path: definition_path,
                    source,
                })?;
        if definition.process_id != process_id {
            return Err(ProcessSnapshotError::InvalidProcessId {
                process: entry.folder_name.clone(),
                expected: process_id,
                actual: definition.process_id,
            });
        }
        if definition.folder_name != entry.folder_name {
            return Err(ProcessSnapshotError::ProcessNameMismatch {
                expected: entry.folder_name.clone(),
                actual: definition.folder_name,
            });
        }
        if definition.generation_type != entry.kind {
            return Err(ProcessSnapshotError::ProcessKindMismatch {
                process: entry.folder_name.clone(),
                expected: entry.kind,
                actual: definition.generation_type,
            });
        }
        let mut collection = match entry.kind {
            GenerationType::Amplitude => ProcessCollection::new_amplitude(),
            GenerationType::CrossSection => ProcessCollection::new_cross_section(),
        };
        for child in &entry.children {
            if selected.is_some_and(|names| !names.contains(&child.name)) {
                continue;
            }
            let child_path = OutputPath::new(path)?.join(&child.name)?;
            match &mut collection {
                ProcessCollection::Amplitudes(amplitudes) => {
                    let amplitude = Amplitude::load(child_path.as_path(), context, child.generated)
                        .map_err(|source| {
                            ProcessSnapshotError::payload(
                                "load amplitude payload",
                                &entry.folder_name,
                                child_path.as_path(),
                                source,
                            )
                        })?;
                    if amplitude.name != child.name {
                        return Err(ProcessSnapshotError::ChildNameMismatch {
                            process: entry.folder_name.clone(),
                            expected: child.name.clone(),
                            actual: amplitude.name,
                        });
                    }
                    amplitudes.insert(child.name.clone(), amplitude);
                }
                ProcessCollection::CrossSections(cross_sections) => {
                    let cross_section =
                        CrossSection::load(child_path.as_path(), context, child.generated)
                            .map_err(|source| {
                                ProcessSnapshotError::payload(
                                    "load cross-section payload",
                                    &entry.folder_name,
                                    child_path.as_path(),
                                    source,
                                )
                            })?;
                    if cross_section.name != child.name {
                        return Err(ProcessSnapshotError::ChildNameMismatch {
                            process: entry.folder_name.clone(),
                            expected: child.name.clone(),
                            actual: cross_section.name,
                        });
                    }
                    cross_sections.insert(child.name.clone(), cross_section);
                }
            }
        }
        Ok(Self {
            definition,
            collection,
            settings_history: if entry.settings_history {
                Some(Self::load_settings_history(path)?)
            } else {
                None
            },
        })
    }

    pub(super) fn write_snapshot(&mut self, processes_root: &Path) -> Result<()> {
        let entry = self.snapshot_entry()?;
        let path = OutputPath::new(processes_root)?
            .join(entry.kind_directory())?
            .join(&entry.folder_name)?;
        path.create_directory()?;
        let path = path.as_path();
        let binary = bincode::encode_to_vec(&self.definition, bincode::config::standard())?;
        fs::write(path.join("def.bin"), binary)?;
        if let Some(settings) = &self.settings_history {
            fs::write(
                path.join(SETTINGS_HISTORY_TOML),
                toml::to_string_pretty(settings)?,
            )?;
        }
        match &mut self.collection {
            ProcessCollection::Amplitudes(amplitudes) => {
                for amplitude in amplitudes.values_mut() {
                    amplitude.write_snapshot(path)?;
                }
            }
            ProcessCollection::CrossSections(cross_sections) => {
                for cross_section in cross_sections.values_mut() {
                    cross_section.write_snapshot(path)?;
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
                OutputPath::new(&p)?.ensure_directory()?;
                let p = OutputPath::new(&p)?
                    .join(&self.definition.folder_name)?
                    .as_path()
                    .to_owned();

                let r = OutputPath::new(&p)?.ensure_directory().with_context(|| {
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
                OutputPath::new(&p)?.ensure_directory()?;
                let p = OutputPath::new(&p)?
                    .join(&self.definition.folder_name)?
                    .as_path()
                    .to_owned();

                let r = OutputPath::new(&p)?.ensure_directory().with_context(|| {
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
                let path = OutputPath::new(&p)?
                    .join(&self.definition.folder_name)?
                    .as_path()
                    .to_owned();
                OutputPath::new(&path)?.ensure_directory()?;
                for amp in a.values() {
                    // Create a folder for each amplitude
                    let amp_path = OutputPath::new(&path)?
                        .join(&amp.name)?
                        .as_path()
                        .to_owned();
                    OutputPath::new(&amp_path)?
                        .ensure_directory()
                        .with_context(|| {
                            format!(
                                "Trying to create directory for amplitude {}",
                                amp_path.display()
                            )
                        })?;

                    OutputPath::new(&amp_path)?.check_tree()?;
                    amp.export_standalone(&amp_path, settings)?;
                }
            }
            ProcessCollection::CrossSections(cs) => {
                let p = path.as_ref().join("cross_sections");
                let path = OutputPath::new(&p)?
                    .join(&self.definition.folder_name)?
                    .as_path()
                    .to_owned();
                OutputPath::new(&path)?.ensure_directory()?;
                for cs in cs.values() {
                    // Create a folder for each cross section
                    let cs_path = OutputPath::new(&path)?.join(&cs.name)?.as_path().to_owned();
                    OutputPath::new(&cs_path)?
                        .ensure_directory()
                        .with_context(|| {
                            format!(
                                "Trying to create directory for cross section {}",
                                cs_path.display()
                            )
                        })?;

                    OutputPath::new(&cs_path)?.check_tree()?;
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
                let path = OutputPath::new(&p)?
                    .join(&self.definition.folder_name)?
                    .as_path()
                    .to_owned();
                OutputPath::new(&path)?.ensure_directory()?;
                for (amp_name, amp) in a {
                    // Create a folder for each amplitude
                    let amp_path = OutputPath::new(&path)?
                        .join(&amp.name)?
                        .as_path()
                        .to_owned();
                    OutputPath::new(&amp_path)?
                        .ensure_directory()
                        .with_context(|| {
                            format!(
                                "Trying to create directory for amplitude {}",
                                amp_path.display()
                            )
                        })?;

                    if settings.combine_diagrams {
                        // Save all graphs combined in one file
                        OutputPath::validate_component(amp_name)?;
                        let output_path = amp_path.join(format!("{amp_name}_graphs.dot"));
                        let mut dot = create_overwriting_file(&output_path, "amplitude graph")?;
                        for graph in amp.graphs.iter() {
                            graph.graph.dot_serialize_io(&mut dot, settings)?;
                        }
                    } else {
                        // Save each graph in its own file
                        for graph in amp.graphs.iter() {
                            OutputPath::validate_component(&graph.graph.name)?;
                            let output_path = amp_path.join(format!("{}.dot", graph.graph.name));
                            let mut dot = create_overwriting_file(&output_path, "amplitude graph")?;
                            graph.graph.dot_serialize_io(&mut dot, settings)?;
                        }
                    }
                }
            }
            ProcessCollection::CrossSections(cs) => {
                let p = path.as_ref().join("cross_sections");
                let path = OutputPath::new(&p)?
                    .join(&self.definition.folder_name)?
                    .as_path()
                    .to_owned();
                OutputPath::new(&path)?.ensure_directory()?;
                for (xs_name, cs) in cs {
                    // Create a folder for each cross section
                    let cs_path = OutputPath::new(&path)?.join(&cs.name)?.as_path().to_owned();
                    OutputPath::new(&cs_path)?
                        .ensure_directory()
                        .with_context(|| {
                            format!(
                                "Trying to create directory for cross section {}",
                                cs_path.display()
                            )
                        })?;

                    if settings.combine_diagrams {
                        // Save all graphs combined in one file
                        OutputPath::validate_component(xs_name)?;
                        let output_path = cs_path.join(format!("{xs_name}_graphs.dot"));
                        let mut dot = create_overwriting_file(&output_path, "cross section graph")?;
                        for graph in cs.supergraphs.iter() {
                            graph.graph.dot_serialize_io(&mut dot, settings)?;
                        }
                    } else {
                        // Save each supergraph in its own file
                        for graph in cs.supergraphs.iter() {
                            OutputPath::validate_component(&graph.graph.name)?;
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
        model: &Model,
        path: impl AsRef<Path>,
        integrand_name: &str,
        graph_ids: &[usize],
        settings: &UVForestExportSettings,
    ) -> Result<()> {
        let resolved = self.get_integrand(integrand_name)?;
        resolved.require_generated()?;
        let directory = match &self.collection {
            ProcessCollection::Amplitudes(_) => "amplitudes",
            ProcessCollection::CrossSections(_) => "cross_sections",
        };
        let integrand_path = OutputPath::new(path.as_ref())?
            .join(directory)?
            .join(&self.definition.folder_name)?
            .join(&resolved.canonical_name)?;
        integrand_path.ensure_directory()?;
        let integrand_path = integrand_path.as_path();

        for &graph_id in graph_ids {
            let export = self.export_uv_forest_graph(model, integrand_name, graph_id, settings)?;
            OutputPath::validate_component(&export.graph_name)?;
            let graph_name = &export.graph_name;
            let forest_path = integrand_path.join(format!("{graph_name}.forest.dot"));
            let mut forest_file = create_overwriting_file(&forest_path, "UV forest")?;
            forest_file.write_all(export.forest_dot.as_bytes())?;

            for term in export.node_terms {
                let node_dir = integrand_path
                    .join(format!("{graph_name}_nodes"))
                    .join(format!("forest_{:03}", term.forest_index));
                OutputPath::new(&node_dir)?
                    .ensure_directory()
                    .with_context(|| {
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
        _model: &Model,
    ) -> Result<Self> {
        let mut proc_definition = match definition {
            Some(definition) if definition.generation_type != generation_type => {
                return Err(eyre!(
                    "Process definition kind does not match its graph collection"
                ));
            }
            Some(definition) => definition,
            None => ProcessDefinition {
                generation_type,
                ..ProcessDefinition::default()
            },
        };
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
    use feynkit_generator::GenerationType;
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
    fn computed_uv_forest_process_export_uses_stored_sources_in_all_routes()
    -> color_eyre::Result<()> {
        use crate::{
            graph::{Graph, GroupId},
            initialisation::test_initialise,
            model::Model,
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
        let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build()?;
        let runtime = RuntimeSettings::default();
        let directory = fresh_temp_dir("computed-uv-process-export");
        for (kind, folder, source, model) in [
            (
                GenerationType::Amplitude,
                "amplitudes",
                include_str!(concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/../../tests/resources/graphs/scalar_bubble.dot"
                )),
                load_generic_model("scalars"),
            ),
            (
                GenerationType::CrossSection,
                "cross_sections",
                include_str!(concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/../../tests/resources/graphs/mass_approach_scalar_self_energy.dot"
                )),
                Model::from_json(include_str!(
                    "../../../../assets/models/json/scalars/scalars_2p_3p.json"
                ))?,
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
                    let mut graphs = Graph::import_dot(source, &model)?;
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
                                &model,
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
                            &model,
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
                                for exported in
                                    Graph::from_finalized_runtime_string(&node_dot, &model)?
                                {
                                    *exported_expressions
                                        .entry((forest_index, residue.clone()))
                                        .or_insert(Atom::Zero) +=
                                        feynkit_graph::NumeratorInterface::export(
                                            &exported.global_prefactor.num,
                                        )
                                        .0;
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
                            fs::create_dir_all(&saved)?;
                            let saved_index = processes.write_snapshot(&saved)?;
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
                            processes = ProcessList::load(&saved, &saved_index, context)?;
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
                            &model,
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
                        &model,
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
                            &model,
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

    #[test]
    fn test_proc_definition_encode() {
        let def = crate::processes::ProcessDefinition::default();
        let encoded = bincode::encode_to_vec(&def, bincode::config::standard()).unwrap();
        let model_sm = load_generic_model("sm");

        let mut symbols = Vec::new();
        symbolica::state::State::export(&mut symbols).unwrap();
        let state_map =
            symbolica::state::State::import(&mut std::io::Cursor::new(symbols), None).unwrap();

        let context = GammaLoopContextContainer {
            model: &model_sm,
            state_map: &state_map,
        };

        let (decoded, _): (crate::processes::ProcessDefinition, _) =
            bincode::decode_from_slice_with_context(&encoded, bincode::config::standard(), context)
                .unwrap();
        assert_eq!(def, decoded);
    }
}
