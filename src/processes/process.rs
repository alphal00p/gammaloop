use std::{
    collections::{BTreeMap, HashSet},
    fs::{self, File},
    path::Path,
};

// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::{Help, Result};
use log::debug;

use crate::{
    gammaloop_integrand::NewIntegrand,
    model::ArcParticle,
    settings::{runtime::LockedRuntimeSettings, GlobalSettings},
    GammaLoopContext, GammaLoopContextContainer,
};
use eyre::{eyre, Context};

use crate::{
    feyngen::{FeynGenFilters, GenerationType},
    graph::Graph,
    model::Model,
    settings::global::GenerationSettings,
    settings::RuntimeSettings,
};

use super::{Amplitude, CrossSection};

#[derive(Debug, Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ProcessDefinition {
    pub initial_pdgs: Vec<i64>, // Do we want a pub type Pdg = i64;?
    pub final_pdgs_lists: Vec<Vec<i64>>,
    pub n_unresolved: usize, // we need al this information to know what cuts are considered at runtime
    pub unresolved_cut_content: HashSet<ArcParticle>,
    pub amplitude_filters: FeynGenFilters,
    pub cross_section_filters: FeynGenFilters,
}

impl ProcessDefinition {
    pub(crate) fn new_empty() -> Self {
        Self {
            initial_pdgs: vec![],
            final_pdgs_lists: vec![],
            n_unresolved: 0,
            unresolved_cut_content: HashSet::new(),
            amplitude_filters: FeynGenFilters(vec![]),
            cross_section_filters: FeynGenFilters(vec![]),
        }
    }

    pub(crate) fn folder_name(&self, model: &Model, id: usize) -> String {
        let mut filename = String::new();
        filename.push_str(&id.to_string());
        filename.push_str(&model.name);
        filename.push('_');
        for p in self.initial_pdgs.iter() {
            filename.push_str(&model.get_particle_from_pdg(*p as isize).name);
            filename.push('_');
        }
        filename.push_str("_to_");
        let n = self.final_pdgs_lists.len();

        for (i, list) in self.final_pdgs_lists.iter().enumerate() {
            for p in list {
                filename.push_str(&model.get_particle_from_pdg(*p as isize).name);
                filename.push('_');
            }
            if i < n - 1 {
                filename.push_str("_or_");
            }
        }
        filename
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
    pub(crate) fn warm_up(&mut self) {
        self.collection.warm_up();
    }
    pub(crate) fn preprocess(&mut self, model: &Model, settings: &GlobalSettings) -> Result<()> {
        self.collection
            .preprocess(model, &self.definition, &settings.generation)?;
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

    pub fn save(
        &mut self,
        path: impl AsRef<Path>,
        override_existing: bool,
        id: usize,
        model: &Model,
    ) -> Result<()> {
        match &mut self.collection {
            ProcessCollection::Amplitudes(a) => {
                let p = path.as_ref().join("amplitudes");
                fs::create_dir_all(&p);
                let p = p.join(self.definition.folder_name(model, id));

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
                    serde_yaml::to_writer(
                        File::create(path.as_ref().join("settings_history.yaml"))?,
                        a,
                    )?;
                }

                for (_, amp) in a {
                    amp.save(&p, override_existing)?;
                }
            }
            ProcessCollection::CrossSections(a) => {}
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

    pub(crate) fn export_dot(
        &self,
        path: impl AsRef<Path>,
        model: &Model,
        id: usize,
    ) -> Result<()> {
        let p = path.as_ref().join("amplitudes");
        let path = p.join(self.definition.folder_name(model, id));
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
            ProcessCollection::CrossSections(a) => {}
        }
        Ok(())
    }

    pub fn from_graph_list(
        name: String,
        graphs: Vec<Graph>,
        generation_type: GenerationType,
        definition: ProcessDefinition,
        sub_classes: Option<Vec<Vec<String>>>,
    ) -> Result<Self> {
        match generation_type {
            GenerationType::Amplitude => {
                let mut collection: ProcessCollection = ProcessCollection::new_amplitude();

                if let Some(_sub_classes) = sub_classes {
                    todo!("implement seperation of processes into user defined sub classes");
                } else {
                    let mut amplitude: Amplitude = Amplitude::new(name);

                    for amplitude_graph in graphs {
                        amplitude.add_graph(amplitude_graph)?;
                    }

                    collection.add_amplitude(amplitude);
                    Ok(Self {
                        settings_history: None,
                        definition,
                        collection,
                    })
                }
            }
            GenerationType::CrossSection => {
                let mut collection: ProcessCollection = ProcessCollection::new_cross_section();

                if let Some(_sub_classes) = sub_classes {
                    todo!("implement seperation of processes into user defined sub classes");
                } else {
                    let mut cross_section: CrossSection = CrossSection::new(name);

                    for cross_section_graph in graphs {
                        cross_section.add_supergraph(cross_section_graph)?;
                    }

                    collection.add_cross_section(cross_section);
                    Ok(Self {
                        settings_history: None,
                        definition,
                        collection,
                    })
                }
            }
        }
    }

    pub(super) fn generate_integrands(
        &mut self,
        model: &Model,
        runtime_default: LockedRuntimeSettings,
    ) -> Result<()> {
        self.collection.generate_integrands(model, runtime_default)
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

    fn get_integrand_names(&self) -> Vec<&str> {
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

    pub(crate) fn add_amplitude(&mut self, amplitude: Amplitude) {
        match self {
            Self::Amplitudes(amplitudes) => amplitudes.insert(amplitude.name.clone(), amplitude),
            _ => panic!("Cannot add amplitude to a cross section collection"),
        };
    }

    pub(crate) fn add_cross_section(&mut self, cross_section: CrossSection) {
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
    ) -> Result<()> {
        match self {
            Self::Amplitudes(amplitudes) => {
                for (_, amplitude) in amplitudes {
                    amplitude.preprocess(model, settings)?;
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

    fn warm_up(&mut self) {
        match self {
            Self::Amplitudes(amplitudes) => {
                for (_, amplitude) in amplitudes {
                    amplitude.warm_up();
                }
            }
            Self::CrossSections(cross_sections) => {
                for (_, cross_section) in cross_sections {
                    cross_section.warm_up();
                }
            }
        }
        // Ok(())
    }

    fn generate_integrands(
        &mut self,
        model: &Model,
        runtime_default: LockedRuntimeSettings,
    ) -> Result<()> {
        // let mut result = HashMap::default();
        match self {
            Self::Amplitudes(amplitudes) => {
                for (_, amplitude) in amplitudes {
                    amplitude.build_integrand(model, runtime_default)?;
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
