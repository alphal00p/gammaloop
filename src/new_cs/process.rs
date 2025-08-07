use std::{
    collections::{BTreeMap, HashSet},
    fs::{self, File},
    path::Path,
};

use ahash::HashMap;
// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::{Help, Result};
use indexmap::IndexMap;

use crate::{model::ArcParticle, new_gammaloop_integrand::NewIntegrand, GammaLoopContext};
use eyre::{eyre, Context};

use itertools::Itertools;

use crate::{
    feyngen::{FeynGenFilters, GenerationType},
    integrands::Integrand,
    model::Model,
    new_graph::Graph,
    ProcessSettings, Settings,
};

use super::{Amplitude, AmplitudeState, CrossSection, CrossSectionState, ExportSettings};

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

    pub(crate) fn folder_name(&self, model: &Model) -> String {
        let mut filename = String::new();
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
pub struct Process<A: AmplitudeState = (), C: CrossSectionState = ()> {
    pub definition: ProcessDefinition,
    pub collection: ProcessCollection<A, C>,
}

impl<A: AmplitudeState, C: CrossSectionState> Process<A, C> {
    pub(crate) fn preprocess(&mut self, model: &Model, settings: &ProcessSettings) -> Result<()> {
        self.collection
            .preprocess(model, &self.definition, settings)?;
        Ok(())
    }
}

impl Process {
    pub fn get_integrand(&self, integrand_name: impl AsRef<str>) -> Result<&NewIntegrand> {
        self.collection.get_integrand(integrand_name)
    }

    pub fn get_integrand_mut(
        &mut self,
        integrand_name: impl AsRef<str>,
    ) -> Result<&mut NewIntegrand> {
        self.collection.get_integrand_mut(integrand_name)
    }

    pub(crate) fn export_dot(&self, settings: &ExportSettings, model: &Model) -> Result<()> {
        let path = Path::new(&settings.root_folder).join(self.definition.folder_name(model));

        match &self.collection {
            ProcessCollection::Amplitudes(a) => {
                let p = path.join("amplitudes");
                fs::create_dir_all(&p).with_context(|| {
                    format!(
                        "Trying to create directory to export amplitude dot {}",
                        p.display()
                    )
                })?;
                for (_, amp) in a {
                    let mut dot = File::create_new(p.join(&format!("{}.dot", amp.name)))
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
                let mut collection: ProcessCollection<(), ()> = ProcessCollection::new_amplitude();

                if let Some(_sub_classes) = sub_classes {
                    todo!("implement seperation of processes into user defined sub classes");
                } else {
                    let mut amplitude: Amplitude<()> = Amplitude::new(name);

                    for amplitude_graph in graphs {
                        amplitude.add_graph(amplitude_graph)?;
                    }

                    collection.add_amplitude(amplitude);
                    Ok(Self {
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
                        definition,
                        collection,
                    })
                }
            }
        }
    }

    pub(super) fn generate_integrands(&mut self, settings: &Settings, model: &Model) -> Result<()> {
        self.collection.generate_integrands(settings, model)
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub enum ProcessCollection<A: AmplitudeState = (), C: CrossSectionState = ()> {
    Amplitudes(BTreeMap<String, Amplitude<A>>),
    CrossSections(BTreeMap<String, CrossSection<C>>),
}

impl<A: AmplitudeState, C: CrossSectionState> ProcessCollection<A, C> {
    fn new_amplitude() -> Self {
        Self::Amplitudes(BTreeMap::new())
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

    pub(crate) fn add_amplitude(&mut self, amplitude: Amplitude<A>) {
        match self {
            Self::Amplitudes(amplitudes) => amplitudes.insert(amplitude.name.clone(), amplitude),
            _ => panic!("Cannot add amplitude to a cross section collection"),
        };
    }

    pub(crate) fn add_cross_section(&mut self, cross_section: CrossSection<C>) {
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
        settings: &ProcessSettings,
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

    fn generate_integrands(&mut self, settings: &Settings, model: &Model) -> Result<()> {
        // let mut result = HashMap::default();
        match self {
            Self::Amplitudes(amplitudes) => {
                for (_, amplitude) in amplitudes {
                    amplitude.build_integrand(settings.clone(), model)?;
                }
            }
            Self::CrossSections(cross_sections) => {
                for (_, cross_section) in cross_sections {
                    cross_section.build_integrand(settings.clone(), model)?;
                }
            }
        }
        // result
        Ok(())
    }
}
