use std::{
    cell::RefCell,
    collections::HashSet,
    fmt::{Display, Formatter},
    fs::{self, File},
    io::Write,
    iter,
    marker::PhantomData,
    path::{Path, PathBuf},
};

use ahash::{AHashSet, HashMap};
// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use bitvec::vec::BitVec;
use color_eyre::Result;
use momtrop::SampleGenerator;

use idenso::metric::MS;

use spenso::{
    algebra::{algebraic_traits::IsZero, complex::Complex},
    iterators::Fiber,
    tensors::parametric::SerializableCompiledEvaluator,
};

use crate::{
    cff::{
        cut_expression::SuperGraphOrientationID,
        esurface::{generate_esurface_data, get_existing_esurfaces, EsurfaceDerivedData},
        expression::{
            AmplitudeOrientationID, CFFExpression, OrientationData, SubgraphOrientationID,
        },
        generation::{generate_cff_expression, get_orientations_from_subgraph},
    },
    model::ArcParticle,
    momentum::SignOrZero,
    momentum_sample::ExternalIndex,
    new_gammaloop_integrand::{
        amplitude_integrand::{AmplitudeGraphTerm, AmplitudeIntegrand},
        cross_section_integrand::OrientationEvaluator,
        GenericEvaluator, LmbMultiChannelingSetup, ParamBuilder,
    },
    new_graph::{
        get_cff_inverse_energy_product_impl, Edge, LMBext, LmbIndex, LoopMomentumBasis,
        NumHedgeData, Vertex,
    },
    signature::SignatureLike,
    subtraction::overlap::find_maximal_overlap,
    utils::{external_energy_atom_from_index, f128, ose_atom_from_index, GS, W_},
    uv::UltravioletGraph,
    GammaLoopContext, GammaLoopContextContainer,
};
use eyre::{eyre, Context};
use itertools::Itertools;
use linnet::half_edge::{
    involution::{EdgeVec, Flow, HedgePair, Orientation},
    subgraph::{HedgeNode, Inclusion, InternalSubGraph, OrientedCut},
    HedgeGraph,
};
use log::{debug, warn};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::rational::Rational,
    evaluate::{CompileOptions, FunctionMap, OptimizationSettings},
    function, parse, symbol,
};
use typed_index_collections::TiVec;

use crate::{
    cff::{
        cut_expression::{CFFCutsExpression, CutOrientationData},
        esurface::{Esurface, EsurfaceID},
        generation::generate_cff_with_cuts,
    },
    cross_section::IsPolarizable,
    feyngen::{FeynGenFilters, GenerationType},
    graph::BareGraph,
    integrands::Integrand,
    model::Model,
    momentum::{Rotatable, Rotation, RotationMethod},
    new_gammaloop_integrand::{
        cross_section_integrand::{CrossSectionGraphTerm, CrossSectionIntegrand},
        NewIntegrand,
    },
    new_graph::{ExternalConnection, FeynmanGraph, Graph},
    numerator::{NumeratorState, PythonState},
    utils::F,
    DependentMomentaConstructor, Externals, Polarizations, ProcessSettings, Settings,
};

use super::{Amplitude, AmplitudeState, CrossSection, CrossSectionState, ExportSettings};
use derive_more::{From, Into};

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
                for amp in a {
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

    pub(super) fn generate_integrands(
        &self,
        settings: Settings,
        model: &Model,
    ) -> HashMap<String, Integrand> {
        self.collection.generate_integrands(settings, model)
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub enum ProcessCollection<A: AmplitudeState = (), C: CrossSectionState = ()> {
    Amplitudes(Vec<Amplitude<A>>),
    CrossSections(Vec<CrossSection<C>>),
}

impl<A: AmplitudeState, C: CrossSectionState> ProcessCollection<A, C> {
    fn new_amplitude() -> Self {
        Self::Amplitudes(vec![])
    }

    fn new_cross_section() -> Self {
        Self::CrossSections(vec![])
    }

    pub(crate) fn add_amplitude(&mut self, amplitude: Amplitude<A>) {
        match self {
            Self::Amplitudes(amplitudes) => amplitudes.push(amplitude),
            _ => panic!("Cannot add amplitude to a cross section collection"),
        }
    }

    pub(crate) fn add_cross_section(&mut self, cross_section: CrossSection<C>) {
        match self {
            Self::CrossSections(cross_sections) => cross_sections.push(cross_section),
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
                for amplitude in amplitudes {
                    amplitude.preprocess(model, settings)?;
                }
            }
            Self::CrossSections(cross_sections) => {
                for cross_section in cross_sections {
                    cross_section.preprocess(model, process_definition)?;
                }
            }
        }
        Ok(())
    }

    fn generate_integrands(&self, settings: Settings, model: &Model) -> HashMap<String, Integrand> {
        let mut result = HashMap::default();
        match self {
            Self::Amplitudes(amplitudes) => {
                let name = "default".to_owned();
                for amplitude in amplitudes {
                    let integrand = amplitude.generate_integrand(settings.clone(), model);
                    result.insert(name.clone(), integrand);
                }
            }
            Self::CrossSections(cross_sections) => {
                let name = "default".to_owned();
                for cross_section in cross_sections {
                    let integrand = cross_section.generate_integrand(settings.clone(), model);
                    result.insert(name.clone(), integrand);
                }
            }
        }
        result
    }
}
