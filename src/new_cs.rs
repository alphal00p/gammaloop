use std::{
    cell::RefCell,
    collections::HashSet,
    fmt::{Display, Formatter},
    fs::File,
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
        esurface::{self, generate_esurface_data, EsurfaceDerivedData},
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
    utils::{external_energy_atom_from_index, f128, ose_atom_from_index, GS, W_},
    uv::UltravioletGraph,
    GammaLoopContext, GammaLoopContextContainer,
};
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::{
    involution::{Flow, HedgePair, Orientation},
    subgraph::{HedgeNode, Inclusion, InternalSubGraph, OrientedCut},
    HedgeGraph,
};
use log::{debug, warn};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::rational::Rational,
    evaluate::{CompileOptions, FunctionMap},
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
    pub fn new_empty() -> Self {
        Self {
            initial_pdgs: vec![],
            final_pdgs_lists: vec![],
            n_unresolved: 0,
            unresolved_cut_content: HashSet::new(),
            amplitude_filters: FeynGenFilters(vec![]),
            cross_section_filters: FeynGenFilters(vec![]),
        }
    }

    pub fn folder_name(&self, model: &Model) -> String {
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
pub struct Process<S: NumeratorState = PythonState> {
    pub definition: ProcessDefinition,
    pub collection: ProcessCollection<S>,
}

impl<S: NumeratorState> Process<S> {
    pub fn preprocess(&mut self, model: &Model, settings: &ProcessSettings) -> Result<()> {
        self.collection
            .preprocess(model, &self.definition, settings)?;
        Ok(())
    }
}

impl Process {
    pub fn export_dot(&self, settings: &ExportSettings, model: &Model) -> Result<()> {
        let path = Path::new(&settings.root_folder).join(self.definition.folder_name(model));

        match &self.collection {
            ProcessCollection::Amplitudes(a) => {
                let p = path.join("amplitudes");
                for amp in a {
                    let mut dot = File::create_new(p.join(&format!("{}.dot", amp.name)))?;
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
                let mut collection: ProcessCollection<PythonState> =
                    ProcessCollection::new_amplitude();

                if let Some(_sub_classes) = sub_classes {
                    todo!("implement seperation of processes into user defined sub classes");
                } else {
                    let mut amplitude: Amplitude<PythonState> = Amplitude::new(name);

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
                let mut collection: ProcessCollection<PythonState> =
                    ProcessCollection::new_cross_section();

                if let Some(_sub_classes) = sub_classes {
                    todo!("implement seperation of processes into user defined sub classes");
                } else {
                    let mut cross_section: CrossSection<PythonState> = CrossSection::new(name);

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

    pub fn from_bare_graph_list(
        name: String,
        bare_graphs: Vec<BareGraph>,
        generation_type: GenerationType,
        definition: ProcessDefinition,
        sub_classes: Option<Vec<Vec<String>>>,
    ) -> Result<Self> {
        let graphs = bare_graphs.into_iter().map(Graph::from).collect_vec();
        Self::from_graph_list(name, graphs, generation_type, definition, sub_classes)
    }

    fn generate_integrands(&self, settings: Settings, model: &Model) -> HashMap<String, Integrand> {
        self.collection.generate_integrands(settings, model)
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ProcessList<S: NumeratorState = PythonState> {
    pub processes: Vec<Process<S>>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct ExportSettings {
    pub root_folder: PathBuf,
}

impl Default for ProcessList {
    fn default() -> Self {
        Self::new()
    }
}

// Python Api
impl ProcessList {
    /// Generates a new empty process list
    pub fn new() -> Self {
        ProcessList { processes: vec![] }
    }

    pub fn add_process(&mut self, process: Process) {
        self.processes.push(process);
    }

    /// imports a process list from a folder
    pub fn import(_settings: ExportSettings) -> Result<Self> {
        Ok(Self::new())
    }

    ///preprocesses the process list according to the settings
    pub fn preprocess(&mut self, model: &Model, settings: ProcessSettings) -> Result<()> {
        for process in self.processes.iter_mut() {
            process.preprocess(model, &settings)?;
        }

        Ok(())
    }

    pub fn generate_integrands(
        &self,
        settings: Settings,
        model: &Model,
    ) -> HashMap<String, Integrand> {
        let mut result = HashMap::default();

        for process in self.processes.iter() {
            let integrands = process.generate_integrands(settings.clone(), model);
            result.extend(integrands);
        }

        result
    }

    /// exports a process list to a folder
    pub fn export_amplitudes(
        &self,
        amplitude_names: &[String],
        settings: &ExportSettings,
    ) -> Result<usize> {
        let mut n_exported = 0;
        for process in self.processes.iter() {
            let process_definition_data =
                bincode::encode_to_vec(&process.definition, bincode::config::standard())?;

            let path = settings
                .root_folder
                .join("sources")
                .join("process_definition.bin");

            std::fs::write(path, &process_definition_data)?;

            n_exported += process
                .collection
                .export_amplitudes(amplitude_names, settings)?;
        }
        Ok(n_exported)
    }

    pub fn export_cross_sections(
        &self,
        cross_section_names: &[String],
        settings: &ExportSettings,
    ) -> Result<usize> {
        let mut n_exported = 0;

        for process in self.processes.iter() {
            let process_definition_data =
                bincode::encode_to_vec(&process.definition, bincode::config::standard())?;

            let path = settings
                .root_folder
                .join("sources")
                .join("process_definition.bin");

            std::fs::write(path, &process_definition_data)?;

            n_exported += process
                .collection
                .export_cross_sections(cross_section_names, settings)?;
        }

        Ok(n_exported)
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub enum ProcessCollection<S: NumeratorState = PythonState> {
    Amplitudes(Vec<Amplitude<S>>),
    CrossSections(Vec<CrossSection<S>>),
}

impl<S: NumeratorState> ProcessCollection<S> {
    fn new_amplitude() -> Self {
        Self::Amplitudes(vec![])
    }

    fn new_cross_section() -> Self {
        Self::CrossSections(vec![])
    }

    pub fn add_amplitude(&mut self, amplitude: Amplitude<S>) {
        match self {
            Self::Amplitudes(amplitudes) => amplitudes.push(amplitude),
            _ => panic!("Cannot add amplitude to a cross section collection"),
        }
    }

    pub fn add_cross_section(&mut self, cross_section: CrossSection<S>) {
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

    fn export_amplitudes(
        &self,
        amplitude_names: &[String],
        settings: &ExportSettings,
    ) -> Result<usize> {
        let mut n_exported = 0;
        match self {
            Self::Amplitudes(amplitudes) => {
                for amplitude in amplitudes {
                    if amplitude_names.contains(&amplitude.name.to_string()) {
                        amplitude.export(settings)?;
                        n_exported += 1;
                    }
                }
            }
            _ => {}
        }
        Ok(n_exported)
    }

    fn export_cross_sections(
        &self,
        cross_section_names: &[String],
        settings: &ExportSettings,
    ) -> Result<usize> {
        let mut n_exported = 0;
        match self {
            Self::CrossSections(cross_sections) => {
                for cross_section in cross_sections {
                    if cross_section_names.contains(&cross_section.name.to_string()) {
                        cross_section.export(settings)?;
                        n_exported += 1;
                    }
                }
            }
            _ => {}
        }
        Ok(n_exported)
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct Amplitude<S: NumeratorState = PythonState> {
    pub name: String,
    pub graphs: Vec<AmplitudeGraph<S>>,
    pub external_particles: Vec<ArcParticle>,
    pub external_signature: SignatureLike<ExternalIndex>,
}

impl<S: NumeratorState> Amplitude<S> {
    pub fn preprocess(&mut self, model: &Model, settings: &ProcessSettings) -> Result<()> {
        for amplitude_graph in self.graphs.iter_mut() {
            amplitude_graph.preprocess(model, settings)?;
        }
        Ok(())
    }

    pub fn generate_integrand(&self, settings: Settings, model: &Model) -> Integrand {
        let terms = self
            .graphs
            .iter()
            .map(|graph| graph.generate_term_for_graph(&settings))
            .collect_vec();

        let rotations: Vec<Rotation> = Some(Rotation::new(RotationMethod::Identity))
            .into_iter()
            .chain(
                settings
                    .stability
                    .rotation_axis
                    .iter()
                    .map(|axis| Rotation::new(axis.rotation_method())),
            )
            .collect();

        let orig_polarizations = self.polarizations(&settings.kinematics.externals);

        let polarizations = rotations
            .iter()
            .map(|r| orig_polarizations.rotate(r))
            .collect();

        let amplitude_integrand = AmplitudeIntegrand {
            settings,
            rotations,
            polarizations,
            graph_terms: terms,
            external_signature: self.external_signature.clone(),
            model_parameter_cache: model.generate_values(),
        };

        Integrand::NewIntegrand(NewIntegrand::Amplitude(amplitude_integrand))
    }

    pub fn export(&self, settings: &ExportSettings) -> Result<()> {
        let path = Path::new(&settings.root_folder)
            .join("sources")
            .join("amplitudes")
            .join(self.name.as_str());

        let data = bincode::encode_to_vec(self, bincode::config::standard())?;
        std::fs::write(
            path.clone().join(&format!("amplitude_{}.bin", self.name)),
            &data,
        )?;

        let mut dot = File::create_new(path.clone().join(&format!("amplitude_{}.dot", self.name)))?;
        self.write_dot(&mut dot)?;
        Ok(())
    }

    pub fn write_dot<W: std::io::Write>(&self, writer: &mut W) -> Result<(), std::io::Error> {
        for graph in &self.graphs {
            graph.write_dot(writer)?;
        }
        Ok(())
    }

    pub fn load_from_file(
        name: &str,
        root_folder: &str,
        context: GammaLoopContextContainer,
    ) -> Result<Self> {
        let path = Path::new(root_folder)
            .join("sources")
            .join("amplitudes")
            .join(name);

        let data = std::fs::read(path.join(format!("amplitude_{}.bin", name)))?;
        let (amplitude, _): (Self, _) =
            bincode::decode_from_slice_with_context(&data, bincode::config::standard(), context)?;
        Ok(amplitude)
    }
}

impl<S: NumeratorState> IsPolarizable for Amplitude<S> {
    fn polarizations(&self, externals: &Externals) -> Polarizations {
        externals.generate_polarizations(
            &self.external_particles,
            DependentMomentaConstructor::Amplitude(&self.external_signature),
        )
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait= GammaLoopContext)]
pub struct AmplitudeGraph<S: NumeratorState> {
    pub graph: Graph,
    pub derived_data: AmplitudeDerivedData<S>,
}

impl<S: NumeratorState> AmplitudeGraph<S> {
    pub(crate) fn write_dot<W: std::io::Write>(
        &self,
        writer: &mut W,
    ) -> Result<(), std::io::Error> {
        self.graph.dot_serialize_io(writer)
    }

    pub fn new(graph: Graph) -> Self {
        AmplitudeGraph {
            graph,
            derived_data: AmplitudeDerivedData {
                cff_expression: None,
                bare_cff_evaluator: None,
                bare_cff_orientation_evaluatos: None,
                _temp_numerator: None,
                lmbs: None,
                tropical_sampler: None,
                multi_channeling_setup: None,
                threshold_counterterms: None,
                esurface_data: None,
            },
        }
    }

    fn generate_cff(&mut self) -> Result<()> {
        let shift_rewrite = self
            .graph
            .underlying
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let cff_expression = generate_cff_expression(
            &self.graph.underlying,
            &shift_rewrite,
            &self.graph.underlying.dummy_list(),
        )?;
        self.derived_data.cff_expression = Some(cff_expression);

        Ok(())
    }

    pub fn preprocess(&mut self, model: &Model, settings: &ProcessSettings) -> Result<()> {
        self.generate_cff()?;
        self.build_evaluator(model);
        self.build_evaluator_for_orientations(model)?;
        self.build_tropical_sampler(settings)?;
        self.build_loop_momentum_bases();
        self.build_multi_channeling_channels();
        self.build_esurface_derived_data()?;
        //self.build_threshold_counterterms(model);

        Ok(())
    }

    pub fn build_esurface_derived_data(&mut self) -> Result<()> {
        let lmbs = self.derived_data.lmbs.as_ref().unwrap();
        let esurfaces = &self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .esurface_cache;

        let esurface_data = generate_esurface_data(&self.graph, lmbs, esurfaces)?;
        Ok(self.derived_data.esurface_data = Some(esurface_data))
    }

    fn build_multi_channeling_channels(&mut self) {
        let channels = self
            .graph
            .build_multi_channeling_channels(self.derived_data.lmbs.as_ref().unwrap());

        self.derived_data.multi_channeling_setup = Some(channels)
    }

    fn get_params(&self, model: &Model) -> Vec<Atom> {
        // the float type does not matter here
        let mut param_builder = ParamBuilder::<f64>::new();

        // this is wrong if we allow for vacuum graphs
        param_builder.external_energies_atom(self.graph.underlying.get_external_energy_atoms());

        // spatial components of external momenta
        param_builder.external_spatial_atom(
            self.graph
                .underlying
                .iter_edges()
                .flat_map(|(pair, edge_id, _)| {
                    if let HedgePair::Unpaired { .. } = pair {
                        let i64_id = Into::<usize>::into(edge_id) as i64;
                        vec![
                            function!(GS.external_mom, i64_id, 1),
                            function!(GS.external_mom, i64_id, 2),
                            function!(GS.external_mom, i64_id, 3),
                        ]
                    } else {
                        vec![]
                    }
                })
                .collect(),
        );

        // spatial EMR
        param_builder.emr_spatial_atom(
            self.graph
                .underlying
                .iter_edges()
                .flat_map(|(pair, edge_id, _)| {
                    if let HedgePair::Paired { .. } = pair {
                        let i64_id = Into::<usize>::into(edge_id) as i64;
                        vec![
                            function!(GS.emr_vec, i64_id, 1),
                            function!(GS.emr_vec, i64_id, 2),
                            function!(GS.emr_vec, i64_id, 3),
                        ]
                    } else {
                        vec![]
                    }
                })
                .collect(),
        );

        param_builder.model_parameters_atom(model.generate_params());

        param_builder.m_uv_atom(Atom::var(GS.m_uv));

        param_builder.mu_r_sq_atom(Atom::var(GS.mu_r_sq));
        param_builder.build_params_amplitude().unwrap()
    }

    fn get_function_map(&self) -> FunctionMap {
        let mut fn_map = FunctionMap::new();
        let pi_rational = Rational::from(std::f64::consts::PI);
        fn_map.add_constant(Atom::PI.into(), pi_rational.into());
        fn_map
    }

    fn add_additional_factors_to_cff_atom(&self, cff_atom: &Atom) -> Atom {
        // let inverse_energy_product = self.graph.underlying.get_cff_inverse_energy_product();
        let factors_of_pi =
            (Atom::var(Atom::PI) * 2).npow(3 * self.graph.underlying.get_loop_number() as i64);

        let result = cff_atom * &self.graph.multiplicity / factors_of_pi;
        debug!("result: {}", result);
        result
    }

    pub fn build_all_orientations_integrand_atom(&self) -> Atom {
        let wood = self.graph.wood(&self.graph.underlying.no_dummy());
        let mut forest = wood.unfold(&self.graph, &self.graph.loop_momentum_basis);

        let canonize_esurface = self
            .graph
            .underlying
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let orientations: TiVec<AmplitudeOrientationID, OrientationData> = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .orientations
            .iter()
            .map(|a| a.data.clone())
            .collect();

        forest.compute(
            &self.graph,
            &self.graph.underlying.no_dummy(),
            &orientations,
            &canonize_esurface,
            &[],
        );

        let mut result = Atom::new();

        for (_orientation_id, orientation_data) in orientations.iter_enumerated() {
            let mut orientation_expr = forest.local_expr(
                &self.graph,
                &self.graph.underlying.no_dummy(),
                orientation_data,
                &canonize_esurface,
                &orientations,
                &[],
            );

            orientation_expr = do_replacement_rules(orientation_expr, &self.graph.underlying);
            orientation_expr = self.add_additional_factors_to_cff_atom(&orientation_expr);

            println!("orientation: {:?}", orientation_data.orientation);
            println!("orientation_expr: {}", orientation_expr);
            //panic!("atom test: {}", spenso_lor_atom(3, 20, GS.dim));
            result += orientation_expr;
        }

        result
    }

    fn build_threshold_counterterms(&mut self, model: &Model) {
        let mut counterterms: TiVec<EsurfaceID, Atom> = TiVec::new();
        let canonize_esurface = self
            .graph
            .underlying
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let full_orientation_list = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .orientations
            .iter()
            .map(|o| o.data.clone())
            .collect::<TiVec<AmplitudeOrientationID, _>>();

        for (esurface_id, esurface) in self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .esurface_cache
            .iter_enumerated()
        {
            if esurface.external_shift.is_empty() {
                // these will never satsify the threshold condition
                // so we can skip them
                counterterms.push(Atom::new());
                continue;
            }

            let (circled, complement) = esurface.get_subgraph_components(&self.graph.underlying);
            let orientations = self
                .derived_data
                .cff_expression
                .as_ref()
                .unwrap()
                .get_orientations_with_esurface(esurface_id);

            let first_orientation = &self
                .derived_data
                .cff_expression
                .as_ref()
                .unwrap()
                .orientations[orientations[0]]
                .data
                .orientation;

            let orientation_of_edges_in_esurface = esurface
                .energies
                .iter()
                .map(|e| first_orientation[*e])
                .collect_vec();

            assert!(orientations.iter().all(|o| {
                let or = &self
                    .derived_data
                    .cff_expression
                    .as_ref()
                    .unwrap()
                    .orientations[*o]
                    .data
                    .orientation;

                let orientation_of_esurface_in_this_orientation =
                    esurface.energies.iter().map(|e| or[*e]).collect_vec();

                if orientation_of_edges_in_esurface != orientation_of_esurface_in_this_orientation {
                    println!("{:?}", orientation_of_edges_in_esurface);
                    println!("{:?}", orientation_of_esurface_in_this_orientation);
                    println!("esurface shift: {:?}", esurface.external_shift);
                    false
                } else {
                    true
                }
            }));

            let cut_momentum_basis = self.derived_data.esurface_data.as_ref().unwrap()[esurface_id]
                .as_ref()
                .unwrap()
                .cut_momentum_basis;

            let circled_wood = self.graph.wood(&circled);
            let complement_wood = self.graph.wood(&complement);

            let mut circled_forest =
                circled_wood.unfold(&self.graph, &self.graph.loop_momentum_basis);

            let mut complement_forest =
                complement_wood.unfold(&self.graph, &self.graph.loop_momentum_basis);

            let reverse_dangling = esurface
                .energies
                .iter()
                .zip(orientation_of_edges_in_esurface)
                .filter_map(|(e, o)| {
                    if o == Orientation::Reversed {
                        Some(*e)
                    } else {
                        None
                    }
                })
                .collect_vec();

            let circled_orientations =
                get_orientations_from_subgraph(&self.graph.underlying, &circled, &reverse_dangling)
                    .into_iter()
                    .map(|cff_graph| cff_graph.global_orientation)
                    .collect::<TiVec<SubgraphOrientationID, _>>();

            let complement_orientations = get_orientations_from_subgraph(
                &self.graph.underlying,
                &complement,
                &reverse_dangling,
            )
            .into_iter()
            .map(|cff_graph| cff_graph.global_orientation)
            .collect::<TiVec<SubgraphOrientationID, _>>();

            circled_forest.compute(
                &self.graph,
                &circled,
                &circled_orientations,
                &canonize_esurface,
                &esurface.energies,
            );

            complement_forest.compute(
                &self.graph,
                &complement,
                &complement_orientations,
                &canonize_esurface,
                &esurface.energies,
            );

            let mut counterterm = Atom::new();
            for orientation in orientations {
                let circled_expr = circled_forest.local_expr(
                    &self.graph,
                    &circled,
                    &full_orientation_list[orientation],
                    &canonize_esurface,
                    &circled_orientations,
                    &esurface.energies,
                );

                let complement_expr = complement_forest.local_expr(
                    &self.graph,
                    &complement,
                    &full_orientation_list[orientation],
                    &canonize_esurface,
                    &complement_orientations,
                    &esurface.energies,
                );

                let mut orientation_result = circled_expr * complement_expr;
            }

            counterterms.push(counterterm);
        }
    }

    pub fn build_evaluator(&mut self, model: &Model) {
        let replace_dots = self.build_all_orientations_integrand_atom();
        let params = self.get_params(model);

        let function_map = self.get_function_map();

        let evaluator = GenericEvaluator::new(&replace_dots, &function_map, &params, Some(1));
        self.derived_data.bare_cff_evaluator = Some(evaluator)
    }

    fn build_loop_momentum_bases(&mut self) {
        let lmbs = self
            .graph
            .generate_loop_momentum_bases(&self.graph.underlying.no_dummy());

        self.derived_data.lmbs = Some(lmbs)
    }

    fn build_tropical_sampler(&mut self, process_settings: &ProcessSettings) -> Result<()> {
        let num_virtual_loop_edges = self.graph.iter_loop_edges().count();
        let num_loops = self.graph.loop_momentum_basis.loop_edges.len();
        let target_omega = process_settings
            .tropical_subgraph_table_settings
            .target_omega;

        let weight = (target_omega + (3 * num_loops) as f64 / 2.) / num_virtual_loop_edges as f64;

        debug!(
            "Building tropical subgraph table with all edge weights set to: {}",
            weight
        );

        let tropical_edges = self
            .graph
            .iter_loop_edges()
            .map(|(pair, _edge_id, edge)| {
                let is_massive = match edge.data.particle.0.mass.value {
                    Some(complex_mass) => complex_mass.is_non_zero(),
                    None => false,
                };

                let vertices = match pair {
                    HedgePair::Paired { source, sink } => (
                        self.graph.underlying.node_id(source).0 as u8,
                        self.graph.underlying.node_id(sink).0 as u8,
                    ),
                    _ => unreachable!(),
                };

                momtrop::Edge {
                    is_massive,
                    weight,
                    vertices,
                }
            })
            .collect_vec();

        let mut external_vertices_pool = AHashSet::new();

        for (pair, _, _) in self.graph.iter_non_loop_edges() {
            match pair {
                HedgePair::Paired { source, sink } => {
                    let source_id = self.graph.underlying.node_id(source).0 as u8;
                    let sink_id = self.graph.underlying.node_id(sink).0 as u8;

                    external_vertices_pool.insert(source_id);
                    external_vertices_pool.insert(sink_id);
                }
                HedgePair::Unpaired { hedge, .. } => {
                    let id = self.graph.underlying.node_id(hedge).0 as u8;
                    external_vertices_pool.insert(id);
                }
                _ => unreachable!(),
            }
        }

        let mut external_vertices = vec![];

        for tropical_edge in &tropical_edges {
            if external_vertices_pool.contains(&tropical_edge.vertices.0) {
                external_vertices.push(tropical_edge.vertices.0);
            }

            if external_vertices_pool.contains(&tropical_edge.vertices.1) {
                external_vertices.push(tropical_edge.vertices.1);
            }
        }

        let tropical_graph = momtrop::Graph {
            edges: tropical_edges,
            externals: external_vertices,
        };

        let loop_part = self
            .graph
            .iter_loop_edges()
            .map(|(_, edge_id, _edge)| {
                self.graph.loop_momentum_basis.edge_signatures[edge_id]
                    .internal
                    .clone()
                    .to_momtrop_format()
            })
            .collect_vec();

        let sampler = tropical_graph
            .build_sampler(loop_part)
            .map_err(|e| eyre!(e))?;

        Ok(self.derived_data.tropical_sampler = Some(sampler))
    }

    fn build_evaluator_for_orientations(&mut self, model: &Model) -> Result<()> {
        let params = self.get_params(model);
        let function_map = self.get_function_map();

        let evaluators = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .get_orientation_atoms()
            .iter()
            .map(|orientation_atom_unsubstituted| {
                let atom_no_prefactor = self
                    .derived_data
                    .cff_expression
                    .as_ref()
                    .unwrap()
                    .surfaces
                    .substitute_energies(orientation_atom_unsubstituted, &[]);

                let atom = self.add_additional_factors_to_cff_atom(&atom_no_prefactor);
                let replacements = self.graph.underlying.get_ose_replacements();
                let replaced_atom = atom.replace_multiple(&replacements);
                let replace_dots = replaced_atom
                    .replace(function!(
                        MS.dot,
                        function!(GS.emr_vec, W_.x_),
                        function!(GS.emr_vec, W_.y_)
                    ))
                    .with(
                        -(function!(GS.emr_vec, W_.x_, 1) * function!(GS.emr_vec, W_.y_, 1)
                            + function!(GS.emr_vec, W_.x_, 2) * function!(GS.emr_vec, W_.y_, 2)
                            + function!(GS.emr_vec, W_.x_, 3) * function!(GS.emr_vec, W_.y_, 3)),
                    );

                GenericEvaluator::new(&replace_dots, &function_map, &params, Some(1))
            })
            .collect();

        Ok(self.derived_data.bare_cff_orientation_evaluatos = Some(evaluators))
    }

    fn generate_term_for_graph(&self, settings: &Settings) -> AmplitudeGraphTerm {
        let estimated_scale = self
            .graph
            .underlying
            .expected_scale(settings.kinematics.e_cm);

        AmplitudeGraphTerm {
            bare_cff_evaluator: self.derived_data.bare_cff_evaluator.clone().unwrap(),
            bare_cff_orientation_evaluators: self
                .derived_data
                .bare_cff_orientation_evaluatos
                .clone()
                .unwrap(),
            tropical_sampler: self.derived_data.tropical_sampler.clone().unwrap(),
            graph: self.graph.clone(),
            multi_channeling_setup: self.derived_data.multi_channeling_setup.clone().unwrap(),
            lmbs: self.derived_data.lmbs.clone().unwrap(),
            estimated_scale,
        }
    }
}

#[derive(Clone, Encode, Decode)]
pub struct AmplitudeDerivedData<S: NumeratorState> {
    pub cff_expression: Option<CFFExpression<AmplitudeOrientationID>>,
    pub bare_cff_evaluator: Option<GenericEvaluator>,
    pub bare_cff_orientation_evaluatos: Option<TiVec<AmplitudeOrientationID, GenericEvaluator>>,
    pub threshold_counterterms: Option<TiVec<EsurfaceID, GenericEvaluator>>,
    pub _temp_numerator: Option<PhantomData<S>>,
    pub lmbs: Option<TiVec<LmbIndex, LoopMomentumBasis>>,
    pub tropical_sampler: Option<SampleGenerator<3>>,
    pub esurface_data: Option<EsurfaceDerivedData>,
    pub multi_channeling_setup: Option<LmbMultiChannelingSetup>,
}

impl<S: NumeratorState> Amplitude<S> {
    pub fn from_dot_string<Str: AsRef<str>>(s: Str, name: String, model: &Model) -> Result<Self> {
        let graphs = Graph::from_string(s, model)?;

        let mut amp = Amplitude::new(name);
        for g in graphs {
            amp.add_graph(g)?;
        }
        Ok(amp)
    }

    pub fn from_dot_file<'a, P>(p: P, name: String, model: &Model) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let graphs = Graph::from_file(p, model)?;

        let mut amp = Amplitude::new(name);
        for g in graphs {
            amp.add_graph(g)?;
        }
        Ok(amp)
    }

    pub fn new(name: String) -> Self {
        Self {
            name,
            graphs: vec![],
            external_particles: vec![],
            external_signature: SignatureLike::from_iter(iter::empty::<i8>()),
        }
    }

    pub fn add_graph(&mut self, graph: Graph) -> Result<()> {
        let new_external_particels = graph.underlying.get_external_partcles();
        let new_external_signature = graph.underlying.get_external_signature();

        if !self.graphs.is_empty() {
            if self.external_particles != new_external_particels {
                return Err(eyre!("amplitude graph has different number of externals"));
            }

            if self.external_signature != new_external_signature {
                return Err(eyre!("wrong external signature"));
            }
        } else {
            self.external_particles = new_external_particels;
            self.external_signature = new_external_signature;
        }

        self.graphs.push(AmplitudeGraph::new(graph));

        //  TODO: validate that the graph is compatible
        Ok(())
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSection<S: NumeratorState> {
    pub name: String,
    pub supergraphs: Vec<CrossSectionGraph<S>>,
    pub external_particles: Vec<ArcParticle>,
    pub external_connections: Vec<ExternalConnection>,
    pub n_incmoming: usize,
}

impl<S: NumeratorState> CrossSection<S> {
    pub fn new(name: String) -> Self {
        Self {
            name,
            supergraphs: vec![],
            external_connections: vec![],
            external_particles: vec![],
            n_incmoming: 0,
        }
    }

    pub fn add_supergraph(&mut self, supergraph: Graph) -> Result<()> {
        if self.external_particles.is_empty() {
            let external_particles = supergraph.underlying.get_external_partcles();
            if external_particles.len() % 2 != 0 {
                return Err(eyre!(
                    "expected even number of externals for forward scattering graph"
                ));
            }
            self.external_particles = external_particles;
            self.n_incmoming = self.external_particles.len() / 2;
        } else if self.external_particles != supergraph.underlying.get_external_partcles() {
            return Err(eyre!(
                "attempt to add supergraph with differnt external particles"
            ));
        }

        let cross_section_graph = CrossSectionGraph::new(supergraph);
        self.supergraphs.push(cross_section_graph);

        // TODO: validate that the graph is compatible
        Ok(())
    }

    pub fn preprocess(
        &mut self,
        model: &Model,
        process_definition: &ProcessDefinition,
    ) -> Result<()> {
        for supergraph in &mut self.supergraphs {
            supergraph.preprocess(model, process_definition)?;
        }
        Ok(())
    }

    pub fn generate_integrand(&self, settings: Settings, model: &Model) -> Integrand {
        let terms = self
            .supergraphs
            .iter()
            .map(|sg| sg.generate_term_for_graph(&settings))
            .collect_vec();

        let rotations: Vec<Rotation> = Some(Rotation::new(RotationMethod::Identity))
            .into_iter()
            .chain(
                settings
                    .stability
                    .rotation_axis
                    .iter()
                    .map(|axis| Rotation::new(axis.rotation_method())),
            )
            .collect(); // want this to include the identity rotation (i.e the first sample)

        let orig_polarizations = self.polarizations(&settings.kinematics.externals);

        let polarizations = rotations
            .iter()
            .map(|r| orig_polarizations.rotate(r))
            .collect();

        let model_parameter_cache = model.generate_values();

        let cross_section_integrand = CrossSectionIntegrand {
            rotations,
            external_connections: self.external_connections.clone(),
            n_incoming: self.n_incmoming,
            polarizations,
            settings,
            graph_terms: terms,
            model_parameter_cache,
        };

        Integrand::NewIntegrand(NewIntegrand::CrossSection(cross_section_integrand))
    }

    fn export(&self, settings: &ExportSettings) -> Result<()> {
        let path = Path::new(&settings.root_folder)
            .join("sources")
            .join("cross_sections")
            .join(self.name.as_str());

        let data = bincode::encode_to_vec(self, bincode::config::standard())?;
        std::fs::write(
            path.clone()
                .join(&format!("cross_section_{}.bin", self.name)),
            &data,
        )?;

        Ok(())
    }

    pub fn load_from_file(
        name: &str,
        root_folder: &str,
        context: GammaLoopContextContainer,
    ) -> Result<Self> {
        let path = Path::new(root_folder)
            .join("sources")
            .join("cross_sections")
            .join(name);

        let data = std::fs::read(path.join(format!("cross_section_{}.bin", name)))?;
        let (cross_section, _): (Self, _) =
            bincode::decode_from_slice_with_context(&data, bincode::config::standard(), context)?;
        Ok(cross_section)
    }
}

impl<S: NumeratorState> IsPolarizable for CrossSection<S> {
    fn polarizations(&self, externals: &Externals) -> Polarizations {
        externals.generate_polarizations(
            &self.external_particles,
            DependentMomentaConstructor::CrossSection {
                external_connections: &self.external_connections,
            },
        )
    }
}

#[derive(
    Debug, Clone, Serialize, Decode, Deserialize, From, Into, Hash, PartialEq, Copy, Eq, Encode,
)]
pub struct CutId(pub usize);

impl Display for CutId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSectionGraph<S: NumeratorState = PythonState> {
    pub graph: Graph,
    pub source_nodes: HedgeNode,
    pub target_nodes: HedgeNode,
    pub cuts: TiVec<CutId, CrossSectionCut>,
    pub cut_esurface: TiVec<CutId, Esurface>,
    pub cut_esurface_id_map: TiVec<CutId, EsurfaceID>,
    pub derived_data: CrossSectionDerivedData<S>,
}

impl<S: NumeratorState> CrossSectionGraph<S> {
    pub fn new(graph: Graph) -> Self {
        let mut source_nodes = AHashSet::new();
        let mut target_nodes = AHashSet::new();

        for (hedge_pair, _, _) in graph.underlying.iter_edges() {
            match hedge_pair {
                HedgePair::Unpaired { hedge, flow } => {
                    let node_id = graph.underlying.node_id(hedge);
                    match flow {
                        Flow::Source => {
                            target_nodes.insert(node_id);
                        }
                        Flow::Sink => {
                            source_nodes.insert(node_id);
                        }
                    }
                }
                _ => continue,
            }
        }

        assert_eq!(source_nodes.len(), target_nodes.len());

        let source_node_vec = source_nodes.into_iter().collect_vec();
        let target_node_vec = target_nodes.into_iter().collect_vec();

        let source_node = graph
            .underlying
            .combine_to_single_hedgenode(&source_node_vec);

        let target_node = graph
            .underlying
            .combine_to_single_hedgenode(&target_node_vec);

        Self {
            graph,
            source_nodes: source_node,
            target_nodes: target_node,
            cuts: TiVec::new(),
            cut_esurface: TiVec::new(),
            cut_esurface_id_map: TiVec::new(),
            derived_data: CrossSectionDerivedData::<S>::new_empty(),
        }
    }

    pub fn preprocess(
        &mut self,
        model: &Model,
        process_definition: &ProcessDefinition,
    ) -> Result<()> {
        self.generate_cuts(model, process_definition)?;
        self.generate_esurface_cuts();
        self.generate_cff()?;
        self.update_surface_cache();

        //self.build_cut_evaluators(model, None);
        //self.build_orientation_evaluators(model);
        self.build_lmbs();
        self.build_esurface_derived_data()?;
        Ok(self.build_multi_channeling_channels())
    }

    pub fn build_esurface_derived_data(&mut self) -> Result<()> {
        let lmbs = self.derived_data.lmbs.as_ref().unwrap();
        let esurfaces = &self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .esurface_cache;

        let esurface_data = generate_esurface_data(&self.graph, lmbs, esurfaces)?;
        Ok(self.derived_data.esurface_data = Some(esurface_data))
    }

    pub fn update_surface_cache(&mut self) {
        let esurface_cache = &mut self
            .derived_data
            .cff_expression
            .as_mut()
            .unwrap()
            .surfaces
            .esurface_cache;

        // if a cut was not generated during cff, we still add it to the surface cache such that it has an esurface_id
        for esurface in self.cut_esurface.iter() {
            if let Some(esurface_id) = esurface_cache.iter().position(|e| e == esurface) {
                self.cut_esurface_id_map.push(esurface_id.into());
            } else {
                self.cut_esurface_id_map.push(esurface_cache.len().into());
                esurface_cache.push(esurface.clone());
            }
        }
    }

    fn generate_cff(&mut self) -> Result<()> {
        // hardcorde 1 to n for now
        debug!("generating cff");

        let shift_rewrite = self
            .graph
            .underlying
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let cff_cut_expression =
            generate_cff_with_cuts(&self.graph.underlying, &shift_rewrite, &self.cuts)?;

        Ok(self.derived_data.cff_expression = Some(cff_cut_expression))
    }

    fn generate_cuts(
        &mut self,
        model: &Model,
        process_definition: &ProcessDefinition,
    ) -> Result<()> {
        debug!("generatig cuts for graph: {}", self.graph.name);

        let all_st_cuts = self
            .graph
            .underlying
            .all_cuts(self.source_nodes.clone(), self.target_nodes.clone());

        debug!("num s_t cuts: {}", all_st_cuts.len());

        let mut cuts: TiVec<CutId, CrossSectionCut> = all_st_cuts
            .into_iter()
            .map(|(left, cut, right)| CrossSectionCut { cut, left, right })
            .filter_map(
                |cut| match cut.is_valid_for_process(self, process_definition, model) {
                    Ok(true) => Some(Ok(cut)),
                    Ok(false) => None,
                    Err(e) => Some(Err(e)),
                },
            )
            .collect::<Result<_>>()?;

        cuts.sort_by(|a, b| a.cut.cmp(&b.cut));

        self.cuts = cuts;

        debug!(
            "found {} cuts for graph: {}",
            self.cuts.len(),
            self.graph.name
        );

        Ok(())
    }

    fn generate_esurface_cuts(&mut self) {
        debug!("generating esurfaces for cuts");

        let esurfaces: TiVec<CutId, Esurface> = self
            .cuts
            .iter()
            .map(|cut| Esurface::new_from_cut_left(&self.graph.underlying, cut))
            .collect();

        self.cut_esurface = esurfaces;
    }

    fn build_left_right_amplitudes(&self, cut: CutId) -> (Atom, Atom) {
        let (left_amplitude, right_amplitude) = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .to_atom_for_cut(cut);

        let cut_edges = self
            .graph
            .underlying
            .iter_edges_of(&self.cuts[cut].cut)
            .map(|(_, id, _)| id)
            .collect_vec();

        let left_amplitude_energy_sub = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .substitute_energies(&left_amplitude, &cut_edges);

        let right_amplitude_energy_sub = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .substitute_energies(&right_amplitude, &cut_edges);

        let left_ose_product =
            get_cff_inverse_energy_product_impl(&self.graph.underlying, &self.cuts[cut].left, &[]);
        let right_ose_product =
            get_cff_inverse_energy_product_impl(&self.graph.underlying, &self.cuts[cut].right, &[]);

        (
            left_amplitude_energy_sub * left_ose_product,
            right_amplitude_energy_sub * right_ose_product,
        )
    }

    pub fn build_left_right_amplitudes_for_orientation(
        &self,
        cut: CutId,
        orientation: SuperGraphOrientationID,
    ) -> (Atom, Atom) {
        let (left_amplitude, right_amplitude) = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .get_atom_for_orientation_and_cut(orientation, cut);

        let cut_edges = self
            .graph
            .underlying
            .iter_edges_of(&self.cuts[cut].cut)
            .map(|(_, id, _)| id)
            .collect_vec();

        let left_amplitude_energy_sub = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .substitute_energies(&left_amplitude, &cut_edges);

        let right_amplitude_energy_sub = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .substitute_energies(&right_amplitude, &cut_edges);

        let left_ose_product =
            get_cff_inverse_energy_product_impl(&self.graph.underlying, &self.cuts[cut].left, &[]);
        let right_ose_product =
            get_cff_inverse_energy_product_impl(&self.graph.underlying, &self.cuts[cut].right, &[]);

        (
            left_amplitude_energy_sub * left_ose_product,
            right_amplitude_energy_sub * right_ose_product,
        )
    }

    fn build_atom_for_cut(&self, cut_id: CutId) -> Atom {
        let (left_amplitude, right_amplitude) = self.build_left_right_amplitudes(cut_id);

        let product = left_amplitude * right_amplitude;

        let ose_atom = self.add_additional_factors_to_cff_atom(&product, cut_id);
        debug!("ose atom for cut: {}: {}", cut_id, ose_atom);

        let replacements = self.graph.underlying.get_ose_replacements();
        let replaced_atom = ose_atom.replace_multiple(&replacements);
        let replace_dots = replaced_atom
            .replace(function!(
                MS.dot,
                function!(GS.emr_vec, W_.x_),
                function!(GS.emr_vec, W_.y_)
            ))
            .with(
                -(function!(GS.emr_vec, W_.x_, 1) * function!(GS.emr_vec, W_.y_, 1)
                    + function!(GS.emr_vec, W_.x_, 2) * function!(GS.emr_vec, W_.y_, 2)
                    + function!(GS.emr_vec, W_.x_, 3) * function!(GS.emr_vec, W_.y_, 3)),
            );

        // debug!("replaced atom: {}", replace_dots);
        replace_dots
    }

    fn build_atom_for_orientation_and_cut(
        &self,
        cut_id: CutId,
        orientation: SuperGraphOrientationID,
    ) -> Atom {
        let (left_amplitude, right_amplitude) =
            self.build_left_right_amplitudes_for_orientation(cut_id, orientation);

        let product = left_amplitude * right_amplitude;

        let ose_atom = self.add_additional_factors_to_cff_atom(&product, cut_id);
        let replacements = self.graph.underlying.get_ose_replacements();
        let replaced_atom = ose_atom.replace_multiple(&replacements);
        let replace_dots = replaced_atom
            .replace(function!(
                MS.dot,
                function!(GS.emr_vec, W_.x_),
                function!(GS.emr_vec, W_.y_)
            ))
            .with(
                -(function!(GS.emr_vec, W_.x_, 1) * function!(GS.emr_vec, W_.y_, 1)
                    + function!(GS.emr_vec, W_.x_, 2) * function!(GS.emr_vec, W_.y_, 2)
                    + function!(GS.emr_vec, W_.x_, 3) * function!(GS.emr_vec, W_.y_, 3)),
            );

        debug!("replaced atom: {}", replace_dots);
        replace_dots
    }

    pub fn add_additional_factors_to_cff_atom(&self, cut_atom: &Atom, cut_id: CutId) -> Atom {
        let loop_3 = self.graph.underlying.get_loop_number() as i64 * 3;
        let t_star_factor = Atom::var(GS.rescale_star).npow(loop_3);

        let h_function = Atom::var(GS.hfunction);
        let grad_eta = Atom::var(GS.deta);

        let factors_of_pi = (Atom::num(2) * Atom::var(GS.pi)).npow(loop_3);

        let cut_inverse_energy_product = Atom::num(1)
            / self
                .graph
                .underlying
                .iter_edges_of(&self.cuts[cut_id].cut)
                .map(|(_, edge_id, _)| Atom::num(2) * ose_atom_from_index(edge_id))
                .fold(Atom::num(1), |product, factor| product * factor);

        let result = cut_atom
            * cut_inverse_energy_product
            * t_star_factor
            * h_function
            * &self.graph.multiplicity
            / grad_eta
            / factors_of_pi;

        //debug!("result: {}", result);
        result
    }

    fn get_params(&self, model: &Model) -> Vec<Atom> {
        let mut params = vec![];

        // all external energies
        params.extend(self.graph.underlying.get_external_energy_atoms());

        // spatial components of external momenta
        for (pair, edge_id, _) in self.graph.underlying.iter_edges() {
            match pair {
                HedgePair::Unpaired { .. } => {
                    let i64_id = Into::<usize>::into(edge_id) as i64;
                    let external_spatial = [
                        function!(GS.external_mom, i64_id, 1),
                        function!(GS.external_mom, i64_id, 2),
                        function!(GS.external_mom, i64_id, 3),
                    ];
                    params.extend(external_spatial);
                }
                _ => {}
            }
        }

        // spatial EMR
        for (pair, edge_id, _) in self.graph.underlying.iter_edges() {
            match pair {
                HedgePair::Paired { .. } => {
                    let i64_id = Into::<usize>::into(edge_id) as i64;
                    let emr_components = [
                        function!(GS.emr_vec, i64_id, 1),
                        function!(GS.emr_vec, i64_id, 2),
                        function!(GS.emr_vec, i64_id, 3),
                    ];
                    params.extend(emr_components)
                }
                _ => {}
            }
        }

        // add model parameters
        params.extend(model.generate_params());
        // add additional parameters
        params.push(Atom::var(GS.m_uv));
        params.push(Atom::var(GS.mu_r_sq));
        params.push(Atom::var(GS.rescale_star));
        params.push(Atom::var(GS.hfunction));
        params.push(Atom::var(GS.deta));

        params
    }

    fn get_function_map(&self) -> FunctionMap {
        let mut fn_map = FunctionMap::new();
        let pi_rational = Rational::from(std::f64::consts::PI);
        fn_map.add_constant(Atom::PI.into(), pi_rational.into());
        fn_map
    }

    pub fn build_cut_evaluators(
        &mut self,
        model: &Model,
        overwrite_atoms_for_test: Option<TiVec<CutId, Atom>>,
    ) {
        let evaluators = self
            .cuts
            .iter_enumerated()
            .map(|(cut_id, _)| {
                let atom = if let Some(atoms) = &overwrite_atoms_for_test {
                    atoms[cut_id].clone()
                } else {
                    self.build_atom_for_cut(cut_id)
                };

                let params = self.get_params(model);
                let mut tree = atom
                    .to_evaluation_tree(&self.get_function_map(), &params)
                    .unwrap();

                tree.horner_scheme();
                tree.common_subexpression_elimination();

                let tree_double = tree.map_coeff::<F<f64>, _>(&|r| (&r.re).into());
                let tree_quad = tree.map_coeff::<F<f128>, _>(&|r| (&r.re).into());

                let tree_complex_double = tree.map_coeff::<Complex<F<f64>>, _>(&|r| {
                    Complex::new(F(r.re.to_f64()), F(r.im.to_f64()))
                });
                let tree_complex_quad = tree.map_coeff::<Complex<F<f128>>, _>(&|r| {
                    Complex::new(F((&r.re).into()), F((&r.im).into()))
                });

                let lin = tree_double.linearize(None);
                let lin_complex = tree_complex_double.linearize(None);

                let filename = format!("{}_cut_{}.cpp", self.graph.name, cut_id);
                let function_name = format!("{}_cut_{}", self.graph.name, cut_id);
                let lib_name = format!("{}_cut_{}.so", self.graph.name, cut_id);
                let _exp = lin
                    .export_cpp(
                        &filename,
                        &function_name,
                        true,
                        symbolica::evaluate::InlineASM::None,
                    )
                    .unwrap()
                    .compile(&lib_name, CompileOptions::default())
                    .unwrap();

                GenericEvaluator {
                    f64_compiled: Some(RefCell::new(
                        SerializableCompiledEvaluator::load(&lib_name, &function_name).unwrap(),
                    )),
                    f64_eager: RefCell::new(lin_complex.into()),
                    f128: RefCell::new(tree_complex_quad.linearize(None)),
                }
            })
            .collect();

        self.derived_data.bare_cff_evaluators = Some(evaluators)
    }

    fn build_orientation_evaluators(&mut self, model: &Model) {
        let orientation_data = &self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .orientation_data;

        let substituted_energies = orientation_data
            .iter_enumerated()
            .map(|(orientation_id, data)| {
                data.cuts
                    .iter()
                    .map(|cut_id| self.build_atom_for_orientation_and_cut(*cut_id, orientation_id))
                    .collect_vec()
            })
            .collect::<TiVec<SuperGraphOrientationID, _>>();

        let orientation_evaluators = substituted_energies
            .iter()
            .zip(orientation_data)
            .map(|(cut_atoms, orientation_data)| {
                let cut_evaluators = cut_atoms
                    .iter()
                    .map(|cut_atom| {
                        let params = self.get_params(model);
                        let mut tree = cut_atom
                            .to_evaluation_tree(&self.get_function_map(), &params)
                            .unwrap();

                        tree.horner_scheme();
                        tree.common_subexpression_elimination();

                        let tree_double = tree.map_coeff::<Complex<F<f64>>, _>(&|r| {
                            Complex::new(F(r.re.to_f64()), F(r.im.to_f64()))
                        });
                        let tree_quad = tree.map_coeff::<Complex<F<f128>>, _>(&|r| {
                            Complex::new(F((&r.re).into()), F((&r.im).into()))
                        });

                        GenericEvaluator {
                            f64_compiled: None,
                            f64_eager: RefCell::new(tree_double.linearize(Some(1))),
                            f128: RefCell::new(tree_quad.linearize(Some(1))),
                        }
                    })
                    .collect_vec();

                OrientationEvaluator {
                    orientation_data: orientation_data.clone(),
                    evaluators: cut_evaluators,
                }
            })
            .collect();

        self.derived_data.bare_cff_orientation_evaluators = Some(orientation_evaluators);
    }

    fn build_lmbs(&mut self) {
        let lmbs = self
            .graph
            .underlying
            .generate_loop_momentum_bases(&self.graph.underlying.full_filter());

        self.derived_data.lmbs = Some(lmbs)
    }

    fn build_multi_channeling_channels(&mut self) {
        let channels = self
            .graph
            .build_multi_channeling_channels(self.derived_data.lmbs.as_ref().unwrap());

        self.derived_data.multi_channeling_setup = Some(channels)
    }

    fn generate_term_for_graph(&self, settings: &Settings) -> CrossSectionGraphTerm {
        let estimated_scale = self
            .graph
            .underlying
            .expected_scale(settings.kinematics.e_cm);

        CrossSectionGraphTerm {
            multi_channeling_setup: self.derived_data.multi_channeling_setup.clone().unwrap(),
            bare_cff_evaluators: self.derived_data.bare_cff_evaluators.clone().unwrap(),
            bare_cff_orientation_evaluators: self
                .derived_data
                .bare_cff_orientation_evaluators
                .clone()
                .unwrap_or_else(|| vec![].into()),
            graph: self.graph.clone(),
            cuts: self.cuts.clone(),
            cut_esurface: self.cut_esurface.clone(),
            lmbs: self.derived_data.lmbs.clone().unwrap(),
            estimated_scale,
        }
    }
}

#[derive(Clone, Encode, Decode)]
pub struct CrossSectionDerivedData<S: NumeratorState> {
    pub orientations: Option<TiVec<SuperGraphOrientationID, CutOrientationData>>,
    pub bare_cff_evaluators: Option<TiVec<CutId, GenericEvaluator>>,
    pub bare_cff_orientation_evaluators:
        Option<TiVec<SuperGraphOrientationID, OrientationEvaluator>>,
    pub cff_expression: Option<CFFCutsExpression>,
    pub lmbs: Option<TiVec<LmbIndex, LoopMomentumBasis>>,
    pub multi_channeling_setup: Option<LmbMultiChannelingSetup>,
    pub esurface_data: Option<EsurfaceDerivedData>,
    pub _temp_numerator: Option<PhantomData<S>>,
}

impl<S: NumeratorState> CrossSectionDerivedData<S> {
    fn new_empty() -> Self {
        Self {
            orientations: None,
            cff_expression: None,
            _temp_numerator: None,
            bare_cff_evaluators: None,
            bare_cff_orientation_evaluators: None,
            lmbs: None,
            multi_channeling_setup: None,
            esurface_data: None,
        }
    }
}

#[derive(Clone, bincode::Encode, bincode::Decode)]
pub struct CrossSectionCut {
    pub cut: OrientedCut,
    #[bincode(with_serde)]
    pub left: BitVec,
    #[bincode(with_serde)]
    pub right: BitVec,
}

impl CrossSectionCut {
    pub fn is_s_channel<S: NumeratorState>(
        &self,
        cross_section_graph: &CrossSectionGraph<S>,
    ) -> Result<bool> {
        let nodes_of_left_cut: Vec<_> = cross_section_graph
            .graph
            .underlying
            .iter_nodes_of(&self.left)
            .map(|(nid, _, _)| nid)
            .collect();

        let left_node = cross_section_graph
            .graph
            .underlying
            .combine_to_single_hedgenode(&nodes_of_left_cut);
        let res = left_node.includes(&cross_section_graph.source_nodes)
            && cross_section_graph.target_nodes.weakly_disjoint(&left_node);

        if !res {
            warn!("s channel check wrong");
        }

        Ok(true)
    }

    pub fn is_valid_for_process<S: NumeratorState>(
        &self,
        cross_section_graph: &CrossSectionGraph<S>,
        process: &ProcessDefinition,
        model: &Model,
    ) -> Result<bool> {
        if self.is_s_channel(cross_section_graph)? {
            let cut_content_builder = self
                .cut
                .iter_edges(&cross_section_graph.graph.underlying)
                .map(|(orientation, edge_data)| {
                    if orientation == Orientation::Reversed {
                        edge_data.data.particle.0.get_anti_particle(model)
                    } else {
                        edge_data.data.particle.clone()
                    }
                })
                .collect_vec();

            let any_pdg_list_passes = process
                .final_pdgs_lists
                .iter()
                .map(|x| {
                    x.iter()
                        .map(|pdg| model.get_particle_from_pdg(*pdg as isize))
                })
                .any(|particle_content| {
                    let mut cut_content = cut_content_builder.clone();
                    debug!(
                        "cut content: {:?}",
                        cut_content.iter().map(|p| p.0.pdg_code).collect_vec()
                    );

                    for particle in particle_content {
                        if let Some(index) = cut_content.iter().position(|p| p == &particle) {
                            cut_content.remove(index);
                        } else {
                            debug!("wrong particles");
                            return false;
                        }
                    }

                    if cut_content.len() > process.n_unresolved {
                        debug!(" too many unresolved particles");
                        return false;
                    }

                    if !cut_content
                        .iter()
                        .all(|particle| process.unresolved_cut_content.contains(particle))
                    {
                        debug!("wrong unresolved particles");
                        return false;
                    }

                    true
                });

            if !any_pdg_list_passes {
                debug!("wrong pdg list");
                return Ok(false);
            }

            let amplitude_couplings = process.amplitude_filters.get_coupling_orders();
            let amplitude_loop_count = process.amplitude_filters.get_loop_count_range();

            if let Some((min_loop, max_loop)) = amplitude_loop_count {
                let loop_range = min_loop..=max_loop;
                let left_internal_subgraph = InternalSubGraph::cleaned_filter_pessimist(
                    self.left.clone(),
                    &cross_section_graph.graph.underlying,
                );

                let right_internal_subgraph = InternalSubGraph::cleaned_filter_pessimist(
                    self.right.clone(),
                    &cross_section_graph.graph.underlying,
                );

                let left_loop = cross_section_graph
                    .graph
                    .underlying
                    .cyclotomatic_number(&left_internal_subgraph);

                let right_loop = cross_section_graph
                    .graph
                    .underlying
                    .cyclotomatic_number(&right_internal_subgraph);

                let total_loops = left_loop + right_loop;

                if !loop_range.contains(&total_loops) {
                    debug!("incorrect loop count");
                    return Ok(false);
                }
            }

            if amplitude_couplings.is_some() {
                todo!("waiting for update")
            }

            Ok(true)
        } else {
            debug!("cut is not s channel");
            Ok(false)
        }
    }
}

fn do_replacement_rules(
    mut orientation_expr: Atom,
    graph: &HedgeGraph<Edge, Vertex, NumHedgeData>,
) -> Atom {
    // add Feynman rules of external edges
    for (_p, edge_index, d) in graph.iter_edges_of(&graph.external_filter()) {
        let edge_id = usize::from(edge_index) as i64;
        orientation_expr = (orientation_expr * &d.data.spin_num)
            .replace(function!(GS.emr_mom, edge_id, W_.y_))
            .with_map(move |m| {
                let index = m.get(W_.y_).unwrap().to_atom();

                function!(GS.energy, edge_id, index) + function!(GS.emr_vec, edge_id, index)
            });
    }

    let spenso_mink = symbol!("spenso::mink");

    // substitute all OSEs from subgraphs, they are in the form OSE(edge_id, momentum, mass^2, mom.mom + mass^2)
    // the sqrt has already been applied
    // should simplify the expression
    orientation_expr = orientation_expr
        .replace(function!(GS.ose, W_.x_, W_.y_, W_.z_, W_.prop_, W_.a___))
        .with(function!(GS.ose, 100, W_.prop_, W_.a___));

    // contract all dot products, set all cross terms ose.q3 to 0
    // MS.dot is a 4d dot product
    orientation_expr = orientation_expr
        .replace(
            function!(GS.emr_vec, W_.x__, function!(spenso_mink, W_.z__))
                * function!(GS.ose, W_.y__, function!(spenso_mink, W_.z__)),
        )
        .with(Atom::Zero)
        .replace(
            function!(GS.emr_vec, W_.x__, function!(spenso_mink, W_.z__))
                * function!(GS.energy, W_.y__, function!(spenso_mink, W_.z__)),
        )
        .with(Atom::Zero)
        .replace(
            function!(GS.emr_vec, W_.x__, function!(spenso_mink, W_.z__))
                * function!(GS.ose, W_.y__, function!(spenso_mink, W_.z__)).pow(Atom::var(W_.b_)),
        )
        .with(Atom::Zero)
        .expand() // TODO: prevent expansion
        .replace(
            function!(GS.emr_vec, W_.x__, function!(spenso_mink, W_.z__))
                * function!(GS.ose, W_.y__, function!(spenso_mink, W_.z__)),
        )
        .with(Atom::Zero)
        .replace(
            function!(GS.emr_vec, W_.x__, function!(spenso_mink, W_.z__))
                * function!(GS.energy, W_.y__, function!(spenso_mink, W_.z__)),
        )
        .with(Atom::Zero)
        .replace(
            function!(GS.emr_vec, W_.x__, function!(spenso_mink, W_.z__))
                * function!(GS.ose, W_.y__, function!(spenso_mink, W_.z__)).pow(Atom::var(W_.b_)),
        )
        .with(Atom::Zero)
        .replace(function!(
            MS.dot,
            function!(GS.emr_vec, W_.x_),
            function!(GS.ose, W_.y_)
        ))
        .with(Atom::Zero)
        .replace(function!(
            MS.dot,
            function!(GS.emr_vec, W_.x_),
            function!(GS.energy, W_.y_)
        ))
        .with(Atom::Zero)
        .replace(function!(GS.emr_vec, W_.x__, W_.y_).npow(2))
        .with(function!(
            MS.dot,
            function!(GS.emr_vec, W_.x__),
            function!(GS.emr_vec, W_.x__)
        ))
        .replace(function!(GS.emr_vec, W_.x__, W_.a_) * function!(GS.emr_vec, W_.y__, W_.a_))
        .repeat()
        .with(function!(
            MS.dot,
            function!(GS.emr_vec, W_.x__),
            function!(GS.emr_vec, W_.y__)
        ))
        .replace(function!(GS.ose, W_.y__, function!(spenso_mink, W_.z__)).npow(2))
        .with(function!(GS.ose, W_.y__).npow(2))
        .replace(function!(GS.energy, W_.y__, function!(spenso_mink, W_.z__)).npow(2))
        .with(function!(GS.energy, W_.y__).npow(2))
        .replace(
            function!(GS.ose, W_.x__, function!(spenso_mink, W_.z__))
                * function!(GS.ose, W_.y__, function!(spenso_mink, W_.z__)),
        )
        .repeat()
        .with(function!(GS.ose, W_.x__) * function!(GS.ose, W_.y__))
        .replace(
            function!(GS.energy, W_.x__, function!(spenso_mink, W_.z__))
                * function!(GS.energy, W_.y__, function!(spenso_mink, W_.z__)),
        )
        .repeat()
        .with(function!(GS.energy, W_.x__) * function!(GS.energy, W_.y__))
        .replace(
            function!(GS.energy, W_.x__, function!(spenso_mink, W_.z__))
                * function!(GS.ose, W_.y__, function!(spenso_mink, W_.z__)),
        )
        .repeat()
        .with(function!(GS.energy, W_.x__) * function!(GS.ose, W_.y__))
        .replace(
            function!(GS.ose, W_.x__, function!(spenso_mink, W_.z__)).pow(Atom::var(W_.a_))
                * function!(GS.ose, W_.y__, function!(spenso_mink, W_.z__)),
        )
        .repeat()
        .with(function!(GS.ose, W_.x__).pow(Atom::var(W_.a_)) * function!(GS.ose, W_.y__))
        .replace(
            function!(GS.ose, W_.x__, function!(spenso_mink, W_.z__)).pow(Atom::var(W_.a_))
                * function!(GS.ose, W_.y__, function!(spenso_mink, W_.z__)).pow(Atom::var(W_.b_)),
        )
        .repeat()
        .with(
            function!(GS.ose, W_.x__).pow(Atom::var(W_.a_))
                * function!(GS.ose, W_.y__).pow(Atom::var(W_.b_)),
        )
        .replace(
            function!(GS.ose, W_.x__, function!(spenso_mink, W_.z__)).pow(Atom::var(W_.a_))
                * function!(GS.energy, W_.y__, function!(spenso_mink, W_.z__)),
        )
        .repeat()
        .with(function!(GS.ose, W_.x__).pow(Atom::var(W_.a_)) * function!(GS.energy, W_.y__))
        .replace(
            function!(GS.ose, W_.x__, function!(spenso_mink, W_.z__)).pow(Atom::var(W_.a_))
                * function!(GS.ose, W_.y__, function!(spenso_mink, W_.z__)).pow(Atom::var(W_.b_)),
        )
        .repeat()
        .with(
            function!(GS.ose, W_.x__).pow(Atom::var(W_.a_))
                * function!(GS.ose, W_.y__).pow(Atom::var(W_.b_)),
        )
        .replace(function!(
            MS.dot,
            function!(GS.ose, W_.x__),
            function!(GS.ose, W_.y__)
        ))
        .with(function!(GS.ose, W_.x__) * function!(GS.ose, W_.y__))
        .replace(function!(
            MS.dot,
            function!(GS.emr_mom, W_.x__),
            function!(GS.emr_vec, W_.y__)
        ))
        .with(function!(GS.emr_vec, W_.x__) * function!(GS.emr_vec, W_.y__)); // substitute all OSEs from subgraphs, they are in the form OSE(edge_id, momentum, mass^2, mom.mom + mass^2)
                                                                              // the sqrt has already been applied
    orientation_expr = orientation_expr
        .replace(function!(GS.ose, W_.x_, W_.y_, W_.z_, W_.prop_))
        .with(function!(GS.ose, 100, W_.prop_)) // do in two steps to get slightly nicer output
        .replace(function!(GS.ose, 100, W_.prop_))
        .with(Atom::var(W_.prop_).sqrt().npow(2))
        .replace(function!(GS.ose, 100, W_.prop_, W_.x_))
        .with(Atom::var(W_.prop_).sqrt().npow(2)); // it could be that GS.ose(mu)^1/2 fused into GS.ose(mu)^1 which leaves a fake dummy index

    // simplify nested exponents
    orientation_expr = orientation_expr
        .replace(Atom::var(W_.x_).pow(Atom::var(W_.a_)).pow(Atom::var(W_.b_)))
        .repeat()
        .with(Atom::var(W_.x_).pow(Atom::var(W_.a_) * Atom::var(W_.b_)));

    orientation_expr = orientation_expr
        .replace(function!(GS.external_mom, W_.x_, W_.y_))
        .with(function!(GS.energy, W_.x_));

    // set the external energies
    for (_p, edge_index, _d) in graph.iter_edges_of(&graph.external_filter()) {
        let edge_id = usize::from(edge_index) as i64;
        orientation_expr = orientation_expr
            .replace(function!(GS.energy, edge_id))
            .with(external_energy_atom_from_index(edge_index));
    }

    orientation_expr = orientation_expr
        .replace(function!(GS.emr_vec, W_.x_ + W_.y__))
        .repeat()
        .with(function!(GS.emr_vec, W_.x_) + function!(GS.emr_vec, W_.y__));
    orientation_expr = orientation_expr
        .replace(function!(GS.emr_vec, -function!(GS.emr_mom, W_.x_)))
        .with(-function!(GS.emr_vec, W_.x_));

    orientation_expr = orientation_expr
        .replace(function!(GS.emr_vec, function!(GS.emr_mom, W_.x_)))
        .with(function!(GS.emr_vec, W_.x_));

    orientation_expr = orientation_expr.replace(parse!("ZERO")).with(Atom::new());

    let replacements = graph.get_ose_replacements();
    orientation_expr = orientation_expr.replace_multiple(&replacements);

    orientation_expr = orientation_expr
        .replace(function!(
            MS.dot,
            function!(GS.emr_vec, W_.x_),
            function!(GS.emr_vec, W_.y_)
        ))
        .with(
            -(function!(GS.emr_vec, W_.x_, 1) * function!(GS.emr_vec, W_.y_, 1)
                + function!(GS.emr_vec, W_.x_, 2) * function!(GS.emr_vec, W_.y_, 2)
                + function!(GS.emr_vec, W_.x_, 3) * function!(GS.emr_vec, W_.y_, 3)),
        )
        .replace(parse!("ZERO"))
        .with(Atom::new());

    // set the external spatial parts
    for (_p, edge_index, _d) in graph.iter_edges_of(&graph.external_filter()) {
        let edge_id = usize::from(edge_index) as i64;

        orientation_expr = orientation_expr
            .replace(function!(GS.emr_vec, edge_id, W_.x_))
            .with(function!(GS.external_mom, edge_id, W_.x_));
    }

    orientation_expr
}

#[cfg(test)]
mod tests {
    use std::fs::OpenOptions;

    use linnet::half_edge::involution::EdgeIndex;
    use spenso::network::library::TensorLibraryData;
    use symbolica::{atom::Atom, state::State};

    use crate::{
        dot,
        new_graph::{Graph, LoopMomentumBasis},
        numerator::{
            EvaluatorOptions, GammaAlgebraMode, GlobalPrefactor, NumeratorEvaluatorOptions,
            NumeratorParseMode, NumeratorSettings, UnInit,
        },
        signature::LoopExtSignature,
        tests_from_pytest::load_generic_model,
        GammaLoopContextContainer, GammaloopCompileOptions, TropicalSubgraphTableSettings,
    };

    use super::AmplitudeGraph;

    #[test]
    fn test_encode_decode_amplitude_graph() {
        // load the model and hack the masses, go through serializable model since arc is not mutable
        let model = load_generic_model("sm");

        let mut graphs = dot!(
            digraph G{
                e1      [flow=sink]
                e2      [flow=source]
                e3      [flow=source]
                e1 -> n1  [particle=h]
                e2 -> n4    [particle=h]
                n1 -> n2    [particle=h]
                n1 -> n3    [particle=h]
                n2 -> n3    [particle=t]
                n3 -> n4    [particle=t]
                n4 -> n2    [particle=t]
            }
        )
        .unwrap();
        let graph = &mut graphs[0];
        let mut loop_momentum_basis = LoopMomentumBasis {
            tree: None,
            loop_edges: vec![EdgeIndex::from(0), EdgeIndex::from(4)].into(),
            ext_edges: vec![EdgeIndex::from(5), EdgeIndex::from(6)].into(),
            edge_signatures: graph
                .underlying
                .new_edgevec(|_, _, _| LoopExtSignature::from((vec![], vec![]))),
        };

        // loop_momentum_basis
        //     .set_edge_signatures(&graph.underlying)
        //     .unwrap();

        graph.loop_momentum_basis = loop_momentum_basis;

        let mut amplitude: AmplitudeGraph<UnInit> = AmplitudeGraph::new(graph.clone());

        amplitude
            .preprocess(
                &model,
                &crate::ProcessSettings {
                    compile_cff: false,
                    compile_separate_orientations: false,
                    cpe_rounds_cff: None,
                    numerator_settings: NumeratorSettings {
                        eval_settings: NumeratorEvaluatorOptions::Single(EvaluatorOptions {
                            cpe_rounds: None,
                            compile_options: crate::numerator::NumeratorCompileOptions::NotCompiled,
                        }),
                        parse_mode: NumeratorParseMode::Direct,
                        dump_expression: None,
                        global_numerator: None,
                        global_prefactor: GlobalPrefactor {
                            color: Atom::one(),
                            colorless: Atom::one(),
                        },
                        gamma_algebra: GammaAlgebraMode::Concrete,
                    },
                    gammaloop_compile_options: GammaloopCompileOptions {
                        inline_asm: false,
                        fast_math: false,
                        optimization_level: 0,
                        unsafe_math: false,
                        compiler: "g++".into(),
                        custom: Vec::new(),
                    },
                    tropical_subgraph_table_settings: TropicalSubgraphTableSettings {
                        panic_on_fail: false,
                        target_omega: 1.0,
                    },
                },
            )
            .unwrap();

        let mut temp = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open("test.bin")
            .unwrap();

        State::export(&mut temp).unwrap();
        drop(temp);

        let mut temp = OpenOptions::new().read(true).open("test.bin").unwrap();
        let state_map = State::import(&mut temp, None).unwrap();

        let context = GammaLoopContextContainer {
            model: &model,
            state_map: &state_map,
        };

        println!("context created");

        let encoded_amplitude =
            bincode::encode_to_vec(&amplitude, bincode::config::standard()).unwrap();

        let _amplitude: AmplitudeGraph<UnInit> = bincode::decode_from_slice_with_context(
            &encoded_amplitude,
            bincode::config::standard(),
            context,
        )
        .expect("amplitude decode failed")
        .0;

        println!("amplitude graph passed");
    }
}
