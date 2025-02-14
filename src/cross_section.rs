use crate::gammaloop_integrand::GammaLoopIntegrand;
use crate::graph::{BareGraph, Graph, SerializableGraph};
use crate::model::{Model, Particle};
use crate::momentum::Signature;
use crate::numerator::{
    AppliedFeynmanRule, ContractionSettings, Evaluators, GetSingleAtom, Numerator, NumeratorState,
    PythonState, TypedNumeratorState, UnInit,
};
use crate::{utils::*, ExportSettings, Externals, Polarizations, Settings};
use bincode;
use color_eyre::Result;
use color_eyre::{Help, Report};
#[allow(unused_imports)]
use eyre::{eyre, Context};
use log::debug;
use serde::{Deserialize, Serialize};
use serde_yaml::Error;
use smartstring::{LazyCompact, SmartString};
use std::fs;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;
use symbolica::printer::{AtomPrinter, PrintOptions};
use symbolica::state::State;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum OutputType {
    #[serde(rename = "amplitudes")]
    Amplitudes,
    #[serde(rename = "cross_sections")]
    CrossSections,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct OutputMetaData {
    pub model_name: String,
    pub output_type: OutputType,
    pub contents: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SerializableSuperGraphCut {
    pub cut_edges: Vec<SmartString<LazyCompact>>,
    pub forward_scattering_graph: SerializableForwardScatteringGraph,
}

impl SerializableSuperGraphCut {
    pub fn from_supergraph_cut(
        graph: &BareGraph,
        supergraph_cut: &SuperGraphCut,
    ) -> SerializableSuperGraphCut {
        SerializableSuperGraphCut {
            cut_edges: supergraph_cut
                .cut_edges
                .iter()
                .map(|&edge| graph.edges[edge].name.clone())
                .collect(),
            forward_scattering_graph:
                SerializableForwardScatteringGraph::from_forward_scattering_graph(
                    &supergraph_cut.forward_scattering_graph,
                ),
        }
    }
}

#[derive(Debug, Clone)]
pub struct SuperGraphCut {
    pub cut_edges: Vec<usize>,
    pub forward_scattering_graph: ForwardScatteringGraph,
}

impl SuperGraphCut {
    pub fn from_serializable_supergraph_cut(
        model: &Model,
        graph: &BareGraph,
        serializable_supergraph_cut: &SerializableSuperGraphCut,
    ) -> SuperGraphCut {
        SuperGraphCut {
            cut_edges: serializable_supergraph_cut
                .cut_edges
                .iter()
                .map(|edge_name| {
                    graph
                        .get_edge_position(edge_name)
                        .unwrap_or_else(|| panic!("Could not find edge with name {}", edge_name))
                })
                .collect(),
            forward_scattering_graph:
                ForwardScatteringGraph::from_serializable_forward_scattering_graph(
                    model,
                    &serializable_supergraph_cut.forward_scattering_graph,
                ),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SerializableSuperGraph {
    pub sg_id: usize,
    pub graph: SerializableGraph,
    pub multiplicity: String,
    // This identifier of the topology class is mostly a stub for now
    pub topology_class: Vec<usize>,
    pub cuts: Vec<SerializableSuperGraphCut>,
}

impl SerializableSuperGraph {
    pub fn from_supergraph(supergraph: &SuperGraph) -> SerializableSuperGraph {
        SerializableSuperGraph {
            sg_id: supergraph.sg_id,
            graph: SerializableGraph::from_graph(&supergraph.graph.bare_graph),
            multiplicity: supergraph.multiplicity.clone(),
            topology_class: supergraph.topology_class.clone(),
            cuts: supergraph
                .cuts
                .iter()
                .map(|cut| {
                    SerializableSuperGraphCut::from_supergraph_cut(
                        &supergraph.graph.bare_graph,
                        cut,
                    )
                })
                .collect(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct SuperGraph {
    pub sg_id: usize,
    pub graph: Graph<Evaluators>,
    pub multiplicity: String,
    // This identifier of the topology class is mostly a stub for now
    pub topology_class: Vec<usize>,
    pub cuts: Vec<SuperGraphCut>,
}

impl SuperGraph {
    pub fn from_serializable_supergraph(
        model: &Model,
        serializable_supergraph: &SerializableSuperGraph,
    ) -> SuperGraph {
        let _g = Graph::from_serializable_graph(model, &serializable_supergraph.graph);
        // let cuts = serializable_supergraph
        //     .cuts
        //     .iter()
        //     .map(|cut| SuperGraphCut::from_serializable_supergraph_cut(model, &g.bare_graph, cut))
        //     .collect();
        // SuperGraph {
        //     sg_id: serializable_supergraph.sg_id,
        //     graph: g,
        //     multiplicity: serializable_supergraph.multiplicity,
        //     topology_class: serializable_supergraph.topology_class.clone(),
        //     cuts,
        // }
        unimplemented!()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SerializableForwardScatteringGraphCut {
    pub cut_edges: Vec<SmartString<LazyCompact>>,
    pub amplitudes: [SerializableAmplitude; 2],
}

impl SerializableForwardScatteringGraphCut {
    pub fn from_forward_scattering_graph_cut(
        graph: &BareGraph,
        forward_scattering_graph_cut: &ForwardScatteringGraphCut,
    ) -> SerializableForwardScatteringGraphCut {
        SerializableForwardScatteringGraphCut {
            cut_edges: forward_scattering_graph_cut
                .cut_edges
                .iter()
                .map(|&edge| graph.edges[edge].name.clone())
                .collect(),
            amplitudes: [
                SerializableAmplitude::from_amplitude(&forward_scattering_graph_cut.amplitudes[0]),
                SerializableAmplitude::from_amplitude(&forward_scattering_graph_cut.amplitudes[1]),
            ],
        }
    }
}
#[derive(Debug, Clone)]
pub struct ForwardScatteringGraphCut {
    pub cut_edges: Vec<usize>,
    pub amplitudes: [Amplitude<UnInit>; 2],
}

impl ForwardScatteringGraphCut {
    pub fn from_serializable_forward_scattering_graph_cut(
        model: &Model,
        graph: &BareGraph,
        serializable_forward_scattering_graph_cut: &SerializableForwardScatteringGraphCut,
    ) -> ForwardScatteringGraphCut {
        ForwardScatteringGraphCut {
            cut_edges: serializable_forward_scattering_graph_cut
                .cut_edges
                .iter()
                .map(|edge_name| {
                    graph
                        .get_edge_position(edge_name)
                        .unwrap_or_else(|| panic!("Could not find edge with name {}", edge_name))
                })
                .collect(),
            amplitudes: [
                Amplitude::from_serializable_amplitude(
                    model,
                    &serializable_forward_scattering_graph_cut.amplitudes[0],
                ),
                Amplitude::from_serializable_amplitude(
                    model,
                    &serializable_forward_scattering_graph_cut.amplitudes[1],
                ),
            ],
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SerializableForwardScatteringGraph {
    pub sg_id: usize,
    pub sg_cut_id: usize,
    pub graph: SerializableGraph,
    pub multiplicity: String,
    pub cuts: Vec<SerializableForwardScatteringGraphCut>,
}

impl SerializableForwardScatteringGraph {
    pub fn from_forward_scattering_graph(
        forward_scattering_graph: &ForwardScatteringGraph,
    ) -> SerializableForwardScatteringGraph {
        SerializableForwardScatteringGraph {
            sg_id: forward_scattering_graph.sg_id,
            sg_cut_id: forward_scattering_graph.sg_cut_id,
            graph: SerializableGraph::from_graph(&forward_scattering_graph.graph),
            multiplicity: forward_scattering_graph.multiplicity.clone(),
            cuts: forward_scattering_graph
                .cuts
                .iter()
                .map(|cut| {
                    SerializableForwardScatteringGraphCut::from_forward_scattering_graph_cut(
                        &forward_scattering_graph.graph,
                        cut,
                    )
                })
                .collect(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct ForwardScatteringGraph {
    pub sg_id: usize,
    pub sg_cut_id: usize,
    pub graph: BareGraph,
    pub multiplicity: String,
    pub cuts: Vec<ForwardScatteringGraphCut>,
}

impl ForwardScatteringGraph {
    pub fn from_serializable_forward_scattering_graph(
        model: &Model,
        forward_scattering_graph: &SerializableForwardScatteringGraph,
    ) -> ForwardScatteringGraph {
        let g = BareGraph::from_serializable_graph(model, &forward_scattering_graph.graph);
        let cuts = forward_scattering_graph
            .cuts
            .iter()
            .map(|cut| {
                ForwardScatteringGraphCut::from_serializable_forward_scattering_graph_cut(
                    model, &g, cut,
                )
            })
            .collect();
        ForwardScatteringGraph {
            sg_id: forward_scattering_graph.sg_id,
            sg_cut_id: forward_scattering_graph.sg_cut_id,
            graph: g,
            multiplicity: forward_scattering_graph.multiplicity.clone(),
            cuts,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SerializableAmplitudeGraph {
    pub sg_id: usize,
    pub sg_cut_id: usize,
    pub fs_cut_id: usize,
    pub amplitude_side: Side,
    pub graph: SerializableGraph,
    pub multiplicity: String,
    pub multi_channeling_channels: Vec<usize>, // empty list defaults to all channels if multi_channeling is enabled
}

impl SerializableAmplitudeGraph {
    pub fn from_amplitude_graph<S: NumeratorState>(
        amplitude_graph: &AmplitudeGraph<S>,
    ) -> SerializableAmplitudeGraph {
        SerializableAmplitudeGraph {
            sg_id: amplitude_graph.sg_id,
            sg_cut_id: amplitude_graph.sg_cut_id,
            fs_cut_id: amplitude_graph.fs_cut_id,
            amplitude_side: amplitude_graph.amplitude_side.clone(),
            graph: SerializableGraph::from_graph(&amplitude_graph.graph.bare_graph),
            multiplicity: amplitude_graph.multiplicity.clone(),
            multi_channeling_channels: amplitude_graph.multi_channeling_channels.clone(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct AmplitudeGraph<NumState: NumeratorState> {
    pub sg_id: usize,
    pub sg_cut_id: usize,
    pub fs_cut_id: usize,
    pub amplitude_side: Side,
    pub graph: Graph<NumState>,
    pub multiplicity: String,
    pub multi_channeling_channels: Vec<usize>,
}

impl<S: TypedNumeratorState> AmplitudeGraph<S> {
    pub fn try_from_python(ag: AmplitudeGraph<PythonState>) -> Result<Self> {
        ag.map_res(Graph::<S>::try_from_python)
    }
}

impl AmplitudeGraph<PythonState> {
    pub fn sync(&mut self, _model: &Model) {
        // self.graph.sync();
    }
}

impl<S: NumeratorState> AmplitudeGraph<S> {
    pub fn load_derived_data_mut(
        &mut self,
        model: &Model,
        path: &Path,
        settings: &Settings,
    ) -> Result<()> {
        let g = self
            .graph
            .bare_graph
            .clone()
            .load_derived_data::<S>(model, path, settings)?;

        self.graph = g;
        Ok(())
        //  load_derived_data::<S>(path)
    }

    pub fn load_derived_data<T: NumeratorState>(
        self,
        model: &Model,
        path: &Path,
        settings: &Settings,
    ) -> Result<AmplitudeGraph<T>> {
        let g = self
            .graph
            .bare_graph
            .load_derived_data::<T>(model, path, settings)?;

        Ok(AmplitudeGraph {
            sg_id: self.sg_id,
            sg_cut_id: self.sg_cut_id,
            fs_cut_id: self.fs_cut_id,
            amplitude_side: self.amplitude_side,
            graph: g,
            multiplicity: self.multiplicity,
            multi_channeling_channels: self.multi_channeling_channels,
        })
        //  load_derived_data::<S>(path)
    }

    pub fn map<F, U: NumeratorState>(self, mut f: F) -> AmplitudeGraph<U>
    where
        F: FnMut(Graph<S>) -> Graph<U>,
    {
        AmplitudeGraph {
            sg_id: self.sg_id,
            sg_cut_id: self.sg_cut_id,
            fs_cut_id: self.fs_cut_id,
            amplitude_side: self.amplitude_side,
            graph: f(self.graph),
            multiplicity: self.multiplicity,
            multi_channeling_channels: self.multi_channeling_channels,
        }
    }

    pub fn map_mut<F>(&mut self, mut f: F)
    where
        F: FnMut(&mut Graph<S>),
    {
        f(&mut self.graph);
    }

    pub fn map_res<F, U: NumeratorState, E>(self, mut f: F) -> Result<AmplitudeGraph<U>, E>
    where
        F: FnMut(Graph<S>) -> Result<Graph<U>, E>,
    {
        Ok(AmplitudeGraph {
            sg_id: self.sg_id,
            sg_cut_id: self.sg_cut_id,
            fs_cut_id: self.fs_cut_id,
            amplitude_side: self.amplitude_side,
            graph: f(self.graph)?,
            multiplicity: self.multiplicity,
            multi_channeling_channels: self.multi_channeling_channels,
        })
    }
}

impl AmplitudeGraph<UnInit> {
    pub fn from_amplitude_graph(
        model: &Model,
        amplitude_graph: &SerializableAmplitudeGraph,
    ) -> Self {
        Self {
            sg_id: amplitude_graph.sg_id,
            sg_cut_id: amplitude_graph.sg_cut_id,
            fs_cut_id: amplitude_graph.fs_cut_id,
            amplitude_side: amplitude_graph.amplitude_side.clone(),
            graph: Graph::from_serializable_graph(model, &amplitude_graph.graph),
            multiplicity: amplitude_graph.multiplicity.clone(),
            multi_channeling_channels: amplitude_graph.multi_channeling_channels.clone(),
        }
    }

    pub fn apply_feynman_rules(
        self,
        export_settings: &ExportSettings,
    ) -> AmplitudeGraph<AppliedFeynmanRule> {
        let graph = self.graph.apply_feynman_rules(export_settings);
        AmplitudeGraph {
            sg_id: self.sg_id,
            sg_cut_id: self.sg_cut_id,
            fs_cut_id: self.fs_cut_id,
            amplitude_side: self.amplitude_side,
            graph,
            multiplicity: self.multiplicity.clone(),
            multi_channeling_channels: self.multi_channeling_channels,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SerializableCrossSection {
    name: SmartString<LazyCompact>,
    supergraphs: Vec<SerializableSuperGraph>,
}

impl SerializableCrossSection {
    pub fn from_file(file_path: String) -> Result<SerializableCrossSection, Report> {
        let f = File::open(file_path.clone())
            .wrap_err_with(|| format!("Could not open cross-section yaml file {}", file_path))
            .suggestion("Does the path exist?")?;
        serde_yaml::from_reader(f)
            .map_err(|e| eyre!(format!("Error parsing cross-section yaml: {}", e)))
            .suggestion("Is it a correct yaml file")
    }

    pub fn from_cross_section(cross_section: &CrossSection) -> SerializableCrossSection {
        SerializableCrossSection {
            name: cross_section.name.clone(),
            supergraphs: cross_section
                .supergraphs
                .iter()
                .map(SerializableSuperGraph::from_supergraph)
                .collect(),
        }
    }

    pub fn from_yaml_str(yaml_str: String) -> Result<SerializableCrossSection, Report> {
        serde_yaml::from_str(yaml_str.as_str())
            .map_err(|e| eyre!(format!("Error parsing cross-section yaml: {}", e)))
            .suggestion("Is it a correct yaml file")
    }
}

#[derive(Debug, Clone)]
pub struct CrossSection {
    pub name: SmartString<LazyCompact>,
    pub supergraphs: Vec<SuperGraph>,
}

impl CrossSection {
    pub fn from_file(model: &Model, file_path: String) -> Result<CrossSection, Report> {
        SerializableCrossSection::from_file(file_path).map(|serializable_cross_section| {
            CrossSection::from_serializable_cross_section(model, &serializable_cross_section)
        })
    }

    pub fn from_yaml_str(model: &Model, yaml_str: String) -> Result<CrossSection, Report> {
        SerializableCrossSection::from_yaml_str(yaml_str).map(|serializable_cross_section| {
            CrossSection::from_serializable_cross_section(model, &serializable_cross_section)
        })
    }

    pub fn from_serializable_cross_section(
        model: &Model,
        serializable_cross_section: &SerializableCrossSection,
    ) -> CrossSection {
        CrossSection {
            name: serializable_cross_section.name.clone(),
            supergraphs: serializable_cross_section
                .supergraphs
                .iter()
                .map(|sg| SuperGraph::from_serializable_supergraph(model, sg))
                .collect(),
        }
    }

    pub fn to_serializable(&self) -> SerializableCrossSection {
        SerializableCrossSection::from_cross_section(self)
    }

    #[allow(unused)]
    pub fn export(&self, export_root: &str, model: &Model) -> Result<(), Report> {
        // TODO process cross-section by adding lots of additional information necessary for runtime.
        // e.g. generate e-surface, cff expression, counterterms, etc.

        // Then dumped the new yaml representation of the amplitude now containing all that additional information
        let path = Path::new(export_root)
            .join("sources")
            .join("cross_sections")
            .join(self.name.as_str());

        fs::write(
            path.join("cross_section.yaml"),
            serde_yaml::to_string(&self.to_serializable())?,
        )?;

        // Additional files can be written too, e.g. the lengthy cff expressions can be dumped in separate files
        fs::write(path.join("cff_expression.yaml"), "TODO")?;

        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SerializableAmplitude {
    name: SmartString<LazyCompact>,
    amplitude_graphs: Vec<SerializableAmplitudeGraph>,
}

impl SerializableAmplitude {
    pub fn from_amplitude<S: NumeratorState>(amplitude: &Amplitude<S>) -> SerializableAmplitude {
        SerializableAmplitude {
            name: amplitude.name.clone(),
            amplitude_graphs: amplitude
                .amplitude_graphs
                .iter()
                .map(SerializableAmplitudeGraph::from_amplitude_graph)
                .collect(),
        }
    }

    pub fn from_file(file_path: String) -> Result<SerializableAmplitude, Report> {
        let f = File::open(file_path.clone())
            .wrap_err_with(|| format!("Could not open amplitude yaml file {}", file_path))
            .suggestion("Does the path exist?")?;
        serde_yaml::from_reader(f)
            .map_err(|e| eyre!(format!("Error parsing amplitude yaml: {}", e)))
            .suggestion("Is it a correct yaml file")
    }

    pub fn from_yaml_str(yaml_str: String) -> Result<SerializableAmplitude, Report> {
        serde_yaml::from_str(yaml_str.as_str())
            .map_err(|e| eyre!(format!("Error parsing amplitude yaml: {}", e)))
            .suggestion("Is it a correct yaml file")
    }
}

pub trait IsPolarizable {
    fn polarizations(&self, externals: &Externals) -> Polarizations;
}

#[derive(Debug, Clone)]
pub struct Amplitude<NumState: NumeratorState = Evaluators> {
    pub name: SmartString<LazyCompact>,
    pub external_signature: Signature,
    pub external_particles: Vec<Arc<Particle>>,
    pub amplitude_graphs: Vec<AmplitudeGraph<NumState>>,
}

impl<S: TypedNumeratorState> Amplitude<S> {
    pub fn try_from_python(amp: Amplitude<PythonState>) -> Result<Self> {
        amp.map_res(AmplitudeGraph::<S>::try_from_python)
    }
}

impl Amplitude<PythonState> {
    pub fn sync(&mut self, model: &Model) {
        self.amplitude_graphs
            .iter_mut()
            .for_each(|ag| ag.sync(model));
    }
}

impl<S: NumeratorState> Amplitude<S> {
    pub fn external_signature(&self) -> Signature {
        let external_signature = self
            .amplitude_graphs
            .first()
            .unwrap()
            .graph
            .bare_graph
            .external_in_or_out_signature();
        for amplitude_graph in self.amplitude_graphs.iter() {
            assert_eq!(
                amplitude_graph
                    .graph
                    .bare_graph
                    .external_in_or_out_signature(),
                external_signature
            );
        }

        external_signature
    }

    pub fn external_particle_spin_and_masslessness(&self) -> Vec<(isize, bool)> {
        self.amplitude_graphs
            .first()
            .unwrap()
            .graph
            .bare_graph
            .external_particle_spin_and_masslessness()
    }

    pub fn load_derived_data_mut(
        &mut self,
        model: &Model,
        path: &Path,
        settings: &Settings,
    ) -> Result<()> {
        for amplitude_graph in self.amplitude_graphs.iter_mut() {
            amplitude_graph.load_derived_data_mut(model, path, settings)?;
        }
        Ok(())
    }

    pub fn map<F, U: NumeratorState>(self, f: F) -> Amplitude<U>
    where
        F: FnMut(AmplitudeGraph<S>) -> AmplitudeGraph<U>,
    {
        Amplitude {
            name: self.name,
            external_particles: self.external_particles,
            external_signature: self.external_signature,
            amplitude_graphs: self.amplitude_graphs.into_iter().map(f).collect(),
        }
    }

    pub fn map_mut<F>(&mut self, f: F)
    where
        F: FnMut(&mut AmplitudeGraph<S>),
    {
        self.amplitude_graphs.iter_mut().for_each(f);
    }

    pub fn map_res<F, U: NumeratorState, E>(self, f: F) -> Result<Amplitude<U>, E>
    where
        F: FnMut(AmplitudeGraph<S>) -> Result<AmplitudeGraph<U>, E>,
    {
        let new_amp_graphs: Result<_, E> = self.amplitude_graphs.into_iter().map(f).collect();
        Ok(Amplitude {
            external_particles: self.external_particles,
            external_signature: self.external_signature,
            name: self.name,
            amplitude_graphs: new_amp_graphs?,
        })
    }
}

impl<S: NumeratorState> IsPolarizable for Amplitude<S> {
    fn polarizations(&self, externals: &Externals) -> Polarizations {
        externals.generate_polarizations(&self.external_particles, &self.external_signature)
    }
}

impl Amplitude<UnInit> {
    pub fn from_file(model: &Model, file_path: String) -> Result<Self, Report> {
        SerializableAmplitude::from_file(file_path).map(|serializable_amplitude| {
            Amplitude::from_serializable_amplitude(model, &serializable_amplitude)
        })
    }

    pub fn from_yaml_str(model: &Model, yaml_str: String) -> Result<Self, Report> {
        SerializableAmplitude::from_yaml_str(yaml_str).map(|serializable_amplitude| {
            Amplitude::from_serializable_amplitude(model, &serializable_amplitude)
        })
    }

    pub fn from_serializable_amplitude(
        model: &Model,
        serializable_amplitude: &SerializableAmplitude,
    ) -> Self {
        let amplitude_graphs: Vec<_> = serializable_amplitude
            .amplitude_graphs
            .iter()
            .map(|sg| AmplitudeGraph::from_amplitude_graph(model, sg))
            .collect();

        let external_signature = amplitude_graphs
            .first()
            .unwrap()
            .graph
            .bare_graph
            .external_in_or_out_signature();

        let external_particles = amplitude_graphs
            .first()
            .unwrap()
            .graph
            .bare_graph
            .external_particles();

        for ag in &amplitude_graphs {
            let other_signature = ag.graph.bare_graph.external_in_or_out_signature();

            if external_signature != other_signature {
                panic!(
                    "External signature mismatch: {:?} != {:?}",
                    external_signature, other_signature
                );
            }

            let other_particles = ag.graph.bare_graph.external_particles();

            if other_particles != external_particles {
                panic!(
                    "External particles mismatch: {:?} != {:?}",
                    external_particles, other_particles
                );
            }
        }

        Self {
            name: serializable_amplitude.name.clone(),
            amplitude_graphs,
            external_particles,
            external_signature,
        }
    }

    pub fn apply_feynman_rules(
        self,
        export_settings: &ExportSettings,
    ) -> Amplitude<AppliedFeynmanRule> {
        let graphs = self
            .amplitude_graphs
            .into_iter()
            .map(|ag| ag.apply_feynman_rules(export_settings))
            .collect();

        Amplitude {
            external_particles: self.external_particles,
            external_signature: self.external_signature,
            name: self.name,
            amplitude_graphs: graphs,
        }
    }

    pub fn load_derived_data<S: NumeratorState>(
        self,
        model: &Model,
        path: &Path,
        settings: &Settings,
    ) -> Result<Amplitude<S>, Report> {
        self.map_res(|a| a.load_derived_data(model, path, settings))
    }
}

impl Amplitude<PythonState> {
    pub fn export(
        &mut self,
        export_root: &str,
        model: &Model,
        export_settings: &ExportSettings,
        no_evaluators: bool,
    ) -> Result<(), Report> {
        // TODO process amplitude by adding lots of additional information necessary for runtime.
        // e.g. generate e-surface, cff expression, counterterms, etc.

        // Then dumped the new yaml representation of the amplitude now containing all that additional information
        let path = Path::new(export_root)
            .join("sources")
            .join("amplitudes")
            .join(self.name.as_str());

        // generate cff and ltd for each graph in the ampltiudes, ltd also generates lmbs

        self.map_mut(|a| {
            a.map_mut(|g| {
                g.generate_cff();
                if !g.bare_graph.is_tree() {
                    g.generate_ltd();
                } else {
                    g.generate_loop_momentum_bases();
                }
                g.generate_tropical_subgraph_table(
                    &export_settings.tropical_subgraph_table_settings,
                );
                g.generate_esurface_data().unwrap();
                g.build_compiled_expression(path.clone(), export_settings)
                    .unwrap();

                if no_evaluators {
                    g.forgetfull_apply::<_, UnInit>(|d, b| {
                        d.process_numerator_no_eval(
                            b,
                            ContractionSettings::Normal,
                            path.clone(),
                            export_settings,
                        )
                        .unwrap()
                    })
                } else {
                    g.statefull_apply::<_, UnInit, Evaluators>(|d, b| {
                        d.process_numerator(
                            b,
                            model,
                            ContractionSettings::Normal,
                            path.clone(),
                            export_settings,
                        )
                        .unwrap()
                    })
                }
                .unwrap();
            })
        });

        fs::write(
            path.clone().join("amplitude.yaml"),
            serde_yaml::to_string(&self.to_serializable())?,
        )?;

        let mut state = fs::File::create(path.clone().join("state.bin"))?;
        State::export(&mut state)?;

        // dump the derived data in a binary file
        for amplitude_graph in self.amplitude_graphs.iter() {
            debug!("dumping derived data");
            fs::write(
                path.clone().join(format!(
                    "derived_data_{}.bin",
                    amplitude_graph.graph.bare_graph.name
                )),
                &bincode::encode_to_vec(
                    amplitude_graph
                        .graph
                        .derived_data
                        .as_ref()
                        .ok_or(eyre!("Empty derived data"))?,
                    bincode::config::standard(),
                )?,
            )?;
        }

        // Additional files can be written too, e.g. the lengthy cff expressions can be dumped in separate files

        Ok(())
    }
}
impl<S: GetSingleAtom + NumeratorState> Amplitude<S> {
    pub fn export_expressions(
        &self,
        export_root: &str,
        printer_ops: PrintOptions,
        export_settings: &ExportSettings,
    ) -> Result<(), Report> {
        let path = Path::new(export_root)
            .join("sources")
            .join("amplitudes")
            .join(self.name.as_str())
            .join("expressions");
        std::fs::create_dir_all(&path)?;
        for amplitude_graph in self.amplitude_graphs.iter() {
            let num = Numerator::default()
                .from_graph(
                    &amplitude_graph.graph.bare_graph,
                    &export_settings.numerator_settings.global_prefactor,
                )
                .get_single_atom();

            let dens: Vec<(String, String)> = amplitude_graph
                .graph
                .bare_graph
                .denominator_print(printer_ops);
            let rep_rules: Vec<(String, String)> = amplitude_graph
                .graph
                .bare_graph
                .rep_rules_print(printer_ops);

            let out = (
                format!(
                    "{}",
                    AtomPrinter::new_with_options(num.as_ref().unwrap().0.as_view(), printer_ops)
                ),
                rep_rules,
                dens,
            );

            fs::write(
                path.join(format!(
                    "{}_exp.json",
                    amplitude_graph.graph.bare_graph.name
                )),
                serde_json::to_string_pretty(&out).unwrap(),
            )?;
        }
        Ok(())
    }
}
impl<S: NumeratorState> Amplitude<S> {
    pub fn to_serializable(&self) -> SerializableAmplitude {
        SerializableAmplitude::from_amplitude(self)
    }

    pub fn export_denominator(
        &self,
        export_root: &str,
        printer_ops: PrintOptions,
    ) -> Result<(), Report> {
        let path = Path::new(export_root)
            .join("sources")
            .join("amplitudes")
            .join(self.name.as_str())
            .join("denominator");

        for amplitude_graph in self.amplitude_graphs.iter() {
            let dens: Vec<(String, String)> = amplitude_graph
                .graph
                .bare_graph
                .edges
                .iter()
                .map(|e| {
                    let (mom, mass) = e.denominator(&amplitude_graph.graph.bare_graph);
                    (
                        format!(
                            "{}",
                            AtomPrinter::new_with_options(mom.as_view(), printer_ops)
                        ),
                        format!(
                            "{}",
                            AtomPrinter::new_with_options(mass.as_view(), printer_ops)
                        ),
                    )
                })
                .collect();
            fs::write(
                path.join(format!(
                    "{}_den.json",
                    amplitude_graph.graph.bare_graph.name
                )),
                serde_json::to_string_pretty(&dens).unwrap(),
            )?;
        }
        Ok(())
    }
}
impl Amplitude<PythonState> {
    pub fn generate_integrand(
        &self,
        path_to_settings: &Path,
    ) -> Result<GammaLoopIntegrand, Report> {
        let settings_string = fs::read_to_string(path_to_settings)
            .wrap_err_with(|| {
                format!(
                    "Could not open settings yaml file {}",
                    path_to_settings.display()
                )
            })
            .suggestion("does the path exist?")?;

        let mut settings: Settings = serde_yaml::from_str(&settings_string)
            .wrap_err("Could not parse settings yaml content")
            .suggestion("Is it a correct yaml file")?;

        let amp = Amplitude::<Evaluators>::try_from_python(self.clone())?;

        settings.sync_with_amplitude(&amp)?;

        Ok(GammaLoopIntegrand::amplitude_integrand_constructor(
            amp, settings,
        ))
    }
}
#[derive(Debug, Clone, Default)]
pub struct CrossSectionList {
    pub container: Vec<CrossSection>,
}

impl CrossSectionList {
    pub fn sync(&mut self, _model: &Model) {
        // self.container.iter_mut().for_each(|cs| cs.sync(model));
    }
    pub fn from_file(model: &Model, file_path: String) -> Result<CrossSectionList, Report> {
        let f = File::open(file_path.clone())
            .wrap_err_with(|| format!("Could not open cross-section yaml file {}", file_path))
            .suggestion("Does the path exist?")?;
        serde_yaml::from_reader(f)
            .wrap_err("Could not parse cross-section yaml content")
            .suggestion("Is it a correct yaml file")
            .map(
                |serializable_cross_sections: Vec<SerializableCrossSection>| {
                    CrossSectionList::from_serializable_cross_sections(
                        model,
                        &serializable_cross_sections,
                    )
                },
            )
    }

    pub fn from_yaml_str(model: &Model, yaml_str: String) -> Result<CrossSectionList, Report> {
        serde_yaml::from_str(yaml_str.as_str())
            .wrap_err("Could not parse cross-section yaml content")
            .suggestion("Is it a correct yaml file")
            .map(
                |serializable_cross_sections: Vec<SerializableCrossSection>| {
                    CrossSectionList::from_serializable_cross_sections(
                        model,
                        &serializable_cross_sections,
                    )
                },
            )
    }

    pub fn from_serializable_cross_sections(
        model: &Model,
        serializable_cross_sections: &[SerializableCrossSection],
    ) -> CrossSectionList {
        CrossSectionList {
            container: serializable_cross_sections
                .iter()
                .map(|sg| CrossSection::from_serializable_cross_section(model, sg))
                .collect(),
        }
    }

    pub fn to_serializable(&self) -> Vec<SerializableCrossSection> {
        self.container
            .iter()
            .map(|cs| cs.to_serializable())
            .collect()
    }

    pub fn to_yaml(&self) -> Result<String, Error> {
        serde_yaml::to_string(&self.to_serializable())
    }

    pub fn add_cross_section(&mut self, cross_section: CrossSection) {
        self.container.push(cross_section);
    }
}

#[derive(Debug, Clone, Default)]
pub struct AmplitudeList<S: NumeratorState> {
    pub container: Vec<Amplitude<S>>,
}

impl<S: TypedNumeratorState> AmplitudeList<S> {
    pub fn try_from_python(al: AmplitudeList<PythonState>) -> Result<AmplitudeList<S>> {
        al.map_res(Amplitude::<S>::try_from_python)
    }
}

impl AmplitudeList<PythonState> {
    pub fn sync(&mut self, model: &Model) {
        self.container.iter_mut().for_each(|a| a.sync(model));
    }
}

impl<S: NumeratorState> AmplitudeList<S> {
    pub fn load_derived_data_mut(
        &mut self,
        model: &Model,
        path: &Path,
        settings: &Settings,
    ) -> Result<()> {
        for amplitude in self.container.iter_mut() {
            let ampltitude_path = path.join(amplitude.name.as_str());
            amplitude.load_derived_data_mut(model, &ampltitude_path, settings)?;
            for amplitude_graph in amplitude.amplitude_graphs.iter_mut() {
                amplitude_graph.graph.generate_esurface_data()?;
            }
        }
        Ok(())
    }
    pub fn map<F, U: NumeratorState>(self, f: F) -> AmplitudeList<U>
    where
        F: FnMut(Amplitude<S>) -> Amplitude<U>,
    {
        AmplitudeList {
            container: self.container.into_iter().map(f).collect(),
        }
    }

    pub fn map_mut<F>(&mut self, f: F)
    where
        F: FnMut(&mut Amplitude<S>),
    {
        self.container.iter_mut().for_each(f);
    }

    pub fn map_mut_graphs<F>(&mut self, f: F)
    where
        F: FnMut(&mut Graph<S>) + Copy,
    {
        self.container
            .iter_mut()
            .for_each(|a| a.map_mut(|ag| ag.map_mut(f)));
    }

    pub fn map_res<F, U: NumeratorState, E>(self, f: F) -> Result<AmplitudeList<U>, E>
    where
        F: FnMut(Amplitude<S>) -> Result<Amplitude<U>, E>,
    {
        Ok(AmplitudeList {
            container: self
                .container
                .into_iter()
                .map(f)
                .collect::<Result<_, E>>()?,
        })
    }
}

impl AmplitudeList<UnInit> {
    pub fn from_file(model: &Model, file_path: String) -> Result<Self, Report> {
        let f = File::open(file_path.clone())
            .wrap_err_with(|| format!("Could not open amplitude list yaml file {}", file_path))
            .suggestion("Does the path exist?")?;
        serde_yaml::from_reader(f)
            .wrap_err("Could not parse amplitude list yaml content")
            .suggestion("Is it a correct yaml file")
            .map(|serializable_amplitudes: Vec<SerializableAmplitude>| {
                AmplitudeList::from_serializable_amplitudes(model, &serializable_amplitudes)
            })
    }

    pub fn from_yaml_str(model: &Model, yaml_str: String) -> Result<Self, Report> {
        serde_yaml::from_str(yaml_str.as_str())
            .wrap_err("Could not parse amplitude list yaml content")
            .suggestion("Is it a correct yaml file")
            .map(|serializable_amplitudes: Vec<SerializableAmplitude>| {
                AmplitudeList::from_serializable_amplitudes(model, &serializable_amplitudes)
            })
    }

    pub fn from_serializable_amplitudes(
        model: &Model,
        serializable_amplitudes: &[SerializableAmplitude],
    ) -> Self {
        Self {
            container: serializable_amplitudes
                .iter()
                .map(|sg| Amplitude::from_serializable_amplitude(model, sg))
                .collect(),
        }
    }

    pub fn load_derived_data<S: NumeratorState>(
        self,
        model: &Model,
        path: &Path,
        settings: &Settings,
    ) -> Result<AmplitudeList<S>, Report> {
        self.map_res(|a| {
            let a = a.map_res(|g| {
                g.map_res(|mut g| {
                    g.generate_esurface_data()?;
                    Result::<_, Report>::Ok(g)
                })
            })?;
            let ampltitude_path = path.join(a.name.as_str());
            a.load_derived_data::<S>(model, &ampltitude_path, settings)
        })
    }

    pub fn generate_numerator(
        self,
        export_settings: &ExportSettings,
    ) -> AmplitudeList<AppliedFeynmanRule> {
        let container = self
            .container
            .into_iter()
            .map(|a| a.apply_feynman_rules(export_settings))
            .collect();

        AmplitudeList { container }
    }
}

impl<S: NumeratorState> AmplitudeList<S> {
    pub fn to_serializable(&self) -> Vec<SerializableAmplitude> {
        self.container
            .iter()
            .map(|cs| cs.to_serializable())
            .collect()
    }

    pub fn to_yaml(&self) -> Result<String, Error> {
        serde_yaml::to_string(&self.to_serializable())
    }

    pub fn add_amplitude(&mut self, amplitude: Amplitude<S>) {
        self.container.push(amplitude);
    }
}
