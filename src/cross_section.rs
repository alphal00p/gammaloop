use crate::gammaloop_integrand::GammaLoopIntegrand;
use crate::graph::{Graph, SerializableGraph};
use crate::model::Model;
use crate::{utils::*, Settings};
use bincode;
use color_eyre::{Help, Report};
#[allow(unused_imports)]
use eyre::{eyre, Context};
use serde::{Deserialize, Serialize};
use serde_yaml::Error;
use smartstring::{LazyCompact, SmartString};
use std::fs;
use std::fs::File;
use std::path::Path;
use symbolica::printer::{AtomPrinter, PrintOptions};

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
        graph: &Graph,
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
        graph: &Graph,
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
    pub multiplicity: f64,
    // This identifier of the topology class is mostly a stub for now
    pub topology_class: Vec<usize>,
    pub cuts: Vec<SerializableSuperGraphCut>,
}

impl SerializableSuperGraph {
    pub fn from_supergraph(supergraph: &SuperGraph) -> SerializableSuperGraph {
        SerializableSuperGraph {
            sg_id: supergraph.sg_id,
            graph: SerializableGraph::from_graph(&supergraph.graph),
            multiplicity: supergraph.multiplicity,
            topology_class: supergraph.topology_class.clone(),
            cuts: supergraph
                .cuts
                .iter()
                .map(|cut| SerializableSuperGraphCut::from_supergraph_cut(&supergraph.graph, cut))
                .collect(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct SuperGraph {
    pub sg_id: usize,
    pub graph: Graph,
    pub multiplicity: f64,
    // This identifier of the topology class is mostly a stub for now
    pub topology_class: Vec<usize>,
    pub cuts: Vec<SuperGraphCut>,
}

impl SuperGraph {
    pub fn from_serializable_supergraph(
        model: &Model,
        serializable_supergraph: &SerializableSuperGraph,
    ) -> SuperGraph {
        let g = Graph::from_serializable_graph(model, &serializable_supergraph.graph);
        let cuts = serializable_supergraph
            .cuts
            .iter()
            .map(|cut| SuperGraphCut::from_serializable_supergraph_cut(model, &g, cut))
            .collect();
        SuperGraph {
            sg_id: serializable_supergraph.sg_id,
            graph: g,
            multiplicity: serializable_supergraph.multiplicity,
            topology_class: serializable_supergraph.topology_class.clone(),
            cuts,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SerializableForwardScatteringGraphCut {
    pub cut_edges: Vec<SmartString<LazyCompact>>,
    pub amplitudes: [SerializableAmplitude; 2],
}

impl SerializableForwardScatteringGraphCut {
    pub fn from_forward_scattering_graph_cut(
        graph: &Graph,
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
    pub amplitudes: [Amplitude; 2],
}

impl ForwardScatteringGraphCut {
    pub fn from_serializable_forward_scattering_graph_cut(
        model: &Model,
        graph: &Graph,
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
    pub multiplicity: f64,
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
            multiplicity: forward_scattering_graph.multiplicity,
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
    pub graph: Graph,
    pub multiplicity: f64,
    pub cuts: Vec<ForwardScatteringGraphCut>,
}

impl ForwardScatteringGraph {
    pub fn from_serializable_forward_scattering_graph(
        model: &Model,
        forward_scattering_graph: &SerializableForwardScatteringGraph,
    ) -> ForwardScatteringGraph {
        let g = Graph::from_serializable_graph(model, &forward_scattering_graph.graph);
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
            multiplicity: forward_scattering_graph.multiplicity,
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
    pub multi_channeling_channels: Vec<usize>, // empty list defaults to all channels if multi_channeling is enabled
}

impl SerializableAmplitudeGraph {
    pub fn from_amplitude_graph(amplitude_graph: &AmplitudeGraph) -> SerializableAmplitudeGraph {
        SerializableAmplitudeGraph {
            sg_id: amplitude_graph.sg_id,
            sg_cut_id: amplitude_graph.sg_cut_id,
            fs_cut_id: amplitude_graph.fs_cut_id,
            amplitude_side: amplitude_graph.amplitude_side.clone(),
            graph: SerializableGraph::from_graph(&amplitude_graph.graph),
            multi_channeling_channels: amplitude_graph.multi_channeling_channels.clone(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct AmplitudeGraph {
    pub sg_id: usize,
    pub sg_cut_id: usize,
    pub fs_cut_id: usize,
    pub amplitude_side: Side,
    pub graph: Graph,
    pub multi_channeling_channels: Vec<usize>,
}

impl AmplitudeGraph {
    pub fn from_amplitude_graph(
        model: &Model,
        amplitude_graph: &SerializableAmplitudeGraph,
    ) -> AmplitudeGraph {
        AmplitudeGraph {
            sg_id: amplitude_graph.sg_id,
            sg_cut_id: amplitude_graph.sg_cut_id,
            fs_cut_id: amplitude_graph.fs_cut_id,
            amplitude_side: amplitude_graph.amplitude_side.clone(),
            graph: Graph::from_serializable_graph(model, &amplitude_graph.graph),
            multi_channeling_channels: amplitude_graph.multi_channeling_channels.clone(),
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
    pub fn from_amplitude(amplitude: &Amplitude) -> SerializableAmplitude {
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

#[derive(Debug, Clone)]
pub struct Amplitude {
    pub name: SmartString<LazyCompact>,
    pub amplitude_graphs: Vec<AmplitudeGraph>,
}

impl Amplitude {
    pub fn from_file(model: &Model, file_path: String) -> Result<Amplitude, Report> {
        SerializableAmplitude::from_file(file_path).map(|serializable_amplitude| {
            Amplitude::from_serializable_amplitude(model, &serializable_amplitude)
        })
    }

    pub fn from_yaml_str(model: &Model, yaml_str: String) -> Result<Amplitude, Report> {
        SerializableAmplitude::from_yaml_str(yaml_str).map(|serializable_amplitude| {
            Amplitude::from_serializable_amplitude(model, &serializable_amplitude)
        })
    }

    pub fn from_serializable_amplitude(
        model: &Model,
        serializable_amplitude: &SerializableAmplitude,
    ) -> Amplitude {
        Amplitude {
            name: serializable_amplitude.name.clone(),
            amplitude_graphs: serializable_amplitude
                .amplitude_graphs
                .iter()
                .map(|sg| AmplitudeGraph::from_amplitude_graph(model, sg))
                .collect(),
        }
    }

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
                .edges
                .iter()
                .map(|e| {
                    let (mom, mass) = e.denominator(&amplitude_graph.graph);
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
                path.join(format!("{}_den.json", amplitude_graph.graph.name)),
                serde_json::to_string_pretty(&dens).unwrap(),
            )?;
        }
        Ok(())
    }

    pub fn export_expressions(
        &self,
        export_root: &str,
        printer_ops: PrintOptions,
    ) -> Result<(), Report> {
        let path = Path::new(export_root)
            .join("sources")
            .join("amplitudes")
            .join(self.name.as_str())
            .join("expressions");
        for amplitude_graph in self.amplitude_graphs.iter() {
            if let Some(num) = &amplitude_graph.graph.derived_data.numerator {
                let dens: Vec<(String, String)> = amplitude_graph
                    .graph
                    .edges
                    .iter()
                    .map(|e| {
                        let (mom, mass) = e.denominator(&amplitude_graph.graph);
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

                let rep_rules: Vec<(String, String)> = amplitude_graph
                    .graph
                    .generate_lmb_replacement_rules()
                    .iter()
                    .map(|(lhs, rhs)| {
                        (
                            format!(
                                "{}",
                                AtomPrinter::new_with_options(lhs.as_view(), printer_ops)
                            ),
                            format!(
                                "{}",
                                AtomPrinter::new_with_options(rhs.as_view(), printer_ops)
                            ),
                        )
                    })
                    .collect();

                let out = (
                    format!(
                        "{}",
                        AtomPrinter::new_with_options(num.as_view(), printer_ops)
                    ),
                    rep_rules,
                    dens,
                );

                fs::write(
                    path.join(format!("{}_exp.json", amplitude_graph.graph.name)),
                    serde_json::to_string_pretty(&out).unwrap(),
                )?;
            }
        }
        Ok(())
    }

    #[allow(unused)]
    pub fn export(&mut self, export_root: &str, model: &Model) -> Result<(), Report> {
        // TODO process amplitude by adding lots of additional information necessary for runtime.
        // e.g. generate e-surface, cff expression, counterterms, etc.

        // generate cff and ltd for each graph in the ampltiudes, ltd also generates lmbs
        for amplitude_graph in self.amplitude_graphs.iter_mut() {
            amplitude_graph.graph.generate_cff();
            amplitude_graph.graph.generate_ltd();
            amplitude_graph.graph.generate_tropical_subgraph_table();
        }

        // Then dumped the new yaml representation of the amplitude now containing all that additional information
        let path = Path::new(export_root)
            .join("sources")
            .join("amplitudes")
            .join(self.name.as_str());

        fs::write(
            path.join("amplitude.yaml"),
            serde_yaml::to_string(&self.to_serializable())?,
        )?;

        // dump the derived data in a binary file
        for amplitude_graph in self.amplitude_graphs.iter() {
            fs::write(
                path.join(format!("derived_data_{}.bin", amplitude_graph.graph.name)),
                bincode::serialize(&amplitude_graph.graph.derived_data.to_serializable())?,
            );
        }

        // Additional files can be written too, e.g. the lengthy cff expressions can be dumped in separate files

        Ok(())
    }

    pub fn load_derived_data(&mut self, path: &Path) -> Result<(), Report> {
        for ampltitude_graph in self.amplitude_graphs.iter_mut() {
            let graph_path = path.join(format!(
                "derived_data_{}.bin",
                ampltitude_graph.graph.name.as_str()
            ));
            ampltitude_graph.graph.load_derived_data(&graph_path)?;
        }
        Ok(())
    }

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

        let settings: Settings = serde_yaml::from_str(&settings_string)
            .wrap_err("Could not parse settings yaml content")
            .suggestion("Is it a correct yaml file")?;

        Ok(GammaLoopIntegrand::amplitude_integrand_constructor(
            self.clone(),
            settings.clone(),
        ))
    }
}

#[derive(Debug, Clone, Default)]
pub struct CrossSectionList {
    pub container: Vec<CrossSection>,
}

impl CrossSectionList {
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
pub struct AmplitudeList {
    pub container: Vec<Amplitude>,
}

impl AmplitudeList {
    pub fn from_file(model: &Model, file_path: String) -> Result<AmplitudeList, Report> {
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

    pub fn from_yaml_str(model: &Model, yaml_str: String) -> Result<AmplitudeList, Report> {
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
    ) -> AmplitudeList {
        AmplitudeList {
            container: serializable_amplitudes
                .iter()
                .map(|sg| Amplitude::from_serializable_amplitude(model, sg))
                .collect(),
        }
    }

    pub fn to_serializable(&self) -> Vec<SerializableAmplitude> {
        self.container
            .iter()
            .map(|cs| cs.to_serializable())
            .collect()
    }

    pub fn to_yaml(&self) -> Result<String, Error> {
        serde_yaml::to_string(&self.to_serializable())
    }

    pub fn add_amplitude(&mut self, amplitude: Amplitude) {
        self.container.push(amplitude);
    }

    pub fn load_derived_data(&mut self, path: &str) -> Result<(), Report> {
        let path = Path::new(path);
        for amplitude in self.container.iter_mut() {
            let ampltitude_path = path.join(amplitude.name.as_str());
            amplitude.load_derived_data(&ampltitude_path)?;
        }
        Ok(())
    }

    pub fn generate_numerator(&mut self, model: &Model) {
        for amplitude in self.container.iter_mut() {
            for amplitude_graph in amplitude.amplitude_graphs.iter_mut() {
                amplitude_graph.graph.generate_numerator(model);
            }
        }
    }
}
