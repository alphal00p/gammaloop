use std::fmt;

use color_eyre::Result;
use eyre::{eyre, Context};
use gammalooprs::{
    graph::Graph,
    integrands::process::{ActiveF64Backend, LmbMultiChannelingSetup, ParamBuilder},
    model::Model,
    processes::{
        Amplitude, CrossSection, CrossSectionCut, CutId, ProcessCollection,
        ThresholdCountertermAssociation, ThresholdCountertermStatus,
    },
    settings::{global::FrozenCompilationMode, runtime::ParameterizationSettings, RuntimeSettings},
};
use linnet::half_edge::involution::{EdgeVec, Orientation};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::state::State;

#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum IntegrandKind {
    Amplitude,
    CrossSection,
}

impl fmt::Display for IntegrandKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Amplitude => f.write_str("amplitude"),
            Self::CrossSection => f.write_str("cross section"),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct IntegrandGraphInfo {
    pub graph_id: usize,
    pub name: String,
    pub is_master: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct IntegrandOrientationInfo {
    pub orientation_id: usize,
    pub signature: Vec<i8>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct IntegrandLoopMomentumBasisInfo {
    pub basis_id: usize,
    pub channel_id: Option<usize>,
    pub edge_ids: Vec<usize>,
    pub matches_generation_basis: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct IntegrandCutInfo {
    pub cut_id: usize,
    pub edge_ids: Vec<usize>,
    pub raising_power: usize,
    pub left_thresholds: Vec<IntegrandCutThresholdInfo>,
    pub right_thresholds: Vec<IntegrandCutThresholdInfo>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum IntegrandThresholdStatus {
    NoRadialDependence,
    AlwaysPinched,
    CanBecomePinched,
    ProvenNonExisting,
    PotentiallyExisting,
}

impl IntegrandThresholdStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::NoRadialDependence => "no_radial_dependence",
            Self::AlwaysPinched => "always_pinched",
            Self::CanBecomePinched => "can_become_pinched",
            Self::ProvenNonExisting => "proven_non_existing",
            Self::PotentiallyExisting => "potentially_existing",
        }
    }

    pub fn is_currently_viable(self) -> bool {
        matches!(self, Self::CanBecomePinched | Self::PotentiallyExisting)
    }

    pub fn can_become_pinched(self) -> bool {
        self == Self::CanBecomePinched
    }
}

impl From<ThresholdCountertermStatus> for IntegrandThresholdStatus {
    fn from(value: ThresholdCountertermStatus) -> Self {
        match value {
            ThresholdCountertermStatus::NoRadialDependence => Self::NoRadialDependence,
            ThresholdCountertermStatus::AlwaysPinched => Self::AlwaysPinched,
            ThresholdCountertermStatus::CanBecomePinched => Self::CanBecomePinched,
            ThresholdCountertermStatus::ProvenNonExisting => Self::ProvenNonExisting,
            ThresholdCountertermStatus::PotentiallyExisting => Self::PotentiallyExisting,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct IntegrandCutThresholdInfo {
    pub esurface_id: usize,
    pub status: IntegrandThresholdStatus,
    pub cut_boundary_edge_ids: Vec<usize>,
    pub threshold_boundary_edge_ids: Vec<usize>,
    pub invariant_bound_is_applicable: bool,
}

impl IntegrandCutThresholdInfo {
    fn from_association(
        association: &ThresholdCountertermAssociation,
        graph: &Graph,
        cut: &CrossSectionCut,
        model: &Model,
        param_builder: &ParamBuilder,
        settings: &RuntimeSettings,
    ) -> Self {
        Self {
            esurface_id: association.esurface_id.0,
            status: association
                .classify_for_model(
                    graph,
                    cut,
                    model,
                    param_builder,
                    settings,
                    settings.subtraction.esurface_existence_threshold,
                )
                .into(),
            cut_boundary_edge_ids: association
                .cut_boundary_edges
                .iter()
                .map(|edge_id| edge_id.0)
                .collect(),
            threshold_boundary_edge_ids: association
                .threshold_boundary_edges
                .iter()
                .map(|edge_id| edge_id.0)
                .collect(),
            invariant_bound_is_applicable: association.invariant_bound_is_applicable,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct IntegrandActiveThresholdCutInfo {
    pub cut_id: usize,
    pub can_become_pinched: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct IntegrandThresholdEsurfaceInfo {
    pub esurface_id: usize,
    pub edge_ids: Vec<usize>,
    pub active_cuts: Vec<IntegrandActiveThresholdCutInfo>,
}

impl IntegrandCutInfo {
    fn from_cross_section_cut(
        graph_term: &gammalooprs::integrands::process::cross_section::CrossSectionGraphTerm,
        cut_id: CutId,
        cut: &CrossSectionCut,
        raising_power: usize,
        model: &Model,
        param_builder: &ParamBuilder,
        settings: &RuntimeSettings,
    ) -> Self {
        let associations = &graph_term.cut_threshold_associations[cut_id];
        let threshold_info = |association| {
            IntegrandCutThresholdInfo::from_association(
                association,
                &graph_term.graph,
                cut,
                model,
                param_builder,
                settings,
            )
        };

        Self {
            // Retain `CutId` for internal indexing and convert only at the API boundary.
            cut_id: usize::from(cut_id),
            edge_ids: cut_edge_ids(&graph_term.graph, cut),
            raising_power,
            left_thresholds: associations.left.iter().map(&threshold_info).collect(),
            right_thresholds: associations.right.iter().map(threshold_info).collect(),
        }
    }

    fn active_threshold_cut(&self, esurface_id: usize) -> Option<IntegrandActiveThresholdCutInfo> {
        let mut is_currently_viable = false;
        let mut can_become_pinched = false;

        for threshold in self
            .left_thresholds
            .iter()
            .chain(self.right_thresholds.iter())
            .filter(|threshold| threshold.esurface_id == esurface_id)
        {
            is_currently_viable |= threshold.status.is_currently_viable();
            can_become_pinched |= threshold.status.can_become_pinched();
        }

        is_currently_viable.then_some(IntegrandActiveThresholdCutInfo {
            cut_id: self.cut_id,
            can_become_pinched,
        })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct IntegrandGraphGroupInfo {
    pub group_id: usize,
    pub graphs: Vec<IntegrandGraphInfo>,
    pub orientation_edge_ids: Vec<usize>,
    pub orientations: Vec<IntegrandOrientationInfo>,
    pub loop_momentum_bases: Vec<IntegrandLoopMomentumBasisInfo>,
    pub threshold_esurface_ids: Vec<usize>,
    pub threshold_esurfaces: Vec<IntegrandThresholdEsurfaceInfo>,
    pub cuts: Vec<IntegrandCutInfo>,
}

fn threshold_esurface_edge_ids(
    esurfaces: &gammalooprs::cff::esurface::EsurfaceCollection,
    esurface_id: usize,
) -> Vec<usize> {
    esurfaces
        .iter()
        .nth(esurface_id)
        .expect("threshold esurface id should resolve in collection")
        .energies
        .iter()
        .map(|edge_id| edge_id.0)
        .collect()
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct IntegrandInfo {
    pub process_id: usize,
    pub process_name: String,
    pub integrand_name: String,
    pub kind: IntegrandKind,
    pub generation_compilation: FrozenCompilationMode,
    pub active_f64_backend: ActiveF64Backend,
    pub graph_count: usize,
    pub graph_group_count: usize,
    pub record_size_bytes: usize,
    pub graph_groups: Vec<IntegrandGraphGroupInfo>,
}

pub(crate) fn collect_integrand_info(
    state: &State,
    process_id: usize,
    integrand_name: &str,
) -> Result<IntegrandInfo> {
    let process = &state.process_list.processes[process_id];
    let resolved = process.get_integrand(integrand_name)?;
    let generated = resolved.require_generated()?;
    let generation_compilation = generated.frozen_compilation().clone();
    let active_f64_backend = generated.active_f64_backend();

    match (generated, &process.collection) {
        (
            gammalooprs::integrands::process::ProcessIntegrand::Amplitude(integrand),
            ProcessCollection::Amplitudes(amplitudes),
        ) => {
            let amplitude = amplitudes.get(&resolved.canonical_name).ok_or_else(|| {
                eyre!(
                    "Could not resolve amplitude record '{}' in process #{} ({})",
                    resolved.canonical_name,
                    process.definition.process_id,
                    process.definition.folder_name
                )
            })?;
            Ok(IntegrandInfo {
                process_id: process.definition.process_id,
                process_name: process.definition.folder_name.clone(),
                integrand_name: resolved.canonical_name,
                kind: IntegrandKind::Amplitude,
                generation_compilation: generation_compilation.clone(),
                active_f64_backend,
                graph_count: amplitude.graphs.len(),
                graph_group_count: amplitude.graph_group_structure.len(),
                record_size_bytes: integrand_record_size_from_amplitude(amplitude)?,
                graph_groups: amplitude_graph_groups(integrand)?,
            })
        }
        (
            gammalooprs::integrands::process::ProcessIntegrand::CrossSection(integrand),
            ProcessCollection::CrossSections(cross_sections),
        ) => {
            let cross_section = cross_sections
                .get(&resolved.canonical_name)
                .ok_or_else(|| {
                    eyre!(
                        "Could not resolve cross-section record '{}' in process #{} ({})",
                        resolved.canonical_name,
                        process.definition.process_id,
                        process.definition.folder_name
                    )
                })?;
            let model = state.resolve_model_for_integrand(
                process.definition.process_id,
                &resolved.canonical_name,
            )?;
            Ok(IntegrandInfo {
                process_id: process.definition.process_id,
                process_name: process.definition.folder_name.clone(),
                integrand_name: resolved.canonical_name,
                kind: IntegrandKind::CrossSection,
                generation_compilation,
                active_f64_backend,
                graph_count: cross_section.supergraphs.len(),
                graph_group_count: cross_section.graph_group_structure.len(),
                record_size_bytes: integrand_record_size_from_cross_section(cross_section)?,
                graph_groups: cross_section_graph_groups(integrand, &model)?,
            })
        }
        _ => Err(eyre!(
            "Process/integrand type mismatch for '{}' in process #{} ({})",
            resolved.canonical_name,
            process.definition.process_id,
            process.definition.folder_name
        )),
    }
}

fn orientation_signature(orientation: &EdgeVec<Orientation>) -> Vec<i8> {
    orientation
        .iter()
        .map(|(_, orientation)| match *orientation {
            Orientation::Default => 1,
            Orientation::Reversed => -1,
            Orientation::Undirected => 0,
        })
        .collect()
}

fn orientation_edge_ids(orientation: &EdgeVec<Orientation>) -> Vec<usize> {
    orientation.iter().map(|(edge_id, _)| edge_id.0).collect()
}

fn cut_edge_ids(graph: &Graph, cut: &CrossSectionCut) -> Vec<usize> {
    cut.cut
        .iter_edges(&graph.underlying)
        .map(|(_, edge)| {
            graph
                .underlying
                .iter_edges()
                .find_map(|(_, edge_id, edge_data)| {
                    (edge_data.data.name == edge.data.name).then_some(edge_id.0)
                })
                .expect("cut edge should resolve in graph")
        })
        .collect()
}

fn cut_raising_powers(
    graph: &gammalooprs::integrands::process::cross_section::CrossSectionGraphTerm,
) -> typed_index_collections::TiVec<CutId, usize> {
    let mut raising_powers: typed_index_collections::TiVec<CutId, usize> =
        vec![1; graph.cuts.len()].into();
    for cut_group in graph.cut_group_data.cut_groups.iter() {
        let raising_power = cut_group.related_esurface_group.max_occurence;
        for cut_id in &cut_group.cuts {
            raising_powers[*cut_id] = raising_power;
        }
    }
    raising_powers
}

fn lmb_channel_ids(
    lmbs: &typed_index_collections::TiVec<
        gammalooprs::graph::LmbIndex,
        gammalooprs::graph::LoopMomentumBasis,
    >,
    multi_channeling_setup: &LmbMultiChannelingSetup,
    graph_name: &str,
    parameterization_settings: &ParameterizationSettings,
) -> Result<Vec<Option<usize>>> {
    let mut channel_ids = vec![None; lmbs.len()];
    for (channel_id, lmb_index) in multi_channeling_setup
        .effective_channels(graph_name, parameterization_settings)?
        .into_iter()
        .enumerate()
    {
        channel_ids[usize::from(lmb_index)] = Some(channel_id);
    }
    Ok(channel_ids)
}

fn amplitude_graph_groups(
    integrand: &gammalooprs::integrands::process::amplitude::AmplitudeIntegrand,
) -> Result<Vec<IntegrandGraphGroupInfo>> {
    let parameterization_settings = integrand
        .settings
        .sampling
        .get_parameterization_settings()
        .unwrap_or_default();
    integrand
        .data
        .graph_group_structure
        .iter_enumerated()
        .map(|(group_id, group)| {
            let master_graph_id = group
                .into_iter()
                .next()
                .expect("graph group should not be empty");
            let master_graph = &integrand.data.graph_terms[master_graph_id];
            let channel_ids = lmb_channel_ids(
                &master_graph.lmbs,
                &master_graph.multi_channeling_setup,
                &master_graph.graph.name,
                &parameterization_settings,
            )?;
            // Amplitude overlap centers are solved from the canonical union across all members
            // of the graph group. Report that same union rather than the master graph's local
            // raised-surface numbering.
            let threshold_esurfaces = integrand.data.group_derived_data[group_id]
                .esurface_map
                .iter_enumerated()
                .filter_map(|(group_esurface_id, raised_esurface_map)| {
                    let (graph_term, raised_esurface_id) = raised_esurface_map
                        .iter_enumerated()
                        .filter_map(|(graph_group_position, raised_esurface_id)| {
                            raised_esurface_id.map(|raised_esurface_id| {
                                let graph_id = group[graph_group_position];
                                (&integrand.data.graph_terms[graph_id], raised_esurface_id)
                            })
                        })
                        .find(|(graph_term, raised_esurface_id)| {
                            graph_term.threshold_counterterm.generated_mask[*raised_esurface_id]
                        })?;
                    let local_esurface_id =
                        graph_term.threshold_counterterm.raised_data.raised_groups
                            [raised_esurface_id]
                            .esurface_ids[0];

                    Some(IntegrandThresholdEsurfaceInfo {
                        esurface_id: group_esurface_id.0,
                        edge_ids: threshold_esurface_edge_ids(
                            &graph_term.esurfaces,
                            local_esurface_id.0,
                        ),
                        active_cuts: Vec::new(),
                    })
                })
                .collect::<Vec<_>>();
            let threshold_esurface_ids = threshold_esurfaces
                .iter()
                .map(|threshold| threshold.esurface_id)
                .collect();
            Ok(IntegrandGraphGroupInfo {
                group_id: usize::from(group_id),
                graphs: group
                    .into_iter()
                    .map(|graph_id| IntegrandGraphInfo {
                        graph_id,
                        name: integrand.data.graph_terms[graph_id].graph.name.clone(),
                        is_master: graph_id == master_graph_id,
                    })
                    .collect(),
                orientation_edge_ids: master_graph
                    .orientations
                    .first()
                    .map(orientation_edge_ids)
                    .unwrap_or_default(),
                orientations: master_graph
                    .orientations
                    .iter()
                    .enumerate()
                    .map(|(orientation_id, orientation)| IntegrandOrientationInfo {
                        orientation_id,
                        signature: orientation_signature(orientation),
                    })
                    .collect(),
                loop_momentum_bases: master_graph
                    .lmbs
                    .iter_enumerated()
                    .map(|(basis_id, lmb)| IntegrandLoopMomentumBasisInfo {
                        basis_id: usize::from(basis_id),
                        channel_id: channel_ids[usize::from(basis_id)],
                        edge_ids: lmb.loop_edges.iter().map(|edge_id| edge_id.0).collect(),
                        matches_generation_basis: lmb.loop_edges
                            == master_graph.graph.loop_momentum_basis.loop_edges,
                    })
                    .collect(),
                threshold_esurface_ids,
                threshold_esurfaces,
                cuts: Vec::new(),
            })
        })
        .collect()
}

fn cross_section_graph_groups(
    integrand: &gammalooprs::integrands::process::cross_section::CrossSectionIntegrand,
    model: &Model,
) -> Result<Vec<IntegrandGraphGroupInfo>> {
    let parameterization_settings = integrand
        .settings
        .sampling
        .get_parameterization_settings()
        .unwrap_or_default();
    integrand
        .data
        .graph_group_structure
        .iter()
        .enumerate()
        .map(|(group_id, group)| {
            let master_graph_id = group
                .into_iter()
                .next()
                .expect("graph group should not be empty");
            let master_graph = &integrand.data.graph_terms[master_graph_id];
            let mut active_model_param_builder: ParamBuilder =
                master_graph.graph.param_builder.clone();
            active_model_param_builder.update_model_values(model);
            let channel_ids = lmb_channel_ids(
                &master_graph.lmbs,
                &master_graph.multi_channeling_setup,
                &master_graph.graph.name,
                &parameterization_settings,
            )?;
            let cut_raising_powers = cut_raising_powers(master_graph);

            let cuts = master_graph
                .cuts
                .iter_enumerated()
                .map(|(cut_id, cut)| {
                    IntegrandCutInfo::from_cross_section_cut(
                        master_graph,
                        cut_id,
                        cut,
                        cut_raising_powers[cut_id],
                        model,
                        &active_model_param_builder,
                        &integrand.settings,
                    )
                })
                .collect::<Vec<_>>();

            let threshold_esurface_ids = master_graph
                .threshold_candidate_esurface_ids
                .iter()
                .map(|esurface_id| esurface_id.0)
                .collect::<Vec<_>>();

            let threshold_esurfaces = threshold_esurface_ids
                .iter()
                .copied()
                .map(|esurface_id| {
                    let active_cuts = cuts
                        .iter()
                        .filter_map(|cut| cut.active_threshold_cut(esurface_id))
                        .collect();

                    IntegrandThresholdEsurfaceInfo {
                        esurface_id,
                        edge_ids: threshold_esurface_edge_ids(
                            &master_graph.graph.surface_cache.esurface_cache,
                            esurface_id,
                        ),
                        active_cuts,
                    }
                })
                .collect::<Vec<_>>();

            Ok(IntegrandGraphGroupInfo {
                group_id,
                graphs: group
                    .into_iter()
                    .map(|graph_id| IntegrandGraphInfo {
                        graph_id,
                        name: integrand.data.graph_terms[graph_id].graph.name.clone(),
                        is_master: graph_id == master_graph_id,
                    })
                    .collect(),
                orientation_edge_ids: master_graph
                    .orientations
                    .first()
                    .map(orientation_edge_ids)
                    .unwrap_or_default(),
                orientations: master_graph
                    .orientations
                    .iter()
                    .enumerate()
                    .map(|(orientation_id, orientation)| IntegrandOrientationInfo {
                        orientation_id,
                        signature: orientation_signature(orientation),
                    })
                    .collect(),
                loop_momentum_bases: master_graph
                    .lmbs
                    .iter_enumerated()
                    .map(|(basis_id, lmb)| IntegrandLoopMomentumBasisInfo {
                        basis_id: usize::from(basis_id),
                        channel_id: channel_ids[usize::from(basis_id)],
                        edge_ids: lmb.loop_edges.iter().map(|edge_id| edge_id.0).collect(),
                        matches_generation_basis: lmb.loop_edges
                            == master_graph.graph.loop_momentum_basis.loop_edges,
                    })
                    .collect(),
                threshold_esurface_ids,
                threshold_esurfaces,
                cuts,
            })
        })
        .collect()
}

fn integrand_record_size_from_amplitude(amplitude: &Amplitude) -> Result<usize> {
    let mut record = amplitude.clone();
    record.integrand = None;
    let encoded =
        bincode::encode_to_vec(&record, bincode::config::standard()).with_context(|| {
            format!(
                "While serializing amplitude '{}' for integrand info",
                amplitude.name
            )
        })?;
    Ok(encoded.len())
}

fn integrand_record_size_from_cross_section(cross_section: &CrossSection) -> Result<usize> {
    let mut record = cross_section.clone();
    record.integrand = None;
    let encoded =
        bincode::encode_to_vec(&record, bincode::config::standard()).with_context(|| {
            format!(
                "While serializing cross-section '{}' for integrand info",
                cross_section.name
            )
        })?;
    Ok(encoded.len())
}
