use std::{collections::BTreeSet, fmt};

use color_eyre::Result;
use eyre::{eyre, Context};
use gammalooprs::{
    graph::Graph,
    integrands::process::ActiveF64Backend,
    processes::{Amplitude, CrossSection, CrossSectionCut, ProcessCollection, RaisedCutId},
    settings::global::FrozenCompilationMode,
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
    pub left_threshold_esurface_ids: Vec<usize>,
    pub right_threshold_esurface_ids: Vec<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct IntegrandThresholdEsurfaceInfo {
    pub esurface_id: usize,
    pub edge_ids: Vec<usize>,
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
                graph_groups: amplitude_graph_groups(integrand),
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
                graph_groups: cross_section_graph_groups(integrand),
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
) -> Vec<usize> {
    let mut raising_powers = vec![1; graph.cuts.len()];
    for raised_cut_group in graph.raised_data.raised_cut_groups.iter() {
        let raising_power = raised_cut_group.related_esurface_group.max_occurence;
        for cut_id in &raised_cut_group.cuts {
            raising_powers[cut_id.0] = raising_power;
        }
    }
    raising_powers
}

fn lmb_channel_ids(
    graph: &Graph,
    lmbs: &typed_index_collections::TiVec<
        gammalooprs::graph::LmbIndex,
        gammalooprs::graph::LoopMomentumBasis,
    >,
) -> Vec<Option<usize>> {
    let mut channels = Vec::new();

    for (lmb_index, lmb) in lmbs.iter_enumerated() {
        let massless_edges = lmb
            .loop_edges
            .iter()
            .filter(|&&edge_id| graph.underlying[edge_id].particle.mass_atom().is_zero())
            .collect::<Vec<_>>();

        if massless_edges.is_empty() {
            continue;
        }

        if channels
            .iter()
            .any(|included_channel: &gammalooprs::graph::LmbIndex| {
                let basis = &lmbs[*included_channel].loop_edges;
                massless_edges.iter().all(|edge_id| basis.contains(edge_id))
            })
        {
            continue;
        }

        channels.push(lmb_index);
    }

    channels.sort_by_key(|lmb_index| usize::from(*lmb_index));
    if channels.is_empty() {
        if let Some(current_lmb_index) = lmbs
            .iter_enumerated()
            .find(|(_lmb_index, lmb)| lmb.loop_edges == graph.loop_momentum_basis.loop_edges)
            .map(|(lmb_index, _)| lmb_index)
        {
            channels.push(current_lmb_index);
        }
    }

    let mut channel_ids = vec![None; lmbs.len()];
    for (channel_id, lmb_index) in channels.into_iter().enumerate() {
        channel_ids[usize::from(lmb_index)] = Some(channel_id);
    }
    channel_ids
}

fn amplitude_graph_groups(
    integrand: &gammalooprs::integrands::process::amplitude::AmplitudeIntegrand,
) -> Vec<IntegrandGraphGroupInfo> {
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
            let channel_ids = lmb_channel_ids(&master_graph.graph, &master_graph.lmbs);
            let threshold_esurface_ids = master_graph
                .threshold_counterterm
                .generated_mask
                .iter_enumerated()
                .filter_map(|(esurface_id, is_generated)| is_generated.then_some(esurface_id.0))
                .collect::<Vec<_>>();
            let threshold_esurfaces = threshold_esurface_ids
                .iter()
                .copied()
                .map(|esurface_id| IntegrandThresholdEsurfaceInfo {
                    esurface_id,
                    edge_ids: threshold_esurface_edge_ids(&master_graph.esurfaces, esurface_id),
                })
                .collect::<Vec<_>>();
            IntegrandGraphGroupInfo {
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
                cuts: Vec::new(),
            }
        })
        .collect()
}

fn cross_section_graph_groups(
    integrand: &gammalooprs::integrands::process::cross_section::CrossSectionIntegrand,
) -> Vec<IntegrandGraphGroupInfo> {
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
            let channel_ids = lmb_channel_ids(&master_graph.graph, &master_graph.lmbs);
            let cut_raising_powers = cut_raising_powers(master_graph);
            let mut cut_to_raised_cut = vec![None; master_graph.cuts.len()];
            for (raised_cut_id, raised_cut_group) in
                master_graph.raised_data.raised_cut_groups.iter_enumerated()
            {
                for cut_id in &raised_cut_group.cuts {
                    cut_to_raised_cut[cut_id.0] = Some(raised_cut_id);
                }
            }

            let cuts = master_graph
                .cuts
                .iter()
                .enumerate()
                .map(|(cut_id, cut)| {
                    let (left_threshold_esurface_ids, right_threshold_esurface_ids) =
                        cut_to_raised_cut[cut_id].map_or_else(
                            || (Vec::new(), Vec::new()),
                            |raised_cut_id: RaisedCutId| {
                                master_graph.threshold_esurface_ids_for_raised_cut(raised_cut_id)
                            },
                        );

                    IntegrandCutInfo {
                        cut_id,
                        edge_ids: cut_edge_ids(&master_graph.graph, cut),
                        raising_power: cut_raising_powers[cut_id],
                        left_threshold_esurface_ids,
                        right_threshold_esurface_ids,
                    }
                })
                .collect::<Vec<_>>();

            let threshold_esurface_ids = cuts
                .iter()
                .flat_map(|cut| {
                    cut.left_threshold_esurface_ids
                        .iter()
                        .chain(cut.right_threshold_esurface_ids.iter())
                })
                .copied()
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();

            let threshold_esurfaces = threshold_esurface_ids
                .iter()
                .copied()
                .map(|esurface_id| IntegrandThresholdEsurfaceInfo {
                    esurface_id,
                    edge_ids: threshold_esurface_edge_ids(
                        &master_graph.graph.surface_cache.esurface_cache,
                        esurface_id,
                    ),
                })
                .collect::<Vec<_>>();

            IntegrandGraphGroupInfo {
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
            }
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
