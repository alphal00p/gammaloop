use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::{Display, Formatter},
    fs::{self, File},
    io::Write,
    ops::{Index, IndexMut},
    path::{Path, PathBuf},
};

use ahash::HashMap;
// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use itertools::Itertools;
use rayon::{
    ThreadPool,
    iter::{IntoParallelRefMutIterator, ParallelIterator},
};
use spenso::algebra::algebraic_traits::IsZero;
use tracing::info;
use vakint::Vakint;

use crate::{
    DependentMomentaConstructor, GammaLoopContext, GammaLoopContextContainer,
    cff::{
        CutCFFIndex,
        esurface::{RaisedEsurfaceData, RaisedEsurfaceGroup, RaisedEsurfaceId},
        expression::{CFFExpression, OrientationID},
    },
    debug_tags, define_index,
    graph::{
        GraphGroup, GroupId, LMBext, LmbChannelFallback, LmbIndex, LoopMomentumBasis,
        ThresholdPinchStatus,
        cuts::{CutSet, ResidueSelector},
        edge::EdgeMass,
        parse::complete_group_parsing,
    },
    integrands::process::{
        GenericEvaluator, LmbMultiChannelingSetup, ParamBuilder,
        cross_section::CrossSectionIntegrandData, graph_to_group_id_for_group_structure,
    },
    model::ArcParticle,
    momentum::{
        Helicity,
        sample::{ExternalIndex, SubspaceData},
    },
    processes::{
        DotExportSettings, EvaluatorSettings, GraphCutSelectionSubject, GraphGenerationStats,
        GraphGroupSelectionMode, GraphGroupSelectionPlan, GraphGroupSelectionSpec,
        GraphSelectionSignatureInventory, GraphSelectionSubject, NamedGraphGenerationReport,
        RaisedCutSignatureInventory,
    },
    settings::{
        GlobalSettings, RuntimeSettings, global::GenerationSettings, runtime::LockedRuntimeSettings,
    },
    utils::{
        DEFAULT_ESURFACE_EXISTENCE_THRESHOLD, F, GS, W_,
        hyperdual_utils::{shape_from_cut_cff_index, simple_n_deriv_shape},
    },
    uv::{approx::CutStructure, forest::ParametricIntegrands, wood::CutWoods},
};
use eyre::{Context, eyre};
use linnet::half_edge::{
    involution::{EdgeIndex, EdgeVec, Orientation},
    subgraph::{
        HedgeNode, Inclusion, InternalSubGraph, ModifySubSet, OrientedCut, SuBitGraph,
        SubGraphLike, SubSetLike, SubSetOps,
    },
};
use serde::{Deserialize, Serialize};
use symbolica::{domains::dual::HyperDual, prelude::*};
use tracing::{debug, warn};
use typed_index_collections::{TiVec, ti_vec};

use super::generation_progress::{self, GenerationProcessKind, GenerationProgressPhase};

use crate::{
    cff::esurface::{Esurface, EsurfaceID},
    graph::{ExternalConnection, FeynmanGraph, Graph},
    integrands::process::{
        ProcessIntegrand,
        cross_section::{CrossSectionGraphTerm, CrossSectionIntegrand},
    },
    model::Model,
};

use crate::processes::ProcessDefinition;

#[derive(Clone, Debug, Encode, Decode)]
pub struct IteratedCtCollection<T> {
    data: Vec<T>,
    num_right_thresholds: usize,
}

impl<T> IteratedCtCollection<T> {
    pub(crate) fn new(
        data: Vec<T>,
        num_left_thresholds: usize,
        num_right_thresholds: usize,
    ) -> Self {
        let expected_len = num_left_thresholds
            .checked_mul(num_right_thresholds)
            .expect("iterated threshold-counterterm dimensions overflow usize");
        assert_eq!(
            data.len(),
            expected_len,
            "iterated threshold-counterterm data must contain one entry per left/right pair"
        );
        Self {
            data,
            num_right_thresholds,
        }
    }

    pub fn map_ref<U, F>(&self, f: F) -> IteratedCtCollection<U>
    where
        F: Fn(&T) -> U,
    {
        let data = self.data.iter().map(f).collect();
        IteratedCtCollection {
            data,
            num_right_thresholds: self.num_right_thresholds,
        }
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = &T> {
        self.data.iter()
    }

    pub(crate) fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.data.iter_mut()
    }

    pub(crate) fn num_right_thresholds(&self) -> usize {
        self.num_right_thresholds
    }
}

impl<T> Index<(LeftThresholdId, RightThresholdId)> for IteratedCtCollection<T> {
    type Output = T;

    fn index(&self, index: (LeftThresholdId, RightThresholdId)) -> &Self::Output {
        let (left_id, right_id) = index;
        &self.data[left_id.0 * self.num_right_thresholds + right_id.0]
    }
}

impl<T> IndexMut<(LeftThresholdId, RightThresholdId)> for IteratedCtCollection<T> {
    fn index_mut(&mut self, index: (LeftThresholdId, RightThresholdId)) -> &mut Self::Output {
        let (left_id, right_id) = index;
        &mut self.data[left_id.0 * self.num_right_thresholds + right_id.0]
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct LUCounterTermData {
    pub left_thresholds: TiVec<LeftThresholdId, RaisedEsurfaceGroup>,
    pub right_thresholds: TiVec<RightThresholdId, RaisedEsurfaceGroup>,
    pub left_atoms: TiVec<LeftThresholdId, ParametricIntegrands>,
    pub right_atoms: TiVec<RightThresholdId, ParametricIntegrands>,
    pub iterated: IteratedCtCollection<ParametricIntegrands>,
}

fn max_dual_size_for_cut_cff_indices<'a>(
    cut_cff_indices: impl Iterator<Item = &'a CutCFFIndex>,
) -> usize {
    cut_cff_indices
        .map(|cut_cff_index| {
            shape_from_cut_cff_index(cut_cff_index)
                .map(|shape| HyperDual::<F<f64>>::new(shape).values.len())
                .unwrap_or(1)
        })
        .max()
        .unwrap_or(1)
}

define_index! {pub struct RightThresholdId;}
define_index! {pub struct LeftThresholdId;}
define_index! {pub struct CutGroupId;}

/// Eligibility of one threshold E-surface relative to one side of one physical cut.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThresholdCountertermStatus {
    /// The E-surface does not depend on any active loop coordinate of the cut-side subspace.
    NoRadialDependence,
    /// The zero can only occur at a pinched boundary and has no strictly existing region.
    AlwaysPinched,
    /// A strictly existing region is possible, but its boundary includes a pinched limit.
    CanBecomePinched,
    /// Current resolved masses and kinematics rigorously exclude a zero over the full cut phase space.
    ProvenNonExisting,
    /// A zero is possible, or the available metadata is insufficient to exclude one safely.
    PotentiallyExisting,
}

impl ThresholdCountertermStatus {
    pub fn is_eligible_for_generation(self, check_current_model: bool) -> bool {
        match self {
            Self::NoRadialDependence | Self::AlwaysPinched => false,
            Self::ProvenNonExisting => !check_current_model,
            Self::CanBecomePinched | Self::PotentiallyExisting => true,
        }
    }
}

/// Persisted topology and invariant-bound metadata for one generated cut-relative CT association.
#[derive(Clone, Debug, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ThresholdCountertermAssociation {
    pub esurface_id: EsurfaceID,
    pub cut_boundary_edges: Vec<EdgeIndex>,
    pub threshold_boundary_edges: Vec<EdgeIndex>,
    pub invariant_bound_is_applicable: bool,
}

#[derive(Clone, Debug, Default, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CutThresholdCountertermAssociations {
    pub left: Vec<ThresholdCountertermAssociation>,
    pub right: Vec<ThresholdCountertermAssociation>,
}

impl ThresholdCountertermAssociation {
    fn mass_sum(edges: &[EdgeIndex], masses: &EdgeVec<Option<F<f64>>>) -> Option<F<f64>> {
        edges.iter().try_fold(F(0.0), |sum, edge_id| {
            masses[*edge_id].as_ref().map(|mass| sum + mass.abs())
        })
    }

    fn incoming_invariant_mass(
        graph: &Graph,
        runtime_settings: &RuntimeSettings,
    ) -> Option<F<f64>> {
        let external_momenta = runtime_settings
            .kinematics
            .externals
            .get_dependent_externals(DependentMomentaConstructor::CrossSection)
            .ok()?;
        if external_momenta.is_empty()
            || external_momenta.len() != graph.get_edges_in_initial_state_cut().len()
        {
            return None;
        }

        let total = external_momenta
            .into_iter()
            .reduce(|sum, momentum| sum + momentum)?;
        let invariant_squared = total.square();
        if invariant_squared.is_nan()
            || invariant_squared.is_infinite()
            || invariant_squared <= F(0.0)
        {
            None
        } else {
            Some(invariant_squared.sqrt())
        }
    }

    fn classify_from_invariant_bounds(
        &self,
        cut_boundary_mass: &F<f64>,
        threshold_boundary_mass: &F<f64>,
        maximum_cut_invariant: Option<F<f64>>,
        tolerance: &F<f64>,
    ) -> ThresholdCountertermStatus {
        if !self.invariant_bound_is_applicable
            || cut_boundary_mass.is_nan()
            || cut_boundary_mass.is_infinite()
            || threshold_boundary_mass.is_nan()
            || threshold_boundary_mass.is_infinite()
            || tolerance.is_nan()
            || tolerance.is_infinite()
            || maximum_cut_invariant
                .as_ref()
                .is_some_and(|bound| bound.is_nan() || bound.is_infinite())
        {
            return ThresholdCountertermStatus::PotentiallyExisting;
        }

        let minimum_margin = cut_boundary_mass * cut_boundary_mass
            - threshold_boundary_mass * threshold_boundary_mass;
        let masses_match = minimum_margin.abs() <= *tolerance;
        let threshold_minimum_is_reachable = minimum_margin <= *tolerance;

        if masses_match
            && (self.cut_boundary_edges.len() <= 1 || self.threshold_boundary_edges.len() <= 1)
        {
            return ThresholdCountertermStatus::AlwaysPinched;
        }

        if self.cut_boundary_edges.len() == 1 && self.threshold_boundary_edges.len() == 1 {
            return ThresholdCountertermStatus::ProvenNonExisting;
        }

        if self.cut_boundary_edges.len() > 1 && self.threshold_boundary_edges.len() == 1 {
            let lower_margin = threshold_boundary_mass * threshold_boundary_mass
                - cut_boundary_mass * cut_boundary_mass;
            if lower_margin <= *tolerance {
                return ThresholdCountertermStatus::ProvenNonExisting;
            }
        }

        let Some(maximum_cut_invariant) = maximum_cut_invariant else {
            return if masses_match
                && self.cut_boundary_edges.len() > 1
                && self.threshold_boundary_edges.len() > 1
            {
                ThresholdCountertermStatus::CanBecomePinched
            } else {
                ThresholdCountertermStatus::PotentiallyExisting
            };
        };

        let maximum_margin = maximum_cut_invariant * maximum_cut_invariant
            - threshold_boundary_mass * threshold_boundary_mass;
        if maximum_margin.abs() <= *tolerance {
            ThresholdCountertermStatus::AlwaysPinched
        } else if maximum_margin < -*tolerance {
            ThresholdCountertermStatus::ProvenNonExisting
        } else if threshold_minimum_is_reachable
            && self.cut_boundary_edges.len() > 1
            && self.threshold_boundary_edges.len() > 1
        {
            ThresholdCountertermStatus::CanBecomePinched
        } else {
            ThresholdCountertermStatus::PotentiallyExisting
        }
    }

    pub fn classify_for_model(
        &self,
        graph: &Graph,
        cut: &CrossSectionCut,
        model: &Model,
        param_builder: &ParamBuilder,
        runtime_settings: &RuntimeSettings,
        normalized_margin_tolerance: f64,
    ) -> ThresholdCountertermStatus {
        if !self.invariant_bound_is_applicable {
            return ThresholdCountertermStatus::PotentiallyExisting;
        }

        match graph
            .classify_threshold_pinch(&self.cut_boundary_edges, &self.threshold_boundary_edges)
        {
            ThresholdPinchStatus::Always => {
                return ThresholdCountertermStatus::AlwaysPinched;
            }
            ThresholdPinchStatus::CanBecome | ThresholdPinchStatus::NotProven => {}
        }

        // Unlike `get_real_mass_vector`, retain unresolved or complex masses as unknown:
        // replacing either with zero could turn an inconclusive bound into a false rejection.
        let masses = graph.new_edgevec(|edge, _, _| {
            if matches!(edge.mass, EdgeMass::Zero) {
                Some(F(0.0))
            } else {
                edge.mass_value::<f64>(model, param_builder)
                    .and_then(|mass| IsZero::is_zero(&mass.im).then_some(mass.re))
            }
        });
        let Some(cut_boundary_mass) = Self::mass_sum(&self.cut_boundary_edges, &masses) else {
            return ThresholdCountertermStatus::PotentiallyExisting;
        };
        let Some(threshold_boundary_mass) = Self::mass_sum(&self.threshold_boundary_edges, &masses)
        else {
            return ThresholdCountertermStatus::PotentiallyExisting;
        };
        let e_cm = F(runtime_settings.kinematics.e_cm.abs());
        // The classified margins are differences of invariant masses squared, so the
        // dimensionless threshold is normalized by E_cm squared.
        let normalized_margin_tolerance = if normalized_margin_tolerance.is_finite() {
            normalized_margin_tolerance.abs()
        } else {
            DEFAULT_ESURFACE_EXISTENCE_THRESHOLD
        };
        let tolerance = F(normalized_margin_tolerance) * e_cm * e_cm;

        let maximum_cut_invariant = if self.cut_boundary_edges.len() == 1 {
            Some(cut_boundary_mass)
        } else {
            let incoming_invariant = Self::incoming_invariant_mass(graph, runtime_settings);
            let complement_mass = graph
                .iter_edges_of(&cut.cut)
                .map(|(_, edge_id, _)| edge_id)
                .filter(|edge_id| !self.cut_boundary_edges.contains(edge_id))
                .try_fold(F(0.0), |sum, edge_id| {
                    masses[edge_id].as_ref().map(|mass| sum + mass.abs())
                });
            incoming_invariant
                .zip(complement_mass)
                .and_then(|(invariant, complement_mass)| {
                    let maximum = invariant - complement_mass;
                    (maximum > F(0.0)).then_some(maximum)
                })
        };
        self.classify_from_invariant_bounds(
            &cut_boundary_mass,
            &threshold_boundary_mass,
            maximum_cut_invariant,
            &tolerance,
        )
    }
}

#[derive(Clone)]
struct TopologicalThresholdCandidate {
    left: SuBitGraph,
    cut: OrientedCut,
    right: SuBitGraph,
    esurface_id: EsurfaceID,
}

use derive_more::{From, Into};
#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSection {
    pub name: String,
    pub integrand: Option<ProcessIntegrand>,
    pub supergraphs: Vec<CrossSectionGraph>,
    pub external_particles: Vec<ArcParticle>,
    pub external_connections: Vec<ExternalConnection>,
    pub n_incmoming: usize,
    pub graph_group_structure: TiVec<GroupId, GraphGroup>,
}

impl CrossSection {
    pub fn plan_graph_group_selection(
        &self,
        spec: &GraphGroupSelectionSpec,
    ) -> Result<GraphGroupSelectionPlan> {
        if spec.mode() == GraphGroupSelectionMode::CrossSectionAmplitudeGraphs {
            return Err(eyre!(
                "`select --amplitude-graphs` for cross sections requires process and generation settings context."
            ));
        }
        if spec.has_raised_cut_rules() {
            return Err(eyre!(
                "Raised-cut signature selection for cross sections requires process and generation settings context."
            ));
        }
        spec.plan(&self.graph_group_structure, |graph_id| {
            self.supergraphs.get(graph_id).map(|graph| &graph.graph)
        })
    }

    pub fn plan_graph_group_selection_with_context(
        &self,
        spec: &GraphGroupSelectionSpec,
        model: &Model,
        process_definition: &ProcessDefinition,
        generation_settings: &GenerationSettings,
    ) -> Result<GraphGroupSelectionPlan> {
        match spec.mode() {
            GraphGroupSelectionMode::MasterGraphs => {
                spec.plan_with_analysis_contexts(
                    &self.graph_group_structure,
                    |graph_id| self.supergraphs.get(graph_id).map(|graph| &graph.graph),
                    |_master_graph_id, master_graph| {
                        Ok(vec![GraphSelectionSubject::whole_graph(master_graph)])
                    },
                    |master_graph_id, master_graph| {
                        self.cut_selection_subjects(
                            master_graph_id,
                            master_graph,
                            model,
                            process_definition,
                            generation_settings,
                        )
                    },
                    "Graph-group selection structural filters have no graph analysis subjects.",
                    "Raised-cut signature selection found no process-valid Cutkosky cuts to analyse.",
                )
            }
            GraphGroupSelectionMode::CrossSectionAmplitudeGraphs => {
                spec.plan_with_analysis_contexts(
                    &self.graph_group_structure,
                    |graph_id| self.supergraphs.get(graph_id).map(|graph| &graph.graph),
                    |master_graph_id, master_graph| {
                        let cross_section_graph =
                            self.supergraphs.get(master_graph_id).ok_or_else(|| {
                                eyre!(
                                    "Graph group refers to missing master supergraph id {}.",
                                    master_graph_id
                                )
                            })?;
                        let cuts = cross_section_graph.process_valid_cuts(
                            model,
                            process_definition,
                            generation_settings,
                        )?;
                        Ok(cuts
                            .into_iter()
                            .flat_map(|cut| cut.amplitude_side_subjects(master_graph))
                            .collect())
                    },
                    |master_graph_id, master_graph| {
                        self.cut_selection_subjects(
                            master_graph_id,
                            master_graph,
                            model,
                            process_definition,
                            generation_settings,
                        )
                    },
                    "`select --amplitude-graphs` found no process-valid Cutkosky cuts to analyse for structural filters.",
                    "Raised-cut signature selection found no process-valid Cutkosky cuts to analyse.",
                )
            }
        }
    }

    fn cut_selection_subjects<'a>(
        &'a self,
        master_graph_id: usize,
        master_graph: &'a Graph,
        model: &Model,
        process_definition: &ProcessDefinition,
        generation_settings: &GenerationSettings,
    ) -> Result<Vec<GraphCutSelectionSubject<'a>>> {
        let cross_section_graph = self.supergraphs.get(master_graph_id).ok_or_else(|| {
            eyre!(
                "Graph group refers to missing master supergraph id {}.",
                master_graph_id
            )
        })?;
        let cuts = cross_section_graph.process_valid_cuts(
            model,
            process_definition,
            generation_settings,
        )?;
        Ok(cuts
            .into_iter()
            .map(|cut| GraphCutSelectionSubject::new(master_graph, cut.cut.as_subgraph()))
            .collect())
    }

    pub fn amplitude_graph_signature_inventory(
        &self,
        model: &Model,
        process_definition: &ProcessDefinition,
        generation_settings: &GenerationSettings,
    ) -> Result<GraphSelectionSignatureInventory> {
        let subjects = self
            .graph_group_structure
            .iter()
            .flat_map(|group| group.into_iter().next())
            .map(|master_graph_id| {
                let cross_section_graph =
                    self.supergraphs.get(master_graph_id).ok_or_else(|| {
                        eyre!(
                            "Graph group refers to missing master supergraph id {}.",
                            master_graph_id
                        )
                    })?;
                let cuts = cross_section_graph.process_valid_cuts(
                    model,
                    process_definition,
                    generation_settings,
                )?;
                Ok(cuts
                    .into_iter()
                    .flat_map(|cut| cut.amplitude_side_subjects(&cross_section_graph.graph))
                    .collect::<Vec<_>>())
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();

        Ok(GraphSelectionSignatureInventory::from_analysis_subjects(
            subjects,
        ))
    }

    pub fn raised_cut_signature_inventory(
        &self,
        model: &Model,
        process_definition: &ProcessDefinition,
        generation_settings: &GenerationSettings,
    ) -> Result<RaisedCutSignatureInventory> {
        let subjects = self
            .graph_group_structure
            .iter()
            .flat_map(|group| group.into_iter().next())
            .map(|master_graph_id| {
                let cross_section_graph =
                    self.supergraphs.get(master_graph_id).ok_or_else(|| {
                        eyre!(
                            "Graph group refers to missing master supergraph id {}.",
                            master_graph_id
                        )
                    })?;
                let cuts = cross_section_graph.process_valid_cuts(
                    model,
                    process_definition,
                    generation_settings,
                )?;
                Ok(cuts
                    .into_iter()
                    .map(|cut| {
                        GraphCutSelectionSubject::new(
                            &cross_section_graph.graph,
                            cut.cut.as_subgraph(),
                        )
                    })
                    .collect::<Vec<_>>())
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();

        Ok(RaisedCutSignatureInventory::from_cut_subjects(subjects))
    }

    pub fn validate_graph_group_selection_plan(
        &self,
        plan: &GraphGroupSelectionPlan,
    ) -> Result<()> {
        for &old_group_id in plan.retained_group_ids() {
            plan.new_group_id_for_old(old_group_id).ok_or_else(|| {
                eyre!(
                    "Selection plan is missing compact group id for old group {}.",
                    old_group_id.0
                )
            })?;
            if old_group_id.0 >= self.graph_group_structure.len() {
                return Err(eyre!(
                    "Selection plan refers to missing graph group {}.",
                    old_group_id.0
                ));
            }
            let group = &self.graph_group_structure[old_group_id];
            for old_graph_id in group {
                if old_graph_id >= self.supergraphs.len() {
                    return Err(eyre!(
                        "Graph group {} refers to missing graph id {}.",
                        old_group_id.0,
                        old_graph_id
                    ));
                }
            }
            let master = group.master();
            if master >= self.supergraphs.len() {
                return Err(eyre!(
                    "Graph group {} refers to missing master graph id {}.",
                    old_group_id.0,
                    master
                ));
            }
        }

        Ok(())
    }

    pub fn apply_graph_group_selection(&mut self, plan: &GraphGroupSelectionPlan) -> Result<()> {
        self.validate_graph_group_selection_plan(plan)?;

        let mut old_graph_to_new_group = vec![None; self.supergraphs.len()];
        let mut old_graph_is_master = vec![false; self.supergraphs.len()];
        for &old_group_id in plan.retained_group_ids() {
            let new_group_id = plan.new_group_id_for_old(old_group_id).ok_or_else(|| {
                eyre!(
                    "Selection plan is missing compact group id for old group {}.",
                    old_group_id.0
                )
            })?;
            let group = &self.graph_group_structure[old_group_id];
            for old_graph_id in group {
                old_graph_to_new_group[old_graph_id] = Some(new_group_id);
            }
            old_graph_is_master[group.master()] = true;
        }

        let mut new_supergraphs = self
            .supergraphs
            .iter()
            .cloned()
            .enumerate()
            .filter_map(|(old_graph_id, mut graph)| {
                let new_group_id = old_graph_to_new_group[old_graph_id]?;
                graph.graph.group_id = Some(new_group_id);
                graph.graph.is_group_master = old_graph_is_master[old_graph_id];
                if let Some(multi_channeling_setup) = &mut graph.derived_data.multi_channeling_setup
                {
                    multi_channeling_setup.graph.group_id = Some(new_group_id);
                    multi_channeling_setup.graph.is_group_master =
                        old_graph_is_master[old_graph_id];
                }
                Some(graph)
            })
            .collect::<Vec<_>>();

        let mut parsed_graphs = new_supergraphs
            .iter()
            .map(|graph| graph.graph.clone())
            .collect::<Vec<_>>();
        let new_graph_group_structure = complete_group_parsing(&mut parsed_graphs)?;
        for (graph, parsed_graph) in new_supergraphs.iter_mut().zip(parsed_graphs) {
            graph.graph.group_id = parsed_graph.group_id;
            graph.graph.is_group_master = parsed_graph.is_group_master;
            if let Some(multi_channeling_setup) = &mut graph.derived_data.multi_channeling_setup {
                multi_channeling_setup.graph.group_id = graph.graph.group_id;
                multi_channeling_setup.graph.is_group_master = graph.graph.is_group_master;
            }
        }

        self.supergraphs = new_supergraphs;
        self.graph_group_structure = new_graph_group_structure;
        self.integrand = None;
        Ok(())
    }

    fn storage_path(&self, base: &Path) -> PathBuf {
        base.join(&self.name)
    }

    pub fn export_standalone(
        &self,
        path: impl AsRef<Path>,
        settings: &crate::processes::StandaloneExportSettings,
    ) -> Result<()> {
        if let Some(integrand) = &self.integrand {
            integrand.export_standalone(path, settings)?;
        } else {
            return Err(eyre!(
                "Cannot export standalone cross section {} without integrand",
                self.name
            ));
        }

        Ok(())
    }

    #[allow(dead_code)]
    pub(crate) fn write_dot<W: std::io::Write>(
        &self,
        writer: &mut W,
        settings: &DotExportSettings,
    ) -> Result<(), std::io::Error> {
        for graph in &self.supergraphs {
            graph.write_dot(writer, settings)?;
            writeln!(writer)?;
        }
        Ok(())
    }

    pub fn write_dot_fmt<W: std::fmt::Write>(
        &self,
        writer: &mut W,
        settings: &DotExportSettings,
    ) -> Result<(), std::fmt::Error> {
        for graph in &self.supergraphs {
            graph.write_dot_fmt(writer, settings)?;
            writeln!(writer)?;
        }
        Ok(())
    }

    pub(crate) fn new(name: String) -> Self {
        Self {
            name,
            integrand: None,
            supergraphs: vec![],
            external_connections: vec![],
            external_particles: vec![],
            n_incmoming: 0,
            graph_group_structure: TiVec::new(),
        }
    }

    pub fn from_graph_list(name: String, mut graphs: Vec<Graph>, _model: &Model) -> Result<Self> {
        let mut cross_section = CrossSection::new(name);
        cross_section.graph_group_structure = complete_group_parsing(&mut graphs)?;
        // println!("group structure: {:?}", cross_section.graph_group_structure);

        for cross_section_graph in graphs {
            // cross_section_graph.param_builder =  ParamBuilder::new(&cross_section_graph, model);
            cross_section.add_supergraph(cross_section_graph)?;
        }
        Ok(cross_section)
    }

    pub(crate) fn warm_up(&mut self, model: &Model) -> Result<()> {
        if let Some(integrand) = &mut self.integrand {
            integrand.warm_up(model)
        } else {
            Err(eyre!(
                "Cannot warm up amplitude {} without integrand",
                self.name
            ))
        }
    }

    fn add_supergraph(&mut self, supergraph: Graph) -> Result<()> {
        if self.external_particles.is_empty() {
            let external_particles = supergraph.get_external_partcles();
            if !external_particles.len().is_multiple_of(2) {
                return Err(eyre!(
                    "expected even number of externals for forward scattering graph"
                ));
            }
            self.external_particles = external_particles;
            self.n_incmoming = self.external_particles.len() / 2;
        } else if self.external_particles != supergraph.get_external_partcles() {
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
        global_settings: &GenerationSettings,
        runtime_default: LockedRuntimeSettings,
        generation_pool: &ThreadPool,
    ) -> Result<Vec<NamedGraphGenerationReport>> {
        let integrand_name = self.name.clone();
        generation_progress::begin_phase(
            GenerationProgressPhase::GraphPreprocessing,
            GenerationProcessKind::CrossSection,
            &process_definition.folder_name,
            &integrand_name,
            self.supergraphs.len(),
            None,
        );
        generation_pool.install(|| {
            self.supergraphs
                .par_iter_mut()
                .map(|supergraph| {
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    let graph_name = supergraph.graph.name.clone();
                    generation_progress::graph_started(
                        GenerationProcessKind::CrossSection,
                        &integrand_name,
                        &graph_name,
                        None,
                    );
                    let stats = supergraph.preprocess(
                        model,
                        process_definition,
                        global_settings,
                        runtime_default,
                    )?;
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    generation_progress::graph_finished(
                        GenerationProcessKind::CrossSection,
                        &integrand_name,
                        &graph_name,
                        &stats,
                        None,
                    );
                    Ok(NamedGraphGenerationReport {
                        integrand_name: integrand_name.clone(),
                        graph_name,
                        stats,
                    })
                })
                .collect::<Result<Vec<_>>>()
        })
    }

    pub fn build_integrand(
        &mut self,
        model: &Model,
        process_name: &str,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
        generation_pool: &ThreadPool,
    ) -> Result<Vec<NamedGraphGenerationReport>> {
        let started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #graph, #summary;
            stage = "cross_section_build_integrand_start",
            integrand = %self.name,
            graph_count = self.supergraphs.len(),
            "Generation timing milestone"
        );
        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        let integrand_name = self.name.clone();
        let total_cuts = self
            .supergraphs
            .iter()
            .map(|sg| {
                sg.derived_data
                    .cut_group_data
                    .cut_groups
                    .iter()
                    .map(|cut_group| cut_group.cuts.len())
                    .sum::<usize>()
            })
            .sum();
        generation_progress::begin_phase(
            GenerationProgressPhase::GraphGeneration,
            GenerationProcessKind::CrossSection,
            process_name,
            &integrand_name,
            self.supergraphs.len(),
            Some(total_cuts),
        );
        let mut graph_reports = Vec::new();
        let terms = generation_pool.install(|| {
            self.supergraphs
                .par_iter_mut()
                .map(|sg| {
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    let graph_started = std::time::Instant::now();
                    crate::debug_tags!(#generation, #profile, #graph, #summary;
                        stage = "generate_term_for_graph_start",
                        integrand = %integrand_name,
                        graph = %sg.graph.name,
                        "Generation timing milestone"
                    );
                    generation_progress::graph_started(
                        GenerationProcessKind::CrossSection,
                        &integrand_name,
                        &sg.graph.name,
                        Some(sg.cuts.len()),
                    );
                    let _progress_context_guard = generation_progress::enter_progress_context(
                        format!("{} / {} cuts", sg.graph.name, sg.cuts.len()),
                    );
                    let (term, mut stats) = sg.generate_term_for_graph(model, global_settings)?;
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    stats.total_time += graph_started.elapsed();
                    crate::debug_tags!(#generation, #profile, #graph, #summary;
                        stage = "generate_term_for_graph_done",
                        integrand = %integrand_name,
                        graph = %sg.graph.name,
                        elapsed_ms = graph_started.elapsed().as_secs_f64() * 1000.0,
                        "Generation timing milestone"
                    );
                    generation_progress::graph_finished(
                        GenerationProcessKind::CrossSection,
                        &integrand_name,
                        &sg.graph.name,
                        &stats,
                        None,
                    );
                    Ok((
                        term,
                        NamedGraphGenerationReport {
                            integrand_name: integrand_name.clone(),
                            graph_name: sg.graph.name.clone(),
                            stats,
                        },
                    ))
                })
                .collect::<Result<Vec<_>>>()
        })?;
        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        for (_, report) in &terms {
            graph_reports.push(report.clone());
        }
        let mut terms = terms.into_iter().map(|(term, _)| term).collect::<Vec<_>>();

        for group in self.graph_group_structure.iter() {
            let master = group.master();
            let mc_of_master = self.supergraphs[master]
                .derived_data
                .multi_channeling_setup
                .as_ref()
                .unwrap();

            for graph_id in group.into_iter() {
                terms[graph_id].multi_channeling_setup = mc_of_master.clone();
            }
        }

        let backend_started = std::time::Instant::now();
        let graph_count = terms.len();
        crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
            stage = "prepare_runtime_backends_start",
            integrand = %self.name,
            graph_count,
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );
        generation_progress::backend_started(
            GenerationProcessKind::CrossSection,
            &self.name,
            graph_count,
        );
        let mut cross_section_integrand = CrossSectionIntegrand {
            settings: runtime_default.into(),
            data: CrossSectionIntegrandData {
                compilation: global_settings
                    .generation
                    .compile
                    .frozen_mode(&global_settings.generation.evaluator),
                loop_cache_id: 0,
                external_cache_id: 0,
                base_external_cache_id: 0,
                rotations: None,
                name: self.name.clone(),
                external_connections: self.external_connections.clone(),
                n_incoming: self.n_incmoming,
                // polarizations,
                graph_terms: terms,
                graph_group_structure: self.graph_group_structure.clone(),
                graph_to_group_id: graph_to_group_id_for_group_structure(
                    &self.graph_group_structure,
                ),
            },
            event_processing_runtime: Default::default(),
            active_f64_backend: Default::default(),
        };
        let compile_times = cross_section_integrand
            .prepare_runtime_backends_after_generation_with_compile_times()?;
        crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
            stage = "prepare_runtime_backends_done",
            integrand = %self.name,
            graph_count,
            elapsed_ms = backend_started.elapsed().as_secs_f64() * 1000.0,
            total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );
        generation_progress::backend_finished(
            GenerationProcessKind::CrossSection,
            &self.name,
            backend_started.elapsed(),
        );
        for (report, compile_time) in graph_reports.iter_mut().zip(compile_times) {
            report.stats.evaluator_compile_time += compile_time;
            report.stats.total_time += compile_time;
        }

        self.integrand = Some(ProcessIntegrand::CrossSection(cross_section_integrand));
        crate::debug_tags!(#generation, #profile, #graph, #summary;
            stage = "cross_section_build_integrand_done",
            integrand = %self.name,
            graph_count = self.supergraphs.len(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );
        Ok(graph_reports)
    }

    pub fn compile(
        &mut self,
        path: impl AsRef<Path>,
        override_existing: bool,
        thread_pool: &ThreadPool,
    ) -> Result<Vec<NamedGraphGenerationReport>> {
        info!("Compiling cross section {}", self.name);
        let p = self.storage_path(path.as_ref());

        let result = fs::create_dir_all(&p).with_context(|| {
            format!(
                "Trying to create directory to compile cross section {}",
                p.display()
            )
        });

        if override_existing {
            result?;
        }

        if let Some(integrand) = &mut self.integrand {
            let compile_times = integrand.compile(&p, override_existing, thread_pool)?;
            return Ok(compile_times
                .into_iter()
                .map(|(graph_name, duration)| NamedGraphGenerationReport {
                    integrand_name: self.name.clone(),
                    graph_name,
                    stats: GraphGenerationStats {
                        total_time: duration,
                        evaluator_compile_time: duration,
                        ..GraphGenerationStats::default()
                    },
                })
                .collect());
        }
        Ok(Vec::new())
    }

    pub fn save(&mut self, path: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        let p = self.storage_path(path.as_ref());
        let r = fs::create_dir_all(&p).with_context(|| {
            format!(
                "Trying to create directory to save cross section {}",
                p.display()
            )
        });

        if override_existing {
            r?;
        }

        let integrand = self.integrand.take();
        if let Some(integrand) = &integrand {
            integrand.save(&p, override_existing)?;
        }

        let binary = bincode::encode_to_vec(&(*self), bincode::config::standard())?;
        if override_existing {
            fs::write(p.join("cs.bin"), &binary)?;
        } else {
            let mut file = File::create_new(p.join("cs.bin"))?;
            file.write_all(&binary)?;
        }

        self.integrand = integrand;
        Ok(())
    }

    pub(crate) fn load(path: impl AsRef<Path>, context: GammaLoopContextContainer) -> Result<Self> {
        let binary = fs::read(path.as_ref().join("cs.bin"))?;
        let (mut cs, _): (Self, _) =
            bincode::decode_from_slice_with_context(&binary, bincode::config::standard(), context)?;

        if path.as_ref().join("integrand").exists() {
            let integrand = CrossSectionIntegrand::load(path.as_ref().join("integrand"), context)?;
            cs.integrand = Some(ProcessIntegrand::CrossSection(integrand));
        }

        Ok(cs)
    }
}

#[derive(Clone, bincode::Encode, bincode::Decode)]
pub struct CrossSectionCut {
    pub cut: OrientedCut,
    pub left: SuBitGraph,
    pub right: SuBitGraph,
}

impl CrossSectionCut {
    fn amplitude_side_subjects<'a>(&self, graph: &'a Graph) -> [GraphSelectionSubject<'a>; 2] {
        let cut_edges = self.cut.as_subgraph();
        [
            GraphSelectionSubject::cut_side_amplitude_subgraph(
                graph,
                self.left.subtract(&cut_edges),
            ),
            GraphSelectionSubject::cut_side_amplitude_subgraph(
                graph,
                self.right.subtract(&cut_edges),
            ),
        ]
    }

    pub(crate) fn is_s_channel(&self, cross_section_graph: &CrossSectionGraph) -> Result<bool> {
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

    pub(crate) fn is_valid_for_process(
        &self,
        cross_section_graph: &CrossSectionGraph,
        process: &ProcessDefinition,
        model: &Model,
    ) -> Result<bool> {
        if self.is_s_channel(cross_section_graph)? {
            let cut_content_builder = self
                .cut
                .iter_edges(&cross_section_graph.graph.underlying)
                .filter_map(|(orientation, edge_data)| {
                    Some(if orientation == Orientation::Reversed {
                        edge_data.data.particle()?.get_anti_particle(model)
                    } else {
                        edge_data.data.particle()?.clone()
                    })
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
                        cut_content.iter().map(|p| p.name.clone()).collect_vec()
                    );

                    for particle in particle_content {
                        if let Some(index) = cut_content.iter().position(|p| p == &particle) {
                            cut_content.remove(index);
                        } else {
                            debug!("wrong particles");
                            return false;
                        }
                    }

                    let (n_unresolved, unresolved_cut_content) =
                        process.unresolved_cut_content(model);

                    if cut_content.len() > n_unresolved {
                        debug!(" too many unresolved particles");
                        return false;
                    }

                    if !cut_content
                        .iter()
                        .all(|particle| unresolved_cut_content.contains(particle))
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
pub struct CrossSectionGraph {
    pub graph: Graph,
    pub source_nodes: HedgeNode,
    pub target_nodes: HedgeNode,
    pub cuts: TiVec<CutId, CrossSectionCut>,
    pub cut_esurface: TiVec<CutId, Esurface>,
    pub cut_esurface_id_map: TiVec<CutId, EsurfaceID>,
    pub derived_data: CrossSectionDerivedData,
}

impl CrossSectionGraph {
    pub(crate) fn new(graph: Graph) -> Self {
        let (source_node, target_node) = graph.get_source_and_target();

        Self {
            graph,
            source_nodes: source_node,
            target_nodes: target_node,
            cuts: TiVec::new(),
            cut_esurface: TiVec::new(),
            cut_esurface_id_map: TiVec::new(),
            derived_data: CrossSectionDerivedData::new_empty(),
        }
    }

    pub(crate) fn apply_spin_sum(
        &mut self,
        model: &Model,
        generation_settings: &GenerationSettings,
        locked_runtime_settings: &LockedRuntimeSettings,
    ) -> Result<()> {
        for (extid, hel) in locked_runtime_settings.helicities().iter().enumerate() {
            //println!("{extid},{hel:?}");
            let eid = self.graph.loop_momentum_basis.ext_edges[ExternalIndex(extid)];

            let Some(p) = self.graph.underlying[eid].particle() else {
                continue;
            };

            match hel {
                Helicity::Summed => {
                    let Some(p) = p.polarization_sum(
                        eid,
                        false,
                        generation_settings.vector_polarization_sum_gauge,
                    )?
                    else {
                        continue;
                    };
                    self.graph.global_prefactor.projector =
                        self.graph.global_prefactor.projector.replace_multiple(&[p]);
                }
                Helicity::SummedAveraged => {
                    let Some(p) = p.polarization_sum(
                        eid,
                        true,
                        generation_settings.vector_polarization_sum_gauge,
                    )?
                    else {
                        continue;
                    };
                    self.graph.global_prefactor.projector =
                        self.graph.global_prefactor.projector.replace_multiple(&[p]);
                }
                _ => {}
            }
        }

        self.graph.polarizations = self.graph.global_prefactor.polarizations();
        self.graph.param_builder = ParamBuilder::new(
            &self.graph,
            model,
            &self.graph.loop_momentum_basis,
            &self.graph.param_builder.pairs.additional_params.params,
        );
        Ok(())
    }

    pub(crate) fn preprocess(
        &mut self,
        model: &Model,
        process_definition: &ProcessDefinition,
        settings: &GenerationSettings,
        runtime_default: LockedRuntimeSettings,
    ) -> Result<GraphGenerationStats> {
        let preprocess_started = std::time::Instant::now();
        let mut stats = GraphGenerationStats::default();
        self.apply_spin_sum(model, settings, &runtime_default)?;
        debug_tags!(#generation; "generating cuts");
        self.generate_cuts(model, process_definition, settings)?;
        debug_tags!(#generation; "generating esurfaces corresponding to cuts");
        self.generate_esurface_cuts();
        debug_tags!(#generation; "generating cff");
        stats.merge_in_place(&self.generate_cff(settings)?);
        debug_tags!(#generation; "building lmbs");
        self.build_lmbs()?;
        debug_tags!(#generation; "building multi channeling channels");

        if self.graph.is_group_master {
            self.build_multi_channeling_channels(settings.override_lmb_heuristics)?;
        }

        let vk = crate::utils::vakint()?;
        debug_tags!(#generation; "building parametric integrand");
        self.build_parametric_integrand(settings, vk)?;
        //self.build_parametric_integrand_cut_groups(settings)?;

        let threshold_candidates = self.topological_threshold_candidates()?;
        self.derived_data.threshold_candidate_esurface_ids = threshold_candidates
            .iter()
            .map(|candidate| candidate.esurface_id)
            .sorted()
            .dedup()
            .collect();
        self.derived_data.cut_threshold_associations =
            ti_vec![CutThresholdCountertermAssociations::default(); self.cuts.len()];

        if settings.threshold_subtraction.enable_thresholds {
            debug_tags!(#generation, #subtraction; "building threshold counterterm");
            self.build_subspace_data()?;
            let runtime_settings: RuntimeSettings = runtime_default.into();
            self.build_threshold_counterterm(
                model,
                settings,
                &runtime_settings,
                &threshold_candidates,
                vk,
            )?;
        }

        stats.total_time += preprocess_started.elapsed();
        Ok(stats)
    }

    #[allow(dead_code)]
    pub(crate) fn write_dot<W: std::io::Write>(
        &self,
        writer: &mut W,
        settings: &DotExportSettings,
    ) -> Result<(), std::io::Error> {
        self.graph.dot_serialize_io(writer, settings)
    }

    pub(crate) fn write_dot_fmt<W: std::fmt::Write>(
        &self,
        writer: &mut W,
        settings: &DotExportSettings,
    ) -> Result<(), std::fmt::Error> {
        self.graph.dot_serialize_fmt(writer, settings)
    }

    fn generate_cff(&mut self, settings: &GenerationSettings) -> Result<GraphGenerationStats> {
        let canonize_esurface = self
            .graph
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let contract_edges = self
            .graph
            .iter_edges_of(
                &self
                    .graph
                    .tree_edges
                    .subtract(&self.graph.initial_state_cut),
            )
            .map(|x| x.1)
            .collect_vec();

        let global_cff = self.graph.generate_cff(
            &contract_edges,
            &canonize_esurface,
            &settings.orientation_pattern,
        )?;

        let cut_esurface_map = self
            .cut_esurface
            .iter()
            .map(|esurface| {
                if let Some(pos) = self
                    .graph
                    .surface_cache
                    .esurface_cache
                    .iter()
                    .position(|e_sf| e_sf == esurface)
                    .map(Into::<EsurfaceID>::into)
                {
                    pos
                } else {
                    let pos = self.graph.surface_cache.esurface_cache.len();
                    self.graph
                        .surface_cache
                        .esurface_cache
                        .push(esurface.clone());
                    EsurfaceID(pos)
                }
            })
            .collect();

        self.cut_esurface_id_map = cut_esurface_map;

        let esurface_raised_data = self
            .graph
            .determine_raised_esurfaces_from_expression(&global_cff);

        let (cut_group_data, cut_group_stats) = CutGroupData::new_from_esurface(
            &esurface_raised_data,
            &self.cut_esurface_id_map,
            &settings.evaluator,
        );

        self.derived_data.global_cff_expression = Some(global_cff);
        self.derived_data.cut_group_data = cut_group_data;

        Ok(cut_group_stats)
    }

    pub(crate) fn process_valid_cuts(
        &self,
        model: &Model,
        process_definition: &ProcessDefinition,
        settings: &GenerationSettings,
    ) -> Result<TiVec<CutId, CrossSectionCut>> {
        if !self.cuts.is_empty() {
            return Ok(self.cuts.clone());
        }
        self.compute_process_valid_cuts(model, process_definition, settings)
            .map(|(_, cuts)| cuts)
    }

    fn compute_process_valid_cuts(
        &self,
        model: &Model,
        process_definition: &ProcessDefinition,
        settings: &GenerationSettings,
    ) -> Result<(usize, TiVec<CutId, CrossSectionCut>)> {
        let all_st_cuts = self.graph.all_st_cuts_for_cs(
            self.source_nodes.clone(),
            self.target_nodes.clone(),
            &self.graph.get_initial_state_tree().0,
        );
        let num_st_cuts = all_st_cuts.len();

        let mut cuts: TiVec<CutId, CrossSectionCut> = all_st_cuts
            .into_iter()
            .map(|(left, cut, right)| CrossSectionCut { cut, left, right })
            .filter(|cut| cut.cut.nedges(&self.graph) > 1)
            .filter_map(
                |cut| match cut.is_valid_for_process(self, process_definition, model) {
                    Ok(true) => Some(Ok(cut)),
                    Ok(false) => None,
                    Err(e) => Some(Err(e)),
                },
            )
            .collect::<Result<_>>()?;

        cuts.sort_by(|a, b| a.cut.cmp(&b.cut));

        if !settings.force_cuts.is_empty() {
            let force_cuts_sorted = settings
                .force_cuts
                .clone()
                .into_iter()
                .map(|cut_edges| cut_edges.into_iter().sorted().collect_vec())
                .collect_vec();

            cuts.retain(|cut| {
                let edges_in_cut = self
                    .graph
                    .iter_edges_of(&cut.cut)
                    .map(|(_, _, e)| e.data.name.value.clone())
                    .sorted()
                    .collect_vec();

                force_cuts_sorted.contains(&edges_in_cut)
            });
        }

        Ok((num_st_cuts, cuts))
    }

    fn generate_cuts(
        &mut self,
        model: &Model,
        process_definition: &ProcessDefinition,
        settings: &GenerationSettings,
    ) -> Result<()> {
        debug_tags!(#generation, #profile, #graph;
            stage = "cross_section_generate_cuts_start",
            graph = %self.graph.name,
            "Cut discovery timing milestone"
        );
        let started = std::time::Instant::now();
        let (num_st_cuts, cuts) =
            self.compute_process_valid_cuts(model, process_definition, settings)?;
        self.cuts = cuts;
        debug_tags!(#generation, #profile, #graph;
            stage = "cross_section_generate_cuts_done",
            graph = %self.graph.name,
            st_cut_count = num_st_cuts,
            valid_cut_count = self.cuts.len(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Cut discovery timing milestone"
        );
        generation_progress::cuts_discovered("", &self.graph.name, num_st_cuts, self.cuts.len());

        Ok(())
    }

    fn generate_esurface_cuts(&mut self) {
        debug!("generating esurfaces for cuts");

        let esurfaces: TiVec<CutId, Esurface> = self
            .cuts
            .iter()
            .map(|cut| {
                Esurface::new_from_cut_left(
                    &self.graph.underlying,
                    cut,
                    Some(&self.graph.initial_state_cut),
                )
            })
            .collect();

        debug!("generated esurfaces {:?}", esurfaces);

        self.cut_esurface = esurfaces;
    }

    pub(crate) fn build_parametric_integrand(
        &mut self,
        settings: &GenerationSettings,
        vakint: &Vakint,
    ) -> Result<()> {
        self.derived_data.cut_paramatric_integrand = self.build_integrand(settings, vakint)?;
        Ok(())
    }

    fn build_integrand(
        &mut self,
        settings: &GenerationSettings,
        vakint: &Vakint,
    ) -> Result<TiVec<CutGroupId, ParametricIntegrands>> {
        let started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "supergraph_build_integrand_start",
            graph = %self.graph.name,
            subtract_uv = settings.uv.subtract_uv,
            generate_integrated = settings.uv.generate_integrated,
            only_integrated = settings.uv.only_integrated,
            "Generation timing milestone"
        );
        let max_order = self
            .derived_data
            .cut_group_data
            .cut_groups
            .iter()
            .map(|cut_group| cut_group.related_esurface_group.max_occurence)
            .max()
            .unwrap();

        self.graph.param_builder.initialize_duals(max_order);
        crate::debug_tags!(#generation, #profile, #graph, #summary;
            stage = "supergraph_initialize_duals_done",
            graph = %self.graph.name,
            max_order,
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );

        let cuts = self
            .derived_data
            .cut_group_data
            .cut_groups
            .iter()
            .map(|cuts| CutSet {
                residue_selector: ResidueSelector {
                    lu_cut: Some(cuts.related_esurface_group.clone()),
                    left_th_cut: None,
                    right_th_cut: None,
                },
                union: cuts
                    .cuts
                    .iter()
                    .map(|cut_id| self.cuts[*cut_id].cut.as_subgraph())
                    .reduce(|cut_1, cut_2| cut_1.union(&cut_2))
                    .unwrap_or_else(|| self.graph.empty_subgraph()),
                canonicalize_external_shifts: false,
            })
            .collect();

        let cut_structure = CutStructure { cuts };
        crate::debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "supergraph_cutsets_done",
            graph = %self.graph.name,
            cut_count = cut_structure.cuts.len(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );

        let cut_woods_started = std::time::Instant::now();
        let cut_woods = CutWoods::new(cut_structure, &self.graph, &settings.uv);
        crate::debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "supergraph_cut_woods_done",
            graph = %self.graph.name,
            cut_count = cut_woods.cuts.cuts.len(),
            wood_count = cut_woods.woods.len(),
            elapsed_ms = cut_woods_started.elapsed().as_secs_f64() * 1000.0,
            total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );
        let valid_orientations: Vec<_> = self
            .derived_data
            .global_cff_expression
            .as_ref()
            .expect("global_cff_expression should have been created")
            .orientations
            .iter()
            .map(|orientation| orientation.data.orientation.clone())
            .collect();

        let lu_prefactor = self.lu_prefactor_helper();

        let unfold_started = std::time::Instant::now();
        let mut cut_forests = cut_woods.unfold(&self.graph);
        crate::debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "supergraph_cut_forests_unfold_done",
            graph = %self.graph.name,
            forest_count = cut_forests.forests.len(),
            elapsed_ms = unfold_started.elapsed().as_secs_f64() * 1000.0,
            total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );
        let forests_started = std::time::Instant::now();
        cut_forests.compute(
            &mut self.graph,
            vakint,
            &valid_orientations,
            &settings.uv,
            &settings.orientation_pattern,
        )?;
        crate::debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "supergraph_cut_forests_compute_done",
            graph = %self.graph.name,
            elapsed_ms = forests_started.elapsed().as_secs_f64() * 1000.0,
            total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );

        let orientation_started = std::time::Instant::now();
        let parametric_integrands =
            cut_forests.orientation_parametric_exprs(&self.graph, &settings.uv)?;
        crate::debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "supergraph_orientation_parametric_exprs_done",
            graph = %self.graph.name,
            parametric_integrand_count = parametric_integrands.len(),
            elapsed_ms = orientation_started.elapsed().as_secs_f64() * 1000.0,
            total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );

        let finalize_started = std::time::Instant::now();
        let result = parametric_integrands
            .into_iter()
            .map(|integrand| integrand.map(|a| a * &lu_prefactor))
            .collect();
        crate::debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "supergraph_build_integrand_done",
            graph = %self.graph.name,
            elapsed_ms = finalize_started.elapsed().as_secs_f64() * 1000.0,
            total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );
        Ok(result)
    }

    fn lu_prefactor_helper(&self) -> Atom {
        let loop_number = self.graph.cyclotomatic_number(&self.graph.full_filter())
            - self.graph.initial_state_cut.nedges(&self.graph);

        let loop_3 = loop_number as i64 * 3;
        let energy_conservation_delta_factor = Atom::num(2) * Atom::var(GS.pi);

        crate::debug_tags!(#generation, #normalization, #lu, #graph, #summary;
            stage = "cross_section_lu_prefactor",
            graph = %self.graph.name,
            loop_number,
            loop_3,
            "Cross-section LU prefactor normalization"
        );

        let tstar = Atom::var(GS.rescale_star);
        let tsrat_pow = tstar.pow(loop_3);
        let hfunction = Atom::var(GS.hfunction_lu_cut);
        tsrat_pow * hfunction * energy_conservation_delta_factor
    }

    fn single_th_prefactor_helper_atom(
        &self,
        order: u8,
        subspace_loop_count: usize,
        is_on_right: bool,
        include_integrated: bool,
    ) -> Atom {
        let loop_3 = subspace_loop_count as i64 * 3;

        let i = Atom::i();

        let radius = if is_on_right {
            Atom::var(GS.radius_right)
        } else {
            Atom::var(GS.radius_left)
        };

        let radius_star = if is_on_right {
            Atom::var(GS.radius_star_right)
        } else {
            Atom::var(GS.radius_star_left)
        };
        let uv_damp_plus = if is_on_right {
            Atom::var(GS.uv_damp_plus_right)
        } else {
            Atom::var(GS.uv_damp_plus_left)
        };
        let uv_damp_minus = if is_on_right {
            Atom::var(GS.uv_damp_minus_right)
        } else {
            Atom::var(GS.uv_damp_minus_left)
        };
        let hfunction = if is_on_right {
            Atom::var(GS.hfunction_right_th)
        } else {
            Atom::var(GS.hfunction_left_th)
        };

        let laurent_coeff_indices = (1..=order).map(|i| -(i as i8));

        let mut laurent_coeffs = laurent_coeff_indices.map(|laurent_coeff_index| {
            build_derivative_structure_atom(order, laurent_coeff_index)
                .replace(GS.rescale_star)
                .with(radius_star.clone())
        });

        let delta_r_plus = &radius - &radius_star;
        let delta_r_minus = -&radius - &radius_star;

        let jacobian_ratio = (Atom::one() / &radius).pow(loop_3 - 1);

        let local_prefactor =
            &jacobian_ratio * (uv_damp_plus / &delta_r_plus + uv_damp_minus / &delta_r_minus);

        let integrated_prefactor = if include_integrated {
            if is_on_right {
                -i * Atom::var(GS.pi) * &jacobian_ratio * hfunction
            } else {
                i * Atom::var(GS.pi) * &jacobian_ratio * hfunction
            }
        } else {
            Atom::zero()
        };

        let mut result = (local_prefactor + integrated_prefactor) * laurent_coeffs.next().unwrap();

        for pow in 2..=order {
            result += laurent_coeffs.next().unwrap()
                * &jacobian_ratio
                * (Atom::one() / delta_r_plus.pow(pow as i64)
                    + Atom::one() / delta_r_minus.pow(pow as i64));
        }

        debug!(
            "Threshold counterterm helper atom for order {} and loop number {}: {}",
            order, subspace_loop_count, result
        );

        result
    }

    fn iterated_th_prefactor_helper_atom(
        &self,
        left_order: u8,
        right_order: u8,
        left_subspace_loop_count: usize,
        right_subspace_loop_count: usize,
        include_integrated: bool,
    ) -> Atom {
        let left_prefactor = self.single_th_prefactor_helper_atom(
            left_order,
            left_subspace_loop_count,
            false,
            include_integrated,
        );
        let right_prefactor = self.single_th_prefactor_helper_atom(
            right_order,
            right_subspace_loop_count,
            true,
            include_integrated,
        );

        let mut product = left_prefactor * right_prefactor;

        product = product.replace_multiple(Self::fuse_left_right_replacement());
        product
    }

    fn single_th_prefactor_helper_params(
        &self,
        order: u8,
        _subspace_loop_count: usize,
        is_on_right: bool,
    ) -> Vec<Atom> {
        let radius_star = if is_on_right {
            Atom::var(GS.radius_star_right)
        } else {
            Atom::var(GS.radius_star_left)
        };

        let uv_damp_plus = if is_on_right {
            Atom::var(GS.uv_damp_plus_right)
        } else {
            Atom::var(GS.uv_damp_plus_left)
        };
        let uv_damp_minus = if is_on_right {
            Atom::var(GS.uv_damp_minus_right)
        } else {
            Atom::var(GS.uv_damp_minus_left)
        };
        let hfunction = if is_on_right {
            Atom::var(GS.hfunction_right_th)
        } else {
            Atom::var(GS.hfunction_left_th)
        };

        let radius = if is_on_right {
            Atom::var(GS.radius_right)
        } else {
            Atom::var(GS.radius_left)
        };

        let mut params = params_for_derivative_order(order)
            .into_iter()
            .map(|param| param.replace(GS.rescale_star).with(radius_star.clone()))
            .collect_vec();

        params.push(radius);
        params.push(radius_star);
        params.push(uv_damp_plus);
        params.push(uv_damp_minus);
        params.push(hfunction);
        params
    }

    fn iterated_th_prefactor_helper_params(&self, left_order: u8, right_order: u8) -> Vec<Atom> {
        let mut iterated_params = params_for_iterated_threshold_ct(left_order, right_order);
        let left_radius_star = Atom::var(GS.radius_star_left);
        let right_radius_star = Atom::var(GS.radius_star_right);
        let left_radius = Atom::var(GS.radius_left);
        let right_radius = Atom::var(GS.radius_right);
        let left_uv_damp_plus = Atom::var(GS.uv_damp_plus_left);
        let left_uv_damp_minus = Atom::var(GS.uv_damp_minus_left);
        let right_uv_damp_plus = Atom::var(GS.uv_damp_plus_right);
        let right_uv_damp_minus = Atom::var(GS.uv_damp_minus_right);
        let left_hfunction = Atom::var(GS.hfunction_left_th);
        let right_hfunction = Atom::var(GS.hfunction_right_th);

        iterated_params.push(left_radius);
        iterated_params.push(left_radius_star);
        iterated_params.push(left_uv_damp_plus);
        iterated_params.push(left_uv_damp_minus);
        iterated_params.push(left_hfunction);
        iterated_params.push(right_radius);
        iterated_params.push(right_radius_star);
        iterated_params.push(right_uv_damp_plus);
        iterated_params.push(right_uv_damp_minus);
        iterated_params.push(right_hfunction);
        iterated_params
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn single_th_helper(
        &self,
        order: u8,
        subspace_loop_count: usize,
        is_on_right: bool,
        include_integrated: bool,
        dual_shape: Option<Vec<Vec<usize>>>,
        optimization_settings: OptimizationSettings,
        evaluator_settings: &EvaluatorSettings,
    ) -> Result<GenericEvaluator> {
        let atom = self.single_th_prefactor_helper_atom(
            order,
            subspace_loop_count,
            is_on_right,
            include_integrated,
        );
        let params =
            self.single_th_prefactor_helper_params(order, subspace_loop_count, is_on_right);

        let mut fn_map = FunctionMap::new();
        fn_map
            .add_aliases([(
                GS.pi.into(),
                Atom::num(Rational::try_from(std::f64::consts::PI).unwrap()),
            )])
            .unwrap();

        let evaluator = GenericEvaluator::new_from_raw_params(
            [atom],
            &params,
            &fn_map,
            vec![],
            optimization_settings,
            dual_shape,
            evaluator_settings,
        )?
        .into_eager_only();

        Ok(evaluator)
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn iterated_th_helper(
        &self,
        left_order: u8,
        right_order: u8,
        left_subspace_loop_count: usize,
        right_subspace_loop_count: usize,
        include_integrated: bool,
        dual_shape: Option<Vec<Vec<usize>>>,
        optimization_settings: OptimizationSettings,
        evaluator_settings: &EvaluatorSettings,
    ) -> Result<GenericEvaluator> {
        let atom = self.iterated_th_prefactor_helper_atom(
            left_order,
            right_order,
            left_subspace_loop_count,
            right_subspace_loop_count,
            include_integrated,
        );

        let params = self.iterated_th_prefactor_helper_params(left_order, right_order);

        let mut fn_map = FunctionMap::new();
        fn_map
            .add_aliases([(
                GS.pi.into(),
                Atom::num(Rational::try_from(std::f64::consts::PI).unwrap()),
            )])
            .unwrap();

        let evaluator = GenericEvaluator::new_from_raw_params(
            [atom],
            &params,
            &fn_map,
            vec![],
            optimization_settings,
            dual_shape,
            evaluator_settings,
        )?
        .into_eager_only();

        Ok(evaluator)
    }

    fn fuse_left_right_replacement() -> Vec<Replacement> {
        let f = symbol!("f");

        vec![
            Replacement::new(
                (function!(f, GS.radius_star_left) * function!(f, GS.radius_star_right))
                    .to_pattern(),
                function!(f, GS.radius_star_left, GS.radius_star_right),
            ),
            Replacement::new(
                (function!(f, GS.radius_star_left)
                    * function!(Symbol::DERIVATIVE, W_.x_, f, GS.radius_star_right))
                .to_pattern(),
                function!(
                    Symbol::DERIVATIVE,
                    0,
                    W_.x_,
                    f,
                    GS.radius_star_left,
                    GS.radius_star_right
                ),
            ),
            Replacement::new(
                (function!(Symbol::DERIVATIVE, W_.x_, f, GS.radius_star_left)
                    * function!(f, GS.radius_star_right))
                .to_pattern(),
                function!(
                    Symbol::DERIVATIVE,
                    W_.x_,
                    0,
                    f,
                    GS.radius_star_left,
                    GS.radius_star_right
                ),
            ),
            Replacement::new(
                (function!(Symbol::DERIVATIVE, W_.x_, f, GS.radius_star_left)
                    * function!(Symbol::DERIVATIVE, W_.y_, f, GS.radius_star_right))
                .to_pattern(),
                function!(
                    Symbol::DERIVATIVE,
                    W_.x_,
                    W_.y_,
                    f,
                    GS.radius_star_left,
                    GS.radius_star_right
                ),
            ),
        ]
    }
    //fn th_prefactor_helper(
    //    &self,
    //    subspace_loop_count: usize,
    //    is_on_right: bool,
    //    include_integrated: bool,
    //) -> Atom {
    //    let radius = if is_on_right {
    //        Atom::var(GS.radius_right)
    //    } else {
    //        Atom::var(GS.radius_left)
    //    };

    //    let radius_star = if is_on_right {
    //        Atom::var(GS.radius_star_right)
    //    } else {
    //        Atom::var(GS.radius_star_left)
    //    };

    //    let grad_eta = if is_on_right {
    //        Atom::var(GS.deta_right_th)
    //    } else {
    //        Atom::var(GS.deta_left_th)
    //    };

    //    let uv_damp_plus = if is_on_right {
    //        Atom::var(GS.uv_damp_plus_right)
    //    } else {
    //        Atom::var(GS.uv_damp_plus_left)
    //    };

    //    let uv_damp_minus = if is_on_right {
    //        Atom::var(GS.uv_damp_minus_right)
    //    } else {
    //        Atom::var(GS.uv_damp_minus_left)
    //    };

    //    let h_function = if is_on_right {
    //        Atom::var(GS.hfunction_right_th)
    //    } else {
    //        Atom::var(GS.hfunction_left_th)
    //    };

    //    let i = if is_on_right { -Atom::i() } else { Atom::i() };
    //    let pi = Atom::var(GS.pi);

    //    let jacobian_ratio = (&radius_star / &radius).pow(subspace_loop_count as i64 * 3 - 1);

    //    let local_prefactor = (&uv_damp_plus / (&radius - &radius_star)
    //        + &uv_damp_minus / (-&radius - &radius_star))
    //        / &grad_eta
    //        * &jacobian_ratio;

    //    let integrated_prefactor = if include_integrated {
    //        &h_function * &i * &pi * &jacobian_ratio / &grad_eta
    //    } else {
    //        Atom::new()
    //    };

    //    debug!(
    //        "th prefactor local: {}, integrated: {}",
    //        local_prefactor, integrated_prefactor
    //    );

    //    local_prefactor + integrated_prefactor
    //}

    fn build_lmbs(&mut self) -> Result<()> {
        let mut lmbs: TiVec<LmbIndex, LoopMomentumBasis> = vec![].into();

        let externals: SuBitGraph = self.graph.empty_subgraph();
        let full_filter = self.graph.full_filter();
        let cut_graph = full_filter.subtract(&self.graph.initial_state_cut.right);

        for s in self.graph.all_spanning_forests_of(&cut_graph) {
            let mut lmb = self.graph.lmb_impl(&full_filter, &s, externals.clone())?;
            let mut exts = vec![];

            for i in lmb.loop_edges.iter() {
                let (_, p) = &self.graph[i];

                if self.graph.initial_state_cut.intersects(p) {
                    exts.push(*i);
                }
            }

            exts.sort();

            for e in exts {
                let mut loopid = None;
                for (l, s) in lmb.edge_signatures[e].internal.iter_enumerated() {
                    if s.is_non_zero() {
                        if loopid.is_none() {
                            loopid = Some(l);
                        } else {
                            panic!("external edge has multiple loop momenta")
                        }
                    }
                }
                lmb.put_loop_to_ext(loopid.unwrap());
            }
            let external_momentum_edge_order = self.graph.external_momentum_edge_order();
            lmb.canonicalize_external_order(&external_momentum_edge_order);
            //lmbs.push(self.graph.lmb_impl(&full_filter, &s, externals.clone()));
            lmbs.push(lmb);
        }

        let sorted_graph_loop_edges = self
            .graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .copied()
            .sorted()
            .collect_vec();

        let matching_lmb_index = lmbs.iter_enumerated().find_map(|(lmb_index, lmb)| {
            let sorted_loop_edges = lmb.loop_edges.iter().copied().sorted().collect_vec();
            (sorted_loop_edges == sorted_graph_loop_edges).then_some(lmb_index)
        });

        if let Some(matching_lmb_index) = matching_lmb_index {
            {
                let matching_lmb = &mut lmbs[matching_lmb_index];
                for (target_pos, target_edge) in
                    self.graph.loop_momentum_basis.loop_edges.iter_enumerated()
                {
                    let current_pos = matching_lmb
                        .loop_edges
                        .iter()
                        .find_position(|edge| *edge == target_edge)
                        .map(|(pos, _)| pos)
                        .unwrap();

                    if current_pos != target_pos.0 {
                        matching_lmb.swap_loops(current_pos.into(), target_pos);
                    }
                }
            }

            let matching_lmb_pos: usize = matching_lmb_index.into();
            if matching_lmb_pos != 0 {
                lmbs.raw[..=matching_lmb_pos].rotate_right(1);
            }
        } else {
            warn!(
                "Could not match current graph LMB against generated LMBs for graph {}",
                self.graph.name
            );
        }

        let sorted_graph_loop_edges = self
            .graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .copied()
            .sorted()
            .collect_vec();

        let matching_lmb_index = lmbs.iter_enumerated().find_map(|(lmb_index, lmb)| {
            let sorted_loop_edges = lmb.loop_edges.iter().copied().sorted().collect_vec();
            (sorted_loop_edges == sorted_graph_loop_edges).then_some(lmb_index)
        });

        if let Some(matching_lmb_index) = matching_lmb_index {
            {
                let matching_lmb = &mut lmbs[matching_lmb_index];
                for (target_pos, target_edge) in
                    self.graph.loop_momentum_basis.loop_edges.iter_enumerated()
                {
                    let current_pos = matching_lmb
                        .loop_edges
                        .iter()
                        .find_position(|edge| *edge == target_edge)
                        .map(|(pos, _)| pos)
                        .unwrap();

                    if current_pos != target_pos.0 {
                        matching_lmb.swap_loops(current_pos.into(), target_pos);
                    }
                }
            }

            let matching_lmb_pos: usize = matching_lmb_index.into();
            if matching_lmb_pos != 0 {
                lmbs.raw[..=matching_lmb_pos].rotate_right(1);
            }
        } else {
            warn!(
                "Could not match current graph LMB against generated LMBs for graph {}",
                self.graph.name
            );
        }

        self.derived_data.lmbs = Some(lmbs);
        Ok(())
    }

    fn build_multi_channeling_channels(&mut self, override_lmb_heuristics: bool) -> Result<()> {
        let lmbs = self.derived_data.lmbs.as_ref().unwrap();
        let channels = if override_lmb_heuristics {
            self.graph
                .build_multi_channeling_channels(lmbs, override_lmb_heuristics)
        } else {
            self.build_cross_section_multi_channeling_channels(lmbs)?
        };

        self.derived_data.multi_channeling_setup = Some(channels);
        Ok(())
    }

    fn build_cross_section_multi_channeling_channels(
        &self,
        lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
    ) -> Result<LmbMultiChannelingSetup> {
        let channels = self.select_cross_section_lmb_channel_indices(lmbs)?;
        debug!(
            "number of lmbs: {}, number of cross-section channels: {}",
            lmbs.len(),
            channels.len()
        );
        Ok(LmbMultiChannelingSetup {
            channels,
            graph: self.graph.clone(),
            all_bases: lmbs.clone(),
        })
    }

    fn select_cross_section_lmb_channel_indices(
        &self,
        lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
    ) -> Result<TiVec<crate::integrands::process::ChannelIndex, LmbIndex>> {
        if self.cuts.is_empty() {
            return Ok(self.graph.select_amplitude_lmb_channel_indices(
                lmbs,
                false,
                LmbChannelFallback::CurrentGraphBasis,
            ));
        }

        let mut lmb_index_by_loop_edges = BTreeMap::<Vec<EdgeIndex>, LmbIndex>::new();
        for (lmb_index, lmb) in lmbs.iter_enumerated() {
            lmb_index_by_loop_edges
                .entry(lmb.loop_edges.iter().copied().sorted().collect())
                .or_insert(lmb_index);
        }

        let mut selected_loop_edge_sets = BTreeSet::<Vec<EdgeIndex>>::new();
        for cut in &self.cuts {
            let cut_edges = cut.cut.as_subgraph();
            let cut_edge_ids = self
                .graph
                .underlying
                .iter_edges_of(&cut_edges)
                .map(|(_, edge_id, _)| edge_id)
                .sorted()
                .collect_vec();
            let Some(excluded_cut_edge) = self.excluded_cut_edge_for_lmb_channel(&cut_edge_ids)
            else {
                continue;
            };

            let fixed_cut_edges = cut_edge_ids
                .iter()
                .copied()
                .filter(|edge_id| *edge_id != excluded_cut_edge)
                .collect::<BTreeSet<_>>();
            let left_basis_edge_sets =
                self.selected_cut_side_lmb_edge_sets(cut.left.subtract(&cut_edges))?;
            let right_basis_edge_sets =
                self.selected_cut_side_lmb_edge_sets(cut.right.subtract(&cut_edges))?;

            for left_basis_edges in &left_basis_edge_sets {
                for right_basis_edges in &right_basis_edge_sets {
                    let loop_edges = fixed_cut_edges
                        .iter()
                        .copied()
                        .chain(left_basis_edges.iter().copied())
                        .chain(right_basis_edges.iter().copied())
                        .collect::<BTreeSet<_>>()
                        .into_iter()
                        .collect_vec();
                    selected_loop_edge_sets.insert(loop_edges);
                }
            }
        }

        let mut channels = Vec::<LmbIndex>::new();
        for loop_edges in selected_loop_edge_sets {
            let lmb_index = lmb_index_by_loop_edges.get(&loop_edges).ok_or_else(|| {
                eyre!(
                    "Could not find a generated cross-section LMB with loop edges [{}] for graph '{}'.",
                    loop_edges.iter().map(|edge| edge.to_string()).join(", "),
                    self.graph.name
                )
            })?;
            channels.push(*lmb_index);
        }

        if channels.is_empty() {
            Ok(self.graph.select_amplitude_lmb_channel_indices(
                lmbs,
                false,
                LmbChannelFallback::CurrentGraphBasis,
            ))
        } else {
            Ok(channels.into_iter().sorted().dedup().collect())
        }
    }

    fn excluded_cut_edge_for_lmb_channel(&self, cut_edge_ids: &[EdgeIndex]) -> Option<EdgeIndex> {
        Self::excluded_cut_edge_for_lmb_channel_in(&self.graph, cut_edge_ids)
    }

    fn excluded_cut_edge_for_lmb_channel_in(
        graph: &Graph,
        cut_edge_ids: &[EdgeIndex],
    ) -> Option<EdgeIndex> {
        cut_edge_ids
            .iter()
            .copied()
            .filter(|edge_id| graph[*edge_id].particle.is_massive())
            .min()
            .or_else(|| {
                cut_edge_ids
                    .iter()
                    .copied()
                    .filter(|edge_id| graph[*edge_id].particle.is_fermion())
                    .min()
            })
            .or_else(|| cut_edge_ids.iter().copied().min())
    }

    fn selected_cut_side_lmb_edge_sets(
        &self,
        mut side_subgraph: SuBitGraph,
    ) -> Result<Vec<Vec<EdgeIndex>>> {
        for (pair, _, edge) in self.graph.underlying.iter_edges() {
            if edge.data.is_dummy {
                side_subgraph.sub(pair);
            }
        }

        if self.graph.underlying.cyclotomatic_number(&side_subgraph) == 0 {
            return Ok(vec![Vec::new()]);
        }

        let side_lmbs = self.graph.generate_loop_momentum_bases_of(&side_subgraph);
        if side_lmbs.is_empty() {
            return Err(eyre!(
                "Could not generate cut-side LMBs for a non-tree side of graph '{}'.",
                self.graph.name
            ));
        }

        let selected_side_lmbs = self.graph.select_amplitude_lmb_channel_indices(
            &side_lmbs,
            false,
            LmbChannelFallback::FirstBasis,
        );
        let mut edge_sets = BTreeSet::<Vec<EdgeIndex>>::new();
        for lmb_index in selected_side_lmbs {
            edge_sets.insert(
                side_lmbs[lmb_index]
                    .loop_edges
                    .iter()
                    .copied()
                    .sorted()
                    .collect(),
            );
        }

        if edge_sets.is_empty() {
            Ok(vec![Vec::new()])
        } else {
            Ok(edge_sets.into_iter().collect())
        }
    }

    fn threshold_counterterm_association(
        &self,
        threshold_candidate: &TopologicalThresholdCandidate,
        sandwich: &SuBitGraph,
        cut: &OrientedCut,
    ) -> ThresholdCountertermAssociation {
        let boundary_edges = |boundary: &OrientedCut| {
            let boundary_filter = boundary.left.union(&boundary.right).intersection(sandwich);
            let edges = self
                .graph
                .iter_edges_of(&boundary_filter)
                .map(|(_, edge_id, _)| edge_id)
                .sorted()
                .collect();
            let has_left_orientation =
                boundary_filter.intersection(&boundary.left).n_included() > 0;
            let has_right_orientation =
                boundary_filter.intersection(&boundary.right).n_included() > 0;
            (edges, has_left_orientation != has_right_orientation)
        };

        let (cut_boundary_edges, cut_bound_is_applicable) = boundary_edges(cut);
        let (threshold_boundary_edges, threshold_bound_is_applicable) =
            boundary_edges(&threshold_candidate.cut);

        ThresholdCountertermAssociation {
            esurface_id: threshold_candidate.esurface_id,
            cut_boundary_edges,
            threshold_boundary_edges,
            invariant_bound_is_applicable: cut_bound_is_applicable && threshold_bound_is_applicable,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn classify_threshold_counterterm_association(
        &self,
        association: &ThresholdCountertermAssociation,
        cut_id: CutId,
        cut: &CrossSectionCut,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        model: &Model,
        runtime_settings: &RuntimeSettings,
        settings: &GenerationSettings,
    ) -> ThresholdCountertermStatus {
        let esurface = &self.graph.surface_cache.esurface_cache[association.esurface_id];
        if !esurface.has_radial_dependence_in_subspace(subspace, all_lmbs, &self.graph) {
            debug!(
                "Skipping graph '{}' cut {} threshold E-surface {} after cut-relative classification {:?}",
                self.graph.name,
                cut_id.0,
                association.esurface_id.0,
                ThresholdCountertermStatus::NoRadialDependence,
            );
            return ThresholdCountertermStatus::NoRadialDependence;
        }

        let structural_status = if association.invariant_bound_is_applicable {
            self.graph.classify_threshold_pinch(
                &association.cut_boundary_edges,
                &association.threshold_boundary_edges,
            )
        } else {
            ThresholdPinchStatus::NotProven
        };
        if structural_status == ThresholdPinchStatus::Always {
            debug!(
                "Skipping graph '{}' cut {} threshold E-surface {} after cut-relative classification {:?}",
                self.graph.name,
                cut_id.0,
                association.esurface_id.0,
                ThresholdCountertermStatus::AlwaysPinched,
            );
            return ThresholdCountertermStatus::AlwaysPinched;
        }

        let status = if settings.threshold_subtraction.check_esurface_at_generation {
            association.classify_for_model(
                &self.graph,
                cut,
                model,
                &self.graph.param_builder,
                runtime_settings,
                settings.threshold_subtraction.esurface_existence_threshold,
            )
        } else {
            match structural_status {
                ThresholdPinchStatus::CanBecome => ThresholdCountertermStatus::CanBecomePinched,
                ThresholdPinchStatus::NotProven => ThresholdCountertermStatus::PotentiallyExisting,
                ThresholdPinchStatus::Always => unreachable!(
                    "identically pinched threshold associations return before classification"
                ),
            }
        };
        if !status
            .is_eligible_for_generation(settings.threshold_subtraction.check_esurface_at_generation)
        {
            debug!(
                "Skipping graph '{}' cut {} threshold E-surface {} after cut-relative classification {:?}",
                self.graph.name, cut_id.0, association.esurface_id.0, status,
            );
        }
        status
    }

    fn topological_threshold_candidates(&self) -> Result<Vec<TopologicalThresholdCandidate>> {
        let mut candidates = self.graph.all_st_cuts_for_cs(
            self.source_nodes.clone(),
            self.target_nodes.clone(),
            &self.graph.get_initial_state_tree().0,
        );
        candidates.retain(|(_left, cut, _right)| cut.nedges(&self.graph) > 1);
        candidates.sort_by(|a, b| a.1.cmp(&b.1));

        candidates
            .into_iter()
            .map(|(left, cut, right)| {
                let threshold_esurface = Esurface::new_from_cut_left(
                    &self.graph.underlying,
                    &CrossSectionCut {
                        cut: cut.clone(),
                        left: left.clone(),
                        right: right.clone(),
                    },
                    Some(&self.graph.initial_state_cut),
                );
                let esurface_id = self
                    .graph
                    .surface_cache
                    .esurface_cache
                    .position(|esurface| esurface == &threshold_esurface)
                    .ok_or_else(|| {
                        eyre!(
                            "Topology-discovered threshold surface {:?} is missing from graph '{}' CFF surface cache",
                            threshold_esurface.energies,
                            self.graph.name,
                        )
                    })?;
                Ok(TopologicalThresholdCandidate {
                    left,
                    cut,
                    right,
                    esurface_id,
                })
            })
            .collect()
    }

    fn build_threshold_counterterm(
        &mut self,
        model: &Model,
        settings: &GenerationSettings,
        runtime_settings: &RuntimeSettings,
        all_possible_thresholds: &[TopologicalThresholdCandidate],
        vakint: &Vakint,
    ) -> Result<()> {
        // Thresholds are topology-discovered independently of whether a CT association survives
        // the structural and optional current-model checks below.
        let mut cut_threshold_associations: TiVec<CutId, CutThresholdCountertermAssociations> =
            ti_vec![CutThresholdCountertermAssociations::default(); self.cuts.len()];

        let mut cut_group_by_cut = vec![None; self.cuts.len()];
        for (cut_group_id, cut_group) in self
            .derived_data
            .cut_group_data
            .cut_groups
            .iter_enumerated()
        {
            for cut_id in &cut_group.cuts {
                cut_group_by_cut[cut_id.0] = Some(cut_group_id);
            }
        }
        let all_lmbs = self
            .derived_data
            .lmbs
            .as_ref()
            .expect("threshold generation requires loop-momentum bases");

        // Keep topology-discovered E-surfaces in the graph inventory, but do not generate a
        // threshold CT association for any E-surface that is also one of this graph's physical
        // Cutkosky cuts when the corresponding generation setting is enabled. This comparison
        // must cover every physical cut, not only the cut currently being processed below.
        let subtraction_threshold_candidates = all_possible_thresholds
            .iter()
            .filter(|threshold_candidate| {
                !settings.threshold_subtraction.skip_thresholds_that_are_cuts
                    || !self
                        .cuts
                        .iter()
                        .any(|physical_cut| physical_cut.cut == threshold_candidate.cut)
            })
            .collect_vec();

        for (cut_id, cut) in self.cuts.iter_enumerated() {
            let cut_group_id = cut_group_by_cut[cut_id.0].ok_or_else(|| {
                eyre!(
                    "Physical cut {} in graph '{}' is missing from the cut-group partition",
                    cut_id.0,
                    self.graph.name,
                )
            })?;
            let (left_subspace, right_subspace) = &self.derived_data.subspace_data[cut_group_id];

            for threshold_candidate in &subtraction_threshold_candidates {
                // if the subgraph on the left of the threshold cut is a subgraph of the left amplitude, then the threshold is on the left of the cut
                if cut.left.includes(&threshold_candidate.left) {
                    let sandwich = cut.left.intersection(&threshold_candidate.right);

                    // The shared projected-signature check below determines whether the threshold
                    // cuts an active loop. Connectivity remains a separate structural requirement.
                    if self.graph.underlying.is_connected(&sandwich) {
                        let association = self.threshold_counterterm_association(
                            threshold_candidate,
                            &sandwich,
                            &cut.cut,
                        );
                        let status = self.classify_threshold_counterterm_association(
                            &association,
                            cut_id,
                            cut,
                            left_subspace,
                            all_lmbs,
                            model,
                            runtime_settings,
                            settings,
                        );
                        if status.is_eligible_for_generation(
                            settings.threshold_subtraction.check_esurface_at_generation,
                        ) {
                            cut_threshold_associations[cut_id].left.push(association);
                        }
                    }
                } else if cut.right.includes(&threshold_candidate.right) {
                    let sandwich = cut.right.intersection(&threshold_candidate.left);
                    if self.graph.underlying.is_connected(&sandwich) {
                        // threshold_esurface.subspace_graph =
                        //     InternalSubGraph::cleaned_filter_pessimist(
                        //         cut.right.clone(),
                        //         &self.graph,
                        //     );
                        let association = self.threshold_counterterm_association(
                            threshold_candidate,
                            &sandwich,
                            &cut.cut,
                        );
                        let status = self.classify_threshold_counterterm_association(
                            &association,
                            cut_id,
                            cut,
                            right_subspace,
                            all_lmbs,
                            model,
                            runtime_settings,
                            settings,
                        );
                        if status.is_eligible_for_generation(
                            settings.threshold_subtraction.check_esurface_at_generation,
                        ) {
                            cut_threshold_associations[cut_id].right.push(association);
                        }
                    }
                }
            }
        }

        let left_cut_threshold_data: TiVec<CutId, Vec<EsurfaceID>> = cut_threshold_associations
            .iter()
            .map(|associations| {
                associations
                    .left
                    .iter()
                    .map(|association| association.esurface_id)
                    .collect()
            })
            .collect();
        let right_cut_threshold_data: TiVec<CutId, Vec<EsurfaceID>> = cut_threshold_associations
            .iter()
            .map(|associations| {
                associations
                    .right
                    .iter()
                    .map(|association| association.esurface_id)
                    .collect()
            })
            .collect();
        self.derived_data.cut_threshold_associations = cut_threshold_associations;

        let threshold_raised_data = self.graph.determine_raised_esurfaces_from_expression(
            self.derived_data
                .global_cff_expression
                .as_ref()
                .expect("global_cff_expression should have been created"),
        );
        let mut raised_threshold_ids: TiVec<EsurfaceID, Option<RaisedEsurfaceId>> =
            ti_vec![None; self.graph.surface_cache.esurface_cache.len()];

        for (raised_threshold_id, raised_group) in
            threshold_raised_data.raised_groups.iter_enumerated()
        {
            for &esurface_id in &raised_group.esurface_ids {
                raised_threshold_ids[esurface_id] = Some(raised_threshold_id);
            }
        }

        let raised_threshold_ids: TiVec<EsurfaceID, RaisedEsurfaceId> = raised_threshold_ids
            .into_iter()
            .map(|raised_threshold_id| {
                raised_threshold_id
                    .expect("every esurface should belong to exactly one raised threshold group")
            })
            .collect();

        let collect_raised_threshold_groups = |threshold_ids: Vec<EsurfaceID>| {
            let mut groups = Vec::new();
            for esurface_id in threshold_ids.into_iter().sorted().dedup() {
                let raised_threshold_id = raised_threshold_ids[esurface_id];
                let raised_group = threshold_raised_data.raised_groups[raised_threshold_id].clone();
                if !groups.contains(&raised_group) {
                    groups.push(raised_group);
                }
            }
            groups
        };

        let mut left_cut_group_threshold_data: TiVec<
            CutGroupId,
            TiVec<LeftThresholdId, RaisedEsurfaceGroup>,
        > = TiVec::new();

        let mut right_cut_group_threshold_data: TiVec<
            CutGroupId,
            TiVec<RightThresholdId, RaisedEsurfaceGroup>,
        > = TiVec::new();

        for cut_group in self.derived_data.cut_group_data.cut_groups.iter() {
            let left_thresholds = collect_raised_threshold_groups(
                cut_group
                    .cuts
                    .iter()
                    .flat_map(|cut_id| left_cut_threshold_data[*cut_id].iter().copied())
                    .collect(),
            );

            let mut right_thresholds = collect_raised_threshold_groups(
                cut_group
                    .cuts
                    .iter()
                    .flat_map(|cut_id| right_cut_threshold_data[*cut_id].iter().copied())
                    .collect(),
            );

            right_thresholds.retain(|raised_group| !left_thresholds.contains(raised_group));

            left_cut_group_threshold_data.push(left_thresholds.into());
            right_cut_group_threshold_data.push(right_thresholds.into());
        }

        let mut cut_structure = vec![];

        for (cut_group_id, cut_group) in self
            .derived_data
            .cut_group_data
            .cut_groups
            .iter_enumerated()
        {
            let left_thresholds = &left_cut_group_threshold_data[cut_group_id];
            let right_thresholds = &right_cut_group_threshold_data[cut_group_id];

            let cutkosky_cut_union = cut_group
                .cuts
                .iter()
                .map(|cut_id| self.cuts[*cut_id].cut.as_subgraph())
                .reduce(|a, b| a.union(&b))
                .unwrap_or(self.graph.empty_subgraph());

            let add_threshold_group_to_union =
                |base: SuBitGraph, raised_group: &RaisedEsurfaceGroup| {
                    let representative_esurface =
                        &self.graph.surface_cache.esurface_cache[raised_group.esurface_ids[0]];

                    representative_esurface
                        .energies
                        .iter()
                        .map(|edge_id| self.graph.get_edge_subgraph(*edge_id))
                        .fold(base, |acc, subgraph| acc.union(&subgraph))
                };

            for raised_esurface_group in left_thresholds {
                let esurface_cut_union =
                    add_threshold_group_to_union(cutkosky_cut_union.clone(), raised_esurface_group);

                cut_structure.push(CutSet {
                    residue_selector: ResidueSelector {
                        lu_cut: Some(cut_group.related_esurface_group.clone()),
                        left_th_cut: Some(raised_esurface_group.clone()),
                        right_th_cut: None,
                    },
                    union: esurface_cut_union,
                    canonicalize_external_shifts: false,
                });
            }

            for raised_esurface_group in right_thresholds {
                let esurface_cut_union =
                    add_threshold_group_to_union(cutkosky_cut_union.clone(), raised_esurface_group);

                cut_structure.push(CutSet {
                    residue_selector: ResidueSelector {
                        lu_cut: Some(cut_group.related_esurface_group.clone()),
                        left_th_cut: None,
                        right_th_cut: Some(raised_esurface_group.clone()),
                    },
                    union: esurface_cut_union,
                    canonicalize_external_shifts: false,
                });
            }

            for (left_raised_esurface_group, right_raised_esurface_group) in left_thresholds
                .iter()
                .cartesian_product(right_thresholds.iter())
            {
                let esurface_cut_union = add_threshold_group_to_union(
                    add_threshold_group_to_union(
                        cutkosky_cut_union.clone(),
                        left_raised_esurface_group,
                    ),
                    right_raised_esurface_group,
                );

                cut_structure.push(CutSet {
                    residue_selector: ResidueSelector {
                        lu_cut: Some(cut_group.related_esurface_group.clone()),
                        left_th_cut: Some(left_raised_esurface_group.clone()),
                        right_th_cut: Some(right_raised_esurface_group.clone()),
                    },
                    union: esurface_cut_union,
                    canonicalize_external_shifts: false,
                });
            }
        }

        let cut_structure = CutStructure {
            cuts: cut_structure,
        };

        let cut_woods = CutWoods::new(cut_structure, &self.graph, &settings.uv);
        let mut cut_forests = cut_woods.unfold(&self.graph);
        let valid_orientations: Vec<_> = self
            .derived_data
            .global_cff_expression
            .as_ref()
            .expect("global_cff_expression should have been created")
            .orientations
            .iter()
            .map(|orientation| orientation.data.orientation.clone())
            .collect();

        cut_forests.compute(
            &mut self.graph,
            vakint,
            &valid_orientations,
            &settings.uv,
            &settings.orientation_pattern,
        )?;

        let mut threshold_counterterms = cut_forests
            .orientation_parametric_exprs(&self.graph, &settings.uv)?
            .into_iter();

        let lu_prefactor = self.lu_prefactor_helper();

        let mut result = TiVec::<CutGroupId, LUCounterTermData>::new();
        for (cut_group_id, _cut_group) in self
            .derived_data
            .cut_group_data
            .cut_groups
            .iter_enumerated()
        {
            let (left_subspace, right_subspace) = &self.derived_data.subspace_data[cut_group_id];

            let left_rstar_pow =
                Atom::var(GS.radius_star_left).pow(left_subspace.loopcount() as i32 * 3 - 1);

            let right_rstar_pow =
                Atom::var(GS.radius_star_right).pow(right_subspace.loopcount() as i32 * 3 - 1);

            let mut left_atoms = TiVec::<LeftThresholdId, _>::new();
            let mut right_atoms = TiVec::<RightThresholdId, _>::new();
            let mut iterated_atoms = vec![];

            for _ in 0..left_cut_group_threshold_data[cut_group_id].len() {
                left_atoms.push(
                    threshold_counterterms
                        .next()
                        .unwrap()
                        .map(|x| x * &lu_prefactor * &left_rstar_pow),
                );
            }

            for _ in 0..right_cut_group_threshold_data[cut_group_id].len() {
                right_atoms.push(
                    threshold_counterterms
                        .next()
                        .unwrap()
                        .map(|x| x * &lu_prefactor * &right_rstar_pow),
                );
            }

            for _ in 0..(left_cut_group_threshold_data[cut_group_id].len()
                * right_cut_group_threshold_data[cut_group_id].len())
            {
                iterated_atoms.push(
                    threshold_counterterms
                        .next()
                        .unwrap()
                        .map(|x| x * &lu_prefactor * &left_rstar_pow * &right_rstar_pow),
                );
            }

            let iterated_collection =
                IteratedCtCollection::new(iterated_atoms, left_atoms.len(), right_atoms.len());

            let counterterm_data = LUCounterTermData {
                left_thresholds: left_cut_group_threshold_data[cut_group_id].clone(),
                right_thresholds: right_cut_group_threshold_data[cut_group_id].clone(),
                left_atoms,
                right_atoms,
                iterated: iterated_collection,
            };
            result.push(counterterm_data);
        }

        let max_dual_size =
            max_dual_size_for_cut_cff_indices(result.iter().flat_map(|counterterm_data| {
                counterterm_data
                    .left_atoms
                    .iter()
                    .chain(counterterm_data.right_atoms.iter())
                    .chain(counterterm_data.iterated.iter())
                    .flat_map(|integrands| integrands.integrands.keys())
            }));
        self.graph.param_builder.initialize_duals(max_dual_size);

        self.derived_data.threshold_counterterms = result;
        Ok(())
    }

    fn build_subspace_data(&mut self) -> Result<()> {
        let all_lmbs = self.derived_data.lmbs.as_ref().unwrap();

        let subspace_data = self
            .derived_data
            .cut_group_data
            .cut_groups
            .iter_enumerated()
            .map(|(cut_group_id, cut_group)| {
                let representative_cut_id = cut_group.cuts.first().copied().ok_or_else(|| {
                    eyre!(
                        "Graph '{}' has an empty cut group {} while building threshold-counterterm subspaces",
                        self.graph.name,
                        cut_group_id.0,
                    )
                })?;
                let valid_subspace_lmbs = all_lmbs
                    .iter_enumerated()
                    .filter_map(|(index, lmb)| {
                        let mut edges_in_cut = self
                            .graph
                            .underlying
                            .iter_edges_of(&self.cuts[representative_cut_id].cut)
                            .map(|(_, e, _)| e)
                            .collect_vec();

                        edges_in_cut.retain(|e| !lmb.loop_edges.contains(e));
                        if edges_in_cut.len() == 1 {
                            Some(index)
                        } else {
                            None
                        }
                    })
                    .collect_vec();

                let left_subgraphs = cut_group
                    .cuts
                    .iter()
                    .map(|cut_id| self.cuts[*cut_id].left.clone())
                    .sorted_by(|subgraph_a, subgraph_b| {
                        let n_edges_a = subgraph_a.nedges(&self.graph);
                        let n_edges_b = subgraph_b.nedges(&self.graph);
                        n_edges_a.cmp(&n_edges_b)
                    })
                    .collect_vec();

                let right_subgraphs = cut_group
                    .cuts
                    .iter()
                    .map(|cut_id| self.cuts[*cut_id].right.clone())
                    .sorted_by(|subgraph_a, subgraph_b| {
                        let n_edges_a = subgraph_a.nedges(&self.graph);
                        let n_edges_b = subgraph_b.nedges(&self.graph);
                        n_edges_a.cmp(&n_edges_b)
                    })
                    .collect_vec();

                let smallest_left_subgraph = left_subgraphs.first().cloned().ok_or_else(|| {
                    eyre!(
                        "Graph '{}' cut group {} has no left subgraph",
                        self.graph.name,
                        cut_group_id.0,
                    )
                })?;
                let smallest_right_subgraph = right_subgraphs.first().cloned().ok_or_else(|| {
                    eyre!(
                        "Graph '{}' cut group {} has no right subgraph",
                        self.graph.name,
                        cut_group_id.0,
                    )
                })?;

                let mut possible_subspaces = Vec::new();
                let mut rejected_lmbs = Vec::new();
                for lmb_index in valid_subspace_lmbs {
                    let left = SubspaceData::new_with_user_selected_lmb(
                        smallest_left_subgraph.clone(),
                        lmb_index,
                        &self.graph,
                        all_lmbs,
                    );
                    let right = SubspaceData::new_with_user_selected_lmb(
                        smallest_right_subgraph.clone(),
                        lmb_index,
                        &self.graph,
                        all_lmbs,
                    );

                    match (left, right) {
                        (Ok(left), Ok(right)) if left.is_mergable_with(&right) => {
                            possible_subspaces.push((left, right));
                        }
                        (Ok(left), Ok(right)) => rejected_lmbs.push(format!(
                            "LMB {} produced non-disjoint subspaces: left={:?}, right={:?}",
                            usize::from(lmb_index),
                            left.iter_lmb_indices().collect_vec(),
                            right.iter_lmb_indices().collect_vec(),
                        )),
                        (left, right) => rejected_lmbs.push(format!(
                            "LMB {}: left={}, right={}",
                            usize::from(lmb_index),
                            left.err()
                                .map(|error| format!("{error:#}"))
                                .unwrap_or_else(|| "compatible".to_string()),
                            right
                                .err()
                                .map(|error| format!("{error:#}"))
                                .unwrap_or_else(|| "compatible".to_string()),
                        )),
                    }
                }

                possible_subspaces.sort_by_key(|(left, right)| {
                    (
                        left.iter_basis_edges(all_lmbs).collect_vec(),
                        right.iter_basis_edges(all_lmbs).collect_vec(),
                    )
                });

                possible_subspaces.first().cloned().ok_or_else(|| {
                    eyre!(
                        "No topology-compatible parent LMB found for graph '{}' cut group {}. Rejections:\n{}",
                        self.graph.name,
                        cut_group_id.0,
                        rejected_lmbs.join("\n"),
                    )
                })
            })
            .collect::<Result<_>>()?;

        let _: () = self.derived_data.subspace_data = subspace_data;
        Ok(())
    }

    fn generate_term_for_graph(
        &self,
        _model: &Model,
        settings: &GlobalSettings,
    ) -> Result<(CrossSectionGraphTerm, GraphGenerationStats)> {
        CrossSectionGraphTerm::from_cross_section_graph(self, settings)
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSectionDerivedData {
    pub orientations: Option<TiVec<OrientationID, EdgeVec<Orientation>>>,
    pub cut_paramatric_integrand: TiVec<CutGroupId, ParametricIntegrands>,
    pub global_cff_expression: Option<CFFExpression<OrientationID>>,
    pub lmbs: Option<TiVec<LmbIndex, LoopMomentumBasis>>,
    pub multi_channeling_setup: Option<LmbMultiChannelingSetup>,
    pub threshold_counterterms: TiVec<CutGroupId, LUCounterTermData>,
    /// Graph-level inventory of every topology-discovered threshold E-surface. This remains
    /// independent of whether a counterterm is generated for any particular physical cut.
    pub threshold_candidate_esurface_ids: Vec<EsurfaceID>,
    /// Exact generated left/right threshold associations for each physical cut. Runtime
    /// evaluators aggregate these into cut groups, but display and model reclassification
    /// must retain the physical-cut-relative eligibility information.
    pub cut_threshold_associations: TiVec<CutId, CutThresholdCountertermAssociations>,
    pub subspace_data: TiVec<CutGroupId, (SubspaceData, SubspaceData)>,
    pub cut_group_data: CutGroupData,
}

/// Partition and residue-derivative metadata for all Cutkosky cuts, including order-one cuts.
#[derive(Clone, Encode, Decode, Debug)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CutGroupData {
    pub cut_groups: TiVec<CutGroupId, CutGroup>,
    pub dual_shapes: Vec<Vec<Vec<usize>>>,
    pub pass_two_evaluators: Vec<GenericEvaluator>,
}

/// Physical cuts that share one normalized E-surface residue, whether ordinary or raised.
#[derive(Clone, Encode, Decode, Debug)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CutGroup {
    pub cuts: Vec<CutId>,
    pub related_esurface_group: RaisedEsurfaceGroup,
}

impl Default for CutGroupData {
    fn default() -> Self {
        Self::new()
    }
}

impl CutGroupData {
    pub fn new() -> Self {
        CutGroupData {
            cut_groups: TiVec::new(),
            dual_shapes: vec![],
            pass_two_evaluators: vec![],
        }
    }

    pub fn new_from_esurface(
        raised_esurface_data: &RaisedEsurfaceData,
        cut_esurface_map: &TiVec<CutId, EsurfaceID>,
        evaluator_settings: &EvaluatorSettings,
    ) -> (Self, GraphGenerationStats) {
        let mut stats = GraphGenerationStats::default();
        let reversed_map = cut_esurface_map
            .iter_enumerated()
            .map(|(cut_id, &esurface_id)| (esurface_id, cut_id))
            .collect::<HashMap<EsurfaceID, CutId>>();

        let mut groups = TiVec::new();

        for (_raised_esurface_id, raised_esurface_group) in
            raised_esurface_data.raised_groups.iter_enumerated()
        {
            if cut_esurface_map.contains(&raised_esurface_group.esurface_ids[0]) {
                let cuts = raised_esurface_group
                    .esurface_ids
                    .iter()
                    .map(|esurface_id| reversed_map[esurface_id])
                    .collect::<Vec<_>>();

                let cut_group = CutGroup {
                    cuts,
                    related_esurface_group: raised_esurface_group.clone(),
                };

                groups.push(cut_group);
            } else {
                continue;
            }
        }

        let global_max_occurence = groups
            .iter()
            .map(|group| group.related_esurface_group.max_occurence)
            .max()
            .unwrap_or_else(|| {
                println!("corrupted groups");
                panic!();
            });

        let dual_shapes = (1..global_max_occurence)
            .map(simple_n_deriv_shape)
            .collect();

        let pass_two_evaluators = (1..=global_max_occurence)
            .map(|i| {
                let evaluator_started = std::time::Instant::now();
                let evaluator = build_derivative_structure(i as u8, -1, evaluator_settings);
                stats.evaluator_symbolica_time += evaluator_started.elapsed();
                stats.evaluator_count += 1;
                evaluator
            })
            .collect();

        (
            Self {
                cut_groups: groups,
                dual_shapes,
                pass_two_evaluators,
            },
            stats,
        )
    }
}

impl CrossSectionDerivedData {
    fn new_empty() -> Self {
        Self {
            orientations: None,
            global_cff_expression: None,
            cut_paramatric_integrand: TiVec::new(),
            lmbs: None,
            multi_channeling_setup: None,
            threshold_counterterms: TiVec::new(),
            threshold_candidate_esurface_ids: Vec::new(),
            cut_threshold_associations: TiVec::new(),
            subspace_data: TiVec::new(),
            cut_group_data: CutGroupData::new(),
        }
    }
}

pub(crate) fn build_derivative_structure_atom(
    singularity_order: u8,
    laurent_coefficient: i8,
) -> Atom {
    assert!(
        laurent_coefficient <= -1,
        "only laurent coefficients up to -1 are supported"
    );

    assert!(
        singularity_order >= 1,
        "eta order must be at least 1, got {singularity_order}"
    );
    assert!(
        singularity_order >= -laurent_coefficient as u8,
        "eta order must be at least the negative of the laurent coefficient, got {singularity_order} for laurent coefficient {laurent_coefficient}"
    );

    let order = singularity_order as i32;
    let laurent_coefficient = laurent_coefficient as i32;
    let f = symbol!("f");

    let expansion = parse!("η(t)")
        .series(GS.rescale, Atom::var(GS.rescale_star), (order, 1))
        .unwrap()
        .to_atom()
        .replace(function!(symbol!("η"), GS.rescale_star))
        .level_range((0, Some(0)))
        .with(0);

    let mut expression_to_derive = function!(f, GS.rescale)
        * expansion.pow(-order)
        * (GS.rescale - GS.rescale_star).pow(order);

    for _ in 1..=(order + laurent_coefficient) {
        expression_to_derive = expression_to_derive.derivative(GS.rescale);
    }

    expression_to_derive = expression_to_derive
        .replace(GS.rescale - GS.rescale_star)
        .with(parse!("delta_t"));

    let polynomial_in_delta_t = expression_to_derive
        .series(symbol!("delta_t"), Atom::num(0), (0, 1))
        .unwrap();

    let factorial_prefactor = (2..=(order + laurent_coefficient)).product::<i32>();
    debug!("factorial prefactor: {}", factorial_prefactor);
    let mut expression_to_derive = polynomial_in_delta_t.to_atom() / Atom::num(factorial_prefactor);

    expression_to_derive = expression_to_derive
        .replace(GS.rescale)
        .with(GS.rescale_star);

    expression_to_derive
}

pub(crate) fn build_derivative_structure(
    singularity_order: u8,
    laurent_coefficient: i8,
    evaluator_settings: &EvaluatorSettings,
) -> GenericEvaluator {
    let expression_to_derive =
        build_derivative_structure_atom(singularity_order, laurent_coefficient);

    let params = params_for_derivative_order(singularity_order);

    GenericEvaluator::new_from_raw_params(
        [expression_to_derive],
        &params,
        &FunctionMap::default(),
        vec![],
        evaluator_settings.optimization_settings(),
        None,
        evaluator_settings,
    )
    .unwrap()
    .into_eager_only()
}

fn ordered_f_derivative_params(
    base: Atom,
    derivative_vars: &[Symbol],
    derivative_shape: Option<&Vec<Vec<usize>>>,
) -> Vec<Atom> {
    derivative_shape
        .cloned()
        .unwrap_or_else(|| vec![vec![0; derivative_vars.len()]])
        .into_iter()
        .map(|orders| {
            debug_assert_eq!(orders.len(), derivative_vars.len());

            let mut param = base.clone();
            for (derivative_var, order) in derivative_vars.iter().zip(orders) {
                for _ in 0..order {
                    param = param.derivative(*derivative_var);
                }
            }

            param
        })
        .collect()
}

pub(crate) fn params_for_derivative_order(singularity_order: u8) -> Vec<Atom> {
    let f = symbol!("f");
    let eta = symbol!("η");

    let f_0 = function!(f, GS.rescale_star);
    let eta_1 = function!(eta, GS.rescale_star).derivative(GS.rescale_star);

    let f_derivative_shape = shape_from_cut_cff_index(&CutCFFIndex {
        left_threshold_order: None,
        right_threshold_order: None,
        lu_cut_order: Some(singularity_order as usize),
    });

    let f_parameters =
        ordered_f_derivative_params(f_0, &[GS.rescale_star], f_derivative_shape.as_ref());
    let mut eta_params = vec![eta_1.clone()];

    for _ in 2..=singularity_order {
        let next_eta = eta_params.last().unwrap().derivative(GS.rescale_star);
        eta_params.push(next_eta);
    }

    let mut result = vec![];
    result.extend(f_parameters);
    result.extend(eta_params);
    result
}

pub(crate) fn params_for_iterated_threshold_ct(
    left_singularit_order: u8,
    right_singularity_order: u8,
) -> Vec<Atom> {
    let f = symbol!("f");
    let eta_left = symbol!("η_left");
    let eta_right = symbol!("η_right");

    let cut_cff_index = CutCFFIndex {
        left_threshold_order: Some(left_singularit_order as usize),
        right_threshold_order: Some(right_singularity_order as usize),
        lu_cut_order: None,
    };
    let f_derivative_shape = shape_from_cut_cff_index(&cut_cff_index);

    let f_base = function!(f, GS.radius_star_left, GS.radius_star_right);
    let derivative_vars = match (left_singularit_order > 1, right_singularity_order > 1) {
        (true, true) => vec![GS.radius_star_left, GS.radius_star_right],
        (true, false) => vec![GS.radius_star_left],
        (false, true) => vec![GS.radius_star_right],
        (false, false) => vec![],
    };

    let eta_left_d1 = function!(eta_left, GS.radius_star_left).derivative(GS.radius_star_left);
    let eta_right_d1 = function!(eta_right, GS.radius_star_right).derivative(GS.radius_star_right);

    let mut eta_left_params = vec![eta_left_d1.clone()];
    let mut eta_right_params = vec![eta_right_d1.clone()];

    for _ in 2..=left_singularit_order {
        let next_eta_left = eta_left_params
            .last()
            .unwrap()
            .derivative(GS.radius_star_left);

        eta_left_params.push(next_eta_left);
    }

    for _ in 2..=right_singularity_order {
        let next_eta_right = eta_right_params
            .last()
            .unwrap()
            .derivative(GS.radius_star_right);

        eta_right_params.push(next_eta_right);
    }

    let f_product =
        ordered_f_derivative_params(f_base, &derivative_vars, f_derivative_shape.as_ref());

    let mut result = vec![];
    result.extend(f_product);
    result.extend(eta_left_params);
    result.extend(eta_right_params);
    result
}

#[cfg(test)]
mod tests {
    use std::{
        fs,
        path::PathBuf,
        time::{SystemTime, UNIX_EPOCH},
    };

    use symbolica::{atom::AtomCore, function, symbol};

    use crate::{
        cff::CutCFFIndex, dot, graph::parse::from_dot::IntoGraph, initialisation::test_initialise,
        utils::GS,
    };
    use linnet::half_edge::{
        involution::EdgeIndex,
        subgraph::{OrientedCut, SuBitGraph},
    };

    fn threshold_association(
        cut_boundary_size: usize,
        threshold_boundary_size: usize,
    ) -> super::ThresholdCountertermAssociation {
        super::ThresholdCountertermAssociation {
            esurface_id: crate::cff::esurface::EsurfaceID(0),
            cut_boundary_edges: (0..cut_boundary_size).map(EdgeIndex::from).collect(),
            threshold_boundary_edges: (0..threshold_boundary_size)
                .map(|index| EdgeIndex::from(cut_boundary_size + index))
                .collect(),
            invariant_bound_is_applicable: true,
        }
    }

    #[test]
    fn singleton_boundary_classification_uses_generic_mass_hierarchy() {
        use super::ThresholdCountertermStatus;
        use crate::utils::F;

        let association = threshold_association(1, 2);
        let tolerance = F(1.0e-12);

        assert_eq!(
            association.classify_from_invariant_bounds(&F(3.0), &F(4.0), Some(F(3.0)), &tolerance,),
            ThresholdCountertermStatus::ProvenNonExisting,
        );
        assert_eq!(
            association.classify_from_invariant_bounds(&F(4.0), &F(4.0), Some(F(4.0)), &tolerance,),
            ThresholdCountertermStatus::AlwaysPinched,
        );
        assert_eq!(
            association.classify_from_invariant_bounds(&F(5.0), &F(4.0), Some(F(5.0)), &tolerance,),
            ThresholdCountertermStatus::PotentiallyExisting,
        );

        let fixed_association = threshold_association(1, 1);
        assert_eq!(
            fixed_association.classify_from_invariant_bounds(
                &F(3.0),
                &F(4.0),
                Some(F(3.0)),
                &tolerance,
            ),
            ThresholdCountertermStatus::ProvenNonExisting,
        );
        assert_eq!(
            fixed_association.classify_from_invariant_bounds(
                &F(5.0),
                &F(4.0),
                Some(F(5.0)),
                &tolerance,
            ),
            ThresholdCountertermStatus::ProvenNonExisting,
        );
    }

    #[test]
    fn current_model_filter_is_controlled_only_by_generation_setting() {
        use super::ThresholdCountertermStatus;

        assert!(ThresholdCountertermStatus::ProvenNonExisting.is_eligible_for_generation(false));
        assert!(!ThresholdCountertermStatus::ProvenNonExisting.is_eligible_for_generation(true));
        assert!(!ThresholdCountertermStatus::AlwaysPinched.is_eligible_for_generation(false));
        assert!(!ThresholdCountertermStatus::NoRadialDependence.is_eligible_for_generation(false));
        assert!(ThresholdCountertermStatus::PotentiallyExisting.is_eligible_for_generation(true));
    }

    #[test]
    fn multiparticle_boundary_classification_is_conservative() {
        use super::ThresholdCountertermStatus;
        use crate::utils::F;

        let association = threshold_association(2, 2);
        let tolerance = F(1.0e-12);

        assert_eq!(
            association.classify_from_invariant_bounds(&F(4.0), &F(4.0), Some(F(8.0)), &tolerance,),
            ThresholdCountertermStatus::CanBecomePinched,
        );
        assert_eq!(
            association.classify_from_invariant_bounds(&F(3.0), &F(5.0), Some(F(4.0)), &tolerance,),
            ThresholdCountertermStatus::ProvenNonExisting,
        );
        assert_eq!(
            association.classify_from_invariant_bounds(&F(3.0), &F(5.0), Some(F(5.0)), &tolerance,),
            ThresholdCountertermStatus::AlwaysPinched,
        );
        assert_eq!(
            association.classify_from_invariant_bounds(&F(3.0), &F(5.0), Some(F(8.0)), &tolerance,),
            ThresholdCountertermStatus::CanBecomePinched,
        );
        assert_eq!(
            association.classify_from_invariant_bounds(&F(5.0), &F(3.0), Some(F(8.0)), &tolerance,),
            ThresholdCountertermStatus::PotentiallyExisting,
        );
        assert_eq!(
            association.classify_from_invariant_bounds(&F(3.0), &F(5.0), None, &tolerance,),
            ThresholdCountertermStatus::PotentiallyExisting,
        );
    }

    #[test]
    fn invariant_bound_rejects_mixed_or_unresolved_inputs_conservatively() {
        use super::ThresholdCountertermStatus;
        use crate::utils::F;

        let mut association = threshold_association(2, 2);
        let tolerance = F(1.0e-12);
        association.invariant_bound_is_applicable = false;
        assert_eq!(
            association.classify_from_invariant_bounds(&F(3.0), &F(5.0), Some(F(4.0)), &tolerance,),
            ThresholdCountertermStatus::PotentiallyExisting,
        );

        association.invariant_bound_is_applicable = true;
        assert_eq!(
            association.classify_from_invariant_bounds(
                &F(f64::NAN),
                &F(5.0),
                Some(F(4.0)),
                &tolerance,
            ),
            ThresholdCountertermStatus::PotentiallyExisting,
        );

        let masses: linnet::half_edge::involution::EdgeVec<Option<F<f64>>> =
            vec![Some(F(3.0)), None].into();
        assert_eq!(
            super::ThresholdCountertermAssociation::mass_sum(
                &[EdgeIndex::from(0), EdgeIndex::from(1)],
                &masses,
            ),
            None,
        );
    }

    #[test]
    fn singleton_threshold_must_lie_strictly_inside_available_interval() {
        use super::ThresholdCountertermStatus;
        use crate::utils::F;

        let association = threshold_association(2, 1);
        let tolerance = F(1.0e-12);

        assert_eq!(
            association.classify_from_invariant_bounds(&F(5.0), &F(4.0), Some(F(9.0)), &tolerance,),
            ThresholdCountertermStatus::ProvenNonExisting,
        );
        assert_eq!(
            association.classify_from_invariant_bounds(&F(3.0), &F(5.0), Some(F(9.0)), &tolerance,),
            ThresholdCountertermStatus::PotentiallyExisting,
        );
        assert_eq!(
            association.classify_from_invariant_bounds(&F(3.0), &F(9.0), Some(F(9.0)), &tolerance,),
            ThresholdCountertermStatus::AlwaysPinched,
        );
    }

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
    fn cross_section_storage_path_stays_inside_process_folder() {
        let temp = fresh_temp_dir("cross-section-storage-path");
        let cross_section = super::CrossSection::new("NLO".to_string());

        assert_eq!(cross_section.storage_path(&temp), temp.join("NLO"));
        fs::remove_dir_all(temp).unwrap();
    }

    #[test]
    fn cross_section_lmb_cut_edge_exclusion_prefers_massive_then_fermion_then_id() {
        test_initialise().unwrap();
        let graph = dot!(
            digraph cut_edge_priority {
                edge [num=1]
                node [num=1]
                A -> B [id=0 particle="g"]
                A -> B [id=1 particle="d"]
                A -> B [id=2 particle="t" mass=1]
                A -> B [id=3 particle="a"]
                A -> B [id=4 mass=0]
            }
        )
        .unwrap();

        assert_eq!(
            super::CrossSectionGraph::excluded_cut_edge_for_lmb_channel_in(
                &graph,
                &[EdgeIndex::from(0), EdgeIndex::from(1), EdgeIndex::from(2)]
            ),
            Some(EdgeIndex::from(2))
        );
        assert_eq!(
            super::CrossSectionGraph::excluded_cut_edge_for_lmb_channel_in(
                &graph,
                &[EdgeIndex::from(0), EdgeIndex::from(1), EdgeIndex::from(3)]
            ),
            Some(EdgeIndex::from(1))
        );
        assert_eq!(
            super::CrossSectionGraph::excluded_cut_edge_for_lmb_channel_in(
                &graph,
                &[EdgeIndex::from(0), EdgeIndex::from(3)]
            ),
            Some(EdgeIndex::from(0))
        );

        // A structural zero and a DOT-overridden mass are resolved information for current-model
        // invariant bounds. Keep genuinely unresolved and complex masses conservative, but do not
        // let a massless edge disable a proof independent of model parameters.
        assert_eq!(
            graph.underlying[EdgeIndex::from(2)]
                .particle
                .mass_atom()
                .to_string(),
            "1"
        );
        let empty: SuBitGraph = graph.underlying.empty_subgraph();
        let cut = super::CrossSectionCut {
            cut: OrientedCut {
                left: empty.clone(),
                right: empty.clone(),
            },
            left: empty.clone(),
            right: empty,
        };
        let association = super::ThresholdCountertermAssociation {
            esurface_id: crate::cff::esurface::EsurfaceID(0),
            cut_boundary_edges: vec![EdgeIndex::from(4)],
            threshold_boundary_edges: vec![EdgeIndex::from(2)],
            invariant_bound_is_applicable: true,
        };
        assert_eq!(
            association.classify_for_model(
                &graph,
                &cut,
                &crate::model::Model::default(),
                &graph.param_builder,
                &crate::settings::RuntimeSettings::default(),
                1.0e-7,
            ),
            super::ThresholdCountertermStatus::ProvenNonExisting,
        );
    }

    #[test]
    fn iterated_counterterms_use_the_right_threshold_count_as_row_stride() {
        let two_by_three = super::IteratedCtCollection::new((0..6).collect(), 2, 3);
        assert_eq!(
            two_by_three[(
                super::LeftThresholdId::from(0),
                super::RightThresholdId::from(0)
            )],
            0
        );
        assert_eq!(
            two_by_three[(
                super::LeftThresholdId::from(0),
                super::RightThresholdId::from(2)
            )],
            2
        );
        assert_eq!(
            two_by_three[(
                super::LeftThresholdId::from(1),
                super::RightThresholdId::from(0)
            )],
            3
        );
        assert_eq!(
            two_by_three[(
                super::LeftThresholdId::from(1),
                super::RightThresholdId::from(2)
            )],
            5
        );

        let three_by_two = super::IteratedCtCollection::new((0..6).collect(), 3, 2);
        assert_eq!(
            three_by_two[(
                super::LeftThresholdId::from(0),
                super::RightThresholdId::from(1)
            )],
            1
        );
        assert_eq!(
            three_by_two[(
                super::LeftThresholdId::from(1),
                super::RightThresholdId::from(0)
            )],
            2
        );
        assert_eq!(
            three_by_two[(
                super::LeftThresholdId::from(2),
                super::RightThresholdId::from(1)
            )],
            5
        );
    }

    #[test]
    fn max_dual_size_for_cut_cff_indices_tracks_mixed_threshold_shapes() {
        let cut_cff_indices = [
            CutCFFIndex {
                left_threshold_order: None,
                right_threshold_order: None,
                lu_cut_order: Some(2),
            },
            CutCFFIndex {
                left_threshold_order: Some(2),
                right_threshold_order: Some(2),
                lu_cut_order: Some(2),
            },
        ];

        assert_eq!(
            super::max_dual_size_for_cut_cff_indices(cut_cff_indices.iter()),
            8
        );
    }

    #[test]
    fn single_threshold_params_follow_effective_f_then_eta_contract() {
        let params = super::params_for_derivative_order(3);
        let f = symbol!("f");
        let eta = symbol!("η");
        let f_base = function!(f, GS.rescale_star);
        let eta_base = function!(eta, GS.rescale_star);

        let expected = [
            f_base.clone(),
            f_base.clone().derivative(GS.rescale_star),
            f_base
                .derivative(GS.rescale_star)
                .derivative(GS.rescale_star),
            eta_base.clone().derivative(GS.rescale_star),
            eta_base
                .clone()
                .derivative(GS.rescale_star)
                .derivative(GS.rescale_star),
            eta_base
                .derivative(GS.rescale_star)
                .derivative(GS.rescale_star)
                .derivative(GS.rescale_star),
        ]
        .into_iter()
        .map(|atom| atom.to_string())
        .collect::<Vec<_>>();

        let actual = params
            .iter()
            .map(|atom| atom.to_string())
            .collect::<Vec<_>>();

        assert_eq!(actual, expected);
    }

    #[test]
    fn iterated_threshold_f_params_follow_mixed_left_right_shape_order() {
        let params = super::params_for_iterated_threshold_ct(2, 2);
        let f = symbol!("f");
        let base = function!(f, GS.radius_star_left, GS.radius_star_right);
        let expected = [
            base.clone(),
            base.clone().derivative(GS.radius_star_left),
            base.clone().derivative(GS.radius_star_right),
            base.derivative(GS.radius_star_left)
                .derivative(GS.radius_star_right),
        ]
        .into_iter()
        .map(|atom| atom.to_string())
        .collect::<Vec<_>>();
        let actual = params[..expected.len()]
            .iter()
            .map(|atom| atom.to_string())
            .collect::<Vec<_>>();

        assert_eq!(actual, expected);
    }

    #[test]
    fn iterated_threshold_params_append_left_then_right_eta_families() {
        let params = super::params_for_iterated_threshold_ct(2, 2);
        let eta_left = symbol!("η_left");
        let eta_right = symbol!("η_right");

        let expected = [
            function!(eta_left, GS.radius_star_left).derivative(GS.radius_star_left),
            function!(eta_left, GS.radius_star_left)
                .derivative(GS.radius_star_left)
                .derivative(GS.radius_star_left),
            function!(eta_right, GS.radius_star_right).derivative(GS.radius_star_right),
            function!(eta_right, GS.radius_star_right)
                .derivative(GS.radius_star_right)
                .derivative(GS.radius_star_right),
        ]
        .into_iter()
        .map(|atom| atom.to_string())
        .collect::<Vec<_>>();

        let actual = params[4..]
            .iter()
            .map(|atom| atom.to_string())
            .collect::<Vec<_>>();

        assert_eq!(actual, expected);
    }

    #[test]
    fn iterated_threshold_f_params_follow_active_single_axis_shape_order() {
        let params = super::params_for_iterated_threshold_ct(1, 3);
        let f = symbol!("f");
        let base = function!(f, GS.radius_star_left, GS.radius_star_right);
        let expected = [
            base.clone(),
            base.clone().derivative(GS.radius_star_right),
            base.derivative(GS.radius_star_right)
                .derivative(GS.radius_star_right),
        ]
        .into_iter()
        .map(|atom| atom.to_string())
        .collect::<Vec<_>>();
        let actual = params[..expected.len()]
            .iter()
            .map(|atom| atom.to_string())
            .collect::<Vec<_>>();

        assert_eq!(actual, expected);
    }
}
