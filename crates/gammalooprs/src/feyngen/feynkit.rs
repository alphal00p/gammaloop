//! GammaLoop policy adapter for standalone FeynKit generation.
//!
//! FeynKit owns model-aware diagram generation. This module is the single
//! upward-facing boundary that translates GammaLoop process/settings types and
//! enriches the resulting structural diagrams with GammaLoop graph state.

use std::{collections::BTreeMap, sync::Arc};

use ahash::HashMap;
use feynkit_generator::{
    GenerationControl, GenerationFilter, GenerationOptions, GenerationReport, GenerationResult,
    Generator, NumeratorGrouping, Process,
    SelfEnergyFilterOptions as FeynkitSelfEnergyFilterOptions,
    SewnFilterOptions as FeynkitSewnFilterOptions, SnailFilterOptions as FeynkitSnailFilterOptions,
    TadpoleFilterOptions as FeynkitTadpoleFilterOptions,
};
use feynkit_graph::{ExternalState, FeynmanDiagram};
use idenso::color::{ColorSimplifier, ColorSimplifySettings};
use itertools::Itertools;
use linnet::{
    half_edge::{
        HedgeGraphError,
        involution::EdgeIndex,
        subgraph::{ModifySubSet, SubSetLike, SubSetOps},
        tree::SimpleTraversalTree,
    },
    permutation::Permutation,
};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::finite_field::PrimeIteratorU64,
    function,
    graph::Graph as SymbolicaGraph,
    id::Replacement,
    symbol,
};
use thiserror::Error;
use tracing::debug;

use super::{
    FeynGenError, FeynGenFilter, GenerationType, NumeratorAwareGraphGroupingOption,
    diagram_generator::{
        EdgeColor, NodeColorWithVertexRule, NodeColorWithoutVertexRule,
        ProcessedNumeratorForComparison,
    },
};
use crate::graph::parse::StripParse;
use crate::{
    graph::{FeynmanGraph, Graph, LMBext, LmbError, parse::ParseGraph},
    integrands::process::ParamBuilder,
    model::{ArcVertexRule, ColorStructure, Model, VertexRule},
    momentum::sample::LoopIndex,
    numerator::{NumeratorStateError, symbolica_ext::NumeratorAtomExt},
    processes::ProcessDefinition,
    settings::GlobalSettings,
    uv::UltravioletGraph,
};

/// Errors raised while crossing from GammaLoop application policy into
/// standalone FeynKit generation.
#[derive(Debug, Error)]
pub enum FeynkitAdapterError {
    #[error("generation interrupted by user")]
    Interrupted,
    #[error("failed to serialize GammaLoop model '{model}' as normalized JSON: {source}")]
    ModelSerialization {
        model: String,
        #[source]
        source: serde_json::Error,
    },
    #[error("normalized GammaLoop model is not a valid FeynKit model: {0}")]
    Model(#[from] feynkit_model::ModelError),
    #[error(transparent)]
    Process(#[from] feynkit_generator::ProcessError),
    #[error(transparent)]
    Generation(#[from] feynkit_generator::GenerationError),
    #[error("{generation_type} process has no final-state specification")]
    MissingFinalState { generation_type: GenerationType },
    #[error("{generation_type} process has inactive {filter_set} filters")]
    InactiveFilterSet {
        generation_type: GenerationType,
        filter_set: &'static str,
    },
    #[error(
        "amplitude symmetrization combination initial={initial}, final={final_state}, left_right={left_right} is invalid"
    )]
    InvalidAmplitudeSymmetrization {
        initial: bool,
        final_state: bool,
        left_right: bool,
    },
    #[error("GammaLoop filter '{filter}' is not meaningful as a cut-amplitude filter")]
    UnsupportedCutAmplitudeFilter { filter: String },
    #[error("diagram '{diagram}' references unknown GammaLoop interaction '{interaction}'")]
    UnknownInteraction {
        diagram: String,
        interaction: String,
    },
    #[error("process external PDG code {pdg} cannot fit in GammaLoop's particle lookup")]
    ExternalPdgOutOfRange { pdg: i64 },
    #[error("process references unknown external particle with PDG code {pdg}: {source}")]
    UnknownExternalParticle {
        pdg: i64,
        #[source]
        source: color_eyre::Report,
    },
    #[error("diagram '{diagram}' vertex {vertex} has neither interaction nor external-leg data")]
    UnclassifiedVertex { diagram: String, vertex: usize },
    #[error(
        "diagram '{diagram}' external leg {index} is {actual:?}, but the GammaLoop process expects {expected:?}"
    )]
    InvalidExternalState {
        diagram: String,
        index: usize,
        actual: ExternalState,
        expected: ExternalState,
    },
    #[error("diagram '{diagram}' is missing external leg {index}")]
    MissingExternalLeg { diagram: String, index: usize },
    #[error("diagram '{diagram}' contains unexpected external leg {index}")]
    UnexpectedExternalLeg { diagram: String, index: usize },
    #[error("diagram '{diagram}' external index {index} cannot fit in GammaLoop's i32 tag")]
    ExternalIndexOutOfRange { diagram: String, index: usize },
    #[error("diagram '{diagram}' edge {edge} has PDG code {pdg}, which cannot fit in isize")]
    PdgOutOfRange {
        diagram: String,
        edge: usize,
        pdg: i64,
    },
    #[error("diagram '{diagram}' edge {edge} references an unknown GammaLoop particle: {source}")]
    UnknownParticle {
        diagram: String,
        edge: usize,
        #[source]
        source: color_eyre::Report,
    },
    #[error("failed to construct Symbolica edge {edge} for diagram '{diagram}': {message}")]
    SymbolicaGraph {
        diagram: String,
        edge: usize,
        message: String,
    },
    #[error("failed to parse overall factor for diagram '{diagram}': {source}")]
    OverallFactor {
        diagram: String,
        #[source]
        source: color_eyre::Report,
    },
    #[error("failed to convert FeynKit diagram '{diagram}' into a GammaLoop graph: {source}")]
    GraphConversion {
        diagram: String,
        #[source]
        source: color_eyre::Report,
    },
    #[error("GammaLoop post-generation policy failed for diagram '{diagram}': {message}")]
    GammaLoopPolicy { diagram: String, message: String },
    #[error("failed to construct the numerator for diagram '{diagram}': {source}")]
    Numerator {
        diagram: String,
        #[source]
        source: NumeratorStateError,
    },
    #[error("explicit loop-momentum basis for graph '{graph}' references unknown edge {edge}")]
    UnknownLoopMomentumEdge { graph: String, edge: usize },
    #[error("explicit loop-momentum basis for graph '{graph}' leaves no spanning-tree edge")]
    EmptyLoopMomentumTree { graph: String },
    #[error("failed to traverse the spanning tree for graph '{graph}': {source}")]
    LoopMomentumTreeTraversal {
        graph: String,
        #[source]
        source: HedgeGraphError,
    },
    #[error("failed to build loop-momentum basis for graph '{graph}': {source}")]
    LoopMomentumBasis {
        graph: String,
        #[source]
        source: LmbError,
    },
    #[error("initial-state edge {edge} is not a loop edge in graph '{graph}'")]
    InitialStateLoopEdge { graph: String, edge: usize },
    #[error("requested edge {edge} is not a loop edge in graph '{graph}'")]
    RequestedLoopEdge { graph: String, edge: usize },
    #[error("diagram index {diagram} cannot fit in GammaLoop's symbolic integer representation")]
    GroupDiagramIndexOutOfRange { diagram: usize },
}

struct PreparedDiagram {
    graph_id: usize,
    topology: SymbolicaGraph<NodeColorWithoutVertexRule, String>,
    external_orbit_key: Option<SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>>,
    external_ordering: Vec<usize>,
    ordering_key: SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
    comparison_graph: Graph,
    representative_graph: Graph,
}

type ExternalCanonicalizationPolicy = (HashMap<i32, i32>, Option<(bool, bool)>);

fn combine_prepared_factors(representative: &mut PreparedDiagram, diagram: &PreparedDiagram) {
    let symmetry_pattern = function!(
        symbol!("NumeratorIndependentSymmetryGrouping"),
        Atom::var(symbol!("x_"))
    )
    .to_pattern();
    let denominator = representative
        .representative_graph
        .overall_factor
        .replace(&symmetry_pattern)
        .with(Atom::num(1).to_pattern());
    let ratio = (super::diagram_generator::evaluate_overall_factor(
        diagram.representative_graph.overall_factor.as_view(),
    ) / super::diagram_generator::evaluate_overall_factor(denominator.as_view()))
    .expand();
    let factor = if representative
        .representative_graph
        .overall_factor
        .pattern_match(&symmetry_pattern, None, None)
        .next()
        .is_some()
    {
        representative
            .representative_graph
            .overall_factor
            .replace(&symmetry_pattern)
            .with(
                function!(
                    symbol!("NumeratorIndependentSymmetryGrouping"),
                    (Atom::var(symbol!("x_")) + ratio).expand()
                )
                .to_pattern(),
            )
    } else {
        &representative.representative_graph.overall_factor
            * function!(
                symbol!("NumeratorIndependentSymmetryGrouping"),
                (Atom::num(1) + ratio).expand()
            )
    };
    representative.comparison_graph.overall_factor = factor.clone();
    representative.representative_graph.overall_factor = factor;
}

struct PooledGraphData {
    graph_id: usize,
    numerator_data: Option<ProcessedNumeratorForComparison>,
    ratio: Atom,
    graph: Graph,
}

/// GammaLoop graphs together with the standalone generation report.
pub struct FeynkitAdapterResult {
    pub graphs: Vec<Graph>,
    pub report: GenerationReport,
}

/// Typed, downward-only bridge from GammaLoop policy to FeynKit.
pub struct FeynkitGeneratorAdapter<'a> {
    definition: &'a ProcessDefinition,
    model: &'a Model,
    settings: &'a GlobalSettings,
}

impl ProcessDefinition {
    /// Generate GammaLoop graphs through the standalone structural generator,
    /// then apply GammaLoop-only numerator and loop-momentum policy locally.
    #[tracing::instrument(skip_all)]
    pub fn generate(
        &self,
        model: &Model,
        settings: &GlobalSettings,
    ) -> Result<Vec<Graph>, FeynGenError> {
        match FeynkitGeneratorAdapter::new(self, model, settings).generate() {
            Ok(result) => Ok(result.graphs),
            Err(FeynkitAdapterError::Interrupted) => Err(FeynGenError::Interrupted),
            Err(error) => Err(FeynGenError::Feynkit(error)),
        }
    }
}

impl<'a> FeynkitGeneratorAdapter<'a> {
    pub fn new(
        definition: &'a ProcessDefinition,
        model: &'a Model,
        settings: &'a GlobalSettings,
    ) -> Self {
        Self {
            definition,
            model,
            settings,
        }
    }

    /// Convert GammaLoop's normalized serializable model at the application
    /// boundary. Executable UFO handling remains outside this adapter.
    pub fn normalized_model(&self) -> Result<feynkit_model::Model, FeynkitAdapterError> {
        let json = serde_json::to_string(&self.model.to_serializable()).map_err(|source| {
            FeynkitAdapterError::ModelSerialization {
                model: self.model.name.to_string(),
                source,
            }
        })?;
        Ok(feynkit_model::Model::from_json(&json)?)
    }

    pub fn process(&self) -> Result<Process, FeynkitAdapterError> {
        let Some(first_final_state) = self.definition.final_pdgs_lists.first() else {
            return Err(FeynkitAdapterError::MissingFinalState {
                generation_type: self.definition.generation_type,
            });
        };
        self.validate_amplitude_symmetrization()?;
        let process = match self.definition.generation_type {
            GenerationType::Amplitude => Process::amplitude(
                self.definition.initial_pdgs.clone(),
                first_final_state.clone(),
            )
            .with_final_state_alternatives(self.definition.final_pdgs_lists.clone())?,
            GenerationType::CrossSection => Process::cross_section(
                self.definition.initial_pdgs.clone(),
                first_final_state.clone(),
            )
            .with_final_state_alternatives(self.definition.final_pdgs_lists.clone())?,
        }
        .with_loop_count(
            self.definition.loop_count_range.0,
            self.definition.loop_count_range.1,
        )?
        .symmetrize_initial(self.definition.symmetrize_initial_states)
        .symmetrize_final(self.definition.symmetrize_final_states)
        .symmetrize_left_right(self.definition.symmetrize_left_right_states)
        .symmetrize_external_fermions(
            self.definition
                .allow_symmetrization_of_external_fermions_in_amplitudes,
        );
        Ok(process)
    }

    fn validate_amplitude_symmetrization(&self) -> Result<(), FeynkitAdapterError> {
        if self.definition.generation_type != GenerationType::Amplitude {
            return Ok(());
        }
        let symmetrization = (
            self.definition.symmetrize_initial_states,
            self.definition.symmetrize_final_states,
            self.definition.symmetrize_left_right_states,
        );
        if !matches!(
            symmetrization,
            (false, false, false)
                | (true, false, false)
                | (false, true, false)
                | (true, true, false)
                | (true, true, true)
        ) {
            return Err(FeynkitAdapterError::InvalidAmplitudeSymmetrization {
                initial: symmetrization.0,
                final_state: symmetrization.1,
                left_right: symmetrization.2,
            });
        }
        Ok(())
    }

    pub fn options(&self) -> Result<GenerationOptions, FeynkitAdapterError> {
        let mut options = GenerationOptions::default()
            .threads(self.settings.n_cores.feyngen)
            .allow_self_loops(!self.definition.filter_self_loop)
            .allow_zero_flow_edges(!self.definition.filter_zero_flow_edges)
            // Tensor-aware numerator policy remains in GammaLoop. FeynKit is
            // responsible only for structural diagram generation here.
            .numerator_grouping(NumeratorGrouping::None)
            .graph_prefix(&self.definition.graph_prefix)
            .cancellation_check(crate::is_interrupt_requested)
            .progress(|_| {
                if crate::is_interrupt_requested() {
                    GenerationControl::Cancel
                } else {
                    GenerationControl::Continue
                }
            });

        match self.definition.generation_type {
            GenerationType::Amplitude => {
                if !self.definition.cross_section_filters.0.is_empty() {
                    return Err(FeynkitAdapterError::InactiveFilterSet {
                        generation_type: self.definition.generation_type,
                        filter_set: "cross-section",
                    });
                }
                for filter in &self.definition.amplitude_filters.0 {
                    options = options.with_graph_filter(Self::filter(filter)?);
                }
            }
            GenerationType::CrossSection => {
                for filter in &self.definition.cross_section_filters.0 {
                    options = options.with_graph_filter(Self::filter(filter)?);
                }
                for filter in &self.definition.amplitude_filters.0 {
                    let filter = match filter {
                        FeynGenFilter::CouplingOrders(_) | FeynGenFilter::LoopCountRange(_) => {
                            Self::filter(filter)?
                        }
                        filter => {
                            return Err(FeynkitAdapterError::UnsupportedCutAmplitudeFilter {
                                filter: filter.to_string(),
                            });
                        }
                    };
                    options = options.with_cut_amplitude_filter(filter);
                }
            }
        }
        Ok(options)
    }

    pub fn generate(&self) -> Result<FeynkitAdapterResult, FeynkitAdapterError> {
        let result = (|| {
            let normalized_model = self.normalized_model()?;
            // GammaLoop applies its external-state orbit policy after all
            // structural filters. Generate labeled assignments here so paired
            // cross-section permutations and left-right swaps retain their
            // complete multiplicities.
            let process = self
                .process()?
                .symmetrize_initial(false)
                .symmetrize_final(false)
                .symmetrize_left_right(false)
                .symmetrize_external_fermions(false);
            let options = self.options()?;
            let generated = Generator::new(normalized_model).generate(&process, &options)?;
            let GenerationResult {
                diagrams,
                mut report,
                ..
            } = generated;
            if !report.completed {
                return Err(FeynkitAdapterError::Interrupted);
            }
            let external_connections = self.external_connections();
            let mut prepared = Vec::with_capacity(diagrams.len());

            for (graph_id, diagram) in diagrams.iter().enumerate() {
                if crate::is_interrupt_requested() {
                    return Err(FeynkitAdapterError::Interrupted);
                }
                prepared.push(self.prepare_diagram(
                    graph_id,
                    diagram,
                    diagram.name(),
                    &external_connections,
                )?);
            }

            // Preserve the legacy amplitude pregrouping key, which includes
            // complete vertex-rule colors. Skeleton canonicalization alone can
            // otherwise conflate automorphic assignments before the external
            // orbit representative is selected.
            let mut external_orbits: HashMap<_, Vec<_>> = HashMap::default();
            let mut externally_grouped = Vec::with_capacity(prepared.len());
            for mut diagram in prepared {
                if crate::is_interrupt_requested() {
                    return Err(FeynkitAdapterError::Interrupted);
                }
                if let Some(key) = diagram.external_orbit_key.take() {
                    external_orbits.entry(key).or_default().push(diagram);
                } else {
                    externally_grouped.push(diagram);
                }
            }
            for mut orbit in external_orbits.into_values() {
                if crate::is_interrupt_requested() {
                    return Err(FeynkitAdapterError::Interrupted);
                }
                orbit.sort_by(|left, right| left.external_ordering.cmp(&right.external_ordering));
                let mut orbit = orbit.into_iter();
                if let Some(mut representative) = orbit.next() {
                    for diagram in orbit {
                        if crate::is_interrupt_requested() {
                            return Err(FeynkitAdapterError::Interrupted);
                        }
                        combine_prepared_factors(&mut representative, &diagram);
                    }
                    externally_grouped.push(representative);
                }
            }
            let mut prepared = externally_grouped;

            // Legacy GammaLoop assigns graph IDs only after canonical flow
            // normalization. Preserve that boundary so numerator-group masters,
            // graph selection, and explicit LMB names remain stable even when
            // FeynKit's structural generation order differs.
            prepared.sort_by(|left, right| left.ordering_key.cmp(&right.ordering_key));

            // Combine duplicates after canonical flow normalization. They can
            // arise when externally symmetrized assignments select the same
            // representative, and their numerator-independent factors must be
            // accumulated before stable GammaLoop graph IDs are assigned.
            let mut combined: Vec<PreparedDiagram> = Vec::with_capacity(prepared.len());
            for diagram in prepared {
                if crate::is_interrupt_requested() {
                    return Err(FeynkitAdapterError::Interrupted);
                }
                if let Some(representative) = combined
                    .last_mut()
                    .filter(|representative| representative.ordering_key == diagram.ordering_key)
                {
                    combine_prepared_factors(representative, &diagram);
                } else {
                    combined.push(diagram);
                }
            }
            let mut prepared = combined;
            let padding_width = prepared.len().saturating_sub(1).to_string().len().max(1);
            for (graph_id, diagram) in prepared.iter_mut().enumerate() {
                if crate::is_interrupt_requested() {
                    return Err(FeynkitAdapterError::Interrupted);
                }
                let name = format!(
                    "{}{:0padding_width$}",
                    self.definition.graph_prefix, graph_id
                );
                diagram.graph_id = graph_id;
                diagram.comparison_graph.name.clone_from(&name);
                diagram.representative_graph.name = name;
            }
            prepared.retain(|diagram| {
                let name = &diagram.representative_graph.name;
                !self
                    .definition
                    .selected_graphs
                    .as_ref()
                    .is_some_and(|selected| !selected.contains(name))
                    && !self
                        .definition
                        .vetoed_graphs
                        .as_ref()
                        .is_some_and(|vetoed| vetoed.contains(name))
            });

            let (mut graphs, zero_numerator_count) = self.postprocess_numerators(prepared)?;
            report.retained_count = graphs.len();
            report.zero_numerator_count = zero_numerator_count;
            for graph in &mut graphs {
                if crate::is_interrupt_requested() {
                    return Err(FeynkitAdapterError::Interrupted);
                }
                if let Some(edges) = self
                    .definition
                    .loop_momentum_bases
                    .as_ref()
                    .and_then(|bases| bases.get(&graph.name))
                {
                    self.apply_loop_momentum_basis(graph, edges)?;
                }
            }
            if crate::is_interrupt_requested() {
                return Err(FeynkitAdapterError::Interrupted);
            }

            Ok(FeynkitAdapterResult { graphs, report })
        })();
        crate::clear_interrupt_request();
        result
    }

    fn postprocess_numerators(
        &self,
        diagrams: Vec<PreparedDiagram>,
    ) -> Result<(Vec<Graph>, usize), FeynkitAdapterError> {
        if matches!(
            self.definition.numerator_grouping,
            NumeratorAwareGraphGroupingOption::NoGrouping
        ) {
            let mut graphs = Vec::with_capacity(diagrams.len());
            for diagram in diagrams {
                if crate::is_interrupt_requested() {
                    return Err(FeynkitAdapterError::Interrupted);
                }
                if !super::diagram_generator::evaluate_overall_factor(
                    diagram.representative_graph.overall_factor.as_view(),
                )
                .expand()
                .is_zero()
                {
                    graphs.push(diagram.representative_graph);
                }
            }
            return Ok((graphs, 0));
        }

        let samples = if let Some(options) = self.definition.numerator_grouping.get_options() {
            let mut sample_iterator = PrimeIteratorU64::new(1);
            sample_iterator.nth(options.numerical_sample_seed as usize);
            (0..options.number_of_numerical_samples)
                .map(|_| {
                    self.definition.sample_lib(
                        &mut sample_iterator,
                        options.fully_numerical_substitution_when_comparing_numerators,
                        options.symmetric_polarizations,
                        self.model,
                    )
                })
                .collect::<Vec<_>>()
        } else {
            Vec::new()
        };
        let coupling_replacements = self
            .model
            .couplings
            .values()
            .map(|coupling| {
                let [left, right] = coupling.rep_rule();
                Replacement::new(left.to_pattern(), right)
            })
            .collect::<Vec<_>>();
        let mut pools: HashMap<
            SymbolicaGraph<NodeColorWithoutVertexRule, String>,
            Vec<Vec<PooledGraphData>>,
        > = HashMap::default();
        let mut zero_numerator_count = 0;

        for diagram in diagrams {
            if crate::is_interrupt_requested() {
                return Err(FeynkitAdapterError::Interrupted);
            }
            let graph_name = diagram.comparison_graph.name.clone();
            let mut numerator = diagram.comparison_graph.numerator(
                &diagram.comparison_graph.no_dummy(),
                &diagram.comparison_graph.empty_subgraph(),
            );
            // The graph symmetry factor is deliberately excluded. Grouping
            // compares numerators and combines each member's factor afterwards.
            numerator.state.expr *= &diagram.comparison_graph.global_prefactor.num
                * &diagram.comparison_graph.global_prefactor.projector;
            numerator.state.expr = numerator
                .state
                .expr
                .replace_multiple(&coupling_replacements);

            let numerator_color_simplified = numerator
                .get_single_atom()
                .map_err(|source| FeynkitAdapterError::Numerator {
                    diagram: graph_name.clone(),
                    source,
                })?
                .to_param_color()
                .simplify_color_with(
                    ColorSimplifySettings::default().with_cof_dimension_invariants(),
                );

            // The color part must be expanded to identify exact zeroes.
            if numerator_color_simplified
                .expand_color()
                .iter()
                .fold(Atom::Zero, |sum, (color, lorentz)| color * lorentz + sum)
                .is_zero()
            {
                zero_numerator_count += 1;
                continue;
            }

            let topology_groups = pools.entry(diagram.topology).or_default();
            if matches!(
                self.definition.numerator_grouping,
                NumeratorAwareGraphGroupingOption::OnlyDetectZeroes
            ) {
                topology_groups.push(vec![PooledGraphData {
                    graph_id: diagram.graph_id,
                    numerator_data: None,
                    ratio: Atom::num(1),
                    graph: diagram.representative_graph,
                }]);
                continue;
            }

            // Numerator comparison uses the sorted graph before flow
            // normalization. The retained representative uses normalized flow.
            let numerator_data =
                ProcessedNumeratorForComparison::from_numerator_symbolic_expression(
                    diagram.graph_id,
                    &diagram.comparison_graph,
                    numerator_color_simplified,
                    &samples,
                    self.settings,
                    &self.definition.numerator_grouping,
                )
                .map_err(|source| FeynkitAdapterError::GammaLoopPolicy {
                    diagram: graph_name,
                    message: source.to_string(),
                })?;
            if numerator_data.evaluated_samples_are_zero() {
                zero_numerator_count += 1;
                continue;
            }

            let matching_group =
                topology_groups
                    .iter()
                    .enumerate()
                    .find_map(|(group_index, group)| {
                        ProcessDefinition::compare_numerator_tensors(
                            &self.definition.numerator_grouping,
                            &numerator_data,
                            group[0].numerator_data.as_ref()?,
                        )
                        .map(|ratio| (group_index, ratio))
                    });
            if let Some((group_index, ratio)) = matching_group {
                topology_groups[group_index].push(PooledGraphData {
                    graph_id: diagram.graph_id,
                    numerator_data: None,
                    ratio,
                    graph: diagram.representative_graph,
                });
            } else {
                topology_groups.push(vec![PooledGraphData {
                    graph_id: diagram.graph_id,
                    numerator_data: Some(numerator_data),
                    ratio: Atom::num(1),
                    graph: diagram.representative_graph,
                }]);
            }
        }

        Ok((self.collapse_groups(pools)?, zero_numerator_count))
    }

    fn collapse_groups(
        &self,
        pools: HashMap<
            SymbolicaGraph<NodeColorWithoutVertexRule, String>,
            Vec<Vec<PooledGraphData>>,
        >,
    ) -> Result<Vec<Graph>, FeynkitAdapterError> {
        let mut collapsed = Vec::new();
        for topology_groups in pools.into_values() {
            for mut group in topology_groups {
                if crate::is_interrupt_requested() {
                    return Err(FeynkitAdapterError::Interrupted);
                }
                group.sort_by_key(|diagram| diagram.graph_id);
                let reference_ratio = group[0].ratio.clone();
                let reference_id = group[0].graph_id;
                let mut representative = group[0].graph.clone();
                if group.len() > 1 {
                    let mut combined_factor = Atom::num(0);
                    for diagram in group {
                        let diagram_id = i64::try_from(diagram.graph_id).map_err(|_| {
                            FeynkitAdapterError::GroupDiagramIndexOutOfRange {
                                diagram: diagram.graph_id,
                            }
                        })?;
                        // Ratios are measured against the insertion-order
                        // reference. Dividing by the retained reference ratio
                        // keeps this correct if those orders ever differ.
                        combined_factor += function!(
                            symbol!("NumeratorDependentGrouping"),
                            Atom::num(diagram_id),
                            (&diagram.ratio / &reference_ratio).expand(),
                            diagram.graph.overall_factor
                        );
                    }
                    representative.overall_factor = combined_factor;
                }
                if super::diagram_generator::evaluate_overall_factor(
                    representative.overall_factor.as_view(),
                )
                .expand()
                .is_zero()
                {
                    debug!(
                        overall_factor = %representative.overall_factor,
                        "Numerator-aware grouping produced an exact cancellation"
                    );
                } else {
                    collapsed.push((reference_id, representative));
                }
            }
        }
        collapsed.sort_by_key(|(graph_id, _)| *graph_id);
        Ok(collapsed.into_iter().map(|(_, graph)| graph).collect())
    }

    fn filter(filter: &FeynGenFilter) -> Result<GenerationFilter, FeynkitAdapterError> {
        Ok(match filter {
            FeynGenFilter::SelfEnergyFilter(options) => {
                GenerationFilter::SelfEnergy(FeynkitSelfEnergyFilterOptions {
                    veto_massive: options.veto_self_energy_of_massive_lines,
                    veto_massless: options.veto_self_energy_of_massless_lines,
                    only_scaleless: options.veto_only_scaleless_self_energy,
                })
            }
            FeynGenFilter::TadpolesFilter(options) => {
                GenerationFilter::Tadpoles(FeynkitTadpoleFilterOptions {
                    veto_attached_to_massive: options.veto_tadpoles_attached_to_massive_lines,
                    veto_attached_to_massless: options.veto_tadpoles_attached_to_massless_lines,
                    only_scaleless: options.veto_only_scaleless_tadpoles,
                })
            }
            FeynGenFilter::ZeroSnailsFilter(options) => {
                GenerationFilter::ZeroSnails(FeynkitSnailFilterOptions {
                    veto_attached_to_massive: options.veto_snails_attached_to_massive_lines,
                    veto_attached_to_massless: options.veto_snails_attached_to_massless_lines,
                    only_scaleless: options.veto_only_scaleless_snails,
                })
            }
            FeynGenFilter::ParticleVeto(pdgs) => GenerationFilter::ParticleVeto(pdgs.clone()),
            FeynGenFilter::VertexAllow(vertices) => GenerationFilter::VertexAllow(vertices.clone()),
            FeynGenFilter::VertexVeto(vertices) => GenerationFilter::VertexVeto(vertices.clone()),
            FeynGenFilter::MaxNumberOfBridges(maximum) => {
                GenerationFilter::MaxNumberOfBridges(*maximum)
            }
            FeynGenFilter::CouplingOrders(orders) => GenerationFilter::CouplingOrders(
                orders
                    .iter()
                    .map(|(name, range)| (name.clone(), *range))
                    .collect::<BTreeMap<_, _>>(),
            ),
            FeynGenFilter::LoopCountRange(range) => GenerationFilter::LoopCountRange(*range),
            FeynGenFilter::BlobRange(range) => GenerationFilter::BlobRange(range.clone()),
            FeynGenFilter::SpectatorRange(range) => GenerationFilter::SpectatorRange(range.clone()),
            FeynGenFilter::PerturbativeOrders(orders) => GenerationFilter::PerturbativeOrders(
                orders
                    .iter()
                    .map(|(name, order)| (name.clone(), *order))
                    .collect::<BTreeMap<_, _>>(),
            ),
            FeynGenFilter::FermionLoopCountRange(range) => {
                GenerationFilter::FermionLoopCountRange(*range)
            }
            FeynGenFilter::FactorizedLoopTopologiesCountRange(range) => {
                GenerationFilter::FactorizedLoopTopologiesCountRange(*range)
            }
            FeynGenFilter::SewedFilter(options) => {
                GenerationFilter::Sewn(FeynkitSewnFilterOptions {
                    filter_tadpoles: options.filter_tadpoles,
                })
            }
        })
    }

    fn external_connections(&self) -> Vec<(Option<usize>, Option<usize>)> {
        match self.definition.generation_type {
            GenerationType::Amplitude => {
                self.definition
                    .initial_pdgs
                    .iter()
                    .enumerate()
                    .map(|(index, _)| (Some(index + 1), None))
                    .chain(self.definition.final_pdgs_lists[0].iter().enumerate().map(
                        |(index, _)| (None, Some(self.definition.initial_pdgs.len() + index + 1)),
                    ))
                    .collect()
            }
            GenerationType::CrossSection => (0..self.definition.initial_pdgs.len())
                .map(|index| {
                    (
                        Some(index + 1),
                        Some(self.definition.initial_pdgs.len() + index + 1),
                    )
                })
                .collect(),
        }
    }

    fn external_canonicalization_policy(
        &self,
    ) -> Result<ExternalCanonicalizationPolicy, FeynkitAdapterError> {
        let mut external_colors = HashMap::default();
        let manual_ordering = match self.definition.generation_type {
            GenerationType::CrossSection => {
                let manual_ordering = (self.definition.symmetrize_initial_states
                    || self.definition.symmetrize_left_right_states)
                    .then_some((
                        self.definition.symmetrize_initial_states,
                        self.definition.symmetrize_left_right_states,
                    ));
                if manual_ordering.is_some() {
                    for initial_color in 1..=self.definition.initial_pdgs.len() {
                        external_colors.insert(initial_color as i32, -2000 - initial_color as i32);
                    }
                    for final_color in self.definition.initial_pdgs.len() + 1
                        ..=2 * self.definition.initial_pdgs.len()
                    {
                        external_colors.insert(
                            final_color as i32,
                            -3000 - (final_color - self.definition.initial_pdgs.len()) as i32,
                        );
                    }
                }
                manual_ordering
            }
            GenerationType::Amplitude => {
                let (initial_color, final_color) = match (
                    self.definition.symmetrize_initial_states,
                    self.definition.symmetrize_final_states,
                    self.definition.symmetrize_left_right_states,
                ) {
                    (false, false, false) => (None, None),
                    (true, false, false) => (Some(-2), None),
                    (false, true, false) => (None, Some(-3)),
                    (true, true, false) => (Some(-2), Some(-3)),
                    (true, true, true) => (Some(-1), Some(-1)),
                    (initial, final_state, left_right) => {
                        return Err(FeynkitAdapterError::InvalidAmplitudeSymmetrization {
                            initial,
                            final_state,
                            left_right,
                        });
                    }
                };
                let can_symmetrize = |pdg: i64| -> Result<bool, FeynkitAdapterError> {
                    let lookup = isize::try_from(pdg)
                        .map_err(|_| FeynkitAdapterError::ExternalPdgOutOfRange { pdg })?;
                    let particle =
                        self.model
                            .try_get_particle_from_pdg(lookup)
                            .map_err(|source| FeynkitAdapterError::UnknownExternalParticle {
                                pdg,
                                source,
                            })?;
                    Ok(!particle.0.is_fermion()
                        || self
                            .definition
                            .allow_symmetrization_of_external_fermions_in_amplitudes)
                };
                if let Some(color) = initial_color {
                    for (index, pdg) in self.definition.initial_pdgs.iter().copied().enumerate() {
                        if can_symmetrize(pdg)? {
                            external_colors.insert(index as i32 + 1, color);
                        }
                    }
                }
                if let Some(color) = final_color {
                    for (offset, pdg) in self.definition.final_pdgs_lists[0]
                        .iter()
                        .copied()
                        .enumerate()
                    {
                        if can_symmetrize(pdg)? {
                            external_colors.insert(
                                (self.definition.initial_pdgs.len() + offset) as i32 + 1,
                                color,
                            );
                        }
                    }
                }
                None
            }
        };
        Ok((external_colors, manual_ordering))
    }

    fn prepare_diagram(
        &self,
        graph_id: usize,
        diagram: &FeynmanDiagram,
        name: &str,
        external_connections: &[(Option<usize>, Option<usize>)],
    ) -> Result<PreparedDiagram, FeynkitAdapterError> {
        self.validate_external_legs(diagram)?;
        let symbolica_graph = self.symbolica_graph(diagram)?;
        // GammaLoop selects the representative of each labeled external-state
        // orbit here. The first pass makes that choice; the second canonizes
        // vertex ordering for it, preserving legacy graph IDs and
        // numerator-group masters.
        let (external_colors, manual_ordering) = self.external_canonicalization_policy()?;
        let (external_orbit_key, external_ordering) = if self.definition.generation_type
            == GenerationType::Amplitude
            && !external_colors.is_empty()
        {
            let mut recolored = symbolica_graph.clone();
            for (index, node) in symbolica_graph.nodes().iter().enumerate() {
                if let Some(color) = external_colors.get(&node.data.external_tag) {
                    let mut node = node.data.clone();
                    node.external_tag = *color;
                    recolored.set_node_data(index, node);
                }
            }
            let external_ordering = symbolica_graph
                .nodes()
                .iter()
                .filter(|node| node.data.external_tag > 0)
                .filter_map(|node| node.edges.first().copied())
                .collect();
            (Some(recolored.canonize().graph), external_ordering)
        } else {
            (None, Vec::new())
        };
        let (mut topology, mut comparison_symbolica_graph) = self
            .definition
            .canonicalize_edge_and_vertex_ordering(
                self.model,
                &symbolica_graph,
                &external_colors,
                &self.definition.numerator_grouping,
                manual_ordering,
            )
            .map_err(|source| FeynkitAdapterError::GammaLoopPolicy {
                diagram: name.to_owned(),
                message: source.to_string(),
            })?;
        if manual_ordering.is_some() {
            (topology, comparison_symbolica_graph) = self
                .definition
                .canonicalize_edge_and_vertex_ordering(
                    self.model,
                    &comparison_symbolica_graph,
                    &external_colors,
                    &self.definition.numerator_grouping,
                    None,
                )
                .map_err(|source| FeynkitAdapterError::GammaLoopPolicy {
                    diagram: name.to_owned(),
                    message: source.to_string(),
                })?;
        }
        // Normalize fermion directions only for the retained representative.
        // Flipping those edges in the comparison graph would misalign the
        // signs and orientations of its loop-momentum basis during numerator
        // comparison and diagram grouping.
        let (representative_symbolica_graph, _) = self
            .definition
            .normalize_flows(&comparison_symbolica_graph, self.model)
            .map_err(|source| FeynkitAdapterError::GammaLoopPolicy {
                diagram: name.to_owned(),
                message: source.to_string(),
            })?;
        let mut overall_factor =
            diagram
                .overall_factor()
                .strip_parse::<Atom>()
                .map_err(|source| FeynkitAdapterError::OverallFactor {
                    diagram: name.to_owned(),
                    source,
                })?;
        if self.definition.generation_type == GenerationType::Amplitude
            && self
                .definition
                .allow_symmetrization_of_external_fermions_in_amplitudes
            && (self.definition.symmetrize_initial_states
                || self.definition.symmetrize_final_states)
        {
            let external_sign = function!(
                symbol!("ExternalFermionOrderingSign"),
                Atom::var(symbol!("x_"))
            )
            .to_pattern();
            overall_factor = overall_factor
                .replace(&external_sign)
                .with(Atom::num(1).to_pattern());
        }
        let mut comparison_parsed = ParseGraph::from_symbolica_graph(
            self.model,
            name,
            &comparison_symbolica_graph,
            overall_factor.clone(),
            external_connections,
        )
        .map_err(|source| FeynkitAdapterError::GraphConversion {
            diagram: name.to_owned(),
            source,
        })?;

        // GammaLoop projectors and user prefactors are application policy and
        // intentionally do not live in the standalone diagram IR.
        comparison_parsed.global_data.num = self.definition.prefactor.num.clone();
        comparison_parsed.global_data.projectors = (self.definition.prefactor.projector
            != Atom::num(1))
        .then(|| self.definition.prefactor.projector.clone());
        let mut representative_parsed = ParseGraph::from_symbolica_graph(
            self.model,
            name,
            &representative_symbolica_graph,
            overall_factor,
            external_connections,
        )
        .map_err(|source| FeynkitAdapterError::GraphConversion {
            diagram: name.to_owned(),
            source,
        })?;
        representative_parsed.global_data = comparison_parsed.global_data.clone();
        let comparison_graph =
            Graph::from_parsed(comparison_parsed, self.model).map_err(|source| {
                FeynkitAdapterError::GraphConversion {
                    diagram: name.to_owned(),
                    source,
                }
            })?;
        let representative_graph =
            Graph::from_parsed(representative_parsed, self.model).map_err(|source| {
                FeynkitAdapterError::GraphConversion {
                    diagram: name.to_owned(),
                    source,
                }
            })?;
        Ok(PreparedDiagram {
            graph_id,
            topology,
            external_orbit_key,
            external_ordering,
            ordering_key: representative_symbolica_graph,
            comparison_graph,
            representative_graph,
        })
    }

    fn validate_external_legs(&self, diagram: &FeynmanDiagram) -> Result<(), FeynkitAdapterError> {
        let expected = match self.definition.generation_type {
            GenerationType::Amplitude => {
                self.definition.initial_pdgs.len() + self.definition.final_pdgs_lists[0].len()
            }
            GenerationType::CrossSection => 2 * self.definition.initial_pdgs.len(),
        };
        let mut states = vec![None; expected];
        for (_, vertex) in diagram.vertices() {
            let Some(external) = &vertex.external else {
                continue;
            };
            let Some(slot) = states.get_mut(external.index) else {
                return Err(FeynkitAdapterError::UnexpectedExternalLeg {
                    diagram: diagram.name().to_owned(),
                    index: external.index,
                });
            };
            let expected_state = if external.index < self.definition.initial_pdgs.len() {
                ExternalState::Incoming
            } else {
                ExternalState::Outgoing
            };
            if external.state != expected_state {
                return Err(FeynkitAdapterError::InvalidExternalState {
                    diagram: diagram.name().to_owned(),
                    index: external.index,
                    actual: external.state,
                    expected: expected_state,
                });
            }
            *slot = Some(external.state);
        }
        if let Some(index) = states.iter().position(Option::is_none) {
            return Err(FeynkitAdapterError::MissingExternalLeg {
                diagram: diagram.name().to_owned(),
                index,
            });
        }
        Ok(())
    }

    fn symbolica_graph(
        &self,
        diagram: &FeynmanDiagram,
    ) -> Result<SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>, FeynkitAdapterError> {
        let external_rule = ArcVertexRule(Arc::new(VertexRule {
            name: "external".into(),
            couplings: Vec::new(),
            lorentz_structures: Vec::new(),
            particles: Vec::new(),
            color_structures: ColorStructure {
                color_structure: Vec::new(),
            },
            dod: 0,
        }));
        let mut graph = SymbolicaGraph::new();
        for (vertex_id, vertex) in diagram.vertices() {
            let (external_tag, vertex_rule) = if let Some(external) = &vertex.external {
                let external_tag = i32::try_from(external.index)
                    .ok()
                    .and_then(|index| index.checked_add(1))
                    .ok_or_else(|| FeynkitAdapterError::ExternalIndexOutOfRange {
                        diagram: diagram.name().to_owned(),
                        index: external.index,
                    })?;
                (external_tag, external_rule.clone())
            } else if let Some(interaction) = &vertex.interaction {
                let Some(vertex_rule) = self
                    .model
                    .vertex_rules
                    .iter()
                    .find(|rule| rule.name.as_str() == interaction)
                    .cloned()
                else {
                    return Err(FeynkitAdapterError::UnknownInteraction {
                        diagram: diagram.name().to_owned(),
                        interaction: interaction.clone(),
                    });
                };
                (0, vertex_rule)
            } else {
                return Err(FeynkitAdapterError::UnclassifiedVertex {
                    diagram: diagram.name().to_owned(),
                    vertex: vertex_id.0,
                });
            };
            graph.add_node(NodeColorWithVertexRule {
                external_tag,
                vertex_rule,
            });
        }
        for (edge_id, endpoints, edge) in diagram.edges() {
            let pdg = isize::try_from(edge.particle.pdg).map_err(|_| {
                FeynkitAdapterError::PdgOutOfRange {
                    diagram: diagram.name().to_owned(),
                    edge: edge_id.0,
                    pdg: edge.particle.pdg,
                }
            })?;
            self.model
                .try_get_particle_from_pdg(pdg)
                .map_err(|source| FeynkitAdapterError::UnknownParticle {
                    diagram: diagram.name().to_owned(),
                    edge: edge_id.0,
                    source,
                })?;
            graph
                .add_edge(
                    endpoints.source.0,
                    endpoints.target.0,
                    edge.directed,
                    EdgeColor { pdg },
                )
                .map_err(|error| FeynkitAdapterError::SymbolicaGraph {
                    diagram: diagram.name().to_owned(),
                    edge: edge_id.0,
                    message: error.to_string(),
                })?;
        }
        Ok(graph)
    }

    fn apply_loop_momentum_basis(
        &self,
        graph: &mut Graph,
        requested_edges: &[usize],
    ) -> Result<(), FeynkitAdapterError> {
        for &edge in requested_edges {
            if edge >= graph.underlying.n_edges() {
                return Err(FeynkitAdapterError::UnknownLoopMomentumEdge {
                    graph: graph.name.clone(),
                    edge,
                });
            }
        }

        let full_filter = graph.full_filter();
        let mut cut_graph = full_filter.subtract(&graph.initial_state_cut.right);
        for &edge in requested_edges {
            cut_graph.sub(graph[&EdgeIndex(edge)].1);
        }
        let Some(root_edge) = cut_graph.included_iter().next() else {
            return Err(FeynkitAdapterError::EmptyLoopMomentumTree {
                graph: graph.name.clone(),
            });
        };
        let tree = SimpleTraversalTree::depth_first_traverse(
            graph,
            &cut_graph,
            &graph.node_id(root_edge),
            None,
        )
        .map_err(|source| FeynkitAdapterError::LoopMomentumTreeTraversal {
            graph: graph.name.clone(),
            source,
        })?;
        let external = graph.internal_crown(&full_filter);
        let mut basis = graph
            .lmb_impl(&full_filter, &tree.tree_subgraph, external)
            .map_err(|source| FeynkitAdapterError::LoopMomentumBasis {
                graph: graph.name.clone(),
                source,
            })?;

        for edge in graph.initial_state_cut.left.included_iter() {
            let edge = graph[&edge];
            let Some((position, _)) = basis
                .loop_edges
                .iter()
                .find_position(|candidate| **candidate == edge)
            else {
                return Err(FeynkitAdapterError::InitialStateLoopEdge {
                    graph: graph.name.clone(),
                    edge: edge.0,
                });
            };
            basis.put_loop_to_ext(LoopIndex(position));
        }
        graph.canonicalize_lmb_external_order(&mut basis);

        let mut requested_positions = Vec::with_capacity(requested_edges.len());
        for &edge in requested_edges {
            let Some((position, _)) = basis
                .loop_edges
                .iter()
                .find_position(|candidate| candidate.0 == edge)
            else {
                return Err(FeynkitAdapterError::RequestedLoopEdge {
                    graph: graph.name.clone(),
                    edge,
                });
            };
            requested_positions.push(LoopIndex(position));
        }
        for (left, right) in Permutation::sort(&requested_positions)
            .inverse()
            .transpositions()
        {
            basis.swap_loops(LoopIndex(left), LoopIndex(right));
        }

        graph.loop_momentum_basis = basis;
        graph.param_builder = ParamBuilder::new::<_, Vec<Atom>, _>(
            graph,
            self.model,
            &graph.loop_momentum_basis,
            Vec::new(),
        );
        Ok(())
    }
}
