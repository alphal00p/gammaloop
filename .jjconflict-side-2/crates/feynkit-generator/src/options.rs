use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
    ops::RangeInclusive,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};

use feynkit_graph::{DiagramId, EdgeId};
use feynkit_model::Model;
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore},
    parser::ParseSettings,
};

use crate::{ParticleSelector, SelectorError, VertexSelector};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GraphGroupingOptions {
    pub numerical_sample_seed: u16,
    pub number_of_numerical_samples: usize,
    pub differentiate_particle_masses_only: bool,
    pub fully_numerical_substitution: bool,
    pub check_canonical_numerator: bool,
}

impl Default for GraphGroupingOptions {
    fn default() -> Self {
        Self {
            numerical_sample_seed: 3,
            number_of_numerical_samples: 5,
            differentiate_particle_masses_only: true,
            fully_numerical_substitution: false,
            check_canonical_numerator: false,
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NumeratorGrouping {
    #[default]
    None,
    OnlyDetectZeroes,
    Identical(GraphGroupingOptions),
    UpToSign(GraphGroupingOptions),
    UpToScalar(GraphGroupingOptions),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct SelfEnergyFilterOptions {
    pub veto_massive: bool,
    pub veto_massless: bool,
    pub only_scaleless: bool,
}

impl Default for SelfEnergyFilterOptions {
    fn default() -> Self {
        Self {
            veto_massive: true,
            veto_massless: true,
            only_scaleless: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct SnailFilterOptions {
    pub veto_attached_to_massive: bool,
    pub veto_attached_to_massless: bool,
    pub only_scaleless: bool,
}

impl Default for SnailFilterOptions {
    fn default() -> Self {
        Self {
            veto_attached_to_massive: false,
            veto_attached_to_massless: true,
            only_scaleless: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct TadpoleFilterOptions {
    pub veto_attached_to_massive: bool,
    pub veto_attached_to_massless: bool,
    pub only_scaleless: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct SewnFilterOptions {
    pub filter_tadpoles: bool,
}

impl Default for SewnFilterOptions {
    fn default() -> Self {
        Self {
            filter_tadpoles: true,
        }
    }
}

impl Default for TadpoleFilterOptions {
    fn default() -> Self {
        Self {
            veto_attached_to_massive: true,
            veto_attached_to_massless: true,
            only_scaleless: false,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GenerationFilter {
    SelfEnergy(SelfEnergyFilterOptions),
    Tadpoles(TadpoleFilterOptions),
    ZeroSnails(SnailFilterOptions),
    Sewn(SewnFilterOptions),
    ParticleVeto(Vec<ParticleSelector>),
    VertexAllow(Vec<VertexSelector>),
    VertexVeto(Vec<VertexSelector>),
    MaxNumberOfBridges(usize),
    CouplingOrders(BTreeMap<String, (usize, Option<usize>)>),
    LoopCountRange((usize, usize)),
    BlobRange(RangeInclusive<usize>),
    SpectatorRange(RangeInclusive<usize>),
    PerturbativeOrders(BTreeMap<String, usize>),
    FermionLoopCountRange((usize, usize)),
    FactorizedLoopTopologiesCountRange((usize, usize)),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FilterScope {
    Graph,
    CutAmplitude,
}

impl fmt::Display for FilterScope {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Graph => "graph",
            Self::CutAmplitude => "cut amplitude",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GenerationFilterKind {
    SelfEnergy,
    Tadpoles,
    ZeroSnails,
    Sewn,
    ParticleVeto,
    VertexAllow,
    VertexVeto,
    MaxNumberOfBridges,
    CouplingOrders,
    LoopCountRange,
    BlobRange,
    SpectatorRange,
    PerturbativeOrders,
    FermionLoopCountRange,
    FactorizedLoopTopologiesCountRange,
}

impl fmt::Display for GenerationFilterKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::SelfEnergy => "self_energy",
            Self::Tadpoles => "tadpoles",
            Self::ZeroSnails => "zero_snails",
            Self::Sewn => "sewn",
            Self::ParticleVeto => "particle_veto",
            Self::VertexAllow => "vertex_allow",
            Self::VertexVeto => "vertex_veto",
            Self::MaxNumberOfBridges => "max_number_of_bridges",
            Self::CouplingOrders => "coupling_orders",
            Self::LoopCountRange => "loop_count_range",
            Self::BlobRange => "blob_range",
            Self::SpectatorRange => "spectator_range",
            Self::PerturbativeOrders => "perturbative_orders",
            Self::FermionLoopCountRange => "fermion_loop_count_range",
            Self::FactorizedLoopTopologiesCountRange => "factorized_loop_topologies_count_range",
        })
    }
}

impl GenerationFilter {
    pub const fn kind(&self) -> GenerationFilterKind {
        match self {
            Self::SelfEnergy(_) => GenerationFilterKind::SelfEnergy,
            Self::Tadpoles(_) => GenerationFilterKind::Tadpoles,
            Self::ZeroSnails(_) => GenerationFilterKind::ZeroSnails,
            Self::Sewn(_) => GenerationFilterKind::Sewn,
            Self::ParticleVeto(_) => GenerationFilterKind::ParticleVeto,
            Self::VertexAllow(_) => GenerationFilterKind::VertexAllow,
            Self::VertexVeto(_) => GenerationFilterKind::VertexVeto,
            Self::MaxNumberOfBridges(_) => GenerationFilterKind::MaxNumberOfBridges,
            Self::CouplingOrders(_) => GenerationFilterKind::CouplingOrders,
            Self::LoopCountRange(_) => GenerationFilterKind::LoopCountRange,
            Self::BlobRange(_) => GenerationFilterKind::BlobRange,
            Self::SpectatorRange(_) => GenerationFilterKind::SpectatorRange,
            Self::PerturbativeOrders(_) => GenerationFilterKind::PerturbativeOrders,
            Self::FermionLoopCountRange(_) => GenerationFilterKind::FermionLoopCountRange,
            Self::FactorizedLoopTopologiesCountRange(_) => {
                GenerationFilterKind::FactorizedLoopTopologiesCountRange
            }
        }
    }
}

impl fmt::Display for GenerationFilter {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GenerationProgress {
    pub generated_graphs: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenerationControl {
    Continue,
    Cancel,
}

/// Cloneable cancellation state that can be shared with another thread.
#[derive(Debug, Clone, Default)]
pub struct CancellationToken(Arc<AtomicBool>);

impl CancellationToken {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn cancel(&self) {
        self.0.store(true, Ordering::Release);
    }

    pub fn is_cancelled(&self) -> bool {
        self.0.load(Ordering::Acquire)
    }
}

type ProgressCallback =
    Arc<dyn Fn(GenerationProgress) -> GenerationControl + Send + Sync + 'static>;
type CancellationCheck = Arc<dyn Fn() -> bool + Send + Sync + 'static>;

#[derive(Clone)]
pub struct GenerationOptions {
    pub(crate) threads: Option<usize>,
    pub(crate) max_vertices: Option<usize>,
    pub(crate) allow_self_loops: bool,
    pub(crate) allow_zero_flow_edges: bool,
    pub(crate) graph_filters: Vec<GenerationFilter>,
    pub(crate) cut_amplitude_filters: Vec<GenerationFilter>,
    pub(crate) numerator_grouping: NumeratorGrouping,
    pub(crate) graph_prefix: String,
    pub(crate) selected_diagram_ids: Option<BTreeSet<DiagramId>>,
    pub(crate) selected_diagram_names: Option<BTreeSet<String>>,
    pub(crate) vetoed_diagram_ids: BTreeSet<DiagramId>,
    pub(crate) vetoed_diagram_names: BTreeSet<String>,
    pub(crate) loop_momentum_bases: BTreeMap<DiagramId, Vec<EdgeId>>,
    pub(crate) named_loop_momentum_bases: BTreeMap<String, Vec<EdgeId>>,
    pub(crate) forced_cuts: Vec<BTreeSet<EdgeId>>,
    pub(crate) numerator_prefactor: Atom,
    pub(crate) projector: Option<Atom>,
    pub(crate) cancellation: CancellationToken,
    pub(crate) cancellation_check: Option<CancellationCheck>,
    pub(crate) progress: Option<ProgressCallback>,
}

/// Stable persistent portion of [`GenerationOptions`]. Callbacks and
/// cancellation state are deliberately runtime-only and are restored to their
/// defaults after deserialization.
#[derive(Serialize, Deserialize, PartialEq, Eq)]
struct GenerationOptionsSerde {
    threads: Option<usize>,
    max_vertices: Option<usize>,
    allow_self_loops: bool,
    allow_zero_flow_edges: bool,
    graph_filters: Vec<GenerationFilter>,
    cut_amplitude_filters: Vec<GenerationFilter>,
    numerator_grouping: NumeratorGrouping,
    graph_prefix: String,
    selected_diagram_ids: Option<BTreeSet<DiagramId>>,
    selected_diagram_names: Option<BTreeSet<String>>,
    vetoed_diagram_ids: BTreeSet<DiagramId>,
    vetoed_diagram_names: BTreeSet<String>,
    loop_momentum_bases: BTreeMap<DiagramId, Vec<EdgeId>>,
    named_loop_momentum_bases: BTreeMap<String, Vec<EdgeId>>,
    forced_cuts: Vec<BTreeSet<EdgeId>>,
    numerator_prefactor: String,
    projector: Option<String>,
}

impl From<&GenerationOptions> for GenerationOptionsSerde {
    fn from(options: &GenerationOptions) -> Self {
        Self {
            threads: options.threads,
            max_vertices: options.max_vertices,
            allow_self_loops: options.allow_self_loops,
            allow_zero_flow_edges: options.allow_zero_flow_edges,
            graph_filters: options.graph_filters.clone(),
            cut_amplitude_filters: options.cut_amplitude_filters.clone(),
            numerator_grouping: options.numerator_grouping.clone(),
            graph_prefix: options.graph_prefix.clone(),
            selected_diagram_ids: options.selected_diagram_ids.clone(),
            selected_diagram_names: options.selected_diagram_names.clone(),
            vetoed_diagram_ids: options.vetoed_diagram_ids.clone(),
            vetoed_diagram_names: options.vetoed_diagram_names.clone(),
            loop_momentum_bases: options.loop_momentum_bases.clone(),
            named_loop_momentum_bases: options.named_loop_momentum_bases.clone(),
            forced_cuts: options.forced_cuts.clone(),
            numerator_prefactor: options.numerator_prefactor.to_canonical_string(),
            projector: options
                .projector
                .as_ref()
                .map(AtomCore::to_canonical_string),
        }
    }
}

impl Serialize for GenerationOptions {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        GenerationOptionsSerde::from(self).serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for GenerationOptions {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let persistent = GenerationOptionsSerde::deserialize(deserializer)?;
        let parse = |expression: String, field: &str| {
            Atom::parse(
                &expression,
                "feynkit_generation_options",
                ParseSettings::default(),
            )
            .map_err(|message| {
                serde::de::Error::custom(format!(
                    "failed to parse {field} '{expression}': {message}"
                ))
            })
        };
        Ok(Self {
            threads: persistent.threads,
            max_vertices: persistent.max_vertices,
            allow_self_loops: persistent.allow_self_loops,
            allow_zero_flow_edges: persistent.allow_zero_flow_edges,
            graph_filters: persistent.graph_filters,
            cut_amplitude_filters: persistent.cut_amplitude_filters,
            numerator_grouping: persistent.numerator_grouping,
            graph_prefix: persistent.graph_prefix,
            selected_diagram_ids: persistent.selected_diagram_ids,
            selected_diagram_names: persistent.selected_diagram_names,
            vetoed_diagram_ids: persistent.vetoed_diagram_ids,
            vetoed_diagram_names: persistent.vetoed_diagram_names,
            loop_momentum_bases: persistent.loop_momentum_bases,
            named_loop_momentum_bases: persistent.named_loop_momentum_bases,
            forced_cuts: persistent.forced_cuts,
            numerator_prefactor: parse(persistent.numerator_prefactor, "numerator prefactor")?,
            projector: persistent
                .projector
                .map(|projector| parse(projector, "projector"))
                .transpose()?,
            cancellation: CancellationToken::new(),
            cancellation_check: None,
            progress: None,
        })
    }
}

impl PartialEq for GenerationOptions {
    fn eq(&self, other: &Self) -> bool {
        GenerationOptionsSerde::from(self) == GenerationOptionsSerde::from(other)
    }
}

impl Eq for GenerationOptions {}

impl Default for GenerationOptions {
    fn default() -> Self {
        Self {
            threads: None,
            max_vertices: None,
            allow_self_loops: false,
            allow_zero_flow_edges: false,
            graph_filters: Vec::new(),
            cut_amplitude_filters: Vec::new(),
            numerator_grouping: NumeratorGrouping::None,
            graph_prefix: "FK".to_owned(),
            selected_diagram_ids: None,
            selected_diagram_names: None,
            vetoed_diagram_ids: BTreeSet::new(),
            vetoed_diagram_names: BTreeSet::new(),
            loop_momentum_bases: BTreeMap::new(),
            named_loop_momentum_bases: BTreeMap::new(),
            forced_cuts: Vec::new(),
            numerator_prefactor: Atom::one(),
            projector: None,
            cancellation: CancellationToken::new(),
            cancellation_check: None,
            progress: None,
        }
    }
}

impl GenerationOptions {
    /// Resolve all boundary selectors to stable model IDs exactly once.
    pub(crate) fn resolve_selectors(&self, model: &Model) -> Result<Self, SelectorError> {
        fn resolve_filter(
            filter: &GenerationFilter,
            model: &Model,
        ) -> Result<GenerationFilter, SelectorError> {
            Ok(match filter {
                GenerationFilter::ParticleVeto(selectors) => GenerationFilter::ParticleVeto(
                    selectors
                        .iter()
                        .map(|selector| ParticleSelector::by_id(model, selector.resolve(model)?))
                        .collect::<Result<_, _>>()?,
                ),
                GenerationFilter::VertexAllow(selectors) => GenerationFilter::VertexAllow(
                    selectors
                        .iter()
                        .map(|selector| VertexSelector::by_id(model, selector.resolve(model)?))
                        .collect::<Result<_, _>>()?,
                ),
                GenerationFilter::VertexVeto(selectors) => GenerationFilter::VertexVeto(
                    selectors
                        .iter()
                        .map(|selector| VertexSelector::by_id(model, selector.resolve(model)?))
                        .collect::<Result<_, _>>()?,
                ),
                other => other.clone(),
            })
        }

        let mut resolved = self.clone();
        resolved.graph_filters = self
            .graph_filters
            .iter()
            .map(|filter| resolve_filter(filter, model))
            .collect::<Result<_, _>>()?;
        resolved.cut_amplitude_filters = self
            .cut_amplitude_filters
            .iter()
            .map(|filter| resolve_filter(filter, model))
            .collect::<Result<_, _>>()?;
        Ok(resolved)
    }
}

impl fmt::Debug for GenerationOptions {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("GenerationOptions")
            .field("threads", &self.threads)
            .field("max_vertices", &self.max_vertices)
            .field("allow_self_loops", &self.allow_self_loops)
            .field("allow_zero_flow_edges", &self.allow_zero_flow_edges)
            .field("graph_filters", &self.graph_filters)
            .field("cut_amplitude_filters", &self.cut_amplitude_filters)
            .field("numerator_grouping", &self.numerator_grouping)
            .field("graph_prefix", &self.graph_prefix)
            .field("selected_diagram_ids", &self.selected_diagram_ids)
            .field("selected_diagram_names", &self.selected_diagram_names)
            .field("vetoed_diagram_ids", &self.vetoed_diagram_ids)
            .field("vetoed_diagram_names", &self.vetoed_diagram_names)
            .field("loop_momentum_bases", &self.loop_momentum_bases)
            .field("named_loop_momentum_bases", &self.named_loop_momentum_bases)
            .field("forced_cuts", &self.forced_cuts)
            .field("numerator_prefactor", &self.numerator_prefactor)
            .field("projector", &self.projector)
            .field("cancelled", &self.cancellation.is_cancelled())
            .field("cancellation_check", &self.cancellation_check.is_some())
            .field("progress", &self.progress.is_some())
            .finish()
    }
}

impl GenerationOptions {
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = Some(threads);
        self
    }

    pub fn max_vertices(mut self, max_vertices: usize) -> Self {
        self.max_vertices = Some(max_vertices);
        self
    }

    pub fn allow_self_loops(mut self, allow: bool) -> Self {
        self.allow_self_loops = allow;
        self
    }

    pub fn allow_zero_flow_edges(mut self, allow: bool) -> Self {
        self.allow_zero_flow_edges = allow;
        self
    }

    pub fn with_graph_filter(mut self, filter: GenerationFilter) -> Self {
        self.graph_filters.push(filter);
        self
    }

    pub fn with_cut_amplitude_filter(mut self, filter: GenerationFilter) -> Self {
        self.cut_amplitude_filters.push(filter);
        self
    }

    pub fn with_scoped_filter(self, scope: FilterScope, filter: GenerationFilter) -> Self {
        match scope {
            FilterScope::Graph => self.with_graph_filter(filter),
            FilterScope::CutAmplitude => self.with_cut_amplitude_filter(filter),
        }
    }

    pub fn filters(&self, scope: FilterScope) -> &[GenerationFilter] {
        match scope {
            FilterScope::Graph => &self.graph_filters,
            FilterScope::CutAmplitude => &self.cut_amplitude_filters,
        }
    }

    pub fn numerator_grouping(mut self, grouping: NumeratorGrouping) -> Self {
        self.numerator_grouping = grouping;
        self
    }

    pub fn graph_prefix(mut self, prefix: impl Into<String>) -> Self {
        self.graph_prefix = prefix.into();
        self
    }

    /// Retain only diagrams with one of the finalized deterministic names.
    pub fn select_diagrams<I, S>(mut self, names: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.selected_diagram_names = Some(names.into_iter().map(Into::into).collect());
        self
    }

    /// Retain only diagrams with one of the stable content-derived IDs.
    pub fn select_diagram_ids<I>(mut self, ids: I) -> Self
    where
        I: IntoIterator<Item = DiagramId>,
    {
        self.selected_diagram_ids = Some(ids.into_iter().collect());
        self
    }

    /// Remove diagrams with one of the finalized deterministic names.
    pub fn veto_diagrams<I, S>(mut self, names: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.vetoed_diagram_names = names.into_iter().map(Into::into).collect();
        self
    }

    /// Remove diagrams with one of the stable content-derived IDs.
    pub fn veto_diagram_ids<I>(mut self, ids: I) -> Self
    where
        I: IntoIterator<Item = DiagramId>,
    {
        self.vetoed_diagram_ids = ids.into_iter().collect();
        self
    }

    /// Select the ordered independent edges for one finalized diagram.
    pub fn with_loop_momentum_basis<I>(mut self, diagram: DiagramId, edges: I) -> Self
    where
        I: IntoIterator<Item = EdgeId>,
    {
        self.loop_momentum_bases
            .insert(diagram, edges.into_iter().collect());
        self
    }

    /// Select a routing using the display name convenience key.
    pub fn with_named_loop_momentum_basis<I>(mut self, diagram: impl Into<String>, edges: I) -> Self
    where
        I: IntoIterator<Item = EdgeId>,
    {
        self.named_loop_momentum_bases
            .insert(diagram.into(), edges.into_iter().collect());
        self
    }

    /// Retain only the listed physical cut-edge sets before diagram grouping.
    pub fn forced_cuts<I, C>(mut self, cuts: I) -> Self
    where
        I: IntoIterator<Item = C>,
        C: IntoIterator<Item = EdgeId>,
    {
        self.forced_cuts = cuts
            .into_iter()
            .map(|cut| cut.into_iter().collect())
            .collect();
        self
    }

    /// Multiply every generated numerator by an application-independent factor.
    pub fn numerator_prefactor(mut self, prefactor: Atom) -> Self {
        self.numerator_prefactor = prefactor;
        self
    }

    /// Override the automatically generated external-state projector.
    ///
    /// Passing `1` explicitly disables external wavefunctions; leaving the
    /// option unset generates the physical wavefunction for every external leg.
    pub fn projector(mut self, projector: Atom) -> Self {
        self.projector = Some(projector);
        self
    }

    pub fn cancellation_token(mut self, token: CancellationToken) -> Self {
        self.cancellation = token;
        self
    }

    /// Poll an application-owned cancellation source during every generation stage.
    pub fn cancellation_check(mut self, check: impl Fn() -> bool + Send + Sync + 'static) -> Self {
        self.cancellation_check = Some(Arc::new(check));
        self
    }

    pub fn progress(
        mut self,
        callback: impl Fn(GenerationProgress) -> GenerationControl + Send + Sync + 'static,
    ) -> Self {
        self.progress = Some(Arc::new(callback));
        self
    }

    pub(crate) fn max_bridges(&self) -> Option<usize> {
        self.graph_filters.iter().find_map(|filter| match filter {
            GenerationFilter::MaxNumberOfBridges(maximum) => Some(*maximum),
            _ => None,
        })
    }

    pub(crate) fn cancellation_requested(&self) -> bool {
        self.cancellation.is_cancelled()
            || self
                .cancellation_check
                .as_ref()
                .is_some_and(|check| check())
    }
}

#[cfg(test)]
mod tests {
    use symbolica::atom::Atom;

    use super::GenerationOptions;

    #[test]
    fn projector_serde_distinguishes_automatic_from_explicit_one() {
        let automatic = GenerationOptions::default();
        let automatic_json = serde_json::to_string(&automatic).unwrap();
        assert!(automatic_json.contains(r#""projector":null"#));
        let automatic_roundtrip: GenerationOptions = serde_json::from_str(&automatic_json).unwrap();
        assert!(automatic_roundtrip.projector.is_none());

        let explicit = GenerationOptions::default().projector(Atom::one());
        let explicit_json = serde_json::to_string(&explicit).unwrap();
        assert!(explicit_json.contains(r#""projector":"1""#));
        let explicit_roundtrip: GenerationOptions = serde_json::from_str(&explicit_json).unwrap();
        assert_eq!(explicit_roundtrip.projector, Some(Atom::one()));
    }
}
