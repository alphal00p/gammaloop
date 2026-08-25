use std::{
    collections::BTreeMap,
    fmt,
    ops::RangeInclusive,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};

use serde::{Deserialize, Serialize};

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
    ParticleVeto(Vec<i64>),
    VertexAllow(Vec<String>),
    VertexVeto(Vec<String>),
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
    pub(crate) cancellation: CancellationToken,
    pub(crate) cancellation_check: Option<CancellationCheck>,
    pub(crate) progress: Option<ProgressCallback>,
}

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
            cancellation: CancellationToken::new(),
            cancellation_check: None,
            progress: None,
        }
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
