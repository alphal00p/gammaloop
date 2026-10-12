use crate::{
    ExternalLabel, GenerationDiagramView, GenerationDressingView, GenerationError,
    GenerationSkeletonView, GenerationType, ParticleSelector, Process, VertexSelector,
};
use feynkit_graph::{DiagramId, EdgeId};
use feynkit_model::{Model, ParticleId, VertexRuleId};
pub use linnet::half_edge::generate::Bridges;
use linnet::half_edge::generate::batch;
use num_bigint::BigUint;
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
    ops::{BitOr, ControlFlow, RangeInclusive},
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use symbolica::{
    atom::{Atom, AtomCore},
    parser::ParseSettings,
};
#[derive(
    Debug,
    Clone,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
pub struct GraphGroupingOptions {
    pub numerical_sample_seed: u16,
    pub number_of_numerical_samples: usize,
    pub differentiate_particle_masses_only: bool,
    pub fully_numerical_substitution: bool,
    pub check_canonical_numerator: bool,
    /// Reuse external wavefunction samples for the two sides of a sewn state.
    pub symmetric_polarizations: bool,
}

impl Default for GraphGroupingOptions {
    fn default() -> Self {
        Self {
            numerical_sample_seed: 3,
            number_of_numerical_samples: 5,
            differentiate_particle_masses_only: true,
            fully_numerical_substitution: false,
            check_canonical_numerator: false,
            symmetric_polarizations: false,
        }
    }
}

#[derive(
    Debug,
    Clone,
    Default,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
#[serde(rename_all = "snake_case")]
pub enum Grouping {
    #[default]
    None,
    Identical(GraphGroupingOptions),
    UpToSign(GraphGroupingOptions),
    UpToScalar(GraphGroupingOptions),
}
bincode::impl_borrow_decode!(Grouping);

/// An explicit unrestricted or inclusive interval. A missing upper bound means infinity.
#[derive(
    Debug,
    Clone,
    Copy,
    Default,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
#[serde(rename_all = "snake_case")]
pub enum Range {
    #[default]
    Unrestricted,
    Inclusive {
        min: usize,
        max: Option<usize>,
    },
}
impl Range {
    pub const fn exact(value: usize) -> Self {
        Self::Inclusive {
            min: value,
            max: Some(value),
        }
    }
    pub const fn at_least(min: usize) -> Self {
        Self::Inclusive { min, max: None }
    }
    pub const fn at_most(max: usize) -> Self {
        Self::Inclusive {
            min: 0,
            max: Some(max),
        }
    }
    pub fn contains(&self, value: usize) -> bool {
        let (min, max) = self.bounds();
        value >= min && max.is_none_or(|max| value <= max)
    }
    pub fn bounds(&self) -> (usize, Option<usize>) {
        match *self {
            Self::Unrestricted => (0, None),
            Self::Inclusive { min, max } => (min, max),
        }
    }
    pub fn is_empty(&self) -> bool {
        let (min, max) = self.bounds();
        max.is_some_and(|max| min > max)
    }
    pub fn intersect(self, other: Self) -> Self {
        let (a, b) = self.bounds();
        let (c, d) = other.bounds();
        Self::Inclusive {
            min: a.max(c),
            max: match (b, d) {
                (Some(x), Some(y)) => Some(x.min(y)),
                (x, None) | (None, x) => x,
            },
        }
    }
}
impl From<RangeInclusive<usize>> for Range {
    fn from(v: RangeInclusive<usize>) -> Self {
        Self::Inclusive {
            min: *v.start(),
            max: Some(*v.end()),
        }
    }
}

/// Total amplitude loops or loops of the sewn forward graph. Every request supplies one.
#[derive(
    Debug,
    Clone,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
#[serde(rename_all = "snake_case")]
pub enum Loops {
    Total(RangeInclusive<usize>),
    Forward(RangeInclusive<usize>),
}
impl Loops {
    pub fn total(range: RangeInclusive<usize>) -> Self {
        Self::Total(range)
    }
    pub fn forward(range: RangeInclusive<usize>) -> Self {
        Self::Forward(range)
    }
    pub fn range(&self) -> RangeInclusive<usize> {
        match self {
            Self::Total(r) | Self::Forward(r) => r.clone(),
        }
    }
}

/// Plain bounds constrain one graph; squared bounds constrain an allowed pair of amplitudes.
#[derive(
    Debug,
    Clone,
    Default,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
pub struct Orders {
    pub plain: BTreeMap<String, Range>,
    pub squared: BTreeMap<String, Range>,
    pub perturbative: BTreeMap<String, usize>,
}
impl Orders {
    pub fn new() -> Self {
        Self::default()
    }
    fn add(mut self, name: impl Into<String>, range: Range, squared: bool) -> Self {
        let table = if squared {
            &mut self.squared
        } else {
            &mut self.plain
        };
        table
            .entry(name.into())
            .and_modify(|old| *old = old.intersect(range))
            .or_insert(range);
        self
    }
    pub fn max(self, name: impl Into<String>, max: usize) -> Self {
        self.add(name, Range::at_most(max), false)
    }
    pub fn min(self, name: impl Into<String>, min: usize) -> Self {
        self.add(name, Range::at_least(min), false)
    }
    pub fn exact(self, name: impl Into<String>, value: usize) -> Self {
        self.add(name, Range::exact(value), false)
    }
    pub fn squared_max(self, name: impl Into<String>, max: usize) -> Self {
        self.add(name, Range::at_most(max), true)
    }
    pub fn squared_min(self, name: impl Into<String>, min: usize) -> Self {
        self.add(name, Range::at_least(min), true)
    }
    pub fn squared_exact(self, name: impl Into<String>, value: usize) -> Self {
        self.add(name, Range::exact(value), true)
    }
    pub fn perturbative(mut self, name: impl Into<String>, power: usize) -> Self {
        self.perturbative.insert(name.into(), power);
        self
    }
    pub fn merge(mut self, other: Self) -> Self {
        for (n, r) in other.plain {
            self = self.add(n, r, false)
        }
        for (n, r) in other.squared {
            self = self.add(n, r, true)
        }
        for (n, p) in other.perturbative {
            self.perturbative
                .entry(n)
                .and_modify(|old| *old = (*old).min(p))
                .or_insert(p);
        }
        self
    }
    pub fn is_empty(&self) -> bool {
        self.plain
            .values()
            .chain(self.squared.values())
            .any(Range::is_empty)
    }
}

/// Allow clauses intersect after model resolution; exclusions unite. An empty clause allows none.
#[derive(
    Debug,
    Clone,
    Default,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
pub struct Restrict {
    pub particle_allow: Vec<BTreeSet<ParticleSelector>>,
    pub particle_exclude: BTreeSet<ParticleSelector>,
    pub vertex_allow: Vec<BTreeSet<VertexSelector>>,
    pub vertex_exclude: BTreeSet<VertexSelector>,
}
impl Restrict {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn only_particles<I, P>(mut self, items: I) -> Self
    where
        I: IntoIterator<Item = P>,
        P: Into<ParticleSelector>,
    {
        self.particle_allow
            .push(items.into_iter().map(Into::into).collect());
        self
    }
    pub fn exclude_particles<I, P>(mut self, items: I) -> Self
    where
        I: IntoIterator<Item = P>,
        P: Into<ParticleSelector>,
    {
        self.particle_exclude
            .extend(items.into_iter().map(Into::into));
        self
    }
    pub fn only_vertices<I, V>(mut self, items: I) -> Self
    where
        I: IntoIterator<Item = V>,
        V: Into<VertexSelector>,
    {
        self.vertex_allow
            .push(items.into_iter().map(Into::into).collect());
        self
    }
    pub fn exclude_vertices<I, V>(mut self, items: I) -> Self
    where
        I: IntoIterator<Item = V>,
        V: Into<VertexSelector>,
    {
        self.vertex_exclude
            .extend(items.into_iter().map(Into::into));
        self
    }
    pub fn merge(mut self, other: Self) -> Self {
        self.particle_allow.extend(other.particle_allow);
        self.particle_exclude.extend(other.particle_exclude);
        self.vertex_allow.extend(other.vertex_allow);
        self.vertex_exclude.extend(other.vertex_exclude);
        self
    }
    pub(crate) fn resolve(&self, model: &Model) -> Result<ResolvedRestrictions, GenerationError> {
        let species = |selectors: &BTreeSet<ParticleSelector>| -> Result<BTreeSet<ParticleId>, GenerationError> {
            let mut ids = BTreeSet::new();
            for selector in selectors {
                let id = selector.resolve(model)?;
                ids.insert(id);
                ids.insert(model.particle_by_id(id)?.antiparticle);
            }
            Ok(ids)
        };
        let mut particles = None;
        for clause in &self.particle_allow {
            let set = species(clause)?;
            particles = Some(match particles {
                None => set,
                Some(old) => set.intersection(&old).copied().collect(),
            });
        }
        let mut vertices = None;
        for clause in &self.vertex_allow {
            let set = clause
                .iter()
                .map(|s| s.resolve(model))
                .collect::<Result<BTreeSet<_>, _>>()?;
            vertices = Some(match vertices {
                None => set,
                Some(old) => set.intersection(&old).copied().collect(),
            });
        }
        Ok(ResolvedRestrictions {
            particles,
            particle_veto: species(&self.particle_exclude)?,
            vertices,
            vertex_veto: self
                .vertex_exclude
                .iter()
                .map(|s| s.resolve(model))
                .collect::<Result<_, _>>()?,
        })
    }
}
#[derive(Debug, Clone)]
pub(crate) struct ResolvedRestrictions {
    pub particles: Option<BTreeSet<ParticleId>>,
    pub particle_veto: BTreeSet<ParticleId>,
    pub vertices: Option<BTreeSet<VertexRuleId>>,
    pub vertex_veto: BTreeSet<VertexRuleId>,
}
impl ResolvedRestrictions {
    pub fn particle(&self, p: ParticleId) -> bool {
        self.particles.as_ref().is_none_or(|set| set.contains(&p))
            && !self.particle_veto.contains(&p)
    }
    pub fn vertex(&self, v: VertexRuleId) -> bool {
        self.vertices.as_ref().is_none_or(|set| set.contains(&v)) && !self.vertex_veto.contains(&v)
    }
}

/// Mass qualifiers for an opt-in topology veto. No veto is the shared default.
#[derive(
    Debug,
    Clone,
    Copy,
    Default,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
pub struct Veto {
    pub massive: bool,
    pub massless: bool,
}
impl Veto {
    pub const fn all() -> Self {
        Self {
            massive: true,
            massless: true,
        }
    }
    pub const fn massive() -> Self {
        Self {
            massive: true,
            massless: false,
        }
    }
    pub const fn massless() -> Self {
        Self {
            massive: false,
            massless: true,
        }
    }
    pub fn enabled(self) -> bool {
        self.massive || self.massless
    }
}
#[derive(
    Debug,
    Clone,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
pub struct TopologyRules {
    pub bridges: Bridges,
    pub self_loops: bool,
    pub zero_flow_edges: bool,
    pub max_vertices: Option<usize>,
    pub self_energies: Veto,
    pub tadpoles: Veto,
    pub snails: Veto,
    pub fermion_loops: Range,
    pub factorized_loops: Range,
}
impl Default for TopologyRules {
    fn default() -> Self {
        Self {
            bridges: Bridges::Any,
            self_loops: true,
            zero_flow_edges: false,
            max_vertices: None,
            self_energies: Veto::default(),
            tadpoles: Veto::default(),
            snails: Veto::default(),
            fermion_loops: Range::Unrestricted,
            factorized_loops: Range::Unrestricted,
        }
    }
}
impl TopologyRules {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn with_bridges(mut self, v: Bridges) -> Self {
        self.bridges = v;
        self
    }
    pub fn with_self_loops(mut self, v: bool) -> Self {
        self.self_loops = v;
        self
    }
    pub fn with_zero_flow_edges(mut self, v: bool) -> Self {
        self.zero_flow_edges = v;
        self
    }
    pub fn with_max_vertices(mut self, v: usize) -> Self {
        self.max_vertices = Some(v);
        self
    }
    pub fn with_self_energies(mut self, v: Veto) -> Self {
        self.self_energies = v;
        self
    }
    pub fn with_tadpoles(mut self, v: Veto) -> Self {
        self.tadpoles = v;
        self
    }
    pub fn with_snails(mut self, v: Veto) -> Self {
        self.snails = v;
        self
    }
    pub fn with_fermion_loops(mut self, v: impl Into<Range>) -> Self {
        self.fermion_loops = v.into();
        self
    }
    pub fn with_factorized_loops(mut self, v: impl Into<Range>) -> Self {
        self.factorized_loops = v.into();
        self
    }
}
#[derive(
    Debug,
    Clone,
    Default,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
pub struct CutAmplitude {
    pub orders: Orders,
    pub loops: Range,
}
#[derive(
    Debug,
    Clone,
    Default,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
pub struct Cuts {
    pub amplitude: CutAmplitude,
    pub blobs: Range,
    pub spectators: Range,
    pub sewn_tadpoles: bool,
    pub symmetrize_left_right: bool,
}
impl Cuts {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn with_cut_amplitude(mut self, v: CutAmplitude) -> Self {
        self.amplitude = v;
        self
    }
    pub fn with_blobs(mut self, v: impl Into<Range>) -> Self {
        self.blobs = v.into();
        self
    }
    pub fn with_spectators(mut self, v: impl Into<Range>) -> Self {
        self.spectators = v.into();
        self
    }
    pub fn with_sewn_tadpoles(mut self, v: bool) -> Self {
        self.sewn_tadpoles = v;
        self
    }
    pub fn with_symmetrize_left_right(mut self, v: bool) -> Self {
        self.symmetrize_left_right = v;
        self
    }
}
#[derive(
    Debug,
    Clone,
    Copy,
    Default,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
pub struct Symmetrize {
    pub initial: bool,
    pub final_state: bool,
    pub external_fermions: bool,
}
impl Symmetrize {
    pub const INITIAL: Self = Self {
        initial: true,
        final_state: false,
        external_fermions: false,
    };
    pub const FINAL: Self = Self {
        initial: false,
        final_state: true,
        external_fermions: false,
    };
    pub const EXTERNAL_FERMIONS: Self = Self {
        initial: false,
        final_state: false,
        external_fermions: true,
    };
}
impl BitOr for Symmetrize {
    type Output = Self;
    fn bitor(self, rhs: Self) -> Self {
        Self {
            initial: self.initial || rhs.initial,
            final_state: self.final_state || rhs.final_state,
            external_fermions: self.external_fermions || rhs.external_fermions,
        }
    }
}

type SkeletonCallback = Arc<dyn Fn(&GenerationSkeletonView<'_>) -> bool + Send + Sync + 'static>;
type DressingCallback = Arc<dyn Fn(&GenerationDressingView<'_>) -> bool + Send + Sync + 'static>;
type SelectionCallback = Arc<dyn Fn(&GenerationDiagramView<'_>) -> bool + Send + Sync + 'static>;
type DeferredDispatch = Arc<
    dyn Fn(Vec<batch::Request<ParticleId, ExternalLabel>>) -> ControlFlow<(), Vec<batch::Reply>>
        + Send
        + Sync
        + 'static,
>;
/// Persisted physical request. Callable hooks are deliberately omitted from serialization.
#[derive(Clone, Serialize, Deserialize)]
pub struct Request {
    pub(crate) process: Process,
    pub(crate) kind: GenerationType,
    pub(crate) loops: Loops,
    pub(crate) restrictions: Restrict,
    pub(crate) orders: Orders,
    pub(crate) topology: TopologyRules,
    pub(crate) cuts: Cuts,
    pub(crate) symmetrization: Symmetrize,
    #[serde(skip)]
    pub(crate) prune_skeletons: Option<SkeletonCallback>,
    #[serde(skip)]
    pub(crate) prune_dressings: Option<DressingCallback>,
    #[serde(skip)]
    pub(crate) select: Option<SelectionCallback>,
    #[serde(skip)]
    pub(crate) deferred_hooks: Option<(batch::HookMask, DeferredDispatch)>,
}
impl Request {
    fn new(process: Process, kind: GenerationType, loops: Loops) -> Self {
        Self {
            process,
            kind,
            loops,
            restrictions: Restrict::new(),
            orders: Orders::new(),
            topology: TopologyRules::default(),
            cuts: Cuts::new(),
            symmetrization: Symmetrize::default(),
            prune_skeletons: None,
            prune_dressings: None,
            select: None,
            deferred_hooks: None,
        }
    }
    pub fn amplitude(process: Process, loops: Loops) -> Self {
        Self::new(process, GenerationType::Amplitude, loops)
    }
    pub fn cross_section(process: Process, loops: Loops) -> Self {
        Self::new(process, GenerationType::CrossSection, loops)
    }
    pub fn process(&self) -> &Process {
        &self.process
    }
    pub fn kind(&self) -> GenerationType {
        self.kind
    }
    pub fn loops(&self) -> &Loops {
        &self.loops
    }
    pub fn orders(&self) -> &Orders {
        &self.orders
    }
    pub fn restrictions(&self) -> &Restrict {
        &self.restrictions
    }
    pub fn topology(&self) -> &TopologyRules {
        &self.topology
    }
    pub fn cuts(&self) -> &Cuts {
        &self.cuts
    }
    pub fn symmetrization(&self) -> Symmetrize {
        self.symmetrization
    }
    pub fn with_orders(mut self, v: Orders) -> Self {
        self.orders = self.orders.merge(v);
        self
    }
    pub fn with_restrict(mut self, v: Restrict) -> Self {
        self.restrictions = self.restrictions.merge(v);
        self
    }
    pub fn with_topology(mut self, v: TopologyRules) -> Self {
        self.topology = v;
        self
    }
    pub fn with_cuts(mut self, v: Cuts) -> Self {
        self.cuts = v;
        self
    }
    pub fn with_symmetrize(mut self, v: Symmetrize) -> Self {
        self.symmetrization = v;
        self
    }
    pub fn prune_skeletons(
        mut self,
        v: impl Fn(&GenerationSkeletonView<'_>) -> bool + Send + Sync + 'static,
    ) -> Self {
        self.prune_skeletons = Some(Arc::new(v));
        self
    }
    pub fn prune_dressings(
        mut self,
        v: impl Fn(&GenerationDressingView<'_>) -> bool + Send + Sync + 'static,
    ) -> Self {
        self.prune_dressings = Some(Arc::new(v));
        self
    }
    pub fn select(
        mut self,
        v: impl Fn(&GenerationDiagramView<'_>) -> bool + Send + Sync + 'static,
    ) -> Self {
        self.select = Some(Arc::new(v));
        self
    }
    #[doc(hidden)]
    pub fn deferred_hooks(
        mut self,
        mask: batch::HookMask,
        dispatch: impl Fn(
            Vec<batch::Request<ParticleId, ExternalLabel>>,
        ) -> ControlFlow<(), Vec<batch::Reply>>
        + Send
        + Sync
        + 'static,
    ) -> Self {
        self.deferred_hooks = Some((mask, Arc::new(dispatch)));
        self
    }
    /// Validate persisted semantics and model selectors without compiling a search plan.
    pub fn validate_in(&self, model: &Model) -> Result<(), GenerationError> {
        self.validate()?;
        self.process.validate_in(model)?;
        if self.kind == GenerationType::Amplitude && !self.orders.perturbative.is_empty() {
            return Err(GenerationError::InvalidRequest(
                "perturbative unresolved-cut orders require a cross section".into(),
            ));
        }
        self.restrictions.resolve(model)?;
        for name in self
            .orders
            .plain
            .keys()
            .chain(self.orders.squared.keys())
            .chain(self.orders.perturbative.keys())
            .chain(self.cuts.amplitude.orders.plain.keys())
        {
            model.order(name)?;
        }
        if !self.cuts.amplitude.orders.squared.is_empty()
            || !self.cuts.amplitude.orders.perturbative.is_empty()
        {
            return Err(GenerationError::InvalidRequest(
                "cut-amplitude orders accept plain bounds only".into(),
            ));
        }
        if self.kind == GenerationType::Amplitude && self.cuts != Cuts::default() {
            return Err(GenerationError::InvalidRequest(
                "physical cut options require a cross section".into(),
            ));
        }
        self.process
            .validate_covariant_cut_restrictions(model, &self.restrictions, self.kind)
    }
    pub(crate) fn validate(&self) -> Result<(), GenerationError> {
        self.process.validate()?;
        let correct = matches!(
            (self.kind, &self.loops),
            (GenerationType::Amplitude, Loops::Total(_))
                | (GenerationType::CrossSection, Loops::Forward(_))
        );
        if !correct {
            return Err(GenerationError::InvalidRequest(
                "loop domain must be total for amplitudes and forward for cross sections".into(),
            ));
        }
        if self.loops.range().is_empty() {
            return Err(GenerationError::InvalidRequest(
                "loop range is empty".into(),
            ));
        }
        if self.kind == GenerationType::CrossSection && self.symmetrization.final_state {
            return Err(GenerationError::CrossSectionFinalStateSymmetrization);
        }
        if self.kind == GenerationType::Amplitude && self.process.outgoing_alternatives().len() != 1
        {
            return Err(crate::ProcessError::MultipleAmplitudeFinalStates.into());
        }
        Ok(())
    }
}
impl fmt::Debug for Request {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Request")
            .field("process", &self.process)
            .field("kind", &self.kind)
            .field("loops", &self.loops)
            .field("restrictions", &self.restrictions)
            .field("orders", &self.orders)
            .field("topology", &self.topology)
            .field("cuts", &self.cuts)
            .field("symmetrization", &self.symmetrization)
            .field("prune_skeletons", &self.prune_skeletons.is_some())
            .field("prune_dressings", &self.prune_dressings.is_some())
            .field("select", &self.select.is_some())
            .finish()
    }
}
impl PartialEq for Request {
    fn eq(&self, other: &Self) -> bool {
        self.process == other.process
            && self.kind == other.kind
            && self.loops == other.loops
            && self.restrictions == other.restrictions
            && self.orders == other.orders
            && self.topology == other.topology
            && self.cuts == other.cuts
            && self.symmetrization == other.symmetrization
    }
}
impl Eq for Request {}
// Native field encoding deliberately excludes callable handles, including deferred transport.
impl bincode::Encode for Request {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        e: &mut E,
    ) -> Result<(), bincode::error::EncodeError> {
        bincode::Encode::encode(
            &(
                &self.process,
                &self.kind,
                &self.loops,
                &self.restrictions,
                &self.orders,
                &self.topology,
                &self.cuts,
                &self.symmetrization,
            ),
            e,
        )
    }
}
impl<C> bincode::Decode<C> for Request {
    fn decode<D: bincode::de::Decoder<Context = C>>(
        d: &mut D,
    ) -> Result<Self, bincode::error::DecodeError> {
        let (process, kind, loops, restrictions, orders, topology, cuts, symmetrization) =
            bincode::Decode::decode(d)?;
        Ok(Self {
            process,
            kind,
            loops,
            restrictions,
            orders,
            topology,
            cuts,
            symmetrization,
            prune_skeletons: None,
            prune_dressings: None,
            select: None,
            deferred_hooks: None,
        })
    }
}
bincode::impl_borrow_decode!(Request);

/// Explicit naming demand; iteration and compact generation never introduce names.
#[derive(
    Debug,
    Clone,
    Default,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
#[serde(rename_all = "snake_case")]
pub enum Naming {
    #[default]
    None,
    Named {
        prefix: String,
    },
}
impl Naming {
    pub fn prefix(prefix: impl Into<String>) -> Self {
        Self::Named {
            prefix: prefix.into(),
        }
    }
    pub fn named() -> Self {
        Self::prefix("GL")
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NumeratorOptions {
    pub prefactor: Atom,
    pub projector: Option<Atom>,
    pub naming: Naming,
    pub filter_zero_color: bool,
    pub loop_momentum_bases: BTreeMap<DiagramId, Vec<EdgeId>>,
    pub named_loop_momentum_bases: BTreeMap<String, Vec<EdgeId>>,
}
impl Default for NumeratorOptions {
    fn default() -> Self {
        Self {
            prefactor: Atom::one(),
            projector: None,
            naming: Naming::None,
            filter_zero_color: false,
            loop_momentum_bases: BTreeMap::new(),
            named_loop_momentum_bases: BTreeMap::new(),
        }
    }
}
impl NumeratorOptions {
    pub fn with_prefactor(mut self, v: Atom) -> Self {
        self.prefactor = v;
        self
    }
    pub fn with_projector(mut self, v: Atom) -> Self {
        self.projector = Some(v);
        self
    }
    pub fn with_naming(mut self, v: Naming) -> Self {
        self.naming = v;
        self
    }
    pub fn with_zero_color_filter(mut self, v: bool) -> Self {
        self.filter_zero_color = v;
        self
    }
    pub fn with_loop_momentum_basis(
        mut self,
        id: DiagramId,
        edges: impl IntoIterator<Item = EdgeId>,
    ) -> Self {
        self.loop_momentum_bases
            .insert(id, edges.into_iter().collect());
        self
    }
    pub fn with_named_loop_momentum_basis(
        mut self,
        name: impl Into<String>,
        edges: impl IntoIterator<Item = EdgeId>,
    ) -> Self {
        self.named_loop_momentum_bases
            .insert(name.into(), edges.into_iter().collect());
        self
    }
}
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Selection {
    pub selected_diagram_ids: Option<BTreeSet<DiagramId>>,
    pub selected_diagram_names: Option<BTreeSet<String>>,
    pub vetoed_diagram_ids: BTreeSet<DiagramId>,
    pub vetoed_diagram_names: BTreeSet<String>,
    pub forced_cuts: Vec<BTreeSet<EdgeId>>,
}
impl Selection {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn select_names<S: Into<String>>(mut self, v: impl IntoIterator<Item = S>) -> Self {
        self.selected_diagram_names = Some(v.into_iter().map(Into::into).collect());
        self
    }
    pub fn select_ids(mut self, v: impl IntoIterator<Item = DiagramId>) -> Self {
        self.selected_diagram_ids = Some(v.into_iter().collect());
        self
    }
    pub fn veto_names<S: Into<String>>(mut self, v: impl IntoIterator<Item = S>) -> Self {
        self.vetoed_diagram_names
            .extend(v.into_iter().map(Into::into));
        self
    }
    pub fn veto_ids(mut self, v: impl IntoIterator<Item = DiagramId>) -> Self {
        self.vetoed_diagram_ids.extend(v);
        self
    }
    pub fn with_forced_cuts<I, C>(mut self, v: I) -> Self
    where
        I: IntoIterator<Item = C>,
        C: IntoIterator<Item = EdgeId>,
    {
        self.forced_cuts = v.into_iter().map(|cut| cut.into_iter().collect()).collect();
        self
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Stage {
    Skeletons,
    Dressing,
    Prepare,
    Numerators,
    Selection,
    ZeroColor,
    GroupingPreparation,
    GroupingSamples,
    GroupingComparison,
    Grouping,
    Done,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenerationProgress {
    pub stage: Stage,
    pub completed: BigUint,
    pub total: Option<BigUint>,
    pub zero_numerator_count: BigUint,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenerationControl {
    Continue,
    Cancel,
}
#[derive(Debug, Clone, Default)]
pub struct CancellationToken(Arc<AtomicBool>);
impl CancellationToken {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn cancel(&self) {
        self.0.store(true, Ordering::Release)
    }
    pub fn is_cancelled(&self) -> bool {
        self.0.load(Ordering::Acquire)
    }
}
type ProgressCallback =
    Arc<dyn Fn(GenerationProgress) -> GenerationControl + Send + Sync + 'static>;
type CancellationCheck = Arc<dyn Fn() -> bool + Send + Sync + 'static>;
/// Execution state is never serialized. Clones share cancellation, not progress snapshots.
#[derive(Clone)]
pub struct Run {
    pub(crate) threads: Option<usize>,
    pub(crate) parallel: bool,
    pub(crate) check_heredity: bool,
    pub(crate) cancellation: CancellationToken,
    pub(crate) cancellation_check: Option<CancellationCheck>,
    pub(crate) progress: Option<ProgressCallback>,
}
impl Default for Run {
    fn default() -> Self {
        Self {
            threads: None,
            parallel: true,
            check_heredity: false,
            cancellation: CancellationToken::new(),
            cancellation_check: None,
            progress: None,
        }
    }
}
impl fmt::Debug for Run {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Run")
            .field("threads", &self.threads)
            .field("parallel", &self.parallel)
            .field("check_heredity", &self.check_heredity)
            .field("cancelled", &self.cancellation.is_cancelled())
            .finish()
    }
}
impl Run {
    pub fn threads(mut self, v: usize) -> Self {
        self.threads = Some(v);
        self
    }
    pub fn parallel(mut self, v: bool) -> Self {
        self.parallel = v;
        self
    }
    pub fn check_heredity(mut self, v: bool) -> Self {
        self.check_heredity = v;
        self
    }
    pub fn cancellation(&self) -> &CancellationToken {
        &self.cancellation
    }
    pub fn cancellation_token(mut self, v: CancellationToken) -> Self {
        self.cancellation = v;
        self
    }
    pub fn cancellation_check(mut self, v: impl Fn() -> bool + Send + Sync + 'static) -> Self {
        self.cancellation_check = Some(Arc::new(v));
        self
    }
    pub fn progress(
        mut self,
        v: impl Fn(GenerationProgress) -> GenerationControl + Send + Sync + 'static,
    ) -> Self {
        self.progress = Some(Arc::new(v));
        self
    }
    pub fn cancellation_requested(&self) -> bool {
        if !self.cancellation.is_cancelled()
            && self.cancellation_check.as_ref().is_some_and(|v| v())
        {
            self.cancellation.cancel()
        }
        self.cancellation.is_cancelled()
    }
    pub(crate) fn report_progress(
        &self,
        stage: Stage,
        completed: impl Into<BigUint>,
        total: Option<usize>,
    ) {
        self.report_progress_with_zero_count(stage, completed, total, 0usize)
    }
    pub(crate) fn report_progress_with_zero_count(
        &self,
        stage: Stage,
        completed: impl Into<BigUint>,
        total: Option<usize>,
        zero: impl Into<BigUint>,
    ) {
        if self.progress.as_ref().is_some_and(|cb| {
            cb(GenerationProgress {
                stage,
                completed: completed.into(),
                total: total.map(Into::into),
                zero_numerator_count: zero.into(),
            }) == GenerationControl::Cancel
        }) {
            self.cancellation.cancel()
        }
    }
}

impl bincode::Encode for Selection {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        e: &mut E,
    ) -> Result<(), bincode::error::EncodeError> {
        let selected = self
            .selected_diagram_ids
            .as_ref()
            .map(|s| s.iter().map(|i| i.0).collect::<BTreeSet<_>>());
        let vetoed = self
            .vetoed_diagram_ids
            .iter()
            .map(|i| i.0)
            .collect::<BTreeSet<_>>();
        let cuts = self
            .forced_cuts
            .iter()
            .map(|s| s.iter().map(|i| i.0).collect::<BTreeSet<_>>())
            .collect::<Vec<_>>();
        bincode::Encode::encode(
            &(
                selected,
                &self.selected_diagram_names,
                vetoed,
                &self.vetoed_diagram_names,
                cuts,
            ),
            e,
        )
    }
}
impl<C> bincode::Decode<C> for Selection {
    fn decode<D: bincode::de::Decoder<Context = C>>(
        d: &mut D,
    ) -> Result<Self, bincode::error::DecodeError> {
        let (selected, selected_diagram_names, vetoed, vetoed_diagram_names, cuts): (
            Option<BTreeSet<u128>>,
            _,
            BTreeSet<u128>,
            _,
            Vec<BTreeSet<usize>>,
        ) = bincode::Decode::decode(d)?;
        Ok(Self {
            selected_diagram_ids: selected.map(|s| s.into_iter().map(DiagramId).collect()),
            selected_diagram_names,
            vetoed_diagram_ids: vetoed.into_iter().map(DiagramId).collect(),
            vetoed_diagram_names,
            forced_cuts: cuts
                .into_iter()
                .map(|s| s.into_iter().map(EdgeId).collect())
                .collect(),
        })
    }
}
bincode::impl_borrow_decode!(Selection);

#[derive(Serialize, Deserialize, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
struct NumeratorData {
    prefactor: String,
    projector: Option<String>,
    naming: Naming,
    filter_zero_color: bool,
    loop_momentum_bases: BTreeMap<u128, Vec<usize>>,
    named_loop_momentum_bases: BTreeMap<String, Vec<usize>>,
}
impl From<&NumeratorOptions> for NumeratorData {
    fn from(v: &NumeratorOptions) -> Self {
        Self {
            prefactor: v.prefactor.to_canonical_string(),
            projector: v.projector.as_ref().map(AtomCore::to_canonical_string),
            naming: v.naming.clone(),
            filter_zero_color: v.filter_zero_color,
            loop_momentum_bases: v
                .loop_momentum_bases
                .iter()
                .map(|(id, edges)| (id.0, edges.iter().map(|e| e.0).collect()))
                .collect(),
            named_loop_momentum_bases: v
                .named_loop_momentum_bases
                .iter()
                .map(|(name, edges)| (name.clone(), edges.iter().map(|e| e.0).collect()))
                .collect(),
        }
    }
}
impl TryFrom<NumeratorData> for NumeratorOptions {
    type Error = String;
    fn try_from(v: NumeratorData) -> Result<Self, Self::Error> {
        let parse = |text: &str| {
            Atom::parse(text, "feynkit_generation_options", ParseSettings::default())
                .map_err(|e| e.to_string())
        };
        Ok(Self {
            prefactor: parse(&v.prefactor)?,
            projector: v.projector.as_deref().map(parse).transpose()?,
            naming: v.naming,
            filter_zero_color: v.filter_zero_color,
            loop_momentum_bases: v
                .loop_momentum_bases
                .into_iter()
                .map(|(id, edges)| (DiagramId(id), edges.into_iter().map(EdgeId).collect()))
                .collect(),
            named_loop_momentum_bases: v
                .named_loop_momentum_bases
                .into_iter()
                .map(|(name, edges)| (name, edges.into_iter().map(EdgeId).collect()))
                .collect(),
        })
    }
}
impl Serialize for NumeratorOptions {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        NumeratorData::from(self).serialize(s)
    }
}
impl<'de> Deserialize<'de> for NumeratorOptions {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        NumeratorData::deserialize(d)?
            .try_into()
            .map_err(serde::de::Error::custom)
    }
}
impl bincode::Encode for NumeratorOptions {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        e: &mut E,
    ) -> Result<(), bincode::error::EncodeError> {
        bincode::Encode::encode(&NumeratorData::from(self), e)
    }
}
impl<C> bincode::Decode<C> for NumeratorOptions {
    fn decode<D: bincode::de::Decoder<Context = C>>(
        d: &mut D,
    ) -> Result<Self, bincode::error::DecodeError> {
        NumeratorData::decode(d)?
            .try_into()
            .map_err(bincode::error::DecodeError::OtherString)
    }
}
bincode::impl_borrow_decode!(NumeratorOptions);

#[cfg(test)]
mod tests {
    use super::*;
    fn model() -> Model {
        Model::from_json(include_str!(
            "../../feynkit-model/tests/fixtures/scalars_2p_3p.json"
        ))
        .unwrap()
    }
    #[test]
    fn grouped_bounds_and_restrictions_intersect_after_alias_resolution() {
        let model = model();
        let id = model.particle_id("scalar_0").unwrap();
        let restriction = Restrict::new()
            .only_particles(["scalar_0"])
            .merge(Restrict::new().only_particles([ParticleSelector::by_id(&model, id).unwrap()]));
        assert!(restriction.resolve(&model).unwrap().particle(id));
        assert!(
            !restriction
                .merge(Restrict::new().exclude_particles([1000i64]))
                .resolve(&model)
                .unwrap()
                .particle(id)
        );
        assert!(
            !Restrict::new()
                .only_particles(std::iter::empty::<&str>())
                .resolve(&model)
                .unwrap()
                .particle(id)
        );
        let orders = Orders::new()
            .min("QCD", 3)
            .merge(Orders::new().max("QCD", 1));
        assert!(orders.plain["QCD"].is_empty());
        assert!(Range::Unrestricted.contains(usize::MAX));
    }
    #[test]
    fn persisted_request_and_build_options_drop_execution_hooks_only() {
        let request = Process::new(["scalar_0"], ["scalar_0"])
            .amplitude(Loops::total(0..=2))
            .with_restrict(Restrict::new().only_vertices(["V_3_SCALAR_000"]))
            .with_topology(TopologyRules::default().with_bridges(Bridges::OnePIInLoops))
            .prune_skeletons(|_| true)
            .select(|_| false);
        let decoded: Request =
            serde_json::from_str(&serde_json::to_string(&request).unwrap()).unwrap();
        assert_eq!(decoded, request);
        assert!(decoded.prune_skeletons.is_none() && decoded.select.is_none());
        let encoded = bincode::encode_to_vec(&request, bincode::config::standard()).unwrap();
        let (decoded, n): (Request, usize) =
            bincode::decode_from_slice(&encoded, bincode::config::standard()).unwrap();
        assert_eq!(n, encoded.len());
        assert_eq!(decoded, request);
        assert!(decoded.prune_skeletons.is_none());
        let atom = Atom::parse("(x+y)^3", "P7_wire", ParseSettings::default()).unwrap();
        let options = NumeratorOptions::default()
            .with_prefactor(atom.clone())
            .with_projector(atom)
            .with_naming(Naming::named())
            .with_loop_momentum_basis(DiagramId(u128::MAX), [EdgeId(9)]);
        let bytes = bincode::encode_to_vec(&options, bincode::config::standard()).unwrap();
        let (decoded, n): (NumeratorOptions, usize) =
            bincode::decode_from_slice(&bytes, bincode::config::standard()).unwrap();
        assert_eq!(n, bytes.len());
        assert_eq!(decoded, options);
        let decoded: NumeratorOptions =
            serde_json::from_str(&serde_json::to_string(&options).unwrap()).unwrap();
        assert_eq!(decoded, options);
    }
    #[test]
    fn domains_are_explicit_and_perturbative_orders_remain_cross_section_only() {
        let model = model();
        let process = Process::new(["scalar_0"], ["scalar_0"]);
        assert!(
            process
                .clone()
                .amplitude(Loops::forward(1..=1))
                .validate_in(&model)
                .is_err()
        );
        assert!(
            process
                .clone()
                .cross_section(Loops::total(1..=1))
                .validate_in(&model)
                .is_err()
        );
        let request = process
            .amplitude(Loops::total(0..=0))
            .with_orders(Orders::new().perturbative("QCD", 2));
        assert!(request.validate_in(&model).is_err());
    }
}
