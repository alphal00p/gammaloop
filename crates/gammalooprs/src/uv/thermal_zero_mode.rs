//! Hosts of the finite-temperature Matsubara zero-mode subtraction.
//!
//! At finite temperature, a set Σ of massless bosonic loop lines without a
//! chemical potential can carry their static Matsubara mode while the rest
//! of the graph stays hard. With L_Σ soft loops, Σ has the static infrared
//! degree `ω(Σ) = 3 L_Σ + Σ_e dod(e) + Σ_v dod(v)`
//! ([`UltravioletGraph::static_ir_degree`]), and `ω(Σ) ≤ 0` makes that
//! region divergent. The complementary spinney γ hosts its subtraction: the
//! zero-mode operation `S0_γ` projects the Matsubara zero mode of the
//! frequency classes of Σ and Taylor-expands the gapped hard part γ in its
//! soft momenta, to the relative order `J = -ω(Σ)` for the static mode.
//! Tree lines carry no loop momentum and belong to no soft set.
//!
//! [`assign_zero_mode_hosts`] selects the hosts of one cut subgraph and
//! records their [`ThermalZeroModeRegion`]s on the classified spinneys:
//! - a connected `IR` spinney hosts the region on its own wood node, whose
//!   operation is `S0 + U - U S0`;
//! - a disjoint union of `IR` spinneys hosts the joint region of both
//!   factors on a copy of the union, a separate wood node above their
//!   factorized join, whose operation is `S0` alone;
//! - every other `IR` spinney becomes a plain `MUV` U.
//!
//! The pass runs only where completed thermal soft branches are frozen
//! ([`UVgenerationSettings::freezes_thermal_soft_branches`]), independently of
//! vacuum subtraction, which only adds the vacuum-soft order of a class.
//!
//! The host policy of a cut subgraph runs in numbered passes:
//! 1. candidates: connected `IR` spinneys of positive degree, and disjoint
//!    unions of such spinneys;
//! 2. the soft set Σ of each candidate, which must be eligible;
//! 3. the frequency classes of Σ;
//! 4. the static infrared degree `ω(Σ)`;
//! 5. the order tuples the candidate needs;
//! 6. Taylor safety and zero-bin coverage in the candidate's canonical chart,
//!    inner candidates first: the accepted hosts;
//! 7. pruning of the inner hosts an outer host subsumes;
//! 8. the nesting of the hosts and their single uncut structure;
//! 9. a host for every divergent soft set;
//! 10. installing the regions of the hosts;
//! 11. demoting every other `IR` spinney to `MUV`.
//!
//! Passes 12 and 13 restrict the full observable under vacuum subtraction
//! and check that a host keeps its enclosing whole-graph U
//! ([`Graph::restrict_full_observable_spinney`]).

use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::mem;

use itertools::Itertools;
use linnet::half_edge::involution::EdgeIndex;
use linnet::half_edge::subgraph::{Inclusion, InternalSubGraph, SuBitGraph, SubSetLike, SubSetOps};
use thiserror::Error;
use three_dimensional_reps::utils::extend_independent_rows;

use crate::debug_tags;
use crate::graph::{Graph, LMBext, LmbError, LoopMomentumBasis};
use crate::settings::global::MediumSettings;
use crate::uv::approx::CutStructure;
use crate::uv::uv_graph::UVE;
use crate::uv::{ApproximationType, Spinney, UVgenerationSettings, UltravioletGraph};

/// Also the expected regions and policy helpers of the hedge-poset wood tests.
#[cfg(test)]
#[path = "thermal_zero_mode_tests.rs"]
pub(in crate::uv) mod tests;

/// How a zero-mode host enters the hedge-poset wood.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) enum ZeroModeHostKind {
    /// A connected spinney. Its own node subtracts the zero mode and the UV
    /// of the hard part: `S0 + U - U S0`.
    Connected,
    /// A disjoint union of divergent spinneys, expanded jointly in the soft
    /// momenta. A copy of the union above its factorized join applies `S0`
    /// alone; the U of each factor precedes it in the forest.
    Union,
}

/// The distributions a zero-mode coefficient is read from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) enum ZeroModeSector {
    /// The sector as handed to the operation. It serves the static order
    /// `2 nu` of every class.
    Complete,
    /// Every distribution outside the projected classes at its vacuum value.
    /// Under vacuum subtraction it serves the vacuum 4D-soft order `2 nu - 1`
    /// of a single class: that coefficient keeps its sum here, while its CFF
    /// variants stay bounded where those of the complete sector cancel only
    /// in the sum. No other order is read from this sector. A summand that a
    /// class has already projected at this order is read off unchanged, so it
    /// must come from an earlier hard-vacuum projection; one that keeps the
    /// thermal distribution of a line outside the classes is an internal
    /// error of the projector.
    HardVacuum,
}

/// A frequency class of a soft set: raised edges that carry one Matsubara
/// frequency at one on-shell energy. Its static propagator power `nu` is its
/// number of edges.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct ZeroModeClassSpec {
    edges: BTreeSet<EdgeIndex>,
    /// The smallest edge.
    representative: EdgeIndex,
    power: u32,
}

impl ZeroModeClassSpec {
    /// The class of the raised `edges`, or `None` when it has no edge.
    pub(crate) fn new(edges: BTreeSet<EdgeIndex>) -> Option<Self> {
        let representative = *edges.first()?;
        let power = u32::try_from(edges.len())
            .expect("a frequency class of one graph has fewer than 2^32 edges");
        Some(Self {
            edges,
            representative,
            power,
        })
    }

    /// The static propagator power `nu`.
    pub(crate) fn power(&self) -> u32 {
        self.power
    }

    /// The smallest edge, whose energy stands for the class.
    pub(crate) fn representative(&self) -> EdgeIndex {
        self.representative
    }

    pub(crate) fn edges(&self) -> &BTreeSet<EdgeIndex> {
        &self.edges
    }
}

impl Display for ZeroModeClassSpec {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "{{{}}}^{}", self.edges.iter().join(","), self.power)
    }
}

/// One Laurent order of a zero-mode operation and the Taylor order of the
/// hard part that accompanies it.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct ZeroModeTuple {
    /// The order `a_c` of the coefficient `prod_c x_c^-a_c`, per class.
    orders: Vec<u32>,
    sector: ZeroModeSector,
    /// The relative Taylor order `J = -ω(Σ) - Σ_c (2 nu_c - a_c)`.
    taylor_order: i64,
}

impl ZeroModeTuple {
    pub(crate) fn orders(&self) -> &[u32] {
        &self.orders
    }

    pub(crate) fn sector(&self) -> ZeroModeSector {
        self.sector
    }

    pub(crate) fn taylor_order(&self) -> i64 {
        self.taylor_order
    }

    /// The tuples a host needs, for classes of the static powers `powers`, or
    /// `None` when the static mode converges (`J = -ω < 0`).
    ///
    /// The static tuple `a = 2 nu` reads the complete sector at `J = -ω`.
    /// Under vacuum subtraction a single class also has the vacuum-soft pole
    /// `2 nu - 1`, read from the hard-vacuum sector, at `J = -ω - 1`. A tuple
    /// with static and vacuum orders vanishes: the thermal part is even in
    /// each class energy, and the vacuum part has the pole in every class at
    /// once.
    ///
    /// # Errors
    ///
    /// A needed order `a ≤ 2 nu - 2`, and the joint vacuum pole of several
    /// classes, which the projection cannot supply exactly.
    fn required(
        powers: &[u32],
        omega: i32,
        vacuum_subtraction: bool,
    ) -> Result<Option<Vec<Self>>, UnsupportedOrders> {
        // J(a) = -ω - Σ_c (2 nu_c - a_c).
        let taylor_order = |orders: &[u32]| {
            powers
                .iter()
                .zip(orders)
                .fold(-i64::from(omega), |order, (power, laurent)| {
                    order - 2 * i64::from(*power) + i64::from(*laurent)
                })
        };
        let static_orders = powers.iter().map(|power| 2 * power).collect_vec();
        let static_order = taylor_order(&static_orders);
        if static_order < 0 {
            return Ok(None);
        }
        // The highest lower order: one class two below its static order.
        if let Some(class) = powers.iter().position(|power| *power >= 2) {
            let mut orders = static_orders.clone();
            orders[class] -= 2;
            let order = taylor_order(&orders);
            if order >= 0 {
                return Err(UnsupportedOrders::BelowVacuumPole {
                    orders,
                    taylor_order: order,
                });
            }
        }
        let mut tuples = vec![Self {
            orders: static_orders.clone(),
            sector: ZeroModeSector::Complete,
            taylor_order: static_order,
        }];
        if vacuum_subtraction {
            let orders = static_orders.iter().map(|order| order - 1).collect_vec();
            let order = taylor_order(&orders);
            if order >= 0 {
                if powers.len() > 1 {
                    return Err(UnsupportedOrders::MultiClassVacuumPole {
                        orders,
                        taylor_order: order,
                    });
                }
                tuples.push(Self {
                    orders,
                    sector: ZeroModeSector::HardVacuum,
                    taylor_order: order,
                });
            }
        }
        Ok(Some(tuples))
    }
}

impl Display for ZeroModeTuple {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(
            f,
            "({}) J{} {:?}",
            self.orders.iter().join(","),
            self.taylor_order,
            self.sector
        )
    }
}

/// The finite-temperature zero-mode region a spinney hosts.
///
/// The operation `S0` of the host projects every class of [`Self::classes`]
/// onto the orders of [`Self::tuples`], in the canonical coordinate basis of
/// [`Self::coordinate_loop_edges`]. Lines of Σ outside the classes are
/// static conditions implied by the classes: they keep their energies and
/// get no zero-mode temperature.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct ThermalZeroModeRegion {
    kind: ZeroModeHostKind,
    sigma: BTreeSet<EdgeIndex>,
    classes: Vec<ZeroModeClassSpec>,
    omega: i32,
    soft_loops: usize,
    tuples: Vec<ZeroModeTuple>,
    zero_bins: Vec<SuBitGraph>,
    coordinate_loop_edges: Vec<EdgeIndex>,
}

impl ThermalZeroModeRegion {
    pub(crate) fn kind(&self) -> ZeroModeHostKind {
        self.kind
    }

    /// The soft set Σ: the loop lines of the cut subgraph outside the host.
    pub(crate) fn sigma(&self) -> &BTreeSet<EdgeIndex> {
        &self.sigma
    }

    /// The Taylor order `-ω` of the static mode.
    fn static_taylor_order(&self) -> i64 {
        -i64::from(self.omega)
    }

    /// The projected frequency classes, by decreasing power, then
    /// representative.
    pub(crate) fn classes(&self) -> &[ZeroModeClassSpec] {
        &self.classes
    }

    /// The requested orders: the static order first.
    pub(crate) fn tuples(&self) -> &[ZeroModeTuple] {
        &self.tuples
    }

    /// The loop edges of the host's canonical coordinate basis
    /// ([`Graph::canonical_component_lmb`]), the one chart of every operation
    /// of the host.
    pub(crate) fn coordinate_loop_edges(&self) -> &[EdgeIndex] {
        &self.coordinate_loop_edges
    }

    /// The first frequency class of which the integrated lines `contracted`
    /// hold some edges but not all, if any. A wholly contracted class leaves
    /// the zero mode no class energy, so the zero-mode operation vanishes
    /// exactly there; a partly contracted one would leave its Laurent
    /// projection incomplete.
    pub(crate) fn partly_contracted_class(
        &self,
        contracted: &BTreeSet<EdgeIndex>,
    ) -> Option<&ZeroModeClassSpec> {
        self.classes.iter().find(|class| {
            let held = class.edges().intersection(contracted).count();
            held != 0 && held != class.edges().len()
        })
    }
}

impl Display for ThermalZeroModeRegion {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(
            f,
            "{:?} host, Σ {{{}}}, classes [{}], ω {}, L {}, tuples [{}], zero bins [{}], LMB [{}]",
            self.kind,
            self.sigma.iter().join(","),
            self.classes.iter().join(", "),
            self.omega,
            self.soft_loops,
            self.tuples.iter().join(", "),
            self.zero_bins
                .iter()
                .map(SuBitGraph::string_label)
                .join(", "),
            self.coordinate_loop_edges.iter().join(","),
        )
    }
}

/// Why zero-mode hosts could not be assigned.
#[derive(Debug, Error)]
pub(crate) enum ZeroModeHostError {
    #[error(
        "the soft set {{{}}} (ω = {omega}) is infrared divergent at finite temperature, but no \
         zero-mode host subtracts it; finite-temperature soft subtraction needs \
         `renormalization_prescription.massless_power_divergent = IR` and a supported \
         zero-mode host",
        .sigma.iter().join(",")
    )]
    Hostless {
        sigma: BTreeSet<EdgeIndex>,
        omega: i32,
    },
    #[error(
        "zero-mode hosts {{{}}} and {{{}}} are not nested",
        .first.iter().join(","),
        .second.iter().join(",")
    )]
    IncompatibleHosts {
        first: BTreeSet<EdgeIndex>,
        second: BTreeSet<EdgeIndex>,
    },
    #[error(
        "zero-mode host {{{}}} contains host {{{}}}; only a union host may contain another \
         host, and only as one of its components",
        .outer.iter().join(","),
        .inner.iter().join(",")
    )]
    UnsupportedNesting {
        outer: BTreeSet<EdgeIndex>,
        inner: BTreeSet<EdgeIndex>,
    },
    #[error(
        "zero-mode host {{{}}} needs orders the zero-mode projection cannot supply exactly",
        .host.iter().join(",")
    )]
    UnsupportedOrders {
        host: BTreeSet<EdgeIndex>,
        #[source]
        reason: UnsupportedOrders,
    },
    #[error(
        "zero-mode host {{{}}} lies on a cut structure with nonempty or several cuts, whose \
         spinneys are merged across cuts; zero-mode hosts need a single uncut structure",
        .host.iter().join(",")
    )]
    NonTrivialCut { host: BTreeSet<EdgeIndex> },
    #[error(
        "zero-mode host {{{}}} needs the enclosing whole-graph U (its static cograph is \
         log-divergent at the top Taylor order), but the full graph was not retained \
         (dod {dod})",
        .host.iter().join(",")
    )]
    MissingEnclosingSpinney { host: BTreeSet<EdgeIndex>, dod: i32 },
    #[error(
        "incompatible zero-mode hosts in a factorized join: the factors {} each wait for \
         another one, because a connected host must replay before every factor its soft set \
         meets",
        .factors
            .iter()
            .map(|factor| format!("{{{}}}", factor.iter().join(",")))
            .join(", ")
    )]
    MutuallyCoveringJoin { factors: Vec<BTreeSet<EdgeIndex>> },
    /// A bug: the components of a union host are classified spinneys, so its
    /// plain twin is a factorized join.
    #[error(
        "internal error: the zero-mode union host {{{}}} has the scope of a collective plain \
         union, whose trace operations cannot be told apart from its own",
        .host.iter().join(",")
    )]
    CollectiveUnionScope { host: BTreeSet<EdgeIndex> },
    /// A bug: an integrated prefix reaches the zero-mode operation of a host
    /// either inside the host, which leaves its soft set alone, or in a
    /// replay that keeps the host's plain U alone.
    #[error(
        "internal error: the integrated prefix {{{}}} contracts part of the frequency class {} \
         of the zero-mode host {{{}}}",
        .prefix.iter().join(","),
        .class,
        .host.iter().join(",")
    )]
    PartlyContractedClass {
        host: BTreeSet<EdgeIndex>,
        class: ZeroModeClassSpec,
        prefix: BTreeSet<EdgeIndex>,
    },
    /// A bug: the integrated prefixes that the replay of a factorized join
    /// sends through a union host are joins of its components' prefixes, so
    /// they lie inside the host, apart from its soft set. A union host has no
    /// U of its own to replay on a prefix in its soft set.
    #[error(
        "internal error: the integrated prefix {{{}}} lies in the soft set of the union \
         zero-mode host {{{}}}",
        .prefix.iter().join(","),
        .host.iter().join(",")
    )]
    PrefixInUnionSoftSet {
        host: BTreeSet<EdgeIndex>,
        prefix: BTreeSet<EdgeIndex>,
    },
    /// A bug: the frequency classes of a soft set always span its loops.
    #[error(
        "internal error: the frequency classes of the soft set {{{}}} span {rank} of its \
         {loops} loops",
        .sigma.iter().join(",")
    )]
    ClassRank {
        sigma: BTreeSet<EdgeIndex>,
        rank: usize,
        loops: usize,
    },
    #[error(
        "failed to induce the canonical coordinate basis of the zero-mode candidate {{{}}}",
        .host.iter().join(",")
    )]
    CoordinateLmb {
        host: BTreeSet<EdgeIndex>,
        source: LmbError,
    },
    #[error(
        "failed to build a loop basis of the soft set {{{}}} with its hard part contracted",
        .sigma.iter().join(",")
    )]
    QuotientLmb {
        sigma: BTreeSet<EdgeIndex>,
        source: LmbError,
    },
}

/// A needed tuple that the zero-mode projection cannot supply exactly. It
/// is an error only once its candidate becomes a host.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub(crate) enum UnsupportedOrders {
    /// Some class needs an order `a ≤ 2 nu - 2`.
    #[error(
        "the Laurent orders ({}) at Taylor order {taylor_order} lie below the vacuum pole of a \
         class; only the orders 2nu and 2nu-1 are exact",
        .orders.iter().join(",")
    )]
    BelowVacuumPole { orders: Vec<u32>, taylor_order: i64 },
    /// Several classes need their vacuum pole `2 nu - 1` at once.
    #[error(
        "the vacuum pole ({}) of several classes is needed at Taylor order {taylor_order}",
        .orders.iter().join(",")
    )]
    MultiClassVacuumPole { orders: Vec<u32>, taylor_order: i64 },
}

/// A soft set that can carry its static mode: nonempty and eligible, with
/// at least one loop once its hard part is hard.
struct SoftSet {
    edges: BTreeSet<EdgeIndex>,
    loops: usize,
    /// The static infrared degree `ω`.
    omega: i32,
}

/// The soft sets of one cut subgraph.
struct SoftSets<'a> {
    graph: &'a Graph,
    /// The loop lines of the cut subgraph: its paired half-edges less its
    /// bridges, which carry no loop momentum.
    internal: SuBitGraph,
    loops: usize,
    raised_groups: Vec<Vec<EdgeIndex>>,
}

impl<'a> SoftSets<'a> {
    fn new(graph: &'a Graph, cut_subgraph: &SuBitGraph) -> Self {
        let internal = graph.loop_lines(cut_subgraph);
        Self {
            graph,
            loops: graph.n_loops(&internal),
            internal,
            raised_groups: graph.get_raised_edge_groups(),
        }
    }

    fn edges(&self, subgraph: &SuBitGraph) -> BTreeSet<EdgeIndex> {
        self.graph
            .underlying
            .paired_edges(subgraph)
            .into_iter()
            .collect()
    }

    /// Σ for the hard subgraph `hard`, if it can carry its static mode.
    fn soft_set(&self, hard: &SuBitGraph) -> Option<SoftSet> {
        let subgraph = self.internal.subtract(hard);
        let edges = self.edges(&subgraph);
        if !self.is_eligible(&edges) {
            return None;
        }
        let loops = self
            .loops
            .checked_sub(self.graph.n_loops(hard))
            .filter(|loops| *loops >= 1)?;
        let omega = self.graph.static_ir_degree(&subgraph, loops);
        Some(SoftSet {
            edges,
            loops,
            omega,
        })
    }

    /// A massless boson without a chemical potential has a static Matsubara
    /// mode that is not gapped.
    fn has_static_mode(&self, edge: EdgeIndex) -> bool {
        let data = &self.graph[edge];
        data.mass_atom().is_zero() && data.zero_mode().is_ok()
    }

    /// Σ can be soft as a whole: it is nonempty, and each line has a static
    /// mode and carries loop momentum but no graph-external momentum.
    fn is_eligible(&self, sigma: &BTreeSet<EdgeIndex>) -> bool {
        !sigma.is_empty()
            && sigma.iter().all(|edge| {
                let signature = &self.graph.loop_momentum_basis.edge_signatures[*edge];
                self.has_static_mode(*edge)
                    && signature.internal.iter().any(|sign| sign.is_sign())
                    && signature.external.iter().all(|sign| sign.is_zero())
            })
    }

    /// A line of a host that a soft Σ would render massless at small
    /// momentum: a static-mode line carrying crown momentum in `lmb`. A
    /// Taylor expansion of the host in that crown momentum is not safe.
    fn is_unsafe(&self, edge: EdgeIndex, lmb: &LoopMomentumBasis) -> bool {
        self.has_static_mode(edge)
            && lmb.edge_signatures[edge]
                .external
                .iter()
                .any(|sign| !sign.is_zero())
    }

    /// The frequency classes of Σ that get a zero-mode energy. Classes are
    /// tried by decreasing power, then representative, and kept while their
    /// signatures in a loop basis of the quotient `G/host` are independent.
    fn classes(
        &self,
        host: &InternalSubGraph,
        sigma: &BTreeSet<EdgeIndex>,
        soft_loops: usize,
    ) -> Result<Vec<ZeroModeClassSpec>, ZeroModeHostError> {
        let mut candidates = self
            .raised_groups
            .iter()
            .filter_map(|group| {
                ZeroModeClassSpec::new(
                    group
                        .iter()
                        .copied()
                        .filter(|edge| sigma.contains(edge))
                        .collect(),
                )
            })
            .collect_vec();
        candidates.sort_by_key(|class| (Reverse(class.power()), class.representative()));

        let externals = self.graph.dummy_stripped_external_flows_of(&self.internal);
        let quotient = self
            .graph
            .shrunken_sub_lmb(&self.internal, host, externals, None)
            .map_err(|source| ZeroModeHostError::QuotientLmb {
                sigma: sigma.clone(),
                source,
            })?;
        let rows = candidates.into_iter().map(|class| {
            let row = quotient.edge_signatures[class.representative()]
                .internal
                .iter()
                .map(|sign| i64::from(i8::from(*sign)))
                .collect();
            (class, row)
        });
        let classes = extend_independent_rows(&mut Vec::new(), rows, soft_loops);
        if classes.len() < soft_loops {
            return Err(ZeroModeHostError::ClassRank {
                sigma: sigma.clone(),
                rank: classes.len(),
                loops: soft_loops,
            });
        }
        Ok(classes)
    }
}

/// A spinney that may host a region: passes 1 to 5.
struct Candidate {
    index: usize,
    kind: ZeroModeHostKind,
    filter: SuBitGraph,
    edges: BTreeSet<EdgeIndex>,
    /// The filters of its connected components.
    components: Vec<SuBitGraph>,
    soft_set: SoftSet,
    classes: Vec<ZeroModeClassSpec>,
    tuples: Result<Vec<ZeroModeTuple>, UnsupportedOrders>,
}

/// An accepted host and the data its compatibility checks need.
struct Host {
    index: usize,
    filter: SuBitGraph,
    edges: BTreeSet<EdgeIndex>,
    /// The filters of its connected components.
    components: Vec<SuBitGraph>,
    lmb: LoopMomentumBasis,
    region: ThermalZeroModeRegion,
}

impl Host {
    /// Whether this host is a proper subgraph of `outer`.
    fn lies_strictly_inside(&self, outer: &SuBitGraph) -> bool {
        self.filter != *outer && outer.includes(&self.filter)
    }

    /// Whether this host makes the host `inner` redundant: `inner` lies
    /// inside it, the extra soft lines of `inner` are safe in this host's
    /// coordinates, and the static order of `inner` does not exceed this
    /// host's. Then `S_this S_inner = S_inner`.
    fn subsumes(&self, inner: &Self, soft_sets: &SoftSets<'_>) -> bool {
        inner.lies_strictly_inside(&self.filter)
            && inner
                .region
                .sigma
                .difference(&self.region.sigma)
                .all(|line| !soft_sets.is_unsafe(*line, &self.lmb))
            && inner.region.static_taylor_order() <= self.region.static_taylor_order()
    }

    /// Whether this host's region also subtracts the larger soft set
    /// `soft_set`: the extra lines are safe in this host's coordinates, and
    /// its Taylor order leaves the remainder convergent.
    fn certifies(&self, soft_sets: &SoftSets<'_>, soft_set: &SoftSet) -> bool {
        self.region.sigma.is_subset(&soft_set.edges)
            && self.region.sigma != soft_set.edges
            && soft_set
                .edges
                .difference(&self.region.sigma)
                .all(|line| !soft_sets.is_unsafe(*line, &self.lmb))
            && i64::from(soft_set.omega) + self.region.static_taylor_order() + 1 > 0
    }
}

/// The zero-mode host policy of one cut subgraph: passes 1 to 11.
struct HostPolicy<'a> {
    soft_sets: SoftSets<'a>,
    cut_subgraph: &'a SuBitGraph,
    vacuum_subtraction: bool,
}

impl<'a> HostPolicy<'a> {
    fn new(graph: &'a Graph, cut_subgraph: &'a SuBitGraph, medium: &MediumSettings) -> Self {
        Self {
            soft_sets: SoftSets::new(graph, cut_subgraph),
            cut_subgraph,
            vacuum_subtraction: medium.vacuum_subtraction,
        }
    }

    /// Passes 1 to 5 for every spinney.
    fn candidates(&self, spinneys: &[Spinney]) -> Result<Vec<Candidate>, ZeroModeHostError> {
        spinneys
            .iter()
            .enumerate()
            .filter_map(|(index, spinney)| self.candidate(index, spinney, spinneys).transpose())
            .collect()
    }

    /// Passes 1 to 5: a candidate host, its eligible soft set, its classes,
    /// its degree and the tuples it needs.
    fn candidate(
        &self,
        index: usize,
        spinney: &Spinney,
        spinneys: &[Spinney],
    ) -> Result<Option<Candidate>, ZeroModeHostError> {
        let components = spinney
            .components()
            .iter()
            .map(|component| {
                InternalSubGraph::cleaned_filter_optimist(
                    component.clone(),
                    self.soft_sets.graph.as_ref(),
                )
                .filter
            })
            .collect_vec();
        let Some(kind) = Self::candidate_kind(spinney, &components, spinneys) else {
            return Ok(None);
        };
        let Some(soft_set) = self.soft_sets.soft_set(spinney.filter()) else {
            return Ok(None);
        };
        // A convergent static mode needs no tuple; skip its class basis.
        if soft_set.omega > 0 {
            return Ok(None);
        }
        let classes = self
            .soft_sets
            .classes(&spinney.subgraph, &soft_set.edges, soft_set.loops)?;
        let powers = classes.iter().map(ZeroModeClassSpec::power).collect_vec();
        let Some(tuples) =
            ZeroModeTuple::required(&powers, soft_set.omega, self.vacuum_subtraction).transpose()
        else {
            return Ok(None);
        };
        Ok(Some(Candidate {
            index,
            kind,
            filter: spinney.filter().clone(),
            edges: self.soft_sets.edges(spinney.filter()),
            components,
            soft_set,
            classes,
            tuples,
        }))
    }

    /// Pass 1: a connected `IR` spinney of positive degree, or a disjoint
    /// union whose every component, of the cleaned filters `components`, is
    /// one. `classified_spinneys` tags such unions `MUV`, so their components
    /// decide.
    fn candidate_kind(
        spinney: &Spinney,
        components: &[SuBitGraph],
        spinneys: &[Spinney],
    ) -> Option<ZeroModeHostKind> {
        let divergent_ir = |spinney: &Spinney| {
            spinney.n_components() == 1
                && spinney.renormalization_scheme == ApproximationType::IR
                && spinney.dod > 0
        };
        if spinney.thermal_zero_mode().is_some() {
            return None;
        }
        match spinney.n_components() {
            0 => None,
            1 => divergent_ir(spinney).then_some(ZeroModeHostKind::Connected),
            _ => components
                .iter()
                .all(|factor| {
                    spinneys
                        .iter()
                        .any(|spinney| divergent_ir(spinney) && spinney.filter() == factor)
                })
                .then_some(ZeroModeHostKind::Union),
        }
    }

    /// Pass 6, inner candidates first: Taylor safety and zero-bin coverage.
    fn accept_hosts(&self, mut candidates: Vec<Candidate>) -> Result<Vec<Host>, ZeroModeHostError> {
        candidates.sort_by_key(|candidate| (candidate.edges.len(), candidate.kind));
        let mut hosts = Vec::new();
        for candidate in candidates {
            if let Some(host) = self.accept(candidate, &hosts)? {
                hosts.push(host);
            }
        }
        Ok(hosts)
    }

    /// Pass 6 for one candidate, given the hosts accepted inside it so far.
    /// A line of the candidate that a soft set would render massless at
    /// small momentum must lie in the soft set of an inner host, its zero
    /// bin; otherwise the candidate gets no region.
    fn accept(
        &self,
        candidate: Candidate,
        hosts: &[Host],
    ) -> Result<Option<Host>, ZeroModeHostError> {
        let graph = self.soft_sets.graph;
        let lmb = graph
            .canonical_component_lmb(&candidate.filter)
            .map_err(|source| ZeroModeHostError::CoordinateLmb {
                host: candidate.edges.clone(),
                source,
            })?;
        let mut zero_bins = BTreeSet::new();
        for line in candidate
            .edges
            .iter()
            .filter(|line| self.soft_sets.is_unsafe(**line, &lmb))
        {
            let Some(inner) = hosts.iter().find(|inner| {
                inner.lies_strictly_inside(&candidate.filter) && inner.region.sigma.contains(line)
            }) else {
                debug_tags!(#uv, #generation, #zero_mode;
                    candidate = %candidate.edges.iter().join(","),
                    line = %line,
                    "Zero-mode candidate has a soft line without a zero bin; it gets no region"
                );
                return Ok(None);
            };
            zero_bins.insert(inner.filter.clone());
        }
        let tuples = candidate
            .tuples
            .map_err(|reason| ZeroModeHostError::UnsupportedOrders {
                host: candidate.edges.clone(),
                reason,
            })?;
        let region = ThermalZeroModeRegion {
            kind: candidate.kind,
            sigma: candidate.soft_set.edges,
            classes: candidate.classes,
            omega: candidate.soft_set.omega,
            soft_loops: candidate.soft_set.loops,
            tuples,
            zero_bins: zero_bins.into_iter().collect(),
            coordinate_loop_edges: lmb.loop_edges.iter().copied().collect(),
        };
        Ok(Some(Host {
            index: candidate.index,
            filter: candidate.filter,
            edges: candidate.edges,
            components: candidate.components,
            lmb,
            region,
        }))
    }

    /// Pass 7: an inner host is redundant when an outer host subsumes it, so
    /// its terms cancel in the forest. A zero bin is never redundant.
    fn prune(&self, hosts: Vec<Host>) -> Vec<Host> {
        let zero_bins = hosts
            .iter()
            .flat_map(|host| host.region.zero_bins.iter().cloned())
            .collect::<BTreeSet<_>>();
        let redundant = hosts
            .iter()
            .map(|inner| {
                !zero_bins.contains(&inner.filter)
                    && hosts
                        .iter()
                        .any(|outer| outer.subsumes(inner, &self.soft_sets))
            })
            .collect_vec();
        hosts
            .into_iter()
            .zip(redundant)
            .filter_map(|(host, redundant)| {
                if redundant {
                    debug_tags!(#uv, #generation, #zero_mode;
                        host = %host.edges.iter().join(","),
                        "Pruned a zero-mode host, which an enclosing host subsumes"
                    );
                }
                (!redundant).then_some(host)
            })
            .collect()
    }

    /// Pass 8: hosts form a chain, only a union host contains another, as
    /// one of its components, and hosts need one uncut structure, since
    /// spinneys are merged across cuts.
    fn check_nesting(hosts: &[Host], cuts: &CutStructure) -> Result<(), ZeroModeHostError> {
        for (position, first) in hosts.iter().enumerate() {
            for second in &hosts[position + 1..] {
                let (inner, outer) = if second.filter.includes(&first.filter) {
                    (first, second)
                } else if first.filter.includes(&second.filter) {
                    (second, first)
                } else {
                    return Err(ZeroModeHostError::IncompatibleHosts {
                        first: first.edges.clone(),
                        second: second.edges.clone(),
                    });
                };
                if outer.region.kind != ZeroModeHostKind::Union
                    || !outer.components.contains(&inner.filter)
                {
                    return Err(ZeroModeHostError::UnsupportedNesting {
                        outer: outer.edges.clone(),
                        inner: inner.edges.clone(),
                    });
                }
            }
        }
        if let Some(host) = hosts.first()
            && !cuts.is_trivial()
        {
            return Err(ZeroModeHostError::NonTrivialCut {
                host: host.edges.clone(),
            });
        }
        Ok(())
    }

    /// Pass 9: every divergent soft set has a host. The spinney walk never
    /// visits soft sets whose hard part converges in the UV, so check every
    /// cycle union as a hard part. The first hostless soft set, in edge
    /// order, is the error.
    fn check_hostless(&self, hosts: &[Host]) -> Result<(), ZeroModeHostError> {
        // A cut subgraph without a loop has no soft set, and no cycle union.
        if self.soft_sets.loops == 0 {
            return Ok(());
        }
        let graph = self.soft_sets.graph;
        let hostless = graph
            .all_cycle_unions(self.cut_subgraph)
            .into_iter()
            .map(|cycles| cycles.filter)
            .chain([graph.underlying.empty_subgraph::<SuBitGraph>()])
            .filter_map(|hard| self.soft_sets.soft_set(&hard))
            .filter(|soft_set| {
                soft_set.omega <= 0
                    && !hosts.iter().any(|host| {
                        host.region.sigma == soft_set.edges
                            || host.certifies(&self.soft_sets, soft_set)
                    })
            })
            .map(|soft_set| (soft_set.edges, soft_set.omega))
            .min();
        match hostless {
            Some((sigma, omega)) => Err(ZeroModeHostError::Hostless { sigma, omega }),
            None => Ok(()),
        }
    }

    /// Passes 10 and 11: install the regions, a union host as a copy of its
    /// plain union, and demote every other `IR` spinney to a plain U. The
    /// empty root is no operation.
    fn install(&self, hosts: Vec<Host>, spinneys: &mut Vec<Spinney>) {
        let mut connected = BTreeMap::new();
        let mut unions = Vec::new();
        for host in hosts {
            match host.region.kind {
                ZeroModeHostKind::Connected => {
                    connected.insert(host.index, host.region);
                }
                ZeroModeHostKind::Union => {
                    unions.push(
                        spinneys[host.index]
                            .clone()
                            .with_thermal_zero_mode(host.region),
                    );
                }
            }
        }
        *spinneys = mem::take(spinneys)
            .into_iter()
            .enumerate()
            .map(|(index, mut spinney)| {
                if let Some(region) = connected.remove(&index) {
                    return spinney.with_thermal_zero_mode(region);
                }
                if spinney.renormalization_scheme == ApproximationType::IR
                    && spinney.thermal_zero_mode().is_none()
                    && !spinney.subgraph.is_empty()
                {
                    debug_tags!(#uv, #generation, #zero_mode;
                        spinney = %self.soft_sets.edges(spinney.filter()).iter().join(","),
                        dod = spinney.dod,
                        "Demoted an IR spinney without a zero-mode region to MUV"
                    );
                    spinney.demote_to_muv();
                }
                spinney
            })
            .chain(unions)
            .collect();
    }

    /// Log the classification the policy leaves.
    fn log_classification(&self, spinneys: &[Spinney]) {
        for spinney in spinneys {
            debug_tags!(#uv, #generation, #zero_mode;
                spinney = %self.soft_sets.edges(spinney.filter()).iter().join(","),
                scheme = %spinney.renormalization_scheme,
                dod = spinney.dod,
                region = %spinney
                    .thermal_zero_mode()
                    .map(ToString::to_string)
                    .unwrap_or_default(),
                "Zero-mode host policy classification"
            );
        }
    }
}

/// Assign the thermal zero-mode hosts of one cut subgraph to its classified
/// `spinneys`, and demote every other `IR` spinney to `MUV`.
///
/// The pass runs iff UV subtraction is on and the route freezes completed
/// thermal soft branches, in both vacuum-subtraction modes; otherwise it
/// leaves `spinneys` unchanged. A connected host gets its region on its own
/// spinney. A union host is pushed as an `IR` copy of the union, which stays
/// a plain `MUV` join.
///
/// # Errors
///
/// - a soft set of non-positive degree that no host subtracts;
/// - hosts that are not nested, or nested other than as a component of a
///   union host;
/// - a host that needs a Laurent order the projection cannot supply;
/// - a host on a structure of `cuts` other than one uncut subgraph, whose
///   spinneys are merged across cuts;
/// - internal and loop-momentum-basis errors.
pub(crate) fn assign_zero_mode_hosts(
    graph: &Graph,
    cut_subgraph: &SuBitGraph,
    spinneys: &mut Vec<Spinney>,
    settings: &UVgenerationSettings,
    medium: &MediumSettings,
    cuts: &CutStructure,
) -> Result<(), ZeroModeHostError> {
    if !(settings.subtract_uv && settings.freezes_thermal_soft_branches(medium)) {
        return Ok(());
    }
    let policy = HostPolicy::new(graph, cut_subgraph, medium);
    let candidates = policy.candidates(spinneys)?;
    let hosts = policy.prune(policy.accept_hosts(candidates)?);
    HostPolicy::check_nesting(&hosts, cuts)?;
    policy.check_hostless(&hosts)?;
    policy.install(hosts, spinneys);
    policy.log_classification(spinneys);
    Ok(())
}
