<<<<<<< qrspmutx ee257326 "add feynkit ufo and Python bindings" (rebase destination)
use std::{
    collections::{BTreeMap, BTreeSet, VecDeque},
    hash::Hash,
    mem::{size_of, size_of_val},
    sync::Arc,
    time::{Duration, Instant},
};
||||||| qrspmutx 7a048a39 "add feynkit ufo and Python bindings" (parents of rebased revision)
use std::collections::HashMap;
=======
//! GammaLoop runtime extensions for FeynKit-owned CFF generation.
>>>>>>> xzuonylm 3e5f0ddf (rebased revision)

<<<<<<< qrspmutx ee257326 "add feynkit ufo and Python bindings" (rebase destination)
use crate::utils::GS;
use crate::{
    cff::{
        VertexSet,
        hsurface::{Hsurface, HsurfaceID},
        surface::{
            HybridSurfaceID, LinearEnergyExpr, LinearSurface, LinearSurfaceID, LinearSurfaceKind,
        },
    },
    graph::{
        ExactUvSubLmbFrame, FeynmanGraph, FourDDenominator, Graph, GraphThreeDSource,
        LoopMomentumBasis,
    },
    numerator::energy_degree::EnergyPowerAssignmentPlan,
    settings::global::{GenerationSettings, UniformNumeratorSamplingScale},
    uv::approx::{
        local_4d::{CanonicalUvDenominatorClass, CanonicalUvSector, FourDSector},
        projected_4d::Local4dProjectionContext,
    },
};
use ahash::{AHashMap, HashSet};
||||||| qrspmutx 7a048a39 "add feynkit ufo and Python bindings" (parents of rebased revision)
use crate::{
    cff::{
        expression::OrientationData,
        hsurface::{Hsurface, HsurfaceID},
        surface::{HybridSurface, HybridSurfaceID, InfiniteSurface},
        tree::Tree,
    },
    graph::{Graph, LoopMomentumBasis, get_cff_inverse_energy_product_impl},
    processes::{CrossSectionCut, CutId},
    settings::global::OrientationPattern,
};
use ahash::HashSet;
use bincode::{Decode, Encode};
use color_eyre::Report;
=======
>>>>>>> xzuonylm 3e5f0ddf (rebased revision)
use color_eyre::Result;
<<<<<<< qrspmutx ee257326 "add feynkit ufo and Python bindings" (rebase destination)
use itertools::Itertools;
use linnet::half_edge::{
    involution::{EdgeIndex, Hedge},
    subgraph::{SuBitGraph, SubSetLike},
};
use linnet::num_traits::SignOrZero;
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::domains::rational::Rational;
use three_dimensional_reps::{
    CffEnergyFactorOwnership, CffGenerationContext, EnergyEdgeIndexMap,
    Generate3DExpressionOptions, GeneratedThreeDExpression, NumeratorSamplingScaleMode,
    ParsedGraph, RepresentationMode, ThreeDGraphSource,
    tree::{NodeId, Tree},
};

use tracing::debug;

use super::{
    CffEnergyBoundSourceKind, CffEnergyDegreeBoundReport, PlannedExactSourceNumerator,
    esurface::{Esurface, EsurfaceID, ExternalShift},
    expression::CFFExpression,
};

#[derive(Debug, Clone)]
pub struct ShiftRewrite {
    pub dependent_momentum: EdgeIndex,
    pub dependent_momentum_expr: ExternalShift,
}

#[derive(Debug, Clone, Copy)]
struct SurfaceMapEntry {
    surface_id: HybridSurfaceID,
    sign: i64,
}

struct RemappedDenominatorTree {
    prefactor: i64,
    surface_signs: BTreeMap<HybridSurfaceID, i64>,
    denominator: Tree<HybridSurfaceID>,
}

/// Generation-local, deterministic LRU retention. Callers charge owned keys,
/// containers and symbolic payload at insertion; this is not a process RSS cap.
/// Cache state changes work only, never the admitted plans or their semantics.
pub(crate) struct GenerationCache<K, V> {
    entries: AHashMap<K, GenerationCacheEntry<V>>,
    recency: BTreeMap<u64, K>,
    clock: u64,
    max_bytes: usize,
    max_entries: usize,
    retained_bytes: usize,
    pub(crate) hits: usize,
    pub(crate) misses: usize,
    pub(crate) evictions: usize,
    // Callers include key construction and accounting, excluding value work.
    pub(crate) cache_time: Duration,
}

struct GenerationCacheEntry<V> {
    value: V,
    bytes: usize,
    used: u64,
}

impl<K: Clone + Eq + Hash, V> GenerationCache<K, V> {
    pub(crate) fn new(max_bytes: usize, max_entries: usize) -> Self {
        Self {
            entries: AHashMap::new(),
            recency: BTreeMap::new(),
            clock: 0,
            max_bytes,
            max_entries,
            retained_bytes: 0,
            hits: 0,
            misses: 0,
            evictions: 0,
            cache_time: Duration::ZERO,
        }
    }

    pub(crate) fn get(&mut self, key: &K) -> Option<&V> {
        let Some(entry) = self.entries.get_mut(key) else {
            self.misses += 1;
            return None;
        };
        self.hits += 1;
        self.clock = self
            .clock
            .checked_add(1)
            .expect("generation cache access counter exhausted");
        let retained_key = self
            .recency
            .remove(&entry.used)
            .expect("LRU index and cache entries agree");
        entry.used = self.clock;
        self.recency.insert(self.clock, retained_key);
        Some(&entry.value)
    }

    pub(crate) fn insert(&mut self, key: K, value: V, accounted_bytes: usize) {
        if let Some(previous) = self.entries.remove(&key) {
            self.retained_bytes -= previous.bytes;
            self.recency.remove(&previous.used);
        }
        if self.max_entries == 0 || self.max_bytes == 0 || accounted_bytes > self.max_bytes {
            return;
        }
        while self.entries.len() >= self.max_entries
            || accounted_bytes > self.max_bytes - self.retained_bytes
        {
            let Some((_, oldest)) = self.recency.pop_first() else {
                break;
            };
            let removed = self
                .entries
                .remove(&oldest)
                .expect("LRU index and cache entries agree");
            self.retained_bytes -= removed.bytes;
            self.evictions += 1;
        }
        self.clock = self
            .clock
            .checked_add(1)
            .expect("generation cache access counter exhausted");
        self.retained_bytes += accounted_bytes;
        self.recency.insert(self.clock, key.clone());
        self.entries.insert(
            key,
            GenerationCacheEntry {
                value,
                bytes: accounted_bytes,
                used: self.clock,
            },
        );
    }

    pub(crate) fn clear(&mut self) {
        let started = Instant::now();
        self.entries.clear();
        self.recency.clear();
        self.retained_bytes = 0;
        self.cache_time += started.elapsed();
    }

    pub(crate) fn len(&self) -> usize {
        self.entries.len()
    }
    pub(crate) fn can_retain(&self, accounted_bytes: usize) -> bool {
        self.max_entries != 0 && self.max_bytes != 0 && accounted_bytes <= self.max_bytes
    }
    pub(crate) fn retained_bytes(&self) -> usize {
        self.retained_bytes
    }
}

#[derive(Clone, PartialEq, Eq, Hash)]
struct ExactCffGenerationKey {
    topology: ParsedGraph,
    internal_energy_edges: BTreeMap<usize, usize>,
    external_energy_edges: BTreeMap<usize, usize>,
    orientation_edge_count: usize,
    options: Generate3DExpressionOptions,
}

struct ExactCffGenerationEntry {
    count: usize,
    payload: Option<GeneratedThreeDExpression>,
}

impl ExactCffGenerationKey {
    fn accounted_bytes(&self) -> usize {
        let graph = &self.topology;
        size_of::<Self>()
            + ExactCffGenerationCache::vector_bytes(&graph.internal_edges)
            + graph
                .internal_edges
                .iter()
                .map(|edge| {
                    edge.label.capacity()
                        + edge.mass_key.as_ref().map_or(0, String::capacity)
                        + ExactCffGenerationCache::vector_bytes(&edge.signature.loop_signature)
                        + ExactCffGenerationCache::vector_bytes(&edge.signature.external_signature)
                })
                .sum::<usize>()
            + ExactCffGenerationCache::vector_bytes(&graph.external_edges)
            + graph
                .external_edges
                .iter()
                .map(|edge| {
                    edge.label.capacity()
                        + ExactCffGenerationCache::vector_bytes(&edge.external_coefficients)
                })
                .sum::<usize>()
            + ExactCffGenerationCache::vector_bytes(&graph.initial_state_cut_edges)
            + ExactCffGenerationCache::vector_bytes(&graph.loop_names)
            + ExactCffGenerationCache::vector_bytes(&graph.external_names)
            + graph
                .loop_names
                .iter()
                .chain(&graph.external_names)
                .map(String::capacity)
                .sum::<usize>()
            + graph
                .node_name_to_internal
                .keys()
                .map(|name| size_of::<String>() + name.capacity() + size_of::<usize>() + 32)
                .sum::<usize>()
            + (self.internal_energy_edges.len() + self.external_energy_edges.len())
                * (2 * size_of::<usize>() + 32)
            + self
                .options
                .energy_degree_bounds
                .as_ref()
                .map_or(0, ExactCffGenerationCache::vector_bytes)
            + ExactCffGenerationCache::vector_bytes(
                &self.options.preserve_internal_edges_as_four_d_denominators,
            )
    }
}

pub(crate) struct ExactCffGenerationCache {
    // Keep each occurrence's requested capacity, even when source provenance
    // identifies equal physical energies. Reuse requires canonical topology
    // and identical per-edge bounds; independent terms must not combine into
    // a larger Cartesian capacity or redistribute their numerator ownership.
    // Losing trials keep only their count, so repeated terms can compare the
    // same complete source capacity without retaining or rebuilding its trees.
    entries: GenerationCache<ExactCffGenerationKey, ExactCffGenerationEntry>,
    pub(crate) native_generations: usize,
}

impl Default for ExactCffGenerationCache {
    fn default() -> Self {
        Self {
            entries: GenerationCache::new(64 * 1024 * 1024, 4096),
            native_generations: 0,
        }
    }
}

impl ExactCffGenerationCache {
    pub(crate) fn statistics(&self) -> (usize, usize, usize, usize) {
        (
            self.entries.hits,
            self.entries.misses,
            self.entries.evictions,
            self.entries.retained_bytes(),
        )
    }

    pub(crate) fn len(&self) -> usize {
        self.entries
            .entries
            .values()
            .filter(|entry| entry.value.payload.is_some())
            .count()
    }

    fn count(&mut self, key: &ExactCffGenerationKey) -> Option<usize> {
        self.entries.get(key).map(|entry| entry.count)
    }

    fn record_count(&mut self, key: ExactCffGenerationKey, count: usize) {
        let bytes = key.accounted_bytes() * 2 + size_of::<ExactCffGenerationEntry>() + 64;
        self.entries.insert(
            key,
            ExactCffGenerationEntry {
                count,
                payload: None,
            },
            bytes,
        );
    }

    fn record_payload(&mut self, key: ExactCffGenerationKey, payload: &GeneratedThreeDExpression) {
        let count = payload.expression.orientations.len();
        let bytes = key.accounted_bytes() * 2
            + size_of::<ExactCffGenerationEntry>()
            + 64
            + Self::payload_bytes(payload);
        if bytes > self.entries.max_bytes {
            self.record_count(key, count);
        } else {
            self.entries.insert(
                key,
                ExactCffGenerationEntry {
                    count,
                    payload: Some(payload.clone()),
                },
                bytes,
            );
        }
    }

    fn rational_bytes(value: &Rational) -> usize {
        [value.numerator_ref(), value.denominator_ref()]
            .into_iter()
            .map(|integer| {
                if let symbolica::domains::integer::Integer::Large(value) = integer {
                    value.significant_bits().div_ceil(8) as usize
                } else {
                    0
                }
            })
            .sum()
    }

    fn vector_bytes<T>(values: &Vec<T>) -> usize {
        values.capacity() * size_of::<T>()
    }

    fn energy_bytes(energy: &LinearEnergyExpr) -> usize {
        Self::vector_bytes(&energy.internal_terms)
            + Self::vector_bytes(&energy.external_terms)
            + energy
                .internal_terms
                .iter()
                .chain(&energy.external_terms)
                .map(|(_, coefficient)| Self::rational_bytes(coefficient))
                .sum::<usize>()
            + Self::rational_bytes(&energy.constant)
            + Self::rational_bytes(&energy.uniform_scale_coeff)
    }

    fn payload_bytes(payload: &GeneratedThreeDExpression) -> usize {
        // Native payloads contain rational data and owned trees, not mapped
        // numerator Atoms. Count owned vector/string capacity once at insertion;
        // the private tree/EdgeVec containers are cloned to their live length.
        // No printing, serialization or per-hit traversal is necessary.
        let expression = &payload.expression;
        size_of_val(payload)
            + Self::vector_bytes(&payload.source_energy_degree_bounds)
            + Self::vector_bytes(&payload.energy_factor_components)
            + payload
                .energy_factor_components
                .iter()
                .map(|component| Self::vector_bytes(&component.internal_edge_ids))
                .sum::<usize>()
            + expression.surfaces.linear_surface_cache.capacity() * size_of::<LinearSurface>()
            + expression
                .surfaces
                .linear_surface_cache
                .iter()
                .map(|surface| Self::energy_bytes(&surface.expression))
                .sum::<usize>()
            + Self::vector_bytes(&expression.residual_denominators)
            + expression
                .residual_denominators
                .iter()
                .map(|denominator| denominator.origin.as_ref().map_or(0, String::capacity))
                .sum::<usize>()
            + expression.orientations.capacity()
                * size_of::<three_dimensional_reps::expression::OrientationExpression>()
            + expression
                .orientations
                .iter()
                .map(|orientation| {
                    orientation.data.orientation.iter().count()
                        * size_of::<linnet::half_edge::involution::Orientation>()
                        + orientation.data.label.as_ref().map_or(0, String::capacity)
                        + Self::vector_bytes(&orientation.loop_energy_map)
                        + Self::vector_bytes(&orientation.edge_energy_map)
                        + orientation
                            .loop_energy_map
                            .iter()
                            .chain(&orientation.edge_energy_map)
                            .map(Self::energy_bytes)
                            .sum::<usize>()
                        + Self::vector_bytes(&orientation.variants)
                        + orientation
                            .variants
                            .iter()
                            .map(|variant| {
                                variant.origin.as_ref().map_or(0, String::capacity)
                                    + Self::rational_bytes(&variant.prefactor)
                                    + Self::vector_bytes(&variant.half_edges)
                                    + Self::vector_bytes(&variant.denominator_edges)
                                    + Self::vector_bytes(&variant.numerator_surfaces)
                                    + variant.denominator_surface_signs.len()
                                        * (size_of::<(HybridSurfaceID, i64)>() + 32)
                                    + variant
                                        .denominator_edge_support_signs
                                        .keys()
                                        .map(|support| {
                                            size_of::<(Vec<EdgeIndex>, i64)>()
                                                + Self::vector_bytes(support)
                                                + 32
                                        })
                                        .sum::<usize>()
                                    + variant
                                        .denominator
                                        .iter_nodes()
                                        .map(|node| {
                                            size_of_val(node) + Self::vector_bytes(&node.children)
                                        })
                                        .sum::<usize>()
                            })
                            .sum::<usize>()
                })
                .sum::<usize>()
    }

    fn generation_key(
        parsed: &ParsedGraph,
        energy_edges: &EnergyEdgeIndexMap,
        options: &Generate3DExpressionOptions,
    ) -> ExactCffGenerationKey {
        let mut topology = parsed.clone();
        topology.node_name_to_internal = topology
            .internal_edges
            .iter()
            .flat_map(|edge| [edge.tail, edge.head])
            .chain(
                topology
                    .external_edges
                    .iter()
                    .flat_map(|edge| edge.source.into_iter().chain(edge.destination)),
            )
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .map(|node| (format!("__gammaloop_exact_node_{node}"), node))
            .collect();
        ExactCffGenerationKey {
            topology,
            internal_energy_edges: energy_edges.internal.clone(),
            external_energy_edges: energy_edges.external.clone(),
            orientation_edge_count: energy_edges.orientation_edge_count,
            options: options.clone(),
        }
    }
}

/// These disjoint identities share one preparation budget and LRU order.
#[derive(Clone, PartialEq, Eq, Hash)]
pub(crate) enum PreparationKey {
    Canonical(Arc<FourDSector>),
    Source(Arc<ExactCffPreparationKey>),
    Numerator(super::ExactNumeratorTemplateKey),
}

pub(crate) enum PreparationValue {
    Canonical(Arc<CanonicalUvSector>),
    Source(Arc<ExactCffGenerationPreparation>),
    Numerator(Arc<super::PlannedExactSourceNumerator>),
}

#[derive(PartialEq, Eq, Hash)]
pub(crate) struct ExactCffPreparationKey {
    pub(super) denominators: Vec<FourDDenominator>,
    pub(super) uv_edges: Vec<EdgeIndex>,
    pub(super) boundary_hedges: Vec<Hedge>,
    pub(super) coordinates: Option<(LoopMomentumBasis, ExactUvSubLmbFrame)>,
    pub(super) classes: Vec<CanonicalUvDenominatorClass>,
    pub(super) numerator: Atom,
    pub(super) options: Generate3DExpressionOptions,
}

impl ExactCffPreparationKey {
    pub(super) fn accounted_bytes(&self) -> usize {
        size_of::<Self>()
            + self.denominators.capacity() * size_of::<FourDDenominator>()
            + self
                .denominators
                .iter()
                .map(FourDDenominator::accounted_bytes)
                .sum::<usize>()
            + self.uv_edges.capacity() * size_of::<EdgeIndex>()
            + self.boundary_hedges.capacity() * size_of::<Hedge>()
            + self
                .coordinates
                .as_ref()
                .map_or(0, |(lmb, _)| lmb.accounted_bytes())
            + self.classes.capacity() * size_of::<CanonicalUvDenominatorClass>()
            + self
                .classes
                .iter()
                .map(CanonicalUvDenominatorClass::accounted_bytes)
                .sum::<usize>()
            + self.numerator.as_view().get_byte_size()
            + self
                .options
                .energy_degree_bounds
                .as_ref()
                .map_or(0, ExactCffGenerationCache::vector_bytes)
            + ExactCffGenerationCache::vector_bytes(
                &self.options.preserve_internal_edges_as_four_d_denominators,
            )
    }
}

pub(crate) struct ExactCffGenerationPreparation {
    parsed: ParsedGraph,
    energy_edges: EnergyEdgeIndexMap,
    source_options: Generate3DExpressionOptions,
    exact_source_energy_mapper: Arc<crate::graph::three_d_source::ExactSourceEnergyMapper>,
    energy_assignment_plans: Vec<Arc<EnergyPowerAssignmentPlan>>,
    physical_energy_degree_bounds: Vec<(usize, usize)>,
    pub(super) physical_energy_edges: EnergyEdgeIndexMap,
    pub(super) physical_cut_support_edges: BTreeMap<usize, Vec<EdgeIndex>>,
    pub(super) physical_surface_edges: BTreeSet<usize>,
    pub(super) inverse_energy_product: Atom,
    pub(super) active_loop_count: usize,
    pub(super) contract_subgraph: SuBitGraph,
}

impl ThreeDGraphSource for ExactCffGenerationPreparation {
    fn to_three_d_parsed_graph(&self) -> three_dimensional_reps::graph_io::Result<ParsedGraph> {
        Ok(self.parsed.clone())
    }

    fn energy_edge_index_map(&self, _: &ParsedGraph) -> Option<EnergyEdgeIndexMap> {
        Some(self.energy_edges.clone())
    }
}

impl ExactCffGenerationPreparation {
    pub(super) fn accounted_bytes(&self) -> Result<usize> {
        // Charge every Arc payload reachable from this entry independently.
        // Another preparation may retain the same payload after this is evicted.
        Ok(size_of::<Self>()
            + 128
            + ExactCffGenerationKey {
                topology: self.parsed.clone(),
                internal_energy_edges: self.energy_edges.internal.clone(),
                external_energy_edges: self.energy_edges.external.clone(),
                orientation_edge_count: self.energy_edges.orientation_edge_count,
                options: self.source_options.clone(),
            }
            .accounted_bytes()
            + self.exact_source_energy_mapper.accounted_bytes()?
            + self.energy_assignment_plans.capacity() * size_of::<Arc<EnergyPowerAssignmentPlan>>()
            + self
                .energy_assignment_plans
                .iter()
                .map(|plan| 32 + plan.accounted_bytes())
                .sum::<usize>()
            + self.physical_energy_degree_bounds.capacity() * size_of::<(usize, usize)>()
            + (self.physical_energy_edges.internal.len()
                + self.physical_energy_edges.external.len())
                * 64
            + self
                .physical_cut_support_edges
                .values()
                .map(|edges| 64 + edges.capacity() * size_of::<EdgeIndex>())
                .sum::<usize>()
            + self.physical_surface_edges.len() * 48
            + self.inverse_energy_product.as_view().get_byte_size()
            + self.contract_subgraph.size().div_ceil(8)
            + 32)
    }
}

impl Graph {
    pub(crate) fn prepare_3d_expression_for_4d_term(
        &self,
        source: &GraphThreeDSource<'_>,
        options: &Generate3DExpressionOptions,
        analysis_numerator: &Atom,
        classes: &[CanonicalUvDenominatorClass],
    ) -> Result<ExactCffGenerationPreparation> {
        let preparation_started = Instant::now();
        let source_options = options.clone();
        let initial_state_cut_edges = self
            .iter_edges_of(&self.initial_state_cut)
            .map(|(_, edge_id, _)| edge_id)
            .collect::<HashSet<_>>();
        let bridge_edges = self
            .iter_edges_of(&self.tree_edges)
            .map(|(_, edge_id, _)| edge_id)
            .filter(|edge_id| !initial_state_cut_edges.contains(edge_id))
            .collect::<HashSet<_>>();
        let excluded_numerator_edges = initial_state_cut_edges
            .iter()
            .chain(&bridge_edges)
            .copied()
            .chain(source.factorized_external_emr_edges())
            .collect::<HashSet<_>>();
        let contracted_source = source.contract_subgraph();
        for (_, edge_id, _) in self.iter_edges_of(&contracted_source) {
            if !bridge_edges.contains(&edge_id) {
                continue;
            }
            let Some(coordinates) = source.reconstructible_outer_loop_coordinates(edge_id) else {
                return Err(eyre::eyre!(
                    "contracted bridge edge {} retains an inner-loop energy coordinate and cannot be projected independently of the exact CFF source",
                    usize::from(edge_id),
                ));
            };
            if coordinates.iter().any(|coordinate| *coordinate != 0) {
                return Err(eyre::eyre!(
                    "contracted bridge edge {} has nonzero outer-loop coordinates {:?}; a CFF-external tree edge must have a purely external affine energy",
                    usize::from(edge_id),
                    coordinates,
                ));
            }
        }
        // This report lives in the physical parent's namespace. Canonical
        // classes are expanded only for degree reporting, never for dispatch;
        // remote component classes carry no energy of this source.
        let degree_analysis_started = Instant::now();
        let (prepared_numerator, fixed_affine_blocks) =
            source.prepare_affine_numerator(analysis_numerator)?;
        let analysis_numerator = &prepared_numerator;
        let physical_report_numerator = analysis_numerator.replace_map(|view, _, output| {
            if let AtomView::Fun(denominator) = view
                && denominator.get_symbol() == GS.den
                && denominator.get_nargs() == 4
            {
                **output = denominator.get(3).to_owned();
            }
        });
        let physical_report_numerator = physical_report_numerator.replace_map(|view, _, output| {
            let AtomView::Fun(momentum) = view else {
                return;
            };
            if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() < 2 {
                return;
            }
            let Some(id) = GS.uv_class_data(momentum.get(0)) else {
                return;
            };
            **output = classes
                .iter()
                .find(|class| class.id == id)
                .map_or(Atom::Zero, |class| {
                    class.momentum_with_indices(
                        &momentum
                            .iter()
                            .skip(1)
                            .map(|index| index.to_owned())
                            .collect::<Vec<_>>(),
                    )
                });
        });
        let physical_energy_degree_bounds = self
            .automatic_numerator_energy_degree_bounds_in_atoms_excluding_with_min_degree(
                [&physical_report_numerator],
                excluded_numerator_edges.iter().copied(),
                1,
            )
            .map_err(|error| {
                eyre::eyre!("could not analyze numerator in physical EMR energy variables: {error}")
            })?;
        let degree_analysis_ms = degree_analysis_started.elapsed().as_secs_f64() * 1000.0;
        // Parse first so a malformed exact rational source returns its
        // structural error instead of being hidden behind the mapper's
        // optional convenience API.
        let source_mapping_started = Instant::now();
        let parsed = source.to_three_d_parsed_graph()?;
        let energy_edges = source
            .energy_edge_index_map(&parsed)
            .expect("exact 4D source has an occurrence-local energy map");
        let exact_source_energy_mapper = source.exact_source_energy_mapper(classes)?;
        // Exact sources have occurrence-local denominator IDs. Physical
        // originals retain their base occurrence and their derived factors may
        // use serial copies; canonical UV classes use their certified pools.
        // Analyze all physical active edges first so unused, unrelated
        // candidate groups cannot reject a constant numerator.
        let mut candidates = exact_source_energy_mapper
            .equivalent_energy_candidates(
                physical_energy_degree_bounds
                    .iter()
                    .map(|(edge, _)| EdgeIndex(*edge)),
            )
            .map_err(|error| {
                eyre::eyre!(
                    "could not certify exact 4D CFF numerator energies for graph `{}`: {error}",
                    self.name,
                )
            })?;
        candidates.fixed_affine_blocks = fixed_affine_blocks;
        let source_mapping_ms = source_mapping_started.elapsed().as_secs_f64() * 1000.0;
        // Each immutable factor-local plan owns both its exact bounds and the
        // later numerator substitutions. This keeps the numerator factorized
        // and prevents generation from understating the expression actually
        // sampled in a residue or contact sector. Rank proposes a bounded set
        // of plans; the real source map count chooses between them below.
        let allocation_started = Instant::now();
        let energy_assignment_plans = self
            .plan_numerator_energy_assignment_proposals_in_atom_excluding(
                analysis_numerator,
                excluded_numerator_edges.iter().copied(),
                &candidates,
            )
            .map_err(|error| {
                eyre::eyre!(
                    "could not plan exact 4D CFF numerator energies for graph `{}`: {error}",
                    self.name,
                )
            })?;
        let allocation_ms = allocation_started.elapsed().as_secs_f64() * 1000.0;
        debug!(
            graph = %self.name,
            physical_energy_degree_bounds = ?physical_energy_degree_bounds,
            equivalent_energy_candidates = ?candidates,
            candidate_bounds = ?energy_assignment_plans.iter().map(|plan| plan.energy_degree_bounds()).collect::<Vec<_>>(),
            "planned factorized exact-CFF numerator energy assignment proposals"
        );
        let preparation = ExactCffGenerationPreparation {
            parsed,
            energy_edges,
            source_options,
            exact_source_energy_mapper: Arc::new(exact_source_energy_mapper),
            energy_assignment_plans: energy_assignment_plans.into_iter().map(Arc::new).collect(),
            physical_energy_degree_bounds,
            physical_energy_edges: source
                .physical_energy_edge_index_map()
                .expect("exact source has a physical energy projection"),
            physical_cut_support_edges: source
                .physical_cut_support_edge_index_map()
                .expect("exact source has a cut-support projection"),
            physical_surface_edges: source.physical_surface_energy_edges(),
            inverse_energy_product: source
                .exact_inverse_energy_product()
                .expect("exact source has an inverse energy product"),
            active_loop_count: source.active_loop_count(),
            contract_subgraph: source.contract_subgraph(),
        };
        crate::debug_tags!(#generation, #uv, #local, #four_d, #cff, #profile;
            stage = "exact_cff_preparation",
            graph = %self.name,
            elapsed_ms = preparation_started.elapsed().as_secs_f64() * 1000.0,
            degree_analysis_ms,
            source_mapping_ms,
            allocation_ms,
            proposal_count = preparation.energy_assignment_plans.len(),
            "Prepared factorized energy assignments and certified source mapping"
        );
        Ok(preparation)
    }

    pub(crate) fn generate_3d_expression_for_4d_term(
        &self,
        preparation: &ExactCffGenerationPreparation,
        mut context: Option<&mut Local4dProjectionContext>,
    ) -> Result<(
        GeneratedThreeDExpression,
        Arc<PlannedExactSourceNumerator>,
        Arc<EnergyPowerAssignmentPlan>,
        CffEnergyDegreeBoundReport,
    )> {
        let ExactCffGenerationPreparation {
            parsed,
            energy_edges,
            source_options,
            exact_source_energy_mapper,
            energy_assignment_plans,
            physical_energy_degree_bounds,
            ..
        } = preparation;
        let mut source_options = source_options.clone();
        let selection_started = std::time::Instant::now();
        let mut certificate_time = std::time::Duration::ZERO;
        let mut native_times = BTreeMap::new();
        let mut template_build_times = BTreeMap::new();
        let mut template_time = Duration::ZERO;
        let generate = |source_options: &Generate3DExpressionOptions| {
            crate::debug_tags!(#generation, #uv, #local, #four_d, #cff, #profile;
                graph = %self.name,
                term_local_bounds = ?source_options.energy_degree_bounds,
                file.parsed_source = ?parsed,
                "Generating exact CFF at its term-local capacity"
            );
            three_dimensional_reps::generate_3d_expression(preparation, source_options).map_err(
                |error| {
                    eyre::eyre!(
                        "generalized CFF expression generation failed for exact 4D source in graph `{}` with physical EMR bounds {:?} and term-local exact-occurrence bounds {:?}: {error}\n{}",
                        self.name,
                        physical_energy_degree_bounds,
                        source_options.energy_degree_bounds,
                        three_d_source_summary(parsed),
                    )
                },
            )
        };
        type Candidate = (
            usize,
            Vec<usize>,
            Arc<PlannedExactSourceNumerator>,
            Option<GeneratedThreeDExpression>,
        );
        let mut selected: Option<Candidate> = None;
        let mut pending = VecDeque::from(energy_assignment_plans.clone());
        let mut seen_bounds = Vec::new();
        let mut challenged = false;
        loop {
            if pending.is_empty() && !challenged {
                challenged = true;
                if let Some((_, _, best, _)) = &selected
                    && let Some(challenger) = best
                        .binding
                        .binding
                        .assignment
                        .placement_challenger(&seen_bounds)?
                {
                    pending.push_back(Arc::new(challenger));
                }
            }
            let Some(plan) = pending.pop_front() else {
                break;
            };
            let proposal = seen_bounds.len();
            let certificate_started = std::time::Instant::now();
            exact_source_energy_mapper.certify_assignment(&plan)?;
            certificate_time += certificate_started.elapsed();
            // Preparing the immutable mapper also certifies its exact signed
            // diagonal, including affine shifts, before any native CFF trial.
            // A complete binding cache hit reuses that same certified template.
            let (numerator, build_time, cache_time) = PlannedExactSourceNumerator::prepare(
                Arc::clone(exact_source_energy_mapper),
                Arc::clone(&plan),
                context.as_deref_mut(),
            )?;
            template_time += build_time + cache_time;
            template_build_times.insert(plan.energy_degree_bounds().to_vec(), build_time);
            seen_bounds.push(plan.energy_degree_bounds().to_vec());
            let mut rank_envelope = plan
                .energy_degree_bounds()
                .iter()
                .map(|(_, degree)| *degree)
                .collect::<Vec<_>>();
            rank_envelope.sort_unstable_by(|left, right| right.cmp(left));
            source_options.energy_degree_bounds = Some(plan.energy_degree_bounds().to_vec());
            let key =
                ExactCffGenerationCache::generation_key(parsed, energy_edges, &source_options);
            let known_count = context
                .as_deref_mut()
                .map(|context| &mut context.generation_cache)
                .and_then(|cache| cache.count(&key));
            let started = std::time::Instant::now();
            // A known contender needs no expression until it wins. Generated
            // expressions clone their tree containers; count-only loser records
            // avoid retaining or rebuilding those trees on repeated requests.
            let generated = if known_count.is_some() {
                None
            } else {
                if let Some(cache) = context
                    .as_deref_mut()
                    .map(|context| &mut context.generation_cache)
                {
                    cache.native_generations += 1;
                }
                let native_started = std::time::Instant::now();
                let generated = generate(&source_options)?;
                native_times.insert(
                    plan.energy_degree_bounds().to_vec(),
                    native_started.elapsed(),
                );
                Some(generated)
            };
            let map_count = known_count.unwrap_or_else(|| {
                generated
                    .as_ref()
                    .expect("an unknown candidate was freshly generated")
                    .expression
                    .orientations
                    .len()
            });
            if known_count.is_none()
                && let Some(cache) = context
                    .as_deref_mut()
                    .map(|context| &mut context.generation_cache)
            {
                cache.record_count(key, map_count);
            }
            crate::debug_tags!(#generation, #uv, #local, #four_d, #cff, #profile;
                graph = %self.name,
                proposal,
                native_source_maps = map_count,
                bounds = ?source_options.energy_degree_bounds,
                count_memo_hit = known_count.is_some(),
                elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
                "Scored bounded exact-CFF assignment proposal"
            );
            // Native rows are primary; the descending envelope starts with
            // maximum rank and supplies the remaining deterministic tie break.
            // Equal scores retain the earlier admitted candidate. This bounded
            // search does not claim a global minimum.
            if selected
                .as_ref()
                .is_none_or(|(best_count, best_envelope, _, _)| {
                    (map_count, &rank_envelope) < (*best_count, best_envelope)
                })
            {
                selected = Some((map_count, rank_envelope, numerator, generated));
            }
        }
        let (selected_count, _, numerator, generated) =
            selected.expect("rank planning always provides a baseline assignment");
        let energy_assignment_plan = Arc::clone(&numerator.binding.binding.assignment);
        source_options.energy_degree_bounds =
            Some(energy_assignment_plan.energy_degree_bounds().to_vec());
        let key = ExactCffGenerationCache::generation_key(parsed, energy_edges, &source_options);
        let cached = context
            .as_deref_mut()
            .map(|context| &mut context.generation_cache)
            .and_then(|cache| cache.entries.get(&key))
            .and_then(|entry| entry.payload.as_ref());
        let cache_hit = cached.is_some();
        let generated = if let Some(generated) = generated {
            generated
        } else if let Some(cached) = cached {
            cached.clone()
        } else {
            // A previous losing key can win against a different proposal set.
            // Its count was enough to select it; obtain its payload only now.
            if let Some(cache) = context
                .as_deref_mut()
                .map(|context| &mut context.generation_cache)
            {
                cache.native_generations += 1;
            }
            let native_started = std::time::Instant::now();
            let generated = generate(&source_options)?;
            *native_times
                .entry(energy_assignment_plan.energy_degree_bounds().to_vec())
                .or_default() += native_started.elapsed();
            generated
        };
        debug_assert_eq!(generated.expression.orientations.len(), selected_count);
        if !cache_hit && let Some(context) = context {
            let cache = &mut context.generation_cache;
            // Reuse requires both canonical topology and identical occurrence
            // capacity; keep the term's assignment plan unchanged. Retain
            // only the winning native payload. Losing counts and certified
            // templates have their own bounded retention.
            cache.record_payload(key, &generated);
        }
        let winning_native_time = native_times
            .get(energy_assignment_plan.energy_degree_bounds())
            .copied()
            .unwrap_or_default();
        let native_time = native_times.values().copied().sum::<std::time::Duration>();
        let winning_template_build_time = template_build_times
            .get(energy_assignment_plan.energy_degree_bounds())
            .copied()
            .unwrap_or_default();
        let template_build_time = template_build_times.values().copied().sum::<Duration>();
        let selection_time = selection_started.elapsed();
        crate::debug_tags!(#generation, #uv, #local, #four_d, #cff, #profile;
            stage = "exact_cff_assignment_selection",
            graph = %self.name,
            admitted_candidates = seen_bounds.len(),
            native_rows = selected_count,
            candidate_certificate_ms = certificate_time.as_secs_f64() * 1000.0,
            winning_native_ms = winning_native_time.as_secs_f64() * 1000.0,
            losing_native_ms = (native_time - winning_native_time).as_secs_f64() * 1000.0,
            winning_template_build_ms = winning_template_build_time.as_secs_f64() * 1000.0,
            losing_template_build_ms = (template_build_time - winning_template_build_time).as_secs_f64() * 1000.0,
            selection_overhead_ms = selection_time.saturating_sub(winning_native_time).saturating_sub(certificate_time).saturating_sub(template_time).as_secs_f64() * 1000.0,
            native_payload_cache_hit = cache_hit,
            "Selected a certified exact CFF assignment; losing trials count as selection overhead"
        );
        let energy_degree_bound_report = CffEnergyDegreeBoundReport {
            source_kind: CffEnergyBoundSourceKind::ExactFourD,
            physical_parent_bounds: physical_energy_degree_bounds.clone(),
            assigned_cff_source_bounds: energy_assignment_plan.energy_degree_bounds().to_vec(),
        };
        Ok((
            generated,
            numerator,
            energy_assignment_plan,
            energy_degree_bound_report,
        ))
    }

    pub(crate) fn generate_3d_expression_for_integrand(
||||||| qrspmutx 7a048a39 "add feynkit ufo and Python bindings" (parents of rebased revision)
use itertools::Itertools;
use linnet::half_edge::{
    HedgeGraph,
    involution::{EdgeVec, HedgePair},
    subgraph::{OrientedCut, SubGraphLike, SubSetOps},
};
use linnet::half_edge::{
    involution::{EdgeIndex, Orientation},
    subgraph::InternalSubGraph,
};
use symbolica::{
    atom::{Atom, AtomCore},
    id::{Pattern, Replacement},
};
use typed_index_collections::TiVec;

use serde::{Deserialize, Serialize};

use tracing::debug;

use super::{
    cff_graph::CFFGenerationGraph,
    esurface::{Esurface, EsurfaceCollection, EsurfaceID, ExternalShift},
    expression::{CFFExpression, OrientationID},
    hsurface::HsurfaceCollection,
    surface::{HybridSurfaceRef, UnitSurface},
};

#[derive(Debug, Clone)]
struct GenerationData {
    graph: CFFGenerationGraph,
    surface_id: Option<HybridSurfaceID>,
}

#[derive(Debug, Clone)]
pub struct ShiftRewrite {
    pub dependent_momentum: EdgeIndex,
    pub dependent_momentum_expr: ExternalShift,
}

fn forget_graphs(data: GenerationData) -> HybridSurfaceID {
    data.surface_id.expect("corrupted expression tree")
}

impl GenerationData {
    fn insert_esurface(&mut self, surface_id: HybridSurfaceID) {
        self.surface_id = Some(surface_id);
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[allow(dead_code)]
struct CFFTreeNodePointer {
    term_id: usize,
    node_id: usize,
}

// Orientation is a bitset that represents the orientation of the edges in a graph.
// 0 means +, 1 means -, in this representation the original graph is represented by the number 0
#[derive(Debug, Clone, Copy)]
struct OrientationGenerator {
    identifier: usize,
    num_edges: usize,
}

impl OrientationGenerator {
    #[allow(unused)]
    fn default(num_edges: usize) -> Self {
        Self {
            identifier: 0,
            num_edges,
        }
    }
}

impl IntoIterator for OrientationGenerator {
    type Item = Orientation;
    type IntoIter = OrientationIterator;

    fn into_iter(self) -> Self::IntoIter {
        OrientationIterator {
            identifier: self.identifier,
            current_location: 0,
            num_edges: self.num_edges,
        }
    }
}

// OrientationIterator allows us to iterate over the edges in a graph, and
// view their orientation as a boolean
struct OrientationIterator {
    identifier: usize,
    current_location: usize,
    num_edges: usize,
}

impl Iterator for OrientationIterator {
    type Item = Orientation;
    fn next(&mut self) -> Option<Self::Item> {
        if self.current_location < self.num_edges {
            let result_bool = self.identifier & (1 << self.current_location) == 0;
            let result = match result_bool {
                true => Orientation::Default,
                false => Orientation::Reversed,
            };

            self.current_location += 1;
            Some(result)
        } else {
            None
        }
    }
}

// This function returns an iterator over all possible orientations of a graph
fn iterate_possible_orientations(num_edges: usize) -> impl Iterator<Item = OrientationGenerator> {
    if num_edges > 64 {
        panic!("Maximum number of edges supported is currently 64")
    }

    let max_size = 2_usize.pow(num_edges as u32);
    (0..max_size).map(move |x| OrientationGenerator {
        identifier: x,
        num_edges,
    })
}

#[cfg(test)]
fn get_orientations<E, V, H>(
    graph: &HedgeGraph<E, V, H>,
    dummy_edges: &[EdgeIndex],
) -> Vec<CFFGenerationGraph> {
    let internal_subgraph = InternalSubGraph::cleaned_filter_pessimist(graph.full_filter(), graph);
    let num_virtual_edges = graph.count_internal_edges(&internal_subgraph);
    let virtual_possible_orientations = iterate_possible_orientations(num_virtual_edges);

    virtual_possible_orientations
        .map(|orientation_of_virtuals| {
            let mut orientation_of_virtuals = orientation_of_virtuals.into_iter();

            let global_orientation = graph.new_edgevec(|_, __, hedge_pair| match hedge_pair {
                HedgePair::Unpaired { .. } => Orientation::Default,
                HedgePair::Paired { .. } => orientation_of_virtuals
                    .next()
                    .expect(" unable to reconstruct orientation"),
                HedgePair::Split { .. } => todo!(),
            });

            assert!(
                orientation_of_virtuals.next().is_none(),
                "did not saturate virtual orientations when constructing global orientation"
            );

            CFFGenerationGraph::new(graph, global_orientation, dummy_edges)
        })
        .collect_vec()
}

pub(crate) fn get_orientations_from_subgraph<E, V, H, S: SubGraphLike>(
    graph: &HedgeGraph<E, V, H>,
    subgraph: &S,
    reversed_dangling: &[EdgeIndex],
) -> Vec<CFFGenerationGraph> {
    let num_virtual_edges = graph.count_internal_edges(subgraph);
    let virtual_possible_orientations = iterate_possible_orientations(num_virtual_edges);

    virtual_possible_orientations
        .map(|orientation_of_virtuals| {
            let mut orientation_of_virtuals = orientation_of_virtuals.into_iter();

            let global_orientation = graph.new_edgevec(|_, edge_id, _| {
                if let Some((pair, _, _)) = graph
                    .iter_edges_of(subgraph)
                    .find(|(_pair, id, _)| *id == edge_id)
                {
                    match pair {
                        HedgePair::Paired { .. } => orientation_of_virtuals
                            .next()
                            .expect("orientation generation corrupted, not enough edges"),
                        HedgePair::Unpaired { .. } => Orientation::Default,
                        HedgePair::Split { .. } => {
                            if reversed_dangling.contains(&edge_id) {
                                Orientation::Reversed
                            } else {
                                Orientation::Default
                            }
                        }
                    }
                } else {
                    Orientation::Undirected
                }
            });

            CFFGenerationGraph::new_from_subgraph(graph, global_orientation, subgraph).unwrap()
        })
        .filter(|cff_graph| !cff_graph.has_directed_cycle_initial())
        .collect()
}

#[allow(unused)]
fn get_orientations_with_cut<E, V, H>(
    graph: &HedgeGraph<E, V, H>,
    oriented_cut: &OrientedCut,
) -> Vec<EdgeVec<Orientation>> {
    let internal_subgraph = InternalSubGraph::cleaned_filter_pessimist(graph.full_filter(), graph);
    let num_virtual_edges = graph.count_internal_edges(&internal_subgraph);

    let virtual_possible_orientations = iterate_possible_orientations(num_virtual_edges);

    let orientations_consistent_with_cut = virtual_possible_orientations
        .map(|orientation_of_virtuals| {
            // pad a virtual orientation with orientations of externals.
            let mut orientation_of_virtuals = orientation_of_virtuals.into_iter();

            let global_orientation = graph.new_edgevec(|_, __, hedge_pair| match hedge_pair {
                HedgePair::Unpaired { .. } => Orientation::Default,
                HedgePair::Paired { .. } => orientation_of_virtuals
                    .next()
                    .expect(" unable to reconstruct orientation"),
                HedgePair::Split { .. } => todo!(),
            });

            assert!(
                orientation_of_virtuals.next().is_none(),
                "did not saturate virtual orientations when constructing global orientation"
            );

            global_orientation
        })
        .filter(|global_orientation| {
            // filter out orientations that are not consistent with the cut
            let edges_in_cut = graph.iter_edges_of(oriented_cut).map(|(_, id, _)| id);
            let orientation_of_edges_in_cut = oriented_cut.iter_edges(graph).map(|(or, _)| or);

            edges_in_cut
                .zip(orientation_of_edges_in_cut)
                .all(|(edge_id, orientation)| global_orientation[edge_id] == orientation)
        })
        .filter(|global_orientation| {
            // filter out orientations that have a directed cycle
            let graph = CFFGenerationGraph::new(graph, global_orientation.clone(), &[]);
            !graph.has_directed_cycle_initial()
        });

    orientations_consistent_with_cut.collect_vec()
}

#[cfg(test)]
fn generate_cff_expression<E, V, H>(
    graph: &HedgeGraph<E, V, H>,
    canonize_esurface: &Option<ShiftRewrite>,
    edges_in_initial_state_cut: &[EdgeIndex],
    dummy_edges: &[EdgeIndex],
) -> Result<CFFExpression<OrientationID>> {
    let graphs = get_orientations(graph, dummy_edges);
    debug!("number of orientations: {}", graphs.len());
    let mut surface_cache = SurfaceCache {
        esurface_cache: EsurfaceCollection::from_iter(std::iter::empty()),
        hsurface_cache: HsurfaceCollection::from_iter(std::iter::empty()),
    };
    let graph_cff = generate_cff_from_orientations(
        graphs,
        &mut surface_cache,
        edges_in_initial_state_cut,
        canonize_esurface,
    )?;

    // patch the surface cache
    Ok(graph_cff)
}

impl Graph {
    pub(crate) fn generate_cff(
=======
use feynkit_cff::{CffResult, ShiftRewrite, SurfaceCache};
use linnet::half_edge::involution::EdgeIndex;
use symbolica::{
    atom::{Atom, AtomCore},
    id::Replacement,
};

use crate::{graph::LoopMomentumBasis, settings::global::OrientationPattern};

/// GammaLoop entry points layered on the canonical FeynKit CFF generator.
pub(crate) trait GammaLoopGraphCffExt {
    fn generate_cff(
>>>>>>> xzuonylm 3e5f0ddf (rebased revision)
        &mut self,
        contract_edges: &[EdgeIndex],
        canonize_esurface: &Option<ShiftRewrite>,
<<<<<<< qrspmutx ee257326 "add feynkit ufo and Python bindings" (rebase destination)
        options: &Generate3DExpressionOptions,
        analysis_numerator: Option<&Atom>,
    ) -> Result<GeneratedThreeDExpression<Esurface, Hsurface>> {
        let generated = self.generate_raw_3d_expression_for_integrand(
            contract_edges,
            options,
            analysis_numerator,
        )?;
        let initial_state_cut_edges = self
            .iter_edges_of(&self.initial_state_cut)
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        self.convert_generated_expression_surfaces(
            generated,
            canonize_esurface,
            &initial_state_cut_edges,
||||||| qrspmutx 7a048a39 "add feynkit ufo and Python bindings" (parents of rebased revision)
        orientation_pattern: &OrientationPattern,
    ) -> Result<CFFExpression<OrientationID>> {
        let mut seed_graph = CFFGenerationGraph::new_from_graph(self);

        for edge in contract_edges {
            seed_graph = seed_graph.contract_edge(*edge);
        }

        let edges_in_initial_state_cut = self
            .iter_edges_of(&self.initial_state_cut)
            .map(|x| x.1)
            .collect_vec();

        let virtual_edges_of_contracted_graph = seed_graph.num_virtual_edges();

        let orientations = iterate_possible_orientations(virtual_edges_of_contracted_graph);

        let mut oriented_acyclic_graphs = vec![];

        for orientation in orientations {
            let mut orientation_iterator = orientation.into_iter();

            let global_orientation = self.new_edgevec(|_, edge_id, hedge_pair| {
                if hedge_pair.is_unpaired() || contract_edges.contains(&edge_id) {
                    Orientation::Undirected
                } else if edges_in_initial_state_cut.contains(&edge_id) {
                    Orientation::Default
                } else {
                    orientation_iterator
                        .next()
                        .expect("orientation generation corrupted, not enough edges")
                }
            });

            if orientation_pattern.filter(&global_orientation) {
                let mut cff_graph = seed_graph.clone();
                cff_graph.apply_orientation(global_orientation)?;

                if !cff_graph.has_directed_cycle_initial() {
                    oriented_acyclic_graphs.push(cff_graph);
                }
            }
        }

        generate_cff_from_orientations(
            oriented_acyclic_graphs,
            &mut self.surface_cache,
            &edges_in_initial_state_cut,
            canonize_esurface,
=======
        orientation_pattern: &OrientationPattern,
    ) -> Result<CffResult>;
}

/// GammaLoop symbolic lowering over the canonical FeynKit surface arena.
pub trait GammaLoopSurfaceCacheExt {
    fn substitute_energies(&self, atom: &Atom, cut_edges: &[EdgeIndex]) -> Atom;
    fn get_all_replacements(&self, cut_edges: &[EdgeIndex]) -> Vec<Replacement>;
    fn get_all_replacements_in_lmb(
        &self,
        cut_edges: &[EdgeIndex],
        lmb: &LoopMomentumBasis,
    ) -> Vec<Replacement>;
}

impl GammaLoopSurfaceCacheExt for SurfaceCache {
    fn substitute_energies(&self, atom: &Atom, cut_edges: &[EdgeIndex]) -> Atom {
        atom.replace_multiple(self.get_all_replacements(cut_edges))
    }

    fn get_all_replacements(&self, cut_edges: &[EdgeIndex]) -> Vec<Replacement> {
        self.replacements_with(
            |edge| {
                if cut_edges.contains(&edge) {
                    crate::utils::cut_energy(edge)
                } else {
                    crate::utils::ose_atom_from_index(edge)
                }
            },
            crate::utils::external_energy_atom_from_index,
>>>>>>> xzuonylm 3e5f0ddf (rebased revision)
        )
    }
<<<<<<< qrspmutx ee257326 "add feynkit ufo and Python bindings" (rebase destination)

    /// Generate source-owned maps before interning any physical graph surface.
    /// Bounded dispatch can discard a raw proposal without leaving graph state.
    pub(crate) fn generate_raw_3d_expression_for_integrand(
        &self,
        contract_edges: &[EdgeIndex],
        options: &Generate3DExpressionOptions,
        analysis_numerator: Option<&Atom>,
    ) -> Result<GeneratedThreeDExpression> {
        let request_started = Instant::now();
        let initial_state_cut_edges = self
            .iter_edges_of(&self.initial_state_cut)
            .map(|(_, edge_id, _)| edge_id)
            .collect_vec();
        let source_contract_edges =
            source_contract_edges_for_3d_expression(contract_edges, &initial_state_cut_edges);
        if !source_contract_edges.is_empty() {
            debug!(
                "contracting {} internal edges before generalized 3D expression generation",
                source_contract_edges.len()
            );
        }
        let reconstruction_started = Instant::now();
        let source = GraphThreeDSource::new(self, &source_contract_edges)?;
        let reconstruction_time = reconstruction_started.elapsed();
        let bridge_edges = self
            .iter_edges_of(&self.tree_edges)
            .map(|(_, edge_id, _)| edge_id)
            .filter(|edge_id| !initial_state_cut_edges.contains(edge_id))
            .collect::<HashSet<_>>();
        for edge_id in source_contract_edges
            .iter()
            .filter(|edge_id| bridge_edges.contains(edge_id))
        {
            let Some(coordinates) = source.reconstructible_outer_loop_coordinates(*edge_id) else {
                return Err(eyre::eyre!(
                    "contracted bridge edge {} retains an inner-loop energy coordinate and cannot be projected independently of the CFF source",
                    usize::from(*edge_id),
                ));
            };
            if coordinates.iter().any(|coordinate| *coordinate != 0) {
                return Err(eyre::eyre!(
                    "contracted bridge edge {} has nonzero outer-loop coordinates {:?}; a CFF-external tree edge must have a purely external affine energy",
                    usize::from(*edge_id),
                    coordinates,
                ));
            }
        }
        let mut source_options = options.clone();
        if let Some(numerator) = analysis_numerator {
            source_options.energy_degree_bounds = Some(
                self.automatic_numerator_energy_degree_bounds_in_atoms_excluding_with_min_degree(
                    [numerator],
                    initial_state_cut_edges.iter().chain(&bridge_edges).copied(),
                    1,
                )?,
            );
            debug!(
                graph = %self.name,
                bounds = ?source_options.energy_degree_bounds,
                "using source-edge numerator energy-degree bounds"
            );
        }
        if let Some(bounds) = &mut source_options.energy_degree_bounds {
            // Cut aliases and tree denominators are external to the CFF graph,
            // but their energies remain in the complete numerator. Removing
            // their EMR bounds here therefore changes only CFF capacity; in
            // particular, an explicitly bounded numerator class remains
            // `Some([])` when every bound belongs to this external sector.
            bounds.retain(|(edge_id, _)| {
                let edge_id = EdgeIndex(*edge_id);
                !initial_state_cut_edges.contains(&edge_id) && !bridge_edges.contains(&edge_id)
            });
        }
        for &(edge_id, degree) in source_options
            .energy_degree_bounds
            .as_deref()
            .unwrap_or(&[])
        {
            let edge_id = EdgeIndex(edge_id);
            if !source_contract_edges.contains(&edge_id) {
                continue;
            }
            if source
                .reconstructible_outer_loop_coordinates(edge_id)
                .is_some()
            {
                return Err(eyre::eyre!(
                    "numerator energy-degree bound {degree} for shrunken EMR edge {} cannot yet be supplied to generalized CFF generation: its exact outer affine energy is retained for numerator evaluation, but the standalone generator has no numerator-only/shrunken-edge bound channel",
                    usize::from(edge_id),
                ));
            }
            return Err(eyre::eyre!(
                "numerator energy-degree bound {degree} for shrunken EMR edge {} retains contracted inner-loop dependence and cannot be represented in the outer CFF source",
                usize::from(edge_id),
            ));
        }
        let preparation_time = request_started
            .elapsed()
            .saturating_sub(reconstruction_time);
        let native_started = Instant::now();
        let mut generated = {
            let result = three_dimensional_reps::generate_3d_expression(&source, &source_options);
            result.map_err(|error| {
                    let source_summary = source
                        .to_three_d_parsed_graph()
                        .map(|parsed| three_d_source_summary(&parsed))
                        .unwrap_or_else(|source_error| {
                            format!("failed to rebuild 3D source summary: {source_error}")
                    });
                    eyre::eyre!(
                        "generalized CFF expression generation failed for graph `{}` with source-edge numerator energy-degree bounds {:?}: {error}\n{source_summary}",
                        self.name,
                        source_options.energy_degree_bounds,
                    )
                })
        }?;
        let native_time = native_started.elapsed();

        // Generic edge-index remapping leaves omitted parent edges at zero. A
        // contracted edge can still have an exact outer-energy map when its
        // source coordinates contain no inner-loop component; restore only
        // that source-owned information and leave genuinely inner-dependent
        // contracted edges unavailable.
        for edge_id in &source_contract_edges {
            let Some(coordinates) = source.reconstructible_outer_loop_coordinates(*edge_id) else {
                continue;
            };
            let signature = &self.loop_momentum_basis.edge_signatures[*edge_id];
            for orientation in generated.expression.orientations.iter_mut() {
                if coordinates.len() != orientation.loop_energy_map.len() {
                    return Err(eyre::eyre!(
                        "contracted edge {} has {} outer coordinates for {} generated loop-energy maps",
                        usize::from(*edge_id),
                        coordinates.len(),
                        orientation.loop_energy_map.len(),
                    ));
                }
                let mut edge_energy = coordinates.iter().zip(&orientation.loop_energy_map).fold(
                    LinearEnergyExpr::zero(),
                    |sum, (coefficient, loop_energy)| {
                        sum + loop_energy.clone().scale_rational(coefficient.clone())
                    },
                );
                for (external_edge, sign) in self
                    .loop_momentum_basis
                    .ext_edges
                    .iter()
                    .zip(&signature.external)
                {
                    let coefficient = match sign {
                        SignOrZero::Zero => 0,
                        SignOrZero::Plus => 1,
                        SignOrZero::Minus => -1,
                    };
                    edge_energy =
                        edge_energy + LinearEnergyExpr::external(*external_edge, coefficient);
                }
                let edge_energy_slot = orientation
                    .edge_energy_map
                    .get_mut(usize::from(*edge_id))
                    .ok_or_else(|| {
                    eyre::eyre!(
                        "contracted edge {} is outside the generated parent energy map",
                        usize::from(*edge_id),
                    )
                })?;
                *edge_energy_slot = edge_energy;
            }
        }

        let elapsed = request_started.elapsed();
        crate::debug_tags!(#generation, #cff, #profile;
            stage = "raw_cff_generation",
            graph = %self.name,
            native_source_maps = generated.expression.orientations.len(),
            elapsed_ms = elapsed.as_secs_f64() * 1000.0,
            source_reconstruction_ms = reconstruction_time.as_secs_f64() * 1000.0,
            preparation_ms = preparation_time.as_secs_f64() * 1000.0,
            native_generation_ms = native_time.as_secs_f64() * 1000.0,
            postprocessing_ms = elapsed.saturating_sub(reconstruction_time + preparation_time + native_time).as_secs_f64() * 1000.0,
            "Generated a raw CFF request before surface conversion"
        );
        Ok(generated)
    }

    pub(crate) fn production_cff_3d_expression_options(
        &self,
        settings: &GenerationSettings,
    ) -> Result<Generate3DExpressionOptions> {
        self.cff_3d_expression_options(numerator_sampling_scale_mode(
            settings.uniform_numerator_sampling_scale,
        ))
    }

    pub fn cff_3d_expression_options(
        &self,
        numerator_sampling_scale: NumeratorSamplingScaleMode,
    ) -> Result<Generate3DExpressionOptions> {
        let numerator = self.production_numerator_atom_for_full_3d_expression();
        let cff_external_edges = self
            .iter_edges_of(&self.initial_state_cut)
            .chain(self.iter_edges_of(&self.tree_edges))
            .map(|(_, edge_id, _)| edge_id);
        // Analyze the factorized numerator compositionally while keeping every
        // active energy degree attached to its EMR source edge. Initial-cut
        // aliases and structural bridges stay in their separately projected
        // external sector. The standalone CFF generator remaps active source
        // IDs only after extracting its local graph, so a loop-momentum basis
        // never defines numerator ownership.
        let energy_degree_bounds = self
            .automatic_numerator_energy_degree_bounds_in_atoms_excluding_with_min_degree(
                [&numerator],
                cff_external_edges,
                1,
            )?;
        debug!(
            graph = %self.name,
            bounds = ?energy_degree_bounds,
            "using production CFF source-edge numerator energy-degree bounds"
        );
        Ok(Generate3DExpressionOptions {
            representation: RepresentationMode::Cff,
            cff_generation_context: CffGenerationContext::Standalone,
            energy_degree_bounds: Some(energy_degree_bounds),
            numerator_sampling_scale,
            preserve_internal_edges_as_four_d_denominators: Vec::new(),
        })
    }

    #[cfg(test)]
    pub(crate) fn denominator_only_cff_3d_expression_options(&self) -> Generate3DExpressionOptions {
        Generate3DExpressionOptions {
            representation: RepresentationMode::Cff,
            cff_generation_context: CffGenerationContext::Standalone,
            energy_degree_bounds: Some(Vec::new()),
            numerator_sampling_scale: NumeratorSamplingScaleMode::None,
            preserve_internal_edges_as_four_d_denominators: Vec::new(),
        }
    }

    pub(crate) fn convert_generated_expression_surfaces(
        &mut self,
        generated: GeneratedThreeDExpression,
        canonize_esurface: &Option<ShiftRewrite>,
        initial_state_cut_edges: &[EdgeIndex],
    ) -> Result<GeneratedThreeDExpression<Esurface, Hsurface>> {
        Ok(self
            .convert_generated_expression_surfaces_impl(
                generated,
                canonize_esurface,
                initial_state_cut_edges,
                None,
            )?
            .0)
    }

    pub(crate) fn convert_4d_expression_surfaces(
        &mut self,
        generated: GeneratedThreeDExpression,
        physical_surfaces: &[Option<LinearSurface>],
    ) -> Result<(
        GeneratedThreeDExpression<Esurface, Hsurface>,
        BTreeMap<EsurfaceID, CffEnergyFactorOwnership>,
    )> {
        let canonize_esurface = self.get_esurface_canonization(&self.loop_momentum_basis);
        let initial_state_cut_edges = self
            .iter_edges_of(&self.initial_state_cut)
            .map(|(_, edge_id, _)| edge_id)
            .collect_vec();
        self.convert_generated_expression_surfaces_impl(
            generated,
            &canonize_esurface,
            &initial_state_cut_edges,
            Some(physical_surfaces),
        )
    }

    fn convert_generated_expression_surfaces_impl(
        &mut self,
        generated: GeneratedThreeDExpression,
        canonize_esurface: &Option<ShiftRewrite>,
        initial_state_cut_edges: &[EdgeIndex],
        physical_surfaces: Option<&[Option<LinearSurface>]>,
    ) -> Result<(
        GeneratedThreeDExpression<Esurface, Hsurface>,
        BTreeMap<EsurfaceID, CffEnergyFactorOwnership>,
    )> {
        let GeneratedThreeDExpression {
            mut expression,
            energy_factor_ownership,
            energy_factor_components,
            source_energy_degree_bounds,
            denominator_only_global_prefactor_sign,
            core_global_prefactor_sign,
        } = generated;
        if !expression.residual_denominators.is_empty() {
            return Err(eyre::eyre!(
                "GammaLoop production CFF conversion received residual four-dimensional denominator edges [{}]; production callers must contract tree edges and attach their projected denominators separately",
                expression
                    .residual_denominators
                    .iter()
                    .map(|denominator| denominator.edge_id.0)
                    .join(", "),
            ));
        }
        let mut linear_surface_components = BTreeMap::new();
        for (surface_id, surface) in expression.surfaces.linear_surface_cache.iter_enumerated() {
            let mut surface_component = None;
            for (edge_id, _) in &surface.expression.internal_terms {
                for (component_index, component) in energy_factor_components.iter().enumerate() {
                    if !component.internal_edge_ids.contains(&usize::from(*edge_id)) {
                        continue;
                    }
                    if let Some(previous) = surface_component
                        && previous != component_index
                    {
                        return Err(eyre::eyre!(
                            "generated causal surface {surface_id:?} spans disconnected energy-factor components"
                        ));
                    }
                    surface_component = Some(component_index);
                }
            }
            if let Some(component_index) = surface_component {
                linear_surface_components.insert(surface_id, component_index);
            }
        }
        let mut linear_surface_map = BTreeMap::<LinearSurfaceID, SurfaceMapEntry>::new();
        let mut retained_linear_surfaces = Vec::new();
        for (linear_surface_id, surface) in
            expression.surfaces.linear_surface_cache.iter_enumerated()
        {
            let physical = physical_surfaces
                .and_then(|surfaces| surfaces.get(usize::from(linear_surface_id)))
                .and_then(Option::as_ref);
            if physical_surfaces.is_some() && physical.is_none() {
                retained_linear_surfaces.push((linear_surface_id, surface.clone()));
            } else {
                let converted = self.intern_generated_linear_surface(
                    physical.unwrap_or(surface),
                    canonize_esurface,
                    initial_state_cut_edges,
                )?;
                linear_surface_map.insert(linear_surface_id, converted);
            }
        }
        let mut surface_cache = self.surface_cache.clone();
        for (source_id, surface) in retained_linear_surfaces {
            let target_id = LinearSurfaceID(surface_cache.linear_surface_cache.len());
            surface_cache.linear_surface_cache.push(surface);
            linear_surface_map.insert(
                source_id,
                SurfaceMapEntry {
                    surface_id: HybridSurfaceID::Linear(target_id),
                    sign: 1,
                },
            );
        }
        let mut mapped_surface_components = BTreeMap::<EsurfaceID, usize>::new();
        for (source_id, component_index) in linear_surface_components {
            let Some(SurfaceMapEntry {
                surface_id: HybridSurfaceID::Esurface(target_id),
                ..
            }) = linear_surface_map.get(&source_id)
            else {
                continue;
            };
            if let Some(previous) = mapped_surface_components.insert(*target_id, component_index)
                && previous != component_index
            {
                return Err(eyre::eyre!(
                    "physical causal surface {target_id:?} merges disconnected energy-factor components"
                ));
            }
        }
        let mapped_surface_ownership = mapped_surface_components
            .into_iter()
            .map(|(surface_id, component_index)| {
                (
                    surface_id,
                    energy_factor_components[component_index].ownership,
                )
            })
            .collect();
        for orientation in expression.orientations.iter_mut() {
            let mut remapped_variants = Vec::with_capacity(orientation.variants.len());
            for mut variant in std::mem::take(&mut orientation.variants) {
                let signed_denominators =
                    remap_denominator_tree_surface_ids(&variant.denominator, &linear_surface_map);
                for surface_id in &mut variant.numerator_surfaces {
                    let numerator_sign =
                        remap_generated_surface_id(surface_id, &linear_surface_map);
                    if numerator_sign < 0 {
                        variant.prefactor = -variant.prefactor;
                    }
                }
                match energy_factor_ownership {
                    CffEnergyFactorOwnership::GlobalSourceProduct => {
                        // GammaLoop's CFF evaluator convention keeps the
                        // on-shell-energy factors as one global product
                        // 1/prod(-2E_i) for the ordinary CFF denominator sector.
                        variant.half_edges.clear();
                    }
                    // Generalized variants retain their local positive-energy
                    // factors. Their branch-relative denominator convention is
                    // already consumed while each powered channel is lowered.
                    CffEnergyFactorOwnership::VariantLocal => {}
                }
                for remapped_denominator in signed_denominators {
                    let mut signed_variant = variant.clone();
                    signed_variant.denominator = remapped_denominator.denominator;
                    signed_variant.denominator_surface_signs = remapped_denominator.surface_signs;
                    signed_variant.prefactor *= Rational::from(remapped_denominator.prefactor);
                    remapped_variants.push(signed_variant);
                }
            }
            orientation.variants = remapped_variants;
        }

        Ok((
            GeneratedThreeDExpression {
                expression: CFFExpression {
                    orientations: expression.orientations,
                    surfaces: surface_cache,
                    residual_denominators: Vec::new(),
                },
                energy_factor_ownership,
                energy_factor_components,
                source_energy_degree_bounds,
                denominator_only_global_prefactor_sign,
                core_global_prefactor_sign,
            },
            mapped_surface_ownership,
        ))
    }

    fn intern_generated_linear_surface(
        &mut self,
        surface: &LinearSurface,
        canonize_esurface: &Option<ShiftRewrite>,
        initial_state_cut_edges: &[EdgeIndex],
    ) -> Result<SurfaceMapEntry> {
        if surface.numerator_only
            && (!surface.expression.uniform_scale_coeff.is_zero()
                || !surface.expression.constant.is_zero())
        {
            let surface_id =
                Into::<LinearSurfaceID>::into(self.surface_cache.linear_surface_cache.len());
            self.surface_cache
                .linear_surface_cache
                .push(surface.clone());
            return Ok(SurfaceMapEntry {
                surface_id: HybridSurfaceID::Linear(surface_id),
                sign: 1,
            });
        }

        if !surface.expression.uniform_scale_coeff.is_zero()
            || !surface.expression.constant.is_zero()
        {
            return Err(eyre::eyre!(
                "generalized CFF production cannot convert non-homogeneous linear surface {:?}",
                surface
            ));
        }

        let mut positive_energies = Vec::new();
        let mut negative_energies = Vec::new();
        let mut external_shift = Vec::new();
        collect_linear_surface_terms(
            &surface.expression.internal_terms,
            initial_state_cut_edges,
            &mut positive_energies,
            &mut negative_energies,
            &mut external_shift,
        )?;
        for (edge_id, coeff) in &surface.expression.external_terms {
            let coeff = integer_coeff(coeff)?;
            external_shift.push((*edge_id, coeff));
        }
        positive_energies.sort();
        negative_energies.sort();
        external_shift.sort_by_key(|(edge_id, _)| *edge_id);

        let (surface_id, sign) = match surface.kind {
            LinearSurfaceKind::Esurface => {
                let mut sign = 1;
                if positive_energies.is_empty() && !negative_energies.is_empty() {
                    sign = -1;
                    positive_energies = negative_energies;
                    negative_energies = Vec::new();
                    for (_, coeff) in &mut external_shift {
                        *coeff = -*coeff;
                    }
                }
                if !negative_energies.is_empty() {
                    return Err(eyre::eyre!(
                        "generalized 3D production cannot convert E-surface with mixed-sign internal terms {:?}",
                        surface
                    ));
                }
                let mut esurface = Esurface {
                    energies: positive_energies,
                    external_shift,
                    vertex_set: VertexSet::dummy(),
                };
                if let Some(shift_rewrite) = canonize_esurface {
                    esurface.canonicalize_shift(shift_rewrite);
                }

                let esurface_id = self
                    .surface_cache
                    .esurface_cache
                    .position(|existing| *existing == esurface)
                    .unwrap_or_else(|| {
                        self.surface_cache.esurface_cache.push(esurface);
                        Into::<EsurfaceID>::into(self.surface_cache.esurface_cache.len() - 1)
                    });

                (HybridSurfaceID::Esurface(esurface_id), sign)
            }
            LinearSurfaceKind::Hsurface => {
                let hsurface = Hsurface {
                    positive_energies,
                    negative_energies,
                    external_shift,
                    vertex_set: VertexSet::dummy(),
                };
                let hsurface_id = self
                    .surface_cache
                    .hsurface_cache
                    .position(|existing| existing == &hsurface)
                    .unwrap_or_else(|| {
                        self.surface_cache.hsurface_cache.push(hsurface);
                        Into::<HsurfaceID>::into(self.surface_cache.hsurface_cache.len() - 1)
                    });

                (HybridSurfaceID::Hsurface(hsurface_id), 1)
            }
        };
        Ok(SurfaceMapEntry { surface_id, sign })
    }
}

fn three_d_source_summary(parsed: &three_dimensional_reps::ParsedGraph) -> String {
    let internal_edges = parsed
        .internal_edges
        .iter()
        .map(|edge| {
            format!(
                "edge {}: {} -> {}, loops={:?}, externals={:?}",
                edge.edge_id,
                edge.tail,
                edge.head,
                edge.signature.loop_signature,
                edge.signature.external_signature
            )
        })
        .join("; ");
    format!(
        "3D source has {} loop names {:?}, {} internal edges [{}]",
        parsed.loop_names.len(),
        parsed.loop_names,
        parsed.internal_edges.len(),
        internal_edges
    )
}

fn source_contract_edges_for_3d_expression(
    contract_edges: &[EdgeIndex],
    initial_state_cut_edges: &[EdgeIndex],
) -> Vec<EdgeIndex> {
    let non_contractible_edges = initial_state_cut_edges
        .iter()
        .copied()
        .collect::<HashSet<_>>();
    contract_edges
        .iter()
        .copied()
        .filter(|edge_id| !non_contractible_edges.contains(edge_id))
        .collect_vec()
}

fn remap_generated_surface_id(
    surface_id: &mut HybridSurfaceID,
    linear_surface_map: &BTreeMap<LinearSurfaceID, SurfaceMapEntry>,
) -> i64 {
    if let HybridSurfaceID::Linear(linear_surface_id) = surface_id {
        let entry = *linear_surface_map
            .get(linear_surface_id)
            .expect("all generated linear surfaces should have been interned");
        *surface_id = entry.surface_id;
        entry.sign
    } else {
        1
    }
}

fn remap_denominator_tree_surface_ids(
    denominator: &Tree<HybridSurfaceID>,
    linear_surface_map: &BTreeMap<LinearSurfaceID, SurfaceMapEntry>,
) -> Vec<RemappedDenominatorTree> {
    let node_signs = denominator
        .iter_nodes()
        .map(|node| {
            let sign = if let HybridSurfaceID::Linear(linear_surface_id) = node.data {
                linear_surface_map
                    .get(&linear_surface_id)
                    .map(|entry| entry.sign)
                    .unwrap_or(1)
            } else {
                1
            };
            (node.node_id, sign)
        })
        .collect::<BTreeMap<_, _>>();

    let mut chains_by_signature =
        BTreeMap::<(i64, Vec<(HybridSurfaceID, i64)>), BTreeMap<Vec<HybridSurfaceID>, usize>>::new(
        );
    for leaf in denominator.get_bottom_layer() {
        let mut sign = 1;
        let mut surface_signs = BTreeMap::<HybridSurfaceID, i64>::new();
        let mut chain = Vec::new();
        let mut current = Some(leaf);
        while let Some(node_id) = current {
            let mut surface_id = denominator.get_node(node_id).data;
            let node_sign = node_signs.get(&node_id).copied().unwrap_or(1);
            sign *= node_sign;
            remap_generated_surface_id(&mut surface_id, linear_surface_map);
            if surface_id != HybridSurfaceID::Unit {
                if node_sign < 0 {
                    let entry = surface_signs.entry(surface_id).or_insert(1);
                    *entry *= node_sign;
                    if *entry > 0 {
                        surface_signs.remove(&surface_id);
                    }
                }
                chain.push(surface_id);
            }
            current = denominator.get_node(node_id).parent;
        }
        chain.reverse();
        chains_by_signature
            .entry((sign, surface_signs.into_iter().collect()))
            .or_default()
            .entry(chain)
            .and_modify(|multiplicity| *multiplicity += 1)
            .or_insert(1);
    }

    if chains_by_signature.is_empty() {
        return vec![RemappedDenominatorTree {
            prefactor: 1,
            surface_signs: BTreeMap::new(),
            denominator: Tree::from_root(HybridSurfaceID::Unit),
        }];
    }

    chains_by_signature
        .into_iter()
        .flat_map(|((prefactor_sign, surface_signs), chains)| {
            // Surface projection is not injective: distinct occurrence-local
            // chains can become the same physical chain. The tree shares
            // equal branches structurally, so retain their additive count in
            // the variant prefactor and still factor chains with equal counts.
            let mut chains_by_multiplicity = BTreeMap::<usize, Vec<Vec<HybridSurfaceID>>>::new();
            for (chain, multiplicity) in chains {
                chains_by_multiplicity
                    .entry(multiplicity)
                    .or_default()
                    .push(chain);
            }
            chains_by_multiplicity
                .into_iter()
                .map(move |(multiplicity, chains)| RemappedDenominatorTree {
                    prefactor: prefactor_sign
                        * i64::try_from(multiplicity)
                            .expect("denominator-chain multiplicity should fit in i64"),
                    surface_signs: surface_signs.iter().copied().collect(),
                    denominator: denominator_tree_from_chains(&chains),
                })
        })
        .collect()
}

fn denominator_tree_from_chains(chains: &[Vec<HybridSurfaceID>]) -> Tree<HybridSurfaceID> {
    if chains.is_empty() || chains.iter().all(Vec::is_empty) {
        return Tree::from_root(HybridSurfaceID::Unit);
    }

    let mut tree = Tree::from_root(HybridSurfaceID::Unit);
    for chain in chains {
        if chain.is_empty() {
            insert_terminal_unit_if_missing(&mut tree, NodeId::root());
            continue;
        }
        let mut parent = NodeId::root();
        for surface_id in chain {
            let existing_child = tree
                .get_node(parent)
                .children
                .iter()
                .copied()
                .find(|child| tree.get_node(*child).data == *surface_id);
            if let Some(child) = existing_child {
                parent = child;
                continue;
            }
            let child = NodeId(tree.get_num_nodes());
            tree.insert_node(parent, *surface_id);
            parent = child;
        }
        insert_terminal_unit_if_missing(&mut tree, parent);
    }
    tree
}

fn insert_terminal_unit_if_missing(tree: &mut Tree<HybridSurfaceID>, parent: NodeId) {
    let has_terminal_unit = tree
        .get_node(parent)
        .children
        .iter()
        .any(|child| tree.get_node(*child).data == HybridSurfaceID::Unit);
    if !has_terminal_unit {
        tree.insert_node(parent, HybridSurfaceID::Unit);
    }
}

fn collect_linear_surface_terms(
    terms: &[(EdgeIndex, Rational)],
    initial_state_cut_edges: &[EdgeIndex],
    positive_energies: &mut Vec<EdgeIndex>,
    negative_energies: &mut Vec<EdgeIndex>,
    external_shift: &mut Vec<(EdgeIndex, i64)>,
) -> Result<()> {
    for (edge_id, coeff) in terms {
        let coeff = integer_coeff(coeff)?;
        if coeff == 0 {
            continue;
        }
        if initial_state_cut_edges.contains(edge_id) {
            // GammaLoop cut E-surfaces store initial-state energies on the
            // external-shift side with the opposite sign.
            external_shift.push((*edge_id, -coeff));
            continue;
        }

        let target = if coeff > 0 {
            &mut *positive_energies
        } else {
            &mut *negative_energies
        };
        for _ in 0..coeff.unsigned_abs() {
            target.push(*edge_id);
        }
    }
    Ok(())
}

fn integer_coeff(coeff: &Rational) -> Result<i64> {
    coeff
        .is_integer()
        .then(|| coeff.numerator_ref().to_i64())
        .flatten()
        .ok_or_else(|| eyre::eyre!("expected integer linear-surface coefficient, found {coeff}"))
}

fn numerator_sampling_scale_mode(
    setting: UniformNumeratorSamplingScale,
) -> NumeratorSamplingScaleMode {
    match setting {
        UniformNumeratorSamplingScale::None => NumeratorSamplingScaleMode::None,
        UniformNumeratorSamplingScale::BeyondQuadratic => {
            NumeratorSamplingScaleMode::BeyondQuadratic
        }
        UniformNumeratorSamplingScale::All => NumeratorSamplingScaleMode::All,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        cff::{expression::GammaLoopOrientationExpression, surface::GammaLoopSurfaceCache},
        dot,
        graph::{
            ExactUvSubLmbFrame, FeynmanGraph, FourDDenominator,
            cuts::{CutSet, LuCutSelection},
            parse::from_dot::IntoGraph,
        },
        initialisation::test_initialise,
        settings::global::{GenerationSettings, OrientationPattern},
        utils::GS,
        uv::uv_graph::UVE,
    };
    use linnet::half_edge::subgraph::SubSetOps;
    use symbolica::atom::{AtomCore, FunctionBuilder};

    #[test]
    fn generation_cache_obeys_byte_entry_and_access_order_limits() {
        let mut cache = GenerationCache::new(10, 2);
        cache.insert(1, "first", 4);
        cache.insert(2, "second", 4);
        assert_eq!(cache.get(&1), Some(&"first"));
        cache.insert(3, "third", 4);
        assert_eq!(cache.get(&2), None);
        assert_eq!(cache.get(&1), Some(&"first"));
        assert_eq!(cache.get(&3), Some(&"third"));

        cache.insert(4, "oversized", 11);
        assert_eq!(cache.get(&4), None);
        cache.insert(3, "replacement", 7);
        assert_eq!(cache.get(&1), None);
        assert_eq!(cache.get(&3), Some(&"replacement"));
        cache.clear();
        assert_eq!(cache.get(&3), None);

        let mut disabled = GenerationCache::new(0, 10);
        disabled.insert(1, "uncached", 0);
        assert_eq!(disabled.get(&1), None);
    }

    #[test]
    fn canonical_sunset_allocator_matches_native_small_oracle() -> Result<()> {
        use crate::{
            cff::surface::GammaLoopLinearEnergyExpr,
            numerator::energy_degree::{EnergyPowerAnalyzer, EquivalentEnergyCandidates},
            utils::symbols::UvDenominatorClassId,
        };
        use three_dimensional_reps::{
            CffGlobalPrefactorSign, MomentumSignature, graph_io::ParsedGraphInternalEdge,
        };

        test_initialise()?;
        // Fixed source incidence: two short channels and five serial copies of
        // the third. This polynomial diagnostic is independent of reconstruction.
        let parsed = ParsedGraph {
            internal_edges: [(0, 1, vec![1, 0]), (0, 1, vec![-1, 1])]
                .into_iter()
                .chain((1..6).map(|tail| (tail, (tail + 1) % 6, vec![0, 1])))
                .enumerate()
                .map(
                    |(edge_id, (tail, head, loop_signature))| ParsedGraphInternalEdge {
                        edge_id,
                        tail,
                        head,
                        label: format!("q{edge_id}"),
                        mass_key: Some("m_uv".into()),
                        signature: MomentumSignature {
                            loop_signature,
                            external_signature: Vec::new(),
                        },
                        had_pow: false,
                    },
                )
                .collect(),
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["K0".into(), "K1".into()],
            external_names: Vec::new(),
            node_name_to_internal: (0..6).map(|node| (format!("n{node}"), node)).collect(),
        };
        let mut options = Generate3DExpressionOptions {
            cff_generation_context: CffGenerationContext::EmbeddedCffFactor,
            ..Default::default()
        };
        let mut requests = Vec::new();
        for a in 0..=5 {
            for b in 0..=5 - a {
                for c in 0..=5 - a - b {
                    for d in 0..=5 - a - b - c {
                        let loads = [a, b, c, d, 5 - a - b - c - d];
                        let bounds = std::iter::once((0, 1))
                            .chain(loads.iter().enumerate().filter_map(|(index, degree)| {
                                (*degree > 0).then_some((index + 2, *degree))
                            }))
                            .collect::<Vec<_>>();
                        requests.push((5, bounds));
                    }
                }
            }
        }
        requests.push((7, vec![(0, 1), (2, 7)]));

        let classes = [UvDenominatorClassId(0), UvDenominatorClassId(1)];
        let mut candidates = EquivalentEnergyCandidates::try_from_source_occurrences([])?;
        candidates.add_uv_classes([(classes[0], vec![0]), (classes[1], (2..7).collect())])?;
        let q0 = FunctionBuilder::new(GS.emr_mom)
            .add_arg(GS.uv_class_ref(classes[0]))
            .add_arg(GS.cind(0))
            .finish();
        let q1 = FunctionBuilder::new(GS.emr_mom)
            .add_arg(GS.uv_class_ref(classes[1]))
            .add_arg(GS.cind(0))
            .finish();
        let analyzer = EnergyPowerAnalyzer::for_physical_emr_edges([]);
        for degree in [5, 7] {
            let plans =
                analyzer.plan_atom_assignment_proposals(&(&q0 * q1.pow(degree)), &candidates)?;
            assert!(!plans.is_empty());
            let seen = plans
                .iter()
                .map(|plan| plan.energy_degree_bounds().to_vec())
                .collect::<Vec<_>>();
            for plan in plans {
                requests.push((degree, plan.energy_degree_bounds().to_vec()));
                if let Some(challenger) = plan.placement_challenger(&seen)? {
                    requests.push((degree, challenger.energy_degree_bounds().to_vec()));
                }
            }
        }
        // Compare the complete signed contour for every degree distribution,
        // including the bounded allocator's proposals. Native row counts and
        // proposal ordering may change without changing this public value.
        let energy_replacements = (0..7)
            .map(|edge| {
                symbolica::id::Replacement::new(
                    GS.ose(EdgeIndex(edge)).to_pattern(),
                    Atom::num(match edge {
                        0 => 2,
                        1 => 3,
                        _ => 7,
                    })
                    .to_pattern(),
                )
            })
            .collect::<Vec<_>>();
        // Integrating k0 first gives -k1/[2 E1 (k1²-(E0+E1)²)].
        // The remaining Below contour is the negative sum of the simple pole
        // at E0+E1 and the fifth-order pole at E2. At (E0,E1,E2)=(2,3,7),
        // these independent residues give the exact values below.
        for (degree, bounds) in requests {
            options.energy_degree_bounds = Some(bounds.clone());
            let generated = three_dimensional_reps::generate_3d_expression(&parsed, &options)?;
            let surfaces = generated.expression.surfaces.get_all_replacements_gs(&[]);
            let mut contour = Atom::Zero;
            for orientation in &generated.expression.orientations {
                let numerator = bounds
                    .iter()
                    .map(|(edge, exponent)| {
                        orientation.edge_energy_map[*edge]
                            .to_atom_gs(&[])
                            .replace_multiple(&energy_replacements)
                            .pow(*exponent as u64)
                    })
                    .product::<Atom>();
                contour += orientation
                    .to_atom_gs()
                    .replace_multiple(&surfaces)
                    .replace_multiple(&energy_replacements)
                    * numerator;
            }
            let source_sign =
                generated
                    .energy_factor_components
                    .iter()
                    .fold(1, |sign, component| {
                        let frame = match component.ownership {
                            CffEnergyFactorOwnership::GlobalSourceProduct => {
                                component.core_global_prefactor_sign
                            }
                            CffEnergyFactorOwnership::VariantLocal => {
                                component.denominator_only_global_prefactor_sign
                            }
                        };
                        sign * CffGlobalPrefactorSign::from_exponent(
                            component.internal_edge_ids.len(),
                        )
                        .product(frame)
                        .factor()
                    });
            contour *= Atom::num(source_sign);
            let expected = match degree {
                5 => Atom::num(1213) / Atom::num(16_387_080_192i64),
                7 => -Atom::num(365) / Atom::num(47_775_744),
                _ => unreachable!(),
            };
            assert_eq!(contour, expected, "degree {degree}, bounds {bounds:?}");
        }
        Ok(())
    }

    #[test]
    fn projected_denominator_chains_preserve_occurrence_multiplicity() {
        let local_surfaces = [0, 1, 2, 3].map(LinearSurfaceID);
        let physical_surfaces = [
            HybridSurfaceID::Esurface(EsurfaceID(0)),
            HybridSurfaceID::Esurface(EsurfaceID(1)),
            HybridSurfaceID::Esurface(EsurfaceID(2)),
        ];
        let denominator = denominator_tree_from_chains(&[
            vec![HybridSurfaceID::Linear(local_surfaces[0])],
            vec![HybridSurfaceID::Linear(local_surfaces[1])],
            vec![HybridSurfaceID::Linear(local_surfaces[2])],
            vec![HybridSurfaceID::Linear(local_surfaces[3])],
        ]);
        let surface_map = BTreeMap::from([
            (
                local_surfaces[0],
                SurfaceMapEntry {
                    surface_id: physical_surfaces[0],
                    sign: 1,
                },
            ),
            (
                local_surfaces[1],
                SurfaceMapEntry {
                    surface_id: physical_surfaces[0],
                    sign: 1,
                },
            ),
            (
                local_surfaces[2],
                SurfaceMapEntry {
                    surface_id: physical_surfaces[1],
                    sign: 1,
                },
            ),
            (
                local_surfaces[3],
                SurfaceMapEntry {
                    surface_id: physical_surfaces[2],
                    sign: 1,
                },
            ),
        ]);

        let remapped = remap_denominator_tree_surface_ids(&denominator, &surface_map);
        let actual = remapped.iter().fold(Atom::new(), |sum, term| {
            sum + Atom::num(term.prefactor) * term.denominator.to_atom_inv()
        });
        let expected = Atom::num(2) / Atom::from(physical_surfaces[0])
            + Atom::num(1) / Atom::from(physical_surfaces[1])
            + Atom::num(1) / Atom::from(physical_surfaces[2]);

        assert!(
            (actual.clone() - expected).expand().is_zero(),
            "projecting distinct local chains onto one physical chain lost multiplicity: {actual}"
        );

        let negative_map = BTreeMap::from([
            (
                local_surfaces[0],
                SurfaceMapEntry {
                    surface_id: physical_surfaces[0],
                    sign: -1,
                },
            ),
            (
                local_surfaces[1],
                SurfaceMapEntry {
                    surface_id: physical_surfaces[0],
                    sign: -1,
                },
            ),
        ]);
        let negative_denominator = denominator_tree_from_chains(&[
            vec![HybridSurfaceID::Linear(local_surfaces[0])],
            vec![HybridSurfaceID::Linear(local_surfaces[1])],
        ]);
        let negative = remap_denominator_tree_surface_ids(&negative_denominator, &negative_map);
        let value = negative.iter().fold(Atom::Zero, |sum, term| {
            sum + Atom::num(term.prefactor) * term.denominator.to_atom_inv()
        });
        assert!(
            (value + Atom::num(2) / Atom::from(physical_surfaces[0]))
                .together()
                .is_zero()
        );
    }

    #[test]
    fn exact_cff_cache_preserves_occurrence_capacities() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_occurrence_cache {
            edge [num=1 mass=1]
            node [num=1]
            a -> b [id=0 lmb_id=0]
            b -> a [id=1]
        })?;
        let owners = [EdgeIndex(0), EdgeIndex(1)];
        let denominators = owners.map(|source_edge| FourDDenominator {
            source_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(source_edge))
                .finish(),
            mass_squared: Atom::one(),
            full_expr: Atom::one(),
        });
        let source = GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
            &graph,
            &denominators,
            owners,
            [],
            &graph.loop_momentum_basis,
            ExactUvSubLmbFrame::TaylorVacuum,
        )?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let left = GS.emr_mom(owners[0], GS.cind(0));
        let right = GS.emr_mom(owners[1], GS.cind(0));
        let numerators = [left.clone().pow(2), right.clone().pow(2), left * right];
        let mut cache = Local4dProjectionContext::default();
        for numerator in numerators.iter().cycle().take(2 * numerators.len()) {
            let mut values = Vec::new();
            for generation_cache in [Some(&mut cache), None] {
                let preparation =
                    graph.prepare_3d_expression_for_4d_term(&source, &options, numerator, &[])?;
                let (generated, mapper, plan, _) =
                    graph.generate_3d_expression_for_4d_term(&preparation, generation_cache)?;
                let mapper = &mapper.mapper;
                let surfaces = generated.expression.surfaces.get_all_replacements_gs(&[]);
                let mut value = Atom::Zero;
                for orientation in &generated.expression.orientations {
                    value += orientation.to_atom_gs().replace_multiple(&surfaces)
                        * mapper.map_planned_numerator(
                            &orientation.loop_energy_map,
                            &orientation.edge_energy_map,
                            &plan,
                        )?;
                }
                values.push(value.replace_multiple(mapper.exact_ose_replacements()));
            }
            assert!(
                !values[0].collect_factors().is_zero(),
                "the occurrence-capacity requests must exercise nonzero residues"
            );
            assert!(
                (values[0].collect_factors() - values[1].collect_factors())
                    .collect_factors()
                    .is_zero(),
                "cached and uncached complete residues must agree for each independent numerator capacity: {numerator}"
            );
        }
        Ok(())
    }

    #[test]
    fn exact_cff_batch_reuses_owner_relabelled_sub_lmb_topology() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_owner_relabelled_cache {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
        })?;
        let denominators = [EdgeIndex(0), EdgeIndex(1)].map(|source_edge| FourDDenominator {
            source_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(source_edge))
                .finish(),
            mass_squared: graph.underlying[source_edge].mass_atom().pow(2),
            full_expr: Atom::one(),
        });
        let mut relabelled = denominators.clone();
        relabelled[1].source_edge = EdgeIndex(0);
        let options = graph.denominator_only_cff_3d_expression_options();

        let mut cache = crate::uv::approx::projected_4d::Local4dProjectionContext::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let mut values = Vec::new();
        for denominators in [&denominators, &relabelled] {
            for generation_cache in [Some(&mut cache), None] {
                let (cff, _) = graph.clone().cff_from_4d_denominators_in_uv_sub_lmb(
                    denominators,
                    [EdgeIndex(0), EdgeIndex(1)],
                    [],
                    &graph.loop_momentum_basis,
                    ExactUvSubLmbFrame::TaylorVacuum,
                    &cutset,
                    &options,
                    &Atom::one(),
                    generation_cache,
                    &[],
                )?;
                let mut value = Atom::Zero;
                for term in cff.terms.values() {
                    for orientation in &term.orientations {
                        value += &orientation.expression
                            * term.map_exact_source_numerator(&orientation.orientation, None)?;
                    }
                }
                values.push(value * Atom::num(cff.production_prefactor_factor()));
            }
        }
        // Unchanged denominators retain physical OSE symbols, while a relabelled
        // owner may retain a literal square root. Compare both at the same two
        // positive-energy loop points, after production has resolved surfaces,
        // occurrence energies, ownership and convention factors.
        for qx in [Atom::num(3) / Atom::num(4), Atom::num(4) / Atom::num(3)] {
            let energy = (Atom::one() + qx.clone().pow(2)).sqrt();
            let at_point = |value: &Atom| {
                let mut value = value.clone();
                for (edge, spatial_x) in [(EdgeIndex(0), qx.clone()), (EdgeIndex(1), -&qx)] {
                    value = value
                        .replace(GS.emr_mom(edge, GS.cind(1)))
                        .with(spatial_x)
                        .replace(GS.emr_mom(edge, GS.cind(2)))
                        .with(Atom::Zero)
                        .replace(GS.emr_mom(edge, GS.cind(3)))
                        .with(Atom::Zero)
                        .replace(GS.ose(edge))
                        .with(energy.clone());
                }
                value.collect_factors()
            };
            let expected = at_point(&values[0]);
            assert!(
                !expected.is_zero(),
                "the scalar bubble residue must be nonzero"
            );
            assert!(
                values
                    .iter()
                    .skip(1)
                    .all(|value| (at_point(value) - &expected).together().is_zero()),
                "cached and uncached complete residues must agree under compatible owner relabelling at qx={qx}"
            );
        }

        Ok(())
    }

    #[test]
    fn bounded_exact_cff_dispatch_preserves_owner_contour() -> Result<()> {
        use crate::cff::{
            expression::GammaLoopOrientationExpression, surface::GammaLoopSurfaceCache,
        };
        use three_dimensional_reps::CffGlobalPrefactorSign;

        test_initialise()?;
        let graph: Graph = dot!(digraph bounded_exact_dispatch {
            edge [num=1 mass=2]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            b -> a [id=1]
            a -> a [id=2 lmb_id=1]
        })?;
        let owners = [EdgeIndex(0), EdgeIndex(1), EdgeIndex(2)];
        let base = owners.map(|source_edge| FourDDenominator {
            source_edge,
            momentum: graph.loop_momentum_basis.loop_atom::<Atom>(
                source_edge,
                GS.emr_mom,
                &[],
                true,
            ),
            mass_squared: graph.underlying[source_edge].mass_atom().pow(2),
            full_expr: Atom::one(),
        });
        let options = Generate3DExpressionOptions {
            cff_generation_context: CffGenerationContext::EmbeddedCffFactor,
            numerator_sampling_scale: NumeratorSamplingScaleMode::None,
            ..graph.denominator_only_cff_3d_expression_options()
        };
        for dotted_owner in [0, 1] {
            let mut denominators = base.to_vec();
            denominators.extend([base[dotted_owner].clone(), base[dotted_owner].clone()]);
            let source = GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
                &graph,
                &denominators,
                owners,
                [],
                &graph.loop_momentum_basis,
                ExactUvSubLmbFrame::TaylorVacuum,
            )?;
            let tagged = |derived| {
                FunctionBuilder::new(GS.emr_mom)
                    .add_arg(GS.uv_momentum_provenance_tag(
                        Atom::num(dotted_owner as i64).as_view(),
                        derived,
                        base[dotted_owner].momentum.as_view(),
                    ))
                    .add_arg(GS.cind(0))
                    .finish()
            };
            let fixed = tagged(false);
            let derived = tagged(true);
            // Keep the original quartic factor on its retained occurrence; only
            // the two new denominator-derived powers may use the serial copies.
            let numerator = fixed.clone().pow(4) * derived.clone().pow(2);
            let physical_edges = source.physical_energy_edge_index_map().unwrap();
            // Copy assignments can tie when neither copy owns the loop basis.
            // Exercise both source topologies without constraining tie resolution.
            // Compare their complete contours against the independent analytic
            // residues below, rather than mock costs or an envelope formula.
            let mut cache = Local4dProjectionContext::default();
            for phase in ["fresh", "cached", "uncached", "disabled", "evicted"] {
                if phase == "disabled" {
                    cache.generation_cache.entries = GenerationCache::new(0, 0);
                } else if phase == "evicted" {
                    cache.generation_cache.entries = GenerationCache::new(64 * 1024 * 1024, 1);
                    let parsed = source.to_three_d_parsed_graph()?;
                    let key = ExactCffGenerationCache::generation_key(
                        &parsed,
                        &source.energy_edge_index_map(&parsed).unwrap(),
                        &options,
                    );
                    // A different, denominator-only request occupies the sole
                    // retention slot; the tested bounded request must evict it.
                    cache.generation_cache.record_count(key, 0);
                }
                let preparation =
                    graph.prepare_3d_expression_for_4d_term(&source, &options, &numerator, &[])?;
                let (generated, mapper, plan, _) = graph.generate_3d_expression_for_4d_term(
                    &preparation,
                    (phase != "uncached").then_some(&mut cache),
                )?;
                let mapper = &mapper.mapper;

                let replacements = generated.expression.surfaces.get_all_replacements_gs(&[]);
                let mut contour = Atom::Zero;
                for orientation in &generated.expression.orientations {
                    let mapped = mapper.map_planned_numerator(
                        &orientation.loop_energy_map,
                        &orientation.edge_energy_map,
                        &plan,
                    )?;
                    contour += orientation
                        .to_atom_gs()
                        .replace_multiple(&replacements)
                        .replace_multiple(mapper.exact_ose_replacements())
                        * mapped;
                }
                let source_sign =
                    generated
                        .energy_factor_components
                        .iter()
                        .fold(1, |sign, component| {
                            let frame = match component.ownership {
                                CffEnergyFactorOwnership::GlobalSourceProduct => {
                                    component.core_global_prefactor_sign
                                }
                                CffEnergyFactorOwnership::VariantLocal => {
                                    component.denominator_only_global_prefactor_sign
                                }
                            };
                            sign * CffGlobalPrefactorSign::from_exponent(
                                component.internal_edge_ids.len(),
                            )
                            .product(frame)
                            .factor()
                        });
                contour *= Atom::num(source_sign);
                let energy = |physical| {
                    let occurrence = physical_edges
                        .internal
                        .iter()
                        .find_map(|(occurrence, owner)| (*owner == physical).then_some(*occurrence))
                        .unwrap();
                    crate::utils::ose_atom_from_index(EdgeIndex(occurrence))
                        .replace_multiple(mapper.exact_ose_replacements())
                };
                // Independent clockwise Below contours are -5/(32E_b)
                // for q^6/D_b^4 and -1/(2E_t) for the attached tadpole.
                // The original quartic factor stays on its original occurrence.
                let expected = Atom::num(5) / (Atom::num(64) * energy(2) * energy(dotted_owner));
                // The exact mapper keeps an unchanged propagator as E(owner),
                // while a rewritten equal-momentum denominator uses its literal
                // square root. The repeated-pole oracle requires the same positive
                // physical mass shell for both dialects.
                let physical_shell = |mut expression: Atom| {
                    for physical in owners {
                        let on_shell_energy = (1..=3)
                            .fold(
                                graph.underlying[physical].mass_atom().pow(2),
                                |square, spatial_index| {
                                    square
                                        + graph
                                            .loop_momentum_basis
                                            .loop_atom::<Atom>(
                                                physical,
                                                GS.emr_mom,
                                                &[GS.cind(spatial_index)],
                                                true,
                                            )
                                            .pow(2)
                                },
                            )
                            .sqrt();
                        expression = expression.replace(GS.ose(physical)).with(on_shell_energy);
                    }
                    expression
                };
                let difference = (physical_shell(contour) - physical_shell(expected)).together();
                assert!(
                    difference.is_zero(),
                    "dotted owner {dotted_owner}, {phase}: dispatch changed the analytic contour: {difference}"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn direct_cff_matches_signed_bubble_and_contracted_unit_contours() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph direct_duplicate_sign {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
        })?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let cutset = CutSet::empty(graph.n_hedges());
        let contract: linnet::half_edge::subgraph::SuBitGraph = graph.empty_subgraph();
        let cff = graph.cff(
            &contract,
            &cutset,
            &OrientationPattern::default(),
            &options,
            None,
        )?;
        let mut value = cff
            .terms
            .values()
            .flat_map(|term| &term.orientations)
            .fold(Atom::Zero, |sum, orientation| sum + &orientation.expression)
            * Atom::num(cff.production_prefactor_factor());
        for edge in 0..2 {
            value = value.replace(GS.ose(EdgeIndex(edge))).with(1);
        }
        // At unit on-shell energies, Wick rotation reduces the scalar contour
        // to the independently integrated Euclidean propagator product.
        let signed_contour = value * (Atom::num(2) * Atom::var(GS.pi)).pow(3) / Atom::i().pow(1);
        assert!(
            (signed_contour - Atom::one() / Atom::num(4))
                .together()
                .is_zero()
        );

        let contracted_loop = graph
            .get_edge_subgraph(EdgeIndex(0))
            .union(&graph.get_edge_subgraph(EdgeIndex(1)));
        assert_eq!(graph.cyclotomatic_number(&contracted_loop), 1);
        let reduced = graph.cff(
            &contracted_loop,
            &cutset,
            &OrientationPattern::default(),
            &options,
            None,
        )?;
        let reduced_value = reduced
            .terms
            .values()
            .flat_map(|term| &term.orientations)
            .fold(Atom::Zero, |sum, orientation| sum + &orientation.expression)
            * Atom::num(reduced.production_prefactor_factor());
        assert_eq!(
            reduced_value,
            Atom::one(),
            "contracting the only loop must leave the unit zero-loop source"
        );
        Ok(())
    }

    #[test]
    fn two_loop_cff_matches_signed_theta_contour() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph direct_two_loop_sign {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            a -> b [id=1 lmb_id=1]
            b -> a [id=2]
        })?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let cutset = CutSet::empty(graph.n_hedges());
        let contract: linnet::half_edge::subgraph::SuBitGraph = graph.empty_subgraph();
        let cff = graph.cff(
            &contract,
            &cutset,
            &OrientationPattern::default(),
            &options,
            None,
        )?;
        let mut value = cff
            .terms
            .values()
            .flat_map(|term| &term.orientations)
            .fold(Atom::Zero, |sum, orientation| sum + &orientation.expression)
            * Atom::num(cff.production_prefactor_factor());
        for edge in 0..3 {
            value = value.replace(GS.ose(EdgeIndex(edge))).with(1);
        }
        // At unit on-shell energies, Wick rotation reduces the scalar contour
        // to the independently integrated Euclidean propagator product.
        let signed_contour = value * (Atom::num(2) * Atom::var(GS.pi)).pow(6) / Atom::i().pow(2);
        assert!(
            (signed_contour - -Atom::one() / Atom::num(12))
                .together()
                .is_zero()
        );

        Ok(())
    }

    #[test]
    fn contracted_raised_generation_rejects_absent_emr_bound_without_aliasing() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(
            digraph factorized_source_bounds {
                edge [pdg=1000 num=1 mass=0]
                node [num=1]
                ext_in [style=invis]
                ext_out [style=invis]
                ext_in -> A [id=4]
                C -> ext_out [id=5]
                A -> B [id=0 num="Q(0,spenso::cind(0))+1"]
                B -> C [id=1 num="Q(1,spenso::cind(0))+2"]
                C -> D [id=2]
                D -> A [id=3]
            },
            "scalars"
        )?;
        let options = graph.production_cff_3d_expression_options(&GenerationSettings::default())?;
        assert_eq!(options.energy_degree_bounds, Some(vec![(0, 1), (1, 1)]));

        let numerator = graph.production_numerator_atom_for_full_3d_expression();
        let mut contracted_graph = graph.clone();
        let error = contracted_graph
            .generate_3d_expression_for_integrand(
                &[EdgeIndex(0)],
                &None,
                &options,
                Some(&numerator),
            )
            .expect_err("a bound owned by a shrunken EMR edge must not be silently reassigned");
        let error = format!("{error:#}");
        assert!(
            error.contains("shrunken EMR edge 0")
                && error.contains("no numerator-only/shrunken-edge bound channel"),
            "unexpected error: {error}"
        );
        Ok(())
    }

    #[test]
    fn selected_raw_cff_matches_normal_generation_after_discarding_a_trial() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph selected_raw_cff {
            edge [num=1 mass=1]
            node [num=1]
            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
            b -> b [id=2 lmb_id=1]
        })?;
        let mut ordinary_graph = graph.clone();
        let contract = graph.get_edge_subgraph(EdgeIndex(2));
        let options = Generate3DExpressionOptions {
            cff_generation_context: CffGenerationContext::EmbeddedCffFactor,
            ..graph.denominator_only_cff_3d_expression_options()
        };
        let numerator = GS.emr_mom(EdgeIndex(0), GS.cind(0)).pow(2);
        let discarded = graph.generate_raw_3d_expression_for_integrand(
            &[EdgeIndex(2)],
            &options,
            Some(&GS.emr_mom(EdgeIndex(0), GS.cind(0)).pow(3)),
        )?;
        drop(discarded);
        let selected = graph.generate_raw_3d_expression_for_integrand(
            &[EdgeIndex(2)],
            &options,
            Some(&numerator),
        )?;
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let selected = graph.convert_generated_expression_surfaces(selected, &canonization, &[])?;
        let cutset = CutSet::empty(graph.n_hedges());
        let pattern = OrientationPattern::default();
        let selected = graph
            .cff_from_generated_expression(selected, &contract, &cutset, &pattern, &options)?;
        let ordinary =
            ordinary_graph.cff(&contract, &cutset, &pattern, &options, Some(&numerator))?;
        assert_eq!(
            selected.energy_degree_bound_report,
            ordinary.energy_degree_bound_report
        );
        assert_eq!(
            selected.terms.keys().collect::<Vec<_>>(),
            ordinary.terms.keys().collect::<Vec<_>>()
        );
        // Discarded raw proposals must not affect any selected physical residue;
        // persistent surface numbering and interpolation layout remain internal.
        for (index, selected_term) in &selected.terms {
            let ordinary_term = &ordinary.terms[index];
            let values = [
                (
                    selected_term,
                    selected.production_prefactor_factor(),
                    &graph,
                ),
                (
                    ordinary_term,
                    ordinary.production_prefactor_factor(),
                    &ordinary_graph,
                ),
            ]
            .map(|(term, prefactor, graph)| {
                term.orientations
                    .iter()
                    .fold(Atom::Zero, |sum, orientation| {
                        sum + &orientation.expression
                            * numerator.replace_multiple(
                                orientation.orientation.energy_replacements_gs(graph),
                            )
                    })
                    * Atom::num(prefactor)
            });
            assert!(
                (values[0].collect_factors() - values[1].collect_factors())
                    .collect_factors()
                    .is_zero(),
                "discarded trials must leave the complete physical residue unchanged at {index}"
            );
        }
        Ok(())
    }

    #[test]
    fn raised_lu_cff_preserves_complete_production_residues() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph raised_lu {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> v1 [id=0]
            v1 -> v2 [id=1 lmb_id=0]
            v2 -> v3 [id=2]
            v1 -> v3 [id=3]
            v3 -> outgoing [id=4]
        })?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let numerator = GS.emr_mom(EdgeIndex(1), GS.cind(0)).pow(2);
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let generated = graph.generate_3d_expression_for_integrand(
            &[],
            &canonization,
            &options,
            Some(&numerator),
        )?;
        let production = &generated.expression;
        let lu_cut = graph
            .determine_raised_esurfaces_from_expression(production)
            .raised_groups
            .into_iter()
            .find(|group| {
                group.max_occurence > 1
                    && group.esurface_ids.iter().any(|esurface_id| {
                        !production.surfaces.esurface_cache[*esurface_id]
                            .external_shift
                            .is_empty()
                    })
            })
            .expect("raised production CFF should contain a physical repeated LU surface");
        let mut cutset = CutSet::empty(graph.n_hedges());
        cutset.residue_selector.lu = Some(LuCutSelection {
            raised_group: lu_cut.clone(),
            cut_edge_alternatives: lu_cut
                .esurface_ids
                .iter()
                .map(|esurface_id| {
                    production.surfaces.esurface_cache[*esurface_id]
                        .energies
                        .clone()
                })
                .collect(),
        });
        let direct = graph.cff(
            &cutset.union,
            &cutset,
            &OrientationPattern::default(),
            &options,
            Some(&numerator),
        )?;

        assert_eq!(direct.terms.len(), lu_cut.max_occurence);
        assert!(direct.terms.keys().all(|index| {
            index.lu_cut_order.is_some()
                && index.left_threshold_order.is_none()
                && index.right_threshold_order.is_none()
        }));
        assert!(direct.terms.iter().any(|(index, term)| {
            index.lu_cut_order == Some(2) && !term.orientations.is_empty()
        }));
        let stored = graph.cff_from_production_expression(
            &generated,
            &cutset,
            &OrientationPattern::default(),
        )?;
        for (index, term) in &direct.terms {
            let value = |term: &crate::cff::CFFTerm, prefactor| {
                term.orientations
                    .iter()
                    .fold(Atom::Zero, |sum, orientation| {
                        sum + &orientation.expression
                            * numerator.replace_multiple(
                                orientation.orientation.energy_replacements_gs(&graph),
                            )
                    })
                    * Atom::num(prefactor)
            };
            assert!(
                (value(term, direct.production_prefactor_factor()).collect_factors()
                    - value(&stored.terms[index], stored.production_prefactor_factor())
                        .collect_factors())
                .is_zero(),
                "selected LU generation must preserve the complete stored production residue"
            );
        }
        Ok(())
    }

    #[test]
    fn exact_massless_raised_lu_cff_retains_selected_residue() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_massless_raised_lu {
            edge [num=1 mass="UFO::ZERO"]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> v1 [id=0]
            v1 -> v2 [id=1 lmb_id=0]
            v2 -> v3 [id=2]
            v1 -> v3 [id=3]
            v3 -> outgoing [id=4]
        })?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let numerator = GS.emr_mom(EdgeIndex(1), GS.cind(0)).pow(2);
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph
            .generate_3d_expression_for_integrand(&[], &canonization, &options, Some(&numerator))?
            .expression;
        let lu_cut = graph
            .determine_raised_esurfaces_from_expression(&production)
            .raised_groups
            .into_iter()
            .find(|group| {
                group.max_occurence > 1
                    && group.esurface_ids.iter().any(|esurface_id| {
                        !production.surfaces.esurface_cache[*esurface_id]
                            .external_shift
                            .is_empty()
                    })
            })
            .expect("raised production CFF should contain a physical repeated LU surface");
        let mut cutset = CutSet::empty(graph.n_hedges());
        cutset.residue_selector.lu = Some(LuCutSelection {
            raised_group: lu_cut.clone(),
            cut_edge_alternatives: lu_cut
                .esurface_ids
                .iter()
                .map(|esurface_id| {
                    production.surfaces.esurface_cache[*esurface_id]
                        .energies
                        .clone()
                })
                .collect(),
        });
        let denominators = [1, 2, 3].map(|edge| {
            let edge = EdgeIndex(edge);
            FourDDenominator {
                source_edge: edge,
                momentum: FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(edge))
                    .finish(),
                mass_squared: graph.underlying[edge].mass_atom().pow(2),
                full_expr: Atom::one(),
            }
        });
        assert!(
            denominators
                .iter()
                .all(|denominator| denominator.mass_squared.is_zero()),
            "UFO::ZERO masses must be normalized before exact surface projection"
        );
        let source = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        let mapper = source
            .exact_source_energy_mapper(&[])
            .expect("exact source has a physical EMR mapper");
        let candidates = mapper.equivalent_energy_candidates([EdgeIndex(1)])?;
        let plan = graph.plan_numerator_energy_assignment_in_atom_excluding(
            &numerator,
            std::iter::empty(),
            &candidates,
        )?;
        let assigned_occurrence = source
            .physical_energy_edge_index_map()
            .expect("exact source has a physical occurrence map")
            .internal
            .into_iter()
            .find_map(|(occurrence, owner)| (owner == 1).then_some(occurrence))
            .expect("the physical numerator carrier has one exact occurrence");
        assert_eq!(
            plan.energy_degree_bounds(),
            &[(assigned_occurrence, 2)],
            "canonical exact-edge IDs may change, but the quadratic bound must remain on the occurrence owned by physical edge 1",
        );

        let (exact, _) =
            graph.cff_from_4d_denominators(&denominators, &cutset, &options, &numerator)?;
        assert_eq!(
            exact.energy_degree_bound_report,
            CffEnergyDegreeBoundReport {
                source_kind: CffEnergyBoundSourceKind::ExactFourD,
                physical_parent_bounds: vec![(1, 2)],
                assigned_cff_source_bounds: plan.energy_degree_bounds().to_vec(),
            },
            "exact CFF diagnostics must keep physical-parent and assigned occurrence bounds distinct",
        );
        assert_eq!(exact.terms.len(), lu_cut.max_occurence);
        let selected = exact
            .terms
            .iter()
            .find(|(index, _)| index.lu_cut_order == Some(2))
            .expect("the exact raised CFF should retain its second LU residue");
        let selected_value = selected
            .1
            .orientations
            .iter()
            .map(|orientation| {
                Ok(&orientation.expression
                    * selected
                        .1
                        .map_exact_source_numerator(&orientation.orientation, None)?)
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .sum::<Atom>()
            * Atom::num(exact.production_prefactor_factor());
        assert!(
            !selected_value.collect_factors().is_zero(),
            "the complete mapped exact second-order LU residue must remain nonzero"
        );
        Ok(())
||||||| qrspmutx 7a048a39 "add feynkit ufo and Python bindings" (parents of rebased revision)
}

pub fn generate_cff_expression_from_subgraph<E, V, H, S: SubGraphLike>(
    graph: &HedgeGraph<E, V, H>,
    subgraph: &S,
    canonize_esurface: &Option<ShiftRewrite>,
    reversed_dangling: &[EdgeIndex],
    edges_in_initial_state_cut: &[EdgeIndex],
    surface_cache: &mut SurfaceCache,
) -> Result<CFFExpression<OrientationID>> {
    let graphs = get_orientations_from_subgraph(graph, subgraph, reversed_dangling);
    let cff = generate_cff_from_orientations(
        graphs,
        surface_cache,
        edges_in_initial_state_cut,
        canonize_esurface,
    )?;
    Ok(cff)
}

#[derive(Copy, Clone, Debug)]
pub struct ConstraintData<'a> {
    pub constraints: &'a [&'a Esurface],
    pub illegal_esurfaces: &'a [&'a Esurface],
}

#[derive(Copy, Clone, Debug)]
pub struct UvCffTopology<'a> {
    pub contract_edges: &'a [EdgeIndex],
    pub edges_in_initial_state_cut: &'a [EdgeIndex],
    pub orientation: &'a EdgeVec<Orientation>,
    pub cut_edges: &'a [EdgeIndex],
}

pub fn generate_uv_cff<E, V, H, S: SubGraphLike>(
    graph: &HedgeGraph<E, V, H>,
    subgraph: &S,
    canonize_esurface: &Option<ShiftRewrite>,
    topology: UvCffTopology<'_>,
    setup: PostProcessingSetup<'_>,
) -> Result<Atom> {
    let mut generation_graph =
        CFFGenerationGraph::new_from_subgraph(graph, topology.orientation.clone(), subgraph)?;

    for contracted_edge in topology.contract_edges {
        generation_graph = generation_graph.contract_edge(*contracted_edge);
    }

    generation_graph.remove_self_edges();

    if generation_graph.has_directed_cycle_initial() {
        return Ok(Atom::new());
    }

    let mut surface_cache = SurfaceCache {
        esurface_cache: EsurfaceCollection::from_iter(std::iter::empty()),
        hsurface_cache: HsurfaceCollection::from_iter(std::iter::empty()),
    };

    let generate_tree_for_orientation = generate_tree_for_orientation(
        generation_graph,
        &mut surface_cache,
        topology.edges_in_initial_state_cut,
        canonize_esurface,
    );

    let mut tree = generate_tree_for_orientation.map(forget_graphs);

    post_process(
        &mut tree,
        topology.orientation,
        subgraph,
        &surface_cache,
        setup,
    );

    let surface_cache_to_use = setup
        .rewrite_esurfaces
        .map_or(&surface_cache, |rewrite| rewrite.allowed_targets);

    let atom_tree = tree.to_atom_inv();
    let atom_tree_substituted =
        surface_cache_to_use.substitute_energies(&atom_tree, topology.cut_edges);
    let inverse_energies =
        get_cff_inverse_energy_product_impl(graph, subgraph, topology.contract_edges);

    Ok(atom_tree_substituted * &inverse_energies)
}

#[derive(Clone, Copy)]
pub struct PostProcessingSetup<'a> {
    pub constraint_data: Option<ConstraintData<'a>>,
    pub rewrite_esurfaces: Option<EsurfaceRewritingInstructions<'a>>,
}

#[derive(Clone, Copy)]
pub struct EsurfaceRewritingInstructions<'a> {
    pub allowed_targets: &'a SurfaceCache,
    pub graph: &'a Graph,
    pub cuts: &'a TiVec<CutId, CrossSectionCut>,
    pub subgraph_location: (Option<CutId>, Option<CutId>),
}

fn post_process<S: SubGraphLike>(
    tree: &mut Tree<HybridSurfaceID>,
    orientation: &EdgeVec<Orientation>,
    subgraph: &S,
    surface_cache: &SurfaceCache,
    setup: PostProcessingSetup<'_>,
) {
    if let Some(constraint_data) = setup.constraint_data {
        tree.map_mut(|surface_id| {
            let esurface_is_allowed = match surface_id {
                HybridSurfaceID::Esurface(esurface_id) => {
                    let esurface_to_compare = &surface_cache.esurface_cache[*esurface_id];
                    constraint_data
                        .illegal_esurfaces
                        .iter()
                        .all(|illegal_esurface| esurface_to_compare != *illegal_esurface)
                }
                HybridSurfaceID::Hsurface(hsurface_id) => {
                    let hsurface_to_compare = &surface_cache.hsurface_cache[*hsurface_id];
                    constraint_data
                        .illegal_esurfaces
                        .iter()
                        .all(|illegal_esurface| {
                            !hsurface_to_compare
                                .equality_under_energy_conservation(
                                    illegal_esurface,
                                    constraint_data.constraints,
                                )
                                .unwrap_or(
                                    hsurface_to_compare.equality_by_try_convert(illegal_esurface),
                                )
                        })
                }
                HybridSurfaceID::Unit => true,
                HybridSurfaceID::Infinite => true,
            };

            if !esurface_is_allowed {
                *surface_id = HybridSurfaceID::Infinite
            }
        });
    }

    if let Some(rewrite_esurfaces) = setup.rewrite_esurfaces {
        let hashset_of_appearing_ids = tree
            .iter_nodes()
            .map(|node| node.data)
            .collect::<HashSet<HybridSurfaceID>>();

        let mut id_map = HashMap::<HybridSurfaceID, HybridSurfaceID>::new();
        id_map.insert(HybridSurfaceID::Unit, HybridSurfaceID::Unit);
        id_map.insert(HybridSurfaceID::Infinite, HybridSurfaceID::Infinite);

        for appearing_id in hashset_of_appearing_ids.iter() {
            let surface_to_rewrite = surface_cache.get_surface(*appearing_id);

            match surface_to_rewrite {
                HybridSurfaceRef::Unit(_) => continue,
                HybridSurfaceRef::Infinite(_) => continue,
                HybridSurfaceRef::Esurface(esurface) => {
                    if let Some(esurface_id) = rewrite_esurfaces
                        .allowed_targets
                        .esurface_cache
                        .position(|allowed_esurface| allowed_esurface == esurface)
                    {
                        let new_id = HybridSurfaceID::Esurface(esurface_id);
                        id_map.insert(*appearing_id, new_id);
                    } else {
                        let complete_to_right =
                            if let Some(cut_id) = rewrite_esurfaces.subgraph_location.1 {
                                let edges_in_cut = rewrite_esurfaces
                                    .graph
                                    .iter_edges_of(&rewrite_esurfaces.cuts[cut_id].cut)
                                    .map(|(_, edge_id, _)| edge_id)
                                    .collect_vec();

                                edges_in_cut
                                    .iter()
                                    .all(|edge_id| esurface.energies.contains(edge_id))
                            } else {
                                false
                            };

                        let complete_to_left =
                            if let Some(cut_id) = rewrite_esurfaces.subgraph_location.0 {
                                let edges_in_cut = rewrite_esurfaces
                                    .graph
                                    .iter_edges_of(&rewrite_esurfaces.cuts[cut_id].cut)
                                    .map(|(_, edge_id, _)| edge_id)
                                    .collect_vec();

                                edges_in_cut
                                    .iter()
                                    .all(|edge_id| esurface.energies.contains(edge_id))
                            } else {
                                false
                            };

                        if complete_to_left && complete_to_right {
                            panic!("esurface has no connected component");
                        }

                        if !complete_to_left && !complete_to_right {
                            println!("esurface: {:#?}", esurface);
                            panic!("esurface cannot be rewritten to any allowed target");
                        }

                        let vertices_to_add = if complete_to_left {
                            let cut_id = rewrite_esurfaces.subgraph_location.0.unwrap();
                            &rewrite_esurfaces.cuts[cut_id].left
                        } else if complete_to_right {
                            let cut_id = rewrite_esurfaces.subgraph_location.1.unwrap();
                            &rewrite_esurfaces.cuts[cut_id].right
                        } else {
                            unreachable!()
                        };

                        let new_esurface_subgraph = esurface
                            .vertex_set
                            .subgraph(rewrite_esurfaces.graph)
                            .union(vertices_to_add);

                        let new_esurface = Esurface::new_from_subgraph(
                            &new_esurface_subgraph,
                            rewrite_esurfaces.graph,
                            orientation,
                        );

                        let new_esurface_id = rewrite_esurfaces
                            .allowed_targets
                            .esurface_cache
                            .position(|allowed_esurface| allowed_esurface == &new_esurface)
                            .expect("constructed esurface not in allowed targets");

                        let new_id = HybridSurfaceID::Esurface(new_esurface_id);
                        id_map.insert(*appearing_id, new_id);
                    }
                }
                HybridSurfaceRef::Hsurface(hsurface) => {
                    let complete_to_left =
                        if let Some(cut_id) = rewrite_esurfaces.subgraph_location.0 {
                            let edges_in_cut = rewrite_esurfaces
                                .graph
                                .iter_edges_of(&rewrite_esurfaces.cuts[cut_id].cut)
                                .map(|(_, edge_id, _)| edge_id)
                                .collect_vec();

                            hsurface
                                .negative_energies
                                .iter()
                                .all(|edge_id| edges_in_cut.contains(edge_id))
                        } else {
                            false
                        };

                    let complete_to_right =
                        if let Some(cut_id) = rewrite_esurfaces.subgraph_location.1 {
                            let edges_in_cut = rewrite_esurfaces
                                .graph
                                .iter_edges_of(&rewrite_esurfaces.cuts[cut_id].cut)
                                .map(|(_, edge_id, _)| edge_id)
                                .collect_vec();

                            hsurface
                                .negative_energies
                                .iter()
                                .all(|edge_id| edges_in_cut.contains(edge_id))
                        } else {
                            false
                        };

                    if complete_to_left && complete_to_right {
                        panic!(
                            "hsurface has no connected component supergraph, it cannot exist, but it does"
                        );
                    }

                    if !complete_to_left && !complete_to_right {
                        println!("hsurface: {:#?}", hsurface);
                        panic!("hsurface cannot be rewritten to any allowed target");
                    }

                    let vertices_to_add = if complete_to_left {
                        let cut_id = rewrite_esurfaces.subgraph_location.0.unwrap();
                        &rewrite_esurfaces.cuts[cut_id].left
                    } else if complete_to_right {
                        let cut_id = rewrite_esurfaces.subgraph_location.1.unwrap();
                        &rewrite_esurfaces.cuts[cut_id].right
                    } else {
                        unreachable!()
                    };

                    let new_esurface_subgraph = hsurface
                        .vertex_set
                        .subgraph(rewrite_esurfaces.graph)
                        .union(vertices_to_add);

                    let new_esurface = Esurface::new_from_subgraph(
                        &new_esurface_subgraph,
                        rewrite_esurfaces.graph,
                        orientation,
                    );

                    let new_esurface_id = rewrite_esurfaces
                        .allowed_targets
                        .esurface_cache
                        .position(|allowed_esurface| allowed_esurface == &new_esurface)
                        .unwrap_or_else(|| {
                            println!("for graph: {}", rewrite_esurfaces.graph.name.clone());
                            println!("dot: \n {}", rewrite_esurfaces.graph.debug_dot());
                            println!("subgraph: \n {}", rewrite_esurfaces.graph.dot(subgraph));

                            println!("from hsurface: {:?}", hsurface);
                            println!("constructed esurface: {:?}", new_esurface);
                            panic!("constructed esurface not in allowed targets");
                        });

                    let new_id = HybridSurfaceID::Esurface(new_esurface_id);
                    id_map.insert(*appearing_id, new_id);
                }
            }
        }

        tree.map_mut(|surface_id| *surface_id = id_map[surface_id]);
    }
}

fn generate_cff_from_orientations<O: From<usize> + Into<usize>>(
    orientations_and_graphs: Vec<CFFGenerationGraph>,
    generator_cache: &mut SurfaceCache,
    edges_in_initial_state_cut: &[EdgeIndex],
    canonize_esurface: &Option<ShiftRewrite>,
) -> Result<CFFExpression<O>, Report> {
    // filter cyclic orientations beforehand
    let acyclic_orientations_and_graphs = orientations_and_graphs
        .into_iter()
        .filter(|graph| !graph.has_directed_cycle_initial())
        .collect_vec();

    debug!(
        "number of acyclic orientations: {}",
        acyclic_orientations_and_graphs.len()
    );

    let terms = acyclic_orientations_and_graphs
        .into_iter()
        .map(|graph| {
            let global_orientation = graph.global_orientation.clone();
            let tree = generate_tree_for_orientation(
                graph.clone(),
                generator_cache,
                edges_in_initial_state_cut,
                canonize_esurface,
            );
            let expression = tree.map(forget_graphs);

            crate::cff::expression::OrientationExpression {
                expression,
                data: OrientationData {
                    orientation: global_orientation,
                },
            }
        })
        .collect_vec();

    Ok(CFFExpression {
        orientations: terms.into(),
        surfaces: generator_cache.clone(),
    })
}

#[derive(Clone, Debug, Serialize, Deserialize, Encode, Decode)]
pub struct SurfaceCache {
    #[bincode(with_serde)]
    pub esurface_cache: EsurfaceCollection, // Esurfaces of the supergraph
    #[bincode(with_serde)]
    pub hsurface_cache: HsurfaceCollection, // Anything else.
}

impl SurfaceCache {
    pub fn substitute_energies(&self, atom: &Atom, cut_edges: &[EdgeIndex]) -> Atom {
        let replacement_rules = self.get_all_replacements(cut_edges);
        atom.replace_multiple(&replacement_rules)
    }

    pub(crate) fn iter_all_surfaces(
        &'_ self,
    ) -> impl Iterator<Item = (HybridSurfaceID, HybridSurfaceRef<'_>)> + '_ {
        let esurface_id_iter = self.esurface_cache.iter_enumerated().map(|(id, esurface)| {
            (
                HybridSurfaceID::Esurface(id),
                HybridSurfaceRef::Esurface(esurface),
            )
        });

        let hsurface_id_iter = self.hsurface_cache.iter_enumerated().map(|(id, hsurface)| {
            (
                HybridSurfaceID::Hsurface(id),
                HybridSurfaceRef::Hsurface(hsurface),
            )
        });

        esurface_id_iter.chain(hsurface_id_iter)
    }

    pub(crate) fn get_all_replacements(&self, cut_edges: &[EdgeIndex]) -> Vec<Replacement> {
        self.iter_all_surfaces()
            .map(|(id, surface)| {
                let id_atom = Pattern::from(Atom::from(id));
                let surface_atom = Pattern::from(surface.to_atom(cut_edges));
                Replacement::new(id_atom, surface_atom)
            })
            .collect()
    }

    pub(crate) fn get_all_replacements_in_lmb(
        &self,
        cut_edges: &[EdgeIndex],
        lmb: &LoopMomentumBasis,
    ) -> Vec<Replacement> {
        self.iter_all_surfaces()
            .map(|(id, surface)| {
                let id_atom = Pattern::from(Atom::from(id));
                let surface_atom = Pattern::from(surface.to_atom_in_lmb(cut_edges, lmb));
                Replacement::new(id_atom, surface_atom)
            })
            .collect()
    }

    #[allow(dead_code)]
    pub(crate) fn get_surface(&self, surface_id: HybridSurfaceID) -> HybridSurfaceRef<'_> {
        match surface_id {
            HybridSurfaceID::Esurface(id) => HybridSurfaceRef::Esurface(&self.esurface_cache[id]),
            HybridSurfaceID::Hsurface(id) => HybridSurfaceRef::Hsurface(&self.hsurface_cache[id]),
            HybridSurfaceID::Unit => HybridSurfaceRef::Unit(UnitSurface {}),
            HybridSurfaceID::Infinite => HybridSurfaceRef::Infinite(InfiniteSurface {}),
        }
    }

    #[allow(dead_code)]
    pub(crate) fn new() -> Self {
        Self {
            esurface_cache: EsurfaceCollection::from_iter(std::iter::empty()),
            hsurface_cache: HsurfaceCollection::from_iter(std::iter::empty()),
        }
    }
}

fn generate_tree_for_orientation(
    graph: CFFGenerationGraph,
    generator_cache: &mut SurfaceCache,
    edges_in_initial_state_cut: &[EdgeIndex],
    canonize_esurface: &Option<ShiftRewrite>,
) -> Tree<GenerationData> {
    let mut tree = Tree::from_root(GenerationData {
        graph,
        surface_id: None,
    });

    while let Some(()) = advance_tree(
        &mut tree,
        generator_cache,
        edges_in_initial_state_cut,
        canonize_esurface,
    ) {}

    tree
}

fn advance_tree(
    tree: &mut Tree<GenerationData>,
    generator_cache: &mut SurfaceCache,
    edges_in_initial_state_cut: &[EdgeIndex],
    canonize_esurface: &Option<ShiftRewrite>,
) -> Option<()> {
    let bottom_layer = tree.get_bottom_layer();

    let (children_optional, new_surfaces_for_tree): (
        Vec<Option<Vec<CFFGenerationGraph>>>,
        Vec<HybridSurfaceID>,
    ) = bottom_layer
        .iter()
        .map(|&node_id| {
            let node = &tree.get_node(node_id);
            let graph = &node.data.graph;

            let (option_children, surface) = graph.generate_children();

            // treat the edges in the initial state cut as true externals
            let surface = match surface {
                HybridSurface::Esurface(esurface) => {
                    let energies_to_be_moved = esurface
                        .energies
                        .iter()
                        .filter(|edge_id| edges_in_initial_state_cut.contains(edge_id))
                        .copied()
                        .collect_vec();

                    if energies_to_be_moved.is_empty() {
                        HybridSurface::Esurface(esurface)
                    } else {
                        let new_energies = esurface
                            .energies
                            .iter()
                            .filter(|edge_id| !energies_to_be_moved.contains(edge_id))
                            .copied()
                            .collect_vec();

                        let mut new_shift = esurface.external_shift.clone();
                        for energy_to_move in energies_to_be_moved.iter() {
                            new_shift.push((*energy_to_move, 1));
                        }

                        new_shift.sort_by_key(|(edge_id, _)| *edge_id);

                        HybridSurface::Esurface(Esurface {
                            energies: new_energies,
                            external_shift: new_shift,
                            vertex_set: esurface.vertex_set,
                        })
                    }
                }
                HybridSurface::Unit(unit) => HybridSurface::Unit(unit),
                HybridSurface::Infinite(infinite) => HybridSurface::Infinite(infinite),
                HybridSurface::Hsurface(hsurface) => {
                    let positive_energies_to_be_moved = hsurface
                        .positive_energies
                        .iter()
                        .filter(|edge_id| edges_in_initial_state_cut.contains(edge_id))
                        .copied()
                        .collect_vec();

                    let negative_energies_to_be_moved = hsurface
                        .negative_energies
                        .iter()
                        .filter(|edge_id| edges_in_initial_state_cut.contains(edge_id))
                        .copied()
                        .collect_vec();

                    if positive_energies_to_be_moved.is_empty()
                        && negative_energies_to_be_moved.is_empty()
                    {
                        HybridSurface::Hsurface(hsurface)
                    } else if !positive_energies_to_be_moved.is_empty()
                        && negative_energies_to_be_moved.is_empty()
                    {
                        let new_positive_energies = hsurface
                            .positive_energies
                            .iter()
                            .filter(|edge_id| !positive_energies_to_be_moved.contains(edge_id))
                            .copied()
                            .collect_vec();

                        let mut new_shift = hsurface.external_shift.clone();

                        for positive_energy_to_move in positive_energies_to_be_moved.iter() {
                            new_shift.push((*positive_energy_to_move, 1));
                        }

                        new_shift.sort_by_key(|(edge_id, _)| *edge_id);

                        HybridSurface::Hsurface(Hsurface {
                            positive_energies: new_positive_energies,
                            negative_energies: hsurface.negative_energies.clone(),
                            external_shift: new_shift,
                            vertex_set: hsurface.vertex_set,
                        })
                    } else if !negative_energies_to_be_moved.is_empty()
                        && positive_energies_to_be_moved.is_empty()
                    {
                        let new_negative_energies = hsurface
                            .negative_energies
                            .iter()
                            .filter(|edge_id| !negative_energies_to_be_moved.contains(edge_id))
                            .copied()
                            .collect_vec();

                        let mut new_shift = hsurface.external_shift.clone();

                        for negative_energy_to_move in negative_energies_to_be_moved.iter() {
                            new_shift.push((*negative_energy_to_move, -1));
                        }

                        new_shift.sort_by_key(|(edge_id, _)| *edge_id);

                        if new_negative_energies.is_empty() {
                            HybridSurface::Esurface(Esurface {
                                energies: hsurface.positive_energies.clone(),
                                external_shift: new_shift,
                                vertex_set: hsurface.vertex_set,
                            })
                        } else {
                            HybridSurface::Hsurface(Hsurface {
                                positive_energies: hsurface.positive_energies.clone(),
                                negative_energies: new_negative_energies,
                                external_shift: new_shift,
                                vertex_set: hsurface.vertex_set,
                            })
                        }
                    } else {
                        unreachable!()
                    }
                }
            };

            let surface_id = match surface {
                HybridSurface::Esurface(mut esurface) => {
                    if let Some(shift_rewrite) = canonize_esurface {
                        esurface.canonicalize_shift(shift_rewrite);
                    }
                    let option_esurface_id = generator_cache
                        .esurface_cache
                        .position(|val| *val == esurface);

                    let esurface_id = match option_esurface_id {
                        Some(esurface_id) => esurface_id,
                        None => {
                            generator_cache.esurface_cache.push(esurface);
                            Into::<EsurfaceID>::into(generator_cache.esurface_cache.len() - 1)
                        }
                    };

                    HybridSurfaceID::Esurface(esurface_id)
                }
                HybridSurface::Hsurface(hsurface) => {
                    let option_hsurface_id = generator_cache
                        .hsurface_cache
                        .position(|val| val == &hsurface);

                    let hsurface_id = match option_hsurface_id {
                        Some(hsurface_id) => hsurface_id,
                        None => {
                            generator_cache.hsurface_cache.push(hsurface);
                            Into::<HsurfaceID>::into(generator_cache.hsurface_cache.len() - 1)
                        }
                    };

                    HybridSurfaceID::Hsurface(hsurface_id)
                }
                HybridSurface::Unit(_) => HybridSurfaceID::Unit,
                HybridSurface::Infinite(_) => HybridSurfaceID::Infinite,
            };

            (option_children, surface_id)
        })
        .unzip();

    bottom_layer
        .iter()
        .zip(new_surfaces_for_tree)
        .for_each(|(&node_id, esurface_id)| {
            tree.apply_mut_closure(node_id, |data| data.insert_esurface(esurface_id))
        });

    let all_some = children_optional.iter().all(Option::is_some);
    let all_none = children_optional.iter().all(Option::is_none);

    assert!(
        all_some || all_none,
        "Some cff branches have finished earlier than others"
    );

    let children = if all_some && !all_none {
        children_optional
            .into_iter()
            .map(Option::unwrap)
            .collect_vec()
    } else {
        return None;
    };

    bottom_layer
        .iter()
        .zip(children)
        .for_each(|(&node_id, children)| {
            children.into_iter().for_each(|child| {
                let child_node = GenerationData {
                    graph: child,
                    surface_id: None,
                };

                tree.insert_node(node_id, child_node);
            });
        });
    Some(())
}

#[cfg(test)]
mod tests_cff {
    use std::{ops::Range, vec};

    use ahash::HashMap;

    use linnet::half_edge::{
        builder::HedgeGraphBuilder, involution::Flow, nodestore::NodeStorageVec,
    };
    use symbolica::{
        evaluate::{ExpressionEvaluator, FunctionMap, OptimizationSettings},
        parse, symbol,
    };
    use utils::FloatLike;

    use crate::{
        cff::cff_graph::CFFEdgeType,
        momentum::{FourMomentum, ThreeMomentum},
        settings::global::OrientationPattern,
        utils::{
            self, F, RefDefault, external_energy_atom_from_index, ose_atom_from_index,
            test_utils::dummy_hedge_graph,
        },
    };

    use super::*;

    // helper function to make a symbolica evaluator
    impl CFFExpression<OrientationID> {
        fn quick_symbolica_evaluator(
            &self,
            external_range: Range<usize>,
            virtual_range: Range<usize>,
        ) -> ExpressionEvaluator<F<f64>> {
            let expression_atom_no_energy_sub = self.to_atom(OrientationPattern::default());
            let num_energies = external_range.end.max(virtual_range.end);
            let mut expression_atom = self
                .surfaces
                .substitute_energies(&expression_atom_no_energy_sub, &[]);
            for edge_id in 0..num_energies {
                let edge_id = EdgeIndex::from(edge_id);
                expression_atom = expression_atom
                    .replace(ose_atom_from_index(edge_id))
                    .with(external_energy_atom_from_index(edge_id));
            }

            let params = (0..num_energies)
                .map(|i| external_energy_atom_from_index(EdgeIndex::from(i)))
                .collect_vec();

            let function_map = FunctionMap::new();

            let mut tree = expression_atom
                .as_view()
                .to_evaluation_tree(&function_map, &params)
                .unwrap();

            tree.horner_scheme();
            tree.common_subexpression_elimination();
            tree.linearize(&OptimizationSettings::default())
                .map_coeff(&|c| (&c.re).into())
        }
    }

    // helper function to do some quick tests
    #[allow(unused)]
    fn generate_orientations_for_testing(
        edges: Vec<(usize, usize)>,
        incoming_vertices: Vec<usize>,
    ) -> Vec<CFFGenerationGraph> {
        let num_edges = edges.len();
        let incoming_vertices = incoming_vertices
            .into_iter()
            .map(|v| (v, CFFEdgeType::External))
            .collect_vec();

        iterate_possible_orientations(num_edges)
            .map(|or| {
                let orientation_vector = or.into_iter().collect_vec();
                let mut new_edges = edges.clone();
                for (edge_id, edge_orientation) in orientation_vector.iter().enumerate() {
                    match edge_orientation {
                        Orientation::Default => {
                            new_edges[edge_id] = edges[edge_id];
                        }
                        Orientation::Reversed => {
                            let rotated_edge = (edges[edge_id].1, edges[edge_id].0);
                            new_edges[edge_id] = rotated_edge;
                        }
                        Orientation::Undirected => {
                            unreachable!("unexpected orientation")
                        }
                    }
                }

                CFFGenerationGraph::from_vec(new_edges, incoming_vertices.clone(), None)
            })
            .filter(|graph| !graph.has_directed_cycle_initial())
            .collect_vec()
    }

    #[allow(unused)]
    fn compute_one_loop_energy<T: FloatLike>(
        k: ThreeMomentum<F<T>>,
        p: ThreeMomentum<F<T>>,
        m: F<T>,
    ) -> F<T> {
        ((k + p).norm_squared() + &m * &m).sqrt()
    }

    #[test]
    fn test_orientation_struct() {
        let orientations = iterate_possible_orientations(3).collect_vec();
        assert_eq!(orientations.len(), 8);

        let orientation1 = orientations[0].into_iter().collect_vec();
        assert_eq!(
            orientation1,
            vec![
                Orientation::Default,
                Orientation::Default,
                Orientation::Default
            ]
        );

        let orientation2 = orientations[1].into_iter().collect_vec();
        assert_eq!(
            orientation2,
            vec![
                Orientation::Reversed,
                Orientation::Default,
                Orientation::Default
            ]
        );

        let orientation3 = orientations[2].into_iter().collect_vec();
        assert_eq!(
            orientation3,
            vec![
                Orientation::Default,
                Orientation::Reversed,
                Orientation::Default
            ]
        );

        let orientation4 = orientations[3].into_iter().collect_vec();
        assert_eq!(
            orientation4,
            vec![
                Orientation::Reversed,
                Orientation::Reversed,
                Orientation::Default
            ]
        );

        let orientation5 = orientations[4].into_iter().collect_vec();
        assert_eq!(
            orientation5,
            vec![
                Orientation::Default,
                Orientation::Default,
                Orientation::Reversed
            ]
        );

        let orientation6 = orientations[5].into_iter().collect_vec();
        assert_eq!(
            orientation6,
            vec![
                Orientation::Reversed,
                Orientation::Default,
                Orientation::Reversed
            ]
        );

        let orientation7 = orientations[6].into_iter().collect_vec();
        assert_eq!(
            orientation7,
            vec![
                Orientation::Default,
                Orientation::Reversed,
                Orientation::Reversed
            ]
        );

        let orientation8 = orientations[7].into_iter().collect_vec();
        assert_eq!(
            orientation8,
            vec![
                Orientation::Reversed,
                Orientation::Reversed,
                Orientation::Reversed
            ]
        );
    }

    #[test]
    fn fishnet2b2() {
        let edges = vec![
            (0, 1),
            (1, 2),
            (3, 4),
            (4, 5),
            (6, 7),
            (7, 8),
            (0, 3),
            (1, 4),
            (2, 5),
            (3, 6),
            (4, 7),
            (5, 8),
        ];

        let incoming_vertices = vec![0, 2, 6, 8];

        let dep_mom = EdgeIndex::from(3);
        let dep_mom_expr = vec![
            (EdgeIndex::from(0), -1),
            (EdgeIndex::from(1), -1),
            (EdgeIndex::from(2), -1),
        ];

        let shift_rewrite = ShiftRewrite {
            dependent_momentum: dep_mom,
            dependent_momentum_expr: dep_mom_expr,
        };

        let orientations = generate_orientations_for_testing(edges, incoming_vertices);

        // get time before cff generation
        let start = std::time::Instant::now();

        let mut surface_cache = SurfaceCache::new();

        let _cff = generate_cff_from_orientations::<OrientationID>(
            orientations,
            &mut surface_cache,
            &[],
            &Some(shift_rewrite),
        )
        .unwrap();

        let finish = std::time::Instant::now();
        println!("time to generate cff: {:?}", finish - start);
    }

    #[test]
    fn cube() {
        let edges = vec![
            (0, 1),
            (1, 3),
            (3, 2),
            (2, 0),
            (4, 5),
            (5, 7),
            (7, 6),
            (6, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];

        let mut external_data = HashMap::default();
        for v in 0..8 {
            external_data.insert(v, vec![12 + v]);
        }

        let mut position_map = HashMap::default();
        for i in 0..edges.len() {
            position_map.insert(i, i);
        }

        let dep_mom = EdgeIndex::from(7);
        let dep_mom_expr = (0..7).map(|i| (EdgeIndex::from(i), -1)).collect();

        let shift_rewrite = ShiftRewrite {
            dependent_momentum: dep_mom,
            dependent_momentum_expr: dep_mom_expr,
        };

        let incoming_vertices = vec![0, 1, 2, 3, 4, 5, 6, 7];

        let orientations = generate_orientations_for_testing(edges, incoming_vertices);

        // get time before cff generation
        let _start = std::time::Instant::now();

        let mut surface_cache = SurfaceCache::new();

        let _cff = generate_cff_from_orientations::<OrientationID>(
            orientations,
            &mut surface_cache,
            &[],
            &Some(shift_rewrite),
        )
        .unwrap();

        let _finish = std::time::Instant::now();
    }

    fn proper_atom(graph: &HedgeGraph<(), ()>) -> Atom {
        let cff = generate_cff_expression(graph, &None, &[], &[]).unwrap();

        let mut cff_atom = cff.to_atom(OrientationPattern::default());
        cff_atom = cff.surfaces.substitute_energies(&cff_atom, &[]);
        let inverse_energy_product =
            get_cff_inverse_energy_product_impl(graph, &graph.full_graph(), &[]);

        cff_atom *= inverse_energy_product;
        cff_atom
    }

    #[test]
    fn test_dot_trick_bubble() {
        let mut dotted_topology_builder = HedgeGraphBuilder::new();
        let dotted_nodes = (0..3)
            .map(|_| dotted_topology_builder.add_node(()))
            .collect_vec();

        dotted_topology_builder.add_edge(dotted_nodes[0], dotted_nodes[1], (), false);
        dotted_topology_builder.add_edge(dotted_nodes[1], dotted_nodes[2], (), false);
        dotted_topology_builder.add_edge(dotted_nodes[2], dotted_nodes[0], (), false);
        let dotted_topology = dotted_topology_builder.build();

        let mut dotted_cff_atom = proper_atom(&dotted_topology);
        dotted_cff_atom = dotted_cff_atom
            .replace(parse!("OSE(2)"))
            .with(parse!("OSE(1)"));

        let mut topology_builder = HedgeGraphBuilder::new();
        let nodes = (0..2).map(|_| topology_builder.add_node(())).collect_vec();

        topology_builder.add_edge(nodes[0], nodes[1], (), false);
        topology_builder.add_edge(nodes[0], nodes[1], (), false);
        let toplogy = topology_builder.build();
        let mut cff_atom = proper_atom(&toplogy);

        cff_atom = cff_atom.replace(parse!("OSE(1)")).with(parse!("OSE1"));
        cff_atom = cff_atom.derivative(symbol!("OSE1"));
        cff_atom = cff_atom.replace(parse!("OSE1")).with(parse!("OSE(1)")) / parse!("2*OSE(1)");

        let diff = (&cff_atom - &dotted_cff_atom).expand();

        println!("cff_atom: {}", cff_atom.expand());
        println!("dotted_cff_atom: {}", dotted_cff_atom.expand());

        println!("diff: {}", diff);
    }

    #[test]
    fn test_dot_trick_amg() {
        let mut dotted_topology_builder = HedgeGraphBuilder::new();
        let dotted_nodes = (0..5)
            .map(|_| dotted_topology_builder.add_node(()))
            .collect_vec();

        dotted_topology_builder.add_edge(dotted_nodes[0], dotted_nodes[3], (), false);
        dotted_topology_builder.add_edge(dotted_nodes[0], dotted_nodes[2], (), false);
        dotted_topology_builder.add_edge(dotted_nodes[0], dotted_nodes[1], (), false);
        dotted_topology_builder.add_edge(dotted_nodes[3], dotted_nodes[1], (), false);
        dotted_topology_builder.add_edge(dotted_nodes[1], dotted_nodes[2], (), false);
        dotted_topology_builder.add_edge(dotted_nodes[3], dotted_nodes[4], (), false);
        dotted_topology_builder.add_edge(dotted_nodes[4], dotted_nodes[2], (), false);
        let dotted_topology = dotted_topology_builder.build();

        let mut dotted_cff_atom = proper_atom(&dotted_topology);

        dotted_cff_atom = dotted_cff_atom
            .replace(parse!("OSE(6)"))
            .with(parse!("OSE(5)"));

        let mut topology_builder = HedgeGraphBuilder::new();
        let _nodes = (0..4).map(|_| topology_builder.add_node(())).collect_vec();

        topology_builder.add_edge(dotted_nodes[0], dotted_nodes[3], (), false);
        topology_builder.add_edge(dotted_nodes[0], dotted_nodes[2], (), false);
        topology_builder.add_edge(dotted_nodes[0], dotted_nodes[1], (), false);
        topology_builder.add_edge(dotted_nodes[3], dotted_nodes[1], (), false);
        topology_builder.add_edge(dotted_nodes[1], dotted_nodes[2], (), false);
        topology_builder.add_edge(dotted_nodes[3], dotted_nodes[2], (), false);
        let topology = topology_builder.build();

        let mut cff_atom = proper_atom(&topology);
        cff_atom = cff_atom.replace(parse!("OSE(5)")).with(parse!("OSE5"));
        cff_atom = cff_atom.derivative(symbol!("OSE5"));
        cff_atom = cff_atom.replace(parse!("OSE5")).with(parse!("OSE(5)")) / parse!("2*OSE(5)");

        let diff = (cff_atom + dotted_cff_atom).expand();
        //.replace(function!(GS.ose, W_.x_))
        //.with(parse!("E"))
        //.expand();

        println!("diff: {}", diff);
    }

    #[test]
    fn test_cff_generation_triangle() {
        let triangle = vec![(2, 0), (0, 1), (1, 2)];

        let incoming_vertices = vec![0, 1, 2];
        let orientations = generate_orientations_for_testing(triangle, incoming_vertices);
        assert_eq!(orientations.len(), 6);

        let dep_mom = EdgeIndex::from(2);
        let dep_mom_expr = vec![(EdgeIndex::from(0), -1), (EdgeIndex::from(1), -1)];

        let shift_rewrite = Some(ShiftRewrite {
            dependent_momentum: dep_mom,
            dependent_momentum_expr: dep_mom_expr,
        });

        let mut surface_cache = SurfaceCache::new();

        let cff = generate_cff_from_orientations(
            orientations,
            &mut surface_cache,
            &[],
            &shift_rewrite.clone(),
        )
        .unwrap();
        assert_eq!(
            cff.surfaces.esurface_cache.len(),
            6,
            "too many esurfaces: {:#?}",
            cff.surfaces.esurface_cache,
        );

        let p1 = FourMomentum::from_args(F(1.), F(3.), F(4.), F(5.));
        let p2 = FourMomentum::from_args(F(1.), F(6.), F(7.), F(8.));
        let p3 = -p1 - p2;
        let zero = FourMomentum::from_args(F(0.), F(0.), F(0.), F(0.));
        let m = F(0.);

        let k = ThreeMomentum::new(F(1.), F(2.), F(3.));

        let virtual_energy_cache = [
            compute_one_loop_energy(k, zero.spatial, m),
            compute_one_loop_energy(k, p1.spatial, m),
            compute_one_loop_energy(k, p1.spatial + p2.spatial, m),
        ];

        let external_energy_cache = [p1.temporal.value, p2.temporal.value, p3.temporal.value];

        // combine the virtual and external energies
        let mut energy_cache = external_energy_cache.to_vec();
        energy_cache.extend(virtual_energy_cache);

        let energy_cache = dummy_hedge_graph(6)
            .new_edgevec_from_iter(energy_cache)
            .unwrap();

        let energy_prefactor = virtual_energy_cache
            .iter()
            .map(|e| (F(2.) * e).inv())
            .reduce(|acc, x| acc * x)
            .unwrap();

        let mut evaluator = cff.quick_symbolica_evaluator(0..3, 3..6);

        let cff_res: F<f64> = energy_prefactor
            * evaluator.evaluate_single(energy_cache.clone().as_ref())
            * F((2. * std::f64::consts::PI).powi(-3));

        let target_res = F(6.333_549_225_536_17e-9_f64);
        let absolute_error = cff_res - target_res;
        let relative_error = absolute_error.abs() / cff_res.abs();

        assert!(
            relative_error.abs() < F(1.0e-15),
            "relative error: {:+e} (ground truth: {:+e} vs reproduced: {:+e})",
            relative_error,
            target_res,
            cff_res
        );

        // test cff from hedge graph
        let mut triangle_hedge_graph_builder = HedgeGraphBuilder::new();

        let nodes = (0..3)
            .map(|_| triangle_hedge_graph_builder.add_node(()))
            .collect_vec();

        for node in nodes.clone() {
            triangle_hedge_graph_builder.add_external_edge(
                node,
                (),
                Orientation::Undirected,
                Flow::Sink,
            );
        }

        triangle_hedge_graph_builder.add_edge(nodes[2], nodes[0], (), Orientation::Undirected);
        triangle_hedge_graph_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);
        triangle_hedge_graph_builder.add_edge(nodes[1], nodes[2], (), Orientation::Undirected);

        let triangle_hedge_graph: HedgeGraph<(), (), ()> =
            triangle_hedge_graph_builder.build::<NodeStorageVec<()>>();

        let cff_hedge =
            generate_cff_expression(&triangle_hedge_graph, &shift_rewrite, &[], &[]).unwrap();
        let mut cff_hedge_evaluator = cff_hedge.quick_symbolica_evaluator(0..3, 3..6);

        let cff_res: F<f64> = energy_prefactor
            * cff_hedge_evaluator.evaluate_single(energy_cache.as_ref())
            * F((2. * std::f64::consts::PI).powi(-3));

        let target_res = F(6.333_549_225_536_17e-9_f64);
        let absolute_error = cff_res - target_res;
        let relative_error = absolute_error.abs() / cff_res.abs();

        assert!(
            relative_error.abs() < F(1.0e-15),
            "relative error: {:+e} (ground truth: {:+e} vs reproduced: {:+e})",
            relative_error,
            target_res,
            cff_res
        );
    }

    mod failing {
        use super::*;

        #[test]
        fn test_cff_test_double_triangle() {
            let double_triangle_edges = vec![(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)];
            let incoming_vertices = vec![0, 3];

            let orientations =
                generate_orientations_for_testing(double_triangle_edges, incoming_vertices);

            let dep_mom = EdgeIndex::from(1);
            let dep_mom_expr = vec![(EdgeIndex::from(0), -1)];

            let shift_rewrite = Some(ShiftRewrite {
                dependent_momentum: dep_mom,
                dependent_momentum_expr: dep_mom_expr,
            });

            let mut surface_cache = SurfaceCache::new();

            let cff = generate_cff_from_orientations(
                orientations,
                &mut surface_cache,
                &[],
                &shift_rewrite,
            )
            .unwrap();

            let q = FourMomentum::from_args(F(1.), F(2.), F(3.), F(4.));
            let zero = FourMomentum::from_args(F(0.), F(0.), F(0.), F(0.));

            let k = ThreeMomentum::new(F(6.), F(23.), F(9.));
            let l = ThreeMomentum::new(F(3.), F(12.), F(34.));

            let m = F::from_f64(0.);

            let virtual_energy_cache = [
                compute_one_loop_energy(k, zero.spatial, m),
                compute_one_loop_energy(q.spatial - k, zero.spatial, m),
                compute_one_loop_energy(k - l, zero.spatial, m),
                compute_one_loop_energy(l, zero.spatial, m),
                compute_one_loop_energy(q.spatial - l, zero.spatial, m),
            ];

            let external_energy_cache = [q.temporal.value, -q.temporal.value];

            let mut energy_cache = external_energy_cache.to_vec();
            energy_cache.extend(virtual_energy_cache);

            let energy_cache = dummy_hedge_graph(energy_cache.len())
                .new_edgevec_from_iter(energy_cache)
                .unwrap();

            let energy_prefactor = virtual_energy_cache
                .iter()
                .map(|e| (F(2.) * e).inv())
                .reduce(|acc, x| acc * x)
                .unwrap();

            let mut evaluator = cff.quick_symbolica_evaluator(0..2, 2..7);

            let cff_res =
                energy_prefactor * evaluator.evaluate_single(energy_cache.clone().as_ref());

            let target = F(1.0794792137096797e-13);
            let absolute_error = cff_res - target;
            let relative_error = absolute_error / cff_res;

            assert!(
                relative_error.abs() < F(1.0e-15),
                "relative error: {:+e}, target: {:+e}, result: {:+e}",
                relative_error,
                target,
                cff_res
            );

            let mut hedge_double_triangle_builder = HedgeGraphBuilder::new();
            let nodes = (0..4)
                .map(|_| hedge_double_triangle_builder.add_node(()))
                .collect_vec();

            hedge_double_triangle_builder.add_external_edge(
                nodes[0],
                (),
                Orientation::Undirected,
                Flow::Sink,
            );
            hedge_double_triangle_builder.add_external_edge(
                nodes[3],
                (),
                Orientation::Undirected,
                Flow::Sink,
            );

            hedge_double_triangle_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);
            hedge_double_triangle_builder.add_edge(nodes[0], nodes[2], (), Orientation::Undirected);
            hedge_double_triangle_builder.add_edge(nodes[1], nodes[2], (), Orientation::Undirected);
            hedge_double_triangle_builder.add_edge(nodes[1], nodes[3], (), Orientation::Undirected);
            hedge_double_triangle_builder.add_edge(nodes[2], nodes[3], (), Orientation::Undirected);

            let hedge_double_traingle: HedgeGraph<(), (), ()> =
                hedge_double_triangle_builder.build::<NodeStorageVec<()>>();
            let cff_hedge =
                generate_cff_expression(&hedge_double_traingle, &shift_rewrite, &[], &[]).unwrap();
            let mut cff_hedge_evaluator = cff_hedge.quick_symbolica_evaluator(0..2, 2..7);
            let cff_res =
                energy_prefactor * cff_hedge_evaluator.evaluate_single(energy_cache.as_ref());

            let target = F(1.0794792137096797e-13);
            let absolute_error = cff_res - target;
            let relative_error = absolute_error / cff_res;

            assert!(
                relative_error.abs() < F(1.0e-15),
                "relative error: {:+e}, target: {:+e}, result: {:+e}",
                relative_error,
                target,
                cff_res
            );

            let node_3 = hedge_double_traingle.iter_crown(nodes[3]).into();
            let node_0 = hedge_double_traingle.iter_crown(nodes[0]).into();

            let cuts = hedge_double_traingle.all_cuts(node_3, node_0);
            let mut num_with_6_ors = 0;
            let mut num_with_4_ors = 0;
            assert_eq!(cuts.len(), 4);
            for (_, cut, _) in &cuts {
                let orientations = get_orientations_with_cut(&hedge_double_traingle, cut);
                if orientations.len() == 4 {
                    num_with_4_ors += 1
                } else if orientations.len() == 6 {
                    num_with_6_ors += 1
                }
            }

            assert_eq!(num_with_4_ors, 2);
            assert_eq!(num_with_6_ors, 2);
        }

        #[test]
        fn test_cff_tbt() {
            let tbt_edges = vec![
                (0, 1),
                (2, 0),
                (1, 2),
                (1, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (5, 4),
            ];

            let incoming_vertices = vec![0, 5];

            let dep_mom = EdgeIndex::from(1);
            let dep_mom_expr = vec![(EdgeIndex::from(0), -1)];

            let shift_rewrite = Some(ShiftRewrite {
                dependent_momentum: dep_mom,
                dependent_momentum_expr: dep_mom_expr,
            });
            let mut surface_cache = SurfaceCache::new();
            let orientataions = generate_orientations_for_testing(tbt_edges, incoming_vertices);
            let cff = generate_cff_from_orientations(
                orientataions,
                &mut surface_cache,
                &[],
                &shift_rewrite,
            )
            .unwrap();

            let q = FourMomentum::from_args(F(1.0), F(2.0), F(3.0), F(4.0));
            let zero_vector = q.default();

            let p0 = q;
            let p5 = -q;

            let k = ThreeMomentum::new(F(6.), F(23.), F(9.));
            let l = ThreeMomentum::new(F(3.), F(12.), F(34.));
            let m = ThreeMomentum::new(F(7.), F(24.), F(1.));

            let mass = F(0.);

            let energies_cache = [
                p0.temporal.value,
                p5.temporal.value,
                compute_one_loop_energy(k, zero_vector.spatial, mass),
                compute_one_loop_energy(k - q.spatial, zero_vector.spatial, mass),
                compute_one_loop_energy(k - l, zero_vector.spatial, mass),
                compute_one_loop_energy(l, zero_vector.spatial, mass),
                compute_one_loop_energy(q.spatial - l, zero_vector.spatial, mass),
                compute_one_loop_energy(l - m, zero_vector.spatial, mass),
                compute_one_loop_energy(m, zero_vector.spatial, mass),
                compute_one_loop_energy(m - q.spatial, zero_vector.spatial, mass),
            ];

            let virtual_energy_cache = energies_cache[2..].to_vec();

            let energy_prefactor = virtual_energy_cache
                .iter()
                .map(|e| (F(2.) * e).inv())
                .reduce(|acc, x| acc * x)
                .unwrap();

            let energies_cache = dummy_hedge_graph(energies_cache.len())
                .new_edgevec_from_iter(energies_cache)
                .unwrap();

            let mut evaluator = cff.quick_symbolica_evaluator(0..2, 2..10);

            let res = evaluator.evaluate_single(energies_cache.clone().as_ref()) * energy_prefactor;

            let absolute_error = res - F(1.2625322619777278e-21);
            let relative_error = absolute_error / res;
            assert!(
                relative_error.abs() < F(1.0e-15),
                "relative error: {:+e}",
                relative_error
            );

            let mut tbt_hedge_builder = HedgeGraphBuilder::new();
            let nodes = (0..6).map(|_| tbt_hedge_builder.add_node(())).collect_vec();
            tbt_hedge_builder.add_external_edge(nodes[0], (), Orientation::Undirected, Flow::Sink);
            tbt_hedge_builder.add_external_edge(nodes[5], (), Orientation::Undirected, Flow::Sink);

            tbt_hedge_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);
            tbt_hedge_builder.add_edge(nodes[2], nodes[0], (), Orientation::Undirected);
            tbt_hedge_builder.add_edge(nodes[1], nodes[2], (), Orientation::Undirected);
            tbt_hedge_builder.add_edge(nodes[1], nodes[3], (), Orientation::Undirected);
            tbt_hedge_builder.add_edge(nodes[2], nodes[4], (), Orientation::Undirected);
            tbt_hedge_builder.add_edge(nodes[3], nodes[4], (), Orientation::Undirected);
            tbt_hedge_builder.add_edge(nodes[3], nodes[5], (), Orientation::Undirected);
            tbt_hedge_builder.add_edge(nodes[5], nodes[4], (), Orientation::Undirected);

            let tbt_hedge: HedgeGraph<(), (), ()> = tbt_hedge_builder.build::<NodeStorageVec<()>>();
            let cff_hedge = generate_cff_expression(&tbt_hedge, &shift_rewrite, &[], &[]).unwrap();

            let mut cff_hedge_evaluator = cff_hedge.quick_symbolica_evaluator(0..2, 2..10);
            let res =
                cff_hedge_evaluator.evaluate_single(energies_cache.as_ref()) * energy_prefactor;

            let absolute_error = res - F(1.2625322619777278e-21);
            let relative_error = absolute_error / res;
            assert!(
                relative_error.abs() < F(1.0e-15),
                "relative error: {:+e}",
                relative_error
            );

            let node_0 = tbt_hedge.iter_crown(nodes[0]).into();
            let node_5 = tbt_hedge.iter_crown(nodes[5]).into();

            let cuts = tbt_hedge.all_cuts(node_0, node_5).clone();
            assert_eq!(cuts.len(), 9);
            let mut num_with_24 = 0;
            let mut num_with_16 = 0;
            let mut num_with_42 = 0;
            let mut num_with_36 = 0;
            for (_, cut, _) in cuts.iter() {
                let orientations = get_orientations_with_cut(&tbt_hedge, cut);
                if orientations.len() == 24 {
                    num_with_24 += 1;
                }
                if orientations.len() == 16 {
                    num_with_16 += 1;
                }
                if orientations.len() == 42 {
                    num_with_42 += 1
                }
                if orientations.len() == 36 {
                    num_with_36 += 1;
                }
            }

            assert_eq!(num_with_24, 4);
            assert_eq!(num_with_16, 2);
            assert_eq!(num_with_42, 2);
            assert_eq!(num_with_36, 1);
        }
=======

    fn get_all_replacements_in_lmb(
        &self,
        cut_edges: &[EdgeIndex],
        lmb: &LoopMomentumBasis,
    ) -> Vec<Replacement> {
        self.replacements_with(
            |edge| {
                if cut_edges.contains(&edge) {
                    crate::utils::cut_energy(edge)
                } else {
                    crate::utils::ose_atom_from_index(edge)
                }
            },
            |edge| {
                lmb.edge_signatures[edge].external.iter_enumerated().fold(
                    Atom::Zero,
                    |sum, (external_index, sign)| {
                        let energy = crate::utils::external_energy_atom_from_index(
                            lmb.ext_edges[external_index],
                        );
                        match sign {
                            crate::momentum::SignOrZero::Zero => sum,
                            crate::momentum::SignOrZero::Plus => sum + energy,
                            crate::momentum::SignOrZero::Minus => sum - energy,
                        }
                    },
                )
            },
        )
>>>>>>> xzuonylm 3e5f0ddf (rebased revision)
    }
}
