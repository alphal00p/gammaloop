//! Resolution of user-facing graph sampling channel selections.
//!
//! Settings deliberately keep selectors as strings so that TOML and Python
//! have the same compact interface.  This module turns those strings into a
//! typed, deterministic selection after the graph name is known.  Numerical
//! channel construction is intentionally kept out of this layer.

use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
    str::FromStr,
};

use crate::settings::runtime::{
    HFunctionSettings, SamplingChannelDefinition, SamplingChannelSelection, SamplingRadialProfile,
};
use color_eyre::eyre::{Result, WrapErr, eyre};
use serde::{Deserialize, Serialize};
use symbolica::{atom::Atom, symbol, try_parse};

use super::sampling_context::{
    PreparedLUHost, SamplingMapContext, SamplingProposalKey, SamplingProposalPolicies,
};
use super::sampling_maps::{SamplingEvaluationError, combine_contracts};
use super::{
    ImplicitSurfaceRadialMap, SamplingExpressionEvaluator, SamplingMapAffine, SamplingMapComponent,
    SamplingMapComposition, SamplingMapContract, SamplingMapDefinition, SamplingMapEmbedding,
    SamplingMapEvaluation, SamplingMapKernel, SamplingPartition, SamplingPartitionMode,
    SamplingScoreFunction, SamplingSupport, SharedEnergyJointMap, SurfaceRadialMap,
};
use crate::integrands::evaluation::EvaluationMetaData;
use crate::momentum::sample::{LoopMomenta, MomentumSample};
use crate::settings::runtime::ParameterizationSettings;
use crate::settings::runtime::kinematic::Externals;
use crate::utils::newton_solver::RadialRootDiagnostics;
use crate::utils::{F, FloatLike};
use crate::{DependentMomentaConstructor, momentum::ThreeMomentum};

/// Precision-independent proxy, radial-profile and joint programs in canonical
/// order. Entries reuse the one neutral joint program compiled for the catalogue.
pub(crate) type SamplingChannelPrograms = (
    Option<SamplingExpressionEvaluator>,
    Option<SamplingExpressionEvaluator>,
    Option<SamplingExpressionEvaluator>,
);

/// Built-in selectors understood by the channel catalogue.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SamplingChannelPreset {
    /// Every admissible loop-momentum basis channel.
    Lmb,
    /// The existing optimized LMB heuristic and its soft-coverage audit.
    OptimizedLmb,
    /// Surface-aware channels together with the ordinary coverage channels.
    Surfaces,
}

impl SamplingChannelPreset {
    pub const LMB: &'static str = "auto:lmb";
    pub const OPTIMIZED_LMB: &'static str = "auto:optimized_lmb";
    pub const SURFACES: &'static str = "auto:surfaces";

    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Lmb => Self::LMB,
            Self::OptimizedLmb => Self::OPTIMIZED_LMB,
            Self::Surfaces => Self::SURFACES,
        }
    }
}

/// One normalized selector in a graph's channel selection.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum SamplingChannelSelector {
    Preset(SamplingChannelPreset),
    Named(String),
}

impl SamplingChannelSelector {
    pub fn parse(source: &str) -> Result<Self, SamplingSelectionError> {
        let source = source.trim();
        if source.is_empty() {
            return Err(SamplingSelectionError::EmptySelector);
        }
        match source {
            SamplingChannelPreset::LMB => Ok(Self::Preset(SamplingChannelPreset::Lmb)),
            SamplingChannelPreset::OPTIMIZED_LMB => {
                Ok(Self::Preset(SamplingChannelPreset::OptimizedLmb))
            }
            SamplingChannelPreset::SURFACES => Ok(Self::Preset(SamplingChannelPreset::Surfaces)),
            value if value.starts_with("auto:") => {
                Err(SamplingSelectionError::UnknownPreset(value.to_owned()))
            }
            value => Ok(Self::Named(value.to_owned())),
        }
    }

    pub const fn preset(&self) -> Option<SamplingChannelPreset> {
        match self {
            Self::Preset(preset) => Some(*preset),
            Self::Named(_) => None,
        }
    }

    pub fn name(&self) -> Option<&str> {
        match self {
            Self::Named(name) => Some(name),
            Self::Preset(_) => None,
        }
    }
}

impl fmt::Display for SamplingChannelSelector {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Preset(preset) => formatter.write_str(preset.as_str()),
            Self::Named(name) => formatter.write_str(name),
        }
    }
}

impl FromStr for SamplingChannelSelector {
    type Err = SamplingSelectionError;

    fn from_str(source: &str) -> Result<Self, Self::Err> {
        Self::parse(source)
    }
}

/// A named channel together with its parsed Symbolica map definition.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedNamedSamplingChannel {
    pub name: String,
    pub definition: SamplingChannelDefinition,
    pub map: SamplingMapDefinition,
    pub singularity_proxy: Option<Atom>,
    /// One shared declaration-order dependency plan for host and compiler.
    pub blocks: Vec<ResolvedSamplingBlock>,
}

/// One coordinate block in declaration order. Host and compiler consume this
/// same plan; physical edge labels never stand in for coordinate positions.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResolvedSamplingBlock {
    pub target: SamplingMapDefinition,
    pub active_lmb: Vec<usize>,
    pub preceding_lmb: Vec<usize>,
    pub remaining_lmb: Vec<usize>,
}

impl ResolvedSamplingBlock {
    pub fn geometry_key(&self, parent_lmb: &[usize]) -> SamplingGeometryKey {
        (
            self.target.clone(),
            parent_lmb.to_vec(),
            self.active_lmb.clone(),
            self.preceding_lmb.clone(),
        )
    }
}

impl ResolvedNamedSamplingChannel {
    fn compile_map<T: FloatLike>(
        &self,
        context: &SamplingChannelCompileContext<T>,
        radial: Option<&SamplingExpressionEvaluator>,
    ) -> Result<CompiledSamplingMap<T>> {
        use SamplingMapDefinition as Map;
        let invalid = |error: String| {
            color_eyre::Report::new(SamplingChannelCompileError::InvalidChannel {
                channel: self.name.clone(),
                error,
            })
        };
        if !self.definition.on_cut.is_empty()
            && !self
                .blocks
                .iter()
                .any(|block| block.target.host_cut().is_some())
        {
            return Err(invalid(
                "on_cut requires an explicit physical host".to_owned(),
            ));
        }
        if let Map::Lmb(edges) = &self.map {
            return compile_lmb_map(&self.name, None, edges, context).map_err(Into::into);
        }
        let parent = &self.definition.parent_lmb;
        if parent.len() != context.n_loop_momenta || self.blocks.is_empty() {
            return Err(invalid(
                "named channel needs a complete parent and resolved coordinate blocks".to_owned(),
            ));
        }
        let mut top = &self.map;
        while let Map::AtCut { map, .. } = top {
            top = map;
        }
        let ordered = !matches!(top, Map::Product(_));
        let mut used = BTreeSet::new();
        let mut children: Vec<Box<dyn SamplingMapComponent<T>>> = Vec::new();
        let mut output_indices = Vec::new();
        let mut only_map = None;
        for block in &self.blocks {
            for edge in &block.active_lmb {
                let Some(position) = parent.iter().position(|candidate| candidate == edge) else {
                    return Err(invalid(format!(
                        "block {:?} contains edge {edge} outside parent LMB {parent:?}",
                        block.active_lmb
                    )));
                };
                if !used.insert(*edge) {
                    return Err(invalid(format!(
                        "overlapping child blocks on edge {edge}; blocks must be disjoint in parent LMB {parent:?}"
                    )));
                }
                output_indices.extend([3 * position, 3 * position + 1, 3 * position + 2]);
            }
            let mut map = match &block.target {
                Map::Lmb(_) | Map::Complement(_) => SamplingMapKernel::new(
                    block.target.clone(),
                    context.parameterization_settings.clone(),
                    context.e_cm,
                    block.active_lmb.len(),
                )
                .map(CompiledSamplingMap::Lmb)
                .map_err(|error| invalid(error.to_string()))?,
                target
                    if !target.energy_edge_sets().is_empty() && !matches!(target, Map::Cut(_)) =>
                {
                    if let Some(host) = target.host_cut() {
                        let ids = context.physical_cut_ids.get(host).ok_or_else(|| {
                            invalid(format!(
                                "target {target:?} has no graph-resolved physical host cut {host:?}"
                            ))
                        })?;
                        if !self.definition.on_cut.is_empty()
                            && !ids.iter().any(|id| self.definition.on_cut.contains(id))
                        {
                            return Err(invalid(format!(
                                "on_cut {:?} disagrees with physical host {host:?} IDs {ids:?}",
                                self.definition.on_cut
                            )));
                        }
                    } else if !self.definition.on_cut.is_empty() {
                        return Err(invalid(format!(
                            "on_cut {:?} requires an explicit physical host",
                            self.definition.on_cut
                        )));
                    }
                    context
                        .geometry_maps
                        .get(&block.geometry_key(parent))
                        .cloned()
                        .ok_or_else(|| SamplingChannelCompileError::MissingGeometry {
                            channel: self.name.clone(),
                            target: target.clone(),
                            subspace_lmb: block.active_lmb.clone(),
                        })?
                }
                unsupported => {
                    return Err(unsupported_map_error(
                        &self.name,
                        format!("{unsupported:?}"),
                        unsupported,
                    )
                    .into());
                }
            };
            if map.dimensions() != 3 * block.active_lmb.len() {
                return Err(invalid(format!(
                    "target {:?} map dimension {} disagrees with active block {:?}",
                    block.target,
                    map.dimensions(),
                    block.active_lmb
                )));
            }
            if map.contract().requires_context && (!ordered || block.preceding_lmb.is_empty()) {
                return Err(invalid(format!(
                    "conditional surface {:?} requires then(...) with declared preceding coordinates, received {:?}",
                    block.target, block.preceding_lmb
                )));
            }
            if block.target.is_phase_space() {
                match (&self.definition.radial_profile, radial) {
                    (Some(profile), Some(program)) => {
                        let host = block.target.host_cut().unwrap();
                        let order =
                            *context
                                .physical_cut_max_occurrences
                                .get(host)
                                .ok_or_else(|| {
                                    invalid(
                                        "LU h profile requires actual maximum cut residue order"
                                            .to_owned(),
                                    )
                                })?;
                        map = map
                            .with_lu_h_profile(program.clone(), profile, order)
                            .wrap_err_with(|| {
                                format!("channel {}: native LU h profile binding", self.name)
                            })?;
                    }
                    (None, None) => {}
                    _ => {
                        return Err(invalid(
                            "compiled radial-profile presence disagrees with channel metadata"
                                .to_owned(),
                        ));
                    }
                }
            }
            if self.blocks.len() == 1 {
                only_map = Some(map);
            } else {
                children.push(Box::new(map));
            }
        }
        if used.len() != parent.len() {
            return Err(invalid(format!(
                "child blocks leave parent edges {:?} uncovered",
                parent
                    .iter()
                    .filter(|edge| !used.contains(edge))
                    .collect::<Vec<_>>()
            )));
        }
        let mut map = if let Some(map) = only_map {
            map
        } else {
            let composition = if ordered {
                SamplingMapComposition::then(children)
            } else {
                SamplingMapComposition::product(children)
            };
            composition
                .and_then(|map| SamplingMapEmbedding::from_composition(map, output_indices))
                .map(CompiledSamplingMap::Embedded)
                .map_err(|error| invalid(error.to_string()))?
        };
        if parent != &context.parent_lmb {
            let frame = context.lmb_frame_maps_by_edges.get(parent).ok_or_else(|| {
                SamplingChannelCompileError::MissingLmbFrameMap {
                    channel: self.name.clone(),
                    basis_id: None,
                    edges: parent.clone(),
                    parent_lmb: context.parent_lmb.clone(),
                }
            })?;
            if frame.dimension() != map.dimensions() {
                return Err(invalid(
                    "native-parent affine routing has the wrong dimension".to_owned(),
                ));
            }
            map = CompiledSamplingMap::Affine {
                map: Box::new(map),
                frame: frame.clone(),
            };
        }
        Ok(map)
    }

    /// Resolve the ordered block partition once alongside the Symbolica AST.
    /// Inferred complements are sampled before a physical fiber. Explicit
    /// products retain independent children; a host must prove that their
    /// equations do not depend on the other blocks before binding them.
    fn resolve_blocks(&mut self) -> Result<()> {
        use SamplingMapDefinition as Map;
        #[allow(clippy::too_many_arguments)]
        fn visit(
            map: &Map,
            parent: &[usize],
            default_active: &[usize],
            active: Option<&[usize]>,
            host: Option<&[usize]>,
            side: Option<super::SamplingCutSide>,
            ordered: bool,
            composition_depth: usize,
            blocks: &mut Vec<ResolvedSamplingBlock>,
        ) -> Result<()> {
            match map {
                Map::Block { edges, map } => {
                    if active.is_some_and(|outer| {
                        outer.len() != edges.len() || outer.iter().any(|edge| !edges.contains(edge))
                    }) {
                        return Err(eyre!(
                            "nested block descriptors select different active edges"
                        ));
                    }
                    visit(
                        map,
                        parent,
                        default_active,
                        Some(edges),
                        host,
                        side,
                        ordered,
                        composition_depth,
                        blocks,
                    )
                }
                Map::AtCut { cut, map } => {
                    let start = blocks.len();
                    visit(
                        map,
                        parent,
                        default_active,
                        active,
                        Some(cut),
                        side,
                        ordered,
                        composition_depth,
                        blocks,
                    )?;
                    if !blocks[start..]
                        .iter()
                        .any(|block| !block.target.energy_edge_sets().is_empty())
                    {
                        return Err(eyre!(
                            "at_cut requires a physical target; an ordinary LMB has no cut host"
                        ));
                    }
                    Ok(())
                }
                Map::Left(inner) | Map::Right(inner) => {
                    if side.is_some() {
                        return Err(eyre!("nested side qualifiers are ambiguous"));
                    }
                    if matches!(inner.as_ref(), Map::PhaseSpace(_)) {
                        return Err(eyre!("a phase-space host has no amplitude side"));
                    }
                    let side = if matches!(map, Map::Left(_)) {
                        super::SamplingCutSide::Left
                    } else {
                        super::SamplingCutSide::Right
                    };
                    // The qualifier is carried by the outer constructor, not
                    // inferred from the target's energy edges.
                    visit(
                        inner,
                        parent,
                        default_active,
                        active,
                        host,
                        Some(side),
                        ordered,
                        composition_depth,
                        blocks,
                    )
                }
                Map::Product(children) | Map::Then(children) => {
                    if composition_depth != 0 {
                        return Err(eyre!(
                            "nested block compositions are not yet supported, including through qualifiers"
                        ));
                    }
                    if active.is_some() || side.is_some() {
                        return Err(eyre!(
                            "block/side qualifiers must identify one physical child"
                        ));
                    }
                    let is_ordered = matches!(map, Map::Then(_));
                    let mut inherited = host.map(<[usize]>::to_vec);
                    for child in children {
                        visit(
                            child,
                            parent,
                            default_active,
                            None,
                            inherited.as_deref(),
                            None,
                            is_ordered,
                            composition_depth + 1,
                            blocks,
                        )?;
                        if is_ordered
                            && let Some(block) = blocks.last()
                            && block.target.is_phase_space()
                        {
                            if inherited.is_some()
                                && inherited.as_deref() != block.target.host_cut()
                            {
                                return Err(eyre!(
                                    "multiple phase-space hosts require an explicit at_cut on every dependent target"
                                ));
                            }
                            inherited = block.target.host_cut().map(<[usize]>::to_vec);
                        }
                    }
                    Ok(())
                }
                _ => {
                    if matches!(map, Map::Intersect(_)) && map.energy_edge_sets().len() != 2 {
                        return Err(eyre!(
                            "intersect currently requires exactly two distinct surface(...) targets; qualify the whole joint block with its host or side"
                        ));
                    }
                    let physical = !map.energy_edge_sets().is_empty();
                    let requested = match map {
                        Map::Lmb(edges) | Map::Complement(edges) => {
                            if active.is_some_and(|active| active != edges) {
                                return Err(eyre!(
                                    "block LMB descriptor disagrees with its ordinary map"
                                ));
                            }
                            edges.as_slice()
                        }
                        _ => active.unwrap_or(if default_active.is_empty() && !physical {
                            parent
                        } else {
                            default_active
                        }),
                    };
                    if requested.is_empty() {
                        return Err(eyre!(
                            "physical map requires subspace_lmb or an explicit block(lmb(...),map)"
                        ));
                    }
                    let active_lmb = if physical {
                        if requested.iter().any(|edge| !parent.contains(edge))
                            || requested.iter().copied().collect::<BTreeSet<_>>().len()
                                != requested.len()
                        {
                            return Err(eyre!(
                                "active block {requested:?} must be unique and contained in parent LMB {parent:?}"
                            ));
                        }
                        parent
                            .iter()
                            .copied()
                            .filter(|edge| requested.contains(edge))
                            .collect()
                    } else {
                        requested.to_vec()
                    };
                    if matches!(map, Map::Intersect(_)) && active_lmb.len() != 1 {
                        return Err(eyre!(
                            "intersect requires one three-dimensional active LMB cycle, received {active_lmb:?}"
                        ));
                    }
                    let mut target = map.clone();
                    if let Some(side) = side {
                        target = match side {
                            super::SamplingCutSide::Left => Map::Left(Box::new(target)),
                            super::SamplingCutSide::Right => Map::Right(Box::new(target)),
                        };
                    }
                    if let Some(host) = host.filter(|_| physical) {
                        if map.is_phase_space() {
                            if map.host_cut() != Some(host) {
                                return Err(eyre!("phase-space target and explicit host disagree"));
                            }
                        } else {
                            target = Map::AtCut {
                                cut: host.to_vec(),
                                map: Box::new(target),
                            };
                        }
                    } else if side.is_some() {
                        return Err(eyre!(
                            "left/right target requires an explicit or preceding phase-space host"
                        ));
                    }
                    let preceding_lmb = if ordered {
                        blocks
                            .iter()
                            .flat_map(|block| block.active_lmb.iter().copied())
                            .collect()
                    } else {
                        Vec::new()
                    };
                    blocks.push(ResolvedSamplingBlock {
                        target,
                        active_lmb,
                        preceding_lmb,
                        remaining_lmb: Vec::new(),
                    });
                    Ok(())
                }
            }
        }
        let mut blocks = Vec::new();
        visit(
            &self.map,
            &self.definition.parent_lmb,
            &self.definition.subspace_lmb,
            None,
            None,
            None,
            true,
            0,
            &mut blocks,
        )?;
        let parent = &self.definition.parent_lmb;
        if blocks.len() == 1
            && !blocks[0].target.energy_edge_sets().is_empty()
            && blocks[0].active_lmb.len() < parent.len()
        {
            let complement = parent
                .iter()
                .copied()
                .filter(|edge| !blocks[0].active_lmb.contains(edge))
                .collect::<Vec<_>>();
            blocks[0].preceding_lmb = complement.clone();
            blocks.insert(
                0,
                ResolvedSamplingBlock {
                    target: Map::Complement(complement.clone()),
                    active_lmb: complement,
                    preceding_lmb: Vec::new(),
                    remaining_lmb: Vec::new(),
                },
            );
        }
        for block in &mut blocks {
            block.remaining_lmb = parent
                .iter()
                .copied()
                .filter(|edge| {
                    !block.active_lmb.contains(edge) && !block.preceding_lmb.contains(edge)
                })
                .collect();
        }
        self.blocks = blocks;
        Ok(())
    }
}

/// The complete selection for one graph after applying the settings' fallback.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedSamplingChannelSelection {
    pub graph_name: String,
    pub selectors: Vec<SamplingChannelSelector>,
    pub named_channels: Vec<ResolvedNamedSamplingChannel>,
}

/// One entry in the graph-local sampling catalogue. The catalogue is the
/// migration target for grid construction and evaluation. The existing setup
/// only supplies generated entries until that runtime migration is complete;
/// it is not a second production channel universe.
#[derive(Clone, Debug, PartialEq)]
pub enum SamplingCatalogueEntry {
    Lmb {
        basis_id: usize,
        edges: Vec<usize>,
        preset: SamplingChannelPreset,
    },
    /// An automatically enumerated E-surface candidate. Its geometry is
    /// prepared later from the current cut/orientation kinematics.
    Surface {
        edges: Vec<usize>,
        parent_lmb: Vec<usize>,
    },
    Named(ResolvedNamedSamplingChannel),
}

/// Graph-scoped, deterministically ordered channel catalogue.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingChannelCatalogue {
    pub graph_name: String,
    pub selectors: Vec<SamplingChannelSelector>,
    pub entries: Vec<SamplingCatalogueEntry>,
}

/// Conservative coverage facts for a resolved catalogue.
///
/// This report only makes claims that can be established from the canonical
/// catalogue itself.  In particular, an LMB entry is known to provide a
/// full-domain ordinary map, while a named graph-aware map is not assumed to
/// have full support until its compiled contract has been checked.  The soft
/// audit is deliberately elementary: it records massless loop edges that are
/// present in at least one generated LMB and whether an ordinary selected
/// channel exposes each one.  It does not claim coverage of compound soft or
/// collinear strata.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SamplingCoverageReport {
    pub has_full_domain_lmb: bool,
    pub covered_massless_loop_edges: Vec<usize>,
    pub missing_massless_loop_edges: Vec<usize>,
}

/// Stable, side-effect-free view of the resolved catalogue for diagnostics.
///
/// This report is deliberately built from [`SamplingChannelCatalogue`] rather
/// than enumerating channels independently.  CLI and Python inspection can
/// therefore show exactly the catalogue that the production sampling path
/// will compile, without constructing kinematics or running a sample.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SamplingChannelInspection {
    pub graph_name: String,
    pub selectors: Vec<String>,
    pub entries: Vec<String>,
}

/// One graph-scoped geometry identity: qualified target, complete native parent,
/// active cycles in that parent's order, and the declared prerequisite order.
/// A callback's input layout is part of its identity: the same physical target
/// may follow different blocks in different channels. Host cuts are physical
/// edge sets, never sampling channel IDs or numeric CutIds.
pub type SamplingGeometryKey = (SamplingMapDefinition, Vec<usize>, Vec<usize>, Vec<usize>);

/// Input context for graph-independent compilation of the canonical catalogue.
/// Geometry and native frame transforms are bound from the same improved
/// externals. Per-draw LU data are prepared only inside a certified ordered map;
/// no partially sampled point is presented as a complete MomentumSample.
#[derive(Clone, Debug)]
pub struct SamplingChannelCompileContext<T: FloatLike = f64> {
    pub master_graph: String,
    /// Complete output frame shared by all channel densities.
    pub parent_lmb: Vec<usize>,
    pub parameterization_settings: ParameterizationSettings,
    pub e_cm: f64,
    pub n_loop_momenta: usize,
    pub geometry_maps: BTreeMap<SamplingGeometryKey, CompiledSamplingMap<T>>,
    /// Physical cut identities, independently of the sampling catalogue.
    pub physical_cut_ids: BTreeMap<Vec<usize>, Vec<usize>>,
    /// Highest physical residue order across all equivalent active cut groups.
    pub physical_cut_max_occurrences: BTreeMap<Vec<usize>, usize>,
    pub orientation: Option<usize>,
    /// Exact generated-parent to master routing, including external shifts.
    pub lmb_frame_maps: BTreeMap<usize, SamplingMapAffine<T>>,
    pub lmb_frame_maps_by_edges: BTreeMap<Vec<usize>, SamplingMapAffine<T>>,
}

impl<T: FloatLike> SamplingChannelCompileContext<T> {
    pub fn new(
        master_graph: impl Into<String>,
        parent_lmb: Vec<usize>,
        parameterization_settings: ParameterizationSettings,
        e_cm: f64,
        n_loop_momenta: usize,
    ) -> Self {
        Self {
            master_graph: master_graph.into(),
            parent_lmb,
            parameterization_settings,
            e_cm,
            n_loop_momenta,
            geometry_maps: BTreeMap::new(),
            physical_cut_ids: BTreeMap::new(),
            physical_cut_max_occurrences: BTreeMap::new(),
            orientation: None,
            lmb_frame_maps: BTreeMap::new(),
            lmb_frame_maps_by_edges: BTreeMap::new(),
        }
    }

    /// Register one native physical block transactionally. Distinct host, parent
    /// and active-cycle choices cannot overwrite each other's geometry.
    pub fn insert_geometry_map(
        &mut self,
        target: SamplingMapDefinition,
        parent_lmb: Vec<usize>,
        subspace_lmb: Vec<usize>,
        preceding_lmb: Vec<usize>,
        map: CompiledSamplingMap<T>,
    ) -> Result<()> {
        let unique =
            |edges: &[usize]| edges.iter().copied().collect::<BTreeSet<_>>().len() == edges.len();
        let energy_sets = target.energy_edge_sets();
        if energy_sets.is_empty()
            || energy_sets
                .iter()
                .any(|edges| edges.is_empty() || !unique(edges))
            || (energy_sets.len() == 2 && subspace_lmb.len() != 1)
            || parent_lmb.len() != self.n_loop_momenta
            || !unique(&parent_lmb)
            || subspace_lmb.is_empty()
            || !unique(&subspace_lmb)
            || subspace_lmb.iter().any(|edge| !parent_lmb.contains(edge))
            || !unique(&preceding_lmb)
            || preceding_lmb
                .iter()
                .any(|edge| !parent_lmb.contains(edge) || subspace_lmb.contains(edge))
        {
            return Err(eyre!(
                "geometry map needs a supported physical target, complete unique parent and nonempty active subset (one cycle for intersect): target={target:?}, parent={parent_lmb:?}, active={subspace_lmb:?}"
            ));
        }
        if map.dimensions() != 3 * subspace_lmb.len() {
            return Err(eyre!(
                "geometry map has dimension {}, expected {} for active {:?}",
                map.dimensions(),
                3 * subspace_lmb.len(),
                subspace_lmb
            ));
        }
        self.geometry_maps
            .insert((target, parent_lmb, subspace_lmb, preceding_lmb), map);
        Ok(())
    }
}

/// The map kernels currently compilable without a graph-specific implicit
/// solver.  Physical conditional/cut maps stay in the typed catalogue until
/// their prepared kinematic context is supplied by the process layer; the
/// bounded direct-product route is compiled there when its blocks are explicit.
#[derive(Clone, Debug)]
pub enum CompiledSamplingMap<T: FloatLike = f64> {
    Lmb(SamplingMapKernel),
    Affine {
        map: Box<CompiledSamplingMap<T>>,
        frame: SamplingMapAffine<T>,
    },
    Surface(SurfaceRadialMap<T>),
    ImplicitSurface(ImplicitSurfaceRadialMap<T>),
    Joint(SharedEnergyJointMap<T>),
    Embedded(SamplingMapEmbedding<T>),
}

impl<T: FloatLike> CompiledSamplingMap<T> {
    /// Bind a per-name proposal only after cloning shared physical geometry.
    /// Transparent frame transforms preserve access to the zero-centered cut
    /// primitive; side maps never inherit the host's LU h profile.
    pub fn with_lu_h_profile(
        self,
        program: SamplingExpressionEvaluator,
        profile: &SamplingRadialProfile,
        max_occurrence: usize,
    ) -> Result<Self> {
        match self {
            Self::ImplicitSurface(map) => map
                .with_lu_h_profile(program, profile, max_occurrence)
                .map(Self::ImplicitSurface),
            Self::Affine { map, frame } => Ok(Self::Affine {
                map: Box::new(map.with_lu_h_profile(program, profile, max_occurrence)?),
                frame,
            }),
            _ => Err(eyre!(
                "LU h profile requires a graph-resolved implicit cut primitive"
            )),
        }
    }

    pub fn contract(&self) -> SamplingMapContract {
        match self {
            Self::Lmb(map) => map.contract(),
            Self::Affine { map, frame } => combine_contracts(map.contract(), frame.contract()),
            Self::Surface(map) => map.contract(),
            Self::ImplicitSurface(map) => map.contract(),
            Self::Joint(map) => map.contract(),
            Self::Embedded(map) => map.contract(),
        }
    }

    pub fn dimensions(&self) -> usize {
        match self {
            Self::Lmb(map) => map.dimensions(),
            Self::Affine { map, .. } => map.dimensions(),
            Self::Surface(map) => map.dimension(),
            Self::ImplicitSurface(map) => map.dimension(),
            Self::Joint(map) => map.dimensions(),
            Self::Embedded(map) => map.dimensions(),
        }
    }

    pub fn as_component(&self) -> &dyn SamplingMapComponent<T> {
        self
    }

    pub fn forward(&self, coordinates: &[T]) -> Result<SamplingMapEvaluation<T>> {
        self.forward_with_context(coordinates, &mut SamplingMapContext::detached(&[]))
    }

    pub fn forward_with_context(
        &self,
        coordinates: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<SamplingMapEvaluation<T>> {
        <Self as SamplingMapComponent<T>>::forward(self, coordinates, context)
    }

    pub fn inverse(&self, point: &[T]) -> Result<Option<SamplingMapEvaluation<T>>> {
        self.inverse_with_context(point, &mut SamplingMapContext::detached(&[]))
    }

    pub fn inverse_with_context(
        &self,
        point: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<SamplingMapEvaluation<T>>> {
        <Self as SamplingMapComponent<T>>::inverse(self, point, context)
    }
}

impl<T: FloatLike> SamplingMapComponent<T> for CompiledSamplingMap<T> {
    fn dimensions(&self) -> usize {
        self.dimensions()
    }

    fn output_dimensions(&self) -> usize {
        self.dimensions()
    }

    fn contract(&self) -> SamplingMapContract {
        self.contract()
    }

    fn name(&self) -> &'static str {
        match self {
            Self::Lmb(_) => "lmb",
            Self::Affine { .. } => "affine",
            Self::Surface(_) => "surface",
            Self::ImplicitSurface(_) => "implicit_surface",
            Self::Joint(_) => "joint",
            Self::Embedded(_) => "embedded",
        }
    }

    fn forward(
        &self,
        coordinates: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<SamplingMapEvaluation<T>> {
        match self {
            Self::Lmb(map) => SamplingMapComponent::forward(map, coordinates, context),
            Self::Surface(map) => SamplingMapComponent::forward(map, coordinates, context),
            Self::ImplicitSurface(map) => SamplingMapComponent::forward(map, coordinates, context),
            Self::Joint(map) => SamplingMapComponent::forward(map, coordinates, context),
            Self::Embedded(map) => SamplingMapComponent::forward(map, coordinates, context),
            Self::Affine { map, frame } => {
                let previous = context.previous;
                let lmb_evaluation = SamplingMapComponent::forward(
                    map.as_ref(),
                    coordinates,
                    &mut context.reborrow(previous, Some(0)),
                )?;
                let frame_evaluation = SamplingMapComponent::forward(
                    frame,
                    &lmb_evaluation.point,
                    &mut context.reborrow(previous, Some(1)),
                )?;
                SamplingMapEvaluation {
                    coordinates: lmb_evaluation.coordinates,
                    point: frame_evaluation.point,
                    jacobian: (F(lmb_evaluation.jacobian) * F(frame_evaluation.jacobian)).0,
                    inverse_jacobian: (F(lmb_evaluation.inverse_jacobian)
                        * F(frame_evaluation.inverse_jacobian))
                    .0,
                    residual: F(lmb_evaluation.residual)
                        .max(F(frame_evaluation.residual))
                        .0,
                    support: combine_contracts(map.contract(), frame.contract()).support,
                    diagnostics: lmb_evaluation
                        .diagnostics
                        .into_iter()
                        .chain(frame_evaluation.diagnostics)
                        .collect(),
                }
                .validate("compiled affine map")
            }
        }
    }

    fn inverse(
        &self,
        point: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<SamplingMapEvaluation<T>>> {
        match self {
            Self::Lmb(map) => SamplingMapComponent::inverse(map, point, context),
            Self::Surface(map) => SamplingMapComponent::inverse(map, point, context),
            Self::ImplicitSurface(map) => SamplingMapComponent::inverse(map, point, context),
            Self::Joint(map) => SamplingMapComponent::inverse(map, point, context),
            Self::Embedded(map) => SamplingMapComponent::inverse(map, point, context),
            Self::Affine { map, frame } => {
                let previous = context.previous;
                let frame_evaluation = frame.inverse(point, previous)?;
                let Some(lmb_evaluation) = SamplingMapComponent::inverse(
                    map.as_ref(),
                    &frame_evaluation.coordinates,
                    &mut context.reborrow(previous, Some(0)),
                )?
                else {
                    return Ok(None);
                };
                SamplingMapEvaluation {
                    coordinates: lmb_evaluation.coordinates,
                    point: point.to_vec(),
                    jacobian: (F(lmb_evaluation.jacobian) * F(frame_evaluation.jacobian)).0,
                    inverse_jacobian: (F(lmb_evaluation.inverse_jacobian)
                        * F(frame_evaluation.inverse_jacobian))
                    .0,
                    residual: F(lmb_evaluation.residual)
                        .max(F(frame_evaluation.residual))
                        .0,
                    support: combine_contracts(map.contract(), frame.contract()).support,
                    diagnostics: lmb_evaluation
                        .diagnostics
                        .into_iter()
                        .chain(frame_evaluation.diagnostics)
                        .collect(),
                }
                .validate("compiled affine inverse")
                .map(Some)
            }
        }
    }
}

/// A compiled graph channel and the master-graph raw-frame block it occupies.
#[derive(Clone, Debug)]
pub struct CompiledSamplingChannel<T: FloatLike = f64> {
    pub name: String,
    pub master_graph: String,
    pub basis_id: Option<usize>,
    pub definition: SamplingMapDefinition,
    /// Ordered master-graph edge ids for this raw coordinate block.
    pub embedded_edges: Vec<usize>,
    pub map: CompiledSamplingMap<T>,
    /// Optional eager positive score used by the singularity-proxy partition.
    /// It is compiled once from the channel metadata during warmup.
    pub singularity_proxy: Option<SamplingScoreFunction<T>>,
}

/// Canonical identifier for a compiled sampling channel.
///
/// The same id is used by the graph evaluator and by the sampling catalogue;
/// no second legacy channel-index domain is maintained.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Ord, PartialOrd, Serialize, Deserialize)]
pub struct SamplingChannelId(pub usize);

impl SamplingChannelId {
    pub const fn index(self) -> usize {
        self.0
    }
}

impl From<usize> for SamplingChannelId {
    fn from(index: usize) -> Self {
        Self(index)
    }
}

impl From<SamplingChannelId> for usize {
    fn from(index: SamplingChannelId) -> Self {
        index.0
    }
}

/// A bounded bridge from compiled graph channels to the existing graph
/// evaluator.  The bridge deliberately deals in the complete master raw
/// frame: it never silently embeds a lower-dimensional surface block or
/// reinterprets a channel in another graph's loop basis.
#[derive(Clone, Debug)]
pub struct SamplingChannelBridge<T: FloatLike = f64> {
    channels: Vec<CompiledSamplingChannel<T>>,
    dimensions: usize,
    partition_mode: SamplingPartitionMode,
    /// Explicit runtime accuracy budget; fresh numerical bridges use native sqrt(epsilon).
    relative_density_tolerance: Option<f64>,
}

/// Per-sample context supplied to conditional channel maps and to every
/// inverse-density score in a common partition.  Contexts are indexed by the
/// canonical catalogue position; they are never compacted into a second
/// channel list.  An absent context means the ordinary empty context and is
/// valid for maps which do not depend on earlier blocks. Production also borrows
/// the sole original-draw policy owner; detached numerical fixtures opt out
/// explicitly, without installing a mutable cache on a warmed bridge.
#[derive(Debug)]
pub struct SamplingChannelRuntimeContexts<'a, T: FloatLike = f64> {
    contexts: Vec<Option<Vec<T>>>,
    policies: Option<&'a mut SamplingProposalPolicies>,
    graph_id: usize,
    generating_channel: SamplingChannelId,
    prepared_lu_hosts: Vec<PreparedLUHost<T>>,
    /// None only during the initial selected operation; later inverse work
    /// cannot extend the frozen source authority exported to the physical body.
    selected_host_count: Option<usize>,
    radial_root_diagnostics: Option<&'a mut RadialRootDiagnostics>,
    detached_root_diagnostics: RadialRootDiagnostics,
}

impl<'a, T: FloatLike> SamplingChannelRuntimeContexts<'a, T> {
    /// Detached represented-geometry evaluation, without production decisions.
    pub fn new(channel_count: usize) -> Self {
        Self {
            contexts: vec![None; channel_count],
            policies: None,
            graph_id: 0,
            generating_channel: SamplingChannelId(0),
            prepared_lu_hosts: Vec::new(),
            selected_host_count: Some(0),
            radial_root_diagnostics: None,
            detached_root_diagnostics: RadialRootDiagnostics::default(),
        }
    }

    /// Bind one generating row to collecting or sealed original-draw records.
    pub(crate) fn for_draw(
        channel_count: usize,
        graph_id: usize,
        generating_channel: SamplingChannelId,
        metadata: &'a mut EvaluationMetaData,
    ) -> Self {
        Self {
            contexts: vec![None; channel_count],
            policies: Some(&mut metadata.sampling_proposal_policies),
            graph_id,
            generating_channel,
            prepared_lu_hosts: Vec::new(),
            selected_host_count: Some(0),
            radial_root_diagnostics: Some(&mut metadata.radial_root_diagnostics),
            detached_root_diagnostics: RadialRootDiagnostics::default(),
        }
    }

    fn begin_selected(&mut self, channel_id: SamplingChannelId) -> Result<()> {
        if self.policies.is_some() && self.generating_channel != channel_id {
            return Err(eyre!(
                "selected channel {channel_id:?} does not match generating row {:?}",
                self.generating_channel
            ));
        }
        self.generating_channel = channel_id;
        self.prepared_lu_hosts.clear();
        self.selected_host_count = None;
        Ok(())
    }

    fn seal_selected(&mut self) {
        self.selected_host_count = Some(self.prepared_lu_hosts.len());
    }

    fn selected_hosts(&self) -> Vec<PreparedLUHost<T>> {
        self.prepared_lu_hosts[..self.selected_host_count.unwrap_or(0)].to_vec()
    }

    pub fn channel_count(&self) -> usize {
        self.contexts.len()
    }

    pub fn set(&mut self, channel_id: SamplingChannelId, context: Vec<T>) -> Result<()> {
        let Some(slot) = self.contexts.get_mut(channel_id.0) else {
            return Err(eyre!(
                "sampling channel context id {} is out of range for {} channels",
                channel_id.0,
                self.contexts.len()
            ));
        };
        if context.iter().any(|value| !value.is_finite()) {
            return Err(eyre!(
                "sampling channel context {} contains a non-finite value",
                channel_id.0
            ));
        }
        *slot = Some(context);
        Ok(())
    }

    pub(crate) fn for_channel(
        &mut self,
        channel_id: SamplingChannelId,
    ) -> Result<SamplingMapContext<'_, T>> {
        let previous = self
            .contexts
            .get(channel_id.0)
            .ok_or_else(|| {
                eyre!(
                    "sampling channel context id {} is out of range for {} channels",
                    channel_id.0,
                    self.contexts.len()
                )
            })
            .map(|context| context.as_deref().unwrap_or(&[]))?;
        Ok(SamplingMapContext {
            previous,
            policies: self.policies.as_deref_mut(),
            key: SamplingProposalKey {
                graph_id: self.graph_id,
                generating_channel: self.generating_channel,
                target_channel: channel_id,
                block_path: Vec::new(),
            },
            prepared_lu_hosts: Some(&mut self.prepared_lu_hosts),
            radial_root_diagnostics: Some(
                self.radial_root_diagnostics
                    .as_deref_mut()
                    .unwrap_or(&mut self.detached_root_diagnostics),
            ),
            selecting_lu_hosts: self.selected_host_count.is_none(),
        })
    }
}

/// One push-forward/inverse result together with the common raw-frame
/// multichannel partition.
#[derive(Clone, Debug)]
pub struct SamplingChannelBridgeEvaluation<T: FloatLike = f64> {
    pub channel_id: SamplingChannelId,
    pub channel_name: String,
    /// The complete master raw coordinate frame passed to graph evaluation.
    pub raw_coordinates: Vec<T>,
    pub map: SamplingMapEvaluation<T>,
    pub partition: SamplingPartition<T>,
    pub(crate) prepared_lu_hosts: Vec<PreparedLUHost<T>>,
}

/// Summary of a deterministic acceptance run through a complete channel
/// bridge.  The reference integrand is a normalized isotropic Gaussian in
/// the bridge's raw master frame.  Every canonical channel is sampled with
/// equal probability and the estimator includes its exact map determinant
/// and the common inverse-density partition.  This makes the report useful
/// for acceptance tests of arbitrary graph-resolved channel catalogues,
/// including products and ordered surface compositions.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingChannelBridgeAcceptanceReport {
    pub sample_count: usize,
    pub channel_count: usize,
    pub finite_sample_count: usize,
    pub normalization: f64,
    pub normalization_stderr: f64,
    pub second_moment: f64,
    pub expected_second_moment: f64,
    pub second_moment_stderr: f64,
    pub partition_min: f64,
    pub partition_max: f64,
    pub jacobian_min: f64,
    pub jacobian_max: f64,
    /// Maximum forward/inverse coordinate residual seen by the acceptance run.
    /// A finite small value is required before a channel is trusted in a
    /// graph-level acceptance test.
    pub round_trip_residual_max: f64,
}

impl SamplingChannelBridgeAcceptanceReport {
    /// Integrate a normalized Gaussian over all channels in `bridge`.
    ///
    /// The returned normalization converges to one as `sample_count` grows.
    /// A sample count applies independently to each channel; no channel is
    /// silently omitted from the acceptance run.  The deterministic Halton
    /// points keep this suitable for reproducible unit and acceptance tests.
    pub fn normalized_gaussian(
        bridge: &SamplingChannelBridge,
        sample_count: usize,
        width: f64,
        center: &[f64],
    ) -> Result<Self> {
        if sample_count == 0 {
            return Err(eyre!(
                "sampling bridge acceptance harness needs at least one sample"
            ));
        }
        if !width.is_finite() || width <= 0.0 {
            return Err(eyre!(
                "sampling bridge acceptance Gaussian width must be positive and finite"
            ));
        }
        if center.len() != bridge.dimensions {
            return Err(eyre!(
                "sampling bridge acceptance Gaussian has dimension {}, expected {}",
                center.len(),
                bridge.dimensions
            ));
        }
        if center.iter().any(|component| !component.is_finite()) {
            return Err(eyre!(
                "sampling bridge acceptance Gaussian centre must be finite"
            ));
        }

        let channel_count = bridge.channels.len();
        if channel_count == 0 {
            return Err(eyre!(
                "sampling bridge acceptance harness needs at least one channel"
            ));
        }
        let dimension = bridge.dimensions as f64;
        let gaussian_normalization =
            (2.0 * std::f64::consts::PI * width * width).powf(-0.5 * dimension);
        let expected_second_moment = center
            .iter()
            .map(|component| component * component)
            .sum::<f64>()
            + dimension * width.powi(2);
        let mut report = Self {
            sample_count,
            channel_count,
            finite_sample_count: 0,
            normalization: 0.0,
            normalization_stderr: 0.0,
            second_moment: 0.0,
            expected_second_moment,
            second_moment_stderr: 0.0,
            partition_min: f64::INFINITY,
            partition_max: f64::NEG_INFINITY,
            jacobian_min: f64::INFINITY,
            jacobian_max: f64::NEG_INFINITY,
            round_trip_residual_max: 0.0,
        };
        let mut square_sum = 0.0;
        let mut second_moment_square_sum = 0.0;
        let total_samples = sample_count * channel_count;
        for channel_index in 0..channel_count {
            let channel_id = SamplingChannelId::from(channel_index);
            for sample in 1..=sample_count {
                let coordinates = (0..bridge.channels[channel_index].dimensions())
                    .map(|axis| bridge_halton(sample, bridge_prime(axis)))
                    .collect::<Vec<_>>();
                let evaluation = bridge.forward(channel_id, &coordinates)?;
                // Exercise the inverse on every accepted point as well.  This
                // catches a channel whose forward determinant is finite while
                // its inverse uses a different branch or frame, which would
                // invalidate the graph-level map-density estimator.
                let inverse = bridge
                    .inverse(channel_id, &evaluation.raw_coordinates)?
                    .ok_or_else(|| eyre!("accepted forward point has no selected inverse"))?;
                let coordinate_residual = inverse
                    .map
                    .coordinates
                    .iter()
                    .zip(&coordinates)
                    .map(|(recovered, original)| (recovered - original).abs())
                    .fold(inverse.map.residual, f64::max);
                if !coordinate_residual.is_finite() {
                    return Err(eyre!(
                        "sampling bridge acceptance inverse residual is non-finite for channel {} sample {}",
                        channel_index,
                        sample
                    ));
                }
                report.round_trip_residual_max =
                    report.round_trip_residual_max.max(coordinate_residual);
                let partition_sum = evaluation.partition.weights.iter().sum::<f64>();
                if !partition_sum.is_finite() {
                    return Err(eyre!(
                        "sampling bridge acceptance partition is non-finite for channel {} sample {}",
                        channel_index,
                        sample
                    ));
                }
                report.partition_min = report.partition_min.min(partition_sum);
                report.partition_max = report.partition_max.max(partition_sum);
                let radius_squared = evaluation
                    .raw_coordinates
                    .iter()
                    .zip(center)
                    .map(|(point, centre)| (point - centre).powi(2))
                    .sum::<f64>();
                let weight = gaussian_normalization
                    * (-0.5 * radius_squared / width.powi(2)).exp()
                    * evaluation.selected_factor()?
                    * channel_count as f64;
                report.jacobian_min = report.jacobian_min.min(evaluation.map.jacobian);
                report.jacobian_max = report.jacobian_max.max(evaluation.map.jacobian);
                if !weight.is_finite() {
                    continue;
                }
                report.finite_sample_count += 1;
                report.normalization += weight;
                square_sum += weight * weight;
                let raw_radius_squared = evaluation
                    .raw_coordinates
                    .iter()
                    .map(|component| component.powi(2))
                    .sum::<f64>();
                let second_moment = weight * raw_radius_squared;
                report.second_moment += second_moment;
                second_moment_square_sum += second_moment * second_moment;
            }
        }
        if report.finite_sample_count > 0 {
            report.normalization /= total_samples as f64;
            report.normalization_stderr =
                ((square_sum / total_samples as f64 - report.normalization.powi(2)).max(0.0)
                    / total_samples as f64)
                    .sqrt();
            report.second_moment /= total_samples as f64;
            report.second_moment_stderr = ((second_moment_square_sum / total_samples as f64
                - report.second_moment.powi(2))
            .max(0.0)
                / total_samples as f64)
                .sqrt();
        }
        Ok(report)
    }
}

fn bridge_prime(index: usize) -> u64 {
    const PRIMES: [u64; 16] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];
    PRIMES.get(index).copied().unwrap_or(59 + 2 * index as u64)
}

fn bridge_halton(mut index: usize, base: u64) -> f64 {
    let mut fraction = 1.0;
    let mut value = 0.0;
    while index > 0 {
        fraction /= base as f64;
        value += fraction * (index as u64 % base) as f64;
        index /= base as usize;
    }
    value
}

/// Runtime context needed to turn one bridge point into the sample consumed by
/// the graph evaluator.  The bridge itself is graph-frame aware, while this
/// small context supplies the cache and external-momentum ownership that are
/// deliberately process-local.
#[derive(Clone, Copy, Debug)]
pub struct SamplingMomentumSampleContext<'a> {
    pub loop_mom_cache_id: usize,
    pub external_moms: &'a Externals,
    pub external_mom_cache_id: usize,
    pub dependent_momenta_constructor: DependentMomentaConstructor<'a>,
    pub orientation: Option<usize>,
}

impl<T: FloatLike> SamplingChannelBridgeEvaluation<T> {
    /// Exact positive map determinant times the selected partition weight.
    /// A numerical zero here is never a physical cancellation or selector zero.
    pub fn selected_factor(&self) -> Result<T> {
        let weight = self
            .partition
            .weight(self.channel_id.index())
            .ok_or_else(|| {
                eyre!(
                    "sampling partition has no channel {}",
                    self.channel_id.index()
                )
            })?;
        let jacobian = F(self.map.jacobian.clone());
        let weight = F(weight);
        let factor = &jacobian * &weight;
        if [&jacobian, &weight, &factor]
            .iter()
            .any(|value| !value.0.is_finite() || **value <= value.zero())
        {
            return Err(
                super::sampling_maps::SamplingEvaluationError::Unrepresentable {
                    operation: "selected sampling factor",
                    detail: format!(
                        "channel '{}' has a nonpositive or unrepresentable J*w product {factor}",
                        self.channel_name
                    ),
                }
                .into(),
            );
        }
        Ok(factor.0)
    }

    /// Materialize this exact full-frame bridge point as a graph-evaluator
    /// momentum sample. The production sampler uses the same conversion after
    /// graph-level channel and cut preparation.
    pub fn to_momentum_sample(
        &self,
        context: SamplingMomentumSampleContext<'_>,
    ) -> Result<MomentumSample<T>> {
        if self.raw_coordinates.len() != self.map.point.len() {
            return Err(eyre!(
                "sampling bridge point has {} coordinates but its raw frame has {}",
                self.map.point.len(),
                self.raw_coordinates.len()
            ));
        }
        if self.raw_coordinates.is_empty() || self.raw_coordinates.len() % 3 != 0 {
            return Err(eyre!(
                "sampling bridge raw frame has {} coordinates; expected a non-empty multiple of 3",
                self.raw_coordinates.len()
            ));
        }
        if self.raw_coordinates.iter().any(|value| !value.is_finite()) {
            return Err(eyre!(
                "sampling bridge raw frame contains a non-finite momentum component"
            ));
        }
        if !self.map.jacobian.is_finite() || self.map.jacobian <= self.map.jacobian.zero() {
            return Err(eyre!(
                "sampling bridge map Jacobian must be finite and positive, got {}",
                self.map.jacobian
            ));
        }

        let loop_momenta =
            LoopMomenta::from_iter(self.raw_coordinates.chunks_exact(3).map(|components| {
                ThreeMomentum::new(
                    F(components[0].clone()),
                    F(components[1].clone()),
                    F(components[2].clone()),
                )
            }));
        MomentumSample::new(
            loop_momenta,
            context.loop_mom_cache_id,
            context.external_moms,
            context.external_mom_cache_id,
            F(self.map.jacobian.clone()),
            context.dependent_momenta_constructor,
            context.orientation,
        )
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SamplingChannelBridgeError {
    Empty,
    DimensionMismatch {
        channel: String,
        dimensions: usize,
        expected: usize,
    },
    IncompatibleFrames {
        channel: String,
        edges: Vec<usize>,
        expected: Vec<usize>,
    },
    MasterGraphMismatch {
        channel: String,
        graph: String,
        expected: String,
    },
    UnpreparedContext {
        channel: String,
    },
    MissingFullSupport,
    InexactJacobian {
        channel: String,
        jacobian: super::SamplingJacobian,
    },
    InvalidChannel {
        channel: String,
        error: String,
    },
    UnknownChannel {
        channel: SamplingChannelId,
    },
}

impl fmt::Display for SamplingChannelBridgeError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => {
                formatter.write_str("sampling channel bridge requires at least one channel")
            }
            Self::DimensionMismatch {
                channel,
                dimensions,
                expected,
            } => write!(
                formatter,
                "sampling channel `{channel}` has {dimensions} raw dimensions; the bridge requires {expected}"
            ),
            Self::IncompatibleFrames {
                channel,
                edges,
                expected,
            } => write!(
                formatter,
                "sampling channel `{channel}` uses master-frame edges {edges:?}, incompatible with bridge frame {expected:?}"
            ),
            Self::MasterGraphMismatch {
                channel,
                graph,
                expected,
            } => write!(
                formatter,
                "sampling channel `{channel}` belongs to master graph `{graph}`, but the bridge frame is `{expected}`"
            ),
            Self::UnpreparedContext { channel } => write!(
                formatter,
                "sampling channel `{channel}` requires unprovided context; conditional maps must be completed before graph evaluation"
            ),
            Self::MissingFullSupport => formatter.write_str(
                "restricted sampling channels require an explicitly selected full-support sibling; compact-only selections omit raw volume",
            ),
            Self::InexactJacobian { channel, jacobian } => write!(
                formatter,
                "sampling channel `{channel}` supplies {jacobian:?}; the evaluator bridge requires an exact map Jacobian"
            ),
            Self::InvalidChannel { channel, error } => {
                write!(
                    formatter,
                    "sampling channel `{channel}` is invalid: {error}"
                )
            }
            Self::UnknownChannel { channel } => {
                write!(
                    formatter,
                    "sampling channel index {} is out of range",
                    channel.0
                )
            }
        }
    }
}

impl std::error::Error for SamplingChannelBridgeError {}

impl<T: FloatLike> SamplingChannelBridge<T> {
    /// Construct an exact map-density bridge. Every map must span the complete
    /// raw frame with discharged prerequisites and an exact determinant. The
    /// selected catalogue must include a full-support map so restricted charts
    /// cannot silently omit raw volume.
    pub fn new(
        channels: Vec<CompiledSamplingChannel<T>>,
    ) -> std::result::Result<Self, SamplingChannelBridgeError> {
        Self::new_with_partition_mode(channels, SamplingPartitionMode::MapDensity)
    }

    /// Construct a bridge with an explicit common partition score.  The
    /// singularity-proxy mode requires every compiled channel to carry an
    /// explicit proxy; missing metadata is diagnosed by the partition rather
    /// than silently falling back to map densities.
    pub fn new_with_partition_mode(
        channels: Vec<CompiledSamplingChannel<T>>,
        partition_mode: SamplingPartitionMode,
    ) -> std::result::Result<Self, SamplingChannelBridgeError> {
        let Some(first) = channels.first() else {
            return Err(SamplingChannelBridgeError::Empty);
        };
        let dimensions = first.dimensions();
        let frame = &first.embedded_edges;
        let master_graph = &first.master_graph;
        for channel in &channels {
            if channel.master_graph != *master_graph {
                return Err(SamplingChannelBridgeError::MasterGraphMismatch {
                    channel: channel.name.clone(),
                    graph: channel.master_graph.clone(),
                    expected: master_graph.clone(),
                });
            }
            if channel.embedded_edges != *frame {
                return Err(SamplingChannelBridgeError::IncompatibleFrames {
                    channel: channel.name.clone(),
                    edges: channel.embedded_edges.clone(),
                    expected: frame.clone(),
                });
            }
            if channel.dimensions() != dimensions {
                return Err(SamplingChannelBridgeError::DimensionMismatch {
                    channel: channel.name.clone(),
                    dimensions: channel.dimensions(),
                    expected: dimensions,
                });
            }
            let contract = channel.contract();
            if contract.requires_context {
                return Err(SamplingChannelBridgeError::UnpreparedContext {
                    channel: channel.name.clone(),
                });
            }
            if matches!(contract.jacobian, super::SamplingJacobian::ProxyOnly) {
                return Err(SamplingChannelBridgeError::InexactJacobian {
                    channel: channel.name.clone(),
                    jacobian: contract.jacobian,
                });
            }
        }
        if !channels
            .iter()
            .any(|channel| channel.contract().support == SamplingSupport::Full)
        {
            return Err(SamplingChannelBridgeError::MissingFullSupport);
        }
        Ok(Self {
            channels,
            dimensions,
            partition_mode,
            relative_density_tolerance: None,
        })
    }

    /// Reserve a fraction of the configured precision target for consistency
    /// between the selected forward determinant and its actual inverse density.
    pub(crate) fn with_relative_density_tolerance(mut self, tolerance: f64) -> Result<Self> {
        if !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(eyre!(
                "sampling inverse density tolerance must be positive and finite"
            ));
        }
        self.relative_density_tolerance = Some(tolerance);
        Ok(self)
    }

    pub(crate) fn relative_density_tolerance(&self, precision: &F<T>) -> F<T> {
        self.relative_density_tolerance
            .map(F::<T>::from_f64)
            .unwrap_or_else(|| precision.epsilon().sqrt())
    }

    pub fn channels(&self) -> &[CompiledSamplingChannel<T>] {
        &self.channels
    }

    pub fn dimensions(&self) -> usize {
        self.dimensions
    }

    pub fn partition(&self, raw_coordinates: &[T]) -> Result<SamplingPartition<T>> {
        self.partition_with_contexts(
            raw_coordinates,
            &mut SamplingChannelRuntimeContexts::new(self.channels.len()),
        )
    }

    /// Build the common inverse-density partition using the context belonging
    /// to every canonical channel.  A selected channel's context alone is
    /// insufficient: all denominator scores must describe the same raw point
    /// and their own map branches.
    pub fn partition_with_contexts(
        &self,
        raw_coordinates: &[T],
        contexts: &mut SamplingChannelRuntimeContexts<'_, T>,
    ) -> Result<SamplingPartition<T>> {
        if contexts.selected_host_count.is_none() {
            return Err(eyre!(
                "sampling partition cannot run before selected host authority is sealed"
            ));
        }
        if raw_coordinates.len() != self.dimensions {
            return Err(color_eyre::eyre::eyre!(
                "raw sampling frame has {}, expected {} dimensions",
                raw_coordinates.len(),
                self.dimensions
            ));
        }
        if contexts.channel_count() != self.channels.len() {
            return Err(eyre!(
                "sampling channel runtime context has {} channels, expected {}",
                contexts.channel_count(),
                self.channels.len()
            ));
        }
        SamplingPartition::from_log_scores(
            self.partition_mode,
            raw_coordinates,
            self.channels.iter().map(|channel| channel.name.as_str()),
            |index| {
                let channel = &self.channels[index];
                let mut context = contexts.for_channel(SamplingChannelId(index))?;
                match self.partition_mode {
                    SamplingPartitionMode::MapDensity => {
                        let Some(evaluation) =
                            channel.inverse_with_context(raw_coordinates, &mut context)?
                        else {
                            return Ok(None);
                        };
                        let density = F(evaluation.inverse_jacobian);
                        if !density.0.is_finite() || density <= density.zero() {
                            return Err(SamplingEvaluationError::Unrepresentable {
                                operation: "inverse map density",
                                detail: format!(
                                    "supported density is not finite and positive: {density}"
                                ),
                            }
                            .into());
                        }
                        Ok(Some(density.ln().0))
                    }
                    SamplingPartitionMode::SingularityProxy => {
                        let proxy = channel.singularity_proxy.as_ref().ok_or_else(|| {
                            eyre!(
                                "channel '{}' has no singularity_proxy metadata; provide a positive expression for every channel when sampling_channel_weight = 'singularity_proxy'",
                                channel.name
                            )
                        })?;
                        // Proxy scores cannot give a compact chart weight at a
                        // point with no preimage. Full maps keep the cheap proxy
                        // path; each restricted inverse uses its own context.
                        if channel.contract().support == SamplingSupport::Restricted
                            && channel
                                .inverse_with_context(raw_coordinates, &mut context)?
                                .is_none()
                        {
                            return Ok(None);
                        }
                        proxy.evaluate(raw_coordinates)
                    }
                }
            },
        )
    }

    pub fn forward(
        &self,
        channel_id: SamplingChannelId,
        coordinates: &[T],
    ) -> Result<SamplingChannelBridgeEvaluation<T>> {
        self.forward_with_context(channel_id, coordinates, &[])
    }

    /// Forward a channel while supplying the output of an earlier conditional
    /// map (for example a sampled complement block) to the selected map.
    pub fn forward_with_context(
        &self,
        channel_id: SamplingChannelId,
        coordinates: &[T],
        context: &[T],
    ) -> Result<SamplingChannelBridgeEvaluation<T>> {
        let mut contexts = SamplingChannelRuntimeContexts::new(self.channels.len());
        contexts.set(channel_id, context.to_vec())?;
        self.forward_with_runtime_contexts(channel_id, coordinates, &mut contexts)
    }

    pub fn forward_with_runtime_contexts(
        &self,
        channel_id: SamplingChannelId,
        coordinates: &[T],
        contexts: &mut SamplingChannelRuntimeContexts<'_, T>,
    ) -> Result<SamplingChannelBridgeEvaluation<T>> {
        let channel_index = channel_id.0;
        let channel =
            self.channels
                .get(channel_index)
                .ok_or(SamplingChannelBridgeError::UnknownChannel {
                    channel: channel_id,
                })?;
        contexts.begin_selected(channel_id)?;
        let map =
            channel.forward_with_context(coordinates, &mut contexts.for_channel(channel_id)?)?;
        contexts.seal_selected();
        if map.point.len() != self.dimensions {
            return Err(SamplingChannelBridgeError::DimensionMismatch {
                channel: channel.name.clone(),
                dimensions: map.point.len(),
                expected: self.dimensions,
            }
            .into());
        }
        let partition = self.partition_with_contexts(&map.point, contexts)?;
        // A small Cartesian round-trip residual does not certify the density
        // near a focused shell. Compare the selected determinant to the inverse
        // density at the actual generated point, reusing the partition score
        // when it already evaluates that inverse. Proxy partitions additionally
        // invert restricted foreign maps to establish their actual support.
        let log_density = match self.partition_mode {
            SamplingPartitionMode::MapDensity => F(partition.log_scores[channel_index]
                .clone()
                .ok_or_else(|| SamplingEvaluationError::UncertainGeometry {
                    detail: format!(
                        "forward point is outside selected channel '{}' inverse support",
                        channel.name
                    ),
                })?),
            SamplingPartitionMode::SingularityProxy => {
                let inverse = channel
                    .inverse_with_context(&map.point, &mut contexts.for_channel(channel_id)?)?
                    .ok_or_else(|| SamplingEvaluationError::UncertainGeometry {
                        detail: format!(
                            "forward point is outside selected channel '{}' inverse support",
                            channel.name
                        ),
                    })?;
                F(inverse.inverse_jacobian).ln()
            }
        };
        let jacobian = F(map.jacobian.clone());
        let one = jacobian.one();
        let tolerance = self.relative_density_tolerance(&one);
        let log_ratio = jacobian.ln() + log_density;
        let two = &one + &one;
        // log(1 +/- tolerance) must retain a sub-epsilon budget. The atanh
        // identity avoids forming 1 +/- tolerance in that small-budget regime.
        let upper = if tolerance < &one / &two {
            &two * F((&tolerance / (&two + &tolerance)).0.atanh())
        } else {
            (&one + &tolerance).ln()
        };
        let lower = if tolerance < &one / &two {
            Some(-&two * F((&tolerance / (&two - &tolerance)).0.atanh()))
        } else if tolerance < one {
            Some((&one - &tolerance).ln())
        } else {
            None
        };
        if !jacobian.0.is_finite()
            || jacobian <= jacobian.zero()
            || !log_ratio.0.is_finite()
            || log_ratio > upper
            || lower.is_some_and(|lower| log_ratio < lower)
        {
            return Err(super::sampling_maps::SamplingEvaluationError::Unrepresentable {
                operation: "sampling inverse density consistency",
                detail: format!(
                    "channel '{}' has log(J_forward * q_inverse)={log_ratio}, outside relative tolerance {tolerance}",
                    channel.name
                ),
            }.into());
        }
        Ok(SamplingChannelBridgeEvaluation {
            channel_id,
            channel_name: channel.name.clone(),
            raw_coordinates: map.point.clone(),
            map,
            partition,
            prepared_lu_hosts: contexts.selected_hosts(),
        })
    }

    pub fn inverse(
        &self,
        channel_id: SamplingChannelId,
        raw_coordinates: &[T],
    ) -> Result<Option<SamplingChannelBridgeEvaluation<T>>> {
        self.inverse_with_context(channel_id, raw_coordinates, &[])
    }

    pub fn inverse_with_context(
        &self,
        channel_id: SamplingChannelId,
        raw_coordinates: &[T],
        context: &[T],
    ) -> Result<Option<SamplingChannelBridgeEvaluation<T>>> {
        let mut contexts = SamplingChannelRuntimeContexts::new(self.channels.len());
        contexts.set(channel_id, context.to_vec())?;
        self.inverse_with_runtime_contexts(channel_id, raw_coordinates, &mut contexts)
    }

    pub fn inverse_with_runtime_contexts(
        &self,
        channel_id: SamplingChannelId,
        raw_coordinates: &[T],
        contexts: &mut SamplingChannelRuntimeContexts<'_, T>,
    ) -> Result<Option<SamplingChannelBridgeEvaluation<T>>> {
        let channel_index = channel_id.0;
        let channel =
            self.channels
                .get(channel_index)
                .ok_or(SamplingChannelBridgeError::UnknownChannel {
                    channel: channel_id,
                })?;
        contexts.begin_selected(channel_id)?;
        let Some(map) = channel
            .inverse_with_context(raw_coordinates, &mut contexts.for_channel(channel_id)?)?
        else {
            return Ok(None);
        };
        contexts.seal_selected();
        let partition = self.partition_with_contexts(raw_coordinates, contexts)?;
        Ok(Some(SamplingChannelBridgeEvaluation {
            channel_id,
            channel_name: channel.name.clone(),
            raw_coordinates: raw_coordinates.to_vec(),
            map,
            partition,
            prepared_lu_hosts: contexts.selected_hosts(),
        }))
    }
}

impl<T: FloatLike> CompiledSamplingChannel<T> {
    pub fn contract(&self) -> SamplingMapContract {
        self.map.contract()
    }

    pub fn dimensions(&self) -> usize {
        self.map.dimensions()
    }

    pub fn forward(&self, coordinates: &[T]) -> Result<SamplingMapEvaluation<T>> {
        self.forward_with_context(coordinates, &mut SamplingMapContext::detached(&[]))
    }

    pub fn forward_with_context(
        &self,
        coordinates: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<SamplingMapEvaluation<T>> {
        self.map
            .forward_with_context(coordinates, context)
            .wrap_err_with(|| {
                format!(
                    "sampling graph '{}' channel '{}' target {:?}, master-frame edges {:?}",
                    self.master_graph, self.name, self.definition, self.embedded_edges,
                )
            })
    }

    pub fn inverse(&self, point: &[T]) -> Result<Option<SamplingMapEvaluation<T>>> {
        self.inverse_with_context(point, &mut SamplingMapContext::detached(&[]))
    }

    pub fn inverse_with_context(
        &self,
        point: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<SamplingMapEvaluation<T>>> {
        self.map
            .inverse_with_context(point, context)
            .wrap_err_with(|| {
                format!(
                    "sampling graph '{}' channel '{}' target {:?}, master-frame edges {:?}",
                    self.master_graph, self.name, self.definition, self.embedded_edges,
                )
            })
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SamplingChannelCompileError {
    EmptyMasterGraph,
    UnsupportedMap {
        channel: String,
        map: String,
    },
    /// A soft or collinear selector was parsed successfully but cannot yet be
    /// turned into a map-density channel from edge ids alone.  These targets
    /// need a graph-resolved frame and an explicit radial/angular profile;
    /// retaining this distinction avoids accidentally treating a proxy score
    /// as the selected channel's integration Jacobian.
    UnsupportedSingularPrimitive {
        channel: String,
        primitive: String,
        reason: String,
    },
    MissingGeometry {
        channel: String,
        target: SamplingMapDefinition,
        subspace_lmb: Vec<usize>,
    },
    MissingLmbFrameMap {
        channel: String,
        basis_id: Option<usize>,
        edges: Vec<usize>,
        parent_lmb: Vec<usize>,
    },
    InvalidLmbFrameMap {
        channel: String,
        basis_id: usize,
        error: String,
    },
    InvalidChannel {
        channel: String,
        error: String,
    },
}

impl fmt::Display for SamplingChannelCompileError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyMasterGraph => {
                formatter.write_str("sampling channel compilation requires a master graph")
            }
            Self::UnsupportedMap { channel, map } => write!(
                formatter,
                "sampling channel `{channel}` uses map `{map}`, which needs prepared graph context before it can be compiled"
            ),
            Self::UnsupportedSingularPrimitive {
                channel,
                primitive,
                reason,
            } => write!(
                formatter,
                "sampling channel `{channel}` uses singular primitive `{primitive}`, which is not compiled: {reason}"
            ),
            Self::MissingGeometry {
                channel,
                target,
                subspace_lmb,
            } => write!(
                formatter,
                "sampling channel `{channel}` has no prepared geometry for target {target:?} in subspace_lmb {subspace_lmb:?}"
            ),
            Self::MissingLmbFrameMap {
                channel,
                basis_id,
                edges,
                parent_lmb,
            } => write!(
                formatter,
                "sampling channel `{channel}` uses non-master LMB edges {edges:?} (parent LMB {parent_lmb:?}) but no affine frame map was supplied for basis {basis_id:?}"
            ),
            Self::InvalidLmbFrameMap {
                channel,
                basis_id,
                error,
            } => write!(
                formatter,
                "sampling channel `{channel}` has invalid affine frame map for basis {basis_id}: {error}"
            ),
            Self::InvalidChannel { channel, error } => write!(
                formatter,
                "sampling channel `{channel}` could not be compiled: {error}"
            ),
        }
    }
}

impl std::error::Error for SamplingChannelCompileError {}

/// Keep parsed soft/collinear syntax from becoming a silently approximate
/// channel.  A useful singular proposal must be resolved against the actual
/// routed vectors (and, for collinearity, a relative angular frame) and carry
/// a normalized profile.  Neither can be reconstructed from the edge labels
/// in a graph-independent compiler.  Until those ingredients are supplied,
/// report a capability error rather than installing a proxy-only map in the
/// exact-density catalogue.
fn unsupported_map_error(
    channel: &str,
    map_label: String,
    map: &SamplingMapDefinition,
) -> SamplingChannelCompileError {
    match singular_primitive(map) {
        Some(SamplingMapDefinition::Soft(edge)) => {
            SamplingChannelCompileError::UnsupportedSingularPrimitive {
                channel: channel.to_owned(),
                primitive: format!("soft({edge})"),
                reason: "a graph-resolved routed three-momentum frame and a normalized radial profile are required; edge ids alone do not define an exact push-forward density".to_owned(),
            }
        }
        Some(SamplingMapDefinition::Collinear(a, b)) => {
            SamplingChannelCompileError::UnsupportedSingularPrimitive {
                channel: channel.to_owned(),
                primitive: format!("collinear({a},{b})"),
                reason: "a graph-resolved relative angular frame, branch/support rules and a normalized angular profile are required; edge ids alone do not define an exact push-forward density".to_owned(),
            }
        }
        Some(other) => SamplingChannelCompileError::UnsupportedMap {
            channel: channel.to_owned(),
            map: format!("{}: unsupported primitive {other:?}", map_label),
        },
        None => SamplingChannelCompileError::UnsupportedMap {
            channel: channel.to_owned(),
            map: map_label,
        },
    }
}

fn singular_primitive(map: &SamplingMapDefinition) -> Option<&SamplingMapDefinition> {
    match map {
        SamplingMapDefinition::Soft(_) | SamplingMapDefinition::Collinear(_, _) => Some(map),
        SamplingMapDefinition::Product(maps)
        | SamplingMapDefinition::Intersect(maps)
        | SamplingMapDefinition::Then(maps) => maps.iter().find_map(singular_primitive),
        SamplingMapDefinition::PhaseSpace(map)
        | SamplingMapDefinition::Left(map)
        | SamplingMapDefinition::Right(map)
        | SamplingMapDefinition::AtCut { map, .. }
        | SamplingMapDefinition::Block { map, .. } => singular_primitive(map),
        _ => None,
    }
}

fn compile_lmb_map<T: FloatLike>(
    channel: &str,
    basis_id: Option<usize>,
    edges: &[usize],
    context: &SamplingChannelCompileContext<T>,
) -> Result<CompiledSamplingMap<T>, SamplingChannelCompileError> {
    let definition = SamplingMapDefinition::Lmb(edges.to_vec());
    let lmb = SamplingMapKernel::new(
        definition,
        context.parameterization_settings.clone(),
        context.e_cm,
        context.n_loop_momenta,
    )
    .map_err(|error| SamplingChannelCompileError::InvalidChannel {
        channel: channel.to_owned(),
        error: error.to_string(),
    })?;
    if edges == context.parent_lmb {
        return Ok(CompiledSamplingMap::Lmb(lmb));
    }
    let frame = match basis_id {
        Some(id) => context.lmb_frame_maps.get(&id),
        None => context.lmb_frame_maps_by_edges.get(edges),
    };
    let Some(frame) = frame else {
        return Err(SamplingChannelCompileError::MissingLmbFrameMap {
            channel: channel.to_owned(),
            basis_id,
            edges: edges.to_vec(),
            parent_lmb: context.parent_lmb.clone(),
        });
    };
    let expected_dimension = 3 * context.n_loop_momenta;
    if frame.dimension() != expected_dimension {
        return Err(SamplingChannelCompileError::InvalidLmbFrameMap {
            channel: channel.to_owned(),
            basis_id: basis_id.unwrap_or(usize::MAX),
            error: format!(
                "affine frame dimension {} does not match parent raw dimension {expected_dimension}",
                frame.dimension()
            ),
        });
    }
    Ok(CompiledSamplingMap::Affine {
        map: Box::new(CompiledSamplingMap::Lmb(lmb)),
        frame: frame.clone(),
    })
}

impl SamplingChannelCatalogue {
    /// Return generated-LMB metadata carried by canonical catalogue entries.
    ///
    /// The returned `basis_id` is an implementation detail used to prepare
    /// affine frame maps; it is deliberately not a channel index.  Callers
    /// that need the production channel axis must enumerate `entries` and
    /// assign [`SamplingChannelId`] from the entry position.
    pub fn lmb_basis_entries(&self) -> impl Iterator<Item = (usize, &[usize])> {
        self.entries.iter().filter_map(|entry| match entry {
            SamplingCatalogueEntry::Lmb {
                basis_id, edges, ..
            } => Some((*basis_id, edges.as_slice())),
            SamplingCatalogueEntry::Named(_) => None,
            SamplingCatalogueEntry::Surface { .. } => None,
        })
    }

    pub fn named_entries(&self) -> impl Iterator<Item = &ResolvedNamedSamplingChannel> {
        self.entries.iter().filter_map(|entry| match entry {
            SamplingCatalogueEntry::Named(channel) => Some(channel),
            SamplingCatalogueEntry::Lmb { .. } => None,
            SamplingCatalogueEntry::Surface { .. } => None,
        })
    }

    /// Audit ordinary full-domain and elementary soft coverage using the same
    /// catalogue entries that drive production channel enumeration.
    pub fn coverage_report(
        &self,
        all_lmbs: &[(usize, Vec<usize>)],
        massless_edges: &[usize],
    ) -> SamplingCoverageReport {
        let selected_lmbs = self
            .entries
            .iter()
            .filter_map(|entry| match entry {
                SamplingCatalogueEntry::Lmb {
                    basis_id, edges, ..
                } => Some((*basis_id, edges.as_slice())),
                _ => None,
            })
            .collect::<Vec<_>>();
        // A zero-loop graph has a zero-dimensional momentum domain and needs
        // no ordinary LMB coordinate block to cover it.
        let has_full_domain_lmb = !selected_lmbs.is_empty() || all_lmbs.is_empty();
        let relevant_massless = massless_edges
            .iter()
            .copied()
            .filter(|edge| all_lmbs.iter().any(|(_, lmb)| lmb.contains(edge)))
            .collect::<BTreeSet<_>>();
        let covered_massless_loop_edges = relevant_massless
            .iter()
            .copied()
            .filter(|edge| selected_lmbs.iter().any(|(_, lmb)| lmb.contains(edge)))
            .collect::<Vec<_>>();
        let covered_set = covered_massless_loop_edges
            .iter()
            .copied()
            .collect::<BTreeSet<_>>();
        let missing_massless_loop_edges = relevant_massless
            .difference(&covered_set)
            .copied()
            .collect();
        SamplingCoverageReport {
            has_full_domain_lmb,
            covered_massless_loop_edges,
            missing_massless_loop_edges,
        }
    }

    /// Compile user supplied proxies, radial profiles and the neutral joint map
    /// once in canonical channel order. Native geometry bindings borrow these
    /// same programs; distinct joint channels clone the compiled expression.
    /// Proxy expressions use the complete master raw frame's x0, x1, ...;
    /// LU profiles inherit the actual physical runtime h-function settings.
    pub fn compile_programs(
        &self,
        dimensions: usize,
        lu_h: &HFunctionSettings,
    ) -> Result<Vec<SamplingChannelPrograms>> {
        let joint_program = self
            .named_entries()
            .any(|channel| {
                channel
                    .blocks
                    .iter()
                    .any(|block| block.target.energy_edge_sets().len() == 2)
            })
            .then(SharedEnergyJointMap::<f64>::compile_program)
            .transpose()?;
        self.entries.iter().map(|entry| {
            let SamplingCatalogueEntry::Named(channel) = entry else { return Ok((None,None,None)); };
            let proxy = channel.singularity_proxy.as_ref().map(|expression| {
                if dimensions == 0 {
                    return Err(eyre!("sampling singularity proxy requires a non-empty raw coordinate frame"));
                }
                let parameters = (0..dimensions)
                    .map(|index| Atom::var(symbol!(format!("x{index}"))))
                    .collect::<Vec<_>>();
                SamplingExpressionEvaluator::new_positive_proxy(expression.clone(), parameters)
            }).transpose()?;
            let radial = channel.definition.radial_profile.as_ref().map(|profile| {
                if channel.blocks.iter().filter(|block| block.target.is_phase_space()).count() != 1 {
                    return Err(eyre!("channel {}: radial_profile=lu_h requires exactly one phase_space(cut(...)) block", channel.name));
                }
                SamplingExpressionEvaluator::new_lu_h_profile(profile, lu_h)
            }).transpose()?;
            let joint = if channel.blocks.iter().any(|block| block.target.energy_edge_sets().len() == 2) {
                joint_program.clone()
            } else { None };
            Ok((proxy,radial,joint))
        }).collect()
    }

    /// Compile the resolved block catalogue against native graph geometry.
    /// Qualified target, parent, active and preceding edges identify each
    /// binding. Cut and side maps prepare their auxiliary LU scale and physical
    /// prerequisites at the actual point. A joint target is one registered 3D
    /// block; bare `cut` and missing kinematics remain rejected, avoiding
    /// plausible maps in an incorrect frame. Production joint binding remains
    /// gated by its graph matcher and retained proposal policy.
    /// Structural diagnostics remain typed; native binding failures propagate
    /// unchanged so warmup and stability can retry another precision.
    pub fn compile<T: FloatLike>(
        &self,
        context: &SamplingChannelCompileContext<T>,
        programs: &[SamplingChannelPrograms],
    ) -> Result<Vec<CompiledSamplingChannel<T>>> {
        if programs.len() != self.entries.len() {
            return Err(SamplingChannelCompileError::InvalidChannel {
                channel: self.graph_name.clone(),
                error: format!(
                    "compiled sampling program vector has {} entries, expected {} canonical channels",
                    programs.len(),
                    self.entries.len()
                ),
            }.into());
        }
        if context.master_graph.trim().is_empty() {
            return Err(SamplingChannelCompileError::EmptyMasterGraph.into());
        }
        let mut parent_edges = BTreeSet::new();
        if context.parent_lmb.len() != context.n_loop_momenta
            || context
                .parent_lmb
                .iter()
                .any(|edge| !parent_edges.insert(*edge))
        {
            return Err(SamplingChannelCompileError::InvalidChannel {
                channel: context.master_graph.clone(),
                error: format!(
                    "parent LMB {:?} must contain one unique edge per loop momentum (expected {})",
                    context.parent_lmb, context.n_loop_momenta
                ),
            }
            .into());
        }
        let mut compiled = Vec::with_capacity(self.entries.len());
        for (entry, (proxy, radial, joint)) in self.entries.iter().zip(programs) {
            let (name, basis_id, definition, map, singularity_proxy) = match entry {
                SamplingCatalogueEntry::Lmb {
                    basis_id, edges, ..
                } => {
                    let definition = SamplingMapDefinition::Lmb(edges.clone());
                    let map = compile_lmb_map(
                        &format!("lmb[{basis_id}]"),
                        Some(*basis_id),
                        edges,
                        context,
                    )?;
                    (
                        format!("lmb[{basis_id}]"),
                        Some(*basis_id),
                        definition,
                        map,
                        None,
                    )
                }
                SamplingCatalogueEntry::Surface { edges, parent_lmb } => {
                    if parent_lmb != &context.parent_lmb {
                        return Err(SamplingChannelCompileError::InvalidChannel {
                            channel: format!("surface:{edges:?}"),
                            error: format!(
                                "parent LMB {:?} does not match master context {:?}",
                                parent_lmb, context.parent_lmb
                            ),
                        }
                        .into());
                    }
                    let definition = SamplingMapDefinition::Surface(edges.clone());
                    let key = (
                        definition.clone(),
                        parent_lmb.clone(),
                        parent_lmb.clone(),
                        Vec::new(),
                    );
                    let map = context.geometry_maps.get(&key).cloned().ok_or_else(|| {
                        SamplingChannelCompileError::MissingGeometry {
                            channel: format!("surface:{edges:?}"),
                            target: definition.clone(),
                            subspace_lmb: parent_lmb.clone(),
                        }
                    })?;
                    (format!("surface:{edges:?}"), None, definition, map, None)
                }
                SamplingCatalogueEntry::Named(channel) => {
                    let definition = channel.map.clone();
                    if joint.is_some()
                        != channel
                            .blocks
                            .iter()
                            .any(|block| block.target.energy_edge_sets().len() == 2)
                    {
                        return Err(SamplingChannelCompileError::InvalidChannel {
                            channel: channel.name.clone(),
                            error: "compiled joint program presence does not match the canonical definition".into(),
                        }.into());
                    }
                    let map = channel.compile_map(context, radial.as_ref())?;
                    if proxy.is_some() != channel.singularity_proxy.is_some() {
                        return Err(SamplingChannelCompileError::InvalidChannel {
                            channel: channel.name.clone(),
                            error:
                                "compiled proxy presence does not match the canonical definition"
                                    .into(),
                        }
                        .into());
                    }
                    let singularity_proxy = proxy
                        .as_ref()
                        .map(|program| -> Result<_> {
                            if program.parameter_count() != 3 * context.n_loop_momenta {
                                return Err(SamplingChannelCompileError::InvalidChannel {
                                    channel: channel.name.clone(),
                                    error:
                                        "compiled proxy does not use the complete master raw frame"
                                            .into(),
                                }
                                .into());
                            }
                            SamplingScoreFunction::from_compiled_symbolica_positive_expression(
                                program.clone(),
                            )
                            .map_err(|error| {
                                SamplingChannelCompileError::InvalidChannel {
                                    channel: channel.name.clone(),
                                    error: format!("invalid compiled singularity_proxy: {error}"),
                                }
                                .into()
                            })
                        })
                        .transpose()?;
                    (
                        channel.name.clone(),
                        None,
                        definition,
                        map,
                        singularity_proxy,
                    )
                }
            };
            // Physical energy labels identify constraints, not output axes.
            // Every compiled complete channel now returns the same master frame,
            // including native-parent affine routing and explicit complements.
            let embedded_edges = context.parent_lmb.clone();
            compiled.push(CompiledSamplingChannel {
                name,
                master_graph: context.master_graph.clone(),
                basis_id,
                definition,
                embedded_edges,
                map,
                singularity_proxy,
            });
        }
        Ok(compiled)
    }

    /// Human-readable inspection rows, stable across runs and suitable for
    /// CLI/API diagnostics.
    pub fn inspection_rows(&self) -> Vec<String> {
        self.entries
            .iter()
            .enumerate()
            .map(|(index, entry)| match entry {
                SamplingCatalogueEntry::Lmb {
                    basis_id,
                    edges,
                    preset,
                } => format!(
                    "{index}: lmb basis={basis_id} edges={edges:?} source={}",
                    preset.as_str()
                ),
                SamplingCatalogueEntry::Surface { edges, parent_lmb } => {
                    format!("{index}: surface edges={edges:?} parent_lmb={parent_lmb:?}")
                }
                SamplingCatalogueEntry::Named(channel) => format!(
                    "{index}: {} around={} subspace_lmb={:?} parent_lmb={:?} on_cut={:?}",
                    channel.name,
                    channel.definition.around,
                    channel.definition.subspace_lmb,
                    channel.definition.parent_lmb,
                    channel.definition.on_cut
                ),
            })
            .collect()
    }

    /// Return the canonical resolved selectors and catalogue rows for
    /// read-only CLI/API inspection.
    pub fn inspection(&self) -> SamplingChannelInspection {
        SamplingChannelInspection {
            graph_name: self.graph_name.clone(),
            selectors: self.selectors.iter().map(ToString::to_string).collect(),
            entries: self.inspection_rows(),
        }
    }
}

/// Expand a resolved graph selection against the generated LMB catalogue.
/// Entries are deduplicated by their resolved identity while preserving the
/// first selector's order.  A surface preset includes the optimized LMB
/// fallback so that a surface-aware selection remains defined where a surface
/// is absent; numerical use of these entries is implemented by the process
/// sampler in a later phase.
pub fn build_sampling_channel_catalogue(
    resolved: &ResolvedSamplingChannelSelection,
    all_lmbs: &[(usize, Vec<usize>)],
    optimized_lmbs: &[usize],
) -> SamplingChannelCatalogue {
    build_sampling_channel_catalogue_with_surfaces(resolved, all_lmbs, optimized_lmbs, &[], &[])
}

/// Expand a selection and append automatically enumerated surface candidates.
/// `surface_edges` and `parent_lmb` are already resolved in the master graph
/// frame; no graph mapping is inferred here.
pub fn build_sampling_channel_catalogue_with_surfaces(
    resolved: &ResolvedSamplingChannelSelection,
    all_lmbs: &[(usize, Vec<usize>)],
    optimized_lmbs: &[usize],
    surface_edges: &[Vec<usize>],
    parent_lmb: &[usize],
) -> SamplingChannelCatalogue {
    build_sampling_channel_catalogue_with_surfaces_and_coverage(
        resolved,
        all_lmbs,
        optimized_lmbs,
        surface_edges,
        parent_lmb,
        &[],
    )
}

/// Expand a selection while applying the automatic elementary soft-coverage
/// audit.  Physical E-surface candidates remain an explicit input prepared by
/// the process layer; this helper never invents them from edge ids.
pub fn build_sampling_channel_catalogue_with_surfaces_and_coverage(
    resolved: &ResolvedSamplingChannelSelection,
    all_lmbs: &[(usize, Vec<usize>)],
    optimized_lmbs: &[usize],
    surface_edges: &[Vec<usize>],
    parent_lmb: &[usize],
    massless_edges: &[usize],
) -> SamplingChannelCatalogue {
    let mut entries = Vec::new();
    for selector in &resolved.selectors {
        let Some(preset) = selector.preset() else {
            if let SamplingChannelSelector::Named(name) = selector {
                if let Some(channel) = resolved.named(name) {
                    let entry = SamplingCatalogueEntry::Named(channel.clone());
                    if !entries.contains(&entry) {
                        entries.push(entry);
                    }
                }
            }
            continue;
        };
        let mut basis_ids: Vec<usize> = match preset {
            SamplingChannelPreset::Lmb => all_lmbs.iter().map(|(id, _)| *id).collect(),
            SamplingChannelPreset::OptimizedLmb | SamplingChannelPreset::Surfaces => optimized_lmbs
                .iter()
                .copied()
                .filter(|id| all_lmbs.iter().any(|(candidate, _)| candidate == id))
                .collect(),
        };
        if matches!(
            preset,
            SamplingChannelPreset::OptimizedLmb | SamplingChannelPreset::Surfaces
        ) {
            // A surface-aware or optimized preset must retain at least one
            // ordinary full-domain channel when generated LMBs exist.  This
            // is a deterministic fallback, not a hidden second catalogue.
            if basis_ids.is_empty() && optimized_lmbs.is_empty() {
                if let Some((basis_id, _)) = all_lmbs.first() {
                    basis_ids.push(*basis_id);
                }
            }
            // Add the smallest available ordinary channel exposing each
            // elementary massless loop edge omitted by the heuristic.  This
            // is intentionally weaker than a compound-soft guarantee.
            if !basis_ids.is_empty() {
                for edge in massless_edges {
                    if basis_ids.iter().any(|basis_id| {
                        all_lmbs
                            .iter()
                            .find(|(candidate, _)| candidate == basis_id)
                            .is_some_and(|(_, lmb)| lmb.contains(edge))
                    }) {
                        continue;
                    }
                    if let Some((basis_id, _)) = all_lmbs.iter().find(|(_, lmb)| lmb.contains(edge))
                    {
                        basis_ids.push(*basis_id);
                    }
                }
            }
        }
        for basis_id in basis_ids {
            let Some((_, edges)) = all_lmbs.iter().find(|(id, _)| *id == basis_id) else {
                continue;
            };
            let entry = SamplingCatalogueEntry::Lmb {
                basis_id,
                edges: edges.clone(),
                preset,
            };
            let duplicate = entries.iter().any(|existing| {
                matches!(
                    existing,
                    SamplingCatalogueEntry::Lmb {
                        basis_id: existing_basis,
                        edges: existing_edges,
                        ..
                    } if *existing_basis == basis_id && existing_edges == edges
                )
            });
            if !duplicate {
                entries.push(entry);
            }
        }
        if preset == SamplingChannelPreset::Surfaces {
            for edges in surface_edges {
                if !edges.is_empty()
                    && !entries.iter().any(|entry| {
                        matches!(entry, SamplingCatalogueEntry::Surface {
                            edges: existing, parent_lmb: parent
                        } if existing == edges && parent == parent_lmb)
                    })
                {
                    entries.push(SamplingCatalogueEntry::Surface {
                        edges: edges.clone(),
                        parent_lmb: parent_lmb.to_vec(),
                    });
                }
            }
        }
    }
    SamplingChannelCatalogue {
        graph_name: resolved.graph_name.clone(),
        selectors: resolved.selectors.clone(),
        entries,
    }
}

impl ResolvedSamplingChannelSelection {
    pub fn has_preset(&self, preset: SamplingChannelPreset) -> bool {
        self.selectors
            .iter()
            .any(|selector| selector.preset() == Some(preset))
    }

    pub fn named(&self, name: &str) -> Option<&ResolvedNamedSamplingChannel> {
        self.named_channels
            .iter()
            .find(|channel| channel.name == name)
    }
}

/// Errors from selector parsing or graph-scoped resolution.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SamplingSelectionError {
    EmptyGraphName,
    EmptySelector,
    UnknownPreset(String),
    MissingChannelDefinition {
        graph: String,
        channel: String,
        available: Vec<String>,
    },
    MissingParentLmb {
        graph: String,
        channel: String,
    },
    MissingSubspaceLmb {
        graph: String,
        channel: String,
    },
    InvalidChannelDefinition {
        graph: String,
        channel: String,
        error: String,
    },
}

impl fmt::Display for SamplingSelectionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyGraphName => {
                formatter.write_str("sampling channel selection requires a graph name")
            }
            Self::EmptySelector => {
                formatter.write_str("sampling channel selection contains an empty selector")
            }
            Self::UnknownPreset(preset) => write!(
                formatter,
                "unknown sampling channel preset `{preset}`; expected auto:lmb, auto:optimized_lmb or auto:surfaces"
            ),
            Self::MissingChannelDefinition {
                graph,
                channel,
                available,
            } => write!(
                formatter,
                "sampling channel `{graph}.{channel}` was selected but has no channel definition; available definitions: {available:?}"
            ),
            Self::MissingParentLmb { graph, channel } => write!(
                formatter,
                "sampling channel `{graph}.{channel}` must declare a non-empty parent_lmb"
            ),
            Self::MissingSubspaceLmb { graph, channel } => write!(
                formatter,
                "sampling channel `{graph}.{channel}` contains a surface/cut map and must declare a non-empty subspace_lmb"
            ),
            Self::InvalidChannelDefinition {
                graph,
                channel,
                error,
            } => write!(
                formatter,
                "sampling channel `{graph}.{channel}` has an invalid around expression: {error}"
            ),
        }
    }
}

impl std::error::Error for SamplingSelectionError {}

/// Resolve the selectors that apply to `graph_name`.
///
/// A graph-specific entry replaces the default list. If no entry exists for
/// this graph, the default list applies. Entries within the selected list are
/// additive and deduplicated in declaration order. An omitted default with no
/// graph override selects `auto:optimized_lmb`; an explicit empty graph list
/// stays empty. Counting and compiling channels share this default resolution.
pub fn resolve_sampling_channel_selection(
    graph_name: &str,
    selection: &SamplingChannelSelection,
) -> Result<ResolvedSamplingChannelSelection, SamplingSelectionError> {
    resolve_selection(graph_name, selection, true)
}

/// Resolve with an additive default plus graph-specific list. This is useful
/// for callers that intentionally construct a shared baseline catalogue; the
/// runtime TOML contract uses replacement semantics above.
pub fn resolve_sampling_channel_selection_with_defaults(
    graph_name: &str,
    selection: &SamplingChannelSelection,
) -> Result<ResolvedSamplingChannelSelection, SamplingSelectionError> {
    resolve_selection(graph_name, selection, false)
}

/// Resolve using the explicit graph entry as a replacement for the defaults.
///
/// This is retained as a named compatibility entry point for callers that
/// want to state the replacement policy explicitly.
pub fn resolve_sampling_channel_selection_replacing_default(
    graph_name: &str,
    selection: &SamplingChannelSelection,
) -> Result<ResolvedSamplingChannelSelection, SamplingSelectionError> {
    resolve_selection(graph_name, selection, true)
}

fn contains_surface_map(map: &SamplingMapDefinition) -> bool {
    match map {
        SamplingMapDefinition::Surface(_) | SamplingMapDefinition::Cut(_) => true,
        SamplingMapDefinition::Product(maps)
        | SamplingMapDefinition::Intersect(maps)
        | SamplingMapDefinition::Then(maps) => maps.iter().any(contains_surface_map),
        SamplingMapDefinition::PhaseSpace(map)
        | SamplingMapDefinition::Left(map)
        | SamplingMapDefinition::Right(map)
        | SamplingMapDefinition::AtCut { map, .. } => contains_surface_map(map),
        SamplingMapDefinition::Block { .. } => false,
        _ => false,
    }
}

fn resolve_selection(
    graph_name: &str,
    selection: &SamplingChannelSelection,
    replace_default: bool,
) -> Result<ResolvedSamplingChannelSelection, SamplingSelectionError> {
    let graph_name = graph_name.trim();
    if graph_name.is_empty() {
        return Err(SamplingSelectionError::EmptyGraphName);
    }

    let graph_selectors = selection.channel_selection.get(graph_name);
    let raw_selectors = if replace_default {
        graph_selectors
            .map_or(
                selection.default_channel_selection.as_slice(),
                Vec::as_slice,
            )
            .iter()
            .collect::<Vec<_>>()
    } else {
        selection
            .default_channel_selection
            .iter()
            .chain(graph_selectors.into_iter().flatten())
            .collect::<Vec<_>>()
    };
    let mut selectors = Vec::with_capacity(raw_selectors.len());
    if raw_selectors.is_empty() && graph_selectors.is_none() {
        selectors.push(SamplingChannelSelector::Preset(
            SamplingChannelPreset::OptimizedLmb,
        ));
    }
    for raw in raw_selectors {
        let selector = SamplingChannelSelector::parse(raw)?;
        if !selectors.contains(&selector) {
            selectors.push(selector);
        }
    }

    let definitions = selection
        .channel_definitions
        .get(graph_name)
        .cloned()
        .unwrap_or_default();
    let available = definitions.keys().cloned().collect::<Vec<_>>();
    let mut named_channels = Vec::new();
    for selector in &selectors {
        let SamplingChannelSelector::Named(name) = selector else {
            continue;
        };
        let Some(definition) = definitions.get(name) else {
            return Err(SamplingSelectionError::MissingChannelDefinition {
                graph: graph_name.to_owned(),
                channel: name.clone(),
                available: available.clone(),
            });
        };
        if definition.parent_lmb.is_empty() {
            return Err(SamplingSelectionError::MissingParentLmb {
                graph: graph_name.to_owned(),
                channel: name.clone(),
            });
        }
        let mut parent_seen = BTreeSet::new();
        if definition
            .parent_lmb
            .iter()
            .any(|edge| !parent_seen.insert(*edge))
        {
            return Err(SamplingSelectionError::InvalidChannelDefinition {
                graph: graph_name.to_owned(),
                channel: name.clone(),
                error: format!(
                    "parent_lmb {:?} must contain unique ordered edge ids",
                    definition.parent_lmb
                ),
            });
        }
        let map = SamplingMapDefinition::parse(&definition.around).map_err(|error| {
            SamplingSelectionError::InvalidChannelDefinition {
                graph: graph_name.to_owned(),
                channel: name.clone(),
                error: error.to_string(),
            }
        })?;
        let singularity_proxy = definition
            .singularity_proxy
            .as_deref()
            .map(|proxy| {
                try_parse!(proxy.trim()).map_err(|error| {
                    SamplingSelectionError::InvalidChannelDefinition {
                        graph: graph_name.to_owned(),
                        channel: name.clone(),
                        error: format!(
                            "failed to parse singularity_proxy `{proxy}` with Symbolica: {error}"
                        ),
                    }
                })
            })
            .transpose()?;
        if contains_surface_map(&map) {
            if definition.subspace_lmb.is_empty() {
                return Err(SamplingSelectionError::MissingSubspaceLmb {
                    graph: graph_name.to_owned(),
                    channel: name.clone(),
                });
            }
            let mut seen = BTreeSet::new();
            if definition
                .subspace_lmb
                .iter()
                .any(|edge| !seen.insert(*edge) || !definition.parent_lmb.contains(edge))
            {
                return Err(SamplingSelectionError::InvalidChannelDefinition {
                    graph: graph_name.to_owned(),
                    channel: name.clone(),
                    error: format!(
                        "subspace_lmb {:?} must be unique and contained in parent_lmb {:?}",
                        definition.subspace_lmb, definition.parent_lmb
                    ),
                });
            }
        }
        let mut channel = ResolvedNamedSamplingChannel {
            name: name.clone(),
            definition: definition.clone(),
            map,
            singularity_proxy,
            blocks: Vec::new(),
        };
        channel.resolve_blocks().map_err(|error| {
            SamplingSelectionError::InvalidChannelDefinition {
                graph: graph_name.to_owned(),
                channel: name.clone(),
                error: error.to_string(),
            }
        })?;
        named_channels.push(channel);
    }

    Ok(ResolvedSamplingChannelSelection {
        graph_name: graph_name.to_owned(),
        selectors,
        named_channels,
    })
}

/// Return graph names having explicit selection entries, in deterministic order.
pub fn explicitly_selected_graphs(selection: &SamplingChannelSelection) -> Vec<&str> {
    selection
        .channel_selection
        .keys()
        .map(String::as_str)
        .collect()
}

/// Return a copy of graph definitions for diagnostics and catalogue builders.
pub fn graph_channel_definitions(
    selection: &SamplingChannelSelection,
    graph_name: &str,
) -> BTreeMap<String, SamplingChannelDefinition> {
    selection
        .channel_definitions
        .get(graph_name)
        .cloned()
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::initialisation::test_initialise;
    use crate::integrands::process::{SamplingCutSide, SamplingSupport};
    use crate::settings::runtime::ParameterizationSettings;

    #[test]
    fn native_lu_host_reuses_exact_sources_without_replacing_selected_authority() -> Result<()> {
        use super::super::sampling_context::SamplingLUHostPlan;
        use crate::{
            cff::{VertexSet, esurface::Esurface},
            dot,
            graph::{Graph, parse::from_dot::IntoGraph},
            momentum::{FourMomentum, sample::ExternalFourMomenta},
            processes::{CutGroupId, CutId},
            utils::{ArbPrec, QuadFloat},
        };
        use linnet::half_edge::involution::EdgeIndex;
        use std::cell::Cell;

        test_initialise()?;
        let graph: Graph = dot!(digraph host_source {
            ext [style=invis]
            node [num=1]
            edge [num=1 mass=1]
            ext -> a:0 [id=0]
            a -> b [id=1 lmb_id=0]
            a -> b [id=2]
            b:1 -> ext [id=3]
        })?;
        let surface = Esurface {
            energies: vec![EdgeIndex(1), EdgeIndex(2)],
            external_shift: vec![(EdgeIndex(0), -1)],
            vertex_set: VertexSet::dummy(),
        };
        fn check<T: FloatLike>(graph: &Graph, surface: &Esurface) -> Result<()> {
            let one = F::<T>::default().one();
            let zero = one.zero();
            let e_cm = one.from_usize(4);
            let masses = graph.underlying.new_edgevec_from_iter([
                zero.clone(),
                one.clone(),
                one.clone(),
                zero.clone(),
            ])?;
            let externals = ExternalFourMomenta::from_iter((0..2).map(|_| {
                FourMomentum::from_args(e_cm.clone(), zero.clone(), zero.clone(), zero.clone())
            }));
            let plan = Arc::new(SamplingLUHostPlan {
                graph_name: "host_source".into(),
                cut_group_id: CutGroupId::from(0),
                representative_cut_id: CutId(0),
                parent_lmb: graph
                    .loop_momentum_basis
                    .loop_edges
                    .iter()
                    .map(|edge| edge.0)
                    .collect(),
                required_prior_lmb: graph
                    .loop_momentum_basis
                    .loop_edges
                    .iter()
                    .map(|edge| edge.0)
                    .collect(),
            });
            let calls = Cell::new(0);
            let identities = std::cell::RefCell::new(Vec::new());
            let prepare =
                |prior: &[T],
                 diagnostics: &mut RadialRootDiagnostics,
                 identity: &crate::utils::newton_solver::RadialRootIdentity| {
                    calls.set(calls.get() + 1);
                    identities.borrow_mut().push(identity.to_string());
                    let loops = LoopMomenta::from_iter([ThreeMomentum::new(
                        F(prior[0].clone()),
                        F(prior[1].clone()),
                        F(prior[2].clone()),
                    )]);
                    surface
                        .solve_lu_cut(
                            &loops,
                            &externals,
                            &masses,
                            &graph.loop_momentum_basis,
                            &e_cm,
                            diagnostics,
                            identity,
                        )
                        .map_err(|error| eyre!("{error:?}"))
                };
            let prior = vec![one.0.clone(), zero.0.clone(), zero.0.clone()];
            let mut metadata = EvaluationMetaData::new_empty();
            let mut row =
                SamplingChannelRuntimeContexts::for_draw(2, 7, SamplingChannelId(0), &mut metadata);
            row.begin_selected(SamplingChannelId(0))?;
            {
                let mut context = row.for_channel(SamplingChannelId(0))?;
                let host =
                    context.prepare_lu_host(&plan, prior.clone(), |d, i| prepare(&prior, d, i))?;
                assert!(
                    (&host.solution.solution - one.from_usize(3).sqrt()).abs()
                        < one.epsilon() * one.from_usize(64)
                );
                assert_eq!(host.source.graph_id, 7);
                assert_eq!(host.prior, prior);
                context
                    .reborrow(&[], Some(1))
                    .prepare_lu_host(&plan, prior.clone(), |d, i| prepare(&prior, d, i))?;
            }
            assert_eq!(calls.get(), 1);
            let different = vec![
                (&one + one.epsilon() * one.from_usize(8)).0,
                zero.0.clone(),
                zero.0.clone(),
            ];
            assert_ne!(different, prior);
            let error = row
                .for_channel(SamplingChannelId(0))?
                .prepare_lu_host(&plan, different.clone(), |d, i| prepare(&different, d, i))
                .unwrap_err();
            assert!(matches!(
                error.downcast_ref::<SamplingEvaluationError>(),
                Some(SamplingEvaluationError::UncertainGeometry { .. })
            ));
            assert_eq!(calls.get(), 1);
            row.seal_selected();
            // A supplied-point inverse has distinct represented prerequisites.
            // Its actual equation is solved, while selected source authority stays fixed.
            row.for_channel(SamplingChannelId(0))?.prepare_lu_host(
                &plan,
                different.clone(),
                |d, i| prepare(&different, d, i),
            )?;
            row.for_channel(SamplingChannelId(1))?.prepare_lu_host(
                &plan,
                different.clone(),
                |d, i| prepare(&different, d, i),
            )?;
            assert_eq!(calls.get(), 2);
            assert_eq!(row.prepared_lu_hosts.len(), 2);
            let selected = row.selected_hosts();
            assert_eq!(selected.len(), 1);
            assert_eq!(selected[0].prior, prior);
            assert_eq!(row.prepared_lu_hosts[1].prior, different);
            let inverse_expected = one.from_usize(3).sqrt() / F(different[0].clone());
            assert!(
                (&row.prepared_lu_hosts[1].solution.solution - inverse_expected).abs()
                    < one.epsilon() * one.from_usize(64)
            );
            assert!(identities.borrow()[0].contains("native selected"));
            assert!(identities.borrow()[1].contains("native inverse"));
            // A standalone partition owns no selected prefix; a new direct raw
            // initial inverse establishes its own source, even in this same row object.
            let mut partition_row = SamplingChannelRuntimeContexts::new(2);
            partition_row
                .for_channel(SamplingChannelId(0))?
                .prepare_lu_host(&plan, prior.clone(), |d, i| prepare(&prior, d, i))?;
            assert!(partition_row.selected_hosts().is_empty());
            row.begin_selected(SamplingChannelId(0))?;
            row.for_channel(SamplingChannelId(0))?.prepare_lu_host(
                &plan,
                different.clone(),
                |d, i| prepare(&different, d, i),
            )?;
            row.seal_selected();
            assert_eq!(row.selected_hosts()[0].prior, different);
            assert_eq!(calls.get(), 4);
            Ok(())
        }
        check::<f64>(&graph, &surface)?;
        check::<QuadFloat>(&graph, &surface)?;
        check::<ArbPrec>(&graph, &surface)
    }

    #[test]
    fn joint_bridge_requires_full_coverage_and_gates_actual_inverse_support() -> Result<()> {
        let program = std::thread::Builder::new()
            .stack_size(64 * 1024 * 1024)
            .spawn(|| {
                test_initialise().unwrap();
                SharedEnergyJointMap::<f64>::compile_program()
            })?
            .join()
            .unwrap()?;
        let compact = SharedEnergyJointMap::new(
            Arc::new(|_| {
                Ok(super::super::SharedEnergyJointGeometry {
                    shifts: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                    masses: [1.0; 3],
                    energy_sums: [4.0; 2],
                })
            }),
            0,
            0.125,
            1.0,
            1.0,
            program,
        )?;
        let channel = |name: &str, map| CompiledSamplingChannel {
            name: name.into(),
            master_graph: "shared-energy-pair".into(),
            basis_id: None,
            definition: SamplingMapDefinition::Intersect(vec![
                SamplingMapDefinition::Surface(vec![0, 1]),
                SamplingMapDefinition::Surface(vec![0, 2]),
            ]),
            embedded_edges: vec![0],
            map,
            singularity_proxy: Some(SamplingScoreFunction::from_positive_function(|_| {
                Ok(Some(1.0))
            })),
        };
        let compact = channel("compact", CompiledSamplingMap::Joint(compact));
        assert!(matches!(
            SamplingChannelBridge::new(vec![compact.clone()]),
            Err(SamplingChannelBridgeError::MissingFullSupport)
        ));
        let ordinary = channel(
            "ordinary",
            CompiledSamplingMap::Surface(SurfaceRadialMap::absent(3, vec![0.0; 3], 1.0, 1.0)?),
        );
        for mode in [
            SamplingPartitionMode::MapDensity,
            SamplingPartitionMode::SingularityProxy,
        ] {
            let bridge = SamplingChannelBridge::new_with_partition_mode(
                vec![compact.clone(), ordinary.clone()],
                mode,
            )?;
            let mapped = bridge.forward(SamplingChannelId(0), &[0.43, 0.31, 0.37])?;
            assert!(
                mapped
                    .map
                    .diagnostics
                    .iter()
                    .any(|value| value.contains("certified normal disk"))
            );
            assert!(mapped.partition.weights.iter().all(|weight| *weight > 0.0));
            assert!((mapped.partition.weight_sum() - 1.0).abs() < 1.0e-14);
            assert!(
                bridge
                    .inverse(SamplingChannelId(0), &mapped.raw_coordinates)?
                    .is_some()
            );

            // Original energies put this regular raw point far outside the
            // small normal disk; this is geometric exclusion, not a failed root.
            let outside: [f64; 3] = [0.2, 0.3, 0.4];
            let e0 = (1.0 + outside.iter().map(|x| x * x).sum::<f64>()).sqrt();
            let h = e0 + (1.0 + 1.2_f64.powi(2) + 0.3_f64.powi(2) + 0.4_f64.powi(2)).sqrt() - 4.0;
            let z = e0 + (1.0 + 0.2_f64.powi(2) + 1.3_f64.powi(2) + 0.4_f64.powi(2)).sqrt() - 4.0;
            assert!(h.hypot(z) > 0.125);
            assert!(bridge.inverse(SamplingChannelId(0), &outside)?.is_none());
            assert_eq!(bridge.partition(&outside)?.weights, vec![0.0, 1.0]);
            assert!(bridge.inverse(SamplingChannelId(1), &outside)?.is_some());
            if mode == SamplingPartitionMode::MapDensity {
                let report = SamplingChannelBridgeAcceptanceReport::normalized_gaussian(
                    &bridge,
                    4096,
                    0.9,
                    &[0.4, -0.3, 0.2],
                )?;
                assert_eq!(report.finite_sample_count, 8192);
                assert!((report.normalization - 1.0).abs() < 0.02, "{report:?}");
                assert!(
                    (report.second_moment / report.expected_second_moment - 1.0).abs() < 0.04,
                    "{report:?}"
                );
            }
            // The constant proxy checks support only; its corner power is not
            // sufficient evidence for the density-based damping oracle below.
        }
        Ok(())
    }

    #[test]
    fn joint_density_partition_bounds_inverse_radius_from_either_channel() -> Result<()> {
        let program = std::thread::Builder::new()
            .stack_size(64 * 1024 * 1024)
            .spawn(|| {
                test_initialise().unwrap();
                SharedEnergyJointMap::<f64>::compile_program()
            })?
            .join()
            .unwrap()?;
        fn check<T: FloatLike>(program: SamplingExpressionEvaluator) -> Result<()> {
            let one = F::<T>::default().one();
            let zero = one.zero();
            let compact = SharedEnergyJointMap::new(
                Arc::new(|_| {
                    let one = F::<T>::default().one();
                    let zero = one.zero();
                    Ok(super::super::SharedEnergyJointGeometry {
                        shifts: [
                            [one.0.clone(), zero.0.clone(), zero.0.clone()],
                            [zero.0.clone(), one.0.clone(), zero.0],
                        ],
                        masses: [one.0.clone(), one.0.clone(), one.0.clone()],
                        energy_sums: [one.from_usize(4).0, one.from_usize(4).0],
                    })
                }),
                0,
                (&one / one.from_usize(8)).0,
                one.0.clone(),
                1.0,
                program,
            )?;
            let channel = |name: &str, map| CompiledSamplingChannel {
                name: name.into(),
                master_graph: "shared-energy-pair".into(),
                basis_id: None,
                definition: SamplingMapDefinition::Intersect(vec![
                    SamplingMapDefinition::Surface(vec![0, 1]),
                    SamplingMapDefinition::Surface(vec![0, 2]),
                ]),
                embedded_edges: vec![0],
                map,
                singularity_proxy: None,
            };
            let bridge = SamplingChannelBridge::new(vec![
                channel("compact", CompiledSamplingMap::Joint(compact)),
                channel(
                    "ordinary",
                    CompiledSamplingMap::Surface(SurfaceRadialMap::absent(
                        3,
                        vec![zero.0.clone(); 3],
                        1.0,
                        1.0,
                    )?),
                ),
            ])?;
            // Independently evaluate the original two energy equations at the
            // supplied raw point. Never infer R from a Jacobian or inverse q.
            let radius = |point: &[T]| {
                let x = point.iter().cloned().map(F).collect::<Vec<_>>();
                let e0 = (&one + x.iter().fold(zero.clone(), |sum, x| sum + x.square())).sqrt();
                let e1 = (&one + (&x[0] + &one).square() + x[1].square() + x[2].square()).sqrt();
                let e2 = (&one + x[0].square() + (&x[1] + &one).square() + x[2].square()).sqrt();
                ((&e0 + e1 - one.from_usize(4)).square() + (e0 + e2 - one.from_usize(4)).square())
                    .sqrt()
            };
            let theta = one.from_usize(31) / one.from_usize(100);
            let phi = one.from_usize(37) / one.from_usize(100);
            let probe = bridge.forward(
                SamplingChannelId(0),
                &[(&one / one.from_usize(2)).0, theta.0.clone(), phi.0.clone()],
            )?;
            assert!(
                probe
                    .map
                    .diagnostics
                    .iter()
                    .any(|value| value.contains("certified normal disk"))
            );
            let rho = radius(&probe.raw_coordinates) * one.from_usize(2);
            // At h=z=0 the independent circle has k=31, D=163/2,
            // E0=60/31 + sqrt(D)/31*cos(phi). The limiting product is
            // 4*pi^2*rho * E0*(4-E0)^2/sqrt(31), for this alpha=1 frame.
            let pi = F(one.0.pi());
            let angle = one.from_usize(2) * &pi * &phi;
            let e0 = (one.from_usize(60)
                + (one.from_usize(163) / one.from_usize(2)).sqrt() * F(angle.0.cos()))
                / one.from_usize(31);
            let limit =
                one.from_usize(4) * pi.square() * &rho * &e0 * (one.from_usize(4) - &e0).square()
                    / one.from_usize(31).sqrt();
            let mut errors = Vec::new();
            for decade in 1..=4 {
                let cube = [
                    (&one / one.from_usize(10).powi(decade)).0,
                    theta.0.clone(),
                    phi.0.clone(),
                ];
                let point = bridge.forward(SamplingChannelId(0), &cube)?.raw_coordinates;
                let mut factors = Vec::new();
                for selected in 0..2 {
                    let inverse = bridge
                        .inverse(SamplingChannelId(selected), &point)?
                        .expect("both regular inverses exist");
                    let replay =
                        bridge.forward(SamplingChannelId(selected), &inverse.map.coordinates)?;
                    let factor = F(replay.selected_factor()?) / radius(&replay.raw_coordinates);
                    assert!(factor.0.is_finite() && factor > zero);
                    factors.push(factor);
                }
                assert!(
                    (&factors[0] / &factors[1] - &one).abs()
                        < one.epsilon().sqrt() * one.from_usize(1024)
                );
                errors.push((&factors[0] / &limit - &one).abs());
            }
            assert!(errors[3] < &errors[0] / one.from_usize(20), "{errors:?}");
            assert!(errors[3] < &one / one.from_usize(1000), "{errors:?}");
            Ok(())
        }
        check::<crate::utils::QuadFloat>(program.clone())?;
        check::<crate::utils::ArbPrec>(program)?;
        Ok(())
    }

    #[test]
    fn conditional_joint_bridge_uses_foreign_context_and_affine_inverse_support() -> Result<()> {
        let program = std::thread::Builder::new()
            .stack_size(64 * 1024 * 1024)
            .spawn(|| {
                test_initialise().unwrap();
                SharedEnergyJointMap::<f64>::compile_program()
            })?
            .join()
            .unwrap()?;
        let conditional = |scale: f64| -> Result<CompiledSamplingMap> {
            let joint = SharedEnergyJointMap::new(
                Arc::new(|prior| {
                    Ok(super::super::SharedEnergyJointGeometry {
                        shifts: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                        masses: [1.0; 3],
                        energy_sums: [4.0 + prior[0] / 5.0, 4.0 + prior[1] / 5.0],
                    })
                }),
                3,
                0.125,
                1.0,
                1.0,
                program.clone(),
            )?;
            let active = SamplingMapEmbedding::from_composition(
                SamplingMapComposition::then(vec![Box::new(joint)])?,
                vec![0, 1, 2],
            )?
            .with_context_transform(Arc::new(move |context| {
                let prior = context.previous;
                if prior.len() != 3 {
                    return Err(eyre!(
                        "joint fixture needs its three declared prerequisites"
                    ));
                }
                let physical = prior.iter().map(|x| scale * x).collect::<Vec<_>>();
                let frame = SamplingMapAffine::new(
                    vec![
                        vec![scale, scale / 4.0, 0.0],
                        vec![0.0, 1.5 * scale, 0.0],
                        vec![0.0, 0.0, 0.75 * scale],
                    ],
                    physical.clone(),
                )?;
                Ok((physical, frame))
            }));
            assert!(active.contract().requires_context);
            let map = SamplingMapEmbedding::from_composition(
                SamplingMapComposition::then(vec![
                    Box::new(SurfaceRadialMap::absent(3, vec![0.0; 3], 1.0, 1.0)?),
                    Box::new(active),
                ])?,
                vec![3, 4, 5, 0, 1, 2],
            )?;
            assert_eq!(map.contract().support, SamplingSupport::Restricted);
            assert!(!map.contract().requires_context);
            Ok(CompiledSamplingMap::Embedded(map))
        };
        let channel = |name: &str, map| CompiledSamplingChannel {
            name: name.into(),
            master_graph: "conditional-shared-energy-pair".into(),
            basis_id: None,
            definition: SamplingMapDefinition::Intersect(vec![
                SamplingMapDefinition::Surface(vec![0, 1]),
                SamplingMapDefinition::Surface(vec![0, 2]),
            ]),
            embedded_edges: vec![0, 1],
            map,
            singularity_proxy: Some(SamplingScoreFunction::from_positive_function(|_| {
                Ok(Some(1.0))
            })),
        };
        let channels = vec![
            channel("scale-one", conditional(1.0)?),
            channel("scale-two", conditional(2.0)?),
            channel(
                "ordinary",
                CompiledSamplingMap::Surface(SurfaceRadialMap::absent(6, vec![0.0; 6], 1.0, 1.0)?),
            ),
        ];
        for mode in [
            SamplingPartitionMode::MapDensity,
            SamplingPartitionMode::SingularityProxy,
        ] {
            let bridge = SamplingChannelBridge::new_with_partition_mode(channels.clone(), mode)?;
            for selected in 0..2 {
                let cube = [0.1, 0.31, 0.43, 0.47, 0.29, 0.37];
                let mapped = bridge.forward(SamplingChannelId(selected), &cube)?;
                assert!(
                    mapped
                        .map
                        .diagnostics
                        .iter()
                        .any(|d| d.contains("certified normal disk"))
                );
                let inverse = bridge
                    .inverse(SamplingChannelId(selected), &mapped.raw_coordinates)?
                    .expect("selected regular joint inverse");
                for (actual, expected) in inverse.map.coordinates.iter().zip(cube) {
                    assert!((actual - expected).abs() < 1.0e-8);
                }
                assert!(
                    bridge
                        .inverse(SamplingChannelId(1 - selected), &mapped.raw_coordinates)?
                        .is_none()
                );
                assert_eq!(mapped.partition.weights[1 - selected], 0.0);
                assert!(
                    mapped.partition.weights[selected] > 0.0 && mapped.partition.weights[2] > 0.0
                );
                assert!((mapped.partition.weight_sum() - 1.0).abs() < 1.0e-14);

                // Undo the known shear and permutation independently. Each
                // chart must use its own scaled raw complement, including in
                // a foreign inverse; a shared selected context gives the wrong support.
                let scale = (selected + 1) as f64;
                let raw = &mapped.raw_coordinates;
                let y = (raw[1] / scale - raw[4]) / 1.5;
                let x = raw[0] / scale - raw[3] - y / 4.0;
                let z = (raw[2] / scale - raw[5]) / 0.75;
                let e0 = (1.0 + x * x + y * y + z * z).sqrt();
                let h = e0 + (1.0 + (x + 1.0).powi(2) + y * y + z * z).sqrt()
                    - (4.0 + scale * raw[3] / 5.0);
                let zeta = e0 + (1.0 + x * x + (y + 1.0).powi(2) + z * z).sqrt()
                    - (4.0 + scale * raw[4] / 5.0);
                assert!(h.hypot(zeta) < 0.125);
            }
            // Generating rows remain distinct while each selected/foreign map
            // revisits the same compiled path in forward and inverse. Sealing
            // forbids an incidental later lane from filling a missing record.
            let cube = [0.1, 0.31, 0.43, 0.47, 0.29, 0.37];
            let mut metadata = EvaluationMetaData::new_empty();
            metadata.sampling_proposal_policies.begin_collection();
            for generator in 0..2 {
                let mut contexts = SamplingChannelRuntimeContexts::for_draw(
                    3,
                    7,
                    SamplingChannelId(generator),
                    &mut metadata,
                );
                bridge.forward_with_runtime_contexts(
                    SamplingChannelId(generator),
                    &cube,
                    &mut contexts,
                )?;
            }
            assert_eq!(metadata.sampling_proposal_policies.len(), 4);
            metadata.sampling_proposal_policies.seal();
            let sealed = metadata.sampling_proposal_policies.clone();
            for generator in 0..2 {
                let mut contexts = SamplingChannelRuntimeContexts::for_draw(
                    3,
                    7,
                    SamplingChannelId(generator),
                    &mut metadata,
                );
                let mapped = bridge.forward_with_runtime_contexts(
                    SamplingChannelId(generator),
                    &cube,
                    &mut contexts,
                )?;
                assert!(
                    bridge
                        .inverse_with_runtime_contexts(
                            SamplingChannelId(generator),
                            &mapped.raw_coordinates,
                            &mut contexts
                        )?
                        .is_some()
                );
            }
            assert_eq!(metadata.sampling_proposal_policies, sealed);
            let mut empty = EvaluationMetaData::new_empty();
            let mut contexts =
                SamplingChannelRuntimeContexts::for_draw(3, 7, SamplingChannelId(0), &mut empty);
            let error = bridge
                .forward_with_runtime_contexts(SamplingChannelId(0), &cube, &mut contexts)
                .unwrap_err();
            assert!(matches!(
                error.downcast_ref::<SamplingEvaluationError>(),
                Some(SamplingEvaluationError::UncertainGeometry { .. })
            ));
            assert!(format!("{error:#}").contains("conditional-shared-energy-pair"));
            assert!(empty.sampling_proposal_policies.is_empty());
        }
        Ok(())
    }

    #[test]
    fn bridge_certifies_selected_inverse_density_in_both_partition_modes() {
        use super::super::sampling_maps::SamplingEvaluationError;
        use crate::utils::{ArbPrec, QuadFloat};
        use std::sync::atomic::{AtomicUsize, Ordering};

        // Isolate a faulty inverse determinant from coordinate round trips:
        // the underlying affine map remains exactly invertible with zero residual.
        #[derive(Debug)]
        struct InverseDensityProbe<T: FloatLike> {
            map: SamplingMapAffine<T>,
            multiplier: T,
            calls: Arc<AtomicUsize>,
            unavailable: bool,
        }
        impl<T: FloatLike> SamplingMapComponent<T> for InverseDensityProbe<T> {
            fn dimensions(&self) -> usize {
                self.map.dimensions()
            }
            fn output_dimensions(&self) -> usize {
                self.map.output_dimensions()
            }
            fn contract(&self) -> SamplingMapContract {
                self.map.contract()
            }
            fn name(&self) -> &'static str {
                "inverse_density_probe"
            }
            fn forward(
                &self,
                coordinates: &[T],
                context: &mut SamplingMapContext<'_, T>,
            ) -> Result<SamplingMapEvaluation<T>> {
                self.map.forward(coordinates, context)
            }
            fn inverse(
                &self,
                point: &[T],
                context: &mut SamplingMapContext<'_, T>,
            ) -> Result<Option<SamplingMapEvaluation<T>>> {
                self.calls.fetch_add(1, Ordering::Relaxed);
                if self.unavailable {
                    return Err(eyre!("foreign inverse must not be evaluated"));
                }
                let mut inverse = self.map.inverse(point, context.previous)?;
                inverse.inverse_jacobian =
                    (F(inverse.inverse_jacobian) * F(self.multiplier.clone())).0;
                Ok(Some(inverse))
            }
        }
        fn check<T: FloatLike>(exponent: i32) {
            let one = F::<T>::default().one();
            let delta = one.from_i64(10).powi(exponent);
            let coordinates = vec![(&one / one.from_i64(3)).0; 3];
            let selected_calls = Arc::new(AtomicUsize::new(0));
            let foreign_calls = Arc::new(AtomicUsize::new(0));
            let channel = |name: &str, multiplier: F<T>, calls, unavailable| {
                let map = SamplingMapAffine::new(
                    (0..3)
                        .map(|i| {
                            (0..3)
                                .map(|j| if i == j { one.0.clone() } else { one.zero().0 })
                                .collect()
                        })
                        .collect(),
                    vec![one.zero().0; 3],
                )
                .unwrap();
                let map = SamplingMapComposition::product(vec![Box::new(InverseDensityProbe {
                    map,
                    multiplier: multiplier.0,
                    calls,
                    unavailable,
                })])
                .unwrap();
                CompiledSamplingChannel {
                    name: name.into(),
                    master_graph: "G".into(),
                    basis_id: None,
                    definition: SamplingMapDefinition::Lmb(vec![1]),
                    embedded_edges: vec![1],
                    map: CompiledSamplingMap::Embedded(
                        SamplingMapEmbedding::from_composition(map, vec![0, 1, 2]).unwrap(),
                    ),
                    singularity_proxy: Some(SamplingScoreFunction::from_log_function(|_| {
                        Ok(Some(F::<T>::default().zero().0))
                    })),
                }
            };
            for mode in [
                SamplingPartitionMode::MapDensity,
                SamplingPartitionMode::SingularityProxy,
            ] {
                let channels = vec![
                    channel("selected", &one + &delta, selected_calls.clone(), false),
                    channel(
                        "foreign",
                        one.clone(),
                        foreign_calls.clone(),
                        mode == SamplingPartitionMode::SingularityProxy,
                    ),
                ];
                let bridge =
                    SamplingChannelBridge::new_with_partition_mode(channels, mode).unwrap();
                for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
                    assert!(
                        bridge
                            .clone()
                            .with_relative_density_tolerance(invalid)
                            .is_err()
                    );
                }
                selected_calls.store(0, Ordering::Relaxed);
                foreign_calls.store(0, Ordering::Relaxed);
                let loose = 10.0_f64.powi(exponent + 1);
                let evaluation = bridge
                    .clone()
                    .with_relative_density_tolerance(loose)
                    .unwrap()
                    .forward(SamplingChannelId(0), &coordinates)
                    .unwrap();
                assert_eq!(evaluation.map.residual, one.zero().0);
                assert_eq!(
                    selected_calls.load(Ordering::Relaxed),
                    1,
                    "selected inverse must be reused from map-density partition"
                );
                assert_eq!(
                    foreign_calls.load(Ordering::Relaxed),
                    usize::from(mode == SamplingPartitionMode::MapDensity)
                );
                let strict = 10.0_f64.powi(exponent - 1);
                let error = bridge
                    .clone()
                    .with_relative_density_tolerance(strict)
                    .unwrap()
                    .forward(SamplingChannelId(0), &coordinates)
                    .unwrap_err();
                assert!(matches!(
                    error.downcast_ref::<SamplingEvaluationError>(),
                    Some(SamplingEvaluationError::Unrepresentable {
                        operation: "sampling inverse density consistency",
                        ..
                    })
                ));
                assert_eq!(
                    bridge.forward(SamplingChannelId(0), &coordinates).is_ok(),
                    delta < one.epsilon().sqrt()
                );
            }
        }
        check::<f64>(-6);
        check::<QuadFloat>(-20);
        check::<ArbPrec>(-20);
    }

    #[test]
    fn selected_sampling_factor_rejects_underflow_before_physical_cancellation() {
        use super::super::sampling_maps::{SamplingEvaluationError, SamplingSupport};
        use crate::utils::ArbPrec;
        let evaluation = SamplingChannelBridgeEvaluation {
            prepared_lu_hosts: Vec::new(),
            channel_id: SamplingChannelId(0),
            channel_name: "focused".into(),
            raw_coordinates: vec![1.0; 3],
            map: SamplingMapEvaluation {
                coordinates: vec![0.5; 3],
                point: vec![1.0; 3],
                jacobian: 1.0e-200,
                inverse_jacobian: 1.0e200,
                residual: 0.0,
                support: SamplingSupport::Full,
                diagnostics: Vec::new(),
            },
            partition: SamplingPartition {
                mode: SamplingPartitionMode::SingularityProxy,
                log_scores: vec![Some(1.0e-200_f64.ln()), Some(0.0)],
                weights: vec![1.0e-200, 1.0],
                log_denominator: 0.0,
            },
        };
        assert!(evaluation.map.jacobian.is_finite() && evaluation.map.jacobian > 0.0);
        assert!(
            evaluation.partition.weights[0].is_finite() && evaluation.partition.weights[0] > 0.0
        );
        let error = evaluation.selected_factor().unwrap_err();
        assert!(matches!(
            error.downcast_ref::<SamplingEvaluationError>(),
            Some(SamplingEvaluationError::Unrepresentable { .. })
        ));
        let native = SamplingChannelBridgeEvaluation {
            prepared_lu_hosts: Vec::new(),
            channel_id: evaluation.channel_id,
            channel_name: evaluation.channel_name,
            raw_coordinates: evaluation
                .raw_coordinates
                .into_iter()
                .map(|value| F::<ArbPrec>::from_f64(value).0)
                .collect(),
            map: SamplingMapEvaluation {
                coordinates: evaluation
                    .map
                    .coordinates
                    .into_iter()
                    .map(|value| F::<ArbPrec>::from_f64(value).0)
                    .collect(),
                point: evaluation
                    .map
                    .point
                    .into_iter()
                    .map(|value| F::<ArbPrec>::from_f64(value).0)
                    .collect(),
                jacobian: F::<ArbPrec>::from_f64(evaluation.map.jacobian).0,
                inverse_jacobian: F::<ArbPrec>::from_f64(evaluation.map.inverse_jacobian).0,
                residual: F::<ArbPrec>::from_f64(0.0).0,
                support: evaluation.map.support,
                diagnostics: Vec::new(),
            },
            partition: SamplingPartition {
                mode: evaluation.partition.mode,
                log_scores: evaluation
                    .partition
                    .log_scores
                    .into_iter()
                    .map(|value| value.map(|value| F::<ArbPrec>::from_f64(value).0))
                    .collect(),
                weights: evaluation
                    .partition
                    .weights
                    .into_iter()
                    .map(|value| F::<ArbPrec>::from_f64(value).0)
                    .collect(),
                log_denominator: F::<ArbPrec>::from_f64(0.0).0,
            },
        };
        let factor = F(native.selected_factor().unwrap());
        assert!(factor > factor.zero());
        let expected = factor.one().from_usize(10).powi(-400);
        assert!((&factor / expected - factor.one()).abs() < factor.one().from_usize(10).powi(-15));
    }

    #[test]
    fn native_bridge_preserves_composed_points_and_foreign_inverse_densities() {
        fn check<T: FloatLike>() {
            let one = F::<T>::default().one();
            let zero = one.zero();
            let three = one.from_i64(3);
            let two = one.from_i64(2);
            let active =
                SurfaceRadialMap::new(3, vec![zero.0.clone(); 3], Some(three.0.clone()), 2.0, 2.0)
                    .unwrap();
            let complement =
                SurfaceRadialMap::absent(3, vec![zero.0.clone(); 3], 2.0, 1.0).unwrap();
            let product =
                SamplingMapComposition::product(vec![Box::new(active), Box::new(complement)])
                    .unwrap();
            let product =
                SamplingMapEmbedding::from_composition(product, vec![3, 4, 5, 0, 1, 2]).unwrap();
            let target = three.clone();
            let implicit = ImplicitSurfaceRadialMap::new(
                6,
                vec![zero.0.clone(); 6],
                2.0,
                2.0,
                Arc::new(move |_, radius: T| {
                    let radius = F(radius);
                    Ok((
                        (radius.square() - target.square()).0,
                        (radius.from_i64(2) * radius).0,
                    ))
                }),
            )
            .unwrap();
            let first =
                SurfaceRadialMap::new(3, vec![zero.0.clone(); 3], Some(three.0.clone()), 2.0, 1.0)
                    .unwrap();
            let conditional = ImplicitSurfaceRadialMap::new(
                3,
                vec![zero.0.clone(); 3],
                2.0,
                1.0,
                Arc::new(|_, radius: T| {
                    let radius = F(radius);
                    Ok(((&radius - radius.from_i64(2)).0, radius.one().0))
                }),
            )
            .unwrap()
            .with_context_preparer(Arc::new(|context: &[T]| {
                Ok((
                    context[..3].to_vec(),
                    crate::integrands::process::PreparedSurfaceStatus::existing(None)?,
                ))
            }));
            let ordered = SamplingMapEmbedding::from_composition(
                SamplingMapComposition::then(vec![Box::new(first), Box::new(conditional)]).unwrap(),
                (0..6).collect(),
            )
            .unwrap();
            let channels = [
                CompiledSamplingMap::Embedded(product),
                CompiledSamplingMap::ImplicitSurface(implicit),
                CompiledSamplingMap::Embedded(ordered),
            ]
            .into_iter()
            .enumerate()
            .map(|(index, map)| CompiledSamplingChannel {
                name: format!("native_{index}"),
                master_graph: "G".to_owned(),
                basis_id: None,
                definition: SamplingMapDefinition::Surface(vec![1, 2]),
                embedded_edges: vec![1, 2],
                map,
                singularity_proxy: None,
            })
            .collect();
            let bridge = SamplingChannelBridge::new(channels).unwrap();
            let coordinates =
                [13, 23, 31, 43, 59, 71].map(|value| (one.from_i64(value) / one.from_i64(100)).0);
            let tolerance = one.epsilon() * one.from_i64(1000000);
            for selected in 0..3 {
                let evaluation = bridge
                    .forward(SamplingChannelId(selected), &coordinates)
                    .unwrap();
                let densities = bridge
                    .channels()
                    .iter()
                    .map(|channel| {
                        channel
                            .inverse(&evaluation.raw_coordinates)
                            .unwrap()
                            .expect("full-support inverse")
                            .inverse_jacobian
                    })
                    .map(F)
                    .collect::<Vec<_>>();
                let density_sum = densities
                    .iter()
                    .fold(zero.clone(), |sum, value| sum + value);
                for (index, density) in densities.iter().enumerate() {
                    let expected = density / &density_sum;
                    let actual = F(evaluation.partition.weights[index].clone());
                    assert!((actual - expected).abs() < tolerance);
                }
                let inverse = bridge
                    .inverse(SamplingChannelId(selected), &evaluation.raw_coordinates)
                    .unwrap()
                    .expect("full-support inverse");
                for (expected, actual) in coordinates.iter().zip(&inverse.map.coordinates) {
                    assert!((F(actual.clone()) - F(expected.clone())).abs() < tolerance);
                }
                assert!(
                    (F(evaluation.map.jacobian) * F(inverse.map.inverse_jacobian) - &one).abs()
                        < tolerance
                );
                assert!(
                    (evaluation
                        .partition
                        .weights
                        .iter()
                        .cloned()
                        .map(F)
                        .fold(zero.clone(), |sum, value| sum + value)
                        - &one)
                        .abs()
                        < tolerance
                );
            }
            // A perturbation smaller than binary64 resolution survives both
            // native radial geometry and the common foreign inverse scores.
            let tiny = &one / two.powi(80);
            if tiny > one.epsilon() * one.from_i64(128) {
                let radius = &three + &tiny;
                let native = SurfaceRadialMap::new(
                    3,
                    vec![zero.0.clone(); 3],
                    Some(radius.0.clone()),
                    2.0,
                    1.0,
                )
                .unwrap();
                assert!(F(native.threshold_radius().unwrap()) > three);
                assert_eq!(F(native.threshold_radius().unwrap()), radius);
            }
        }
        test_initialise().unwrap();
        check::<f64>();
        check::<crate::utils::QuadFloat>();
        check::<crate::utils::ArbPrec>();
    }

    fn definition(around: &str) -> SamplingChannelDefinition {
        SamplingChannelDefinition {
            radial_profile: None,
            around: around.to_owned(),
            subspace_lmb: vec![1, 2],
            parent_lmb: vec![1, 2],
            on_cut: vec![],
            singularity_proxy: None,
        }
    }

    #[test]
    fn resolved_host_blocks_preserve_exact_prerequisite_order_and_reject_wrapped_nesting() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["both".into(), "direct".into()],
            ..Default::default()
        };
        let mut both = definition(
            "then(block(lmb(2),phase_space(cut(9,7))),block(lmb(1),left(surface(5,3))),block(lmb(4),right(surface(8,6))))",
        );
        both.parent_lmb = vec![1, 2, 4];
        both.subspace_lmb.clear();
        let mut direct = definition("at_cut(cut(9,7),surface(5,3))");
        direct.parent_lmb = vec![1, 2, 4];
        direct.subspace_lmb = vec![1];
        let definitions = selection.channel_definitions.entry("G".into()).or_default();
        definitions.insert("both".into(), both);
        definitions.insert("direct".into(), direct);
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let both = resolved.named("both").unwrap();
        assert_eq!(both.blocks.len(), 3);
        assert_eq!(both.blocks[0].active_lmb, vec![2]);
        assert_eq!(both.blocks[0].remaining_lmb, vec![1, 4]);
        assert_eq!(both.blocks[1].preceding_lmb, vec![2]);
        assert_eq!(both.blocks[1].remaining_lmb, vec![4]);
        assert_eq!(both.blocks[2].preceding_lmb, vec![2, 1]);
        assert_eq!(both.blocks[1].target.host_cut(), Some([7, 9].as_slice()));
        assert_eq!(
            both.blocks[2].target.cut_side(),
            Some(SamplingCutSide::Right)
        );
        let direct = resolved.named("direct").unwrap();
        assert_eq!(direct.blocks[0].active_lmb, vec![2, 4]);
        assert_eq!(direct.blocks[1].preceding_lmb, vec![2, 4]);
        assert_eq!(direct.blocks[1].target.cut_side(), None);
        assert_eq!(
            direct
                .blocks
                .iter()
                .map(|block| block.active_lmb.len())
                .sum::<usize>(),
            3
        );
        selection.default_channel_selection = vec!["bad".into()];
        let mut bad =
            definition("then(complement(2),at_cut(cut(7,9),product(block(lmb(1),surface(3,5)))))");
        bad.subspace_lmb.clear();
        selection
            .channel_definitions
            .get_mut("G")
            .unwrap()
            .insert("bad".into(), bad);
        let error = resolve_sampling_channel_selection("G", &selection).unwrap_err();
        assert!(error.to_string().contains("nested block compositions"));
        let definitions = selection.channel_definitions.get_mut("G").unwrap();
        definitions.insert(
            "bad".into(),
            definition("block(lmb(1),block(lmb(2),surface(3)))"),
        );
        let error = resolve_sampling_channel_selection("G", &selection).unwrap_err();
        assert!(error.to_string().contains("different active edges"));
        selection.channel_definitions.get_mut("G").unwrap().insert(
            "bad".into(),
            definition("block(lmb(2,1),block(lmb(1,2),surface(3)))"),
        );
        let repeated = resolve_sampling_channel_selection("G", &selection).unwrap();
        assert_eq!(
            repeated.named("bad").unwrap().blocks[0].active_lmb,
            vec![1, 2]
        );
    }

    #[test]
    fn joint_targets_resolve_as_one_ordered_three_dimensional_block() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["joint".into()],
            ..Default::default()
        };
        for around in [
            "intersect(surface(6,2),surface(6,3))",
            "at_cut(cut(9,7),intersect(surface(6,2),surface(6,3)))",
            "then(complement(4),block(lmb(6),at_cut(cut(9,7),left(intersect(surface(6,2),surface(6,3))))))",
        ] {
            selection
                .channel_definitions
                .entry("G".into())
                .or_default()
                .insert(
                    "joint".into(),
                    SamplingChannelDefinition {
                        parent_lmb: vec![6, 4],
                        subspace_lmb: vec![6],
                        ..definition(around)
                    },
                );
            let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
            let blocks = &resolved.named("joint").unwrap().blocks;
            assert_eq!(blocks.len(), 2);
            assert_eq!(blocks[0].active_lmb, vec![4]);
            assert_eq!(blocks[1].active_lmb, vec![6]);
            assert_eq!(blocks[1].preceding_lmb, vec![4]);
            assert!(blocks[1].remaining_lmb.is_empty());
            assert_eq!(
                blocks[1].target.energy_edge_sets(),
                vec![[2, 6].as_slice(), [3, 6].as_slice()]
            );
            assert_eq!(
                blocks.iter().map(|b| 3 * b.active_lmb.len()).sum::<usize>(),
                6
            );
            if around.contains("cut") {
                assert_eq!(blocks[1].target.host_cut(), Some([7, 9].as_slice()));
            }
        }
        for (around, active, expected) in [
            ("intersect(surface(2))", vec![6], "exactly two distinct"),
            (
                "intersect(surface(2),surface(2))",
                vec![6],
                "exactly two distinct",
            ),
            (
                "intersect(surface(2),surface(3),surface(5))",
                vec![6],
                "exactly two distinct",
            ),
            (
                "intersect(surface(2),soft(3))",
                vec![6],
                "exactly two distinct",
            ),
            (
                "intersect(surface(2),at_cut(cut(9),surface(3)))",
                vec![6],
                "qualify the whole",
            ),
            (
                "intersect(surface(2),product(surface(3)))",
                vec![6],
                "exactly two distinct",
            ),
            (
                "intersect(surface(2),surface(3))",
                vec![6, 4],
                "one three-dimensional",
            ),
        ] {
            selection.channel_definitions.get_mut("G").unwrap().insert(
                "joint".into(),
                SamplingChannelDefinition {
                    parent_lmb: vec![6, 4],
                    subspace_lmb: active,
                    ..definition(around)
                },
            );
            let error = resolve_sampling_channel_selection("G", &selection).unwrap_err();
            assert!(error.to_string().contains(expected), "{around}: {error}");
        }
    }

    #[test]
    fn joint_catalogue_programs_and_registry_preserve_qualified_pair_identity() -> Result<()> {
        test_initialise()?;
        let pair = "intersect(surface(2,6),surface(3,6))";
        let selection = SamplingChannelSelection {
            default_channel_selection: vec!["joint".into(), "hosted".into(), "ordinary".into()],
            channel_definitions: BTreeMap::from([(
                "G".into(),
                BTreeMap::from([
                    (
                        "joint".into(),
                        SamplingChannelDefinition {
                            parent_lmb: vec![4, 6],
                            subspace_lmb: vec![6],
                            ..definition(pair)
                        },
                    ),
                    (
                        "hosted".into(),
                        SamplingChannelDefinition {
                            parent_lmb: vec![4, 6],
                            subspace_lmb: vec![6],
                            ..definition(&format!("at_cut(cut(9),{pair})"))
                        },
                    ),
                    (
                        "ordinary".into(),
                        SamplingChannelDefinition {
                            parent_lmb: vec![4, 6],
                            subspace_lmb: vec![],
                            ..definition("lmb(4,6)")
                        },
                    ),
                ]),
            )]),
            ..Default::default()
        };
        let resolved = resolve_sampling_channel_selection("G", &selection)?;
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let programs = std::thread::Builder::new()
            .stack_size(64 * 1024 * 1024)
            .spawn({
                let catalogue = catalogue.clone();
                move || catalogue.compile_programs(6, &HFunctionSettings::default())
            })?
            .join()
            .unwrap()?;
        assert_eq!(
            programs
                .iter()
                .filter(|(_, _, joint)| joint.is_some())
                .count(),
            2
        );
        assert!(
            programs
                .iter()
                .all(|(proxy, radial, _)| proxy.is_none() && radial.is_none())
        );
        let program = programs
            .iter()
            .find_map(|(_, _, joint)| joint.clone())
            .unwrap();
        let map = SharedEnergyJointMap::new(
            Arc::new(|_| {
                Ok(super::super::SharedEnergyJointGeometry {
                    shifts: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                    masses: [1.0; 3],
                    energy_sums: [4.0; 2],
                })
            }),
            3,
            0.125,
            1.0,
            1.0,
            program,
        )?;
        let mut context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![4, 6],
            ParameterizationSettings::default(),
            5.0,
            2,
        );
        context.physical_cut_ids.insert(vec![9], vec![0]);
        let error = catalogue.compile(&context, &programs).unwrap_err();
        let Some(SamplingChannelCompileError::MissingGeometry { target, .. }) =
            error.downcast_ref::<SamplingChannelCompileError>()
        else {
            panic!("{error}")
        };
        assert_eq!(
            target.energy_edge_sets(),
            vec![[2, 6].as_slice(), [3, 6].as_slice()]
        );
        for channel in catalogue.named_entries() {
            for block in &channel.blocks {
                if block.target.energy_edge_sets().len() != 2 {
                    continue;
                }
                context.insert_geometry_map(
                    block.target.clone(),
                    channel.definition.parent_lmb.clone(),
                    block.active_lmb.clone(),
                    block.preceding_lmb.clone(),
                    CompiledSamplingMap::Joint(map.clone()),
                )?;
                assert!(
                    context
                        .geometry_maps
                        .contains_key(&block.geometry_key(&channel.definition.parent_lmb))
                );
            }
        }
        assert_eq!(context.geometry_maps.len(), 2);
        let compiled = catalogue.compile(&context, &programs)?;
        assert!(compiled.iter().all(|channel| channel.dimensions() == 6));
        let restricted = compiled
            .iter()
            .filter(|channel| channel.map.contract().support == SamplingSupport::Restricted)
            .cloned()
            .collect();
        assert!(matches!(
            SamplingChannelBridge::new(restricted),
            Err(SamplingChannelBridgeError::MissingFullSupport)
        ));
        let bridge = SamplingChannelBridge::new(compiled)?;
        let selected = bridge
            .channels()
            .iter()
            .position(|channel| channel.name == "joint")
            .unwrap();
        let mapped = bridge.forward(
            SamplingChannelId(selected),
            &[0.3, 0.31, 0.37, 0.43, 0.31, 0.37],
        )?;
        assert!(
            mapped
                .map
                .diagnostics
                .iter()
                .any(|d| d.contains("certified normal disk"))
        );
        assert!(
            bridge
                .inverse(SamplingChannelId(selected), &mapped.raw_coordinates)?
                .is_some()
        );
        assert!((mapped.partition.weight_sum() - 1.0).abs() < 1.0e-14);
        Ok(())
    }

    #[test]
    fn ordinary_lmb_rejects_cut_metadata_and_meaningless_host_qualifiers() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["ordinary".into()],
            ..Default::default()
        };
        for around in [
            "at_cut(cut(3),lmb(1,2))",
            "at_cut(cut(3),then(lmb(1),lmb(2)))",
        ] {
            selection
                .channel_definitions
                .entry("G".into())
                .or_default()
                .insert("ordinary".into(), definition(around));
            assert!(
                resolve_sampling_channel_selection("G", &selection)
                    .unwrap_err()
                    .to_string()
                    .contains("at_cut requires a physical target")
            );
        }
        let mut ordinary = definition("lmb(1,2)");
        ordinary.on_cut = vec![1];
        selection
            .channel_definitions
            .get_mut("G")
            .unwrap()
            .insert("ordinary".into(), ordinary);
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            5.0,
            2,
        );
        let error = catalogue
            .compile(&context, &[(None, None, None)])
            .unwrap_err();
        assert!(matches!(
            error.downcast_ref::<SamplingChannelCompileError>(),
            Some(SamplingChannelCompileError::InvalidChannel { .. })
        ));
        assert!(
            error
                .to_string()
                .contains("on_cut requires an explicit physical host")
        );
    }

    #[test]
    fn native_lu_profile_binding_preserves_retryable_error_through_catalogue() {
        // Exact compiled constants isolate the numerical binding boundary from
        // physical root solving. Binary64 cannot represent this positive shape.
        std::thread::Builder::new()
            .stack_size(64 * 1024 * 1024)
            .spawn(|| {
                test_initialise().unwrap();
                let parameters = (0..5)
                    .map(|i| Atom::var(symbol!(format!("fit{i}"))))
                    .collect::<Vec<_>>();
                let program = SamplingExpressionEvaluator::new(
                    ["0", "1/10^400", "0", "0", "0"].map(|x| try_parse!(x).unwrap()),
                    parameters,
                    true,
                )
                .unwrap();
                let mut selection = SamplingChannelSelection {
                    default_channel_selection: vec!["cut".into()],
                    ..Default::default()
                };
                let mut channel = definition("phase_space(cut(3))");
                channel.parent_lmb = vec![1];
                channel.subspace_lmb = vec![1];
                channel.radial_profile = Some(SamplingRadialProfile::default());
                selection
                    .channel_definitions
                    .entry("G".into())
                    .or_default()
                    .insert("cut".into(), channel);
                let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
                let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
                let mut context = SamplingChannelCompileContext::<f64>::new(
                    "G",
                    vec![1],
                    ParameterizationSettings::default(),
                    5.0,
                    1,
                );
                context.physical_cut_ids.insert(vec![3], vec![0]);
                context.physical_cut_max_occurrences.insert(vec![3], 1);
                context
                    .insert_geometry_map(
                        SamplingMapDefinition::PhaseSpace(Box::new(SamplingMapDefinition::Cut(
                            vec![3],
                        ))),
                        vec![1],
                        vec![1],
                        vec![],
                        CompiledSamplingMap::ImplicitSurface(
                            ImplicitSurfaceRadialMap::new(
                                3,
                                vec![0.0; 3],
                                2.0,
                                1.0,
                                Arc::new(|_, r| Ok((r - 2.0, 1.0))),
                            )
                            .unwrap(),
                        ),
                    )
                    .unwrap();
                let error = catalogue
                    .compile(&context, &[(None, Some(program.clone()), None)])
                    .unwrap_err();
                assert!(matches!(
                    error.downcast_ref::<super::super::sampling_maps::SamplingEvaluationError>(),
                    Some(
                        super::super::sampling_maps::SamplingEvaluationError::Unrepresentable { .. }
                    )
                ));
                let zero = F::<crate::utils::ArbPrec>::default().zero();
                assert!(
                    ImplicitSurfaceRadialMap::new(
                        3,
                        vec![zero.0; 3],
                        2.0,
                        1.0,
                        Arc::new(|_, r: crate::utils::ArbPrec| {
                            let r = F(r);
                            Ok(((r.clone() - r.from_i64(2)).0, r.one().0))
                        })
                    )
                    .unwrap()
                    .with_lu_h_profile(program, &SamplingRadialProfile::default(), 1)
                    .is_ok()
                );
            })
            .unwrap()
            .join()
            .unwrap();
    }

    #[test]
    fn geometry_registry_keeps_distinct_prerequisite_layouts_and_foreign_inverses() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["first".into(), "reordered".into()],
            ..Default::default()
        };
        for (name, prior) in [("first", "2,4"), ("reordered", "4,2")] {
            let mut channel = definition(&format!("then(lmb({prior}),block(lmb(1),surface(3)))"));
            channel.parent_lmb = vec![1, 2, 4];
            channel.subspace_lmb.clear();
            selection
                .channel_definitions
                .entry("G".into())
                .or_default()
                .insert(name.into(), channel);
        }
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let mut context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1, 2, 4],
            ParameterizationSettings::default(),
            5.0,
            3,
        );
        for channel in catalogue.named_entries() {
            let block = &channel.blocks[1];
            let position = block
                .preceding_lmb
                .iter()
                .position(|edge| *edge == 2)
                .unwrap()
                * 3;
            let map = ImplicitSurfaceRadialMap::new(
                3,
                vec![0.0; 3],
                2.0,
                1.0,
                Arc::new(|_, r| Ok((r - 2.0, 1.0))),
            )
            .unwrap()
            .with_context_preparer(Arc::new(move |prior| {
                if prior.len() != 6 {
                    return Err(eyre!("wrong declared prerequisite dimension"));
                }
                Ok((
                    prior[position..position + 3]
                        .iter()
                        .map(|value| value * 0.25)
                        .collect(),
                    super::super::PreparedSurfaceStatus::existing(None)?,
                ))
            }));
            context
                .insert_geometry_map(
                    block.target.clone(),
                    channel.definition.parent_lmb.clone(),
                    block.active_lmb.clone(),
                    block.preceding_lmb.clone(),
                    CompiledSamplingMap::ImplicitSurface(map),
                )
                .unwrap();
        }
        assert_eq!(context.geometry_maps.len(), 2);
        let keys = context.geometry_maps.keys().cloned().collect::<Vec<_>>();
        assert_eq!(
            (&keys[0].0, &keys[0].1, &keys[0].2),
            (&keys[1].0, &keys[1].1, &keys[1].2)
        );
        assert_ne!(keys[0].3, keys[1].3);
        let programs = catalogue
            .compile_programs(9, &HFunctionSettings::default())
            .unwrap();
        let bridge =
            SamplingChannelBridge::new(catalogue.compile(&context, &programs).unwrap()).unwrap();
        let cube = [0.17, 0.29, 0.61, 0.23, 0.37, 0.73, 0.31, 0.43, 0.67];
        let forward = bridge.forward(SamplingChannelId(0), &cube).unwrap();
        assert!((forward.partition.weight_sum() - 1.0).abs() < 1.0e-13);
        for id in [0, 1] {
            let inverse = bridge
                .inverse(SamplingChannelId(id), &forward.raw_coordinates)
                .unwrap()
                .expect("full-support inverse");
            let again = bridge
                .forward(SamplingChannelId(id), &inverse.map.coordinates)
                .unwrap();
            assert!((inverse.map.inverse_jacobian * again.map.jacobian - 1.0).abs() < 1.0e-10);
            for (actual, expected) in again.raw_coordinates.iter().zip(&forward.raw_coordinates) {
                assert!((actual - expected).abs() < 1.0e-9);
            }
        }
        let map = context.geometry_maps.values().next().unwrap().clone();
        assert!(
            context
                .insert_geometry_map(
                    SamplingMapDefinition::Surface(vec![3]),
                    vec![1, 2, 4],
                    vec![1],
                    vec![2, 2],
                    map
                )
                .is_err()
        );
        assert_eq!(
            context.geometry_maps.len(),
            2,
            "invalid registration is transactional"
        );
    }

    fn proxy_selection(first: Option<&str>, second: Option<&str>) -> SamplingChannelSelection {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["first".into(), "second".into()];
        let definitions = selection.channel_definitions.entry("G".into()).or_default();
        for (name, proxy) in [("first", first), ("second", second)] {
            let channel = SamplingChannelDefinition {
                radial_profile: None,
                around: "lmb(1)".into(),
                subspace_lmb: Vec::new(),
                parent_lmb: vec![1],
                on_cut: Vec::new(),
                singularity_proxy: proxy.map(str::to_owned),
            };
            definitions.insert(name.into(), channel);
        }
        selection
    }

    fn proxy_bridge(selection: &SamplingChannelSelection) -> Result<SamplingChannelBridge> {
        test_initialise()?;
        let resolved = resolve_sampling_channel_selection("G", selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1],
            ParameterizationSettings::default(),
            100.0,
            1,
        );
        let proxies = catalogue
            .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())?;
        let channels = catalogue.compile(&context, &proxies)?;
        SamplingChannelBridge::new_with_partition_mode(
            channels,
            SamplingPartitionMode::SingularityProxy,
        )
        .map_err(Into::into)
    }

    #[test]
    fn named_proxy_is_eagerly_compiled_and_uses_complete_frame_coordinates() {
        let selection = proxy_selection(Some("1 + x0^2 + x1^2 + x2^2"), Some("2 + x0^2"));
        let bridge = proxy_bridge(&selection).unwrap();
        assert_eq!(bridge.dimensions(), 3);
        let partition = bridge.partition(&[0.2, 0.3, 0.4]).unwrap();
        let first = 1.0 + 0.2_f64.powi(2) + 0.3_f64.powi(2) + 0.4_f64.powi(2);
        let second = 2.0 + 0.2_f64.powi(2);
        assert!((partition.weights[0] - first / (first + second)).abs() < 1.0e-13);
        assert!((partition.weights[1] - second / (first + second)).abs() < 1.0e-13);
    }

    #[test]
    fn singularity_proxy_mode_rejects_missing_channel_metadata() {
        let selection = proxy_selection(Some("1 + x0^2"), None);
        let error = proxy_bridge(&selection)
            .unwrap()
            .partition(&[0.2, 0.3, 0.4]);
        let error = error.unwrap_err();
        assert!(
            error
                .to_string()
                .contains("has no singularity_proxy metadata")
        );
    }

    #[test]
    fn singularity_proxy_reports_symbolica_parse_and_positivity_errors() {
        let invalid = proxy_selection(Some("x0 + ("), Some("1"));
        let error = resolve_sampling_channel_selection("G", &invalid).unwrap_err();
        assert!(
            error
                .to_string()
                .contains("failed to parse singularity_proxy")
        );

        let non_positive = proxy_selection(Some("x0 - x0"), Some("1"));
        let error = proxy_bridge(&non_positive).unwrap_err();
        assert!(error.to_string().contains("strictly positive"));
    }

    #[test]
    fn graph_entry_replaces_default_and_deduplicates_entries() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["auto:lmb".into(), "common".into()],
            ..Default::default()
        };
        selection.channel_selection.insert(
            "G".into(),
            vec!["auto:surfaces".into(), "named".into(), "named".into()],
        );
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("named".into(), definition("surface(1,2)"));

        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        assert_eq!(resolved.graph_name, "G");
        assert_eq!(
            resolved.selectors,
            vec![
                SamplingChannelSelector::Preset(SamplingChannelPreset::Surfaces),
                SamplingChannelSelector::Named("named".into())
            ]
        );
        assert!(resolved.named("named").is_some());
    }

    #[test]
    fn additive_resolver_keeps_default_and_graph_entries() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["auto:lmb".into(), "common".into()],
            ..Default::default()
        };
        selection.channel_selection.insert(
            "G".into(),
            vec!["auto:surfaces".into(), "common".into(), "named".into()],
        );
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("common".into(), definition("lmb(1,2)"));
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("named".into(), definition("surface(1,2)"));

        let resolved = resolve_sampling_channel_selection_with_defaults("G", &selection).unwrap();
        assert_eq!(
            resolved.selectors,
            vec![
                SamplingChannelSelector::Preset(SamplingChannelPreset::Lmb),
                SamplingChannelSelector::Named("common".into()),
                SamplingChannelSelector::Preset(SamplingChannelPreset::Surfaces),
                SamplingChannelSelector::Named("named".into())
            ]
        );
    }

    #[test]
    fn replacement_variant_ignores_defaults() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["auto:lmb".into()],
            ..Default::default()
        };
        selection
            .channel_selection
            .insert("G".into(), vec!["auto:surfaces".into()]);
        let resolved =
            resolve_sampling_channel_selection_replacing_default("G", &selection).unwrap();
        assert_eq!(
            resolved.selectors,
            vec![SamplingChannelSelector::Preset(
                SamplingChannelPreset::Surfaces
            )]
        );
    }

    #[test]
    fn implicit_default_preserves_empty_and_named_graph_selections() {
        let mut selection = SamplingChannelSelection::default();
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        assert_eq!(
            resolved.selectors,
            vec![SamplingChannelSelector::Preset(
                SamplingChannelPreset::OptimizedLmb
            )]
        );
        selection.channel_selection.insert("G".into(), vec![]);
        assert!(
            resolve_sampling_channel_selection("G", &selection)
                .unwrap()
                .selectors
                .is_empty()
        );
        selection
            .channel_selection
            .insert("G".into(), vec!["custom".into()]);
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("custom".into(), definition("lmb(1,2)"));
        assert_eq!(
            resolve_sampling_channel_selection("G", &selection)
                .unwrap()
                .selectors,
            vec![SamplingChannelSelector::Named("custom".into())]
        );
    }

    #[test]
    fn default_is_used_when_graph_has_no_override() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["auto:optimized_lmb".into(), "named".into()],
            ..Default::default()
        };
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("named".into(), definition("surface(1)"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        assert!(resolved.has_preset(SamplingChannelPreset::OptimizedLmb));
        assert_eq!(resolved.named_channels.len(), 1);
    }

    #[test]
    fn missing_named_definition_reports_graph_and_available_names() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["missing".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("available".into(), definition("surface(1)"));
        let error = resolve_sampling_channel_selection("G", &selection).unwrap_err();
        assert_eq!(
            error,
            SamplingSelectionError::MissingChannelDefinition {
                graph: "G".into(),
                channel: "missing".into(),
                available: vec!["available".into()]
            }
        );
    }

    #[test]
    fn named_channel_requires_explicit_parent_lmb() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["named".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert(
                "named".into(),
                SamplingChannelDefinition {
                    radial_profile: None,
                    around: "lmb(1,2)".into(),
                    subspace_lmb: Vec::new(),
                    parent_lmb: Vec::new(),
                    on_cut: Vec::new(),
                    singularity_proxy: None,
                },
            );
        assert!(matches!(
            resolve_sampling_channel_selection("G", &selection),
            Err(SamplingSelectionError::MissingParentLmb { .. })
        ));
    }

    #[test]
    fn named_surface_requires_explicit_active_subspace() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["threshold".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert(
                "threshold".into(),
                SamplingChannelDefinition {
                    radial_profile: None,
                    around: "surface(2,4)".into(),
                    subspace_lmb: Vec::new(),
                    parent_lmb: vec![1, 2, 4],
                    on_cut: Vec::new(),
                    singularity_proxy: None,
                },
            );
        assert!(matches!(
            resolve_sampling_channel_selection("G", &selection),
            Err(SamplingSelectionError::MissingSubspaceLmb { .. })
        ));
    }

    #[test]
    fn unknown_auto_preset_is_rejected() {
        assert_eq!(
            SamplingChannelSelector::parse("auto:not_a_mode").unwrap_err(),
            SamplingSelectionError::UnknownPreset("auto:not_a_mode".into())
        );
    }

    #[test]
    fn catalogue_expands_presets_and_keeps_named_diagnostics() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["auto:optimized_lmb".into(), "surface_hz".into()],
            ..Default::default()
        };
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("surface_hz".into(), definition("surface(2,4)"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue =
            build_sampling_channel_catalogue(&resolved, &[(0, vec![1, 2]), (1, vec![2, 4])], &[1]);
        assert_eq!(catalogue.entries.len(), 2);
        assert_eq!(catalogue.lmb_basis_entries().next(), Some((1, &[2, 4][..])));
        assert_eq!(catalogue.named_entries().next().unwrap().name, "surface_hz");
        assert!(catalogue.inspection_rows()[0].contains("basis=1"));
        assert!(catalogue.inspection_rows()[1].contains("parent_lmb=[1, 2]"));
    }

    #[test]
    fn catalogue_keeps_one_id_axis_when_named_channel_precedes_lmb() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["threshold".into(), "auto:lmb".into()],
            ..Default::default()
        };
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("threshold".into(), definition("surface(1,2)"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[(7, vec![1, 2])], &[7]);

        // The LMB basis number is metadata of the second catalogue entry. It
        // must not become a compacted positional channel number after the
        // named entry, otherwise LMB-only consumers would enumerate a second
        // channel axis and disagree with the canonical SamplingChannelId.
        assert!(matches!(
            catalogue.entries.first(),
            Some(SamplingCatalogueEntry::Named(channel)) if channel.name == "threshold"
        ));
        assert!(matches!(
            catalogue.entries.get(1),
            Some(SamplingCatalogueEntry::Lmb { basis_id: 7, .. })
        ));
        assert_eq!(
            catalogue
                .entries
                .iter()
                .enumerate()
                .map(|(index, _)| SamplingChannelId::from(index))
                .collect::<Vec<_>>(),
            vec![SamplingChannelId::from(0), SamplingChannelId::from(1)]
        );
    }

    #[test]
    fn inspection_reports_resolved_selectors_and_canonical_rows() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["auto:optimized_lmb".into()];
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[(0, vec![1, 2])], &[0]);

        let report = catalogue.inspection();
        assert_eq!(report.graph_name, "G");
        assert_eq!(report.selectors, vec!["auto:optimized_lmb"]);
        assert_eq!(report.entries, catalogue.inspection_rows());
        assert_eq!(
            report.entries,
            vec!["0: lmb basis=0 edges=[1, 2] source=auto:optimized_lmb"]
        );
    }

    #[test]
    fn surface_preset_has_optimized_lmb_fallback() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["auto:surfaces".into(), "auto:lmb".into()];
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue =
            build_sampling_channel_catalogue(&resolved, &[(0, vec![1]), (1, vec![2])], &[0, 1]);
        assert_eq!(catalogue.lmb_basis_entries().count(), 2);
    }

    #[test]
    fn surface_preset_appends_master_frame_surface_candidates() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["auto:surfaces".into()];
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue_with_surfaces(
            &resolved,
            &[(0, vec![1, 2])],
            &[0],
            &[vec![1, 2], vec![2, 3]],
            &[1, 2],
        );
        assert!(catalogue.entries.iter().any(|entry| matches!(
            entry,
            SamplingCatalogueEntry::Surface { edges, parent_lmb }
                if edges == &vec![1, 2] && parent_lmb == &vec![1, 2]
        )));
    }

    #[test]
    fn optimized_and_surface_presets_retain_full_domain_and_elementary_soft_coverage() {
        let all_lmbs = vec![(0, vec![4, 5]), (1, vec![1, 5]), (2, vec![2, 4])];
        for preset in ["auto:optimized_lmb", "auto:surfaces"] {
            let mut selection = SamplingChannelSelection::default();
            selection.default_channel_selection = vec![preset.into()];
            let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
            // The heuristic selected only basis 0.  The automatic audit adds
            // the smallest available basis exposing each omitted
            // loop-dependent soft edge, while leaving the catalogue's single
            // canonical axis intact.
            let catalogue = build_sampling_channel_catalogue_with_surfaces_and_coverage(
                &resolved,
                &all_lmbs,
                &[0],
                &[],
                &[4, 5],
                &[1, 2, 99],
            );
            assert_eq!(
                catalogue
                    .lmb_basis_entries()
                    .map(|(basis, _)| basis)
                    .collect::<Vec<_>>(),
                vec![0, 1, 2]
            );
            let report = catalogue.coverage_report(&all_lmbs, &[1, 2, 99]);
            assert!(report.has_full_domain_lmb);
            assert_eq!(report.covered_massless_loop_edges, vec![1, 2]);
            assert!(report.missing_massless_loop_edges.is_empty());
        }
    }

    #[test]
    fn optimized_preset_falls_back_to_one_full_domain_lmb_when_heuristic_is_empty() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["auto:optimized_lmb".into()];
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue_with_surfaces_and_coverage(
            &resolved,
            &[(7, vec![1, 2])],
            &[],
            &[],
            &[1, 2],
            &[],
        );
        assert_eq!(
            catalogue.lmb_basis_entries().collect::<Vec<_>>(),
            vec![(7, &[1, 2][..])]
        );
        assert!(
            catalogue
                .coverage_report(&[(7, vec![1, 2])], &[])
                .has_full_domain_lmb
        );
    }

    #[test]
    fn catalogue_compiles_lmb_and_surface_with_master_embedding() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["auto:lmb".into(), "threshold".into()],
            ..Default::default()
        };
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("threshold".into(), definition("surface(1,2)"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[(0, vec![1, 2])], &[0]);
        let mut context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context
            .insert_geometry_map(
                SamplingMapDefinition::Surface(vec![1, 2]),
                context.parent_lmb.clone(),
                vec![1, 2],
                vec![],
                CompiledSamplingMap::Surface(
                    SurfaceRadialMap::new(
                        3 * (vec![1, 2]).len(),
                        vec![0.0; 6],
                        Some(3.0),
                        2.0,
                        1.0,
                    )
                    .unwrap(),
                ),
            )
            .unwrap();
        let compiled = catalogue
            .compile(
                &context,
                &catalogue
                    .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                    .unwrap(),
            )
            .unwrap();
        assert_eq!(compiled.len(), 2);
        assert!(matches!(compiled[0].map, CompiledSamplingMap::Lmb(_)));
        assert_eq!(compiled[0].embedded_edges, vec![1, 2]);
        assert!(matches!(compiled[1].map, CompiledSamplingMap::Surface(_)));
        assert_eq!(compiled[1].master_graph, "G");
        assert_eq!(compiled[1].dimensions(), 6);
    }

    #[test]
    fn singular_primitives_fail_with_actionable_exact_density_diagnostics() {
        for (name, around, expected_primitive, expected_requirement) in [
            (
                "soft_target",
                "soft(1)",
                "soft(1)",
                "routed three-momentum frame",
            ),
            (
                "collinear_target",
                "collinear(1,2)",
                "collinear(1,2)",
                "relative angular frame",
            ),
            (
                "nested_soft_target",
                "then(lmb(2), block(lmb(1),soft(1)))",
                "soft(1)",
                "normalized radial profile",
            ),
        ] {
            let mut selection = SamplingChannelSelection::default();
            selection.default_channel_selection = vec![name.to_owned()];
            selection
                .channel_definitions
                .entry("G".to_owned())
                .or_default()
                .insert(name.to_owned(), definition(around));
            let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
            let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
            let context = SamplingChannelCompileContext::<f64>::new(
                "G",
                vec![1, 2],
                ParameterizationSettings::default(),
                100.0,
                2,
            );
            let error = catalogue
                .compile(
                    &context,
                    &catalogue
                        .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                        .unwrap(),
                )
                .unwrap_err();
            match error.downcast_ref::<SamplingChannelCompileError>() {
                Some(SamplingChannelCompileError::UnsupportedSingularPrimitive {
                    primitive,
                    reason,
                    ..
                }) => {
                    assert_eq!(primitive, expected_primitive);
                    assert!(reason.contains(expected_requirement), "{reason}");
                    assert!(reason.contains("edge ids alone"), "{reason}");
                }
                other => panic!("expected singular primitive diagnostic, got {other:?}"),
            }
        }
    }

    #[test]
    fn catalogue_accepts_prepared_exact_implicit_surface_map() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["threshold".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("threshold".into(), definition("surface(1,2)"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let mut context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context
            .insert_geometry_map(
                SamplingMapDefinition::Surface(vec![1, 2]),
                context.parent_lmb.clone(),
                vec![1, 2],
                vec![],
                CompiledSamplingMap::ImplicitSurface(
                    ImplicitSurfaceRadialMap::new(
                        6,
                        vec![0.0; 6],
                        2.0,
                        2.0,
                        Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
                    )
                    .unwrap(),
                ),
            )
            .unwrap();
        let compiled = catalogue
            .compile(
                &context,
                &catalogue
                    .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                    .unwrap(),
            )
            .unwrap();
        assert!(matches!(
            compiled[0].map,
            CompiledSamplingMap::ImplicitSurface(_)
        ));
        let bridge = SamplingChannelBridge::new(compiled).unwrap();
        let mapped = bridge.forward(
            SamplingChannelId::from(0),
            &[0.31, 0.42, 0.57, 0.23, 0.68, 0.81],
        );
        assert!(mapped.is_ok());
    }

    #[test]
    fn catalogue_compiles_then_with_context_surface_and_bridges_full_composition() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["conditional".into()];
        let mut channel = definition("then(lmb(1),surface(2,4))");
        channel.subspace_lmb = vec![2];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("conditional".into(), channel);
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let mut context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context
            .insert_geometry_map(
                SamplingMapDefinition::Surface(vec![2, 4]),
                context.parent_lmb.clone(),
                vec![2],
                vec![1],
                CompiledSamplingMap::ImplicitSurface(
                    ImplicitSurfaceRadialMap::new(
                        3,
                        vec![0.0; 3],
                        1.0,
                        1.0,
                        Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
                    )
                    .unwrap()
                    .with_context_evaluator(Arc::new(
                        |_, radius, _, context| {
                            let shift = context.first().copied().unwrap_or_default().abs().min(0.2);
                            Ok((radius - (1.0 + shift), 1.0))
                        },
                    )),
                ),
            )
            .unwrap();
        let compiled = catalogue
            .compile(
                &context,
                &catalogue
                    .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                    .unwrap(),
            )
            .unwrap();
        assert_eq!(compiled.len(), 1);
        let channel = &compiled[0];
        assert_eq!(channel.embedded_edges, vec![1, 2]);
        assert_eq!(channel.map.contract().support, SamplingSupport::Full);
        let coordinates = [0.31, 0.42, 0.57, 0.23, 0.68, 0.81];
        let mapped = channel.map.forward(&coordinates).unwrap();
        assert_eq!(mapped.support, SamplingSupport::Full);
        let inverse = channel
            .map
            .inverse(&mapped.point)
            .unwrap()
            .expect("full-support inverse");
        assert!(inverse.residual < 1.0e-9, "{}", inverse.residual);
        assert!(
            inverse
                .coordinates
                .iter()
                .zip(coordinates)
                .all(|(actual, expected)| (actual - expected).abs() < 1.0e-9)
        );
        let bridge = SamplingChannelBridge::new(compiled).unwrap();
        let bridged = bridge
            .forward(SamplingChannelId::from(0), &coordinates)
            .unwrap();
        assert!((bridged.partition.weight_sum() - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn runtime_context_partition_uses_each_canonical_channel_context() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["left_chart".into(), "right_chart".into()];
        for name in ["left_chart", "right_chart"] {
            let mut channel = definition("then(lmb(1),surface(2,4))");
            channel.subspace_lmb = vec![2];
            selection
                .channel_definitions
                .entry("G".into())
                .or_default()
                .insert(name.into(), channel);
        }
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let mut context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context
            .insert_geometry_map(
                SamplingMapDefinition::Surface(vec![2, 4]),
                context.parent_lmb.clone(),
                vec![2],
                vec![1],
                CompiledSamplingMap::ImplicitSurface(
                    ImplicitSurfaceRadialMap::new(
                        3,
                        vec![0.0; 3],
                        1.0,
                        1.0,
                        Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
                    )
                    .unwrap()
                    .with_context_evaluator(Arc::new(
                        |_, radius, _, context| {
                            let shift = context.first().copied().unwrap_or_default();
                            Ok((radius - (1.0 + shift), 1.0))
                        },
                    )),
                ),
            )
            .unwrap();
        let bridge = SamplingChannelBridge::new(
            catalogue
                .compile(
                    &context,
                    &catalogue
                        .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                        .unwrap(),
                )
                .unwrap(),
        )
        .unwrap();
        let mut contexts = SamplingChannelRuntimeContexts::new(bridge.channels().len());
        contexts.set(SamplingChannelId::from(0), vec![0.0]).unwrap();
        contexts.set(SamplingChannelId::from(1), vec![0.2]).unwrap();
        let coordinates = [0.31, 0.42, 0.57, 0.23, 0.68, 0.81];
        let mapped = bridge
            .forward_with_runtime_contexts(SamplingChannelId::from(0), &coordinates, &mut contexts)
            .unwrap();
        assert!((mapped.partition.weight_sum() - 1.0).abs() < 1.0e-12);
        assert!(mapped.partition.weights.iter().all(|weight| *weight > 0.0));
        let inverse = bridge
            .inverse_with_runtime_contexts(
                SamplingChannelId::from(0),
                &mapped.raw_coordinates,
                &mut contexts,
            )
            .unwrap()
            .expect("full-support inverse");
        assert!(inverse.map.residual < 1.0e-9, "{}", inverse.map.residual);
        assert!(
            inverse
                .map
                .coordinates
                .iter()
                .zip(coordinates)
                .all(|(actual, expected)| (actual - expected).abs() < 1.0e-9)
        );
    }

    #[test]
    fn distinct_physical_surfaces_share_one_active_subspace_without_overwrite() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["h".into(), "z".into()];
        let mut h = definition("surface(2,4)");
        h.subspace_lmb = vec![1, 2];
        let mut z = definition("surface(3,10)");
        z.subspace_lmb = vec![1, 2];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .extend([(String::from("h"), h), (String::from("z"), z)]);
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let mut context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        for edges in [vec![2, 4], vec![3, 10]] {
            context
                .insert_geometry_map(
                    SamplingMapDefinition::Surface(edges),
                    context.parent_lmb.clone(),
                    vec![1, 2],
                    vec![],
                    CompiledSamplingMap::ImplicitSurface(
                        ImplicitSurfaceRadialMap::new(
                            6,
                            vec![0.0; 6],
                            2.0,
                            2.0,
                            Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
                        )
                        .unwrap(),
                    ),
                )
                .unwrap();
        }
        let compiled = catalogue
            .compile(
                &context,
                &catalogue
                    .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                    .unwrap(),
            )
            .unwrap();
        assert_eq!(compiled.len(), 2);
        assert_eq!(compiled[0].name, "h");
        assert_eq!(compiled[1].name, "z");
        assert!(
            compiled
                .iter()
                .all(|channel| channel.embedded_edges == [1, 2])
        );
        assert_eq!(
            compiled[0].definition,
            SamplingMapDefinition::Surface(vec![2, 4])
        );
        assert_eq!(
            compiled[1].definition,
            SamplingMapDefinition::Surface(vec![3, 10])
        );
        let bridge = SamplingChannelBridge::new(compiled).unwrap();
        let point = bridge
            .forward(SamplingChannelId(0), &[0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
            .unwrap();
        assert!((point.partition.weight_sum() - 1.0).abs() < 1.0e-13);
    }

    #[test]
    fn catalogue_rejects_non_parent_lmb_until_affine_routing_is_compiled() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["auto:lmb".into()];
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue =
            build_sampling_channel_catalogue(&resolved, &[(0, vec![1, 2]), (1, vec![2, 4])], &[0]);
        let context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        let error = catalogue
            .compile(
                &context,
                &catalogue
                    .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                    .unwrap(),
            )
            .unwrap_err();
        assert!(matches!(
            error.downcast_ref::<SamplingChannelCompileError>(),
            Some(SamplingChannelCompileError::MissingLmbFrameMap { .. })
        ));
        assert!(error.to_string().contains("no affine frame map"));
    }

    #[test]
    fn catalogue_wraps_non_parent_lmb_in_supplied_affine_frame() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["auto:lmb".into()];
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[(0, vec![2, 4])], &[0]);
        let mut context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context.lmb_frame_maps.insert(
            0,
            SamplingMapAffine::new(
                vec![
                    vec![2.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                    vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                    vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                    vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                    vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                ],
                vec![0.5; 6],
            )
            .unwrap(),
        );
        let compiled = catalogue
            .compile(
                &context,
                &catalogue
                    .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                    .unwrap(),
            )
            .unwrap();
        assert!(matches!(
            compiled[0].map,
            CompiledSamplingMap::Affine { .. }
        ));
        assert_eq!(compiled[0].embedded_edges, vec![1, 2]);
        let point = compiled[0]
            .forward(&[0.31, 0.42, 0.57, 0.23, 0.68, 0.81])
            .unwrap();
        assert!(point.residual < 1.0e-10);
        assert!((point.point[0] - 0.5).abs() > 1.0e-6);
        let inverse = compiled[0]
            .inverse(&point.point)
            .unwrap()
            .expect("full-support inverse");
        assert!(inverse.residual < 1.0e-10);
        for (recovered, original) in inverse.coordinates.iter().zip(&point.coordinates) {
            assert!((recovered - original).abs() < 1.0e-10);
        }
        assert!((point.jacobian * inverse.inverse_jacobian - 1.0).abs() < 1.0e-10);
        let bridge = SamplingChannelBridge::new(compiled.clone()).unwrap();
        let bridged = bridge
            .forward(
                SamplingChannelId::from(0),
                &[0.31, 0.42, 0.57, 0.23, 0.68, 0.81],
            )
            .unwrap();
        assert_eq!(bridged.raw_coordinates, point.point);
        assert!((bridged.partition.weight_sum() - 1.0).abs() < 1.0e-12);

        // The same catalogue binds a distinct native numerical frame without
        // re-enumerating channels or rounding the supplied routing to f64.
        let one = F::<crate::utils::f128>::from_f64(1.0);
        let displacement = one.from_i64(2).powi(-80);
        assert_eq!((one + displacement).into_f64(), 1.0);
        let identity = (0..6)
            .map(|row| {
                (0..6)
                    .map(|column| if row == column { one.0 } else { one.zero().0 })
                    .collect()
            })
            .collect::<Vec<Vec<_>>>();
        let mut native_context = SamplingChannelCompileContext::<crate::utils::QuadFloat>::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        native_context.lmb_frame_maps.insert(
            0,
            SamplingMapAffine::new(identity.clone(), vec![one.zero().0; 6]).unwrap(),
        );
        let coordinates = [0.31, 0.42, 0.57, 0.23, 0.68, 0.81]
            .map(|value| F::<crate::utils::f128>::from_f64(value).0);
        let native = catalogue
            .compile(
                &native_context,
                &catalogue
                    .compile_programs(
                        3 * native_context.n_loop_momenta,
                        &HFunctionSettings::default(),
                    )
                    .unwrap(),
            )
            .unwrap();
        let original = native[0].forward(&coordinates).unwrap();
        let mut matrix = identity;
        matrix[0][0] = (one + displacement).0;
        native_context.lmb_frame_maps.insert(
            0,
            SamplingMapAffine::new(matrix, vec![displacement.0; 6]).unwrap(),
        );
        let shifted = catalogue
            .compile(
                &native_context,
                &catalogue
                    .compile_programs(
                        3 * native_context.n_loop_momenta,
                        &HFunctionSettings::default(),
                    )
                    .unwrap(),
            )
            .unwrap();
        assert_eq!(shifted[0].name, native[0].name);
        assert_eq!(shifted[0].embedded_edges, native[0].embedded_edges);
        let mapped = shifted[0].forward(&coordinates).unwrap();
        assert_ne!(mapped.point[0], original.point[0]);
        assert_eq!(
            F(mapped.point[0]),
            F(original.point[0]) * (one + displacement) + displacement
        );
        assert_eq!(
            F(mapped.jacobian),
            F(original.jacobian) * (one + displacement)
        );
        let inverse = shifted[0]
            .inverse(&mapped.point)
            .unwrap()
            .expect("full-support inverse");
        for (actual, expected) in inverse.coordinates.iter().zip(coordinates) {
            assert!((F(*actual) - F(expected)).abs() < one.from_i64(10).powi(-25));
        }
    }

    #[test]
    fn named_non_parent_lmb_uses_edge_keyed_affine_frame() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["named".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert(
                "named".into(),
                SamplingChannelDefinition {
                    radial_profile: None,
                    around: "lmb(2,4)".into(),
                    subspace_lmb: Vec::new(),
                    parent_lmb: vec![1, 2],
                    on_cut: Vec::new(),
                    singularity_proxy: None,
                },
            );
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let mut context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context.lmb_frame_maps_by_edges.insert(
            vec![2, 4],
            SamplingMapAffine::new(
                vec![
                    vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                    vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                    vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                    vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                ],
                vec![0.0; 6],
            )
            .unwrap(),
        );
        let compiled = catalogue
            .compile(
                &context,
                &catalogue
                    .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                    .unwrap(),
            )
            .unwrap();
        assert!(matches!(
            compiled[0].map,
            CompiledSamplingMap::Affine { .. }
        ));
        assert_eq!(compiled[0].embedded_edges, vec![1, 2]);
    }

    #[test]
    fn catalogue_embeds_partial_surface_with_complement_lmb() {
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["threshold".into()],
            ..Default::default()
        };
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("threshold".into(), {
                let mut channel = definition("surface(2)");
                channel.subspace_lmb = vec![2];
                channel
            });
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let mut context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context
            .insert_geometry_map(
                SamplingMapDefinition::Surface(vec![2]),
                context.parent_lmb.clone(),
                vec![2],
                vec![1],
                CompiledSamplingMap::ImplicitSurface(
                    ImplicitSurfaceRadialMap::new(
                        3,
                        vec![0.0; 3],
                        2.0,
                        1.0,
                        Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
                    )
                    .unwrap(),
                ),
            )
            .unwrap();
        let compiled = catalogue
            .compile(
                &context,
                &catalogue
                    .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                    .unwrap(),
            )
            .unwrap();
        assert_eq!(compiled.len(), 1);
        assert!(matches!(compiled[0].map, CompiledSamplingMap::Embedded(_)));
        assert_eq!(compiled[0].embedded_edges, vec![1, 2]);
        let point = compiled[0]
            .forward(&[0.31, 0.42, 0.57, 0.23, 0.68, 0.81])
            .unwrap();
        let inverse = compiled[0]
            .inverse(&point.point)
            .unwrap()
            .expect("full-support inverse");
        assert!(inverse.residual < 1.0e-10);
    }

    #[test]
    fn phase_space_channel_requires_registered_physical_host() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["cut_chart".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("cut_chart".into(), definition("phase_space(cut(4,7))"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        let error = catalogue
            .compile(
                &context,
                &catalogue
                    .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                    .unwrap(),
            )
            .unwrap_err();
        assert!(error.to_string().contains("graph-resolved physical host"));
    }

    #[test]
    fn standalone_cut_charts_share_the_canonical_partition_without_side_context() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["first".into(), "second".into()];
        let definitions = selection.channel_definitions.entry("G".into()).or_default();
        for (name, around, cut) in [
            ("first", "phase_space(cut(4,7))", 3),
            ("second", "phase_space(cut(5,8))", 6),
        ] {
            let mut channel = definition(around);
            channel.parent_lmb = vec![1];
            channel.subspace_lmb = vec![1];
            channel.on_cut = vec![cut];
            definitions.insert(name.into(), channel);
        }
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let mut context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1],
            ParameterizationSettings::default(),
            2.0,
            1,
        );
        for (edges, cut, radius) in [(vec![4, 7], 3, 2.0), (vec![5, 8], 6, 3.0)] {
            context.physical_cut_ids.insert(edges.clone(), vec![cut]);
            context
                .insert_geometry_map(
                    SamplingMapDefinition::PhaseSpace(Box::new(SamplingMapDefinition::Cut(edges))),
                    context.parent_lmb.clone(),
                    vec![1],
                    vec![],
                    CompiledSamplingMap::ImplicitSurface(
                        ImplicitSurfaceRadialMap::new(
                            3,
                            vec![0.0; 3],
                            2.0,
                            1.0,
                            std::sync::Arc::new(move |_: &[f64], r: f64| {
                                Ok((r * r - radius * radius, 2.0 * r))
                            }),
                        )
                        .unwrap(),
                    ),
                )
                .unwrap();
        }
        let bridge = SamplingChannelBridge::new(
            catalogue
                .compile(
                    &context,
                    &catalogue
                        .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                        .unwrap(),
                )
                .unwrap(),
        )
        .unwrap();
        assert_eq!(bridge.channels().len(), 2);
        for (id, channel) in bridge.channels().iter().enumerate() {
            assert_eq!(channel.embedded_edges, vec![1]);

            let mapped = bridge
                .forward(SamplingChannelId::from(id), &[0.31, 0.47, 0.68])
                .unwrap();
            assert!((mapped.partition.weight_sum() - 1.0).abs() < 1.0e-13);
            assert!(
                bridge
                    .inverse(SamplingChannelId::from(id), &mapped.raw_coordinates)
                    .unwrap()
                    .expect("full-support inverse")
                    .map
                    .residual
                    < 1.0e-9
            );
        }
        let report = SamplingChannelBridgeAcceptanceReport::normalized_gaussian(
            &bridge, 4096, 1.0, &[0.0; 3],
        )
        .unwrap();
        assert!((report.normalization - 1.0).abs() < 0.015, "{report:?}");
        assert!((report.second_moment - 3.0).abs() < 0.06, "{report:?}");
        // Neither a cut label without a map nor a map without a graph cut
        // registration is sufficient to enable the physical-cut syntax.
        context.physical_cut_ids.remove(&vec![4, 7]);
        assert!(
            catalogue
                .compile(
                    &context,
                    &catalogue
                        .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                        .unwrap()
                )
                .is_err()
        );
        context.physical_cut_ids.insert(vec![4, 7], vec![9]);
        assert!(
            catalogue
                .compile(
                    &context,
                    &catalogue
                        .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                        .unwrap()
                )
                .unwrap_err()
                .to_string()
                .contains("on_cut [3]")
        );
        // Hosted children are now compiled from declared prerequisites; no
        // external complete-sample handoff is required after the map is bound.
    }

    #[test]
    fn surface_compilation_requires_prepared_geometry() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["threshold".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("threshold".into(), definition("surface(1,2)"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        assert!(matches!(
            catalogue
                .compile(
                    &context,
                    &catalogue
                        .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                        .unwrap()
                )
                .unwrap_err()
                .downcast_ref::<SamplingChannelCompileError>(),
            Some(SamplingChannelCompileError::MissingGeometry { .. })
        ));
    }

    #[test]
    fn composite_map_with_overlapping_blocks_is_rejected_with_context_diagnostic() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["joint".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("joint".into(), definition("product(surface(1), lmb(1,2))"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        let mut context = context;
        context
            .insert_geometry_map(
                SamplingMapDefinition::Surface(vec![1]),
                context.parent_lmb.clone(),
                vec![1, 2],
                vec![],
                CompiledSamplingMap::Surface(
                    SurfaceRadialMap::new(
                        3 * (vec![1, 2]).len(),
                        vec![0.0; 6],
                        Some(3.0),
                        2.0,
                        1.0,
                    )
                    .unwrap(),
                ),
            )
            .unwrap();
        let error = catalogue
            .compile(
                &context,
                &catalogue
                    .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                    .unwrap(),
            )
            .unwrap_err()
            .to_string();
        assert!(error.contains("overlapping child blocks on edge 1"));
        assert!(error.contains("parent LMB [1, 2]"));
    }

    #[test]
    fn product_surface_and_complement_compile_as_one_master_frame_channel() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["joint".into()];
        let mut channel = definition("product(surface(2), complement(1))");
        channel.subspace_lmb = vec![2];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("joint".into(), channel);
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let mut context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context
            .insert_geometry_map(
                SamplingMapDefinition::Surface(vec![2]),
                context.parent_lmb.clone(),
                vec![2],
                vec![],
                CompiledSamplingMap::ImplicitSurface(
                    ImplicitSurfaceRadialMap::new(
                        3,
                        vec![0.0; 3],
                        2.0,
                        1.0,
                        Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
                    )
                    .unwrap(),
                ),
            )
            .unwrap();
        let compiled = catalogue
            .compile(
                &context,
                &catalogue
                    .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                    .unwrap(),
            )
            .unwrap();
        assert!(matches!(compiled[0].map, CompiledSamplingMap::Embedded(_)));
        assert_eq!(compiled[0].embedded_edges, vec![1, 2]);
        assert_eq!(
            compiled[0].map.contract().support,
            crate::integrands::process::SamplingSupport::Full
        );
        let point = compiled[0]
            .forward(&[0.31, 0.42, 0.57, 0.23, 0.68, 0.81])
            .unwrap();
        let inverse = compiled[0]
            .inverse(&point.point)
            .unwrap()
            .expect("full-support inverse");
        assert!(point.residual < 1.0e-10);
        assert!(inverse.residual < 1.0e-10);
        let bridge = SamplingChannelBridge::new(compiled).unwrap();
        let bridged = bridge
            .forward(
                SamplingChannelId::from(0),
                &[0.31, 0.42, 0.57, 0.23, 0.68, 0.81],
            )
            .unwrap();
        assert!((bridged.partition.weight_sum() - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn product_rejects_overlapping_child_blocks_with_parent_diagnostic() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["joint".into()];
        let mut channel = definition("product(surface(2), complement(2))");
        channel.subspace_lmb = vec![2];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("joint".into(), channel);
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let mut context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context
            .insert_geometry_map(
                SamplingMapDefinition::Surface(vec![2]),
                context.parent_lmb.clone(),
                vec![2],
                vec![],
                CompiledSamplingMap::Surface(
                    SurfaceRadialMap::new(3 * (vec![2]).len(), vec![0.0; 3], Some(3.0), 2.0, 1.0)
                        .unwrap(),
                ),
            )
            .unwrap();
        let error = catalogue
            .compile(
                &context,
                &catalogue
                    .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                    .unwrap(),
            )
            .unwrap_err();
        let diagnostic = error.to_string();
        assert!(diagnostic.contains("overlapping child blocks on edge 2"));
        assert!(diagnostic.contains("parent LMB [1, 2]"));
    }

    #[test]
    fn bridge_pushes_selected_channel_and_partitions_the_raw_frame() {
        let settings = ParameterizationSettings::default();
        let map = SamplingMapKernel::new(
            SamplingMapDefinition::Lmb(vec![1]),
            settings.clone(),
            100.0,
            1,
        )
        .unwrap();
        let second =
            SamplingMapKernel::new(SamplingMapDefinition::Lmb(vec![1]), settings, 100.0, 1)
                .unwrap();
        let bridge = SamplingChannelBridge::new(vec![
            CompiledSamplingChannel {
                name: "left".into(),
                master_graph: "G".into(),
                basis_id: Some(0),
                definition: SamplingMapDefinition::Lmb(vec![1]),
                embedded_edges: vec![1],
                map: CompiledSamplingMap::Lmb(map),
                singularity_proxy: None,
            },
            CompiledSamplingChannel {
                name: "right".into(),
                master_graph: "G".into(),
                basis_id: Some(1),
                definition: SamplingMapDefinition::Lmb(vec![1]),
                embedded_edges: vec![1],
                map: CompiledSamplingMap::Lmb(second),
                singularity_proxy: None,
            },
        ])
        .unwrap();
        let evaluation = bridge
            .forward(SamplingChannelId::from(1), &[0.31, 0.42, 0.57])
            .unwrap();
        assert_eq!(evaluation.raw_coordinates, evaluation.map.point);
        assert_eq!(evaluation.partition.weights.len(), 2);
        assert!((evaluation.partition.weight_sum() - 1.0).abs() < 1.0e-12);
        assert!(
            evaluation
                .partition
                .weights
                .iter()
                .all(|weight| *weight > 0.0)
        );
        let inverse = bridge
            .inverse(SamplingChannelId::from(1), &evaluation.raw_coordinates)
            .unwrap()
            .expect("full-support inverse");
        assert!(inverse.map.residual < 1.0e-10);

        let externals = Externals::default();
        let sample = evaluation
            .to_momentum_sample(SamplingMomentumSampleContext {
                loop_mom_cache_id: 4,
                external_moms: &externals,
                external_mom_cache_id: 7,
                dependent_momenta_constructor: DependentMomentaConstructor::CrossSection,
                orientation: Some(2),
            })
            .unwrap();
        assert_eq!(sample.loop_moms().0.len(), 1);
        assert_eq!(sample.sample.orientation, Some(2));
        assert!((sample.jacobian().0 - evaluation.map.jacobian).abs() < 1.0e-12);
    }

    #[test]
    fn canonical_lmb_bridge_integrates_normalized_gaussian_with_map_partition() {
        // Build both channels through the canonical catalogue.  The selected
        // map is sampled uniformly in channel space; multiplying by
        // `n_channels * J_i * w_i` converts that sample to the exact
        // map-density mixture estimator, where w_i = rho_i / sum_j rho_j.
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["auto:lmb".into()];
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue =
            build_sampling_channel_catalogue(&resolved, &[(0, vec![1]), (1, vec![1])], &[0, 1]);
        let context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1],
            ParameterizationSettings::default(),
            2.0,
            1,
        );
        let bridge = SamplingChannelBridge::new(
            catalogue
                .compile(
                    &context,
                    &catalogue
                        .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                        .unwrap(),
                )
                .unwrap(),
        )
        .unwrap();

        let bins = 24usize;
        let normalisation = (2.0 * std::f64::consts::PI).powf(-1.5);
        let channel_count = bridge.channels().len() as f64;
        let mut integral = 0.0;
        for i in 0..bins {
            for j in 0..bins {
                for k in 0..bins {
                    let coordinates = [
                        (i as f64 + 0.5) / bins as f64,
                        (j as f64 + 0.5) / bins as f64,
                        (k as f64 + 0.5) / bins as f64,
                    ];
                    let channel = SamplingChannelId::from((i + j + k) % bridge.channels().len());
                    let evaluation = bridge.forward(channel, &coordinates).unwrap();
                    let radius_squared = evaluation
                        .raw_coordinates
                        .iter()
                        .map(|component| component.powi(2))
                        .sum::<f64>();
                    let target = normalisation * (-0.5 * radius_squared).exp();
                    let partition_weight = evaluation.partition.weight(channel.index()).unwrap();
                    integral += target * evaluation.map.jacobian * channel_count * partition_weight;
                }
            }
        }
        integral /= bins.pow(3) as f64;
        assert!(
            (integral - 1.0).abs() < 2.0e-2,
            "canonical map-density mixture integral = {integral}"
        );
    }

    #[test]
    fn bridge_acceptance_report_integrates_every_canonical_channel() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["auto:lmb".into()];
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue =
            build_sampling_channel_catalogue(&resolved, &[(0, vec![1]), (1, vec![1])], &[0, 1]);
        let context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1],
            ParameterizationSettings::default(),
            2.0,
            1,
        );
        let bridge = SamplingChannelBridge::new(
            catalogue
                .compile(
                    &context,
                    &catalogue
                        .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                        .unwrap(),
                )
                .unwrap(),
        )
        .unwrap();
        let report = SamplingChannelBridgeAcceptanceReport::normalized_gaussian(
            &bridge,
            1024,
            1.0,
            &[0.0, 0.0, 0.0],
        )
        .unwrap();
        assert_eq!(report.sample_count, 1024);
        assert_eq!(report.channel_count, 2);
        assert_eq!(report.finite_sample_count, 2048);
        assert!((report.partition_min - 1.0).abs() < 1.0e-12);
        assert!((report.partition_max - 1.0).abs() < 1.0e-12);
        assert!((report.normalization - 1.0).abs() < 5.0e-2);
        assert!(report.normalization_stderr.is_finite());
        assert!((report.second_moment - report.expected_second_moment).abs() < 2.0e-1);
        assert!(report.second_moment_stderr.is_finite());
        assert!(report.round_trip_residual_max < 1.0e-10);
    }

    #[test]
    fn bridge_acceptance_report_integrates_absent_surface_and_lmb_catalogue() {
        // A named rootless map and an automatically generated LMB map share the same
        // master frame but remain two distinct canonical catalogue entries.
        // The acceptance harness must visit both IDs; silently compacting the
        // generated basis would make the mixed selection under-sample one map.
        let mut selection = SamplingChannelSelection {
            default_channel_selection: vec!["auto:lmb".into(), "absent".into()],
            ..Default::default()
        };
        let mut absent = definition("surface(2,3)");
        absent.subspace_lmb = vec![1];
        absent.parent_lmb = vec![1];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("absent".into(), absent);
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[(0, vec![1])], &[0]);
        assert_eq!(catalogue.entries.len(), 2);
        let mut context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1],
            ParameterizationSettings::default(),
            2.0,
            1,
        );
        context
            .insert_geometry_map(
                SamplingMapDefinition::Surface(vec![2, 3]),
                context.parent_lmb.clone(),
                vec![1],
                vec![],
                CompiledSamplingMap::Surface(
                    SurfaceRadialMap::new(
                        3 * (vec![1]).len(),
                        vec![0.2, -0.3, 0.1],
                        None,
                        2.0,
                        1.0,
                    )
                    .unwrap(),
                ),
            )
            .unwrap();
        let bridge = SamplingChannelBridge::new(
            catalogue
                .compile(
                    &context,
                    &catalogue
                        .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                        .unwrap(),
                )
                .unwrap(),
        )
        .unwrap();
        let report = SamplingChannelBridgeAcceptanceReport::normalized_gaussian(
            &bridge,
            512,
            1.0,
            &[0.0, 0.0, 0.0],
        )
        .unwrap();
        assert_eq!(report.channel_count, 2);
        assert_eq!(report.finite_sample_count, 1024);
        assert!((report.normalization - 1.0).abs() < 5.0e-2);
        assert!((report.partition_min - 1.0).abs() < 1.0e-12);
        assert!((report.partition_max - 1.0).abs() < 1.0e-12);
        assert!(report.round_trip_residual_max < 1.0e-10);
    }

    #[test]
    fn bridge_acceptance_report_supports_a_single_canonical_channel() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["auto:lmb".into()];
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[(0, vec![1])], &[0]);
        let context = SamplingChannelCompileContext::<f64>::new(
            "G",
            vec![1],
            ParameterizationSettings::default(),
            2.0,
            1,
        );
        let bridge = SamplingChannelBridge::new(
            catalogue
                .compile(
                    &context,
                    &catalogue
                        .compile_programs(3 * context.n_loop_momenta, &HFunctionSettings::default())
                        .unwrap(),
                )
                .unwrap(),
        )
        .unwrap();
        let report = SamplingChannelBridgeAcceptanceReport::normalized_gaussian(
            &bridge,
            512,
            1.0,
            &[0.0, 0.0, 0.0],
        )
        .unwrap();
        assert_eq!(report.channel_count, 1);
        assert_eq!(report.finite_sample_count, report.sample_count);
        assert!((report.normalization - 1.0).abs() < 5.0e-2);
        assert!((report.partition_min - 1.0).abs() < 1.0e-12);
        assert!((report.partition_max - 1.0).abs() < 1.0e-12);
        assert!(report.round_trip_residual_max < 1.0e-10);
        // The density is centered, but the reported observable is raw |K|^2.
        // This displacement separates that moment from the centered variance.
        let shifted = SamplingChannelBridgeAcceptanceReport::normalized_gaussian(
            &bridge,
            8192,
            1.0,
            &[1.0, -0.5, 0.75],
        )
        .unwrap();
        assert_eq!(shifted.expected_second_moment, 4.8125);
        assert!((shifted.normalization - 1.0).abs() < 0.025, "{shifted:?}");
        assert!((shifted.second_moment - 4.8125).abs() < 0.12, "{shifted:?}");
    }

    #[test]
    fn bridge_rejects_partial_or_mismatched_channels() {
        let settings = ParameterizationSettings::default();
        let map = SamplingMapKernel::new(SamplingMapDefinition::Lmb(vec![1]), settings, 100.0, 1)
            .unwrap();
        let channel = CompiledSamplingChannel {
            name: "one".into(),
            master_graph: "G".into(),
            basis_id: Some(0),
            definition: SamplingMapDefinition::Lmb(vec![1]),
            embedded_edges: vec![1],
            map: CompiledSamplingMap::Lmb(map),
            singularity_proxy: None,
        };
        let bridge = SamplingChannelBridge::new(vec![channel.clone()]).unwrap();
        assert!(
            bridge
                .forward(SamplingChannelId::from(1), &[0.2, 0.3, 0.4])
                .is_err()
        );
        let mut malformed = channel;
        malformed.name = "wrong".into();
        // A map with a different loop count cannot be put into the common raw frame.
        let other = SamplingMapKernel::new(
            SamplingMapDefinition::Lmb(vec![1, 2]),
            ParameterizationSettings::default(),
            100.0,
            2,
        )
        .unwrap();
        malformed.map = CompiledSamplingMap::Lmb(other);
        assert!(matches!(
            SamplingChannelBridge::new(vec![malformed, bridge.channels()[0].clone()]),
            Err(SamplingChannelBridgeError::DimensionMismatch { .. })
        ));

        let different_frame = CompiledSamplingChannel {
            name: "different-frame".into(),
            master_graph: "G".into(),
            basis_id: Some(2),
            definition: SamplingMapDefinition::Lmb(vec![2]),
            embedded_edges: vec![2],
            map: CompiledSamplingMap::Lmb(
                SamplingMapKernel::new(
                    SamplingMapDefinition::Lmb(vec![2]),
                    ParameterizationSettings::default(),
                    100.0,
                    1,
                )
                .unwrap(),
            ),
            singularity_proxy: None,
        };
        assert!(matches!(
            SamplingChannelBridge::new(vec![bridge.channels()[0].clone(), different_frame]),
            Err(SamplingChannelBridgeError::IncompatibleFrames { .. })
        ));
    }
}
