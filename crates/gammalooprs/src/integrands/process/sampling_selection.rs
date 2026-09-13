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

use crate::settings::runtime::{SamplingChannelDefinition, SamplingChannelSelection};
use color_eyre::eyre::Result;

use super::{
    SamplingMapComponent, SamplingMapContract, SamplingMapDefinition, SamplingMapEvaluation,
    SamplingMapKernel, SurfaceRadialMap,
};
use crate::settings::runtime::ParameterizationSettings;

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
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResolvedNamedSamplingChannel {
    pub name: String,
    pub definition: SamplingChannelDefinition,
    pub map: SamplingMapDefinition,
}

/// The complete selection for one graph after applying the settings' fallback.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResolvedSamplingChannelSelection {
    pub graph_name: String,
    pub selectors: Vec<SamplingChannelSelector>,
    pub named_channels: Vec<ResolvedNamedSamplingChannel>,
}

/// One entry in the graph-local sampling catalogue.  The catalogue is only
/// an inspection/resolution product at this stage; ordinary LMB sampling keeps
/// using `LmbMultiChannelingSetup` unchanged.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SamplingCatalogueEntry {
    Lmb {
        basis_id: usize,
        edges: Vec<usize>,
        preset: SamplingChannelPreset,
    },
    Named(ResolvedNamedSamplingChannel),
}

/// Graph-scoped, deterministically ordered channel catalogue.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SamplingChannelCatalogue {
    pub graph_name: String,
    pub selectors: Vec<SamplingChannelSelector>,
    pub entries: Vec<SamplingCatalogueEntry>,
}

/// Kinematic data needed when compiling a graph-local surface channel.
///
/// The centre and threshold radius are intentionally supplied by the process
/// layer: they depend on the prepared external/cut kinematics and cannot be
/// inferred from a symbolic `surface(...)` expression alone.  `loop_edges`
/// are master-graph edge ids and define the local raw-frame embedding of the
/// compiled block.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingSurfaceGeometry {
    pub center: Vec<f64>,
    pub threshold_radius: Option<f64>,
    pub beta: f64,
    pub power: f64,
}

/// Input context for graph-independent compilation of catalogue entries.
///
/// A caller must identify the master graph explicitly.  This prevents a
/// surface block from being mistaken for a coordinate block in a different
/// graph when channel definitions are shared across a graph group.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingChannelCompileContext {
    pub master_graph: String,
    /// Complete ordered parent LMB in the master graph frame. Every named
    /// channel is resolved against this exact list before compilation.
    pub parent_lmb: Vec<usize>,
    pub parameterization_settings: ParameterizationSettings,
    pub e_cm: f64,
    pub n_loop_momenta: usize,
    pub surfaces: BTreeMap<Vec<usize>, SamplingSurfaceGeometry>,
}

impl SamplingChannelCompileContext {
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
            surfaces: BTreeMap::new(),
        }
    }
}

/// The map kernels currently compilable without a graph-specific implicit
/// solver.  Composite/cut maps stay in the typed catalogue until their
/// prepared kinematic context is supplied by the process layer.
#[derive(Clone, Debug)]
pub enum CompiledSamplingMap {
    Lmb(SamplingMapKernel),
    Surface(SurfaceRadialMap),
}

impl CompiledSamplingMap {
    pub fn contract(&self) -> SamplingMapContract {
        match self {
            Self::Lmb(map) => map.contract(),
            Self::Surface(map) => map.contract(),
        }
    }

    pub fn dimensions(&self) -> usize {
        match self {
            Self::Lmb(map) => map.dimensions(),
            Self::Surface(map) => map.dimension(),
        }
    }

    pub fn as_component(&self) -> &dyn SamplingMapComponent {
        match self {
            Self::Lmb(map) => map,
            Self::Surface(map) => map,
        }
    }

    pub fn forward(&self, coordinates: &[f64]) -> Result<SamplingMapEvaluation> {
        self.as_component().forward(coordinates, &[])
    }

    pub fn inverse(&self, point: &[f64]) -> Result<SamplingMapEvaluation> {
        self.as_component().inverse(point, &[])
    }
}

/// A compiled graph channel and the master-graph raw-frame block it occupies.
#[derive(Clone, Debug)]
pub struct CompiledSamplingChannel {
    pub name: String,
    pub master_graph: String,
    pub basis_id: Option<usize>,
    pub definition: SamplingMapDefinition,
    /// Ordered master-graph edge ids for this raw coordinate block.
    pub embedded_edges: Vec<usize>,
    pub map: CompiledSamplingMap,
}

impl CompiledSamplingChannel {
    pub fn contract(&self) -> SamplingMapContract {
        self.map.contract()
    }

    pub fn dimensions(&self) -> usize {
        self.map.dimensions()
    }

    pub fn forward(&self, coordinates: &[f64]) -> Result<SamplingMapEvaluation> {
        self.map.forward(coordinates)
    }

    pub fn inverse(&self, point: &[f64]) -> Result<SamplingMapEvaluation> {
        self.map.inverse(point)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SamplingChannelCompileError {
    EmptyMasterGraph,
    UnsupportedMap { channel: String, map: String },
    MissingSurfaceGeometry { channel: String, edges: Vec<usize> },
    InvalidChannel { channel: String, error: String },
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
            Self::MissingSurfaceGeometry { channel, edges } => write!(
                formatter,
                "sampling channel `{channel}` has no prepared surface geometry for master-graph edges {edges:?}"
            ),
            Self::InvalidChannel { channel, error } => write!(
                formatter,
                "sampling channel `{channel}` could not be compiled: {error}"
            ),
        }
    }
}

impl std::error::Error for SamplingChannelCompileError {}

impl SamplingChannelCatalogue {
    pub fn lmb_entries(&self) -> impl Iterator<Item = (usize, &[usize])> {
        self.entries.iter().filter_map(|entry| match entry {
            SamplingCatalogueEntry::Lmb {
                basis_id, edges, ..
            } => Some((*basis_id, edges.as_slice())),
            SamplingCatalogueEntry::Named(_) => None,
        })
    }

    pub fn named_entries(&self) -> impl Iterator<Item = &ResolvedNamedSamplingChannel> {
        self.entries.iter().filter_map(|entry| match entry {
            SamplingCatalogueEntry::Named(channel) => Some(channel),
            SamplingCatalogueEntry::Lmb { .. } => None,
        })
    }

    /// Compile the ordinary LMB and radial surface entries in this catalogue.
    ///
    /// Surface geometry is looked up by its canonical master-graph edge list.
    /// Maps such as `cut`, `intersect` and `then` are deliberately rejected
    /// here until the process has prepared their conditional kinematics; this
    /// avoids creating a numerically plausible map with an incorrect frame.
    pub fn compile(
        &self,
        context: &SamplingChannelCompileContext,
    ) -> Result<Vec<CompiledSamplingChannel>, SamplingChannelCompileError> {
        if context.master_graph.trim().is_empty() {
            return Err(SamplingChannelCompileError::EmptyMasterGraph);
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
            });
        }
        let mut compiled = Vec::with_capacity(self.entries.len());
        for entry in &self.entries {
            let (name, basis_id, definition, map) = match entry {
                SamplingCatalogueEntry::Lmb {
                    basis_id, edges, ..
                } => {
                    let definition = SamplingMapDefinition::Lmb(edges.clone());
                    let map = SamplingMapKernel::new(
                        definition.clone(),
                        context.parameterization_settings.clone(),
                        context.e_cm,
                        context.n_loop_momenta,
                    )
                    .map_err(|error| {
                        SamplingChannelCompileError::InvalidChannel {
                            channel: format!("lmb[{basis_id}]"),
                            error: error.to_string(),
                        }
                    })?;
                    (
                        format!("lmb[{basis_id}]"),
                        Some(*basis_id),
                        definition,
                        CompiledSamplingMap::Lmb(map),
                    )
                }
                SamplingCatalogueEntry::Named(channel) => {
                    if channel.definition.parent_lmb != context.parent_lmb {
                        return Err(SamplingChannelCompileError::InvalidChannel {
                            channel: channel.name.clone(),
                            error: format!(
                                "parent LMB {:?} does not match master context {:?}",
                                channel.definition.parent_lmb, context.parent_lmb
                            ),
                        });
                    }
                    let definition = channel.map.clone();
                    let map = match &definition {
                        SamplingMapDefinition::Lmb(_edges) => SamplingMapKernel::new(
                            definition.clone(),
                            context.parameterization_settings.clone(),
                            context.e_cm,
                            context.n_loop_momenta,
                        )
                        .map(CompiledSamplingMap::Lmb)
                        .map_err(|error| SamplingChannelCompileError::InvalidChannel {
                            channel: channel.name.clone(),
                            error: error.to_string(),
                        })?,
                        SamplingMapDefinition::Surface(edges) => {
                            if edges.len() != context.n_loop_momenta {
                                return Err(SamplingChannelCompileError::UnsupportedMap {
                                    channel: channel.name.clone(),
                                    map: format!(
                                        "surface({edges:?}) is only a full-frame map; supply an explicit product(surface(...), complement(...)) with prepared embedding for a proper subspace"
                                    ),
                                });
                            }
                            let Some(geometry) = context.surfaces.get(edges) else {
                                return Err(SamplingChannelCompileError::MissingSurfaceGeometry {
                                    channel: channel.name.clone(),
                                    edges: edges.clone(),
                                });
                            };
                            let expected_dimension = 3 * edges.len();
                            if geometry.center.len() != expected_dimension {
                                return Err(SamplingChannelCompileError::InvalidChannel {
                                    channel: channel.name.clone(),
                                    error: format!(
                                        "surface centre has dimension {}, expected {} for edges {edges:?}",
                                        geometry.center.len(),
                                        expected_dimension
                                    ),
                                });
                            }
                            SurfaceRadialMap::new(
                                expected_dimension,
                                geometry.center.clone(),
                                geometry.threshold_radius,
                                geometry.beta,
                                geometry.power,
                            )
                            .map(CompiledSamplingMap::Surface)
                            .map_err(|error| {
                                SamplingChannelCompileError::InvalidChannel {
                                    channel: channel.name.clone(),
                                    error: error.to_string(),
                                }
                            })?
                        }
                        unsupported => {
                            return Err(SamplingChannelCompileError::UnsupportedMap {
                                channel: channel.name.clone(),
                                map: format!("{unsupported:?}"),
                            });
                        }
                    };
                    (channel.name.clone(), None, definition, map)
                }
            };
            let embedded_edges = match &definition {
                SamplingMapDefinition::Lmb(edges)
                | SamplingMapDefinition::Surface(edges)
                | SamplingMapDefinition::Complement(edges)
                | SamplingMapDefinition::Cut(edges) => edges.clone(),
                _ => Vec::new(),
            };
            compiled.push(CompiledSamplingChannel {
                name,
                master_graph: context.master_graph.clone(),
                basis_id,
                definition,
                embedded_edges,
                map,
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
                SamplingCatalogueEntry::Named(channel) => format!(
                    "{index}: {} around={} parent_lmb={:?} on_cut={:?}",
                    channel.name,
                    channel.definition.around,
                    channel.definition.parent_lmb,
                    channel.definition.on_cut
                ),
            })
            .collect()
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
        let basis_ids: Vec<usize> = match preset {
            SamplingChannelPreset::Lmb => all_lmbs.iter().map(|(id, _)| *id).collect(),
            SamplingChannelPreset::OptimizedLmb | SamplingChannelPreset::Surfaces => {
                optimized_lmbs.to_vec()
            }
        };
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
/// The default and graph-specific entries form an additive union.  Duplicates
/// are removed while preserving the first occurrence's order.  This is useful
/// to callers constructing a catalogue for several graphs at once: a common
/// baseline can be supplied by the default entry and graph-specific channels
/// can be appended without repeating it.
pub fn resolve_sampling_channel_selection(
    graph_name: &str,
    selection: &SamplingChannelSelection,
) -> Result<ResolvedSamplingChannelSelection, SamplingSelectionError> {
    resolve_selection(graph_name, selection, false)
}

/// Resolve using the explicit graph entry as a replacement for the defaults.
///
/// The public TOML contract currently uses the additive resolver above.  This
/// variant is retained for catalogue builders that need strict per-graph
/// replacement semantics without reimplementing validation.
pub fn resolve_sampling_channel_selection_replacing_default(
    graph_name: &str,
    selection: &SamplingChannelSelection,
) -> Result<ResolvedSamplingChannelSelection, SamplingSelectionError> {
    resolve_selection(graph_name, selection, true)
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
    let raw_selectors = (!replace_default)
        .then_some(selection.default_channel_selection.as_slice())
        .into_iter()
        .flatten()
        .chain(graph_selectors.into_iter().flatten())
        .collect::<Vec<_>>();
    let mut selectors = Vec::with_capacity(raw_selectors.len());
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
        let map = SamplingMapDefinition::parse(&definition.around).map_err(|error| {
            SamplingSelectionError::InvalidChannelDefinition {
                graph: graph_name.to_owned(),
                channel: name.clone(),
                error: error.to_string(),
            }
        })?;
        named_channels.push(ResolvedNamedSamplingChannel {
            name: name.clone(),
            definition: definition.clone(),
            map,
        });
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
    use super::*;
    use crate::settings::runtime::ParameterizationSettings;

    fn definition(around: &str) -> SamplingChannelDefinition {
        SamplingChannelDefinition {
            around: around.to_owned(),
            parent_lmb: vec![1, 2],
            on_cut: vec![],
        }
    }

    #[test]
    fn graph_entry_unions_default_and_deduplicates_entries() {
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
                SamplingChannelSelector::Preset(SamplingChannelPreset::Lmb),
                SamplingChannelSelector::Named("common".into()),
                SamplingChannelSelector::Preset(SamplingChannelPreset::Surfaces),
                SamplingChannelSelector::Named("named".into())
            ]
        );
        assert!(resolved.named("named").is_some());
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
        assert_eq!(catalogue.lmb_entries().next(), Some((1, &[2, 4][..])));
        assert_eq!(catalogue.named_entries().next().unwrap().name, "surface_hz");
        assert!(catalogue.inspection_rows()[0].contains("basis=1"));
        assert!(catalogue.inspection_rows()[1].contains("parent_lmb=[1, 2]"));
    }

    #[test]
    fn surface_preset_has_optimized_lmb_fallback() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["auto:surfaces".into(), "auto:lmb".into()];
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue =
            build_sampling_channel_catalogue(&resolved, &[(0, vec![1]), (1, vec![2])], &[0, 1]);
        assert_eq!(catalogue.lmb_entries().count(), 2);
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
        let catalogue =
            build_sampling_channel_catalogue(&resolved, &[(0, vec![1, 2]), (1, vec![2, 4])], &[0]);
        let mut context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        context.surfaces.insert(
            vec![1, 2],
            SamplingSurfaceGeometry {
                center: vec![0.0; 6],
                threshold_radius: Some(3.0),
                beta: 2.0,
                power: 1.0,
            },
        );
        let compiled = catalogue.compile(&context).unwrap();
        assert_eq!(compiled.len(), 3);
        assert!(matches!(compiled[0].map, CompiledSamplingMap::Lmb(_)));
        assert_eq!(compiled[0].embedded_edges, vec![1, 2]);
        assert!(matches!(compiled[2].map, CompiledSamplingMap::Surface(_)));
        assert_eq!(compiled[2].master_graph, "G");
        assert_eq!(compiled[2].dimensions(), 6);
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
        let context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        assert!(matches!(
            catalogue.compile(&context),
            Err(SamplingChannelCompileError::MissingSurfaceGeometry { .. })
        ));
    }

    #[test]
    fn unresolved_composite_map_is_rejected_without_context() {
        let mut selection = SamplingChannelSelection::default();
        selection.default_channel_selection = vec!["joint".into()];
        selection
            .channel_definitions
            .entry("G".into())
            .or_default()
            .insert("joint".into(), definition("product(surface(1), lmb(1,2))"));
        let resolved = resolve_sampling_channel_selection("G", &selection).unwrap();
        let catalogue = build_sampling_channel_catalogue(&resolved, &[], &[]);
        let context = SamplingChannelCompileContext::new(
            "G",
            vec![1, 2],
            ParameterizationSettings::default(),
            100.0,
            2,
        );
        assert!(matches!(
            catalogue.compile(&context),
            Err(SamplingChannelCompileError::UnsupportedMap { .. })
        ));
    }
}
