//! Resolution of user-facing graph sampling channel selections.
//!
//! Settings deliberately keep selectors as strings so that TOML and Python
//! have the same compact interface.  This module turns those strings into a
//! typed, deterministic selection after the graph name is known.  Numerical
//! channel construction is intentionally kept out of this layer.

use std::{collections::BTreeMap, fmt, str::FromStr};

use crate::settings::runtime::{SamplingChannelDefinition, SamplingChannelSelection};

use super::SamplingMapDefinition;

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
}
