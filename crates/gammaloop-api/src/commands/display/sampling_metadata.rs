use std::{collections::BTreeMap, fmt::Write};

use color_eyre::Result;
use colored::Colorize;
use eyre::{eyre, Context};
use gammalooprs::{
    integrands::process::{
        ProcessIntegrand, SamplingCatalogueEntry, SamplingChannelCatalogue,
        SamplingChannelSelector, SamplingMapDefinition,
    },
    processes::GraphGroupSelectionSpec,
    settings::runtime::{SamplingChannelWeight, SamplingSettingsParser, SumMode},
    utils::{serde_utils::ShowDefaultsGuard, symbolica_ext::LOGPRINTOPTS},
};
use symbolica::atom::AtomCore;
use tabled::{
    builder::Builder,
    settings::{object::Columns, Modify, Style, Width},
};

use crate::{settings_tree::serialize_settings_with_defaults, state::State};

use super::{format_json_value, MetadataDisplayFormat};

pub(super) fn render(
    state: &State,
    process_id: usize,
    integrand_name: &str,
    requested_graphs: &[String],
    format: MetadataDisplayFormat,
) -> Result<String> {
    let resolved = state
        .process_list
        .get_integrand(process_id, integrand_name)?;
    let integrand = resolved.require_generated()?;
    let metadata = SamplingMetadata::collect(integrand, requested_graphs)?;
    match format {
        MetadataDisplayFormat::Pretty => {
            let heading = format!(
                "Sampling for '{}' in process #{}",
                resolved.canonical_name, process_id
            );
            Ok(format!("{}\n{}", heading.bold().blue(), metadata.pretty()?))
        }
        MetadataDisplayFormat::Toml => {
            let sampling = metadata.export_settings()?;
            metadata.validate_export(integrand, &sampling)?;
            SamplingMetadata::toml(&sampling)
        }
    }
}

struct SamplingMetadata {
    sampling: SamplingSettingsParser,
    graphs: Vec<GraphSamplingMetadata>,
}

struct GraphSamplingMetadata {
    name: String,
    graph_id: usize,
    group_id: usize,
    /// Index in the integration view after applying sampling.graph_names.
    graph_index: Option<usize>,
    orientations: usize,
    parent_lmb: Vec<usize>,
    catalogue: Option<SamplingChannelCatalogue>,
}

impl SamplingMetadata {
    fn collect(integrand: &ProcessIntegrand, requested_graphs: &[String]) -> Result<Self> {
        for name in requested_graphs {
            integrand.resolve_group_id_by_master_name(name)?;
        }
        let settings = &integrand.get_settings().sampling;
        let sampling = settings.as_parser();
        let parameterization = settings.get_parameterization_settings();
        // Reuse the same subset plan as integration, without cloning or warming an integrand.
        let plan = if settings.selected_graph_names().is_empty() {
            None
        } else {
            let selection = GraphGroupSelectionSpec::from_master_graph_names(
                settings.selected_graph_names().to_vec(),
            );
            Some(match integrand {
                ProcessIntegrand::Amplitude(amplitude) => selection
                    .plan(&amplitude.data.graph_group_structure, |id| {
                        amplitude.data.graph_terms.get(id).map(|term| &term.graph)
                    })?,
                ProcessIntegrand::CrossSection(cross_section) => {
                    selection.plan(&cross_section.data.graph_group_structure, |id| {
                        cross_section
                            .data
                            .graph_terms
                            .get(id)
                            .map(|term| &term.graph)
                    })?
                }
            })
        };
        let mut graphs = Vec::new();
        for name in integrand.graph_group_master_names() {
            if !requested_graphs.is_empty() && !requested_graphs.iter().any(|item| item == name) {
                continue;
            }
            let group_id = integrand.resolve_group_id_by_master_name(name)?;
            let graph_id = integrand
                .find_graph_id_by_name(name)
                .expect("known master graph");
            let setup = match integrand {
                ProcessIntegrand::Amplitude(amplitude) => {
                    &amplitude.data.graph_terms[graph_id].multi_channeling_setup
                }
                ProcessIntegrand::CrossSection(cross_section) => {
                    &cross_section.data.graph_terms[graph_id].multi_channeling_setup
                }
            };
            let catalogue = parameterization
                .as_ref()
                .map(|parameters| setup.canonical_sampling_catalogue(name, parameters))
                .transpose()
                .with_context(|| format!("While resolving sampling for master graph '{name}'"))?;
            if !sampling.sampling_multichanneling {
                if let Some(parameters) = &parameterization {
                    setup.selected_lmb_basis_id(name, parameters)?;
                }
            }
            graphs.push(GraphSamplingMetadata {
                name: name.to_owned(),
                graph_id,
                group_id: group_id.0,
                graph_index: plan.as_ref().map_or(Some(group_id.0), |plan| {
                    plan.new_group_id_for_old(group_id).map(|id| id.0)
                }),
                orientations: integrand.group_orientation_count(group_id).unwrap_or(0),
                parent_lmb: setup
                    .graph
                    .loop_momentum_basis
                    .loop_edges
                    .iter()
                    .map(|edge| edge.0)
                    .collect(),
                catalogue,
            });
        }
        Ok(Self { sampling, graphs })
    }

    fn export_settings(&self) -> Result<SamplingSettingsParser> {
        let mut sampling = self.sampling.clone();
        sampling.default_channel_selection.clear();
        sampling.channel_selection.clear();
        sampling.channel_definitions.clear();
        sampling
            .lmb_basis_ids
            .retain(|name, _| self.graphs.iter().any(|graph| &graph.name == name));
        for graph in &self.graphs {
            let Some(catalogue) = &graph.catalogue else {
                continue;
            };
            let mut selectors = Vec::with_capacity(catalogue.entries.len());
            let mut definitions = BTreeMap::new();
            for entry in &catalogue.entries {
                match entry {
                    SamplingCatalogueEntry::Lmb { edges, .. } => {
                        selectors.push(SamplingChannelSelector::Lmb(edges.clone()).to_string());
                    }
                    SamplingCatalogueEntry::Named(channel) => {
                        selectors.push(channel.name.clone());
                        definitions.insert(channel.name.clone(), channel.definition.clone());
                    }
                    SamplingCatalogueEntry::Surface { .. } => {
                        return Err(eyre!(
                            "Cannot export an automatic surface channel for '{}': the native selector syntax cannot preserve its generated identity.",
                            graph.name
                        ));
                    }
                }
            }
            sampling
                .channel_selection
                .insert(graph.name.clone(), selectors);
            if !definitions.is_empty() {
                sampling
                    .channel_definitions
                    .insert(graph.name.clone(), definitions);
            }
        }
        // Graph-specific arrays replace defaults on merge; unselected graph keys must be absent.
        // An empty shared default is omitted by its native serializer, preserving the live default.
        Ok(sampling)
    }

    fn validate_export(
        &self,
        integrand: &ProcessIntegrand,
        sampling: &SamplingSettingsParser,
    ) -> Result<()> {
        let Some(mut parameters) = integrand
            .get_settings()
            .sampling
            .get_parameterization_settings()
        else {
            return Ok(());
        };
        parameters.lmb_basis_ids = sampling.lmb_basis_ids.clone();
        parameters.sampling_channels.channel_selection = sampling.channel_selection.clone();
        parameters.sampling_channels.channel_definitions = sampling.channel_definitions.clone();
        for graph in &self.graphs {
            let setup = match integrand {
                ProcessIntegrand::Amplitude(amplitude) => {
                    &amplitude.data.graph_terms[graph.graph_id].multi_channeling_setup
                }
                ProcessIntegrand::CrossSection(cross_section) => {
                    &cross_section.data.graph_terms[graph.graph_id].multi_channeling_setup
                }
            };
            let exported = setup.canonical_sampling_catalogue(&graph.name, &parameters)
                .with_context(|| format!("Cannot export sampling for '{}': explicit selectors must reproduce its canonical channel IDs", graph.name))?;
            let original = graph
                .catalogue
                .as_ref()
                .expect("parameterized graphs have a catalogue");
            if original.entries.len() != exported.entries.len()
                || !original.entries.iter().zip(&exported.entries).all(
                    |(original, exported)| match (original, exported) {
                        (
                            SamplingCatalogueEntry::Lmb {
                                basis_id: old_id,
                                edges: old_edges,
                                ..
                            },
                            SamplingCatalogueEntry::Lmb {
                                basis_id, edges, ..
                            },
                        ) => old_id == basis_id && old_edges == edges,
                        _ => original == exported,
                    },
                )
            {
                return Err(eyre!("Cannot export sampling for '{}': explicit selectors changed the canonical channel IDs or definitions", graph.name));
            }
        }
        Ok(())
    }

    fn toml(sampling: &SamplingSettingsParser) -> Result<String> {
        let document = BTreeMap::from([("sampling", sampling)]);
        let _defaults = ShowDefaultsGuard::new(true);
        toml::to_string(&document).context("While serializing explicit sampling settings")
    }

    fn pretty(&self) -> Result<String> {
        let mut output = String::new();
        let discrete_graphs = self.sampling.graphs == SumMode::MonteCarlo;
        let sampled_orientations = self.sampling.orientations == SumMode::MonteCarlo;
        let sampled_channels = self.sampling.sampling_multichanneling
            && self.sampling.sampling_channels == SumMode::MonteCarlo;
        let mut axes = Vec::new();
        if discrete_graphs {
            axes.push("graph index");
            if sampled_orientations {
                axes.push("orientation ID");
            }
            if sampled_channels {
                axes.push("channel ID");
            }
        }
        writeln!(
            output,
            "{} {}",
            "Discrete coordinates:".bold(),
            if axes.is_empty() {
                "none".to_owned()
            } else {
                format!("[{}]", axes.join(", "))
            }
        )?;
        writeln!(output, "{}", "Channel IDs, generated basis IDs, physical cut IDs, and stored graph/group IDs are separate zero-based indices.".dimmed())?;
        let values =
            serialize_settings_with_defaults(&self.sampling, "effective sampling settings")?;
        let mut summary = Builder::new();
        summary.push_record([
            "setting".bold().blue().to_string(),
            "effective value".bold().blue().to_string(),
        ]);
        for (name, value) in values
            .as_object()
            .expect("sampling settings serialize as an object")
        {
            if [
                "channel_selection",
                "channel_definitions",
                "default_channel_selection",
                "lmb_basis_ids",
            ]
            .contains(&name.as_str())
            {
                continue;
            }
            summary.push_record([
                name.cyan().to_string(),
                format_json_value(value).yellow().to_string(),
            ]);
        }
        let mut summary = summary.build();
        summary.with(Style::rounded().remove_horizontals());
        summary.with(Modify::new(Columns::one(1)).with(Width::wrap(76)));
        writeln!(output, "{summary}")?;
        for graph in &self.graphs {
            writeln!(
                output,
                "\n{}",
                format!(
                    "{}  ·  graph #{}  ·  group #{}",
                    graph.name, graph.graph_id, graph.group_id
                )
                .bold()
                .blue()
            )?;
            let state = if graph.graph_index.is_none() {
                "excluded by sampling.graph_names".yellow().to_string()
            } else if graph.catalogue.is_none() {
                "tropical map (no channel axis)".green().to_string()
            } else if !self.sampling.sampling_multichanneling {
                "default map: channel #0 only".green().to_string()
            } else if sampled_channels {
                "channels sampled".green().to_string()
            } else {
                "channels summed".green().to_string()
            };
            let index = if discrete_graphs {
                graph
                    .graph_index
                    .map_or("—".to_owned(), |id| id.to_string())
            } else {
                "none (graphs summed)".to_owned()
            };
            writeln!(
                output,
                "{state}; graph index: {}; runtime orientation entries: {} ({}); parent LMB: {}",
                index.yellow(),
                graph.orientations.to_string().yellow(),
                if sampled_orientations {
                    "sampled"
                } else {
                    "summed"
                },
                format!("{:?}", graph.parent_lmb).cyan()
            )?;
            let Some(catalogue) = &graph.catalogue else {
                writeln!(
                    output,
                    "{}",
                    "Tropical coordinates do not use the sampling-channel catalogue.".dimmed()
                )?;
                continue;
            };
            let original = self
                .sampling
                .channel_selection
                .get(&graph.name)
                .unwrap_or(&self.sampling.default_channel_selection);
            let mut table = Builder::new();
            table.push_record([
                "ID / basis".bold().blue().to_string(),
                "channel / source".bold().blue().to_string(),
                "map and settings".bold().blue().to_string(),
            ]);
            let graph_override = self.sampling.channel_selection.contains_key(&graph.name);
            table.push_record([
                String::new(),
                format!(
                    "selectors ({})",
                    if graph_override { "graph" } else { "default" }
                )
                .blue()
                .to_string(),
                if original.is_empty() {
                    if graph_override {
                        "[]"
                    } else {
                        "(implicit auto:optimized_lmb)"
                    }
                    .to_owned()
                } else {
                    original.join(", ").cyan().to_string()
                },
            ]);
            table.push_record([
                String::new(),
                "expanded channels".blue().to_string(),
                catalogue
                    .entries
                    .iter()
                    .map(|entry| match entry {
                        SamplingCatalogueEntry::Lmb { edges, .. } => {
                            SamplingChannelSelector::Lmb(edges.clone())
                                .to_string()
                                .magenta()
                                .to_string()
                        }
                        SamplingCatalogueEntry::Named(channel) => channel.name.cyan().to_string(),
                        SamplingCatalogueEntry::Surface { edges, .. } => {
                            SamplingMapDefinition::Surface(edges.clone())
                                .to_atom()
                                .printer(LOGPRINTOPTS.clone())
                                .to_string()
                        }
                    })
                    .collect::<Vec<_>>()
                    .join(", "),
            ]);
            if let Some(restrictions) = self.sampling.lmb_basis_ids.get(&graph.name) {
                table.push_record([
                    String::new(),
                    "admissible basis IDs".blue().to_string(),
                    format!("{restrictions:?}").yellow().to_string(),
                ]);
            }
            if catalogue.entries.is_empty() {
                table.push_record([
                    String::new(),
                    "selection".blue().to_string(),
                    "No selected sampling channels.".yellow().to_string(),
                ]);
            }
            if !self.sampling.sampling_multichanneling {
                table.push_record([
                    String::new(),
                    "usage".blue().to_string(),
                    "Default sampling uses channel #0; remaining catalogue entries and multichannel weights are inactive.".to_owned(),
                ]);
            }
            let mut profiles = Vec::new();
            for (id, entry) in catalogue.entries.iter().enumerate() {
                let (basis, label, mut detail, weight) = match entry {
                    SamplingCatalogueEntry::Lmb {
                        basis_id,
                        edges,
                        source,
                    } => (
                        format!("\nbasis #{basis_id}"),
                        source.to_string(),
                        SamplingMapDefinition::Lmb(edges.clone())
                            .to_atom()
                            .printer(LOGPRINTOPTS.clone())
                            .to_string()
                            .magenta()
                            .to_string(),
                        self.sampling.sampling_channel_weight,
                    ),
                    SamplingCatalogueEntry::Surface { edges, parent_lmb } => (
                        String::new(),
                        "auto:surfaces".to_owned(),
                        format!(
                            "{}\nparent: {parent_lmb:?}",
                            SamplingMapDefinition::Surface(edges.clone())
                                .to_atom()
                                .printer(LOGPRINTOPTS.clone())
                                .to_string()
                                .magenta()
                        ),
                        self.sampling.sampling_channel_weight,
                    ),
                    SamplingCatalogueEntry::Named(channel) => {
                        let definition = &channel.definition;
                        let mut detail = definition.around.magenta().to_string();
                        write!(
                            detail,
                            "\nparent: {}; subspace: {}; physical cut IDs: {}",
                            format!("{:?}", definition.parent_lmb).cyan(),
                            format!("{:?}", definition.subspace_lmb).yellow(),
                            format!("{:?}", definition.on_cut).yellow()
                        )?;
                        if let Some(profile) = &definition.radial_profile {
                            let profile_id = profiles
                                .iter()
                                .position(|existing| *existing == profile)
                                .unwrap_or_else(|| {
                                    profiles.push(profile);
                                    profiles.len() - 1
                                });
                            write!(detail, "\nprofile: {}", format!("#{profile_id}").cyan())?;
                        }
                        if let Some(proxy) = &definition.singularity_proxy {
                            write!(detail, "\nproxy: {}", proxy.magenta())?;
                        }
                        for (index, block) in channel.blocks.iter().enumerate() {
                            let target = block
                                .target
                                .to_atom()
                                .printer(LOGPRINTOPTS.clone())
                                .to_string();
                            write!(
                                detail,
                                "\n{} {}\n  active: {}; preceding: {}; remaining: {}",
                                format!("block {index}:").bold().blue(),
                                target.magenta(),
                                format!("{:?}", block.active_lmb).yellow(),
                                format!("{:?}", block.preceding_lmb).cyan(),
                                format!("{:?}", block.remaining_lmb).dimmed()
                            )?;
                        }
                        (
                            String::new(),
                            channel.name.clone(),
                            detail,
                            definition
                                .channel_weight
                                .unwrap_or(self.sampling.sampling_channel_weight),
                        )
                    }
                };
                let weight = match weight {
                    SamplingChannelWeight::MapDensity => "map_density",
                    SamplingChannelWeight::InverseJacobian => "inverse_jacobian",
                    SamplingChannelWeight::Ose => "ose",
                    SamplingChannelWeight::SingularityProxy => "singularity_proxy",
                };
                write!(detail, "\nweight: {}", weight.green())?;
                let inactive = graph.graph_index.is_none()
                    || (!self.sampling.sampling_multichanneling && id != 0);
                let row = [format!("#{id}{basis}"), label, detail];
                table.push_record(if inactive {
                    row.map(|cell| {
                        console::strip_ansi_codes(&cell)
                            .as_ref()
                            .dimmed()
                            .to_string()
                    })
                } else {
                    let [id, label, detail] = row;
                    [id.yellow().to_string(), label.cyan().to_string(), detail]
                });
            }
            let mut table = table.build();
            table.with(Style::rounded().remove_horizontals());
            table.with(Modify::new(Columns::one(1)).with(Width::wrap(24)));
            table.with(Modify::new(Columns::one(2)).with(Width::wrap(76)));
            writeln!(output, "{table}")?;
            if !profiles.is_empty() {
                writeln!(
                    output,
                    "{}",
                    "Radial profiles (display IDs for this graph)".bold().blue()
                )?;
                let mut table = Builder::new();
                table.push_record([
                    "profile".bold().blue().to_string(),
                    "settings".bold().blue().to_string(),
                ]);
                for (id, profile) in profiles.iter().enumerate() {
                    let settings =
                        serialize_settings_with_defaults(profile, "sampling radial profile")?;
                    table.push_record([
                        format!("#{id}").cyan().to_string(),
                        serde_json::to_string(&settings)?.cyan().to_string(),
                    ]);
                }
                let mut table = table.build();
                table.with(Style::rounded().remove_horizontals());
                table.with(Modify::new(Columns::one(1)).with(Width::wrap(96)));
                writeln!(output, "{table}")?;
            }
        }
        writeln!(output, "{}", "Map recipes and cut hosts are static; roots, centers, and support depend on the sampled point.".dimmed())?;
        Ok(output)
    }
}

#[cfg(test)]
mod tests {
    use super::{GraphSamplingMetadata, SamplingMetadata};
    use std::collections::BTreeMap;

    use figment::{providers::Serialized, Figment};
    use gammalooprs::{
        integrands::process::{
            build_sampling_channel_catalogue, resolve_sampling_channel_selection,
            SamplingCatalogueEntry,
        },
        settings::{
            runtime::{
                SamplingChannelDefinition, SamplingRadialProfile, SamplingSettingsParser, SumMode,
            },
            RuntimeSettings,
        },
    };

    use crate::commands::set::SetArgs;

    const JOINT: &str = "then(block(lmb(6,10),phase_space(cut(2,6,10))),complement(7),block(lmb(3),at_cut(cut(2,6,10),intersect(surface(2,4,12),surface(3,10,13)))))";

    fn fixture(names: &[&str]) -> (RuntimeSettings, SamplingMetadata) {
        let parent = vec![3, 6, 7, 10];
        let definition = SamplingChannelDefinition {
            around: JOINT.to_owned(),
            parent_lmb: parent.clone(),
            on_cut: vec![1],
            radial_profile: Some(SamplingRadialProfile {
                scale: Some(1.3),
                shape: Some(2.5),
                ..Default::default()
            }),
            singularity_proxy: Some("1 + x0^2".to_owned()),
            ..Default::default()
        };
        let sampling = SamplingSettingsParser {
            graphs: SumMode::MonteCarlo,
            graph_names: vec!["B".to_owned(), "A".to_owned()],
            orientations: SumMode::MonteCarlo,
            sampling_multichanneling: true,
            sampling_channels: SumMode::MonteCarlo,
            b: 0.3,
            default_channel_selection: vec!["auto:lmb".to_owned()],
            lmb_basis_ids: ["A", "B", "C"]
                .into_iter()
                .map(|name| (name.to_owned(), vec![7]))
                .collect(),
            channel_selection: ["A", "B", "C"]
                .into_iter()
                .map(|name| {
                    (
                        name.to_owned(),
                        vec!["joint".to_owned(), "auto:lmb".to_owned()],
                    )
                })
                .collect(),
            channel_definitions: ["A", "B", "C"]
                .into_iter()
                .map(|name| {
                    (
                        name.to_owned(),
                        BTreeMap::from([("joint".to_owned(), definition.clone())]),
                    )
                })
                .collect(),
            ..Default::default()
        };
        let source: RuntimeSettings =
            toml::from_str(&toml::to_string(&BTreeMap::from([("sampling", &sampling)])).unwrap())
                .unwrap();
        let parameterization = source.sampling.get_parameterization_settings().unwrap();
        let graphs = names
            .iter()
            .enumerate()
            .map(|(index, name)| {
                let resolved =
                    resolve_sampling_channel_selection(name, &parameterization.sampling_channels)
                        .unwrap();
                let catalogue =
                    build_sampling_channel_catalogue(&resolved, &[(7, parent.clone())], &[7])
                        .unwrap();
                GraphSamplingMetadata {
                    name: (*name).to_owned(),
                    graph_id: 17 + index,
                    group_id: 4 + index,
                    graph_index: if *name == "C" { None } else { Some(index) },
                    orientations: 3,
                    parent_lmb: parent.clone(),
                    catalogue: Some(catalogue),
                }
            })
            .collect();
        (source, SamplingMetadata { sampling, graphs })
    }

    #[test]
    fn sampling_toml_preserves_named_joint_metadata_and_expands_lmb_edges() {
        let (source, metadata) = fixture(&["A"]);
        let output = SamplingMetadata::toml(&metadata.export_settings().unwrap()).unwrap();
        let document: BTreeMap<String, SamplingSettingsParser> = toml::from_str(&output).unwrap();
        let sampling = &document["sampling"];
        assert_eq!(sampling.channel_selection["A"], ["joint", "lmb(3,6,7,10)"]);
        assert_eq!(sampling.channel_selection.len(), 1);
        assert_eq!(sampling.channel_definitions.len(), 1);
        assert_eq!(
            sampling.lmb_basis_ids,
            BTreeMap::from([("A".to_owned(), vec![7])])
        );
        assert_eq!(
            sampling.channel_definitions["A"]["joint"],
            metadata.sampling.channel_definitions["A"]["joint"]
        );
        assert_eq!(sampling.graph_names, metadata.sampling.graph_names);
        assert_eq!(sampling.b, 0.3);
        assert!(!output.contains("auto:"));
        assert!(!output.contains("default_channel_selection"));
        assert!(!output.contains('\u{1b}'));
        assert!(output.contains("parent_lmb = [3, 6, 7, 10]"));
        assert_eq!(source.sampling.as_parser(), metadata.sampling);
    }

    #[test]
    fn sampling_fragments_merge_in_either_order_without_resetting_other_graphs() {
        let (source, combined) = fixture(&["A", "B"]);
        let (_, a) = fixture(&["A"]);
        let (_, b) = fixture(&["B"]);
        let a = SamplingMetadata::toml(&a.export_settings().unwrap()).unwrap();
        let b = SamplingMetadata::toml(&b.export_settings().unwrap()).unwrap();
        let together = SamplingMetadata::toml(&combined.export_settings().unwrap()).unwrap();
        let mut results = Vec::new();
        for fragments in [vec![&a, &b], vec![&b, &a], vec![&together]] {
            let mut figment = Figment::from(Serialized::defaults(&source));
            for fragment in fragments {
                figment = SetArgs::String {
                    string: fragment.to_owned(),
                }
                .merge_figment(figment)
                .unwrap();
            }
            let result: RuntimeSettings = figment.extract().unwrap();
            results.push(result);
        }
        assert_eq!(results[0], results[1]);
        assert_eq!(results[0], results[2]);
        let expected = source.sampling.as_parser();
        let actual = results[0].sampling.as_parser();
        assert_eq!(
            actual.default_channel_selection,
            expected.default_channel_selection
        );
        assert_eq!(
            actual.channel_selection["C"],
            expected.channel_selection["C"]
        );
        assert_eq!(actual.channel_definitions, expected.channel_definitions);
        assert_eq!(actual.lmb_basis_ids, expected.lmb_basis_ids);
        assert_eq!(actual.graph_names, expected.graph_names);
        assert_eq!(actual.b, expected.b);
        let parameters = results[0].sampling.get_parameterization_settings().unwrap();
        for graph in &combined.graphs {
            let resolved =
                resolve_sampling_channel_selection(&graph.name, &parameters.sampling_channels)
                    .unwrap();
            let catalogue =
                build_sampling_channel_catalogue(&resolved, &[(7, vec![3, 6, 7, 10])], &[7])
                    .unwrap();
            let original = graph.catalogue.as_ref().unwrap();
            assert_eq!(catalogue.entries.len(), original.entries.len());
            assert_eq!(catalogue.entries[0], original.entries[0]);
            assert!(
                matches!(&catalogue.entries[1], SamplingCatalogueEntry::Lmb { basis_id: 7, edges, .. } if edges == &[3,6,7,10])
            );
        }
    }

    #[test]
    fn sampling_pretty_distinguishes_ids_and_retains_joint_blocks_and_profiles() {
        let (_, metadata) = fixture(&["A", "C"]);
        let output = metadata.pretty().unwrap();
        let plain = console::strip_ansi_codes(&output);
        for expected in [
            "[graph index, orientation ID, channel ID]",
            "graph #17",
            "group #4",
            "graph index: 0",
            "basis #7",
            "expanded channels",
            "lmb(3,6,7,10)",
            "physical cut IDs: [1]",
            "profile: #0",
            "Radial profiles (display IDs for this graph)",
            "proxy: 1 + x0^2",
            "block 0:",
            "active:",
            "preceding:",
            "remaining:",
            "phase_space",
            "intersect",
            "excluded by sampling.graph_names",
            "Map recipes and cut hosts are static",
        ] {
            assert!(plain.contains(expected), "missing {expected:?}: {plain}");
        }
        assert!(
            plain
                .lines()
                .all(|line| console::measure_text_width(line) <= 145),
            "{plain}"
        );
    }

    #[test]
    fn sampling_pretty_shares_equal_profiles_and_retains_distinct_overrides() {
        let (_, mut metadata) = fixture(&["A"]);
        let catalogue = metadata.graphs[0].catalogue.as_mut().unwrap();
        let SamplingCatalogueEntry::Named(channel) = &catalogue.entries[0] else {
            panic!("fixture starts with a named channel");
        };
        let mut shared = channel.clone();
        shared.name = "shared".to_owned();
        let mut distinct = channel.clone();
        distinct.name = "distinct".to_owned();
        distinct.definition.radial_profile.as_mut().unwrap().scale = Some(4.5);
        catalogue.entries.extend([
            SamplingCatalogueEntry::Named(shared),
            SamplingCatalogueEntry::Named(distinct),
        ]);
        let output = metadata.pretty().unwrap();
        let plain = console::strip_ansi_codes(&output);
        assert_eq!(plain.matches("profile: #0").count(), 2);
        assert_eq!(plain.matches("profile: #1").count(), 1);
        assert_eq!(
            plain.matches("\"approximation\":\"log_logistic\"").count(),
            2
        );
        assert_eq!(plain.matches("\"scale\":1.3").count(), 1);
        assert_eq!(plain.matches("\"scale\":4.5").count(), 1);
        assert_eq!(plain.matches("\"shape\":2.5").count(), 2);
    }

    #[test]
    fn sampling_empty_graph_selection_stays_explicit() {
        let (_, mut metadata) = fixture(&["A"]);
        metadata.graphs[0]
            .catalogue
            .as_mut()
            .unwrap()
            .entries
            .clear();
        let output = SamplingMetadata::toml(&metadata.export_settings().unwrap()).unwrap();
        let document: BTreeMap<String, SamplingSettingsParser> = toml::from_str(&output).unwrap();
        assert_eq!(
            document["sampling"].channel_selection["A"],
            Vec::<String>::new()
        );
    }

    #[test]
    fn sampling_pretty_explains_default_first_channel_without_a_discrete_axis() {
        let (_, mut metadata) = fixture(&["A"]);
        metadata.sampling.graphs = SumMode::Summed;
        metadata.sampling.orientations = SumMode::Summed;
        metadata.sampling.sampling_multichanneling = false;
        metadata.sampling.sampling_channels = SumMode::Summed;
        metadata
            .sampling
            .channel_selection
            .insert("A".to_owned(), vec!["auto:lmb".to_owned()]);
        let catalogue = metadata.graphs[0].catalogue.as_mut().unwrap();
        catalogue
            .entries
            .retain(|entry| matches!(entry, SamplingCatalogueEntry::Lmb { .. }));
        catalogue.selectors = vec![
            gammalooprs::integrands::process::SamplingChannelSelector::Preset(
                gammalooprs::integrands::process::SamplingChannelPreset::Lmb,
            ),
        ];
        let output = metadata.pretty().unwrap();
        let plain = console::strip_ansi_codes(&output);
        assert!(plain.contains("Discrete coordinates: none"));
        assert!(plain.contains("default map: channel #0 only"));
        assert!(plain.contains("remaining catalogue entries"));
        assert!(plain.contains("weights are inactive"));
    }
}
