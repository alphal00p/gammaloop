use std::fmt::Write;

use color_eyre::Result;
use colored::Colorize;
use eyre::{eyre, Context};
use gammalooprs::{
    graph::{threshold_counterterms::ThresholdCountertermSpec, Graph, LmbIndex, LoopMomentumBasis},
    momentum::SignOrZero,
    processes::{
        ProcessCollection, ResolvedThresholdCounterterms, ThresholdCountertermOrigin,
        ThresholdCountertermSide,
    },
    utils::serde_utils::ShowDefaultsGuard,
};
use itertools::Itertools;
use linnet::{half_edge::involution::EdgeIndex, parser::GlobalData};
use tabled::{
    builder::Builder,
    settings::{object::Columns, Alignment, Style, Width},
};
use typed_index_collections::TiVec;

use crate::state::State;

use super::MetadataDisplayFormat;

pub(super) fn render(
    state: &State,
    process_id: usize,
    integrand_name: &str,
    requested_graphs: &[String],
    format: MetadataDisplayFormat,
) -> Result<String> {
    if format == MetadataDisplayFormat::Toml && requested_graphs.len() != 1 {
        return Err(eyre!(
            "Threshold TOML output requires exactly one master graph selected with --graph"
        ));
    }
    let process = state
        .process_list
        .processes
        .get(process_id)
        .ok_or_else(|| eyre!("Unknown process #{process_id}"))?;
    let integrand = process.get_integrand(integrand_name)?;
    integrand.require_generated()?;
    let graphs = match &process.collection {
        ProcessCollection::Amplitudes(amplitudes) => {
            let amplitude = &amplitudes[&integrand.canonical_name];
            amplitude
                .graph_group_structure
                .iter()
                .map(|group| {
                    let graph = &amplitude.graphs[group
                        .into_iter()
                        .next()
                        .expect("graph groups have a master")];
                    (
                        &graph.graph,
                        graph.derived_data.resolved_threshold_counterterms.as_ref(),
                        graph
                            .derived_data
                            .threshold_topology
                            .as_ref()
                            .map(|(_, lmbs)| lmbs)
                            .or(graph.derived_data.lmbs.as_ref()),
                    )
                })
                .collect_vec()
        }
        ProcessCollection::CrossSections(cross_sections) => {
            let cross_section = &cross_sections[&integrand.canonical_name];
            cross_section
                .graph_group_structure
                .iter()
                .map(|group| {
                    let graph = &cross_section.supergraphs[group
                        .into_iter()
                        .next()
                        .expect("graph groups have a master")];
                    (
                        &graph.graph,
                        graph.derived_data.resolved_threshold_counterterms.as_ref(),
                        graph.derived_data.lmbs.as_ref(),
                    )
                })
                .collect_vec()
        }
    };
    for requested in requested_graphs {
        if !graphs.iter().any(|(graph, _, _)| graph.name == *requested) {
            return Err(eyre!(
                "Unknown master graph '{}' for integrand '{}'. Available master graphs: {}",
                requested,
                integrand.canonical_name,
                graphs
                    .iter()
                    .map(|(graph, _, _)| graph.name.as_str())
                    .join(", "),
            ));
        }
    }
    graphs
        .into_iter()
        .filter(|(graph, _, _)| {
            requested_graphs.is_empty() || requested_graphs.contains(&graph.name)
        })
        .map(|(graph, resolved, lmbs)| render_graph(graph, resolved, lmbs, format))
        .collect::<Result<Vec<_>>>()
        .map(|reports| reports.join("\n\n"))
}

fn render_graph(
    graph: &Graph,
    resolved: Option<&ResolvedThresholdCounterterms>,
    lmbs: Option<&TiVec<LmbIndex, LoopMomentumBasis>>,
    format: MetadataDisplayFormat,
) -> Result<String> {
    let generation = resolved
        .map(|resolved| {
            let lmbs = lmbs.ok_or_else(|| {
                eyre!(
                    "Graph '{}' has resolved threshold metadata without its generation LMBs",
                    graph.name,
                )
            })?;
            Ok::<_, eyre::Report>((resolved, lmbs))
        })
        .transpose()?;
    let spec = generation.map_or_else(
        || graph.threshold_counterterms.value.clone(),
        |(resolved, lmbs)| resolved.materialized_spec(&graph.threshold_counterterms, lmbs),
    );
    match format {
        MetadataDisplayFormat::Toml => {
            if generation.is_none() {
                return Err(eyre!(
                    "Threshold subtraction was disabled during generation of graph '{}'. Its generation toggle cannot be represented by a threshold_counterterms DOT attribute",
                    graph.name,
                ));
            }
            render_toml(&spec)
        }
        MetadataDisplayFormat::Pretty => render_pretty(&graph.name, &spec, generation),
    }
}

fn render_toml(spec: &ThresholdCountertermSpec) -> Result<String> {
    let toml = {
        let _defaults = ShowDefaultsGuard::new(true);
        toml::to_string(spec).context("Serializing threshold counterterm settings")?
    };
    let mut attributes = GlobalData::from(());
    attributes
        .statements
        .insert("threshold_counterterms".to_string(), format!("\n{toml}"));
    Ok(format!("{attributes:width$}", width = 0))
}

fn render_pretty(
    graph_name: &str,
    spec: &ThresholdCountertermSpec,
    generation: Option<(
        &ResolvedThresholdCounterterms,
        &TiVec<LmbIndex, LoopMomentumBasis>,
    )>,
) -> Result<String> {
    let mut output = format!(
        "{}\n",
        format!("Threshold subtraction · {graph_name}")
            .bold()
            .blue()
    );
    let generated = generation.map_or(0, |(resolved, _)| {
        resolved
            .variants
            .iter()
            .flat_map(|variant| &variant.associations)
            .filter(|association| association.eligible)
            .count()
    });
    if generation.is_none() {
        writeln!(output, "{}", "Threshold subtraction was disabled during generation. Supplied declarations are shown below.".yellow())?;
    } else {
        writeln!(
            output,
            "{generated} generated threshold associations; static generation metadata"
        )?;
    }

    let edges = |edges: &[EdgeIndex]| format!("[{}]", edges.iter().map(|edge| edge.0).join(", "));
    let mut table = Builder::new();
    table.push_record(
        [
            "Cut / threshold",
            "Variant / generation",
            "Solve basis / signed cycles",
            "Multiplier",
        ]
        .map(|heading| heading.bold().blue().to_string()),
    );
    let mut row_count = 0;
    for cut in &spec.cuts {
        for threshold in &cut.thresholds {
            if threshold.counterterms.is_empty() {
                table.push_record([
                    format!(
                        "cut {}\nthreshold {}",
                        edges(&cut.edges).cyan(),
                        edges(&threshold.edges).yellow()
                    ),
                    "default\nskipped / unresolved".yellow().to_string(),
                    "automatic\nsubspace: unresolved\nparent: unresolved"
                        .dimmed()
                        .to_string(),
                    "1".dimmed().to_string(),
                ]);
                row_count += 1;
            }
            for declaration in &threshold.counterterms {
                let name = declaration.name.as_deref().unwrap_or("default");
                let resolved = generation.and_then(|(resolved, lmbs)| {
                    resolved
                        .variants
                        .iter()
                        .filter(|variant| variant.name == name)
                        .find_map(|variant| {
                            variant
                                .associations
                                .iter()
                                .find(|association| {
                                    association.cut_edges == cut.edges
                                        && association.threshold_edges == threshold.edges
                                })
                                .map(|association| (association, variant.side, lmbs))
                        })
                });
                let status = if declaration.disable {
                    "disabled".red()
                } else if let Some((association, _, _)) = resolved {
                    if association.eligible {
                        match association.origin {
                            ThresholdCountertermOrigin::Explicit => "generated · explicit".green(),
                            ThresholdCountertermOrigin::Autogenerated => {
                                "generated · default".green()
                            }
                        }
                    } else {
                        "skipped".yellow()
                    }
                } else {
                    "skipped / unresolved".yellow()
                };
                let side = resolved
                    .map(|(_, side, _)| match side {
                        ThresholdCountertermSide::Amplitude => "amplitude".cyan(),
                        ThresholdCountertermSide::Left => "left".green(),
                        ThresholdCountertermSide::Right => "right".yellow(),
                    })
                    .unwrap_or_else(|| "unresolved".dimmed());
                let mut solve = format!(
                    "{}\nsubspace: {}\nparent: {}",
                    declaration
                        .group_id
                        .map(|id| format!("group_id = {id}").yellow().to_string())
                        .unwrap_or_else(|| "automatic group".dimmed().to_string()),
                    declaration
                        .subspace
                        .as_deref()
                        .map(|value| edges(value).yellow().to_string())
                        .unwrap_or_else(|| "maximal (unresolved)".dimmed().to_string()),
                    declaration
                        .parent_lmb
                        .as_deref()
                        .map(|value| edges(value).cyan().to_string())
                        .unwrap_or_else(|| "unresolved".dimmed().to_string()),
                );
                if let Some((association, _, lmbs)) = resolved {
                    for (edge, cycle) in association.subspace.solve_signature(lmbs) {
                        let cycle = cycle
                            .iter()
                            .map(|(edge, sign)| {
                                let sign = match sign {
                                    SignOrZero::Plus => "+".green(),
                                    SignOrZero::Minus => "-".red(),
                                    SignOrZero::Zero => "0".dimmed(),
                                };
                                format!("{sign}e{}", edge.0)
                            })
                            .join(" ");
                        write!(
                            solve,
                            "\n{} → {cycle}",
                            format!("e{}", edge.0).bold().cyan()
                        )?;
                    }
                }
                let multiplier = if let Some(multiplier) = &declaration.multiplier {
                    let mut text = format!(
                        "{}\nsymmetrize = {}\nopaque_derivatives = {}",
                        multiplier.expression.magenta(),
                        multiplier.symmetrize.to_string().yellow(),
                        multiplier.opaque_derivatives.to_string().yellow(),
                    );
                    for (name, definition) in &multiplier.function_map {
                        write!(
                            text,
                            "\nlocal {} = {}",
                            name.bold().cyan(),
                            definition.magenta()
                        )?;
                    }
                    text
                } else {
                    "1".dimmed().to_string()
                };
                table.push_record([
                    format!(
                        "cut {}\nthreshold {}\n{side}",
                        edges(&cut.edges).cyan(),
                        edges(&threshold.edges).yellow()
                    ),
                    format!("{}\n{status}", name.bold().cyan()),
                    solve,
                    multiplier,
                ]);
                row_count += 1;
            }
        }
    }
    if row_count == 0 {
        writeln!(output, "No threshold declarations.")?;
    } else {
        let mut table = table.build();
        table
            .with(Style::rounded())
            .with(Alignment::left())
            .with(Alignment::top());
        for (column, width) in [26, 24, 38, 48].into_iter().enumerate() {
            table.modify(Columns::one(column), Width::wrap(width).keep_words(true));
        }
        writeln!(output, "{table}")?;
        writeln!(
            output,
            "{} {} → {}e3 {}e14 means varying basis momentum e3 adds to edge 3 and subtracts from edge 14.",
            "Cycles:".bold().blue(),
            "e3".bold().cyan(),
            "+".green(),
            "-".red()
        )?;
        writeln!(output, "{}", "Signs follow edge orientation. Selected edges + signed cycles define automatic solve groups; group_id can separate them.".dimmed())?;
    }
    if !spec.function_map.is_empty() {
        writeln!(
            output,
            "{}",
            "Shared function definitions (local definitions override these)"
                .bold()
                .blue()
        )?;
        let mut definitions = Builder::new();
        definitions.push_record(
            ["Function", "Definition"].map(|heading| heading.bold().blue().to_string()),
        );
        for (name, definition) in &spec.function_map {
            definitions.push_record([
                name.bold().cyan().to_string(),
                definition.magenta().to_string(),
            ]);
        }
        let mut definitions = definitions.build();
        definitions
            .with(Style::rounded())
            .with(Alignment::left())
            .modify(Columns::one(0), Width::wrap(32).keep_words(true))
            .modify(Columns::one(1), Width::wrap(100).keep_words(true));
        writeln!(output, "{definitions}")?;
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use linnet::parser::DotGraph;

    #[test]
    fn gl638_dot_metadata_roundtrips_with_verbose_flags_and_function_scopes() {
        let graph: DotGraph = DotGraph::from_string(include_str!(
            "../../../../../examples/cli/epem_a_ttxh/NNLO/graphs/GL638.dot"
        ))
        .unwrap();
        let mut spec: ThresholdCountertermSpec =
            toml::from_str(&graph.global_data.statements["threshold_counterterms"]).unwrap();
        let variant = &mut spec.cuts[0].thresholds[0].counterterms[0];
        variant.name = Some("quoted \"variant\" λ\\line\nend\\".into());
        variant
            .multiplier
            .as_mut()
            .unwrap()
            .function_map
            .insert("wQ()".into(), "local_override()".into());
        let attribute = render_toml(&spec).unwrap();
        assert!(attribute.starts_with("threshold_counterterms = \"\n"));
        assert!(!attribute.contains('\u{1b}'));
        assert!(attribute.contains("disable = false"));
        assert!(attribute.contains("symmetrize = false"));
        assert!(attribute.contains("opaque_derivatives = true"));
        assert!(attribute.contains("parent_lmb = [3, 6, 7, 10]"));
        let reparsed: DotGraph =
            DotGraph::from_string(format!("digraph G {{\n{attribute}\n}}")).unwrap();
        let roundtrip: ThresholdCountertermSpec =
            toml::from_str(&reparsed.global_data.statements["threshold_counterterms"]).unwrap();
        assert_eq!(roundtrip, spec);
        assert_eq!(roundtrip.function_map["wQ()"], "P(0,cind(0))+P(1,cind(0))");
    }

    #[test]
    fn threshold_pretty_preserves_gl638_definitions_and_dormant_variants() {
        let graph: DotGraph = DotGraph::from_string(include_str!(
            "../../../../../examples/cli/epem_a_ttxh/NNLO/graphs/GL638.dot"
        ))
        .unwrap();
        let spec: ThresholdCountertermSpec =
            toml::from_str(&graph.global_data.statements["threshold_counterterms"]).unwrap();
        let report = render_pretty("GL638", &spec, None).unwrap();
        let report = console::strip_ansi_codes(&report);
        assert!(report.contains("GL638"));
        assert!(report.contains("disabled during generation"));
        assert!(report.contains("shared_1l"));
        assert!(report.contains("native_2l"));
        assert!(report.contains("parent: [3, 6, 7, 10]"));
        assert!(report.contains("skipped / unresolved"));
        assert!(report.contains("WHc()"));
        assert!(report.contains("Solve basis / signed cycles"));
        assert!(report.contains("varying basis momentum e3"));
        assert!(report.contains("Signs follow edge orientation"));
        assert_eq!(report.matches("Shared function definitions").count(), 1);
        assert!(
            report
                .lines()
                .all(|line| console::measure_text_width(line) <= 150),
            "{report}"
        );
    }

    #[test]
    fn empty_enabled_thresholds_are_distinct_from_disabled_generation() -> Result<()> {
        use gammalooprs::{initialisation::test_initialise, utils::load_generic_model};

        test_initialise()?;
        let spec = ThresholdCountertermSpec::default();
        let resolved = ResolvedThresholdCounterterms {
            legacy_equivalent: true,
            variants: TiVec::new(),
            cross_section_cut_groups: TiVec::new(),
        };
        let lmbs = TiVec::new();
        let enabled = render_pretty("empty", &spec, Some((&resolved, &lmbs))).unwrap();
        assert!(enabled.contains("0 generated threshold associations"));
        assert!(!enabled.contains("disabled during generation"));
        let disabled = render_pretty("empty", &spec, None).unwrap();
        assert!(disabled.contains("disabled during generation"));
        let attribute = render_toml(&spec).unwrap();
        assert!(attribute.contains("cuts = []"));

        let model = load_generic_model("scalars");
        let graph = Graph::from_string(
            include_str!("../../../../../tests/resources/graphs/scalar_bubble.dot"),
            &model,
        )?
        .remove(0);
        let lmbs = TiVec::from(vec![graph.loop_momentum_basis.clone()]);
        let before = graph.threshold_counterterms.clone();
        let enabled = render_graph(
            &graph,
            Some(&resolved),
            Some(&lmbs),
            MetadataDisplayFormat::Toml,
        )?;
        assert_eq!(enabled, attribute);
        let disabled =
            render_graph(&graph, None, Some(&lmbs), MetadataDisplayFormat::Toml).unwrap_err();
        assert!(disabled.to_string().contains("disabled during generation"));
        let missing_lmbs =
            render_graph(&graph, Some(&resolved), None, MetadataDisplayFormat::Toml).unwrap_err();
        assert!(missing_lmbs
            .to_string()
            .contains("without its generation LMBs"));
        assert_eq!(graph.threshold_counterterms, before);
        Ok(())
    }
}
