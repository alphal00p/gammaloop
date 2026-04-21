use std::{
    collections::{BTreeMap, HashSet},
    fmt::Display,
};

use color_eyre::eyre::Result;
use colored::Colorize;
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, Orientation};
use rand::Rng;
use symbolica::numerical_integration::MonteCarloRng;
use tabled::{builder::Builder, settings::Style};
use tracing::warn;
use typed_index_collections::TiVec;

use crate::{
    DependentMomentaConstructor,
    cff::esurface::{EsurfaceID, ExistingEsurfaceId, ExistingEsurfaces, GroupEsurfaceId},
    graph::{FeynmanGraph, GraphGroupPosition, LmbError, lmb::LMBwithEdges},
    integrands::{
        evaluation::PreciseEvaluationResult,
        process::{
            GraphTerm, OrientationProfileMode, ProcessIntegrandImpl,
            amplitude::{AmplitudeGraphTerm, AmplitudeIntegrand},
            cross_section::{CrossSectionGraphTerm, CrossSectionIntegrand},
            evaluate_profile_momentum_point_precise, orientation_labels_for_graph,
        },
    },
    model::Model,
    momentum::{
        ThreeMomentum,
        sample::{LoopIndex, LoopMomenta, MomentumSample},
    },
    observables::events::AdditionalWeightKey,
    settings::{
        RuntimeSettings, SamplingSettings,
        runtime::{
            DiscreteGraphSamplingSettings, DiscreteGraphSamplingType, ParameterizationMapping,
            ParameterizationMode, ParameterizationSettings,
        },
    },
    subtraction::amplitude_counterterm::OverlapStructureWithKinematics,
    utils::{
        ArbPrec, F, FloatLike, box_muller,
        fitting::{constant_dropped_fit_points, log_log_slope_constant_dropped},
    },
};

/// The range is from 10^start to 10^end.
pub struct IRProfileSetting {
    pub lambda_exp_start: f64,
    pub lambda_exp_end: f64,
    pub steps: usize,
    pub seed: u64,
    pub select_limits_and_graphs: Option<String>,
    pub orientation_mode: OrientationProfileMode,
    pub show_per_cut_info: bool,
}

impl AmplitudeGraphTerm {
    fn enumerate_ir_limits(&self) -> Vec<IrLimit> {
        let mut limits: HashSet<IrLimit> = HashSet::new();

        let massless_edges: Vec<EdgeIndex> = self
            .graph
            .iter_edges_of(&!self.graph.tree_edges.clone())
            .filter_map(|(_a, b, c)| {
                if c.data.particle.is_massless() {
                    Some(b)
                } else {
                    None
                }
            })
            .collect();

        for subset in massless_edges.iter().powerset() {
            if subset.is_empty() {
                continue;
            }
            let _lmb = match self.lmb_with_loop_edges(subset.as_slice()) {
                Ok(lmb) => lmb,
                Err(err) => match err {
                    LmbError::NotLoopEdges { .. } => {
                        // warn!("{loop_edges} is not a valid loop edge subset");
                        continue;
                    }
                    a => panic!("Failed to build IR loop momentum basis for cut graph:\n{a}"),
                },
            };

            let ir_limit = IrLimit::new_pure_soft(subset.into_iter().copied().collect());

            limits.insert(ir_limit);
        }

        limits.into_iter().sorted().collect()
    }
}

impl CrossSectionGraphTerm {
    fn enumerate_ir_limits(&self) -> Vec<IrLimit> {
        let mut limits: HashSet<IrLimit> = HashSet::new();
        let loop_count = self.graph.loop_momentum_basis.loop_edges.len();

        for cut_group in self.raised_data.raised_cut_groups.iter() {
            let mut limits_of_cut: HashSet<IrLimit> = HashSet::new();

            let representative_cut_esurface = &self.cut_esurface[*cut_group.cuts.first().unwrap()];
            let massless_edges_in_cut = representative_cut_esurface
                .energies
                .iter()
                .filter(|edge_id| self.graph[**edge_id].particle.is_massless())
                .copied()
                .collect_vec();

            if massless_edges_in_cut.len() >= 2 {
                let subsets = massless_edges_in_cut
                    .iter()
                    .powerset()
                    .filter(|subset| {
                        subset.len() >= 2
                            && subset.len() <= loop_count
                            && subset.len() < representative_cut_esurface.energies.len()
                    })
                    .collect_vec();

                for subset in subsets {
                    let ir_limit =
                        IrLimit::new_pure_colinear(subset.into_iter().copied().collect());
                    limits_of_cut.insert(ir_limit);
                }
            }

            if !massless_edges_in_cut.is_empty() {
                let subsets = massless_edges_in_cut
                    .iter()
                    .powerset()
                    .filter(|subset| {
                        !subset.is_empty()
                            && subset.len() <= loop_count
                            && subset.len() < representative_cut_esurface.energies.len()
                    })
                    .collect_vec();
                for subset in subsets {
                    let ir_limit = IrLimit::new_pure_soft(subset.into_iter().copied().collect());
                    limits_of_cut.insert(ir_limit);
                }
            }

            for limit in limits_of_cut.drain() {
                limits.insert(limit);
            }
        }

        limits.into_iter().sorted().collect()
    }
}

pub struct IrLimitTestReport {
    pub all_passed: bool,
    pub results_per_graph: Vec<GraphIRLimitReport>,
}

pub struct GraphIRLimitReport {
    pub graph_name: String,
    pub all_limits_passed: bool,
    pub cut_definitions: Vec<GraphCutDefinition>,
    pub single_limit_reports: Vec<SingleLimitReport>,
}

#[derive(Debug, Clone)]
pub struct GraphCutDefinition {
    pub cut_id: usize,
    pub edges: Vec<EdgeIndex>,
}

pub struct SingleLimitReport {
    pub limit_name: String,
    pub orientation_label: Option<String>,
    pub passed: bool,
    pub power_law_fit: PowerLawFit,
    pub scaling: f64,
    pub per_cut_reports: Vec<CutLimitReport>,
    pub display_only_reports: Vec<DisplayOnlyLimitReport>,
    num_soft: usize,
}

pub struct CutLimitReport {
    pub cut_id: usize,
    pub power_law_fit: Option<PowerLawFit>,
    pub scaling: Option<f64>,
    pub fit_error: Option<String>,
}

pub struct DisplayOnlyLimitReport {
    pub label: String,
    pub power_law_fit: Option<PowerLawFit>,
    pub scaling: Option<f64>,
    pub fit_error: Option<String>,
}

impl Display for IrLimitTestReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let overall_status = if self.all_passed {
            "PASS".green().bold()
        } else {
            "FAIL".red().bold()
        };

        let passed_graphs = self
            .results_per_graph
            .iter()
            .filter(|graph_report| graph_report.all_limits_passed)
            .count();

        writeln!(
            f,
            "IR limit tests: {} ({}/{})",
            overall_status,
            passed_graphs,
            self.results_per_graph.len()
        )?;

        let mut graph_summary_table = Builder::new();
        graph_summary_table.push_record(["graph", "status", "passed", "total"]);

        for graph_report in &self.results_per_graph {
            let graph_status = if graph_report.all_limits_passed {
                "PASS".green().bold().to_string()
            } else {
                "FAIL".red().bold().to_string()
            };

            let passed_limits = graph_report
                .single_limit_reports
                .iter()
                .filter(|report| report.passed)
                .count();

            graph_summary_table.push_record([
                graph_report.graph_name.clone(),
                graph_status,
                passed_limits.to_string(),
                graph_report.single_limit_reports.len().to_string(),
            ]);
        }

        writeln!(f, "{}", graph_summary_table.build().with(Style::rounded()))?;

        for graph_report in &self.results_per_graph {
            writeln!(f)?;
            writeln!(f, "{graph_report}")?;
        }

        Ok(())
    }
}

impl Display for GraphIRLimitReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let graph_status = if self.all_limits_passed {
            "PASS".green().bold()
        } else {
            "FAIL".red().bold()
        };

        let passed_limits = self
            .single_limit_reports
            .iter()
            .filter(|report| report.passed)
            .count();

        writeln!(
            f,
            "  {} {} ({}/{})",
            graph_status,
            self.graph_name.bold(),
            passed_limits,
            self.single_limit_reports.len()
        )?;

        if !self.cut_definitions.is_empty() {
            render_graph_cut_definitions(f, &self.cut_definitions)?;
            writeln!(f)?;
        }

        let mut limit_table = Builder::new();
        let mut separators_after_data_rows = Vec::new();
        let mut data_row_count = 0;
        limit_table.push_record([
            "status",
            "limit",
            "orientation",
            "item",
            "scaling",
            "p",
            "r_squared",
            "n_soft",
        ]);

        for (report_index, report) in self.single_limit_reports.iter().enumerate() {
            let status = if report.passed {
                "PASS".green().bold().to_string()
            } else {
                "FAIL".red().bold().to_string()
            };

            limit_table.push_record([
                status,
                report.limit_name.clone(),
                report
                    .orientation_label
                    .clone()
                    .unwrap_or_else(|| "sum".to_string()),
                "sum".to_string(),
                format!("{:+.4}", report.scaling),
                format!("{:+.4}", report.power_law_fit.exponent),
                format!("{:.4}", report.power_law_fit.r_squared),
                report.num_soft.to_string(),
            ]);
            data_row_count += 1;

            if !report.per_cut_reports.is_empty() || !report.display_only_reports.is_empty() {
                separators_after_data_rows.push(data_row_count);
            }

            for cut_report in &report.per_cut_reports {
                let [item, scaling, exponent, r_squared] = cut_report_display_row(cut_report);
                limit_table.push_record([
                    "INFO".cyan().bold().to_string(),
                    String::new(),
                    String::new(),
                    item,
                    scaling,
                    exponent,
                    r_squared,
                    String::new(),
                ]);
                data_row_count += 1;
            }

            for display_only_report in &report.display_only_reports {
                let [item, scaling, exponent, r_squared] =
                    display_only_report_display_row(display_only_report);
                limit_table.push_record([
                    "INFO".cyan().bold().to_string(),
                    String::new(),
                    String::new(),
                    item,
                    scaling,
                    exponent,
                    r_squared,
                    String::new(),
                ]);
                data_row_count += 1;
            }

            if report_index + 1 < self.single_limit_reports.len() {
                separators_after_data_rows.push(data_row_count);
            }
        }

        let mut table = limit_table.build();
        table.with(Style::rounded());
        let rendered =
            insert_limit_table_separators(table.to_string(), &separators_after_data_rows);

        write!(f, "{rendered}")?;

        Ok(())
    }
}

impl Display for SingleLimitReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let status = if self.passed {
            "PASS".green().bold()
        } else {
            "FAIL".red().bold()
        };

        write!(
            f,
            "{} {}{} | scaling={:+.4} | p={:+.4} | R²={:.4} | n_soft={}",
            status,
            self.limit_name,
            self.orientation_label
                .as_ref()
                .map(|label| format!(" @ {label}"))
                .unwrap_or_default(),
            self.scaling,
            self.power_law_fit.exponent,
            self.power_law_fit.r_squared,
            self.num_soft
        )?;

        render_per_cut_reports(f, self, "")?;
        render_display_only_reports(f, self, "")
    }
}

fn render_per_cut_reports(
    f: &mut std::fmt::Formatter<'_>,
    report: &SingleLimitReport,
    indent: &str,
) -> std::fmt::Result {
    if report.per_cut_reports.is_empty() {
        return Ok(());
    }

    writeln!(f)?;
    writeln!(
        f,
        "\n{indent}per-cut fits for {}{}",
        report.limit_name,
        report
            .orientation_label
            .as_ref()
            .map(|label| format!(" @ {label}"))
            .unwrap_or_default(),
    )?;

    let mut cut_table = Builder::new();
    cut_table.push_record(["cut", "scaling", "p", "r_squared"]);

    for cut_report in &report.per_cut_reports {
        cut_table.push_record(cut_report_display_row(cut_report));
    }

    write!(f, "{}", cut_table.build().with(Style::rounded()))
}

fn render_graph_cut_definitions(
    f: &mut std::fmt::Formatter<'_>,
    cut_definitions: &[GraphCutDefinition],
) -> std::fmt::Result {
    writeln!(f, "  cut definitions")?;

    let mut cut_table = Builder::new();
    cut_table.push_record(["cut", "edges"]);

    for cut_definition in cut_definitions {
        cut_table.push_record([
            cut_definition.cut_id.to_string(),
            cut_definition
                .edges
                .iter()
                .map(ToString::to_string)
                .join(", "),
        ]);
    }

    writeln!(f, "{}", cut_table.build().with(Style::rounded()))
}

fn render_display_only_reports(
    f: &mut std::fmt::Formatter<'_>,
    report: &SingleLimitReport,
    indent: &str,
) -> std::fmt::Result {
    if report.display_only_reports.is_empty() {
        return Ok(());
    }

    writeln!(f)?;
    writeln!(
        f,
        "\n{indent}display-only fits for {}{}",
        report.limit_name,
        report
            .orientation_label
            .as_ref()
            .map(|label| format!(" @ {label}"))
            .unwrap_or_default(),
    )?;

    let mut display_only_table = Builder::new();
    display_only_table.push_record(["component", "scaling", "p", "r_squared"]);

    for display_only_report in &report.display_only_reports {
        display_only_table.push_record(display_only_report_display_row(display_only_report));
    }

    write!(f, "{}", display_only_table.build().with(Style::rounded()))
}

fn cut_report_display_row(cut_report: &CutLimitReport) -> [String; 4] {
    let (scaling, exponent, r_squared) = match (&cut_report.power_law_fit, &cut_report.fit_error) {
        (Some(fit), None) => (
            format!("{:+.4}", cut_report.scaling.unwrap_or(fit.exponent)),
            format!("{:+.4}", fit.exponent),
            format!("{:.4}", fit.r_squared),
        ),
        (None, Some(_)) => ("-".to_string(), "-".to_string(), "-".to_string()),
        _ => ("-".to_string(), "-".to_string(), "-".to_string()),
    };

    [
        format!("cut {}", cut_report.cut_id),
        scaling,
        exponent,
        r_squared,
    ]
}

fn display_only_report_display_row(display_only_report: &DisplayOnlyLimitReport) -> [String; 4] {
    let (scaling, exponent, r_squared) = match (
        &display_only_report.power_law_fit,
        &display_only_report.fit_error,
    ) {
        (Some(fit), None) => (
            format!(
                "{:+.4}",
                display_only_report.scaling.unwrap_or(fit.exponent)
            ),
            format!("{:+.4}", fit.exponent),
            format!("{:.4}", fit.r_squared),
        ),
        (None, Some(_)) => ("-".to_string(), "-".to_string(), "-".to_string()),
        _ => ("-".to_string(), "-".to_string(), "-".to_string()),
    };

    [
        display_only_report.label.clone(),
        scaling,
        exponent,
        r_squared,
    ]
}

fn insert_limit_table_separators(rendered: String, separators_after_data_rows: &[usize]) -> String {
    if separators_after_data_rows.is_empty() {
        return rendered;
    }

    let lines = rendered.lines().collect_vec();
    if lines.len() < 4 {
        return rendered;
    }

    let line_count = lines.len();
    let separator_line = lines[2].to_string();
    let mut next_separator = separators_after_data_rows.iter().copied().peekable();
    let mut output = Vec::with_capacity(lines.len() + separators_after_data_rows.len());
    let mut seen_data_rows = 0;

    for (line_index, line) in lines.into_iter().enumerate() {
        output.push(line.to_string());

        if line_index >= 3 && line_index + 1 < line_count {
            seen_data_rows += 1;
            while next_separator.peek().copied() == Some(seen_data_rows) {
                output.push(separator_line.clone());
                next_separator.next();
            }
        }
    }

    output.join("\n")
}

fn ir_profile_completion_entries(
    limits: Vec<(String, Vec<ProfileLimit>)>,
) -> Vec<(String, Vec<String>)> {
    limits
        .into_iter()
        .map(|(graph_name, limits)| {
            (
                graph_name,
                limits.into_iter().map(|limit| limit.to_string()).collect(),
            )
        })
        .collect()
}

fn graph_cut_definitions_for_cross_section_term(
    graph_term: &CrossSectionGraphTerm,
) -> Vec<GraphCutDefinition> {
    graph_term
        .cuts
        .iter_enumerated()
        .map(|(cut_id, cut)| GraphCutDefinition {
            cut_id: cut_id.into(),
            edges: graph_term
                .graph
                .underlying
                .iter_edges_of(&cut.cut)
                .map(|(_, edge_id, _)| edge_id)
                .sorted()
                .collect(),
        })
        .collect()
}

fn graph_id_by_name<I: ProcessIntegrandImpl>(integrand: &I, graph_name: &str) -> Option<usize> {
    (0..integrand.graph_count())
        .find(|graph_id| integrand.get_graph(*graph_id).name() == graph_name)
}

fn parse_select_limits_and_graphs<I: ProcessIntegrandImpl>(
    integrand: &I,
    input: &str,
) -> Result<Vec<(String, Vec<ProfileLimit>)>> {
    input
        .split(';')
        .map(|graph_info_string| {
            let mut parts = graph_info_string.split(' ');
            let graph_name = parts
                .next()
                .ok_or_else(|| eyre!("Expected graph name in select_limits_and_graphs"))?
                .to_string();

            let Some(graph_id) = graph_id_by_name(integrand, &graph_name) else {
                return Err(eyre!(
                    "Graph name '{}' in select_limits_and_graphs does not match any graph in the integrand",
                    graph_name
                ));
            };

            let profile_limits = parts
                .map(ProfileLimit::parse_limit)
                .collect::<Result<Vec<_>, _>>()?;

            if profile_limits.is_empty() {
                return Err(eyre!(
                    "No limits specified for graph '{}' in select_limits_and_graphs",
                    graph_name
                ));
            }

            let loop_number = integrand
                .get_graph(graph_id)
                .get_graph()
                .loop_momentum_basis
                .loop_edges
                .len();
            if profile_limits
                .iter()
                .any(|limit| !limit.is_valid(loop_number).is_ok())
            {
                return Err(eyre!(
                    "One or more limits specified for graph '{}' in select_limits_and_graphs are not valid",
                    graph_name
                ));
            }

            Ok((graph_name, profile_limits))
        })
        .collect::<Result<Vec<_>, _>>()
}

fn requested_orientations<I: ProcessIntegrandImpl>(
    integrand: &I,
    graph_id: usize,
    settings: &IRProfileSetting,
) -> Result<Vec<(Option<usize>, Option<String>)>> {
    if settings.orientation_mode.profiles_per_orientation() {
        Ok(orientation_labels_for_graph(integrand, graph_id)?
            .into_iter()
            .enumerate()
            .map(|(orientation_id, label)| (Some(orientation_id), Some(label)))
            .collect())
    } else {
        Ok(vec![(None, None)])
    }
}

fn build_single_limit_report(
    ir_limit: &IrLimit,
    orientation_label: Option<String>,
    slope: PowerLawFit,
    per_cut_reports: Vec<CutLimitReport>,
) -> SingleLimitReport {
    let num_soft = ir_limit.num_soft();
    let scaling = slope.exponent + ((num_soft * 3) as f64);
    SingleLimitReport {
        limit_name: format!("{}", ir_limit),
        orientation_label,
        passed: scaling > 0.0,
        power_law_fit: slope,
        scaling,
        per_cut_reports,
        display_only_reports: Vec::new(),
        num_soft,
    }
}

fn build_threshold_limit_report(
    threshold_limit: &ThresholdLimit,
    orientation_label: Option<String>,
    slope: PowerLawFit,
    per_cut_reports: Vec<CutLimitReport>,
) -> SingleLimitReport {
    let scaling = slope.exponent;
    SingleLimitReport {
        limit_name: format!("{}", threshold_limit),
        orientation_label,
        passed: scaling > 0.0,
        power_law_fit: slope,
        scaling,
        per_cut_reports,
        display_only_reports: Vec::new(),
        num_soft: 0,
    }
}

fn build_cut_limit_reports(
    num_soft: usize,
    cut_fits: Vec<(usize, Result<PowerLawFit>)>,
) -> Vec<CutLimitReport> {
    cut_fits
        .into_iter()
        .map(|(cut_id, fit)| match fit {
            Ok(power_law_fit) => CutLimitReport {
                cut_id,
                scaling: Some(power_law_fit.exponent + ((num_soft * 3) as f64)),
                power_law_fit: Some(power_law_fit),
                fit_error: None,
            },
            Err(error) => CutLimitReport {
                cut_id,
                scaling: None,
                power_law_fit: None,
                fit_error: Some(error.to_string()),
            },
        })
        .collect()
}

fn build_display_only_limit_reports(
    component_fits: Vec<(AdditionalWeightKey, Result<PowerLawFit>)>,
) -> Vec<DisplayOnlyLimitReport> {
    component_fits
        .into_iter()
        .map(|(key, fit)| match fit {
            Ok(power_law_fit) => DisplayOnlyLimitReport {
                label: display_only_limit_label(key),
                scaling: Some(power_law_fit.exponent),
                power_law_fit: Some(power_law_fit),
                fit_error: None,
            },
            Err(error) => DisplayOnlyLimitReport {
                label: display_only_limit_label(key),
                scaling: None,
                power_law_fit: None,
                fit_error: Some(error.to_string()),
            },
        })
        .collect()
}

fn display_only_limit_label(key: AdditionalWeightKey) -> String {
    match key {
        AdditionalWeightKey::Original => "original".to_string(),
        AdditionalWeightKey::ThresholdCounterterm { subset_index } => {
            format!("ct_{subset_index}")
        }
        AdditionalWeightKey::FullMultiplicativeFactor => "full multiplicative factor".to_string(),
    }
}

fn threshold_approach_loop_momenta<T: FloatLike>(
    overlap_group_center: &LoopMomenta<F<f64>>,
    threshold_point: &MomentumSample<T>,
    lambda: &F<T>,
) -> LoopMomenta<F<T>> {
    let overlap_group_center = overlap_group_center.cast::<T>();
    let threshold_loop_momenta = threshold_point.loop_moms();
    let offset_towards_center =
        (&overlap_group_center - threshold_loop_momenta).rescale(lambda, None);

    threshold_loop_momenta + &offset_towards_center
}

fn run_ir_profile<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    ir_profile_settings: &IRProfileSetting,
    model: &Model,
    enumerate_limits: impl Fn(&I) -> Vec<(String, Vec<ProfileLimit>)>,
    graph_cut_definitions: impl Fn(&I, usize) -> Vec<GraphCutDefinition>,
    points_on_threshold: &[OverlapStructureWithKinematics<ArbPrec>],
    mut test_single_limit: impl FnMut(
        &mut I,
        usize,
        &ProfileLimit,
        &mut MonteCarloRng,
        &IRProfileSetting,
        &Model,
        &[OverlapStructureWithKinematics<ArbPrec>],
    ) -> Result<Vec<SingleLimitReport>>,
) -> Result<IrLimitTestReport> {
    let mut rng = MonteCarloRng::new(ir_profile_settings.seed, 0);
    let limits_to_check =
        if let Some(select_limits_and_graphs) = &ir_profile_settings.select_limits_and_graphs {
            parse_select_limits_and_graphs(integrand, select_limits_and_graphs)?
        } else {
            enumerate_limits(integrand)
        };

    let mut result = IrLimitTestReport {
        all_passed: false,
        results_per_graph: Vec::new(),
    };

    for (graph_name, limits) in limits_to_check {
        let Some(graph_id) = graph_id_by_name(integrand, &graph_name) else {
            return Err(eyre!("Graph name '{}' not found in integrand", graph_name));
        };

        let mut graph_report = GraphIRLimitReport {
            graph_name: graph_name.clone(),
            all_limits_passed: false,
            cut_definitions: graph_cut_definitions(integrand, graph_id),
            single_limit_reports: Vec::new(),
        };

        for limit in limits {
            graph_report.single_limit_reports.extend(test_single_limit(
                integrand,
                graph_id,
                &limit,
                &mut rng,
                ir_profile_settings,
                model,
                points_on_threshold,
            )?);
        }

        graph_report.all_limits_passed = graph_report
            .single_limit_reports
            .iter()
            .all(|report| report.passed);

        result.results_per_graph.push(graph_report);
    }

    result.all_passed = result
        .results_per_graph
        .iter()
        .all(|graph_report| graph_report.all_limits_passed);

    Ok(result)
}

impl AmplitudeIntegrand {
    pub fn ir_profile_completion_entries(&self) -> Vec<(String, Vec<String>)> {
        ir_profile_completion_entries(self.enumerate_ir_limits())
    }

    pub fn test_ir(
        &mut self,
        ir_profile_settings: &IRProfileSetting,
        model: &Model,
    ) -> Result<IrLimitTestReport> {
        // override the sampling to be in momentum space
        self.settings.sampling = SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings {
            sample_orientations: false,
            sampling_type: DiscreteGraphSamplingType::Default(ParameterizationSettings {
                mode: ParameterizationMode::MomentumSpace,
                mapping: ParameterizationMapping::default(),
                b: 10.0,
            }),
        });

        let previous_generate_events = self.settings.general.generate_events;
        let previous_store_additional_weights =
            self.settings.general.store_additional_weights_in_event;
        self.settings.general.generate_events = true;
        self.settings.general.store_additional_weights_in_event = true;

        let result = (|| {
            self.warm_up(model)?;

            let dependent_momenta_constructor = DependentMomentaConstructor::Amplitude(
                &self.data.graph_terms[0].graph.get_external_signature(),
            );

            let mut rng = MonteCarloRng::new(ir_profile_settings.seed, 0);
            let random_loop_momenta = (0..self.data.graph_terms[0]
                .graph
                .loop_momentum_basis
                .loop_edges
                .len())
                .map(|_| sample_random_unit_vector(&mut rng))
                .collect();

            let momentum_sample = MomentumSample::new(
                random_loop_momenta,
                0,
                &self.settings.kinematics.externals,
                0,
                F::from_f64(1.0),
                dependent_momenta_constructor,
                None,
            )?;

            let points_on_threshold =
                self.kinematics_for_threshold_approach(&momentum_sample, model)?;

            run_ir_profile(
                self,
                ir_profile_settings,
                model,
                Self::enumerate_ir_limits,
                Self::graph_cut_definitions,
                &points_on_threshold,
                Self::test_single_ir_limit_impl,
            )
        })();

        if self.settings.general.generate_events != previous_generate_events {
            self.settings.general.generate_events = previous_generate_events;
        }
        if self.settings.general.store_additional_weights_in_event
            != previous_store_additional_weights
        {
            self.settings.general.store_additional_weights_in_event =
                previous_store_additional_weights;
        }

        result
    }

    fn enumerate_ir_limits(&self) -> Vec<(String, Vec<ProfileLimit>)> {
        self.data
            .graph_terms
            .iter()
            .map(|term| {
                let graph_name = term.graph.name.clone();
                let mut limits = term
                    .enumerate_ir_limits()
                    .into_iter()
                    .map(ProfileLimit::Ir)
                    .collect_vec();
                limits.extend(
                    ThresholdLimit::enumerate_from_overlap_structure(
                        &term.threshold_counterterm.overlap.existing_esurfaces,
                        &term.threshold_counterterm.esurface_map,
                        term.threshold_counterterm.own_group_position,
                    )
                    .into_iter()
                    .map(ProfileLimit::Threshold),
                );
                limits.sort();
                limits.dedup();
                (graph_name, limits)
            })
            .collect()
    }

    fn graph_cut_definitions(&self, _graph_id: usize) -> Vec<GraphCutDefinition> {
        Vec::new()
    }

    fn test_single_ir_limit_impl(
        &mut self,
        graph_id: usize,
        profile_limit: &ProfileLimit,
        rng: &mut MonteCarloRng,
        approach_settings: &IRProfileSetting,
        model: &Model,
        points_on_threshold: &[OverlapStructureWithKinematics<ArbPrec>],
    ) -> Result<Vec<SingleLimitReport>> {
        match profile_limit {
            ProfileLimit::Ir(ir_limit) => {
                let edges_in_limit = ir_limit.get_all_edges()?;
                let lmb = self.data.graph_terms[graph_id]
                    .lmb_with_loop_edges(edges_in_limit.as_slice())?;
                let momenta = ir_limit.get_momenta(rng, &self.settings, approach_settings)?;
                let non_limit_loops = lmb
                    .loop_edges
                    .iter_enumerated()
                    .filter_map(|(loop_id, edge_id)| {
                        if !edges_in_limit.contains(edge_id) {
                            Some(loop_id)
                        } else {
                            None
                        }
                    })
                    .collect_vec();

                let non_limit_momenta = non_limit_loops
                    .iter()
                    .map(|loop_id| (*loop_id, sample_random_unit_vector(rng)))
                    .collect_vec();

                let externals = self.data.graph_terms[graph_id]
                    .graph
                    .get_external_signature();

                let dependent_momenta_constructor =
                    DependentMomentaConstructor::Amplitude(&externals);

                let loop_number = lmb.loop_edges.len();
                let orientations = requested_orientations(self, graph_id, approach_settings)?;
                let mut reports = Vec::with_capacity(orientations.len());

                for (orientation, orientation_label) in orientations {
                    let mut limit_data = LimitData { data: Vec::new() };

                    for (loop_mom_id, lambda_point) in momenta.iter().cloned().enumerate() {
                        let mut loop_moms: LoopMomenta<F<_>> = (0..loop_number)
                            .map(|_| {
                                ThreeMomentum::new(
                                    F::from_f64(0.0),
                                    F::from_f64(0.0),
                                    F::from_f64(0.0),
                                )
                            })
                            .collect();

                        for (loop_id, momentum) in non_limit_momenta.iter() {
                            loop_moms[*loop_id] = *momentum;
                        }

                        for tagged_momenta in &lambda_point.momenta {
                            let edge_id = tagged_momenta.tag;
                            let loop_id = lmb
                                .loop_edges
                                .iter()
                                .position(|loop_edge| loop_edge == &edge_id)
                                .unwrap_or_else(|| {
                                    unreachable!("corrupted lmb and ir limit: {}", ir_limit);
                                });

                            loop_moms[LoopIndex(loop_id)] = tagged_momenta.momentum;
                        }

                        let sample_in_cmb = MomentumSample::new(
                            loop_moms,
                            loop_mom_id,
                            &self.settings.kinematics.externals,
                            0,
                            F::from_f64(1.0),
                            dependent_momenta_constructor,
                            orientation,
                        )?;

                        let sample = sample_in_cmb.lmb_transform(
                            &lmb,
                            &self.data.graph_terms[graph_id].graph.loop_momentum_basis,
                        );

                        limit_data.data.push(LambdaPointEval {
                            lambda: lambda_point.lambda,
                            value: evaluate_profile_momentum_point_arb(
                                self,
                                model,
                                graph_id,
                                orientation,
                                sample.loop_moms().iter().cloned().collect_vec(),
                                approach_settings.show_per_cut_info,
                            )?,
                        });
                    }

                    let (power_law_fit, cut_fits, component_fits) = limit_data.extract_power()?;
                    let mut report = build_single_limit_report(
                        ir_limit,
                        orientation_label,
                        power_law_fit,
                        build_cut_limit_reports(ir_limit.num_soft(), cut_fits),
                    );
                    report.display_only_reports = build_display_only_limit_reports(component_fits);
                    reports.push(report);
                }

                Ok(reports)
            }
            ProfileLimit::Threshold(threshold_limit) => {
                let overlap_structure = points_on_threshold.get(graph_id).ok_or_else(|| {
                    eyre!(
                        "Missing threshold-approach kinematics for amplitude graph {}",
                        graph_id
                    )
                })?;

                let existing_esurface_id = threshold_limit.resolve_existing_esurface_id(
                    &self.data.graph_terms[graph_id]
                        .threshold_counterterm
                        .esurface_map,
                    self.data.graph_terms[graph_id]
                        .threshold_counterterm
                        .own_group_position,
                    &overlap_structure.existing_esurfaces,
                )?;

                let momenta_per_overlap_group = threshold_limit.get_momenta_per_overlap_group(
                    overlap_structure,
                    existing_esurface_id,
                    approach_settings,
                )?;
                let orientations = requested_orientations(self, graph_id, approach_settings)?;
                let mut reports =
                    Vec::with_capacity(orientations.len() * momenta_per_overlap_group.len());

                for (overlap_group_label, momenta) in momenta_per_overlap_group {
                    for (orientation, orientation_label) in &orientations {
                        let mut limit_data = LimitData { data: Vec::new() };

                        for lambda_point in &momenta {
                            limit_data.data.push(LambdaPointEval {
                                lambda: lambda_point.lambda.clone(),
                                value: evaluate_profile_momentum_point_arb(
                                    self,
                                    model,
                                    graph_id,
                                    *orientation,
                                    lambda_point
                                        .loop_momenta
                                        .iter()
                                        .map(|momentum| {
                                            ThreeMomentum::new(
                                                F::from_f64(momentum.px.into_f64()),
                                                F::from_f64(momentum.py.into_f64()),
                                                F::from_f64(momentum.pz.into_f64()),
                                            )
                                        })
                                        .collect_vec(),
                                    approach_settings.show_per_cut_info,
                                )?,
                            });
                        }

                        let (power_law_fit, cut_fits, component_fits) =
                            limit_data.extract_power()?;
                        let context_label = Some(match orientation_label.as_deref() {
                            Some(orientation_label) => {
                                format!("{overlap_group_label} / {orientation_label}")
                            }
                            None => overlap_group_label.clone(),
                        });

                        let mut report = build_threshold_limit_report(
                            threshold_limit,
                            context_label,
                            power_law_fit,
                            build_cut_limit_reports(0, cut_fits),
                        );
                        report.display_only_reports =
                            build_display_only_limit_reports(component_fits);
                        reports.push(report);
                    }
                }

                Ok(reports)
            }
        }
    }
}

impl CrossSectionIntegrand {
    pub fn ir_profile_completion_entries(&self) -> Vec<(String, Vec<String>)> {
        ir_profile_completion_entries(self.enumerate_ir_limits())
    }

    pub fn test_ir(
        &mut self,
        ir_profile_settings: &IRProfileSetting,
        model: &Model,
    ) -> Result<IrLimitTestReport> {
        // override the sampling to be in momentum space
        self.settings.sampling = SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings {
            sample_orientations: false,
            sampling_type: DiscreteGraphSamplingType::Default(ParameterizationSettings {
                mode: ParameterizationMode::MomentumSpace,
                mapping: ParameterizationMapping::default(),
                b: 10.0,
            }),
        });

        let previous_generate_events = self.settings.general.generate_events;
        if ir_profile_settings.show_per_cut_info {
            self.settings.general.generate_events = true;
        }

        let points_on_thrshold = vec![]; // threshold limits are not yet supported for cross-section IR profiling

        let result = (|| {
            self.warm_up(model)?;
            run_ir_profile(
                self,
                ir_profile_settings,
                model,
                Self::enumerate_ir_limits,
                Self::graph_cut_definitions,
                &points_on_thrshold,
                Self::test_single_ir_limit_impl,
            )
        })();

        if self.settings.general.generate_events != previous_generate_events {
            self.settings.general.generate_events = previous_generate_events;
        }

        result
    }

    fn enumerate_ir_limits(&self) -> Vec<(String, Vec<ProfileLimit>)> {
        self.data
            .graph_terms
            .iter()
            .map(|term| {
                let graph_name = term.graph.name.clone();
                let limits = term
                    .enumerate_ir_limits()
                    .into_iter()
                    .map(ProfileLimit::Ir)
                    .collect();
                (graph_name, limits)
            })
            .collect()
    }

    fn graph_cut_definitions(&self, graph_id: usize) -> Vec<GraphCutDefinition> {
        graph_cut_definitions_for_cross_section_term(&self.data.graph_terms[graph_id])
    }

    fn test_single_ir_limit_impl(
        &mut self,
        graph_id: usize,
        profile_limit: &ProfileLimit,
        rng: &mut MonteCarloRng,
        approach_settings: &IRProfileSetting,
        model: &Model,
        _points_on_threshold: &[OverlapStructureWithKinematics<ArbPrec>],
    ) -> Result<Vec<SingleLimitReport>> {
        let ir_limit = match profile_limit {
            ProfileLimit::Ir(ir_limit) => ir_limit,
            ProfileLimit::Threshold(threshold_limit) => {
                return Err(eyre!(
                    "Threshold limit '{}' is not yet supported in cross-section IR profiling",
                    threshold_limit
                ));
            }
        };

        let edges_in_limit = ir_limit.get_all_edges()?;

        // find cut that for that as all edges of the limit
        let (cut_id, _esurface) = self.data.graph_terms[graph_id]
            .cut_esurface
            .iter_enumerated()
            .find(|(_cut_id, esurface)| {
                edges_in_limit
                    .iter()
                    .all(|edge| esurface.energies.contains(edge))
            })
            .ok_or(eyre!(
                "could not find cut with all edges of the limit: {}",
                ir_limit
            ))?;

        let cs_cut = &self.data.graph_terms[graph_id].cuts[cut_id];

        let edges_to_flip = cs_cut
            .cut
            .iter_edges(&self.data.graph_terms[graph_id].graph.underlying)
            .map(|(or, _)| or)
            .zip(
                self.data.graph_terms[graph_id]
                    .graph
                    .underlying
                    .iter_edges_of(&cs_cut.cut)
                    .map(|x| x.1),
            )
            .filter_map(|(orientation, edge_id)| {
                if edges_in_limit.contains(&edge_id) && matches!(orientation, Orientation::Reversed)
                {
                    Some(edge_id)
                } else {
                    None
                }
            })
            .collect_vec();

        let lmb = self.data.graph_terms[graph_id].lmb_with_loop_edges(edges_in_limit.as_slice())?;
        let momenta = ir_limit.get_momenta(rng, &self.settings, approach_settings)?;
        let non_limit_loops = lmb
            .loop_edges
            .iter_enumerated()
            .filter_map(|(loop_id, edge_id)| {
                if !edges_in_limit.contains(edge_id) {
                    Some(loop_id)
                } else {
                    None
                }
            })
            .collect_vec();

        let non_limit_momenta = non_limit_loops
            .iter()
            .map(|loop_id| (*loop_id, sample_random_unit_vector(rng)))
            .collect_vec();

        let dependent_momenta_constructor = DependentMomentaConstructor::CrossSection;

        let loop_number = lmb.loop_edges.len();
        let orientations = requested_orientations(self, graph_id, approach_settings)?;
        let mut reports = Vec::with_capacity(orientations.len());

        for (orientation, orientation_label) in orientations {
            let mut limit_data = LimitData { data: Vec::new() };

            for (loop_mom_id, lambda_point) in momenta.iter().cloned().enumerate() {
                let mut loop_moms: LoopMomenta<F<_>> = (0..loop_number)
                    .map(|_| {
                        ThreeMomentum::new(F::from_f64(0.0), F::from_f64(0.0), F::from_f64(0.0))
                    })
                    .collect();

                for (loop_id, momentum) in non_limit_momenta.iter() {
                    loop_moms[*loop_id] = *momentum;
                }

                for tagged_momenta in &lambda_point.momenta {
                    let edge_id = tagged_momenta.tag;
                    let loop_id = lmb
                        .loop_edges
                        .iter()
                        .position(|loop_edge| loop_edge == &edge_id)
                        .unwrap_or_else(|| {
                            unreachable!("corrupted lmb and ir limit: {}", ir_limit);
                        });

                    loop_moms[LoopIndex(loop_id)] = if edges_to_flip.contains(&edge_id) {
                        -tagged_momenta.momentum
                    } else {
                        tagged_momenta.momentum
                    };
                }

                let sample_in_cmb = MomentumSample::new(
                    loop_moms,
                    loop_mom_id,
                    &self.settings.kinematics.externals,
                    0,
                    F::from_f64(1.0),
                    dependent_momenta_constructor,
                    orientation,
                )?;

                let sample = sample_in_cmb.lmb_transform(
                    &lmb,
                    &self.data.graph_terms[graph_id].graph.loop_momentum_basis,
                );

                limit_data.data.push(LambdaPointEval {
                    lambda: lambda_point.lambda,
                    value: evaluate_profile_momentum_point_arb(
                        self,
                        model,
                        graph_id,
                        orientation,
                        sample.loop_moms().iter().cloned().collect_vec(),
                        approach_settings.show_per_cut_info,
                    )?,
                });
            }

            let (power_law_fit, cut_fits, _component_fits) = limit_data.extract_power()?;
            reports.push(build_single_limit_report(
                ir_limit,
                orientation_label,
                power_law_fit,
                build_cut_limit_reports(ir_limit.num_soft(), cut_fits),
            ));
        }

        Ok(reports)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct IrLimit {
    colinear: Vec<Vec<HardOrSoft>>,
    soft: Vec<EdgeIndex>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct ThresholdLimit {
    esurface_id: EsurfaceID,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum ProfileLimit {
    Ir(IrLimit),
    Threshold(ThresholdLimit),
}

enum MomentumBuilder<T: FloatLike> {
    Colinear {
        edge_id: EdgeIndex,
        x: F<T>,
        colinear_direction: ThreeMomentum<F<T>>,
        perpendicular_direction: ThreeMomentum<F<T>>,
        is_soft: bool,
    },
    Soft {
        edge_id: EdgeIndex,
        direction: ThreeMomentum<F<T>>,
    },
}

impl Display for IrLimit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for colinear_set in &self.colinear {
            write!(f, "C[")?;

            let mut iter = colinear_set.iter();
            if let Some(last) = iter.next_back() {
                for item in iter {
                    write!(f, "{},", item)?;
                }
                write!(f, "{}]", last)?;
            } else {
                write!(f, "]")?;
            }
        }

        for soft in self.soft.iter() {
            write!(f, "S({})", soft)?;
        }

        Ok(())
    }
}

impl Display for ThresholdLimit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "T(t{})", self.esurface_id.0)
    }
}

impl Display for ProfileLimit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ProfileLimit::Ir(ir_limit) => write!(f, "{}", ir_limit),
            ProfileLimit::Threshold(threshold_limit) => write!(f, "{}", threshold_limit),
        }
    }
}

impl ThresholdLimit {
    fn enumerate_from_overlap_structure(
        existing_esurfaces: &ExistingEsurfaces,
        esurface_map: &TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<EsurfaceID>>>,
        own_group_position: GraphGroupPosition,
    ) -> Vec<Self> {
        existing_esurfaces
            .iter()
            .filter_map(|group_esurface_id| {
                esurface_map[*group_esurface_id][own_group_position]
                    .map(|esurface_id| Self { esurface_id })
            })
            .collect::<HashSet<_>>()
            .into_iter()
            .sorted()
            .collect()
    }

    fn parse_threshold(threshold: &str) -> Result<Self> {
        let mut threshold = String::from(threshold);
        threshold = String::from(threshold.trim());

        if threshold.len() < 2 {
            return Err(eyre!("Threshold must be at least two characters long"));
        }

        if threshold.remove(0) != 't' {
            return Err(eyre!("Threshold must start with 't'"));
        }

        let threshold_id: usize = threshold
            .parse()
            .map_err(|_| eyre!("Threshold must be a valid integer, got: {}", threshold))?;

        Ok(Self {
            esurface_id: EsurfaceID::from(threshold_id),
        })
    }

    fn resolve_existing_esurface_id(
        &self,
        esurface_map: &TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<EsurfaceID>>>,
        own_group_position: GraphGroupPosition,
        existing_esurfaces: &ExistingEsurfaces,
    ) -> Result<ExistingEsurfaceId> {
        existing_esurfaces
            .iter_enumerated()
            .find_map(|(existing_esurface_id, group_esurface_id)| {
                esurface_map.get(*group_esurface_id).and_then(|graph_map| {
                    (graph_map[own_group_position] == Some(self.esurface_id))
                        .then_some(existing_esurface_id)
                })
            })
            .ok_or_else(|| {
                eyre!(
                    "Threshold '{}' does not exist in the selected overlap structure",
                    self
                )
            })
    }

    fn get_momenta_per_overlap_group<T: FloatLike>(
        &self,
        overlap_structure: &OverlapStructureWithKinematics<T>,
        existing_esurface_id: ExistingEsurfaceId,
        approach_settings: &IRProfileSetting,
    ) -> Result<Vec<(String, Vec<LambdaLoopMomentaPoint<T>>)>> {
        let lambda_values = constant_dropped_fit_points(
            &F::from_f64(10.0_f64.powf(approach_settings.lambda_exp_start)),
            &F::from_f64(10.0_f64.powf(approach_settings.lambda_exp_end)),
            approach_settings.steps,
        )?;

        let mut momenta_per_overlap_group = Vec::new();

        for (overlap_group_id, overlap_group_with_kinematics) in overlap_structure
            .overlap_groups_with_kinematics
            .iter()
            .enumerate()
        {
            let mut contains_existing_esurface = false;
            let mut threshold_point = None;

            for (group_existing_esurface_id, maybe_threshold_point) in overlap_group_with_kinematics
                .overlap_group
                .existing_esurfaces
                .iter()
                .copied()
                .zip(
                    overlap_group_with_kinematics
                        .loop_momenta_at_esurface
                        .iter(),
                )
            {
                if group_existing_esurface_id != existing_esurface_id {
                    continue;
                }

                contains_existing_esurface = true;
                threshold_point = maybe_threshold_point.as_ref();
                break;
            }

            if !contains_existing_esurface {
                continue;
            }

            let threshold_point = threshold_point.ok_or_else(|| {
                eyre!(
                    "Threshold '{}' is missing stored approach kinematics for overlap group {}",
                    self,
                    overlap_group_id
                )
            })?;

            momenta_per_overlap_group.push((
                format!("overlap group {}", overlap_group_id),
                lambda_values
                    .iter()
                    .cloned()
                    .map(|lambda| LambdaLoopMomentaPoint {
                        loop_momenta: threshold_approach_loop_momenta(
                            &overlap_group_with_kinematics.overlap_group.center,
                            threshold_point,
                            &lambda,
                        ),
                        lambda,
                    })
                    .collect(),
            ));
        }

        if momenta_per_overlap_group.is_empty() {
            return Err(eyre!(
                "Threshold '{}' does not appear in any overlap group for the selected graph term",
                self
            ));
        }

        Ok(momenta_per_overlap_group)
    }
}

impl ProfileLimit {
    fn parse_limit(limit: &str) -> Result<Self> {
        let mut colinear_sets = Vec::new();
        let mut soft_edges = Vec::new();
        let mut threshold_limit = None;

        let mut char_iter = limit.chars().enumerate();

        while let Some((char_position, char)) = char_iter.next() {
            match char {
                'C' => {
                    let (_opening_bracket_position, opening_bracket) =
                        char_iter.next().ok_or_else(|| {
                            eyre!(
                                "Expected opening bracket after 'C' at position {}",
                                char_position
                            )
                        })?;

                    if opening_bracket != '[' {
                        return Err(eyre!(
                            "Expected '[' after 'C' at position {}, found '{}'",
                            char_position,
                            opening_bracket
                        ));
                    }

                    let mut colinear_set_str = String::new();

                    let mut closing_bracket_found = false;
                    for (_next_char_position, next_char) in char_iter.by_ref() {
                        if next_char == ']' {
                            closing_bracket_found = true;
                            break;
                        }
                        colinear_set_str.push(next_char);
                    }

                    if !closing_bracket_found {
                        return Err(eyre!(
                            "Expected closing bracket ']' for colinear set at position {}",
                            char_position
                        ));
                    }

                    let edges = colinear_set_str.trim().split(',');

                    let mut colinear_set = Vec::new();

                    for edge in edges {
                        let trimmed_edge = edge.trim();
                        if trimmed_edge.is_empty() {
                            return Err(eyre!(
                                "Empty edge found in colinear set at position {}",
                                char_position
                            ));
                        }

                        if trimmed_edge.starts_with('S') {
                            let mut trimmed_edge_iter = trimmed_edge.chars().skip(1);
                            let opening_bracket = trimmed_edge_iter
                                .next()
                                .ok_or(eyre!("Expected '(' after 'S' in soft edge at position"))?;

                            if opening_bracket != '(' {
                                return Err(eyre!(
                                    "Expected '(' after 'S' in soft edge at position , found ''",
                                ));
                            }

                            let mut edge_str = String::new();
                            let mut closing_bracket_found = false;

                            for next_char in trimmed_edge_iter {
                                if next_char == ')' {
                                    closing_bracket_found = true;
                                    break;
                                }
                                edge_str.push(next_char);
                            }

                            if !closing_bracket_found {
                                return Err(eyre!(
                                    "Expected closing bracket ')' for soft edge at position {}",
                                    char_position
                                ));
                            }

                            let edge_index = IrLimit::parse_edge(&edge_str)?;
                            colinear_set.push(HardOrSoft::Soft(edge_index));
                        } else {
                            let edge_index = IrLimit::parse_edge(trimmed_edge)?;
                            colinear_set.push(HardOrSoft::Hard(edge_index));
                        }
                    }
                    colinear_sets.push(colinear_set);
                }
                'S' => {
                    let (_opening_bracket_position, opening_bracket) =
                        char_iter.next().ok_or_else(|| {
                            eyre!(
                                "Expected opening bracket '(' after 'S' at position {}",
                                char_position
                            )
                        })?;

                    if opening_bracket != '(' {
                        return Err(eyre!(
                            "Expected '(' after 'S' at position {}, found '{}'",
                            char_position,
                            opening_bracket
                        ));
                    }

                    let mut edge_str = String::new();
                    let mut closing_bracket_found = false;

                    for (_next_char_position, next_char) in char_iter.by_ref() {
                        if next_char == ')' {
                            closing_bracket_found = true;
                            break;
                        }
                        edge_str.push(next_char);
                    }

                    if !closing_bracket_found {
                        return Err(eyre!(
                            "Expected closing bracket ')' for soft edge at position {}",
                            char_position
                        ));
                    }

                    let edge_index = IrLimit::parse_edge(&edge_str)?;
                    soft_edges.push(edge_index);
                }
                'T' => {
                    if threshold_limit.is_some() {
                        return Err(eyre!("Only one threshold limit can be specified"));
                    }

                    let (_opening_bracket_position, opening_bracket) =
                        char_iter.next().ok_or_else(|| {
                            eyre!(
                                "Expected opening bracket '(' after 'T' at position {}",
                                char_position
                            )
                        })?;

                    if opening_bracket != '(' {
                        return Err(eyre!(
                            "Expected '(' after 'T' at position {}, found '{}'",
                            char_position,
                            opening_bracket
                        ));
                    }

                    let mut threshold_str = String::new();
                    let mut closing_bracket_found = false;

                    for (_next_char_position, next_char) in char_iter.by_ref() {
                        if next_char == ')' {
                            closing_bracket_found = true;
                            break;
                        }
                        threshold_str.push(next_char);
                    }

                    if !closing_bracket_found {
                        return Err(eyre!(
                            "Expected closing bracket ')' for threshold at position {}",
                            char_position
                        ));
                    }

                    threshold_limit = Some(ThresholdLimit::parse_threshold(&threshold_str)?);
                }
                _ => {
                    return Err(eyre!(
                        "Unexpected character '{}' at position {}",
                        char,
                        char_position
                    ));
                }
            }
        }

        if let Some(threshold_limit) = threshold_limit {
            if !colinear_sets.is_empty() || !soft_edges.is_empty() {
                return Err(eyre!(
                    "Threshold limits cannot be combined with soft or colinear limits"
                ));
            }

            return Ok(ProfileLimit::Threshold(threshold_limit));
        }

        let mut ir_limit = IrLimit {
            colinear: colinear_sets,
            soft: soft_edges,
        };

        ir_limit.canonize();

        Ok(ProfileLimit::Ir(ir_limit))
    }

    fn is_valid(&self, loop_number: usize) -> Result<()> {
        match self {
            ProfileLimit::Ir(ir_limit) => ir_limit.is_valid(loop_number),
            ProfileLimit::Threshold(_) => Ok(()),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum HardOrSoft {
    Hard(EdgeIndex),
    Soft(EdgeIndex),
}

impl HardOrSoft {
    fn index(&self) -> EdgeIndex {
        match self {
            HardOrSoft::Hard(index) => *index,
            HardOrSoft::Soft(index) => *index,
        }
    }
}

impl Display for HardOrSoft {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HardOrSoft::Hard(index) => write!(f, "{}", index),
            HardOrSoft::Soft(index) => write!(f, "S({})", index),
        }
    }
}

impl IrLimit {
    fn canonize(&mut self) {
        for colinear_set in &mut self.colinear.iter_mut() {
            colinear_set.sort();
        }

        self.colinear.sort();
        self.soft.sort();
    }

    fn new_pure_colinear(colinear_edges: Vec<EdgeIndex>) -> Self {
        let colinear = vec![colinear_edges.into_iter().map(HardOrSoft::Hard).collect()];

        let mut result = IrLimit {
            colinear,
            soft: Vec::new(),
        };
        result.canonize();
        result
    }

    fn new_pure_soft(soft_edges: Vec<EdgeIndex>) -> Self {
        let mut result = IrLimit {
            colinear: Vec::new(),
            soft: soft_edges,
        };
        result.canonize();
        result
    }

    fn num_soft(&self) -> usize {
        self.colinear
            .iter()
            .flatten()
            .filter(|edge| matches!(edge, HardOrSoft::Soft(_)))
            .count()
            + self.soft.len()
    }

    fn check_min_colinear_size(&self) -> bool {
        self.colinear
            .iter()
            .all(|colinear_set| colinear_set.len() >= 2)
    }

    fn is_valid(&self, loop_number: usize) -> Result<()> {
        if !self.check_min_colinear_size() {
            return Err(eyre!("colinear sets must have at least two edges"));
        }

        let all_edges = self.get_all_edges()?;

        if all_edges.len() > loop_number {
            return Err(eyre!("not enough degrees of freedom to setup IR limit"));
        }

        Ok(())
    }

    fn get_all_edges(&self) -> Result<Vec<EdgeIndex>> {
        let colinear_edges = self
            .colinear
            .iter()
            .flatten()
            .map(HardOrSoft::index)
            .collect_vec();
        let soft_edges = self.soft.iter().copied().collect_vec();

        let all_edges: Vec<EdgeIndex> = colinear_edges
            .into_iter()
            .chain(soft_edges)
            .sorted()
            .collect();

        // check for duplicates
        let mut unique_edges = all_edges.clone();
        unique_edges.dedup();

        if unique_edges.len() != all_edges.len() {
            return Err(eyre!("Edges specified in ir limit must be unique")); // duplicates found
        }

        Ok(all_edges)
    }

    fn parse_edge(edge: &str) -> Result<EdgeIndex> {
        let mut edge = String::from(edge);
        edge = String::from(edge.trim());

        if edge.len() < 2 {
            return Err(eyre!("Edge must be at least two characters long"));
        }

        if edge.remove(0) != 'e' {
            return Err(eyre!("Edge must start with 'e'"));
        }

        let edge_id: usize = edge
            .parse()
            .map_err(|_| eyre!("Edge must be a valid integer, got: {}", edge))?;

        Ok(EdgeIndex::from(edge_id))
    }

    fn get_momentum_builders(&self, rng: &mut MonteCarloRng) -> Vec<MomentumBuilder<f64>> {
        let mut momentum_builder = Vec::new();

        for colinear_set in &self.colinear {
            let direction_for_set: ThreeMomentum<F<f64>> = sample_random_unit_vector(rng);

            let x_variables: Vec<F<f64>> = (0..colinear_set.len())
                .map(|_| F::from_f64(rng.random::<f64>() * 0.8 + 0.1))
                .sorted_by(|a, b| a.partial_cmp(b).unwrap())
                .collect_vec();

            for (edge, x) in colinear_set.iter().zip(x_variables) {
                let edge_id = edge.index();
                let direction: ThreeMomentum<F<f64>> = sample_random_unit_vector(rng);

                let perpendicular = direction - direction * (direction * direction_for_set);

                let perpendicular_norm = perpendicular.norm();
                let perpendicular = perpendicular * perpendicular_norm.inv();

                let is_soft = matches!(edge, HardOrSoft::Soft(_));

                momentum_builder.push(MomentumBuilder::Colinear {
                    edge_id,
                    x,
                    colinear_direction: direction_for_set,
                    perpendicular_direction: perpendicular,
                    is_soft,
                });
            }
        }

        for soft_edge in &self.soft {
            let direction = sample_random_unit_vector(rng);
            momentum_builder.push(MomentumBuilder::Soft {
                edge_id: *soft_edge,
                direction,
            });
        }

        momentum_builder
    }

    fn get_momenta(
        &self,
        rng: &mut MonteCarloRng,
        settings: &RuntimeSettings,
        approach_settings: &IRProfileSetting,
    ) -> Result<Vec<LambdaPoint<f64>>> {
        let momentum_builders = self.get_momentum_builders(rng);

        let lambda_values = constant_dropped_fit_points(
            &F::from_f64(10.0_f64.powf(approach_settings.lambda_exp_start)),
            &F::from_f64(10.0_f64.powf(approach_settings.lambda_exp_end)),
            approach_settings.steps,
        )?;

        Ok(lambda_values
            .into_iter()
            .map(|lambda| LambdaPoint {
                momenta: momentum_builders
                    .iter()
                    .map(|builder| match builder {
                        MomentumBuilder::Colinear {
                            edge_id,
                            x,
                            colinear_direction,
                            perpendicular_direction,
                            is_soft,
                        } => {
                            let momentum = if *is_soft {
                                (colinear_direction * x + perpendicular_direction * lambda)
                                    * F::from_f64(settings.kinematics.e_cm)
                                    * lambda
                            } else {
                                (colinear_direction * x + perpendicular_direction * lambda)
                                    * F::from_f64(settings.kinematics.e_cm)
                            };
                            TaggedMomenta {
                                momentum,
                                tag: *edge_id,
                            }
                        }
                        MomentumBuilder::Soft { edge_id, direction } => {
                            let momentum =
                                direction * lambda * F::from_f64(settings.kinematics.e_cm);
                            TaggedMomenta {
                                momentum,
                                tag: *edge_id,
                            }
                        }
                    })
                    .collect(),
                lambda,
            })
            .collect())
    }
}

fn evaluate_profile_momentum_point_arb<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    model: &Model,
    graph_id: usize,
    orientation: Option<usize>,
    loop_momenta: Vec<ThreeMomentum<F<f64>>>,
    show_per_cut_info: bool,
) -> Result<ProfilePointValue> {
    match evaluate_profile_momentum_point_precise(
        integrand,
        model,
        graph_id,
        orientation,
        loop_momenta,
        true,
    )? {
        PreciseEvaluationResult::Arb(result) => {
            let zero = result.integrand_result.re.zero();
            let per_cut = if show_per_cut_info {
                let mut per_cut = BTreeMap::new();
                for event_group in result.event_groups.iter() {
                    for event in event_group.iter() {
                        let entry = per_cut
                            .entry(event.cut_info.cut_id)
                            .or_insert_with(|| zero.clone());
                        *entry += event.weight.re.clone();
                    }
                }
                per_cut
            } else {
                BTreeMap::new()
            };

            let mut display_only_components = BTreeMap::new();
            for event_group in result.event_groups.iter() {
                for event in event_group.iter() {
                    for (key, weight) in &event.additional_weights.weights {
                        if !matches!(
                            key,
                            AdditionalWeightKey::Original
                                | AdditionalWeightKey::ThresholdCounterterm { .. }
                        ) {
                            continue;
                        }

                        let entry = display_only_components
                            .entry(*key)
                            .or_insert_with(|| zero.clone());
                        *entry += weight.re.clone();
                    }
                }
            }

            Ok(ProfilePointValue {
                total: result.integrand_result.re,
                per_cut,
                display_only_components,
            })
        }
        PreciseEvaluationResult::Double(_) => Err(eyre!(
            "IR profiling requested arbitrary precision but received a double-precision result"
        )),
        PreciseEvaluationResult::Quad(_) => Err(eyre!(
            "IR profiling requested arbitrary precision but received a quad-precision result"
        )),
    }
}

fn sample_random_unit_vector<T: FloatLike>(rng: &mut MonteCarloRng) -> ThreeMomentum<F<T>> {
    let x_1 = F::<T>::from_f64(rng.random::<f64>());
    let x_2 = F::<T>::from_f64(rng.random::<f64>());
    let x_3 = F::<T>::from_f64(rng.random::<f64>());
    let x_4 = F::<T>::from_f64(rng.random::<f64>());

    let (k_x, k_y) = box_muller(x_1, x_2);
    let k_z = box_muller(x_3, x_4).0;

    let unnormalized_momentum = ThreeMomentum::new(k_x, k_y, k_z);
    let norm = unnormalized_momentum.norm();
    unnormalized_momentum * norm.inv()
}

#[derive(Clone)]
struct TaggedMomenta<T> {
    momentum: ThreeMomentum<T>,
    tag: EdgeIndex,
}

#[derive(Clone)]
struct LambdaPoint<T: FloatLike> {
    lambda: F<T>,
    momenta: Vec<TaggedMomenta<F<T>>>,
}

#[derive(Clone)]
struct LambdaLoopMomentaPoint<T: FloatLike> {
    lambda: F<T>,
    loop_momenta: LoopMomenta<F<T>>,
}

struct LambdaPointEval<T: FloatLike> {
    lambda: F<T>,
    value: ProfilePointValue,
}

struct ProfilePointValue {
    total: F<ArbPrec>,
    per_cut: BTreeMap<usize, F<ArbPrec>>,
    display_only_components: BTreeMap<AdditionalWeightKey, F<ArbPrec>>,
}

struct LimitData<T: FloatLike> {
    data: Vec<LambdaPointEval<T>>,
}

impl<T: FloatLike> LimitData<T> {
    fn extract_power(
        &self,
    ) -> Result<(
        PowerLawFit,
        Vec<(usize, Result<PowerLawFit>)>,
        Vec<(AdditionalWeightKey, Result<PowerLawFit>)>,
    )> {
        let x = self
            .data
            .iter()
            .map(|point_eval| F::<ArbPrec>::from_ff64(point_eval.lambda.into_ff64()))
            .collect_vec();

        let y = self
            .data
            .iter()
            .map(|point_eval| point_eval.value.total.clone())
            .collect_vec();

        let result = fit_power_law(x.clone(), y.clone())?;

        let zero = x
            .first()
            .map(|value| value.zero())
            .unwrap_or_else(|| F::<ArbPrec>::from_f64(0.0));
        let cut_fits = self
            .data
            .iter()
            .flat_map(|point_eval| point_eval.value.per_cut.keys().copied())
            .unique()
            .sorted()
            .map(|cut_id| {
                let y = self
                    .data
                    .iter()
                    .map(|point_eval| {
                        point_eval
                            .value
                            .per_cut
                            .get(&cut_id)
                            .cloned()
                            .unwrap_or_else(|| zero.clone())
                    })
                    .collect_vec();
                (cut_id, fit_power_law(x.clone(), y))
            })
            .collect_vec();

        let component_fits = self
            .data
            .iter()
            .flat_map(|point_eval| point_eval.value.display_only_components.keys().copied())
            .filter(|key| {
                matches!(
                    key,
                    AdditionalWeightKey::Original
                        | AdditionalWeightKey::ThresholdCounterterm { .. }
                )
            })
            .unique()
            .sorted()
            .map(|key| {
                let y = self
                    .data
                    .iter()
                    .map(|point_eval| {
                        point_eval
                            .value
                            .display_only_components
                            .get(&key)
                            .cloned()
                            .unwrap_or_else(|| zero.clone())
                    })
                    .collect_vec();
                (key, fit_power_law(x.clone(), y))
            })
            .collect_vec();

        if result.r_squared < 0.9 {
            warn!("low r^2 value found for input data");
            warn!(
                "x: {:?}",
                x.iter().map(|value| format!("{}", value)).collect_vec()
            );
            warn!(
                "y: {:?}",
                y.iter().map(|value| format!("{}", value)).collect_vec()
            );
        }

        Ok((result, cut_fits, component_fits))
    }
}

#[derive(Debug, Clone)]
pub struct PowerLawFit {
    exponent: f64,
    r_squared: f64,
}

fn fit_power_law(x: Vec<F<ArbPrec>>, y: Vec<F<ArbPrec>>) -> Result<PowerLawFit> {
    if x.len() != y.len() {
        return Err(eyre!(
            "fit_power_law requires x and y to have the same length"
        ));
    }
    if x.len() < 3 {
        return Err(eyre!("fit_power_law requires at least three observations"));
    }
    if x.iter()
        .any(|value| value.is_nan() || value.is_infinite() || value <= &value.zero())
    {
        return Err(eyre!(
            "fit_power_law requires strictly positive, finite x values"
        ));
    }
    if y.iter().any(|value| value.is_nan() || value.is_infinite()) {
        return Err(eyre!("fit_power_law requires finite y values"));
    }

    let fit = log_log_slope_constant_dropped(&x, &y)?;

    Ok(PowerLawFit {
        exponent: fit.slope.into_f64(),
        r_squared: fit.r_squared().into_f64(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use typed_index_collections::ti_vec;

    #[test]
    fn test_display() {
        let ir_limit = IrLimit {
            colinear: vec![
                vec![
                    HardOrSoft::Hard(EdgeIndex::from(1)),
                    HardOrSoft::Hard(EdgeIndex::from(2)),
                    HardOrSoft::Hard(EdgeIndex::from(3)),
                ],
                vec![
                    HardOrSoft::Hard(EdgeIndex::from(4)),
                    HardOrSoft::Soft(EdgeIndex::from(5)),
                ],
            ],
            soft: vec![EdgeIndex::from(6), EdgeIndex::from(7)],
        };

        let display = ir_limit.to_string();
        let expected = "C[e1,e2,e3]C[e4,S(e5)]S(e6)S(e7)";

        assert_eq!(display, expected);
    }

    #[test]
    fn test_threshold_display() {
        let threshold_limit = ThresholdLimit {
            esurface_id: EsurfaceID::from(8usize),
        };

        let display = threshold_limit.to_string();
        let expected = "T(t8)";

        assert_eq!(display, expected);
    }

    #[test]
    fn parse_edge() {
        let edge_str = "e5";
        let edge_index = IrLimit::parse_edge(edge_str).unwrap();
        assert_eq!(edge_index, EdgeIndex::from(5));

        let invalid_edge_str = "5"; // missing 'e'
        assert!(IrLimit::parse_edge(invalid_edge_str).is_err());

        let invalid_edge_str2 = "e"; // too short
        assert!(IrLimit::parse_edge(invalid_edge_str2).is_err());

        let invalid_edge_str3 = "e5a"; // not a valid integer
        assert!(IrLimit::parse_edge(invalid_edge_str3).is_err());
    }

    #[test]
    fn parse_threshold() {
        let threshold_str = "t5";
        let threshold_limit = ThresholdLimit::parse_threshold(threshold_str).unwrap();
        assert_eq!(threshold_limit.esurface_id, EsurfaceID::from(5usize));

        let invalid_threshold_str = "5"; // missing 't'
        assert!(ThresholdLimit::parse_threshold(invalid_threshold_str).is_err());

        let invalid_threshold_str2 = "t"; // too short
        assert!(ThresholdLimit::parse_threshold(invalid_threshold_str2).is_err());

        let invalid_threshold_str3 = "t5a"; // not a valid integer
        assert!(ThresholdLimit::parse_threshold(invalid_threshold_str3).is_err());
    }

    #[test]
    fn parse_limit() {
        let limit_str = "C[e1,e2,e3]C[e4,S(e5)]S(e6)S(e7)";
        let limit = ProfileLimit::parse_limit(limit_str).unwrap();
        let ir_limit = match limit {
            ProfileLimit::Ir(ir_limit) => ir_limit,
            ProfileLimit::Threshold(_) => panic!("Expected an IR limit"),
        };

        assert_eq!(ir_limit.colinear.len(), 2, "Expected two colinear sets");
        assert_eq!(ir_limit.soft.len(), 2, "Expected two soft edges");

        assert_eq!(
            ir_limit.colinear[0],
            vec![
                HardOrSoft::Hard(EdgeIndex::from(1)),
                HardOrSoft::Hard(EdgeIndex::from(2)),
                HardOrSoft::Hard(EdgeIndex::from(3))
            ],
            "First colinear set does not match"
        );
        assert_eq!(
            ir_limit.colinear[1],
            vec![
                HardOrSoft::Hard(EdgeIndex::from(4)),
                HardOrSoft::Soft(EdgeIndex::from(5))
            ],
            "Second colinear set does not match"
        );

        assert_eq!(
            ir_limit.soft,
            vec![EdgeIndex::from(6), EdgeIndex::from(7)],
            "Soft edges do not match"
        );

        let threshold_limit = ProfileLimit::parse_limit("T(t8)").unwrap();
        assert_eq!(
            threshold_limit,
            ProfileLimit::Threshold(ThresholdLimit {
                esurface_id: EsurfaceID::from(8usize),
            }),
            "Threshold limit does not match"
        );

        let invalid_limit_str = "C[e1,e2,e3]C[e4,S(e5)]S(e6)S(e7, e8)";
        assert!(
            ProfileLimit::parse_limit(invalid_limit_str).is_err(),
            "Expected error"
        );

        let invalid_limit_str2 = "C[e1,e2,e3C[e4,S(e5)]S(e6)";
        assert!(
            ProfileLimit::parse_limit(invalid_limit_str2).is_err(),
            "Expected error for unmatched brackets"
        );

        let invalid_limit_str3 = "C[e1,e2,e3]T(e8)";
        assert!(
            ProfileLimit::parse_limit(invalid_limit_str3).is_err(),
            "Expected error for invalid threshold syntax"
        );

        let invalid_limit_str4 = "C[e1,e2,e3]T(t8)";
        assert!(
            ProfileLimit::parse_limit(invalid_limit_str4).is_err(),
            "Expected error for mixed threshold and IR limit syntax"
        );

        let invalid_limit_str5 = "T(t8)T(t9)";
        assert!(
            ProfileLimit::parse_limit(invalid_limit_str5).is_err(),
            "Expected error for multiple threshold limits"
        );
    }

    #[test]
    fn resolve_existing_esurface_id_for_threshold_limit() {
        let threshold_limit = ThresholdLimit {
            esurface_id: EsurfaceID::from(7usize),
        };
        let esurface_map = ti_vec![
            ti_vec![
                Some(EsurfaceID::from(5usize)),
                Some(EsurfaceID::from(6usize))
            ],
            ti_vec![Some(EsurfaceID::from(7usize)), None],
        ];
        let existing_esurfaces =
            ti_vec![GroupEsurfaceId::from(0usize), GroupEsurfaceId::from(1usize)];

        let existing_esurface_id = threshold_limit
            .resolve_existing_esurface_id(
                &esurface_map,
                GraphGroupPosition::from(0usize),
                &existing_esurfaces,
            )
            .unwrap();

        assert_eq!(existing_esurface_id, ExistingEsurfaceId::from(1usize));
    }

    #[test]
    fn resolve_existing_esurface_id_rejects_threshold_missing_from_graph() {
        let threshold_limit = ThresholdLimit {
            esurface_id: EsurfaceID::from(9usize),
        };
        let esurface_map = ti_vec![ti_vec![
            Some(EsurfaceID::from(5usize)),
            Some(EsurfaceID::from(6usize))
        ]];
        let existing_esurfaces = ti_vec![GroupEsurfaceId::from(0usize)];

        assert!(
            threshold_limit
                .resolve_existing_esurface_id(
                    &esurface_map,
                    GraphGroupPosition::from(0usize),
                    &existing_esurfaces,
                )
                .is_err()
        );
    }

    #[test]
    fn resolve_existing_esurface_id_rejects_threshold_missing_from_overlap() {
        let threshold_limit = ThresholdLimit {
            esurface_id: EsurfaceID::from(7usize),
        };
        let esurface_map = ti_vec![ti_vec![Some(EsurfaceID::from(7usize))]];
        let existing_esurfaces = ti_vec![GroupEsurfaceId::from(1usize)];

        assert!(
            threshold_limit
                .resolve_existing_esurface_id(
                    &esurface_map,
                    GraphGroupPosition::from(0usize),
                    &existing_esurfaces,
                )
                .is_err()
        );
    }

    #[test]
    fn enumerate_threshold_limits_from_overlap_structure() {
        let esurface_map = ti_vec![
            ti_vec![
                Some(EsurfaceID::from(5usize)),
                Some(EsurfaceID::from(8usize))
            ],
            ti_vec![Some(EsurfaceID::from(7usize)), None],
            ti_vec![
                Some(EsurfaceID::from(5usize)),
                Some(EsurfaceID::from(9usize))
            ],
            ti_vec![None, Some(EsurfaceID::from(3usize))],
        ];
        let existing_esurfaces = ti_vec![
            GroupEsurfaceId::from(2usize),
            GroupEsurfaceId::from(0usize),
            GroupEsurfaceId::from(1usize),
            GroupEsurfaceId::from(3usize),
        ];

        let threshold_limits = ThresholdLimit::enumerate_from_overlap_structure(
            &existing_esurfaces,
            &esurface_map,
            GraphGroupPosition::from(0usize),
        );
        let threshold_limits_for_other_group = ThresholdLimit::enumerate_from_overlap_structure(
            &existing_esurfaces,
            &esurface_map,
            GraphGroupPosition::from(1usize),
        );

        assert_eq!(
            threshold_limits,
            vec![
                ThresholdLimit {
                    esurface_id: EsurfaceID::from(5usize),
                },
                ThresholdLimit {
                    esurface_id: EsurfaceID::from(7usize),
                },
            ]
        );
        assert_eq!(
            threshold_limits_for_other_group,
            vec![
                ThresholdLimit {
                    esurface_id: EsurfaceID::from(3usize),
                },
                ThresholdLimit {
                    esurface_id: EsurfaceID::from(8usize),
                },
                ThresholdLimit {
                    esurface_id: EsurfaceID::from(9usize),
                },
            ]
        );
    }

    fn test_ir_profile_settings(steps: usize) -> IRProfileSetting {
        IRProfileSetting {
            lambda_exp_start: -3.0,
            lambda_exp_end: -1.0,
            steps,
            seed: 0,
            select_limits_and_graphs: None,
            orientation_mode: OrientationProfileMode::Summed,
            show_per_cut_info: false,
        }
    }

    fn test_momentum_sample(loop_momenta: Vec<ThreeMomentum<F<f64>>>) -> MomentumSample<f64> {
        MomentumSample {
            sample: crate::momentum::sample::BareMomentumSample {
                loop_moms: loop_momenta.into_iter().collect(),
                dual_loop_moms: None,
                loop_mom_cache_id: 0,
                loop_mom_base_cache_id: 0,
                external_moms: ti_vec![],
                external_mom_cache_id: 0,
                external_mom_base_cache_id: 0,
                jacobian: F::from_f64(1.0),
                orientation: None,
            },
        }
    }

    #[test]
    fn threshold_approach_loop_momenta_starts_at_stored_threshold_point() {
        let overlap_group_center: LoopMomenta<F<f64>> = vec![ThreeMomentum::new(
            F::from_f64(5.0),
            F::from_f64(7.0),
            F::from_f64(9.0),
        )]
        .into_iter()
        .collect();
        let threshold_point = test_momentum_sample(vec![ThreeMomentum::new(
            F::from_f64(1.0),
            F::from_f64(3.0),
            F::from_f64(5.0),
        )]);

        let at_threshold = threshold_approach_loop_momenta(
            &overlap_group_center,
            &threshold_point,
            &F::from_f64(0.0),
        );
        let halfway = threshold_approach_loop_momenta(
            &overlap_group_center,
            &threshold_point,
            &F::from_f64(0.5),
        );
        let at_center = threshold_approach_loop_momenta(
            &overlap_group_center,
            &threshold_point,
            &F::from_f64(1.0),
        );

        assert_eq!(at_threshold, threshold_point.loop_moms().clone());
        assert_eq!(
            halfway,
            vec![ThreeMomentum::new(
                F::from_f64(3.0),
                F::from_f64(5.0),
                F::from_f64(7.0),
            )]
            .into_iter()
            .collect()
        );
        assert_eq!(at_center, overlap_group_center);
    }

    #[test]
    fn threshold_limit_builds_group_trajectories_for_matching_overlap_groups() {
        let threshold_limit = ThresholdLimit {
            esurface_id: EsurfaceID::from(7usize),
        };
        let threshold_point = test_momentum_sample(vec![ThreeMomentum::new(
            F::from_f64(1.0),
            F::from_f64(2.0),
            F::from_f64(3.0),
        )]);
        let overlap_structure = OverlapStructureWithKinematics {
            existing_esurfaces: ti_vec![
                GroupEsurfaceId::from(0usize),
                GroupEsurfaceId::from(1usize)
            ],
            overlap_groups_with_kinematics: vec![
                crate::subtraction::amplitude_counterterm::OverlapGroupWithKinematics {
                    overlap_group: crate::subtraction::overlap::OverlapGroup {
                        existing_esurfaces: vec![ExistingEsurfaceId::from(1usize)],
                        complement: vec![],
                        center: vec![ThreeMomentum::new(
                            F::from_f64(5.0),
                            F::from_f64(6.0),
                            F::from_f64(7.0),
                        )]
                        .into_iter()
                        .collect(),
                        prefactor_evaluator: None,
                    },
                    loop_momenta_at_esurface: ti_vec![Some(threshold_point.clone())],
                },
                crate::subtraction::amplitude_counterterm::OverlapGroupWithKinematics {
                    overlap_group: crate::subtraction::overlap::OverlapGroup {
                        existing_esurfaces: vec![ExistingEsurfaceId::from(0usize)],
                        complement: vec![],
                        center: vec![ThreeMomentum::new(
                            F::from_f64(8.0),
                            F::from_f64(9.0),
                            F::from_f64(10.0),
                        )]
                        .into_iter()
                        .collect(),
                        prefactor_evaluator: None,
                    },
                    loop_momenta_at_esurface: ti_vec![Some(test_momentum_sample(vec![
                        ThreeMomentum::new(F::from_f64(2.0), F::from_f64(3.0), F::from_f64(4.0),),
                    ]))],
                },
            ],
        };

        let trajectories = threshold_limit
            .get_momenta_per_overlap_group(
                &overlap_structure,
                ExistingEsurfaceId::from(1usize),
                &test_ir_profile_settings(4),
            )
            .unwrap();

        assert_eq!(trajectories.len(), 1);
        assert_eq!(trajectories[0].0, "overlap group 0");
        assert_eq!(trajectories[0].1.len(), 4);

        for lambda_point in &trajectories[0].1 {
            assert_eq!(
                lambda_point.loop_momenta,
                threshold_approach_loop_momenta(
                    &overlap_structure.overlap_groups_with_kinematics[0]
                        .overlap_group
                        .center,
                    &threshold_point,
                    &lambda_point.lambda,
                )
            );
        }
    }

    #[test]
    fn threshold_limit_rejects_group_missing_threshold_kinematics() {
        let threshold_limit = ThresholdLimit {
            esurface_id: EsurfaceID::from(7usize),
        };
        let overlap_structure: OverlapStructureWithKinematics<f64> =
            OverlapStructureWithKinematics {
                existing_esurfaces: ti_vec![GroupEsurfaceId::from(1usize)],
                overlap_groups_with_kinematics: vec![
                    crate::subtraction::amplitude_counterterm::OverlapGroupWithKinematics {
                        overlap_group: crate::subtraction::overlap::OverlapGroup {
                            existing_esurfaces: vec![ExistingEsurfaceId::from(0usize)],
                            complement: vec![],
                            center: vec![ThreeMomentum::new(
                                F::from_f64(5.0),
                                F::from_f64(6.0),
                                F::from_f64(7.0),
                            )]
                            .into_iter()
                            .collect(),
                            prefactor_evaluator: None,
                        },
                        loop_momenta_at_esurface: ti_vec![None],
                    },
                ],
            };

        assert!(
            threshold_limit
                .get_momenta_per_overlap_group(
                    &overlap_structure,
                    ExistingEsurfaceId::from(0usize),
                    &test_ir_profile_settings(3),
                )
                .is_err()
        );
    }

    #[test]
    fn fit_power_law_recovers_known_parameters() {
        let exponent = -1.75_f64;
        let coefficient = 3.2_f64;
        let offset = 0.6_f64;

        let x = constant_dropped_fit_points(
            &F::<ArbPrec>::from_f64(0.2_f64 * 1.6_f64.powi(7)),
            &F::<ArbPrec>::from_f64(0.2_f64),
            8,
        )
        .expect("geometric fit points should be generated");
        let y = x
            .iter()
            .map(|xv| F::<ArbPrec>::from_f64(coefficient * xv.into_f64().powf(exponent) + offset))
            .collect::<Vec<_>>();

        let fit = fit_power_law(x, y).expect("power-law fit should succeed");

        assert!((fit.exponent - exponent).abs() < 1e-10);
        assert!(fit.r_squared > 0.999_999);
    }

    #[test]
    fn fit_power_law_rejects_non_geometric_grid() {
        let exponent = -1.75_f64;
        let coefficient = 3.2_f64;
        let offset = 0.6_f64;

        let x = [
            0.2_f64, 0.37_f64, 0.55_f64, 0.92_f64, 1.3_f64, 1.85_f64, 2.75_f64, 3.6_f64,
        ]
        .into_iter()
        .map(F::<ArbPrec>::from_f64)
        .collect::<Vec<_>>();
        let y = x
            .iter()
            .map(|xv| F::<ArbPrec>::from_f64(coefficient * xv.into_f64().powf(exponent) + offset))
            .collect::<Vec<_>>();

        let fit = fit_power_law(x, y);

        assert!(fit.is_err());
    }

    #[test]
    fn negative_cut_scaling_does_not_fail_limit() {
        let ir_limit = IrLimit::new_pure_soft(vec![EdgeIndex::from(1)]);
        let total_fit = PowerLawFit {
            exponent: -2.5,
            r_squared: 0.999,
        };
        let cut_reports = build_cut_limit_reports(
            ir_limit.num_soft(),
            vec![
                (
                    0,
                    Ok(PowerLawFit {
                        exponent: -4.5,
                        r_squared: 0.999,
                    }),
                ),
                (
                    1,
                    Ok(PowerLawFit {
                        exponent: -2.5,
                        r_squared: 0.999,
                    }),
                ),
            ],
        );

        let report = build_single_limit_report(&ir_limit, None, total_fit, cut_reports);

        assert!(report.passed);
        assert_eq!(report.per_cut_reports.len(), 2);
        assert!(report.per_cut_reports[0].scaling.unwrap() < 0.0);
        assert!(report.per_cut_reports[1].scaling.unwrap() > 0.0);
    }

    #[test]
    fn single_limit_display_includes_per_cut_table() {
        let ir_limit = IrLimit::new_pure_soft(vec![EdgeIndex::from(1)]);
        let total_fit = PowerLawFit {
            exponent: -2.5,
            r_squared: 0.999,
        };
        let cut_reports = build_cut_limit_reports(
            ir_limit.num_soft(),
            vec![(
                0,
                Ok(PowerLawFit {
                    exponent: -4.5,
                    r_squared: 0.999,
                }),
            )],
        );

        let report = build_single_limit_report(&ir_limit, None, total_fit, cut_reports);
        let rendered = format!("{report}");

        assert!(rendered.contains("per-cut fits for"));
        assert!(rendered.contains("cut"));
        assert!(rendered.contains("r_squared"));
    }

    #[test]
    fn single_limit_display_includes_display_only_table() {
        let ir_limit = IrLimit::new_pure_soft(vec![EdgeIndex::from(1)]);
        let total_fit = PowerLawFit {
            exponent: -2.5,
            r_squared: 0.999,
        };
        let mut report = build_single_limit_report(&ir_limit, None, total_fit, Vec::new());
        report.display_only_reports = build_display_only_limit_reports(vec![
            (
                AdditionalWeightKey::Original,
                Ok(PowerLawFit {
                    exponent: -1.5,
                    r_squared: 0.995,
                }),
            ),
            (
                AdditionalWeightKey::ThresholdCounterterm { subset_index: 0 },
                Ok(PowerLawFit {
                    exponent: -0.5,
                    r_squared: 0.991,
                }),
            ),
        ]);

        let rendered = format!("{report}");

        assert!(rendered.contains("display-only fits for"));
        assert!(rendered.contains("original"));
        assert!(rendered.contains("ct_0"));
        assert!(!rendered.contains("note"));
    }

    #[test]
    fn graph_limit_display_weaves_per_cut_rows_into_limit_table() {
        let ir_limit = IrLimit::new_pure_soft(vec![EdgeIndex::from(1)]);
        let total_fit = PowerLawFit {
            exponent: -2.5,
            r_squared: 0.999,
        };
        let cut_reports = build_cut_limit_reports(
            ir_limit.num_soft(),
            vec![
                (
                    0,
                    Ok(PowerLawFit {
                        exponent: -4.5,
                        r_squared: 0.999,
                    }),
                ),
                (
                    1,
                    Ok(PowerLawFit {
                        exponent: -2.5,
                        r_squared: 0.999,
                    }),
                ),
            ],
        );

        let mut report = build_single_limit_report(&ir_limit, None, total_fit, cut_reports);
        report.display_only_reports = build_display_only_limit_reports(vec![
            (
                AdditionalWeightKey::Original,
                Ok(PowerLawFit {
                    exponent: -1.5,
                    r_squared: 0.995,
                }),
            ),
            (
                AdditionalWeightKey::ThresholdCounterterm { subset_index: 0 },
                Ok(PowerLawFit {
                    exponent: -0.5,
                    r_squared: 0.991,
                }),
            ),
        ]);
        let second_report = build_single_limit_report(
            &ir_limit,
            Some("ori-1".to_string()),
            PowerLawFit {
                exponent: -2.5,
                r_squared: 0.999,
            },
            build_cut_limit_reports(
                ir_limit.num_soft(),
                vec![(
                    2,
                    Ok(PowerLawFit {
                        exponent: -1.5,
                        r_squared: 0.995,
                    }),
                )],
            ),
        );
        let graph_report = GraphIRLimitReport {
            graph_name: "GL0".to_string(),
            all_limits_passed: true,
            cut_definitions: vec![
                GraphCutDefinition {
                    cut_id: 0,
                    edges: vec![EdgeIndex::from(1), EdgeIndex::from(3)],
                },
                GraphCutDefinition {
                    cut_id: 1,
                    edges: vec![EdgeIndex::from(2), EdgeIndex::from(4)],
                },
            ],
            single_limit_reports: vec![report, second_report],
        };
        let rendered = format!("{graph_report}");

        assert!(rendered.contains("cut definitions"));
        assert!(rendered.contains("edges"));
        assert!(rendered.contains("e1, e3"));
        assert!(rendered.find("cut definitions").unwrap() < rendered.find("status").unwrap());
        assert!(rendered.contains("item"));
        assert!(rendered.contains("cut 0"));
        assert!(rendered.contains("original"));
        assert!(rendered.contains("ct_0"));
        assert!(!rendered.contains("note"));
        assert!(rendered.contains("sum"));
        assert!(rendered.contains("INFO"));
        assert!(!rendered.contains("per-cut fits for"));
        assert!(rendered.matches('├').count() >= 3);
    }

    #[test]
    fn limit_data_extract_power_includes_display_only_component_fits() {
        let lambdas = constant_dropped_fit_points(
            &F::<f64>::from_f64(1.0e-3),
            &F::<f64>::from_f64(1.0e-1),
            5,
        )
        .unwrap();

        let data = lambdas
            .into_iter()
            .map(|lambda| {
                let lambda_f64 = lambda.into_ff64().0;
                let mut display_only_components = BTreeMap::new();
                display_only_components.insert(
                    AdditionalWeightKey::Original,
                    F::<ArbPrec>::from_f64(3.0 * lambda_f64.powf(-1.5) + 0.5),
                );
                display_only_components.insert(
                    AdditionalWeightKey::ThresholdCounterterm { subset_index: 0 },
                    F::<ArbPrec>::from_f64(5.0 * lambda_f64.powf(-0.5) + 1.0),
                );
                display_only_components.insert(
                    AdditionalWeightKey::FullMultiplicativeFactor,
                    F::<ArbPrec>::from_f64(7.0 * lambda_f64.powf(-2.5) + 1.0),
                );

                LambdaPointEval {
                    lambda,
                    value: ProfilePointValue {
                        total: F::<ArbPrec>::from_f64(2.0 * lambda_f64.powf(-2.0) + 1.0),
                        per_cut: BTreeMap::new(),
                        display_only_components,
                    },
                }
            })
            .collect();

        let (total_fit, _cut_fits, component_fits) = LimitData { data }.extract_power().unwrap();

        assert!((total_fit.exponent + 2.0).abs() < 1.0e-10);
        assert_eq!(component_fits.len(), 2);

        let component_fits = component_fits.into_iter().collect::<BTreeMap<_, _>>();
        let original_fit = component_fits
            .get(&AdditionalWeightKey::Original)
            .unwrap()
            .as_ref()
            .unwrap();
        let ct_fit = component_fits
            .get(&AdditionalWeightKey::ThresholdCounterterm { subset_index: 0 })
            .unwrap()
            .as_ref()
            .unwrap();

        assert!((original_fit.exponent + 1.5).abs() < 1.0e-10);
        assert!((ct_fit.exponent + 0.5).abs() < 1.0e-10);
    }
}
