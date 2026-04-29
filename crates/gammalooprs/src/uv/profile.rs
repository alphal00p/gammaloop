//! UV Profile Analysis Module
//!
//! This module provides functionality for analyzing ultraviolet behavior of loop integrands
//! by evaluating them at different momentum scalings and computing degrees of divergence.

use std::collections::BTreeMap;
use std::path::Path;
use std::sync::{Arc, Mutex};

use crate::DependentMomentaConstructor;
use crate::cff::expression::{GammaLoopGraphOrientation, OrientationData};
use crate::graph::Graph;
use crate::graph::parse::string_utils::ToOrderedSimple;
use crate::integrands::evaluation::EvaluationResult;
use crate::model::Model;
use crate::momentum::ThreeMomentum;
use crate::momentum::sample::{ExternalIndex, LoopIndex};
use crate::processes::{Amplitude, AmplitudeGraph};
use crate::settings::RuntimeSettings;
use crate::utils::F;
use crate::uv::UltravioletGraph;
use crate::{
    graph::LoopMomentumBasis,
    integrands::process::{
        OrientationProfileMode, ProcessIntegrand, evaluate_profile_momentum_point,
        orientation_labels_for_graph,
    },
};
use color_eyre::{Result, eyre::Context};
use colored::Colorize;
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::PowersetIterator;
use linnet::half_edge::involution::{EdgeIndex, SignOrZero};
use linnet::half_edge::subgraph::subset::SubSet;
use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike, SubSetOps};
use linnet::half_edge::tree::SimpleTraversalTree;
use rand::Rng;
use rayon::prelude::*;
use serde::Serialize;
use symbolica::domains::atom::AtomField;
use symbolica::numerical_integration::MonteCarloRng;
use symbolica::poly::series::Series;
use symbolica::symbol;
use tabled::{
    Table, Tabled,
    builder::Builder,
    settings::{Modify, Span, Style},
};
use tracing::{info_span, instrument};
use tracing_indicatif::{span_ext::IndicatifSpanExt, style::ProgressStyle};
use typed_index_collections::TiVec;

type ExternalMomenta = TiVec<ExternalIndex, ThreeMomentum<F<f64>>>;
type LoopMomentumSample = TiVec<LoopIndex, ThreeMomentum<F<f64>>>;
const UV_PROFILE_RETRY_MAX_DOD: f64 = -0.9;

pub struct ProfileSettings {
    pub n_points: usize,
    pub min_scale_exponent: f64,
    pub max_scale_exponent: f64,
    pub seed: u64,
    pub use_f128: bool,
    pub analyse_analytically: bool,
    pub orientation_mode: OrientationProfileMode,
}

impl Default for ProfileSettings {
    fn default() -> Self {
        ProfileSettings {
            n_points: 15,
            min_scale_exponent: 3.0,
            max_scale_exponent: 6.0,
            seed: 42,
            analyse_analytically: false,
            use_f128: false,
            orientation_mode: OrientationProfileMode::Summed,
        }
    }
}

pub fn logspace(start: f64, stop: f64, num: usize, base: f64) -> Vec<f64> {
    let log_start = start;
    let log_stop = stop;
    let step = (log_stop - log_start) / (num - 1) as f64;

    (0..num)
        .map(|i| {
            let exponent = log_start + step * i as f64;
            base.powf(exponent)
        })
        .collect()
}

fn lmb_seed(base_seed: u64, graph_id: usize, lmb_index: usize) -> u64 {
    base_seed
        .wrapping_add((graph_id as u64).wrapping_mul(0x9E3779B97F4A7C15))
        .wrapping_add(lmb_index as u64)
}

struct UVProfileRunner<'a> {
    integrand: &'a Arc<Mutex<ProcessIntegrand>>,
    scales: &'a [f64],
    externals: &'a ExternalMomenta,
    model: &'a Model,
    settings: &'a RuntimeSettings,
    profile_settings: &'a ProfileSettings,
    base_seed: u64,
}

pub trait UVProfileable {
    fn profile(
        &mut self,
        model: &Model,
        // settings: &RuntimeSettings,
        profile_settings: &ProfileSettings,
    ) -> Result<UVProfile>;
}

impl UVProfileable for Amplitude {
    #[instrument(skip_all)]
    fn profile(
        &mut self,
        model: &Model,
        // settings: &RuntimeSettings,
        profile_settings: &ProfileSettings,
    ) -> Result<UVProfile> {
        let scales = logspace(
            profile_settings.min_scale_exponent,
            profile_settings.max_scale_exponent,
            profile_settings.n_points,
            10.0,
        );

        let settings = self
            .integrand
            .as_ref()
            .ok_or(eyre!("Integrand Not built yet"))?
            .get_settings()
            .clone();
        let externals: ExternalMomenta = settings
            .kinematics
            .externals
            .get_dependent_externals(DependentMomentaConstructor::Amplitude(
                &self.external_signature,
            ))
            .unwrap()
            .into_iter()
            .map(|a| a.spatial)
            .collect();

        let base_seed = profile_settings.seed;
        let integrand = Arc::new(Mutex::new(self.integrand.take().unwrap()));

        let profile_span = info_span!("Profiling graphs", indicatif.pb_show = true);
        profile_span.pb_set_style(&ProgressStyle::with_template(
            "{wide_bar} {pos}/{len} {msg}",
        )?);
        profile_span.pb_set_length(self.graphs.len() as u64);
        profile_span.pb_set_message("Profiling graphs");
        profile_span.pb_set_finish_message("all graphs profiled");
        let _profile_span_enter = profile_span.enter();

        let runner = UVProfileRunner {
            integrand: &integrand,
            scales: &scales,
            externals: &externals,
            model,
            settings: &settings,
            profile_settings,
            base_seed,
        };

        let per_graph = self
            .graphs
            .par_iter()
            .enumerate()
            .map(|(i, g)| {
                let res = runner.sample_graph(i, g)?;
                profile_span.pb_inc(1);
                Ok(res)
            })
            .collect::<Result<Vec<_>>>()?;

        drop(_profile_span_enter);
        drop(profile_span);

        let integrand = Arc::try_unwrap(integrand)
            .map_err(|_| color_eyre::eyre::eyre!("integrand still shared"))?
            .into_inner()
            .expect("integrand mutex poisoned");
        self.integrand = Some(integrand);

        Ok(UVProfile { per_graph, scales })
        // results.push((inspect_res, analytic_res));
    }
}

pub struct UVProfile {
    pub per_graph: Vec<UVSamplingResult>,
    pub scales: Vec<f64>,
}

impl UVProfile {
    pub fn analyse(&self) -> UVProfileAnalysis {
        let graphs = self
            .per_graph
            .iter()
            .enumerate()
            .map(|(graph_index, graph)| {
                let lmbs: Vec<UVProfileLmbAnalysis> = graph
                    .per_lmb
                    .iter()
                    .enumerate()
                    .map(|(lmb_index, lmb)| {
                        let lmb_label = lmb_label(&lmb.lmb);
                        let subsets: Vec<UVProfileSubsetAnalysis> = lmb
                            .per_subsets
                            .iter()
                            .enumerate()
                            .map(|(subset_index, (subset, subset_result))| {
                                let mut not_included: SubSet<LoopIndex> =
                                    SubSet::full(subset.size());
                                not_included.subtract_with(subset);
                                let free: Vec<EdgeIndex> = subset
                                    .included_iter()
                                    .map(|loop_index| lmb.lmb.loop_edges[loop_index])
                                    .collect();

                                let fixed: Vec<EdgeIndex> = not_included
                                    .included_iter()
                                    .map(|loop_index| lmb.lmb.loop_edges[loop_index])
                                    .collect();
                                let analysis = subset_result.analyse(&self.scales);
                                let per_orientation_inspect_entries =
                                    analysis.per_orientation_inspect_entries();
                                let analytic_entries =
                                    analysis.analytic.as_ref().and_then(|analytic| {
                                        let entries = analytic
                                            .per_orientation
                                            .iter()
                                            .map(|(orientation, orientation_analysis)| {
                                                let (orientation_edges, orientation_signs) =
                                                    orientation_signs(orientation);
                                                UVProfileAnalyticEntry {
                                                    graph_index,
                                                    lmb_index,
                                                    subset_index,
                                                    fixed: fixed.clone(),
                                                    free: free.clone(),
                                                    orientation_edges,
                                                    orientation_signs,
                                                    is_constant: orientation_analysis.is_constant,
                                                    leading_coef: orientation_analysis
                                                        .leading_coef
                                                        .to_string(),
                                                }
                                            })
                                            .collect::<Vec<_>>();
                                        if entries.is_empty() {
                                            None
                                        } else {
                                            Some(entries)
                                        }
                                    });
                                UVProfileSubsetAnalysis {
                                    subset_index,
                                    fixed,
                                    free,
                                    initial_dod: subset_result.initial_dod,
                                    analysis,
                                    per_orientation_inspect_entries,
                                    analytic_entries,
                                }
                            })
                            .collect();

                        UVProfileLmbAnalysis {
                            lmb_index,
                            lmb_label,
                            subsets,
                        }
                    })
                    .collect();

                UVProfileGraphAnalysis { graph_index, lmbs }
            })
            .collect();

        UVProfileAnalysis {
            scales: self.scales.clone(),
            graphs,
        }
    }

    pub fn write_profile_data<P: AsRef<Path>>(
        &self,
        _settings: &ProfileSettings,
        out_dir: P,
    ) -> Result<()> {
        self.analyse().write_profile_data(out_dir)
    }

    pub fn write_typst_bundle<P: AsRef<Path>>(
        &self,
        settings: &ProfileSettings,
        out_dir: P,
    ) -> Result<()> {
        self.write_profile_data(settings, out_dir)
    }

    pub fn pass_fail(&self, max_dod: f64, _settings: &ProfileSettings) -> UVProfilePassFail {
        self.analyse().pass_fail(max_dod)
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct UVProfileAnalysis {
    pub scales: Vec<f64>,
    pub graphs: Vec<UVProfileGraphAnalysis>,
}

#[derive(Debug, Clone, Serialize)]
pub struct UVProfileGraphAnalysis {
    pub graph_index: usize,
    pub lmbs: Vec<UVProfileLmbAnalysis>,
}

#[derive(Debug, Clone, Serialize)]
pub struct UVProfileLmbAnalysis {
    pub lmb_index: usize,
    pub lmb_label: String,
    pub subsets: Vec<UVProfileSubsetAnalysis>,
}

#[derive(Debug, Clone, Serialize)]
pub struct UVProfileSubsetAnalysis {
    pub subset_index: usize,
    pub fixed: Vec<EdgeIndex>,
    pub free: Vec<EdgeIndex>,
    pub initial_dod: i32,
    // pub subset_label: String,
    pub analysis: Analysis,
    pub per_orientation_inspect_entries: Option<Vec<UVProfileOrientationInspectEntry>>,
    pub analytic_entries: Option<Vec<UVProfileAnalyticEntry>>,
}

impl UVProfileSubsetAnalysis {
    pub fn estimated_dod(&self) -> Option<i64> {
        self.analysis
            .inspect_level
            .as_ref()
            .map(|analysis| analysis.estimated_dod)
    }

    pub fn bare_dod_matches_estimate(&self) -> bool {
        self.estimated_dod() == Some(i64::from(self.initial_dod))
    }
}

// #[derive(Debug, Clone, Serialize)]
// pub struct UVProfilePoint {
//     norm:f64,
// }

#[derive(Debug, Clone, Serialize)]
pub struct UVProfileSummary {
    pub subset_count: usize,
    pub slope_min: Option<f64>,
    pub slope_max: Option<f64>,
    pub dod_min: Option<i64>,
    pub dod_max: Option<i64>,
    pub stable_min: Option<f64>,
    pub stable_max: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct UVProfilePassFail {
    pub max_dod: f64,
    pub total: usize,
    pub failed: usize,
    pub failures: Vec<UVProfileFailure>,
}

#[derive(Debug, Clone, Serialize)]
pub struct UVProfileFailure {
    pub graph_index: usize,
    pub lmb_index: usize,
    pub fixed: Vec<EdgeIndex>,
    pub free: Vec<EdgeIndex>,
    pub orientation_label: Option<String>,
    pub reason: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct UVProfileAnalyticEntry {
    pub graph_index: usize,
    pub lmb_index: usize,
    pub subset_index: usize,
    pub fixed: Vec<EdgeIndex>,
    pub free: Vec<EdgeIndex>,
    pub orientation_edges: Vec<String>,
    pub orientation_signs: Vec<String>,
    pub is_constant: bool,
    pub leading_coef: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct UVProfileOrientationInspectEntry {
    pub orientation_label: String,
    pub analysis: Option<InspectAnalysis>,
}

#[derive(Tabled)]
struct UVProfileSubsetRow {
    #[tabled(rename = "fixed")]
    fixed: String,
    #[tabled(rename = "→ ∞")]
    free: String,
    #[tabled(rename = "slope")]
    slope: String,
    #[tabled(rename = "r2")]
    r_squared: String,
    #[tabled(rename = "DOD")]
    estimated_dod: String,
    #[tabled(rename = "bare DOD")]
    initial_dod: String,
    #[tabled(rename = "inspect")]
    inspect: String,
}

#[derive(Tabled)]
struct UVProfileOrientationSubsetRow {
    #[tabled(rename = "fixed")]
    fixed: String,
    #[tabled(rename = "→ ∞")]
    free: String,
    #[tabled(rename = "orientation")]
    orientation_label: String,
    #[tabled(rename = "slope")]
    slope: String,
    #[tabled(rename = "r2")]
    r_squared: String,
    #[tabled(rename = "DOD")]
    estimated_dod: String,
    #[tabled(rename = "bare DOD")]
    initial_dod: String,
    #[tabled(rename = "inspect")]
    inspect: String,
}

#[derive(Debug, Clone, Serialize)]
struct FitResult {
    slope: f64,
    points: Vec<f64>,
    points_detail: Vec<EvaluationResult>,
    intercept: f64,
    r_squared: f64,
}

impl UVProfileAnalysis {
    pub fn tables_per_graph(&self, max_dod: f64) -> Vec<Table> {
        self.graphs
            .iter()
            .map(|graph| {
                let rows = graph
                    .lmbs
                    .iter()
                    .flat_map(|lmb| {
                        lmb.subsets.iter().map(|subset| {
                            let (slope, r_squared, estimated_dod) =
                                match &subset.analysis.inspect_level {
                                    Some(analysis) => {
                                        let r2_text = format!("{:.3}", analysis.result.r_squared);
                                        let r2_text = if analysis.result.r_squared >= 0.99 {
                                            r2_text.green()
                                        } else {
                                            r2_text.red()
                                        }
                                        .to_string();

                                        let dod = analysis.estimated_dod;
                                        let dod_text = dod.to_string();
                                        let dod_text = if (dod as f64) <= max_dod {
                                            dod_text.green()
                                        } else {
                                            dod_text.red()
                                        }
                                        .to_string();

                                        (format!("{:.6}", analysis.result.slope), r2_text, dod_text)
                                    }
                                    None => ("-".to_string(), "-".to_string(), "-".to_string()),
                                };

                            UVProfileSubsetRow {
                                fixed: format!(
                                    "{{{}}}",
                                    subset.fixed.iter().map(ToString::to_string).join(",")
                                ),
                                free: format!(
                                    "{{{}}}",
                                    subset.free.iter().map(ToString::to_string).join(",")
                                ),
                                slope,
                                r_squared,
                                estimated_dod,
                                initial_dod: if subset.initial_dod >= 0 {
                                    subset.initial_dod.to_string().red().to_string()
                                } else {
                                    subset.initial_dod.to_string().green().to_string()
                                },
                                inspect: inspect_retry_label(
                                    subset.analysis.inspect_level.as_ref(),
                                ),
                            }
                        })
                    })
                    .collect::<Vec<_>>();

                let mut table = Table::new(rows);
                table.with(Style::rounded());
                table
            })
            .collect()
    }

    pub fn analytic_tables_per_graph(&self) -> Vec<Option<Table>> {
        self.graphs
            .iter()
            .map(|graph| {
                let mut groups: Vec<Vec<&UVProfileAnalyticEntry>> = Vec::new();
                let mut orientation_headers: Option<Vec<String>> = None;

                for lmb in &graph.lmbs {
                    for subset in &lmb.subsets {
                        if let Some(entries) = &subset.analytic_entries {
                            if orientation_headers.is_none() {
                                orientation_headers =
                                    entries.first().map(|entry| entry.orientation_edges.clone());
                            }
                            groups.push(entries.iter().collect());
                        }
                    }
                }

                if groups.is_empty() {
                    return None;
                }

                let mut builder = Builder::new();
                let mut header = vec!["fixed".to_string(), "→ ∞".to_string()];
                header.extend(orientation_headers.unwrap_or_default());
                header.extend(["const".to_string(), "leading coef".to_string()]);
                builder.push_record(header);

                let mut span_ops: Vec<(usize, usize, usize)> = Vec::new();
                let mut row_index = 1;

                for group in groups {
                    let span_len = group.len();
                    let start_row = row_index;

                    for (i, entry) in group.into_iter().enumerate() {
                        let fixed = if i == 0 {
                            format!(
                                "{{{}}}",
                                entry.fixed.iter().map(ToString::to_string).join(",")
                            )
                        } else {
                            String::new()
                        };
                        let free = if i == 0 {
                            format!(
                                "{{{}}}",
                                entry.free.iter().map(ToString::to_string).join(",")
                            )
                        } else {
                            String::new()
                        };
                        let mut row = vec![fixed, free];
                        row.extend(entry.orientation_signs.iter().cloned());
                        row.extend([entry.is_constant.to_string(), "".to_string()]);
                        builder.push_record(row);
                        row_index += 1;
                    }

                    if span_len > 1 {
                        span_ops.push((start_row, 0, span_len));
                        span_ops.push((start_row, 1, span_len));
                    }
                }

                let mut table = builder.build();
                for (row, col, span_len) in span_ops {
                    table.with(Modify::new((row, col)).with(Span::row(span_len as isize)));
                }
                table.with(Style::rounded());
                Some(table)
            })
            .collect()
    }

    pub fn per_orientation_tables_per_graph(&self, max_dod: f64) -> Vec<Option<Table>> {
        self.graphs
            .iter()
            .map(|graph| {
                let rows = graph
                    .lmbs
                    .iter()
                    .flat_map(|lmb| {
                        lmb.subsets.iter().flat_map(|subset| {
                            subset
                                .per_orientation_inspect_entries
                                .iter()
                                .flatten()
                                .map(|entry| {
                                    let (slope, r_squared, estimated_dod) = match &entry.analysis {
                                        Some(analysis) => {
                                            let r2_text =
                                                format!("{:.3}", analysis.result.r_squared);
                                            let r2_text = if analysis.result.r_squared >= 0.99 {
                                                r2_text.green()
                                            } else {
                                                r2_text.red()
                                            }
                                            .to_string();

                                            let dod = analysis.estimated_dod;
                                            let dod_text = dod.to_string();
                                            let dod_text = if (dod as f64) <= max_dod {
                                                dod_text.green()
                                            } else {
                                                dod_text.red()
                                            }
                                            .to_string();

                                            (
                                                format!("{:.6}", analysis.result.slope),
                                                r2_text,
                                                dod_text,
                                            )
                                        }
                                        None => ("-".to_string(), "-".to_string(), "-".to_string()),
                                    };

                                    UVProfileOrientationSubsetRow {
                                        fixed: format!(
                                            "{{{}}}",
                                            subset.fixed.iter().map(ToString::to_string).join(",")
                                        ),
                                        free: format!(
                                            "{{{}}}",
                                            subset.free.iter().map(ToString::to_string).join(",")
                                        ),
                                        orientation_label: entry.orientation_label.clone(),
                                        slope,
                                        r_squared,
                                        estimated_dod,
                                        initial_dod: if subset.initial_dod >= 0 {
                                            subset.initial_dod.to_string().red().to_string()
                                        } else {
                                            subset.initial_dod.to_string().green().to_string()
                                        },
                                        inspect: inspect_retry_label(entry.analysis.as_ref()),
                                    }
                                })
                                .collect::<Vec<_>>()
                        })
                    })
                    .collect::<Vec<_>>();

                if rows.is_empty() {
                    None
                } else {
                    let mut table = Table::new(rows);
                    table.with(Style::rounded());
                    Some(table)
                }
            })
            .collect()
    }

    pub fn write_profile_data<P: AsRef<Path>>(&self, out_dir: P) -> Result<()> {
        let out_dir = out_dir.as_ref();
        std::fs::create_dir_all(out_dir).context("failed to create UV profile output directory")?;

        let json_path = out_dir.join("uv_profile.json");
        let json = serde_json::to_string_pretty(self)
            .context("failed to serialize UV profile analysis to JSON")?;
        std::fs::write(&json_path, json).context("failed to write UV profile JSON output")?;

        Ok(())
    }

    pub fn pass_fail(&self, max_dod: f64) -> UVProfilePassFail {
        let mut failures = Vec::new();
        let mut total = 0;

        for graph in &self.graphs {
            for lmb in &graph.lmbs {
                for subset in &lmb.subsets {
                    let per_orientation_reason_count = subset
                        .per_orientation_inspect_entries
                        .as_ref()
                        .map(|entries| {
                            entries
                                .iter()
                                .filter(|entry| match &entry.analysis {
                                    None => true,
                                    Some(analysis) => {
                                        analysis.result.slope > max_dod
                                            || analysis.result.slope.is_nan()
                                    }
                                })
                                .count()
                        })
                        .unwrap_or(0);
                    let all_orientations_pass = subset.per_orientation_inspect_entries.is_some()
                        && per_orientation_reason_count == 0;

                    total += 1;
                    let reason = match &subset.analysis.inspect_level {
                        None => Some("missing_fit"),
                        Some(analysis)
                            if analysis.result.slope > max_dod
                                || analysis.result.slope.is_nan() =>
                        {
                            Some("dod_exceeds_threshold")
                        }
                        _ => None,
                    };

                    if let Some(reason) = reason.filter(|_| !all_orientations_pass) {
                        failures.push(UVProfileFailure {
                            graph_index: graph.graph_index,
                            lmb_index: lmb.lmb_index,
                            fixed: subset.fixed.clone(),
                            free: subset.free.clone(),
                            orientation_label: None,
                            reason: reason.to_string(),
                        });
                    }

                    for entry in subset.per_orientation_inspect_entries.iter().flatten() {
                        total += 1;
                        let reason = match &entry.analysis {
                            None => Some("missing_fit"),
                            Some(analysis)
                                if analysis.result.slope > max_dod
                                    || analysis.result.slope.is_nan() =>
                            {
                                Some("dod_exceeds_threshold")
                            }
                            _ => None,
                        };

                        if let Some(reason) = reason {
                            failures.push(UVProfileFailure {
                                graph_index: graph.graph_index,
                                lmb_index: lmb.lmb_index,
                                fixed: subset.fixed.clone(),
                                free: subset.free.clone(),
                                orientation_label: Some(entry.orientation_label.clone()),
                                reason: reason.to_string(),
                            });
                        }
                    }
                }
            }
        }

        UVProfilePassFail {
            max_dod,
            total,
            failed: failures.len(),
            failures,
        }
    }
}

fn lmb_label(lmb: &LoopMomentumBasis) -> String {
    let edges: Vec<String> = lmb.loop_edges.iter().map(|edge| edge.to_string()).collect();
    format!("loop_edges=[{}]", edges.join(","))
}

fn orientation_signs(orientation: &OrientationData) -> (Vec<String>, Vec<String>) {
    let edges = orientation
        .orientation
        .iter()
        .map(|(edge, _)| edge.to_string())
        .collect::<Vec<_>>();
    let signs = orientation
        .orientation
        .iter()
        .map(|(_, sign)| SignOrZero::from(*sign).to_string())
        .collect::<Vec<_>>();
    (edges, signs)
}

pub struct UVSamplingResult {
    pub per_lmb: Vec<LMBResult>,
}

impl<'a> UVProfileRunner<'a> {
    fn sample_graph(&self, graph_id: usize, g: &AmplitudeGraph) -> Result<UVSamplingResult> {
        let lmbs = g.derived_data.lmbs.as_ref().unwrap();
        let integrand_expr = &g.derived_data.all_mighty_integrand;
        let analytic_orientations: Vec<_> = g
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .orientations
            .iter()
            .collect();
        let orientation_labels = if self
            .profile_settings
            .orientation_mode
            .profiles_per_orientation()
        {
            let integrand = self.integrand.lock().expect("integrand mutex poisoned");
            match &*integrand {
                ProcessIntegrand::Amplitude(amplitude) => {
                    Some(orientation_labels_for_graph(amplitude, graph_id)?)
                }
                ProcessIntegrand::CrossSection(_) => {
                    unreachable!("UV profiling expects amplitudes")
                }
            }
        } else {
            None
        };
        let lmb_refs: Vec<_> = lmbs.iter().enumerate().collect();

        let lmb_span = info_span!(
            "Profiling loop momentum bases",
            indicatif.pb_show = true,
            graph_id = graph_id
        );
        lmb_span.pb_set_style(
            &ProgressStyle::with_template("{wide_bar} {pos}/{len} {msg}")
                .expect("invalid progress bar template"),
        );
        lmb_span.pb_set_length(lmb_refs.len() as u64);
        lmb_span.pb_set_message("Profiling loop momentum bases");
        lmb_span.pb_set_finish_message("all loop momentum bases profiled");
        let _lmb_span_enter = lmb_span.enter();

        let per_lmb = lmb_refs
            .par_iter()
            .map(|(lmb_index, lmb)| {
                let mut res = self.sample_lmb(
                    graph_id,
                    &g.graph,
                    *lmb_index,
                    lmb,
                    orientation_labels.as_deref(),
                )?;

                if self.profile_settings.analyse_analytically {
                    let orientation_limits: Vec<(
                        SubSet<LoopIndex>,
                        OrientationData,
                        Series<AtomField>,
                    )> = analytic_orientations
                        .par_iter()
                        .map(|o| {
                            let oatom = o.data.orientation.select_gs(integrand_expr);
                            g.graph
                                .all_limits(&g.graph.full_filter(), &oatom, symbol!("lambd"), lmb)
                                .into_iter()
                                .map(|(l, v)| (l, o.data.clone(), v))
                                .collect::<Vec<_>>()
                        })
                        .reduce(Vec::new, |mut acc, mut v| {
                            acc.append(&mut v);
                            acc
                        });

                    for (l, odata, v) in orientation_limits {
                        let subset = res
                            .per_subsets
                            .get_mut(&l)
                            .expect("subset missing for orientation limits");
                        let analytic = subset.analytic.get_or_insert_with(|| AnalyticResult {
                            per_orientations: BTreeMap::new(),
                        });
                        analytic.per_orientations.insert(odata, v);
                    }
                }
                lmb_span.pb_inc(1);
                Ok(res)
            })
            .collect::<Result<Vec<_>>>()?;

        drop(_lmb_span_enter);
        drop(lmb_span);

        Ok(UVSamplingResult { per_lmb })
    }
}

pub struct LMBResult {
    pub(crate) lmb: LoopMomentumBasis,
    pub(crate) per_subsets: BTreeMap<SubSet<LoopIndex>, SubSetResult>,
}

impl<'a> UVProfileRunner<'a> {
    fn sample_lmb(
        &self,
        graph_id: usize,
        graph: &Graph,
        lmb_index: usize,
        lmb: &LoopMomentumBasis,
        orientation_labels: Option<&[String]>,
    ) -> Result<LMBResult> {
        let mut rng = MonteCarloRng::new(lmb_seed(self.base_seed, graph_id, lmb_index), 0);
        let sample: LoopMomentumSample = lmb
            .loop_edges
            .iter()
            .map(|_| ThreeMomentum {
                px: F(
                    rng.random_range(-self.settings.kinematics.e_cm..self.settings.kinematics.e_cm)
                ),
                py: F(
                    rng.random_range(-self.settings.kinematics.e_cm..self.settings.kinematics.e_cm)
                ),
                pz: F(
                    rng.random_range(-self.settings.kinematics.e_cm..self.settings.kinematics.e_cm)
                ),
            })
            .collect();

        let mut loops = PowersetIterator::<LoopIndex>::new(lmb.loop_edges.len() as u8);
        loops.next();

        let subsets: Vec<_> = loops.collect();
        let subset_span = info_span!(
            "Profiling subsets",
            indicatif.pb_show = true,
            graph_id = graph_id,
            lmb_index = lmb_index
        );
        subset_span.pb_set_style(
            &ProgressStyle::with_template("{wide_bar} {pos}/{len} {msg}")
                .expect("invalid progress bar template"),
        );
        subset_span.pb_set_length(subsets.len() as u64);
        subset_span.pb_set_message("Profiling subsets");
        let _subset_span_enter = subset_span.enter();

        let per_subsets_vec: Vec<(SubSet<LoopIndex>, SubSetResult)> = subsets
            .into_par_iter()
            .map_init(
                || {
                    self.integrand
                        .lock()
                        .expect("integrand mutex poisoned")
                        .clone()
                },
                |integrand, ls| {
                    let res = self.sample_subset(
                        integrand,
                        SubsetSampleInput {
                            graph_id,
                            graph,
                            subset: &ls,
                            lmb,
                            sample: &sample,
                            orientation_labels,
                        },
                    )?;
                    subset_span.pb_inc(1);
                    Ok((ls, res))
                },
            )
            .collect::<Result<Vec<_>>>()?;
        let per_subsets = per_subsets_vec.into_iter().collect();

        drop(_subset_span_enter);
        drop(subset_span);

        Ok(LMBResult {
            lmb: lmb.clone(),
            per_subsets,
        })
    }
}

pub struct SubSetResult {
    pub(crate) initial_dod: i32,
    pub(crate) inspect: InspectSamples,
    pub(crate) analytic: Option<AnalyticResult>,
}

#[derive(Debug, Clone, Default)]
pub(crate) struct InspectSamples {
    pub(crate) summed: Vec<InspectResult>,
    pub(crate) summed_used_arb_prec_retry: bool,
    pub(crate) per_orientation: Vec<OrientationInspectSamples>,
}

#[derive(Debug, Clone)]
pub(crate) struct OrientationInspectSamples {
    pub(crate) label: String,
    pub(crate) inspect: Vec<InspectResult>,
    pub(crate) used_arb_prec_retry: bool,
}

#[derive(Debug, Clone)]
struct InspectRun {
    inspect: Vec<InspectResult>,
    used_arb_prec_retry: bool,
}

#[derive(Clone, Copy)]
struct SubsetOrientationInput<'a> {
    graph_id: usize,
    subset: &'a SubSet<LoopIndex>,
    lmb: &'a LoopMomentumBasis,
    sample: &'a LoopMomentumSample,
    orientation: Option<usize>,
}

struct SubsetSampleInput<'a> {
    graph_id: usize,
    graph: &'a Graph,
    subset: &'a SubSet<LoopIndex>,
    lmb: &'a LoopMomentumBasis,
    sample: &'a LoopMomentumSample,
    orientation_labels: Option<&'a [String]>,
}

impl<'a> UVProfileRunner<'a> {
    fn sample_subset(
        &self,
        integrand: &mut ProcessIntegrand,
        input: SubsetSampleInput<'_>,
    ) -> Result<SubSetResult> {
        let mut subgraph: SuBitGraph = input.graph.empty_subgraph();
        for l in input.subset.included_iter() {
            let eid = input.lmb.loop_edges[l];
            let cut = input.graph[&eid].1.any_hedge();
            let root_node = input.graph.node_id(cut);

            let tree = SimpleTraversalTree::depth_first_traverse(
                input.graph,
                &input.lmb.tree,
                &root_node,
                None,
            )
            .unwrap();
            subgraph.union_with(
                &tree
                    .get_cycle(cut, input.graph.underlying.as_ref())
                    .unwrap()
                    .filter,
            );
        }

        let initial_dod = input.graph.dod(&subgraph);
        let inspect = if let Some(orientation_labels) = input.orientation_labels {
            let per_orientation = orientation_labels
                .iter()
                .enumerate()
                .map(|(orientation_id, label)| {
                    let inspect = self.sample_subset_orientation(
                        integrand,
                        SubsetOrientationInput {
                            graph_id: input.graph_id,
                            subset: input.subset,
                            lmb: input.lmb,
                            sample: input.sample,
                            orientation: Some(orientation_id),
                        },
                    )?;
                    Ok(OrientationInspectSamples {
                        label: label.clone(),
                        inspect: inspect.inspect,
                        used_arb_prec_retry: inspect.used_arb_prec_retry,
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            let summed = sum_orientation_inspect_samples(&per_orientation);
            InspectSamples {
                summed,
                summed_used_arb_prec_retry: per_orientation
                    .iter()
                    .any(|orientation| orientation.used_arb_prec_retry),
                per_orientation,
            }
        } else {
            let inspect = self.sample_subset_orientation(
                integrand,
                SubsetOrientationInput {
                    graph_id: input.graph_id,
                    subset: input.subset,
                    lmb: input.lmb,
                    sample: input.sample,
                    orientation: None,
                },
            )?;
            InspectSamples {
                summed: inspect.inspect,
                summed_used_arb_prec_retry: inspect.used_arb_prec_retry,
                per_orientation: Vec::new(),
            }
        };
        let analytic = None;

        Ok(SubSetResult {
            inspect,
            initial_dod,
            analytic,
        })
    }

    fn sample_subset_orientation(
        &self,
        integrand: &mut ProcessIntegrand,
        input: SubsetOrientationInput<'_>,
    ) -> Result<InspectRun> {
        let inspect = self.sample_subset_orientation_with_precision(
            integrand,
            input,
            self.profile_settings.use_f128,
        )?;
        if !self.profile_settings.use_f128
            && inspect_results_need_arbprec_retry(&inspect, self.scales)
        {
            return Ok(InspectRun {
                inspect: self.sample_subset_orientation_with_precision(integrand, input, true)?,
                used_arb_prec_retry: true,
            });
        }
        Ok(InspectRun {
            inspect,
            used_arb_prec_retry: false,
        })
    }

    fn sample_subset_orientation_with_precision(
        &self,
        integrand: &mut ProcessIntegrand,
        input: SubsetOrientationInput<'_>,
        use_arb_prec: bool,
    ) -> Result<Vec<InspectResult>> {
        let n_included = input.subset.n_included() as i32;
        self.scales
            .iter()
            .map(|s| {
                let prefactor = s.powi(3 * n_included);
                let mut scaled_sample = input.sample.clone();
                for l in input.subset.included_iter() {
                    scaled_sample[l] = scaled_sample[l].map_ref(&|a| a * F(*s));
                }
                let loop_momenta = input
                    .lmb
                    .loop_edges
                    .iter()
                    .map(|edge| {
                        input.lmb.edge_signatures[*edge]
                            .compute_momentum(&scaled_sample, self.externals)
                    })
                    .collect::<Vec<_>>();

                let inspect_res_eval = evaluate_momentum_space_point(
                    integrand,
                    self.model,
                    loop_momenta,
                    input.graph_id,
                    input.orientation,
                    use_arb_prec,
                )?;

                Ok(InspectResult {
                    result: inspect_res_eval,
                    prefactor,
                })
            })
            .collect()
    }
}

impl SubSetResult {
    pub fn analyse_inspect(&self, scales: &[f64]) -> Option<InspectAnalysis> {
        analyse_inspect_results(
            &self.inspect.summed,
            scales,
            self.inspect.summed_used_arb_prec_retry,
        )
    }

    pub fn analyse(&self, scales: &[f64]) -> Analysis {
        Analysis {
            inspect_level: self.analyse_inspect(scales),
            per_orientation_inspect: (!self.inspect.per_orientation.is_empty()).then(|| {
                self.inspect
                    .per_orientation
                    .iter()
                    .map(|orientation| OrientationInspectAnalysis {
                        orientation_label: orientation.label.clone(),
                        analysis: analyse_inspect_results(
                            &orientation.inspect,
                            scales,
                            orientation.used_arb_prec_retry,
                        ),
                    })
                    .collect()
            }),
            analytic: self.analyse_analytic(),
        }
    }

    pub fn analyse_analytic(&self) -> Option<AnalyticAnalysis> {
        //             .derived_data
        //             .cff_expression
        //             .as_ref()
        //             .unwrap()
        //             .orientations
        //             .iter()
        //             .enumerate()
        //         {
        //             let expansion = symbol!("lambd");

        //             for (ls, res) in &analytic_res[i_lmb][i] {
        //                 // print!("res:{res}");
        //                 let l = res.coefficient_list::<i8>(&[Atom::var(expansion)]);

        //                 println!(
        //                     "In the limit of {:?} going to infinity for orientation \n{}:",
        //                     ls.included_iter()
        //                         .map(|l| lmb.loop_edges[l].to_string())
        //                         .collect::<Vec<_>>(),
        //                     o.data
        //                 );
        //                 if l.is_empty() {
        //                     println!("\tFull cancellation to order 1");
        //                 }
        //                 for (t, a) in l {
        //                     println!("\t{}: {}", t, a);
        //                 }
        //             }
        //         }

        Some(AnalyticAnalysis {
            per_orientation: self
                .analytic
                .as_ref()?
                .per_orientations
                .par_iter()
                .map(|(k, v)| {
                    (
                        k.clone(),
                        OrientationAnalyticAnalysis {
                            is_constant: v.is_constant(),
                            leading_coef: v.lcoeff().to_ordered_simple(),
                        },
                    )
                })
                .collect(),
        })
    }
}

fn analyse_inspect_results(
    inspect: &[InspectResult],
    scales: &[f64],
    used_arb_prec_retry: bool,
) -> Option<InspectAnalysis> {
    let result = log_log_slope(inspect, scales)?;
    let dod = result.slope.round() as i64;
    Some(InspectAnalysis {
        result,
        estimated_dod: dod,
        used_arb_prec_retry,
    })
}

fn inspect_results_need_arbprec_retry(inspect: &[InspectResult], scales: &[f64]) -> bool {
    match analyse_inspect_results(inspect, scales, false) {
        None => true,
        Some(analysis) => {
            analysis.result.slope.is_nan() || analysis.result.slope > UV_PROFILE_RETRY_MAX_DOD
        }
    }
}

fn log_log_slope(inspect: &[InspectResult], scales: &[f64]) -> Option<FitResult> {
    let mut valid_samples = Vec::new();
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_xy = 0.0;
    let mut sum_x2 = 0.0;
    let mut points = vec![];
    let mut points_detail = vec![];

    for (x, s) in inspect.iter().zip(scales) {
        let norm = x.magnitude();
        if norm <= 0.0 {
            println!("{s}:\t{}", x.result.evaluation_metadata);
            continue;
        }
        if !norm.is_finite() {
            continue;
        }
        points_detail.push(x.result.clone());
        points.push(norm);
        let y = (norm).log10();
        let x = s.log10();
        if !y.is_finite() || !x.is_finite() {
            continue;
        }
        valid_samples.push((x, y));
        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_x2 += x * x;
    }

    if valid_samples.len() < 2 {
        return None;
    }

    let n = valid_samples.len() as f64;
    let denominator = n * sum_x2 - sum_x * sum_x;
    if denominator.abs() < 1e-15 {
        return None;
    }

    let slope = (n * sum_xy - sum_x * sum_y) / denominator;
    let intercept = (sum_y - slope * sum_x) / n;

    let y_mean = sum_y / n;
    let mut ss_tot = 0.0;
    let mut ss_res = 0.0;
    for (x, y) in &valid_samples {
        let y_pred = intercept + slope * x;
        ss_tot += (y - y_mean).powi(2);
        ss_res += (y - y_pred).powi(2);
    }

    let r_squared = if ss_tot > 1e-15 {
        1.0 - ss_res / ss_tot
    } else {
        0.0
    };

    Some(FitResult {
        points,
        points_detail,
        slope,
        intercept,
        r_squared,
    })
}

fn sum_orientation_inspect_samples(
    per_orientation: &[OrientationInspectSamples],
) -> Vec<InspectResult> {
    let Some((first, rest)) = per_orientation.split_first() else {
        return Vec::new();
    };

    (0..first.inspect.len())
        .map(|point_index| {
            let mut summed = first.inspect[point_index].clone();
            for orientation in rest {
                summed.result.integrand_result +=
                    orientation.inspect[point_index].result.integrand_result;
            }
            summed
        })
        .collect()
}

fn inspect_retry_label(analysis: Option<&InspectAnalysis>) -> String {
    match analysis {
        Some(analysis) if analysis.used_arb_prec_retry => "arb retry".yellow().to_string(),
        Some(_) => String::new(),
        None => "-".to_string(),
    }
}

fn evaluate_momentum_space_point(
    integrand: &mut ProcessIntegrand,
    model: &Model,
    loop_momenta: Vec<ThreeMomentum<F<f64>>>,
    graph_id: usize,
    orientation: Option<usize>,
    use_arb_prec: bool,
) -> Result<EvaluationResult> {
    match integrand {
        ProcessIntegrand::Amplitude(amplitude) => evaluate_profile_momentum_point(
            amplitude,
            model,
            graph_id,
            orientation,
            loop_momenta,
            use_arb_prec,
        ),
        ProcessIntegrand::CrossSection(cross_section) => evaluate_profile_momentum_point(
            cross_section,
            model,
            graph_id,
            orientation,
            loop_momenta,
            use_arb_prec,
        ),
    }
}

#[derive(Debug, Clone)]
pub struct InspectResult {
    pub(crate) result: EvaluationResult,
    pub(crate) prefactor: f64,
}

impl InspectResult {
    fn magnitude(&self) -> f64 {
        self.result.integrand_result.norm_squared().sqrt().0 * self.prefactor
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct Analysis {
    ///Is None if the fit hasn't worked
    inspect_level: Option<InspectAnalysis>,
    #[serde(skip_serializing)]
    per_orientation_inspect: Option<Vec<OrientationInspectAnalysis>>,
    ///Is None if the analytic analysis is disabled
    #[serde(skip_serializing)]
    analytic: Option<AnalyticAnalysis>,
}

impl Analysis {
    fn per_orientation_inspect_entries(&self) -> Option<Vec<UVProfileOrientationInspectEntry>> {
        self.per_orientation_inspect.as_ref().map(|entries| {
            entries
                .iter()
                .map(|entry| UVProfileOrientationInspectEntry {
                    orientation_label: entry.orientation_label.clone(),
                    analysis: entry.analysis.clone(),
                })
                .collect()
        })
    }
}

#[derive(Debug, Clone)]
struct OrientationInspectAnalysis {
    orientation_label: String,
    analysis: Option<InspectAnalysis>,
}

#[derive(Debug, Clone, Serialize)]
pub struct AnalyticAnalysis {
    per_orientation: BTreeMap<OrientationData, OrientationAnalyticAnalysis>,
}

#[derive(Debug, Clone, Serialize)]
pub struct OrientationAnalyticAnalysis {
    is_constant: bool,
    leading_coef: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct InspectAnalysis {
    result: FitResult,
    estimated_dod: i64,
    used_arb_prec_retry: bool,
}

pub struct AnalyticResult {
    pub(crate) per_orientations: BTreeMap<OrientationData, Series<AtomField>>,
}
