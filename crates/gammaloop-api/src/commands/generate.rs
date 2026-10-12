// file: src/cmd/generate.rs
#![allow(clippy::too_many_arguments)]

use color_eyre::owo_colors::OwoColorize;
use tracing::info;

use std::collections::BTreeMap;
use std::ffi::OsStr;
use std::fs;
use std::num::ParseIntError;
use std::path::Path;
use std::str::FromStr;

use clap::{Args, Parser, Subcommand, ValueEnum};
use color_eyre::Result;
#[cfg(test)]
use regex::Regex;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use symbolica::parse;
use tabled::{
    builder::Builder,
    settings::{style::HorizontalLine, themes::Theme, Style},
};
use thiserror::Error;
use walkdir::WalkDir;

use eyre::{eyre, Context};
use feynkit_generator::{
    Bridges, GenerationError, GenerationType, GraphGroupingOptions, Grouping, Loops, Naming,
    ParticleSelector, Range, Request, RequestParseError, Symmetrize, TopologyRules, VertexSelector,
    Veto,
};
use feynkit_graph::EdgeId;
use gammalooprs::model::Model;
use gammalooprs::processes::amplitude::Amplitude;
use gammalooprs::processes::{
    merge_generated_graph_reports, CrossSection, GeneratedGraphReport, GraphGenerationStats,
    Process, ProcessDefinition, ProcessList,
};
use gammalooprs::settings::{GlobalSettings, RuntimeSettings};

use crate::commands::set::KvPair;
use crate::completion::CompletionArgExt;
use crate::state::{GenerationResourceSummary, ProcessRef, State};

// =================== CLI containers (kept close to your structure) ===================

#[derive(Debug, Parser, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
/// Generate cross-section or amplitude integrands from a process specification.
pub struct Generate {
    /// Keep generated C++ source files after external compilation
    #[arg(long = "keep-sources", default_value_t = false, global = true)]
    #[serde(default)]
    pub keep_sources: bool,
    /// Selects cross-section, amplitude, or existing-process generation; when
    /// omitted, generates all integrands in the active state.
    #[command(subcommand)]
    pub mode: Option<GenerateCmd>,
}

#[derive(Debug, Subcommand, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
/// Select cross-section, amplitude, or existing-process integrand generation.
pub enum GenerateCmd {
    /// Generate a cross-section integrand through the forward-scattering construction.
    Xs(SpecArgs),

    /// Generate an amplitude integrand for the supplied process specification.
    #[command(alias = "amplitude")]
    Amp(SpecArgs),

    /// Create another integrand from an existing process without regenerating its graphs.
    Existing(ProcessArgs),
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq, ValueEnum)]
pub enum GroupingChoice {
    #[clap(name = "no_grouping")]
    NoGrouping,
    #[clap(name = "group_identical_graphs_up_to_sign")]
    GroupIdenticalGraphsUpToSign,
    #[clap(name = "group_identical_graphs_up_to_scalar_rescaling")]
    GroupIdenticalGraphsUpToScalarRescaling,
}

impl GroupingChoice {
    fn to_strategy(
        self,
        seed: Option<u16>,
        num_samples: Option<usize>,
        differentiate_particle_masses_only: Option<bool>,
        fully_numerical_substitution_when_comparing_numerators: Option<bool>,
        test_canonized_numerator: Option<bool>,
        symmetric_polarizations: Option<bool>,
    ) -> Grouping {
        let mut graph_grouping_options = GraphGroupingOptions::default();
        if let Some(seed) = seed {
            graph_grouping_options.numerical_sample_seed = seed;
        }
        if let Some(num_samples) = num_samples {
            graph_grouping_options.number_of_numerical_samples = num_samples;
        }
        if let Some(differentiate_particle_masses_only) = differentiate_particle_masses_only {
            graph_grouping_options.differentiate_particle_masses_only =
                differentiate_particle_masses_only;
        }
        if let Some(test_canonized_numerator) = test_canonized_numerator {
            graph_grouping_options.check_canonical_numerator = test_canonized_numerator;
        }
        if let Some(fully_numerical_substitution_when_comparing_numerators) =
            fully_numerical_substitution_when_comparing_numerators
        {
            graph_grouping_options.fully_numerical_substitution =
                fully_numerical_substitution_when_comparing_numerators;
        }
        if let Some(symmetric_polarizations) = symmetric_polarizations {
            graph_grouping_options.symmetric_polarizations = symmetric_polarizations;
        }
        match self {
            GroupingChoice::GroupIdenticalGraphsUpToScalarRescaling => {
                Grouping::UpToScalar(graph_grouping_options)
            }
            GroupingChoice::GroupIdenticalGraphsUpToSign => {
                Grouping::UpToSign(graph_grouping_options)
            }
            GroupingChoice::NoGrouping => Grouping::None,
        }
    }
}

#[derive(Args, Debug, Clone, Default, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct SpecArgs {
    /// Quote-free spec.
    ///
    /// Grammar (unquoted): `INITIAL to FINAL [process_options] [generation_options]`
    /// Also accepts: `INITIAL > FINAL [ ... ]`
    ///
    /// Examples:
    ///   e+ e- to d d~ g [ {1} ] veto u d QED^2==4 QCD^2<=4
    ///   e+ e- > { z z, a a } [ {1} {{2}} QCD=2 QED=1 ] / u d | g ghG QED==2 QCD>=2 QCD<=4
    ///
    /// Notes:
    /// - `veto` and `only` are synonyms for `/` and `|`.
    /// - The `[ ... ]` block holds perturbative spec: `{n}`, `{{n}}`, and orders `QCD=n` with shorthand `QCD`≡1.
    #[arg(value_name = "TOKENS", num_args = 1..)]
    pub tokens: Vec<String>,

    /// Append the generated process instead of replacing the current process selection.
    #[arg(short = 'a', default_value_t = false)]
    pub append: bool,

    /// Clear pre-existing processes
    #[arg(
        long = "clear-existing-processes",
        short = 'c',
        alias = "clear",
        default_value_t = false
    )]
    pub clear_existing_processes: bool,

    /// Optional human-readable process name
    #[arg(long = "process-name", short = 'p')]
    pub process_name: Option<String>,

    /// Name assigned to the generated integrand; defaults to the generated process name.
    #[arg(
        long = "integrand-name",
        short = 'i',
        completion_disable_special_value()
    )]
    pub integrand_name: Option<String>,

    /// Stop after diagram generation without constructing an evaluable integrand.
    #[arg(long = "only-diagrams", short = 'o', default_value_t = false)]
    pub only_diagrams: bool,

    /// Veto self-energy topologies; omission keeps them.
    #[arg(long = "filter-selfenergies")]
    pub filter_selfenergies: Option<bool>,
    /// Veto snail topologies; omission keeps them.
    #[arg(long = "filter-snails")]
    pub filter_snails: Option<bool>,
    /// Veto tadpole topologies; omission keeps them.
    #[arg(long = "filter-tadpoles")]
    pub filter_tadpoles: Option<bool>,

    /// Exclude graphs containing any listed vertex-interaction name
    #[arg(long = "veto-vertex-interactions", num_args = 0..)]
    pub veto_vertex_interactions: Option<Vec<String>>,

    /// Restrict generation to an allow-list of vertex interactions
    #[arg(long = "allowed-vertex-interactions", num_args = 0..)]
    pub allowed_vertex_interactions: Option<Vec<String>>,

    /// Filter tadpoles formed while sewing amplitudes into a cross-section graph.
    #[arg(long = "filter-cross-section-tadpoles")]
    pub filter_cross_section_tadpoles: Option<bool>,

    /// Veto tadpoles attached to massive propagators.
    #[arg(long = "veto-tadpoles-attached-to-massive-lines")]
    pub veto_tadpoles_attached_to_massive_lines: Option<bool>,
    /// Veto tadpoles attached to massless propagators.
    #[arg(long = "veto-tadpoles-attached-to-massless-lines")]
    pub veto_tadpoles_attached_to_massless_lines: Option<bool>,

    /// Veto snail insertions attached to massive propagators.
    #[arg(long = "veto-snails-attached-to-massive-lines")]
    pub veto_snails_attached_to_massive_lines: Option<bool>,
    /// Veto snail insertions attached to massless propagators.
    #[arg(long = "veto-snails-attached-to-massless-lines")]
    pub veto_snails_attached_to_massless_lines: Option<bool>,

    /// Veto self-energy insertions on massive propagators.
    #[arg(long = "veto-self-energy-of-massive-lines")]
    pub veto_self_energy_of_massive_lines: Option<bool>,
    /// Veto self-energy insertions on massless propagators.
    #[arg(long = "veto-self-energy-of-massless-lines")]
    pub veto_self_energy_of_massless_lines: Option<bool>,

    /// Maximum number of bridges allowed in a generated graph; omission is unrestricted.
    #[arg(long = "max-n-bridges", short = 'b')]
    pub max_n_bridges: Option<usize>,
    /// Inclusive minimum and maximum number of factorized loop subtopologies.
    #[arg(
        long = "number-of-factorized-loop-subtopologies",
        short = 'f',
        alias = "nfactl",
        num_args = 2
    )]
    pub number_of_factorized_loop_subtopologies: Option<Vec<usize>>,
    /// Inclusive minimum and maximum number of closed fermion loops, excluding ghosts.
    #[arg(long = "number-of-fermion-loops", short = 'L', num_args = 2)]
    pub number_of_fermion_loops: Option<Vec<usize>>,

    /// Inclusive minimum and maximum number of cut blobs on either side of a cross-section cut.
    #[arg(long = "n-cut-blobs", short = 'B', num_args = 2)]
    pub n_cut_blobs: Option<Vec<usize>>,
    /// Inclusive minimum and maximum number of spectators crossing a cross-section cut.
    #[arg(long = "n-cut-spectators", short = 'S', num_args = 2)]
    pub n_cut_spectators: Option<Vec<usize>>,

    /// Permit external fermion permutations when identifying equivalent amplitudes.
    #[arg(
        long = "allow-symmetrization-of-external-fermions-in-amplitudes",
        alias = "symferm"
    )]
    pub allow_symmetrization_of_external_fermions_in_amplitudes: Option<bool>,
    /// Identify graphs related by permutations of identical initial-state particles. A cross
    /// section keeps one graph per exchange of its two identical beams, weighted by the number
    /// of graphs it stands for; this assumes observables, selectors and beam helicities
    /// symmetric under the exchange, and can raise the integration error at equal sample counts.
    /// With left/right symmetrization the weights inherit its deviations and can differ from
    /// its total too.
    #[arg(long = "symmetrize-initial-states")]
    pub symmetrize_initial_states: Option<bool>,
    /// Identify amplitude graphs related by permutations of identical final-state particles;
    /// cross sections reject it, their identical final states being cut edges.
    #[arg(long = "symmetrize-final-states")]
    pub symmetrize_final_states: Option<bool>,
    /// Identify cross-section graphs under CP-based left/right exchange; requires CP symmetry.
    /// From two loops on the physical cuts can depend on which side is searched, so the total
    /// can differ from the unsymmetrized one.
    #[arg(long = "symmetrize-left-right-states")]
    pub symmetrize_left_right_states: Option<bool>,

    /// Numerator-aware grouping choice
    #[arg(long = "numerator-grouping", short = 'G', value_enum)]
    pub numerator_aware_isomorphism_grouping: Option<GroupingChoice>,

    /// Remove diagrams whose numerators vanish under color algebra before grouping.
    #[arg(long = "filter-zero-color", default_value_t = false)]
    #[serde(default)]
    pub filter_zero_color: bool,

    /// Deterministic random seed used for numerical numerator comparisons.
    #[arg(long = "numerical-samples-seed")]
    pub numerical_samples_seed: Option<u16>,
    /// Number of phase-space samples used to decide whether two numerators are equivalent.
    #[arg(long = "number-of-samples-for-numerator-comparisons")]
    pub number_of_samples_for_numerator_comparisons: Option<usize>,
    /// Distinguish mass parameters while ignoring other model-parameter labels during grouping.
    #[arg(long = "consider-internal-masses-only-in-numerator-isomorphisms")]
    pub consider_internal_masses_only_in_numerator_isomorphisms: Option<bool>,
    /// Numerically substitute every parameter when comparing candidate numerator isomorphisms.
    #[arg(long = "fully-numerical-substitution-when-comparing-numerators")]
    pub fully_numerical_substitution_when_comparing_numerators: Option<bool>,
    /// Compare canonicalized numerators in addition to graph topology when grouping.
    #[arg(long = "compare-canonized-numerator")]
    pub compare_canonized_numerator: Option<bool>,
    /// Treat left- and right-side polarization factors symmetrically during numerator grouping.
    #[arg(long = "symmetric-left-right-polarizations")]
    pub symmetric_left_right_polarizations: Option<bool>,

    /// Graph processing toggles
    ///
    /// Format:
    ///   --loop-momentum-bases "GL_12=7,10 GL_77=4,2"
    ///   --select-graphs "GL_12 GL_13"
    ///   --veto-graphs "GL_11 GL_15"
    #[arg(long = "loop-momentum-bases", value_name = "KEY=VALUE", num_args = 0.., value_parser = KvPair::from_str)]
    pub loop_momentum_bases: Option<Vec<KvPair>>,
    /// Generate only graphs whose names appear in this allow-list.
    #[arg(long = "select-graphs", num_args = 0..)]
    pub select_graphs: Option<Vec<String>>,
    /// Exclude graphs whose names appear in this deny-list.
    #[arg(long = "veto-graphs", num_args = 0..)]
    pub veto_graphs: Option<Vec<String>>,
    /// Prefix assigned to generated graph names; defaults to `GL`
    #[arg(long = "graph-prefix", short = 'g')]
    pub graph_prefix: Option<String>,

    /// Global prefactor projector (Symbolica atom string)
    #[arg(
        long = "global-prefactor-projector",
        value_name = "ATOM",
        allow_hyphen_values = true
    )]
    pub global_prefactor_projector: Option<String>,

    /// Global prefactor numerator (Symbolica atom string)
    #[arg(
        long = "global-prefactor-num",
        value_name = "ATOM",
        allow_hyphen_values = true
    )]
    pub global_prefactor_num: Option<String>,

    /// Filter graph-theoretic self-loops explicitly; omission keeps them.
    #[arg(long = "filter-self-loop")]
    pub filter_self_loop: Option<bool>,

    /// Filter edges that carry zero momentum flow in the selected routing.
    #[arg(long = "filter-zero-flow-edges")]
    pub filter_zero_flow_edges: Option<bool>,
}

impl SpecArgs {
    pub(crate) fn from_process_spec_string(spec: impl Into<String>) -> Self {
        Self {
            tokens: vec![spec.into()],
            ..Self::default()
        }
    }
}

// =================== Runner ===================

fn format_generation_duration(duration: std::time::Duration) -> String {
    if duration.as_secs() >= 60 {
        let minutes = duration.as_secs() / 60;
        let seconds = duration.as_secs_f64() - (minutes * 60) as f64;
        format!("{minutes}m {seconds:.1}s")
    } else if duration.as_secs_f64() >= 1.0 {
        format!("{:.2}s", duration.as_secs_f64())
    } else {
        format!("{}ms", duration.as_millis())
    }
}

fn format_generation_fraction(
    numerator: std::time::Duration,
    total: std::time::Duration,
) -> String {
    let percent = if total.is_zero() {
        0.0
    } else {
        (numerator.as_secs_f64() / total.as_secs_f64()) * 100.0
    };

    if percent >= 10.0 {
        format!("{percent:.0}%")
    } else if percent >= 1.0 {
        format!("{percent:.1}%")
    } else {
        format!("{percent:.2}%")
    }
}

fn format_generation_memory(bytes: u64) -> String {
    const KIB: f64 = 1024.0;
    const MIB: f64 = KIB * 1024.0;
    const GIB: f64 = MIB * 1024.0;
    const TIB: f64 = GIB * 1024.0;

    let value = bytes as f64;
    if value < KIB {
        format!("{bytes} B")
    } else if value < MIB {
        format!("{:.2} KiB", value / KIB)
    } else if value < GIB {
        format!("{:.2} MiB", value / MIB)
    } else if value < TIB {
        format!("{:.2} GiB", value / GIB)
    } else {
        format!("{:.2} TiB", value / TIB)
    }
}

pub(crate) fn render_generation_summary(
    reports: &[GeneratedGraphReport],
    peak_ram_bytes: u64,
    generation_cores: Option<usize>,
    title: Option<&str>,
) -> Option<String> {
    if reports.is_empty() {
        return None;
    }

    let mut sorted_reports = reports.to_vec();
    sorted_reports.sort_by(|left, right| {
        (
            left.process_id,
            left.integrand_name.as_str(),
            left.graph_name.as_str(),
        )
            .cmp(&(
                right.process_id,
                right.integrand_name.as_str(),
                right.graph_name.as_str(),
            ))
    });

    let mut builder = Builder::new();
    builder.push_record([
        "integrand".bold().blue().to_string(),
        "graph".bold().blue().to_string(),
        "# evals".bold().blue().to_string(),
        "expr build".bold().blue().to_string(),
        "spenso".bold().blue().to_string(),
        "symbolica eval".bold().blue().to_string(),
        "compile".bold().blue().to_string(),
    ]);

    let mut total_stats = GraphGenerationStats::default();
    for report in &sorted_reports {
        total_stats.merge_in_place(&report.stats);
        let expr_time = report.stats.expression_build_time();
        let total_time = report.stats.total_time;
        let expr_value = format!(
            "{} ({})",
            format_generation_duration(expr_time).magenta(),
            format_generation_fraction(expr_time, total_time).cyan()
        );
        let spenso_value = format!(
            "{} ({})",
            format_generation_duration(report.stats.evaluator_spenso_time).magenta(),
            format_generation_fraction(report.stats.evaluator_spenso_time, total_time).cyan()
        );
        let symbolica_value = format!(
            "{} ({})",
            format_generation_duration(report.stats.evaluator_symbolica_time).magenta(),
            format_generation_fraction(report.stats.evaluator_symbolica_time, total_time).cyan()
        );
        let compile_value = format!(
            "{} ({})",
            format_generation_duration(report.stats.evaluator_compile_time).magenta(),
            format_generation_fraction(report.stats.evaluator_compile_time, total_time).cyan()
        );

        builder.push_record([
            report.integrand_name.yellow().to_string(),
            report.graph_name.yellow().to_string(),
            report
                .stats
                .evaluator_count
                .to_string()
                .yellow()
                .to_string(),
            expr_value,
            spenso_value,
            symbolica_value,
            compile_value,
        ]);
    }

    let total_time = total_stats.total_time;
    let total_expr_time = total_stats.expression_build_time();
    let total_expr_value = format!(
        "{} ({})",
        format_generation_duration(total_expr_time).magenta(),
        format_generation_fraction(total_expr_time, total_time).cyan()
    );
    let total_spenso_value = format!(
        "{} ({})",
        format_generation_duration(total_stats.evaluator_spenso_time).magenta(),
        format_generation_fraction(total_stats.evaluator_spenso_time, total_time).cyan()
    );
    let total_symbolica_value = format!(
        "{} ({})",
        format_generation_duration(total_stats.evaluator_symbolica_time).magenta(),
        format_generation_fraction(total_stats.evaluator_symbolica_time, total_time).cyan()
    );
    let total_compile_value = format!(
        "{} ({})",
        format_generation_duration(total_stats.evaluator_compile_time).magenta(),
        format_generation_fraction(total_stats.evaluator_compile_time, total_time).cyan()
    );
    builder.push_record([
        "Total".bold().yellow().to_string(),
        String::new(),
        total_stats
            .evaluator_count
            .to_string()
            .bold()
            .yellow()
            .to_string(),
        total_expr_value,
        total_spenso_value,
        total_symbolica_value,
        total_compile_value,
    ]);

    let mut table = builder.build();
    let mut style = Theme::from_style(Style::rounded());
    style.insert_horizontal_line(
        sorted_reports.len() + 1,
        HorizontalLine::inherit(Style::modern()),
    );
    table.with(style);
    let mut sections = Vec::new();
    if let Some(title) = title {
        sections.push(title.bold().blue().to_string());
    }
    let resources = match generation_cores {
        Some(cores) => {
            format!(
                "{} {} | {} {}",
                "peak RAM".bold().blue(),
                format_generation_memory(peak_ram_bytes).yellow(),
                "cores".bold().blue(),
                cores.to_string().yellow(),
            )
        }
        None => format!(
            "{} {}",
            "peak RAM".bold().blue(),
            format_generation_memory(peak_ram_bytes).yellow(),
        ),
    };
    sections.push(resources);
    sections.push(table.to_string());
    Some(sections.join("\n"))
}

fn log_generation_summary_table(
    reports: &[GeneratedGraphReport],
    resources: GenerationResourceSummary,
) {
    if let Some(summary) = render_generation_summary(
        reports,
        resources.peak_ram_bytes,
        Some(resources.generation_cores),
        Some("Integrand generation summary"),
    ) {
        info!("\n{summary}");
    }
}

fn remove_compiled_cpp_sources(root: &Path) -> Result<usize> {
    let mut removed_count = 0usize;
    for entry in WalkDir::new(root) {
        let entry = entry?;
        if !entry.file_type().is_file() {
            continue;
        }

        let source_path = entry.path();
        if source_path.extension() != Some(OsStr::new("cpp")) {
            continue;
        }

        let library_path = source_path.with_extension("so");
        if !library_path.is_file() {
            continue;
        }

        fs::remove_file(source_path).with_context(|| {
            format!(
                "Trying to remove generated C++ source {} after compiling {}",
                source_path.display(),
                library_path.display()
            )
        })?;
        removed_count += 1;
    }

    Ok(removed_count)
}

fn compile_integrands_for_generation(
    state: &mut State,
    compile_folder: &Path,
    override_existing_compiled: bool,
    global_settings: &GlobalSettings,
    process_id: Option<usize>,
    integrand_name: Option<String>,
    keep_sources: bool,
) -> Result<Vec<GeneratedGraphReport>> {
    let reports = state.compile_integrands(
        compile_folder,
        override_existing_compiled,
        global_settings,
        process_id,
        integrand_name,
    )?;

    if !keep_sources {
        let removed_count = remove_compiled_cpp_sources(compile_folder)?;
        if removed_count > 0 {
            info!(
                "Removed {} generated C++ source file(s). Use '{}' to keep them.",
                removed_count,
                "--keep-sources".green()
            );
        }
    }

    Ok(reports)
}

fn finish_generation(
    state: &mut State,
    reports: &[GeneratedGraphReport],
    resources: GenerationResourceSummary,
) {
    state.record_generation_summary(reports, resources);
    log_generation_summary_table(reports, resources);
}

impl Generate {
    pub fn run(
        &self,
        state: &mut State,
        compile_folder: impl AsRef<Path>,
        override_existing_compiled: bool,
        global_settings: &GlobalSettings,
        runtime_settings: &RuntimeSettings,
    ) -> Result<()> {
        let compile_folder = compile_folder.as_ref();
        let generation_mode = match &self.mode {
            Some(GenerateCmd::Xs(a)) => Some((GenerationType::CrossSection, a)),
            Some(GenerateCmd::Amp(a)) => Some((GenerationType::Amplitude, a)),
            _ => None,
        };
        if let Some((_, args)) = generation_mode.as_ref() {
            if !state.process_list.processes.is_empty() && args.clear_existing_processes {
                info!(
                    "Clearing all {} existing processes as requested.",
                    state.process_list.processes.len()
                );
                state.process_list = ProcessList::default();
            }
        }
        let generation_info = if let Some((gen_mode, args)) = generation_mode {
            let mut spec = parse_spec_with_model(args, gen_mode, &state.model)?;
            spec.process_definition.process_id = state.process_list.processes.len();

            let mut existing_process = None;
            if let Some(ep) = state
                .process_list
                .processes
                .iter_mut()
                .find(|p| p.definition.folder_name == spec.process_definition.folder_name)
            {
                if ep.definition != spec.process_definition {
                    if !args.append {
                        return Err(eyre!(
                            "Process with name '{}' already exists.\n> Use 'existing' subcommand to continue generation of this process.\n> Use '--clear-existing-processes' to remove all existing ones.\n> Or specify a different process name with '--process-name <chosen_process_name>'.",
                            spec.process_definition.folder_name
                        ));
                    }
                } else {
                    info!(
                        "Identical process definition, with name '{}', already exists. Gammaloop will recycle it.",
                        spec.process_definition.folder_name
                    );
                    return Ok(());
                }
                spec.process_definition.process_id = ep.definition.process_id;
                existing_process = Some(ep);
            }
            Some((spec, existing_process))
        } else {
            None
        };
        match &self.mode {
            Some(GenerateCmd::Amp(args)) | Some(GenerateCmd::Xs(args)) => {
                let generation_type = generation_mode.as_ref().unwrap().0;
                let model: &Model = &state.model;
                let (spec, existing_process) = generation_info.unwrap();
                let this_process_id = spec.process_definition.process_id;
                // TODO handle existing process and continue
                let graphs = spec.process_definition.generate(model, global_settings)?;
                info!(
                    "Generated {} {} graphs.",
                    if matches!(self.mode, Some(GenerateCmd::Amp(_))) {
                        "amplitude"
                    } else {
                        "cross-section"
                    },
                    graphs.len()
                );
                // Keep the possibility of changing default name for the two modes
                let integrand_base_name = matches!(self.mode, Some(GenerateCmd::Amp(_)))
                    .then(|| args.integrand_name.clone().unwrap_or("default".to_string()))
                    .unwrap_or_else(|| {
                        args.integrand_name.clone().unwrap_or("default".to_string())
                    });
                let generated_integrand_name = if let Some(p) = existing_process {
                    let existing_names = p.collection.get_integrand_names();
                    let integrand_name = if existing_names.contains(&integrand_base_name.as_str()) {
                        let mut integrand_i = 0;
                        while existing_names
                            .iter()
                            .any(|ce| *ce == format!("{}_{}", integrand_base_name, integrand_i))
                        {
                            integrand_i += 1;
                        }
                        format!("{}_{}", integrand_base_name, integrand_i)
                    } else {
                        integrand_base_name
                    };
                    match &self.mode {
                        Some(GenerateCmd::Amp(_)) => {
                            p.collection.add_amplitude(Amplitude::from_graph_list(
                                integrand_name.clone(),
                                graphs,
                            )?);
                        }
                        Some(GenerateCmd::Xs(_)) => {
                            p.collection
                                .add_cross_section(CrossSection::from_graph_list(
                                    integrand_name.clone(),
                                    graphs,
                                )?);
                        }
                        _ => unreachable!(),
                    }
                    integrand_name
                } else {
                    let process = Process::from_graph_list(
                        spec.process_definition.folder_name.clone(),
                        integrand_base_name.clone(),
                        graphs,
                        generation_type,
                        Some(spec.process_definition),
                        None,
                        model,
                    )?;
                    state.process_list.add_process(process);
                    integrand_base_name
                };
                if !args.only_diagrams {
                    let generated_integrand_name_for_compile = generated_integrand_name.clone();
                    let generation = state.generate_integrand(
                        global_settings,
                        runtime_settings.into(),
                        this_process_id,
                        Some(generated_integrand_name),
                    )?;
                    let mut reports = generation.reports;
                    if global_settings.generation.evaluator.compile
                        && global_settings
                            .generation
                            .compile
                            .requires_external_compilation()
                    {
                        merge_generated_graph_reports(
                            &mut reports,
                            compile_integrands_for_generation(
                                state,
                                compile_folder,
                                override_existing_compiled,
                                global_settings,
                                Some(this_process_id),
                                Some(generated_integrand_name_for_compile),
                                self.keep_sources,
                            )?,
                        );
                    }
                    finish_generation(state, &reports, generation.resources);
                    Ok(())
                } else {
                    info!(
                        "Only diagram generation was requested, skipping integrand generation. You can generate integrands later using the '{}' command.",
                        "generate existing <options>".green()
                    );
                    Ok(())
                }
            }
            Some(GenerateCmd::Existing(process_args)) => {
                let process_id = state.resolve_process_ref(process_args.process.as_ref())?;

                let generation = state.generate_integrand(
                    global_settings,
                    runtime_settings.into(),
                    process_id,
                    process_args.integrand_name.clone(),
                )?;
                let mut reports = generation.reports;
                if global_settings.generation.evaluator.compile
                    && global_settings
                        .generation
                        .compile
                        .requires_external_compilation()
                {
                    merge_generated_graph_reports(
                        &mut reports,
                        compile_integrands_for_generation(
                            state,
                            compile_folder,
                            override_existing_compiled,
                            global_settings,
                            Some(process_id),
                            process_args.integrand_name.clone(),
                            self.keep_sources,
                        )?,
                    );
                }
                finish_generation(state, &reports, generation.resources);
                Ok(())
            }
            None => {
                let generation =
                    state.generate_integrands(global_settings, runtime_settings.into())?;
                let mut reports = generation.reports;
                if global_settings.generation.evaluator.compile
                    && global_settings
                        .generation
                        .compile
                        .requires_external_compilation()
                {
                    merge_generated_graph_reports(
                        &mut reports,
                        compile_integrands_for_generation(
                            state,
                            compile_folder,
                            override_existing_compiled,
                            global_settings,
                            None,
                            None,
                            self.keep_sources,
                        )?,
                    );
                }
                finish_generation(state, &reports, generation.resources);
                Ok(())
            }
        }
    }
}

#[derive(Args, Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct ProcessArgs {
    /// Process reference: `#<id>`, `name:<name>`, or `<id>/<name>`
    #[arg(
        long = "process",
        short = 'p',
        value_name = "PROCESS",
        completion_process_selector(crate::completion::SelectorKind::Any)
    )]
    pub process: Option<ProcessRef>,

    /// Optional human-readable integrand name used to disambiguate the process
    #[arg(
        long = "integrand-name",
        short = 'i',
        completion_integrand_selector(crate::completion::SelectorKind::Any)
    )]
    pub integrand_name: Option<String>,
}

/// A parsed CLI specification with phase settings owned by the saved definition.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct ProcessSpec {
    #[schemars(with = "serde_json::Value")]
    pub process_definition: ProcessDefinition,
}

impl ProcessSpec {
    fn particle_name(model: &Model, selector: &ParticleSelector) -> String {
        selector
            .resolve(model)
            .and_then(|id| model.particle_by_id(id).map_err(Into::into))
            .map(|particle| particle.name.to_ascii_lowercase())
            .unwrap_or_else(|_| selector.to_string())
    }

    fn order_terms(orders: &BTreeMap<String, Range>, squared: bool, shell: bool) -> Vec<String> {
        let mut terms = Vec::new();
        for (name, range) in orders {
            let name = if squared {
                format!("{name}{}", if shell { "sq" } else { "^2" })
            } else {
                name.clone()
            };
            let (min, max) = range.bounds();
            if max == Some(min) {
                terms.push(format!("{name}{}{min}", if shell { "_eq_" } else { "==" }));
            } else {
                if min > 0 {
                    terms.push(format!("{name}{}{min}", if shell { "_ge_" } else { ">=" }));
                }
                if let Some(max) = max {
                    terms.push(format!("{name}{}{max}", if shell { "_le_" } else { "<=" }));
                }
            }
        }
        terms
    }

    fn plain_orders(&self) -> &BTreeMap<String, Range> {
        let request = &self.process_definition.request;
        match request.kind() {
            GenerationType::Amplitude => &request.orders().plain,
            GenerationType::CrossSection => &request.cuts().amplitude.orders.plain,
        }
    }

    /// Deterministic shell name derived from the resolved request.
    pub fn process_shell_name(&self, model: &Model, short: bool) -> String {
        let request = &self.process_definition.request;
        let sanitize = |s: &str| s.replace('~', "x").replace('+', "p").replace('-', "m");
        let states = |particles: &[ParticleSelector]| {
            if particles.is_empty() {
                "empty".to_owned()
            } else {
                particles
                    .iter()
                    .map(|p| Self::particle_name(model, p))
                    .collect::<String>()
            }
        };
        let finals = request
            .process()
            .outgoing_alternatives()
            .iter()
            .map(|s| states(s))
            .collect::<Vec<_>>()
            .join("_or_");
        let mut pieces = vec![sanitize(&format!(
            "{}_{finals}",
            states(request.process().incoming())
        ))];
        if !short {
            let mut veto = request
                .restrictions()
                .particle_exclude
                .iter()
                .map(|p| Self::particle_name(model, p))
                .collect::<Vec<_>>();
            veto.sort();
            if !veto.is_empty() {
                pieces.push(sanitize(&format!("no__{}", veto.join("_"))));
            }
            for terms in [
                Self::order_terms(self.plain_orders(), false, true),
                Self::order_terms(&request.orders().squared, true, true),
            ] {
                if !terms.is_empty() {
                    pieces.push(terms.join("__"));
                }
            }
            if !request.orders().perturbative.is_empty() {
                pieces.push(
                    request
                        .orders()
                        .perturbative
                        .iter()
                        .map(|(name, n)| format!("{name}loop_eq_{n}"))
                        .collect::<Vec<_>>()
                        .join("__"),
                );
            }
        }
        pieces.join("-")
    }

    /// Canonical text for the shared process grammar.
    pub fn repr_str(&self, model: &Model) -> String {
        let request = &self.process_definition.request;
        let states = |particles: &[ParticleSelector]| {
            if particles.is_empty() {
                "{}".to_owned()
            } else {
                particles
                    .iter()
                    .map(|p| Self::particle_name(model, p))
                    .collect::<Vec<_>>()
                    .join(" ")
            }
        };
        let finals = request.process().outgoing_alternatives();
        let final_text = if finals.len() == 1 {
            states(&finals[0])
        } else {
            format!(
                "{{ {} }}",
                finals
                    .iter()
                    .map(|s| states(s))
                    .collect::<Vec<_>>()
                    .join(", ")
            )
        };
        let mut out = vec![
            states(request.process().incoming()),
            ">".to_owned(),
            final_text,
        ];
        let mut veto = request
            .restrictions()
            .particle_exclude
            .iter()
            .map(|p| Self::particle_name(model, p))
            .collect::<Vec<_>>();
        veto.sort();
        if !veto.is_empty() {
            out.push("/".to_owned());
            out.extend(veto);
        }
        for clause in &request.restrictions().particle_allow {
            let mut names = clause
                .iter()
                .map(|p| Self::particle_name(model, p))
                .collect::<Vec<_>>();
            names.sort();
            out.push("|".to_owned());
            out.extend(names);
        }
        out.extend(Self::order_terms(self.plain_orders(), false, false));
        out.push("[".to_owned());
        if request.kind() == GenerationType::CrossSection {
            let (min, max) = request.cuts().amplitude.loops.bounds();
            if max == Some(min) {
                out.push(format!("{{{min}}}"));
            }
        }
        match request.loops() {
            Loops::Total(range) => out.push(format!("{{{}}}", range.start())),
            Loops::Forward(range) => out.push(format!("{{{{{}}}}}", range.start())),
        }
        for (name, n) in &request.orders().perturbative {
            out.push(if *n == 1 {
                name.clone()
            } else {
                format!("{name}={n}")
            });
        }
        out.push("]".to_owned());
        out.extend(Self::order_terms(&request.orders().squared, true, false));
        out.join(" ")
    }
}

#[derive(Debug, Error)]
pub enum ParseError {
    #[error(transparent)]
    Request(#[from] RequestParseError),
    #[error(transparent)]
    Generation(#[from] GenerationError),
    #[error("invalid lmb specification: {0}")]
    InvalidLmbSpec(ParseIntError),
    #[error("{0} requires an inclusive minimum and maximum")]
    InvalidRange(&'static str),
}

impl SpecArgs {
    fn range(values: &Option<Vec<usize>>, name: &'static str) -> Result<Range, ParseError> {
        match values.as_deref() {
            None => Ok(Range::Unrestricted),
            Some([min, max]) if min <= max => Ok((*min..=*max).into()),
            Some(_) => Err(ParseError::InvalidRange(name)),
        }
    }

    fn veto(
        enabled: Option<bool>,
        massive: Option<bool>,
        massless: Option<bool>,
        defaults: Veto,
    ) -> Veto {
        if enabled == Some(false) {
            return Veto::default();
        }
        let defaults = if enabled == Some(true) {
            defaults
        } else {
            Veto::default()
        };
        Veto {
            massive: massive.unwrap_or(defaults.massive),
            massless: massless.unwrap_or(defaults.massless),
        }
    }
}

pub fn parse_spec_with_model(
    args: &SpecArgs,
    kind: GenerationType,
    model: &Model,
) -> Result<ProcessSpec, ParseError> {
    let mut request = Request::parse(model, kind, &args.tokens.join(" "), None)?;
    let mut restrictions = feynkit_generator::Restrict::new();
    if let Some(allowed) = &args.allowed_vertex_interactions {
        restrictions =
            restrictions.only_vertices(allowed.iter().cloned().map(VertexSelector::Name));
    }
    if let Some(vetoed) = &args.veto_vertex_interactions {
        restrictions =
            restrictions.exclude_vertices(vetoed.iter().cloned().map(VertexSelector::Name));
    }
    let topology = TopologyRules {
        bridges: args.max_n_bridges.map_or(Bridges::Any, Bridges::AtMost),
        self_loops: !args.filter_self_loop.unwrap_or(false),
        zero_flow_edges: !args.filter_zero_flow_edges.unwrap_or(true),
        self_energies: SpecArgs::veto(
            args.filter_selfenergies,
            args.veto_self_energy_of_massive_lines,
            args.veto_self_energy_of_massless_lines,
            Veto::all(),
        ),
        tadpoles: SpecArgs::veto(
            args.filter_tadpoles,
            args.veto_tadpoles_attached_to_massive_lines,
            args.veto_tadpoles_attached_to_massless_lines,
            Veto::all(),
        ),
        snails: SpecArgs::veto(
            args.filter_snails,
            args.veto_snails_attached_to_massive_lines,
            args.veto_snails_attached_to_massless_lines,
            Veto::massless(),
        ),
        fermion_loops: SpecArgs::range(&args.number_of_fermion_loops, "--number-of-fermion-loops")?,
        factorized_loops: SpecArgs::range(
            &args.number_of_factorized_loop_subtopologies,
            "--number-of-factorized-loop-subtopologies",
        )?,
        ..TopologyRules::default()
    };
    let amplitude_lr =
        kind == GenerationType::Amplitude && args.symmetrize_left_right_states.unwrap_or(false);
    let symmetrize = Symmetrize {
        initial: args.symmetrize_initial_states.unwrap_or(amplitude_lr),
        final_state: args.symmetrize_final_states.unwrap_or(amplitude_lr),
        external_fermions: args
            .allow_symmetrization_of_external_fermions_in_amplitudes
            .unwrap_or(false),
    };
    let mut cuts = request.cuts().clone();
    cuts.blobs = SpecArgs::range(&args.n_cut_blobs, "--n-cut-blobs")?;
    cuts.spectators = SpecArgs::range(&args.n_cut_spectators, "--n-cut-spectators")?;
    cuts.sewn_tadpoles = args.filter_cross_section_tadpoles.unwrap_or(false);
    cuts.symmetrize_left_right =
        kind == GenerationType::CrossSection && args.symmetrize_left_right_states.unwrap_or(false);
    request = request
        .with_restrict(restrictions)
        .with_topology(topology)
        .with_cuts(cuts)
        .with_symmetrize(symmetrize);
    request.validate_in(model)?;
    let mut definition = ProcessDefinition::new(request);
    definition.grouping = build_grouping_option(args);
    definition.numerators.naming =
        Naming::prefix(args.graph_prefix.clone().unwrap_or_else(|| "GL".to_owned()));
    definition.numerators.filter_zero_color = args.filter_zero_color;
    if let Some(selected) = &args.select_graphs {
        definition.selection = definition.selection.select_names(selected.clone());
    }
    if let Some(vetoed) = &args.veto_graphs {
        definition.selection = definition.selection.veto_names(vetoed.clone());
    }
    if let Some(bases) = &args.loop_momentum_bases {
        for basis in bases {
            let edges = basis
                .value
                .split(',')
                .flat_map(str::split_whitespace)
                .map(|edge| edge.parse::<usize>().map(EdgeId))
                .collect::<Result<Vec<_>, _>>()
                .map_err(ParseError::InvalidLmbSpec)?;
            definition.numerators = definition
                .numerators
                .with_named_loop_momentum_basis(&basis.key, edges);
        }
    }
    if let Some(projector) = &args.global_prefactor_projector {
        definition.numerators.projector = Some(parse!(projector));
    }
    if let Some(prefactor) = &args.global_prefactor_num {
        definition.numerators.prefactor = parse!(prefactor);
    }
    let mut spec = ProcessSpec {
        process_definition: definition,
    };
    spec.process_definition.folder_name = args
        .process_name
        .clone()
        .unwrap_or_else(|| spec.process_shell_name(model, true));
    Ok(spec)
}

pub(crate) fn parse_process_spec_string(
    spec: &str,
    kind: GenerationType,
    model: &Model,
) -> Result<ProcessSpec, ParseError> {
    parse_spec_with_model(&SpecArgs::from_process_spec_string(spec), kind, model)
}

fn build_grouping_option(args: &SpecArgs) -> Grouping {
    args.numerator_aware_isomorphism_grouping
        .unwrap_or(GroupingChoice::NoGrouping)
        .to_strategy(
            args.numerical_samples_seed,
            args.number_of_samples_for_numerator_comparisons,
            args.consider_internal_masses_only_in_numerator_isomorphisms,
            args.fully_numerical_substitution_when_comparing_numerators,
            args.compare_canonized_numerator,
            args.symmetric_left_right_polarizations,
        )
}

// =================== Tests ===================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{commands::Commands, Repl};
    use clap::Parser;
    use gammalooprs::initialisation::test_initialise;
    use gammalooprs::utils::load_generic_model;
    use std::time::Duration;

    #[test]
    fn remove_compiled_cpp_sources_only_removes_sources_with_matching_libraries() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let compiled_source = temp.path().join("integrand.cpp");
        let compiled_library = temp.path().join("integrand.so");
        let uncompiled_source = temp.path().join("orphan.cpp");
        std::fs::write(&compiled_source, "compiled source")?;
        std::fs::write(&compiled_library, "compiled library")?;
        std::fs::write(&uncompiled_source, "orphan source")?;

        let removed_count = remove_compiled_cpp_sources(temp.path())?;

        assert_eq!(removed_count, 1);
        assert!(!compiled_source.exists());
        assert!(compiled_library.is_file());
        assert!(uncompiled_source.is_file());
        Ok(())
    }

    #[test]
    fn parse_generate_keep_sources_after_subcommand() {
        let repl =
            Repl::try_parse_from(["gammaloop", "generate", "existing", "--keep-sources"]).unwrap();

        match repl.command {
            Commands::Generate(generate) => {
                assert!(generate.keep_sources);
                assert!(matches!(generate.mode, Some(GenerateCmd::Existing(_))));
            }
            other => panic!("Expected generate command, got {other:?}"),
        }
    }

    #[test]
    fn generation_summary_includes_separator_and_total_row() {
        let reports = vec![
            GeneratedGraphReport {
                process_id: 0,
                integrand_name: "itg".to_string(),
                graph_name: "GL01".to_string(),
                stats: GraphGenerationStats {
                    evaluator_count: 2,
                    total_time: Duration::from_secs(4),
                    evaluator_spenso_time: Duration::from_secs(1),
                    evaluator_symbolica_time: Duration::from_secs(1),
                    evaluator_compile_time: Duration::from_secs(1),
                    ..GraphGenerationStats::default()
                },
            },
            GeneratedGraphReport {
                process_id: 0,
                integrand_name: "itg".to_string(),
                graph_name: "GL02".to_string(),
                stats: GraphGenerationStats {
                    evaluator_count: 3,
                    total_time: Duration::from_secs(6),
                    evaluator_spenso_time: Duration::from_secs(2),
                    evaluator_symbolica_time: Duration::ZERO,
                    evaluator_compile_time: Duration::from_secs(3),
                    ..GraphGenerationStats::default()
                },
            },
        ];

        let summary = render_generation_summary(&reports, 0, Some(4), None).unwrap();
        let plain = Regex::new(r"\x1b\[[0-9;]*m")
            .unwrap()
            .replace_all(&summary, "")
            .into_owned();

        assert!(plain.contains("Total"));
        assert!(plain.contains("5"));
        assert!(plain.contains("20%"));
        assert!(plain.contains("30%"));
        assert!(plain.contains("10%"));
        assert!(plain.contains("40%"));

        let lines = plain.lines().collect::<Vec<_>>();
        let total_line_index = lines
            .iter()
            .position(|line| line.contains("Total"))
            .expect("missing Total row");
        let separator_line = lines[..total_line_index]
            .iter()
            .rev()
            .find(|line| !line.trim().is_empty())
            .expect("missing line before Total row");
        assert!(
            separator_line.chars().any(|ch| "┼╪╫┿├┤".contains(ch)),
            "missing horizontal separator before Total row: {separator_line}"
        );
    }

    fn base_args(tokens: &str) -> SpecArgs {
        SpecArgs {
            tokens: tokens.split_whitespace().map(str::to_owned).collect(),
            only_diagrams: true,
            ..SpecArgs::default()
        }
    }

    // Every fixture states its loop domain explicitly, just like a real CLI request.
    fn parse_ok_amp(s: &str) -> ProcessSpec {
        test_initialise().unwrap();
        parse_spec_with_model(
            &base_args(s),
            GenerationType::Amplitude,
            &load_generic_model("sm"),
        )
        .unwrap()
    }
    fn parse_ok_xs(s: &str) -> ProcessSpec {
        test_initialise().unwrap();
        parse_spec_with_model(
            &base_args(s),
            GenerationType::CrossSection,
            &load_generic_model("sm"),
        )
        .unwrap()
    }
    fn pdgs(
        selectors: impl IntoIterator<Item = impl std::borrow::Borrow<ParticleSelector>>,
    ) -> std::collections::BTreeSet<i64> {
        let model = load_generic_model("sm");
        selectors
            .into_iter()
            .map(|p| {
                model
                    .particle_by_id(p.borrow().resolve(&model).unwrap())
                    .unwrap()
                    .pdg_code
            })
            .collect()
    }

    #[test]
    fn zero_color_filter_is_independent_of_grouping_and_forwarded_in_both_modes() {
        let model = load_generic_model("sm");
        for (mode, kind, loops) in [
            ("amp", GenerationType::Amplitude, "[{0}]"),
            ("xs", GenerationType::CrossSection, "[{{1}}]"),
        ] {
            for enabled in [false, true] {
                let mut command = vec!["gammaloop", "generate", mode, "a", ">", "d", "d~", loops];
                if enabled {
                    command.push("--filter-zero-color");
                }
                let repl = Repl::try_parse_from(command).unwrap();
                let Commands::Generate(Generate {
                    mode: Some(GenerateCmd::Amp(args) | GenerateCmd::Xs(args)),
                    ..
                }) = repl.command
                else {
                    panic!("expected diagram generation")
                };
                assert_eq!(args.filter_zero_color, enabled);
                assert!(matches!(build_grouping_option(&args), Grouping::None));
                let decoded: SpecArgs =
                    serde_json::from_value(serde_json::to_value(&args).unwrap()).unwrap();
                assert_eq!(decoded, args);
                let spec = parse_spec_with_model(&decoded, kind, &model).unwrap();
                assert_eq!(
                    spec.process_definition.numerators.filter_zero_color,
                    enabled
                );
            }
        }
    }

    #[test]
    fn cp_symmetrization_option_defaults_false_and_roundtrips() {
        let model = load_generic_model("sm");
        for option in [None, Some(false), Some(true)] {
            let mut command = vec![
                "gammaloop",
                "generate",
                "xs",
                "a",
                ">",
                "d",
                "d~",
                "[{{1}}]",
            ];
            if let Some(enabled) = option {
                command.extend([
                    "--symmetrize-left-right-states",
                    if enabled { "true" } else { "false" },
                ]);
            }
            let repl = Repl::try_parse_from(command).unwrap();
            let Commands::Generate(Generate {
                mode: Some(GenerateCmd::Xs(args)),
                ..
            }) = repl.command
            else {
                panic!("expected cross-section generation")
            };
            assert_eq!(args.symmetrize_left_right_states, option);
            let serialized = serde_json::to_value(&args).unwrap();
            assert_eq!(
                serialized["symmetrize_left_right_states"],
                serde_json::json!(option)
            );
            let decoded: SpecArgs = serde_json::from_value(serialized).unwrap();
            assert_eq!(decoded, args);
            let spec =
                parse_spec_with_model(&decoded, GenerationType::CrossSection, &model).unwrap();
            assert_eq!(
                spec.process_definition.request.cuts().symmetrize_left_right,
                option.unwrap_or(false)
            );
        }
    }

    #[test]
    fn cross_sections_default_to_unsymmetrized_external_states() {
        let spec = parse_ok_xs("e+ e- > d d~ [{{1}}]");
        assert_eq!(
            spec.process_definition.request.symmetrization(),
            Symmetrize::default()
        );
    }

    #[test]
    fn amplitude_left_right_shorthand_respects_explicit_side_options() {
        let model = load_generic_model("sm");
        let mut args = base_args("e+ e- > d d~ [{1}]");
        for (lr, initial, final_state, expected) in [
            (None, None, None, (false, false)),
            (Some(true), None, None, (true, true)),
            (Some(true), Some(false), None, (false, true)),
            (Some(true), None, Some(false), (true, false)),
        ] {
            args.symmetrize_left_right_states = lr;
            args.symmetrize_initial_states = initial;
            args.symmetrize_final_states = final_state;
            let request = parse_spec_with_model(&args, GenerationType::Amplitude, &model)
                .unwrap()
                .process_definition
                .request;
            assert_eq!(
                (
                    request.symmetrization().initial,
                    request.symmetrization().final_state
                ),
                expected
            );
            // The amplitude shorthand is represented by the two external-state settings.
            assert!(!request.cuts().symmetrize_left_right);
        }
    }

    #[test]
    fn fermion_loop_filter_is_added_in_both_generation_modes() {
        let model = load_generic_model("sm");
        for (mode, kind, loops) in [
            ("amp", GenerationType::Amplitude, "[{1}]"),
            ("xs", GenerationType::CrossSection, "[{{1}}]"),
        ] {
            let repl = Repl::try_parse_from([
                "gammaloop",
                "generate",
                mode,
                "e+",
                "e-",
                ">",
                "d",
                "d~",
                loops,
                "--number-of-fermion-loops",
                "1",
                "2",
            ])
            .unwrap();
            let Commands::Generate(Generate {
                mode: Some(GenerateCmd::Amp(args) | GenerateCmd::Xs(args)),
                ..
            }) = repl.command
            else {
                panic!("expected diagram generation")
            };
            let request = parse_spec_with_model(&args, kind, &model)
                .unwrap()
                .process_definition
                .request;
            assert_eq!(request.topology().fermion_loops, Range::from(1..=2));
        }
    }

    #[test]
    fn allowed_vertex_interactions_are_added_to_generation_filters() {
        let model = load_generic_model("sm");
        let mut args = base_args("e+ e- > d d~ [{1}]");
        args.allowed_vertex_interactions = Some(vec!["V_6".to_owned(), "V_9".to_owned()]);
        let request = parse_spec_with_model(&args, GenerationType::Amplitude, &model)
            .unwrap()
            .process_definition
            .request;
        assert_eq!(
            request.restrictions().vertex_allow,
            vec![[
                VertexSelector::Name("V_6".to_owned()),
                VertexSelector::Name("V_9".to_owned())
            ]
            .into_iter()
            .collect::<std::collections::BTreeSet<_>>()]
        );
    }

    #[test]
    fn basic_process_list_arrow() {
        let ps = parse_ok_amp("e+ e- > d d~ g [{1}]");
        let p = ps.process_definition.request.process();
        assert_eq!(
            p.incoming(),
            &[ParticleSelector::from("e+"), ParticleSelector::from("e-")]
        );
        assert_eq!(
            pdgs(p.incoming()),
            std::collections::BTreeSet::from([-11, 11])
        );
        assert_eq!(
            p.outgoing_pdgs(&load_generic_model("sm")).unwrap(),
            vec![vec![1, -1, 21]]
        );
    }
    #[test]
    fn basic_process_list_to() {
        let ps = parse_ok_amp("e+ e- to d d~ g [{1}]");
        assert_eq!(
            ps.process_definition.request,
            parse_ok_amp("e+ e- > d d~ g [{1}]")
                .process_definition
                .request
        );
    }
    #[test]
    fn empty_sets_allowed() {
        let ps = parse_ok_xs("{} to {} [{{1}}]");
        assert!(ps
            .process_definition
            .request
            .process()
            .incoming()
            .is_empty());
        assert_eq!(
            ps.process_definition
                .request
                .process()
                .outgoing_alternatives(),
            &[Vec::<ParticleSelector>::new()]
        );
    }
    #[test]
    fn alternatives_in_final_states() {
        let ps = parse_ok_xs("e+ e- > { Z Z, a a } [{{1}}]");
        assert_eq!(
            ps.process_definition
                .request
                .process()
                .outgoing_pdgs(&load_generic_model("sm"))
                .unwrap(),
            vec![vec![23; 2], vec![22; 2]]
        );
    }
    #[test]
    fn veto_and_only() {
        let ps = parse_ok_xs("e+ e- > mu+ mu- / u d c g ghG e- | u g ghG [{{1}}]");
        let r = ps.process_definition.request.restrictions();
        let veto = pdgs(&r.particle_exclude);
        assert!(veto.contains(&2));
        assert!(veto.contains(&11));
        let only = pdgs(&r.particle_allow[0]);
        for pdg in [2, 21, 9000005] {
            assert!(only.contains(&pdg));
        }
    }
    #[test]
    fn amp_order_constraints() {
        let ps = parse_ok_amp("e+ e- > mu+ mu- QED==2 QCD>=2 QCD<=4 [{1}]");
        assert_eq!(
            ps.process_definition.request.orders().plain["QCD"],
            Range::from(2..=4)
        );
        assert_eq!(
            ps.process_definition.request.orders().plain["QED"],
            Range::exact(2)
        );
    }
    #[test]
    fn xs_order_constraints_powered() {
        let ps = parse_ok_xs("e+ e- > mu+ mu- QED^2==2 QCD^2>=2 QCD^2<=4 [{{1}}]");
        assert_eq!(
            ps.process_definition.request.orders().squared["QCD"],
            Range::from(2..=4)
        );
        assert_eq!(
            ps.process_definition.request.orders().squared["QED"],
            Range::exact(2)
        );
    }
    #[test]
    fn perturbative_block_all_forms() {
        let ps = parse_ok_xs("e+ e- > z [ {1} {{2}} QCD=2 QED=1 ]");
        let r = &ps.process_definition.request;
        assert_eq!(r.cuts().amplitude.loops, Range::exact(1));
        assert_eq!(r.loops(), &Loops::forward(2..=2));
        assert_eq!(r.orders().perturbative.get("QCD"), Some(&2));
        assert_eq!(r.orders().perturbative.get("QED"), Some(&1));
    }
    #[test]
    fn perturbative_block_shorthand() {
        let ps = parse_ok_xs("e+ e- > z [ {{1}} QCD ]");
        assert_eq!(
            ps.process_definition
                .request
                .orders()
                .perturbative
                .get("QCD"),
            Some(&1)
        );
    }
    #[test]
    fn amplitude_generation_type_is_set() {
        assert_eq!(
            parse_ok_amp("e+ e- > mu+ mu- [{1}]")
                .process_definition
                .request
                .kind(),
            GenerationType::Amplitude
        );
    }
    #[test]
    fn rejects_malformed_tokens() {
        let err = parse_spec_with_model(
            &base_args("e+ > z [ {{}} ]"),
            GenerationType::CrossSection,
            &load_generic_model("sm"),
        )
        .unwrap_err();
        assert!(matches!(err,ParseError::Request(RequestParseError::InvalidToken(s)) if s=="{{}}"));
    }
    #[test]
    fn missing_arrow_detected() {
        let err = parse_spec_with_model(
            &base_args("e+ e- z z"),
            GenerationType::CrossSection,
            &load_generic_model("sm"),
        )
        .unwrap_err();
        assert!(matches!(
            err,
            ParseError::Request(RequestParseError::MissingArrow)
        ));
    }
    #[test]
    fn vacuum_defaults_and_filters_xs() {
        let ps = parse_ok_xs("{} to {} [{{1}}]");
        let r = &ps.process_definition.request;
        assert_eq!(r.kind(), GenerationType::CrossSection);
        // Vacuum requests use the same minimal defaults as other requests.
        assert_eq!(r.topology(), &TopologyRules::default());
        assert_eq!(r.cuts(), &feynkit_generator::Cuts::default());
        assert!(matches!(ps.process_definition.grouping, Grouping::None));
    }
    #[test]
    fn particle_veto_resolution_cross_section() {
        // Veto both particles and antiparticles, plus a charged lepton.
        let ps = parse_ok_xs("e+ e- > mu+ mu- / u d u~ e+ [{{1}}]");
        let set = pdgs(
            &ps.process_definition
                .request
                .restrictions()
                .particle_exclude,
        );
        for pdg in [2, 1, -2, -11] {
            assert!(set.contains(&pdg), "missing PDG {pdg} in veto");
        }
    }
    #[test]
    fn particle_inclusion_resolution_cross_section() {
        let ps = parse_ok_xs("e+ e- > mu+ mu- | u d u~ e+ [{{1}}]");
        let r = ps.process_definition.request.restrictions();
        let set = pdgs(&r.particle_allow[0]);
        assert_eq!(set, std::collections::BTreeSet::from([-11, -2, 1, 2]));
        // The allow clause is kept explicitly instead of materializing its complement.
        for pdg in [
            3, 4, 5, 6, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 250, 251, 9000001, 9000002,
            9000003, 9000004, 9000005,
        ] {
            assert!(!set.contains(&pdg), "unexpected allowed PDG {pdg}");
        }
        assert!(r.particle_exclude.is_empty());
    }
    #[test]
    fn final_sets_and_orders_and_loops_cross_section() {
        let ps = parse_ok_xs("e+ e- > { Z Z, a a, H H } [ {{3}} QCD=2 QED=1 ]");
        let r = &ps.process_definition.request;
        assert_eq!(r.process().outgoing_alternatives().len(), 3);
        assert_eq!(r.kind(), GenerationType::CrossSection);
        assert_eq!(r.loops(), &Loops::forward(3..=3));
        assert_eq!(r.orders().perturbative.get("QCD"), Some(&2));
        assert_eq!(r.orders().perturbative.get("QED"), Some(&1));
    }
    #[test]
    fn amplitude_with_powered_xs_constraints_and_only() {
        let ps = parse_ok_amp("e+ e- > mu+ mu- | g ghG QED^2==2 QCD^2>=2 QCD^2<=4 [ {1} ]");
        let r = &ps.process_definition.request;
        assert_eq!(r.loops(), &Loops::total(1..=1));
        assert!(r.orders().perturbative.is_empty());
        assert_eq!(r.orders().squared["QED"], Range::exact(2));
        assert_eq!(r.orders().squared["QCD"], Range::from(2..=4));
        let only = pdgs(&r.restrictions().particle_allow[0]);
        assert!(only.contains(&21));
        assert!(only.contains(&9000005));
    }
    #[test]
    fn amplitude_case_insensitive_and_anti_in_veto() {
        // Mixed case names and antiparticles resolve to the same model entries.
        let ps = parse_ok_amp("E+ e- to Z Z / u~ d~ A [{1}]");
        let set = pdgs(
            &ps.process_definition
                .request
                .restrictions()
                .particle_exclude,
        );
        for pdg in [-2, -1, 22] {
            assert!(set.contains(&pdg), "missing PDG {pdg} in AMP veto");
        }
    }
    #[test]
    fn shell_name_long_with_sets_veto_orders_pert_and_xs() {
        let model = load_generic_model("sm");
        // Alternatives belong to cross sections; amplitudes require a single final state.
        let ps = parse_ok_xs(
            "e+ e- > { Z Z, a a } / u d QED==2 QCD>=1 QED^2<=4 QCD^2>=2 [ {{1}} QED=1 QCD=2 ]",
        );
        assert_eq!(
            ps.process_shell_name(&model, false),
            "epem_zz_or_aa-no__d_u-QCD_ge_1__QED_eq_2-QCDsq_ge_2__QEDsq_le_4-QCDloop_eq_2__QEDloop_eq_1"
        );
    }
    #[test]
    fn repr_str_canonical_with_amp_xs_veto_loops_orders() {
        let model = load_generic_model("sm");
        let ps = parse_ok_xs(
            "mu+ mu- > g g / ghg gha QED==1 QCD<=3 [ {{2}} QED=0 QCD=1 ] QED^2==0 QCD^2>=2",
        );
        // Explicit zero orders remain in the canonical representation.
        assert_eq!(
            ps.repr_str(&model),
            "mu+ mu- > g g / gha ghg QCD<=3 QED==1 [ {{2}} QCD QED=0 ] QCD^2>=2 QED^2==0"
        );
    }
    #[test]
    fn repr_str_with_both_loop_counts_and_sorted_orders() {
        let model = load_generic_model("sm");
        let ps = parse_ok_xs("e+ e- > z z [ {{2}} QED=2 QCD=1 {1} ]");
        assert_eq!(ps.repr_str(&model), "e+ e- > z z [ {1} {{2}} QCD QED=2 ]");
    }
    #[test]
    fn shell_name_short_canonical_and_sanitization() {
        let model = load_generic_model("sm");
        let ps = parse_ok_xs("W+ W- > {} [{{1}}]");
        // Short form uses only the base slug and sanitizes +/- and ~.
        assert_eq!(ps.process_shell_name(&model, true), "wpwm_empty");
    }
    #[test]
    fn all_cli_requests_require_loops_and_reject_the_wrong_domain() {
        let model = load_generic_model("sm");
        for kind in [GenerationType::Amplitude, GenerationType::CrossSection] {
            for text in ["e+ e- > mu+ mu-", "{} > {}"] {
                assert!(matches!(
                    parse_spec_with_model(&base_args(text), kind, &model),
                    Err(ParseError::Request(RequestParseError::MissingLoops))
                ));
            }
        }
        assert!(matches!(
            parse_spec_with_model(
                &base_args("{} > {} [{{1}}]"),
                GenerationType::Amplitude,
                &model
            ),
            Err(ParseError::Request(RequestParseError::WrongLoopDomain))
        ));
    }
    #[test]
    fn explicit_empty_allow_lists_and_veto_qualifiers_are_preserved() {
        let model = load_generic_model("sm");
        let mut args = base_args("e+ e- > d d~ [{1}] |");
        args.allowed_vertex_interactions = Some(vec![]);
        args.veto_self_energy_of_massive_lines = Some(true);
        let r = parse_spec_with_model(&args, GenerationType::Amplitude, &model)
            .unwrap()
            .process_definition
            .request;
        assert_eq!(
            r.restrictions().vertex_allow,
            vec![std::collections::BTreeSet::<VertexSelector>::new()]
        );
        assert_eq!(
            r.restrictions().particle_allow,
            vec![std::collections::BTreeSet::<ParticleSelector>::new()]
        );
        assert_eq!(r.topology().self_energies, Veto::massive());
    }
}
