use std::path::PathBuf;

use crate::{
    completion::CompletionArgExt,
    state::{ProcessRef, State},
    CLISettings,
};
use color_eyre::Result;
use eyre::eyre;
use gammalooprs::{
    integrands::process::{
        ir::{IRProfileSetting, IrLimitTestReport},
        OrientationProfileMode, ProcessIntegrand,
    },
    processes::ProcessCollection,
    uv::{
        profile::{ProfileSettings, UVLimitSelection, UVProfileFixedRay, UVProfileable},
        UVProfileAnalysis,
    },
};
use linnet::half_edge::involution::EdgeIndex;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use tracing::{info, instrument};

use clap::{Args, Subcommand};

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Profile {
    /// Ultraviolet profile analysis
    UltraViolet(#[command(flatten)] UltraVioletProfile),
    /// Bulk profile analysis
    #[command(name = "bulk")]
    InfraRed(#[command(flatten)] InfraRedProfile),
}

#[derive(Args, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct UltraVioletProfile {
    /// Process reference: #<id>, name:<name>, or <id>/<name>
    #[arg(
        short = 'p',
        long = "process",
        value_name = "PROCESS",
        completion_process_selector(crate::completion::SelectorKind::Any)
    )]
    pub process: Option<ProcessRef>,

    /// The integrand name to inspect
    #[arg(
        short = 'i',
        long = "integrand-name",
        value_name = "NAME",
        completion_integrand_selector(crate::completion::SelectorKind::Any)
    )]
    pub integrand_name: Option<String>,

    /// Restrict profiling to this graph
    #[arg(
        short = 'g',
        long = "graph",
        value_name = "GRAPH",
        completion_selected_master_graph()
    )]
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub graph: Option<String>,

    /// Restrict a cross section to the Cutkosky cut with these edge IDs
    #[arg(
        long = "cutkosky-cut",
        visible_alias = "cut-edges",
        value_name = "EDGE",
        num_args = 1..,
        value_delimiter = ',',
        requires = "graph"
    )]
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub cutkosky_cut: Vec<usize>,

    /// UV limits to profile
    #[arg(long = "selected-limits", value_enum, default_value = "only-divergent")]
    #[serde(default)]
    pub selected_limits: UVLimitSelection,

    /// Number of scaling points to sample
    #[arg(long = "n-points", default_value_t = 20)]
    pub n_points: usize,

    /// Minimum base-10 scaling exponent
    #[arg(long = "min-scaling", default_value_t = 3.0)]
    pub min_scale_exponent: f64,

    /// Maximum base-10 scaling exponent
    #[arg(long = "max-scaling", default_value_t = 6.0)]
    pub max_scale_exponent: f64,

    /// Use f128 precision for evaluation
    #[arg(long = "use_f128")]
    pub use_f128: bool,

    /// Also derive analytic UV series (amplitudes only)
    #[arg(long = "analyse_analytically")]
    pub analyse_analytically: bool,

    /// Profile each visible orientation separately and include per-orientation results
    #[arg(long = "per-orientation")]
    pub per_orientation: bool,

    /// Random seed for momentum sampling
    #[arg(long = "seed")]
    pub seed: Option<u64>,

    /// Fixed UV ray directions as flattened 3-vectors, repeated for all rays if one direction is supplied
    #[arg(
        long = "uv-ray-directions",
        num_args = 1..,
        value_delimiter = ',',
        allow_negative_numbers = true
    )]
    pub uv_ray_directions: Vec<f64>,

    /// Fixed UV ray starting norms, repeated for all rays if one norm is supplied
    #[arg(
        long = "uv-ray-norms",
        num_args = 1..,
        value_delimiter = ',',
        allow_negative_numbers = true,
        requires = "uv_ray_directions"
    )]
    pub uv_ray_norms: Vec<f64>,

    /// Output directory for uv_profile.json (optional)
    #[arg(short = 'o', long = "output", value_hint = clap::ValueHint::DirPath)]
    pub output_file: Option<PathBuf>,
}

#[derive(Args, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct InfraRedProfile {
    /// Process reference: #<id>, name:<name>, or <id>/<name>
    #[arg(
        short = 'p',
        long = "process",
        value_name = "PROCESS",
        completion_process_selector(crate::completion::SelectorKind::CrossSection)
    )]
    pub process: Option<ProcessRef>,

    /// The cross-section name to inspect
    #[arg(
        short = 'i',
        long = "integrand-name",
        value_name = "NAME",
        completion_integrand_selector(crate::completion::SelectorKind::CrossSection)
    )]
    pub integrand_name: Option<String>,

    /// Number of scaling points to sample
    #[arg(long = "n-points", default_value_t = 20)]
    pub n_points: usize,

    /// Minimum scaling factor
    #[arg(long = "min-scaling", default_value_t = -2.0)]
    pub min_scale_exponent: f64,

    /// Maximum scaling factor
    #[arg(long = "max-scaling", default_value_t = -3.0)]
    pub max_scale_exponent: f64,

    /// Random seed for momentum sampling
    #[arg(long = "seed")]
    pub seed: Option<u64>,

    /// Output file for results (optional)
    #[arg(short = 'o', long = "output", value_hint = clap::ValueHint::FilePath)]
    pub output_file: Option<PathBuf>,

    /// restrict test to particular graphs or limits
    #[arg(short = 's', long = "select")]
    pub select: Option<String>,

    /// Profile each visible orientation separately and report one row per orientation
    #[arg(long = "per-orientation")]
    pub per_orientation: bool,

    /// Retain generated event data so per-cut IR fits can be computed and displayed
    #[arg(long = "show-per-cut-info")]
    pub show_per_cut_info: bool,
}

impl Default for InfraRedProfile {
    fn default() -> Self {
        Self {
            process: None,
            integrand_name: None,
            n_points: 20,
            min_scale_exponent: -2.0,
            max_scale_exponent: -3.0,
            seed: None,
            output_file: None,
            select: None,
            per_orientation: false,
            show_per_cut_info: false,
        }
    }
}
impl Default for UltraVioletProfile {
    fn default() -> Self {
        Self {
            process: None,
            integrand_name: None,
            graph: None,
            cutkosky_cut: Vec::new(),
            selected_limits: UVLimitSelection::OnlyDivergent,
            n_points: 20,
            min_scale_exponent: 3.0,
            max_scale_exponent: 6.0,
            use_f128: false,
            analyse_analytically: false,
            per_orientation: false,
            seed: None,
            uv_ray_directions: Vec::new(),
            uv_ray_norms: Vec::new(),
            output_file: None,
        }
    }
}

pub enum ProfileResult {
    UltraViolet(UVProfileAnalysis),
    InfraRed(IrLimitTestReport),
}

impl ProfileResult {
    pub fn unwrap_uv(self) -> UVProfileAnalysis {
        match self {
            ProfileResult::UltraViolet(uv_analyis) => uv_analyis,
            _ => panic!("result does not contain uv profilel analysis data"),
        }
    }

    pub fn unwrap_ir(self) -> IrLimitTestReport {
        match self {
            ProfileResult::InfraRed(ir_analyis) => ir_analyis,
            _ => panic!("result does not contain ir profile analysis data"),
        }
    }
}

impl Profile {
    #[instrument(skip_all)]
    #[allow(clippy::needless_update)]
    pub fn run(
        &self,
        state: &mut State,
        global_cli_settings: &CLISettings,
    ) -> Result<ProfileResult> {
        match self {
            Profile::UltraViolet(UltraVioletProfile {
                process,
                integrand_name,
                graph,
                cutkosky_cut,
                selected_limits,
                n_points,
                min_scale_exponent,
                max_scale_exponent,
                use_f128,
                seed,
                uv_ray_directions,
                uv_ray_norms,
                analyse_analytically,
                per_orientation,
                output_file,
            }) => {
                let (process_id, integrand_name) =
                    state.find_integrand_ref(process.as_ref(), integrand_name.as_ref())?;
                let model = state.resolve_model_for_integrand(process_id, &integrand_name)?;
                let (default_uv_ray_norm, graph_id) = {
                    let integrand = state
                        .process_list
                        .get_integrand_mut(process_id, &integrand_name)?;
                    integrand.warm_up(&model)?;
                    let graph_id = graph
                        .as_ref()
                        .map(|graph_selector| {
                            let numeric_selector =
                                graph_selector.strip_prefix('#').unwrap_or(graph_selector);
                            if let Some(graph_id) = numeric_selector
                                .parse::<usize>()
                                .ok()
                                .filter(|graph_id| *graph_id < integrand.graph_count())
                            {
                                return Ok(graph_id);
                            }
                            if let Some(graph_id) = integrand.find_graph_id_by_name(graph_selector)
                            {
                                return Ok(graph_id);
                            }

                            let available = (0..integrand.graph_count())
                                .filter_map(|graph_id| {
                                    integrand
                                        .graph_name_by_id(graph_id)
                                        .map(|name| format!("{graph_id}:{name}"))
                                })
                                .collect::<Vec<_>>()
                                .join(", ");
                            Err(eyre!(
                                "No graph '{}' exists in integrand '{}'. Available graphs: {}",
                                graph_selector,
                                integrand_name,
                                available
                            ))
                        })
                        .transpose()?;
                    if !cutkosky_cut.is_empty()
                        && matches!(integrand, ProcessIntegrand::Amplitude(_))
                    {
                        return Err(eyre!(
                            "Cutkosky-cut selection is only supported for cross sections"
                        ));
                    }
                    (integrand.get_settings().kinematics.e_cm, graph_id)
                };

                let fixed_uv_ray = if uv_ray_directions.is_empty() {
                    None
                } else {
                    let uv_ray_norms = if uv_ray_norms.is_empty() {
                        vec![default_uv_ray_norm]
                    } else {
                        uv_ray_norms.clone()
                    };
                    Some(UVProfileFixedRay::from_flat_components(
                        uv_ray_directions,
                        &uv_ray_norms,
                    )?)
                };

                let profile_settings = ProfileSettings {
                    n_points: *n_points,
                    min_scale_exponent: *min_scale_exponent,
                    max_scale_exponent: *max_scale_exponent,
                    seed: (*seed).unwrap_or(42),
                    use_f128: *use_f128,
                    analyse_analytically: *analyse_analytically,
                    orientation_mode: if *per_orientation {
                        OrientationProfileMode::PerOrientation
                    } else {
                        OrientationProfileMode::Summed
                    },
                    fixed_uv_ray,
                    graph_id,
                    cutkosky_cut: (!cutkosky_cut.is_empty())
                        .then(|| cutkosky_cut.iter().copied().map(EdgeIndex::from).collect()),
                    selected_limits: *selected_limits,
                    ..Default::default()
                };
                let profile_res = {
                    let process = &mut state.process_list.processes[process_id];
                    match &mut process.collection {
                        ProcessCollection::Amplitudes(amplitudes) => amplitudes
                            .get_mut(&integrand_name)
                            .ok_or_else(|| {
                                eyre!(
                                    "No amplitude named '{}' in process '{}'",
                                    integrand_name,
                                    process.definition.folder_name
                                )
                            })?
                            .profile(&model, &profile_settings)?,
                        ProcessCollection::CrossSections(cross_sections) => cross_sections
                            .get_mut(&integrand_name)
                            .ok_or_else(|| {
                                eyre!(
                                    "No cross section named '{}' in process '{}'",
                                    integrand_name,
                                    process.definition.folder_name
                                )
                            })?
                            .profile(&model, &profile_settings)?,
                    }
                }
                .analyse();

                for t in profile_res.tables_per_graph(-0.9) {
                    info!("\n{}", t);
                }

                for t in profile_res.analytic_tables_per_graph() {
                    let Some(t) = t else {
                        continue;
                    };
                    info!("\n{}", t);
                }

                for t in profile_res.per_orientation_tables_per_graph(-0.9) {
                    let Some(t) = t else {
                        continue;
                    };
                    info!("\n{}", t);
                }

                info!("\n{}", profile_res.pass_fail(-0.9));

                if let Some(file) = output_file {
                    global_cli_settings
                        .ensure_write_target_outside_active_state(file, "write profile output")?;
                    profile_res.write_profile_data(file)?
                }

                Ok(ProfileResult::UltraViolet(profile_res))
            }
            Profile::InfraRed(InfraRedProfile {
                process,
                integrand_name,
                n_points,
                min_scale_exponent,
                max_scale_exponent,
                seed,
                output_file: _,
                select,
                per_orientation,
                show_per_cut_info,
            }) => {
                let ir_profile_settings = IRProfileSetting {
                    lambda_exp_start: *min_scale_exponent,
                    lambda_exp_end: *max_scale_exponent,
                    steps: *n_points,
                    seed: seed.unwrap_or(420),
                    select_limits_and_graphs: select.clone(),
                    orientation_mode: if *per_orientation {
                        OrientationProfileMode::PerOrientation
                    } else {
                        OrientationProfileMode::Summed
                    },
                    show_per_cut_info: *show_per_cut_info,
                };

                let (process_id, integrand_name) =
                    state.find_integrand_ref(process.as_ref(), integrand_name.as_ref())?;

                let model = state.resolve_model_for_integrand(process_id, &integrand_name)?;

                let process = &mut state.process_list.processes[process_id];
                let integrand = match &mut process.collection {
                    ProcessCollection::Amplitudes(amplitudes) => amplitudes
                        .get_mut(&integrand_name)
                        .ok_or_else(|| {
                            eyre!(
                                "No amplitude named '{}' in process '{}'",
                                integrand_name,
                                process.definition.folder_name
                            )
                        })?
                        .integrand
                        .as_mut(),
                    ProcessCollection::CrossSections(xs) => xs
                        .get_mut(&integrand_name)
                        .ok_or_else(|| {
                            eyre!(
                                "No xs named '{}' in process '{}'",
                                integrand_name,
                                process.definition.folder_name
                            )
                        })?
                        .integrand
                        .as_mut(),
                };

                let profile_result = match integrand.ok_or(eyre!(
                    "Integrand {} has not yet been generated",
                    integrand_name
                ))? {
                    ProcessIntegrand::CrossSection(cross_section_integrand) => {
                        cross_section_integrand.test_ir(&ir_profile_settings, &model)?
                    }
                    ProcessIntegrand::Amplitude(amplitude_integrand) => {
                        amplitude_integrand.test_ir(&ir_profile_settings, &model)?
                    }
                };

                info!("\n{}", profile_result);
                Ok(ProfileResult::InfraRed(profile_result))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use clap::Parser;
    use gammalooprs::uv::profile::UVLimitSelection;

    use crate::{commands::Commands, Repl};

    use super::Profile;

    #[test]
    fn uv_profile_cli_defaults_to_only_divergent_limits() {
        let repl = Repl::try_parse_from(["gammaloop", "profile", "ultra-violet"]).unwrap();
        let Commands::Profile(Profile::UltraViolet(profile)) = repl.command else {
            panic!("expected an ultraviolet profile command");
        };

        assert_eq!(profile.selected_limits, UVLimitSelection::OnlyDivergent);
        assert_eq!(profile.graph, None);
        assert!(profile.cutkosky_cut.is_empty());
    }

    #[test]
    fn uv_profile_cli_parses_graph_cut_and_all_limits() {
        let repl = Repl::try_parse_from([
            "gammaloop",
            "profile",
            "ultra-violet",
            "--graph",
            "GL2",
            "--cutkosky-cut",
            "5,2",
            "--selected-limits",
            "all",
        ])
        .unwrap();
        let Commands::Profile(Profile::UltraViolet(profile)) = repl.command else {
            panic!("expected an ultraviolet profile command");
        };

        assert_eq!(profile.selected_limits, UVLimitSelection::All);
        assert_eq!(profile.graph.as_deref(), Some("GL2"));
        assert_eq!(profile.cutkosky_cut, [5, 2]);

        let alias = Repl::try_parse_from([
            "gammaloop",
            "profile",
            "ultra-violet",
            "--graph",
            "GL2",
            "--cut-edges",
            "5,2",
        ])
        .unwrap();
        let Commands::Profile(Profile::UltraViolet(alias)) = alias.command else {
            panic!("expected an ultraviolet profile command");
        };
        assert_eq!(alias.cutkosky_cut, [5, 2]);
    }

    #[test]
    fn uv_profile_cli_requires_a_graph_for_a_cut() {
        let result = Repl::try_parse_from([
            "gammaloop",
            "profile",
            "ultra-violet",
            "--cutkosky-cut",
            "5,2",
        ]);

        assert!(result.is_err());
    }
}
