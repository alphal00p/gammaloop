use std::path::PathBuf;

use crate::{
    completion::CompletionArgExt,
    state::{ProcessListExt, ProcessRef, State},
    CLISettings,
};
use color_eyre::Result;
use eyre::eyre;
use gammalooprs::{
    integrands::process::ir::{IRProfileSetting, IrLimitTestReport},
    integrands::process::ProcessIntegrand,
    uv::{
        profile::{ProfileSettings, UVProfileable},
        UVProfileAnalysis,
    },
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use tracing::{info, instrument};

use clap::{Args, Subcommand};

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Profile {
    /// Ultraviolet profile analysis
    UltraViolet(#[command(flatten)] UltraVioletProfile),
    /// Infrared profile analysis
    InfraRed(#[command(flatten)] InfraRedProfile),
}

#[derive(Args, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct UltraVioletProfile {
    /// Process reference: #<id>, name:<name>, or <id>/<name>
    #[arg(
        short = 'p',
        long = "process",
        value_name = "PROCESS",
        completion_process_selector(crate::completion::SelectorKind::Amplitude)
    )]
    pub process: Option<ProcessRef>,

    /// The amplitude name to inspect
    #[arg(
        short = 'i',
        long = "integrand-name",
        value_name = "NAME",
        completion_integrand_selector(crate::completion::SelectorKind::Amplitude)
    )]
    pub integrand_name: Option<String>,

    /// Number of scaling points to sample
    #[arg(long = "n-points", default_value_t = 20)]
    pub n_points: usize,

    /// Minimum scaling factor
    #[arg(long = "min-scaling", default_value_t = 3.0)]
    pub min_scale_exponent: f64,

    /// Maximum scaling factor
    #[arg(long = "max-scaling", default_value_t = 5.0)]
    pub max_scale_exponent: f64,

    /// Use f128 precision for evaluation
    #[arg(long = "use_f128")]
    pub use_f128: bool,

    #[arg(long = "analyse_analytically")]
    pub analyse_analytically: bool,

    /// Random seed for momentum sampling
    #[arg(long = "seed")]
    pub seed: Option<u64>,

    /// Output file for results (optional)
    #[arg(short = 'o', long = "output", value_hint = clap::ValueHint::FilePath)]
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
}

impl Default for UltraVioletProfile {
    fn default() -> Self {
        Self {
            process: None,
            integrand_name: None,
            n_points: 20,
            min_scale_exponent: 1.0,
            max_scale_exponent: 2.0,
            use_f128: false,
            analyse_analytically: false,
            seed: None,
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
}

impl Profile {
    #[instrument(skip_all)]
    pub fn run(
        &self,
        state: &mut State,
        _global_cli_settings: &CLISettings,
    ) -> Result<ProfileResult> {
        match self {
            Profile::UltraViolet(UltraVioletProfile {
                process,
                integrand_name,
                n_points,
                min_scale_exponent,
                max_scale_exponent,
                use_f128,
                seed,
                analyse_analytically,
                output_file,
            }) => {
                let (process_id, integrand_name) =
                    state.find_integrand_ref(process.as_ref(), integrand_name.as_ref())?;
                let model = state.resolve_model_for_integrand(process_id, &integrand_name)?;
                let process_ref = ProcessRef::Id(process_id);
                let amplitude = state
                    .process_list
                    .get_amplitude_mut_ref(Some(&process_ref), Some(&integrand_name))?;

                amplitude
                    .integrand
                    .as_mut()
                    .ok_or(eyre!(
                        "Integrand {} has not yet been generated, but exists",
                        amplitude.name
                    ))?
                    .warm_up(&model)?;

                let profile_settings = ProfileSettings {
                    n_points: *n_points,
                    min_scale_exponent: *min_scale_exponent,
                    max_scale_exponent: *max_scale_exponent,
                    seed: (*seed).unwrap_or(42),
                    use_f128: *use_f128,
                    analyse_analytically: *analyse_analytically,
                    ..Default::default()
                };
                let profile_res = amplitude.profile(&model, &profile_settings)?.analyse();

                for t in profile_res.tables_per_graph(-0.9) {
                    info!("\n{}", t);
                }

                for t in profile_res.analytic_tables_per_graph() {
                    let Some(t) = t else {
                        continue;
                    };
                    info!("\n{}", t);
                }

                if let Some(file) = output_file {
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
            }) => {
                let ir_profile_settings = IRProfileSetting {
                    lambda_exp_start: *min_scale_exponent,
                    lambda_exp_end: *max_scale_exponent,
                    steps: *n_points,
                    seed: seed.unwrap_or(420),
                    select_limits_and_graphs: select.clone(),
                };

                let (process_id, integrand_name) =
                    state.find_integrand_ref(process.as_ref(), integrand_name.as_ref())?;
                let model = state.resolve_model_for_integrand(process_id, &integrand_name)?;
                let process_ref = ProcessRef::Id(process_id);
                let cross_section = state
                    .process_list
                    .get_cross_section_mut_ref(Some(&process_ref), Some(&integrand_name))?;

                let profile_result = match cross_section.integrand.as_mut().ok_or(eyre!(
                    "Integrand {} has not yet been generated",
                    cross_section.name
                ))? {
                    ProcessIntegrand::CrossSection(cross_section_integrand) => {
                        cross_section_integrand.test_ir(&ir_profile_settings, &model)?
                    }
                    _ => unreachable!(),
                };

                info!("\n{}", profile_result);
                Ok(ProfileResult::InfraRed(profile_result))
            }
        }
    }
}
