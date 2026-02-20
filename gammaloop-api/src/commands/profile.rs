use std::path::PathBuf;

use crate::{
    state::{ProcessListExt, ProcessRef, State},
    CLISettings,
};
use color_eyre::Result;
use eyre::eyre;
use gammalooprs::uv::{
    profile::{ProfileSettings, UVProfileable},
    UVProfileAnalysis,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use tracing::{info, instrument};

use clap::{Args, Subcommand};

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Profile {
    /// Ultraviolet profile analysis
    UltraViolet(#[command(flatten)] UltraVioletProfile),
}

#[derive(Args, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct UltraVioletProfile {
    /// Process reference: #<id>, name:<name>, or <id>/<name>
    #[arg(short = 'p', long = "process", value_name = "PROCESS")]
    pub process: Option<ProcessRef>,

    /// The name of the process to inspect
    #[arg(short = 'n', long = "name", value_name = "NAME")]
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
    #[arg(short = 'o', long = "output")]
    pub output_file: Option<PathBuf>,
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

impl Profile {
    #[instrument(skip_all)]
    pub fn run(
        &self,
        state: &mut State,
        _global_cli_settings: &CLISettings,
    ) -> Result<UVProfileAnalysis> {
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
                let amplitude = state
                    .process_list
                    .get_amplitude_mut_ref(process.as_ref(), integrand_name.as_ref())?;

                amplitude
                    .integrand
                    .as_mut()
                    .ok_or(eyre!(
                        "Integrand {} has not yet been generated, but exists",
                        amplitude.name
                    ))?
                    .warm_up(&state.model)?;

                let profile_settings = ProfileSettings {
                    n_points: *n_points,
                    min_scale_exponent: *min_scale_exponent,
                    max_scale_exponent: *max_scale_exponent,
                    seed: (*seed).unwrap_or(42),
                    use_f128: *use_f128,
                    analyse_analytically: *analyse_analytically,
                    ..Default::default()
                };
                let profile_res = amplitude
                    .profile(&state.model, &profile_settings)?
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

                if let Some(file) = output_file {
                    profile_res.write_profile_data(file)?
                }

                Ok(profile_res)
            }
        }
    }
}
