use crate::{state::State, CLISettings};
use clap::Subcommand;
use color_eyre::Result;
use gammalooprs::{
    settings::RuntimeSettings,
    uv::{display_results_summary, run_uv_profile, write_results, UVProfileConfig},
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Profile {
    /// Ultraviolet profile analysis
    UltraViolet {
        /// The process id to inspect
        #[arg(short = 'i', long = "process-id", value_name = "ID")]
        process_id: Option<usize>,

        /// The name of the process to inspect
        #[arg(short = 'n', long = "name", value_name = "NAME")]
        integrand_name: Option<String>,

        /// Number of scaling points to sample
        #[arg(long = "n-points", default_value = "50")]
        n_points: usize,

        /// Minimum scaling factor
        #[arg(long = "min-scaling", default_value = "1e-6")]
        min_scaling: f64,

        /// Maximum scaling factor
        #[arg(long = "max-scaling", default_value = "1e6")]
        max_scaling: f64,

        /// Loop momentum basis to use (e.g., "defining", "LU")
        #[arg(long = "lmb", default_value = "defining")]
        lmb: String,

        /// UV indices to analyze (comma-separated list, e.g., "0,1,2")
        #[arg(long = "uv-indices")]
        uv_indices: Option<String>,

        /// Use f128 precision for evaluation
        #[arg(long = "f128")]
        use_f128: bool,

        /// Disable f128 precision (force f64)
        #[arg(long = "no-f128")]
        no_f128: bool,

        /// Generate plots of the profile results
        #[arg(long = "plots")]
        plots: bool,

        /// Scale cuts by the same factor as momenta
        #[arg(long = "scale-cuts")]
        scale_cuts: bool,

        /// Random seed for momentum sampling
        #[arg(long = "seed")]
        seed: Option<u64>,

        /// Target scaling factor for specific analysis
        #[arg(long = "target-scaling")]
        target_scaling: Option<f64>,

        /// Verbose output
        #[arg(short = 'v', long = "verbose")]
        verbose: bool,

        /// Maximum number of evaluations per scaling point
        #[arg(long = "n-max", default_value = "1000")]
        n_max: usize,

        /// Stability threshold for numerical evaluation
        #[arg(long = "stability-threshold", default_value = "1e-12")]
        stability_threshold: f64,

        /// Output file for results (optional)
        #[arg(short = 'o', long = "output")]
        output_file: Option<String>,

        /// Format for output results
        #[arg(long = "format", default_value = "yaml")]
        output_format: String,
    },
}

impl Profile {
    pub fn run(
        &self,
        state: &mut State,
        _global_cli_settings: &CLISettings,
        runtime_settings: &RuntimeSettings,
    ) -> Result<()> {
        match self {
            Profile::UltraViolet {
                process_id,
                integrand_name,
                n_points,
                min_scaling,
                max_scaling,
                lmb,
                uv_indices,
                use_f128,
                no_f128,
                plots,
                scale_cuts,
                seed,
                target_scaling,
                verbose,
                n_max,
                stability_threshold,
                output_file,
                output_format,
            } => {
                // Parse UV indices if provided
                let parsed_uv_indices = if let Some(indices_str) = uv_indices {
                    let indices: Result<Vec<usize>, _> = indices_str
                        .split(',')
                        .map(|s| s.trim().parse::<usize>())
                        .collect();
                    Some(indices.map_err(|e| {
                        color_eyre::eyre::eyre!(
                            "Failed to parse UV indices '{}': {}",
                            indices_str,
                            e
                        )
                    })?)
                } else {
                    None
                };

                // Create configuration
                let config = UVProfileConfig {
                    n_points: *n_points,
                    min_scaling: *min_scaling,
                    max_scaling: *max_scaling,
                    lmb: lmb.clone(),
                    uv_indices: parsed_uv_indices,
                    use_f128: *use_f128,
                    no_f128: *no_f128,
                    plots: *plots,
                    scale_cuts: *scale_cuts,
                    seed: *seed,
                    target_scaling: *target_scaling,
                    verbose: *verbose,
                    n_max: *n_max,
                    stability_threshold: *stability_threshold,
                    output_file: output_file.clone(),
                    output_format: output_format.clone(),
                };

                // Determine which process/integrand to use
                let (proc_id, integrand_name) = match (process_id, integrand_name) {
                    (Some(id), Some(name)) => (*id, name.clone()),
                    (Some(id), None) => {
                        // Use first available integrand in the process
                        let process = state.process_list.processes.get(*id).ok_or_else(|| {
                            color_eyre::eyre::eyre!("Process ID {} not found", id)
                        })?;

                        let first_name = match &process.collection {
                            gammalooprs::processes::ProcessCollection::Amplitudes(amplitudes) => {
                                amplitudes
                                    .keys()
                                    .next()
                                    .ok_or_else(|| {
                                        color_eyre::eyre::eyre!(
                                            "No amplitudes found in process {}",
                                            id
                                        )
                                    })?
                                    .clone()
                            }
                            gammalooprs::processes::ProcessCollection::CrossSections(
                                cross_sections,
                            ) => cross_sections
                                .keys()
                                .next()
                                .ok_or_else(|| {
                                    color_eyre::eyre::eyre!(
                                        "No cross sections found in process {}",
                                        id
                                    )
                                })?
                                .clone(),
                        };
                        (*id, first_name)
                    }
                    (None, Some(name)) => {
                        // Find the process containing this integrand
                        let mut found_id = None;
                        for (id, process) in state.process_list.processes.iter().enumerate() {
                            let contains_integrand = match &process.collection {
                                gammalooprs::processes::ProcessCollection::Amplitudes(
                                    amplitudes,
                                ) => amplitudes.contains_key(name),
                                gammalooprs::processes::ProcessCollection::CrossSections(
                                    cross_sections,
                                ) => cross_sections.contains_key(name),
                            };
                            if contains_integrand {
                                found_id = Some(id);
                                break;
                            }
                        }
                        let id = found_id.ok_or_else(|| {
                            color_eyre::eyre::eyre!("Integrand '{}' not found in any process", name)
                        })?;
                        (id, name.clone())
                    }
                    (None, None) => {
                        return Err(color_eyre::eyre::eyre!(
                            "Must specify either process-id, integrand name, or both"
                        ));
                    }
                };

                // Follow inspect pattern for integrand access
                state.process_list.processes[proc_id].warm_up(&state.model)?;

                let integrand = state
                    .process_list
                    .get_integrand_mut(proc_id, integrand_name.clone())?;

                println!(
                    "Running UV profile analysis on integrand '{}'",
                    integrand_name
                );
                if config.verbose {
                    println!("Configuration: {:?}", config);
                }

                // Run the UV profile analysis
                let result = run_uv_profile(config, integrand, &state.model, runtime_settings)?;

                // Display results
                display_results_summary(&result);

                // Write results to file if requested
                write_results(&result)?;

                Ok(())
            }
        }
    }
}
