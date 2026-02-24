use std::{fs, path::Path};

// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use eyre::Context;
use rayon::ThreadPool;
use schemars::JsonSchema;
use symbolica::evaluate::OptimizationSettings;
use tracing::debug;

use crate::{
    GammaLoopContext, GammaLoopContextContainer,
    settings::{GlobalSettings, runtime::LockedRuntimeSettings},
    utils::serde_utils::{IsDefault, is_false, is_usize},
};
use serde::{Deserialize, Serialize};

use crate::model::Model;

#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
pub struct EvaluatorSettings {
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub iterative_orientation_optimization: bool,
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub summed: bool,
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub summed_function_map: bool,
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub compile: bool,
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub store_atom: bool,
    #[serde(default, skip_serializing_if = "is_usize::<10>")]
    pub horner_iterations: usize,
    #[serde(default, skip_serializing_if = "is_usize::<10>")]
    pub n_cores: usize,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub cpe_iterations: Option<usize>,
    #[serde(default, skip_serializing_if = "is_usize::<10>")]
    pub abort_level: usize,

    #[serde(default, skip_serializing_if = "is_usize::<10>")]
    pub max_horner_scheme_variables: usize,

    #[serde(default, skip_serializing_if = "is_usize::<10>")]
    pub max_common_pair_cache_entries: usize,

    #[serde(default, skip_serializing_if = "is_usize::<10>")]
    pub max_common_pair_distance: usize,
    #[serde(default, skip_serializing_if = "is_false")]
    pub verbose: bool,
}

impl Default for EvaluatorSettings {
    fn default() -> Self {
        Self {
            iterative_orientation_optimization: true,
            summed: false,
            summed_function_map: false,
            compile: false,
            store_atom: false,
            horner_iterations: 10,
            n_cores: 1,
            cpe_iterations: None,
            abort_level: 0,
            max_horner_scheme_variables: 500,
            max_common_pair_cache_entries: 1_000_000,
            max_common_pair_distance: 1000,
            verbose: false,
        }
    }
}

impl EvaluatorSettings {
    pub fn optimization_settings(&self) -> OptimizationSettings {
        OptimizationSettings {
            horner_iterations: self.horner_iterations,
            n_cores: self.n_cores,
            cpe_iterations: self.cpe_iterations,
            hot_start: None,
            abort_check: None,
            abort_level: self.abort_level,
            max_horner_scheme_variables: self.max_horner_scheme_variables,
            max_common_pair_cache_entries: self.max_common_pair_cache_entries,
            max_common_pair_distance: self.max_common_pair_distance,
            verbose: self.verbose,
        }
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ProcessList {
    // pub amplitude_process: Vec<Process<Amplitude>>
    // pub cross_section_process: Vec<Process<CrossSection>>
    pub processes: Vec<Process>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
pub struct DotExportSettings {
    // pub root_folder: PathBuf,
    pub combine_diagrams: bool,
    pub with_uv: bool,
    pub output_full_numerator: bool,
    pub split_xs_by_initial_states: bool,
    pub do_gamma_algebra: bool,
    pub do_color_algebra: bool,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
pub struct StandaloneExportSettings {
    #[serde(default)]
    pub mode: StandaloneExportMode,
}

#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[derive(
    Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq, Default, Encode, Decode, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum StandaloneExportMode {
    #[default]
    Rust,
    Python,
}

#[cfg(feature = "python_api")]
#[cfg_attr(feature = "python_api", pyo3::pymethods)]
impl DotExportSettings {
    #[new]
    fn new() -> Self {
        DotExportSettings::default()
    }
}

impl Default for DotExportSettings {
    fn default() -> Self {
        Self {
            // root_folder: PathBuf::from("."),
            split_xs_by_initial_states: false,
            output_full_numerator: false,
            do_gamma_algebra: false,
            do_color_algebra: true,
            combine_diagrams: false,
            with_uv: false,
        }
    }
}

impl Default for StandaloneExportSettings {
    fn default() -> Self {
        Self {
            mode: StandaloneExportMode::Rust,
        }
    }
}

impl Default for ProcessList {
    fn default() -> Self {
        Self::new()
    }
}

// Python Api
impl ProcessList {
    /// Generates a new empty process list
    pub fn new() -> Self {
        ProcessList { processes: vec![] }
    }

    pub fn warm_up_all_processes(&mut self, model: &Model) -> Result<()> {
        for process in &mut self.processes.iter_mut() {
            process.warm_up(model)?;
        }

        Ok(())
    }

    pub fn load(path: impl AsRef<Path>, context: GammaLoopContextContainer) -> Result<Self> {
        let mut process_list = Self::new();

        let path = path.as_ref().join("processes");
        let amplitudes_path = path.join("amplitudes");
        if amplitudes_path.exists() {
            debug!("Looking for amplitudes in {}", amplitudes_path.display());

            for entry in fs::read_dir(amplitudes_path).context("Error reading dir")? {
                let Ok(entry) = entry else {
                    debug!("Skipping invalid entry");

                    continue;
                };
                let path = entry.path();
                process_list.processes.push(
                    Process::load_amplitude(path, context).context("Error loading amplitude")?,
                );
            }
        }

        let cross_sections_path = path.join("cross_sections");
        if cross_sections_path.exists() {
            debug!(
                "Looking for cross sections in {}",
                cross_sections_path.display()
            );

            for entry in fs::read_dir(cross_sections_path)? {
                let entry = entry?;
                let path = entry.path();
                process_list
                    .processes
                    .push(Process::load_cross_section(path, context)?);
            }
        }
        process_list
            .processes
            .sort_by_key(|p| p.definition.process_id);

        Ok(process_list)
    }

    pub fn save(&mut self, folder: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        let path = folder.as_ref().join("processes");

        let r = fs::create_dir_all(&path);
        if !override_existing {
            r?;
        }

        for p in self.processes.iter_mut() {
            p.save(&path, override_existing)?;
        }

        Ok(())
    }

    pub fn compile(
        &mut self,
        folder: impl AsRef<Path>,
        override_existing: bool,
        settings: &GlobalSettings,
        process_id: Option<usize>,
        integrand_name: Option<String>,
        thread_pool: &ThreadPool,
    ) -> Result<()> {
        let path = folder.as_ref().join("processes");

        let r = fs::create_dir_all(&path);
        if !override_existing {
            r?;
        }

        for p in self.processes.iter_mut() {
            if let Some(id) = process_id
                && p.definition.process_id != id
            {
                continue;
            }
            p.compile(
                &path,
                override_existing,
                settings,
                integrand_name.clone(),
                thread_pool,
            )?;
        }

        Ok(())
    }

    pub fn get_integrand(
        &self,
        process_id: usize,
        integrand_name: impl AsRef<str>,
    ) -> Result<&crate::gammaloop_integrand::GLIntegrand> {
        let process = &self.processes[process_id];
        process.get_integrand(integrand_name)
    }

    pub fn get_integrand_mut(
        &mut self,
        process_id: usize,
        integrand_name: impl AsRef<str>,
    ) -> Result<&mut crate::gammaloop_integrand::GLIntegrand> {
        let process = &mut self.processes[process_id];
        process.get_integrand_mut(integrand_name)
    }

    pub fn export_dot(&self, path: impl AsRef<Path>, settings: &DotExportSettings) -> Result<()> {
        let path = path.as_ref().join("processes");
        fs::create_dir_all(&path)?;

        for p in self.processes.iter() {
            p.export_dot(&path, settings)?;
        }
        Ok(())
    }

    pub fn export_standalone(
        &self,
        path: impl AsRef<Path>,
        settings: &StandaloneExportSettings,
    ) -> Result<()> {
        let path = path.as_ref().join("processes");
        fs::create_dir_all(&path)?;

        for p in self.processes.iter() {
            p.export_standalone(&path, settings)?;
        }
        Ok(())
    }

    pub fn add_process(&mut self, process: Process) {
        self.processes.push(process);
    }

    ///preprocesses the process list according to the settings
    pub fn preprocess(
        &mut self,
        model: &Model,
        settings: &GlobalSettings,
        locked_runtime_settings: &LockedRuntimeSettings,
        thread_pool: &ThreadPool,
    ) -> Result<()> {
        for process in self.processes.iter_mut() {
            process.preprocess(model, settings, locked_runtime_settings, thread_pool)?;
        }
        Ok(())
    }

    pub fn generate_integrands(
        &mut self,
        model: &Model,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
        thread_pool: &ThreadPool,
    ) -> Result<()> {
        for process in &mut self.processes {
            process.generate_integrands(model, global_settings, runtime_default, thread_pool)?;
        }
        Ok(())
    }

    pub fn find_process(&self, process_id: Option<usize>) -> Result<usize> {
        if let Some(id) = process_id {
            if id >= self.processes.len() {
                return Err(color_eyre::eyre::eyre!(
                    "Invalid process id {}. Number of processes: {}",
                    id,
                    self.processes.len()
                ));
            }
            Ok(id)
        } else if self.processes.is_empty() {
            Err(color_eyre::eyre::eyre!("No processes generated yet."))
        } else if self.processes.len() > 1 {
            Err(color_eyre::eyre::eyre!(
                "There are {} processes available. Please specify a process.",
                self.processes.len()
            ))
        } else {
            Ok(0)
        }
    }

    pub fn find_integrand(
        &self,
        process_id: Option<usize>,
        integrand_name: Option<&String>,
    ) -> Result<(usize, String)> {
        let p_id = self.find_process(process_id)?;
        let integrand_name = self.processes[p_id]
            .collection
            .find_integrand(integrand_name.cloned())?;
        Ok((p_id, integrand_name))
    }
}

pub mod process;
pub use process::*;
pub mod amplitude;
pub use amplitude::*;
pub mod cross_section;
pub use cross_section::*;

#[cfg(test)]
mod tests {
    use std::fs::OpenOptions;

    use linnet::half_edge::{
        involution::EdgeIndex,
        subgraph::{SuBitGraph, SubSetLike},
    };

    use symbolica::state::State;

    use crate::{
        GammaLoopContextContainer, dot,
        graph::{Graph, LoopMomentumBasis, parse::IntoGraph},
        settings::{
            RuntimeSettings,
            global::{
                CompilationOptimizationLevel, GammaloopCompileOptions, GenerationSettings,
                TropicalSubgraphTableSettings,
            },
            runtime::LockedRuntimeSettings,
        },
        signature::LoopExtSignature,
        utils::test_utils::load_generic_model,
    };

    use super::AmplitudeGraph;

    #[test]
    fn test_encode_decode_amplitude_graph() {
        // load the model and hack the masses, go through serializable model since arc is not mutable
        let model = load_generic_model("sm");

        let mut graph: Graph = dot!(
            digraph G{
                e1      [flow=sink]
                e2      [flow=source]
                e3      [flow=source]
                e1 -> n1  [particle=h]
                e2 -> n4    [particle=h]
                n1 -> n2    [particle=h]
                n1 -> n3    [particle=h]
                n2 -> n3    [particle=t]
                n3 -> n4    [particle=t]
                n4 -> n2    [particle=t]
            }
        )
        .unwrap();
        let loop_momentum_basis = LoopMomentumBasis {
            tree: SuBitGraph::empty(0),
            loop_edges: vec![EdgeIndex::from(0), EdgeIndex::from(4)].into(),
            ext_edges: vec![EdgeIndex::from(5), EdgeIndex::from(6)].into(),
            edge_signatures: graph
                .underlying
                .new_edgevec(|_, _, _| LoopExtSignature::from((vec![], vec![]))),
        };

        // loop_momentum_basis
        //     .set_edge_signatures(&graph.underlying)
        //     .unwrap();

        graph.loop_momentum_basis = loop_momentum_basis;

        let mut amplitude: AmplitudeGraph = AmplitudeGraph::new(graph.clone());

        amplitude
            .preprocess(
                &model,
                &GenerationSettings {
                    compile: GammaloopCompileOptions {
                        inline_asm: false,
                        fast_math: false,
                        optimization_level: CompilationOptimizationLevel::O0,
                        unsafe_math: false,
                        compiler: "g++".into(),
                        custom: Vec::new(),
                    },
                    tropical_subgraph_table: TropicalSubgraphTableSettings {
                        panic_on_fail: false,
                        target_omega: 1.0,
                        ..Default::default()
                    },
                    ..Default::default()
                },
                &LockedRuntimeSettings::from(&RuntimeSettings::default()),
            )
            .unwrap();

        let mut temp = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open("test.bin")
            .unwrap();

        State::export(&mut temp).unwrap();
        drop(temp);

        let mut temp = OpenOptions::new().read(true).open("test.bin").unwrap();
        let state_map = State::import(&mut temp, None).unwrap();

        let context = GammaLoopContextContainer {
            model: &model,
            state_map: &state_map,
        };

        println!("context created");

        let encoded_amplitude =
            bincode::encode_to_vec(&amplitude, bincode::config::standard()).unwrap();

        let _amplitude: AmplitudeGraph = bincode::decode_from_slice_with_context(
            &encoded_amplitude,
            bincode::config::standard(),
            context,
        )
        .expect("amplitude decode failed")
        .0;

        println!("amplitude graph passed");
    }
}
