use std::path::PathBuf;

use ahash::HashMap;
// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;

use crate::GammaLoopContext;
use serde::{Deserialize, Serialize};

use crate::{integrands::Integrand, model::Model, ProcessSettings, Settings};

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ProcessList<A: AmplitudeState = (), C: CrossSectionState = ()> {
    pub processes: Vec<Process<A, C>>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct ExportSettings {
    pub root_folder: PathBuf,
}

impl Default for ProcessList {
    fn default() -> Self {
        Self::new()
    }
}

// Python Api
impl ProcessList {
    /// Generates a new empty process list
    pub(crate) fn new() -> Self {
        ProcessList { processes: vec![] }
    }

    pub fn get_integrand(
        &self,
        process_id: usize,
        integrand_name: impl AsRef<str>,
    ) -> Result<&crate::new_gammaloop_integrand::NewIntegrand> {
        let process = &self.processes[process_id];
        process.get_integrand(integrand_name)
    }

    pub fn get_integrand_mut(
        &mut self,
        process_id: usize,
        integrand_name: impl AsRef<str>,
    ) -> Result<&mut crate::new_gammaloop_integrand::NewIntegrand> {
        let process = &mut self.processes[process_id];
        process.get_integrand_mut(integrand_name)
    }

    pub fn export_dot(&self, settings: &ExportSettings, model: &Model) -> Result<()> {
        for p in &self.processes {
            p.export_dot(settings, model)?;
        }
        Ok(())
    }

    pub fn add_process(&mut self, process: Process) {
        self.processes.push(process);
    }

    /// imports a process list from a folder
    pub(crate) fn import(_settings: ExportSettings) -> Result<Self> {
        Ok(Self::new())
    }

    ///preprocesses the process list according to the settings
    pub fn preprocess(&mut self, model: &Model, settings: ProcessSettings) -> Result<()> {
        for process in self.processes.iter_mut() {
            process.preprocess(model, &settings)?;
        }
        Ok(())
    }

    pub fn generate_integrands(&mut self, settings: &Settings, model: &Model) -> Result<()> {
        for process in &mut self.processes {
            process.generate_integrands(&settings, model)?;
        }
        Ok(())
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

    use linnet::half_edge::involution::EdgeIndex;
    use spenso::network::library::TensorLibraryData;
    use symbolica::{atom::Atom, state::State};

    use crate::{
        dot,
        new_graph::{parse::IntoGraph, Graph, LoopMomentumBasis},
        numerator::{
            EvaluatorOptions, GammaAlgebraMode, GlobalPrefactor, NumeratorEvaluatorOptions,
            NumeratorParseMode, NumeratorSettings, UnInit,
        },
        signature::LoopExtSignature,
        tests_from_pytest::load_generic_model,
        GammaLoopContextContainer, GammaloopCompileOptions, TropicalSubgraphTableSettings,
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
        let mut loop_momentum_basis = LoopMomentumBasis {
            tree: None,
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
                &crate::ProcessSettings {
                    compile_cff: false,
                    compile_separate_orientations: false,
                    cpe_rounds_cff: None,
                    numerator_settings: NumeratorSettings {
                        eval_settings: NumeratorEvaluatorOptions::Single(EvaluatorOptions {
                            cpe_rounds: None,
                            compile_options: crate::numerator::NumeratorCompileOptions::NotCompiled,
                        }),
                        parse_mode: NumeratorParseMode::Direct,
                        dump_expression: None,
                        // global_numerator: None,
                        // global_prefactor: GlobalPrefactor::default(),
                        gamma_algebra: GammaAlgebraMode::Concrete,
                    },
                    gammaloop_compile_options: GammaloopCompileOptions {
                        inline_asm: false,
                        fast_math: false,
                        optimization_level: 0,
                        unsafe_math: false,
                        compiler: "g++".into(),
                        custom: Vec::new(),
                    },
                    tropical_subgraph_table_settings: TropicalSubgraphTableSettings {
                        panic_on_fail: false,
                        target_omega: 1.0,
                    },
                },
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
