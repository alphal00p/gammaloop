use crate::{
    cli_functions::cli,
    cross_section::{Amplitude, AmplitudeList, CrossSection, CrossSectionList},
    inspect,
    integrands::Integrand,
    integrate::{
        havana_integrate, print_integral_result, MasterNode, SerializableBatchResult,
        SerializableIntegrationState,
    },
    model::Model,
    HasIntegrand, Settings,
};
use ahash::HashMap;
use colored::Colorize;
use git_version::git_version;
use log::{info, warn};
use std::{
    fs,
    path::{Path, PathBuf},
};
use symbolica::printer::PrintOptions;

const GIT_VERSION: &str = git_version!();

#[allow(unused)]
use pyo3::{
    exceptions,
    prelude::*,
    pyclass,
    pyclass::CompareOp,
    pyfunction, pymethods, pymodule,
    types::{PyModule, PyTuple, PyType},
    wrap_pyfunction, FromPyObject, IntoPy, PyObject, PyRef, PyResult, Python,
};

#[pyfunction]
#[pyo3(name = "rust_cli")]
fn cli_wrapper(py: Python) -> PyResult<()> {
    /*
    // This one includes python and the name of the wrapper script itself, e.g.
    // `["/home/ferris/.venv/bin/python", "/home/ferris/.venv/bin/print_cli_args", "a", "b", "c"]`
    println!("{:?}", env::args().collect::<Vec<_>>());
    // This one includes only the name of the wrapper script itself, e.g.
    // `["/home/ferris/.venv/bin/print_cli_args", "a", "b", "c"])`
    println!(
        "{:?}",
        py.import("sys")?
            .getattr("argv")?
            .extract::<Vec<String>>()?
    );
    Ok(())
    */
    crate::set_interrupt_handler();
    cli(&py
        .import("sys")?
        .getattr("argv")?
        .extract::<Vec<String>>()?)
    .map_err(|e| exceptions::PyException::new_err(e.to_string()))
}

#[pymodule]
#[pyo3(name = "_gammaloop")]
fn gammalooprs(_py: Python, m: &PyModule) -> PyResult<()> {
    // TODO: Verify that indeed Python logger level is used in that case.
    pyo3_log::init();
    crate::set_interrupt_handler();
    m.add_class::<PythonWorker>()?;
    m.add("git_version", GIT_VERSION)?;
    m.add_wrapped(wrap_pyfunction!(cli_wrapper))?;
    Ok(())
}

pub struct OutputOptions {}

#[pyclass(name = "Worker")]
pub struct PythonWorker {
    pub model: Model,
    pub cross_sections: CrossSectionList,
    pub amplitudes: AmplitudeList,
    pub integrands: HashMap<String, Integrand>,
    pub master_node: Option<MasterNode>,
}

impl Clone for PythonWorker {
    fn clone(&self) -> PythonWorker {
        PythonWorker {
            model: self.model.clone(),
            cross_sections: self.cross_sections.clone(),
            amplitudes: self.amplitudes.clone(),
            integrands: self.integrands.clone(),
            master_node: self.master_node.clone(),
        }
    }
}

// TODO: Improve error broadcasting to Python so as to show rust backtrace
#[pymethods]
impl PythonWorker {
    #[classmethod]
    pub fn new(_cls: &PyType) -> PyResult<PythonWorker> {
        crate::set_interrupt_handler();
        Ok(PythonWorker {
            model: Model::default(),
            cross_sections: CrossSectionList::default(),
            amplitudes: AmplitudeList::default(),
            integrands: HashMap::default(),
            master_node: None,
        })
    }

    pub fn load_model(&mut self, file_path: &str) -> PyResult<()> {
        Model::from_file(String::from(file_path))
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))
            .map(|m| self.model = m)
    }

    pub fn load_model_from_yaml_str(&mut self, yaml_str: &str) -> PyResult<()> {
        Model::from_yaml_str(String::from(yaml_str))
            .map_err(|e| exceptions::PyException::new_err(e.root_cause().to_string()))
            .map(|m| self.model = m)
    }

    // Note: one could consider returning a PyModel class containing the serialisable model as well,
    // but since python already has its native class for this, it is better for now to pass a yaml representation
    // which will be deserialize in said native class.
    pub fn get_model(&self) -> PyResult<String> {
        self.model
            .to_yaml()
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))
    }

    pub fn add_cross_section_from_yaml_str(&mut self, yaml_str: &str) -> PyResult<()> {
        if self.model.is_empty() {
            return Err(exceptions::PyException::new_err(
                "Model must be loaded before cross sections",
            ));
        }
        match CrossSection::from_yaml_str(&self.model, String::from(yaml_str)) {
            Ok(cs) => {
                self.cross_sections.add_cross_section(cs);
                Ok(())
            }
            Err(e) => Err(exceptions::PyException::new_err(e.to_string())),
        }
    }

    pub fn load_cross_sections(&mut self, file_path: &str) -> PyResult<()> {
        if self.model.is_empty() {
            return Err(exceptions::PyException::new_err(
                "Model must be loaded before cross sections",
            ));
        }
        CrossSectionList::from_file(&self.model, String::from(file_path))
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))
            .map(|cs| self.cross_sections = cs)
    }

    pub fn load_cross_sections_from_yaml_str(&mut self, yaml_str: &str) -> PyResult<()> {
        if self.model.is_empty() {
            return Err(exceptions::PyException::new_err(
                "Model must be loaded before cross sections",
            ));
        }
        CrossSectionList::from_yaml_str(&self.model, String::from(yaml_str))
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))
            .map(|cs| self.cross_sections = cs)
    }

    // Note: one could consider returning a PyCrossSectionList class containing the serialisable model as well,
    // but since python already has its native class for this, it is better for now to pass a yaml representation
    // which will be deserialize in said native class.
    pub fn get_cross_sections(&self) -> PyResult<String> {
        self.cross_sections
            .to_yaml()
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))
    }

    pub fn reset_cross_sections(&mut self) -> PyResult<()> {
        self.cross_sections = CrossSectionList::default();
        Ok(())
    }

    pub fn add_amplitude_from_yaml_str(&mut self, yaml_str: &str) -> PyResult<()> {
        if self.model.is_empty() {
            return Err(exceptions::PyException::new_err(
                "Model must be loaded before cross sections",
            ));
        }
        match Amplitude::from_yaml_str(&self.model, String::from(yaml_str)) {
            Ok(amp) => {
                self.amplitudes.add_amplitude(amp);
                Ok(())
            }
            Err(e) => Err(exceptions::PyException::new_err(e.to_string())),
        }
    }

    pub fn load_amplitudes(&mut self, file_path: &str) -> PyResult<()> {
        if self.model.is_empty() {
            return Err(exceptions::PyException::new_err(
                "Model must be loaded before amplitudes",
            ));
        }
        AmplitudeList::from_file(&self.model, String::from(file_path))
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))
            .map(|a| self.amplitudes = a)
    }

    pub fn load_amplitudes_from_yaml_str(&mut self, yaml_str: &str) -> PyResult<()> {
        if self.model.is_empty() {
            return Err(exceptions::PyException::new_err(
                "Model must be loaded before amplitudes",
            ));
        }
        AmplitudeList::from_yaml_str(&self.model, String::from(yaml_str))
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))
            .map(|a| self.amplitudes = a)
    }

    pub fn load_amplitudes_derived_data(&mut self, path: &str) -> PyResult<()> {
        self.amplitudes
            .load_derived_data(path)
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))
    }

    pub fn generate_numerators(&mut self) {
        self.amplitudes.generate_numerator(&self.model);
    }

    // Note: one could consider returning a PyAmpltiudeList class containing the serialisable model as well,
    // but since python already has its native class for this, it is better for now to pass a yaml representation
    // which will be deserialize in said native class.
    pub fn get_amplitudes(&self) -> PyResult<String> {
        self.amplitudes
            .to_yaml()
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))
    }

    pub fn reset_amplitudes(&mut self) -> PyResult<()> {
        self.amplitudes = AmplitudeList::default();
        Ok(())
    }

    pub fn export_cross_sections(
        &mut self,
        export_root: &str,
        cross_section_names: Vec<&str>,
    ) -> PyResult<String> {
        let mut n_exported: usize = 0;
        for cross_section in &self.cross_sections.container {
            if cross_section_names.contains(&cross_section.name.as_str()) {
                n_exported += 1;
                let res = cross_section.export(export_root, &self.model);
                if let Err(err) = res {
                    return Err(exceptions::PyException::new_err(err.to_string()));
                }
            }
        }
        if n_exported != cross_section_names.len() {
            return Err(exceptions::PyException::new_err(format!(
                "Could not find all cross sections to export: {:?}",
                cross_section_names
            )));
        }
        Ok("Successful export".to_string())
    }

    pub fn export_amplitudes(
        &mut self,
        export_root: &str,
        amplitude_names: Vec<&str>,
    ) -> PyResult<String> {
        let mut n_exported: usize = 0;
        for amplitude in self.amplitudes.container.iter_mut() {
            if amplitude_names.contains(&amplitude.name.as_str()) {
                n_exported += 1;
                let res = amplitude.export(export_root, &self.model);
                if let Err(err) = res {
                    return Err(exceptions::PyException::new_err(err.to_string()));
                }
            }
        }
        if n_exported != amplitude_names.len() {
            return Err(exceptions::PyException::new_err(format!(
                "Could not find all amplitudes to export: {:?}",
                amplitude_names
            )));
        }
        Ok("Successful export".to_string())
    }

    pub fn load_amplitude_integrands(&mut self, path_to_settings: &str) -> PyResult<String> {
        self.integrands.clear();

        let mut integrand_counter = 0;
        for amplitude in &self.amplitudes.container {
            let integrand = match amplitude.generate_integrand(Path::new(path_to_settings)) {
                Ok(integrand) => integrand,
                Err(err) => return Err(exceptions::PyException::new_err(err.to_string())),
            };
            self.integrands.insert(
                amplitude.name.to_string(),
                Integrand::GammaLoopIntegrand(integrand),
            );
            integrand_counter += 1;
        }

        log::info!("Loaded integrands {:?}", self.integrands.keys());

        Ok(format!(
            "Loaded {} integrands from {} amplitudes",
            integrand_counter,
            self.amplitudes.container.len()
        ))
    }

    pub fn export_expressions(&mut self, export_root: &str, format: &str) -> PyResult<String> {
        self.generate_numerators();

        for amplitude in self.amplitudes.container.iter_mut() {
            amplitude
                .export_expressions(export_root, Self::printer_options(format))
                .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;
        }
        Ok("Exported expressions".to_string())
    }

    pub fn export_coupling_replacement_rules(
        &self,
        export_root: &str,
        format: &str,
    ) -> PyResult<String> {
        self.model
            .export_coupling_replacement_rules(export_root, Self::printer_options(format))
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;
        Ok("Exported coupling substitutions".to_string())
    }

    pub fn inspect_integrand(
        &mut self,
        integrand: &str,
        pt: Vec<f64>,
        term: Vec<usize>,
        force_radius: bool,
        is_momentum_space: bool,
        use_f128: bool,
    ) -> PyResult<String> {
        match self.integrands.get_mut(integrand) {
            Some(integrand) => {
                let settings = match integrand {
                    Integrand::GammaLoopIntegrand(integrand) => integrand.settings.clone(),
                    _ => todo!(),
                };

                inspect::inspect(
                    &settings,
                    integrand,
                    pt,
                    &term,
                    force_radius,
                    is_momentum_space,
                    use_f128,
                );
            }
            None => {
                return Err(exceptions::PyException::new_err(format!(
                    "Could not find integrand {}",
                    integrand
                )))
            }
        };

        Ok(format!("Inspected integrand: {:?}", integrand))
    }

    pub fn integrate_integrand(
        &mut self,
        integrand: &str,
        num_cores: usize,
        result_path: &str,
        workspace_path: &str,
        target: Option<(f64, f64)>,
    ) -> PyResult<String> {
        match self.integrands.get_mut(integrand) {
            Some(integrand_enum) => match integrand_enum {
                Integrand::GammaLoopIntegrand(gloop_integrand) => {
                    let target = match target {
                        Some((re, im)) => Some(num::Complex::new(re, im)),
                        _ => None,
                    };

                    info!("Gammaloop now integrates {}", integrand.green().bold());

                    let workspace_path = PathBuf::from(workspace_path);

                    let path_to_state = workspace_path.join("integration_state");

                    let integration_state = match fs::read(path_to_state) {
                        Ok(state_bytes) => {
                            info!(
                                "{}",
                                "Found integration state, result of previous integration:".yellow()
                            );
                            info!("");

                            let serializable_state: SerializableIntegrationState =
                                bincode::deserialize::<SerializableIntegrationState>(&state_bytes)
                                    .unwrap();

                            let path_to_workspace_settings = workspace_path.join("settings.yaml");
                            let workspace_settings_string =
                                fs::read_to_string(path_to_workspace_settings)
                                    .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

                            let workspace_settings: Settings =
                                serde_yaml::from_str(&workspace_settings_string)
                                    .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

                            // force the settings to be the same as the ones used in the previous integration
                            gloop_integrand.settings = workspace_settings.clone();

                            let state =
                                serializable_state.into_integration_state(&workspace_settings);

                            print_integral_result(
                                &state.all_integrals[0],
                                1,
                                state.iter,
                                "re",
                                target.map(|c| c.re),
                            );

                            print_integral_result(
                                &state.all_integrals[1],
                                2,
                                state.iter,
                                "im",
                                target.map(|c| c.im),
                            );
                            info!("");
                            warn!("Any changes to the settings will be ignored, integrate with the {} option for changes to take effect","--restart".blue());
                            info!("{}", "Resuming integration".yellow());

                            Some(state)
                        }

                        Err(_) => {
                            info!("No integration state found, starting new integration");
                            None
                        }
                    };

                    let settings = gloop_integrand.settings.clone();
                    let result = havana_integrate(
                        &settings,
                        |set| gloop_integrand.user_data_generator(num_cores, set),
                        target,
                        integration_state,
                        Some(workspace_path),
                    );

                    fs::write(
                        result_path,
                        serde_yaml::to_string(&result)
                            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?,
                    )?;

                    Ok(format!("Integrated integrand {}", integrand))
                }
                _ => unimplemented!("unsupported integrand type"),
            },
            None => Err(exceptions::PyException::new_err(format!(
                "Could not find integrand {}",
                integrand
            ))),
        }
    }

    pub fn load_master_node(&mut self, integrand: &str) -> PyResult<String> {
        let selected_integrand = if let Some(selected_integrand) = self.integrands.get(integrand) {
            selected_integrand
        } else {
            return Err(exceptions::PyException::new_err(format!(
                "could not find integrand {}",
                integrand
            )));
        };

        let grid = selected_integrand.create_grid();
        let integrator_settings = selected_integrand.get_integrator_settings();

        let master_node = MasterNode::new(grid, integrator_settings);
        self.master_node = Some(master_node);

        Ok(format!("Initialized master grid for {}", integrand))
    }

    pub fn write_batch_input(
        &mut self,
        num_cores: usize,
        num_samples: usize,
        export_grid: bool,
        output_accumulator: bool,
        workspace_path: &str,
        job_id: usize,
    ) -> PyResult<String> {
        let master_node = self.master_node.as_mut().unwrap();

        // extract the integrated phase in a hacky way
        match master_node
            .write_batch_input(
                num_cores,
                num_samples,
                export_grid,
                output_accumulator,
                workspace_path,
                job_id,
            )
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))
        {
            Ok(_) => Ok(format!("Wrote batch input for job {}", job_id)),
            Err(e) => Err(e),
        }
    }

    pub fn process_batch_output(
        &mut self,
        workspace_path: &str,
        job_id: usize,
    ) -> PyResult<String> {
        let master_node = self.master_node.as_mut().unwrap();

        let job_out_name = format!("job_{}_out", job_id);
        let job_out_path = Path::new(workspace_path).join(job_out_name);

        let output_file = std::fs::read(job_out_path)?;
        let batch_result: SerializableBatchResult = bincode::deserialize(&output_file).unwrap();

        master_node
            .process_batch_output(batch_result.into_batch_result())
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        Ok(format!("Processed job {}", job_id))
    }

    pub fn write_default_settings(&self, path: &str) -> PyResult<String> {
        let default = Settings::default();
        let default_string = serde_yaml::to_string(&default)
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        let path = Path::new(path).join("cards").join("run_card.yaml");

        fs::write(path, default_string)
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        Ok("Wrote default settings file".to_string())
    }

    pub fn display_master_node_status(&self) {
        if let Some(master_node) = &self.master_node {
            master_node.display_status();
        }
    }

    pub fn update_iter(&mut self) {
        if let Some(master_node) = &mut self.master_node {
            master_node.update_iter();
        }
    }
}

impl PythonWorker {
    fn printer_options(format: &str) -> PrintOptions {
        match format {
            "file" => PrintOptions::file(),
            "mathematica" => PrintOptions {
                terms_on_new_line: false,
                color_top_level_sum: false,
                color_builtin_symbols: false,
                print_finite_field: true,
                symmetric_representation_for_finite_field: false,
                explicit_rational_polynomial: false,
                number_thousands_separator: None,
                multiplication_operator: ' ',
                square_brackets_for_function: true,
                num_exp_as_superscript: false,
                latex: false,
            },
            "latex" => PrintOptions {
                terms_on_new_line: false,
                color_top_level_sum: false,
                color_builtin_symbols: false,
                print_finite_field: true,
                symmetric_representation_for_finite_field: false,
                explicit_rational_polynomial: false,
                number_thousands_separator: None,
                multiplication_operator: ' ',
                square_brackets_for_function: false,
                num_exp_as_superscript: false,
                latex: true,
            },
            _ => PrintOptions {
                terms_on_new_line: false,
                color_top_level_sum: true,
                color_builtin_symbols: true,
                print_finite_field: true,
                symmetric_representation_for_finite_field: false,
                explicit_rational_polynomial: false,
                number_thousands_separator: None,
                multiplication_operator: '*',
                square_brackets_for_function: false,
                num_exp_as_superscript: false,
                latex: false,
            },
        }
    }
}
