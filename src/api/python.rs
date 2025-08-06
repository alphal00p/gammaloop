use crate::{
    cli::Cli,
    cross_section::{Amplitude, AmplitudeList, CrossSection, CrossSectionList},
    feyngen::{
        self, diagram_generator::FeynGen, FeynGenError, FeynGenFilters, FeynGenOptions,
        NumeratorAwareGraphGroupingOption, SewedFilterOptions,
    },
    graph::SerializableGraph,
    inspect,
    integrands::Integrand,
    integrate::{
        havana_integrate, print_integral_result, BatchResult, IntegrationState, MasterNode,
    },
    model::Model,
    new_cs::{self, ExportSettings, Process, ProcessCollection, ProcessDefinition, ProcessList},
    new_graph::Graph,
    numerator::{GlobalPrefactor, Numerator, PythonState},
    utils::F,
    GammaLoopContextContainer, HasIntegrand, IntegratedPhase, OutputMetadata, ProcessSettings,
    Settings,
};
use ahash::HashMap;
use chrono::{Datelike, Local, Timelike};
use clap::Parser;
use colored::{ColoredString, Colorize};
use feyngen::{
    FeynGenFilter, GenerationType, SelfEnergyFilterOptions, SnailFilterOptions,
    TadpolesFilterOptions,
};
use git_version::git_version;
use itertools::{self, Itertools};
use log::{debug, info, warn, LevelFilter};
use pyo3::{import_exception, types::PyDict};
use spenso::algebra::complex::Complex;
use std::cmp::Ordering;
use std::{
    fs,
    path::{Path, PathBuf},
    str::FromStr,
    sync::{LazyLock, Mutex},
};

use symbolica::{
    atom::{Atom, AtomCore},
    parse,
    printer::PrintOptions,
    state::State,
};
const GIT_VERSION: &str = git_version!(fallback = "unavailable");

#[allow(unused)]
use pyo3::{
    exceptions,
    prelude::*,
    pyclass,
    pyclass::CompareOp,
    pyfunction, pymethods, pymodule,
    types::{PyComplex, PyModule, PyTuple, PyType},
    wrap_pyfunction, FromPyObject, IntoPy, PyObject, PyRef, PyResult, Python,
};

//use pyo3_log;

pub enum LogFormat {
    Long,
    Short,
    Min,
    None,
}
impl FromStr for LogFormat {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "long" => Ok(LogFormat::Long),
            "short" => Ok(LogFormat::Short),
            "min" => Ok(LogFormat::Min),
            "none" => Ok(LogFormat::None),
            _ => Err(format!("Unknown log format: {}", s)),
        }
    }
}

pub static LOG_LEVEL: LazyLock<Mutex<Option<LevelFilter>>> = LazyLock::new(|| Mutex::new(None));
pub static LOG_FORMAT: LazyLock<Mutex<LogFormat>> = LazyLock::new(|| Mutex::new(LogFormat::Long));

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

    match Cli::parse_from(
        &py.import_bound("sys")?
            .getattr("argv")?
            .extract::<Vec<String>>()?,
    )
    .run()
    {
        Err(e) => Err(exceptions::PyException::new_err(e.to_string())),
        Ok(_) => Ok(()),
    }
}

pub fn format_level(level: log::Level) -> ColoredString {
    match level {
        log::Level::Error => format!("{:<8}", "ERROR").red(),
        log::Level::Warn => format!("{:<8}", "WARNING").yellow(),
        log::Level::Info => format!("{:<8}", "INFO").into(),
        log::Level::Debug => format!("{:<8}", "DEBUG").bright_black(),
        log::Level::Trace => format!("{:<8}", "TRACE").into(),
    }
}

pub fn format_target(target: String, level: log::Level) -> ColoredString {
    let split_targets = target.split("::").collect::<Vec<_>>();
    //[-2..].iter().join("::");
    let start = split_targets.len().saturating_sub(2);
    let mut shortened_path = split_targets[start..].join("::");
    if level < log::Level::Debug && shortened_path.len() > 20 {
        shortened_path = format!("{}...", shortened_path.chars().take(17).collect::<String>());
    }
    format!("{:<20}", shortened_path).bright_blue()
}

#[pyfunction]
#[pyo3(name = "evaluate_graph_overall_factor")]
pub fn evaluate_graph_overall_factor(overall_factor: &str) -> PyResult<String> {
    let overall_factor = parse!(overall_factor);
    let overall_factor_evaluated = FeynGen::evaluate_overall_factor(overall_factor.as_view());
    Ok(overall_factor_evaluated.to_canonical_string())
}

#[pyfunction]
#[pyo3(name = "atom_to_canonical_string")]
pub fn atom_to_canonical_string(atom_str: &str) -> PyResult<String> {
    Ok(parse!(atom_str).to_canonical_string())
}

#[pyfunction]
#[pyo3(name = "setup_rust_logging")]
fn setup_logging(level: String, format: String) -> PyResult<()> {
    let log_format_guard = &mut *LOG_FORMAT.lock().unwrap();
    *log_format_guard = LogFormat::from_str(&format).map_err(exceptions::PyException::new_err)?;

    let level = convert_log_level(&level);
    if (*LOG_LEVEL.lock().unwrap()).is_none() {
        // Configure logger at runtime
        fern::Dispatch::new()
            .filter(|metadata| metadata.level() <= (*LOG_LEVEL.lock().unwrap()).unwrap())
            // Perform allocation-free log formatting
            .format(|out, message, record| {
                let now = Local::now();
                match *LOG_FORMAT.lock().unwrap() {
                    LogFormat::Long => out.finish(format_args!(
                        "[{}] @{} {}: {}",
                        format!(
                            "{:04}-{:02}-{:02} {:02}:{:02}:{:02}.{:03}",
                            now.year(),
                            now.month(),
                            now.day(),
                            now.hour(),
                            now.minute(),
                            now.second(),
                            now.timestamp_subsec_millis()
                        )
                        .bright_green(),
                        format_target(record.target().into(), record.level()),
                        format_level(record.level()),
                        message
                    )),
                    LogFormat::Short => out.finish(format_args!(
                        "[{}] {}: {}",
                        format!("{:02}:{:02}:{:02}", now.hour(), now.minute(), now.second())
                            .bright_green(),
                        format_level(record.level()),
                        message
                    )),
                    LogFormat::Min => out.finish(format_args!(
                        "{}: {}",
                        format_level(record.level()),
                        message
                    )),
                    LogFormat::None => out.finish(format_args!("{}", message)),
                }
            })
            // Output to stdout, files, and other Dispatch configurations
            .chain(std::io::stdout())
            .chain(fern::log_file("gammaloop_rust_output.log")?)
            // Apply globally
            .apply()
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;
    }
    let log_level_guard = &mut *LOG_LEVEL.lock().unwrap();
    *log_level_guard = Some(level);

    Ok(())
}

pub fn get_python_log_level() -> Result<String, pyo3::PyErr> {
    Python::with_gil(|py| {
        let py_code = r#"
def get_gamma_loop_log_level():
    import logging
    return logging.getLevelName(logging.getLogger('GammaLoop').level)
"#;
        let locals = PyDict::new_bound(py);
        py.run_bound(py_code, None, Some(&locals))?;
        let log_level: String = locals
            .get_item("get_gamma_loop_log_level")?
            .unwrap()
            .call0()?
            .extract()?;
        Ok(log_level)
    })
}

pub fn convert_log_level(level: &str) -> LevelFilter {
    match level {
        "DEBUG" => LevelFilter::Debug,
        "INFO" => LevelFilter::Info,
        "WARNING" => LevelFilter::Warn,
        "ERROR" => LevelFilter::Error,
        "CRITICAL" => LevelFilter::Error, // No direct mapping for CRITICAL, map to ERROR
        "TRACE" => LevelFilter::Trace,    // Python does not have TRACE, but we allow it
        _ => LevelFilter::Trace,          // Default to TRACE if unknown
    }
}

#[pymodule]
#[pyo3(name = "_gammaloop")]
fn gammalooprs(_py: Python, m: &Bound<PyModule>) -> PyResult<()> {
    pyo3_pylogger::register("GammaLoopRust");

    crate::initialisation::initialise().expect("initialization failed");
    crate::set_interrupt_handler();
    m.add_class::<PythonWorker>()?;
    m.add_class::<PyFeynGenFilters>()?;
    m.add_class::<PySnailFilterOptions>()?;
    m.add_class::<PySewedFilterOptions>()?;
    m.add_class::<PyTadpolesFilterOptions>()?;
    m.add_class::<PySelfEnergyFilterOptions>()?;
    m.add_class::<PyFeynGenOptions>()?;
    m.add_class::<PyNumeratorAwareGroupingOption>()?;
    m.add("git_version", GIT_VERSION)?;
    m.add_wrapped(wrap_pyfunction!(cli_wrapper))?;
    m.add_wrapped(wrap_pyfunction!(setup_logging))?;
    m.add_wrapped(wrap_pyfunction!(atom_to_canonical_string))?;
    m.add_wrapped(wrap_pyfunction!(evaluate_graph_overall_factor))?;
    Ok(())
}

pub struct OutputOptions {}

/// There should only be 1 worker per instance of gammaloop.
#[pyclass(name = "Worker", unsendable)]
pub struct PythonWorker {
    pub model: Model,
    pub process_list: ProcessList,
    pub integrands: HashMap<String, Integrand>,
    pub master_node: Option<MasterNode>,
}

impl Clone for PythonWorker {
    fn clone(&self) -> PythonWorker {
        PythonWorker {
            model: self.model.clone(),
            process_list: self.process_list.clone(),
            integrands: self.integrands.clone(),
            master_node: self.master_node.clone(),
        }
    }
}

#[pyclass(name = "SnailFilterOptions")]
pub struct PySnailFilterOptions {
    pub filter_options: SnailFilterOptions,
}

#[pymethods]
impl PySnailFilterOptions {
    #[new]
    pub fn __new__(
        veto_snails_attached_to_massive_lines: Option<bool>,
        veto_snails_attached_to_massless_lines: Option<bool>,
        veto_only_scaleless_snails: Option<bool>,
    ) -> PyResult<PySnailFilterOptions> {
        Ok(PySnailFilterOptions {
            filter_options: SnailFilterOptions {
                veto_snails_attached_to_massive_lines: veto_snails_attached_to_massive_lines
                    .unwrap_or(false),
                veto_snails_attached_to_massless_lines: veto_snails_attached_to_massless_lines
                    .unwrap_or(true),
                veto_only_scaleless_snails: veto_only_scaleless_snails.unwrap_or(false),
            },
        })
    }

    pub fn __str__(&self) -> PyResult<String> {
        Ok(format!("{}", self.filter_options))
    }
}

#[pyclass(name = "SelfEnergyFilterOptions")]
pub struct PySelfEnergyFilterOptions {
    pub filter_options: SelfEnergyFilterOptions,
}

#[pymethods]
impl PySelfEnergyFilterOptions {
    #[new]
    pub fn __new__(
        veto_self_energy_of_massive_lines: Option<bool>,
        veto_self_energy_of_massless_lines: Option<bool>,
        veto_only_scaleless_self_energy: Option<bool>,
    ) -> PyResult<PySelfEnergyFilterOptions> {
        Ok(PySelfEnergyFilterOptions {
            filter_options: SelfEnergyFilterOptions {
                veto_self_energy_of_massive_lines: veto_self_energy_of_massive_lines
                    .unwrap_or(true),
                veto_self_energy_of_massless_lines: veto_self_energy_of_massless_lines
                    .unwrap_or(true),
                veto_only_scaleless_self_energy: veto_only_scaleless_self_energy.unwrap_or(false),
            },
        })
    }

    pub fn __str__(&self) -> PyResult<String> {
        Ok(format!("{}", self.filter_options))
    }
}

#[pyclass(name = "SewedFilterOptions")]
pub struct PySewedFilterOptions {
    pub filter_options: SewedFilterOptions,
}

#[pymethods]
impl PySewedFilterOptions {
    #[new]
    pub fn __new__(filter_tadpoles: Option<bool>) -> PyResult<PySewedFilterOptions> {
        Ok(PySewedFilterOptions {
            filter_options: SewedFilterOptions {
                filter_tadpoles: filter_tadpoles.unwrap_or(false),
            },
        })
    }
}

#[pyclass(name = "TadpolesFilterOptions")]
pub struct PyTadpolesFilterOptions {
    pub filter_options: TadpolesFilterOptions,
}

#[pymethods]
impl PyTadpolesFilterOptions {
    #[new]
    pub fn __new__(
        veto_tadpoles_attached_to_massive_lines: Option<bool>,
        veto_tadpoles_attached_to_massless_lines: Option<bool>,
        veto_only_scaleless_tadpoles: Option<bool>,
    ) -> PyResult<PyTadpolesFilterOptions> {
        Ok(PyTadpolesFilterOptions {
            filter_options: TadpolesFilterOptions {
                veto_tadpoles_attached_to_massive_lines: veto_tadpoles_attached_to_massive_lines
                    .unwrap_or(true),
                veto_tadpoles_attached_to_massless_lines: veto_tadpoles_attached_to_massless_lines
                    .unwrap_or(true),
                veto_only_scaleless_tadpoles: veto_only_scaleless_tadpoles.unwrap_or(false),
            },
        })
    }

    pub fn __str__(&self) -> PyResult<String> {
        Ok(format!("{}", self.filter_options))
    }
}

#[pyclass(name = "FeynGenFilters")]
pub struct PyFeynGenFilters {
    pub filters: Vec<FeynGenFilter>,
}
impl<'a> FromPyObject<'a> for PyFeynGenFilters {
    fn extract_bound(ob: &Bound<'a, PyAny>) -> PyResult<Self> {
        if let Ok(a) = ob.extract::<PyFeynGenFilters>() {
            Ok(PyFeynGenFilters { filters: a.filters })
        } else {
            Err(exceptions::PyValueError::new_err(
                "Not a valid Feynman generation filter",
            ))
        }
    }
}

#[pymethods]
impl PyFeynGenFilters {
    pub fn __str__(&self) -> PyResult<String> {
        Ok(self.filters.iter().map(|f| format!(" > {}", f)).join("\n"))
    }

    #[new]
    #[allow(clippy::too_many_arguments)]
    pub fn __new__(
        particle_veto: Option<Vec<i64>>,
        max_number_of_bridges: Option<usize>,
        sewed_filter: Option<PyRef<PySewedFilterOptions>>,
        self_energy_filter: Option<PyRef<PySelfEnergyFilterOptions>>,
        tadpoles_filter: Option<PyRef<PyTadpolesFilterOptions>>,
        zero_snails_filter: Option<PyRef<PySnailFilterOptions>>,
        perturbative_orders: Option<HashMap<String, usize>>,
        coupling_orders: Option<HashMap<String, (usize, Option<usize>)>>,
        loop_count_range: Option<(usize, usize)>,
        fermion_loop_count_range: Option<(usize, usize)>,
        factorized_loop_topologies_count_range: Option<(usize, usize)>,
    ) -> PyResult<PyFeynGenFilters> {
        let mut filters = Vec::new();
        if let Some(sewed_filter) = sewed_filter {
            filters.push(FeynGenFilter::SewedFilter(sewed_filter.filter_options));
        }
        if let Some(self_energy_filter) = self_energy_filter {
            filters.push(FeynGenFilter::SelfEnergyFilter(
                self_energy_filter.filter_options,
            ));
        }
        if let Some(particle_veto) = particle_veto {
            filters.push(FeynGenFilter::ParticleVeto(particle_veto));
        }
        if let Some(max_number_of_bridges) = max_number_of_bridges {
            filters.push(FeynGenFilter::MaxNumberOfBridges(max_number_of_bridges));
        }
        if let Some(tadpoles_filter) = tadpoles_filter {
            filters.push(FeynGenFilter::TadpolesFilter(
                tadpoles_filter.filter_options,
            ));
        }
        if let Some(zero_snails_filter) = zero_snails_filter {
            filters.push(FeynGenFilter::ZeroSnailsFilter(
                zero_snails_filter.filter_options,
            ));
        }
        if let Some(perturbative_orders) = perturbative_orders {
            filters.push(FeynGenFilter::PerturbativeOrders(perturbative_orders));
        }
        if let Some(coupling_orders) = coupling_orders {
            filters.push(FeynGenFilter::CouplingOrders(coupling_orders));
        }
        if let Some(loop_count_range) = loop_count_range {
            filters.push(FeynGenFilter::LoopCountRange(loop_count_range));
        }
        if let Some(fermion_loop_count_range) = fermion_loop_count_range {
            filters.push(FeynGenFilter::FermionLoopCountRange(
                fermion_loop_count_range,
            ));
        }
        if let Some(factorized_loop_topologies_count_range) = factorized_loop_topologies_count_range
        {
            filters.push(FeynGenFilter::FactorizedLoopTopologiesCountRange(
                factorized_loop_topologies_count_range,
            ));
        }
        Ok(PyFeynGenFilters { filters })
    }
}

#[pyclass(name = "FeynGenOptions")]
pub struct PyFeynGenOptions {
    pub options: FeynGenOptions,
}
impl<'a> FromPyObject<'a> for PyFeynGenOptions {
    fn extract_bound(ob: &Bound<'a, PyAny>) -> PyResult<Self> {
        if let Ok(a) = ob.extract::<PyFeynGenOptions>() {
            Ok(PyFeynGenOptions { options: a.options })
        } else {
            Err(exceptions::PyValueError::new_err(
                "Not a valid Feynman generation option structure",
            ))
        }
    }
}

fn feyngen_to_python_error(error: FeynGenError) -> PyErr {
    exceptions::PyValueError::new_err(format!("Feynam diagram generator error | {error}"))
}

#[pymethods]
impl PyFeynGenOptions {
    pub fn __str__(&self) -> PyResult<String> {
        Ok(format!("{}", self.options))
    }
    #[allow(clippy::too_many_arguments)]
    #[new]
    pub fn __new__(
        generation_type: String,
        initial_particles: Vec<i64>,
        final_particles_lists: Vec<Vec<i64>>,
        loop_count_range: (usize, usize),
        n_cut_blobs: (usize, usize),
        n_cut_spectators: (usize, usize),
        symmetrize_initial_states: bool,
        symmetrize_final_states: bool,
        symmetrize_left_right_states: bool,
        allow_symmetrization_of_external_fermions_in_amplitudes: bool,
        max_multiplicity_for_fast_cut_filter: usize,
        amplitude_filters: Option<PyRef<PyFeynGenFilters>>,
        cross_section_filters: Option<PyRef<PyFeynGenFilters>>,
    ) -> PyResult<PyFeynGenOptions> {
        let amplitude_filters = FeynGenFilters(
            amplitude_filters
                .map(|f| f.filters.clone())
                .unwrap_or_default(),
        );

        let mut cross_section_filters = FeynGenFilters(
            cross_section_filters
                .map(|f| f.filters.clone())
                .unwrap_or_default(),
        );

        cross_section_filters
            .0
            .push(FeynGenFilter::BlobRange(n_cut_blobs.0..=n_cut_blobs.1));
        cross_section_filters.0.push(FeynGenFilter::SpectatorRange(
            n_cut_spectators.0..=n_cut_spectators.1,
        ));

        let feyngen_options = FeynGenOptions {
            generation_type: GenerationType::from_str(&generation_type)
                .map_err(feyngen_to_python_error)?,
            initial_pdgs: initial_particles,
            final_pdgs_lists: final_particles_lists,
            loop_count_range,
            symmetrize_initial_states,
            symmetrize_final_states,
            symmetrize_left_right_states,
            allow_symmetrization_of_external_fermions_in_amplitudes,
            max_multiplicity_for_fast_cut_filter,
            amplitude_filters,
            cross_section_filters,
        };
        if feyngen_options.generation_type == GenerationType::Amplitude {
            if feyngen_options.final_pdgs_lists.len() > 1 {
                return Err(exceptions::PyValueError::new_err(
                    "Multiple set of final states are not allowed for amplitude generation",
                ));
            }
        } else if feyngen_options.final_pdgs_lists.len() > 1
            && feyngen_options
                .final_pdgs_lists
                .iter()
                .any(|l| l.is_empty())
        {
            return Err(exceptions::PyValueError::new_err(
                    "When specifying multiple set of final states options, each must contain at least one particle",
                ));
        }
        Ok(PyFeynGenOptions {
            options: feyngen_options,
        })
    }
}

#[pyclass(name = "NumeratorAwareGroupingOption")]
pub struct PyNumeratorAwareGroupingOption {
    pub grouping_options: NumeratorAwareGraphGroupingOption,
}

#[pymethods]
impl PyNumeratorAwareGroupingOption {
    #[new]
    pub fn __new__(
        numerator_aware_grouping_option: Option<String>,
        compare_canonized_numerator: Option<bool>,
        number_of_samples_for_numerator_comparisons: Option<usize>,
        consider_internal_masses_only_in_numerator_isomorphisms: Option<bool>,
        fully_numerical_substitution_when_comparing_numerators: Option<bool>,
        numerical_samples_seed: Option<u16>,
    ) -> PyResult<PyNumeratorAwareGroupingOption> {
        Ok(PyNumeratorAwareGroupingOption {
            grouping_options: NumeratorAwareGraphGroupingOption::new_with_attributes(
                numerator_aware_grouping_option
                    .unwrap_or("group_identical_graphs_up_to_scalar_rescaling".into())
                    .as_str(),
                numerical_samples_seed,
                number_of_samples_for_numerator_comparisons,
                consider_internal_masses_only_in_numerator_isomorphisms,
                fully_numerical_substitution_when_comparing_numerators,
                compare_canonized_numerator,
            )
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?,
        })
    }

    pub fn __str__(&self) -> PyResult<String> {
        Ok(format!("{}", self.grouping_options))
    }
}

// TODO: Improve error broadcasting to Python so as to show rust backtrace
#[pymethods]
impl PythonWorker {
    #[new]
    pub fn new() -> PyResult<PythonWorker> {
        crate::set_interrupt_handler();
        Ok(PythonWorker {
            model: Model::default(),
            process_list: ProcessList::new(),
            integrands: HashMap::default(),
            master_node: None,
        })
    }

    #[pyo3(signature = (path,name=None))]
    pub fn import_amplitude(&mut self, path: PathBuf, name: Option<String>) -> PyResult<()> {
        // println!("")
        let graphs = Graph::from_file(&path, &self.model)
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;
        let name = name.unwrap_or_else(|| path.file_name().unwrap().to_string_lossy().into_owned());
        let process = Process::from_graph_list(
            name,
            graphs,
            GenerationType::Amplitude,
            ProcessDefinition::new_empty(),
            None,
        )
        .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        self.process_list.add_process(process);
        Ok(())
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

    pub fn preprocess(&mut self, export_yaml_str: &str) -> PyResult<()> {
        let process_settings: ProcessSettings = serde_yaml::from_str(export_yaml_str)
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))
            .unwrap();

        if self.process_list.processes.is_empty() {
            warn!("No processes to preprocess");
        }

        self.process_list
            .preprocess(&self.model, process_settings)
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))
    }

    #[allow(clippy::too_many_arguments)]
    pub fn generate_diagrams(
        &mut self,
        generation_options: PyRef<PyFeynGenOptions>,
        proccess_name: String,
        numerator_aware_isomorphism_grouping: Option<PyRef<PyNumeratorAwareGroupingOption>>,
        filter_self_loop: Option<bool>,
        graph_prefix: Option<String>,
        selected_graphs: Option<Vec<String>>,
        vetoed_graphs: Option<Vec<String>>,
        loop_momentum_bases: Option<HashMap<String, Vec<String>>>,
        global_prefactor_color: Option<String>,
        global_prefactor_colorless: Option<String>,
        num_threads: Option<usize>,
    ) -> PyResult<Vec<String>> {
        if self.model.is_empty() {
            return Err(exceptions::PyException::new_err(
                "A physics model must be loaded before generating diagrams",
            ));
        }
        let feyngen_options = generation_options.options.clone();

        // clone some of the options that will be used in the process definition
        let initial_pdgs = feyngen_options.initial_pdgs.clone();
        let final_pdgs_lists = feyngen_options.final_pdgs_lists.clone();
        let generation_type = feyngen_options.generation_type.clone();
        let amplitude_filters = feyngen_options.amplitude_filters.clone();
        let cross_section_filters = feyngen_options.cross_section_filters.clone();

        let diagram_generator = FeynGen::new(feyngen_options);

        let mut global_prefactor = GlobalPrefactor::default();
        if let Some(global_prefactor_color) = global_prefactor_color {
            global_prefactor.num = parse!(&global_prefactor_color);
        }
        if let Some(global_prefactor_colorless) = global_prefactor_colorless {
            global_prefactor.projector = parse!(&global_prefactor_colorless);
        }

        let diagrams = diagram_generator
            .generate(
                &self.model,
                &numerator_aware_isomorphism_grouping
                    .map(|o| o.grouping_options.clone())
                    .unwrap_or(NumeratorAwareGraphGroupingOption::NoGrouping),
                filter_self_loop.unwrap_or(false),
                graph_prefix.unwrap_or("GL".to_string()),
                selected_graphs,
                vetoed_graphs,
                loop_momentum_bases,
                global_prefactor,
                num_threads,
            )
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        let res = Ok(diagrams
            .iter()
            .map(|d| serde_yaml::to_string(&SerializableGraph::from_graph(d)).unwrap())
            .collect());

        // load everything into processlist
        let (n_unresolved, unresolved_cut_content) =
            diagram_generator.unresolved_cut_content(&self.model);

        let process_definition = ProcessDefinition {
            initial_pdgs,
            final_pdgs_lists,
            n_unresolved,
            unresolved_cut_content: unresolved_cut_content.into_iter().collect(),
            amplitude_filters,
            cross_section_filters,
        };

        let process = Process::from_bare_graph_list(
            proccess_name.into(),
            diagrams,
            generation_type,
            process_definition,
            None,
        )
        .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        self.process_list.add_process(process);

        res
    }

    pub fn export_expressions(
        &mut self,
        export_root: &str,
        amplitued_list: Vec<String>,
        format: &str,
        export_yaml_str: &str,
    ) -> PyResult<String> {
        let export_settings: ProcessSettings = serde_yaml::from_str(export_yaml_str)
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))
            .unwrap();

        for amplitude in amplitued_list.into_iter() {
            match Amplitude::from_yaml_str(&self.model, amplitude) {
                Ok(amp) => {
                    amp.map(|a| a.map(|ag| ag.forget_type()))
                        .export_expressions(
                            export_root,
                            Self::printer_options(format),
                            &export_settings,
                        )
                        .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;
                }
                Err(e) => return Err(exceptions::PyException::new_err(e.to_string())),
            }
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
    ) -> PyResult<(f64, f64)> {
        let pt = pt.iter().map(|&x| F(x)).collect::<Vec<F<f64>>>();
        match self.integrands.get_mut(integrand) {
            Some(integrand) => {
                let settings = match integrand {
                    Integrand::GammaLoopIntegrand(integrand) => {
                        integrand.global_data.settings.clone()
                    }
                    Integrand::NewIntegrand(integrand) => integrand.get_settings().clone(),
                    _ => todo!(),
                };

                let res = inspect::inspect(
                    &settings,
                    integrand,
                    pt,
                    &term,
                    force_radius,
                    is_momentum_space,
                    use_f128,
                );

                Ok((res.re.0, res.im.0))
            }
            None => Err(exceptions::PyException::new_err(format!(
                "Could not find integrand {} available integrands: {:?}",
                integrand,
                self.integrands.keys().collect::<Vec<&String>>()
            ))),
        }
    }

    pub fn inspect_lmw_integrand(
        &mut self,
        integrand: &str,
        workspace_path: &str,
        use_f128: bool,
    ) -> PyResult<()> {
        match self.integrands.get_mut(integrand) {
            Some(integrand_struct) => {
                let new_settings = match integrand_struct {
                    Integrand::GammaLoopIntegrand(integrand_struct) => {
                        integrand_struct.global_data.settings.clone()
                    }
                    _ => todo!(),
                };

                let workspace_path = PathBuf::from(workspace_path);
                let path_to_state = workspace_path.join("integration_state");

                match fs::read(path_to_state) {
                    Ok(state_bytes) => {
                        let integration_state: IntegrationState =
                            bincode::decode_from_slice::<IntegrationState, _>(
                                &state_bytes,
                                bincode::config::standard(),
                            )
                            .expect("failed to obtain state")
                            .0;
                        let path_to_workspace_settings = workspace_path.join("settings.yaml");
                        let workspace_settings_string =
                            fs::read_to_string(path_to_workspace_settings)
                                .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

                        let mut workspace_settings: Settings =
                            serde_yaml::from_str(&workspace_settings_string)
                                .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

                        workspace_settings.general.debug = new_settings.general.debug;

                        let max_weight_samples = vec![
                            integration_state.integral.re.max_eval_positive_xs,
                            integration_state.integral.re.max_eval_negative_xs,
                            integration_state.integral.im.max_eval_positive_xs,
                            integration_state.integral.im.max_eval_negative_xs,
                        ];

                        for max_weight_sample in max_weight_samples
                            .into_iter()
                            .filter_map(std::convert::identity)
                        {
                            // bypass inspect function as it does not take a symbolica sample as input
                            let eval_result = integrand_struct.evaluate_sample(
                                &max_weight_sample,
                                F(0.0),
                                1,
                                use_f128,
                                Complex::new_zero(),
                            );

                            let eval = eval_result.integrand_result;

                            info!(
                            "\nFor input point xs: \n\n{}\n\nThe evaluation of integrand '{}' is:\n\n{}\n",
                            format!(
                                "( {:?} )",
                                max_weight_sample,
                            )
                            .blue(),
                            integrand,
                            format!("( {:+.16e}, {:+.16e} i)", eval.re, eval.im).blue(),
                        );
                        }

                        Ok(())
                    }
                    Err(_) => Err(exceptions::PyException::new_err(
                        "No previous run to extract max weight from".to_string(),
                    )),
                }
            }
            None => Err(exceptions::PyException::new_err(format!(
                "Could not find integrand {}",
                integrand
            ))),
        }
    }

    pub fn integrate_integrand(
        &mut self,
        integrand: &str,
        num_cores: usize,
        result_path: &str,
        workspace_path: &str,
        target: Option<(f64, f64)>,
    ) -> PyResult<Vec<(f64, f64)>> {
        let target = target.map(|(re, im)| (F(re), F(im)));
        match self.integrands.get_mut(integrand) {
            Some(integrand_enum) => match integrand_enum {
                Integrand::NewIntegrand(gloop_integrand) => {
                    let target = match target {
                        Some((re, im)) => Some(Complex::new(re, im)),
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

                            let state: IntegrationState =
                                bincode::decode_from_slice::<IntegrationState, _>(
                                    &state_bytes,
                                    bincode::config::standard(),
                                )
                                .expect("Could not deserialize state")
                                .0;

                            let path_to_workspace_settings = workspace_path.join("settings.yaml");
                            let workspace_settings_string =
                                fs::read_to_string(path_to_workspace_settings)
                                    .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

                            let workspace_settings: Settings =
                                serde_yaml::from_str(&workspace_settings_string)
                                    .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

                            // force the settings to be the same as the ones used in the previous integration
                            *gloop_integrand.get_mut_settings() = workspace_settings.clone();

                            print_integral_result(
                                &state.integral.re,
                                1,
                                state.iter,
                                "re",
                                target.map(|c| c.re),
                            );

                            print_integral_result(
                                &state.integral.im,
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

                    let settings = gloop_integrand.get_settings().clone();
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

                    Ok(result
                        .result
                        .iter()
                        .tuple_windows()
                        .map(|(re, im)| (re.0, im.0))
                        .collect())
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
        let master_node = self
            .master_node
            .as_mut()
            .expect("Could not get master node");

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
        let master_node = self
            .master_node
            .as_mut()
            .expect("could not get master node");

        let job_out_name = format!("job_{}_out", job_id);
        let job_out_path = Path::new(workspace_path).join(job_out_name);

        let output_file = std::fs::read(job_out_path)?;
        let batch_result: BatchResult =
            bincode::decode_from_slice(&output_file, bincode::config::standard())
                .expect("Could not deserialize batch")
                .0;

        master_node
            .process_batch_output(batch_result)
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        Ok(format!("Processed job {}", job_id))
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

    pub fn generate_integrands(&mut self, settings_yaml_str: &str) {
        let settings =
            serde_yaml::from_str::<Settings>(settings_yaml_str).expect("Could not parse settings");
        println!("calling generate_integrands");

        let integrands = self.process_list.generate_integrands(settings, &self.model);

        self.integrands = integrands;
    }

    pub fn export(&mut self, export_root: &str) -> PyResult<()> {
        let export_settings = ExportSettings {
            root_folder: PathBuf::from_str(export_root)?,
        };

        // check if the export root exists, if not create it, if it does return error
        if !export_settings.root_folder.exists() {
            fs::create_dir_all(&export_settings.root_folder)
                .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;
        } else {
            return Err(exceptions::PyException::new_err(
            "Export root already exists, please choose a different path or remove the existing directory".to_string(),
            ));
        }

        let mut state_file =
        // info!("Hi");
            fs::File::create(PathBuf::from(export_root).join("symbolica_state.bin"))?;

        State::export(&mut state_file)?;

        self.process_list
            .export_dot(&export_settings, &self.model)
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        let binary = bincode::encode_to_vec(&self.process_list, bincode::config::standard())
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        fs::write(PathBuf::from(export_root).join("process_list.bin"), binary)
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        let model_yaml = serde_yaml::to_string(&self.model.to_serializable())
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        fs::write(PathBuf::from(export_root).join("model.yaml"), model_yaml)
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        Ok(())
    }

    pub fn load(&mut self, import_root: &str) -> PyResult<()> {
        let import_root = PathBuf::from(import_root);

        let model_dir = import_root.join("model.yaml");
        self.load_model(model_dir.to_str().unwrap())
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        let state = State::import(
            &mut fs::File::open(import_root.join("symbolica_state.bin"))
                .map_err(|e| exceptions::PyException::new_err(e.to_string()))?,
            None,
        )
        .map_err(|e| {
            exceptions::PyException::new_err(format!(
                "Could not load state from {}: {}",
                import_root.join("symbolica_state.bin").display(),
                e
            ))
        })?;

        let context = GammaLoopContextContainer {
            state_map: &state,
            model: &self.model,
        };

        let process_list_data = fs::read(import_root.join("process_list.bin"))
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        let (process_list, _) = bincode::decode_from_slice_with_context(
            &process_list_data,
            bincode::config::standard(),
            context,
        )
        .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        self.process_list = process_list;

        Ok(())
    }
}

impl PythonWorker {
    fn printer_options(format: &str) -> PrintOptions {
        match format {
            "file" => PrintOptions::file(),
            "mathematica" => PrintOptions::mathematica(),
            "latex" => PrintOptions::latex(),
            _ => PrintOptions::default(),
        }
    }
}
