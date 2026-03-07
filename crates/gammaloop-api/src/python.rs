use gammalooprs::{
    cff::generation::{generate_cff_expression_from_subgraph, SurfaceCache},
    feyngen::diagram_generator::evaluate_overall_factor,
    graph::{self, FeynmanGraph, Graph, LMBext},
    initialisation::initialise,
    processes::{DotExportSettings, ProcessCollection},
    settings::{global::OrientationPattern, RuntimeSettings},
    utils::tracing::LogLevel,
};
use linnet::half_edge::{
    involution::{EdgeIndex, Orientation},
    subgraph::{ModifySubSet, SuBitGraph},
};
use numpy::{Complex64, IntoPyArray, PyArray1, PyReadonlyArray2};
use typed_index_collections::TiVec;

use crate::{
    commands::{
        import::model::ImportModel,
        inspect::{BatchedInspect, Inspect},
        Commands, Evaluate,
    },
    state::{ProcessRef, RunHistory, State},
    CLISettings, OneShot,
};
use ahash::{HashMap, HashMapExt};

use color_eyre::Result;
use eyre::eyre;

/*
use gammalooprs::feyngen::{
    FeynGenError, FeynGenFilter, FeynGenFilters, GenerationType, SelfEnergyFilterOptions,
    SnailFilterOptions, TadpolesFilterOptions,
};
*/

use git_version::git_version;
use itertools::{self, Itertools};
use pyo3::types::PyFloat;
use std::{path::PathBuf, str::FromStr};

use symbolica::{activate_oem_license, atom::AtomCore, parse};
const GIT_VERSION: &str = git_version!(fallback = "unavailable");

#[allow(unused)]
use pyo3::{
    exceptions,
    prelude::*,
    pyclass,
    pyclass::CompareOp,
    pyfunction, pymethods, pymodule,
    types::{PyComplex, PyModule, PyTuple, PyType},
    wrap_pyfunction, FromPyObject, PyObject, PyRef, Python,
};

#[pyfunction]
#[pyo3(name = "evaluate_graph_overall_factor")]
pub(crate) fn evaluate_graph_overall_factor(overall_factor: &str) -> Result<String> {
    let overall_factor = parse!(overall_factor);
    let overall_factor_evaluated = evaluate_overall_factor(overall_factor.as_view());
    Ok(overall_factor_evaluated.to_canonical_string())
}

#[pyfunction]
#[pyo3(name = "atom_to_canonical_string")]
pub(crate) fn atom_to_canonical_string(atom_str: &str) -> Result<String> {
    Ok(parse!(atom_str).to_canonical_string())
}

#[pymodule]
fn _gammaloop(_py: Python, m: &Bound<PyModule>) -> PyResult<()> {
    gammalooprs::initialisation::initialise().expect("initialization failed");
    gammalooprs::set_interrupt_handler();
    m.add_class::<GammaLoopAPI>()?;
    m.add_class::<LogLevel>()?;
    m.add_class::<DotExportSettings>()?;
    /*
    m.add_class::<PyFeynGenFilters>()?;
    m.add_class::<PySnailFilterOptions>()?;
    m.add_class::<PySewedFilterOptions>()?;
    m.add_class::<PyTadpolesFilterOptions>()?;
    m.add_class::<PySelfEnergyFilterOptions>()?;
    m.add_class::<GlobalSettings>()?;
    m.add_class::<RuntimeSettings>()?;

    // m.add_class::<PyFeynGenOptions>()?;
    m.add_class::<PyNumeratorAwareGroupingOption>()?;
    */
    m.add("git_version", GIT_VERSION)?;
    m.add_wrapped(wrap_pyfunction!(atom_to_canonical_string))?;
    m.add_wrapped(wrap_pyfunction!(evaluate_graph_overall_factor))?;
    Ok(())
}

/*
pub struct OutputOptions {}

/// There should only be 1 worker per instance of gammaloop.
#[pyclass(name = "Worker", unsendable)]
pub struct PythonWorker {
    pub model: Model,
    pub process_list: ProcessList,
    // pub integrands: HashMap<String, Integrand>,
    pub master_node: Option<MasterNode>,
}

impl Clone for PythonWorker {
    fn clone(&self) -> PythonWorker {
        PythonWorker {
            model: self.model.clone(),
            process_list: self.process_list.clone(),
            // integrands: self.integrands.clone(),
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
    pub(crate) fn __new__(
        veto_snails_attached_to_massive_lines: Option<bool>,
        veto_snails_attached_to_massless_lines: Option<bool>,
        veto_only_scaleless_snails: Option<bool>,
    ) -> Result<PySnailFilterOptions> {
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

    pub(crate) fn __str__(&self) -> Result<String> {
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
    pub(crate) fn __new__(
        veto_self_energy_of_massive_lines: Option<bool>,
        veto_self_energy_of_massless_lines: Option<bool>,
        veto_only_scaleless_self_energy: Option<bool>,
    ) -> Result<PySelfEnergyFilterOptions> {
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

    pub(crate) fn __str__(&self) -> Result<String> {
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
    pub(crate) fn __new__(filter_tadpoles: Option<bool>) -> Result<PySewedFilterOptions> {
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
    pub(crate) fn __new__(
        veto_tadpoles_attached_to_massive_lines: Option<bool>,
        veto_tadpoles_attached_to_massless_lines: Option<bool>,
        veto_only_scaleless_tadpoles: Option<bool>,
    ) -> Result<PyTadpolesFilterOptions> {
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

    pub(crate) fn __str__(&self) -> Result<String> {
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
    pub(crate) fn __str__(&self) -> Result<String> {
        Ok(self.filters.iter().map(|f| format!(" > {}", f)).join("\n"))
    }

    #[new]
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn __new__(
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
    ) -> Result<PyFeynGenFilters> {
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

fn feyngen_to_python_error(error: FeynGenError) -> PyErr {
    exceptions::PyValueError::new_err(format!("Feynman diagram generator error | {error}"))
}

#[pyclass(name = "NumeratorAwareGroupingOption")]
pub struct PyNumeratorAwareGroupingOption {
    pub grouping_options: NumeratorAwareGraphGroupingOption,
}

#[pymethods]
impl PyNumeratorAwareGroupingOption {
    #[new]
    pub(crate) fn __new__(
        numerator_aware_grouping_option: Option<String>,
        compare_canonized_numerator: Option<bool>,
        number_of_samples_for_numerator_comparisons: Option<usize>,
        consider_internal_masses_only_in_numerator_isomorphisms: Option<bool>,
        fully_numerical_substitution_when_comparing_numerators: Option<bool>,
        numerical_samples_seed: Option<u16>,
        symmetric_polarizations: Option<bool>,
    ) -> Result<PyNumeratorAwareGroupingOption> {
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
                symmetric_polarizations,
            )
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?,
        })
    }

    pub(crate) fn __str__(&self) -> Result<String> {
        Ok(format!("{}", self.grouping_options))
    }
}
*/

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(unsendable, name = "GammaLoopAPI")
)]
struct GammaLoopAPI {
    gammaloop_state: State,
    cli_settings: CLISettings,
    #[allow(unused)]
    run_history: RunHistory,
    default_runtime_settings: RuntimeSettings,
}

// TODO: Improve error broadcasting to Python everywhere so as to show rust backtrace
#[pymethods]
impl GammaLoopAPI {
    #[new]
    #[pyo3(signature = (state_folder=PathBuf::from("./gammaloop_state"), log_file_name=None, log_level=None))]
    pub fn new_python(
        state_folder: PathBuf,
        log_file_name: Option<String>,
        log_level: Option<LogLevel>,
    ) -> Result<Self> {
        if option_env!("NO_SYMBOLICA_OEM_LICENSE").is_none() {
            activate_oem_license!("SYMBOLICA_OEM_KEY_23177b25");
        };
        initialise().map_err(|e| {
            exceptions::PyException::new_err(format!("Could not initialize GammaLoop API: {}", e))
        })?;
        let mut one_shot = OneShot {
            // Path to the a run file to execute
            run_history: None,
            run_block_names: vec![],
            // Don't actually run anything, just build up run card
            dry_run: false,
            // Path to the state folder
            state_folder,
            // Python API takes an explicit state folder argument.
            state_folder_explicitly_set: true,
            // Path to the model file
            model_file: None,
            // Skip saving state on exit
            no_save_state: false,
            // Save state to file after each call
            override_state: false,
            // Set the name of the file containing all traces from gammaloop (logs) for current session
            trace_logs_filename: log_file_name,
            // Set log level for current session
            level: log_level,
            // Try to serialize using strings when saving run history
            no_try_strings: false,
            command: None,
            fresh_state: false,
        };
        let (state, run_history, cli_settings, default_runtime_settings) =
            one_shot.load().map_err(|e| {
                exceptions::PyException::new_err(format!(
                    "Could not load or create GammaLoop state: {}",
                    e
                ))
            })?;
        Ok(GammaLoopAPI {
            gammaloop_state: state,
            run_history,
            cli_settings,
            default_runtime_settings,
        })
    }

    pub fn inspect<'py>(
        &mut self,
        py: Python<'py>,
        process_id: Option<usize>,
        integrand_name: Option<String>,
        point: Vec<f64>,
        use_arb_prec: bool,
        force_radius: bool,
        momentum_space: bool,
        discrete_dim: Vec<usize>,
    ) -> Result<(Bound<'py, PyComplex>, Option<Bound<'py, PyFloat>>)> {
        let res = Inspect {
            process: process_id.map(ProcessRef::Id),
            integrand_name,
            point,
            use_arb_prec,
            force_radius,
            momentum_space,
            discrete_dim,
        }
        .run(&mut self.gammaloop_state)?;
        Ok((
            PyComplex::from_doubles(py, res.1.re, res.1.im),
            res.0.map(|j| j.into_pyobject(py).unwrap()),
        ))
    }

    pub fn batched_inspect<'py>(
        &mut self,
        py: Python<'py>,
        points: PyReadonlyArray2<'py, f64>,
        discrete_dims: PyReadonlyArray2<'py, usize>,
        process_id: Option<usize>,
        integrand_name: Option<String>,
        use_arb_prec: bool,
        momentum_space: bool,
    ) -> Result<(
        Bound<'py, PyArray1<Complex64>>,
        Option<Bound<'py, PyArray1<f64>>>,
    )> {
        let points_rust = points.as_array();
        let discrete_dims_rust = discrete_dims.as_array();

        let batched_inspect = BatchedInspect {
            process_id,
            integrand_name,
            use_arb_prec,
            momentum_space,
            points: points_rust,
            discrete_dims: discrete_dims_rust,
        };

        let (res, res_jac) = batched_inspect.run(&mut self.gammaloop_state).unwrap();
        let res_map = res
            .mapv_into_any(|c| Complex64::new(c.re, c.im))
            .into_pyarray(py);

        let res_jac_map = res_jac.map(|r| r.into_pyarray(py));

        Ok((res_map, res_jac_map))
    }

    #[pyo3(name="import_graphs", signature = (graphs, process_name=None, process_id=None, integrand_name=None, format="dot".into(), overwrite=false, append=false))]
    pub(crate) fn import_graphs_python(
        &mut self,
        graphs: String,
        process_name: Option<String>,
        process_id: Option<usize>,
        integrand_name: Option<String>,
        format: String,
        overwrite: bool,
        append: bool,
    ) -> Result<()> {
        if overwrite && append {
            return Err(eyre!(
                "Cannot use both overwrite and append options when importing graphs"
            ));
        }
        let (graphs, process_name) = match format.as_str() {
            "dot" => {
                let dot_path = PathBuf::from(&graphs);
                if !dot_path.exists() {
                    return Err(eyre!("Path does not exist: {}", dot_path.display()));
                }

                let process_name = process_name
                    .unwrap_or(dot_path.file_stem().unwrap().to_string_lossy().into_owned());

                let graphs = Graph::from_path(&dot_path, &self.gammaloop_state.model)
                    .map_err(|e| eyre!("Could not parse graphs from path: {}", e.to_string()))?;
                (graphs, Some(process_name))
            }
            "string" => {
                if process_id.is_none() && process_name.is_none() {
                    return Err(eyre!(
                        "When importing graphs from string, either a process id or a process name must be provided"
                    ));
                }

                let graphs = Graph::from_string(&graphs, &self.gammaloop_state.model)
                    .map_err(|e| eyre!("Could not parse graphs from string: {}", e.to_string()))?;
                (graphs, process_name)
            }
            other => {
                return Err(eyre!(
                    "Unknown graph import format: {}. Supported formats are 'yaml' and 'json'",
                    other
                ))
            }
        };

        self.gammaloop_state.import_graphs(
            graphs,
            process_name,
            process_id,
            integrand_name,
            overwrite,
            append,
        )
    }

    #[pyo3(name="get_lmbs", signature = (graphs, format="dot".into()))]
    pub(crate) fn get_lmbs(
        &self,
        graphs: String,
        format: String,
    ) -> Result<
        Vec<
            Vec<(
                Vec<usize>, // loop edges
                Vec<usize>, // external edges
                HashMap<usize, (Vec<i8>, Vec<i8>)>,
            )>,
        >,
    > {
        let graphs = match format.as_str() {
            "dot" => {
                let dot_path = PathBuf::from(&graphs);
                if !dot_path.exists() {
                    return Err(eyre!("Path does not exist: {}", dot_path.display()));
                }

                let graphs = Graph::from_path(&dot_path, &self.gammaloop_state.model)
                    .map_err(|e| eyre!("Could not parse graphs from path: {}", e.to_string()))?;
                graphs
            }
            "string" => {
                let graphs = Graph::from_string(&graphs, &self.gammaloop_state.model)
                    .map_err(|e| eyre!("Could not parse graphs from string: {}", e.to_string()))?;
                graphs
            }
            other => {
                return Err(eyre!(
                    "Unknown graph import format: {}. Supported formats are sting and dot",
                    other
                ))
            }
        };

        let lmbs = graphs
            .into_iter()
            .map(|g| g.generate_loop_momentum_bases_of(&g.no_dummy()))
            .collect_vec();

        Ok(lmbs
            .into_iter()
            .map(|lmb_list| {
                lmb_list
                    .into_iter()
                    .map(|lmb| {
                        let loop_edges = lmb.loop_edges.into_iter().map(|e| e.0).collect_vec();
                        let external_edges = lmb.ext_edges.into_iter().map(|e| e.0).collect_vec();
                        let mut edge_signatures = HashMap::new();

                        for (edge_id, signature) in lmb.edge_signatures.into_iter() {
                            let loop_part = signature
                                .internal
                                .into_iter()
                                .map(|e| e as i8)
                                .collect_vec();

                            let ext_part = signature
                                .external
                                .into_iter()
                                .map(|e| e as i8)
                                .collect_vec();

                            edge_signatures.insert(edge_id.0, (loop_part, ext_part));
                        }
                        (loop_edges, external_edges, edge_signatures)
                    })
                    .collect()
            })
            .collect())
    }

    #[pyo3(name="get_orientations", signature = (graph_name, process_id=None, integrand_name=None))]
    pub(crate) fn get_orientations(
        &self,
        graph_name: String,
        process_id: Option<usize>,
        integrand_name: Option<String>,
    ) -> Result<Vec<HashMap<usize, i8>>> {
        let (pid, name) = self
            .gammaloop_state
            .process_list
            .find_integrand(process_id, integrand_name.as_ref())
            .map_err(|e| {
                exceptions::PyException::new_err(format!("Could not find integrand: {}", e))
            })?;

        let orientations = match &self.gammaloop_state.process_list.processes[pid].collection {
            ProcessCollection::Amplitudes(amplitudes) => {
                let cff = amplitudes
                    .get(&name)
                    .unwrap()
                    .graphs
                    .iter()
                    .find(|g| g.graph.name == graph_name)
                    .as_ref()
                    .unwrap()
                    .derived_data
                    .cff_expression
                    .as_ref()
                    .unwrap();

                cff.orientations
                    .iter()
                    .map(|or_data| or_data.data.orientation.clone())
                    .collect_vec()
            }

            ProcessCollection::CrossSections(cross_sections) => {
                let cff = cross_sections
                    .get(&name)
                    .unwrap()
                    .supergraphs
                    .iter()
                    .find(|g| g.graph.name == graph_name)
                    .as_ref()
                    .unwrap()
                    .derived_data
                    .cff_expression
                    .as_ref()
                    .unwrap();

                cff.orientation_data
                    .iter()
                    .map(|or_data| or_data.orientation.clone())
                    .collect_vec()
            }
        };

        Ok(orientations
            .into_iter()
            .map(|orientation| {
                let mut result = HashMap::new();
                for (edge_id, direction) in orientation.into_iter() {
                    let direction = match direction {
                        Orientation::Default => 1,
                        Orientation::Reversed => -1,
                        Orientation::Undirected => 0,
                    };
                    result.insert(edge_id.0, direction);
                }
                result
            })
            .collect())
    }

    #[pyo3(name = "get_model")]
    pub(crate) fn get_model(&self) -> PyResult<String> {
        let serializable_model = self.gammaloop_state.model.to_serializable();
        serde_json::to_string(&serializable_model).map_err(|e| {
            exceptions::PyException::new_err(format!("Could not serialize model: {}", e))
        })
    }

    #[pyo3(name="evaluate", signature = (process_id=None, graphs_group_name=None, result_path=None, numerical=true, number_of_terms_in_epsilon_expansion=None))]
    pub(crate) fn evaluate_python(
        &mut self,
        process_id: Option<usize>,
        graphs_group_name: Option<String>,
        result_path: Option<PathBuf>,
        numerical: bool,
        number_of_terms_in_epsilon_expansion: Option<usize>,
    ) -> PyResult<String> {
        Evaluate {
            process: process_id.map(ProcessRef::Id),
            graphs_group_name,
            result_path,
            numerical,
            number_of_terms_in_epsilon_expansion,
        }
        .run(
            &mut self.gammaloop_state,
            &self.cli_settings,
            &self.default_runtime_settings,
        )
        .map_err(|e| {
            exceptions::PyException::new_err(format!(
                "Could not load or create GammaLoop state: {}",
                e
            ))
        })
        .map(|res| res.to_canonical_string())
    }

    #[pyo3(name="import_model", signature = (model_specifier, simplify_model=true))]
    pub(crate) fn import_model_python(
        &mut self,
        model_specifier: PathBuf,
        simplify_model: bool,
    ) -> PyResult<()> {
        ImportModel {
            path: model_specifier,
            simplify_model,
        }
        .run(&mut self.gammaloop_state)
        .map_err(|e| exceptions::PyException::new_err(format!("Could not import model: {}", e)))
    }

    #[pyo3(name="list_outputs", signature = ())]
    pub(crate) fn list_outputs(
        &mut self,
    ) -> PyResult<(HashMap<String, usize>, HashMap<String, usize>)> {
        let mut all_amplitudes = HashMap::new();
        let mut all_cross_sections = HashMap::new();
        for (p_id, p) in self
            .gammaloop_state
            .process_list
            .processes
            .iter()
            .enumerate()
        {
            match &p.collection {
                ProcessCollection::Amplitudes(amplitudes) => {
                    amplitudes.keys().cloned().for_each(|amp| {
                        all_amplitudes.insert(amp, p_id);
                    });
                }
                ProcessCollection::CrossSections(cross_sections) => {
                    cross_sections.keys().cloned().for_each(|cs| {
                        all_cross_sections.insert(cs, p_id);
                    });
                }
            }
        }
        Ok((all_amplitudes, all_cross_sections))
    }

    #[pyo3(name="get_integrand_settings", signature = (process_id=None, integrand_name=None))]
    pub(crate) fn get_integrand_settings(
        &mut self,
        process_id: Option<usize>,
        integrand_name: Option<String>,
    ) -> PyResult<RuntimeSettings> {
        let (pid, name) = self
            .gammaloop_state
            .process_list
            .find_integrand(process_id, integrand_name.as_ref())
            .map_err(|e| {
                exceptions::PyException::new_err(format!("Could not find integrand: {}", e))
            })?;
        match &self.gammaloop_state.process_list.processes[pid].collection {
            ProcessCollection::Amplitudes(amplitudes) => {
                match &amplitudes.get(&name).unwrap().integrand {
                    Some(integrand) => Ok::<_, PyErr>(integrand.get_settings().clone()),
                    None => Err(exceptions::PyException::new_err(
                        "Integrand for selected amplitude not yet generated",
                    )),
                }
            }
            ProcessCollection::CrossSections(cross_sections) => {
                match &cross_sections.get(&name).unwrap().integrand {
                    Some(integrand) => Ok::<_, PyErr>(integrand.get_settings().clone()),
                    None => Err(exceptions::PyException::new_err(
                        "Integrand for selected cross-section not yet generated",
                    )),
                }
            }
        }
    }

    #[pyo3(name="get_dot_files", signature = (process_id=None, integrand_name=None,settings=DotExportSettings::default()))]
    pub(crate) fn get_dot_files(
        &mut self,
        process_id: Option<usize>,
        integrand_name: Option<String>,
        settings: DotExportSettings,
    ) -> PyResult<String> {
        let (pid, name) = self
            .gammaloop_state
            .process_list
            .find_integrand(process_id, integrand_name.as_ref())
            .map_err(|e| {
                exceptions::PyException::new_err(format!("Could not find integrand: {}", e))
            })?;
        let mut dot_output = String::new();
        match &self.gammaloop_state.process_list.processes[pid].collection {
            ProcessCollection::Amplitudes(amplitudes) => match &amplitudes.get(&name) {
                Some(amplitude) => {
                    amplitude
                        .write_dot_fmt(&mut dot_output, &settings)
                        .map_err(|e| {
                            exceptions::PyException::new_err(format!(
                                "Could not write DOT format for amplitude {}: {}",
                                name, e
                            ))
                        })
                }
                None => Err(exceptions::PyException::new_err(format!(
                    "Could not find amplitude named {}",
                    name
                ))),
            },
            ProcessCollection::CrossSections(cross_sections) => match &cross_sections.get(&name) {
                Some(cross_section) => cross_section
                    .write_dot_fmt(&mut dot_output, &settings)
                    .map_err(|e| {
                        exceptions::PyException::new_err(format!(
                            "Could not write DOT format for amplitude {}: {}",
                            name, e
                        ))
                    }),
                None => Err(exceptions::PyException::new_err(format!(
                    "Could not find cross-section named {}",
                    name
                ))),
            },
        }?;
        Ok(dot_output)
    }

    #[pyo3(name="run", signature = (command,))]
    pub(crate) fn run_command(&mut self, command: String) -> PyResult<()> {
        let cmd = Commands::from_str(&command)?;

        cmd.run(
            &mut self.gammaloop_state,
            &mut self.run_history,
            &mut self.cli_settings,
            &mut self.default_runtime_settings,
        )
        .map_err(|e| {
            exceptions::PyException::new_err(format!(
                "Execution of command '{}' failed: {}",
                command, e
            ))
        })
        .map(|_| ())
    }

    #[pyo3(name = "generate_cff", signature = (dot_string, subgraph_nodes, reverse_dangling,orientation_pattern=None))]
    pub(crate) fn generate_cff(
        &self,
        dot_string: String,
        subgraph_nodes: Vec<String>,
        reverse_dangling: Vec<usize>,
        orientation_pattern: Option<String>,
    ) -> PyResult<Vec<(HashMap<usize, i32>, String)>> {
        let graph = Graph::from_string(dot_string, &self.gammaloop_state.model)
            .unwrap()
            .pop()
            .unwrap();

        let reverse_dangling = reverse_dangling
            .into_iter()
            .map(EdgeIndex::from)
            .collect_vec();

        let subgraph: SuBitGraph = if subgraph_nodes.is_empty() {
            graph.full_filter()
        } else {
            let mut result: SuBitGraph = graph.empty_subgraph();
            for (_node_id, neighbors, vertex) in graph.iter_nodes() {
                if subgraph_nodes.contains(&vertex.name.to_string()) {
                    neighbors.for_each(|hedge| result.add(hedge));
                }
            }
            result
        };

        let mut surface_cache = SurfaceCache {
            esurface_cache: TiVec::new(),
            hsurface_cache: TiVec::new(),
        };

        let cff = generate_cff_expression_from_subgraph(
            &graph.underlying,
            &subgraph,
            &None,
            &reverse_dangling,
            &graph.get_edges_in_initial_state_cut(),
            &mut surface_cache,
        )
        .map_err(|e| {
            exceptions::PyException::new_err(format!("Could not generate CFF expression: {}", e))
        })?;

        let or_pattern: OrientationPattern = match orientation_pattern {
            Some(pattern) => {
                serde_json::from_str(format!("{{\"pat\": \"{}\"}}", pattern).as_str()).unwrap()
            }
            None => OrientationPattern { pat: None },
        };

        let atoms = cff.get_orientation_atoms_with_data(or_pattern);
        let inverse_energies = graph::get_cff_inverse_energy_product_impl(&graph, &subgraph, &[]);

        let result = atoms
            .into_iter()
            .map(|(atom, orientation_data)| {
                let energy_sub = cff.surfaces.substitute_energies(&atom, &[]) * &inverse_energies;
                let string_atom = energy_sub.to_string();

                let mut orientation_as_hashmap = HashMap::new();
                for (edge_id, direction) in orientation_data.orientation.into_iter() {
                    let direction = match direction {
                        Orientation::Default => 1,
                        Orientation::Reversed => -1,
                        Orientation::Undirected => 0,
                    };
                    orientation_as_hashmap.insert(edge_id.0, direction);
                }
                (orientation_as_hashmap, string_atom)
            })
            .collect_vec();

        Ok(result)
    }

    #[pyo3(
        name = "generate_cff_as_json_string",
        signature = (dot_string, subgraph_nodes, reverse_dangling, orientation_pattern = None)
    )]
    pub(crate) fn generate_cff_as_json_string(
        &self,
        dot_string: String,
        subgraph_nodes: Vec<String>,
        reverse_dangling: Vec<usize>,
        orientation_pattern: Option<String>,
    ) -> PyResult<String> {
        let _ = orientation_pattern;
        let graph = Graph::from_string(dot_string, &self.gammaloop_state.model)
            .unwrap()
            .pop()
            .unwrap();

        let reverse_dangling = reverse_dangling
            .into_iter()
            .map(EdgeIndex::from)
            .collect_vec();

        let subgraph: SuBitGraph = if subgraph_nodes.is_empty() {
            graph.full_filter()
        } else {
            let mut result: SuBitGraph = graph.empty_subgraph();
            for (_node_id, neighbors, vertex) in graph.iter_nodes() {
                if subgraph_nodes.contains(&vertex.name.to_string()) {
                    neighbors.for_each(|hedge| result.add(hedge));
                }
            }
            result
        };

        let mut surface_cache = SurfaceCache {
            esurface_cache: TiVec::new(),
            hsurface_cache: TiVec::new(),
        };

        let cff = generate_cff_expression_from_subgraph(
            &graph.underlying,
            &subgraph,
            &None,
            &reverse_dangling,
            &graph.get_edges_in_initial_state_cut(),
            &mut surface_cache,
        )
        .map_err(|e| {
            exceptions::PyException::new_err(format!("Could not generate CFF expression: {}", e))
        })?;

        let json_string = serde_json::to_string(&cff).map_err(|e| {
            exceptions::PyException::new_err(format!("Could not serialize cff to json: {}", e))
        })?;

        Ok(json_string)
    }

    /*

    #[allow(clippy::too_many_arguments)]
    #[pyo3(name="generate", signature = (process_name, numerator_aware_isomorphism_grouping=None, filter_self_loop=None, graph_prefix=None, selected_graphs=None, vetoed_graphs=None, loop_momentum_bases=None, global_prefactor_color=None, global_prefactor_colorless=None, num_threads=None))]
    pub(crate) fn generate_diagrams(
        &mut self,
        //generation_options: PyRef<PyFeynGenOptions>,
        process_name: String,
        numerator_aware_isomorphism_grouping: Option<PyRef<PyNumeratorAwareGroupingOption>>,
        filter_self_loop: Option<bool>,
        graph_prefix: Option<String>,
        selected_graphs: Option<Vec<String>>,
        vetoed_graphs: Option<Vec<String>>,
        loop_momentum_bases: Option<HashMap<String, Vec<String>>>,
        global_prefactor_color: Option<String>,
        global_prefactor_colorless: Option<String>,
        num_threads: Option<usize>,
    ) -> Result<Vec<String>> {
        todo!();
        //if self.model.is_empty() {
        //    return Err(eyre!(
        //        "A physics model must be loaded before generating diagrams"
        //    ));
        //}
        //let feyngen_options = generation_options.options.clone();

        //// clone some of the options that will be used in the process definition
        //let initial_pdgs = feyngen_options.initial_pdgs.clone();
        //let final_pdgs_lists = feyngen_options.final_pdgs_lists.clone();
        //let generation_type = feyngen_options.generation_type.clone();
        //let amplitude_filters = feyngen_options.amplitude_filters.clone();
        //let cross_section_filters = feyngen_options.cross_section_filters.clone();

        //let diagram_generator = FeynGen::new(feyngen_options);

        //let mut global_prefactor = GlobalPrefactor::default();
        //if let Some(global_prefactor_color) = global_prefactor_color {
        //    global_prefactor.num = parse!(&global_prefactor_color);
        //}
        //if let Some(global_prefactor_colorless) = global_prefactor_colorless {
        //    global_prefactor.projector = parse!(&global_prefactor_colorless);
        //}

        //let _diagrams = diagram_generator
        //    .generate(
        //        &self.model,
        //        &numerator_aware_isomorphism_grouping
        //            .map(|o| o.grouping_options.clone())
        //            .unwrap_or(NumeratorAwareGraphGroupingOption::NoGrouping),
        //        filter_self_loop.unwrap_or(false),
        //        graph_prefix.unwrap_or("GL".to_string()),
        //        selected_graphs,
        //        vetoed_graphs,
        //        loop_momentum_bases,
        //        global_prefactor,
        //        num_threads,
        //    )
        //    .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        // let res = Ok(diagrams
        //     .iter()
        //     .map(|d| serde_yaml::to_string(&SerializableGraph::from_graph(d)).unwrap())
        //     .collect());

        // load everything into processlist
        //let (n_unresolved, unresolved_cut_content) =
        //    diagram_generator.unresolved_cut_content(&self.model);

        //let process_definition = ProcessDefinition {
        //    initial_pdgs,
        //    final_pdgs_lists,
        //    n_unresolved,
        //    unresolved_cut_content: unresolved_cut_content.into_iter().collect(),
        //    amplitude_filters,
        //    cross_section_filters,
        //};

        //let process = Process::from_graph_list(
        //    proccess_name.into(),
        //    vec![],
        //    generation_type,
        //    process_definition,
        //    None,
        //)
        //.map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        //self.process_list.add_process(process);

        //Ok(Vec::new())
    }

    // pub(crate) fn export_expressions(
    //     &mut self,
    //     export_root: &str,
    //     amplitued_list: Vec<String>,
    //     format: &str,
    //     export_yaml_str: &str,
    // ) -> Result<String> {
    //     let export_settings: ProcessSettings = serde_yaml::from_str(export_yaml_str)
    //         .map_err(|e| exceptions::PyException::new_err(e.to_string()))
    //         .unwrap();

    //     for amplitude in amplitued_list.into_iter() {
    //         match Amplitude::from_yaml_str(&self.model, amplitude) {
    //             Ok(amp) => {
    //                 amp.map(|a| a.map(|ag| ag.forget_type()))
    //                     .export_expressions(
    //                         export_root,
    //                         Self::printer_options(format),
    //                         &export_settings,
    //                     )
    //                     .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;
    //             }
    //             Err(e) => return Err(exceptions::PyException::new_err(e.to_string())),
    //         }
    //     }
    //     Ok("Exported expressions".to_string())
    // }

    // pub(crate) fn export_coupling_replacement_rules(
    //     &self,
    //     export_root: &str,
    //     format: &str,
    // ) -> Result<String> {
    //     self.model
    //         .export_coupling_replacement_rules(export_root, Self::printer_options(format))
    //         .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;
    //     Ok("Exported coupling substitutions".to_string())
    // }

    pub fn integrate<'py>(
        &mut self,
        py: Python<'py>,
        process_id: Option<usize>,
        integrand_name: Option<String>,
        result_path: Option<PathBuf>,
        n_cores: usize,
        workspace_path: PathBuf,
        target: Option<Vec<f64>>,
        restart: bool,
    ) -> Result<Vec<Bound<'py, PyComplex>>> {
        let a = Integrate {
            process: process_id.map(ProcessRef::Id),
            integrand_name,
            result_path,
            workspace_path: Some(workspace_path),
            n_cores: Some(n_cores),
            target,
            restart,
        }
        .run(
            &mut self.gammaloop_state,
            &CLISettings {
                ..Default::default()
            },
        )?;

        Ok(vec![PyComplex::from_doubles(
            py,
            a.result.re.0,
            a.result.im.0,
        )])
    }

    // pub(crate) fn inspect_lmw_integrand(
    //     &mut self,
    //     integrand_name: &str,
    //     workspace_path: &str,
    //     use_f128: bool,
    // ) -> Result<()> {
    //     let integrand = self.process_list.get_integrand_mut(0, integrand_name)?;

    //     let settings = integrand.get_settings().clone();
    //     let workspace_path = PathBuf::from(workspace_path);
    //     let path_to_state = workspace_path.join("integration_state");

    //     // match fs::read(path_to_state) {
    //     //     Ok(state_bytes) => {
    //     //         let integration_state: IntegrationState =
    //     //             bincode::decode_from_slice::<IntegrationState, _>(
    //     //                 &state_bytes,
    //     //                 bincode::config::standard(),
    //     //             )
    //     //             .expect("failed to obtain state")
    //     //             .0;
    //     //         let path_to_workspace_settings = workspace_path.join("settings.yaml");
    //     //         let workspace_settings_string = fs::read_to_string(path_to_workspace_settings)
    //     //             .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

    //     //         let mut workspace_settings: Settings =
    //     //             serde_yaml::from_str(&workspace_settings_string)
    //     //                 .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

    //     //         workspace_settings.general.debug = new_settings.general.debug;

    //     //         let max_weight_samples = vec![
    //     //             integration_state.integral.re.max_eval_positive_xs,
    //     //             integration_state.integral.re.max_eval_negative_xs,
    //     //             integration_state.integral.im.max_eval_positive_xs,
    //     //             integration_state.integral.im.max_eval_negative_xs,
    //     //         ];

    //     //         for max_weight_sample in max_weight_samples
    //     //             .into_iter()
    //     //             .filter_map(std::convert::identity)
    //     //         {
    //     //             // bypass inspect function as it does not take a symbolica sample as input
    //     //             let eval_result = integrand_struct.evaluate_sample(
    //     //                 &max_weight_sample,
    //     //                 F(0.0),
    //     //                 1,
    //     //                 use_f128,
    //     //                 Complex::new_zero(),
    //     //             );

    //     //             let eval = eval_result.integrand_result;

    //     //             info!(
    //     //                     "\nFor input point xs: \n\n{}\n\nThe evaluation of integrand '{}' is:\n\n{}\n",
    //     //                     format!(
    //     //                         "( {:?} )",
    //     //                         max_weight_sample,
    //     //                     )
    //     //                     .blue(),
    //     //                     integrand_name,
    //     //                     format!("( {:+.16e}, {:+.16e} i)", eval.re, eval.im).blue(),
    //     //                 );
    //     //         }

    //     //         Ok(())
    //     //     }
    //     // }
    //     //
    //     Ok(())
    // }

    // pub(crate) fn load_master_node(&mut self, integrand: &str) -> Result<String> {
    //     let selected_integrand = self.process_list.get_integrand_mut(0, integrand)?;

    //     let grid = selected_integrand.create_grid();
    //     let integrator_settings = selected_integrand.get_integrator_settings();

    //     let master_node = MasterNode::new(grid, integrator_settings);
    //     self.master_node = Some(master_node);

    //     Ok(format!("Initialized master grid for {}", integrand))
    // }

    // pub(crate) fn write_batch_input(
    //     &mut self,
    //     num_cores: usize,
    //     num_samples: usize,
    //     export_grid: bool,
    //     output_accumulator: bool,
    //     workspace_path: &str,
    //     job_id: usize,
    // ) -> Result<String> {
    //     let master_node = self
    //         .master_node
    //         .as_mut()
    //         .expect("Could not get master node");

    //     // extract the integrated phase in a hacky way
    //     match master_node
    //         .write_batch_input(
    //             num_cores,
    //             num_samples,
    //             export_grid,
    //             output_accumulator,
    //             workspace_path,
    //             job_id,
    //         )
    //         .map_err(|e| exceptions::PyException::new_err(e.to_string()))
    //     {
    //         Ok(_) => Ok(format!("Wrote batch input for job {}", job_id)),
    //         Err(e) => Err(e),
    //     }
    // }

    // pub(crate) fn process_batch_output(
    //     &mut self,
    //     workspace_path: &str,
    //     job_id: usize,
    // ) -> Result<String> {
    //     let master_node = self
    //         .master_node
    //         .as_mut()
    //         .expect("could not get master node");

    //     let job_out_name = format!("job_{}_out", job_id);
    //     let job_out_path = Path::new(workspace_path).join(job_out_name);

    //     let output_file = std::fs::read(job_out_path)?;
    //     let batch_result: BatchResult =
    //         bincode::decode_from_slice(&output_file, bincode::config::standard())
    //             .expect("Could not deserialize batch")
    //             .0;

    //     master_node
    //         .process_batch_output(batch_result)
    //         .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

    //     Ok(format!("Processed job {}", job_id))
    // }

    // pub(crate) fn display_master_node_status(&self) {
    //     if let Some(master_node) = &self.master_node {
    //         master_node.display_status();
    //     }
    // }

    // pub(crate) fn update_iter(&mut self) {
    //     if let Some(master_node) = &mut self.master_node {
    //         master_node.update_iter();
    //     }
    // }


    */
}
