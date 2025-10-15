use gammalooprs::{
    feyngen::{
        //  self, diagram_generator::FeynGen, FeynGenError, FeynGenFilters, FeynGenOptions,
        diagram_generator::evaluate_overall_factor,
        NumeratorAwareGraphGroupingOption,
        SewedFilterOptions,
    },
    initialisation::initialise,
    integrate::MasterNode,
    model::{InputParamCard, Model},
    numerator::GlobalPrefactor,
    processes::{Process, ProcessDefinition, ProcessList},
    settings::{GlobalSettings, RuntimeSettings},
    utils::{FloatLike, F},
};
use numpy::{
    Complex64, IntoPyArray, PyArray, PyArray1, PyArray2, PyReadonlyArray2, PyReadonlyArrayDyn,
};

use crate::{
    commands::{
        inspect::{BatchedInspect, Inspect},
        integrate::Integrate,
    },
    state::State,
    CLISettings,
};
use ahash::HashMap;

use color_eyre::Result;
use eyre::eyre;
use gammalooprs::feyngen::{
    FeynGenError, FeynGenFilter, FeynGenFilters, GenerationType, SelfEnergyFilterOptions,
    SnailFilterOptions, TadpolesFilterOptions,
};
use git_version::git_version;
use itertools::{self, Itertools};
use pyo3::types::PyFloat;
use std::{convert::Infallible, default, path::PathBuf, str::FromStr};

use symbolica::{atom::AtomCore, parse};
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
    m.add_class::<State>()?;
    m.add_class::<PyFeynGenFilters>()?;
    m.add_class::<PySnailFilterOptions>()?;
    m.add_class::<PySewedFilterOptions>()?;
    m.add_class::<PyTadpolesFilterOptions>()?;
    m.add_class::<PySelfEnergyFilterOptions>()?;
    m.add_class::<GlobalSettings>()?;
    m.add_class::<RuntimeSettings>()?;
    // m.add_class::<PyFeynGenOptions>()?;
    m.add_class::<PyNumeratorAwareGroupingOption>()?;
    m.add("git_version", GIT_VERSION)?;
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

//#[pyclass(name = "FeynGenOptions")]
//pub struct PyFeynGenOptions {
//    pub options: FeynGenOptions,
//}
//impl<'a> FromPyObject<'a> for PyFeynGenOptions {
//    fn extract_bound(ob: &Bound<'a, PyAny>) -> PyResult<Self> {
//        if let Ok(a) = ob.extract::<PyFeynGenOptions>() {
//            Ok(PyFeynGenOptions { options: a.options })
//        } else {
//            Err(exceptions::PyValueError::new_err(
//                "Not a valid Feynman generation option structure",
//            ))
//        }
//    }
//}

fn feyngen_to_python_error(error: FeynGenError) -> PyErr {
    exceptions::PyValueError::new_err(format!("Feynam diagram generator error | {error}"))
}

//#[pymethods]
//impl PyFeynGenOptions {
//    pub(crate) fn __str__(&self) -> Result<String> {
//        Ok(format!("{}", self.options))
//    }
//    #[allow(clippy::too_many_arguments)]
//    #[new]
//    pub(crate) fn __new__(
//        generation_type: String,
//        initial_particles: Vec<i64>,
//        final_particles_lists: Vec<Vec<i64>>,
//        loop_count_range: (usize, usize),
//        n_cut_blobs: (usize, usize),
//        n_cut_spectators: (usize, usize),
//        symmetrize_initial_states: bool,
//        symmetrize_final_states: bool,
//        symmetrize_left_right_states: bool,
//        allow_symmetrization_of_external_fermions_in_amplitudes: bool,
//        max_multiplicity_for_fast_cut_filter: usize,
//        amplitude_filters: Option<PyRef<PyFeynGenFilters>>,
//        cross_section_filters: Option<PyRef<PyFeynGenFilters>>,
//    ) -> Result<PyFeynGenOptions> {
//        let amplitude_filters = FeynGenFilters(
//            amplitude_filters
//                .map(|f| f.filters.clone())
//                .unwrap_or_default(),
//        );
//
//        let mut cross_section_filters = FeynGenFilters(
//            cross_section_filters
//                .map(|f| f.filters.clone())
//                .unwrap_or_default(),
//        );
//
//        cross_section_filters
//            .0
//            .push(FeynGenFilter::BlobRange(n_cut_blobs.0..=n_cut_blobs.1));
//        cross_section_filters.0.push(FeynGenFilter::SpectatorRange(
//            n_cut_spectators.0..=n_cut_spectators.1,
//        ));
//
//        let feyngen_options = FeynGenOptions {
//            generation_type: GenerationType::from_str(&generation_type)
//                .map_err(feyngen_to_python_error)?,
//            initial_pdgs: initial_particles,
//            final_pdgs_lists: final_particles_lists,
//            loop_count_range,
//            symmetrize_initial_states,
//            symmetrize_final_states,
//            symmetrize_left_right_states,
//            allow_symmetrization_of_external_fermions_in_amplitudes,
//            max_multiplicity_for_fast_cut_filter,
//            amplitude_filters,
//            cross_section_filters,
//        };
//        if feyngen_options.generation_type == GenerationType::Amplitude {
//            if feyngen_options.final_pdgs_lists.len() > 1 {
//                return Err(eyre!(
//                    "Multiple set of final states are not allowed for amplitude generation",
//                ));
//            }
//        } else if feyngen_options.final_pdgs_lists.len() > 1
//            && feyngen_options
//                .final_pdgs_lists
//                .iter()
//                .any(|l| l.is_empty())
//        {
//            return Err(eyre!(
//                    "When specifying multiple set of final states options, each must contain at least one particle",
//                ));
//        }
//        Ok(PyFeynGenOptions {
//            options: feyngen_options,
//        })
//    }
//}

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

// TODO: Improve error broadcasting to Python so as to show rust backtrace
#[pymethods]
impl State {
    #[new]
    pub fn new_python(state_folder: PathBuf) -> Self {
        initialise().unwrap();

        let a =
            State::load(state_folder.clone(), None, Some("python".to_string())).expect(&format!(
                "Could not load or create state from path '{}'",
                state_folder.display()
            ));

        a
    }

    pub fn inspect<'py>(
        &mut self,
        py: Python<'py>,
        process_id: Option<usize>,
        integrand_name: Option<String>,
        point: Vec<f64>,
        use_f128: bool,
        force_radius: bool,
        momentum_space: bool,
        discrete_dim: Vec<usize>,
    ) -> Result<(Bound<'py, PyComplex>, Option<Bound<'py, PyFloat>>)> {
        let res = Inspect {
            process_id,
            integrand_name,
            point,
            use_f128,
            force_radius,
            momentum_space,
            discrete_dim,
        }
        .run(self)?;
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
        use_f128: bool,
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
            use_f128,
            momentum_space,
            points: points_rust,
            discrete_dims: discrete_dims_rust,
        };

        let (res, res_jac) = batched_inspect.run(self).unwrap();
        let res_map = res
            .mapv_into_any(|c| Complex64::new(c.re, c.im))
            .into_pyarray(py);

        let res_jac_map = res_jac.map(|r| r.into_pyarray(py));

        Ok((res_map, res_jac_map))
    }

    #[pyo3(signature = (path, process_name=None, process_id=None, integrand_name=None))]
    pub(crate) fn import_amplitude_python(
        &mut self,
        path: PathBuf,
        process_name: Option<String>,
        process_id: Option<usize>,
        integrand_name: Option<String>,
    ) -> Result<()> {
        self.import_graphs(path, process_name, process_id, integrand_name)
    }

    pub(crate) fn import_model_python(&mut self, file_path: PathBuf) -> Result<()> {
        self.import_model(file_path)
    }

    //    pub(crate) fn load_model_from_yaml_str(&mut self, yaml_str: &str) -> Result<()> {
    //        self.model = Model::from_str(String::from(yaml_str))?;
    //        Ok(())
    //    }

    pub(crate) fn generate_integrands_python(
        &mut self,
        global_settings: &GlobalSettings,
        runtime_default: &RuntimeSettings,
    ) -> Result<()> {
        self.generate_integrands(global_settings, runtime_default.into())
    }

    pub(crate) fn compile_integrands_python(
        &mut self,
        folder: PathBuf,
        override_existing: bool,
        global_settings: &GlobalSettings,
        process_id: Option<usize>,
        integrand_name: Option<String>,
    ) -> Result<()> {
        self.compile_integrands(
            folder,
            override_existing,
            global_settings,
            process_id,
            integrand_name,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn generate_diagrams(
        &mut self,
        //generation_options: PyRef<PyFeynGenOptions>,
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
            process_id,
            integrand_name,
            result_path,
            workspace_path: Some(workspace_path),
            n_cores: Some(n_cores),
            target,
            restart,
        }
        .run(
            self,
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
}
