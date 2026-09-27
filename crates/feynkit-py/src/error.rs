use pyo3::{PyErr, create_exception, exceptions::PyException, types::PyModuleMethods};

#[cfg(feature = "python_stubgen")]
macro_rules! register_exception_docs {
    ($name:ident, $base:expr, $doc:expr) => {
        pyo3_stub_gen::inventory::submit! {
            pyo3_stub_gen::type_info::PyClassInfo {
                pyclass_name: stringify!($name),
                struct_id: std::any::TypeId::of::<$name>,
                getters: &[],
                setters: &[],
                module: Some("symbolica.community.feynkit"),
                doc: $doc,
                bases: &[|| $base],
                has_eq: false,
                has_ord: false,
                has_hash: false,
                has_str: false,
                subclass: true,
            }
        }
    };
}

macro_rules! define_exception {
    ($name:ident, $base:ty, $stub_base:expr, $doc:literal) => {
        create_exception!(symbolica.community.feynkit, $name, $base, $doc);
        #[cfg(feature = "python_stubgen")]
        register_exception_docs!($name, $stub_base, $doc);
    };
}

define_exception!(
    FeynkitError,
    PyException,
    pyo3_stub_gen::TypeInfo::builtin("Exception"),
    "Base exception for native HEP operations.\n\nExamples\n--------\n>>> from symbolica import S, E\n>>> from symbolica.community import hep\n>>> model = hep.Model.standard_model()\n>>> try:\n...     model.particle_by_pdg(999999)\n... except hep.FeynkitError as error:\n...     message = str(error)"
);
define_exception!(
    ModelError,
    FeynkitError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError"),
    "Invalid particle-model data or model operation.\n\nExamples\n--------\n>>> from symbolica import S, E\n>>> from symbolica.community import hep\n>>> model = hep.Model.standard_model()\n>>> try:\n...     model.particle_by_pdg(999999)\n... except hep.ModelError as error:\n...     message = str(error)"
);
define_exception!(
    DiagramError,
    FeynkitError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError"),
    "Invalid Feynman-diagram topology, annotation, or model reference.\n\nExamples\n--------\n>>> from symbolica import S, E\n>>> from symbolica.community import hep\n>>> model = hep.Model.phi4()\n>>> process = model.process([\"phi\", \"phi\"], [\"phi\", \"phi\"])\n>>> result = process.generate_diagrams(loops=1)\n>>> diagram = result.diagrams[0]\n>>> try:\n...     diagram.with_loop_momentum_edges([])\n... except hep.DiagramError as error:\n...     message = str(error)"
);
define_exception!(
    GenerationError,
    FeynkitError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError"),
    "Invalid process configuration or Feynman-diagram generation failure.\n\nExamples\n--------\n>>> from symbolica import S, E\n>>> from symbolica.community import hep\n>>> model = hep.Model.phi4()\n>>> process = model.process([\"phi\", \"phi\"], [\"phi\", \"phi\"])\n>>> try:\n...     process.generate_diagrams(self_energy=hep.SelfEnergyFilterOptions(only_scaleless=True))\n... except hep.GenerationError as error:\n...     message = str(error)"
);
define_exception!(
    CffError,
    FeynkitError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError"),
    "Failure while constructing a Cross-Free Family representation.\n\nExamples\n--------\n>>> from symbolica import S, E\n>>> from symbolica.community import hep\n>>> model = hep.Model.phi4()\n>>> process = model.process([\"phi\", \"phi\"], [\"phi\", \"phi\"])\n>>> result = process.generate_diagrams(loops=1)\n>>> diagram = result.diagrams[0]\n>>> try:\n...     diagram.build_cff(contracted_edges=[999999])\n... except hep.CffError as error:\n...     message = str(error)"
);
define_exception!(
    KinematicsError,
    FeynkitError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError"),
    "Invalid Lorentz transformation, momentum, or jet-clustering request.\n\nExamples\n--------\n>>> from symbolica import S, E\n>>> from symbolica.community import hep\n>>> try:\n...     hep.JetDefinition.anti_kt(-0.4)\n... except hep.KinematicsError as error:\n...     message = str(error)"
);
define_exception!(
    TensorReductionError,
    FeynkitError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError"),
    "Failure while parsing or reducing a Lorentz tensor.\n\nExamples\n--------\n>>> from symbolica import S, E\n>>> from symbolica.community import hep\n>>> D, k, mu = S(\"D\", \"k\", \"mu\")\n>>> mink = S(\"spenso::mink\")\n>>> reducer = hep.TensorReducer(D).with_integrated_vector(k(mink(D)))\n>>> try:\n...     scalar = reducer.reduce(k(mink(D, mu))**2)\n... except hep.TensorReductionError as error:\n...     message = str(error)"
);
#[cfg(feature = "ufo")]
define_exception!(
    UfoLoadError,
    FeynkitError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError"),
    "Failure while importing or normalizing a UFO model.\n\nExamples\n--------\n>>> from symbolica import S, E\n>>> from symbolica.community import hep\n>>> try:\n...     hep.UfoLoader().load(\"/path/to/missing-model\")\n... except hep.UfoLoadError as error:\n...     message = str(error)"
);

define_exception!(
    AmplitudeError,
    FeynkitError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError"),
    "Invalid amplitude, conjugation, or external-state sum.\n\nExamples\n--------\n>>> from symbolica import S, E\n>>> from symbolica.community import hep\n>>> try:\n...     hep.Amplitude([])\n... except hep.AmplitudeError as error:\n...     message = str(error)"
);

pub(crate) fn amplitude(error: feynkit_amplitude::AmplitudeError) -> PyErr {
    AmplitudeError::new_err(error.to_string())
}

pub(crate) fn model(error: feynkit_model::ModelError) -> PyErr {
    ModelError::new_err(error.to_string())
}

pub(crate) fn recompute(error: feynkit_model::RecomputeError<PyErr>) -> PyErr {
    match error {
        feynkit_model::RecomputeError::Model(error) => model(error),
        feynkit_model::RecomputeError::Evaluator(error) => error,
        error => ModelError::new_err(error.to_string()),
    }
}

pub(crate) fn diagram(error: feynkit_graph::DiagramError) -> PyErr {
    DiagramError::new_err(error.to_string())
}

pub(crate) fn generation(error: feynkit_generator::GenerationError) -> PyErr {
    GenerationError::new_err(error.to_string())
}

pub(crate) fn process(error: feynkit_generator::ProcessError) -> PyErr {
    GenerationError::new_err(error.to_string())
}

pub(crate) fn cff(error: feynkit_cff::CffError) -> PyErr {
    CffError::new_err(error.to_string())
}

pub(crate) fn kinematics(error: impl std::fmt::Display) -> PyErr {
    KinematicsError::new_err(error.to_string())
}

define_exception!(
    IntegralFamilyError,
    FeynkitError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError"),
    "Invalid loop-integral family, dependent propagators, or incomplete scalar-product basis.\n\nExamples\n--------\n>>> from symbolica import S, E\n>>> from symbolica.community import hep\n>>> D, k, p, s = S(\"D\", \"k\", \"p\", \"s\")\n>>> d1, d2, x1, x2 = S(\"d1\", \"d2\", \"x1\", \"x2\")\n>>> kin = hep.Kinematics(D, momenta=[k, p]).with_scalar_product(p, p, s)\n>>> denominators = [kin.scalar_product(k, k), kin.scalar_product(k-p, k-p)]\n>>> family = hep.IntegralFamily([k], [p], denominators, kinematics=kin)\n>>> try:\n...     hep.IntegralFamily([k], [p], [denominators[0], denominators[0]], kinematics=kin).complete()\n... except hep.IntegralFamilyError as error:\n...     message = str(error)"
);

pub(crate) fn integral_family(error: feynkit_graph::IntegralFamilyError) -> PyErr {
    IntegralFamilyError::new_err(error.to_string())
}

pub(crate) fn tensor(error: feynkit_tensor::TensorReductionError) -> PyErr {
    TensorReductionError::new_err(error.to_string())
}

#[cfg(feature = "ufo")]
pub(crate) fn ufo(error: feynkit_ufo::UfoLoadError) -> PyErr {
    UfoLoadError::new_err(error.to_string())
}

pub(crate) fn register(module: &pyo3::Bound<'_, pyo3::types::PyModule>) -> pyo3::PyResult<()> {
    let py = module.py();
    module.add("FeynkitError", py.get_type::<FeynkitError>())?;
    module.add("AmplitudeError", py.get_type::<AmplitudeError>())?;
    module.add("ModelError", py.get_type::<ModelError>())?;
    module.add("DiagramError", py.get_type::<DiagramError>())?;
    module.add("GenerationError", py.get_type::<GenerationError>())?;
    module.add("CffError", py.get_type::<CffError>())?;
    module.add("KinematicsError", py.get_type::<KinematicsError>())?;
    module.add("IntegralFamilyError", py.get_type::<IntegralFamilyError>())?;
    module.add(
        "TensorReductionError",
        py.get_type::<TensorReductionError>(),
    )?;
    #[cfg(feature = "ufo")]
    module.add("UfoLoadError", py.get_type::<UfoLoadError>())?;
    Ok(())
}
