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
    "Base exception for native FeynKit operations.\n\nExamples\n--------\nCatch any model, diagram, generation, CFF, tensor-reduction, or kinematics failure:\n\n>>> try:\n...     result = model.generate_diagrams(incoming, outgoing)\n... except fk.FeynkitError as error:\n...     print(error)"
);
define_exception!(
    ModelError,
    FeynkitError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError"),
    "Invalid particle-model data or model operation.\n\nExamples\n--------\nA missing particle is reported as a model error:\n\n>>> try:\n...     model.particle_by_pdg(999999)\n... except fk.ModelError:\n...     pass"
);
define_exception!(
    DiagramError,
    FeynkitError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError"),
    "Invalid Feynman-diagram topology, annotation, or model reference.\n\nExamples\n--------\nValidate imported diagrams before further physics operations:\n\n>>> try:\n...     diagram.validate()\n... except fk.DiagramError as error:\n...     print(error)"
);
define_exception!(
    GenerationError,
    FeynkitError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError"),
    "Invalid process configuration or Feynman-diagram generation failure.\n\nExamples\n--------\nProcess and topology failures share one public exception type:\n\n>>> try:\n...     result = model.generate_diagrams(incoming, outgoing, loops=1)\n... except fk.GenerationError as error:\n...     print(error)"
);
define_exception!(
    CffError,
    FeynkitError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError"),
    "Failure while constructing a Cross-Free Family representation.\n\nExamples\n--------\nCatch incompatible edge constraints at the CFF boundary:\n\n>>> try:\n...     cff = diagram.build_cff(contracted_edges=[edge_id])\n... except fk.CffError as error:\n...     print(error)"
);
define_exception!(
    KinematicsError,
    FeynkitError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError"),
    "Invalid Lorentz transformation, momentum, or jet-clustering request.\n\nExamples\n--------\nKinematic-domain failures remain distinct from model errors:\n\n>>> try:\n...     jets = fk.JetDefinition.anti_kt(-0.4).cluster(momenta)\n... except fk.KinematicsError as error:\n...     print(error)"
);
define_exception!(
    TensorReductionError,
    FeynkitError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError"),
    "Failure while parsing or reducing a Lorentz tensor.\n\nExamples\n--------\nCatch unsupported tensor structures or an exceeded expansion budget:\n\n>>> try:\n...     scalar = reducer.reduce(numerator)\n... except fk.TensorReductionError as error:\n...     print(error)"
);
#[cfg(feature = "ufo")]
define_exception!(
    UfoLoadError,
    FeynkitError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError"),
    "Failure while importing or normalizing a UFO model.\n\nExamples\n--------\nReport a missing or malformed UFO directory cleanly:\n\n>>> try:\n...     loaded = fk.UfoLoader().load(\"models/sm\")\n... except fk.UfoLoadError as error:\n...     print(error)"
);

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
    module.add("ModelError", py.get_type::<ModelError>())?;
    module.add("DiagramError", py.get_type::<DiagramError>())?;
    module.add("GenerationError", py.get_type::<GenerationError>())?;
    module.add("CffError", py.get_type::<CffError>())?;
    module.add("KinematicsError", py.get_type::<KinematicsError>())?;
    module.add(
        "TensorReductionError",
        py.get_type::<TensorReductionError>(),
    )?;
    #[cfg(feature = "ufo")]
    module.add("UfoLoadError", py.get_type::<UfoLoadError>())?;
    Ok(())
}
