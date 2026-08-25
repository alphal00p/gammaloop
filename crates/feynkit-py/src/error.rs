use pyo3::{PyErr, create_exception, exceptions::PyException, types::PyModuleMethods};

create_exception!(symbolica.community.feynkit, FeynkitError, PyException);
create_exception!(symbolica.community.feynkit, ModelError, FeynkitError);
create_exception!(symbolica.community.feynkit, DiagramError, FeynkitError);
create_exception!(symbolica.community.feynkit, GenerationError, FeynkitError);
create_exception!(symbolica.community.feynkit, CffError, FeynkitError);
create_exception!(symbolica.community.feynkit, KinematicsError, FeynkitError);
#[cfg(feature = "ufo")]
create_exception!(symbolica.community.feynkit, UfoLoadError, FeynkitError);

#[cfg(feature = "python_stubgen")]
macro_rules! register_exception_stub {
    ($name:ident, $base:expr) => {
        pyo3_stub_gen::inventory::submit! {
            pyo3_stub_gen::type_info::PyClassInfo {
                pyclass_name: stringify!($name),
                struct_id: std::any::TypeId::of::<$name>,
                getters: &[],
                setters: &[],
                module: Some("symbolica.community.feynkit"),
                doc: "",
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

#[cfg(feature = "python_stubgen")]
register_exception_stub!(FeynkitError, pyo3_stub_gen::TypeInfo::builtin("Exception"));
#[cfg(feature = "python_stubgen")]
register_exception_stub!(
    ModelError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError")
);
#[cfg(feature = "python_stubgen")]
register_exception_stub!(
    DiagramError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError")
);
#[cfg(feature = "python_stubgen")]
register_exception_stub!(
    GenerationError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError")
);
#[cfg(feature = "python_stubgen")]
register_exception_stub!(
    CffError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError")
);
#[cfg(feature = "python_stubgen")]
register_exception_stub!(
    KinematicsError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError")
);
#[cfg(all(feature = "python_stubgen", feature = "ufo"))]
register_exception_stub!(
    UfoLoadError,
    pyo3_stub_gen::TypeInfo::unqualified("FeynkitError")
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
    #[cfg(feature = "ufo")]
    module.add("UfoLoadError", py.get_type::<UfoLoadError>())?;
    Ok(())
}
