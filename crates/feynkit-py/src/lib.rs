//! Python bindings installed as `symbolica.community.feynkit`.

mod cff;
mod error;
mod generation;
mod graph;
mod kinematics;
mod model;
#[cfg(feature = "ufo")]
mod ufo;

use pyo3::{prelude::*, types::PyModule};
use symbolica::api::python::SymbolicaCommunityModule;

pub use cff::{PyCffGenerator, PyCffOrientation, PyCffReport, PyCffResult, PyCffSurface};
pub use generation::{
    PyCancellationToken, PyDiagramGroup, PyGenerationOptions, PyGenerationReport,
    PyGenerationResult, PyGenerationType, PyGenerator, PyGroupMember, PyParticleSelector,
    PyProcess,
};
pub use graph::{
    PyDiagramEdge, PyDiagramVertex, PyFeynmanDiagram, PyLoopMomentumBasis, PyMomentumSignature,
};
pub use kinematics::{
    PyAxis, PyBoost, PyClusteringResult, PyFourMomentum, PyHelicity, PyJet, PyJetAlgorithm,
    PyJetDefinition, PyRotation, PyThreeMomentum,
};
pub use model::{
    PyCoupling, PyEvaluatedValues, PyEvaluationRequest, PyFormFactor, PyLorentzStructure, PyModel,
    PyModelExpression, PyModelFunction, PyParameter, PyParameterCard, PyParameterNature,
    PyParameterType, PyParticle, PyPropagator, PyVertexRule,
};
#[cfg(feature = "ufo")]
pub use ufo::{PyLoadedModel, PyUfoLoadDiagnostics, PyUfoLoader};

pub struct FeynkitModule;

impl SymbolicaCommunityModule for FeynkitModule {
    fn get_name() -> String {
        "feynkit".to_owned()
    }

    fn register_module(module: &Bound<'_, PyModule>) -> PyResult<()> {
        initialize_feynkit(module)
    }

    fn initialize(_py: Python<'_>) -> PyResult<()> {
        Ok(())
    }
}

/// Register FeynKit classes in an existing Symbolica community module.
pub fn initialize_feynkit(module: &Bound<'_, PyModule>) -> PyResult<()> {
    error::register(module)?;
    model::register(module)?;
    graph::register(module)?;
    generation::register(module)?;
    kinematics::register(module)?;
    cff::register(module)?;
    #[cfg(feature = "ufo")]
    ufo::register(module)?;
    Ok(())
}

/// Gather the FeynKit stub inventory without declaring an independent wheel.
#[cfg(feature = "python_stubgen")]
pub fn stub_info() -> pyo3_stub_gen::Result<pyo3_stub_gen::StubInfo> {
    pyo3_stub_gen::StubInfo::from_project_root(
        "symbolica.community.feynkit".to_owned(),
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("python"),
    )
}

#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use pyo3::types::PyDict;
    use symbolica::api::python::PythonExpression;

    use super::*;

    const NORMALIZED_SCALAR_MODEL: &str = include_str!("../tests/fixtures/scalars_2p_3p.json");

    fn registered_module<'py>(py: Python<'py>) -> Bound<'py, PyModule> {
        let module = PyModule::new(py, "symbolica.community.feynkit").unwrap();
        FeynkitModule::register_module(&module).unwrap();
        FeynkitModule::initialize(py).unwrap();
        module
    }

    #[test]
    fn advertises_the_community_module_name() {
        assert_eq!(FeynkitModule::get_name(), "feynkit");
    }

    #[test]
    fn primary_wrappers_are_send() {
        fn assert_send<T: Send>() {}
        assert_send::<PyModel>();
        assert_send::<PyGenerator>();
        assert_send::<PyFeynmanDiagram>();
    }

    #[test]
    fn registers_as_a_standalone_community_submodule() {
        Python::initialize();
        Python::attach(|py| {
            let module = registered_module(py);
            for class in [
                "Model",
                "Generator",
                "Process",
                "FeynmanDiagram",
                "CffGenerator",
                "FourMomentum",
                "JetDefinition",
                "Helicity",
            ] {
                assert!(
                    module.getattr(class).is_ok(),
                    "missing Python class {class}"
                );
            }
        });
    }

    #[test]
    fn standalone_python_pipeline_covers_generation_graphs_cff_and_kinematics() {
        Python::initialize();
        Python::attach(|py| {
            let module = registered_module(py);
            let locals = PyDict::new(py);
            locals.set_item("fk", &module).unwrap();
            locals
                .set_item("MODEL_JSON", NORMALIZED_SCALAR_MODEL)
                .unwrap();
            let code = CString::new(
                r#"
import sys

assert "_gammaloop" not in sys.modules
assert fk.__name__ == "symbolica.community.feynkit"

model = fk.Model.from_json(MODEL_JSON)
assert model.name == "scalars"
assert model.particle("scalar_0").pdg_code == 1000
model = fk.Model.from_json(model.to_json(pretty=False))

process = fk.Process.amplitude(["scalar_0"], [1000, "scalar_0"])
process = process.with_loop_count(0, 1)
assert isinstance(process, fk.Process)
assert process.generation_type == "amplitude"
assert process.loop_count == (0, 1)

options = fk.GenerationOptions(max_vertices=3, allow_self_loops=True)
options.add_vertex_allow(["V_3_SCALAR_000"])
options.set_coupling_orders({"QCD": (0, None)})
options.set_loop_count_range(0, 1)
options.set_fermion_loop_count_range(0, 0)
options.set_factorized_loop_topologies_count_range(0, 1)
assert options.disable_numerator_grouping() is None
generated = fk.Generator(model).generate(process, options)
assert len(generated) > 0
assert generated.report.completed

loop_diagram = next(
    diagram
    for diagram in generated.diagrams
    if diagram.loop_count == 1
    and all(edge.source != edge.target for edge in diagram.edges)
)
loop_diagram.validate(model)

json_diagram = fk.FeynmanDiagram.from_json(loop_diagram.to_json())
dot_diagram = fk.FeynmanDiagram.from_dot(loop_diagram.to_dot())
json_diagram.validate(model)
dot_diagram.validate(model)
assert json_diagram.loop_count == dot_diagram.loop_count == 1

bases = json_diagram.loop_momentum_bases()
assert bases
assert len(bases[0].loop_edges) == 1
assert len(bases[0].edge_signatures) == len(json_diagram.edges)

cff = fk.CffGenerator().generate(json_diagram)
assert len(cff) > 0
expression = cff.to_expression()

momentum = fk.ThreeMomentum(3.0, 4.0, 0.0)
assert momentum.on_shell().components() == (5.0, 3.0, 4.0, 0.0)
rotated = fk.Rotation.quarter_turn(fk.Axis.Z).apply_three(
    fk.ThreeMomentum(1.0, 0.0, 0.0)
)
assert abs(rotated.px) < 1.0e-12
assert abs(rotated.py - 1.0) < 1.0e-12

clustered = fk.JetDefinition.anti_kt(0.4).cluster([
    fk.FourMomentum(10.0, 10.0, 0.0, 0.0),
    fk.FourMomentum(5.0, 5.0, 0.0, 0.0),
])
assert len(clustered) == 1
assert clustered.jets[0].constituent_indices == [0, 1]
assert clustered.jets[0].momentum.components() == (15.0, 15.0, 0.0, 0.0)
assert "_gammaloop" not in sys.modules
"#,
            )
            .unwrap();

            py.run(&code, Some(&locals), Some(&locals)).unwrap();
            let expression = locals.get_item("expression").unwrap().unwrap();
            assert!(expression.is_instance_of::<PythonExpression>());
        });
    }

    #[test]
    fn standalone_python_errors_keep_their_public_types_and_module_isolation() {
        Python::initialize();
        Python::attach(|py| {
            let module = registered_module(py);
            let locals = PyDict::new(py);
            locals.set_item("fk", &module).unwrap();
            let code = CString::new(
                r#"
import sys

assert "_gammaloop" not in sys.modules

def assert_feynkit_error(error_type, operation):
    try:
        operation()
    except error_type as error:
        assert isinstance(error, fk.FeynkitError)
        assert type(error).__module__ == "symbolica.community.feynkit"
    else:
        raise AssertionError(f"{error_type.__name__} was not raised")

assert_feynkit_error(fk.ModelError, lambda: fk.Model.from_json("{}"))
assert_feynkit_error(
    fk.GenerationError,
    lambda: fk.Process.amplitude([], []).with_loop_count(2, 1),
)
assert_feynkit_error(fk.DiagramError, lambda: fk.FeynmanDiagram.from_json("{}"))
assert_feynkit_error(
    fk.KinematicsError,
    lambda: fk.Boost(fk.ThreeMomentum(1.0, 0.0, 0.0)),
)

model = fk.Model.from_json(MODEL_JSON)
generator = fk.Generator(model)

filter_options = fk.GenerationOptions()
assert filter_options.set_self_energy_filter(
    veto_massive=False,
    veto_massless=True,
    only_scaleless=False,
) is None
assert filter_options.set_tadpole_filter(
    veto_attached_to_massive=False,
    veto_attached_to_massless=True,
    only_scaleless=False,
) is None
assert filter_options.set_zero_snail_filter(
    veto_attached_to_massive=True,
    veto_attached_to_massless=False,
    only_scaleless=False,
) is None
assert filter_options.set_coupling_orders({"QCD": (0, None)}) is None
assert filter_options.set_loop_count_range(0, 2) is None
assert filter_options.set_fermion_loop_count_range(0, 1) is None
assert filter_options.set_factorized_loop_topologies_count_range(0, 1) is None
assert filter_options.set_blob_range(0, 1) is None
assert filter_options.set_spectator_range(0, 1) is None
assert filter_options.set_cut_amplitude_coupling_orders(
    {"QCD": (0, None)}
) is None
assert filter_options.set_cut_amplitude_loop_count_range(0, 1) is None

grouping_options = fk.GenerationOptions()
assert grouping_options.disable_numerator_grouping() is None
assert grouping_options.detect_zero_numerators() is None
grouping_arguments = dict(
    numerical_sample_seed=7,
    number_of_numerical_samples=11,
    differentiate_particle_masses_only=False,
    fully_numerical_substitution=True,
    check_canonical_numerator=True,
)
assert grouping_options.group_identical_numerators(**grouping_arguments) is None
assert grouping_options.group_numerators_up_to_sign(**grouping_arguments) is None
assert grouping_options.group_numerators_up_to_scalar(**grouping_arguments) is None
assert not hasattr(fk.GenerationOptions, "set_numerator_grouping")

amplitude = fk.Process.amplitude(
    ["scalar_0"],
    ["scalar_0", "scalar_0"],
)

def assert_generation_configuration_error(configure):
    invalid_options = fk.GenerationOptions(max_vertices=3)
    configure(invalid_options)
    assert_feynkit_error(
        fk.GenerationError,
        lambda: generator.generate(amplitude, invalid_options),
    )

assert_generation_configuration_error(
    lambda value: value.set_self_energy_filter(only_scaleless=True)
)
assert_generation_configuration_error(
    lambda value: value.set_tadpole_filter(only_scaleless=True)
)
assert_generation_configuration_error(
    lambda value: value.set_zero_snail_filter(only_scaleless=True)
)
assert_generation_configuration_error(
    lambda value: value.set_coupling_orders({"QCD": (2, 1)})
)
assert_generation_configuration_error(
    lambda value: value.set_loop_count_range(2, 1)
)
assert_generation_configuration_error(
    lambda value: value.set_fermion_loop_count_range(2, 1)
)
assert_generation_configuration_error(
    lambda value: value.set_factorized_loop_topologies_count_range(2, 1)
)
assert_generation_configuration_error(
    lambda value: value.set_blob_range(0, 1)
)
assert_generation_configuration_error(
    lambda value: value.set_spectator_range(0, 1)
)
assert_generation_configuration_error(
    lambda value: value.set_cut_amplitude_coupling_orders({"QCD": (0, None)})
)
assert_generation_configuration_error(
    lambda value: value.set_cut_amplitude_loop_count_range(0, 1)
)

valid_grouping = fk.GenerationOptions(max_vertices=3)
valid_grouping.group_identical_numerators(**grouping_arguments)
assert generator.generate(amplitude, valid_grouping).report.completed

options = fk.GenerationOptions(max_vertices=3)
options.add_vertex_allow(["V_3_SCALAR_000"])
diagram = generator.generate(
    fk.Process.amplitude(["scalar_0"], ["scalar_0", "scalar_0"])
    .with_loop_count(1, 1),
    options,
).diagrams[0]
assert_feynkit_error(
    fk.CffError,
    lambda: fk.CffGenerator(max_orientations=0).generate(diagram),
)

assert "_gammaloop" not in sys.modules
"#,
            )
            .unwrap();
            locals
                .set_item("MODEL_JSON", NORMALIZED_SCALAR_MODEL)
                .unwrap();

            py.run(&code, Some(&locals), Some(&locals)).unwrap();
        });
    }
}
