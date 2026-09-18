//! Python bindings installed as `symbolica.community.feynkit`.

mod cff;
mod display;
mod error;
mod generation;
mod graph;
mod kinematics;
mod model;
mod tensor;
#[cfg(feature = "ufo")]
mod ufo;

use pyo3::{prelude::*, types::PyModule};
use symbolica::api::python::SymbolicaCommunityModule;

pub use cff::{PyCffGenerator, PyCffOrientation, PyCffReport, PyCffResult, PyCffSurface};
pub use generation::{
    PyCancellationToken, PyDiagramGroup, PyGenerationProgress, PyGenerationReport,
    PyGenerationResult, PyGenerationType, PyGenerator, PyGroupMember, PyNumeratorGrouping,
    PyParticleSelector, PyProcess, PySelfEnergyFilterOptions, PySnailFilterOptions,
    PyTadpoleFilterOptions,
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
pub use tensor::PyTensorReducer;
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
    tensor::register(module)?;
    #[cfg(feature = "ufo")]
    ufo::register(module)?;
    Ok(())
}

/// Gather the FeynKit stub inventory without declaring an independent wheel.
#[cfg(feature = "python_stubgen")]
pub fn stub_info() -> pyo3_stub_gen::Result<pyo3_stub_gen::StubInfo> {
    // Unqualified dependency classes belong to Symbolica; FeynKit declarations
    // explicitly name their community module.
    let info = pyo3_stub_gen::StubInfo::from_project_root(
        "symbolica".to_owned(),
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("python"),
    )?;
    let module = info
        .modules
        .get("symbolica.community.feynkit")
        .ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "FeynKit did not contribute a symbolica.community.feynkit stub module",
            )
        })?;
    validate_stub_documentation(&module.to_string())?;
    Ok(info)
}

#[cfg(feature = "python_stubgen")]
fn validate_stub_documentation(source: &str) -> PyResult<()> {
    const DOCUMENTATION_AUDIT: &str = r#"
import ast

tree = ast.parse(source)
errors = []

def section_body(doc, heading):
    lines = doc.splitlines()
    try:
        start = lines.index(heading) + 2
    except ValueError:
        return None
    end = len(lines)
    for index in range(start, len(lines) - 1):
        underline = lines[index + 1].strip()
        if lines[index].strip() and len(underline) >= 3 and set(underline) == {"-"}:
            end = index
            break
    return [line.strip() for line in lines[start:end] if line.strip()]

for node in ast.walk(tree):
    if isinstance(node, ast.ClassDef):
        doc = ast.get_docstring(node) or ""
        if "Examples\n--------" not in doc:
            errors.append(f"class {node.name} has no Examples section")
        if "Examples\n--------" in doc and "Parameters\n----------" in doc:
            if doc.index("Examples\n--------") > doc.index("Parameters\n----------"):
                errors.append(f"class {node.name} puts Parameters before Examples")
        parameters = section_body(doc, "Parameters")
        if parameters is not None and (not parameters or parameters == ["None"]):
            errors.append(f"class {node.name} has an empty Parameters section")
        if "isinstance(" in doc or "is None or" in doc:
            errors.append(f"class {node.name} uses a type-check-only example")
        continue
    if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
        continue

    doc = ast.get_docstring(node) or ""
    is_property = any(
        isinstance(decorator, ast.Name) and decorator.id == "property"
        for decorator in node.decorator_list
    )
    arguments = [
        *node.args.posonlyargs,
        *node.args.args,
        *node.args.kwonlyargs,
    ]
    arguments = [argument for argument in arguments if argument.arg not in {"self", "cls"}]
    if node.args.vararg is not None:
        arguments.append(node.args.vararg)
    if node.args.kwarg is not None:
        arguments.append(node.args.kwarg)
    has_examples = "Examples\n--------" in doc
    has_parameters = "Parameters\n----------" in doc
    parameters = section_body(doc, "Parameters")
    documented_parameters = set()
    for line in parameters or []:
        if ":" not in line:
            continue
        names = line.split(":", 1)[0]
        documented_parameters.update(
            name.strip().lstrip("*") for name in names.split(",")
        )

    if not doc:
        errors.append(f"callable {node.name} has no documentation")
    if not is_property and not has_examples:
        errors.append(f"callable {node.name} has no Examples section")
    if not is_property and arguments and not has_parameters:
        errors.append(f"callable {node.name} has undocumented parameters")
    elif arguments:
        missing = [
            argument.arg for argument in arguments
            if argument.arg not in documented_parameters
        ]
        if missing:
            errors.append(
                f"callable {node.name} omits parameters: {', '.join(missing)}"
            )
    if not arguments and has_parameters:
        errors.append(f"callable {node.name} has an empty Parameters section")
    if parameters is not None and (not parameters or parameters == ["None"]):
        errors.append(f"callable {node.name} has an empty Parameters section")
    if has_examples and has_parameters and doc.index("Examples\n--------") > doc.index("Parameters\n----------"):
        errors.append(f"callable {node.name} puts Parameters before Examples")
    if "isinstance(" in doc or "is None or" in doc:
        errors.append(f"callable {node.name} uses a type-check-only example")

if errors:
    raise AssertionError("invalid FeynKit API documentation:\n" + "\n".join(errors))
"#;

    Python::initialize();
    Python::attach(|py| {
        PyModule::import(py, "ast")?
            .getattr("parse")?
            .call1((source,))?;
        let globals = pyo3::types::PyDict::new(py);
        globals.set_item("source", source)?;
        PyModule::import(py, "builtins")?.getattr("exec")?.call1((
            DOCUMENTATION_AUDIT,
            &globals,
            &globals,
        ))?;
        Ok(())
    })
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
        assert_send::<PyTensorReducer>();
    }

    #[cfg(feature = "python_stubgen")]
    #[test]
    fn generated_stub_covers_every_registered_native_class() {
        Python::initialize();
        Python::attach(|py| {
            let module = registered_module(py);
            let info = stub_info().unwrap();
            let source = info
                .modules
                .get("symbolica.community.feynkit")
                .unwrap()
                .to_string();
            let locals = PyDict::new(py);
            locals.set_item("fk", module).unwrap();
            locals.set_item("source", source).unwrap();
            let code = CString::new(
                r#"
import ast

stub_classes = {
    node.name for node in ast.walk(ast.parse(source))
    if isinstance(node, ast.ClassDef)
}
runtime_classes = {
    name for name, value in vars(fk).items()
    if isinstance(value, type) and value.__module__ == fk.__name__
}
missing = sorted(runtime_classes - stub_classes)
assert not missing, f"native classes missing from the generated stub: {missing}"
"#,
            )
            .unwrap();
            py.run(&code, Some(&locals), Some(&locals)).unwrap();
        });
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
                "TensorReducer",
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
            for removed_function in [
                "generate_diagrams",
                "build_cff",
                "cluster_jets",
                "load_ufo_model",
            ] {
                assert!(
                    module.getattr(removed_function).is_err(),
                    "redundant module-level function {removed_function} must not be exported"
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
scalar = model.particle("scalar_0")
assert scalar.pdg_code == 1000
assert scalar.antiparticle.name == "scalar_0"
model = fk.Model.from_json(model.to_json(pretty=False))

scalar = model.particle("scalar_0")
process = fk.Process.amplitude([scalar], [1000, scalar.antiparticle])
process = process.with_loop_count(0, 1)
assert isinstance(process, fk.Process)
assert process.generation_type == "amplitude"
assert process.loop_count == (0, 1)

generation_arguments = dict(
    max_vertices=3,
    allow_self_loops=True,
    vertex_allow=["V_3_SCALAR_000"],
    particle_veto=[model.particle("scalar_1"), 1002],
    coupling_orders={"QCD": (0, None)},
    fermion_loop_count_range=(0, 0),
    factorized_loop_topologies_count_range=(0, 1),
)
generated = model.generate_diagrams(
    [scalar], [1000, scalar.antiparticle], loops=(0, 1), **generation_arguments,
)
assert len(generated) > 0
assert generated[0].name == next(iter(generated)).name
assert generated.report.completed
try:
    generation_html = generated._repr_html_()
except ImportError as error:
    assert "typst-py" in str(error)
else:
    assert "Generation result" in generation_html
assert "retained diagrams" in generated.report._repr_html_()
assert "particles" in model._repr_html_()

loop_diagram = next(
    diagram
    for diagram in generated.diagrams
    if diagram.loop_count == 1
    and all(edge.source != edge.target for edge in diagram.edges)
)
loop_diagram.validate()
assert loop_diagram.superficial_degree_of_divergence() == -2
assert loop_diagram.superficial_degree_of_divergence(dimension=6) == 0
denominator = loop_diagram.denominator_expression()
assert denominator.rank == 0
assert "ZERO" not in str(denominator)
assert denominator != 1
integrand = loop_diagram.numerator_expression() / denominator
try:
    diagram_svg = loop_diagram.to_svg()
except ImportError as error:
    # The minimal Rust test environment does not install optional Python
    # rendering dependencies; the community-venv tests exercise the SVG path.
    assert "typst-py" in str(error)
else:
    assert diagram_svg.startswith("<svg")
    assert loop_diagram.to_html() == loop_diagram._repr_html_()
    assert "<svg" in loop_diagram._repr_svg_()
    assert "Feynman diagram" in loop_diagram._repr_html_()

json_diagram = fk.FeynmanDiagram.from_json(model, loop_diagram.to_json())
dot_diagram = fk.FeynmanDiagram.from_dot(model, loop_diagram.to_dot())
json_diagram.validate()
dot_diagram.validate()
assert json_diagram.loop_count == dot_diagram.loop_count == 1

bases = json_diagram.loop_momentum_bases()
assert bases
assert len(bases[0].loop_edges) == 1
assert len(bases[0].edge_signatures) == len(json_diagram.edges)
assert "Loop-momentum basis" in bases[0]._repr_html_()

external_spatial = [
    fk.ThreeMomentum(0.0, 0.0, 10.0),
    fk.ThreeMomentum(0.0, 0.0, 4.0),
    fk.ThreeMomentum(0.0, 0.0, 6.0),
]
routed = bases[0].route(
    [fk.ThreeMomentum(1.0, 2.0, 3.0)],
    external_spatial,
)
assert set(routed) == {edge.id for edge in json_diagram.edges}

cff = json_diagram.build_cff()
assert len(cff) > 0
expression = cff.to_expression()
assert "Cross-free family" in cff._repr_html_()
assert "candidate orientations" in cff.report._repr_html_()

momentum = fk.ThreeMomentum(3.0, 4.0, 0.0)
assert momentum.on_shell().components() == (5.0, 3.0, 4.0, 0.0)
assert r"\vec{p}" in momentum._repr_latex_()
assert r"p^\mu" in momentum.on_shell()._repr_latex_()
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
assert clustered[0].constituent_indices == [0, 1]
assert next(iter(clustered)).momentum.components() == (15.0, 15.0, 0.0, 0.0)
assert "Clustered jets" in clustered._repr_html_()
assert "constituents" in clustered.jets[0]._repr_html_()
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

model = fk.Model.from_json(MODEL_JSON)

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
assert_feynkit_error(
    fk.DiagramError,
    lambda: fk.FeynmanDiagram.from_json(model, "{}"),
)
assert_feynkit_error(
    fk.KinematicsError,
    lambda: fk.Boost(fk.ThreeMomentum(1.0, 0.0, 0.0)),
)

generator = fk.Generator(model)

filter_arguments = dict(
    self_energy=fk.SelfEnergyFilterOptions(veto_massive=False),
    tadpoles=fk.TadpoleFilterOptions(veto_attached_to_massive=False),
    zero_snails=fk.SnailFilterOptions(
        veto_attached_to_massive=True, veto_attached_to_massless=False,
    ),
    coupling_orders={"QCD": (0, None)},
    fermion_loop_count_range=(0, 1),
    factorized_loop_topologies_count_range=(0, 1),
)
grouping_arguments = dict(
    numerical_sample_seed=7,
    number_of_numerical_samples=11,
    differentiate_particle_masses_only=False,
    fully_numerical_substitution=True,
    check_canonical_numerator=True,
    symmetric_polarizations=True,
)
assert not hasattr(fk, "GenerationOptions")

amplitude = fk.Process.amplitude(
    ["scalar_0"], ["scalar_0", "scalar_0"],
)
assert generator.generate(amplitude, max_vertices=3, **filter_arguments).report.completed
for mode in ("none", "zeroes", "identical", "up_to_sign", "up_to_scalar"):
    grouping = fk.NumeratorGrouping(mode, **grouping_arguments)
    assert generator.generate(
        amplitude, max_vertices=3, numerator_grouping=grouping,
    ).report.completed

for invalid_arguments in (
    dict(self_energy=fk.SelfEnergyFilterOptions(only_scaleless=True)),
    dict(tadpoles=fk.TadpoleFilterOptions(only_scaleless=True)),
    dict(zero_snails=fk.SnailFilterOptions(only_scaleless=True)),
    dict(fermion_loop_count_range=(2, 1)),
    dict(factorized_loop_topologies_count_range=(2, 1)),
    dict(blob_range=(0, 1)),
    dict(spectator_range=(0, 1)),
    dict(cut_amplitude_coupling_orders={"QCD": (0, None)}),
    dict(cut_amplitude_loop_count_range=(0, 1)),
):
    assert_feynkit_error(
        fk.GenerationError,
        lambda: generator.generate(amplitude, max_vertices=3, **invalid_arguments),
    )

for invalid_arguments in (
    dict(coupling_orders={"QCD": (2, 1)}),
    dict(loops=(2, 1)),
    dict(loops=(0, None)),
):
    try:
        model.generate_diagrams(["scalar_0"], ["scalar_0", "scalar_0"], **invalid_arguments)
    except ValueError:
        pass
    else:
        raise AssertionError("invalid generation range was accepted")

# Both entry points share configuration, including reusable exact/ranged orders.
kwargs = dict(max_vertices=3, vertex_allow=["V_3_SCALAR_000"], coupling_orders={"QCD": 1})
exact = generator.generate(amplitude, **kwargs)
ranged = model.generate_diagrams(
    ["scalar_0"], ["scalar_0", "scalar_0"],
    **dict(kwargs, coupling_orders={"QCD": (1, 1)}),
)
assert len(exact) > 0
assert [d.id for d in exact] == [d.id for d in ranged]
assert [d.id for d in generator.generate(amplitude, **kwargs)] == [d.id for d in exact]
assert kwargs["coupling_orders"] == {"QCD": 1}
assert len(generator.generate(amplitude, **dict(kwargs, coupling_orders={"QCD": 0}))) == 0
assert len(generator.generate(amplitude, **dict(kwargs, particle_veto=["scalar_0"]))) == 0
assert len(generator.generate(amplitude, **dict(kwargs, vertex_veto=["V_3_SCALAR_000"]))) == 0
assert len(generator.generate(amplitude, select_diagrams=[exact[0]], **kwargs)) == 1
assert len(generator.generate(amplitude, veto_diagrams=[d.id for d in exact], **kwargs)) == 0

for invalid_arguments in (dict(options=None), dict(coupling_order={"QCD": 0}), dict(coupling_orders={"QCD": True})):
    try:
        generator.generate(amplitude, **invalid_arguments)
    except TypeError:
        pass
    else:
        raise AssertionError("invalid generation keyword/type was accepted")

token = fk.CancellationToken()
token.cancel()
assert not generator.generate(amplitude, cancellation_token=token, **kwargs).report.completed

diagram = generator.generate(
    amplitude.with_loop_count(1, 1), max_vertices=3,
    vertex_allow=["V_3_SCALAR_000"],
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

    #[cfg(feature = "python_stubgen")]
    #[test]
    fn documentation_audit_accepts_documented_classes_methods_and_accessors() {
        validate_stub_documentation(
            r#"
class Diagram:
    """A diagram.

    Examples
    --------
    >>> diagram = model.generate_diagrams([], []).diagrams[0]
    """

    @property
    def loops(self) -> int:
        """Return the loop count."""
        ...

    def render(self, width: int) -> str:
        """Render the diagram.

        Examples
        --------
        >>> diagram.render(640)

        Parameters
        ----------
        width : int
            Target width.
        """
        ...
"#,
        )
        .unwrap();
    }

    #[cfg(feature = "python_stubgen")]
    #[test]
    fn documentation_audit_rejects_undocumented_accessors() {
        let error = validate_stub_documentation(
            r#"
class Diagram:
    """A diagram.

    Examples
    --------
    >>> diagram
    """

    @property
    def loops(self) -> int: ...
"#,
        )
        .unwrap_err();

        assert!(
            error
                .to_string()
                .contains("callable loops has no documentation")
        );
    }

    #[cfg(feature = "python_stubgen")]
    #[test]
    fn documentation_audit_rejects_missing_parameter_entries() {
        let error = validate_stub_documentation(
            r#"
class Generator:
    """A generator.

    Examples
    --------
    >>> generator.generate(diagram, options)
    """

    def generate(self, diagram: object, options: object) -> object:
        """Generate diagrams.

        Examples
        --------
        >>> generator.generate(diagram, options)

        Parameters
        ----------
        diagram : object
            Diagram to process.
        """
        ...
"#,
        )
        .unwrap_err();

        assert!(
            error
                .to_string()
                .contains("callable generate omits parameters: options")
        );
    }
}
