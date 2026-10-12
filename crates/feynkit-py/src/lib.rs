//! Python bindings installed as `symbolica.community.hepkit`.

mod amplitude;
mod cff;
mod display;
mod error;
mod generation;
mod graph;
mod graph_interop;
mod integrals;
mod kinematics;
mod model;
mod progress;
mod render_settings;
mod sector_decomposition;
mod tensor;
#[cfg(feature = "ufo")]
mod ufo;
mod wavefunction;

use pyo3::{prelude::*, types::PyModule};
use symbolica::api::python::{Citation, SymbolicaCommunityModule};

pub use amplitude::{PyAmplitude, PyAmplitudeLeg, PySquaredAmplitude};
pub use cff::{
    PyCffGenerator, PyCffOrientation, PyCffReport, PyCffResult, PyCffSurface, PyCffSurfaceGroup,
    PyCutPropagator,
};
pub use generation::{
    PyAmplitudeOptions, PyBuiltDiagramSet, PyCancellationToken, PyCounts, PyCutAmplitude, PyCuts,
    PyDiagramGroup, PyDiagramSet, PyGenerationProgress, PyGenerationReport, PyGroupMember,
    PyGrouping, PyLoopCounts, PyLoops, PyNaming, PyNumeratorOptions, PyOrders, PyParticleSelector,
    PyPreparedDiagram, PyProcess, PyRange, PyRequest, PyRestrict, PyRun, PySelection, PySkeleton,
    PyStage, PySymmetrize, PyTopologyRules, PyVeto,
};
pub use graph::{
    PyDiagramCut, PyDiagramCutSide, PyDiagramEdge, PyDiagramThresholdCandidate, PyDiagramVertex,
    PyFeynmanDiagram, PyLoopMomentumBasis, PyMomentumSignature,
};
pub use integrals::{PyIntegralFamily, PyIntegralMapping, PyPropagatorMapping};
pub use kinematics::{
    PyAxis, PyBoost, PyClusteringResult, PyFourMomentum, PyHelicity, PyJet, PyJetAlgorithm,
    PyJetDefinition, PyKinematics, PyRotation, PyThreeMomentum,
};
pub use model::{
    PyCoupling, PyEvaluatedValues, PyEvaluationRequest, PyFormFactor, PyLorentzStructure, PyModel,
    PyModelExpression, PyModelFunction, PyParameter, PyParameterCard, PyParameterNature,
    PyParameterType, PyParticle, PyPropagator, PyVertexRule,
};
pub use progress::MarimoProgress;
pub use render_settings::PyDiagramStyle;
pub use tensor::PyTensorReducer;
#[cfg(feature = "ufo")]
pub use ufo::{PyLoadedModel, PyUfoLoadDiagnostics, PyUfoLoader};
pub use wavefunction::PyWavefunction;

mod symbols;

pub struct FeynkitModule;

impl SymbolicaCommunityModule for FeynkitModule {
    fn get_citations() -> Vec<Citation> {
        if !CITATIONS_USED.load(std::sync::atomic::Ordering::Relaxed) {
            return Vec::new();
        }
        // The product registry credits FeynKit contributors collectively.
        let mut citations = vec![
            Citation {
                id: "https://github.com/alphal00p/gammaloop#feynkit".into(),
                reference: "FeynKit contributors. FeynKit (2026).".into(),
                bibtex: r#"@software{feynkit,
  author = {{FeynKit contributors}},
  title = {FeynKit},
  year = {2026},
  url = {https://github.com/alphal00p/gammaloop}
}"#.into(),
                reasons: vec!["Provides the FeynKit functionality in this community module.".into()],
                description: "".into(),
                relevance: None,
            },
            Citation {
                id: "arXiv:0903.5143".into(),
                reference: "Benoît Collins and Sho Matsumoto. On some properties of orthogonal Weingarten functions (2009).".into(),
                bibtex: r#"@article{feynkit_0903.5143,
  author = {Benoît Collins and Sho Matsumoto},
  title = {{On some properties of orthogonal Weingarten functions}},
  year = {2009},
  eprint = {0903.5143},
  archivePrefix = {arXiv}
}"#.into(),
                reasons: vec!["Orthogonal Weingarten functions used by FeynKit tensor reduction.".into()],
                description: "Method and validation source listed in crates/feynkit-tensor/README.md.".into(),
                relevance: None,
            },
            Citation {
                id: "arXiv:1701.04493".into(),
                reference: "Benoît Collins and Sho Matsumoto. Weingarten calculus via orthogonality relations: new applications (2017).".into(),
                bibtex: r#"@article{feynkit_1701.04493,
  author = {Benoît Collins and Sho Matsumoto},
  title = {{Weingarten calculus via orthogonality relations: new applications}},
  year = {2017},
  eprint = {1701.04493},
  archivePrefix = {arXiv}
}"#.into(),
                reasons: vec!["Orthogonality relations underlying FeynKit tensor projectors.".into()],
                description: "Method and validation source listed in crates/feynkit-tensor/README.md.".into(),
                relevance: None,
            },
            Citation {
                id: "arXiv:1801.06084".into(),
                reference: "B. Ruijl and F. Herzog and T. Ueda and J. A. M. Vermaseren and A. Vogt. R*-operation and five-loop calculations (2018).".into(),
                bibtex: r#"@article{feynkit_1801.06084,
  author = {B. Ruijl and F. Herzog and T. Ueda and J. A. M. Vermaseren and A. Vogt},
  title = {{R*-operation and five-loop calculations}},
  year = {2018},
  eprint = {1801.06084},
  archivePrefix = {arXiv}
}"#.into(),
                reasons: vec!["Symmetry-reduced tensor projectors in the R* operation.".into()],
                description: "Method and validation source listed in crates/feynkit-tensor/README.md.".into(),
                relevance: None,
            },
            Citation {
                id: "arXiv:2408.05137".into(),
                reference: "Jae Goode and Franz Herzog and Anthony Kennedy and Sam Teale and Jos Vermaseren. Tensor Reduction for Feynman Integrals with Lorentz and Spinor Indices (2024).".into(),
                bibtex: r#"@article{feynkit_2408.05137,
  author = {Jae Goode and Franz Herzog and Anthony Kennedy and Sam Teale and Jos Vermaseren},
  title = {{Tensor Reduction for Feynman Integrals with Lorentz and Spinor Indices}},
  year = {2024},
  eprint = {2408.05137},
  archivePrefix = {arXiv}
}"#.into(),
                reasons: vec!["Orbit-partition formulation and high-rank tensor projectors.".into()],
                description: "Method and validation source listed in crates/feynkit-tensor/README.md.".into(),
                relevance: None,
            },
            Citation {
                id: "arXiv:2411.02233".into(),
                reference: "Jae Goode and Franz Herzog and Sam Teale. OPITeR: A program for tensor reduction of multi-loop Feynman Integrals (2024).".into(),
                bibtex: r#"@article{feynkit_2411.02233,
  author = {Jae Goode and Franz Herzog and Sam Teale},
  title = {{OPITeR: A program for tensor reduction of multi-loop Feynman Integrals}},
  year = {2024},
  eprint = {2411.02233},
  archivePrefix = {arXiv}
}"#.into(),
                reasons: vec!["Method reference for integrand-symmetry-aware tensor reduction.".into()],
                description: "Method and validation source listed in crates/feynkit-tensor/README.md.".into(),
                relevance: None,
            }
        ];
        if !TENSOR_CITATIONS_USED.load(std::sync::atomic::Ordering::Relaxed) {
            citations.truncate(1);
        }
        citations
    }

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
    module.add_class::<symbols::PySymbols>()?;

    render_settings::register(module)?;
    error::register(module)?;
    amplitude::register(module)?;
    model::register(module)?;
    graph::register(module)?;
    generation::register(module)?;
    kinematics::register(module)?;
    cff::register(module)?;
    tensor::register(module)?;
    integrals::register(module)?;
    #[cfg(feature = "ufo")]
    ufo::register(module)?;
    Ok(())
}

/// Gather the FeynKit stub inventory without declaring an independent wheel.
#[cfg(feature = "python_stubgen")]
pub fn stub_info() -> pyo3_stub_gen::Result<pyo3_stub_gen::StubInfo> {
    // Unqualified dependency classes belong to Symbolica; FeynKit declarations
    // explicitly name their community module.
    let mut info = pyo3_stub_gen::StubInfo::from_project_root(
        "symbolica".to_owned(),
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("python"),
    )?;
    let module = info
        .modules
        .get_mut("symbolica.community.hepkit")
        .ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "FeynKit did not contribute a symbolica.community.hepkit stub module",
            )
        })?;
    // Render empty tensor selector defaults as Python tuples.
    for class in module.class.values_mut() {
        for method in class.methods.values_mut().flatten() {
            if class.name == "TensorReducer" && method.name == "__new__" {
                for parameter in &mut method.parameters.keyword_only {
                    if matches!(parameter.name, "integrated" | "external") {
                        parameter.default =
                            pyo3_stub_gen::generate::ParameterDefault::Expr("()".to_owned());
                    }
                }
            }
        }
    }
    validate_stub_documentation(&module.to_string())?;
    Ok(info)
}

/// Render the same FeynKit surface for the package stub and API documentation.
#[cfg(feature = "python_stubgen")]
pub fn stub_source(module: &pyo3_stub_gen::generate::Module) -> String {
    normalize_stub_source(&module.to_string())
}

#[cfg(feature = "python_stubgen")]
fn normalize_stub_source(source: &str) -> String {
    let mut lines = source
        .lines()
        .map(|line| {
            let line = line.trim_end_matches([' ', '\t']);
            // These shared runtime classes are public reexports, whereas other
            // graph imports are only used in annotations. Same-name aliases
            // make that distinction explicit to Python type checkers.
            if let Some(imports) = line.strip_prefix("from symbolica.community.graph import ") {
                let imports = imports
                    .split(", ")
                    .map(|name| match name {
                        "Bridges" => "Bridges as Bridges",
                        "Classes" => "Classes as Classes",
                        _ => name,
                    })
                    .collect::<Vec<_>>()
                    .join(", ");
                format!("from symbolica.community.graph import {imports}")
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>();
    while lines.last().is_some_and(|line| line.is_empty()) {
        lines.pop();
    }
    lines.join("\n") + "\n"
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
    if not has_examples:
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

static CITATIONS_USED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

#[inline]
pub(crate) fn record_usage() {
    use std::sync::atomic::Ordering;
    if !CITATIONS_USED.load(Ordering::Relaxed) {
        CITATIONS_USED.store(true, Ordering::Relaxed);
    }
}

static TENSOR_CITATIONS_USED: std::sync::atomic::AtomicBool =
    std::sync::atomic::AtomicBool::new(false);

#[inline]
fn record_tensor_usage() {
    record_usage();
    spynso3::SpensoModule::record_usage();
    if !TENSOR_CITATIONS_USED.load(std::sync::atomic::Ordering::Relaxed) {
        TENSOR_CITATIONS_USED.store(true, std::sync::atomic::Ordering::Relaxed);
    }
}

#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use pyo3::types::PyDict;
    use symbolica::api::python::PythonExpression;

    use super::*;

    const NORMALIZED_SCALAR_MODEL: &str = include_str!("../tests/fixtures/scalars_2p_3p.json");

    fn registered_module<'py>(py: Python<'py>) -> Bound<'py, PyModule> {
        let module = PyModule::new(py, "symbolica.community.hepkit").unwrap();
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
        assert_send::<PyProcess>();
        assert_send::<PyFeynmanDiagram>();
        assert_send::<PyTensorReducer>();
    }

    #[cfg(feature = "python_stubgen")]
    #[test]
    fn normalizes_trailing_whitespace_and_end_of_file() {
        assert_eq!(
            normalize_stub_source("class Example:  \n    ...\t\n   \n\n"),
            "class Example:\n    ...\n"
        );
    }

    #[cfg(feature = "python_stubgen")]
    #[test]
    fn shared_graph_reexports_are_explicit_and_idempotent() {
        let source = concat!(
            "from symbolica.community.graph import Bridges, Classes, DiagramRender, RenderSettings  \n",
            "from other import Bridges, Classes\n",
            "class Example:\n",
            "    from symbolica.community.graph import Bridges, Classes\n",
        );
        let expected = source.replace(
            "import Bridges, Classes, DiagramRender, RenderSettings  ",
            "import Bridges as Bridges, Classes as Classes, DiagramRender, RenderSettings",
        );
        assert_eq!(normalize_stub_source(source), expected);
        assert_eq!(normalize_stub_source(&expected), expected);
    }

    #[cfg(feature = "python_stubgen")]
    #[test]
    fn generated_stub_covers_every_registered_native_class() {
        Python::initialize();
        Python::attach(|py| {
            let module = registered_module(py);
            let info = stub_info().unwrap();
            let source = stub_source(info.modules.get("symbolica.community.hepkit").unwrap());
            let graph = PyModule::new(py, "symbolica.community.graph").unwrap();
            linnet_py::register(&graph).unwrap();
            let locals = PyDict::new(py);
            locals.set_item("fk", module).unwrap();
            locals.set_item("graph", graph).unwrap();
            locals.set_item("source", source).unwrap();
            let code = CString::new(
                r#"
import ast

tree = ast.parse(source)
stub_classes = {
    node.name for node in ast.walk(tree)
    if isinstance(node, ast.ClassDef)
}
runtime_classes = {
    name for name, value in vars(fk).items()
    if isinstance(value, type) and value.__module__ == fk.__name__
}
missing = sorted(runtime_classes - stub_classes)
assert not missing, f"native classes missing from the generated stub: {missing}"

shared_imports = {
    alias.name: alias.asname
    for node in tree.body
    if isinstance(node, ast.ImportFrom) and node.module == graph.__name__
    for alias in node.names
    if alias.name in {"Bridges", "Classes"}
}
runtime_shared = {
    name: value for name, value in vars(fk).items()
    if name in {"Bridges", "Classes"} and isinstance(value, type)
}
assert "Bridges" in runtime_shared
assert shared_imports == {name: name for name in runtime_shared}
for name, value in runtime_shared.items():
    assert value is getattr(graph, name), f"{name} is not the shared graph class"
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
            let graph_module = PyModule::new(py, "symbolica.community.graph").unwrap();
            linnet_py::register(&graph_module).unwrap();
            locals.set_item("graph", &graph_module).unwrap();
            locals
                .set_item("MODEL_JSON", NORMALIZED_SCALAR_MODEL)
                .unwrap();
            let code = CString::new(
                r#"
import sys

assert "_gammaloop" not in sys.modules
assert fk.__name__ == "symbolica.community.hepkit"

model = fk.Model.from_json(MODEL_JSON)
assert model.name == "scalars"
scalar = model.particle("scalar_0")
assert scalar.pdg_code == 1000
assert scalar.antiparticle.name == "scalar_0"
model = fk.Model.from_json(model.to_json(pretty=False))

scalar = model.particle("scalar_0")
process = model.process([scalar], [1000, scalar.antiparticle])
assert isinstance(process, fk.Process)
assert not hasattr(process, "loop_count")
restriction = fk.Restrict().only_vertices(["V_3_SCALAR_000"]).exclude_particles([model.particle("scalar_1"),1002])
topology = fk.TopologyRules(max_vertices=3, bridges=graph.Bridges.one_pi_in_loops(), self_loops=True,
    self_energies=fk.Veto.all(), tadpoles=fk.Veto.all(), snails=fk.Veto.massless(),
    fermion_loops=fk.Range.exact(0), factorized_loops=fk.Range.between(0,1))
request = process.amplitude(loops=fk.Loops.total((0,1)), restrict=restriction, orders=fk.Orders().min("QCD",0), topology=topology)
generated = request.generate().build(fk.NumeratorOptions(naming=fk.Naming.named("D")))
assert len(generated) > 0
assert generated[0].name == next(iter(generated)).name
assert generated.report.completed
generation_html = generated._repr_html_()
assert "Built diagrams" in generation_html
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
diagram_svg = loop_diagram.render().to_svg()
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
            let graph_module = PyModule::new(py, "symbolica.community.graph").unwrap();
            linnet_py::register(&graph_module).unwrap();
            locals.set_item("graph", &graph_module).unwrap();
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
        assert type(error).__module__ == "symbolica.community.hepkit"
    else:
        raise AssertionError(f"{error_type.__name__} was not raised")

assert_feynkit_error(fk.ModelError, lambda: fk.Model.from_json("{}"))
assert_feynkit_error(
    fk.GenerationError,
    lambda: model.process(["missing_particle"], []),
)
assert_feynkit_error(
    fk.DiagramError,
    lambda: fk.FeynmanDiagram.from_json(model, "{}"),
)
assert_feynkit_error(
    fk.KinematicsError,
    lambda: fk.Boost(fk.ThreeMomentum(1.0, 0.0, 0.0)),
)


topology = fk.TopologyRules(max_vertices=3, bridges=graph.Bridges.one_pi_in_loops(),
    self_energies=fk.Veto.massless(), tadpoles=fk.Veto.massless(), snails=fk.Veto.massive(),
    fermion_loops=fk.Range.between(0,1), factorized_loops=fk.Range.between(0,1))
grouping_arguments = dict(
    numerical_sample_seed=7, number_of_numerical_samples=11,
    differentiate_particle_masses_only=False, fully_numerical_substitution=True,
    check_canonical_numerator=True, symmetric_polarizations=True,
)
assert not hasattr(fk, "GenerationOptions")
assert not hasattr(fk, "GenerationSettings")
amplitude = model.process(["scalar_0"],["scalar_0","scalar_0"])
restriction = fk.Restrict().only_vertices(["V_3_SCALAR_000"])
request = amplitude.amplitude(loops=fk.Loops.total(0), restrict=restriction, topology=topology, orders=fk.Orders().min("QCD",0))
assert request.generate().report.completed
for mode in ("none","identical","up_to_sign","up_to_scalar"):
    assert request.generate().build().group(fk.Grouping(mode,**grouping_arguments)).report.completed
# Contradictory constraints remain explicit empty requests, rather than clearing filters.
for empty_topology in (
    fk.TopologyRules(max_vertices=3, fermion_loops=fk.Range.between(2,1)),
    fk.TopologyRules(max_vertices=3, factorized_loops=fk.Range.between(2,1)),
):
    assert len(amplitude.amplitude(loops=fk.Loops.total(0),restrict=restriction,topology=empty_topology).generate()) == 0
for cuts in (
    fk.Cuts(blobs=fk.Range.between(0,1)), fk.Cuts(spectators=fk.Range.between(0,1)),
    fk.Cuts(amplitude=fk.CutAmplitude(orders=fk.Orders().min("QCD",0))),
    fk.Cuts(amplitude=fk.CutAmplitude(loops=fk.Range.between(0,1))),
):
    # Amplitude Process has no cuts argument; the typed request prevents this misuse.
    try:
        amplitude.amplitude(loops=fk.Loops.total(0),cuts=cuts)
    except TypeError:
        pass
    else:
        raise AssertionError("amplitude accepted cross-section cut options")
for loop_range in ((2,1),(0,None)):
    try:
        fk.Loops.total(loop_range)
    except ValueError:
        pass
    else:
        raise AssertionError("invalid loop range was accepted")
exact_orders = fk.Orders().exact("QCD",1)
exact = amplitude.amplitude(loops=fk.Loops.total(0),restrict=restriction,orders=exact_orders).generate().build()
ranged = amplitude.amplitude(loops=fk.Loops.total(0),restrict=restriction,orders=fk.Orders(plain={"QCD":fk.Range.between(1,1)})).generate().build()
assert len(exact) > 0
assert [d.id for d in exact] == [d.id for d in ranged]
assert [d.id for d in amplitude.amplitude(loops=fk.Loops.total(0),restrict=restriction,orders=exact_orders).generate().build()] == [d.id for d in exact]
assert len(amplitude.amplitude(loops=fk.Loops.total(0),restrict=restriction,orders=fk.Orders().exact("QCD",0)).generate()) == 0
assert len(amplitude.amplitude(loops=fk.Loops.total(0),restrict=restriction.exclude_particles(["scalar_0"])).generate()) == 0
assert len(amplitude.amplitude(loops=fk.Loops.total(0),restrict=restriction.exclude_vertices(["V_3_SCALAR_000"])).generate()) == 0
assert len(exact.select(fk.Selection(select=[exact[0]]))) == 1
assert len(exact.select(fk.Selection(veto=[d.id for d in exact]))) == 0
for make in (
    lambda: amplitude.amplitude(loops=fk.Loops.total(0),options=None),
    lambda: amplitude.amplitude(loops=fk.Loops.total(0),coupling_order={"QCD":0}),
    lambda: fk.Orders().exact("QCD",True),
    lambda: fk.Range.exact(True),
):
    try:
        make()
    except TypeError:
        pass
    else:
        raise AssertionError("invalid generation keyword/type was accepted")
assert_feynkit_error(fk.GenerationError,lambda: amplitude.amplitude(loops=fk.Loops.total(0),orders=fk.Orders().perturbative("QCD",1)))
token = fk.CancellationToken()
token.cancel()
assert not request.generate(run=fk.Run(cancellation=token)).report.completed
loop_topology=fk.TopologyRules(max_vertices=3,bridges=graph.Bridges.one_pi_in_loops(),self_energies=fk.Veto.all(),tadpoles=fk.Veto.all(),snails=fk.Veto.massless())
diagram = amplitude.amplitude(loops=fk.Loops.total(1),restrict=restriction,topology=loop_topology).generate().build()[0]
assert_feynkit_error(fk.CffError,lambda:fk.CffGenerator(max_orientations=0).generate(diagram))

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
    >>> diagram = model.process([], []).amplitude(loops=hep.Loops.total(2)).generate().build()[0]
    """

    @property
    def loops(self) -> int:
        """Return the loop count.

        Examples
        --------
        >>> diagram.loops
        1
        """
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
