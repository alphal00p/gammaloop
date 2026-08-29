use std::ffi::CString;

use feynkit_graph::{DiagramEdge, DiagramVertex, FeynmanDiagram};
use feynkit_model::Model;
use feynkit_py::{FeynkitModule, PyFeynmanDiagram, PyModel};
use pyo3::{
    prelude::*,
    types::{PyCFunction, PyDict, PyList, PyModule},
};
use spynso3::SpensoModule;
use symbolica::{
    api::python::{SymbolicaCommunityModule, create_symbolica_module},
    atom::Atom,
    parser::ParseSettings,
};

const FEYNKIT_WRAPPER: &str = include_str!("../python/symbolica/community/feynkit/__init__.py");
const MODEL_JSON: &str = include_str!("fixtures/scalars_2p_3p.json");
const SPENSO_WRAPPER: &str = "from ..spenso_native import *\n\ninitialize_module()\n";

fn install_package<'py>(py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyModule>> {
    let package = PyModule::new(py, name)?;
    package.setattr("__package__", name)?;
    package.setattr("__path__", PyList::empty(py))?;
    py.import("sys")?
        .getattr("modules")?
        .set_item(name, &package)?;
    Ok(package)
}

fn register_native<'py, M>(core: &Bound<'py, PyModule>) -> PyResult<Bound<'py, PyModule>>
where
    M: SymbolicaCommunityModule + 'static,
{
    let py = core.py();
    let name = M::get_name();
    let native_name = format!("{name}_native");
    let native = PyModule::new(py, &native_name)?;
    let initialize_module =
        PyCFunction::new_closure(py, Some(c"initialize_module"), None, |args, _kwargs| {
            M::initialize(args.py())
        })?;
    native.add("initialize_module", initialize_module)?;
    M::register_module(&native)?;
    core.add_submodule(&native)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item(format!("symbolica.community.{native_name}"), &native)?;
    Ok(native)
}

fn import_wrapper<'py>(
    py: Python<'py>,
    community: &Bound<'py, PyModule>,
    name: &str,
    source: &str,
) -> PyResult<Bound<'py, PyModule>> {
    let full_name = format!("symbolica.community.{name}");
    let wrapper = PyModule::new(py, &full_name)?;
    wrapper.setattr("__package__", &full_name)?;
    wrapper.setattr("__path__", PyList::empty(py))?;
    py.import("sys")?
        .getattr("modules")?
        .set_item(&full_name, &wrapper)?;
    let source = CString::new(source).expect("community wrapper contains a null byte");
    py.run(&source, Some(&wrapper.dict()), Some(&wrapper.dict()))?;
    community.add(name, &wrapper)?;
    PyModule::import(py, &full_name)
}

fn remove_wrapper(py: Python<'_>, community: &Bound<'_, PyModule>, name: &str) -> PyResult<()> {
    py.import("sys")?
        .getattr("modules")?
        .del_item(format!("symbolica.community.{name}"))?;
    community.delattr(name)
}

fn tensor_diagram(
    model: &Model,
    name: &str,
    numerator: Atom,
    numerator_prefactor: Atom,
    projector: Atom,
) -> FeynmanDiagram {
    let rule = model.vertex_rule_id("V_3_SCALAR_000").unwrap();
    let particle = model.particle_id("scalar_0").unwrap();
    let mut left = DiagramVertex::interaction("left", rule);
    left.numerator = numerator.clone();
    let mut builder = FeynmanDiagram::builder(model.clone(), name)
        .numerator(numerator)
        .numerator_prefactor(numerator_prefactor)
        .projector(projector);
    let left = builder.add_vertex(left);
    let right = builder.add_vertex(DiagramVertex::interaction("right", rule));
    for _ in 0..3 {
        builder
            .add_edge(left, right, DiagramEdge::new(particle, false))
            .unwrap();
    }
    builder.build().unwrap()
}

#[test]
fn feynkit_and_spenso_share_one_symbolica_kernel_in_both_import_orders() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let symbolica = install_package(py, "symbolica")?;
        let core = PyModule::new(py, "symbolica.core")?;
        create_symbolica_module(&core)?;
        py.import("sys")?
            .getattr("modules")?
            .set_item("symbolica.core", &core)?;
        symbolica.add("core", &core)?;

        let community = install_package(py, "symbolica.community")?;
        symbolica.add("community", &community)?;
        register_native::<FeynkitModule>(&core)?;
        register_native::<SpensoModule>(&core)?;

        for order in [["feynkit", "spenso"], ["spenso", "feynkit"]] {
            for name in order {
                let source = match name {
                    "feynkit" => FEYNKIT_WRAPPER,
                    "spenso" => SPENSO_WRAPPER,
                    _ => unreachable!(),
                };
                import_wrapper(py, &community, name, source)?;
            }

            let feynkit = PyModule::import(py, "symbolica.community.feynkit")?;
            let spenso = PyModule::import(py, "symbolica.community.spenso")?;
            let model = Model::from_json(MODEL_JSON).expect("fixture is a valid model");
            let tensor_numerator = Atom::parse(
                "k(spenso::mink(D,mu))*k(spenso::mink(D,nu))",
                "feynkit_py_test",
                ParseSettings::default(),
            )
            .unwrap();
            let tensor_projector = Atom::parse(
                "p(spenso::mink(D,mu))*p(spenso::mink(D,nu))",
                "feynkit_py_test",
                ParseSettings::default(),
            )
            .unwrap();
            let free_tensor_numerator = Atom::parse(
                "k(spenso::mink(D,mu))*k(spenso::mink(D,nu))",
                "feynkit_py_test",
                ParseSettings::default(),
            )
            .unwrap();
            let free_tensor_diagram = tensor_diagram(
                &model,
                "free-tensor",
                free_tensor_numerator,
                Atom::one(),
                Atom::one(),
            );
            let diagram = tensor_diagram(
                &model,
                "shared-kernel",
                tensor_numerator,
                Atom::parse("y + 2", "feynkit_py_test", ParseSettings::default()).unwrap(),
                tensor_projector,
            )
            .with_overall_factor(
                Atom::parse("x + 1", "feynkit_py_test", ParseSettings::default()).unwrap(),
            );
            let diagram = Py::new(py, PyFeynmanDiagram::from(diagram))?;
            let free_tensor_diagram = Py::new(py, PyFeynmanDiagram::from(free_tensor_diagram))?;
            let model = Py::new(py, PyModel::from(model))?;
            let locals = PyDict::new(py);
            locals.set_item("core", &core)?;
            locals.set_item("diagram", diagram)?;
            locals.set_item("free_tensor_diagram", free_tensor_diagram)?;
            locals.set_item("model", model)?;
            locals.set_item("feynkit", feynkit)?;
            locals.set_item("spenso", spenso)?;
            let assertions = CString::new(
                r#"
feynkit_expression = diagram.overall_factor_expression()
spenso_expression = spenso.TensorName("T").to_expression()
named_coupling = core.S("UFO::SCALAR_COUPLING")
expanded_coupling = model.expand_couplings(named_coupling)

tensor_dimension = core.S("feynkit_py_test::D")
tensor_mu = core.S("feynkit_py_test::mu")
tensor_nu = core.S("feynkit_py_test::nu")
mink = core.S("spenso::mink")
dot = core.S("spenso::dot")
loop_momentum = core.S("feynkit_py_test::k")
projector_momentum = core.S("feynkit_py_test::p")
tensor_input = (
    loop_momentum(mink(tensor_dimension, tensor_mu))
    * loop_momentum(mink(tensor_dimension, tensor_nu))
    * projector_momentum(mink(tensor_dimension, tensor_mu))
    * projector_momentum(mink(tensor_dimension, tensor_nu))
)
tensor_reducer = feynkit.TensorReducer(tensor_dimension).with_integrated_head(
    "feynkit_py_test::k"
)
tensor_reduced = tensor_reducer.reduce(tensor_input)
diagram_tensor_reduced = diagram.reduce_tensor_numerator(tensor_reducer)
scalar_graphs = diagram.reduce_tensor_graphs(tensor_reducer)
free_tensor_reduced = free_tensor_diagram.reduce_tensor_numerator(tensor_reducer)
try:
    free_tensor_diagram.reduce_tensor_graphs(tensor_reducer)
except feynkit.TensorReductionError:
    pass
else:
    raise AssertionError("scalar graph splitting accepted residual free Lorentz indices")
tensor_expected = (
    dot(loop_momentum(mink(tensor_dimension)), loop_momentum(mink(tensor_dimension)))
    * dot(
        projector_momentum(mink(tensor_dimension)),
        projector_momentum(mink(tensor_dimension)),
    )
    / tensor_dimension
)
numerator_prefactor = core.S("feynkit_py_test::y") + 2

assert diagram.__class__ is feynkit.FeynmanDiagram
assert type(feynkit_expression) is core.Expression
assert type(spenso_expression) is core.Expression
assert type(feynkit_expression + spenso_expression) is core.Expression
assert type(expanded_coupling) is core.Expression
assert expanded_coupling == model.coupling("SCALAR_COUPLING").expression
assert named_coupling == core.S("UFO::SCALAR_COUPLING")
assert type(tensor_reduced) is core.Expression
assert tensor_reduced == tensor_expected
assert type(diagram_tensor_reduced) is core.Expression
assert diagram_tensor_reduced == tensor_expected
assert len(scalar_graphs) == 1
assert scalar_graphs[0].numerator_expression() == tensor_expected
assert scalar_graphs[0].projector_expression() == 1
assert scalar_graphs[0].numerator_prefactor_expression() == numerator_prefactor
assert scalar_graphs[0].overall_factor_expression() == feynkit_expression
scalar_graphs[0].validate()
assert type(free_tensor_reduced) is core.Expression
"#,
            )
            .unwrap();
            py.run(&assertions, Some(&locals), Some(&locals))?;

            remove_wrapper(py, &community, "feynkit")?;
            remove_wrapper(py, &community, "spenso")?;
        }
        Ok(())
    })
    .unwrap();
}
