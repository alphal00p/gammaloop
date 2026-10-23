use std::ffi::CString;

use feynkit_graph::{DiagramEdge, DiagramVertex, FeynmanDiagram, NumeratorInterface};
use feynkit_model::Model;
use feynkit_py::{FeynkitModule, PyFeynmanDiagram, PyModel};
#[cfg(feature = "linnet-python")]
use linnet_py::GraphModule;
#[cfg(not(feature = "linnet-python"))]
use linnet_render_py::GraphModule;
use pyo3::prelude::{Bound, Py, PyAnyMethods, PyDictMethods, PyModuleMethods, PyResult, Python};
use pyo3::types::{PyCFunction, PyDict, PyList, PyModule};
use spynso3::SpensoModule;
use symbolica::api::python;
use symbolica::api::python::SymbolicaCommunityModule;
use symbolica::atom::Atom;
use symbolica::parser::ParseSettings;

const FEYNKIT_WRAPPER: &str = include_str!("../python/symbolica/community/hepkit/__init__.py");
#[cfg(feature = "linnet-python")]
const GRAPH_WRAPPER: &str =
    include_str!("../../linnet-py/python/symbolica/community/graph/__init__.py");
#[cfg(not(feature = "linnet-python"))]
const GRAPH_WRAPPER: &str =
    include_str!("../../linnet-render-py/python/symbolica/community/graph/__init__.py");
const MODEL_JSON: &str = include_str!("fixtures/scalars_2p_3p.json");
const SPENSO_WRAPPER: &str = "from ..tensor_native import *\n\ninitialize_module()\n";

#[test]
fn shared_numeric_external_states_in_installed_host() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let symbolica = install_package(py, "symbolica")?;
        let core = PyModule::new(py, "symbolica.core")?;
        python::create_symbolica_module(&core)?;
        py.import("sys")?
            .getattr("modules")?
            .set_item("symbolica.core", &core)?;
        symbolica.add("core", &core)?;
        let community = install_package(py, "symbolica.community")?;
        symbolica.add("community", &community)?;
        register_native::<GraphModule>(&core)?;
        register_native::<FeynkitModule>(&core)?;
        import_wrapper(
            py,
            &community,
            "symbolica.community.hepkit",
            FEYNKIT_WRAPPER,
        )?;
        let code = CString::new(include_str!("installed_wavefunctions.py")).unwrap();
        let locals = PyDict::new(py);
        py.run(&code, Some(&locals), Some(&locals))
    })
    .unwrap();
}

#[test]
fn qcd_rich_display_requires_no_python_rendering_packages() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let symbolica = install_package(py, "symbolica")?;
        let core = PyModule::new(py, "symbolica.core")?;
        python::create_symbolica_module(&core)?;
        py.import("sys")?
            .getattr("modules")?
            .set_item("symbolica.core", &core)?;
        symbolica.add("core", &core)?;
        let community = install_package(py, "symbolica.community")?;
        symbolica.add("community", &community)?;
        register_native::<GraphModule>(&core)?;
        register_native::<FeynkitModule>(&core)?;
        register_native::<SpensoModule>(&core)?;
        import_wrapper(py, &community, "symbolica.community.graph", GRAPH_WRAPPER)?;
        #[cfg(not(feature = "linnet-python"))]
        {
            let graph = py.import("symbolica.community.graph")?;
            for name in [
                "DiagramRender", "DiagramView", "RenderSettings", "LayoutSettings", "Stroke",
            ] {
                assert!(graph.hasattr(name)?, "missing shared display class {name}");
            }
            for name in ["Graph", "Subgraph", "node", "edge", "build"] {
                assert!(!graph.hasattr(name)?, "unexpected generic graph API {name}");
            }
        }
        import_wrapper(
            py,
            &community,
            "symbolica.community.hepkit",
            FEYNKIT_WRAPPER,
        )?;
        import_wrapper(py, &community, "symbolica.community.tensor", SPENSO_WRAPPER)?;
        let code = CString::new(include_str!("installed_offline_rendering.py")).unwrap();
        let locals = PyDict::new(py);
        py.run(&code, Some(&locals), Some(&locals))?;
        let network_tests = CString::new(include_str!(
            "../../spynso3/tests/installed_tensor_network_display.py"
        )).unwrap();
        locals.set_item("__name__", "network_display_tests")?;
        py.run(&network_tests, Some(&locals), Some(&locals))?;
        py.run(c"assert unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(NetworkDisplayTests)).wasSuccessful()", Some(&locals), Some(&locals))?;
        symbolica.add("E", core.getattr("E")?)?;
        #[cfg(feature = "linnet-python")]
        {
            let metrics_source = CString::new(include_str!("../../linnet-py/tests/rendering_metrics.py")).unwrap();
            let metrics = PyModule::from_code(py, &metrics_source, c"rendering_metrics.py", c"rendering_metrics")?;
            py.import("sys")?.getattr("modules")?.set_item("rendering_metrics", metrics)?;
            let rendering_tests = CString::new(include_str!(
                "../../linnet-py/tests/test_community_rendering.py"
            )).unwrap();
            locals.set_item("__name__", "community_rendering_tests")?;
            py.run(&rendering_tests, Some(&locals), Some(&locals))?;
            py.run(c"assert unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(CommunityRenderingTests)).wasSuccessful()", Some(&locals), Some(&locals))?;
        }
        // The default host validates every label source without a Rust compiler.
        // The optional typst-py harness owns actual vector-export compilation.
        py.run(
            cr#"
for owner in (network, diagram):
    view = owner.render()
    request = view.measurement_request()
    assert request["labels"]
    snapshot = external_metrics_snapshot(view)
    assert "<math" in snapshot.to_svg()
    assert view.measured().to_svg() == snapshot.to_svg()
assert "typst" not in sys.modules
"#,
            Some(&locals),
            Some(&locals),
        )?;
        Ok(())
    })
    .unwrap();
}

#[test]
fn community_native_forwarding_uses_registered_owners() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let symbolica = install_package(py, "symbolica")?;
        let core = PyModule::new(py, "symbolica.core")?;
        python::create_symbolica_module(&core)?;
        py.import("sys")?
            .getattr("modules")?
            .set_item("symbolica.core", &core)?;
        symbolica.add("core", &core)?;
        let community = install_package(py, "symbolica.community")?;
        symbolica.add("community", &community)?;
        register_native::<GraphModule>(&core)?;
        register_native::<FeynkitModule>(&core)?;
        register_native::<SpensoModule>(&core)?;
        let graph = import_wrapper(py, &community, "symbolica.community.graph", GRAPH_WRAPPER)?;
        let fk = import_wrapper(
            py,
            &community,
            "symbolica.community.hepkit",
            FEYNKIT_WRAPPER,
        )?;
        let tensor = import_wrapper(py, &community, "symbolica.community.tensor", SPENSO_WRAPPER)?;
        let locals = PyDict::new(py);
        locals.set_item("fk", fk)?;
        locals.set_item("tensor", tensor)?;
        locals.set_item("graph", graph)?;
        locals.set_item("core", core)?;
        locals.set_item("MODEL_JSON", MODEL_JSON)?;
        locals.set_item("standalone_host", false)?;
        let code = CString::new(include_str!("native_forwarding.py")).unwrap();
        py.run(&code, Some(&locals), Some(&locals))
    })
    .unwrap();
}

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
    parent: &Bound<'py, PyModule>,
    full_name: &str,
    source: &str,
) -> PyResult<Bound<'py, PyModule>> {
    let wrapper = PyModule::new(py, full_name)?;
    wrapper.setattr("__package__", full_name)?;
    wrapper.setattr("__path__", PyList::empty(py))?;
    py.import("sys")?
        .getattr("modules")?
        .set_item(full_name, &wrapper)?;
    let source = CString::new(source).expect("community wrapper contains a null byte");
    py.run(&source, Some(&wrapper.dict()), Some(&wrapper.dict()))?;
    parent.add(full_name.rsplit('.').next().unwrap(), &wrapper)?;
    PyModule::import(py, full_name)
}

fn remove_wrapper(py: Python<'_>, parent: &Bound<'_, PyModule>, full_name: &str) -> PyResult<()> {
    py.import("sys")?.getattr("modules")?.del_item(full_name)?;
    parent.delattr(full_name.rsplit('.').next().unwrap())
}

fn tensor_diagram(
    model: &Model,
    name: &str,
    numerator: Atom,
    numerator_prefactor: Atom,
    projector: Atom,
) -> FeynmanDiagram {
    let numerator = NumeratorInterface::admit(numerator, None).unwrap();
    let numerator_prefactor = NumeratorInterface::admit(numerator_prefactor, None).unwrap();
    let projector = NumeratorInterface::admit(projector, None).unwrap();
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
        python::create_symbolica_module(&core)?;
        py.import("sys")?
            .getattr("modules")?
            .set_item("symbolica.core", &core)?;
        symbolica.add("core", &core)?;

        let community = install_package(py, "symbolica.community")?;
        symbolica.add("community", &community)?;
        register_native::<GraphModule>(&core)?;
        register_native::<FeynkitModule>(&core)?;
        register_native::<SpensoModule>(&core)?;
        import_wrapper(py, &community, "symbolica.community.graph", GRAPH_WRAPPER)?;

        for order in [["hepkit", "tensor"], ["tensor", "hepkit"]] {
            for name in order {
                let (parent, full_name, source) = match name {
                    "hepkit" => (&community, "symbolica.community.hepkit", FEYNKIT_WRAPPER),
                    "tensor" => (&community, "symbolica.community.tensor", SPENSO_WRAPPER),
                    _ => unreachable!(),
                };
                import_wrapper(py, parent, full_name, source)?;
            }

            let feynkit = PyModule::import(py, "symbolica.community.hepkit")?;
            let spenso = PyModule::import(py, "symbolica.community.tensor")?;
            let model = Model::from_json(MODEL_JSON).expect("fixture is a valid model");
            // Structural ports require registered tensor heads before parsing.
            let registration = PyDict::new(py);
            registration.set_item("spenso", &spenso)?;
            py.run(
                c"spenso.TensorName.vector('feynkit_py_test::k')\nspenso.TensorName.vector('feynkit_py_test::p')",
                Some(&registration),
                Some(&registration),
            )?;
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
            let indexed = Atom::parse(
                "gammalooprs::Q(1,spenso::mink(4,gammalooprs::hedge(7,0)))*gammalooprs::Q(2,spenso::mink(4,gammalooprs::hedge(8,0)))*gammalooprs::Q(3,spenso::mink(4,gammalooprs::edge(7,0)))*gammalooprs::Q(4,spenso::mink(4,gammalooprs::vertex(7,0)))",
                "feynkit_py_test",
                ParseSettings::default(),
            ).unwrap();
            let indexed_diagram = tensor_diagram(&model, "indexed", indexed, Atom::one(), Atom::one());
            let indexed_diagram = Py::new(py, PyFeynmanDiagram::from(indexed_diagram))?;
            let diagram = Py::new(py, PyFeynmanDiagram::from(diagram))?;
            let free_tensor_diagram = Py::new(py, PyFeynmanDiagram::from(free_tensor_diagram))?;
            let model = Py::new(py, PyModel::from(model))?;
            let locals = PyDict::new(py);
            locals.set_item("core", &core)?;
            locals.set_item("diagram", diagram)?;
            locals.set_item("indexed_diagram", indexed_diagram)?;
            locals.set_item("free_tensor_diagram", free_tensor_diagram)?;
            locals.set_item("model", model)?;
            locals.set_item("feynkit", feynkit)?;
            locals.set_item("spenso", spenso)?;
            let assertions = CString::new(
                r#"
def external_metrics_snapshot(view):
    # EXTERNAL METRICS FIXTURE in points, not browser measurements.
    request = view.measurement_request()
    return view.measured(measurements={
        "id": request["id"],
        "metrics": [
            {"id": label["id"], "width": 24.0, "height": 12.0, "baseline": 9.0}
            for label in request["labels"]
        ],
    })


core.set_namespace("community_namespace_test")
assert spenso.TensorName.vector("namespace_vector").to_expression() == core.S(
    "community_namespace_test::namespace_vector"
)
assert spenso.TensorName("explicit_namespace::tensor").to_expression() == core.S(
    "explicit_namespace::tensor"
)
core.set_namespace("python")
feynkit_expression = diagram.overall_factor_expression()
sunrise = feynkit.FeynmanDiagram.from_dot(feynkit.Model.phi4(), '''digraph sunrise {
    ext [style=invis];
    ext -> a [particle="phi"];
    a -> b [particle="phi", lmb_id=0];
    a -> b [particle="phi", lmb_id=1];
    a -> b [particle="phi"];
    b -> ext [particle="phi"];
}''')
family = sunrise.integral_family()
assert family.is_complete and family.is_independent
assert len(family.denominators) == 5
products = [family.kinematics.scalar_product(k, family.external_momenta[0])
            for k in family.loop_momenta]
preferred = sunrise.integral_family(products)
assert preferred.denominators == feynkit.IntegralFamily.from_diagram(sunrise, products).denominators
assert preferred.denominators[:3] == sunrise.propagator_family().denominators
assert preferred.denominators[3:] == products
assert sunrise.integral_family(
    independent_dot_products=products[:1], kinematics=family.kinematics
).is_complete
try:
    feynkit.IntegralFamily.from_diagram(sunrise, products + [products[0]**2])
except feynkit.DiagramError:
    pass
else:
    raise AssertionError("nonlinear auxiliary product was accepted")
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
tensor_reducer = feynkit.TensorReducer(tensor_dimension, integrated=[loop_momentum])
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
assert type(tensor_reduced) is spenso.TensorExpression
assert tensor_reduced.rank == 0
assert tensor_reduced.axes == ()
assert tensor_reduced.to_expression() == tensor_expected
assert type(diagram_tensor_reduced) is spenso.TensorExpression
assert diagram_tensor_reduced.rank == 0
assert diagram_tensor_reduced.axes == ()
assert diagram_tensor_reduced.to_expression() == tensor_expected
assert len(scalar_graphs) == 1
assert scalar_graphs[0].numerator_expression().to_expression() == tensor_expected
assert isinstance(scalar_graphs[0].numerator_expression(), spenso.TensorExpression)
auto_reduced = scalar_graphs[0].tensor_reduce(tensor_dimension)
assert isinstance(auto_reduced, spenso.TensorExpression)
assert auto_reduced == scalar_graphs[0].numerator_expression()
graph_momentum = core.S("gammalooprs::Q")
graph_loop_edge = diagram.loop_momentum_basis.loop_edges[0]
expanded_input = spenso.as_tensor(
    3 + 2
    * graph_momentum(graph_loop_edge, mink(4, tensor_mu))
    * graph_momentum(graph_loop_edge, mink(4, tensor_nu))
    * projector_momentum(mink(4, tensor_mu))
    * projector_momentum(mink(4, tensor_nu))
)
expanded_reduced = diagram.tensor_reduce(tensor_dimension, expression=expanded_input)
graph_loop_vector = graph_momentum(graph_loop_edge, mink(tensor_dimension))
expanded_expected = 3 + 2 * dot(graph_loop_vector, graph_loop_vector) * dot(
    projector_momentum(mink(tensor_dimension)),
    projector_momentum(mink(tensor_dimension)),
) / tensor_dimension
assert isinstance(expanded_reduced, spenso.TensorExpression)
assert (expanded_reduced - expanded_expected).to_expression().expand() == 0
try:
    diagram.tensor_reduce(tensor_dimension, expression=expanded_input, projector=1)
except feynkit.DiagramError:
    pass
else:
    raise AssertionError("tensor reduction accepted both an expression and a projector")
# The parent's native MathML callback renders standalone spaces as family/dimension
# subscripts; ordinary function calls retain parentheses. Both preserve namespace
# options, and neither display changes the stored Symbolica expression.
import re
import xml.etree.ElementTree as ET
ordinary_numerator = scalar_graphs[0].numerator_expression().to_expression()
ordinary_before = bytes(ordinary_numerator.__getstate__())
representation = mink(tensor_dimension)
control_representation = core.S("printer_control::mink")(tensor_dimension)
representation_html = representation._repr_html_()
control_representation_html = control_representation._repr_html_()
representation_math, control_math = [
    ET.fromstring(re.search(r"<math\b[^>]*>.*?</math>", html, re.DOTALL).group())
    for html in (representation_html, control_representation_html)
]
family_dimension = representation_math.find(".//{*}msub")
assert family_dimension is not None and len(family_dimension) == 2, representation_html
assert [node.text for node in family_dimension[0].iter() if node.tag.rsplit("}", 1)[-1] == "mi"] == ["mink"], representation_html
assert [node.text for node in family_dimension[1].iter() if node.tag.rsplit("}", 1)[-1] == "mi"] == ["D"], representation_html
assert representation_math.find(".//{*}mmultiscripts") is None
assert all(node.get("title") is not None for node in control_math.findall(".//{*}msub")), control_representation_html
assert [node.text for node in control_math.findall(".//{*}mi")] == ["mink", "D"], control_representation_html
control_operators = [node.text for node in control_math.findall(".//{*}mo")]
assert "(" in control_operators and ")" in control_operators, control_representation_html
for value, namespace in ((representation, "spenso"), (control_representation, "printer_control")):
    before_display = bytes(value.__getstate__())
    for show_namespaces in (False, True):
        plain = value.format(show_namespaces=show_namespaces, color_builtin_symbols=False, bracket_level_colors=None)
        mathml = value.format(mode=core.PrintMode.MathML, show_namespaces=show_namespaces,
                              color_builtin_symbols=False, bracket_level_colors=None)
        math_root = ET.fromstring(f"<math>{mathml}</math>")
        assert [node.text for node in math_root.findall(".//{*}mi")] == ["mink", "D"], mathml
        head = f"{namespace}::mink" if show_namespaces else "mink"
        dimension = "feynkit_py_test::D" if show_namespaces else "D"
        assert plain == f"{head}({dimension})", plain
        # Native MathML displays namespaces as subordinate badges on identifiers.
        qualified = [node.get("title") for node in math_root.findall(".//{*}msub[@title]")]
        badges = [node.text for node in math_root.findall(".//{*}mtext")]
        if show_namespaces:
            assert qualified == [f"{namespace}::mink", "feynkit_py_test::D"], mathml
            assert badges == [namespace, "feynkit_py_test"], mathml
        else:
            assert "::" not in plain and "::" not in "".join(math_root.itertext())
            assert qualified == badges == [], mathml
    assert bytes(value.__getstate__()) == before_display
assert "spenso::" not in ordinary_numerator._repr_html_()
assert "spenso::mink" in ordinary_numerator.format(show_namespaces=True, color_builtin_symbols=False, bracket_level_colors=None)
assert isinstance(ordinary_numerator, core.Expression)
assert bytes(ordinary_numerator.__getstate__()) == ordinary_before

assert scalar_graphs[0].numerator_expression().rank == 0
denominator = indexed_diagram.denominator_expression()
assert isinstance(denominator, spenso.TensorExpression)
assert denominator.rank == 0
assert "ZERO" not in str(denominator)
assert "Momentum" not in denominator._repr_latex_()
indexed = indexed_diagram.numerator_expression()
inverse_denominator = denominator ** -1
assert indexed / denominator == indexed * inverse_denominator
assert isinstance(indexed, spenso.TensorExpression)
assert isinstance(indexed, core.Expression)
assert indexed.rank == 4
indexed_source = spenso.TensorExpression(
    core.Expression.parse(indexed_diagram.numerator), intern="indices"
)
assert indexed.to_expression() == indexed_source.to_expression()
assert indexed.axes == indexed_source.axes
assert "gammalooprs::Q(" in indexed_diagram.numerator
assert "FeynKit::" not in indexed_diagram.numerator
assert "Momentum" not in str(indexed)
assert "q₁" in str(indexed)
assert "q_{1}" in indexed._repr_latex_()
assert "Momentum" not in indexed.to_typst()
# Raw index notation exposes the graph provenance under the alphabetic display default.
raw_indices = spenso.DisplaySettings(index_style="raw")
for label in ("e₇.₀", "v₇.₀"):
    assert label in indexed.format_tensor(settings=raw_indices)
for label in ("e", "v"):
    assert f"{label}_{{7.0}}" in indexed.to_latex(settings=raw_indices)
    assert f"attach({label},b:(7.0))" in indexed.to_typst(settings=raw_indices)
for vertex in indexed_diagram.vertices:
    assert isinstance(vertex.numerator_expression(), spenso.TensorExpression)
for edge in indexed_diagram.edges:
    assert isinstance(edge.numerator_expression(), spenso.TensorExpression)
# Scalar graphs retain dashed-particle styling in captured native snapshots.
view = diagram.render()
snapshot = external_metrics_snapshot(view)
assert "stroke-dasharray" in snapshot.to_svg()
assert "stroke-dasharray" in snapshot.to_html()
assert view.measured().to_svg() == snapshot.to_svg()
assert not hasattr(diagram, "_repr_html_")
assert "DiagramView" in diagram._repr_mimebundle_(include=["text/plain"])["text/plain"]
# Complete direct MathML keeps normalized inverse factors; paging remains upstream.
rational = indexed / core.Expression.parse("feynkit_py_test::x*feynkit_py_test::y")
before_render = rational.to_expression()
fraction_html = rational.to_html()
fraction_math = ET.fromstring(re.search(r"<math\b[^>]*>.*?</math>", fraction_html, re.DOTALL).group())
assert before_render.get_type() == core.AtomType.Mul
source_factors = list(before_render)
source_inverses = [tuple(factor) for factor in source_factors if factor.get_type() == core.AtomType.Pow]
assert len(source_inverses) == 2
assert {base for base, exponent in source_inverses} == {
    core.S("feynkit_py_test::x"), core.S("feynkit_py_test::y"),
}
assert all(exponent == -1 for base, exponent in source_inverses)
powers = fraction_math.findall(".//{*}msup")
assert len(powers) == 2, fraction_html
rendered_inverses = []
for power in powers:
    assert len(power) == 2, fraction_html
    assert [node.tag.rsplit("}", 1)[-1] for node in power] == ["mrow", "mrow"], fraction_html
    base, exponent = power
    assert len(base) == 1 and base[0].tag.rsplit("}", 1)[-1] == "mi", fraction_html
    assert base[0].text in {"x", "y"}, fraction_html
    assert len(exponent) == 1 and exponent[0].tag.rsplit("}", 1)[-1] == "mrow", fraction_html
    assert [(node.tag.rsplit("}", 1)[-1], node.text) for node in exponent[0]] == [
        ("mo", "−"), ("mn", "1"),
    ], fraction_html
    rendered_inverses.append((core.S(f"feynkit_py_test::{base[0].text}"), -1))
assert rendered_inverses == source_inverses, fraction_html
# The product's immediate children retain the actual normalized source order.
expected_order = [
    ("msup" if factor.get_type() == core.AtomType.Pow else "mmultiscripts", factor[0])
    for factor in source_factors
]
product_rows = [
    node for node in fraction_math.findall(".//{*}mrow")
    if [child.tag.rsplit("}", 1)[-1] for child in node
        if child.tag.rsplit("}", 1)[-1] != "mspace"] == [kind for kind, value in expected_order]
]
assert len(product_rows) == 1, fraction_html
rendered_order = [
    (child.tag.rsplit("}", 1)[-1],
     core.S(f"feynkit_py_test::{child[0][0].text}") if child.tag.rsplit("}", 1)[-1] == "msup"
     else core.Expression.parse("".join(child[1].itertext())))
    for child in product_rows[0] if child.tag.rsplit("}", 1)[-1] != "mspace"
]
assert rendered_order == expected_order, fraction_html
assert rational.to_expression() == before_render
html = indexed.to_html()
# Standalone notebook fragments share a pinned font URL without embedding font files.
assert "data:font" not in html
assert "https://cdn.jsdelivr.net/gh/stipub/stixfonts@2.13b171/" in html
assert "<div data-spenso-math>" in html
assert "Momentum" not in html
# Native MathML uses logical mi text; the math font supplies the italic glyph.
indexed_math = ET.fromstring(re.search(r"<math\b[^>]*>.*?</math>", html, re.DOTALL).group())
# Canonical provenance aliases: edge 7 -> mu, half-edges 7/8 -> nu/rho, vertex 7 -> sigma.
expected_ports = {1: "ν", 2: "ρ", 3: "μ", 4: "σ"}
for math_root, math_html in ((indexed_math, html), (fraction_math, fraction_html)):
    q_labels = [
        node for node in math_root.findall(".//{*}mmultiscripts")
        if [token.text for token in node[0].iter() if token.tag.rsplit("}", 1)[-1] == "mi"] == ["q"]
    ]
    assert len(q_labels) == indexed.rank == 4, math_html
    q_ports = {}
    for label in q_labels:
        assert len(label) >= 5 and len(label) % 2 == 1, math_html
        momentum_id = int("".join(label[1].itertext()))
        assert label[2].tag.rsplit("}", 1)[-1] == "none", math_html
        assert label[-2].tag.rsplit("}", 1)[-1] == "none", math_html
        ports = [token.text for script in list(label)[3:] for token in script.iter()
                 if token.tag.rsplit("}", 1)[-1] == "mi"]
        assert ports == [expected_ports[momentum_id]], math_html
        assert momentum_id not in q_ports, math_html
        q_ports[momentum_id] = ports[0]
    assert q_ports == expected_ports, math_html
# Validate the portable source project without invoking the optional compiler.
project = indexed.to_typst_project()
assert {"main.typ", "render.typ", "notation.typ", "tree.cbor"} <= project.keys()
assert all(isinstance(value, bytes) and value for value in project.values())
for output in (html, indexed.to_typst()):
    assert "SourceIndex" not in output
    assert "SinkIndex" not in output
    assert "EdgeDummy" not in output
    assert "VertexDummy" not in output
assert scalar_graphs[0].projector_expression().to_expression() == 1
assert scalar_graphs[0].numerator_prefactor_expression().to_expression() == numerator_prefactor
assert scalar_graphs[0].overall_factor_expression() == feynkit_expression
scalar_graphs[0].validate()
assert type(free_tensor_reduced) is spenso.TensorExpression
assert free_tensor_reduced.rank == 2
assert free_tensor_reduced.axes == free_tensor_diagram.numerator_expression().axes
"#,
            )
            .unwrap();
            py.run(&assertions, Some(&locals), Some(&locals))?;

            remove_wrapper(py, &community, "symbolica.community.hepkit")?;
            remove_wrapper(py, &community, "symbolica.community.tensor")?;
        }
        Ok(())
    })
    .unwrap();
}
