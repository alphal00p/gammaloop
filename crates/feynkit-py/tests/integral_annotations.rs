use feynkit_py::FeynkitModule;
use pyo3::types::{PyDict, PyDictMethods, PyModule};
use pyo3::{PyResult, Python};
use spynso3::SpensoModule;
use symbolica::api::python;
use symbolica::api::python::SymbolicaCommunityModule;

#[test]
fn propagator_annotations_keep_native_wrapper_and_denominator_semantics() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let core = PyModule::new(py, "annotation_core")?;
        python::create_symbolica_module(&core)?;
        let hep = PyModule::new(py, "annotation_hep")?;
        FeynkitModule::register_module(&hep)?;
        FeynkitModule::initialize(py)?;
        let tensor = PyModule::new(py, "annotation_tensor")?;
        SpensoModule::register_module(&tensor)?;
        SpensoModule::initialize(py)?;
        let locals = PyDict::new(py);
        locals.set_item("core", core)?;
        locals.set_item("hep", hep)?;
        locals.set_item("tensor", tensor)?;
        py.run(
            c"model = hep.Model.phi3()
diagram = hep.FeynmanDiagram.from_dot(model, '''digraph bubble {
  graph [num=\"(x+y)^2*(x-y)^3\"];
  ext [style=invis];
  ext -> a [particle=\"phi\"];
  a -> b [particle=\"phi\", lmb_id=0, power=2];
  a -> b [particle=\"phi\", power=3];
  b -> ext [particle=\"phi\"];
}''')
assert [edge.propagator_power for edge in diagram.internal_edges] == [2, 3]
assert isinstance(diagram.denominator_expression(), tensor.TensorExpression)
product = core.E('1')
for edge in diagram.internal_edges:
    base = edge.denominator_expression(power=1).to_expression()
    recorded = edge.denominator_expression().to_expression()
    assert recorded == base ** edge.propagator_power
    assert edge.denominator_expression(power=-1).to_expression() == base ** -1
    product *= recorded
assert diagram.denominator_expression().to_expression() == product
restored = hep.FeynmanDiagram.from_dot(model, diagram.to_dot())
assert [edge.propagator_power for edge in restored.internal_edges] == [2, 3]
assert restored.loop_momentum_basis.loop_edges == diagram.loop_momentum_basis.loop_edges
assert restored.denominator_expression() == diagram.denominator_expression()
assert '<svg' in diagram.render().to_svg()
for options in ({'method': 'cff'}, {'method': 'cff', 'numerator': core.E('1')}, {'method': 'ltd'}):
    try:
        diagram.integrate_energy(**options)
    except ValueError as error:
        assert 'stored propagator power' in str(error)
    else:
        raise AssertionError('raised powers were integrated as unit powers')
sewn = hep.FeynmanDiagram.from_dot(model, '''digraph cut_bubble {
  graph [final_state=\"phi,phi\"];
  ext [style=invis];
  ext -> a [particle=\"phi\", is_cut=9];
  b -> ext [particle=\"phi\", is_cut=9];
  a -> b [particle=\"phi\", lmb_id=0, power=2];
  b -> a [particle=\"phi\"];
}''')
cut = sewn.cuts[0]
stored = {edge.id: edge.propagator_power for edge in cut.edges}
assert 2 in stored.values()
default_factors = [factor.to_expression() for factor in cut.propagators()]
assert default_factors == [factor.to_expression() for factor in cut.propagators(edge_powers=stored)]
assert default_factors != [factor.to_expression() for factor in cut.propagators(edge_powers={edge: 1 for edge in stored})]
try:
    diagram.internal_edges[0].propagator_power = 5
except AttributeError:
    pass
else:
    raise AssertionError('the graph inspection record is writable')
",
            Some(&locals),
            Some(&locals),
        )
    })
    .unwrap();
}
