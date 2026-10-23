use pyo3::{prelude::*, types::PyDict};
use symbolica::api::python::{SymbolicaCommunityModule, create_symbolica_module};

#[test]
fn propagator_annotations_keep_native_wrapper_and_denominator_semantics() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let core = PyModule::new(py, "annotation_core")?;
        create_symbolica_module(&core)?;
        let hep = PyModule::new(py, "annotation_hep")?;
        feynkit_py::FeynkitModule::register_module(&hep)?;
        feynkit_py::FeynkitModule::initialize(py)?;
        let tensor = PyModule::new(py, "annotation_tensor")?;
        spynso3::SpensoModule::register_module(&tensor)?;
        spynso3::SpensoModule::initialize(py)?;
        let locals = PyDict::new(py);
        locals.set_item("core", core)?;
        locals.set_item("hep", hep)?;
        locals.set_item("tensor", tensor)?;
        py.run(
            c"model = hep.Model.phi3()
diagram = hep.FeynmanDiagram.from_dot(model, '''digraph bubble {
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
