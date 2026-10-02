use std::{ffi::CString, sync::Arc};

use feynkit_graph::{DiagramEdge, DiagramVertex, ExternalLeg, ExternalState, FeynmanDiagram};
use feynkit_model::Model;
use feynkit_py::PyFeynmanDiagram;

use pyo3::{
    prelude::*,
    types::{PyDict, PyModule},
};
use symbolica::api::python::create_symbolica_module;

#[test]
fn generalized_cuts_keep_derivative_actions_and_validate_coordinates() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let core = PyModule::new(py, "core")?;
        create_symbolica_module(&core)?;
        let feynkit = PyModule::new(py, "feynkit")?;
        feynkit_py::initialize_feynkit(&feynkit)?;
        let model = Arc::new(Model::from_json(include_str!("fixtures/scalars_2p_3p.json")).unwrap());
        let rule = model.vertex_rule_id("V_3_SCALAR_000").unwrap();
        let mut builder = FeynmanDiagram::builder(Arc::clone(&model), "cut-algebra-bubble");
        let left = builder.add_vertex(DiagramVertex::interaction("left", rule));
        let right = builder.add_vertex(DiagramVertex::interaction("right", rule));
        let scalar = || DiagramEdge::new(model.particle_id("scalar_0").unwrap(), false);
        let mut incoming = scalar();
        incoming.external = Some(ExternalLeg { name: "p0".into(), index: 0, state: ExternalState::Incoming, connection: 0 });
        let mut outgoing = scalar();
        outgoing.external = Some(ExternalLeg { name: "p1".into(), index: 1, state: ExternalState::Outgoing, connection: 1 });
        builder.add_edge(None, left, incoming).unwrap();
        builder.add_edge(right, None, outgoing).unwrap();
        builder.add_edge(left, right, scalar()).unwrap();
        builder.add_edge(left, right, scalar()).unwrap();
        let diagram = Py::new(py, PyFeynmanDiagram::from(builder.build().unwrap()))?;
        let locals = PyDict::new(py);
        locals.set_item("core", core)?;
        locals.set_item("fk", feynkit)?;
        locals.set_item("diagram", diagram)?;
        let source = CString::new(
            r#"
q = core.Expression.symbol('cut_python::q')
cut = fk.CutPropagator(q, 2, power=3, normalization=1)
assert cut.apply(q**2, q) == -core.Expression.num(1)/128
assert cut.to_expression() != cut.to_expression(covariant=False)
reverse = fk.CutPropagator(q, 2, power=3, orientation=-1, normalization=1)
assert reverse.apply(q**2, q) == cut.apply(q**2, q)
for kwargs in ({'power': 0}, {'orientation': 0}, {'prescription': 2}):
    try:
        fk.CutPropagator(q, 2, **kwargs)
    except fk.CffError:
        pass
    else:
        raise AssertionError('invalid cut accepted')
try:
    fk.CutPropagator(q, q+1).apply(1, q)
except fk.CffError:
    pass
else:
    raise AssertionError('dependent on-shell energy accepted')
result = diagram.build_cff()
other = diagram.build_cff()
surface = result.surfaces[0]
assert result.to_expression() != result.to_expression(expand_surfaces=True)
assert result.to_expression(normalized=True) != result.to_expression(expand_surfaces=True)
assert 'feynkit::E' not in str(result.to_expression())
assert result.surface_expression(surface) != core.Expression.num(0)
try:
    other.surface_expression(surface)
except fk.CffError:
    pass
else:
    raise AssertionError('foreign surface accepted')
group = result.raised_surface_groups()[0]
assert group.max_order == 1
assert len(result.pole_coefficients(group)) == 1
try:
    other.pole_coefficients(group)
except fk.CffError:
    pass
else:
    raise AssertionError('foreign surface group accepted')
assert result.residue(group, variable=q, root=1, surface=q-1, coefficient=q**2) == core.Expression.num(1)
"#,
        )
        .unwrap();
        py.run(&source, Some(&locals), Some(&locals))
    })
    .unwrap();
}
