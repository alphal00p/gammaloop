use std::{ffi::CString, sync::Arc};

use feynkit_graph::{DiagramEdge, DiagramVertex, ExternalLeg, ExternalState, FeynmanDiagram};
use feynkit_model::Model;
use feynkit_py::PyFeynmanDiagram;
use pyo3::{
    prelude::*,
    types::{PyDict, PyList, PyModule},
};
use symbolica::api::python::{SymbolicaCommunityModule, create_symbolica_module};

#[test]
fn cff_explorer_keeps_native_families_and_graph_identity() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let core = PyModule::new(py, "symbolica")?;
        create_symbolica_module(&core)?;
        // The display imports tensor formatters; use this binary's Symbolica types.
        let community = PyModule::new(py, "symbolica.community")?;
        let tensor = PyModule::new(py, "symbolica.community.tensor")?;
        spynso3::SpensoModule::register_module(&tensor)?;
        core.add("community", &community)?;
        community.add("tensor", &tensor)?;
        let modules = py.import("sys")?.getattr("modules")?;
        modules.set_item("symbolica", &core)?;
        modules.set_item("symbolica.community", &community)?;
        modules.set_item("symbolica.community.tensor", &tensor)?;
        let fk = PyModule::new(py, "feynkit")?;
        feynkit_py::initialize_feynkit(&fk)?;
        let model =
            Arc::new(Model::from_json(include_str!("fixtures/scalars_2p_3p.json")).unwrap());
        let rule = model.vertex_rule_id("V_3_SCALAR_000").unwrap();
        let diagrams = PyList::empty(py);
        for loops in 0..=3 {
            let mut builder = FeynmanDiagram::builder(
                Arc::clone(&model),
                if loops == 0 {
                    "ladder-0</script><script>throw 'unsafe'</script>".to_owned()
                } else {
                    format!("{loops}-loop ladder")
                },
            );
            let nodes: Vec<_> = (0..2 * (loops + 1))
                .map(|id| builder.add_vertex(DiagramVertex::interaction(format!("v{id}"), rule)))
                .collect();
            let edge = || DiagramEdge::new(model.particle_id("scalar_0").unwrap(), false);
            for column in 0..=loops {
                builder
                    .add_edge(nodes[column], nodes[column + loops + 1], edge())
                    .unwrap();
                if column < loops {
                    builder
                        .add_edge(nodes[column], nodes[column + 1], edge())
                        .unwrap();
                    builder
                        .add_edge(nodes[column + loops + 1], nodes[column + loops + 2], edge())
                        .unwrap();
                }
            }
            for (index, node) in [0, loops + 1, loops, 2 * loops + 1].into_iter().enumerate() {
                let mut leg = edge();
                leg.external = Some(ExternalLeg {
                    name: format!("p{index}"),
                    index,
                    state: if index < 2 {
                        ExternalState::Incoming
                    } else {
                        ExternalState::Outgoing
                    },
                    connection: index,
                });
                if index < 2 {
                    builder.add_edge(None, nodes[node], leg).unwrap();
                } else {
                    builder.add_edge(nodes[node], None, leg).unwrap();
                }
            }
            diagrams.append(Py::new(
                py,
                PyFeynmanDiagram::from(builder.build().unwrap()),
            )?)?;
        }
        let locals = PyDict::new(py);
        locals.set_item("fk", fk)?;
        locals.set_item("core", core)?;
        locals.set_item("diagrams", diagrams)?;
        py.run(
            &CString::new(include_str!("cff_display.py")).unwrap(),
            Some(&locals),
            Some(&locals),
        )?;
        py.run(
            c"check_cff_displays(fk, diagrams)",
            Some(&locals),
            Some(&locals),
        )
    })
    .unwrap();
}
