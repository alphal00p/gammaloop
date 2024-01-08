use symbolica::{
    representations::{AsAtomView, Atom},
    state::{State, Workspace},
};

use crate::{graph::Graph, model::Model};

pub fn generate_numerator(graph: &Graph, model: &Model, state: &mut State, ws: &Workspace) -> Atom {
    let _ = model;
    let vatoms: Vec<Atom> = graph
        .vertices
        .iter()
        .flat_map(|v| v.apply_vertex_rule(state, ws))
        .collect();
    // graph
    //     .edges
    //     .iter()
    //     .filter(|e| e.edge_type == EdgeType::Virtual)
    //     .map(|e| e.particle);
    let mut builder = vatoms[0].builder(state, ws);

    for v in vatoms[1..].iter() {
        builder = builder * v;
    }

    builder.into_atom()
}
