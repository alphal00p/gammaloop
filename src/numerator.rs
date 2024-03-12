use symbolica::representations::Atom;

use crate::{
    graph::{EdgeType, Graph},
    model::Model,
};

pub fn generate_numerator(graph: &Graph, model: &Model) -> Atom {
    let _ = model;
    let vatoms: Vec<Atom> = graph
        .vertices
        .iter()
        .flat_map(|v| v.apply_vertex_rule(graph))
        .collect();
    // graph
    //     .edges
    //     .iter()
    //     .filter(|e| e.edge_type == EdgeType::Virtual)
    //     .map(|e| e.particle);
    let eatoms: Vec<Atom> = graph
        .edges
        .iter()
        .filter(|e| e.edge_type == EdgeType::Virtual)
        .map(|e| e.numerator(graph))
        .collect();
    let mut builder = Atom::new_num(1);

    for v in vatoms[1..].iter() {
        builder = builder * v;
    }

    for e in eatoms.iter() {
        builder = builder * e;
    }

    builder
}
