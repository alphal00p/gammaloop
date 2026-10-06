
use dot_parser::*;

fn main() {
    let raw = ast::Graph::try_from(
        "digraph { // A simple graph with two nodes and a single edge
A -> /* The edge from A to the subgraph {B -> C} */ subgraph g { B -> C }
}",
    );
    match raw {
        Ok(graph) => {
            println!("{:#?}", graph);
            println!("{:#?}", canonical::Graph::from(graph));
        }
        Err(e) => {
            println!("{}", e);
        }
    }
}
