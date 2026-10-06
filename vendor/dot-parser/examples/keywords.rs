
use dot_parser::*;

fn main() {
    let raw = ast::Graph::try_from(
        "// This graph uses identifiers that begin with keywords.
        digraph { 
node_1 -> node_2
}",
    );
    match raw {
        Ok(graph) => {
            //println!("{:#?}", graph);
            println!("{:#?}", canonical::Graph::from(graph));
        }
        Err(e) => {
            println!("{}", e);
        }
    }
}
