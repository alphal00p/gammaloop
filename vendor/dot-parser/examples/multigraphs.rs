use dot_parser::*;

fn main() {
    let raw = ast::Graphs::try_from(
        "// A  string with two graphs: 
        digraph { 
            A -> B
        }

        graph G {
            G1 -> G2
        }
        ",
    );
    match raw {
        Ok(graphs) => {
            for g in graphs.graphs {
                println!("{:#?}", canonical::Graph::from(g));
            }
        }
        Err(e) => {
            println!("{}", e);
        }
    }
}
