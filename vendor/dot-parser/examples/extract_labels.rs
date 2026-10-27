// The goal of this example is to obtain the list of nodes and the list of edges of a graph, and to
// associate them with the "label" (if defined).
//
// For instance, if we are given the following graph:
//
// graph {
//  A[label="A"];
//  B;
//  C[label="C", color="grey"];
//  A -> B -> C [label="Edge"];
// }
//
// We want to obtain the lists (not correct, just to give an idea):
//  - nodes: Vec<(Node, Option<String>)> = [(A, Some("A")), (B, None), (C, Some("C"))]
//  - edges: Vec<(Edge, Option<String>)> = [(A -> B, Some("Edge")), (B -> C, Some("Edge"))]

use dot_parser::*;

fn extract_label(pair: (&str, &str)) -> Option<String> {
    if let "label" = pair.0 {
        Some(pair.1.to_string())
    } else {
        None
    }
}

fn main() {
    let graph = ast::Graph::try_from(
        "graph {
 A[label=\"label 1 node A\", label=\"label 2 node A\"];
 B;
 C[label=\"label node C\", color=\"grey\"];
 A -> B -> C [label=\"label edges\"];
}",
    )
    .unwrap();
    let graph = graph.filter_map(&extract_label);
    let graph = canonical::Graph::from(graph);

    for node in graph.nodes.set.values() {
        println!("{} -- {}", node.id, node.attr)
    }

    for edge in graph.edges.set {
        println!("{} -> {} -- {}", edge.from, edge.to, edge.attr)
    }
}
