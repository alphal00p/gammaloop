// The goal of this example is to illustrate the difference between ast::Graph
// and canonical::Graph. In practice, unless you need to work with ast::Graph,
// for a very good reason, canonical::Graph will make your life easier.

use dot_parser::ast;
use dot_parser::canonical;

static DOT_STR: &str = "graph {
    A [label=A][color=grey]
    B [label=B][color=grey]
    C [label=C]
    D [label=D]

    A -- B -- subgraph { C D }
}";

fn main() {
    let ast_graph = ast::Graph::try_from(DOT_STR).unwrap();
    let canonical_graph = canonical::Graph::from(ast_graph.clone());

    // Working with ast_graph is painful, because the structure follows closely
    // the syntax: for instance, the `ast::Graph` structure has a `stmts` field,
    // which contains a `NodeStmt` (for A) which contains an `AttrList` field, which
    // contains two `AList`s, each containing a single attribute, because the two attributes were given in separate
    // lists (`[label=A][color=grey]` instead of `[label=A; color=grey]`).
    //
    // Printing all attributes is painful:
    
    println!("Showing all attributes using ast::Graph:");
    for stmt in ast_graph.stmts {
        if let ast::Stmt::NodeStmt(n) = stmt {
            if let Some(attr_list) = n.attr {
                for alist in attr_list {
                    for (attr, val) in alist {
                        println!("({}, {})", attr, val);
                    }
                }
            }
        }
    }

    // Fortunatelly, we can work with a canonical graph, which converts the
    // ast::Graph into something that follows more closely our intuition:
    println!("\nShowing all attributes using canonical::Graph:");
    for node in canonical_graph.nodes.set.values() {
        for (attr, val) in &node.attr {
            println!("({}, {})", attr, val);
        }
    }
    

    // Similarly, the grammar allows to chain edges, and to have subgraphs in
    // chains:
    // (`A -- B -- subgraph { C D }`). 
    //
    // This is equivalent to writing explicitly the three following edges: 
    // ```
    // A -- B
    // B -- C
    // B -- D
    // ```
    //
    // While working with the graph, it would be easier to deal with the list of
    // explicit edges, instead of working the chained edges.
    // Suppose we want to print all the edges, doing it with ast::Graph is so
    // painful that I don't even want to do it here: we first need to split each
    // step of the EdgeStmt into separate edges, then look at the source and the
    // destination, and check whether they are a subgraph or not, and if it is, get
    // all the nodes from the subgraph (for both the source and the destination)
    // and create an edge for each of them... this is super annoying in
    // practice (in particular, think about nested subgraphs, and dealing with
    // edges within subgraphs themselves).

    // Fortunatelly, canonical::Graph deals with this complexity:
    println!("\nShowing all edges using canonical::Graph:");
    for edge in canonical_graph.edges.set {
        println!("{} -- {}", edge.from, edge.to);
    }
}
