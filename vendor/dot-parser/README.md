# dot-parser

This crate provides a parser for [dot files](https://www.graphviz.org/doc/info/lang.html), also known as Graphviz.
The implementation *sticks closely* to the formal grammar given by Graphviz. In
addition, we provide transformation functions (i.e. an equivalent of `flat_map`
for lists) in order to manipulate graphs, as well as a canonical representation
of graphs, which flattens subgraphs and replaces edge statements with multiple right-hand sides with simple edge statements.

## Usage

In `Cargo.toml`:

```
dot-parser = "*"
```

In your Rust source file (to be adapted depending on your needs): 

```
use dot_parser::*;
```

Then, a file can be parsed using the `from_file` method of `dot_parser::ast::Graph`. Parsing a file can fail if there is an issue with opening and reading the file, or if the content of the file is not a valid Graphviz graph.

```
fn main() {
	let res = ast::Graph::from_file("/path/to/graph.dot");
	match res {
		Ok(graph) =>
			println!("{:#?}", graph),
		Err(e) => 
			println!("{}", e),
	}
}
```

Dealing with `dot_parser::ast::Graph`s is quite unconvenient, as the
representation sticks closely to the Graphviz grammar. For instance, we have to
consider subgraphs or lists of list of attributes, etc. Instead,
`dot_parser::canonical::Graph`s are easier to work with. 

```
fn main() {
    let res = ast::Graph::from_file("/path/to/graph.dot").unwrap();
    let canonical = canonical::Graph::from(res);
    println!("{:#?}", canonical);
}
```

## Features

### Outputting a dot graph

With the `display` feature, canonical graphs can be displayed (following the Graphviz grammar), as long as the type of attributes implements `Display`.

```
fn main() {
	let res = ast::Graph::from_file("/path/to/graph.dot").unwrap();
    let canonical = canonical::Graph::from(res).filter_map::<_, String>(&|_| None);
    println!("{}", canonical);
}
```

### Petgraph compatibility (up to 0.4.2)

Canonical graphs can be converted into [Petgraph's graphs](https://crates.io/crates/petgraph) with the `petgraph` feature.

```
fn main() {
	let res = ast::Graph::from_file("/path/to/graph.dot").unwrap();
    let canonical = canonical::Graph::from(res);
    let petgraph: petgraph::graph::Graph<_, _> = canonical.into();
    println!("{:#?}", canonical);
}
```

From version 0.5, this is no longer implemented. The implementation was not super flexible (it was implemented for a single version of petgraph). Instead, reading petgraph from dotfile will be implemented in petgraph directly.

## Compile time import

Graphs can be imported at compile time using the `dot_parser_macros` crate. The macro `from_dot_file!` allows to statically import a graph from a file, while `from_dot_string!` allows you to directly write the graph in the Rust source.

```
use dot_parser_macros::*;
use dot_parser::ast::Graph;

fn main() {
    let graph = from_dot_file!("/tmp/graph.dot");
    println!("{:#?}", graph);
}
```

## Limitations

- Officially, any XML string can be used as an ID, and any HTML string as a label. This is not correctly implemented now.
- Officially, lines begining with a '#' are considered to be preprocessor output and are ignored. This is not implemented now.

## Differences from DOT Formal Grammar

- The rule for `stmt_list` is `stmt_list = ( stmt ~ ";"? )*` instead of `stmt_list : [ stmt [ ';' ] stmt_list ]` due to performances issues (c.f [here](https://codeberg.org/bromind/dot-parser/issues/3)).

## Acknowledgement

This crate was written by 
[Martin Vassor](mailto:martin@vassor.org?subject=[dot-parser]), member of
the [Mobility Reading Group](https://mrg.cs.ox.ac.uk/) at the [University of
Oxford](https://www.ox.ac.uk/) (formerly at [Imperial College,
London](https://www.imperial.ac.uk/)).

Initial development of this work (up to v0.1.3) was supported by EPSRC
[EP/T006544/1](https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/T006544/1),
[EP/K011715/1](https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/K011715/1),
[EP/K034413/1](https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/K034413/1),
[EP/L00058X/1](https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/L00058X/1),
[EP/N027833/1](https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/N027833/1),
[EP/N028201/1](https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/N028201/1),
[EP/T014709/1](https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/T014709/1), 
[EP/V000462/1](https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/V000462/1),
and NCSS/EPSRC VeTSS.
