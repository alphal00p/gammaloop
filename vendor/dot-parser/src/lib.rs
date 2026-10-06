#![warn(missing_docs)]

//! This crate is a parser for [DOT/Graphviz](https://www.graphviz.org) files.
//! The parser strictly sticks to the grammar specified
//! [here](https://www.graphviz.org/doc/info/lang.html). Notice that the
//! official implementation of Graphviz tools *does not* strictly follows this
//! grammar. However all graphs following the grammar are accepted by official
//! tools.
//!
//! The crate contains two modules: [`ast`] and [`canonical`]. Both contain a
//! `Graph` structure that represent a DOT graph. In the [`ast`] module, the
//! [`Graph`][ast::Graph] structure matches the abstract syntax tree of a given
//! DOT graph. However, in practice, such graph is not easy to work with.
//! Therefore, the [`Graph`][canonical::Graph] structure in the [`canonical`]
//! module aims to provide a canonic representation of the graph, much easier to
//! work with.
//!
//! ## [canonical::Graph]
//!
//! Unless you have a good reason to use [ast::Graph], I would recommend using
//! [canonical::Graph]. For instance, since [ast::Graph] follows closely the
//! [Graphviz grammar](https://www.graphviz.org/doc/info/lang.html), iterating
//! over all edges is quite complex: one has to take into account chained edges
//! (e.g. `A -> B -> C` is a single [ast::EdgeStmt], but it represents two
//! edges: `A -> B` and `A -> C`), as well as subgraphs (e.g. `A -> subgraph {B
//! C}` contains two edges: `A -> B` and `A -> C`, and subgraph can be nested).
//! On the other hand, [canonical::Graph] flattens all those edges for you,
//! which makes iterating over the edges extremelly trivial: 
//!
//! ```
//! use dot_parser::*;
//! let graph = canonical::Graph::from(
//!     ast::Graph::try_from("graph { A -> subgraph { B C } }").unwrap()
//!     );
//! for edge in graph.edges.set {
//!     println!("{} -> {}", edge.from, edge.to);
//! }
//! ```
//!
//! ## [ast::Graph]
//! The main structure of this module is [`ast::Graph`]. It implements
//! `TryFrom<&str>`, which allows one to parse a [`&str`]. It also provides a
//! function [`from_file`][ast::Graph::from_file] to read a graph directly from
//! a file. The fields (and subfields) of the structure coincide with rules of
//! the grammar.
//!
//! ## Feature flags
//!
//! This crate offers the following feature flags (all disable by default):
//!  - display: if enabled, most structures implement
//!  [Display][std::fmt::Display], as long as the attribute type implements
//!  [Display][std::fmt::Display]. Graphs can be printed (following the
//!  DOT/Graphviz syntax).
//!  - petgraph: if enabled, [Graph][canonical::Graph] can be converted into
//!  [petgraphs'
//!  graphs](https://docs.rs/petgraph/latest/petgraph/graph/struct.Graph.html).
//!  Essentially, enabling this flag provides the following implementation: 
//!  ```no_compile
//! impl<A> From<Graph<A>> for PetGraph<Node<A>, AList<A>> { ... }
//!  ```
//!  - to_tokens: if enabled, [Graph][ast::Graph]s implement quote's
//!  [ToTokens](https://docs.rs/quote/latest/quote/trait.ToTokens.html) trait,
//!  allowing one to easily generate graphs in macros. This is mainly used in
//!  [dot_parser_macros](https://docs.rs/dot-parser-macros/0.2.0/dot_parser_macros/).

pub mod ast;
pub mod canonical;
mod subgraph_free;
