= Changelog

All notable changes to this project will be documented in this file.

The format is based on #link("https://keepachangelog.com/en/1.0.0/")[Keep a Changelog],
and this project adheres to #link("https://semver.org/spec/v2.0.0.html")[Semantic Versioning].

== Unreleased

=== Fixed

- Cut enumeration now preserves complete node sides and separates both terminal regions
  symmetrically. Every connected component must contain a terminal of its side; dangling
  half-edges remain with their nodes and bridges appear in the boundary only when they cross
  between sides. Disconnected components without terminals and live isolated nodes yield no cuts.

=== Added

- Generation now separates reusable `Theory` data, `Legs`, `Query`, compiled `Plan` and
  execution `Run`. `Plan::count` returns exact totals without constructing output graphs;
  skeletons can be retained and dressed separately. `explain()` reports where constraints act.
- `prune_skeletons` and `prune_dressings` receive typed hereditary views; `select` receives
  only complete diagrams. The Python bindings expose the same staged engine, batch hooks
  on the calling thread, and expire callback views on return. Explicit graph snapshots remain
  independent. Replace the previous generator constructor and graph generation method with
  a theory query and plan execution.

- Canonical forms: `HedgeGraph::canonical_form_by` and `canonical_form` compute the
  `CanonicalForm` of a graph (equal exactly for isomorphic graphs under the same keys) without
  labelling maps, automorphism generators or count. `Canonization::canonical_form` returns the
  same form and generated `Diagram`s compute theirs lazily with `canonical_form()`;
  `CanonicalForm::into_parts` moves out its words and key tables.

- Generation budgets: `Query::weight` bounds sums of integer signature weights (coupling
  and counterterm orders) and `Query::count` with `Countable` the vertices, propagators,
  chords and closed lines of some fields. Field dressing prunes with bounds kept per automaton
  state, the vertex rules are checked exactly, and loop orders that no vertex multiset meets and
  requests whose legs carry a charge or parity that no vertex balances are not searched. Budgets
  that every graph meets cost nothing. New `GenerationError::VertexWeights` and
  `GenerationError::TooManyCounts`.

=== Changed

- Read DOT with Linnet's own dependency-free reader instead of the GPL-2.0-or-later
  `dot-parser` crate. Parsing takes 10 to 16 times fewer instructions, a whole import 5 to 11
  times fewer, and a syntax error reports its line and column. `HedgeParseError` loses its
  lifetime parameter, and its `GraphFromFile` and `ParseError` variants become
  `Syntax(DotSyntaxError)` and `Io(std::io::Error)`. `CompassPt` is now
  `linnet::parser::CompassPt`, with `FromStr` and `Display`.
- `DotGraph::from_string`, `from_file` and `from_string_with_figment` reject input after the
  first graph instead of ignoring it, including a `;` after its closing brace, which the DOT
  grammar does not allow.
- Keywords are case-insensitive, as in Graphviz: `Node [x=1]` is a node default, and
  `Node -> b` is an error.
- Nodes are keyed by their decoded identifier, so `a` and `"a"` are one node, and repeated node
  statements merge their attributes instead of the last statement replacing the others.
  Quoted port names are decoded too: `a:"3"` requests half-edge 3, like `a:3`. Inferred node
  and edge indices follow the decoded identifiers, so they can change in a graph that mixes
  quoted and bare identifiers, or that quotes identifiers holding a character that sorts
  before `"`, such as a space or `!`.
- The statements of an endpoint subgraph, including its attributes and inner edges, are kept,
  and its inner edges come before the edges to and from it. This can change the inferred
  indices of parallel edges.
- `ID = ID` statements and `graph [...]` attributes set graph attributes in source order, so
  the later one wins, as in Graphviz; previously `ID = ID` always won.
- Quoted strings decode `\\` to one backslash and drop a backslash-newline continuation. The
  writer escapes backslashes where needed, so the `__linnet_hex_id_` envelope for strings
  ending in a backslash is gone. DOT written by earlier versions, whose writer escaped only
  quotes, reads differently when a value holds a backslash: a written `\\` becomes one
  backslash, a backslash before a quote (written `\\"`) ends the string early, and envelopes
  stay literal text. No DOT file in this repository is affected.
- Subgraphs nest at most 256 levels deep.
- `DotHedgeData::from(Option<String>)` keeps the statement as given instead of stripping
  quotes from it.

=== Removed

- `DotHedgeData::with_port`, `GlobalData`'s `TryFrom` implementation and the public
  `from_parser` constructors of `DotVertexData` and `DotEdgeData`, which took `dot-parser`
  types.

=== Fixed

- Ports that start with a compass letter, such as `a:ne`, `b:sw` or `a:source`, are no longer
  read as a compass point plus a stray node.
- `write_fmt`, `debug_dot` and `dot_of` no longer corrupt multi-line attribute values.
- Graphviz input that `dot-parser` rejected is accepted: strings ending in `\\`, nested HTML
  strings, non-ASCII bare identifiers, `#` comments, `[]` and uppercase keywords. `--` in a
  digraph, `->` in a graph and keywords used as identifiers are rejected, as in Graphviz.

== #link("https://github.com/alphal00p/linnet/compare/linnet-v0.16.1...linnet-v0.17.0")[0.17.0] - 2026-01-28

=== Other

- force-directed layout
- make the powerset iterator properly size the subset
- upgrade to all spanning _forests_ generation
- Refactor PowersetIterator to support generic ID type
- double ended iterator for subsets
- remove unused crates

== #link("https://github.com/alphal00p/linnet/compare/linnet-v0.16.0...linnet-v0.16.1")[0.16.1] - 2026-01-13

=== Other

- Quote string values in graph and parser outputs
- revert edge, node and hedge index order to be direct (not based on sort)

== #link("https://github.com/alphal00p/linnet/compare/linnet-v0.15.1...linnet-v0.16.0")[0.16.0] - 2026-01-08

=== Other

- Add methods for identifying bridges and non-bridges in subgraphs
- Make `trace_unfold` implementable outside of linnet

== #link("https://github.com/alphal00p/linnet/compare/linnet-v0.15.0...linnet-v0.15.1")[0.15.1] - 2026-01-05

=== Fixed

- fixed ordering problem in merge function and make ordering based on sorting not on direct map

=== Other

- update quote stripping
- Add agent context files and update project tooling
- Add graph algorithms module with topological sort and transitive ops

== #link("https://github.com/alphal00p/linnet/compare/linnet-v0.14.1...linnet-v0.15.0")[0.15.0] - 2025-12-08

=== Other

- Refactor Figment handling and improve data merging logic
- _(clinnet)_ release v0.1.8

== #link("https://github.com/alphal00p/linnet/compare/linnet-v0.14.0...linnet-v0.14.1")[0.14.1] - 2025-12-03

=== Other

- Add parallel processing and directional constraints
- update test
- only swap one edge when connecting identities

== #link("https://github.com/alphal00p/linnet/compare/linnet-v0.13.3...linnet-v0.14.0")[0.14.0] - 2025-11-28

=== Other

- change output name and make it short, and add directional groups

== #link("https://github.com/alphal00p/linnet/compare/linnet-v0.13.2...linnet-v0.13.3")[0.13.3] - 2025-11-24

=== Other

- update symbolica

== #link("https://github.com/alphal00p/linnet/compare/linnet-v0.13.1...linnet-v0.13.2")[0.13.2] - 2025-11-12

=== Other

- update to symbolica 0.20
- better templates
- clinnet for nested dots to pdf
- Add Figment integration for configuration handling

== #link("https://github.com/alphal00p/linnet/compare/linnet-v0.13.0...linnet-v0.13.1")[0.13.1] - 2025-11-05

=== Other

- remove prints

== #link("https://github.com/alphal00p/linnet/compare/v0.11.0...v0.12.0")[0.12.0] - 2025-08-18

=== Other

- all tests pass
- clippy fixes
- write global data outside of `graph[]` scope
- Add ID=ID statements to global data when parsing
- newline after nodes
- dot serialize with name
- write graph name
- correctly parse orientations
- make `iter_edges` iterate with `edge_id` order not involution order

== #link("https://github.com/alphal00p/linnet/compare/v0.10.0...v0.11.0")[0.11.0] - 2025-07-29

=== Other

- serialize global data

== #link("https://github.com/alphal00p/linnet/compare/v0.9.1...v0.10.0")[0.10.0] - 2025-07-29

=== Other

- Add check graph function and extracting by nodes
- Add `is_empty` to Swap trait and update implementations

== #link("https://github.com/alphal00p/linnet/compare/v0.9.0...v0.9.1")[0.9.1] - 2025-07-28

=== Other

- impl asref for newtypes

== #link("https://github.com/alphal00p/linnet/compare/v0.8.0...v0.9.0")[0.9.0] - 2025-07-27

=== Other

- Allow specifying a subgraph in the DOT format.
- Rewrite DOT output to include node, hedge and edge indices and fix bidirectional serialization
- DotGraph has concrete data. Parser now supports global edge and vertex data. Can set edge, vertex and hedge index
- add multi-graph parsing and global graph data

== #link("https://github.com/alphal00p/linnet/compare/v0.7.1...v0.8.0")[0.8.0] - 2025-07-14

=== Other

- clippy fix
- Rename HedgeVec to EdgeVec and introduce HedgeVec
- Add H parameter to HedgeGraph in `iter_edges` and orientation
- add H to SignedCycle
- rename `as_ref` to `to_ref`
- Add `hedge_data` to `hedge_graph`
- add perm by key
- update symbolica
- Fix tree traversal and modify ForestNodeStore API to use Option

== #link("https://github.com/alphal00p/linnet/compare/v0.7.0...v0.7.1")[0.7.1] - 2025-05-24

=== Other

- Fix crown computation by inverting subgraph inclusion check and add inclusion for HedgePair

== #link("https://github.com/alphal00p/linnet/compare/v0.6.3...v0.7.0")[0.7.0] - 2025-05-23

=== Other

- unify iterator signature and remove redundant and misleading functions
- Add zenodo doi

== #link("https://github.com/alphal00p/linnet/compare/v0.6.2...v0.6.3")[0.6.3] - 2025-05-23

=== Other

- improve README

== #link("https://github.com/alphal00p/linnet/compare/v0.6.1...v0.6.2")[0.6.2] - 2025-05-22

=== Other

- remove useless check and set
- Correctly forget identification history and adds subgraph contraction

== #link("https://github.com/alphal00p/linnet/compare/v0.6.0...v0.6.1")[0.6.1] - 2025-05-22

=== Other

- Merge pull request \#11 from `alphal00p/jules_wip_3253120278989361716`

== #link("https://github.com/alphal00p/linnet/compare/v0.5.2...v0.6.0")[0.6.0] - 2025-05-22

=== Fixed

- fix all tests and update symbolica
- fixes `forget_identification_history`, If you delete a subgraph you need

=== Other

- Add crown functions and implement spanning tree enumeration algo based on graph contraction
- PartialEq and Eq for hedgevec
- `from_raw` for hedgevec

== #link("https://github.com/alphal00p/linnet/compare/v0.5.1...v0.5.2")[0.5.2] - 2025-05-09

=== Fixed

- fix involution extraction

=== Other

- Add Permutation iterators, new graph extract test

== #link("https://github.com/alphal00p/linnet/compare/v0.5.0...v0.5.1")[0.5.1] - 2025-05-05

=== Other

- Clean up debug prints and add cycle test

== #link("https://github.com/alphal00p/linnet/compare/v0.4.0...v0.5.0")[0.5.0] - 2025-05-05

=== Fixed

- fixes all around

=== Other

- boundary->crown
- hairs->boundary

== #link("https://github.com/alphal00p/linnet/compare/v0.3.0...v0.4.0")[0.4.0] - 2025-04-30

=== Added

- feature flag bincode and serde

== #link("https://github.com/alphal00p/linnet/compare/v0.2.3...v0.3.0")[0.3.0] - 2025-04-30

=== Fixed

- fixed extraction for childvec forest, all tests pass and fix warnings

=== Other

- adds extract to SmartHedgeVec and involution
- maps now take a FnMut
- derive more (Decode and Encode)

== #link("https://github.com/alphal00p/linnet/compare/v0.2.2...v0.2.3")[0.2.3] - 2025-04-08

=== Other

- unpin serde

== #link("https://github.com/alphal00p/linnet/compare/v0.2.1...v0.2.2")[0.2.2] - 2025-04-01

=== Fixed

- fix test snapshots

=== Other

- raw image

== #link("https://github.com/alphal00p/linnet/compare/v0.2.0...v0.2.1")[0.2.1] - 2025-04-01

=== Other

- Update README.md, logo is a github asset
- put crates badge

== #link("https://github.com/alphal00p/linnet/compare/v0.1.0...v0.2.0")[0.2.0] - 2025-04-01

=== Fixed

- fix cetz orientation

=== Other

- Create release-plz.yml
