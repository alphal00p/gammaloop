#import "../figure.typ": render
#import "../../../../../crates/linnest/typst/src/lib.typ": graph, subgraph

// Cross-section fixture: e+ e- -> t tbar with a physical cut.
// Follow FeynmanDiagram::to_linnest: retain the native graph and flatten its
// authenticated paired-external metadata before opening only the initial state.
#let native = graph.parse(
  read("../../../../../tests/resources/graphs/epemttbar.dot"),
).first()
#let native = graph.map(native, edge: edge => {
  let external = json(bytes(edge.fields.at("external", default: "null")))
  if external == none { return none }
  (
    is_cut: external.connection,
    external-state: external.state,
    external-index: external.index,
    external-name: external.name,
    external-connection: external.connection,
  )
})
#let initial = graph.edges(native).filter(edge => (
  edge.source != none and edge.sink != none
    and type(edge.data) == dictionary and edge.data.at("is_cut", default: none) != none
)).map(edge => edge.edge)
#let opened = graph.cut(
  native,
  left: subgraph.select(native, sink: initial),
  right: subgraph.select(native, source: initial),
)

#render(opened, cross-section-mode: true, momentum-arrows: true)
