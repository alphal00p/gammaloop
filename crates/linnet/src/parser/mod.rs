//! # DOT interchange
//!
//! This module parses the [DOT language](https://graphviz.org/doc/info/lang.html)
//! into Linnet's half-edge model and serializes half-edge graphs back to DOT.
//! [`DotGraph::from_string`] and [`DotGraph::from_file`] parse one graph. A
//! [`DotGraph`] retains [`GlobalData`] alongside a [`HedgeGraph`] whose node,
//! edge, and half-edge payloads are [`DotVertexData`], [`DotEdgeData`], and
//! [`DotHedgeData`]. Use [`set::DotGraphSet`] for input containing multiple graphs.
//!
//! Parsed DOT attributes remain available in those payload types. In particular,
//! external-port flow and edge orientation are preserved in the parsed graph
//! rather than discarded as drawing-only attributes. [`DotGraph`] dereferences to
//! its inner [`HedgeGraph`] for graph algorithms.
//!
//! ## Parse, inspect, and serialize
//!
//! ```rust
//! use linnet::parser::DotGraph;
//!
//! let graph: DotGraph = DotGraph::from_string(
//!     r#"digraph example {
//!         incoming -> interaction [label="particle"];
//!     }"#,
//! )
//! .unwrap();
//!
//! assert_eq!(graph.n_nodes(), 2);
//! assert_eq!(graph.n_edges(), 1);
//! assert_eq!(graph.n_hedges(), 2);
//!
//! let mut output = String::new();
//! graph.write_fmt(&mut output).unwrap();
//! let round_trip: DotGraph = DotGraph::from_string(output).unwrap();
//! assert_eq!(round_trip.n_nodes(), graph.n_nodes());
//! assert_eq!(round_trip.n_edges(), graph.n_edges());
//! ```
//!
//! The [`crate::dot!`] macro is a shorthand for parsing an inline DOT token tree. To
//! serialize a graph with application-specific payloads, use
//! [`HedgeGraph::dot_serialize_io`] or [`HedgeGraph::dot_serialize_fmt`]. Those
//! methods take graph-level metadata plus mappings for half-edge, edge, and node
//! payloads; parsing does not implicitly convert DOT payloads into user types.
//!
//! [`HedgeGraph`]: crate::half_edge::HedgeGraph

use std::{
    collections::BTreeMap,
    fmt::{Debug, Write},
    ops::{Deref, DerefMut},
    path::Path,
};

#[cfg(feature = "serde")]
use figment;

use ahash::{HashSet, HashSetExt};
use dot_parser::ast::CompassPt;
use figment::Figment;
use indenter::CodeFormatter;
use itertools::{Either, Itertools};
use subgraph_free::SubGraphFreeGraph;

use crate::{
    half_edge::{
        builder::HedgeGraphBuilder,
        involution::{EdgeIndex, Hedge},
        nodestore::{DefaultNodeStore, NodeStorage, NodeStorageOps},
        subgraph::{ModifySubSet, SubGraphLike, SubSetLike},
        swap::Swap,
        GVEdgeAttrs, HedgeGraph, NodeIndex,
    },
    permutation::Permutation,
};

/// Strips surrounding quotes from a string if present
pub(crate) fn strip_quotes(s: &str) -> &str {
    if s.len() >= 2 {
        let chars: Vec<char> = s.chars().collect();
        if (chars[0] == '"' && chars[chars.len() - 1] == '"')
            || (chars[0] == '\'' && chars[chars.len() - 1] == '\'')
        {
            return &s[1..s.len() - 1];
        }
    }
    s
}

pub mod set;
pub use set::GraphSet;

#[cfg(feature = "rkyv")]
pub mod archive;
#[cfg(feature = "rkyv")]
pub use archive::{
    ArchivedDotEdgeEndpointsView, ArchivedDotEdgeView, ArchivedDotEndpointView,
    ArchivedDotGraphBytesSetView, ArchivedDotGraphView, ArchivedDotVertexView, DotGraphBytesSet,
};

pub mod global;
pub use global::GlobalData;

pub mod vertex;
pub use vertex::DotVertexData;

pub mod edge;
pub use edge::DotEdgeData;

pub mod hedge;
pub use hedge::DotHedgeData;

#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
#[cfg_attr(feature = "rkyv", archive(check_bytes))]
pub struct DotGraph<N: NodeStorage<NodeData = DotVertexData> = DefaultNodeStore<DotVertexData>> {
    /// Graph-level DOT metadata, including its name, payload, and attribute defaults.
    pub global_data: GlobalData,
    /// Parsed half-edge topology together with its DOT edge, vertex, and half-edge data.
    pub graph: HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData, N>,
}

impl<N: NodeStorage<NodeData = DotVertexData>> Deref for DotGraph<N> {
    type Target = HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData, N>;

    fn deref(&self) -> &Self::Target {
        &self.graph
    }
}

impl<N: NodeStorage<NodeData = DotVertexData>> DerefMut for DotGraph<N> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.graph
    }
}

#[derive(Debug, Clone)]
pub enum NodeIdOrDangling {
    Id(NodeIndex),
    Dangling {
        statements: BTreeMap<String, String>,
    },
}

mod subgraph_free;

impl<S: NodeStorageOps<NodeData = DotVertexData>> DotGraph<S> {
    pub fn compass_subgraph<Sub: ModifySubSet<Hedge> + SubSetLike>(
        &self,
        cps: Option<CompassPt>,
    ) -> Sub {
        let mut a: Sub = self.empty_subgraph();

        for (h, d) in self.graph.iter_hedges() {
            if cps == d.compasspt {
                a.add(h);
            }
        }
        a
    }

    pub fn write_io<W: std::io::Write>(&self, writer: &mut W) -> Result<(), std::io::Error> {
        writeln!(writer, "digraph {}{{", self.global_data.name)?;

        writeln!(writer, "{:4}", self.global_data)?;

        for (n, (_, _, v)) in self.iter_nodes().enumerate() {
            let mut node_data: DotVertexData = v.clone();
            node_data.remove_common(&self.global_data);

            if let Some(name) = &node_data.name {
                write!(writer, "\t{name}")?;
            } else {
                write!(writer, "\t{n}")?;
            }

            let data = node_data.to_string();
            if !data.is_empty() {
                writeln!(writer, " [{data}];")?;
            } else {
                writeln!(writer, ";")?;
            }
        }

        for (hedge_pair, eid, data) in self.iter_edges() {
            let mut edata = data.data.clone();
            edata.remove_common(&self.global_data);
            let attr = GVEdgeAttrs {
                color: None,
                label: None,
                other: Some(edata.to_string()),
            };

            // write!(writer, "  ")?;
            hedge_pair.add_data(self).dot_io(
                writer,
                self,
                eid,
                |h| h.statement.clone(),
                |a| self[a].name.clone().unwrap_or(a.to_string()),
                data.orientation,
                attr,
            )?;
        }
        writeln!(writer, "}}")?;
        Ok(())
    }

    pub fn write_fmt<W: std::fmt::Write>(&self, writer: &mut W) -> Result<(), std::fmt::Error> {
        writeln!(writer, "digraph {}{{", self.global_data.name)?;

        let mut writer = CodeFormatter::new(writer, "  ");
        writer.indent(1);

        write!(writer, "{}", self.global_data)?;
        for (n, (_, _, v)) in self.iter_nodes().enumerate() {
            let mut node_data: DotVertexData = v.clone();
            node_data.remove_common(&self.global_data);
            let data = node_data.to_string();
            write!(
                writer,
                "\n{}{};",
                node_data.name.clone().unwrap_or(n.to_string()),
                if !data.is_empty() {
                    format!("[{data}]")
                } else {
                    String::new()
                }
            )?;
        }
        writer.write_str("\n")?;

        for (hedge_pair, eid, data) in self.iter_edges() {
            let mut edata = data.data.clone();
            edata.remove_common(&self.global_data);
            let attr = GVEdgeAttrs {
                color: None,
                label: None,
                other: Some(edata.to_string()),
            };

            hedge_pair.add_data(self).dot_fmt(
                &mut writer,
                self,
                eid,
                |h| h.statement.clone(),
                |a| self[a].name.clone().unwrap_or(a.to_string()),
                data.orientation,
                attr,
            )?;
        }
        writeln!(writer, "}}")?;
        Ok(())
    }

    #[cfg(feature = "rkyv")]
    pub fn to_rkyv_bytes<const BYTES: usize>(&self) -> Result<rkyv::AlignedVec, String>
    where
        Self: rkyv::Serialize<rkyv::ser::serializers::AllocSerializer<BYTES>>,
    {
        rkyv::to_bytes::<_, BYTES>(self).map_err(|err| err.to_string())
    }

    #[cfg(feature = "rkyv")]
    /// Returns the archived graph root without validating the byte buffer.
    ///
    /// # Safety
    ///
    /// `bytes` must contain a valid rkyv archive produced for this exact
    /// `DotGraph` type, and the returned reference must not outlive `bytes`.
    pub unsafe fn archived_from_bytes(bytes: &[u8]) -> &<Self as rkyv::Archive>::Archived
    where
        Self: rkyv::Archive,
    {
        unsafe { rkyv::archived_root::<Self>(bytes) }
    }

    #[allow(clippy::result_large_err, clippy::type_complexity)]
    pub fn from_file<'a, P>(p: P) -> Result<Self, HedgeParseError<'a, (), (), (), ()>>
    where
        P: AsRef<Path>,
    {
        let ast_graph: SubGraphFreeGraph = dot_parser::ast::Graph::from_file(p)?.into();

        Ok(Self::from((ast_graph, Figment::new())))
    }

    #[allow(clippy::result_large_err, clippy::type_complexity)]
    pub fn from_string<'a, Str: AsRef<str>>(
        s: Str,
    ) -> Result<Self, HedgeParseError<'a, (), (), (), ()>> {
        let ast_graph: SubGraphFreeGraph = dot_parser::ast::Graph::try_from(s.as_ref())?
            .filter_map(&|(k, v)| Some((k.into(), strip_quotes(v).to_string())))
            .into();

        Ok(Self::from((ast_graph, Figment::new())))
    }

    #[cfg(feature = "serde")]
    pub fn from_string_with_figment<'a, Str: AsRef<str>>(
        s: Str,
        figment: figment::Figment,
    ) -> Result<Self, HedgeParseError<'a, (), (), (), ()>> {
        let ast_graph: SubGraphFreeGraph = dot_parser::ast::Graph::try_from(s.as_ref())?
            .filter_map(&|(k, v)| Some((k.into(), strip_quotes(v).to_string())))
            .into();

        let graph = Self::from((ast_graph, figment));
        Ok(graph)
    }

    pub fn back_and_forth_dot(self) -> Self {
        Self::from_string(self.debug_dot()).unwrap_or_else(|_| {
            panic!(
                "Failed to parse back the DOT serialization of the graph: {}",
                self.debug_dot()
            )
        })
    }

    pub fn debug_dot(&self) -> String {
        let mut out = String::new();
        self.write_fmt(&mut out).unwrap();
        // println!("{out}");
        out
    }

    pub fn dot_of<Sub: SubGraphLike>(&self, subgraph: &Sub) -> String {
        let mut output = String::new();
        self.dot_impl_fmt(
            &mut output,
            subgraph,
            self.global_data.to_string(),
            &|s| s.statement.clone(),
            &|s| Some(s.to_string()).filter(|attributes| !attributes.is_empty()),
            &|s| Some(s.to_string()).filter(|attributes| !attributes.is_empty()),
        )
        .unwrap();
        output
    }

    pub fn format_dot(
        self,
        edge_format: impl AsRef<str>,
        vertex_format: impl AsRef<str>,
    ) -> String {
        self.graph.dot_impl(
            &self.graph.full_filter(),
            "",
            &|_| None,
            &|d| Some(format!("{d}label={}", d.format(&edge_format))),
            &|d| Some(format!("{d}label={}", d.format(&vertex_format))),
        )
    }
}

pub mod error;
pub use error::{HedgeParseError, HedgeParseExt};

impl<S: NodeStorageOps<NodeData = DotVertexData>> From<(SubGraphFreeGraph, Figment)>
    for DotGraph<S>
{
    fn from((value, fig): (SubGraphFreeGraph, Figment)) -> Self {
        let is_digraph = value.is_digraph;
        let name = value.name.clone();
        let (attrs, ids, nodes, edges) = value.nodes_and_edges();

        // let can_graph = dot_parser::canonical::Graph::from(ast_graph);
        let mut global_data = GlobalData::try_from((attrs, ids)).unwrap();
        if let Some(name) = name {
            global_data.add_name(name);
        }
        global_data.set_figment(fig);

        let mut g = HedgeGraphBuilder::new();
        let mut map = BTreeMap::new();

        for (id, n) in nodes {
            let idorstatements = match DotVertexData::from_parser(n, &global_data) {
                Either::Left(d) => NodeIdOrDangling::Id(g.add_node(d)),
                Either::Right(statements) => NodeIdOrDangling::Dangling { statements },
            };

            map.insert(id, idorstatements);
        }

        for e in edges
            .set
            .into_iter()
            .sorted_by(|a, b| Ord::cmp(&(&a.from, &a.to), &(&b.from, &b.to)))
        {
            let (data, orientation, source, target) =
                DotEdgeData::from_parser(e, &map, is_digraph, &global_data);
            match target {
                Either::Left(a) => {
                    g.add_edge(source, a, data, orientation);
                }
                Either::Right(flow) => {
                    g.add_external_edge(source, data, orientation, flow);
                }
            }
        }

        let mut g: DotGraph<S> = DotGraph {
            global_data,
            graph: g.build(),
        };

        // println!("Built: {}", g.debug_dot());
        g.apply_explicit_id_ordering();
        g
    }
}

impl<S: NodeStorageOps<NodeData = DotVertexData>> DotGraph<S> {
    pub fn apply_explicit_id_ordering(&mut self) {
        // println!("Built: {}", self.debug_dot());

        let mut used_edges = HashSet::new();
        let n_edges = self.n_edges();

        let mut edge_map = self.new_edgevec(|d, e, _| {
            d.edge_id.inspect(|d| {
                assert!(
                    used_edges.insert(*d),
                    "Duplicate edge ID: {d} for edge {e}:{used_edges:?}"
                );
                assert!(d.0 < n_edges, "Edge {d} out of bounds (len={n_edges})")
            })
        });

        let mut used_hedges = HashSet::new();

        let n_hedges = self.n_hedges();
        let mut hedge_map = self.new_hedgevec(|h, d| {
            d.id.inspect(|d| {
                assert!(
                    used_hedges.insert(*d),
                    "Duplicate hedge ID: {d} for hedge {h}",
                );
                assert!(d.0 < n_hedges, "Hedge {d} out of bounds (len={n_hedges})")
            })
        });

        let mut used_nodes = HashSet::new();
        let n_nodes = self.n_nodes();
        let mut node_map = self.new_nodevec(|ni, _, v| {
            v.index.inspect(|i| {
                assert!(
                    used_nodes.insert(*i),
                    "Duplicate node index: {i} for node {ni}"
                );
                assert!(i.0 < n_nodes, "Node {i} out of bounds (len ={n_nodes})")
            })
        });

        // println!(
        //     "Hedge Map: {}",
        //     hedge_map.display_string(|i| format!("{i:?}"))
        // );
        // println!(
        //     "Node Map: {}",
        //     node_map.display_string(|i| format!("{i:?}"))
        // );
        // println!(
        //     "Edge Map: {}",
        //     edge_map.display_string(|i| format!("{i:?}"))
        // );
        node_map.fill_in(|id| used_nodes.contains(id));
        edge_map.fill_in(|id| used_edges.contains(id));
        hedge_map.fill_in(|id| used_hedges.contains(id));
        // println!(
        //     "Filled Hedge Map: {}",
        //     hedge_map.display_string(|i| format!("{i:?}"))
        // );
        // println!(
        //     "Filled Node Map: {}",
        //     node_map.display_string(|i| format!("{i:?}"))
        // );
        // println!(
        //     "Filled Edge Map: {}",
        //     edge_map.display_string(|i| format!("{i:?}"))
        // );
        let edge_perm: Permutation = edge_map.try_into().unwrap();
        let node_perm: Permutation = node_map.try_into().unwrap();
        let hedge_perm: Permutation = hedge_map.try_into().unwrap();
        // println!("Hedge Perm: {hedge_perm}");
        // println!("Edge Perm: {edge_perm}");
        // println!("Node Perm: {node_perm}");

        <HedgeGraph<_, _, _, _> as Swap<Hedge>>::permute(&mut self.graph, &hedge_perm);
        // println!("Permuted Hedge Graph: {}", g.debug_dot());
        <HedgeGraph<_, _, _, _> as Swap<EdgeIndex>>::permute(&mut self.graph, &edge_perm);
        // println!("Permuted Edge Graph: {}", g.debug_dot());
        <HedgeGraph<_, _, _, _> as Swap<NodeIndex>>::permute(&mut self.graph, &node_perm);
        // println!("Permuted Node Graph:{}", g.debug_dot());
    }
}

#[macro_export]
macro_rules! dot {
    ($($t:tt)*) => {
        $crate::parser::DotGraph::from_string(stringify!($($t)*))
    };
}

impl<E, V, H, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    pub fn dot_serialize_of<S: SubGraphLike>(
        &self,
        subgraph: &S,
        global: impl Into<GlobalData>,
        hedge_map: &impl Fn(&H) -> DotHedgeData,
        edge_map: &impl Fn(&E) -> DotEdgeData,
        node_map: &impl Fn(&V) -> DotVertexData,
    ) -> String {
        let global_data = global.into();

        let g = self.map_data_ref(
            |_, _, v| node_map(v),
            |_, _, _, d| d.map(edge_map),
            |_, h| hedge_map(h),
        );

        let dot_graph = DotGraph {
            graph: g,
            global_data,
        };

        dot_graph.dot_of(subgraph)
    }

    pub fn dot_serialize_io(
        &self,
        writer: &mut impl std::io::Write,
        global: impl Into<GlobalData>,
        hedge_map: &impl Fn(&H) -> DotHedgeData,
        edge_map: &impl Fn(&E) -> DotEdgeData,
        node_map: &impl Fn(&V) -> DotVertexData,
    ) -> Result<(), std::io::Error> {
        let global_data = global.into();

        let g = self.map_data_ref(
            |_, _, v| node_map(v),
            |_, _, _, d| d.map(edge_map),
            |_, h| hedge_map(h),
        );

        let dot_graph = DotGraph {
            graph: g,
            global_data,
        };

        dot_graph.write_io(writer)
    }

    pub fn dot_serialize_fmt(
        &self,
        writer: &mut impl std::fmt::Write,
        global: impl Into<GlobalData>,
        hedge_map: &impl Fn(&H) -> DotHedgeData,
        edge_map: &impl Fn(&E) -> DotEdgeData,
        node_map: &impl Fn(&V) -> DotVertexData,
    ) -> Result<(), std::fmt::Error> {
        let global_data = global.into();

        let g = self.map_data_ref(
            |_, _, v| node_map(v),
            |_, _, _, d| d.map(edge_map),
            |_, h| hedge_map(h),
        );

        let dot_graph = DotGraph {
            graph: g,
            global_data,
        };

        dot_graph.write_fmt(writer)
    }
}

#[cfg(test)]
pub mod test {

    use crate::{
        half_edge::{nodestore::NodeStorageVec, subgraph::SuBitGraph},
        parser::{DotGraph, DotVertexData},
    };

    use super::GraphSet;

    #[test]
    fn orientations() {
        let s: DotGraph = dot!(
        digraph bba1{
            num=1
                        ext    [style=invis]
                        ext -> A:1   [particle=a id=1]
                        ext -> A:2    [particle="b" id=2]
                        A:0  -> ext  [particle="b" id=0]
                    })
        .unwrap();

        let g = s.back_and_forth_dot();
        let gg = g.clone().back_and_forth_dot();
        println!("{}", g.debug_dot());
        assert_eq!(g, gg);
    }

    #[test]
    fn test_from_string() {
        let s = "digraph G {
            A      [style=invis]
            A -> B [label=\"Hello\" sink=\"AAA\"];
            B -> C [label=\"World\" dir=back];
        }";
        let graph: DotGraph = DotGraph::from_string(s).unwrap();

        println!("Parsed as:{}", graph.debug_dot());
        // println!("Inv:{}", graph.graph.as_ref());
        assert_eq!(graph.n_nodes(), 2);
        assert_eq!(graph.n_internals(), 1);
        let g = graph.back_and_forth_dot();
        let gg = g.clone().back_and_forth_dot();
        // println!("{g:?}");
        assert_eq!(g, gg);
    }

    #[test]
    fn dot_of_omits_empty_attribute_lists() {
        let graph: DotGraph = DotGraph::from_string("digraph { A -> B; }").unwrap();
        let serialized = graph.dot_of(&graph.full_filter());

        assert!(!serialized.contains(" []"));
        let reparsed: Result<DotGraph, _> = DotGraph::from_string(serialized);
        assert!(reparsed.is_ok());
    }

    #[test]
    fn multiple_graphs() {
        let s = r#"
            digraph triangle_0 {
            graph [
            overall_factor = 1;
            multiplicity_factor = 1;
            ]
            edge [
            pdg=1000
            dod=-100
            ]
            ext [style=invis]

            }

            digraph tria {
            graph [
            overall_factor = "-1";
            multiplicity_factor = 1;
            ]
            edge [
            pdg=1000
            dod=-100
            ]
            ext [style=invis]
            ext -> v4 ;
            ext -> v5 [name=p2 mom=p2 num="Q(eid,spenso::mink(4,20))" is_dummy=true];
            }
            "#;

        let set = GraphSet::<_, _, _, _, NodeStorageVec<DotVertexData>>::from_string(s).unwrap();
        assert_eq!(set.set.len(), 2);

        for g in set {
            let gg = g.back_and_forth_dot();
            let ggg = gg.clone().back_and_forth_dot();
            assert_eq!(gg, ggg);
        }
    }

    #[test]
    fn test_macro() {
        let _: DotGraph = dot!( digraph {
           node [shape=circle,height=0.1,label=""];  overlap="scale"; layout="neato";
         0 -> 7[ dir=none color="red:blue;0.5",label="a"];
        0 -> 12[ dir=forward color="red:blue;0.5",label="d"];
        1 -> 0[ dir=forward color="red:blue;0.5",label="d"];
        1 -> 3[ dir=none color="red:blue;0.5",label="a"];
        2 -> 1[ dir=forward color="red:blue;0.5",label="d"];
        2 -> 6[ dir=none color="red:blue;0.5",label="a"];
        3 -> 13[ dir=forward color="red:blue;0.5",label="d"];
        4 -> 3[ dir=forward color="red:blue;0.5",label="d"];
        4 -> 5[ dir=none color="red:blue;0.5",label="g"];
        5 -> 2[ dir=forward color="red:blue;0.5",label="d"];
        6 -> 7[ dir=forward color="red:blue;0.5",label="e-"];
        7 -> 11[ dir=forward color="red:blue;0.5",label="e-"];
        8 -> 6[ dir=forward color="red:blue;0.5",label="e-"];
        9 -> 4[ dir=forward color="red:blue;0.5",label="d"];
        10 -> 5[ dir=forward color="red:blue;0.5",label="d"];
        })
        .unwrap();
    }

    #[test]
    fn underlying_alignment() {
        let s = "digraph {
          0 [name=B];
          1 [name=C];
          ext0  [style=invis];
          ext0 -> 0 [dir=forward];
          ext1  [style=invis];
          ext1 -> 0 [dir=forward];
          ext2  [style=invis];
          ext2 -> 1 [dir=back];
          1 -> 0    [dir=back];
          ext5 [style=invis];
          ext5 -> 1 [dir=forward];
        }";
        let graph: DotGraph = DotGraph::from_string(s).unwrap();

        let serialized = graph.debug_dot();

        let colored = graph.dot_of(&graph.full_filter());

        // println!(
        //     "{}",
        //     graph.dot_impl(&graph.full_filter(), "", &|a| None, &|b| Some(format!(
        //         "label={}",
        //         b.id
        //     )))
        // );

        let mut graph2: DotGraph = DotGraph::from_string(serialized.clone()).unwrap();

        let serialized2 = graph.debug_dot();

        let colored2 = graph2.dot_of(&graph2.full_filter());

        assert_eq!(
            serialized, serialized2,
            "{serialized}//not equal to \n{serialized2}",
        );
        assert_eq!(colored, colored2, "{colored}\nnot equal to\n{colored2}");

        println!(
            "{}",
            graph2.dot_impl(&graph.full_filter(), "", &|_| None, &|_| None, &|b| Some(
                format!("label={:?}", b.name)
            ))
        );
        // println!("{}",graph.ed)
        graph2.align_underlying_to_superficial();
        println!(
            "{}",
            graph2.dot_impl(&graph.full_filter(), "", &|_| None, &|_| None, &|b| Some(
                format!("label={:?}", b.name)
            ))
        );

        let serialized2 = graph2.debug_dot();

        println!("{serialized2}");

        let aligned: DotGraph = dot!(
        digraph {
          1 [name=C];
          0:1	-> 1:0	 [id=0 ];
          ext1	 [style=invis];
          1:4	-> ext3	 [id=3 ];
          ext3	 [style=invis];
          ext2	 [style=invis];
          ext4	 [style=invis];
          0 [name=B];
          ext4	-> 1:5	 [id=4 ];
          ext1	-> 0:2	 [id=1 ];
          ext2	-> 0:3	 [id=2 ];
        })
        .unwrap();

        graph2 = graph2.back_and_forth_dot();

        assert_eq!(
            aligned,
            graph2,
            "{}\n//not equal to\n{}",
            aligned.dot_display(&aligned.full_filter()),
            graph2.dot_display(&graph2.full_filter())
        );
        // assert_eq!(graph.n_nodes(), 2);
        // assert_eq!(graph.n_internals(), 1);
    }

    #[test]
    fn subgraph() {
        let aligned: DotGraph = dot!(
        digraph {
          1 [name=C];
          0:1:s	-> 1:0	 [id=0 ];
          ext1	 [style=invis];
          1:4	-> ext3	 [id=3 ];
          ext3	 [style=invis];
          ext2	 [style=invis];
          ext4	 [style=invis];
          0 [name=B];
          ext4	-> 1:5	 [id=4 ];
          ext1	-> 0:2:s	 [id=1 ];
          ext2	-> 0:3:s	 [id=2 ];
        })
        .unwrap();

        let sub: SuBitGraph = aligned.compass_subgraph(Some(dot_parser::ast::CompassPt::S));

        println!("{}", aligned.dot_of(&sub));
    }

    #[test]
    fn test_quote_stripping() {
        let s = r#"digraph G {
            graph [title="Graph Title"];
            A [label="Quoted Label"];
            A -> B [label="Edge Label"];
        }"#;

        let graph: DotGraph = DotGraph::from_string(s).unwrap();

        // Check that quotes are stripped from global graph attributes
        assert_eq!(
            graph.global_data.statements.get("title").unwrap(),
            "Graph Title"
        );

        // Check that quotes are stripped from node attributes
        let node_a = graph
            .iter_nodes()
            .find(|(_, _, data)| data.statements.contains_key("label"))
            .unwrap()
            .2;
        assert_eq!(node_a.statements.get("label").unwrap(), "Quoted Label");

        // Check that quotes are stripped from edge attributes
        let edge = graph.iter_edges().next().unwrap().2.data;
        assert_eq!(edge.statements.get("label").unwrap(), "Edge Label");
    }
}

mod multi {

    #[test]
    fn multiple() {
        assert_eq!(
            dot_parser::ast::Graphs::try_from(
                r#"digraph d{
                overall_factor = 1
                edge [overall_factor=1]
                A
                }

                digraph P{
                overall_factor = 1
                edge [overall_factor=1]
                A
                }
                "#
            )
            .unwrap()
            .graphs
            .len(),
            2
        );
    }
}
