//! # DOT interchange
//!
//! This module parses the [DOT language](https://graphviz.org/doc/info/lang.html)
//! into Linnet's half-edge model and serializes half-edge graphs back to DOT.
//! [`DotGraph::from_string`] and [`DotGraph::from_file`] parse input holding
//! exactly one graph. A [`DotGraph`] retains [`GlobalData`] alongside a
//! [`HedgeGraph`] whose node, edge, and half-edge payloads are [`DotVertexData`],
//! [`DotEdgeData`], and [`DotHedgeData`]. Use [`set::DotGraphSet`] for input
//! containing multiple graphs.
//!
//! Parsed DOT attributes remain available in those payload types. In particular,
//! external-port flow and edge orientation are preserved in the parsed graph
//! rather than discarded as drawing-only attributes. [`DotGraph`] dereferences to
//! its inner [`HedgeGraph`] for graph algorithms. Invalid DOT is reported as a
//! [`DotSyntaxError`] with its line and column; well-formed DOT that does not
//! describe a valid half-edge graph is reported as a [`DotParseError`].
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
//! assert_eq!(graph.n_live(), 2);
//! assert_eq!(graph.n_edges(), 1);
//! assert_eq!(graph.n_hedges(), 2);
//!
//! let mut output = String::new();
//! graph.write_fmt(&mut output).unwrap();
//! let round_trip: DotGraph = DotGraph::from_string(output).unwrap();
//! assert_eq!(round_trip.n_live(), graph.n_live());
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
    borrow::Cow,
    collections::BTreeMap,
    fmt::{self, Debug, Write},
    ops::{Deref, DerefMut},
    path::Path,
};

#[cfg(feature = "serde")]
use figment;

use ahash::{HashSet, HashSetExt};
use figment::Figment;
use itertools::Either;

use crate::{
    half_edge::{
        builder::HedgeGraphBuilder,
        involution::{EdgeIndex, Hedge},
        nodestore::{NodeForest, Topology},
        remap::{HedgeRemap, NodeRemap},
        subgraph::{ModifySubSet, SubGraphLike, SubSetLike},
        swap::Swap,
        GVEdgeAttrs, HedgeGraph, NodeIndex,
    },
    permutation::Permutation,
};

mod syntax;
use syntax::FlatGraph;
pub use syntax::{CompassPt, DotSyntaxError, UnknownCompassPt};

/// Writes an identifier bare when DOT allows it, and quoted otherwise.
pub(crate) fn dot_id(value: &str) -> String {
    let keyword = matches!(
        value.to_ascii_lowercase().as_str(),
        "node" | "edge" | "graph" | "digraph" | "subgraph" | "strict"
    );
    let mut characters = value.chars();
    let plain = characters
        .next()
        .is_some_and(|first| first.is_ascii_alphabetic() || first == '_')
        && characters.all(|character| character.is_ascii_alphanumeric() || character == '_');
    if plain && !keyword {
        value.to_owned()
    } else {
        dot_value(value)
    }
}

/// Writes `value` as a quoted DOT string that this module reads back as
/// `value`. Quotes are escaped, and a backslash is doubled exactly where it
/// would otherwise escape the next character: before a quote, a backslash, a
/// newline, or the closing quote. Other backslashes stay single, so Graphviz
/// escapes such as `\N` and `\l` reach Graphviz unchanged.
pub(crate) fn dot_value(value: &str) -> String {
    let mut quoted = String::with_capacity(value.len() + 2);
    quoted.push('"');
    let mut chars = value.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '"' => quoted.push_str("\\\""),
            '\\' if matches!(chars.peek(), None | Some('"' | '\\' | '\n')) => {
                quoted.push_str("\\\\")
            }
            c => quoted.push(c),
        }
    }
    quoted.push('"');
    quoted
}

/// Indents the statements of the `*_dot_fmt` writers, which emit each
/// statement as a chunk that starts with a newline.
///
/// Such a chunk is written as indented lines, each ending in a newline, after
/// removing the leading whitespace its lines share; any other chunk is written
/// unchanged. Line breaks inside quoted strings belong to the value: they are
/// kept verbatim and do not start a new line.
pub(crate) struct DotIndent<'w, W> {
    writer: &'w mut W,
}

impl<'w, W: Write> DotIndent<'w, W> {
    pub(crate) fn new(writer: &'w mut W) -> Self {
        Self { writer }
    }

    /// Splits `chunk` at the line breaks outside quoted strings.
    fn lines(chunk: &str) -> impl Iterator<Item = &str> + Clone {
        let (mut quoted, mut escaped) = (false, false);
        chunk.split(move |c| {
            let line_break = c == '\n' && !quoted;
            if escaped {
                escaped = false;
            } else if c == '\\' {
                escaped = quoted;
            } else if c == '"' {
                quoted = !quoted;
            }
            line_break
        })
    }
}

impl<W: Write> Write for DotIndent<'_, W> {
    fn write_str(&mut self, chunk: &str) -> fmt::Result {
        let Some(chunk) = chunk.strip_prefix('\n') else {
            return self.writer.write_str(chunk);
        };
        let chunk = chunk.trim_end_matches(|c: char| c.is_ascii_whitespace());
        let indent = |line: &str| {
            line.len()
                - line
                    .trim_start_matches(|c: char| c.is_ascii_whitespace())
                    .len()
        };
        let lines = Self::lines(chunk);
        let common = (lines.clone().map(indent).filter(|&n| n > 0).min()).unwrap_or(0);
        for line in lines {
            let line = &line[common.min(indent(line))..];
            if !line.is_empty() {
                self.writer.write_str("  ")?;
            }
            self.writer.write_str(line)?;
            self.writer.write_char('\n')?;
        }
        Ok(())
    }

    /// Writes a whole `write!` as one chunk.
    fn write_fmt(&mut self, args: fmt::Arguments<'_>) -> fmt::Result {
        match args.as_str() {
            Some(chunk) => self.write_str(chunk),
            None => self.write_str(&args.to_string()),
        }
    }
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
pub struct DotGraph<N: Topology = NodeForest> {
    /// Graph-level DOT metadata, including its name, payload, and attribute defaults.
    pub global_data: GlobalData,
    /// Parsed half-edge topology together with its DOT edge, vertex, and half-edge data.
    pub graph: HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData, N>,
}

impl<N: Topology> Deref for DotGraph<N> {
    type Target = HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData, N>;

    fn deref(&self) -> &Self::Target {
        &self.graph
    }
}

impl<N: Topology> DerefMut for DotGraph<N> {
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

impl<S: Topology> DotGraph<S> {
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
        let name = if self.global_data.name.is_empty() {
            String::new()
        } else {
            format!(" {}", dot_id(&self.global_data.name))
        };
        writeln!(writer, "digraph{name} {{")?;

        writeln!(writer, "{:4}", self.global_data)?;

        for (n, _, v) in self.iter_nodes() {
            let mut node_data: DotVertexData = v.clone();
            node_data.remove_common(&self.global_data);

            if let Some(name) = &node_data.name {
                write!(writer, "\t{}", dot_id(name))?;
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
                |a| {
                    self[a]
                        .name
                        .as_deref()
                        .map(dot_id)
                        .unwrap_or_else(|| a.to_string())
                },
                data.orientation,
                attr,
            )?;
        }
        writeln!(writer, "}}")?;
        Ok(())
    }

    pub fn write_fmt<W: std::fmt::Write>(&self, writer: &mut W) -> Result<(), std::fmt::Error> {
        let name = if self.global_data.name.is_empty() {
            String::new()
        } else {
            format!(" {}", dot_id(&self.global_data.name))
        };
        writeln!(writer, "digraph{name} {{")?;

        let mut writer = DotIndent::new(writer);

        write!(writer, "{}", self.global_data)?;
        for (n, _, v) in self.iter_nodes() {
            let mut node_data: DotVertexData = v.clone();
            node_data.remove_common(&self.global_data);
            let data = node_data.to_string();
            write!(
                writer,
                "\n{}{};",
                node_data
                    .name
                    .as_deref()
                    .map(dot_id)
                    .unwrap_or_else(|| n.to_string()),
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
                |a| {
                    self[a]
                        .name
                        .as_deref()
                        .map(dot_id)
                        .unwrap_or_else(|| a.to_string())
                },
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

    /// Parses a file holding exactly one DOT graph.
    pub fn from_file<P>(p: P) -> Result<Self, HedgeParseError<(), (), (), ()>>
    where
        P: AsRef<Path>,
    {
        let source = std::fs::read_to_string(p)?;
        Self::from_flat(FlatGraph::parse(&source)?, Figment::new())
    }

    /// Parses input holding exactly one DOT graph; trailing input is an error.
    pub fn from_string<Str: AsRef<str>>(s: Str) -> Result<Self, HedgeParseError<(), (), (), ()>> {
        Self::from_flat(FlatGraph::parse(s.as_ref())?, Figment::new())
    }

    #[cfg(feature = "serde")]
    pub fn from_string_with_figment<Str: AsRef<str>>(
        s: Str,
        figment: figment::Figment,
    ) -> Result<Self, HedgeParseError<(), (), (), ()>> {
        Self::from_flat(FlatGraph::parse(s.as_ref())?, figment)
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

    /// Builds the half-edge graph of a flattened DOT graph. Nodes are added in
    /// identifier order and edges in endpoint order before the explicit `id`
    /// and numeric-port reservations are applied.
    fn from_flat(
        graph: FlatGraph<'_>,
        fig: Figment,
    ) -> Result<Self, HedgeParseError<(), (), (), ()>> {
        let FlatGraph {
            directed,
            name,
            graph_attrs,
            node_defaults,
            edge_defaults,
            nodes,
            mut edges,
        } = graph;
        let mut global_data =
            GlobalData::from_statements(graph_attrs, node_defaults, edge_defaults);
        if let Some(name) = name {
            global_data.add_name(name.into_owned());
        }
        global_data.set_figment(fig);

        let mut g = HedgeGraphBuilder::new();
        let map = (nodes.into_iter())
            .map(|(id, attrs)| {
                let node = match DotVertexData::from_parser(&id, attrs, &global_data)? {
                    Either::Left(d) => NodeIdOrDangling::Id(g.add_node(d)),
                    Either::Right(statements) => NodeIdOrDangling::Dangling { statements },
                };
                Ok((id, node))
            })
            .collect::<Result<BTreeMap<Cow<str>, _>, DotParseError>>()?;

        // Stable, so parallel edges with equal endpoints keep their flattened order.
        edges.sort_by(|a, b| (&a.from, &a.to).cmp(&(&b.from, &b.to)));
        for e in edges {
            let (data, orientation, source, target) =
                DotEdgeData::from_parser(e, &map, directed, &global_data)?;
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
        g.apply_explicit_id_ordering()?;
        Ok(g)
    }
}

pub mod error;
pub use error::{DotParseError, ExplicitIdKind, HedgeParseError, HedgeParseExt};

impl<S: Topology> DotGraph<S> {
    pub fn apply_explicit_id_ordering(&mut self) -> Result<(), DotParseError> {
        // println!("Built: {}", self.debug_dot());

        let mut used_edges = HashSet::new();
        let n_edges = self.n_edges();

        let mut edge_map = self.new_edgevec(|d, _, _| d.edge_id);
        for id in edge_map.iter().filter_map(|(_, id)| *id) {
            if id.0 >= n_edges {
                return Err(DotParseError::ExplicitIdOutOfBounds {
                    kind: ExplicitIdKind::Edge,
                    id: id.0,
                    len: n_edges,
                });
            }
            if !used_edges.insert(id) {
                return Err(DotParseError::DuplicateExplicitId {
                    kind: ExplicitIdKind::Edge,
                    id: id.0,
                });
            }
        }

        let mut used_hedges = HashSet::new();

        let n_hedges = self.n_hedges();
        let mut hedge_map = self.new_hedgevec(|_, d| d.id);
        for id in hedge_map.iter().filter_map(|(_, id)| *id) {
            if id.0 >= n_hedges {
                return Err(DotParseError::ExplicitIdOutOfBounds {
                    kind: ExplicitIdKind::Hedge,
                    id: id.0,
                    len: n_hedges,
                });
            }
            if !used_hedges.insert(id) {
                return Err(DotParseError::DuplicateExplicitId {
                    kind: ExplicitIdKind::Hedge,
                    id: id.0,
                });
            }
        }

        let mut used_nodes = HashSet::new();
        let node_bound = self.node_bound();
        let mut node_map = self.new_nodevec(|_, _, v| v.index);
        for id in node_map.iter().filter_map(|(_, id)| *id) {
            if id.0 >= node_bound {
                return Err(DotParseError::ExplicitIdOutOfBounds {
                    kind: ExplicitIdKind::Node,
                    id: id.0,
                    len: node_bound,
                });
            }
            if !used_nodes.insert(id) {
                return Err(DotParseError::DuplicateExplicitId {
                    kind: ExplicitIdKind::Node,
                    id: id.0,
                });
            }
        }

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
        let edge_perm: Permutation =
            edge_map
                .try_into()
                .map_err(|()| DotParseError::InvalidExplicitIdPermutation {
                    kind: ExplicitIdKind::Edge,
                })?;
        let node_perm: Permutation =
            node_map
                .try_into()
                .map_err(|()| DotParseError::InvalidExplicitIdPermutation {
                    kind: ExplicitIdKind::Node,
                })?;
        let hedge_perm: Permutation =
            hedge_map
                .try_into()
                .map_err(|()| DotParseError::InvalidExplicitIdPermutation {
                    kind: ExplicitIdKind::Hedge,
                })?;
        // println!("Hedge Perm: {hedge_perm}");
        // println!("Edge Perm: {edge_perm}");
        // println!("Node Perm: {node_perm}");

        const SETTLED: &str = "parsed graphs have no pending identifications";
        (self.graph)
            .permute_hedges(&HedgeRemap::from_permutation(&hedge_perm))
            .expect(SETTLED);
        <HedgeGraph<_, _, _, _> as Swap<EdgeIndex>>::permute(&mut self.graph, &edge_perm);
        (self.graph)
            .permute_nodes(&NodeRemap::from_permutation(&node_perm))
            .expect(SETTLED);
        Ok(())
    }
}

#[macro_export]
macro_rules! dot {
    ($($t:tt)*) => {
        $crate::parser::DotGraph::from_string(stringify!($($t)*))
    };
}

impl<E, V, H, N: Topology> HedgeGraph<E, V, H, N> {
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
        half_edge::{nodestore::NodeForest, subgraph::SuBitGraph, NodeIndex},
        parser::{CompassPt, DotGraph, DotParseError, ExplicitIdKind, HedgeParseError},
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
        assert_eq!(graph.n_live(), 2);
        assert_eq!(graph.n_internals(), 1);
        let g = graph.back_and_forth_dot();
        let gg = g.clone().back_and_forth_dot();
        // println!("{g:?}");
        assert_eq!(g, gg);
    }

    /// Nodes absorbed by a pending identification are not written: both DOT
    /// writers and the Symbolica conversion show the identified graph, as
    /// compacting first would.
    #[test]
    fn pending_identifications_are_written_as_the_identified_graph() {
        let mut graph: DotGraph = DotGraph::from_string("digraph { a -> b; b -> c; }").unwrap();
        let named = |graph: &DotGraph, name: &str| {
            (graph.iter_nodes())
                .find(|(_, _, node)| node.name.as_deref() == Some(name))
                .map(|(index, _, _)| index)
                .unwrap()
        };
        let (a, c) = (named(&graph, "a"), named(&graph, "c"));
        let merged = graph[a].clone();
        graph.identify_nodes(&[a, c], merged);
        let mut compacted = graph.clone();
        compacted.compact();

        let mut written = String::new();
        graph.write_fmt(&mut written).unwrap();
        let back: DotGraph = DotGraph::from_string(written).unwrap();
        assert_eq!(
            (back.n_live(), back.node_bound(), back.n_edges()),
            (2, 2, 2)
        );
        assert_eq!(
            graph.clone().back_and_forth_dot(),
            compacted.back_and_forth_dot()
        );
        #[cfg(feature = "symbolica")]
        {
            let identified = symbolica::graph::Graph::try_from(graph.graph).unwrap();
            assert_eq!((identified.nodes().len(), identified.edges().len()), (2, 2));
        }
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
    fn canonical_dot_escapes_names_keys_and_quoted_values() {
        let source = r#"digraph "graph \"name" {
            "node \"one" ["node-key"="a\"b\N"];
            "node \"one" -> other ["edge-key"="c\"d\l"];
        }"#;
        let mut graph: DotGraph = DotGraph::from_string(source).unwrap();
        let node_index = graph
            .iter_nodes()
            .find(|(_, _, node)| node.name.as_deref() == Some("node \"one"))
            .map(|(index, _, _)| index)
            .unwrap();
        let edge_index = graph.iter_edges().next().unwrap().1;
        graph
            .global_data
            .statements
            .insert("global-terminal".to_owned(), "graph\\".to_owned());
        graph[node_index]
            .statements
            .insert("node-terminal".to_owned(), "node\\".to_owned());
        graph[edge_index]
            .statements
            .insert("edge-terminal".to_owned(), "edge\\".to_owned());
        let serialized = graph.debug_dot();
        let reparsed: DotGraph = DotGraph::from_string(&serialized).unwrap();

        assert_eq!(reparsed.global_data.name, "graph \"name");
        assert_eq!(
            reparsed
                .global_data
                .statements
                .get("global-terminal")
                .map(String::as_str),
            Some("graph\\")
        );
        assert!(reparsed.iter_nodes().any(|(_, _, node)| {
            node.name.as_deref() == Some("node \"one")
                && node.statements.get("node-key").map(String::as_str) == Some("a\"b\\N")
                && node.statements.get("node-terminal").map(String::as_str) == Some("node\\")
        }));
        assert_eq!(
            reparsed
                .iter_edges()
                .next()
                .and_then(|(_, _, edge)| edge.data.statements.get("edge-key"))
                .map(String::as_str),
            Some("c\"d\\l")
        );
        assert_eq!(
            reparsed
                .iter_edges()
                .next()
                .and_then(|(_, _, edge)| edge.data.statements.get("edge-terminal"))
                .map(String::as_str),
            Some("edge\\")
        );
        assert_eq!(reparsed, reparsed.clone().back_and_forth_dot());
        assert!(serialized.contains("\"graph \\\"name\""));
        assert!(serialized.contains("\"node-key\"=\"a\\\"b\\N\""));
        assert!(serialized.contains("\"edge-key\"=\"c\\\"d\\l\""));
        assert!(serialized.contains(r#""global-terminal" = "graph\\""#));
        assert!(serialized.contains(r#""node-terminal"="node\\""#));
        assert!(serialized.contains(r#""edge-terminal"="edge\\""#));
    }

    /// Every string survives the writer and the parser, including backslashes
    /// next to quotes, newlines or the closing quote, and multi-line values,
    /// whose continuation lines the writer must not indent.
    #[test]
    fn quoted_values_round_trip_through_the_writer() {
        let values = [
            "trailing\\",
            "two\\\\backslashes",
            "graphviz \\N and \\l escapes",
            "a \"quote\" and \\\" an escaped one",
            "backslash\\\nbefore a newline",
            "multi-line\n  indented continuation\n",
            "__linnet_hex_id_5c",
        ];
        let mut graph: DotGraph = DotGraph::from_string("digraph { a -> b; }").unwrap();
        let edge_index = graph.iter_edges().next().unwrap().1;
        for (i, value) in values.into_iter().enumerate() {
            let key = format!("key {i}");
            graph
                .global_data
                .statements
                .insert(key.clone(), value.to_owned());
            graph[NodeIndex(0)]
                .statements
                .insert(key.clone(), value.to_owned());
            graph[edge_index].statements.insert(key, value.to_owned());
        }
        graph.global_data.name = values[0].to_owned();

        let reparsed: DotGraph = DotGraph::from_string(graph.debug_dot()).unwrap();
        assert_eq!(reparsed.global_data, graph.global_data);
        assert_eq!(
            reparsed[NodeIndex(0)].statements,
            graph[NodeIndex(0)].statements
        );
        assert_eq!(
            reparsed.iter_edges().next().unwrap().2.data.statements,
            graph[edge_index].statements
        );
        assert!(graph.debug_dot().contains(r#""key 0" = "trailing\\";"#));
        assert_eq!(reparsed, reparsed.clone().back_and_forth_dot());
    }

    /// The indenting writer used to dedent the continuation lines of a
    /// multi-line value and cut the start of its statement.
    #[test]
    fn multi_line_values_keep_their_statement() {
        let graph: DotGraph =
            DotGraph::from_string("digraph g { a [label=\"x\n  y\"]; a -> b; }").unwrap();
        let written = graph.debug_dot();
        assert!(written.contains("  a[label=\"x\n  y\"];\n"), "{written}");
        let reparsed: DotGraph = DotGraph::from_string(&written).unwrap();
        assert_eq!(reparsed, graph.back_and_forth_dot());
    }

    #[test]
    fn canonical_dot_escapes_half_edge_statements() {
        let source = r#"digraph {
            a -> b [source="source-statement", sink="sink \"quoted\"\N"];
        }"#;
        let mut graph: DotGraph = DotGraph::from_string(source).unwrap();
        let source_hedge = match graph.iter_edges().next().unwrap().0 {
            crate::half_edge::involution::HedgePair::Paired { source, .. } => source,
            _ => unreachable!("the test graph has one paired edge"),
        };
        graph[source_hedge].statement = Some("terminal-backslash\\".to_owned());
        let serialized = graph.debug_dot();
        let reparsed: DotGraph = DotGraph::from_string(&serialized).unwrap();

        assert_eq!(
            reparsed[source_hedge].statement.as_deref(),
            Some("terminal-backslash\\")
        );
        assert_eq!(reparsed, graph.back_and_forth_dot());
        assert!(serialized.contains(r#"source="terminal-backslash\\""#));
        assert!(serialized.contains("sink=\"sink \\\"quoted\\\"\\N\""));
    }

    #[test]
    fn invalid_direction_returns_an_error() {
        for dot in [
            "digraph { a -> b [dir=sideways]; }",
            "digraph { edge [dir=sideways]; a -> b; }",
        ] {
            let parsed: Result<DotGraph, _> = DotGraph::from_string(dot);
            let error = parsed.unwrap_err();
            assert!(matches!(
                error,
                HedgeParseError::DotParseError(DotParseError::InvalidEdgeDirection { value })
                    if value == "sideways"
            ));
        }
    }

    #[test]
    fn graph_sets_return_semantic_errors() {
        let parsed = GraphSet::<_, _, _, _, NodeForest>::from_string(
            "digraph valid { a -> b; } digraph invalid { a -> b [dir=sideways]; }",
        );
        let error = parsed.unwrap_err();

        assert!(matches!(
            error,
            HedgeParseError::DotParseError(DotParseError::InvalidEdgeDirection { value })
                if value == "sideways"
        ));
    }

    #[test]
    fn external_endpoint_data_returns_errors() {
        for (dot, external_source) in [
            ("digraph { a; ext [style=invis]; ext:port -> a; }", true),
            (
                "digraph { a; ext [style=invis]; ext -> a [source=payload]; }",
                true,
            ),
            ("digraph { a; ext [style=invis]; a -> ext:port; }", false),
            (
                "digraph { a; ext [style=invis]; a -> ext [sink=payload]; }",
                false,
            ),
        ] {
            let parsed: Result<DotGraph, _> = DotGraph::from_string(dot);
            let error = parsed.unwrap_err();
            assert!(matches!(
                (error, external_source),
                (
                    HedgeParseError::DotParseError(DotParseError::ExternalSourceEndpointData),
                    true
                ) | (
                    HedgeParseError::DotParseError(DotParseError::ExternalSinkEndpointData),
                    false
                )
            ));
        }
    }

    #[test]
    fn edges_between_external_nodes_return_errors() {
        let parsed: Result<DotGraph, _> = DotGraph::from_string(
            "digraph { source [style=invis]; sink [style=invis]; source -> sink; }",
        );
        let error = parsed.unwrap_err();

        assert!(matches!(
            error,
            HedgeParseError::DotParseError(DotParseError::EdgeBetweenExternalNodes)
        ));
    }

    #[test]
    fn malformed_explicit_ids_return_errors() {
        let cases = [
            ("digraph { a [id=invalid]; }", ExplicitIdKind::Node),
            ("digraph { node [id=invalid]; a; }", ExplicitIdKind::Node),
            ("digraph { a -> b [id=invalid]; }", ExplicitIdKind::Edge),
            (
                "digraph { a:999999999999999999999999999999999999 -> b:0; }",
                ExplicitIdKind::Hedge,
            ),
        ];

        for (dot, expected_kind) in cases {
            let parsed: Result<DotGraph, _> = DotGraph::from_string(dot);
            let error = parsed.unwrap_err();
            assert!(matches!(
                error,
                HedgeParseError::DotParseError(DotParseError::InvalidExplicitId {
                    kind,
                    ..
                }) if kind == expected_kind
            ));
        }
    }

    #[test]
    fn duplicate_explicit_ids_return_errors() {
        let cases = [
            ("digraph { a [id=0]; b [id=0]; }", ExplicitIdKind::Node),
            (
                "digraph { a -> b [id=0]; a -> b [id=0]; }",
                ExplicitIdKind::Edge,
            ),
            ("digraph { a:0 -> b:0; }", ExplicitIdKind::Hedge),
        ];

        for (dot, expected_kind) in cases {
            let parsed: Result<DotGraph, _> = DotGraph::from_string(dot);
            let error = parsed.unwrap_err();
            assert!(matches!(
                error,
                HedgeParseError::DotParseError(DotParseError::DuplicateExplicitId {
                    kind,
                    id: 0,
                }) if kind == expected_kind
            ));
        }
    }

    #[test]
    fn out_of_bounds_explicit_ids_return_errors() {
        let cases = [
            ("digraph { a [id=1]; }", ExplicitIdKind::Node, 1, 1),
            ("digraph { a -> b [id=1]; }", ExplicitIdKind::Edge, 1, 1),
            ("digraph { a:2 -> b:0; }", ExplicitIdKind::Hedge, 2, 2),
        ];

        for (dot, expected_kind, expected_id, expected_len) in cases {
            let parsed: Result<DotGraph, _> = DotGraph::from_string(dot);
            let error = parsed.unwrap_err();
            assert!(matches!(
                error,
                HedgeParseError::DotParseError(DotParseError::ExplicitIdOutOfBounds {
                    kind,
                    id,
                    len,
                }) if kind == expected_kind && id == expected_id && len == expected_len
            ));
        }
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

        let set = GraphSet::<_, _, _, _, NodeForest>::from_string(s).unwrap();
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
        // assert_eq!(graph.n_live(), 2);
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

        let sub: SuBitGraph = aligned.compass_subgraph(Some(CompassPt::S));

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

    /// `ID = ID` and `graph [ID = ID]` set the same graph attributes, and the
    /// later statement wins, as in Graphviz.
    #[test]
    fn later_graph_attributes_win() {
        let graph: DotGraph =
            DotGraph::from_string("digraph { k = 1; graph [k=2, j=1]; j = 2; a }").unwrap();
        let statements = &graph.global_data.statements;
        assert_eq!(
            (statements["k"].as_str(), statements["j"].as_str()),
            ("2", "2")
        );
    }

    /// Edges with equal endpoints keep their statement order. The standard
    /// library sorts up to 20 elements stably even with an unstable sort, so
    /// this needs more edges.
    #[test]
    fn parallel_edges_keep_their_statement_order() {
        let statements: String = (0..64)
            .map(|n| match n % 2 {
                0 => format!("a -> b [n={n}]; "),
                _ => format!("b -> a [n={n}]; "),
            })
            .collect();
        let graph: DotGraph = DotGraph::from_string(format!("digraph {{ {statements} }}")).unwrap();
        let order: Vec<usize> = (graph.iter_edges())
            .map(|(_, _, edge)| edge.data.statements["n"].parse().unwrap())
            .collect();
        let expected: Vec<usize> = (0..64).step_by(2).chain((1..64).step_by(2)).collect();
        assert_eq!(order, expected);
    }
}
