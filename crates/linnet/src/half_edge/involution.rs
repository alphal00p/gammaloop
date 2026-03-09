use std::{
    fmt::Display,
    num::ParseIntError,
    ops::{Index, IndexMut, Mul, Neg},
    str::FromStr,
};

use crate::{define_indexed_vec, num_traits::RefZero};
use rand::{rngs::SmallRng, Rng, SeedableRng};
use thiserror::Error;

use super::{
    builder::HedgeData,
    nodestore::{NodeStorage, NodeStorageOps},
    subgraph::SubSetLike,
    swap::Swap,
    GVEdgeAttrs, HedgeGraph, NodeIndex,
};

define_indexed_vec!(
    /// A type-safe wrapper around a `usize` representing an undirected edge in the graph.
    ///
    /// An `EdgeIndex` typically refers to one of the two `Hedge`s that constitute
    /// the full edge. The specific `Hedge` it points to (e.g., the "source" or "canonical"
    /// half-edge) depends on the conventions used by the `Involution` or `HedgeGraph`.
    /// It's used to retrieve data common to both half-edges of an edge pair.
    pub struct EdgeIndex; // new‑type around usize

    /// Vector whose items are `T`, indexable only by `EntityId`.
    pub struct EdgeVec;
);

impl Display for EdgeIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "e{}", self.0)
    }
}
define_indexed_vec!(
    /// A type-safe wrapper around a `usize` representing a directed half-edge.
    ///
    /// Each undirected edge in a graph is typically composed of two oppositely
    /// directed half-edges. `Hedge` identifies one of these. The actual edge data
    /// and topological information (like its opposite or next half-edge) are
    /// stored within the graph structure (e.g., [`HedgeGraph`]).
    pub struct Hedge; // new‑type around usize

    /// Vector whose items are `T`, indexable only by `EntityId`.
    pub struct HedgeVec;
);

impl<T: Display> Display for HedgeVec<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HedgeVec:{{")?;
        for (i, (_, item)) in self.iter().enumerate() {
            if i > 0 {
                write!(f, ", ")?;
            }
            write!(f, "{item}")?;
        }
        write!(f, "}}")
    }
}

impl Hedge {
    pub fn to_edge_id<E>(self, involution: &Involution<E>) -> HedgePair {
        HedgePair::from_half_edge(self, involution)
    }

    pub fn to_edge_id_with_subgraph<E, S: SubSetLike>(
        self,
        involution: &Involution<E>,
        subgraph: &S,
    ) -> Option<HedgePair> {
        HedgePair::from_half_edge_with_subgraph(self, involution, subgraph)
    }
}

impl FromStr for Hedge {
    type Err = ParseIntError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.parse::<usize>().map(Hedge)
    }
}

pub enum HedgePairWithData<H> {
    Unpaired {
        hedge: Hedge,
        flow: Flow,
        data: H,
        is_in_subgraph: bool,
    },
    Paired {
        source: Hedge,
        source_data: H,
        sink: Hedge,
        sink_data: H,
        is_in_subgraph: bool,
    },
    Split {
        source: Hedge,
        source_data: H,
        sink: Hedge,
        sink_data: H,
        split: Flow,
    },
}

impl HedgeData<(Hedge, Option<String>)> {
    pub fn dot_node(&self) -> String {
        format!(
            "{}:{}{}",
            self.node,
            self.data.0,
            if self.is_in_subgraph { ":s" } else { "" },
        )
    }

    pub fn dot_node_map(&self, map: impl Fn(NodeIndex) -> String) -> String {
        format!(
            "{}:{}{}",
            map(self.node),
            self.data.0,
            if self.is_in_subgraph { ":s" } else { "" },
        )
    }
}

#[allow(clippy::too_many_arguments)]
impl<H> HedgePairWithData<H> {
    pub fn identity_dot_io<W: std::io::Write>(
        writer: &mut W,
        edge_id: EdgeIndex,
        hedge: HedgeData<(Hedge, Option<String>)>,
        map: impl Fn(NodeIndex) -> String,
        attr: &GVEdgeAttrs,
        orientation: Orientation,
        flow: Flow,
    ) -> Result<(), std::io::Error> {
        //we interpret this node as the hidden one. so the flow is reversed.
        match flow {
            Flow::Sink => {
                writeln!(writer, "\text{edge_id}\t [style=invis];",)?;
                write!(
                    writer,
                    "\text{edge_id}\t-> {}\t [",
                    hedge.dot_node_map(&map)
                )?;
            }
            Flow::Source => {
                writeln!(writer, "\text{edge_id}\t [style=invis];")?;
                write!(
                    writer,
                    "\t{}\t-> ext{edge_id}\t [",
                    hedge.dot_node_map(&map)
                )?;
            }
        }

        write!(writer, "id={}", edge_id.0)?;

        match orientation {
            Orientation::Reversed => write!(writer, " dir=back")?,
            Orientation::Undirected => write!(writer, " dir=none")?,
            _ => {}
        }
        if let Some(hedge_label) = hedge.data.1 {
            match flow {
                Flow::Sink => write!(writer, " sink={hedge_label}")?,
                Flow::Source => write!(writer, " source={hedge_label}")?,
            }
        }
        write!(writer, " {attr}")?;

        writeln!(writer, "];")?;
        Ok(())
    }

    pub fn identity_dot_fmt<W: std::fmt::Write>(
        writer: &mut W,
        edge_id: EdgeIndex,
        hedge: HedgeData<(Hedge, Option<String>)>,
        map: impl Fn(NodeIndex) -> String,
        attr: &GVEdgeAttrs,
        orientation: Orientation,
        flow: Flow,
    ) -> Result<(), std::fmt::Error> {
        //we interpret this node as the hidden one. so the flow is reversed.
        let estr = match flow {
            Flow::Sink => {
                write!(writer, "\next{}\t [style=invis];", edge_id.0)?;
                format!("ext{}\t-> {}", edge_id.0, hedge.dot_node_map(&map))
            }
            Flow::Source => {
                write!(writer, "\next{}\t [style=invis];", edge_id.0)?;
                format!("{}\t-> ext{}", hedge.dot_node_map(&map), edge_id.0,)
            }
        };

        write!(
            writer,
            "\n{estr}\t [id={}{}{}{}];",
            edge_id.0,
            match orientation {
                Orientation::Reversed => " dir=back",
                Orientation::Undirected => " dir=none",
                _ => "",
            },
            if let Some(hedge_label) = hedge.data.1 {
                match flow {
                    Flow::Sink => format!(" sink={hedge_label}"),
                    Flow::Source => format!(" source={hedge_label}"),
                }
            } else {
                "".to_string()
            },
            attr
        )
    }

    pub fn pair_dot_fmt<W: std::fmt::Write>(
        writer: &mut W,
        eid: EdgeIndex,
        source_data: HedgeData<(Hedge, Option<String>)>,
        sink_data: HedgeData<(Hedge, Option<String>)>,
        map: impl Fn(NodeIndex) -> String,
        attr: &GVEdgeAttrs,
        orientation: Orientation,
    ) -> Result<(), std::fmt::Error> {
        write!(
            writer,
            "\n{}\t-> {}\t [id={}{}{}{} {}];",
            source_data.dot_node_map(&map),
            sink_data.dot_node_map(&map),
            eid.0,
            match orientation {
                Orientation::Reversed => " dir=back",
                Orientation::Undirected => " dir=none",
                _ => "",
            },
            if let Some(source_label) = source_data.data.1 {
                format!(" source={source_label}")
            } else {
                "".to_string()
            },
            if let Some(sink_label) = sink_data.data.1 {
                format!(" sink={sink_label}")
            } else {
                "".to_string()
            },
            attr
        )
    }

    pub fn pair_dot_io<W: std::io::Write>(
        writer: &mut W,
        eid: EdgeIndex,
        source_data: HedgeData<(Hedge, Option<String>)>,
        sink_data: HedgeData<(Hedge, Option<String>)>,
        map: impl Fn(NodeIndex) -> String,
        attr: &GVEdgeAttrs,
        orientation: Orientation,
    ) -> Result<(), std::io::Error> {
        writeln!(
            writer,
            "\t{}\t-> {}\t [id={}{}{}{} {}];",
            source_data.dot_node_map(&map),
            sink_data.dot_node_map(&map),
            eid.0,
            match orientation {
                Orientation::Reversed => " dir=back",
                Orientation::Undirected => " dir=none",
                _ => "",
            },
            if let Some(source_label) = source_data.data.1 {
                format!(" source={source_label}")
            } else {
                "".to_string()
            },
            if let Some(sink_label) = sink_data.data.1 {
                format!(" sink={sink_label}")
            } else {
                "".to_string()
            },
            attr
        )
    }
}
impl<H> HedgePairWithData<&H> {
    #[allow(clippy::too_many_arguments)]
    pub fn dot_fmt<W: std::fmt::Write, E, V, N: NodeStorageOps<NodeData = V>>(
        &self,
        writer: &mut W,
        graph: &HedgeGraph<E, V, H, N>,
        eid: EdgeIndex,
        labeler: impl Fn(&H) -> Option<String>,
        ider: impl Fn(NodeIndex) -> String,
        orientation: Orientation,
        attr: GVEdgeAttrs,
    ) -> Result<(), std::fmt::Error> {
        match self {
            HedgePairWithData::Unpaired {
                hedge,
                flow,
                data,
                is_in_subgraph,
            } => Self::identity_dot_fmt(
                writer,
                eid,
                graph
                    .node_id(*hedge)
                    .add_data((*hedge, labeler(data)))
                    .set_in_subgraph(*is_in_subgraph),
                &ider,
                &attr,
                orientation,
                *flow,
            ),
            HedgePairWithData::Paired {
                source,
                sink,
                sink_data,
                source_data,
                is_in_subgraph,
            } => Self::pair_dot_fmt(
                writer,
                eid,
                graph
                    .node_id(*source)
                    .add_data((*source, labeler(source_data)))
                    .set_in_subgraph(*is_in_subgraph),
                graph
                    .node_id(*sink)
                    .add_data((*sink, labeler(sink_data)))
                    .set_in_subgraph(*is_in_subgraph),
                &ider,
                &attr,
                orientation,
            ),
            HedgePairWithData::Split {
                source,
                sink,
                sink_data,
                source_data,
                split,
            } => Self::pair_dot_fmt(
                writer,
                eid,
                graph
                    .node_id(*source)
                    .add_data((*source, labeler(source_data)))
                    .set_in_subgraph(*split == Flow::Source),
                graph
                    .node_id(*sink)
                    .add_data((*sink, labeler(sink_data)))
                    .set_in_subgraph(*split == Flow::Sink),
                &ider,
                &attr,
                orientation,
            ),
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn dot_io<W: std::io::Write, E, V, N: NodeStorageOps<NodeData = V>>(
        &self,
        writer: &mut W,
        graph: &HedgeGraph<E, V, H, N>,
        eid: EdgeIndex,
        labeler: impl Fn(&H) -> Option<String>,
        ider: impl Fn(NodeIndex) -> String,
        orientation: Orientation,
        attr: GVEdgeAttrs,
    ) -> Result<(), std::io::Error> {
        match self {
            HedgePairWithData::Unpaired {
                hedge,
                flow,
                data,
                is_in_subgraph,
            } => Self::identity_dot_io(
                writer,
                eid,
                graph
                    .node_id(*hedge)
                    .add_data((*hedge, labeler(data)))
                    .set_in_subgraph(*is_in_subgraph),
                &ider,
                &attr,
                orientation,
                *flow,
            ),
            HedgePairWithData::Paired {
                source,
                sink,
                sink_data,
                source_data,
                is_in_subgraph,
            } => Self::pair_dot_io(
                writer,
                eid,
                graph
                    .node_id(*source)
                    .add_data((*source, labeler(source_data)))
                    .set_in_subgraph(*is_in_subgraph),
                graph
                    .node_id(*sink)
                    .add_data((*sink, labeler(sink_data)))
                    .set_in_subgraph(*is_in_subgraph),
                &ider,
                &attr,
                orientation,
            ),
            HedgePairWithData::Split {
                source,
                sink,
                sink_data,
                source_data,
                split,
            } => Self::pair_dot_io(
                writer,
                eid,
                graph
                    .node_id(*source)
                    .add_data((*source, labeler(source_data)))
                    .set_in_subgraph(*split == Flow::Source),
                graph
                    .node_id(*sink)
                    .add_data((*sink, labeler(sink_data)))
                    .set_in_subgraph(*split == Flow::Source),
                &ider,
                &attr,
                orientation,
            ),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// Describes the pairing status of a half-edge, indicating whether it forms a
/// complete edge with another half-edge, is a dangling/external edge, or is
/// part of a "split" edge in a subgraph context.
pub enum HedgePair {
    /// Represents a half-edge that is unpaired (dangling or external).
    Unpaired {
        /// The half-edge itself.
        hedge: Hedge,
        /// The [`Flow`] (directionality) of this unpaired half-edge.
        flow: Flow,
    },
    /// Represents two half-edges that are paired together to form a complete,
    /// internal, undirected edge.
    Paired {
        /// The half-edge considered the "source" in this pair.
        source: Hedge,
        /// The half-edge considered the "sink" in this pair, opposite to `source`.
        sink: Hedge,
    },
    /// Represents an edge that is "split" in the context of a subgraph.
    /// This typically means one half-edge is inside the subgraph and its
    /// opposite is outside, or vice-versa.
    Split {
        /// The half-edge on one side of the split (often the one inside a subgraph).
        source: Hedge,
        /// The half-edge on the other side of the split (often the one outside a subgraph).
        sink: Hedge,
        /// Indicates which side of the original conceptual edge (`source` or `sink`
        /// of the pair) is considered "split off" or external.
        /// For example, if `split == Flow::Source`, then the `source` half-edge
        /// is effectively dangling from the perspective of the other side of the split.
        split: Flow,
    },
}

impl Display for HedgePair {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HedgePair::Paired { source, sink } => {
                write!(f, "({source} -> {sink})")
            }
            HedgePair::Split {
                source,
                sink,
                split,
            } => match split {
                Flow::Source => {
                    write!(f, "([{source}] -> {sink})")
                }
                Flow::Sink => {
                    write!(f, "({sink} -> [{source}])")
                }
            },
            HedgePair::Unpaired { hedge, flow } => match flow {
                Flow::Source => {
                    write!(f, "([{hedge}] ->)")
                }
                Flow::Sink => {
                    write!(f, "( -> [{hedge}])")
                }
            },
        }
    }
}

impl HedgePair {
    pub fn add_data_of<'a, S: SubSetLike, V, E, H, N: NodeStorage<NodeData = V>>(
        self,
        graph: &'a HedgeGraph<E, V, H, N>,
        subgraph: &S,
    ) -> HedgePairWithData<&'a H> {
        match self {
            HedgePair::Paired { source, sink } => HedgePairWithData::Paired {
                source,
                source_data: &graph[source],
                is_in_subgraph: subgraph.includes(&source) && subgraph.includes(&sink),
                sink,
                sink_data: &graph[sink],
            },

            HedgePair::Split {
                source,
                sink,
                split,
            } => HedgePairWithData::Split {
                source,
                source_data: &graph[source],
                sink,
                sink_data: &graph[sink],
                split,
            },
            HedgePair::Unpaired { hedge, flow } => HedgePairWithData::Unpaired {
                hedge,
                flow,
                data: &graph[hedge],
                is_in_subgraph: subgraph.includes(&hedge),
            },
        }
    }

    pub fn add_data<V, E, H, N: NodeStorage<NodeData = V>>(
        self,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> HedgePairWithData<&H> {
        match self {
            HedgePair::Paired { source, sink } => HedgePairWithData::Paired {
                source,
                source_data: &graph[source],
                is_in_subgraph: false,
                sink,
                sink_data: &graph[sink],
            },

            HedgePair::Split {
                source,
                sink,
                split,
            } => HedgePairWithData::Split {
                source,
                source_data: &graph[source],
                sink,
                sink_data: &graph[sink],
                split,
            },
            HedgePair::Unpaired { hedge, flow } => HedgePairWithData::Unpaired {
                hedge,
                flow,
                data: &graph[hedge],
                is_in_subgraph: false,
            },
        }
    }

    pub fn is_unpaired(&self) -> bool {
        matches!(self, HedgePair::Unpaired { .. })
    }

    pub fn is_paired(&self) -> bool {
        matches!(self, HedgePair::Paired { .. })
    }

    pub fn is_split(&self) -> bool {
        matches!(self, HedgePair::Split { .. })
    }
    pub fn swap(&mut self) {
        match self {
            HedgePair::Unpaired { flow, .. } => {
                *flow = -*flow;
            }
            HedgePair::Paired { source, sink } => {
                std::mem::swap(&mut (*source), &mut (*sink));
            }
            HedgePair::Split {
                source,
                sink,
                split,
            } => {
                std::mem::swap(&mut (*source), &mut (*sink));
                *split = -*split;
            }
        }
    }
    pub fn any_hedge(&self) -> Hedge {
        match self {
            HedgePair::Unpaired { hedge, .. } => *hedge,
            HedgePair::Paired { source, .. } => *source,
            HedgePair::Split {
                source,
                sink,
                split,
            } => {
                if *split == Flow::Source {
                    *source
                } else {
                    *sink
                }
            }
        }
    }

    pub fn from_half_edge<E>(hedge: Hedge, involution: &Involution<E>) -> Self {
        match involution.hedge_data(hedge) {
            InvolutiveMapping::Source { sink_idx, .. } => HedgePair::Paired {
                source: hedge,
                sink: *sink_idx,
            },
            InvolutiveMapping::Sink { source_idx } => HedgePair::Paired {
                source: *source_idx,
                sink: hedge,
            },
            InvolutiveMapping::Identity { underlying, .. } => HedgePair::Unpaired {
                hedge,
                flow: *underlying,
            },
        }
    }

    pub fn from_source<E>(hedge: Hedge, involution: &Involution<E>) -> Option<Self> {
        match involution.hedge_data(hedge) {
            InvolutiveMapping::Source { sink_idx, .. } => Some(HedgePair::Paired {
                source: hedge,
                sink: *sink_idx,
            }),
            InvolutiveMapping::Identity { underlying, .. } => Some(HedgePair::Unpaired {
                hedge,
                flow: *underlying,
            }),
            _ => None,
        }
    }

    pub fn from_source_with_subgraph<E, S: SubSetLike, I: AsRef<Involution<E>>>(
        hedge: Hedge,
        involution: I,
        subgraph: &S,
    ) -> Option<Self> {
        if subgraph.includes(&hedge) {
            match involution.as_ref().hedge_data(hedge) {
                InvolutiveMapping::Source { sink_idx, .. } => {
                    if subgraph.includes(sink_idx) {
                        Some(HedgePair::Paired {
                            source: hedge,
                            sink: *sink_idx,
                        })
                    } else {
                        Some(HedgePair::Split {
                            source: hedge,
                            sink: *sink_idx,
                            split: Flow::Source,
                        })
                    }
                }
                InvolutiveMapping::Sink { source_idx } => {
                    if subgraph.includes(source_idx) {
                        None
                    } else {
                        Some(HedgePair::Split {
                            source: *source_idx,
                            sink: hedge,
                            split: Flow::Sink,
                        })
                    }
                }
                InvolutiveMapping::Identity { underlying, .. } => Some(HedgePair::Unpaired {
                    hedge,
                    flow: *underlying,
                }),
            }
        } else {
            None
        }
    }

    // uses the present color for split edge background
    pub fn fill_color(&self, attr: GVEdgeAttrs) -> GVEdgeAttrs {
        match self {
            HedgePair::Unpaired { flow, .. } => {
                if attr.color.is_some() {
                    attr
                } else {
                    GVEdgeAttrs {
                        color: Some(flow.color().to_owned()),
                        label: attr.label,
                        other: attr.other,
                    }
                }
            }
            HedgePair::Paired { .. } => {
                if attr.color.is_some() {
                    attr
                } else {
                    let color = format!("{}:{};0.5", Flow::Source.color(), Flow::Sink.color());
                    GVEdgeAttrs {
                        color: Some(color),
                        label: attr.label,
                        other: attr.other,
                    }
                }
            }
            HedgePair::Split { split, .. } => {
                let background = if let Some(color) = attr.color {
                    color
                } else {
                    "gray75".to_owned()
                };

                let color = match split {
                    Flow::Source => format!("{}:{};0.5", split.color(), background),
                    Flow::Sink => format!("{}:{};0.5", background, split.color()),
                };
                GVEdgeAttrs {
                    color: Some(color),
                    label: attr.label,
                    other: attr.other,
                }
            }
        }
    }

    pub fn with_subgraph<S: SubSetLike>(self, subgraph: &S) -> Option<Self> {
        match self {
            HedgePair::Paired { source, sink } => {
                match (subgraph.includes(&source), subgraph.includes(&sink)) {
                    (true, true) => Some(HedgePair::Paired { source, sink }),
                    (true, false) => Some(HedgePair::Split {
                        source,
                        sink,
                        split: Flow::Source,
                    }),
                    (false, true) => Some(HedgePair::Split {
                        source,
                        sink,
                        split: Flow::Sink,
                    }),
                    (false, false) => None,
                }
            }
            HedgePair::Split { source, sink, .. } => {
                match (subgraph.includes(&source), subgraph.includes(&sink)) {
                    (true, true) => Some(HedgePair::Paired { source, sink }),
                    (true, false) => Some(HedgePair::Split {
                        source,
                        sink,
                        split: Flow::Source,
                    }),
                    (false, true) => Some(HedgePair::Split {
                        source,
                        sink,
                        split: Flow::Sink,
                    }),
                    (false, false) => None,
                }
            }
            HedgePair::Unpaired { hedge, flow } => {
                if subgraph.includes(&hedge) {
                    Some(HedgePair::Unpaired { hedge, flow })
                } else {
                    None
                }
            }
        }
    }

    pub fn from_half_edge_with_subgraph<E, S: SubSetLike>(
        hedge: Hedge,
        involution: &Involution<E>,
        subgraph: &S,
    ) -> Option<Self> {
        if subgraph.includes(&hedge) {
            Some(match involution.hedge_data(hedge) {
                InvolutiveMapping::Source { sink_idx, .. } => {
                    if subgraph.includes(sink_idx) {
                        HedgePair::Paired {
                            source: hedge,
                            sink: *sink_idx,
                        }
                    } else {
                        HedgePair::Split {
                            source: hedge,
                            sink: *sink_idx,
                            split: Flow::Source,
                        }
                    }
                }
                InvolutiveMapping::Sink { source_idx } => {
                    if subgraph.includes(source_idx) {
                        HedgePair::Paired {
                            source: *source_idx,
                            sink: hedge,
                        }
                    } else {
                        HedgePair::Split {
                            source: *source_idx,
                            sink: hedge,
                            split: Flow::Sink,
                        }
                    }
                }
                InvolutiveMapping::Identity { underlying, .. } => HedgePair::Unpaired {
                    hedge,
                    flow: *underlying,
                },
            })
        } else {
            None
        }
    }
}

impl Display for Hedge {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// Represents the state of a half-edge within an `Involution`.
///
/// An `Involution` maps each half-edge to another, defining the graph's topology.
/// This enum describes whether a half-edge is part of a pair forming a full edge
/// (either as a `Source` or a `Sink`), or if it's an unpaired `Identity` half-edge
/// (representing a dangling or external connection).
pub enum InvolutiveMapping<E> {
    /// The half-edge is an "identity" or "external" half-edge, meaning it is
    /// not paired with another half-edge to form a complete internal edge.
    /// It points back to itself under the involution.
    Identity {
        /// The data associated with this half-edge.
        data: EdgeData<E>,
        /// The underlying [`Flow`] (directionality) of this identity half-edge.
        underlying: Flow,
    },
    /// The half-edge is a "source" half-edge, meaning it is paired with a "sink"
    /// half-edge to form a complete internal edge. It stores the edge data.
    Source {
        /// The data associated with the edge this half-edge belongs to.
        data: EdgeData<E>,
        /// The [`Hedge`] index of the corresponding "sink" half-edge.
        sink_idx: Hedge,
    },
    /// The half-edge is a "sink" half-edge. It does not store data itself,
    /// as the data is stored by its corresponding "source" half-edge.
    Sink {
        /// The [`Hedge`] index of the corresponding "source" half-edge.
        source_idx: Hedge,
    },
}

impl<E> InvolutiveMapping<E> {
    pub fn flow(&self) -> Flow {
        match self {
            InvolutiveMapping::Identity { underlying, .. } => *underlying,
            InvolutiveMapping::Source { .. } => Flow::Source,
            InvolutiveMapping::Sink { .. } => Flow::Sink,
        }
    }

    fn dummy() -> Self {
        InvolutiveMapping::Sink {
            source_idx: Hedge(0),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// Holds the data associated with an edge, including its [`Orientation`]
/// and the custom data of type `E`.
///
/// This struct is typically stored by the "source" half-edge of an edge pair,
/// or by an "identity" (unpaired) half-edge.
pub struct EdgeData<E> {
    /// The superficial [`Orientation`] of the edge (e.g., Default, Reversed, Undirected).
    pub orientation: Orientation,
    /// The custom data of type `E` associated with the edge.
    pub data: E,
}

impl<E: Display> Display for EdgeData<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.orientation {
            Orientation::Default => write!(f, "+{}", self.data),
            Orientation::Reversed => write!(f, "-{}", self.data),
            Orientation::Undirected => write!(f, "{}", self.data),
        }
    }
}

impl<E> From<Orientation> for EdgeData<Option<E>> {
    fn from(orientation: Orientation) -> Self {
        EdgeData::empty(orientation)
    }
}

impl<E> EdgeData<Option<E>> {
    pub fn is_set(&self) -> bool {
        self.data.is_some()
    }
    pub const fn none() -> Self {
        Self {
            orientation: Orientation::Undirected,
            data: None,
        }
    }

    pub fn permute(self) -> Option<EdgeData<E>> {
        Some(EdgeData {
            orientation: self.orientation,
            data: self.data?,
        })
    }
    pub fn empty(orientation: Orientation) -> Self {
        EdgeData {
            orientation,
            data: None,
        }
    }
    pub fn take(&mut self) -> Self {
        EdgeData {
            orientation: self.orientation,
            data: self.data.take(),
        }
    }
    pub fn and_then<U, F>(self, f: F) -> EdgeData<Option<U>>
    where
        F: FnOnce(E) -> U,
    {
        match self.data {
            Some(x) => EdgeData {
                data: Some(f(x)),
                orientation: self.orientation,
            },
            None => EdgeData {
                data: None,
                orientation: self.orientation,
            },
        }
    }
}

impl<E> EdgeData<E> {
    pub fn un_orient(&mut self) {
        self.orientation = Orientation::Undirected;
    }

    pub fn orient(&mut self) {
        if let Orientation::Undirected = self.orientation {
            self.orientation = Orientation::Default;
        }
    }
    pub fn flip_orientation(&mut self) {
        self.orientation = self.orientation.reverse();
    }
    pub fn new(data: E, orientation: Orientation) -> Self {
        EdgeData { data, orientation }
    }
    pub fn map<E2>(self, f: impl FnOnce(E) -> E2) -> EdgeData<E2> {
        EdgeData {
            orientation: self.orientation,
            data: f(self.data),
        }
    }

    pub fn map_option<E2>(self, f: impl FnOnce(E) -> Option<E2>) -> Option<EdgeData<E2>> {
        Some(EdgeData {
            orientation: self.orientation,
            data: f(self.data)?,
        })
    }

    pub fn map_result<E2, O>(self, f: impl FnOnce(E) -> Result<E2, O>) -> Result<EdgeData<E2>, O> {
        Ok(EdgeData {
            orientation: self.orientation,
            data: f(self.data)?,
        })
    }

    pub const fn as_ref(&self) -> EdgeData<&E> {
        EdgeData {
            orientation: self.orientation,
            data: &self.data,
        }
    }
}

#[derive(Clone, Debug, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// Represents the superficial or conventional orientation of an edge.
///
/// This orientation can be distinct from the underlying `Flow` of the half-edges
/// that constitute the edge. It's often used for display or for algorithms
/// where a conventional direction is needed.
pub enum Orientation {
    /// The edge's orientation aligns with the underlying flow of its "source"
    /// half-edge (e.g., if half-edge A is `Flow::Source` and points to B,
    /// a `Default` orientation means the edge is A -> B).
    Default,
    /// The edge's orientation is opposite to the underlying flow of its "source"
    /// half-edge (e.g., if half-edge A is `Flow::Source` and points to B,
    /// a `Reversed` orientation means the edge is B -> A).
    Reversed,
    /// The edge is considered undirected, regardless of the underlying half-edge flows.
    Undirected,
}

impl Mul for Orientation {
    type Output = Orientation;

    fn mul(self, other: Orientation) -> Self::Output {
        match (self, other) {
            (Orientation::Default, _) => other,
            (Orientation::Undirected, _) => Orientation::Undirected,
            (_, Orientation::Default) => self,
            (_, Orientation::Undirected) => Orientation::Undirected,
            (Orientation::Reversed, Orientation::Reversed) => Orientation::Default,
        }
    }
}

#[derive(Clone, Debug, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// Represents the underlying, intrinsic directionality of a half-edge
/// relative to the full edge it is part of.
///
/// Each full edge is composed of two half-edges. One is designated as the
/// `Source` and the other as the `Sink`. This is a fundamental property
/// used by the `Involution` to pair half-edges correctly.
pub enum Flow {
    /// Indicates this half-edge is the "source" part of its full edge.
    /// This is often visualized as the outgoing part of a directed edge
    /// if the edge's `Orientation` is `Default`.
    Source, // outgoing
    /// Indicates this half-edge is the "sink" part of its full edge.
    /// This is often visualized as the incoming part of a directed edge
    /// if the edge's `Orientation` is `Default`.
    Sink, // incoming
}

impl Mul for Flow {
    type Output = Flow;

    fn mul(self, other: Flow) -> Self::Output {
        match (self, other) {
            (Flow::Source, _) => other,
            (Flow::Sink, _) => -other,
        }
    }
}

impl Flow {
    pub fn color(&self) -> &'static str {
        match self {
            Flow::Source => "red",
            Flow::Sink => "blue",
        }
    }
}

impl Neg for Flow {
    type Output = Flow;

    fn neg(self) -> Self::Output {
        match self {
            Flow::Source => Flow::Sink,
            Flow::Sink => Flow::Source,
        }
    }
}

impl From<bool> for Flow {
    fn from(value: bool) -> Self {
        if value {
            Flow::Source
        } else {
            Flow::Sink
        }
    }
}

impl From<Flow> for Orientation {
    fn from(value: Flow) -> Self {
        match value {
            Flow::Source => Orientation::Default,
            Flow::Sink => Orientation::Reversed,
        }
    }
}

impl From<Flow> for SignOrZero {
    fn from(value: Flow) -> Self {
        match value {
            Flow::Source => SignOrZero::Minus,
            Flow::Sink => SignOrZero::Plus,
        }
    }
}

impl TryFrom<Orientation> for Flow {
    type Error = &'static str;

    fn try_from(value: Orientation) -> Result<Self, Self::Error> {
        match value {
            Orientation::Default => Ok(Flow::Source),
            Orientation::Reversed => Ok(Flow::Sink),
            Orientation::Undirected => Err("Cannot convert undirected orientation to flow"),
        }
    }
}

impl Orientation {
    pub fn reverse(self) -> Orientation {
        match self {
            Orientation::Default => Orientation::Reversed,
            Orientation::Reversed => Orientation::Default,
            Orientation::Undirected => Orientation::Undirected,
        }
    }

    pub fn relative_to(&self, other: Flow) -> Orientation {
        match (self, other) {
            (Orientation::Default, Flow::Source) => Orientation::Default,
            (Orientation::Default, Flow::Sink) => Orientation::Reversed,
            (Orientation::Reversed, Flow::Source) => Orientation::Reversed,
            (Orientation::Reversed, Flow::Sink) => Orientation::Default,
            (Orientation::Undirected, _) => Orientation::Undirected,
        }
    }
}

impl From<bool> for Orientation {
    fn from(value: bool) -> Self {
        if value {
            Orientation::Default
        } else {
            Orientation::Undirected
        }
    }
}

#[derive(Error, Debug)]
pub enum SignError {
    #[error("Invalid value for Sign")]
    InvalidValue,
    #[error("Zero is not a valid value for Sign")]
    ZeroValue,
}
#[repr(i8)]
/// Represents a sign (Plus or Minus) or zero.
///
/// Useful for calculations where direction or orientation matters,
/// and a neutral (zero) state is also possible.
pub enum SignOrZero {
    /// Represents a zero or neutral value.
    Zero = 0,
    /// Represents a positive sign or direction.
    Plus = 1,
    /// Represents a negative sign or direction.
    Minus = -1,
}

impl TryFrom<i8> for SignOrZero {
    type Error = SignError;
    fn try_from(value: i8) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(SignOrZero::Zero),
            1 => Ok(SignOrZero::Plus),
            -1 => Ok(SignOrZero::Minus),
            _ => Err(SignError::InvalidValue),
        }
    }
}

impl Display for SignOrZero {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SignOrZero::Zero => write!(f, "."),
            SignOrZero::Plus => write!(f, "+"),
            SignOrZero::Minus => write!(f, "-"),
        }
    }
}

impl From<SignOrZero> for Orientation {
    fn from(value: SignOrZero) -> Self {
        match value {
            SignOrZero::Zero => Orientation::Undirected,
            SignOrZero::Plus => Orientation::Default,
            SignOrZero::Minus => Orientation::Reversed,
        }
    }
}

#[allow(non_upper_case_globals)]
impl SignOrZero {
    pub fn is_zero(&self) -> bool {
        matches!(self, SignOrZero::Zero)
    }

    pub fn is_sign(&self) -> bool {
        matches!(self, SignOrZero::Plus | SignOrZero::Minus)
    }

    pub fn is_positive(&self) -> bool {
        matches!(self, SignOrZero::Plus)
    }

    pub fn is_negative(&self) -> bool {
        matches!(self, SignOrZero::Minus)
    }
}

impl<T: Neg<Output = T> + RefZero> Mul<T> for SignOrZero {
    type Output = T;
    fn mul(self, rhs: T) -> Self::Output {
        match self {
            SignOrZero::Plus => rhs,
            SignOrZero::Minus => -rhs,
            SignOrZero::Zero => rhs.ref_zero(),
        }
    }
}

impl Neg for SignOrZero {
    type Output = Self;
    fn neg(self) -> Self::Output {
        match self {
            SignOrZero::Plus => SignOrZero::Minus,
            SignOrZero::Minus => SignOrZero::Plus,
            SignOrZero::Zero => SignOrZero::Zero,
        }
    }
}

impl From<Orientation> for SignOrZero {
    fn from(value: Orientation) -> Self {
        match value {
            Orientation::Default => SignOrZero::Plus,
            Orientation::Reversed => SignOrZero::Minus,
            Orientation::Undirected => SignOrZero::Zero,
        }
    }
}

impl<E> InvolutiveMapping<Option<E>> {}
impl<E> InvolutiveMapping<E> {
    pub fn is_identity(&self) -> bool {
        matches!(self, InvolutiveMapping::Identity { .. })
    }

    pub fn data(self) -> Option<EdgeData<E>> {
        match self {
            InvolutiveMapping::Identity { data, .. } => Some(data),
            InvolutiveMapping::Source { data, .. } => Some(data),
            _ => None,
        }
    }

    pub fn to_sign(&self) -> SignOrZero {
        match self {
            InvolutiveMapping::Identity { .. } => SignOrZero::Zero,
            InvolutiveMapping::Sink { .. } => SignOrZero::Plus, // incoming
            InvolutiveMapping::Source { .. } => SignOrZero::Minus,
        }
    }

    pub fn as_ref(&self) -> InvolutiveMapping<&E> {
        self.map_data_ref(|e| e)
    }

    pub fn map_data_ref<'a, E2>(&'a self, f: impl FnOnce(&'a E) -> E2) -> InvolutiveMapping<E2> {
        match self {
            InvolutiveMapping::Identity { data, underlying } => InvolutiveMapping::Identity {
                data: data.as_ref().map(f),
                underlying: *underlying,
            },
            InvolutiveMapping::Source { data, sink_idx } => InvolutiveMapping::Source {
                data: data.as_ref().map(f),
                sink_idx: *sink_idx,
            },
            InvolutiveMapping::Sink { source_idx } => InvolutiveMapping::Sink {
                source_idx: *source_idx,
            },
        }
    }

    pub fn map_data<E2>(self, f: impl FnOnce(E) -> E2) -> InvolutiveMapping<E2> {
        match self {
            InvolutiveMapping::Identity { data, underlying } => InvolutiveMapping::Identity {
                data: data.map(f),
                underlying,
            },
            InvolutiveMapping::Source { data, sink_idx } => InvolutiveMapping::Source {
                data: data.map(f),
                sink_idx,
            },
            InvolutiveMapping::Sink { source_idx } => InvolutiveMapping::Sink { source_idx },
        }
    }

    pub fn map_data_option<E2>(
        self,
        f: impl FnOnce(E) -> Option<E2>,
    ) -> Option<InvolutiveMapping<E2>> {
        match self {
            InvolutiveMapping::Identity { data, underlying } => Some(InvolutiveMapping::Identity {
                data: data.map_option(f)?,
                underlying,
            }),
            InvolutiveMapping::Source { data, sink_idx } => Some(InvolutiveMapping::Source {
                data: data.map_option(f)?,
                sink_idx,
            }),
            InvolutiveMapping::Sink { source_idx } => Some(InvolutiveMapping::Sink { source_idx }),
        }
    }

    pub fn map_data_result<E2, O>(
        self,
        f: impl FnOnce(E) -> Result<E2, O>,
    ) -> Result<InvolutiveMapping<E2>, O> {
        match self {
            InvolutiveMapping::Identity { data, underlying } => Ok(InvolutiveMapping::Identity {
                data: data.map_result(f)?,
                underlying,
            }),
            InvolutiveMapping::Source { data, sink_idx } => Ok(InvolutiveMapping::Source {
                data: data.map_result(f)?,
                sink_idx,
            }),
            InvolutiveMapping::Sink { source_idx } => Ok(InvolutiveMapping::Sink { source_idx }),
        }
    }

    pub fn map_empty<E2>(&self) -> InvolutiveMapping<Option<E2>> {
        self.map_data_ref(|_| None)
    }

    pub fn is_internal(&self) -> bool {
        !matches!(self, InvolutiveMapping::Identity { .. })
    }

    pub fn is_source(&self) -> bool {
        matches!(self, InvolutiveMapping::Source { .. })
    }

    pub fn is_sink(&self) -> bool {
        matches!(self, InvolutiveMapping::Sink { .. })
    }

    // pub fn make_source(&mut self, sink_idx: Hedge) -> Option<InvolutiveMapping<E>> {
    //     let data = match self {
    //         InvolutiveMapping::Identity { data, .. } => data.take(),
    //         _ => EdgeData::none(),
    //     };
    //     Some(InvolutiveMapping::Source { data, sink_idx })
    // }

    pub fn new_identity(data: E, orientation: impl Into<Orientation>, underlying: Flow) -> Self {
        let o = orientation.into();
        InvolutiveMapping::Identity {
            data: EdgeData::new(data, o),
            underlying,
        }
    }

    pub fn identity_dot(
        edge_id: Hedge,
        source: NodeIndex,
        attr: Option<&GVEdgeAttrs>,
        orientation: Orientation,
        flow: Flow,
    ) -> String {
        let mut out = "".to_string();

        //we interpret this node as the hidden one. so the flow is reversed.
        match flow {
            Flow::Sink => {
                out.push_str(&format!(
                    "ext{edge_id} [shape=none, label=\"\" flow=source];\n"
                ));
            }
            Flow::Source => {
                out.push_str(&format!(
                    "ext{edge_id} [shape=none, label=\"\" flow=sink];\n"
                ));
            }
        }

        out.push_str(&format!("  ext{edge_id} -> {source}["));
        match (orientation, flow) {
            (Orientation::Default, Flow::Source) => {
                out.push_str("dir=back ");
            }
            (Orientation::Default, Flow::Sink) => {
                out.push_str("dir=forward ");
            }
            (Orientation::Reversed, Flow::Sink) => {
                out.push_str("dir=back ");
            }
            (Orientation::Reversed, Flow::Source) => {
                out.push_str("dir=forward ");
            }
            (Orientation::Undirected, _) => {
                out.push_str("dir=none ");
            }
        }
        if let Some(attr) = attr {
            out.push_str(&format!("{attr}"));
        }
        out.push_str("];\n");
        out
    }

    pub fn pair_dot(
        source: NodeIndex,
        sink: NodeIndex,
        attr: Option<&GVEdgeAttrs>,
        orientation: Orientation,
    ) -> String {
        let mut out = "".to_string();

        out.push_str(&format!("{source} -> {sink}["));
        match orientation {
            Orientation::Default => {
                out.push_str(" dir=forward ");
            }
            Orientation::Reversed => {
                out.push_str(" dir=back ");
            }
            Orientation::Undirected => {
                out.push_str(" dir=none ");
            }
        }
        if let Some(attr) = attr {
            out.push_str(&format!("{attr}];\n"));
        } else {
            out.push_str(" color=\"red:blue;0.5 \" ];\n");
        }
        out
    }
}

#[derive(Clone, Debug, Error)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// Errors that can occur during `Involution` operations.
pub enum InvolutionError {
    #[error("Should have been identity")]
    /// Expected a half-edge to be an identity (unpaired) but it was not.
    NotIdentity,
    #[error("Should have been an paired hedge")]
    /// Expected a half-edge to be part of a pair (source or sink) but it was an identity.
    NotPaired,
    #[error("Involutino is non -involutive for {0}")]
    NonInvolutive(Hedge),
    #[error("Other:{0}")]
    Other(String),
}
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// Manages the topological relationships between half-edges in a graph.
///
/// An involution is a function `f` such that `f(f(x)) = x`. In this context,
/// it maps each half-edge to its unique opposite half-edge. If a half-edge
/// is unpaired (an "identity" or "external" edge), it maps to itself.
///
/// This struct stores a vector of [`InvolutiveMapping<E>`], where each entry
/// corresponds to a half-edge and describes its pairing and associated data.
///
/// # Type Parameters
///
/// - `E`: The type of data associated with each edge. Defaults to `EdgeIndex` if not specified,
///   implying edges primarily store their own index as data, but typically this is
///   instantiated with custom edge data for a `HedgeGraph`.
pub struct Involution<E = EdgeIndex> {
    /// The core data storage: a vector where each element describes the mapping
    /// and data for a specific half-edge. The index in this vector corresponds
    /// to a `Hedge`'s underlying `usize` value.
    pub(super) inv: HedgeVec<InvolutiveMapping<E>>,
}

impl<E> AsRef<Involution<E>> for Involution<E> {
    fn as_ref(&self) -> &Involution<E> {
        self
    }
}

impl<E> Default for Involution<E> {
    fn default() -> Self {
        Involution::new()
    }
}

impl<E> Involution<E> {
    pub fn display(&self) -> String
    where
        E: Display,
    {
        let mut out = "".to_string();
        for (i, e) in self.inv.iter() {
            match e {
                InvolutiveMapping::Identity { underlying, data } => match underlying {
                    Flow::Sink => out.push_str(&format!("{i}<< :{data}\n")),
                    Flow::Source => out.push_str(&format!("{i}>> :{data}\n")),
                },
                InvolutiveMapping::Source { sink_idx, data } => {
                    out.push_str(&format!("{i}->{sink_idx}: {data}\n"));
                }
                InvolutiveMapping::Sink { source_idx } => {
                    out.push_str(&format!("{i}->{source_idx}\n"));
                }
            }
        }
        out
    }

    /// The fundamental function of involutions is to pair edges. This function provides this pairing.
    /// It is a bijective map inv: H -> H, where H is the set of half-edges.
    /// The map is involutive, meaning that inv(inv(h)) = h.
    /// The fixed points of this function correspond to the unpaired hedges, interpreted as dangling or external edges.
    pub fn inv(&self, hedge: Hedge) -> Hedge {
        match self.hedge_data(hedge) {
            InvolutiveMapping::Sink { source_idx } => *source_idx,
            InvolutiveMapping::Source { sink_idx, .. } => *sink_idx,
            _ => hedge,
        }
    }

    // pub fn len(&self) -> usize {
    //     self.inv.len()
    // }

    pub fn is_empty(&self) -> bool {
        self.inv.is_zero_length()
    }

    pub fn new() -> Self {
        Involution {
            inv: HedgeVec::new(),
        }
    }

    /// add a new dangling edge to the involution
    /// returns the index of the new edge
    pub fn add_identity(
        &mut self,
        data: E,
        orientation: impl Into<Orientation>,
        underlying: Flow,
    ) -> Hedge {
        let index = self.inv.len();
        self.inv.push(InvolutiveMapping::new_identity(
            data,
            orientation,
            underlying,
        ));
        index
    }

    /// adds a new hedge that connects to an identity hedge
    pub fn connect_to_identity(&mut self, connect: Hedge) -> Result<Hedge, InvolutionError> {
        let sink = self.inv.len();
        if !self.is_identity(connect) {
            return Err(InvolutionError::NotIdentity);
        }
        self.inv.push(InvolutiveMapping::Sink {
            source_idx: connect,
        });

        let flow = self.flow(connect);
        if let Some(data) = self.get_data_owned(connect) {
            self.set(
                connect,
                InvolutiveMapping::Source {
                    data,
                    sink_idx: sink,
                },
            );
            if let Flow::Sink = flow {
                self.flip_underlying(sink);
            }
            Ok(sink)
        } else {
            panic!("Should have data");
        }
    }

    /// add a connected pair of hedges to the involution
    /// returns the pair of hedges: (source, sink)
    pub fn add_pair(&mut self, data: E, directed: impl Into<Orientation>) -> (Hedge, Hedge) {
        let orientation = directed.into();
        let source = self.add_identity(data, orientation, Flow::Source);
        let sink = self.connect_to_identity(source).unwrap();
        (source, sink)
    }

    pub fn find_from_data(&self, data: &E) -> Option<HedgePair>
    where
        E: PartialEq,
    {
        self.iter_edge_data().find_map(|(i, d)| {
            if &d.data == data {
                Some(self.hedge_pair(i))
            } else {
                None
            }
        })
    }

    pub fn last(&self) -> Option<Hedge> {
        self.len().0.checked_sub(1).map(Hedge)
    }

    fn _validate(&self) -> bool {
        self.iter_idx().all(|h| {
            let pair = self.inv(h);
            if h == pair {
                self.is_identity(h)
            } else {
                self.inv(pair) == h
                    && ((self.is_sink(h) && self.is_source(pair))
                        || (self.is_source(h) && self.is_sink(pair)))
            }
        })
    }

    /// modify the involution to set hedge to sink.
    /// this might make the structure invalid
    fn set(&mut self, hedge: Hedge, mut mapping: InvolutiveMapping<E>) -> Option<EdgeData<E>> {
        std::mem::swap(&mut self.inv[hedge], &mut mapping);

        mapping.data()
    }

    /// extract the data from the hedge replacing it with a dummy mapping
    fn get_data_owned(&mut self, hedge: Hedge) -> Option<EdgeData<E>> {
        self.set(hedge, InvolutiveMapping::dummy())
    }

    /// connect a pair of identity hedges,
    /// taking orientation and data from source_idx
    /// errors if the hedges are not identity
    pub fn connect_identities(
        &mut self,
        source_idx: Hedge,
        sink_idx: Hedge,
        merge_fn: impl FnOnce(Flow, EdgeData<E>, Flow, EdgeData<E>) -> (Flow, EdgeData<E>),
    ) -> Result<HedgePair, InvolutionError>
    where
        E: Clone,
    {
        if self.is_identity(source_idx) && self.is_identity(sink_idx) {
            let source_flow = self.flow(source_idx);
            let sink_flow = self.flow(sink_idx);
            let source_data = self.get_data_owned(source_idx);
            let sink_data = self.get_data_owned(sink_idx);
            if let (Some(so_d), Some(si_d)) = (source_data, sink_data) {
                let (new_flow, new_data) = merge_fn(source_flow, so_d, sink_flow, si_d);

                match new_flow {
                    Flow::Source => {
                        self.set(
                            source_idx,
                            InvolutiveMapping::Source {
                                data: new_data,
                                sink_idx,
                            },
                        );
                        self.set(sink_idx, InvolutiveMapping::Sink { source_idx });
                        Ok(HedgePair::Paired {
                            source: source_idx,
                            sink: sink_idx,
                        })
                    }
                    Flow::Sink => {
                        self.set(
                            source_idx,
                            InvolutiveMapping::Sink {
                                source_idx: sink_idx,
                            },
                        );
                        self.set(
                            sink_idx,
                            InvolutiveMapping::Source {
                                data: new_data,
                                sink_idx: source_idx,
                            },
                        );
                        Ok(HedgePair::Paired {
                            source: sink_idx,
                            sink: source_idx,
                        })
                    }
                }
            } else {
                Err(InvolutionError::Other(format!("Data is None")))
            }
        } else {
            Err(InvolutionError::NotIdentity)
        }
    }

    /// Splits the edge that hedge is a part of into two dangling hedges, adding the data to the side given by hedge.
    /// The underlying orientation of the new edges is the same as the original edge, i.e. the source will now have `Flow::Source` and the sink will have `Flow::Sink`.
    /// The superficial orientation has to be given knowing this.
    pub fn split_edge(&mut self, hedge: Hedge, data: EdgeData<E>) -> Result<(), InvolutionError> {
        if self.is_identity(hedge) {
            Err(InvolutionError::NotPaired)
        } else {
            let invh = self.inv(hedge);
            if let Some(replacing_data) = self.get_data_owned(hedge) {
                //hedge is a source, replace its data with the new data
                self.set(
                    hedge,
                    InvolutiveMapping::Identity {
                        data,
                        underlying: Flow::Source,
                    },
                );
                self.set(
                    invh,
                    InvolutiveMapping::Identity {
                        data: replacing_data,
                        underlying: Flow::Sink,
                    },
                );
            } else {
                // hedge is a sink, give it the new data
                self.set(
                    hedge,
                    InvolutiveMapping::Identity {
                        data,
                        underlying: Flow::Sink,
                    },
                );
                let data = self.get_data_owned(invh).unwrap(); //extract the data from the source
                self.set(
                    invh,
                    InvolutiveMapping::Identity {
                        data,
                        underlying: Flow::Source,
                    },
                );
            }
            Ok(())
        }
    }

    fn set_inv(&mut self, a: Hedge, inva: Hedge) -> Option<()> {
        match &mut self.inv[a] {
            InvolutiveMapping::Source { sink_idx, .. } => {
                *sink_idx = inva;
                Some(())
            }
            InvolutiveMapping::Sink { source_idx } => {
                *source_idx = inva;
                Some(())
            }
            _ => None,
        }
    }

    pub fn check_involutivity(&self) -> Result<(), String> {
        for (h, _) in self.iter() {
            if h != self.inv(self.inv(h)) {
                return Err(format!("Involution is not involutive at half-edge {h}"));
            }
        }
        Ok(())
    }
    // fn put_at_end(&mut self,)

    /// Extracts the edges included in the subgraph into a separate involution. If the involution maps edges across the splitting boundary, i.e. if i in graph but inv(i) not in graph, then the edge involution is modified on each side, to turn this into a dangling edge. The data assigned to a fully owned edge is mapped by the internal_data closure, whilst those that get split, are assigned the data gotten from the split_edge_fn closure. The new dangling edges retain the old underlying orientation.
    pub fn extract<S: SubSetLike, O>(
        &mut self,
        graph: &S,
        mut split_edge_fn: impl FnMut(EdgeData<&E>) -> EdgeData<O>,
        mut internal_data: impl FnMut(EdgeData<E>) -> EdgeData<O>,
    ) -> Involution<O> {
        let (extract, left) = self.delete_impl(graph);

        let extracted_inv = extract
            .into_iter()
            .map(|(_, m)| match m {
                InvolutiveMapping::Identity { data, underlying } => InvolutiveMapping::Identity {
                    data: internal_data(data),
                    underlying,
                },
                InvolutiveMapping::Sink { source_idx } => {
                    if source_idx >= left {
                        // left includes the hedge extracted (it is the length of the left hedges)
                        InvolutiveMapping::Sink {
                            source_idx: Hedge(source_idx.0 - left.0),
                        }
                    } else {
                        let underlying = -self.flow(source_idx);
                        InvolutiveMapping::Identity {
                            data: split_edge_fn(self.edge_data(source_idx).as_ref()),
                            underlying,
                        }
                    }
                }
                InvolutiveMapping::Source { data, sink_idx } => InvolutiveMapping::Source {
                    data: internal_data(data),
                    sink_idx: Hedge(sink_idx.0 - left.0),
                },
            })
            .collect();

        Involution { inv: extracted_inv }
    }

    pub fn delete<S: SubSetLike>(&mut self, graph: &S) {
        self.delete_impl(graph);
    }

    fn delete_impl<S: SubSetLike>(&mut self, graph: &S) -> (HedgeVec<InvolutiveMapping<E>>, Hedge) {
        let mut left = Hedge(0);
        let mut extracted = self.inv.len();
        while left < extracted {
            if !graph.includes(&left) {
                //left is in the right place
                left.0 += 1;
            } else {
                //left needs to be swapped
                extracted.0 -= 1;
                if !graph.includes(&extracted) {
                    //only with an extracted that is in the wrong spot
                    self.swap(left, extracted);
                    left.0 += 1;
                }
            }
        }

        for i in 0..(left.0) {
            if self.inv(Hedge(i)) >= left {
                let flow = if self.set_as_source(Hedge(i)) {
                    Flow::Sink
                } else {
                    Flow::Source
                };

                self.source_to_identity_impl(Hedge(i), flow, flow == Flow::Sink);
            }
        }
        (self.inv.split_off(left), left)
    }

    pub fn hedge_pair(&self, hedge: Hedge) -> HedgePair {
        HedgePair::from_half_edge(hedge, self)
    }

    pub(crate) fn smart_data<S: SubSetLike>(
        &self,
        hedge: Hedge,
        subgraph: &S,
    ) -> Option<&EdgeData<E>> {
        if subgraph.includes(&hedge) {
            match self.hedge_data(hedge) {
                InvolutiveMapping::Identity { .. } => Some(self.edge_data(hedge)),
                InvolutiveMapping::Source { .. } => Some(self.edge_data(hedge)),
                InvolutiveMapping::Sink { source_idx } => {
                    if subgraph.includes(source_idx) {
                        None
                    } else {
                        Some(self.edge_data(*source_idx))
                    }
                }
            }
        } else {
            None
        }
    }

    fn _smart_data_mut<S: SubSetLike>(
        &mut self,
        hedge: Hedge,
        subgraph: &S,
    ) -> Option<&mut EdgeData<E>> {
        if subgraph.includes(&hedge) {
            match self.hedge_data(hedge) {
                InvolutiveMapping::Identity { .. } => Some(self.edge_data_mut(hedge)),
                InvolutiveMapping::Source { .. } => Some(self.edge_data_mut(hedge)),
                InvolutiveMapping::Sink { source_idx } => {
                    if subgraph.includes(source_idx) {
                        None
                    } else {
                        Some(self.edge_data_mut(*source_idx))
                    }
                }
            }
        } else {
            None
        }
    }

    pub fn first_internal<S: SubSetLike>(&self, subgraph: &S) -> Option<Hedge> {
        self.iter_idx().find(|e| self.is_internal(*e, subgraph))
    }

    pub fn n_internals<S: SubSetLike>(&self, subgraph: &S) -> usize {
        subgraph
            .included_iter()
            .filter(|i| self.is_internal(*i, subgraph))
            .count()
            / 2
    }

    /// check if the edge is internal and totally included in the subgraph
    pub fn is_internal<S: SubSetLike>(&self, index: Hedge, subgraph: &S) -> bool {
        if !subgraph.includes(&index) {
            return false;
        }
        match &self.inv[index] {
            InvolutiveMapping::Identity { .. } => false,
            InvolutiveMapping::Source { sink_idx, .. } => subgraph.includes(sink_idx),
            InvolutiveMapping::Sink { source_idx } => subgraph.includes(source_idx),
        }
    }

    pub fn orientation(&self, hedge: Hedge) -> Orientation {
        match self.hedge_data(hedge) {
            InvolutiveMapping::Identity { data, .. } => data.orientation,
            InvolutiveMapping::Source { data, .. } => data.orientation,
            InvolutiveMapping::Sink { source_idx } => self.orientation(*source_idx),
        }
    }

    pub fn set_orientation(&mut self, hedge: Hedge, orientation: Orientation) {
        self.edge_data_mut(hedge).orientation = orientation;
    }

    pub fn flow(&self, hedge: Hedge) -> Flow {
        self.hedge_data(hedge).flow()
    }

    pub fn iter_edge_data(&self) -> impl Iterator<Item = (Hedge, &EdgeData<E>)> {
        self.iter().filter_map(|(e, m)| match m {
            InvolutiveMapping::Source { data, .. } => Some((e, data)),
            InvolutiveMapping::Identity { data, .. } => Some((e, data)),
            _ => None,
        })
    }

    pub fn iter_edge_data_mut(&mut self) -> impl Iterator<Item = &mut EdgeData<E>> {
        self.iter_mut().filter_map(|(_, m)| match m {
            InvolutiveMapping::Source { data, .. } => Some(data),
            InvolutiveMapping::Identity { data, .. } => Some(data),
            _ => None,
        })
    }

    pub fn iter(&self) -> InvolutionIter<'_, E> {
        self.into_iter()
    }

    pub fn iter_mut(&mut self) -> InvolutionIterMut<'_, E> {
        self.into_iter()
    }

    pub fn iter_idx(&self) -> impl Iterator<Item = Hedge> {
        (0..self.inv.len().0).map(Hedge)
    }

    pub fn as_ref(&self) -> Involution<&E> {
        let inv = self.inv.iter().map(|(_, e)| e.as_ref()).collect();
        Involution { inv }
    }

    pub fn map_data_ref<'a, G, E2>(&'a self, g: &G) -> Involution<E2>
    where
        G: FnMut(&'a E) -> E2 + Clone,
    {
        let inv = self
            .inv
            .iter()
            .map(|(h, e)| (h, e.map_data_ref(g.clone())))
            .collect();

        Involution { inv }
    }

    pub fn map_data<G, E2>(self, g: &G) -> Involution<E2>
    where
        G: FnMut(E) -> E2 + Clone,
    {
        let inv = self
            .inv
            .into_iter()
            .map(|(h, e)| (h, e.map_data(g.clone())))
            .collect();

        Involution { inv }
    }

    pub fn map_data_option<F, G, E2>(self, g: &G) -> Option<Involution<E2>>
    where
        G: FnMut(E) -> Option<E2> + Clone,
    {
        let inv = self
            .inv
            .into_iter()
            .map(|(h, e)| e.map_data_option(g.clone()).map(|e| (h, e)))
            .collect::<Option<HedgeVec<_>>>()?;

        Some(Involution { inv })
    }

    pub fn map_data_result<F, G, E2, O>(self, g: &G) -> Result<Involution<E2>, O>
    where
        G: FnMut(E) -> Result<E2, O> + Clone,
    {
        let inv = self
            .inv
            .into_iter()
            .map(|(h, e)| e.map_data_result(g.clone()).map(|e| (h, e)))
            .collect::<Result<HedgeVec<_>, O>>()?;

        Ok(Involution { inv })
    }

    pub fn map_full<E2>(
        self,
        mut g: impl FnMut(HedgePair, EdgeData<E>) -> EdgeData<E2>,
    ) -> Involution<E2> {
        let inv = self
            .into_iter()
            .map(|(i, e)| match e {
                InvolutiveMapping::Identity { data, underlying } => InvolutiveMapping::Identity {
                    data: g(
                        HedgePair::Unpaired {
                            hedge: i,
                            flow: underlying,
                        },
                        data,
                    ),
                    underlying,
                },
                InvolutiveMapping::Source { data, sink_idx } => InvolutiveMapping::Source {
                    data: g(
                        HedgePair::Paired {
                            source: i,
                            sink: sink_idx,
                        },
                        data,
                    ),

                    sink_idx,
                },
                InvolutiveMapping::Sink { source_idx } => InvolutiveMapping::Sink { source_idx },
            })
            .collect();

        Involution { inv }
    }

    pub fn is_sink(&self, hedge: Hedge) -> bool {
        self.hedge_data(hedge).is_sink()
    }

    pub fn is_source(&self, hedge: Hedge) -> bool {
        self.hedge_data(hedge).is_source()
    }

    pub fn is_identity(&self, hedge: Hedge) -> bool {
        self.hedge_data(hedge).is_identity()
    }

    /// If the data at `hedge` was a sink, turn it into a source.
    /// Otherwise do nothing.
    pub fn set_as_source(&mut self, hedge: Hedge) -> bool {
        let is_sink = self.is_sink(hedge);
        if is_sink {
            self.flip_underlying(hedge);
            return true;
        }
        false
    }

    // fn source_to_identity(&mut self, hedge: Hedge, underlying: Flow) {
    //     self.source_to_identity_impl(hedge, underlying, false);
    // }

    pub(crate) fn source_to_identity_impl(
        &mut self,
        hedge: Hedge,
        underlying: Flow,
        reverse_orientation: bool,
    ) {
        if self.is_source(hedge) {
            let mut data = self.get_data_owned(hedge).unwrap();
            if reverse_orientation {
                data.orientation = data.orientation.reverse()
            }
            self.inv[hedge] = InvolutiveMapping::Identity { data, underlying };
        }
    }

    /// If the data at `hedge` was a source, turn it into a sink.
    /// Otherwise do nothing.
    pub fn set_as_sink(&mut self, hedge: Hedge) {
        let is_source = self.is_source(hedge);
        if is_source {
            self.flip_underlying(hedge)
        }
    }

    /// Swap the data carrier in a pair of hedges.
    pub(super) fn flip_underlying(&mut self, hedge: Hedge) {
        let pair = self.inv(hedge);

        self.inv.swap(hedge, pair);

        match self.hedge_data_mut(hedge) {
            InvolutiveMapping::Identity { underlying, data } => {
                *underlying = -*underlying;
                data.orientation = data.orientation.reverse();
            }
            InvolutiveMapping::Source { data, sink_idx } => {
                *sink_idx = pair;
                data.orientation = data.orientation.reverse();
            }
            InvolutiveMapping::Sink { source_idx } => {
                *source_idx = pair;
            }
        }

        match self.hedge_data_mut(pair) {
            InvolutiveMapping::Source { data, sink_idx } => {
                *sink_idx = hedge;
                data.orientation = data.orientation.reverse();
            }
            InvolutiveMapping::Sink { source_idx } => {
                *source_idx = hedge;
            }
            _ => {}
        }
    }

    pub(crate) fn random(len: usize, seed: u64) -> Involution<()> {
        let mut rng = SmallRng::seed_from_u64(seed);

        let mut inv = Involution::new();

        for _ in 0..len {
            let r = rng.gen_bool(0.1);
            if r {
                inv.add_identity((), Orientation::Undirected, Flow::Sink);
            } else {
                inv.add_pair((), false);
            }
        }

        inv
    }

    pub fn edge_data(&self, index: Hedge) -> &EdgeData<E> {
        match &self.inv[index] {
            InvolutiveMapping::Source { data, .. } => data,
            InvolutiveMapping::Identity { data, .. } => data,
            InvolutiveMapping::Sink { source_idx } => self.edge_data(*source_idx),
        }
    }

    #[allow(clippy::needless_return)]
    pub fn edge_data_mut(&mut self, index: Hedge) -> &mut EdgeData<E> {
        if let InvolutiveMapping::Sink { source_idx } = self.inv[index] {
            return self.edge_data_mut(source_idx);
        }
        match &mut self.inv[index] {
            InvolutiveMapping::Source { data, .. } => return data,
            InvolutiveMapping::Identity { data, .. } => return data,
            _ => unreachable!(),
        };
    }

    pub(super) fn hedge_data(&self, hedge: Hedge) -> &InvolutiveMapping<E> {
        &self.inv[hedge]
    }

    pub(crate) fn hedge_data_mut(&mut self, hedge: Hedge) -> &mut InvolutiveMapping<E> {
        &mut self.inv[hedge]
    }

    pub(crate) fn data_inv(&self, hedge: Hedge) -> Hedge {
        match self.hedge_data(hedge) {
            InvolutiveMapping::Sink { source_idx } => *source_idx,
            _ => hedge,
        }
    }

    pub fn print<S: SubSetLike>(
        &self,
        subgraph: &S,
        h_label: &impl Fn(&E) -> Option<String>,
    ) -> String {
        let mut out = "".to_string();
        for (i, e) in self.iter() {
            if !subgraph.includes(&i) {
                continue;
            }
            match e {
                InvolutiveMapping::Identity { .. } => {
                    out.push_str(&format!("{i}\n"));
                }
                InvolutiveMapping::Source { data, sink_idx } => {
                    let d = &data.data;
                    if let Some(l) = h_label(d) {
                        out.push_str(&format!("{i}->{sink_idx}: {l}\n"));
                    } else {
                        out.push_str(&format!("{i}->{sink_idx}\n"));
                    }
                }
                InvolutiveMapping::Sink { source_idx } => {
                    out.push_str(&format!("{i}->{source_idx}\n"));
                }
            }
        }
        out
    }
}

impl<E> FromIterator<InvolutiveMapping<E>> for Involution<E> {
    fn from_iter<I: IntoIterator<Item = InvolutiveMapping<E>>>(iter: I) -> Self {
        Involution {
            inv: iter.into_iter().collect(),
        }
    }
}

impl<E> Index<Hedge> for Involution<E> {
    type Output = E;

    fn index(&self, index: Hedge) -> &Self::Output {
        match &self.inv[index] {
            InvolutiveMapping::Identity { data, .. } => &data.data,
            InvolutiveMapping::Source { data, .. } => &data.data,
            InvolutiveMapping::Sink { source_idx } => &self[*source_idx],
        }
    }
}

impl<E> IndexMut<Hedge> for Involution<E> {
    fn index_mut(&mut self, index: Hedge) -> &mut Self::Output {
        let invh = self.data_inv(index);

        match self.hedge_data_mut(invh) {
            InvolutiveMapping::Identity { data, .. } => &mut data.data,
            InvolutiveMapping::Source { data, .. } => &mut data.data,
            _ => panic!("should have gotten data inv"),
        }
    }
}

impl<E> Index<&HedgePair> for Involution<E> {
    type Output = EdgeData<E>;

    fn index(&self, index: &HedgePair) -> &Self::Output {
        match index {
            HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                if let InvolutiveMapping::Source { data, sink_idx } = &self.inv[*source] {
                    debug_assert_eq!(sink_idx, sink);
                    data
                } else {
                    panic!("{index} should be aligned with involution but here {source} is not a source")
                }
            }
            HedgePair::Unpaired { hedge, flow } => {
                if let InvolutiveMapping::Identity { data, underlying } = &self.inv[*hedge] {
                    debug_assert_eq!(flow, underlying);
                    data
                } else {
                    panic!("{index} should be aligned with involution but here {hedge} is not a identity")
                }
            }
        }
    }
}

// #[derive(Clone, Debug, Serialize, Deserialize)]
// pub struct EdgeVec<E> {
//     inv: Involution<usize>,
//     data: Vec<E>,
// }

pub struct DrainingInvolutionIter<E> {
    into: HedgeVecIntoIter<InvolutiveMapping<E>>,
}

impl<E> Iterator for DrainingInvolutionIter<E> {
    type Item = (Hedge, InvolutiveMapping<E>);
    fn next(&mut self) -> Option<Self::Item> {
        self.into.next()
    }
}

impl<E> IntoIterator for Involution<E> {
    type Item = (Hedge, InvolutiveMapping<E>);
    type IntoIter = DrainingInvolutionIter<E>;
    fn into_iter(self) -> Self::IntoIter {
        DrainingInvolutionIter {
            into: self.inv.into_iter(),
        }
    }
}

pub struct InvolutionIter<'a, E> {
    into: HedgeVecIter<'a, InvolutiveMapping<E>>,
}

impl<'a, E> Iterator for InvolutionIter<'a, E> {
    type Item = (Hedge, &'a InvolutiveMapping<E>);
    fn next(&mut self) -> Option<Self::Item> {
        self.into.next()
    }
}

impl<'a, E> IntoIterator for &'a Involution<E> {
    type Item = (Hedge, &'a InvolutiveMapping<E>);
    type IntoIter = InvolutionIter<'a, E>;
    fn into_iter(self) -> Self::IntoIter {
        InvolutionIter {
            into: self.inv.iter(),
        }
    }
}

pub struct InvolutionIterMut<'a, E> {
    into: HedgeVecIterMut<'a, InvolutiveMapping<E>>,
}

impl<'a, E> Iterator for InvolutionIterMut<'a, E> {
    type Item = (Hedge, &'a mut InvolutiveMapping<E>);
    fn next(&mut self) -> Option<Self::Item> {
        self.into.next()
    }
}

impl<'a, E> IntoIterator for &'a mut Involution<E> {
    type Item = (Hedge, &'a mut InvolutiveMapping<E>);
    type IntoIter = InvolutionIterMut<'a, E>;
    fn into_iter(self) -> Self::IntoIter {
        InvolutionIterMut {
            into: self.inv.iter_mut(),
        }
    }
}

// pub trait Get<H> {
//     type Output<'a>
//     where
//         Self: 'a;

//     type MutOutput<'a>
//     where
//         Self: 'a;
//     fn get(&self, h: H) -> Self::Output<'_>;

//     fn get_mut(&mut self, h: H) -> Self::MutOutput<'_>;
// }

// impl<N, E> Get<Hedge> for Involution<N, E> {
//     type Output<'a>
//         = &'a InvolutiveMapping<E>
//     where
//         Self: 'a;
//     type MutOutput<'a>
//         = &'a mut InvolutiveMapping<E>
//     where
//         Self: 'a;
//     fn get(&self, h: Hedge) -> Self::Output<'_> {
//         &self.inv[h.0]
//     }

//     fn get_mut(&mut self, h: Hedge) -> Self::MutOutput<'_> {
//         &mut self.inv[h.0]
//     }
// }

// impl<N, E> Get<HedgePair> for Involution<N, E> {
//     type Output<'a>
//         = Option<(Flow, &'a EdgeData<E>)>
//     where
//         Self: 'a;

//     type MutOutput<'a>
//         = Option<(Flow, &'a mut EdgeData<E>)>
//     where
//         Self: 'a;

//     fn get(&self, h: HedgePair) -> Self::Output<'_> {
//         let (source, flow) = match h {
//             HedgePair::Paired { source, .. } => (source, Flow::Source),
//             HedgePair::Unpaired { hedge, flow } => (hedge, flow),
//             HedgePair::Split { source, split, .. } => (source, split),
//         };

//         match &self[source] {
//             InvolutiveMapping::Source { data, .. } => Some((flow, data)),
//             InvolutiveMapping::Identity { data, underlying } => Some((*underlying, data)),
//             _ => None,
//         }
//     }

//     fn get_mut(&mut self, h: HedgePair) -> Self::MutOutput<'_> {
//         let (source, flow) = match h {
//             HedgePair::Paired { source, .. } => (source, Flow::Source),
//             HedgePair::Unpaired { hedge, flow } => (hedge, flow),
//             HedgePair::Split { source, split, .. } => (source, split),
//         };

//         match &mut self[source] {
//             InvolutiveMapping::Source { data, .. } => Some((flow, data)),
//             InvolutiveMapping::Identity { data, underlying } => Some((*underlying, data)),
//             _ => None,
//         }
//     }
// }

impl<E> Display for Involution<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut out = "".to_string();
        for (i, e) in self.inv.iter() {
            match e {
                InvolutiveMapping::Identity { underlying, .. } => match underlying {
                    Flow::Sink => out.push_str(&format!("{i}<<\n")),
                    Flow::Source => out.push_str(&format!("{i}>>\n")),
                },
                InvolutiveMapping::Source { sink_idx, .. } => {
                    out.push_str(&format!("{i}->{sink_idx}\n"));
                }
                InvolutiveMapping::Sink { source_idx } => {
                    out.push_str(&format!("{i}<-{source_idx}\n"));
                }
            }
        }
        write!(f, "{out}")
    }
}

impl<E> Swap<Hedge> for Involution<E> {
    fn is_zero_length(&self) -> bool {
        self.inv.is_zero_length()
    }

    // type Item = InvolutiveMapping<E>;

    // fn filter(&self, id: &Hedge, filter: &impl Fn(&Hedge, &Self::Item) -> bool) -> bool {
    //     filter(id, &self.inv[*id])
    // }

    ///Invalidates Subgraphs
    fn swap(&mut self, a: Hedge, b: Hedge) {
        // Determine if the pointer to them needs to be modified
        // Need to do this before set inv otherwise you get into an invalid state
        let inva = self.inv(a);
        let invb = self.inv(b);

        if inva != a {
            self.set_inv(inva, b).unwrap();
        }
        if invb != b {
            self.set_inv(invb, a).unwrap();
        }
        self.inv.swap(a, b);
    }

    fn len(&self) -> Hedge {
        self.inv.len()
    }
}

#[cfg(test)]
pub mod test;
