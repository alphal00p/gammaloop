use std::path::Path;

#[cfg(feature = "serde")]
use figment;
use figment::Figment;

use crate::half_edge::{
    nodestore::{NodeForest, Topology},
    HedgeGraph,
};

use super::{
    syntax::FlatGraph, DotEdgeData, DotGraph, DotHedgeData, DotVertexData, GlobalData,
    HedgeParseError,
};

#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct GraphSet<E, V, H, G, S: Topology> {
    pub global_data: Vec<G>,
    pub set: Vec<HedgeGraph<E, V, H, S>>,
}

/// A collection of DOT graphs parsed from a file or string containing multiple graphs.
pub type DotGraphSet = GraphSet<DotEdgeData, DotVertexData, DotHedgeData, GlobalData, NodeForest>;

impl<S: Topology> GraphSet<DotEdgeData, DotVertexData, DotHedgeData, GlobalData, S> {
    /// Parses a file holding one or more DOT graphs.
    pub fn from_file<P>(p: P) -> Result<Self, HedgeParseError<(), (), (), ()>>
    where
        P: AsRef<Path>,
    {
        Self::from_flat(
            FlatGraph::parse_all(&std::fs::read_to_string(p)?)?,
            Figment::new(),
        )
    }

    /// Parses input holding one or more DOT graphs.
    pub fn from_string<Str: AsRef<str>>(s: Str) -> Result<Self, HedgeParseError<(), (), (), ()>> {
        Self::from_flat(FlatGraph::parse_all(s.as_ref())?, Figment::new())
    }

    #[cfg(feature = "serde")]
    pub fn from_string_with_figment<Str: AsRef<str>>(
        s: Str,
        figment: figment::Figment,
    ) -> Result<Self, HedgeParseError<(), (), (), ()>> {
        Self::from_flat(FlatGraph::parse_all(s.as_ref())?, figment)
    }

    fn from_flat(
        graphs: Vec<FlatGraph<'_>>,
        figment: Figment,
    ) -> Result<Self, HedgeParseError<(), (), (), ()>> {
        let mut set = Vec::with_capacity(graphs.len());
        let mut global_data = Vec::with_capacity(graphs.len());
        for graph in graphs {
            let graph = DotGraph::<S>::from_flat(graph, figment.clone())?;
            set.push(graph.graph);
            global_data.push(graph.global_data);
        }
        Ok(GraphSet { set, global_data })
    }
}

impl<E, V, H, G, S: Topology> GraphSet<E, V, H, G, S> {
    #[cfg(feature = "rkyv")]
    pub fn to_rkyv_bytes<const N: usize>(&self) -> Result<rkyv::AlignedVec, String>
    where
        Self: rkyv::Serialize<rkyv::ser::serializers::AllocSerializer<N>>,
    {
        rkyv::to_bytes::<_, N>(self).map_err(|err| err.to_string())
    }

    #[cfg(feature = "rkyv")]
    /// Returns the archived graph set root without validating the byte buffer.
    ///
    /// # Safety
    ///
    /// `bytes` must contain a valid rkyv archive produced for this exact
    /// `GraphSet` type, and the returned reference must not outlive `bytes`.
    pub unsafe fn archived_from_bytes(bytes: &[u8]) -> &<Self as rkyv::Archive>::Archived
    where
        Self: rkyv::Archive,
    {
        unsafe { rkyv::archived_root::<Self>(bytes) }
    }
}

impl<S: Topology> IntoIterator
    for GraphSet<DotEdgeData, DotVertexData, DotHedgeData, GlobalData, S>
{
    type Item = DotGraph<S>;
    type IntoIter = std::iter::Map<
        std::iter::Zip<
            std::vec::IntoIter<HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData, S>>,
            std::vec::IntoIter<GlobalData>,
        >,
        fn(
            (
                HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData, S>,
                GlobalData,
            ),
        ) -> DotGraph<S>,
    >;

    fn into_iter(self) -> Self::IntoIter {
        self.set
            .into_iter()
            .zip(self.global_data)
            .map(|(graph, global_data)| DotGraph { global_data, graph })
    }
}
