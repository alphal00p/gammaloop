use rkyv::Archived;

use crate::half_edge::{
    involution::{ArchivedFlow, EdgeIndex, Hedge, Orientation},
    nodestore::DefaultNodeStore,
    subgraph::{ModifySubSet, SuBitGraph, SubSetLike},
    HedgeGraph,
};

use super::{
    hedge::compass_pt_to_u8, set::DotGraphSet, DotEdgeData, DotGraph, DotHedgeData, DotVertexData,
    GlobalData,
};

type ArchivedDotGraph = Archived<DotGraph>;
type ArchivedDotGraphBytesSetRoot = Archived<DotGraphBytesSet>;
type ArchivedDotHedgeGraph =
    Archived<HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData, DefaultNodeStore<DotVertexData>>>;

#[derive(Debug, Clone, PartialEq, Eq, rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)]
pub struct DotGraphBytesSet {
    pub graphs: Vec<Vec<u8>>,
}

#[derive(Clone, Copy)]
pub struct ArchivedDotGraphBytesSetView<'a> {
    archived: &'a ArchivedDotGraphBytesSetRoot,
}

#[derive(Clone, Copy)]
pub struct ArchivedDotGraphView<'a> {
    bytes: &'a [u8],
    global_data: &'a Archived<GlobalData>,
    graph: &'a ArchivedDotHedgeGraph,
}

#[derive(Clone, Copy)]
pub struct ArchivedDotVertexView<'a> {
    pub node: crate::half_edge::NodeIndex,
    pub data: &'a Archived<DotVertexData>,
}

#[derive(Clone, Copy)]
pub struct ArchivedDotEdgeView<'a> {
    pub edge: EdgeIndex,
    pub pair: &'a Archived<crate::half_edge::involution::HedgePair>,
    pub data: &'a Archived<DotEdgeData>,
    pub orientation: &'a Archived<Orientation>,
}

#[derive(Clone, Copy)]
pub struct ArchivedDotEndpointView<'a> {
    pub node: crate::half_edge::NodeIndex,
    pub hedge: Hedge,
    pub data: &'a Archived<DotHedgeData>,
}

#[derive(Clone, Copy)]
pub struct ArchivedDotEdgeEndpointsView<'a> {
    pub source: Option<ArchivedDotEndpointView<'a>>,
    pub sink: Option<ArchivedDotEndpointView<'a>>,
}

impl DotGraphSet {
    pub fn into_graph_bytes_set<const N: usize>(self) -> Result<DotGraphBytesSet, String> {
        let mut graphs = Vec::with_capacity(self.set.len());
        for graph in self {
            graphs.push(graph.to_rkyv_bytes::<N>()?.to_vec());
        }
        Ok(DotGraphBytesSet { graphs })
    }

    pub fn to_graph_bytes_set<const N: usize>(&self) -> Result<DotGraphBytesSet, String> {
        self.clone().into_graph_bytes_set::<N>()
    }
}

impl DotGraphBytesSet {
    pub fn archived_view<'a>(bytes: &'a [u8]) -> ArchivedDotGraphBytesSetView<'a> {
        ArchivedDotGraphBytesSetView::from_bytes(bytes)
    }

    pub fn to_rkyv_bytes<const N: usize>(&self) -> Result<rkyv::AlignedVec, String>
    where
        Self: rkyv::Serialize<rkyv::ser::serializers::AllocSerializer<N>>,
    {
        rkyv::to_bytes::<_, N>(self).map_err(|err| err.to_string())
    }

    pub unsafe fn archived_from_bytes<'a>(bytes: &'a [u8]) -> &'a ArchivedDotGraphBytesSetRoot {
        unsafe { rkyv::archived_root::<Self>(bytes) }
    }
}

impl DotGraph {
    pub fn archived_view<'a>(bytes: &'a [u8]) -> ArchivedDotGraphView<'a> {
        ArchivedDotGraphView::from_bytes(bytes)
    }
}

impl<'a> ArchivedDotGraphBytesSetView<'a> {
    pub fn from_bytes(bytes: &'a [u8]) -> Self {
        Self {
            archived: unsafe { DotGraphBytesSet::archived_from_bytes(bytes) },
        }
    }

    pub fn len(&self) -> usize {
        self.archived.graphs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn get_bytes(&self, index: usize) -> Option<&'a [u8]> {
        self.archived
            .graphs
            .as_slice()
            .get(index)
            .map(|bytes| bytes.as_slice())
    }

    pub fn get(&self, index: usize) -> Option<ArchivedDotGraphView<'a>> {
        self.get_bytes(index).map(ArchivedDotGraphView::from_bytes)
    }

    pub fn iter_bytes(&self) -> impl Iterator<Item = &'a [u8]> + 'a {
        self.archived
            .graphs
            .as_slice()
            .iter()
            .map(|bytes| bytes.as_slice())
    }
}

impl<'a> ArchivedDotGraphView<'a> {
    pub(crate) fn from_bytes(bytes: &'a [u8]) -> Self {
        let archived: &'a ArchivedDotGraph = unsafe { DotGraph::archived_from_bytes(bytes) };
        Self {
            bytes,
            global_data: &archived.global_data,
            graph: &archived.graph,
        }
    }

    pub fn bytes(&self) -> &'a [u8] {
        self.bytes
    }

    pub fn global_data(&self) -> &'a Archived<GlobalData> {
        self.global_data
    }

    pub fn n_edges(&self) -> usize {
        self.graph.edge_store.data.0.len()
    }

    pub fn n_hedges(&self) -> usize {
        self.graph.hedge_data.0.len()
    }

    pub fn n_vertices(&self) -> usize {
        self.graph.node_store.node_data.0.len()
    }

    pub fn vertex_data(&self) -> impl Iterator<Item = ArchivedDotVertexView<'_>> + '_ {
        self.graph
            .node_store
            .node_data
            .0
            .as_slice()
            .iter()
            .enumerate()
            .map(|(node, data)| ArchivedDotVertexView {
                node: crate::half_edge::NodeIndex(node),
                data,
            })
    }

    pub fn vertex_data_of<'b>(
        &'b self,
        subgraph: &'b SuBitGraph,
    ) -> impl Iterator<Item = ArchivedDotVertexView<'b>> + 'b {
        self.graph
            .node_store
            .nodes
            .0
            .as_slice()
            .iter()
            .enumerate()
            .filter(move |(_, neighbors)| neighbors.intersects_owned(subgraph))
            .map(move |(node, _)| ArchivedDotVertexView {
                node: crate::half_edge::NodeIndex(node),
                data: &self.graph.node_store.node_data.0.as_slice()[node],
            })
    }

    pub fn edge_data(&self) -> impl Iterator<Item = ArchivedDotEdgeView<'_>> + '_ {
        self.graph
            .edge_store
            .data
            .0
            .as_slice()
            .iter()
            .enumerate()
            .map(move |(edge, (data, pair))| ArchivedDotEdgeView {
                edge: EdgeIndex(edge),
                pair,
                data,
                orientation: self.orientation_of_pair(pair),
            })
    }

    pub fn edge_data_of<'b>(
        &'b self,
        subgraph: &'b SuBitGraph,
    ) -> impl Iterator<Item = ArchivedDotEdgeView<'b>> + 'b {
        self.edge_data()
            .filter(move |edge| edge.pair.intersects(subgraph))
    }

    pub fn endpoints_of_edge(
        &self,
        edge: &ArchivedDotEdgeView<'a>,
    ) -> ArchivedDotEdgeEndpointsView<'a> {
        self.endpoints_of_pair(edge.pair)
    }

    pub fn subgraph_from_bools(
        &self,
        values: impl IntoIterator<Item = bool>,
    ) -> Result<SuBitGraph, String> {
        let values: Vec<bool> = values.into_iter().collect();
        if values.len() != self.n_hedges() {
            return Err(format!(
                "Expected {} subgraph bits, got {}",
                self.n_hedges(),
                values.len()
            ));
        }
        Ok(values.into_iter().collect())
    }

    pub fn subgraph_from_base62(&self, label: &str) -> Result<SuBitGraph, String> {
        SuBitGraph::from_base62(label, self.n_hedges()).ok_or_else(|| {
            format!(
                "Invalid base62 subgraph label for graph with {} hedges",
                self.n_hedges()
            )
        })
    }

    pub fn compass_subgraph(&self, compass_pt: Option<dot_parser::ast::CompassPt>) -> SuBitGraph {
        let target = compass_pt.map(compass_pt_to_u8);
        let mut subgraph = SuBitGraph::empty(self.n_hedges());
        for (hedge, data) in self.graph.hedge_data.0.as_slice().iter().enumerate() {
            if data.compasspt.as_ref().copied() == target {
                subgraph.add(Hedge(hedge));
            }
        }
        subgraph
    }

    fn orientation_of_pair(
        &self,
        pair: &'a Archived<crate::half_edge::involution::HedgePair>,
    ) -> &'a Archived<Orientation> {
        self.orientation_of_hedge(pair.any_hedge())
    }

    fn endpoints_of_pair(
        &self,
        pair: &'a Archived<crate::half_edge::involution::HedgePair>,
    ) -> ArchivedDotEdgeEndpointsView<'a> {
        match pair {
            crate::half_edge::involution::ArchivedHedgePair::Paired { source, sink }
            | crate::half_edge::involution::ArchivedHedgePair::Split { source, sink, .. } => {
                ArchivedDotEdgeEndpointsView {
                    source: Some(self.endpoint_of_hedge(Hedge(source.0.try_into().unwrap()))),
                    sink: Some(self.endpoint_of_hedge(Hedge(sink.0.try_into().unwrap()))),
                }
            }
            crate::half_edge::involution::ArchivedHedgePair::Unpaired { hedge, flow } => {
                let endpoint = self.endpoint_of_hedge(Hedge(hedge.0.try_into().unwrap()));
                match flow {
                    ArchivedFlow::Source => ArchivedDotEdgeEndpointsView {
                        source: Some(endpoint),
                        sink: None,
                    },
                    ArchivedFlow::Sink => ArchivedDotEdgeEndpointsView {
                        source: None,
                        sink: Some(endpoint),
                    },
                }
            }
        }
    }

    fn endpoint_of_hedge(&self, hedge: Hedge) -> ArchivedDotEndpointView<'a> {
        ArchivedDotEndpointView {
            node: self.node_of_hedge(hedge),
            hedge,
            data: &self.graph.hedge_data.0.as_slice()[hedge.0],
        }
    }

    fn node_of_hedge(&self, hedge: Hedge) -> crate::half_edge::NodeIndex {
        let node = &self.graph.node_store.hedge_data.0.as_slice()[hedge.0];
        crate::half_edge::NodeIndex(node.0.try_into().unwrap())
    }

    fn orientation_of_hedge(&self, hedge: Hedge) -> &'a Archived<Orientation> {
        match &self.graph.edge_store.involution.inv.0.as_slice()[hedge.0] {
            crate::half_edge::involution::ArchivedInvolutiveMapping::Identity { data, .. }
            | crate::half_edge::involution::ArchivedInvolutiveMapping::Source { data, .. } => {
                &data.orientation
            }
            crate::half_edge::involution::ArchivedInvolutiveMapping::Sink { source_idx } => {
                self.orientation_of_hedge(Hedge(source_idx.0.try_into().unwrap()))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use dot_parser::ast::CompassPt;

    use super::{DotGraph, DotGraphBytesSet, DotGraphSet};
    use crate::half_edge::{involution::ArchivedOrientation, subgraph::SubSetLike};

    #[test]
    fn archived_graph_bytes_set_supports_indexing_and_iteration() {
        let graphs = DotGraphSet::from_string(
            r#"
            digraph first {
                a [label="A"];
                b [label="B"];
                a:out:n -> b:in:s [label="ab"];
                ext -> a [dir=back, source="src"];
            }
            digraph second { x -> y }
            "#,
        )
        .unwrap();

        let archived_bytes = graphs
            .into_graph_bytes_set::<4096>()
            .unwrap()
            .to_rkyv_bytes::<4096>()
            .unwrap();
        let archived = DotGraphBytesSet::archived_view(&archived_bytes);
        let first = archived.get(0).unwrap();

        assert_eq!(archived.len(), 2);
        assert_eq!(first.global_data().name.as_str(), "first");
        assert_eq!(first.n_vertices(), 3);
        assert_eq!(first.n_edges(), 2);

        let vertices: Vec<_> = first
            .vertex_data()
            .map(|vertex| vertex.data.name.as_ref().unwrap().as_str().to_string())
            .collect();
        assert_eq!(
            vertices,
            vec!["a".to_string(), "b".to_string(), "ext".to_string()]
        );

        let edges: Vec<_> = first
            .edge_data()
            .map(|edge| {
                (
                    edge.data
                        .statements
                        .get("label")
                        .map(|value| value.as_str().to_string()),
                    edge.orientation,
                )
            })
            .collect();
        assert_eq!(edges.len(), 2);
        assert_eq!(edges[0].0.as_deref(), Some("ab"));
        assert!(matches!(edges[0].1, ArchivedOrientation::Default));
        assert!(matches!(edges[1].1, ArchivedOrientation::Reversed));
    }

    #[test]
    fn archived_graph_bytes_support_subgraphs_from_bits_base62_and_compass() {
        let graph: DotGraph = DotGraph::from_string(
            r#"
            digraph first {
                a:out:n -> b:in:s [label="ab"];
                ext -> a [dir=back, source="src"];
            }
            "#,
        )
        .unwrap();

        let bytes = graph.to_rkyv_bytes::<4096>().unwrap();
        let archived = DotGraph::archived_view(&bytes);

        let mut bits = vec![false; archived.n_hedges()];
        bits[0] = true;
        bits[archived.n_hedges() - 1] = true;
        let from_bits = archived.subgraph_from_bools(bits).unwrap();
        assert_eq!(from_bits.n_included(), 2);

        let label = from_bits.string_label();
        let from_label = archived.subgraph_from_base62(&label).unwrap();
        assert_eq!(from_label, from_bits);

        let north = archived.compass_subgraph(Some(CompassPt::N));
        let south = archived.compass_subgraph(Some(CompassPt::S));
        assert_eq!(north.n_included(), 1);
        assert_eq!(south.n_included(), 1);

        let north_edges: Vec<_> = archived.edge_data_of(&north).collect();
        let south_edges: Vec<_> = archived.edge_data_of(&south).collect();
        assert_eq!(north_edges.len(), 1);
        assert_eq!(south_edges.len(), 1);
        assert!(matches!(
            north_edges[0].pair,
            crate::half_edge::involution::ArchivedHedgePair::Paired { .. }
        ));

        let south_vertices: Vec<_> = archived
            .vertex_data_of(&south)
            .map(|vertex| vertex.data.name.as_ref().unwrap().as_str().to_string())
            .collect();
        assert_eq!(south_vertices, vec!["b".to_string()]);
    }
}
