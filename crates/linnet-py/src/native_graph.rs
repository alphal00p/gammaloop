//! Python graph payloads use the native unified node topology.
use crate::graph::{EdgeRecord, HalfEdgeRecord, NodeRecord};
use linnet::half_edge::HedgeGraph;

pub(crate) type PyHedgeGraph = HedgeGraph<EdgeRecord, NodeRecord, HalfEdgeRecord>;
