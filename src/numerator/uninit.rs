use std::sync::atomic::AtomicUsize;

use bitvec::vec::BitVec;
use linnet::half_edge::{involution::HedgePair, subgraph::SubGraph, NodeIndex};
use spenso::network::library::TensorLibraryData;
use symbolica::atom::Atom;
use symbolica::atom::AtomCore;
use tracing::debug;
use tracing::instrument;

use crate::graph::Graph;

use super::{AppliedFeynmanRule, Numerator, UnInit};

static MAXEDGECOUNTER: AtomicUsize = AtomicUsize::new(0);

impl Numerator<UnInit> {
    #[instrument(skip_all, fields(graph=%graph.name,debug_dot=%graph.debug_dot(),subgraph_dot=%graph.dot(subgraph)))]
    pub(crate) fn from_new_graph<S: SubGraph>(
        self,
        graph: &Graph,
        subgraph: &S,
    ) -> Numerator<AppliedFeynmanRule> {
        let mut num = Atom::one();

        let mut seen: BitVec = BitVec::empty(graph.n_nodes());

        for (p, eid, e) in graph.underlying.iter_edges_of(subgraph) {
            let i = MAXEDGECOUNTER.fetch_max(eid.0, std::sync::atomic::Ordering::Relaxed);
            if i == eid.0 {
                // TENSORLIB.write().unwrap().insert_explicit(
                //     ExplicitKey::<Aind>::from_iter(
                //         [Minkowski {}.new_rep(4)],
                //         GS.emr_vec,
                //         Some(vec![Atom::num(eid.0 as i64)]),
                //     )
                //     .map_structure(|s| {
                //         let mut a = ParamTensor::param(DenseTensor::fill(s, Atom::new()).into()).into();
                //         a.set_flat(FlatIndex(0), GS.emr_vec())
                //         a
                //     }),
                // );
            }
            if let HedgePair::Paired { source, sink } = p {
                let source_n = graph.node_id(source);
                if !seen[source_n.0] {
                    seen.set(source_n.0, true);
                    num *= graph[source_n].get_num();
                }
                let sink_n = graph.node_id(sink);
                if !seen[sink_n.0] {
                    seen.set(sink_n.0, true);
                    num *= graph[sink_n].get_num();
                }

                num *= &e.data.num;
            }
        }

        // From all the nodes not yet covered by paired edges, include those included in the subgraph, ignoring dummies
        for n in seen.iter_zeros() {
            let node_id = NodeIndex(n);
            if graph
                .iter_crown(node_id)
                .all(|h| subgraph.includes(&h) || graph[graph[&h]].is_dummy)
            {
                num *= graph[node_id].get_num()
            }
        }

        debug!( numerator = %num.to_canonical_string(),"Numerator constructed",);

        Numerator {
            state: AppliedFeynmanRule {
                expr: num,
                state: Default::default(),
            },
        }
    }
}
