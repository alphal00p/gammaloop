use std::sync::atomic::AtomicUsize;

use linnet::half_edge::subgraph::subset::SubSet;
use linnet::half_edge::subgraph::{ModifySubSet, SubGraphLike, SubSetLike, SubSetOps};
use linnet::half_edge::{NodeIndex, involution::HedgePair};
use spenso::network::library::TensorLibraryData;
use spenso::shadowing::symbolica_utils::AtomCoreExt;
use symbolica::atom::Atom;
use tracing::debug;
use tracing::instrument;

use crate::graph::Graph;

use super::{AppliedFeynmanRule, Numerator, UnInit};

static MAXEDGECOUNTER: AtomicUsize = AtomicUsize::new(0);

impl Numerator<UnInit> {
    #[instrument(skip_all, fields(graph=%graph.name,debug_dot=%graph.debug_dot(),subgraph_dot=%graph.dot(subgraph)))]
    #[allow(clippy::wrong_self_convention)]
    pub(crate) fn from_new_graph<S: SubGraphLike>(
        self,
        graph: &Graph,
        subgraph: &S,
    ) -> Numerator<AppliedFeynmanRule> {
        let mut num = Atom::one();

        let mut seen: SubSet<NodeIndex> = SubSet::empty(graph.n_nodes());

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
                if !seen[source_n] {
                    seen.add(source_n);
                    num *= graph[source_n].get_num();
                }
                let sink_n = graph.node_id(sink);
                if !seen[sink_n] {
                    seen.add(sink_n);
                    num *= graph[sink_n].get_num();
                }

                num *= &e.data.num.value //.kill_color();
            }
        }

        let notseen = !seen;

        // From all the nodes not yet covered by paired edges, include those included in the subgraph, ignoring dummies
        for node_id in notseen.included_iter() {
            if graph
                .iter_crown(node_id)
                .all(|h| subgraph.includes(&h) || graph[graph[&h]].is_dummy)
            {
                num *= graph[node_id].get_num()
            }
        }

        debug!( numerator = %num.to_bare_ordered_string(),"Numerator constructed",);

        Numerator {
            state: AppliedFeynmanRule {
                expr: num,
                state: Default::default(),
            },
        }
    }
    #[instrument(skip_all, fields(graph=%graph.name,debug_dot=%graph.debug_dot(),subgraph_dot=%graph.dot(subgraph)))]
    #[allow(clippy::wrong_self_convention)]
    pub(crate) fn fill_in_reduced<S: SubGraphLike + SubSetOps>(
        self,
        graph: &Graph,
        subgraph: &S,
        ignore: &S,
    ) -> Numerator<AppliedFeynmanRule> {
        let mut num = Atom::one();

        let mut seen: SubSet<NodeIndex> = SubSet::empty(graph.n_nodes());

        for (nid, _, _) in graph.underlying.iter_nodes_of(ignore) {
            seen.add(nid);
        }
        let not_ignored = subgraph.subtract(ignore);

        for (p, _eid, e) in graph.underlying.iter_edges_of(&not_ignored) {
            if let HedgePair::Paired { source, sink } = p {
                let source_n = graph.node_id(source);
                if !seen[source_n] {
                    seen.add(source_n);
                    num *= graph[source_n].get_num();
                }
                let sink_n = graph.node_id(sink);
                if !seen[sink_n] {
                    seen.add(sink_n);
                    num *= graph[sink_n].get_num();
                }

                num *= &e.data.num.value //.kill_color();
            }
        }

        let notseen = !seen;

        // From all the nodes not yet covered by paired edges, include those included in the subgraph, ignoring dummies
        for node_id in notseen.included_iter() {
            if graph
                .iter_crown(node_id)
                .all(|h| subgraph.includes(&h) || graph[graph[&h]].is_dummy)
            {
                num *= graph[node_id].get_num()
            }
        }

        debug!( numerator = %num.to_bare_ordered_string(),"Numerator constructed",);

        Numerator {
            state: AppliedFeynmanRule {
                expr: num,
                state: Default::default(),
            },
        }
    }
}
