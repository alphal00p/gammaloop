use std::sync::atomic::AtomicUsize;

use linnet::half_edge::subgraph::{ModifySubSet, SubGraphLike, SubSetOps};
use symbolica_utils::AtomPrintExt;
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
        for (_, eid, _) in graph.underlying.iter_edges_of(subgraph) {
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
        }
        let mut selected: linnet::half_edge::subgraph::SuBitGraph =
            graph.underlying.empty_subgraph();
        selected.union_with_iter(subgraph.included_iter());
        let num = feynkit_graph::expressions::GraphExpressions::numerator_of(
            &graph.underlying,
            &selected,
            &graph.underlying.empty_subgraph(),
            |_, vertex| vertex.get_num(),
            |edge| edge.num.value.clone(),
            |edge| edge.is_dummy,
        );

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
        let num = feynkit_graph::expressions::GraphExpressions::numerator_of(
            &graph.underlying,
            subgraph,
            ignore,
            |_, vertex| vertex.get_num(),
            |edge| edge.num.value.clone(),
            |edge| edge.is_dummy,
        );

        debug!( numerator = %num.to_bare_ordered_string(),"Numerator constructed",);

        Numerator {
            state: AppliedFeynmanRule {
                expr: num,
                state: Default::default(),
            },
        }
    }
}
