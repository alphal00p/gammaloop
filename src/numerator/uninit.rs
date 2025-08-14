use std::sync::atomic::AtomicUsize;

use linnet::half_edge::subgraph::SubGraph;
use log::debug;
use spenso::{
    network::library::{symbolic::ExplicitKey, TensorLibraryData},
    structure::{
        concrete_index::FlatIndex,
        representation::{LibraryRep, Minkowski, RepName},
        IndexLess, PermutedStructure,
    },
    tensors::{
        data::{DenseTensor, SetTensorData},
        parametric::ParamTensor,
    },
};
use symbolica::atom::Atom;

use crate::{
    new_graph::Graph,
    numerator::aind::Aind,
    utils::{GS, TENSORLIB},
};

use super::{AppliedFeynmanRule, Numerator, UnInit};

static MAXEDGECOUNTER: AtomicUsize = AtomicUsize::new(0);

impl Numerator<UnInit> {
    pub(crate) fn from_new_graph<S: SubGraph>(
        self,
        graph: &Graph,
        subgraph: &S,
        multiply_prefactor: bool,
    ) -> Numerator<AppliedFeynmanRule> {
        debug!("Generating numerator for graph: {}", graph.name);

        debug!("Graph: {}", graph.dot(subgraph));
        let mut num = if multiply_prefactor {
            &graph.global_prefactor.num * &graph.overall_factor
        } else {
            Atom::one()
        };

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
            if p.is_paired() {
                num *= &e.data.num;
            }
        }
        for (_, _, v) in graph.underlying.iter_nodes_of(subgraph) {
            num *= v.get_num();
        }

        Numerator {
            state: AppliedFeynmanRule {
                expr: num,
                state: Default::default(),
            },
        }
    }
}
