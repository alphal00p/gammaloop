use linnet::half_edge::subgraph::SubGraph;
use log::debug;
use spenso::network::library::TensorLibraryData;
use symbolica::atom::Atom;

use crate::new_graph::Graph;

use super::{AppliedFeynmanRule, Global, GlobalPrefactor, Numerator, UnInit};

impl Numerator<UnInit> {
    pub(crate) fn from_new_graph<S: SubGraph>(
        self,
        graph: &Graph,
        subgraph: &S,
        multiply_prefactor: bool,
    ) -> Numerator<AppliedFeynmanRule> {
        debug!("Generating numerator for graph: {}", graph.name);
        let mut num = if multiply_prefactor {
            graph.global_prefactor.num.clone()
                * &graph.overall_factor
                * &graph.global_prefactor.projector
        } else {
            Atom::one()
        };

        for (p, _, e) in graph.underlying.iter_edges_of(subgraph) {
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

    pub(crate) fn from_global(
        self,
        global: Atom,
        // _graph: &BareGraph,
        prefactor: &GlobalPrefactor,
    ) -> Numerator<Global> {
        debug!("Setting global numerator");
        let state = {
            let mut global = global;
            global = global * &prefactor.num * &prefactor.projector;
            Global::new(global.into())
        };
        debug!("Global numerator:\n\t{}", state.expr);
        Numerator { state }
    }

    pub(crate) fn from_global_color(
        self,
        global: Atom,
        // _graph: &BareGraph,
        prefactor: &GlobalPrefactor,
    ) -> Numerator<Global> {
        debug!("Setting global numerator");
        let state = {
            let mut global = global;
            global = global * &prefactor.num * &prefactor.projector;
            Global::new_color(global.into())
        };
        debug!("Global numerator:\n\t{}", state.expr);
        Numerator { state }
    }
}
