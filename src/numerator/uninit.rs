use idenso::metric::MetricSimplifier;
use linnet::half_edge::subgraph::SubGraph;
use log::debug;
use spenso::{
    contraction::Contract,
    network::library::TensorLibraryData,
    structure::{representation::Euclidean, OrderedStructure, ScalarTensor},
    tensors::{data::StorageTensor, parametric::ParamTensor},
};
use symbolica::atom::Atom;

use crate::{graph::BareGraph, numerator::aind::Aind};

use crate::new_graph::Graph;

use super::{AppliedFeynmanRule, Global, GlobalPrefactor, Numerator, UnInit};

impl Numerator<UnInit> {
    pub fn from_graph(
        self,
        graph: &BareGraph,
        prefactor: &GlobalPrefactor,
    ) -> Numerator<AppliedFeynmanRule> {
        debug!("Applying feynman rules");
        let state = AppliedFeynmanRule::from_graph(graph, prefactor);
        debug!(
            "Applied feynman rules:\n\tcolor:{}\n\tcolorless:{}",
            state.color, state.colorless
        );
        Numerator { state }
    }

    pub fn from_new_graph<S: SubGraph>(
        self,
        graph: &Graph,
        subgraph: &S,
        multiply_prefactor: bool,
    ) -> Numerator<AppliedFeynmanRule> {
        debug!("Generating numerator for graph: {}", graph.name);
        let mut colorless_builder = if multiply_prefactor {
            graph.global_prefactor.colorless.clone() * &graph.multiplicity
        } else {
            Atom::one()
        };

        let mut colorful_builder = if multiply_prefactor {
            graph.global_prefactor.color.clone()
        } else {
            Atom::one()
        };

        for (p, _, e) in graph.underlying.iter_edges_of(subgraph) {
            if p.is_paired() {
                colorful_builder = colorful_builder * e.data.color_num.clone();
                colorless_builder = colorless_builder * e.data.spin_num.clone();
            }
        }

        let mut colorless_builder: ParamTensor<OrderedStructure<Euclidean, Aind>> =
            ParamTensor::new_scalar(colorless_builder);
        let mut colorful_builder: ParamTensor<OrderedStructure<Euclidean, Aind>> =
            ParamTensor::new_scalar(colorful_builder);

        for (_, _, v) in graph.underlying.iter_nodes_of(subgraph) {
            colorless_builder = colorless_builder.contract(&v.num_spin).unwrap();
            colorful_builder = colorful_builder.contract(&v.num_color).unwrap();
        }
        Numerator {
            state: AppliedFeynmanRule {
                colorless: colorless_builder.map_data_self(|a| a.simplify_metrics()),

                color: colorful_builder.map_data_self(|a| a.simplify_metrics()),

                state: Default::default(),
            },
        }
    }

    pub fn from_global(
        self,
        global: Atom,
        // _graph: &BareGraph,
        prefactor: &GlobalPrefactor,
    ) -> Numerator<Global> {
        debug!("Setting global numerator");
        let state = {
            let mut global = global;
            global = global * &prefactor.color * &prefactor.colorless;
            Global::new(global.into())
        };
        debug!(
            "Global numerator:\n\tcolor:{}\n\tcolorless:{}",
            state.color, state.colorless
        );
        Numerator { state }
    }

    pub fn from_global_color(
        self,
        global: Atom,
        // _graph: &BareGraph,
        prefactor: &GlobalPrefactor,
    ) -> Numerator<Global> {
        debug!("Setting global numerator");
        let state = {
            let mut global = global;
            global = global * &prefactor.color * &prefactor.colorless;
            Global::new_color(global.into())
        };
        debug!(
            "Global numerator:\n\tcolor:{}\n\tcolorless:{}",
            state.color, state.colorless
        );
        Numerator { state }
    }
}
