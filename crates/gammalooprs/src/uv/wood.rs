use crate::{
    graph::{Graph, LoopMomentumBasis},
    settings::global::GenerationSettings,
    uv::{Spinney, approx::CutStructure, forest::CutForests},
};
use slotmap::SecondaryMap;
use std::collections::VecDeque;
use tracing::instrument;

use linnet::half_edge::{
    HedgeGraph,
    subgraph::{InternalSubGraph, SubSetOps},
};

// use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};

use super::{
    Forest, Poset, UltravioletGraph,
    approx::Approximation,
    poset::{DAG, DagNode, PosetNode},
};

pub struct CutWoods {
    pub cuts: CutStructure,
    pub woods: Vec<Wood>,
    pub settings: Vec<vakint::VakintSettings>,
}

impl CutWoods {
    #[instrument(skip_all)]
    pub(crate) fn new(cuts: CutStructure, graph: &Graph, settings: &GenerationSettings) -> Self {
        let mut woods = vec![];
        let mut vakint_settings = vec![];
        for cut in cuts.cuts.iter() {
            let mut subgraph = graph.full_filter();
            subgraph.subtract_with(&graph.initial_state_cut.left);
            subgraph.subtract_with(&cut.union);

            let spinneys =
                graph.classified_spinneys(&subgraph, settings, &graph.loop_momentum_basis);

            let wood = Wood::from_spinneys(spinneys, graph);

            let mut lvk_settings = settings.uv.vakint.true_settings();
            // Keep the legacy wood path aligned with the hedge-poset path:
            // the downstream integrand builder extracts the epsilon^0 term, so
            // Vakint must provide one term beyond the maximal pole order.
            lvk_settings.number_of_terms_in_epsilon_expansion = wood.max_loops as i64 + 1;
            vakint_settings.push(lvk_settings);
            woods.push(wood);
        }
        CutWoods {
            cuts,
            woods,
            settings: vakint_settings,
        }
    }

    pub(crate) fn unfold(self, graph: &Graph) -> CutForests {
        CutForests {
            cuts: self.cuts,
            forests: self
                .woods
                .iter()
                .map(|a| a.unfold(graph, &graph.loop_momentum_basis))
                .collect(),
            settings: self.settings,
        }
    }
}

pub struct Wood {
    poset: Poset<Spinney, ()>,
    pub max_loops: usize,
    additional_unions: SecondaryMap<PosetNode, Vec<PosetNode>>,
}

impl Wood {
    pub(crate) fn from_spinneys<E, V, H, I: IntoIterator<Item = Spinney>>(
        s: I,
        graph: impl AsRef<HedgeGraph<E, V, H>>,
    ) -> Self {
        let mut poset = Poset::from_iter(s.into_iter().map(|s| (s, ())));
        let ref_graph = graph.as_ref();

        poset.invert();
        poset.compute_topological_order();
        let mut max_loops = 0;

        let mut unions = SecondaryMap::new();

        for (i, sg) in poset.nodes.iter() {
            let cs = ref_graph.connected_components(sg.data.filter());
            max_loops = max_loops.max(sg.data.max_comp_loop_count());

            if cs.len() > 1 {
                // sg is a disjoint union of spinneys (at the level of half-edges) (strongly disjoint)
                let mut union = vec![];

                for &c in sg.parents.iter() {
                    let mut is_in = 0;
                    for comp in &cs {
                        let comp =
                            InternalSubGraph::cleaned_filter_optimist(comp.clone(), ref_graph);
                        if comp == poset.nodes[c].data.subgraph {
                            // find the components in the wood that this union is made of
                            union.push(c);
                            is_in += 1;
                        }
                    }

                    if is_in > 1 {
                        panic!("is in too many components")
                    }
                }

                unions.insert(i, union.clone());
                // sg.children = union;
            }
        }

        // let coverset = poset.to_cover_set();
        Wood {
            max_loops,
            poset,
            additional_unions: unions,
        }
    }

    #[allow(clippy::type_complexity)]
    fn unfold_bfs<E, V, H, G>(
        &self,
        _graph: &G,
        _lmb: &LoopMomentumBasis,
        dag: &mut DAG<Approximation, DagNode, ()>,
        unions: &mut SecondaryMap<PosetNode, Option<Vec<(PosetNode, Option<DagNode>)>>>,
        root: PosetNode,
    ) -> DagNode
    where
        G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>,
    {
        // let graph = graph.as_ref();
        let mut search_front = VecDeque::new();

        let tree_root = dag.add_node(Approximation::new(self.poset.data(root).clone()));
        search_front.push_front((root, tree_root));

        while let Some((node, parent)) = search_front.pop_front() {
            for c in &self.poset.nodes[node].children {
                if let Some(tagged_union) = unions.get_mut(*c) {
                    // Is this node a disjoint union of spinneys
                    if let Some(mut union) = tagged_union.take() {
                        let mut all_supplied = true;
                        for (p, d) in &mut union {
                            if self.poset.nodes[*p].data == self.poset.nodes[node].data {
                                *d = Some(parent);
                            }
                            if d.is_none() {
                                all_supplied = false;
                            }
                        }
                        if all_supplied {
                            let child =
                                dag.add_node(Approximation::new(self.poset.data(*c).clone()));
                            for (_, d) in union {
                                dag.add_edge(d.unwrap(), child);
                            }
                            search_front.push_front((*c, child));
                        } else {
                            *tagged_union = Some(union);
                        }
                    }
                } else {
                    let child = dag.add_node(Approximation::new(self.poset.data(*c).clone()));
                    dag.add_edge(parent, child);
                    search_front.push_front((*c, child));
                }
            }
        }
        tree_root
    }

    pub(crate) fn unfold<E, V, H, G>(&self, graph: &G, lmb: &LoopMomentumBasis) -> Forest
    where
        G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>,
    {
        let mut dag: DAG<Approximation, DagNode, ()> = DAG::new();

        let root = self.poset.minimum().unwrap();

        let mut unions = SecondaryMap::new();

        for (p, u) in self.additional_unions.iter() {
            let union: Vec<(PosetNode, Option<DagNode>)> = u.iter().map(|i| (*i, None)).collect();
            unions.insert(p, Some(union));
        }

        let _ = self.unfold_bfs(graph, lmb, &mut dag, &mut unions, root);

        Forest { dag }
    }
}
