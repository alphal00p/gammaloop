use std::fmt::Display;

use bincode_trait_derive::{Decode, Encode};
use derive_more::{From, Into};
use itertools::Itertools;
use linnet::half_edge::{
    HedgeGraph, NoData,
    involution::{EdgeData, EdgeIndex, EdgeVec, Flow, Hedge, HedgePair, Orientation},
    subgraph::{
        Inclusion, ModifySubSet, SuBitGraph, SubGraphLike, SubGraphOps, SubSetLike, SubSetOps,
        cycle::SignedCycle,
    },
    tree::SimpleTraversalTree,
};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, FunctionBuilder, Symbol},
    function,
    id::Replacement,
    printer::PrintOptions,
    symbol,
};
use tabled::{builder::Builder, settings::Style};
use typed_index_collections::TiVec;

use crate::{
    momentum::SignOrZero,
    momentum::sample::{ExternalIndex, LoopIndex},
    momentum::signature::{LoopExtSignature, SignatureLike},
    utils::{GS, W_, symbolica_ext::CallSymbol},
};

use super::Graph;

#[derive(Debug, Clone, Hash, Eq, PartialEq, Serialize, Deserialize, Encode, Decode)]
pub struct LoopMomentumBasis {
    pub tree: SuBitGraph,
    pub loop_edges: TiVec<LoopIndex, EdgeIndex>,
    pub ext_edges: TiVec<ExternalIndex, EdgeIndex>, //It should have length = to number of externals (not number of independent externals)
    pub edge_signatures: EdgeVec<LoopExtSignature>,
}

impl Display for LoopMomentumBasis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut signature = Builder::new();
        let mut title = vec!["edge".to_string()];
        let mut table = Builder::new();

        table.push_column(["Loop id", "edge id"]);
        for (i, item) in self.loop_edges.iter_enumerated() {
            title.push(format!("L{i} {item}"));
            table.push_column(&[i.to_string(), item.to_string()]);
        }
        table.build().with(Style::rounded()).fmt(f)?;
        writeln!(f)?;
        let mut table = Builder::new();

        table.push_column(["External id", "edge id"]);
        for (i, item) in self.ext_edges.iter_enumerated() {
            title.push(format!("Ext{i} {item}"));
            table.push_column(&[i.to_string(), item.to_string()]);
        }
        table.build().with(Style::rounded()).fmt(f)?;
        writeln!(f)?;
        signature.push_record(&title);
        for (i, s) in &self.edge_signatures {
            let mut signs = s
                .internal
                .iter()
                .chain(s.external.iter())
                .map(|s| s.to_string())
                .collect_vec();
            signs.insert(0, i.to_string());
            signature.push_record(signs);
        }
        signature.build().with(Style::rounded()).fmt(f)?;
        writeln!(f)?;

        Ok(())
    }
}

impl LoopMomentumBasis {
    pub(crate) fn ext_from(&self, eid: EdgeIndex) -> Option<ExternalIndex> {
        self.ext_edges
            .iter()
            .position(|&e| e == eid)
            .map(ExternalIndex)
    }
    pub(crate) fn swap_loops(&mut self, i: LoopIndex, j: LoopIndex) {
        self.loop_edges.swap(i, j);
        self.edge_signatures = self
            .edge_signatures
            .iter()
            .map(|(eid, a)| {
                let mut a = a.clone();
                a.swap_loops(i, j);
                (eid, a)
            })
            .collect();
    }

    pub(crate) fn put_loop_to_ext(&mut self, i: LoopIndex) {
        let a = self.loop_edges.remove(i);
        // let ext_id = ExternalIndex::from(self.ext_edges.len());
        self.ext_edges.push(a);
        self.edge_signatures
            .iter_mut()
            .for_each(|(_, s)| s.put_loop_to_ext(i));
    }
    // pub(crate) fn swap_external(&mut self, i: ExternalIndex, j: ExternalIndex) {
    //     self.ext_edges.swap(i, j);
    //     self.edge_signatures = self
    //         .edge_signatures
    //         .iter()
    //         .map(|(eid, a)| {
    //             let mut a = a.clone();
    //             a.swap_external(i, j);
    //             (eid, a)
    //         })
    //         .collect();
    // }
}

pub trait LMBext {
    fn generate_loop_momentum_bases_of<S: SubGraphLike>(
        &self,
        subgraph: &S,
    ) -> TiVec<LmbIndex, LoopMomentumBasis>
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>;

    fn generate_loop_momentum_bases(&self) -> TiVec<LmbIndex, LoopMomentumBasis>;
    fn uv_wrapped_replacement<'a, S: SubSetLike, I>(
        &self,
        subgraph: &S,
        lmb: &LoopMomentumBasis,
        rep_args: &'a [I],
    ) -> Vec<Replacement>
    where
        &'a I: Into<AtomOrView<'a>>,
    {
        self.replacement_impl(
            |e, a, b| {
                Replacement::new(
                    GS.emr_mom.f(([usize::from(e)], rep_args)).to_pattern(),
                    (a.replace(function!(GS.emr_mom, W_.x_))
                        .allow_new_wildcards_on_rhs(true)
                        .with(
                            FunctionBuilder::new(GS.emr_mom)
                                .add_arg(W_.x_)
                                .add_args(rep_args)
                                .finish(),
                        )
                        + b)
                        .to_pattern(),
                )
            },
            subgraph,
            lmb,
            GS.emr_mom,
            GS.emr_mom,
            &[],
            rep_args,
            HedgePair::is_paired,
            true,
        )
    }

    fn uv_spatial_wrapped_replacement<'a, S: SubSetLike, I>(
        &self,
        subgraph: &S,
        lmb: &LoopMomentumBasis,
        rep_args: &'a [I],
    ) -> Vec<Replacement>
    where
        &'a I: Into<AtomOrView<'a>>,
    {
        self.replacement_impl(
            |e, a, b| {
                Replacement::new(
                    GS.emr_vec.f(([usize::from(e)], rep_args)).to_pattern(),
                    (a.replace(function!(GS.emr_vec, W_.x_))
                        .allow_new_wildcards_on_rhs(true)
                        .with(
                            FunctionBuilder::new(GS.emr_vec)
                                .add_arg(W_.x_)
                                .add_args(rep_args)
                                .finish(),
                        )
                        + b)
                        .to_pattern(),
                )
            },
            subgraph,
            lmb,
            GS.emr_vec,
            GS.emr_vec,
            &[],
            rep_args,
            HedgePair::is_paired,
            true,
        )
    }

    fn normal_emr_replacement<'a, S: SubSetLike, I>(
        &self,
        subgraph: &S,
        lmb: &LoopMomentumBasis,
        rep_args: &'a [I],
        filter_pair: fn(&HedgePair) -> bool,
    ) -> Vec<Replacement>
    where
        &'a I: Into<AtomOrView<'a>>,
    {
        self.replacement_impl(
            |e, a, b| {
                Replacement::new(
                    GS.emr_mom.f(([usize::from(e)], rep_args)).to_pattern(),
                    (a + b).to_pattern(),
                )
            },
            subgraph,
            lmb,
            GS.emr_mom,
            GS.emr_mom,
            rep_args,
            rep_args,
            filter_pair,
            true,
        )
    }

    fn integrand_replacement<'a, S: SubSetLike, I>(
        &self,
        subgraph: &S,
        lmb: &LoopMomentumBasis,
        rep_args: &'a [I],
    ) -> Vec<Replacement>
    where
        &'a I: Into<AtomOrView<'a>>,
    {
        self.replacement_impl(
            |e, a, b| {
                Replacement::new(
                    GS.emr_mom.f(([usize::from(e)], rep_args)).to_pattern(),
                    (a + b).to_pattern(),
                )
            },
            subgraph,
            lmb,
            GS.loop_mom,
            GS.external_mom,
            rep_args,
            rep_args,
            no_filter,
            false,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn replacement_impl<'a, S: SubSetLike, I>(
        &self,
        rep: impl Fn(EdgeIndex, Atom, Atom) -> Replacement,
        subgraph: &S,
        lmb: &LoopMomentumBasis,
        loop_symbol: Symbol,
        ext_symbol: Symbol,
        loop_args: &'a [I],
        ext_args: &'a [I],
        filter_pair: fn(&HedgePair) -> bool,
        emr_id: bool,
    ) -> Vec<Replacement>
    where
        &'a I: Into<AtomOrView<'a>>;

    fn lmb_impl<S: SubGraphLike + SubSetOps + ModifySubSet<HedgePair> + ModifySubSet<Hedge>>(
        &self,
        subgraph: &S,
        tree: &S,
        externals: S,
    ) -> LoopMomentumBasis
    where
        S::Base: ModifySubSet<Hedge> + SubGraphLike;

    fn lmb_of<S: SubGraphLike<Base = SuBitGraph>>(&self, subgraph: &S) -> LoopMomentumBasis;

    fn lmb(&self) -> LoopMomentumBasis;

    fn compatible_sub_lmb<S: SubGraphLike>(
        &self,
        subgraph: &S,
        externals: S::Base,
        lmb: &LoopMomentumBasis,
    ) -> LoopMomentumBasis
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>;

    fn cotree_lmb<
        S: SubGraphLike + SubSetOps + SubGraphOps + ModifySubSet<HedgePair> + ModifySubSet<Hedge>,
    >(
        &self,
        subgraph: &S,
        cotree: &S,
        externals: S,
    ) -> LoopMomentumBasis
    where
        S::Base: ModifySubSet<Hedge> + SubGraphLike,
    {
        let tree = subgraph.subtract(cotree);
        self.lmb_impl(subgraph, &tree, externals)
    }

    fn empty_lmb(&self) -> LoopMomentumBasis;

    fn dot_lmb_of<S: SubGraphLike>(&self, subgraph: &S, lmb: &LoopMomentumBasis) -> String;
}

pub(crate) fn no_filter(_pair: &HedgePair) -> bool {
    true
}

impl<E, V, H> LMBext for HedgeGraph<E, V, H> {
    fn empty_lmb(&self) -> LoopMomentumBasis {
        LoopMomentumBasis {
            tree: SuBitGraph::empty(0),
            loop_edges: vec![].into(),
            ext_edges: vec![].into(),
            edge_signatures: self.new_edgevec(|_, _, _| LoopExtSignature::from((vec![], vec![]))),
        }
    }
    fn lmb(&self) -> LoopMomentumBasis {
        self.lmb_of(&self.full_filter())
    }

    fn dot_lmb_of<S: SubGraphLike>(&self, subgraph: &S, lmb: &LoopMomentumBasis) -> String {
        let reps = self.normal_emr_replacement::<_, Atom>(subgraph, lmb, &[], no_filter);

        let emrgraph = self.map_data_ref(
            |_, _, _| "",
            |_, e, _, _| {
                EdgeData::new(
                    GS.emr_mom
                        .f([usize::from(e)])
                        .replace_multiple(&reps)
                        .printer(PrintOptions {
                            color_builtin_symbols: false,
                            color_top_level_sum: false,
                            bracket_level_colors: None,
                            ..Default::default()
                        })
                        .to_string(),
                    Orientation::Default,
                )
            },
            |_, _| NoData {},
        );
        emrgraph.dot_label(subgraph)
    }

    fn lmb_of<S: SubGraphLike<Base = SuBitGraph>>(&self, subgraph: &S) -> LoopMomentumBasis {
        if subgraph.is_empty() {
            self.empty_lmb()
        } else {
            let external = self.full_crown(subgraph);
            self.lmb_impl(subgraph.included(), subgraph.included(), external)
        }
    }

    fn compatible_sub_lmb<S: SubGraphLike>(
        &self,
        subgraph: &S,
        externals: S::Base,
        lmb: &LoopMomentumBasis,
    ) -> LoopMomentumBasis
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>,
    {
        let n_loops = self.cyclotomatic_number(subgraph);
        if n_loops == 0 {
            return self.empty_lmb();
        }

        // the subgraph may have disconnected components in case the of disjoint graphs in a spinney
        let components = self.count_connected_components(subgraph);

        for v in lmb
            .loop_edges
            .iter()
            .filter(|e| {
                let (_, p) = &self[*e];
                subgraph.includes(p)
            })
            .combinations(n_loops)
        {
            let mut cut_subgraph = subgraph.included().clone();

            for eid in v {
                let (_, p) = &self[eid];
                let HedgePair::Paired { source, sink } = p else {
                    continue;
                };

                //this is a self-loop
                if self.node_id(*source) == self.node_id(*sink) {
                    continue;
                }
                cut_subgraph.sub(*p);
            }

            if self.count_connected_components(&cut_subgraph) == components
                && self.number_of_nodes_in_subgraph(&cut_subgraph)
                    == self.number_of_nodes_in_subgraph(subgraph)
            {
                // let externals = self.full_crown(subgraph);

                return self.lmb_impl(subgraph.included(), &cut_subgraph, externals.clone());
            }

            //
        }

        panic!(
            "No lmb found for {} and lmb \n{}",
            self.dot_lmb_of(subgraph, lmb),
            self.dot_lmb_of(&self.full_filter(), lmb)
        )
    }

    /// The true externals (that will flow through the graph (i.e. not dummy)) are those that are both in the subgraph and in the externals
    fn lmb_impl<S: SubGraphLike + SubSetOps + ModifySubSet<HedgePair> + ModifySubSet<Hedge>>(
        &self,
        subgraph: &S,
        forest_guide: &S, //guide for the forest (can be the full subgraph if no guide necessary), however it must cover the same nodes as subgraph
        mut externals: S, //externals to consider for the flow, cannot contain non-subgraph nodes
    ) -> LoopMomentumBasis
    where
        S::Base: ModifySubSet<Hedge> + SubGraphLike,
    {
        // println!(
        //     "//Lmb of subgraph:\n{}\n//Forest_guide:\n{}//Externals:\n{}",
        //     self.dot(subgraph),
        //     self.dot(forest_guide),
        //     self.dot(&externals),
        // );

        if subgraph.is_empty() {
            return self.empty_lmb();
        };

        let mut not_seen = subgraph.clone();
        let mut forest_edge: SuBitGraph = self.empty_subgraph();

        // The external flows are signed subgraphs (i.e. with only half of the edges to indicate a direction)
        // They always contain the dependent external (except for the flow for the dep ext)
        let mut external_flows: TiVec<ExternalIndex, _> = vec![].into();
        let mut ext_edges: TiVec<ExternalIndex, EdgeIndex> = vec![].into();

        let mut loop_edges: TiVec<LoopIndex, EdgeIndex> = vec![].into();
        let mut cycles = vec![];

        loop {
            let Some(mut root) = not_seen.included_iter().next() else {
                break;
            };

            //we keep removing hedges from not_seen until it is empty
            // we need to get the first root
            // if the externals are not yet empty then take from them
            let tree = if let Some(external_root) = externals.included_iter().next() {
                root = external_root;
                let root_node = self.node_id(root);
                let subgraph_tree =
                    SimpleTraversalTree::depth_first_traverse(self, subgraph, &root_node, None)
                        .unwrap_or_else(|_| {
                            panic!(
                                "Externals \n{}\n contains non subgraph nodes:\n{}\n ",
                                self.dot(&externals),
                                self.dot(subgraph)
                            )
                        });

                let external_cover = subgraph_tree.covers(&externals);

                root = external_cover.included_iter().next_back().unwrap();
                let root_node = self.node_id(root);
                let tree = SimpleTraversalTree::depth_first_traverse(
                    self,
                    forest_guide,
                    &root_node,
                    None,
                )
                .unwrap_or_else(|_| {
                    panic!(
                        "Forest guide \n{}\n,does not cover the same nodes as subgraph \n{}\n",
                        self.dot(forest_guide),
                        self.dot(subgraph)
                    )
                }); //select the last half edge in the external cover of this tree as the dependent one

                debug_assert_eq!(
                    subgraph_tree.covers(subgraph),
                    tree.covers(subgraph),
                    "Forest guide \n{}\n,does not cover the same nodes as subgraph \n{}\n",
                    self.dot(forest_guide),
                    self.dot(subgraph)
                );

                // println!(
                //     "//External cover:\n{}//of \n{}",
                //     self.dot(&external_cover),
                //     self.dot(&tree.tree_subgraph)
                // );

                for (p, e, _) in self.iter_edges_of(&external_cover) {
                    let mut path_to_dep: S = self.empty_subgraph();

                    match p {
                        HedgePair::Split {
                            source,
                            sink,
                            split,
                        } => {
                            let hedge = match split {
                                Flow::Sink => sink,
                                Flow::Source => source,
                            };
                            let ext_sign: SignOrZero = split.into();
                            path_to_dep.add(root);

                            if hedge != root {
                                let ext = tree.hedge_parent(hedge, self.as_ref());
                                if let Some(ext) = ext {
                                    for h in tree.ancestor_iter_hedge(ext, self.as_ref()).step_by(2)
                                    {
                                        path_to_dep.add(h);
                                    }
                                }
                            }
                            external_flows.push((ext_sign, path_to_dep));
                            ext_edges.push(e);
                        }
                        HedgePair::Unpaired { hedge, flow } => {
                            let ext_sign: SignOrZero = flow.into();

                            path_to_dep.add(root);
                            if hedge != root {
                                if self.node_id(hedge) == root_node {
                                } else {
                                    let ext = tree.hedge_parent(hedge, self.as_ref()).unwrap();

                                    for h in tree.ancestor_iter_hedge(ext, self.as_ref()).step_by(2)
                                    {
                                        path_to_dep.add(h);
                                    }
                                }
                            }
                            ext_edges.push(e);
                            external_flows.push((ext_sign, path_to_dep));
                        }
                        HedgePair::Paired { source, .. } => {
                            path_to_dep.add(root);

                            let ext_sign: SignOrZero = Flow::Source.into();
                            if source != root {
                                let ext = tree.hedge_parent(source, self.as_ref());
                                if let Some(ext) = ext {
                                    for h in tree.ancestor_iter_hedge(ext, self.as_ref()).step_by(2)
                                    {
                                        path_to_dep.add(h);
                                    }
                                }
                            }
                            external_flows.push((ext_sign, path_to_dep));
                            ext_edges.push(e);
                        }
                    }
                }

                tree
            } else {
                let root_node = self.node_id(root);
                if forest_guide.is_empty() {
                    SimpleTraversalTree::empty(self)
                } else {
                    SimpleTraversalTree::depth_first_traverse(self, forest_guide, &root_node, None)
                        .unwrap_or_else(|_| panic!(
                            "Forest guide \n{}\n,does not cover the same nodes as subgraph \n{}\n",
                            self.dot(forest_guide),
                            self.dot(subgraph)
                        ))
                }
            };

            forest_edge.union_with(&tree.tree_subgraph);

            let mut cover = tree.covers(subgraph);

            for i in self.iter_crown(self.node_id(root)) {
                if subgraph.includes(&i) {
                    cover.add(i);
                }
            }
            //remove all edges in cover+node_crowns from not_seen and externals
            //if the edge is a non-tree, full internal edge then it is a loop edge
            for (p, e, _) in self.iter_edges_of(&cover) {
                match p {
                    HedgePair::Paired { source, sink } => {
                        for h in self.iter_crown(self.node_id(sink)) {
                            not_seen.sub(h);
                            externals.sub(h);
                        }
                        for h in self.iter_crown(self.node_id(source)) {
                            not_seen.sub(h);
                            externals.sub(h);
                        }
                        if !tree.tree_subgraph.includes(&p) {
                            cycles.push(
                                SignedCycle::from_cycle(
                                    tree.get_cycle(source, self).unwrap(),
                                    source,
                                    self,
                                )
                                .unwrap_or_else(|| {
                                    panic!(
                                        "Failed to get cycle from tree:{}\n{}\n{}",
                                        tree.get_cycle(source, self).unwrap().is_circuit(self),
                                        self.dot(&tree.get_cycle(source, self).unwrap().filter),
                                        self.dot(&cover),
                                    )
                                }),
                            );
                            loop_edges.push(e);
                        }
                    }
                    HedgePair::Split {
                        source,
                        sink,
                        split,
                    } => match split {
                        Flow::Sink => {
                            for h in self.iter_crown(self.node_id(sink)) {
                                not_seen.sub(h);
                                externals.sub(h);
                            }
                        }
                        Flow::Source => {
                            for h in self.iter_crown(self.node_id(source)) {
                                not_seen.sub(h);
                                externals.sub(h);
                            }
                        }
                    },
                    HedgePair::Unpaired { hedge, .. } => {
                        for h in self.iter_crown(self.node_id(hedge)) {
                            not_seen.sub(h);
                            externals.sub(h);
                        }
                    }
                }
            }
        }
        // for (i, e) in external_flows.iter().enumerate() {
        //     println!(
        //         "//Ext flow {} for {}: \n{}",
        //         e.0,
        //         ext_edges[ExternalIndex(i)],
        //         self.dot(&e.1)
        //     );
        // }

        let signature = self.new_edgevec(|_, eid, p| {
            let mut internal = vec![];
            let mut external = vec![];
            // if dep_ext.is_some() {
            // external.push(SignOrZero::Zero);
            // }

            let empty_internal = vec![SignOrZero::Zero; cycles.len()];
            let empty_external = vec![SignOrZero::Zero; external_flows.len()];

            match p {
                HedgePair::Paired { source, sink } => {
                    if subgraph.includes(p) {
                        for l in &cycles {
                            if l.filter.includes(source) {
                                internal.push(SignOrZero::Plus);
                            } else if l.filter.includes(sink) {
                                internal.push(SignOrZero::Minus);
                            } else {
                                internal.push(SignOrZero::Zero);
                            }
                        }
                    } else {
                        internal = empty_internal;
                    }
                    if subgraph.intersects(p) {
                        for (i, (s, e)) in external_flows.iter_enumerated() {
                            if ext_edges[i] == eid {
                                if e.includes(source) || e.includes(sink) {
                                    external.push(SignOrZero::Zero); //This is the dependent momentum
                                } else {
                                    external.push(SignOrZero::Plus);
                                }
                            } else if e.includes(source) {
                                external.push(*s * SignOrZero::Minus);
                            } else if e.includes(sink) {
                                external.push(*s * SignOrZero::Plus);
                            } else {
                                external.push(SignOrZero::Zero);
                            }
                        }
                    } else {
                        external = empty_external;
                    }
                }
                HedgePair::Unpaired { hedge, flow } => {
                    if subgraph.includes(hedge) {
                        for (i, (s, e)) in external_flows.iter_enumerated() {
                            if ext_edges[i] == eid {
                                if e.includes(hedge) {
                                    external.push(SignOrZero::Zero); //This is the dependent momentum
                                } else {
                                    external.push(SignOrZero::Plus);
                                }
                            } else if e.includes(hedge) {
                                match flow {
                                    Flow::Source => external.push(*s * SignOrZero::Minus),
                                    Flow::Sink => external.push(*s * SignOrZero::Plus),
                                }
                            } else {
                                external.push(SignOrZero::Zero);
                            }
                        }
                    } else {
                        external = empty_external;
                        if externals.includes(hedge)
                            && let Some((e, _)) = ext_edges.iter().find_position(|a| *a == &eid)
                        {
                            external[e] = SignOrZero::Plus;
                        };
                    }
                    internal = empty_internal;
                }
                _ => {
                    panic!("Split edge on full graph")
                }
            }

            LoopExtSignature {
                internal: SignatureLike::from_iter(internal),
                external: SignatureLike::from_iter(external),
            }
        });

        LoopMomentumBasis {
            tree: forest_edge,
            edge_signatures: signature,
            ext_edges,
            loop_edges,
        }
    }

    fn generate_loop_momentum_bases_of<S: SubGraphLike>(
        &self,
        subgraph: &S,
    ) -> TiVec<LmbIndex, LoopMomentumBasis>
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>,
    {
        let Some(_) = subgraph.included_iter().next() else {
            return vec![].into();
        };

        let mut lmbs: TiVec<LmbIndex, LoopMomentumBasis> = vec![].into();

        let externals = self.full_crown(subgraph);

        for s in self.all_spanning_forests_of(subgraph) {
            // println!("{}", self.dot(&s));
            lmbs.push(self.lmb_impl(subgraph.included(), &s, externals.clone()));
        }
        lmbs
    }

    fn generate_loop_momentum_bases(&self) -> TiVec<LmbIndex, LoopMomentumBasis> {
        self.generate_loop_momentum_bases_of(&self.full_filter())
    }

    #[allow(clippy::too_many_arguments)]
    fn replacement_impl<'a, S: SubSetLike, I>(
        &self,
        rep: impl Fn(EdgeIndex, Atom, Atom) -> Replacement,
        subgraph: &S,
        lmb: &LoopMomentumBasis,
        loop_symbol: Symbol,
        ext_symbol: Symbol,
        loop_args: &'a [I],
        ext_args: &'a [I],
        filter_pair: fn(&HedgePair) -> bool,
        emr_id: bool,
    ) -> Vec<Replacement>
    where
        &'a I: Into<AtomOrView<'a>>,
    {
        let mut reps = vec![];
        for (p, e, _) in self.iter_edges_of(subgraph) {
            if filter_pair(&p) {
                // println!("{e}");
                let loop_expr = lmb.loop_atom(e, loop_symbol, loop_args, emr_id);
                let external_expr = lmb.ext_atom(e, ext_symbol, ext_args, emr_id);

                // println!("{loop_expr}");

                // println!("{external_expr}");
                reps.push(rep(e, loop_expr, external_expr))
            }
        }

        reps
    }
}

impl LMBext for Graph {
    fn dot_lmb_of<S: SubGraphLike>(&self, subgraph: &S, lmb: &LoopMomentumBasis) -> String {
        self.underlying.dot_lmb_of(subgraph, lmb)
    }

    fn generate_loop_momentum_bases(&self) -> TiVec<LmbIndex, LoopMomentumBasis> {
        self.generate_loop_momentum_bases_of(&self.underlying.full_filter())
    }

    fn lmb(&self) -> LoopMomentumBasis {
        self.lmb_of(&self.underlying.full_filter())
    }

    fn empty_lmb(&self) -> LoopMomentumBasis {
        self.underlying.empty_lmb()
    }
    fn generate_loop_momentum_bases_of<S: SubGraphLike>(
        &self,
        subgraph: &S,
    ) -> TiVec<LmbIndex, LoopMomentumBasis>
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>,
    {
        self.underlying.generate_loop_momentum_bases_of(subgraph)
    }

    fn replacement_impl<'a, S: SubSetLike, I>(
        &self,
        rep: impl Fn(EdgeIndex, Atom, Atom) -> Replacement,
        subgraph: &S,
        lmb: &LoopMomentumBasis,
        loop_symbol: Symbol,
        ext_symbol: Symbol,
        loop_args: &'a [I],
        ext_args: &'a [I],
        filter_pair: fn(&HedgePair) -> bool,
        emr_id: bool,
    ) -> Vec<Replacement>
    where
        &'a I: Into<AtomOrView<'a>>,
    {
        self.underlying.replacement_impl(
            rep,
            subgraph,
            lmb,
            loop_symbol,
            ext_symbol,
            loop_args,
            ext_args,
            filter_pair,
            emr_id,
        )
    }

    fn lmb_impl<S: SubGraphLike + SubSetOps + ModifySubSet<HedgePair> + ModifySubSet<Hedge>>(
        &self,
        subgraph: &S,
        tree: &S,
        externals: S,
    ) -> LoopMomentumBasis
    where
        S::Base: ModifySubSet<Hedge> + SubGraphLike,
    {
        self.underlying.lmb_impl(subgraph, tree, externals)
    }

    fn lmb_of<S: SubGraphLike<Base = SuBitGraph>>(&self, subgraph: &S) -> LoopMomentumBasis {
        self.underlying.lmb_of(subgraph)
    }

    fn compatible_sub_lmb<S: SubGraphLike>(
        &self,
        subgraph: &S,
        externals: S::Base,
        lmb: &LoopMomentumBasis,
    ) -> LoopMomentumBasis
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>,
    {
        self.underlying.compatible_sub_lmb(subgraph, externals, lmb)
    }
}

impl LMBext for &Graph {
    fn dot_lmb_of<S: SubGraphLike>(&self, subgraph: &S, lmb: &LoopMomentumBasis) -> String {
        self.underlying.dot_lmb_of(subgraph, lmb)
    }

    fn lmb(&self) -> LoopMomentumBasis {
        self.lmb_of(&self.underlying.full_filter())
    }

    fn empty_lmb(&self) -> LoopMomentumBasis {
        self.underlying.empty_lmb()
    }
    fn generate_loop_momentum_bases_of<S: SubGraphLike>(
        &self,
        subgraph: &S,
    ) -> TiVec<LmbIndex, LoopMomentumBasis>
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>,
    {
        self.underlying.generate_loop_momentum_bases_of(subgraph)
    }

    fn generate_loop_momentum_bases(&self) -> TiVec<LmbIndex, LoopMomentumBasis> {
        self.generate_loop_momentum_bases_of(&self.underlying.full_filter())
    }

    fn replacement_impl<'a, S: SubSetLike, I>(
        &self,
        rep: impl Fn(EdgeIndex, Atom, Atom) -> Replacement,
        subgraph: &S,
        lmb: &LoopMomentumBasis,
        loop_symbol: Symbol,
        ext_symbol: Symbol,
        loop_args: &'a [I],
        ext_args: &'a [I],
        filter_pair: fn(&HedgePair) -> bool,
        emr_id: bool,
    ) -> Vec<Replacement>
    where
        &'a I: Into<AtomOrView<'a>>,
    {
        self.underlying.replacement_impl(
            rep,
            subgraph,
            lmb,
            loop_symbol,
            ext_symbol,
            loop_args,
            ext_args,
            filter_pair,
            emr_id,
        )
    }

    fn lmb_impl<S: SubGraphLike + SubSetOps + ModifySubSet<HedgePair> + ModifySubSet<Hedge>>(
        &self,
        subgraph: &S,
        tree: &S,
        externals: S,
    ) -> LoopMomentumBasis
    where
        S::Base: ModifySubSet<Hedge> + SubGraphLike,
    {
        self.underlying.lmb_impl(subgraph, tree, externals)
    }

    fn lmb_of<S: SubGraphLike<Base = SuBitGraph>>(&self, subgraph: &S) -> LoopMomentumBasis {
        self.underlying.lmb_of(subgraph)
    }

    fn compatible_sub_lmb<S: SubGraphLike>(
        &self,
        subgraph: &S,
        externals: S::Base,
        lmb: &LoopMomentumBasis,
    ) -> LoopMomentumBasis
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>,
    {
        self.underlying.compatible_sub_lmb(subgraph, externals, lmb)
    }
}

impl LoopMomentumBasis {
    pub fn map_to(&self, other: &Self) -> Vec<Atom> {
        let selfmom = symbol!("K");
        let othermom = symbol!("L");
        let mut sys = vec![];

        for (l, e) in self.loop_edges.iter_enumerated() {
            sys.push(
                other.loop_atom::<Atom>(*e, othermom, &[], false)
                    + other.ext_atom::<Atom>(*e, othermom, &[], false)
                    - selfmom.f(&[l.0]),
            )
        }

        let mut vars = vec![];

        for (l, _) in other.loop_edges.iter_enumerated() {
            vars.push(othermom.f(&[l.0]))
        }

        Atom::solve_linear_system::<u8, _, _>(&sys, &vars).unwrap()
    }
    // pub(crate) fn spatial_emr<T: FloatLike>(
    //     &self,
    //     sample: &BareMomentumSample<T>,
    // ) -> Vec<ThreeMomentum<F<T>>> {
    //     let three_externals: ExternalThreeMomenta<F<T>> = sample
    //         .external_moms
    //         .iter()
    //         .map(|m| m.spatial.clone())
    //         .collect();
    //     self.edge_signatures
    //         .borrow()
    //         .into_iter()
    //         .map(|(_, sig)| sig.compute_momentum(&sample.loop_moms, &three_externals))
    //         .collect()
    // }

    pub(crate) fn loop_atom<'a, I>(
        &self,
        edge_id: EdgeIndex,
        mom_symbol: Symbol,
        additional_args: &'a [I],
        emr_id: bool,
    ) -> Atom
    where
        &'a I: Into<AtomOrView<'a>>,
    {
        self.edge_signatures[edge_id].loop_atom(mom_symbol, additional_args, |l| {
            Atom::num(if emr_id {
                usize::from(self.loop_edges[l])
            } else {
                usize::from(l)
            } as i64)
        })
    }

    pub(crate) fn ext_atom<'a, I>(
        &self,
        edge_id: EdgeIndex,
        mom_symbol: Symbol,
        additional_args: &'a [I],
        emr_id: bool,
    ) -> Atom
    where
        &'a I: Into<AtomOrView<'a>>,
    {
        self.edge_signatures[edge_id].ext_atom(mom_symbol, additional_args, |l| {
            Atom::num(if emr_id {
                usize::from(self.ext_edges[l])
            } else {
                usize::from(l)
            } as i64)
        })
    }

    // pub(crate) fn to_massless_emr<T: FloatLike>(
    //     &self,
    //     sample: &BareMomentumSample<T>,
    // ) -> Vec<FourMomentum<F<T>>> {
    //     self.edge_signatures
    //         .borrow()
    //         .into_iter()
    //         .map(|(_, sig)| {
    //             sig.compute_four_momentum_from_three(&sample.loop_moms, &sample.external_moms)
    //         })
    //         .collect()
    // }

    pub(crate) fn edges_are_raised(&self, edge_1: EdgeIndex, edge_2: EdgeIndex) -> bool {
        let sig_1 = &self.edge_signatures[edge_1];
        let sig_2 = &self.edge_signatures[edge_2];
        sig_1.equality_up_to_sign(sig_2)
    }
}

#[derive(
    Debug,
    Clone,
    Serialize,
    Deserialize,
    bincode::Encode,
    bincode::Decode,
    Copy,
    Hash,
    From,
    Into,
    Eq,
    PartialEq,
    Ord,
    PartialOrd,
)]
pub struct LmbIndex(usize);

#[cfg(test)]
pub mod test {

    use insta::assert_snapshot;
    use linnet::{
        half_edge::{
            involution::Hedge,
            subgraph::{ModifySubSet, SuBitGraph, SubSetOps},
        },
        parser::DotGraph,
    };

    use crate::{
        dot,
        graph::{Graph, LMBext, parse::IntoGraph},
        initialisation::test_initialise,
    };

    #[test]
    fn lmb_for_dummy() {
        test_initialise().unwrap();
        let gs: Vec<Graph> = dot!(
            digraph dxda{
                ext [style=invis]
                node[num=1]
                ext->v1:0[id=0 is_dummy=true]
                ext->v1:1[id=1 ]
                ext->v1:2[id=2 ]
            }

            digraph aa{
                ext [style=invis]
                node[num=1]
                ext->v1:0[id=0 is_dummy=true]
                ext->v1:1[id=1 is_dummy=true]
                v1->v2
                ext->v2:2[id=2 ]
            }
        )
        .unwrap();

        for g in gs {
            insta::with_settings!({
                snapshot_suffix=>format!("{}",g.name),
            }, {
                insta::assert_snapshot!(g.dot_lmb_of(&g.full_filter(), &g.loop_momentum_basis));
            });
        }
    }

    #[test]
    fn complicated() {
        test_initialise().unwrap();
        let g: Graph = dot!(digraph{

            edge[num=1 mass=1]
            node[num=1]

            e[style=invis]

            a->c

            a->e
            a->e
            b->c
            d->c
            d->e
            b->e
            a->b->d->a
            b->b1
            b1->b2
            b1->b2
            b1->b2
            b2->e
        })
        .unwrap();
        assert_snapshot!(g.dot_lmb_of(&g.full_filter(), &g.loop_momentum_basis));
        assert_snapshot!(&g.loop_momentum_basis.to_string());
        let _g = g.generate_loop_momentum_bases_of(&g.full_filter());
    }
    #[test]
    fn disconnected() {
        test_initialise().unwrap();
        let g: Graph = dot!(digraph{
            // layout=neato
            e [style=invis]
            edge[num=1 mass=1]
            node[num=1]
            e->v1
            e->v1
            e->v1
            e->v1
            v1->v1

            e->v2
            e->v2
            e->v2

            v3->v3
            v3->v4
            v4->v4


            e->v5
            e->v5->v6
            v6->v7
            v6->v7
            e->v7
        })
        .unwrap();
        assert_snapshot!(g.dot_lmb_of(&g.full_filter(), &g.loop_momentum_basis));
        assert_snapshot!(&g.loop_momentum_basis.to_string());

        let g: Graph = dot!(digraph{

            edge[num=1 mass=1]
            node[num=1]
            a->b
            a->b
            a->b


            c->e
            e->d
            c->d
            c->d
        })
        .unwrap();

        assert_eq!(g.generate_loop_momentum_bases().len(), 15);
    }
    #[test]
    fn subgraph_with_exts_in_loop() {
        test_initialise().unwrap();
        let g: Graph = dot!(digraph{
            edge[num=1 mass=1]
            node[num=1]
            v3:0->v4:1
            v3:3->v4:2
            v3:4->v4:5
        })
        .unwrap();

        let mut sub = g.full_filter();
        sub.sub(Hedge(0));
        sub.sub(Hedge(1));
        let lmb = g.lmb_of(&sub);
        assert_snapshot!(g.dot_lmb_of(&g.full_filter(), &lmb));
        assert_snapshot!(&lmb.to_string());

        let lmb = g.lmb_impl(&sub, &sub, g.full_crown(&sub));
        assert_snapshot!(g.dot_lmb_of(&g.full_filter(), &lmb));
        assert_snapshot!(&lmb.to_string());
    }

    #[test]
    fn compatible_sub_lmb() {
        let g: DotGraph = linnet::dot!(
        digraph{

                                    node[num=1]

                                    v1->v2
                                    v2->v3
                                    v1->v2
                                    v2->v3
                                    v3:s->v1:s
                                    v1:s->v3:s

                                }
        )
        .unwrap();

        let subgraph: SuBitGraph = g.compass_subgraph(Some(dot_parser::ast::CompassPt::S));

        let lmb = g.lmb_of(&subgraph);
        assert_snapshot!(g.dot_lmb_of(&subgraph, &lmb));
        assert_snapshot!(lmb.to_string());

        let g: DotGraph = linnet::dot!(
            digraph dxda{
                            e1 [style=invis]
                            e2 [style=invis]
                            e3 [style=invis]
                            e4 [style=invis]
                            node[num=1]
                            e1->v1:0:n[id=0]
                            e2->v1:1[id=1 ]
                            v1:s->v2:s
                            v2:s->v3:s
                            v3->v1
                            v1:s->v3:s
                            e4->v3
                            e3->v2:2[id=2 ]
                        }

        )
        .unwrap();

        let subgraph: SuBitGraph = g.compass_subgraph(Some(dot_parser::ast::CompassPt::S));

        let dummy: SuBitGraph = g.compass_subgraph(Some(dot_parser::ast::CompassPt::N));
        let non_dummy = g.full_filter().subtract(&dummy);
        let lmb = g.lmb_of(&non_dummy);
        let non_dummy_sub_ext = g.full_crown(&subgraph).subtract(&dummy);

        let sub_lmb = g.compatible_sub_lmb(&subgraph, non_dummy_sub_ext, &lmb);

        assert_snapshot!(g.dot_lmb_of(&non_dummy, &lmb));
        assert_snapshot!(g.dot_lmb_of(&subgraph, &sub_lmb));

        let g: DotGraph = linnet::dot!(
            digraph {
              0:0:s	-> 0:1:s	    [id=0];
              0:2	-> 1:3	        [id=1];
              1:4:s	-> 1:5:s	    [id=2];
            }

        )
        .unwrap();

        let subgraph: SuBitGraph = g.compass_subgraph(Some(dot_parser::ast::CompassPt::S));
        let non_dummy = g.full_filter();
        let lmb = g.lmb_of(&non_dummy);
        let sub_lmb = g.compatible_sub_lmb(&subgraph, non_dummy, &lmb);
        assert_snapshot!(g.dot_lmb_of(&subgraph, &sub_lmb));
    }
}
