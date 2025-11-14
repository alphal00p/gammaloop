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
    id::{Pattern, Replacement},
    printer::PrintOptions,
    with_default_namespace,
};
use tabled::{builder::Builder, settings::Style};
use typed_index_collections::TiVec;

use crate::{
    GAMMALOOP_NAMESPACE,
    momentum::SignOrZero,
    momentum_sample::{ExternalIndex, LoopIndex},
    signature::{LoopExtSignature, SignatureLike},
    utils::{GS, W_, symbolica_ext::CallSymbol},
};

use super::Graph;

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
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
            .map(|i| ExternalIndex(i))
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
    fn generate_loop_momentum_bases<S: SubGraphLike>(
        &self,
        subgraph: &S,
    ) -> TiVec<LmbIndex, LoopMomentumBasis>
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>;
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
                    (FunctionBuilder::new(GS.emr_mom)
                        .add_arg(usize::from(e))
                        .add_arg(a)
                        .add_args(rep_args)
                        .finish()
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
                                .add_args(&rep_args)
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
        S::Base: ModifySubSet<Hedge>;

    fn lmb<S: SubGraphLike<Base = SuBitGraph>>(&self, subgraph: &S) -> LoopMomentumBasis;

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
        S::Base: ModifySubSet<Hedge>,
    {
        let tree = subgraph.subtract(cotree);
        self.lmb_impl(subgraph, &tree, externals)
    }

    fn empty_lmb(&self) -> LoopMomentumBasis;

    fn dot_lmb<S: SubGraphLike>(&self, subgraph: &S, lmb: &LoopMomentumBasis) -> String;
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

    fn dot_lmb<S: SubGraphLike>(&self, subgraph: &S, lmb: &LoopMomentumBasis) -> String {
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

    fn lmb<S: SubGraphLike<Base = SuBitGraph>>(&self, subgraph: &S) -> LoopMomentumBasis {
        if let Some(i) = subgraph.included_iter().next() {
            let tree =
                SimpleTraversalTree::depth_first_traverse(self, subgraph, &self.node_id(i), None)
                    .unwrap();
            let external = self.full_crown(subgraph);
            // println!("lmb");
            // println!("{}", tree.dot(self));
            return self.lmb_impl(subgraph.included(), &tree.tree_subgraph, external);
        } else {
            panic!("ata")
        };
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
                cut_subgraph.sub(*p);
            }

            if self.count_connected_components(&cut_subgraph) == components {
                // let externals = self.full_crown(subgraph);

                return self.lmb_impl(subgraph.included(), &cut_subgraph, externals.clone());
            }
        }

        panic!("No lmb found for {} and lmb {:?}", self.dot(subgraph), lmb)
    }

    /// The true externals (that will flow through the graph (i.e. not dummy)) are those that are both in the subgraph and in the externals
    fn lmb_impl<S: SubGraphLike + SubSetOps + ModifySubSet<HedgePair> + ModifySubSet<Hedge>>(
        &self,
        subgraph: &S,
        tree: &S,
        true_externals: S,
    ) -> LoopMomentumBasis
    where
        S::Base: ModifySubSet<Hedge>,
    {
        // println!(
        //     "//Lmb of subgraph:\n{}\n//Tree:\n{}\n//Externals\n{}",
        //     self.dot(subgraph),
        //     self.dot(tree),
        //     self.dot(&externals)
        // );
        let Some(h) = subgraph.included_iter().next() else {
            return self.empty_lmb();
        };

        let mut tree_bitvec: SuBitGraph = self.empty_subgraph();
        let mut ext_edges: TiVec<ExternalIndex, EdgeIndex> = vec![].into();

        let mut externals = self.full_crown(subgraph);
        let dep_ext = true_externals.included_iter().next();

        let root_node = self.node_id(dep_ext.unwrap_or(h));

        let cotree = subgraph.subtract(&tree);
        let tree = SimpleTraversalTree::depth_first_traverse(self, tree, &root_node, None).ok();

        if let Some(t) = &tree {
            tree_bitvec = t.tree_subgraph.clone();
        }
        let mut loop_edges: TiVec<LoopIndex, EdgeIndex> = vec![].into();
        let mut cycles = vec![];

        // create all the cycles. The lmb order is in the order of the hedges not in the tree (and their relative order is ultimately due to the relative order of the hedges)
        for (p, e, _) in self.iter_edges_of(&cotree) {
            if let HedgePair::Paired { source, .. } = p {
                if let Some(c) = tree.as_ref().map(|t| t.get_cycle(source, self)).flatten() {
                    if let Some(signed) = SignedCycle::from_cycle(c, source, self) {
                        cycles.push(signed);
                        loop_edges.push(e);
                    }
                }
            }
        }

        // for c in &cycles {
        //     println!("Cycle:{}", self.dot(&c.filter));
        // }

        // The external flows are signed subgraphs (i.e. with only half of the edges to indicate a direction)
        // They always contain the dependent external (except for the flow for the dep ext)
        let mut external_flows: TiVec<ExternalIndex, _> = vec![].into();

        if let Some(ext) = dep_ext {
            ext_edges.push(self[&ext]);
            // print!("{ext}");
            externals.sub(ext);
            external_flows.push((SignOrZero::Zero, self.empty_subgraph()));
        };
        // println!("//Externals {}", self.dot(&externals));

        for (p, e, _) in self.iter_edges_of(&externals) {
            let mut path_to_dep: S = self.empty_subgraph();

            match p {
                HedgePair::Split {
                    source,
                    sink,
                    split,
                } => {
                    if let Some(dep) = dep_ext {
                        path_to_dep.add(dep);
                    }
                    let ext_sign: SignOrZero = split.into();

                    if let Some(tree) = &tree {
                        let ext = match split {
                            Flow::Sink => tree.hedge_parent(sink, self.as_ref()).unwrap(),
                            Flow::Source => tree.hedge_parent(source, self.as_ref()).unwrap(),
                        };
                        for h in tree.ancestor_iter_hedge(ext, self.as_ref()).step_by(2) {
                            path_to_dep.add(h);
                        }
                    }

                    external_flows.push((ext_sign, path_to_dep));
                    ext_edges.push(e);
                }
                HedgePair::Unpaired { hedge, flow } => {
                    let ext_sign: SignOrZero = flow.into();
                    // println!("{hedge}");

                    if true_externals.includes(&hedge) {
                        path_to_dep.add(dep_ext.unwrap());

                        if self.node_id(hedge) == root_node {
                        } else {
                            if let Some(tree) = &tree {
                                let ext = tree.hedge_parent(hedge, self.as_ref()).unwrap();

                                for h in tree.ancestor_iter_hedge(ext, self.as_ref()).step_by(2) {
                                    path_to_dep.add(h);
                                }
                            }
                        }
                    }
                    ext_edges.push(e);
                    external_flows.push((ext_sign, path_to_dep));
                }
                _ => {}
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
            if dep_ext.is_some() {
                external.push(SignOrZero::Zero);
            }

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
                        for (i, (s, e)) in external_flows.iter_enumerated().skip(1) {
                            if ext_edges[i] == eid {
                                external.push(SignOrZero::Plus);
                            } else {
                                if e.includes(source) {
                                    external.push(*s * SignOrZero::Minus);
                                } else if e.includes(sink) {
                                    external.push(*s * SignOrZero::Plus);
                                } else {
                                    external.push(*s * SignOrZero::Zero);
                                }
                            }
                        }
                    } else {
                        external = empty_external;
                    }
                }
                HedgePair::Unpaired { hedge, flow } => {
                    if subgraph.includes(hedge) {
                        for (i, (s, e)) in external_flows.iter_enumerated().skip(1) {
                            if ext_edges[i] == eid {
                                external.push(SignOrZero::Plus);
                            } else {
                                if e.includes(hedge) {
                                    match flow {
                                        Flow::Source => external.push(*s * SignOrZero::Minus),
                                        Flow::Sink => external.push(*s * SignOrZero::Plus),
                                    }
                                } else {
                                    external.push(*s * SignOrZero::Zero);
                                }
                            }
                        }
                    } else {
                        external = empty_external;
                        if externals.includes(hedge) {
                            if let Some((e, _)) = ext_edges.iter().find_position(|a| *a == &eid) {
                                external[e] = SignOrZero::Plus;
                            };
                        }
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
            tree: tree_bitvec,
            edge_signatures: signature,
            ext_edges,
            loop_edges,
        }
    }

    fn generate_loop_momentum_bases<S: SubGraphLike>(
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

        for s in self.all_spanning_trees(subgraph) {
            lmbs.push(self.lmb_impl(subgraph.included(), &s, externals.clone()));
        }
        lmbs
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
    fn dot_lmb<S: SubGraphLike>(&self, subgraph: &S, lmb: &LoopMomentumBasis) -> String {
        self.underlying.dot_lmb(subgraph, lmb)
    }

    fn empty_lmb(&self) -> LoopMomentumBasis {
        self.underlying.empty_lmb()
    }
    fn generate_loop_momentum_bases<S: SubGraphLike>(
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
        self.underlying.generate_loop_momentum_bases(subgraph)
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
        S::Base: ModifySubSet<Hedge>,
    {
        self.underlying.lmb_impl(subgraph, tree, externals)
    }

    fn lmb<S: SubGraphLike<Base = SuBitGraph>>(&self, subgraph: &S) -> LoopMomentumBasis {
        self.underlying.lmb(subgraph)
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
        half_edge::subgraph::{SuBitGraph, SubGraphOps, SubSetOps},
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
                insta::assert_snapshot!(g.dot_lmb(&g.full_filter(), &g.loop_momentum_basis));
            });
        }
    }

    #[test]
    fn compatible_sub_lmb() {
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
        let lmb = g.lmb(&non_dummy);
        let non_dummy_sub_ext = g.full_crown(&subgraph).subtract(&dummy);

        let sub_lmb = g.compatible_sub_lmb(&subgraph, non_dummy_sub_ext, &lmb);

        assert_snapshot!(g.dot_lmb(&non_dummy, &lmb));
        assert_snapshot!(g.dot_lmb(&subgraph, &sub_lmb));
    }
}
