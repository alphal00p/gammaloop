use std::{borrow::Borrow, collections::VecDeque};

use ahash::{HashMap, HashMapExt};
use bincode_trait_derive::{Decode, Encode};
use bitvec::vec::BitVec;
use itertools::Itertools;
use linnet::half_edge::{
    involution::{EdgeData, EdgeIndex, EdgeVec, Flow, Hedge, HedgePair, Orientation},
    subgraph::{cycle::SignedCycle, Inclusion, ModifySubgraph, SubGraph, SubGraphOps},
    tree::SimpleTraversalTree,
    HedgeGraph, NoData, NodeIndex,
};

use color_eyre::Result;
use derive_more::{From, Into};
use eyre::eyre;
use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, FunctionBuilder, Symbol},
    function,
    id::{Pattern, Replacement},
    with_default_namespace,
};
use typed_index_collections::TiVec;

use crate::{
    momentum::{FourMomentum, SignOrZero, ThreeMomentum},
    momentum_sample::{BareMomentumSample, ExternalIndex, ExternalThreeMomenta, LoopIndex},
    signature::{ExternalSignature, LoopExtSignature, LoopSignature, SignatureLike},
    symbolica_ext::CallSymbol,
    utils::{FloatLike, F, GS, W_},
    GAMMALOOP_NAMESPACE,
};

use super::Graph;

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct LoopMomentumBasis {
    pub tree: Option<()>,
    pub loop_edges: TiVec<LoopIndex, EdgeIndex>,
    pub ext_edges: TiVec<ExternalIndex, EdgeIndex>, //It should have length = to number of externals (not number of independent externals)
    pub edge_signatures: EdgeVec<LoopExtSignature>,
}

impl LoopMomentumBasis {
    pub fn swap_loops(&mut self, i: LoopIndex, j: LoopIndex) {
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
    pub fn swap_external(&mut self, i: ExternalIndex, j: ExternalIndex) {
        self.ext_edges.swap(i, j);
        self.edge_signatures = self
            .edge_signatures
            .iter()
            .map(|(eid, a)| {
                let mut a = a.clone();
                a.swap_external(i, j);
                (eid, a)
            })
            .collect();
    }
}

pub trait LMBext {
    fn generate_loop_momentum_bases<S: SubGraph>(
        &self,
        subgraph: &S,
    ) -> TiVec<LmbIndex, LoopMomentumBasis>
    where
        S::Base: SubGraph<Base = S::Base>
            + SubGraphOps
            + Clone
            + ModifySubgraph<HedgePair>
            + ModifySubgraph<Hedge>;
    fn uv_wrapped_replacement<'a, S: SubGraph, I>(
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
                    GS.emr_mom
                        .f(([usize::from(e) as i32], rep_args))
                        .to_pattern(),
                    (FunctionBuilder::new(GS.emr_mom)
                        .add_arg(usize::from(e) as i32)
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

    fn uv_spatial_wrapped_replacement<'a, S: SubGraph, I>(
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
                    GS.emr_vec
                        .f(([usize::from(e) as i32], rep_args))
                        .to_pattern(),
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

    fn normal_emr_replacement<'a, S: SubGraph, I>(
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
                    GS.emr_mom
                        .f(([usize::from(e) as i32], rep_args))
                        .to_pattern(),
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

    fn replacement_impl<'a, S: SubGraph, I>(
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

    fn lmb_impl<S: SubGraph + SubGraphOps + ModifySubgraph<HedgePair> + ModifySubgraph<Hedge>>(
        &self,
        subgraph: &S,
        tree: &S,
        externals: S,
    ) -> LoopMomentumBasis;

    fn lmb<S: SubGraph<Base = BitVec>>(&self, subgraph: &S) -> LoopMomentumBasis;

    fn compatible_sub_lmb<S: SubGraph>(
        &self,
        subgraph: &S,
        lmb: &LoopMomentumBasis,
    ) -> LoopMomentumBasis
    where
        S::Base: SubGraph<Base = S::Base>
            + SubGraphOps
            + Clone
            + ModifySubgraph<HedgePair>
            + ModifySubgraph<Hedge>;

    fn cotree_lmb<S: SubGraph + SubGraphOps + ModifySubgraph<HedgePair> + ModifySubgraph<Hedge>>(
        &self,
        subgraph: &S,
        cotree: &S,
        externals: S,
    ) -> LoopMomentumBasis {
        let tree = subgraph.subtract(cotree);
        self.lmb_impl(subgraph, &tree, externals)
    }

    fn empty_lmb(&self) -> LoopMomentumBasis;

    fn dot_lmb<S: SubGraph>(&self, subgraph: &S, lmb: &LoopMomentumBasis) -> String;
}

pub fn no_filter(pair: &HedgePair) -> bool {
    true
}

impl<E, V, H> LMBext for HedgeGraph<E, V, H> {
    fn empty_lmb(&self) -> LoopMomentumBasis {
        LoopMomentumBasis {
            tree: None,
            loop_edges: vec![].into(),
            ext_edges: vec![].into(),
            edge_signatures: self.new_edgevec(|_, _, _| LoopExtSignature::from((vec![], vec![]))),
        }
    }

    fn dot_lmb<S: SubGraph>(&self, subgraph: &S, lmb: &LoopMomentumBasis) -> String {
        let reps = self.normal_emr_replacement::<_, Atom>(subgraph, lmb, &[], no_filter);

        // for rep in &reps {
        //     println!("{rep}");
        // }
        let emrgraph = self.map_data_ref(
            |_, _, _| "",
            |_, e, _, d| {
                EdgeData::new(
                    GS.emr_mom
                        .f([usize::from(e) as i32])
                        .replace_multiple(&reps),
                    Orientation::Default,
                )
            },
            |_, h| NoData {},
        );
        emrgraph.dot_label(subgraph)
    }

    fn lmb<S: SubGraph<Base = BitVec>>(&self, subgraph: &S) -> LoopMomentumBasis {
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

    fn compatible_sub_lmb<S: SubGraph>(
        &self,
        subgraph: &S,
        lmb: &LoopMomentumBasis,
    ) -> LoopMomentumBasis
    where
        S::Base: SubGraph<Base = S::Base>
            + SubGraphOps
            + Clone
            + ModifySubgraph<HedgePair>
            + ModifySubgraph<Hedge>,
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
                let externals = self.full_crown(subgraph);
                return self.lmb_impl(subgraph.included(), &cut_subgraph, externals);
            }
        }

        panic!("No lmb found for {} and lmb {:?}", self.dot(subgraph), lmb)
    }

    fn lmb_impl<S: SubGraph + SubGraphOps + ModifySubgraph<HedgePair> + ModifySubgraph<Hedge>>(
        &self,
        subgraph: &S,
        tree: &S,
        mut externals: S,
    ) -> LoopMomentumBasis {
        let Some(h) = subgraph.included_iter().next() else {
            return self.empty_lmb();
        };

        let dep_ext = externals.included_iter().next();

        let root_node = self.node_id(dep_ext.unwrap_or(h));

        let cotree = subgraph.subtract(&tree);
        // println!("Tree:{}", self.dot(tree));
        // println!("subgraph:{}", self.dot(subgraph));
        // println!("Cotree:{}", self.dot(&cotree));
        let tree = SimpleTraversalTree::depth_first_traverse(self, tree, &root_node, None).unwrap();

        let mut loop_edges: TiVec<LoopIndex, EdgeIndex> = vec![].into();
        let mut cycles = vec![];
        for (p, e, _) in self.iter_edges_of(&cotree) {
            if let HedgePair::Paired { source, .. } = p {
                if let Some(c) = tree.get_cycle(source, self) {
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

        let mut ext_edges: TiVec<ExternalIndex, EdgeIndex> = vec![].into();
        let mut external_flows: TiVec<ExternalIndex, _> = vec![].into();

        if let Some(ext) = dep_ext {
            ext_edges.push(self[&ext]);
            // print!("{ext}");
            externals.sub(ext);
            external_flows.push((SignOrZero::Zero, self.empty_subgraph()));
        };
        // println!("//Externals {}", self.dot(&externals));

        for (p, e, d) in self.iter_edges_of(&externals) {
            let mut path_to_dep: S = self.empty_subgraph();
            path_to_dep.add(dep_ext.unwrap());

            match p {
                HedgePair::Split {
                    source,
                    sink,
                    split,
                } => {
                    let ext_sign: SignOrZero = split.into();

                    let ext = match split {
                        Flow::Sink => tree.hedge_parent(sink, self.as_ref()).unwrap(),
                        Flow::Source => tree.hedge_parent(source, self.as_ref()).unwrap(),
                    };
                    for h in tree.ancestor_iter_hedge(ext, self.as_ref()).step_by(2) {
                        path_to_dep.add(h);
                    }

                    external_flows.push((ext_sign, path_to_dep));
                    ext_edges.push(e);
                }
                HedgePair::Unpaired { hedge, flow } => {
                    let ext_sign: SignOrZero = flow.into();
                    // println!("{hedge}");

                    if self.node_id(hedge) == root_node {
                        path_to_dep.add(hedge);
                    } else {
                        let ext = tree.hedge_parent(hedge, self.as_ref()).unwrap();

                        for h in tree.ancestor_iter_hedge(ext, self.as_ref()).step_by(2) {
                            path_to_dep.add(h);
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

        let signature = self.new_edgevec(|d, eid, p| {
            let mut internal = vec![];
            let mut external = vec![SignOrZero::Zero];

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
            tree: None,
            edge_signatures: signature,
            ext_edges,
            loop_edges,
        }
    }

    fn generate_loop_momentum_bases<S: SubGraph>(
        &self,
        subgraph: &S,
    ) -> TiVec<LmbIndex, LoopMomentumBasis>
    where
        S::Base: SubGraph<Base = S::Base>
            + SubGraphOps
            + Clone
            + ModifySubgraph<HedgePair>
            + ModifySubgraph<Hedge>,
    {
        let Some(i) = subgraph.included_iter().next() else {
            return vec![].into();
        };

        let mut lmbs: TiVec<LmbIndex, LoopMomentumBasis> = vec![].into();

        let externals = self.full_crown(subgraph);

        for s in self.all_spanning_trees(subgraph) {
            lmbs.push(self.lmb_impl(subgraph.included(), &s, externals.clone()));
        }
        lmbs
    }

    fn replacement_impl<'a, S: SubGraph, I>(
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
    fn dot_lmb<S: SubGraph>(&self, subgraph: &S, lmb: &LoopMomentumBasis) -> String {
        self.underlying.dot_lmb(subgraph, lmb)
    }

    fn empty_lmb(&self) -> LoopMomentumBasis {
        self.underlying.empty_lmb()
    }
    fn generate_loop_momentum_bases<S: SubGraph>(
        &self,
        subgraph: &S,
    ) -> TiVec<LmbIndex, LoopMomentumBasis>
    where
        S::Base: SubGraph<Base = S::Base>
            + SubGraphOps
            + Clone
            + ModifySubgraph<HedgePair>
            + ModifySubgraph<Hedge>,
    {
        self.underlying.generate_loop_momentum_bases(subgraph)
    }

    fn replacement_impl<'a, S: SubGraph, I>(
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

    fn lmb_impl<S: SubGraph + SubGraphOps + ModifySubgraph<HedgePair> + ModifySubgraph<Hedge>>(
        &self,
        subgraph: &S,
        tree: &S,
        externals: S,
    ) -> LoopMomentumBasis {
        self.underlying.lmb_impl(subgraph, tree, externals)
    }

    fn lmb<S: SubGraph<Base = BitVec>>(&self, subgraph: &S) -> LoopMomentumBasis {
        self.underlying.lmb(subgraph)
    }

    fn compatible_sub_lmb<S: SubGraph>(
        &self,
        subgraph: &S,
        lmb: &LoopMomentumBasis,
    ) -> LoopMomentumBasis
    where
        S::Base: SubGraph<Base = S::Base>
            + SubGraphOps
            + Clone
            + ModifySubgraph<HedgePair>
            + ModifySubgraph<Hedge>,
    {
        self.underlying.compatible_sub_lmb(subgraph, lmb)
    }
}

impl LoopMomentumBasis {
    pub(crate) fn spatial_emr<T: FloatLike>(
        &self,
        sample: &BareMomentumSample<T>,
    ) -> Vec<ThreeMomentum<F<T>>> {
        let three_externals: ExternalThreeMomenta<F<T>> = sample
            .external_moms
            .iter()
            .map(|m| m.spatial.clone())
            .collect();
        self.edge_signatures
            .borrow()
            .into_iter()
            .map(|(_, sig)| sig.compute_momentum(&sample.loop_moms, &three_externals))
            .collect()
    }

    pub fn loop_atom<'a, I>(
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

    pub fn ext_atom<'a, I>(
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

    pub fn to_massless_emr<T: FloatLike>(
        &self,
        sample: &BareMomentumSample<T>,
    ) -> Vec<FourMomentum<F<T>>> {
        self.edge_signatures
            .borrow()
            .into_iter()
            .map(|(_, sig)| {
                sig.compute_four_momentum_from_three(&sample.loop_moms, &sample.external_moms)
            })
            .collect()
    }

    pub(crate) fn pattern(&self, edge_id: EdgeIndex) -> Pattern {
        let signature = self.edge_signatures[edge_id].clone();

        let mut atom = Atom::num(0);

        for (i, sign) in signature.internal.into_iter().enumerate() {
            let k = sign
                * Atom::parse(with_default_namespace!(
                    &format!("K({},x_)", i),
                    GAMMALOOP_NAMESPACE
                ))
                .unwrap();

            atom = &atom + &k;
        }

        for (i, sign) in signature.external.into_iter().enumerate() {
            let p = sign
                * Atom::parse(with_default_namespace!(
                    &format!("P({},x_)", i),
                    GAMMALOOP_NAMESPACE
                ))
                .unwrap();
            atom = &atom + &p;
        }

        atom.to_pattern()
    }
}

#[derive(
    Debug,
    Clone,
    Serialize,
    Deserialize,
    Encode,
    bincode_trait_derive::Decode,
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
