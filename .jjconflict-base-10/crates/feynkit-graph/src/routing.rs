//! Momentum routing shared by standalone diagrams and GammaLoop runtime graphs.
//!
//! External carriers are half-edge selections, independent of whether the full
//! graph pairs them (sewn initial states) or leaves them dangling (amplitudes).

use feynkit_kinematics::{MomentumSignature, SignOrZero, Signature};
use itertools::Itertools;
use linnet::half_edge::{
    HedgeGraph, HedgeGraphError, NoData,
    involution::{EdgeData, EdgeIndex, EdgeVec, Flow, Hedge, HedgePair, Orientation},
    subgraph::{
        Inclusion, InternalSubGraph, ModifySubSet, SuBitGraph, SubGraphLike, SubGraphOps,
        SubSetLike, SubSetOps, cycle::SignedCycle,
    },
    tree::SimpleTraversalTree,
};
use symbolica::{
    atom::{Atom, AtomOrView, FunctionBuilder, Symbol},
    id::Replacement,
};
use thiserror::Error;

/// A spanning-forest routing in the underlying graph's edge coordinates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MomentumBasis {
    pub tree: SuBitGraph,
    pub loop_edges: Vec<EdgeIndex>,
    /// Includes the dependent external carrier of each connected component.
    pub ext_edges: Vec<EdgeIndex>,
    pub edge_signatures: EdgeVec<MomentumSignature>,
}

pub type LmbResult<T> = Result<T, LmbError>;

#[derive(Debug, Error)]
pub enum LmbError {
    #[error(
        "loop edges specified are not actual loop edges in the graph:{loop_edges}:\n{loop_edges_dot}"
    )]
    NotLoopEdges {
        loop_edges: String,
        loop_edges_dot: String,
    },
    #[error("externals\n{externals_dot}\ncontain non-subgraph nodes:\n{subgraph_dot}\n")]
    ExternalsOutsideSubgraph {
        externals_dot: String,
        subgraph_dot: String,
    },
    #[error(
        "external cover is empty for externals\n{externals_dot}\nand subgraph\n{subgraph_dot}\n"
    )]
    EmptyExternalCover {
        externals_dot: String,
        subgraph_dot: String,
    },
    #[error(
        "forest guide\n{forest_guide_dot}\ndoes not cover the same nodes as subgraph\n{subgraph_dot}\n"
    )]
    ForestGuideMismatch {
        forest_guide_dot: String,
        subgraph_dot: String,
    },
    #[error(
        "failed to trace external flow from hedge {hedge} to dependent root {root} in tree\n{tree_dot}\n"
    )]
    ExternalFlowPathMissing {
        hedge: Hedge,
        root: Hedge,
        tree_dot: String,
    },
    #[error("failed to get cycle for source hedge {hedge} in tree:\n{tree_dot}\n")]
    MissingCycle { hedge: Hedge, tree_dot: String },
    #[error(
        "no loop-momentum basis compatible with the parent basis was found for subgraph\n{subgraph_dot}\nparent basis\n{parent_lmb_dot}"
    )]
    NoCompatibleSubLmb {
        subgraph_dot: String,
        parent_lmb_dot: String,
    },
    #[error("failed to get cycle from tree:{is_circuit}\n{cycle_dot}\n{cover_dot}")]
    InvalidCycle {
        is_circuit: bool,
        cycle_dot: String,
        cover_dot: String,
    },
    #[error("split edge on full graph")]
    SplitEdgeOnFullGraph,
    #[error("failed to build edge-signature vector")]
    EdgeSignatureVector(#[source] HedgeGraphError),
    #[error(
        "shrunken subgraph is not contained in outer graph\nouter:\n{outer_dot}\nshrunken:\n{shrunken_dot}"
    )]
    ShrunkenOutsideOuter {
        outer_dot: String,
        shrunken_dot: String,
    },
    #[error("invalid shrunken internal subgraph:\n{shrunken_dot}")]
    InvalidShrunkenSubgraph { shrunken_dot: String },
    #[error(
        "failed to build loop momentum basis after shrinking subgraph\nouter:\n{outer_dot}\nshrunken:\n{shrunken_dot}\nremainder:\n{remainder_dot}"
    )]
    NoShrunkenLmb {
        outer_dot: String,
        shrunken_dot: String,
        remainder_dot: String,
        #[source]
        source: Box<LmbError>,
    },
}

impl MomentumBasis {
    /// Expand and reorder external columns without changing any edge momentum.
    pub fn canonicalize_external_order(&mut self, external_edge_order: &[EdgeIndex]) {
        if external_edge_order.is_empty() {
            return;
        }
        let current = self.ext_edges.clone();
        let mut ordered = external_edge_order.to_vec();
        ordered.extend(
            current
                .iter()
                .copied()
                .filter(|edge| !external_edge_order.contains(edge))
                .sorted(),
        );
        if ordered == current {
            return;
        }
        for (_, signature) in self.edge_signatures.iter_mut() {
            let mut external = vec![SignOrZero::Zero; ordered.len()];
            for (old, edge) in current.iter().enumerate() {
                if let Some(new) = ordered.iter().position(|candidate| candidate == edge) {
                    external[new] = signature
                        .external
                        .get(old)
                        .expect("routing signature covers external carriers");
                }
            }
            signature.external = Signature::new(external);
        }
        self.ext_edges = ordered;
    }

    /// Reclassify a sewn initial-state momentum without changing its routing.
    pub fn put_loop_to_ext(&mut self, index: usize) {
        self.ext_edges.push(self.loop_edges.remove(index));
        for (_, signature) in self.edge_signatures.iter_mut() {
            let mut loops: Vec<_> = signature.loops.iter().collect();
            let sign = loops.remove(index);
            signature.loops = Signature::new(loops);
            signature.external = Signature::new(signature.external.iter().chain([sign]));
        }
    }

    pub fn loop_atom<'a, I>(
        &self,
        edge: EdgeIndex,
        symbol: Symbol,
        args: &'a [I],
        edge_ids: bool,
    ) -> Atom
    where
        &'a I: Into<AtomOrView<'a>>,
    {
        Self::signature_atom(
            &self.edge_signatures[edge].loops,
            &self.loop_edges,
            symbol,
            args,
            edge_ids,
        )
    }

    pub fn ext_atom<'a, I>(
        &self,
        edge: EdgeIndex,
        symbol: Symbol,
        args: &'a [I],
        edge_ids: bool,
    ) -> Atom
    where
        &'a I: Into<AtomOrView<'a>>,
    {
        Self::signature_atom(
            &self.edge_signatures[edge].external,
            &self.ext_edges,
            symbol,
            args,
            edge_ids,
        )
    }

    fn signature_atom<'a, I>(
        signature: &Signature,
        edges: &[EdgeIndex],
        symbol: Symbol,
        args: &'a [I],
        edge_ids: bool,
    ) -> Atom
    where
        &'a I: Into<AtomOrView<'a>>,
    {
        signature
            .iter()
            .enumerate()
            .fold(Atom::Zero, |sum, (index, sign)| {
                if sign == SignOrZero::Zero {
                    return sum;
                }
                let index = if edge_ids { edges[index].0 } else { index };
                let term = FunctionBuilder::new(symbol)
                    .add_arg(index)
                    .add_args(args)
                    .finish();
                match sign {
                    SignOrZero::Plus => sum + term,
                    SignOrZero::Minus => sum - term,
                    SignOrZero::Zero => sum,
                }
            })
    }
}

pub trait MomentumRouting {
    /// Enumerate all loop-momentum bases induced by spanning forests of
    /// `subgraph`.
    ///
    /// Each spanning forest covering the same nodes as `subgraph` produces one
    /// basis. Empty subgraphs return an empty list.
    fn generate_loop_momentum_bases_of<S: SubGraphLike>(&self, subgraph: &S) -> Vec<MomentumBasis>
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>;

    /// Enumerate all loop-momentum bases for the full graph.
    fn generate_loop_momentum_bases(&self) -> Vec<MomentumBasis>;

    /// Core implementation shared by the public replacement constructors.
    ///
    /// For each edge of `subgraph` whose [`HedgePair`] passes `filter_pair`, this
    /// computes the loop-dependent and external-flow atoms from `lmb` and passes
    /// them to `rep`.
    ///
    /// `rep` is the final replacement builder: it receives the original edge id,
    /// the loop-dependent atom, and the external-flow atom, and returns the
    /// `Replacement` inserted into the result vector.
    ///
    /// `loop_symbol` and `ext_symbol` select the function symbol used for the
    /// generated loop and external terms, while `loop_args` and `ext_args` are
    /// appended to those function calls. When `emr_id` is `true`, the generated
    /// loop/external indices are the concrete edge ids stored in the basis;
    /// otherwise they are the compact loop/external basis positions.
    #[allow(clippy::too_many_arguments)]
    fn replacement_impl<'a, S: SubSetLike, I>(
        &self,
        rep: impl Fn(EdgeIndex, Atom, Atom) -> Replacement,
        subgraph: &S,
        lmb: &MomentumBasis,
        loop_symbol: Symbol,
        ext_symbol: Symbol,
        loop_args: &'a [I],
        ext_args: &'a [I],
        filter_pair: fn(&HedgePair) -> bool,
        emr_id: bool,
    ) -> Vec<Replacement>
    where
        &'a I: Into<AtomOrView<'a>>;

    /// Build a loop-momentum basis for `subgraph` using `tree` as the spanning
    /// forest guide and `externals` as the external-flow carriers.
    ///
    /// `tree` must cover the same nodes as `subgraph`. `externals` must only
    /// contain nodes from `subgraph`; it chooses which external edges are treated
    /// as true external flows and which one in each connected component becomes
    /// the dependent external.
    fn lmb_impl<S: SubGraphLike + SubSetOps + ModifySubSet<HedgePair> + ModifySubSet<Hedge>>(
        &self,
        subgraph: &S,
        tree: &S,
        externals: S,
    ) -> LmbResult<MomentumBasis>
    where
        S::Base: ModifySubSet<Hedge> + SubGraphLike;

    /// Construct one canonical loop-momentum basis for `subgraph`.
    ///
    /// This uses `subgraph` itself as the forest guide and the full crown of the
    /// subgraph as its external carriers.
    fn lmb_of<S: SubGraphLike<Base = SuBitGraph>>(&self, subgraph: &S) -> MomentumBasis;

    /// Construct the canonical loop-momentum basis for the full graph.
    fn lmb(&self) -> MomentumBasis;

    /// Build the LMB for `outer - shrunken` while each connected component of
    /// `shrunken` acts as a contracted passage node.
    fn shrunken_sub_lmb(
        &self,
        outer: &SuBitGraph,
        shrunken: &InternalSubGraph,
        externals: SuBitGraph,
    ) -> LmbResult<MomentumBasis>;

    /// Construct the canonical shrunken-subgraph LMB using the full crown of
    /// `outer` as external-flow carriers.
    fn shrunken_lmb_of(&self, outer: &SuBitGraph, shrunken: &InternalSubGraph) -> MomentumBasis;

    /// Construct a basis for `subgraph` that reuses loop edges from `lmb`
    /// whenever the induced cut still spans the same connected components.
    ///
    /// This is used when descending into a subgraph while keeping its loop
    /// variables compatible with a parent basis.
    fn compatible_sub_lmb<S: SubGraphLike>(
        &self,
        subgraph: &S,
        externals: S::Base,
        lmb: &MomentumBasis,
    ) -> MomentumBasis
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>;

    /// Fallible form of [`Self::compatible_sub_lmb`] for callers that must handle
    /// an unavailable parent-compatible basis without panicking. The default
    /// preserves compatibility with external trait implementations that only
    /// implement the original infallible method.
    fn try_compatible_sub_lmb<S: SubGraphLike>(
        &self,
        subgraph: &S,
        externals: S::Base,
        lmb: &MomentumBasis,
    ) -> LmbResult<MomentumBasis>
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>,
    {
        Ok(self.compatible_sub_lmb(subgraph, externals, lmb))
    }

    /// Construct a basis from a chosen cotree of `subgraph`.
    ///
    /// The cotree is converted into the corresponding tree by subtracting it
    /// from `subgraph`, then forwarded to [`Self::lmb_impl`].
    fn cotree_lmb<
        S: SubGraphLike + SubSetOps + SubGraphOps + ModifySubSet<HedgePair> + ModifySubSet<Hedge>,
    >(
        &self,
        subgraph: &S,
        cotree: &S,
        externals: S,
    ) -> MomentumBasis
    where
        S::Base: ModifySubSet<Hedge> + SubGraphLike,
    {
        let tree = subgraph.subtract(cotree);
        self.lmb_impl(subgraph, &tree, externals)
            .unwrap_or_else(|err| panic!("Failed to build cotree loop momentum basis:\n{err}"))
    }

    /// Return the empty basis with no loop or external generators.
    fn empty_lmb(&self) -> MomentumBasis;

    /// Render a DOT graph whose edge labels show the explicit momentum carried by
    /// each edge according to `lmb`.
    fn dot_lmb_of<S: SubGraphLike>(&self, subgraph: &S, lmb: &MomentumBasis) -> String;
}

impl<E, V, H> MomentumRouting for HedgeGraph<E, V, H> {
    fn empty_lmb(&self) -> MomentumBasis {
        MomentumBasis {
            tree: SuBitGraph::empty(0),
            loop_edges: vec![],
            ext_edges: vec![],
            edge_signatures: self.new_edgevec(|_, _, _| MomentumSignature::default()),
        }
    }
    fn lmb(&self) -> MomentumBasis {
        self.lmb_of(&self.full_filter())
    }

    fn shrunken_sub_lmb(
        &self,
        outer: &SuBitGraph,
        shrunken: &InternalSubGraph,
        externals: SuBitGraph,
    ) -> LmbResult<MomentumBasis> {
        let graph_size = self.n_hedges();
        let outer_dot = || {
            if outer.size() == graph_size {
                self.dot(outer)
            } else {
                format!(
                    "invalid outer size {}, expected {}; label {}",
                    outer.size(),
                    graph_size,
                    outer.string_label()
                )
            }
        };
        let shrunken_dot = || {
            if shrunken.size() == graph_size {
                self.dot(shrunken)
            } else {
                format!(
                    "invalid shrunken size {}, expected {}; label {}",
                    shrunken.size(),
                    graph_size,
                    shrunken.string_label()
                )
            }
        };

        if shrunken.size() != graph_size || !shrunken.valid(self) {
            return Err(LmbError::InvalidShrunkenSubgraph {
                shrunken_dot: shrunken_dot(),
            });
        }

        if outer.size() != graph_size || !outer.includes(&shrunken.filter) {
            return Err(LmbError::ShrunkenOutsideOuter {
                outer_dot: outer_dot(),
                shrunken_dot: shrunken_dot(),
            });
        }

        if shrunken.is_empty() {
            return self.lmb_impl(outer, outer, externals);
        }

        let remainder = outer.subtract(&shrunken.filter);
        let contracted_externals = externals.subtract(&shrunken.filter);
        let mut contracted = self.to_ref();

        for component in self.connected_components(shrunken) {
            let Some(root) = component.included_iter().next() else {
                continue;
            };
            let node_data = &self[self.node_id(root)];
            contracted.identify_nodes_of_subgraph_without_self_edges::<_, SuBitGraph>(
                &component, node_data,
            );
        }

        contracted
            .lmb_impl(&remainder, &remainder, contracted_externals)
            .map_err(|source| LmbError::NoShrunkenLmb {
                outer_dot: outer_dot(),
                shrunken_dot: shrunken_dot(),
                remainder_dot: self.dot(&remainder),
                source: Box::new(source),
            })
    }

    fn shrunken_lmb_of(&self, outer: &SuBitGraph, shrunken: &InternalSubGraph) -> MomentumBasis {
        let externals = self.full_crown(outer);
        self.shrunken_sub_lmb(outer, shrunken, externals)
            .unwrap_or_else(|err| {
                panic!("Failed to build shrunken-subgraph loop momentum basis:\n{err}")
            })
    }

    fn lmb_of<S: SubGraphLike<Base = SuBitGraph>>(&self, subgraph: &S) -> MomentumBasis {
        if subgraph.is_empty() {
            self.empty_lmb()
        } else {
            let external = self.full_crown(subgraph);
            self.lmb_impl(subgraph.included(), subgraph.included(), external)
                .unwrap_or_else(|err| {
                    panic!("Failed to build loop momentum basis for subgraph:\n{err}")
                })
        }
    }

    fn compatible_sub_lmb<S: SubGraphLike>(
        &self,
        subgraph: &S,
        externals: S::Base,
        lmb: &MomentumBasis,
    ) -> MomentumBasis
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>,
    {
        self.try_compatible_sub_lmb(subgraph, externals, lmb)
            .unwrap_or_else(|err| {
                panic!("Failed to build compatible subgraph loop momentum basis:\n{err}")
            })
    }

    fn try_compatible_sub_lmb<S: SubGraphLike>(
        &self,
        subgraph: &S,
        externals: S::Base,
        lmb: &MomentumBasis,
    ) -> LmbResult<MomentumBasis>
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>,
    {
        let n_loops = self.cyclotomatic_number(subgraph);
        if n_loops == 0 {
            return Ok(self.empty_lmb());
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

        let full_graph = self.full_filter();
        let parent_lmb_has_full_loop_dimension =
            lmb.loop_edges.len() == self.cyclotomatic_number(&full_graph);
        let (subgraph_dot, parent_lmb_dot) = if parent_lmb_has_full_loop_dimension {
            (
                self.dot_lmb_of(subgraph, lmb),
                self.dot_lmb_of(&full_graph, lmb),
            )
        } else {
            // Momentum-label rendering assumes a dimensionally valid parent LMB. Preserve the
            // topology diagnostics without panicking while constructing the fallible error.
            (
                self.dot(subgraph),
                format!(
                    "parent loop edges {:?}; expected {} loops for\n{}",
                    lmb.loop_edges,
                    self.cyclotomatic_number(&full_graph),
                    self.dot(&full_graph),
                ),
            )
        };

        Err(LmbError::NoCompatibleSubLmb {
            subgraph_dot,
            parent_lmb_dot,
        })
    }

    /// The true externals (that will flow through the graph (i.e. not dummy)) are those that are both in the subgraph and in the externals
    fn lmb_impl<S: SubGraphLike + SubSetOps + ModifySubSet<HedgePair> + ModifySubSet<Hedge>>(
        &self,
        subgraph: &S,
        forest_guide: &S, //guide for the forest (can be the full subgraph if no guide necessary), however it must cover the same nodes as subgraph
        mut externals: S, //externals to consider for the flow, cannot contain non-subgraph nodes
    ) -> LmbResult<MomentumBasis>
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
            return Ok(self.empty_lmb());
        };

        let mut not_seen = subgraph.clone();
        let mut forest_edge: SuBitGraph = self.empty_subgraph();

        // The external flows are signed subgraphs (i.e. with only half of the edges to indicate a direction)
        // They always contain the dependent external (except for the flow for the dep ext)
        let external_edge_order = self
            .iter_edges_of(&externals)
            .map(|(_, edge_id, _)| edge_id)
            .unique()
            .collect_vec();

        let mut external_flows: Vec<_> = vec![];
        let mut ext_edges: Vec<EdgeIndex> = vec![];

        let mut loop_edges: Vec<EdgeIndex> = vec![];
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
                        .map_err(|_| LmbError::ExternalsOutsideSubgraph {
                            externals_dot: self.dot(&externals),
                            subgraph_dot: self.dot(subgraph),
                        })?;

                let external_cover = subgraph_tree.covers(&externals);

                root = external_cover.included_iter().next_back().ok_or_else(|| {
                    LmbError::EmptyExternalCover {
                        externals_dot: self.dot(&externals),
                        subgraph_dot: self.dot(subgraph),
                    }
                })?;
                let root_node = self.node_id(root);
                let tree =
                    SimpleTraversalTree::depth_first_traverse(self, forest_guide, &root_node, None)
                        .map_err(|_| LmbError::ForestGuideMismatch {
                            forest_guide_dot: self.dot(forest_guide),
                            subgraph_dot: self.dot(subgraph),
                        })?; //select the last half edge in the external cover of this tree as the dependent one

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
                                    let ext = tree.hedge_parent(hedge, self.as_ref()).ok_or_else(
                                        || LmbError::ExternalFlowPathMissing {
                                            hedge,
                                            root,
                                            tree_dot: self.dot(&tree.tree_subgraph),
                                        },
                                    )?;

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
                        .map_err(|_| LmbError::ForestGuideMismatch {
                            forest_guide_dot: self.dot(forest_guide),
                            subgraph_dot: self.dot(subgraph),
                        })?
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
                            let cycle = tree.get_cycle(source, self).ok_or_else(|| {
                                LmbError::MissingCycle {
                                    hedge: source,
                                    tree_dot: self.dot(&tree.tree_subgraph),
                                }
                            })?;
                            let cycle_is_circuit = cycle.is_circuit(self);
                            let cycle_dot = self.dot(&cycle.filter);
                            cycles.push(SignedCycle::from_cycle(cycle, source, self).ok_or_else(
                                || LmbError::InvalidCycle {
                                    is_circuit: cycle_is_circuit,
                                    cycle_dot,
                                    cover_dot: self.dot(&cover),
                                },
                            )?);
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

        let signature = self
            .new_edgevec_from_iter(
                self.iter_edges()
                    .map(|(p, eid, _)| -> LmbResult<_> {
                        let mut internal = vec![];
                        let mut external = vec![];
                        // if dep_ext.is_some() {
                        // external.push(SignOrZero::Zero);
                        // }

                        let empty_internal = vec![SignOrZero::Zero; cycles.len()];
                        let empty_external = vec![SignOrZero::Zero; external_flows.len()];

                        match p {
                            HedgePair::Paired { source, sink } => {
                                if subgraph.includes(&p) {
                                    for l in &cycles {
                                        if l.filter.includes(&source) {
                                            internal.push(SignOrZero::Plus);
                                        } else if l.filter.includes(&sink) {
                                            internal.push(SignOrZero::Minus);
                                        } else {
                                            internal.push(SignOrZero::Zero);
                                        }
                                    }
                                } else {
                                    internal = empty_internal;
                                }
                                if subgraph.intersects(&p) {
                                    for (i, (s, e)) in external_flows.iter().enumerate() {
                                        if ext_edges[i] == eid {
                                            if e.includes(&source) || e.includes(&sink) {
                                                external.push(SignOrZero::Zero); //This is the dependent momentum
                                            } else {
                                                external.push(SignOrZero::Plus);
                                            }
                                        } else if e.includes(&source) {
                                            external.push(*s * SignOrZero::Minus);
                                        } else if e.includes(&sink) {
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
                                if subgraph.includes(&hedge) {
                                    for (i, (s, e)) in external_flows.iter().enumerate() {
                                        if ext_edges[i] == eid {
                                            if e.includes(&hedge) {
                                                external.push(SignOrZero::Zero); //This is the dependent momentum
                                            } else {
                                                external.push(SignOrZero::Plus);
                                            }
                                        } else if e.includes(&hedge) {
                                            match flow {
                                                Flow::Source => {
                                                    external.push(*s * SignOrZero::Minus)
                                                }
                                                Flow::Sink => external.push(*s * SignOrZero::Plus),
                                            }
                                        } else {
                                            external.push(SignOrZero::Zero);
                                        }
                                    }
                                } else {
                                    external = empty_external;
                                    if externals.includes(&hedge)
                                        && let Some((e, _)) =
                                            ext_edges.iter().find_position(|a| *a == &eid)
                                    {
                                        external[e] = SignOrZero::Plus;
                                    };
                                }
                                internal = empty_internal;
                            }
                            HedgePair::Split { .. } => {
                                return Err(LmbError::SplitEdgeOnFullGraph);
                            }
                        }

                        Ok(MomentumSignature {
                            loops: Signature::from_iter(internal),
                            external: Signature::from_iter(external),
                        })
                    })
                    .collect::<LmbResult<Vec<_>>>()?,
            )
            .map_err(LmbError::EdgeSignatureVector)?;

        let mut lmb = MomentumBasis {
            tree: forest_edge,
            edge_signatures: signature,
            ext_edges,
            loop_edges,
        };
        lmb.canonicalize_external_order(&external_edge_order);

        Ok(lmb)
    }

    fn generate_loop_momentum_bases_of<S: SubGraphLike>(&self, subgraph: &S) -> Vec<MomentumBasis>
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>,
    {
        let Some(_) = subgraph.included_iter().next() else {
            return vec![];
        };

        let mut lmbs: Vec<MomentumBasis> = vec![];

        let externals = self.full_crown(subgraph);

        for s in self.all_spanning_forests_of(subgraph) {
            // println!("{}", self.dot(&s));
            lmbs.push(
                self.lmb_impl(subgraph.included(), &s, externals.clone())
                    .unwrap_or_else(|err| {
                        panic!("Failed to build loop momentum basis from spanning forest:\n{err}")
                    }),
            );
        }
        lmbs
    }

    fn generate_loop_momentum_bases(&self) -> Vec<MomentumBasis> {
        self.generate_loop_momentum_bases_of(&self.full_filter())
    }

    #[allow(clippy::too_many_arguments)]
    fn replacement_impl<'a, S: SubSetLike, I>(
        &self,
        rep: impl Fn(EdgeIndex, Atom, Atom) -> Replacement,
        subgraph: &S,
        lmb: &MomentumBasis,
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
    fn dot_lmb_of<S: SubGraphLike>(&self, subgraph: &S, lmb: &MomentumBasis) -> String {
        self.map_data_ref(
            |_, _, _| "",
            |_, edge, _, _| {
                EdgeData::new(
                    lmb.edge_signatures[edge].format_momentum(),
                    Orientation::Default,
                )
            },
            |_, _| NoData {},
        )
        .dot_label(subgraph)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use linnet::half_edge::builder::HedgeGraphBuilder;

    fn bubble() -> HedgeGraph<(), ()> {
        let mut builder = HedgeGraphBuilder::new();
        let a = builder.add_node(());
        let b = builder.add_node(());
        builder.add_external_edge(a, (), Orientation::Default, Flow::Sink);
        builder.add_edge(a, b, (), Orientation::Default);
        builder.add_edge(a, b, (), Orientation::Default);
        builder.add_external_edge(b, (), Orientation::Default, Flow::Source);
        builder.build()
    }

    #[test]
    fn all_parallel_edge_routings_conserve_signed_momentum() {
        let graph = bubble();
        let bases = graph.generate_loop_momentum_bases();
        assert_eq!(bases.len(), 2);
        for basis in bases {
            assert_eq!(basis.loop_edges.len(), 1);
            assert_eq!(basis.ext_edges, vec![EdgeIndex(0), EdgeIndex(3)]);
            for node in [
                linnet::half_edge::NodeIndex(0),
                linnet::half_edge::NodeIndex(1),
            ] {
                let mut sum = vec![0; 3];
                for hedge in graph.iter_crown(node) {
                    let sign = if graph.flow(hedge) == Flow::Source {
                        1
                    } else {
                        -1
                    };
                    let signature = &basis.edge_signatures[graph[&hedge]];
                    for (total, coefficient) in sum.iter_mut().zip(
                        signature
                            .loops
                            .integer_coefficients()
                            .into_iter()
                            .chain(signature.external.integer_coefficients()),
                    ) {
                        *total += sign * coefficient;
                    }
                }
                assert_eq!(sum, vec![0; 3]);
            }
        }
    }

    #[test]
    fn promotion_and_external_reordering_preserve_edge_momenta() {
        let graph = bubble();
        let mut basis = graph.lmb();
        let original = basis.clone();
        let loop_edge = basis.loop_edges[0];
        basis.put_loop_to_ext(0);
        basis.canonicalize_external_order(&[loop_edge, EdgeIndex(3), EdgeIndex(0)]);
        assert!(basis.loop_edges.is_empty());
        for (edge, signature) in &basis.edge_signatures {
            let previous = &original.edge_signatures[edge];
            assert_eq!(
                signature.external.integer_coefficients(),
                vec![
                    previous.loops.integer_coefficients()[0],
                    previous.external.integer_coefficients()[1],
                    previous.external.integer_coefficients()[0],
                ]
            );
        }
    }
}

impl crate::LoopMomentumBasis {
    /// Convert graph coordinates without changing the selected momentum basis.
    pub fn from_routing(
        graph: &HedgeGraph<crate::DiagramEdge, crate::DiagramVertex>,
        basis: MomentumBasis,
    ) -> Self {
        let dependent_externals = basis
            .ext_edges
            .iter()
            .enumerate()
            .filter_map(|(position, edge)| {
                (!graph[*edge].is_dummy
                    && basis.edge_signatures[*edge].external.get(position)
                        == Some(SignOrZero::Zero))
                .then_some(crate::EdgeId(edge.0))
            })
            .collect();
        Self {
            tree_edges: graph
                .iter_edges_of(&basis.tree)
                .filter_map(|(pair, edge, data)| {
                    (pair.is_paired() && data.data.external.is_none() && !data.data.is_dummy)
                        .then_some(crate::EdgeId(edge.0))
                })
                .collect(),
            loop_edges: basis
                .loop_edges
                .into_iter()
                .map(|edge| crate::EdgeId(edge.0))
                .collect(),
            external_edges: basis
                .ext_edges
                .into_iter()
                .map(|edge| crate::EdgeId(edge.0))
                .collect(),
            dependent_externals,
            edge_signatures: basis
                .edge_signatures
                .iter()
                .map(|(edge, signature)| (crate::EdgeId(edge.0), signature.clone()))
                .collect(),
        }
    }

    pub fn to_routing(
        &self,
        graph: &HedgeGraph<crate::DiagramEdge, crate::DiagramVertex>,
    ) -> MomentumBasis {
        let mut tree: SuBitGraph = graph.empty_subgraph();
        for edge in &self.tree_edges {
            tree.add(graph[&EdgeIndex(edge.0)].1);
        }
        MomentumBasis {
            tree,
            loop_edges: self
                .loop_edges
                .iter()
                .map(|edge| EdgeIndex(edge.0))
                .collect(),
            ext_edges: self
                .external_edges
                .iter()
                .map(|edge| EdgeIndex(edge.0))
                .collect(),
            edge_signatures: graph
                .new_edgevec(|_, edge, _| self.edge_signatures[&crate::EdgeId(edge.0)].clone()),
        }
    }
}

impl crate::FeynmanDiagram {
    /// Momentum-carrying half-edges, excluding dummy attachments.
    pub fn momentum_subgraph(&self) -> SuBitGraph {
        let mut selected = self.graph.full_filter();
        for (pair, _, edge) in self.graph.iter_edges() {
            if edge.data.is_dummy {
                selected.sub(pair);
            }
        }
        selected
    }

    /// Propagator half-edges, excluding both dangling and sewn externals.
    pub fn internal_subgraph(&self) -> SuBitGraph {
        let mut selected = self.momentum_subgraph();
        for (pair, _, edge) in self.graph.iter_edges() {
            if edge.data.external.is_some() {
                selected.sub(pair);
            }
        }
        selected
    }

    fn normalize_routing(&self, mut basis: MomentumBasis) -> crate::LoopMomentumBasis {
        for (pair, edge_id, edge) in self.graph.iter_edges() {
            if pair.is_paired()
                && edge.data.external.is_some()
                && let Some(position) = basis.loop_edges.iter().position(|edge| *edge == edge_id)
            {
                basis.put_loop_to_ext(position);
            }
        }
        let external_order = self
            .graph
            .iter_edges()
            .filter_map(|(_, edge_id, edge)| edge.data.external.as_ref().map(|_| edge_id))
            .collect::<Vec<_>>();
        basis.canonicalize_external_order(&external_order);
        crate::LoopMomentumBasis::from_routing(&self.graph, basis)
    }

    fn routing_externals(&self, selected: &SuBitGraph) -> SuBitGraph {
        let mut externals = self.graph.full_crown(selected);
        for (pair, _, edge) in self.graph.iter_edges() {
            if edge.data.is_dummy {
                externals.sub(pair);
            }
        }
        externals
    }

    pub fn momentum_basis_of(
        &self,
        subgraph: &SuBitGraph,
    ) -> Result<crate::LoopMomentumBasis, crate::DiagramError> {
        let selected = subgraph.intersection(&self.momentum_subgraph());
        let mut forest = selected.intersection(&self.internal_subgraph());
        let externals = self.routing_externals(&selected);
        forest.union_with(&externals);
        self.graph
            .lmb_impl(&selected, &forest, externals)
            .map(|basis| self.normalize_routing(basis))
            .map_err(|error| crate::DiagramError::InvalidLoopMomentumBasis(error.to_string()))
    }

    pub fn loop_momentum_bases_of(
        &self,
        subgraph: &SuBitGraph,
        limit: usize,
    ) -> Result<Vec<crate::LoopMomentumBasis>, crate::DiagramError> {
        if limit == 0 {
            return Ok(Vec::new());
        }
        let selected = subgraph.intersection(&self.momentum_subgraph());
        let internal = selected.intersection(&self.internal_subgraph());
        if limit == 1 || internal.is_empty() {
            return self.momentum_basis_of(&selected).map(|basis| vec![basis]);
        }
        self.graph
            .all_spanning_forests_of(&internal)
            .into_iter()
            .take(limit)
            .map(|mut forest| {
                let externals = self.routing_externals(&selected);
                forest.union_with(&externals);
                self.graph
                    .lmb_impl(&selected, &forest, externals)
                    .map(|basis| self.normalize_routing(basis))
                    .map_err(|error| {
                        crate::DiagramError::InvalidLoopMomentumBasis(error.to_string())
                    })
            })
            .collect()
    }

    pub fn compatible_momentum_basis_of(
        &self,
        subgraph: &SuBitGraph,
        parent: &crate::LoopMomentumBasis,
    ) -> Result<crate::LoopMomentumBasis, crate::DiagramError> {
        let selected = subgraph.intersection(&self.momentum_subgraph());
        self.graph
            .try_compatible_sub_lmb(
                &selected,
                self.routing_externals(&selected),
                &parent.to_routing(&self.graph),
            )
            .map(|basis| self.normalize_routing(basis))
            .map_err(|error| crate::DiagramError::InvalidLoopMomentumBasis(error.to_string()))
    }

    pub fn contracted_momentum_basis_of(
        &self,
        subgraph: &SuBitGraph,
        contracted: &InternalSubGraph,
    ) -> Result<crate::LoopMomentumBasis, crate::DiagramError> {
        let selected = subgraph.intersection(&self.momentum_subgraph());
        self.graph
            .shrunken_sub_lmb(&selected, contracted, self.routing_externals(&selected))
            .map(|basis| self.normalize_routing(basis))
            .map_err(|error| crate::DiagramError::InvalidLoopMomentumBasis(error.to_string()))
    }

    pub(crate) fn basis_from_tree(
        &self,
        requested: &[crate::EdgeId],
    ) -> Result<crate::LoopMomentumBasis, crate::DiagramError> {
        let internal = self.internal_subgraph();
        let mut forest: SuBitGraph = self.graph.empty_subgraph();
        let unique = requested
            .iter()
            .copied()
            .collect::<std::collections::BTreeSet<_>>();
        if unique.len() != requested.len() {
            return Err(crate::DiagramError::InvalidLoopMomentumBasis(
                "tree edges contain duplicates".into(),
            ));
        }
        for edge in requested {
            let Some((pair, _, _)) = self.graph.iter_edges().find(|(_, id, _)| id.0 == edge.0)
            else {
                return Err(crate::DiagramError::InvalidLoopMomentumBasis(format!(
                    "unknown tree edge {}",
                    edge.0
                )));
            };
            if !pair.is_paired() || !internal.includes(&pair) {
                return Err(crate::DiagramError::InvalidLoopMomentumBasis(format!(
                    "tree edge {} is not internal",
                    edge.0
                )));
            }
            forest.add(pair);
        }
        let missing_nodes = self
            .graph
            .number_of_nodes_in_subgraph(&internal)
            .saturating_sub(self.graph.number_of_nodes_in_subgraph(&forest));
        if self.graph.cyclotomatic_number(&forest) != 0
            || self.graph.count_connected_components(&forest) + missing_nodes
                != self.graph.count_connected_components(&internal)
        {
            return Err(crate::DiagramError::InvalidLoopMomentumBasis(format!(
                "requested tree edges {requested:?} do not form a spanning forest"
            )));
        }
        let selected = self.momentum_subgraph();
        let externals = self.routing_externals(&selected);
        forest.union_with(&externals);
        self.graph
            .lmb_impl(&selected, &forest, externals)
            .map(|basis| self.normalize_routing(basis))
            .map_err(|error| crate::DiagramError::InvalidLoopMomentumBasis(error.to_string()))
    }
}

impl crate::LoopMomentumBasis {
    /// Replace every edge momentum by the same loop and external coordinates
    /// used by numerical routing, retaining any tensor indices on that momentum.
    pub fn momentum_replacements(&self) -> Vec<Replacement> {
        use symbolica::{atom::AtomCore, function, symbol};
        let arguments = symbol!("feynkit_graph::routing_arguments___");
        let loop_edges = self
            .loop_edges
            .iter()
            .map(|edge| EdgeIndex(edge.0))
            .collect::<Vec<_>>();
        let external_edges = self
            .external_edges
            .iter()
            .map(|edge| EdgeIndex(edge.0))
            .collect::<Vec<_>>();
        let args = [Atom::var(arguments)];
        self.edge_signatures
            .iter()
            .map(|(edge, signature)| {
                let internal = MomentumBasis::signature_atom(
                    &signature.loops,
                    &loop_edges,
                    crate::symbols::loop_momentum(),
                    &args,
                    false,
                );
                let external = MomentumBasis::signature_atom(
                    &signature.external,
                    &external_edges,
                    crate::symbols::external_momentum(),
                    &args,
                    false,
                );
                Replacement::new(
                    function!(crate::symbols::momentum(), edge.0, arguments).to_pattern(),
                    (internal + external).to_pattern(),
                )
            })
            .collect()
    }

    pub fn route_expression(&self, expression: &Atom) -> Atom {
        use symbolica::atom::AtomCore;
        expression.replace_multiple(&self.momentum_replacements())
    }
}

impl crate::DiagramEdge {
    /// Recover the canonical unoriented endpoints and signed particle species
    /// used by GammaLoop when ordering propagators before momentum selection.
    pub fn canonical_order_key(
        &self,
        model: &feynkit_model::Model,
        edge: crate::EdgeId,
        endpoints: crate::EdgeEndpoints,
    ) -> Result<
        (
            Option<crate::VertexId>,
            Option<crate::VertexId>,
            i64,
            crate::EdgeId,
        ),
        crate::DiagramError,
    > {
        let (source, target, particle) = if endpoints.source <= endpoints.target {
            (endpoints.source, endpoints.target, self.particle)
        } else {
            (
                endpoints.target,
                endpoints.source,
                model.particle_by_id(self.particle)?.antiparticle,
            )
        };
        Ok((
            source,
            target,
            model.particle_by_id(particle)?.pdg_code,
            edge,
        ))
    }
}
