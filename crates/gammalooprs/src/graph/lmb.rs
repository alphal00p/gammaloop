use std::fmt::Display;

use bincode_trait_derive::{Decode, Encode};
use derive_more::{From, Into};
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
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, FunctionBuilder, Symbol},
    function,
    id::Replacement,
    printer::PrintOptions,
    symbol,
};
use tabled::{builder::Builder, settings::Style};
use thiserror::Error;
use typed_index_collections::TiVec;

use crate::{
    integrands::process::{amplitude::AmplitudeGraphTerm, cross_section::CrossSectionGraphTerm},
    momentum::{
        SignOrZero,
        sample::{ExternalIndex, LoopIndex},
        signature::{LoopExtSignature, SignatureLike},
    },
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

pub type LmbResult<T> = std::result::Result<T, LmbError>;

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
    pub(crate) fn canonicalize_external_order(&mut self, external_edge_order: &[EdgeIndex]) {
        if external_edge_order.is_empty() {
            return;
        }

        let current_ext_edges = self.ext_edges.clone();
        let mut ordered_ext_edges = external_edge_order.to_vec();
        ordered_ext_edges.extend(
            current_ext_edges
                .iter()
                .copied()
                .filter(|edge| !external_edge_order.contains(edge))
                .sorted(),
        );

        if ordered_ext_edges == current_ext_edges.raw {
            return;
        }

        for (_, signature) in self.edge_signatures.iter_mut() {
            let mut expanded_external = vec![SignOrZero::Zero; ordered_ext_edges.len()];

            for (old_slot, edge) in current_ext_edges.iter_enumerated() {
                let Some(new_slot) = ordered_ext_edges
                    .iter()
                    .position(|ordered_edge| ordered_edge == edge)
                else {
                    continue;
                };
                expanded_external[new_slot] = signature.external[old_slot];
            }

            signature.external = SignatureLike::from_iter(expanded_external);
        }

        self.ext_edges = ordered_ext_edges.into();
    }
}

/// Helpers for constructing loop-momentum bases and turning them into Symbolica
/// replacement rules.
///
/// The replacement methods decompose an edge momentum into its loop-dependent
/// and external-flow parts using a [`LoopMomentumBasis`]. The basis-building
/// methods pick those signatures from spanning forests of a graph or subgraph.
pub trait LMBext {
    /// Enumerate all loop-momentum bases induced by spanning forests of
    /// `subgraph`.
    ///
    /// Each spanning forest covering the same nodes as `subgraph` produces one
    /// basis. Empty subgraphs return an empty list.
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

    /// Enumerate all loop-momentum bases for the full graph.
    fn generate_loop_momentum_bases(&self) -> TiVec<LmbIndex, LoopMomentumBasis>;

    /// Replace `EMRmom(edge, ..)` by a UV-recursion-friendly decomposition.
    ///
    /// The loop-dependent part stays wrapped in `EMRmom(...)`, but its index is
    /// rewritten from the concrete edge id to the loop-basis edge selected by
    /// `lmb`. The external-flow contribution is added explicitly.
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

    /// Spatial-vector variant of [`Self::uv_wrapped_replacement`].
    ///
    /// This uses `EMRvec` for both the matched pattern and the wrapped loop
    /// contribution.
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

    /// Replace `EMRmom(edge, ..)` by the explicit loop-plus-external momentum
    /// carried by that edge.
    ///
    /// `filter_pair` can restrict which edge kinds are rewritten, for example to
    /// skip split or unpaired half-edge pairs in contexts that only want full
    /// propagators.
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

    /// Replace `EMRmom(edge, ..)` by the integrand momentum variables
    /// `K(...) + P(...)`, i.e. `GS.loop_mom(...) + GS.external_mom(...)`.
    ///
    /// Unlike the UV-wrapped replacements, the generated terms are expressed in
    /// the loop/external variable families used in the integrand rather than in
    /// `EMRmom`.
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
    ) -> LmbResult<LoopMomentumBasis>
    where
        S::Base: ModifySubSet<Hedge> + SubGraphLike;

    /// Construct one canonical loop-momentum basis for `subgraph`.
    ///
    /// This uses `subgraph` itself as the forest guide and the full crown of the
    /// subgraph as its external carriers.
    fn lmb_of<S: SubGraphLike<Base = SuBitGraph>>(&self, subgraph: &S) -> LoopMomentumBasis;

    /// Construct the canonical loop-momentum basis for the full graph.
    fn lmb(&self) -> LoopMomentumBasis;

    /// Build the LMB for `outer - shrunken` while each connected component of
    /// `shrunken` acts as a contracted passage node.
    fn shrunken_sub_lmb(
        &self,
        outer: &SuBitGraph,
        shrunken: &InternalSubGraph,
        externals: SuBitGraph,
    ) -> LmbResult<LoopMomentumBasis>;

    /// Construct the canonical shrunken-subgraph LMB using the full crown of
    /// `outer` as external-flow carriers.
    fn shrunken_lmb_of(&self, outer: &SuBitGraph, shrunken: &InternalSubGraph)
    -> LoopMomentumBasis;

    /// Construct a basis for `subgraph` that reuses loop edges from `lmb`
    /// whenever the induced cut still spans the same connected components.
    ///
    /// This is used when descending into a subgraph while keeping its loop
    /// variables compatible with a parent basis.
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
    ) -> LoopMomentumBasis
    where
        S::Base: ModifySubSet<Hedge> + SubGraphLike,
    {
        let tree = subgraph.subtract(cotree);
        self.lmb_impl(subgraph, &tree, externals)
            .unwrap_or_else(|err| panic!("Failed to build cotree loop momentum basis:\n{err}"))
    }

    /// Return the empty basis with no loop or external generators.
    fn empty_lmb(&self) -> LoopMomentumBasis;

    /// Render a DOT graph whose edge labels show the explicit momentum carried by
    /// each edge according to `lmb`.
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

    fn shrunken_sub_lmb(
        &self,
        outer: &SuBitGraph,
        shrunken: &InternalSubGraph,
        externals: SuBitGraph,
    ) -> LmbResult<LoopMomentumBasis> {
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

    fn shrunken_lmb_of(
        &self,
        outer: &SuBitGraph,
        shrunken: &InternalSubGraph,
    ) -> LoopMomentumBasis {
        let externals = self.full_crown(outer);
        self.shrunken_sub_lmb(outer, shrunken, externals)
            .unwrap_or_else(|err| {
                panic!("Failed to build shrunken-subgraph loop momentum basis:\n{err}")
            })
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
                .unwrap_or_else(|err| {
                    panic!("Failed to build loop momentum basis for subgraph:\n{err}")
                })
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

                return self
                    .lmb_impl(subgraph.included(), &cut_subgraph, externals.clone())
                    .unwrap_or_else(|err| {
                        panic!("Failed to build compatible subgraph loop momentum basis:\n{err}")
                    });
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
    ) -> LmbResult<LoopMomentumBasis>
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
                                    for (i, (s, e)) in external_flows.iter_enumerated() {
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
                                    for (i, (s, e)) in external_flows.iter_enumerated() {
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

                        Ok(LoopExtSignature {
                            internal: SignatureLike::from_iter(internal),
                            external: SignatureLike::from_iter(external),
                        })
                    })
                    .collect::<LmbResult<Vec<_>>>()?,
            )
            .map_err(LmbError::EdgeSignatureVector)?;

        let mut lmb = LoopMomentumBasis {
            tree: forest_edge,
            edge_signatures: signature,
            ext_edges,
            loop_edges,
        };
        lmb.canonicalize_external_order(&external_edge_order);

        Ok(lmb)
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
            lmbs.push(
                self.lmb_impl(subgraph.included(), &s, externals.clone())
                    .unwrap_or_else(|err| {
                        panic!("Failed to build loop momentum basis from spanning forest:\n{err}")
                    }),
            );
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

pub trait LMBwithEdges<E: ?Sized> {
    fn lmb_with_loop_edges(&self, lmb_edges: &E) -> LmbResult<LoopMomentumBasis>;
}

impl LMBwithEdges<SuBitGraph> for Graph {
    fn lmb_with_loop_edges(&self, lmb_edges: &SuBitGraph) -> LmbResult<LoopMomentumBasis> {
        let full_filter = self.full_filter();
        let externals = self.internal_crown(&full_filter);
        let cut_graph = full_filter.subtract(lmb_edges);

        self.lmb_impl(&full_filter, &cut_graph, externals)
    }
}
impl LMBwithEdges<[EdgeIndex]> for Graph {
    fn lmb_with_loop_edges(&self, lmb_edges: &[EdgeIndex]) -> LmbResult<LoopMomentumBasis> {
        let mut lmb_edges_subgraph: SuBitGraph = self.empty_subgraph();

        for e in lmb_edges.iter() {
            lmb_edges_subgraph.add(self[e].1);
        }
        self.lmb_with_loop_edges(&lmb_edges_subgraph)
    }
}

impl LMBwithEdges<[&EdgeIndex]> for Graph {
    fn lmb_with_loop_edges(&self, lmb_edges: &[&EdgeIndex]) -> LmbResult<LoopMomentumBasis> {
        let mut lmb_edges_subgraph: SuBitGraph = self.empty_subgraph();

        for e in lmb_edges.iter() {
            lmb_edges_subgraph.add(self[*e].1);
        }
        self.lmb_with_loop_edges(&lmb_edges_subgraph)
    }
}

impl LMBwithEdges<[&EdgeIndex]> for CrossSectionGraphTerm {
    fn lmb_with_loop_edges(&self, lmb_edges: &[&EdgeIndex]) -> LmbResult<LoopMomentumBasis> {
        Ok(self
            .lmbs
            .iter_enumerated()
            .find(|(_, lmb)| lmb_edges.iter().all(|edge| lmb.loop_edges.contains(edge)))
            .ok_or(LmbError::NotLoopEdges {
                loop_edges: lmb_edges.iter().map(|a| a.to_string()).join(","),
                loop_edges_dot: self.graph.debug_dot(),
            })?
            .1
            .clone())
    }
}

impl LMBwithEdges<[EdgeIndex]> for CrossSectionGraphTerm {
    fn lmb_with_loop_edges(&self, lmb_edges: &[EdgeIndex]) -> LmbResult<LoopMomentumBasis> {
        Ok(self
            .lmbs
            .iter_enumerated()
            .find(|(_, lmb)| lmb_edges.iter().all(|edge| lmb.loop_edges.contains(edge)))
            .ok_or(LmbError::NotLoopEdges {
                loop_edges: lmb_edges.iter().map(|a| a.to_string()).join(","),
                loop_edges_dot: self.graph.debug_dot(),
            })?
            .1
            .clone())
    }
}

impl LMBwithEdges<[&EdgeIndex]> for AmplitudeGraphTerm {
    fn lmb_with_loop_edges(&self, lmb_edges: &[&EdgeIndex]) -> LmbResult<LoopMomentumBasis> {
        Ok(self
            .lmbs
            .iter_enumerated()
            .find(|(_, lmb)| lmb_edges.iter().all(|edge| lmb.loop_edges.contains(edge)))
            .ok_or(LmbError::NotLoopEdges {
                loop_edges: lmb_edges.iter().map(|a| a.to_string()).join(","),
                loop_edges_dot: self.graph.debug_dot(),
            })?
            .1
            .clone())
    }
}

impl LMBwithEdges<[EdgeIndex]> for AmplitudeGraphTerm {
    fn lmb_with_loop_edges(&self, lmb_edges: &[EdgeIndex]) -> LmbResult<LoopMomentumBasis> {
        Ok(self
            .lmbs
            .iter_enumerated()
            .find(|(_, lmb)| lmb_edges.iter().all(|edge| lmb.loop_edges.contains(edge)))
            .ok_or(LmbError::NotLoopEdges {
                loop_edges: lmb_edges.iter().map(|a| a.to_string()).join(","),
                loop_edges_dot: self.graph.debug_dot(),
            })?
            .1
            .clone())
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

    fn shrunken_sub_lmb(
        &self,
        outer: &SuBitGraph,
        shrunken: &InternalSubGraph,
        externals: SuBitGraph,
    ) -> LmbResult<LoopMomentumBasis> {
        let mut lmb = self
            .underlying
            .shrunken_sub_lmb(outer, shrunken, externals)?;
        self.canonicalize_lmb_external_order(&mut lmb);
        Ok(lmb)
    }

    fn shrunken_lmb_of(
        &self,
        outer: &SuBitGraph,
        shrunken: &InternalSubGraph,
    ) -> LoopMomentumBasis {
        let mut lmb = self.underlying.shrunken_lmb_of(outer, shrunken);
        self.canonicalize_lmb_external_order(&mut lmb);
        lmb
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
        let Some(_) = subgraph.included_iter().next() else {
            return vec![].into();
        };

        let externals = self.dummy_stripped_external_flows_of(subgraph);
        let mut lmbs: TiVec<LmbIndex, LoopMomentumBasis> = vec![].into();
        for forest in self.underlying.all_spanning_forests_of(subgraph) {
            let mut lmb = self
                .underlying
                .lmb_impl(subgraph.included(), &forest, externals.clone())
                .unwrap_or_else(|err| {
                    panic!("Failed to build loop momentum basis from spanning forest:\n{err}")
                });
            self.canonicalize_lmb_external_order(&mut lmb);
            lmbs.push(lmb);
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
    ) -> LmbResult<LoopMomentumBasis>
    where
        S::Base: ModifySubSet<Hedge> + SubGraphLike,
    {
        let mut lmb = self.underlying.lmb_impl(subgraph, tree, externals)?;
        self.canonicalize_lmb_external_order(&mut lmb);
        Ok(lmb)
    }

    fn lmb_of<S: SubGraphLike<Base = SuBitGraph>>(&self, subgraph: &S) -> LoopMomentumBasis {
        let mut lmb = self.underlying.lmb_of(subgraph);
        self.canonicalize_lmb_external_order(&mut lmb);
        lmb
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
        let mut sub_lmb = self.underlying.compatible_sub_lmb(subgraph, externals, lmb);
        self.canonicalize_lmb_external_order(&mut sub_lmb);
        sub_lmb
    }
}

impl LMBext for &Graph {
    fn dot_lmb_of<S: SubGraphLike>(&self, subgraph: &S, lmb: &LoopMomentumBasis) -> String {
        self.underlying.dot_lmb_of(subgraph, lmb)
    }

    fn lmb(&self) -> LoopMomentumBasis {
        self.lmb_of(&self.underlying.full_filter())
    }

    fn shrunken_sub_lmb(
        &self,
        outer: &SuBitGraph,
        shrunken: &InternalSubGraph,
        externals: SuBitGraph,
    ) -> LmbResult<LoopMomentumBasis> {
        let mut lmb = self
            .underlying
            .shrunken_sub_lmb(outer, shrunken, externals)?;
        self.canonicalize_lmb_external_order(&mut lmb);
        Ok(lmb)
    }

    fn shrunken_lmb_of(
        &self,
        outer: &SuBitGraph,
        shrunken: &InternalSubGraph,
    ) -> LoopMomentumBasis {
        let mut lmb = self.underlying.shrunken_lmb_of(outer, shrunken);
        self.canonicalize_lmb_external_order(&mut lmb);
        lmb
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
        let Some(_) = subgraph.included_iter().next() else {
            return vec![].into();
        };

        let externals = self.dummy_stripped_external_flows_of(subgraph);
        let mut lmbs: TiVec<LmbIndex, LoopMomentumBasis> = vec![].into();
        for forest in self.underlying.all_spanning_forests_of(subgraph) {
            let mut lmb = self
                .underlying
                .lmb_impl(subgraph.included(), &forest, externals.clone())
                .unwrap_or_else(|err| {
                    panic!("Failed to build loop momentum basis from spanning forest:\n{err}")
                });
            self.canonicalize_lmb_external_order(&mut lmb);
            lmbs.push(lmb);
        }
        lmbs
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
    ) -> LmbResult<LoopMomentumBasis>
    where
        S::Base: ModifySubSet<Hedge> + SubGraphLike,
    {
        let mut lmb = self.underlying.lmb_impl(subgraph, tree, externals)?;
        self.canonicalize_lmb_external_order(&mut lmb);
        Ok(lmb)
    }

    fn lmb_of<S: SubGraphLike<Base = SuBitGraph>>(&self, subgraph: &S) -> LoopMomentumBasis {
        let mut lmb = self.underlying.lmb_of(subgraph);
        self.canonicalize_lmb_external_order(&mut lmb);
        lmb
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
        let mut sub_lmb = self.underlying.compatible_sub_lmb(subgraph, externals, lmb);
        self.canonicalize_lmb_external_order(&mut sub_lmb);
        sub_lmb
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
            involution::{EdgeIndex, Hedge},
            subgraph::{Inclusion, InternalSubGraph, ModifySubSet, SuBitGraph, SubSetOps},
        },
        parser::DotGraph,
    };

    use crate::{
        dot,
        graph::{FeynmanGraph, Graph, LMBext, LmbError, parse::IntoGraph},
        initialisation::test_initialise,
        momentum::SignOrZero,
    };

    static SHRUNKEN_LMB_TEST_INIT: std::sync::Once = std::sync::Once::new();
    static SHRUNKEN_LMB_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

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
    fn generated_lmbs_do_not_use_dummy_external_carriers() {
        test_initialise().unwrap();
        let g: Graph = dot!(digraph{
            ext [style=invis]
            edge[num=1 mass=1]
            node[num=1]
            ext->v1:0[id=0 is_dummy=true]
            ext->v1:1[id=1]
            v1->v2[id=2]
            v2->v1[id=3]
            ext->v2:2[id=4]
        })
        .unwrap();

        let lmbs = g.generate_loop_momentum_bases_of(&g.no_dummy());
        assert!(!lmbs.is_empty());

        for lmb in lmbs {
            assert_eq!(
                lmb.ext_edges[crate::momentum::sample::ExternalIndex(0)],
                EdgeIndex::from(0)
            );

            for edge_id in [1, 2, 3, 4].map(EdgeIndex::from) {
                assert_eq!(
                    lmb.edge_signatures[edge_id].external
                        [crate::momentum::sample::ExternalIndex(0)],
                    SignOrZero::Zero,
                    "non-dummy edge {edge_id} uses the dummy external as a generated LMB carrier"
                );
            }
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

        let lmb = g.lmb_impl(&sub, &sub, g.full_crown(&sub)).unwrap();
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

    #[test]
    fn shrunken_connected_subgraph() {
        let _guard = SHRUNKEN_LMB_TEST_LOCK.lock().unwrap();
        SHRUNKEN_LMB_TEST_INIT.call_once(|| test_initialise().unwrap());
        let g: Graph = dot!(digraph{
            edge[num=1 mass=1]
            node[num=1]

            a:0->b:1[id=0]
            b:2->c:3[id=1]
            c:4->a:5[id=2]
        })
        .unwrap();

        let outer = g.full_filter();
        let mut shrunken_filter: SuBitGraph = g.empty_subgraph();
        let shrunken_edge = EdgeIndex::from(0);
        shrunken_filter.add(g[&shrunken_edge].1);
        let shrunken =
            InternalSubGraph::try_new(shrunken_filter, &g.underlying).expect("valid subgraph");
        let remainder = outer.subtract(&shrunken.filter);

        let lmb = g.shrunken_lmb_of(&outer, &shrunken);

        assert!(
            !lmb.loop_edges
                .iter()
                .any(|edge| shrunken.filter.includes(&g[edge].1))
        );
        assert_snapshot!(g.dot_lmb_of(&remainder, &lmb));
    }

    #[test]
    fn shrunken_disconnected_subgraph() {
        let _guard = SHRUNKEN_LMB_TEST_LOCK.lock().unwrap();
        SHRUNKEN_LMB_TEST_INIT.call_once(|| test_initialise().unwrap());
        let g: Graph = dot!(digraph{
            edge[num=1 mass=1]
            node[num=1]

            a:0->b:1[id=0]
            b:2->c:3[id=1]
            c:4->a:5[id=2]

            d:6->e:7[id=3]
            e:8->f:9[id=4]
            f:10->d:11[id=5]

            c:12->d:13[id=6]
        })
        .unwrap();

        let outer = g.full_filter();
        let mut shrunken_filter: SuBitGraph = g.empty_subgraph();
        for edge in [EdgeIndex::from(0), EdgeIndex::from(3)] {
            shrunken_filter.add(g[&edge].1);
        }
        let shrunken =
            InternalSubGraph::try_new(shrunken_filter, &g.underlying).expect("valid subgraph");
        let remainder = outer.subtract(&shrunken.filter);

        let lmb = g.shrunken_lmb_of(&outer, &shrunken);

        assert!(
            !lmb.loop_edges
                .iter()
                .any(|edge| shrunken.filter.includes(&g[edge].1))
        );
        assert_snapshot!(g.dot_lmb_of(&remainder, &lmb));
    }

    #[test]
    fn shrunken_edge_outside_outer_errors() {
        let _guard = SHRUNKEN_LMB_TEST_LOCK.lock().unwrap();
        SHRUNKEN_LMB_TEST_INIT.call_once(|| test_initialise().unwrap());
        let g: Graph = dot!(digraph{
            edge[num=1 mass=1]
            node[num=1]

            a:0->b:1[id=0]
            b:2->c:3[id=1]
            c:4->a:5[id=2]
        })
        .unwrap();

        let mut shrunken_filter: SuBitGraph = g.empty_subgraph();
        let shrunken_edge = EdgeIndex::from(0);
        shrunken_filter.add(g[&shrunken_edge].1);
        let shrunken =
            InternalSubGraph::try_new(shrunken_filter, &g.underlying).expect("valid subgraph");
        let outer = g.full_filter().subtract(&shrunken.filter);

        let result = g.shrunken_sub_lmb(&outer, &shrunken, g.full_crown(&outer));

        assert!(matches!(result, Err(LmbError::ShrunkenOutsideOuter { .. })));
    }
}
