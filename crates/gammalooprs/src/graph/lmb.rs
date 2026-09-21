use std::fmt::Display;

use bincode_trait_derive::{Decode, Encode};
use derive_more::{From, Into};
use itertools::Itertools;
use linnet::half_edge::{
    HedgeGraph, NoData,
    involution::{EdgeData, EdgeIndex, EdgeVec, Hedge, HedgePair, Orientation},
    subgraph::{
        InternalSubGraph, ModifySubSet, SuBitGraph, SubGraphLike, SubGraphOps, SubSetLike,
        SubSetOps,
    },
};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, FunctionBuilder, Symbol},
    function,
    id::Replacement,
    poly::PolyVariable,
    printer::PrintOptions,
    symbol,
};
use tabled::{builder::Builder, settings::Style};
use typed_index_collections::TiVec;

use crate::{
    integrands::process::{amplitude::AmplitudeGraphTerm, cross_section::CrossSectionGraphTerm},
    momentum::{
        sample::{ExternalIndex, LoopIndex},
        signature::LoopExtSignature,
    },
    utils::{GS, W_},
};

use super::Graph;

#[derive(Debug, Clone, Hash, Eq, PartialEq, Serialize, Deserialize, Encode, Decode)]
pub struct LoopMomentumBasis {
    pub tree: SuBitGraph,
    pub loop_edges: TiVec<LoopIndex, EdgeIndex>,
    pub ext_edges: TiVec<ExternalIndex, EdgeIndex>, //It should have length = to number of externals (not number of independent externals)
    pub edge_signatures: EdgeVec<LoopExtSignature>,
}

pub use feynkit_graph::routing::{LmbError, LmbResult};

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
    pub(crate) fn accounted_bytes(&self) -> usize {
        // Cached LMBs are cloned to their live lengths. Allow a word of slack
        // for the bitset and charge the owned coordinate vectors explicitly.
        std::mem::size_of::<Self>()
            + self.tree.size().div_ceil(8)
            + 32
            + (self.loop_edges.capacity() + self.ext_edges.capacity())
                * std::mem::size_of::<EdgeIndex>()
            + self.edge_signatures.capacity() * std::mem::size_of::<LoopExtSignature>()
            + self
                .edge_signatures
                .iter()
                .map(|(_, signature)| {
                    (signature.internal.len() + signature.external.len())
                        * std::mem::size_of::<SignOrZero>()
                })
                .sum::<usize>()
    }

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
        // The external coordinate is the same selected loop-edge momentum.
        let mut shared: feynkit_graph::routing::MomentumBasis = (&*self).into();
        shared.put_loop_to_ext(i.0);
        *self = shared.into();
    }
    pub(crate) fn canonicalize_external_order(&mut self, external_edge_order: &[EdgeIndex]) {
        let mut shared: feynkit_graph::routing::MomentumBasis = (&*self).into();
        shared.canonicalize_external_order(external_edge_order);
        *self = shared.into();
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
                    FunctionBuilder::new(GS.emr_mom)
                        .add_arg(usize::from(e))
                        .add_args(rep_args)
                        .finish()
                        .to_pattern(),
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
                    FunctionBuilder::new(GS.emr_vec)
                        .add_arg(usize::from(e))
                        .add_args(rep_args)
                        .finish()
                        .to_pattern(),
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
                    FunctionBuilder::new(GS.emr_mom)
                        .add_arg(usize::from(e))
                        .add_args(rep_args)
                        .finish()
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
                    FunctionBuilder::new(GS.emr_mom)
                        .add_arg(usize::from(e))
                        .add_args(rep_args)
                        .finish()
                        .to_pattern(),
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
    /// `shrunken` acts as a contracted passage node. When supplied, the parent
    /// LMB selects a compatible subset of its loop carriers on the contracted
    /// topology instead of independently choosing a canonical basis.
    fn shrunken_sub_lmb(
        &self,
        outer: &SuBitGraph,
        shrunken: &InternalSubGraph,
        externals: SuBitGraph,
        parent_lmb: Option<&LoopMomentumBasis>,
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

    /// Fallible form of [`Self::compatible_sub_lmb`] for callers that must handle
    /// an unavailable parent-compatible basis without panicking. The default
    /// preserves compatibility with external trait implementations that only
    /// implement the original infallible method.
    fn try_compatible_sub_lmb<S: SubGraphLike>(
        &self,
        subgraph: &S,
        externals: S::Base,
        lmb: &LoopMomentumBasis,
    ) -> LmbResult<LoopMomentumBasis>
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

impl From<feynkit_graph::routing::MomentumBasis> for LoopMomentumBasis {
    fn from(basis: feynkit_graph::routing::MomentumBasis) -> Self {
        Self {
            tree: basis.tree,
            loop_edges: basis.loop_edges.into(),
            ext_edges: basis.ext_edges.into(),
            edge_signatures: basis
                .edge_signatures
                .iter()
                .map(|(edge, signature)| {
                    (
                        edge,
                        LoopExtSignature::from(signature.integer_coefficients()),
                    )
                })
                .collect(),
        }
    }
}

impl From<&LoopMomentumBasis> for feynkit_graph::routing::MomentumBasis {
    fn from(basis: &LoopMomentumBasis) -> Self {
        Self {
            tree: basis.tree.clone(),
            loop_edges: basis.loop_edges.raw.clone(),
            ext_edges: basis.ext_edges.raw.clone(),
            edge_signatures: basis
                .edge_signatures
                .iter()
                .map(|(edge, signature)| {
                    (
                        edge,
                        feynkit_kinematics::MomentumSignature::new(
                            feynkit_kinematics::Signature::new(signature.internal.iter().copied()),
                            feynkit_kinematics::Signature::new(signature.external.iter().copied()),
                        ),
                    )
                })
                .collect(),
        }
    }
}

impl<E, V, H> LMBext for HedgeGraph<E, V, H> {
    fn empty_lmb(&self) -> LoopMomentumBasis {
        feynkit_graph::routing::MomentumRouting::empty_lmb(self).into()
    }
    fn lmb(&self) -> LoopMomentumBasis {
        feynkit_graph::routing::MomentumRouting::lmb(self).into()
    }

    fn shrunken_sub_lmb(
        &self,
        outer: &SuBitGraph,
        shrunken: &InternalSubGraph,
        externals: SuBitGraph,
        parent_lmb: Option<&LoopMomentumBasis>,
    ) -> LmbResult<LoopMomentumBasis> {
        feynkit_graph::routing::MomentumRouting::shrunken_sub_lmb(
            self,
            outer,
            shrunken,
            externals,
            parent_lmb
                .map(feynkit_graph::routing::MomentumBasis::from)
                .as_ref(),
        )
        .map(Into::into)
    }

    fn shrunken_lmb_of(
        &self,
        outer: &SuBitGraph,
        shrunken: &InternalSubGraph,
    ) -> LoopMomentumBasis {
        feynkit_graph::routing::MomentumRouting::shrunken_lmb_of(self, outer, shrunken).into()
    }

    fn dot_lmb_of<S: SubGraphLike>(&self, subgraph: &S, lmb: &LoopMomentumBasis) -> String {
        let reps = self.normal_emr_replacement::<_, Atom>(subgraph, lmb, &[], no_filter);

        let emrgraph = self.map_data_ref(
            |_, _, _| "",
            |_, e, _, _| {
                EdgeData::new(
                    GS.emr_mom
                        .call_args([usize::from(e)])
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
        feynkit_graph::routing::MomentumRouting::lmb_of(self, subgraph).into()
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
        feynkit_graph::routing::MomentumRouting::compatible_sub_lmb(
            self,
            subgraph,
            externals,
            &lmb.into(),
        )
        .into()
    }

    fn try_compatible_sub_lmb<S: SubGraphLike>(
        &self,
        subgraph: &S,
        externals: S::Base,
        lmb: &LoopMomentumBasis,
    ) -> LmbResult<LoopMomentumBasis>
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>,
    {
        feynkit_graph::routing::MomentumRouting::try_compatible_sub_lmb(
            self,
            subgraph,
            externals,
            &lmb.into(),
        )
        .map(Into::into)
    }

    /// The true externals (that will flow through the graph (i.e. not dummy)) are those that are both in the subgraph and in the externals
    fn lmb_impl<S: SubGraphLike + SubSetOps + ModifySubSet<HedgePair> + ModifySubSet<Hedge>>(
        &self,
        subgraph: &S,
        forest_guide: &S, //guide for the forest (can be the full subgraph if no guide necessary), however it must cover the same nodes as subgraph
        externals: S,     //externals to consider for the flow, cannot contain non-subgraph nodes
    ) -> LmbResult<LoopMomentumBasis>
    where
        S::Base: ModifySubSet<Hedge> + SubGraphLike,
    {
        feynkit_graph::routing::MomentumRouting::lmb_impl(self, subgraph, forest_guide, externals)
            .map(Into::into)
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
        feynkit_graph::routing::MomentumRouting::generate_loop_momentum_bases_of(self, subgraph)
            .into_iter()
            .map(Into::into)
            .collect()
    }

    fn generate_loop_momentum_bases(&self) -> TiVec<LmbIndex, LoopMomentumBasis> {
        feynkit_graph::routing::MomentumRouting::generate_loop_momentum_bases(self)
            .into_iter()
            .map(Into::into)
            .collect()
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
        feynkit_graph::routing::MomentumRouting::replacement_impl(
            self,
            rep,
            subgraph,
            &lmb.into(),
            loop_symbol,
            ext_symbol,
            loop_args,
            ext_args,
            filter_pair,
            emr_id,
        )
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
        parent_lmb: Option<&LoopMomentumBasis>,
    ) -> LmbResult<LoopMomentumBasis> {
        let mut lmb = self
            .underlying
            .shrunken_sub_lmb(outer, shrunken, externals, parent_lmb)?;
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
        self.try_compatible_sub_lmb(subgraph, externals, lmb)
            .unwrap_or_else(|err| {
                panic!("Failed to build compatible subgraph loop momentum basis:\n{err}")
            })
    }

    fn try_compatible_sub_lmb<S: SubGraphLike>(
        &self,
        subgraph: &S,
        externals: S::Base,
        lmb: &LoopMomentumBasis,
    ) -> LmbResult<LoopMomentumBasis>
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>,
    {
        let mut sub_lmb = self
            .underlying
            .try_compatible_sub_lmb(subgraph, externals, lmb)?;
        self.canonicalize_lmb_external_order(&mut sub_lmb);
        Ok(sub_lmb)
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
        parent_lmb: Option<&LoopMomentumBasis>,
    ) -> LmbResult<LoopMomentumBasis> {
        let mut lmb = self
            .underlying
            .shrunken_sub_lmb(outer, shrunken, externals, parent_lmb)?;
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
        self.try_compatible_sub_lmb(subgraph, externals, lmb)
            .unwrap_or_else(|err| {
                panic!("Failed to build compatible subgraph loop momentum basis:\n{err}")
            })
    }

    fn try_compatible_sub_lmb<S: SubGraphLike>(
        &self,
        subgraph: &S,
        externals: S::Base,
        lmb: &LoopMomentumBasis,
    ) -> LmbResult<LoopMomentumBasis>
    where
        S::Base: SubGraphLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>,
    {
        let mut sub_lmb = self
            .underlying
            .try_compatible_sub_lmb(subgraph, externals, lmb)?;
        self.canonicalize_lmb_external_order(&mut sub_lmb);
        Ok(sub_lmb)
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
                    - selfmom.call_args([l.0]),
            )
        }

        let mut vars = vec![];

        for (l, _) in other.loop_edges.iter_enumerated() {
            vars.push(othermom.call_args([l.0]))
        }

        let solutions = Atom::solve(&sys).wrt_with_exponent::<u8, _>(&vars).unwrap();
        let [solution] = solutions.iter().as_slice() else {
            panic!(
                "expected one loop-momentum basis solution, got {} branches",
                solutions.len()
            );
        };
        assert!(
            solution.free_variables().is_empty(),
            "loop-momentum basis solution is underdetermined"
        );
        vars.iter()
            .map(|variable| {
                let variable = PolyVariable::try_from(variable.clone()).unwrap();
                solution.get(&variable).cloned().unwrap()
            })
            .collect()
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
        finalized_runtime_dot,
        graph::{FeynmanGraph, Graph, LMBext, LmbError, parse::IntoFinalizedRuntimeGraph},
        initialisation::test_initialise,
        momentum::SignOrZero,
    };

    static SHRUNKEN_LMB_TEST_INIT: std::sync::Once = std::sync::Once::new();
    static SHRUNKEN_LMB_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

    #[test]
    fn generated_singleton_lmbs_conserve_external_momentum() {
        test_initialise().unwrap();
        let graphs: Vec<Graph> = dot!(
            digraph contact {
                ext [style=invis]
                node[num=1]
                ext->v:0[id=0]
                ext->v:1[id=1]
                v:2->ext[id=2]
                v:3->ext[id=3]
            }
            digraph tadpole {
                ext [style=invis]
                node[num=1]
                ext->v:0[id=0]
                ext->v:1[id=1]
                v:2->ext[id=2]
                v:3->ext[id=3]
                v->v[id=4 num=1 mass=1]
            }
        )
        .unwrap();

        for (loops, graph) in graphs.iter().enumerate() {
            let lmbs = graph.generate_loop_momentum_bases();
            let lmb = lmbs
                .first()
                .expect("a singleton graph has a valid momentum basis");
            assert_eq!(lmb.loop_edges.len(), loops);
            // The outgoing dependent leg carries p3 = p0 + p1 - p2.
            for (edge, expected) in [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [1, 1, -1, 0]]
                .into_iter()
                .enumerate()
            {
                let signature = &lmb.edge_signatures[EdgeIndex::from(edge)];
                assert!(
                    signature
                        .internal
                        .iter()
                        .all(|sign| *sign == SignOrZero::Zero)
                );
                assert_eq!(
                    signature
                        .external
                        .iter()
                        .map(|sign| *sign * 1_i32)
                        .collect::<Vec<_>>(),
                    expected
                );
            }
            if loops == 1 {
                let signature = &lmb.edge_signatures[EdgeIndex::from(4)];
                assert_eq!(
                    signature.internal.iter().copied().collect::<Vec<_>>(),
                    [SignOrZero::Plus]
                );
                assert!(
                    signature
                        .external
                        .iter()
                        .all(|sign| *sign == SignOrZero::Zero)
                );
            }
        }
    }

    #[test]
    fn lmb_for_dummy() {
        test_initialise().unwrap();
        let gs: Vec<Graph> = finalized_runtime_dot!(
            digraph dxda{
                graph [projector="1"]
                ext [style=invis]
                node[num=1]
                edge[num=1 particle=a dir=none]
                ext->v1:0[id=0 is_dummy=true sink="{ufo_order:0}"]
                ext->v1:1[id=1 sink="{ufo_order:1}"]
                ext->v1:2[id=2 sink="{ufo_order:2}"]
            }
            digraph aa{
                graph [projector="1"]
                ext [style=invis]
                node[num=1]
                edge[num=1 particle=a dir=none]
                ext->v1:0[id=0 is_dummy=true sink="{ufo_order:0}"]
                ext->v1:1[id=1 is_dummy=true sink="{ufo_order:1}"]
                v1->v2[id=3 source="{ufo_order:2}" sink="{ufo_order:0}"]
                ext->v2:2[id=2 sink="{ufo_order:1}"]
            }
        )
        .unwrap();

        for g in gs {
            insta::with_settings!({
                snapshot_suffix=>g.name.to_string(),
            }, {
                insta::assert_snapshot!(g.dot_lmb_of(&g.full_filter(), &g.loop_momentum_basis));
            });
        }
    }

    #[test]
    fn generated_lmbs_do_not_use_dummy_external_carriers() {
        test_initialise().unwrap();
        let g: Graph = finalized_runtime_dot!(digraph{
            graph [projector="1"]
            ext [style=invis]
            edge[num=1 mass=1 particle=a dir=none]
            node[num=1]
            ext->v1:0[id=0 is_dummy=true sink="{ufo_order:0}"]
            ext->v1:1[id=1 sink="{ufo_order:1}"]
            v1->v2[id=2 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:0}"]
            v2->v1[id=3 source="{ufo_order:1}" sink="{ufo_order:3}"]
            ext->v2:2[id=4 sink="{ufo_order:2}"]
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
        let g: DotGraph = linnet::dot!(digraph{

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
        let lmb = g.lmb();
        assert_snapshot!(g.dot_lmb_of(&g.full_filter(), &lmb));
        assert_snapshot!(&lmb.to_string());
        let _g = g.generate_loop_momentum_bases_of(&g.full_filter());
    }
    #[test]
    fn disconnected() {
        test_initialise().unwrap();
        let g: DotGraph = linnet::dot!(digraph{
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
        let lmb = g.lmb();
        assert_snapshot!(g.dot_lmb_of(&g.full_filter(), &lmb));
        assert_snapshot!(&lmb.to_string());

        let g: DotGraph = linnet::dot!(digraph{

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
        let g: DotGraph = linnet::dot!(digraph{
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
        test_initialise().unwrap();
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

        let mut incompatible_parent_lmb = lmb.clone();
        incompatible_parent_lmb.loop_edges.clear();
        let unavailable =
            g.try_compatible_sub_lmb(&subgraph, g.full_crown(&subgraph), &incompatible_parent_lmb);
        assert!(matches!(
            unavailable,
            Err(LmbError::NoCompatibleSubLmb { .. })
        ));

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

        let sub_lmb = g
            .try_compatible_sub_lmb(&subgraph, non_dummy_sub_ext, &lmb)
            .unwrap();

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
        let g: DotGraph = linnet::dot!(digraph{
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
        let shrunken = InternalSubGraph::try_new(shrunken_filter, &g).expect("valid subgraph");
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
        let g: DotGraph = linnet::dot!(digraph{
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
        let shrunken = InternalSubGraph::try_new(shrunken_filter, &g).expect("valid subgraph");
        let remainder = outer.subtract(&shrunken.filter);

        let lmb = g.shrunken_lmb_of(&outer, &shrunken);

        assert!(
            !lmb.loop_edges
                .iter()
                .any(|edge| shrunken.filter.includes(&g[edge].1))
        );
        // Contracting one edge in each triangle preserves two independent loop flows.
        // On the quotient graph, conservation equates each retained edge pair and
        // leaves the bridge with zero momentum, for any valid choice of loop basis.
        let first = &lmb.edge_signatures[EdgeIndex(1)].internal;
        let second = &lmb.edge_signatures[EdgeIndex(4)].internal;
        assert_eq!(first, &lmb.edge_signatures[EdgeIndex(2)].internal);
        assert_eq!(second, &lmb.edge_signatures[EdgeIndex(5)].internal);
        assert!(
            lmb.edge_signatures[EdgeIndex(6)]
                .internal
                .iter()
                .all(|sign| *sign == SignOrZero::Zero)
        );
        for (_, edge, _) in g.iter_edges_of(&remainder) {
            assert!(
                lmb.edge_signatures[edge]
                    .external
                    .iter()
                    .all(|sign| *sign == SignOrZero::Zero)
            );
        }
        let first_norm: i64 = first.iter().map(|sign| (*sign * 1_i64).pow(2)).sum();
        let second_norm: i64 = second.iter().map(|sign| (*sign * 1_i64).pow(2)).sum();
        let overlap: i64 = first
            .iter()
            .zip(second.iter())
            .map(|(left, right)| (*left * 1_i64) * (*right * 1_i64))
            .sum();
        assert!(
            first_norm * second_norm > overlap.pow(2),
            "both surviving cycle flows must be independent"
        );
    }

    #[test]
    fn shrunken_edge_outside_outer_errors() {
        let _guard = SHRUNKEN_LMB_TEST_LOCK.lock().unwrap();
        SHRUNKEN_LMB_TEST_INIT.call_once(|| test_initialise().unwrap());
        let g: DotGraph = linnet::dot!(digraph{
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
        let shrunken = InternalSubGraph::try_new(shrunken_filter, &g).expect("valid subgraph");
        let outer = g.full_filter().subtract(&shrunken.filter);

        let result = g.shrunken_sub_lmb(&outer, &shrunken, g.full_crown(&outer), None);

        assert!(matches!(result, Err(LmbError::ShrunkenOutsideOuter { .. })));
    }

    #[test]
    fn shrunken_whole_component_retains_only_surviving_external_endpoints() {
        use linnet::half_edge::subgraph::SubSetLike;

        let _guard = SHRUNKEN_LMB_TEST_LOCK.lock().unwrap();
        SHRUNKEN_LMB_TEST_INIT.call_once(|| test_initialise().unwrap());
        let g: Graph = dot!(digraph {
            edge[num=1 mass=1]
            node[num=1]
            b:17 -> a:0 [id=0]
            a:1 -> e:2 [id=1]
            a:3 -> f:4 [id=2]
            b:5 -> c:6 [id=3]
            b:7 -> d:8 [id=4]
            c:9 -> d:10 [id=5]
            c:11 -> f:12 [id=6]
            d:13 -> e:14 [id=7]
            e:15 -> f:16 [id=8]
            x:18 -> y:19 [id=9]
        })
        .unwrap();
        let mut outer: SuBitGraph = g.empty_subgraph();
        let mut shrunken_filter: SuBitGraph = g.empty_subgraph();
        for edge in [1, 2, 8].map(EdgeIndex::from) {
            shrunken_filter.add(g[&edge].1);
        }
        outer.union_with(&shrunken_filter);
        for edge in [3, 4, 5].map(EdgeIndex::from) {
            outer.add(g[&edge].1);
        }
        let shrunken = InternalSubGraph::try_new(shrunken_filter, &g.underlying).unwrap();
        let remainder = outer.subtract(&shrunken.filter);
        let externals = g.full_crown(&outer);
        let parent = g.lmb_impl(&outer, &outer, externals.clone()).unwrap();
        assert_eq!(parent.loop_edges.len(), 2);

        let mut contracted = g.underlying.to_ref();
        let root = shrunken.included_iter().next().unwrap();
        contracted.identify_nodes_of_subgraph_without_self_edges::<_, SuBitGraph>(
            &shrunken,
            &g.underlying[g.node_id(root)],
        );
        contracted.forget_identification_history();
        let surviving_crown = contracted.full_crown(&remainder);
        assert_eq!(
            surviving_crown.included_iter().collect::<Vec<_>>(),
            [Hedge(11), Hedge(13), Hedge(17)]
        );
        assert!(matches!(
            contracted.lmb_impl(&remainder, &remainder, externals.clone()),
            Err(LmbError::ExternalsOutsideSubgraph { .. })
        ));

        for parent_lmb in [None, Some(&parent)] {
            let lmb = g
                .shrunken_sub_lmb(&outer, &shrunken, externals.clone(), parent_lmb)
                .unwrap();
            assert_eq!(lmb.loop_edges.len(), 1);
            assert!(
                remainder.includes(&g[&lmb.loop_edges[crate::momentum::sample::LoopIndex(0)]].1)
            );
            let mut external_edges = lmb.ext_edges.iter().copied().collect::<Vec<_>>();
            external_edges.sort();
            assert_eq!(external_edges, [EdgeIndex(0), EdgeIndex(6), EdgeIndex(7)]);
            // The surviving triangle conserves its loop flow at all three
            // vertices; each connector retains exactly its physical endpoint.
            let loop_rows = [3, 4, 5].map(|edge| {
                *lmb.edge_signatures[EdgeIndex(edge)]
                    .internal
                    .iter()
                    .next()
                    .unwrap()
                    * 1_i32
            });
            assert_eq!(loop_rows[0] + loop_rows[1], 0);
            assert_eq!(-loop_rows[0] + loop_rows[2], 0);
            assert_eq!(-loop_rows[1] - loop_rows[2], 0);
        }

        // External-flow carriers can also be internal remainder hedges. After
        // contracting only e1, hedge 3 belongs to the old contracted crown but
        // remains on e2 at a surviving node, so it must retain its flow slot.
        let mut partial_filter: SuBitGraph = g.empty_subgraph();
        partial_filter.add(g[&EdgeIndex(1)].1);
        let partial = InternalSubGraph::try_new(partial_filter, &g.underlying).unwrap();
        let mut internal_externals = externals.clone();
        internal_externals.add(Hedge(3));
        // Both parent carriers must survive this contraction; the earlier
        // parent used e1, which is now contracted despite its loop remaining.
        let mut partial_parent_guide = outer.clone();
        for edge in [2, 3].map(EdgeIndex::from) {
            partial_parent_guide.sub(g[&edge].1);
        }
        let partial_parent = g
            .lmb_impl(&outer, &partial_parent_guide, internal_externals.clone())
            .unwrap();
        let mut partial_carriers = partial_parent
            .loop_edges
            .iter()
            .copied()
            .collect::<Vec<_>>();
        partial_carriers.sort();
        assert_eq!(partial_carriers, [EdgeIndex(2), EdgeIndex(3)]);
        for parent_lmb in [None, Some(&partial_parent)] {
            let lmb = g
                .shrunken_sub_lmb(&outer, &partial, internal_externals.clone(), parent_lmb)
                .unwrap();
            assert_eq!(lmb.loop_edges.len(), 2);
            assert!(lmb.ext_edges.contains(&EdgeIndex(2)));
        }

        // A crown intersection must not silently discard an unrelated invalid
        // external node alongside the legitimately retired component.
        let mut invalid_externals = externals;
        invalid_externals.add(Hedge(18));
        let invalid = g.shrunken_sub_lmb(&outer, &shrunken, invalid_externals, None);
        assert!(matches!(
            invalid,
            Err(LmbError::ExternalsOutsideSubgraph { .. })
        ));
    }
}
