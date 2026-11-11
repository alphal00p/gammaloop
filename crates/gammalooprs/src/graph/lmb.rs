use std::fmt::Display;

use bincode_trait_derive::{Decode, Encode};
use derive_more::{From, Into};
use itertools::Itertools;
use linnet::half_edge::{
    HedgeGraph, NoData,
    involution::{EdgeData, EdgeIndex, EdgeVec, Hedge, HedgePair, Orientation},
    subgraph::{
        Inclusion, InternalSubGraph, ModifySubSet, SuBitGraph, SubGraphLike, SubGraphOps,
        SubSetLike, SubSetOps,
    },
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
    #[error(transparent)]
    Routing(#[from] feynkit_graph::routing::LmbError),
    #[error(
        "no loop-momentum basis compatible with the parent basis was found for subgraph\n{subgraph_dot}\nparent basis\n{parent_lmb_dot}"
    )]
    NoCompatibleSubLmb {
        subgraph_dot: String,
        parent_lmb_dot: String,
    },
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

    fn lmb_impl<S: SubGraphLike + SubSetOps + ModifySubSet<HedgePair> + ModifySubSet<Hedge>>(
        &self,
        subgraph: &S,
        forest_guide: &S,
        externals: S,
    ) -> LmbResult<LoopMomentumBasis>
    where
        S::Base: ModifySubSet<Hedge> + SubGraphLike,
    {
        // Use the shared incidence routing: an omitted chord carries an independent
        // momentum between both endpoints, including when it closes on a UV child.
        let basis = feynkit_graph::routing::MomentumRouting::lmb_impl(
            self,
            subgraph,
            forest_guide,
            externals,
        )?;
        Ok(basis.into())
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

        // Routing coefficients are exact integer graph incidences. The native
        // matrix solve retains the unique-basis requirement; singular systems
        // remain errors rather than choosing values for free loop coordinates.
        let (matrix, rhs) = Atom::system_to_matrix::<u8, _, _>(&sys, &vars).unwrap();
        matrix
            .solve(&rhs)
            .unwrap()
            .into_vec()
            .into_iter()
            .map(|value| value.to_expression())
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
            HedgeGraph, NodeIndex,
            builder::HedgeGraphBuilder,
            involution::{EdgeIndex, Flow, Hedge, Orientation},
            subgraph::{Inclusion, InternalSubGraph, ModifySubSet, SuBitGraph, SubSetOps},
        },
        parser::DotGraph,
    };

    use crate::{
        finalized_runtime_dot,
        graph::{
            FeynmanGraph, Graph, LMBext, LmbError, LoopMomentumBasis,
            parse::IntoFinalizedRuntimeGraph,
        },
        initialisation::test_initialise,
        momentum::SignOrZero,
    };

    static SHRUNKEN_LMB_TEST_INIT: std::sync::Once = std::sync::Once::new();
    static SHRUNKEN_LMB_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

    fn assert_routing_conserved(
        graph: &HedgeGraph<(), ()>,
        basis: &LoopMomentumBasis,
        nodes: &[NodeIndex],
    ) {
        for &node in nodes {
            let mut sum = vec![0; basis.loop_edges.len() + basis.ext_edges.len()];
            for hedge in graph.iter_crown(node) {
                let incidence = if graph.flow(hedge) == Flow::Source {
                    1
                } else {
                    -1
                };
                let signature = &basis.edge_signatures[graph[&hedge]];
                for (total, coefficient) in sum
                    .iter_mut()
                    .zip(signature.internal.iter().chain(signature.external.iter()))
                {
                    let coefficient = match coefficient {
                        SignOrZero::Plus => 1,
                        SignOrZero::Minus => -1,
                        SignOrZero::Zero => 0,
                    };
                    *total += incidence * coefficient;
                }
            }
            assert!(
                sum.iter().all(|coefficient| *coefficient == 0),
                "momentum is not conserved at {node:?}: {sum:?}; {basis}"
            );
        }
    }

    #[test]
    fn omitted_chords_conserve_momentum_for_every_orientation() {
        for orientations in 0..8 {
            let mut builder = HedgeGraphBuilder::new();
            let nodes = (0..4).map(|_| builder.add_node(())).collect::<Vec<_>>();
            for (a, b) in [(0, 1), (1, 2), (2, 3), (3, 0)] {
                builder.add_edge(nodes[a], nodes[b], (), Orientation::Default);
            }
            // Two boundary half-edges at A and C belong to different omitted
            // chords. Each chord transports momentum between both endpoints.
            for (index, (mut a, mut b)) in [(0, 2), (1, 3), (0, 2)].into_iter().enumerate() {
                if orientations & (1 << index) != 0 {
                    std::mem::swap(&mut a, &mut b);
                }
                builder.add_edge(nodes[a], nodes[b], (), Orientation::Default);
            }
            let graph: HedgeGraph<(), ()> = builder.build();
            let mut selected: SuBitGraph = graph.empty_subgraph();
            for edge in 0..4 {
                selected.add(graph[&EdgeIndex(edge)].1);
            }
            let bases = graph.generate_loop_momentum_bases_of(&selected);
            assert_eq!(bases.len(), 4);
            for basis in bases {
                assert_eq!(basis.loop_edges.len(), 1);
                assert_eq!(
                    basis.ext_edges.raw,
                    vec![EdgeIndex(4), EdgeIndex(5), EdgeIndex(6)]
                );
                for (slot, edge) in basis.ext_edges.iter_enumerated() {
                    let signature = &basis.edge_signatures[*edge];
                    assert_eq!(signature.external[slot], SignOrZero::Plus);
                    assert!(
                        signature
                            .internal
                            .iter()
                            .all(|sign| *sign == SignOrZero::Zero)
                    );
                }
                assert_routing_conserved(&graph, &basis, &nodes);
            }
        }
    }

    #[test]
    fn boundary_self_loop_does_not_inject_momentum_into_child() {
        for self_loop_node in 0..4 {
            let mut builder = HedgeGraphBuilder::new();
            let nodes = (0..4).map(|_| builder.add_node(())).collect::<Vec<_>>();
            for (a, b) in [(0, 1), (1, 2), (2, 3), (3, 0)] {
                builder.add_edge(nodes[a], nodes[b], (), Orientation::Default);
            }
            builder.add_edge(
                nodes[self_loop_node],
                nodes[self_loop_node],
                (),
                Orientation::Default,
            );
            builder.add_external_edge(nodes[1], (), Orientation::Default, Flow::Sink);
            builder.add_external_edge(nodes[3], (), Orientation::Default, Flow::Source);
            let graph: HedgeGraph<(), ()> = builder.build();
            let mut selected: SuBitGraph = graph.empty_subgraph();
            for edge in 0..4 {
                selected.add(graph[&EdgeIndex(edge)].1);
            }
            let bases = graph.generate_loop_momentum_bases_of(&selected);
            assert_eq!(bases.len(), 4);
            for basis in bases {
                let self_loop_slot = basis.ext_from(EdgeIndex(4)).unwrap();
                assert_eq!(
                    basis.edge_signatures[EdgeIndex(4)].external[self_loop_slot],
                    SignOrZero::Plus
                );
                for edge in 0..4 {
                    assert_eq!(
                        basis.edge_signatures[EdgeIndex(edge)].external[self_loop_slot],
                        SignOrZero::Zero,
                        "a boundary self-loop cannot supply momentum to the child"
                    );
                }
                assert_routing_conserved(&graph, &basis, &nodes);
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
        // Contracting the disjoint forest edges (a,b) and (d,e) preserves
        // both cycles. The quotient has vertices ab,c,de,f and directed edges
        // ab->c, c->ab, de->f, f->de, c->de. Check its incidence equations
        // independently of the DOT printer and of the selected loop carriers.
        assert_eq!(g.cyclotomatic_number(&outer), 2);
        assert_eq!(lmb.loop_edges.len(), 2);
        assert!(lmb.ext_edges.is_empty());
        let coefficients = [1, 2, 4, 5, 6].map(|edge| {
            let signature = &lmb.edge_signatures[EdgeIndex::from(edge)];
            assert!(
                signature
                    .external
                    .iter()
                    .all(|sign| *sign == SignOrZero::Zero)
            );
            signature
                .internal
                .iter()
                .map(|sign| match sign {
                    SignOrZero::Plus => 1_i8,
                    SignOrZero::Minus => -1,
                    SignOrZero::Zero => 0,
                })
                .collect::<Vec<_>>()
        });
        let [first_out, first_back, second_out, second_back, bridge] = coefficients;
        for index in 0..2 {
            assert_eq!(-first_out[index] + first_back[index], 0, "vertex ab");
            assert_eq!(
                first_out[index] - first_back[index] - bridge[index],
                0,
                "vertex c"
            );
            assert_eq!(
                -second_out[index] + second_back[index] + bridge[index],
                0,
                "vertex de"
            );
            assert_eq!(second_out[index] - second_back[index], 0, "vertex f");
        }
        assert_ne!(
            first_out[0] * second_out[1] - first_out[1] * second_out[0],
            0,
            "the two surviving cycles need independent loop momenta"
        );
        assert_snapshot!(g.dot_lmb_of(&remainder, &lmb));
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

        let result = g.shrunken_sub_lmb(&outer, &shrunken, g.full_crown(&outer));

        assert!(matches!(result, Err(LmbError::ShrunkenOutsideOuter { .. })));
    }
}
