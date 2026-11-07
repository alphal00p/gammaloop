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
        SubSetLike, SubSetOps,
    },
    tree::SimpleTraversalTree,
};
use linnet::union_find::UnionFind;
use std::mem;
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
    #[error("split edge on full graph")]
    SplitEdgeOnFullGraph,
    /// Unreachable for a graph: identifying bridge columns pivots on a network
    /// matrix, whose entries and eliminations stay within -1, 0 and 1.
    #[error(
        "routing bridge {edge:?} between selected components needs a coefficient outside -1, 0, 1"
    )]
    BridgeCoefficientOverflow { edge: EdgeIndex },
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

/// Give each bridge between two selected components one external column.
///
/// Each component routes its crown on its own, so a bridge enters `ext_edges`
/// once per component, and in each its momentum is either that column or,
/// where it is the dependent external, a combination of the component's other
/// columns. Equating the two expressions of every bridge gives a homogeneous
/// linear system over the columns, which always has a solution. Eliminating it
/// leaves one independent momentum per non-bridge external that is not
/// dependent and per independent cycle of the graph of components and
/// bridges.
fn identify_bridge_columns(
    ext_edges: &mut Vec<EdgeIndex>,
    flow_components: &[usize],
    signatures: &mut EdgeVec<MomentumSignature>,
) -> LmbResult<()> {
    let bridges = ext_edges
        .iter()
        .enumerate()
        .filter_map(|(second, edge)| {
            let first = ext_edges[..second].iter().position(|other| other == edge)?;
            Some((*edge, first, second))
        })
        .collect_vec();
    if bridges.is_empty() {
        return Ok(());
    }
    let mut rows = signatures
        .iter()
        .map(|(_, signature)| signature.external.integer_coefficients())
        .collect_vec();
    let overflow = |edge| LmbError::BridgeCoefficientOverflow { edge };
    // Only the component owning a column routes a bridge through it, so the
    // bridge's momentum in that component is its row restricted to them.
    let momentum_in = |row: &[isize], column: usize| {
        row.iter()
            .zip(flow_components)
            .map(|(coefficient, owner)| {
                if *owner == flow_components[column] {
                    *coefficient
                } else {
                    0
                }
            })
            .collect_vec()
    };
    // Each pivot column equals its fully reduced expression in free columns.
    let mut pivots: Vec<(usize, Vec<isize>)> = vec![];
    let substitute = |row: &mut Vec<isize>, pivots: &[(usize, Vec<isize>)]| {
        for (pivot, expression) in pivots {
            let coefficient = mem::take(&mut row[*pivot]);
            for (total, term) in row.iter_mut().zip(expression) {
                *total += coefficient * term;
            }
        }
    };
    for &(edge, first, second) in &bridges {
        let own = &rows[edge.0];
        let mut constraint = momentum_in(own, first);
        for (total, term) in constraint.iter_mut().zip(momentum_in(own, second)) {
            *total -= term;
        }
        substitute(&mut constraint, &pivots);
        if constraint.iter().all(|coefficient| *coefficient == 0) {
            continue;
        }
        // Keep the earliest columns free: a bridge's first column stays its own.
        let pivot = constraint
            .iter()
            .rposition(|coefficient| coefficient.abs() == 1)
            .ok_or_else(|| overflow(edge))?;
        let sign = mem::take(&mut constraint[pivot]);
        let expression = constraint.iter().map(|term| -sign * term).collect_vec();
        for (_, earlier) in &mut pivots {
            let coefficient = mem::take(&mut earlier[pivot]);
            for (total, term) in earlier.iter_mut().zip(&expression) {
                *total += coefficient * term;
            }
        }
        pivots.push((pivot, expression));
    }
    for &(edge, first, _) in &bridges {
        rows[edge.0] = momentum_in(&rows[edge.0], first);
    }
    for row in &mut rows {
        substitute(row, &pivots);
        for &(_, first, second) in &bridges {
            row[first] += mem::take(&mut row[second]);
        }
    }
    let removed = bridges.iter().map(|(_, _, second)| *second).collect_vec();
    for ((edge, signature), row) in signatures.iter_mut().zip(rows) {
        let coefficients = row
            .into_iter()
            .enumerate()
            .filter(|(column, _)| !removed.contains(column))
            .map(|(_, coefficient)| i8::try_from(coefficient).map_err(|_| overflow(edge)))
            .collect::<LmbResult<Vec<_>>>()?;
        signature.external =
            Signature::try_from_integers(coefficients).map_err(|_| overflow(edge))?;
    }
    *ext_edges = ext_edges
        .iter()
        .enumerate()
        .filter(|(column, _)| !removed.contains(column))
        .map(|(_, edge)| *edge)
        .collect();
    Ok(())
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
        parent_lmb: Option<&MomentumBasis>,
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
        parent_lmb: Option<&MomentumBasis>,
    ) -> LmbResult<MomentumBasis> {
        validate_shrunken(self, outer, shrunken, &externals)?;

        if shrunken.is_empty() {
            return match parent_lmb {
                Some(parent_lmb) => self.try_compatible_sub_lmb(outer, externals, parent_lmb),
                None => self.lmb_impl(outer, outer, externals),
            };
        }

        let (contracted, remainder, contracted_externals) =
            contract_shrunken(self, outer, shrunken, externals);
        match parent_lmb {
            Some(parent_lmb) => {
                contracted.try_compatible_sub_lmb(&remainder, contracted_externals, parent_lmb)
            }
            None => contracted.lmb_impl(&remainder, &remainder, contracted_externals),
        }
        .map_err(|source| LmbError::NoShrunkenLmb {
            outer_dot: self.dot(outer),
            shrunken_dot: self.dot(shrunken),
            remainder_dot: self.dot(&remainder),
            source: Box::new(source),
        })
    }

    fn shrunken_lmb_of(&self, outer: &SuBitGraph, shrunken: &InternalSubGraph) -> MomentumBasis {
        let externals = self.full_crown(outer);
        self.shrunken_sub_lmb(outer, shrunken, externals, None)
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
        if self.cyclotomatic_number(subgraph) == 0 {
            return Ok(self.empty_lmb());
        }
        let forest_guide = compatible_forest_guide(self, subgraph, lmb)?;
        self.lmb_impl(subgraph.included(), &forest_guide, externals)
    }

    /// Route the non-dummy external carriers on the selected nodes, including
    /// boundary half-edges absent from the selection. Omitted paired carriers
    /// within one component carry independent soft momenta.
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
        let mut paired_externals = vec![];
        // A bridge between two selected components has a half-edge in each.
        let boundary = externals.clone();
        let mut flow_components: Vec<usize> = vec![];
        let mut component = 0;

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
                let subgraph_cover = subgraph_tree.covers(subgraph);

                // Select the last true boundary half-edge as dependent; omitted chords are
                // independent. A bridge to another selected component keeps its own
                // momentum unless this component has no other boundary.
                let candidates = self
                    .iter_edges_of(&external_cover)
                    .filter_map(|(pair, _, _)| {
                        let hedge = match pair {
                            HedgePair::Unpaired { hedge, .. } => hedge,
                            HedgePair::Split {
                                source,
                                sink,
                                split,
                            } => match split {
                                Flow::Source => source,
                                Flow::Sink => sink,
                            },
                            HedgePair::Paired { source, sink } if subgraph.includes(&pair) => {
                                source.max(sink)
                            }
                            HedgePair::Paired { .. } => return None,
                        };
                        let bridge = matches!(pair, HedgePair::Split { .. })
                            && boundary.includes(&self.inv(hedge));
                        Some((hedge, bridge))
                    })
                    .collect_vec();
                let has_boundary = candidates.iter().any(|(_, bridge)| !bridge);
                let dependent = candidates
                    .iter()
                    .filter(|(_, bridge)| !(has_boundary && *bridge))
                    .map(|(hedge, _)| *hedge)
                    .max();
                root = dependent
                    .or_else(|| external_cover.included_iter().next_back())
                    .ok_or_else(|| LmbError::EmptyExternalCover {
                        externals_dot: self.dot(&externals),
                        subgraph_dot: self.dot(subgraph),
                    })?;
                let root_node = self.node_id(root);
                let tree = if subgraph_tree.tree_subgraph.is_empty() {
                    // A singleton component has an empty spanning forest after external
                    // hairs are stripped. Retain its root, including for contact graphs
                    // and tadpoles, so its external flows still conserve momentum.
                    subgraph_tree
                } else {
                    SimpleTraversalTree::depth_first_traverse(self, forest_guide, &root_node, None)
                        .map_err(|_| LmbError::ForestGuideMismatch {
                            forest_guide_dot: self.dot(forest_guide),
                            subgraph_dot: self.dot(subgraph),
                        })?
                };

                if tree.covers(subgraph) != subgraph_cover {
                    return Err(LmbError::ForestGuideMismatch {
                        forest_guide_dot: self.dot(forest_guide),
                        subgraph_dot: self.dot(subgraph),
                    });
                }

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
                        HedgePair::Paired { source, .. } if !subgraph.includes(&p) => {
                            // An omitted chord has two endpoints in this component.
                            // Its momentum is independent and flows between those
                            // endpoints, rather than to the dependent external root.
                            let signed = tree.get_signed_cycle(source, self).ok_or_else(|| {
                                LmbError::MissingCycle {
                                    hedge: source,
                                    tree_dot: self.dot(&tree.tree_subgraph),
                                }
                            })?;
                            for hedge in signed.filter.included_iter() {
                                path_to_dep.add(hedge);
                            }
                            external_flows.push((SignOrZero::Minus, path_to_dep));
                            paired_externals.push(e);
                            ext_edges.push(e);
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
            flow_components.resize(external_flows.len(), component);
            component += 1;

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
                            let signed_cycle =
                                tree.get_signed_cycle(source, self).ok_or_else(|| {
                                    LmbError::MissingCycle {
                                        hedge: source,
                                        tree_dot: self.dot(&tree.tree_subgraph),
                                    }
                                })?;
                            cycles.push(signed_cycle);
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

        let mut signature = self
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
                                // Boundary carriers still need their own signature when
                                // their half-edges are absent from the selected region.
                                if paired_externals.contains(&eid) {
                                    external = ext_edges
                                        .iter()
                                        .map(|edge| {
                                            if *edge == eid {
                                                SignOrZero::Plus
                                            } else {
                                                SignOrZero::Zero
                                            }
                                        })
                                        .collect();
                                } else if subgraph.intersects(&p) || ext_edges.contains(&eid) {
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
                                if subgraph.includes(&hedge) || ext_edges.contains(&eid) {
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
        identify_bridge_columns(&mut ext_edges, &flow_components, &mut signature)?;

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

/// Remove parent loop edges from `subgraph` until it is a forest guide spanning
/// the same nodes and components, so the induced basis reuses those loops.
fn compatible_forest_guide<E, V, H, S: SubGraphLike>(
    graph: &HedgeGraph<E, V, H>,
    subgraph: &S,
    lmb: &MomentumBasis,
) -> LmbResult<S::Base>
where
    S::Base: SubGraphLike<Base = S::Base>
        + SubSetOps
        + Clone
        + ModifySubSet<HedgePair>
        + ModifySubSet<Hedge>,
{
    let n_loops = graph.cyclotomatic_number(subgraph);
    // the subgraph may have disconnected components in case the of disjoint graphs in a spinney
    let components = graph.count_connected_components(subgraph);

    for v in lmb
        .loop_edges
        .iter()
        .filter(|e| {
            let (_, p) = &graph[*e];
            subgraph.includes(p)
        })
        .combinations(n_loops)
    {
        let mut cut_subgraph = subgraph.included().clone();

        for eid in v {
            let (_, p) = &graph[eid];
            let HedgePair::Paired { source, sink } = p else {
                continue;
            };

            //this is a self-loop
            if graph.node_id(*source) == graph.node_id(*sink) {
                continue;
            }
            cut_subgraph.sub(*p);
        }

        if graph.count_connected_components(&cut_subgraph) == components
            && graph.number_of_nodes_in_subgraph(&cut_subgraph)
                == graph.number_of_nodes_in_subgraph(subgraph)
        {
            return Ok(cut_subgraph);
        }
    }

    let full_graph = graph.full_filter();
    let parent_lmb_has_full_loop_dimension =
        lmb.loop_edges.len() == graph.cyclotomatic_number(&full_graph);
    let (subgraph_dot, parent_lmb_dot) = if parent_lmb_has_full_loop_dimension {
        (
            graph.dot_lmb_of(subgraph, lmb),
            graph.dot_lmb_of(&full_graph, lmb),
        )
    } else {
        // Momentum-label rendering assumes a dimensionally valid parent LMB. Preserve the
        // topology diagnostics without panicking while constructing the fallible error.
        (
            graph.dot(subgraph),
            format!(
                "parent loop edges {:?}; expected {} loops for\n{}",
                lmb.loop_edges,
                graph.cyclotomatic_number(&full_graph),
                graph.dot(&full_graph),
            ),
        )
    };

    Err(LmbError::NoCompatibleSubLmb {
        subgraph_dot,
        parent_lmb_dot,
    })
}

/// Reject a shrunken subgraph, outer region or external carriers that do not
/// fit `graph`.
fn validate_shrunken<E, V, H>(
    graph: &HedgeGraph<E, V, H>,
    outer: &SuBitGraph,
    shrunken: &InternalSubGraph,
    externals: &SuBitGraph,
) -> LmbResult<()> {
    let graph_size = graph.n_hedges();
    let outer_dot = || {
        if outer.size() == graph_size {
            graph.dot(outer)
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
            graph.dot(shrunken)
        } else {
            format!(
                "invalid shrunken size {}, expected {}; label {}",
                shrunken.size(),
                graph_size,
                shrunken.string_label()
            )
        }
    };

    if shrunken.size() != graph_size || !shrunken.valid(graph) {
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

    if !outer.union(&graph.full_crown(outer)).includes(externals) {
        return Err(LmbError::ExternalsOutsideSubgraph {
            externals_dot: graph.dot(externals),
            subgraph_dot: outer_dot(),
        });
    }
    Ok(())
}

/// Contract each connected component of `shrunken` to a passage node, and
/// return the contracted graph with the remainder of `outer` and its external
/// carriers.
fn contract_shrunken<'a, E, V, H>(
    graph: &'a HedgeGraph<E, V, H>,
    outer: &SuBitGraph,
    shrunken: &InternalSubGraph,
    externals: SuBitGraph,
) -> (HedgeGraph<&'a E, &'a V, &'a H>, SuBitGraph, SuBitGraph) {
    let remainder = outer.subtract(&shrunken.filter);
    let mut contracted = graph.to_ref();

    for component in graph.connected_components(shrunken) {
        let Some(root) = component.included_iter().next() else {
            continue;
        };
        let node_data = &graph[graph.node_id(root)];
        contracted
            .identify_nodes_of_subgraph_without_self_edges::<_, SuBitGraph>(&component, node_data);
    }
    contracted.forget_identification_history();

    // A fully contracted component can have no node in the remainder.
    // Remove only its retired crown endpoints; the opposite endpoints at
    // surviving nodes remain external flows. Unrelated invalid externals
    // have already been rejected against the original outer footprint.
    let retired_crown = graph
        .full_crown(shrunken)
        .subtract(&contracted.full_crown(&remainder).union(&remainder));
    let contracted_externals = externals
        .subtract(&shrunken.filter)
        .subtract(&retired_crown);
    (contracted, remainder, contracted_externals)
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
        let routed = basis.ext_edges.clone();
        basis.canonicalize_external_order(&external_order);
        let mut basis = crate::LoopMomentumBasis::from_routing(&self.graph, basis);
        // Diagram externals the region does not touch only add zero columns.
        basis
            .dependent_externals
            .retain(|edge| routed.contains(&EdgeIndex(edge.0)));
        basis
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

    /// Extend a forest of selected propagators by the selection's dangling
    /// hairs. Sewn carriers never enter the forest, so every basis of a region
    /// keeps their momenta external.
    fn routing_forest(&self, selected: &SuBitGraph, mut forest: SuBitGraph) -> SuBitGraph {
        for (pair, _, _) in self.graph.iter_edges_of(selected) {
            if matches!(pair, HedgePair::Unpaired { .. }) {
                forest.add(pair);
            }
        }
        forest
    }

    pub fn momentum_basis_of(
        &self,
        subgraph: &SuBitGraph,
    ) -> Result<crate::LoopMomentumBasis, crate::DiagramError> {
        let selected = subgraph.intersection(&self.momentum_subgraph());
        let forest =
            self.routing_forest(&selected, selected.intersection(&self.internal_subgraph()));
        let externals = self.routing_externals(&selected);
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
                for (pair, _, _) in self.graph.iter_edges_of(&selected) {
                    if matches!(pair, HedgePair::Unpaired { .. }) {
                        forest.add(pair);
                    }
                }
                self.graph
                    .lmb_impl(&selected, &forest, externals)
                    .map(|basis| self.normalize_routing(basis))
                    .map_err(|error| {
                        crate::DiagramError::InvalidLoopMomentumBasis(error.to_string())
                    })
            })
            .collect()
    }

    /// Reuse the parent's loop momenta on the selected propagators, routing
    /// with the same forest constraints as [`Self::momentum_basis_of`].
    pub fn compatible_momentum_basis_of(
        &self,
        subgraph: &SuBitGraph,
        parent: &crate::LoopMomentumBasis,
    ) -> Result<crate::LoopMomentumBasis, crate::DiagramError> {
        let selected = subgraph.intersection(&self.momentum_subgraph());
        let internal = selected.intersection(&self.internal_subgraph());
        compatible_forest_guide(&self.graph, &internal, &parent.to_routing(&self.graph))
            .and_then(|forest| {
                self.graph.lmb_impl(
                    &selected,
                    &self.routing_forest(&selected, forest),
                    self.routing_externals(&selected),
                )
            })
            .map(|basis| self.normalize_routing(basis))
            .map_err(|error| crate::DiagramError::InvalidLoopMomentumBasis(error.to_string()))
    }

    /// Route `subgraph` with each component of `contracted` shrunk to a node,
    /// using the same forest constraints as [`Self::momentum_basis_of`].
    pub fn contracted_momentum_basis_of(
        &self,
        subgraph: &SuBitGraph,
        contracted: &InternalSubGraph,
    ) -> Result<crate::LoopMomentumBasis, crate::DiagramError> {
        let selected = subgraph.intersection(&self.momentum_subgraph());
        let externals = self.routing_externals(&selected);
        validate_shrunken(&self.graph, &selected, contracted, &externals)
            .and_then(|()| {
                let (graph, remainder, externals) =
                    contract_shrunken(&self.graph, &selected, contracted, externals);
                let forest = self.routing_forest(
                    &remainder,
                    remainder.intersection(&self.internal_subgraph()),
                );
                graph.lmb_impl(&remainder, &forest, externals)
            })
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
        let invalid_forest = || {
            crate::DiagramError::InvalidLoopMomentumBasis(format!(
                "requested tree edges {requested:?} do not form a spanning forest"
            ))
        };
        let mut components = UnionFind::new(vec![(); self.graph.n_nodes()]);
        for edge in requested {
            if edge.0 >= self.graph.n_edges() {
                return Err(crate::DiagramError::InvalidLoopMomentumBasis(format!(
                    "unknown tree edge {}",
                    edge.0
                )));
            }
            let pair = self.graph[&EdgeIndex(edge.0)].1;
            if !pair.is_paired() || !internal.includes(&pair) {
                return Err(crate::DiagramError::InvalidLoopMomentumBasis(format!(
                    "tree edge {} is not internal",
                    edge.0
                )));
            }
            let HedgePair::Paired { source, sink } = pair else {
                unreachable!("tree edges are paired");
            };
            let source = Hedge(self.graph.node_id(source).0);
            let sink = Hedge(self.graph.node_id(sink).0);
            if components.find(source) == components.find(sink) {
                return Err(invalid_forest());
            }
            components.union(source, sink, |(), ()| ());
            forest.add(pair);
        }
        // An acyclic subset spans exactly when every internal edge joins nodes
        // already connected by the forest, including singleton components.
        for (pair, _, _) in self.graph.iter_edges_of(&internal) {
            if let HedgePair::Paired { source, sink } = pair
                && components.find(Hedge(self.graph.node_id(source).0))
                    != components.find(Hedge(self.graph.node_id(sink).0))
            {
                return Err(invalid_forest());
            }
        }
        let selected = self.momentum_subgraph();
        let externals = self.routing_externals(&selected);
        for (pair, _, _) in self.graph.iter_edges_of(&selected) {
            if matches!(pair, HedgePair::Unpaired { .. }) {
                forest.add(pair);
            }
        }
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

    /// Route momenta and collect numerical multiples of sums without
    /// distributing products of numerator factors.
    pub fn route_expression(&self, expression: &Atom) -> Atom {
        use symbolica::atom::AtomCore;
        expression
            .replace_multiple(self.momentum_replacements())
            .expand_num()
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

#[cfg(test)]
mod tests {
    use super::*;
    use linnet::half_edge::NodeIndex;
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
    fn omitted_chords_route_between_their_endpoints() {
        let mut builder = HedgeGraphBuilder::new();
        let nodes = (0..4).map(|_| builder.add_node(())).collect::<Vec<_>>();
        for (a, b) in [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 3)] {
            builder.add_edge(nodes[a], nodes[b], (), Orientation::Default);
        }
        let graph: HedgeGraph<(), ()> = builder.build();
        let mut square: SuBitGraph = graph.empty_subgraph();
        for edge in 0..4 {
            square.add(graph[&EdgeIndex(edge)].1);
        }
        let basis = graph.lmb_of(&square);
        assert_eq!(basis.loop_edges.len(), 1);
        assert_eq!(basis.ext_edges, vec![EdgeIndex(4), EdgeIndex(5)]);
        for node in nodes {
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

    /// Net signed momentum leaving `node`, as loop and external coefficients.
    fn outflow(graph: &HedgeGraph<(), ()>, basis: &MomentumBasis, node: NodeIndex) -> Vec<isize> {
        let mut sum = vec![0; basis.loop_edges.len() + basis.ext_edges.len()];
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
        sum
    }

    /// Route the region made of the `bubbles` edges of `graph`, with its full
    /// crown as externals, as `momentum_basis_of` does.
    fn route_bubbles(graph: &HedgeGraph<(), ()>, bubbles: &[EdgeIndex]) -> MomentumBasis {
        let mut selected: SuBitGraph = graph.empty_subgraph();
        for edge in bubbles {
            selected.add(graph[edge].1);
        }
        graph
            .lmb_impl(&selected, &selected, graph.full_crown(&selected))
            .unwrap()
    }

    fn assert_conserved(graph: &HedgeGraph<(), ()>, basis: &MomentumBasis) {
        for node in graph.iter_node_ids() {
            let sum = outflow(graph, basis, node);
            assert!(
                sum.iter().all(|coefficient| *coefficient == 0),
                "momentum is not conserved at {node:?}: {sum:?}\n{}",
                graph.dot_lmb_of(&graph.full_filter(), basis)
            );
        }
    }

    fn external_position(basis: &MomentumBasis, edge: EdgeIndex) -> usize {
        let positions = basis
            .ext_edges
            .iter()
            .positions(|candidate| *candidate == edge)
            .collect_vec();
        assert_eq!(positions.len(), 1, "{edge:?} in {:?}", basis.ext_edges);
        positions[0]
    }

    /// `ext_in -> A`, bubble `A = B`, bridge `B -> C`, bubble `C = D`,
    /// `D -> ext_out`, built in either hedge order. Returns the bridge and
    /// the four bubble edges.
    fn bridged_bubbles(reversed: bool) -> (HedgeGraph<(), ()>, EdgeIndex, Vec<EdgeIndex>) {
        let mut builder = HedgeGraphBuilder::new();
        let mut nodes = (0..4).map(|_| builder.add_node(())).collect_vec();
        if reversed {
            nodes.reverse();
        }
        let [a, b, c, d] = [nodes[0], nodes[1], nodes[2], nodes[3]];
        let mut edges = vec![
            (Some(a), None),
            (Some(a), Some(b)),
            (Some(a), Some(b)),
            (Some(b), Some(c)),
            (Some(c), Some(d)),
            (Some(c), Some(d)),
            (None, Some(d)),
        ];
        if reversed {
            edges.reverse();
        }
        for edge in edges {
            match edge {
                (Some(source), Some(sink)) => {
                    builder.add_edge(source, sink, (), Orientation::Default);
                }
                (Some(node), None) => {
                    builder.add_external_edge(node, (), Orientation::Default, Flow::Sink);
                }
                (None, Some(node)) => {
                    builder.add_external_edge(node, (), Orientation::Default, Flow::Source);
                }
                (None, None) => unreachable!("every edge has an endpoint"),
            }
        }
        let graph: HedgeGraph<(), ()> = builder.build();
        let node_pair = |pair: &HedgePair| match pair {
            HedgePair::Paired { source, sink } => {
                Some((graph.node_id(*source), graph.node_id(*sink)))
            }
            _ => None,
        };
        let bridge = graph
            .iter_edges()
            .find(|(pair, _, _)| node_pair(pair) == Some((b, c)))
            .map(|(_, edge, _)| edge)
            .unwrap();
        let bubbles = graph
            .iter_edges()
            .filter(|(pair, _, _)| {
                matches!(node_pair(pair), Some(ends) if ends == (a, b) || ends == (c, d))
            })
            .map(|(_, edge, _)| edge)
            .collect_vec();
        (graph, bridge, bubbles)
    }

    /// Number of external columns some edge momentum depends on.
    fn independent_momenta(basis: &MomentumBasis) -> usize {
        (0..basis.ext_edges.len())
            .filter(|column| {
                basis
                    .edge_signatures
                    .iter()
                    .any(|(_, signature)| signature.external.get(*column) != Some(SignOrZero::Zero))
            })
            .count()
    }

    /// Double-edge bubbles on a ring of bridges, bubble `i` bridged to bubble
    /// `i + 1` and built in `order`. With `legs`, bubbles `0` and `n / 2` gain
    /// an incoming and an outgoing leg. Returns the graph, its bubble edges and
    /// its bridges.
    fn bubble_ring(
        order: &[usize],
        legs: bool,
    ) -> (HedgeGraph<(), ()>, Vec<EdgeIndex>, Vec<EdgeIndex>) {
        let n = order.len();
        let mut builder = HedgeGraphBuilder::new();
        let mut nodes = vec![[NodeIndex(0); 2]; n];
        for &bubble in order {
            nodes[bubble] = [builder.add_node(()), builder.add_node(())];
        }
        let mut next = 0;
        let mut bubbles = vec![];
        let mut bridges = vec![EdgeIndex(0); n];
        for &bubble in order {
            let [start, end] = nodes[bubble];
            for _ in 0..2 {
                builder.add_edge(start, end, (), Orientation::Default);
                bubbles.push(EdgeIndex(next));
                next += 1;
            }
            builder.add_edge(end, nodes[(bubble + 1) % n][0], (), Orientation::Default);
            bridges[bubble] = EdgeIndex(next);
            next += 1;
            for (leg_bubble, flow) in [(0, Flow::Sink), (n / 2, Flow::Source)] {
                if legs && bubble == leg_bubble {
                    builder.add_external_edge(start, (), Orientation::Default, flow);
                    next += 1;
                }
            }
        }
        (builder.build(), bubbles, bridges)
    }

    #[test]
    fn bubbles_on_a_bridge_ring_route_one_momentum_per_bridge_cycle() {
        let orders = [
            vec![0, 1, 2],
            vec![2, 1, 0],
            vec![1, 2, 0],
            vec![0, 1, 2, 3],
            vec![3, 2, 1, 0],
            vec![1, 3, 0, 2],
            vec![2, 0, 3, 1],
        ];
        for order in orders {
            for legs in [false, true] {
                let context = format!("order={order:?} legs={legs}");
                let (graph, bubbles, bridges) = bubble_ring(&order, legs);
                let basis = route_bubbles(&graph, &bubbles);
                assert_eq!(basis.loop_edges.len(), order.len(), "{context}");
                assert_conserved(&graph, &basis);
                for bridge in &bridges {
                    external_position(&basis, *bridge);
                }
                // The ring is one cycle; the legs add one momentum through it.
                assert_eq!(
                    independent_momenta(&basis),
                    1 + usize::from(legs),
                    "{context}"
                );
                if !legs {
                    let ring = &basis.edge_signatures[bridges[0]];
                    assert!(ring.external.iter().any(|sign| sign != SignOrZero::Zero));
                    for bridge in &bridges {
                        assert_eq!(&basis.edge_signatures[*bridge], ring, "{context}");
                    }
                }
            }
        }
    }

    #[test]
    fn four_selected_bubbles_in_a_ring_have_one_shared_bridge_momentum() {
        let mut builder = HedgeGraphBuilder::new();
        let nodes = (0..8).map(|_| builder.add_node(())).collect::<Vec<_>>();
        for component in 0..4 {
            builder.add_edge(
                nodes[2 * component],
                nodes[2 * component + 1],
                (),
                Orientation::Default,
            );
            builder.add_edge(
                nodes[2 * component],
                nodes[2 * component + 1],
                (),
                Orientation::Default,
            );
        }
        for component in 0..4 {
            builder.add_edge(
                nodes[2 * component + 1],
                nodes[(2 * component + 2) % 8],
                (),
                Orientation::Default,
            );
        }
        let graph: HedgeGraph<(), ()> = builder.build();
        // A closed-form conservation witness: each bubble carries k_i and p-k_i,
        // and all four unselected bridges carry the same independent p.
        let mut witness = vec![vec![0_isize; 5]; 12];
        for component in 0..4 {
            witness[2 * component][component] = 1;
            witness[2 * component + 1][component] = -1;
            witness[2 * component + 1][4] = 1;
            witness[8 + component][4] = 1;
        }
        for node in graph.iter_node_ids() {
            let mut balance = vec![0; 5];
            for hedge in graph.iter_crown(node) {
                let direction = if graph.flow(hedge) == Flow::Source {
                    1
                } else {
                    -1
                };
                for (sum, coefficient) in balance.iter_mut().zip(&witness[graph[&hedge].0]) {
                    *sum += direction * coefficient;
                }
            }
            assert_eq!(balance, [0; 5], "closed-form witness at {node:?}");
        }
        let mut selected: SuBitGraph = graph.empty_subgraph();
        for edge in (0..8).map(EdgeIndex) {
            selected.add(graph[&edge].1);
        }
        let routed = graph.lmb_impl(&selected, &selected, graph.full_crown(&selected));
        assert!(
            routed.is_ok(),
            "valid four-bubble ring was rejected: {routed:?}"
        );
        let basis = routed.unwrap();
        assert_conserved(&graph, &basis);
        for edge in (9..12).map(EdgeIndex) {
            assert_eq!(
                basis.edge_signatures[edge],
                basis.edge_signatures[EdgeIndex(8)]
            );
        }
    }

    #[test]
    fn bridge_between_selected_routing_components_carries_its_own_momentum() {
        for reversed in [false, true] {
            let (graph, bridge, bubbles) = bridged_bubbles(reversed);
            assert_eq!(bubbles.len(), 4);
            let basis = route_bubbles(&graph, &bubbles);
            assert_eq!(basis.loop_edges.len(), 2);
            assert_eq!(basis.ext_edges.len(), 3, "{:?}", basis.ext_edges);
            assert_conserved(&graph, &basis);
            let own = external_position(&basis, bridge);
            let signature = &basis.edge_signatures[bridge];
            assert!(signature.loops.iter().all(|sign| sign == SignOrZero::Zero));
            assert_eq!(
                signature.external.integer_coefficients(),
                (0..3)
                    .map(|position| isize::from(position == own))
                    .collect_vec(),
                "reversed={reversed}"
            );
        }
    }

    #[test]
    fn bridge_to_an_isolated_routing_component_carries_no_momentum() {
        // ext_in -> A = B -> ext_out, with a bubble T = U hung from B by a bridge.
        let mut builder = HedgeGraphBuilder::new();
        let [a, b, t, u] = [(); 4].map(|()| builder.add_node(()));
        builder.add_external_edge(a, (), Orientation::Default, Flow::Sink);
        builder.add_edge(a, b, (), Orientation::Default);
        builder.add_edge(a, b, (), Orientation::Default);
        builder.add_external_edge(b, (), Orientation::Default, Flow::Source);
        builder.add_edge(b, t, (), Orientation::Default);
        builder.add_edge(t, u, (), Orientation::Default);
        builder.add_edge(t, u, (), Orientation::Default);
        let graph: HedgeGraph<(), ()> = builder.build();
        let bridge = EdgeIndex(4);
        let bubbles = [1, 2, 5, 6].map(EdgeIndex);
        let basis = route_bubbles(&graph, &bubbles);
        assert_eq!(basis.loop_edges.len(), 2);
        assert_conserved(&graph, &basis);
        external_position(&basis, bridge);
        let (loops, external) = basis.edge_signatures[bridge].integer_coefficients();
        assert!(loops.into_iter().chain(external).all(|c| c == 0));
    }

    #[test]
    fn bridges_closing_a_loop_between_routing_components_share_one_momentum() {
        // Vacuum graph: bubbles A = B and C = D joined by A -> C and D -> B.
        let mut builder = HedgeGraphBuilder::new();
        let [a, b, c, d] = [(); 4].map(|()| builder.add_node(()));
        builder.add_edge(a, b, (), Orientation::Default);
        builder.add_edge(a, b, (), Orientation::Default);
        builder.add_edge(a, c, (), Orientation::Default);
        builder.add_edge(c, d, (), Orientation::Default);
        builder.add_edge(c, d, (), Orientation::Default);
        builder.add_edge(d, b, (), Orientation::Default);
        let graph: HedgeGraph<(), ()> = builder.build();
        let basis = route_bubbles(&graph, &[0, 1, 3, 4].map(EdgeIndex));
        assert_eq!(basis.loop_edges.len(), 2);
        assert_eq!(basis.ext_edges.len(), 2, "{:?}", basis.ext_edges);
        assert_conserved(&graph, &basis);
        let first = &basis.edge_signatures[EdgeIndex(2)];
        assert_eq!(first, &basis.edge_signatures[EdgeIndex(5)]);
        assert!(first.external.iter().any(|sign| sign != SignOrZero::Zero));
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
