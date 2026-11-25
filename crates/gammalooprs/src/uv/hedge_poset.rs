use std::time::{Duration, Instant};
use std::{
    collections::BTreeMap,
    fmt::{self, Display},
    sync::Arc,
};

#[cfg(test)]
use crate::uv::Integrands;
#[cfg(test)]
use crate::{
    cff::expression::OrientationID,
    uv::approx::direct_3d::{DirectCoordinateFrame, DirectResidueBranches, DirectSector},
};
#[cfg(test)]
use std::cmp::Reverse;

use ahash::AHashMap;
use eyre::{WrapErr, eyre};
use feynkit_uv::{
    ForestNode as UvForestNode, Forests as UvForests, Operation,
    OwnedForestNode as UvOwnedForestNode, Wood as UvWood,
};
use idenso::{
    CookMode, CookSettings,
    color::ColorSimplifySettings,
    representations::{ColorAdjoint, ColorFundamental, ColorSextet},
    tensor::SymbolicTensor,
};
use itertools::Itertools;
use linnet::half_edge::{
    HedgeGraph, NodeIndex,
    algorithms::{DirectionBasis, trace_unfold::UnfoldedTraceGraph},
    involution::HedgePair,
    nodestore::NodeStorageVec,
    subgraph::{ModifySubSet, SuBitGraph, SubSetLike, SubSetOps},
};
use spenso::shadowing::TensorCollectFilter;
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder},
    function,
};
use tracing::debug;
use vakint::{Vakint, VakintSettings};

pub use feynkit_uv::OperationNode;

/// A spinney of the wood at a position of a topological order of the forests.
pub type ForestNode<'a> = UvForestNode<'a, LoopMomentumBasis>;
/// A forest node owning the spinney covered by an operation.
pub type OwnedForestNode = UvOwnedForestNode<LoopMomentumBasis>;

use crate::{
    debug_tags,
    graph::{Graph, LMBext, LoopMomentumBasis, cuts::CutSet, parse::string_utils::ToOrderedSimple},
    model::Model,
    settings::global::GenerationSettings,
    utils::{GS, W_},
    uv::{
        ApproximationType, RenormalizationPart, Spinney, UVgenerationSettings, UltravioletGraph,
        UvGraph,
        approx::{
            CutStructure, ForestNodeLike, OrientationProjection, Rooted, UVCtx,
            direct_3d::{Direct3dApproximation, Direct3dCts},
            final_integrand::{FinalIntegrandBuilder, FinalIntegrands},
            integrated::{Integrated, IntegratedCts},
            local_3d::{Local3DCts, Localizer},
            local_4d::{self, Full4dCts, Local4dCts},
            projected_4d::{Local4dProjectionContext, Projected4dApproximation},
        },
        export::{
            UVForestLocalProvenance, UVForestNodeExpression, UVForestNodeProvenance,
            UVForestProvenance,
        },
        forest::ParametricIntegrands,
        marker::UvMarker,
    },
};
use color_eyre::Result;
use spenso::shadowing::symbolica_utils::{LogPrint, SpensoPrintSettings};

/// The hedge-poset wood of a graph, with the cut structure and Vakint
/// settings that GammaLoop's counterterm computation needs.
pub struct Wood {
    wood: UvWood<LoopMomentumBasis>,
    vakint_settings: VakintSettings,
    cuts: CutStructure,
}

impl Wood {
    pub(crate) fn new(
        cuts: CutStructure,
        graph: &Graph,
        model: &Model,
        settings: &GenerationSettings,
    ) -> Result<Self> {
        let mut subgraph = graph.full_filter();
        subgraph.subtract_with(&graph.initial_state_cut.left);
        let mut spinneys = Vec::new();

        for (cut_index, cut) in cuts.cuts.iter().enumerate() {
            let cut_sub = subgraph.subtract(&cut.union);
            spinneys.extend(
                graph
                    .classified_spinneys(&cut_sub, model, &settings.uv, &graph.loop_momentum_basis)
                    .wrap_err_with(|| {
                        format!(
                            "graph '{}' cut {cut_index}: failed to classify UV counterterm components",
                            graph.name
                        )
                    })?,
            );
        }
        graph.remove_full_observable_spinney(&mut spinneys, &settings.medium);

        Self::from_spinneys(spinneys, graph, cuts, &settings.uv)
    }

    pub(crate) fn from_spinneys<I: IntoIterator<Item = Spinney>>(
        s: I,
        graph: &Graph,
        cuts: CutStructure,
        settings: &UVgenerationSettings,
    ) -> Result<Self> {
        let wood = UvWood::from_spinneys(s, graph, &graph.loop_momentum_basis)?;
        let mut vakint_settings = settings.vakint.true_settings();
        vakint_settings.project_onto_tensor_integrals =
            settings.project_integrated_uv_cts_onto_tensor_integrals;
        // Retain enough positive epsilon powers for finite terms formed when disconnected
        // integrated counterterms are multiplied.
        vakint_settings.number_of_terms_in_epsilon_expansion = wood.max_loop_count() as i64 + 1;

        Ok(Wood {
            wood,
            cuts,
            vakint_settings,
        })
    }

    /// The spinney poset.
    pub fn graph(&self) -> &HedgeGraph<SuBitGraph, Spinney> {
        self.wood.graph()
    }

    /// The combinatorial wood.
    pub fn wood(&self) -> &UvWood<LoopMomentumBasis> {
        &self.wood
    }

    fn unfold_with_cached_node_label_atoms(self, cache_node_label_atoms: bool) -> Forests {
        let forests = self.wood.unfold();

        let mut cuts: Vec<(SuBitGraph, CutSet)> = Vec::new();
        for c in &self.cuts.cuts {
            let mut compatible: SuBitGraph = forests.graph().empty_subgraph();
            for (_, crown, s) in forests.graph().iter_nodes() {
                if s.is_compatible_with(&c.union) {
                    for h in crown {
                        compatible.add(h);
                    }
                }
            }
            cuts.push((compatible, c.clone()));
        }

        let forests = Forests {
            forests,
            cuts,
            vakint_settings: self.vakint_settings,
            cached_node_label_atoms: cache_node_label_atoms.then(AHashMap::new),
            compute_store: ComputeStore::default(),
        };

        if cache_node_label_atoms {
            forests.with_cached_node_label_atoms()
        } else {
            forests
        }
    }

    pub fn unfold(self) -> Forests {
        self.unfold_with_cached_node_label_atoms(true)
    }

    pub fn unfold_uncached(self) -> Forests {
        self.unfold_with_cached_node_label_atoms(false)
    }
}

impl Display for Wood {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.wood.fmt(f)
    }
}

#[derive(Default)]
pub struct ComputeStore {
    entries: AHashMap<OperationNode, ComputeNode>,
}

impl ComputeStore {
    fn get(&self, key: &OperationNode) -> Option<&ComputeNode> {
        self.entries.get(key)
    }

    fn require(&self, key: &OperationNode) -> Result<&ComputeNode> {
        // Union operations, like individual Taylor steps, must be computed
        // before their typed values are consumed by another forest node.
        self.get(key)
            .ok_or_else(|| eyre!("{key} not yet added to compute store"))
    }

    fn entry(
        &mut self,
        key: OperationNode,
    ) -> std::collections::hash_map::Entry<'_, OperationNode, ComputeNode> {
        self.entries.entry(key)
    }
}

/// The unfolded forests of a [`Wood`] with GammaLoop's per-forest computed
/// counterterms.
pub struct Forests {
    forests: UvForests<LoopMomentumBasis>,
    /// Wood subgraph that has compatible
    cuts: Vec<(SuBitGraph, CutSet)>,
    vakint_settings: VakintSettings,
    cached_node_label_atoms: Option<AHashMap<NodeIndex, Atom>>,
    pub compute_store: ComputeStore,
}

#[cfg(test)]
struct LocalLeafOperation {
    op: Operation,
    frontier: NodeIndex,
}

#[cfg(test)]
impl LocalLeafOperation {
    fn new(op: &Operation, frontier: NodeIndex) -> Self {
        Self {
            op: op.clone(),
            frontier,
        }
    }
}

impl ForestNodeLike for ForestNode<'_> {
    fn dod(&self) -> i32 {
        self.spinney.dod()
    }
    fn renormalization_scheme(&self) -> ApproximationType {
        self.spinney.renormalization_scheme()
    }
    fn lmb(&self) -> &LoopMomentumBasis {
        self.spinney.lmb()
    }
    fn reduced_subgraph(&self, given: &Self) -> SuBitGraph {
        UvForestNode::reduced_subgraph(self, given)
    }
    fn subgraph(&self) -> &SuBitGraph {
        self.spinney.filter()
    }
    fn topo_order(&self) -> usize {
        self.topo_order
    }
}

impl ForestNodeLike for OwnedForestNode {
    fn dod(&self) -> i32 {
        self.spinney.dod()
    }

    fn renormalization_scheme(&self) -> ApproximationType {
        self.spinney.renormalization_scheme()
    }

    fn lmb(&self) -> &LoopMomentumBasis {
        self.spinney.lmb()
    }

    fn reduced_subgraph(&self, given: &Self) -> SuBitGraph {
        UvOwnedForestNode::reduced_subgraph(self, given)
    }

    fn subgraph(&self) -> &SuBitGraph {
        self.spinney.filter()
    }

    fn topo_order(&self) -> usize {
        self.topo_order
    }
}

#[derive(Default)]
pub struct ComputeNode {
    local_4d: Option<Local4dCts>,
    integrated: Option<IntegratedCts>,
    cuts: AHashMap<CutSet, CutComputation>,
}

pub struct CutComputation {
    local_3d: Local3DCts,
    final_integrands: FinalIntegrands,
}

impl ComputeNode {
    fn local_4d(&self, operation: &OperationNode) -> Result<&Local4dCts> {
        self.local_4d
            .as_ref()
            .ok_or_else(|| eyre!("{operation} has no computed local 4D counterterm"))
    }

    fn integrated(&self, operation: &OperationNode) -> Result<&IntegratedCts> {
        self.integrated
            .as_ref()
            .ok_or_else(|| eyre!("{operation} has no computed integrated counterterm"))
    }

    fn cut(&self, operation: &OperationNode, cutset: &CutSet) -> Result<&CutComputation> {
        self.cuts.get(cutset).ok_or_else(|| {
            eyre!("{operation} has no computed local counterterms for cut {cutset:?}")
        })
    }
}

impl Forests {
    /// The unfolded trace graph.
    pub fn graph(
        &self,
    ) -> &UnfoldedTraceGraph<UvWood<LoopMomentumBasis>, OperationNode, NodeStorageVec<OperationNode>>
    {
        self.forests.graph()
    }

    /// The combinatorial forests.
    pub fn forests(&self) -> &UvForests<LoopMomentumBasis> {
        &self.forests
    }

    fn source_spinney(&self, node: NodeIndex) -> &Spinney {
        self.forests.source_spinney(node)
    }

    #[cfg(test)]
    fn local_leaf_operations(&self, node: NodeIndex) -> Vec<LocalLeafOperation> {
        let mut leaves = self
            .forests
            .leaf_op_dependency_frontiers(node)
            .map(|(op, frontier)| LocalLeafOperation::new(op, frontier))
            .collect::<Vec<_>>();

        leaves.sort_by_key(|leaf| {
            (
                Reverse(self.forests.graph()[leaf.frontier].key.op_count()),
                leaf.op.order.clone(),
                usize::from(leaf.op.data),
            )
        });
        leaves
    }

    fn recursion_input_4d(&self, node: NodeIndex) -> Result<Full4dCts> {
        // The typed 4D value is the sequential accumulator. It already encodes that an empty
        // dependency frontier starts at the typed root, that other prefixes enter through their
        // reduced branch, and that both local subtraction signs are applied without a raw
        // Foata-level product. Disconnected values likewise retain their factorized local and
        // integrated-prefix cross terms, so no separate active/frozen root-path replay is needed.
        let operation = &self.forests.graph()[node];
        let computed = self.compute_store.require(operation)?;
        if self.forests.graph().is_disjoint_union(node) {
            return Ok(Full4dCts::from_factorized_local(
                computed.local_4d(operation)?,
            ));
        }
        Full4dCts::recursion_input(
            computed.local_4d(operation)?,
            computed.integrated(operation)?,
            self.source_spinney(node).renormalization_scheme(),
            node == self.forests.root(),
            self.source_spinney(node).lmb(),
        )
    }

    fn cached_node_label_atom(&self, node: NodeIndex) -> Option<Atom> {
        self.cached_node_label_atoms
            .as_ref()
            .and_then(|labels| labels.get(&node))
            .cloned()
    }

    fn node_label_atom_factor(&self, frontier: NodeIndex, op: &Operation) -> Atom {
        let approx = FunctionBuilder::new(GS.uv_approx);
        let frontier_atom = self.node_label_atom(frontier);

        let current = function!(op.order.symbol(), usize::from(op.data));
        let argument = if self.forests.graph()[frontier].covers().is_none() {
            current
        } else {
            let previous = self.forests.graph()[frontier]
                .covers()
                .expect("non-empty frontier cover must exist")
                .symbol();
            (current - previous) * frontier_atom
        };
        approx.add_arg(argument).finish()
    }

    fn node_label_atom(&self, node: NodeIndex) -> Atom {
        if self.forests.graph()[node].key.is_empty() {
            return Atom::one();
        }

        if let Some(cached) = self.cached_node_label_atom(node) {
            return cached;
        }

        self.forests
            .leaf_op_dependency_frontiers(node)
            .fold(Atom::one(), |acc, (op, frontier)| {
                acc * self.node_label_atom_factor(frontier, op)
            })
    }

    fn with_cached_node_label_atoms(mut self) -> Self {
        self.cache_node_label_atoms();
        self
    }

    fn cache_node_label_atoms(&mut self) {
        self.cached_node_label_atoms = Some(AHashMap::new());
        for nidx in self
            .forests
            .graph()
            .topo_sort_kahn(DirectionBasis::Underlying)
            .unwrap()
        {
            let atom = self.node_label_atom(nidx);
            self.cached_node_label_atoms
                .as_mut()
                .expect("node-label cache was initialized")
                .insert(nidx, atom);
        }
    }

    fn node_label(&self, node: NodeIndex) -> String {
        let key = &self.forests.graph()[node];
        if key.key.is_empty() {
            return "∅".to_string();
        }

        let mut foata = String::new();
        key.key
            .write_foata_like(&mut foata, |op| op.string_label())
            .expect("writing a trace key into a string must succeed");
        let atom = self.node_label_atom(node);
        format!("{foata}: {}", atom.to_ordered_simple())
    }

    fn dot_serialize_expr_atom() -> Atom {
        Atom::var(GS.expr)
    }

    fn dot_serialize_node_atom_factor(&self, frontier: NodeIndex, op: &Operation) -> Atom {
        let frontier_depth = self.forests.graph()[frontier].key.op_count();
        let (current, given) = self.forests.current_given_pair(op.data, frontier_depth);
        // Structured forest-DOT approximation records use the same operation and subgraph
        // markers as the computed UV expressions.
        function!(
            GS.uv_approx,
            UvMarker::subgraph(current.subgraph(), given.subgraph())
                * self.dot_serialize_node_atom(frontier)
        )
    }

    fn dot_serialize_node_atom(&self, node: NodeIndex) -> Atom {
        if self.forests.graph()[node].key.is_empty() {
            return Self::dot_serialize_expr_atom();
        }

        self.forests
            .leaf_op_dependency_frontiers(node)
            .fold(Atom::one(), |acc, (op, frontier)| {
                acc * self.dot_serialize_node_atom_factor(frontier, op)
            })
    }

    /// The Typst label of a forest in the structured DOT export.
    fn dot_serialize_node_label(&self, node: NodeIndex) -> String {
        self.dot_serialize_node_atom(node)
            .printer(SpensoPrintSettings::typst_options())
            .to_string()
    }

    pub fn dot_serialize(&self) -> String {
        let mut output = String::new();
        self.dot_serialize_fmt(&mut output)
            .expect("writing hedge-poset forest DOT into a string must succeed");
        output
    }

    pub fn dot_serialize_fmt(&self, writer: &mut impl fmt::Write) -> fmt::Result {
        self.forests
            .dot_serialize_with(writer, |node| self.dot_serialize_node_label(node))
    }

    fn compute_4d_for_node(
        &self,
        node: NodeIndex,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
    ) -> Result<(Local4dCts, IntegratedCts)> {
        let operation = &self.forests.graph()[node];
        if operation.key.is_empty() {
            return Ok((Local4dCts::root(), IntegratedCts::root()));
        }

        if self.forests.graph().is_disjoint_union(node) {
            let components = self.forests.disconnected_component_nodes(node)?;
            let mut full_components = Vec::with_capacity(components.len());
            let mut integrated_components = Vec::with_capacity(components.len());
            for component in components {
                let component_operation = &self.forests.graph()[component];
                full_components.push(self.recursion_input_4d(component).wrap_err_with(|| {
                    format!("while loading 4D component {component_operation} for {operation}")
                })?);
                integrated_components.push(
                    self.compute_store
                        .require(component_operation)?
                        .integrated(component_operation)?,
                );
            }

            let depth = graph.n_loops(
                &operation
                    .covers()
                    .expect("a non-root operation has a cover"),
            ) + 1;
            return Ok((
                Local4dCts::from_full_product(full_components)?,
                IntegratedCts::factorized_product(integrated_components, depth)?,
            ));
        }

        let (parent, edge) = self
            .forests
            .graph()
            .unique_parent(node)
            .ok_or_else(|| eyre!("{operation} has no unique 4D parent"))?;
        let full = self
            .recursion_input_4d(parent)
            .wrap_err_with(|| format!("while loading the 4D parent for {operation}"))?;
        let step_order = self.forests.graph()[parent].key.op_count();
        let (current, given) = self.forests.current_given_pair(edge, step_order);
        let ctx = UVCtx::new(graph, model, settings);
        let integrated_approximation = Integrated::new(vakint, &self.vakint_settings);
        let local = local_4d::uv_limit(&full, &ctx, &current, &given, &current, &given)?;
        let integrated = if settings.generate_integrated {
            integrated_approximation.run(&local, &ctx, &current, &given, &current, &given)?
        } else {
            IntegratedCts::root()
        };

        Ok((local, integrated))
    }

    fn local_3d_for_node(
        &self,
        node: NodeIndex,
        graph: &mut Graph,
        cutset: &CutSet,
        localizer: Localizer<'_>,
        settings: &UVgenerationSettings,
        projection_context: &mut Local4dProjectionContext,
    ) -> Result<CutComputation> {
        let operation = &self.forests.graph()[node];
        let forest_node = ForestNode {
            spinney: self.source_spinney(node),
            topo_order: operation.key.op_count(),
        };
        let local_3d = if operation.key.is_empty() {
            Local3DCts::Direct(Direct3dCts::root(graph, localizer)?)
        } else if settings.local_uv_cts_from_expanded_4d_integrands {
            // This is the only route that may produce `Projected4d`: the 4D
            // Taylor coefficient is kept factorized and final assembly later
            // attaches its untouched outer CFF. It never enters any of the
            // `Direct3dApproximation` replay operations below.
            let local_4d = self.compute_store.require(operation)?.local_4d(operation)?;
            Local3DCts::Projected4d(
                Projected4dApproximation::new(localizer, graph, settings)
                    .project_local_4d(local_4d, projection_context)?,
            )
        } else if self.forests.graph().is_disjoint_union(node) {
            // Both direct variants replay the Taylor operators on the complete
            // post-energy-integration CFF. `explicit_orientation_sum_only`
            // changes selector materialization only, not this construction.
            let mut active_sectors = Vec::new();
            for state in self
                .forests
                .union_replay_states(node)?
                .into_iter()
                .filter(|state| !state.local_edges.is_empty())
            {
                let integrated_operation = &self.forests.graph()[state.integrated];
                let mut edges = state.local_edges.into_iter().enumerate();
                let (offset, first_edge) = edges
                    .next()
                    .expect("a proper integrated prefix has a local suffix");
                let step_order = integrated_operation.key.op_count() + offset;
                let (current, given) = self.forests.current_given_pair(first_edge, step_order);

                // An empty integrated prefix starts from the per-cut root integrand. Every
                // other prefix enters through the reduced branch of its first local operation.
                let mut sector = if integrated_operation.key.is_empty() {
                    let root = Direct3dCts::root(graph, localizer)?;
                    Direct3dApproximation::new(localizer, graph, settings)
                        .run_local(&root, &current, &given, &current, &given)?
                } else {
                    // Retain each connected coefficient's owner before the
                    // product enters a local suffix. Equal physical masses in
                    // a disjoint sibling must not acquire the parent's grading.
                    // Keep disconnected vacuum masses owned as well; the
                    // complete union still determines the localization kernel.
                    let components = if self.forests.graph().is_disjoint_union(state.integrated) {
                        self.forests
                            .disconnected_component_nodes(state.integrated)?
                    } else {
                        vec![state.integrated]
                    };
                    let finite_components = components
                        .into_iter()
                        .map(|component| {
                            let operation = &self.forests.graph()[component];
                            let integrated = self
                                .compute_store
                                .require(operation)?
                                .integrated(operation)?;
                            Ok((
                                integrated.physical_finite_counterterm()?,
                                ForestNode {
                                    spinney: self.source_spinney(component),
                                    topo_order: operation.key.op_count(),
                                },
                            ))
                        })
                        .collect::<Result<Vec<_>>>()?;
                    let finite_counterterms = finite_components
                        .iter()
                        .map(|(coefficient, owner)| (coefficient, owner))
                        .collect::<Vec<_>>();
                    // A disabled (or algebraically vanishing) integrated prefix
                    // still owns its typed component-path frame. Replay that
                    // zero coefficient through the local suffix; it remains
                    // algebraically zero without deleting the union sector.
                    let prefix_node = ForestNode {
                        spinney: self.source_spinney(state.integrated),
                        topo_order: integrated_operation.key.op_count(),
                    };
                    Direct3dApproximation::new(localizer, graph, settings).run_integrated(
                        &finite_counterterms,
                        &prefix_node,
                        &current,
                        &given,
                        &current,
                        &given,
                    )?
                };

                for (offset, edge) in edges {
                    let step_order = integrated_operation.key.op_count() + offset;
                    let (current, given) = self.forests.current_given_pair(edge, step_order);
                    sector = Direct3dApproximation::new(localizer, graph, settings)
                        .run_local(&sector, &current, &given, &current, &given)?;
                }

                // Keep each active/frozen split after its root-path replay so any connected
                // descendants rescale only the loop variables still active in that sector.
                active_sectors.extend(sector.sectors()?.iter().cloned());
            }

            Local3DCts::Direct(
                Direct3dCts::from_sectors(active_sectors)
                    .wrap_err_with(|| format!("{operation} has no proper integrated prefixes"))?,
            )
        } else {
            let (parent, edge) = self
                .forests
                .graph()
                .unique_parent(node)
                .ok_or_else(|| eyre!("{operation} has no unique local parent"))?;
            let parent_operation = &self.forests.graph()[parent];
            // An empty dependency frontier starts from the per-cut root integrand;
            // otherwise its typed local result remains the sequential accumulator.
            let root;
            let parent_local = if parent_operation.key.is_empty() {
                root = Direct3dCts::root(graph, localizer)?;
                &root
            } else {
                self.compute_store
                    .require(parent_operation)?
                    .cut(parent_operation, cutset)?
                    .local_3d
                    .direct()?
            };
            let parent_integrated = self
                .compute_store
                .require(parent_operation)?
                .integrated(parent_operation)
                .wrap_err_with(|| {
                    format!("while loading integrated parent {parent_operation} for {operation}")
                })?;
            let step_order = parent_operation.key.op_count();
            let (current, given) = self.forests.current_given_pair(edge, step_order);
            // `run` applies the next Taylor operator to the complete CFF and
            // applies both subtraction signs; no projected 4D coefficient,
            // external sign, or raw Foata-level product enters this route.
            Local3DCts::Direct(Direct3dApproximation::new(localizer, graph, settings).run(
                parent_local,
                parent_integrated,
                &current,
                &given,
                &current,
                &given,
            )?)
        };

        let integrated = self
            .compute_store
            .require(operation)?
            .integrated(operation)?;
        let final_builder = FinalIntegrandBuilder::new(localizer, settings);
        let final_integrands = match &local_3d {
            Local3DCts::Direct(direct) => {
                final_builder.build_direct(graph, &forest_node, direct, integrated)?
            }
            Local3DCts::Projected4d(projected) => {
                final_builder.build_projected(graph, &forest_node, projected, integrated)?
            }
        };

        Ok(CutComputation {
            local_3d,
            final_integrands,
        })
    }

    pub fn integrate(
        &mut self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        for (order, nidx) in self
            .forests
            .graph()
            .topo_sort_kahn(DirectionBasis::Underlying)?
            .into_iter()
            .enumerate()
        {
            debug!(order, nidx=%nidx, key=%self.forests.graph()[nidx], "Computing hedge-poset 4D term");
            let operation = self.forests.graph()[nidx].clone();
            let (local_4d, integrated) =
                self.compute_4d_for_node(nidx, graph, model, vakint, settings)?;
            let cover = operation
                .covers()
                .unwrap_or_else(|| self.forests.graph().empty_subgraph())
                .string_label();
            let source = self.source_spinney(nidx).filter().string_label();
            debug_tags!(#generation, #uv, #integrated, #graph, #term, #inspect;
                stage = "hedge_poset_4d_node_done",
                order,
                node_index = %nidx,
                forest_term = %operation,
                cover = %cover,
                source = %source,
                is_union = self.forests.graph().is_disjoint_union(nidx),
                log.local_4d = local_4d.atom().expression(),
                log.integrated_pole = integrated.physical_pole()?.expression(),
                log.integrated_finite = integrated.physical_finite_counterterm()?.expression(),
                "Computed hedge-poset 4D node"
            );
            let computed = self.compute_store.entry(operation).or_default();
            computed.local_4d = Some(local_4d);
            computed.integrated = Some(integrated);
        }

        Ok(())
    }

    /// Prepare the representation-independent 4D lane once, at the first
    /// selected orientation. The requested model remains the source of rules.
    pub(crate) fn prepare_four_d(
        &mut self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
        projection_context: &mut Local4dProjectionContext,
    ) -> Result<()> {
        let started = Instant::now();
        if settings.generate_integrated || settings.local_uv_cts_from_expanded_4d_integrands {
            // Local-4D projection also needs these coefficients when analytic addbacks are off.
            self.integrate(graph, model, vakint, settings)?;
        } else {
            // Direct local-3D consumes only the zero integrated lane. Seed every
            // operation because a disconnected replay can visit a prefix first.
            let zero = IntegratedCts::root();
            for (_, _, operation) in self.forests.graph().iter_nodes() {
                self.compute_store
                    .entry(operation.clone())
                    .or_default()
                    .integrated = Some(zero.clone());
            }
        }
        projection_context.shared_preparation_time += started.elapsed();
        Ok(())
    }

    /// Return the measured representation-independent part of this computation.
    pub(crate) fn compute(
        &mut self,
        graph: &mut Graph,
        orientation: OrientationProjection<'_>,
        settings: &UVgenerationSettings,
        projection_context: &mut Local4dProjectionContext,
    ) -> Result<Duration> {
        let shared_before = projection_context.shared_preparation_time;
        for (compatible_subset, cutset) in self.cuts.clone() {
            let localizer = Localizer::new(&cutset, orientation);
            for (order, nidx) in self
                .forests
                .compatible_topological_order(&compatible_subset)?
                .into_iter()
                .enumerate()
            {
                debug!(order, nidx=%nidx, key=%self.forests.graph()[nidx], "Computing hedge-poset per-cut term");
                let operation = self.forests.graph()[nidx].clone();
                // Direct local-3D nodes Taylor-expand the complete post-energy-integration CFF.
                // Expanded-4D nodes instead project their typed local coefficients and attach
                // the outer CFF only during final assembly.
                let cut_computation = self.local_3d_for_node(
                    nidx,
                    graph,
                    &cutset,
                    localizer,
                    settings,
                    projection_context,
                )?;
                self.compute_store
                    .entry(operation)
                    .or_default()
                    .cuts
                    .insert(cutset.clone(), cut_computation);
            }
        }

        Ok(projection_context.shared_preparation_time - shared_before)
    }

    pub(crate) fn orientation_parametric_exprs(
        &self,
        graph: &Graph,
        _settings: &UVgenerationSettings,
    ) -> Result<Vec<ParametricIntegrands>> {
        let split_momentum_replacements = graph
            .iter_edges_of(
                &graph
                    .full_filter()
                    .subtract(&graph.initial_state_cut)
                    .subtract(&graph.tree_edges),
            )
            .filter_map(|(pair, edge_index, _)| {
                (!matches!(pair, HedgePair::Unpaired { .. }))
                    .then(|| GS.split_mom_pattern_simple(edge_index))
            })
            .collect::<Vec<_>>();
        let mut expressions = Vec::with_capacity(self.cuts.len());

        for (compatible_subset, cutset) in &self.cuts {
            let mut sum: Option<FinalIntegrands> = None;
            for nidx in self
                .forests
                .compatible_topological_order(compatible_subset)?
            {
                let operation = &self.forests.graph()[nidx];
                let terms = self
                    .compute_store
                    .require(operation)?
                    .cut(operation, cutset)?
                    .final_integrands
                    .collect(TensorCollectFilter::Reps([
                        ColorAdjoint {}.into(),
                        ColorFundamental {}.into(),
                        ColorSextet {}.into(),
                    ]))?;
                sum = Some(match sum {
                    Some(sum) => sum.zip_add([terms]).wrap_err_with(|| {
                        format!("while aggregating hedge-poset term {operation} for cut {cutset:?}")
                    })?,
                    None => terms,
                });
            }
            // A per-node final integrand retains inspect/export provenance but
            // is not independently evaluator-ready.  In particular, temporal
            // components in localized integrated counterterms can cancel only
            // in the completed signed forest.  Preserve that exact cancellation
            // before the final orientation-local momentum split.
            let integrands = sum
                .ok_or_else(|| eyre!("No terms in hedge-poset forest for cut {cutset:?}"))?
                .map_expressions(|integrand| {
                    Ok(integrand
                        .replace_multiple(&split_momentum_replacements)
                        .replace(function!(GS.den, W_.a_, W_.b_, W_.c_, W_.d_))
                        .with(W_.d_))
                })?;
            expressions.push(ParametricIntegrands::from_final(integrands, cutset.clone()));
        }

        Ok(expressions)
    }

    pub(crate) fn export_node_expressions(
        &self,
        forest_index: usize,
        post_process: &mut impl FnMut(Atom) -> Atom,
    ) -> Result<Vec<UVForestNodeExpression>> {
        let (cut_compatible_forest_subset, cutset) = self
            .cuts
            .first()
            .ok_or_else(|| eyre!("No cuts in hedge-poset forest export"))?;
        let mut terms = Vec::new();
        for (node_index, nidx) in self
            .forests
            .compatible_topological_order(cut_compatible_forest_subset)?
            .into_iter()
            .enumerate()
        {
            let operation = &self.forests.graph()[nidx];
            let computed = self.compute_store.require(operation)?;
            let cut_computation = computed.cut(operation, cutset)?;
            let final_integrands = cut_computation
                .final_integrands
                .clone()
                .into_integrands()
                .resolved()?
                .map(|numerator| post_process(numerator.clone()));
            let node_key = operation.to_string();
            let provenance = self.node_export_3d_provenance(nidx, &cut_computation.local_3d);
            for (term_index, (&residue_index, numerator)) in final_integrands.iter().enumerate() {
                terms.push(UVForestNodeExpression {
                    forest_index,
                    node_index,
                    node_key: node_key.clone(),
                    term_index,
                    residue_index,
                    numerator: numerator.clone(),
                    forest_provenance: Some(provenance.clone()),
                });
            }
        }

        Ok(terms)
    }

    pub(crate) fn export_4d_node_expressions(
        &self,
        forest_index: usize,
    ) -> Result<Vec<UVForestNodeExpression>> {
        let (cut_compatible_forest_subset, _) = self
            .cuts
            .first()
            .ok_or_else(|| eyre!("No cuts in hedge-poset forest export"))?;
        self.forests
            .compatible_topological_order(cut_compatible_forest_subset)?
            .into_iter()
            .enumerate()
            .map(|(node_index, nidx)| {
                let operation = &self.forests.graph()[nidx];
                let local_4d = self.compute_store.require(operation)?.local_4d(operation)?;
                Ok(UVForestNodeExpression {
                    forest_index,
                    node_index,
                    node_key: operation.to_string(),
                    term_index: 0,
                    residue_index: crate::cff::CutCFFIndex::new_all_none(),
                    numerator: local_4d.atom().expression().clone(),
                    forest_provenance: Some(self.node_export_4d_provenance(nidx, local_4d)),
                })
            })
            .collect()
    }

    fn node_parent_keys(&self, node: NodeIndex) -> Vec<String> {
        let mut parent_keys = self
            .forests
            .graph()
            .parents(node)
            .map(|(parent, _)| self.forests.graph()[parent].to_string())
            .collect::<Vec<_>>();
        parent_keys.sort();
        parent_keys.dedup();
        parent_keys
    }

    fn node_export_4d_provenance(
        &self,
        node: NodeIndex,
        local_4d: &Local4dCts,
    ) -> Arc<UVForestProvenance> {
        Arc::new(UVForestProvenance {
            node: UVForestNodeProvenance {
                parent_keys: self.node_parent_keys(node),
                local: UVForestLocalProvenance::FourD {
                    components: local_4d.provenance_components().to_vec(),
                    branches: local_4d.provenance_branches().to_vec(),
                },
            },
        })
    }

    fn node_export_3d_provenance(
        &self,
        node: NodeIndex,
        local_3d: &Local3DCts,
    ) -> Arc<UVForestProvenance> {
        Arc::new(UVForestProvenance {
            node: UVForestNodeProvenance {
                parent_keys: self.node_parent_keys(node),
                local: UVForestLocalProvenance::ThreeD {
                    projection_paths: local_3d.projection_paths(),
                },
            },
        })
    }

    pub(crate) fn renormalization_part_of_ends(
        &self,
        graph: &Graph,
        settings: &UVgenerationSettings,
    ) -> Result<RenormalizationPart> {
        let mut sum = Atom::Zero;
        let marker = UvMarker::new(settings);

        let wild = Atom::var(W_.x___);

        let replacements =
            graph.integrand_replacement(&graph.full_filter(), &graph.loop_momentum_basis, &[wild]);
        for (node, mut crown, key) in self.forests.graph().iter_nodes() {
            if crown.any(|r| self.forests.graph().flow(r).is_source()) {
                continue;
            }

            let forest_node = key
                .forest_node(graph, &graph.loop_momentum_basis, key.key.op_count())
                .wrap_err_with(|| format!("graph '{}'", graph.name))?;
            // A disconnected terminal can combine components with different schemes. Select
            // each component's projection before multiplying; the aggregate integrated value
            // only represents the homogeneous all-finite and all-pole projections.
            let components = if self.forests.graph().is_disjoint_union(node) {
                self.forests.disconnected_component_nodes(node)?
            } else {
                vec![node]
            };
            let physical = components.into_iter().try_fold(
                SymbolicTensor::infer_indexed(Atom::one())?,
                |product, component| -> Result<_> {
                    let component_key = &self.forests.graph()[component];
                    let integrated = self
                        .compute_store
                        .require(component_key)?
                        .integrated(component_key)?;
                    let projection = match self.source_spinney(component).renormalization_scheme() {
                        ApproximationType::MUV | ApproximationType::IR => {
                            integrated.physical_finite_counterterm()?
                        }
                        ApproximationType::PolePart => integrated.physical_pole()?,
                        scheme => {
                            return Err(eyre!("No terminal counterterm projection for {scheme}"));
                        }
                    };
                    Ok(product.multiply(projection)?)
                },
            )?;
            let marked = marker.prefix(
                &graph.full_filter(),
                forest_node.subgraph(),
                physical.expression(),
            );
            let physical = physical.with_checked_expression(marked)?;
            debug!(
                key=%key,
               expr = % physical.expression().log_print(None),"Term before simplification"
            );
            let overall = SymbolicTensor::infer_indexed(graph.overall_factor.clone())?;
            let prefactor = graph
                .global_prefactor
                .projector
                .clone()
                .multiply(graph.global_prefactor.num.clone())?
                .multiply(overall)?;
            let cooking = CookSettings::indices()
                .with_mode(CookMode::ReversibleEncoding)
                .with_representation_payloads(true, true);
            let reduced = physical
                .multiply(prefactor)?
                .simplify_algebra(&idenso::tensor::AlgebraSettings {
                    color: Some(ColorSimplifySettings::default()),
                    ..Default::default()
                })?
                .contract(idenso::tensor::ContractSettings {
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                })?
                .to_dots()?;
            let atom = cooking.uncook(reduced.expression().as_view());

            debug!(
                key=%key,
               expr = % atom.log_print(None),"Term"
            );
            sum += atom;
        }

        Ok(RenormalizationPart::new(
            sum.replace_multiple(&replacements)
                .replace(GS.m_uv_expansion)
                .with(GS.m_uv_vacuum),
            self.forests.graph().n_nodes(),
        ))
    }

    pub fn debug_walk(&self) {
        let mut cover_groups: BTreeMap<SuBitGraph, Vec<NodeIndex>> = BTreeMap::new();

        self.forests
            .graph()
            .topo_sort_kahn(DirectionBasis::Underlying)
            .unwrap()
            .iter()
            .for_each(|nidx| {
                let trace_key = &self.forests.graph()[*nidx];
                cover_groups
                    .entry(
                        trace_key
                            .covers()
                            .unwrap_or(self.forests.graph().empty_subgraph()),
                    )
                    .and_modify(|e| e.push(*nidx))
                    .or_insert_with(|| vec![*nidx]);

                println!("Node {}:{}", nidx, self.node_label(*nidx));
            });

        println!("edge [constraint=true style=invis];");
        for (a, b) in cover_groups.values().tuple_windows() {
            println!("{}->{}", a.first().unwrap(), b.first().unwrap())
        }
        println!("edge [style=solid];");

        for (s, g) in cover_groups.iter() {
            println!("subgraph group_{} {{rank=same; ", s.string_label());
            for n in g {
                println!("{};", n.0);
            }
            println!("}}");
        }
    }
}

impl Display for Forests {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let labels: AHashMap<_, _> = self
            .forests
            .graph()
            .iter_nodes()
            .map(|(node, _, key)| (key.clone(), self.node_label(node)))
            .collect();
        self.forests.graph().dot_impl_fmt(
            f,
            &self.forests.graph().full_filter(),
            "start=2;\n",
            &|_| None,
            &|_| None,
            &|v| Some(format!("label=\"{}\"", labels[v])),
        )
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use feynkit_generator::FeynmanDiagramRulesExt;
    use feynkit_graph::FeynmanDiagram;
    use feynkit_model::Model as FeynkitModel;
    use linnet::half_edge::algorithms::trace_unfold::TraceUnfold;
    use linnet::half_edge::involution::{EdgeIndex, Flow};
    use linnet::half_edge::subgraph::{Inclusion, InternalSubGraph};

    use crate::feyngen::feynkit::FeynmanDiagramGammaLoopExt;

    use crate::initialisation;
    use crate::utils;
    use crate::uv::uv_graph::UVE;
    use crate::{
        finalized_runtime_dot,
        graph::{FeynmanGraph, Graph},
        initialisation::test_initialise,
        processes::DotExportSettings,
        uv::{
            UltravioletGraph, Wood as OldWood, approx::local_4d::Local4dBranch,
            settings::RenormalizationPrescriptionSettings,
        },
    };
    use idenso::color::ColorSimplifySettings;
    use idenso::tensor::{AlgebraContraction, AlgebraSettings, ContractSettings};

    use super::*;
    use color_eyre::Result;

    fn scalar_model() -> Model {
        crate::utils::load_generic_model("scalars")
    }

    fn standard_model() -> Model {
        crate::utils::load_generic_model("sm")
    }

    impl Forests {
        fn normalized_node_label(&self, node: NodeIndex) -> String {
            let label = self.node_label(node);
            let mut normalized = String::with_capacity(label.len());
            let mut chars = label.chars().peekable();

            while let Some(ch) = chars.next() {
                normalized.push(ch);
                if ch != '(' {
                    continue;
                }

                let mut digits = String::new();
                while chars.peek().is_some_and(|next| next.is_ascii_digit()) {
                    digits.push(chars.next().expect("peeked digit must exist"));
                }

                if digits.is_empty() || !matches!(chars.peek(), Some(')')) {
                    normalized.push_str(&digits);
                } else {
                    normalized.push('_');
                }
            }

            normalized
        }

        fn normalized_node_labels_with_cover(&self, cover_label: &str) -> Vec<String> {
            let mut labels = self
                .graph()
                .iter_nodes()
                .filter(|(_, _, operation)| {
                    operation
                        .covers()
                        .is_some_and(|cover| cover.string_label() == cover_label)
                })
                .map(|(node, _, _)| self.normalized_node_label(node))
                .collect::<Vec<_>>();
            labels.sort();
            labels
        }
    }

    #[test]
    fn three_d_compute_skips_unused_four_d_atoms_when_integration_is_disabled() -> Result<()> {
        test_initialise().unwrap();
        for project_local_4d in [false, true] {
            let mut graph: Graph = Graph::import_dot(
                stringify!(digraph scalar_self_energy {
                    edge [particle="scalar_1"];
                    node [num=1];
                    ext [style=invis];
                    ext -> a;
                    b -> ext;
                    a -> b;
                    a -> b;
                }),
                &crate::utils::load_generic_model("scalars"),
            )
            .map(|mut graphs| {
                assert_eq!(graphs.len(), 1, "one fixture graph expected");
                graphs.remove(0)
            })?;
            let settings = GenerationSettings {
                uv: UVgenerationSettings {
                    generate_integrated: false,
                    local_uv_cts_from_expanded_4d_integrands: project_local_4d,
                    ..Default::default()
                },
                ..Default::default()
            };
            let orientation_pattern = crate::settings::global::OrientationPattern::default();
            let options = graph.denominator_only_cff_3d_expression_options();
            let production = graph.generate_3d_expression_for_integrand(
                &[],
                &graph.get_esurface_canonization(&graph.loop_momentum_basis),
                &options,
                Some(&Atom::one()),
            )?;
            let mut forests =
                Wood::new(CutStructure::empty(&graph), &graph, &graph.model, &settings)?.unfold();

            let model = graph.model.clone();
            let mut context = Local4dProjectionContext::default();
            forests.prepare_four_d(
                &graph,
                &model,
                crate::utils::vakint()?,
                &settings.uv,
                &mut context,
            )?;
            forests.compute(
                &mut graph,
                OrientationProjection::exact_expression(
                    &production,
                    &options,
                    &orientation_pattern,
                    true,
                ),
                &settings.uv,
                &mut context,
            )?;

            if project_local_4d {
                let mut legacy = crate::uv::wood::CutWoods::new(
                    CutStructure::empty(&graph),
                    &graph,
                    &model,
                    &settings,
                )?
                .unfold(&graph);
                let mut legacy_context = Local4dProjectionContext::default();
                legacy.compute(
                    &mut graph,
                    crate::utils::vakint()?,
                    OrientationProjection::exact_expression(
                        &production,
                        &options,
                        &orientation_pattern,
                        true,
                    ),
                    &settings.uv,
                    true,
                    &mut legacy_context,
                )?;
                let canonical_builds = [
                    context.canonical_preparation_builds,
                    legacy_context.canonical_preparation_builds,
                ];
                assert!(canonical_builds.iter().all(|count| *count > 0));
                for representation in [
                    three_dimensional_reps::RepresentationMode::Ltd,
                    three_dimensional_reps::RepresentationMode::Cff,
                ] {
                    let options = three_dimensional_reps::Generate3DExpressionOptions {
                        representation,
                        ..options.clone()
                    };
                    let production = graph.generate_3d_expression_for_integrand(
                        &[],
                        &graph.get_esurface_canonization(&graph.loop_momentum_basis),
                        &options,
                        Some(&Atom::one()),
                    )?;
                    assert_eq!(production.representation, representation);
                    let orientation = OrientationProjection::exact_expression(
                        &production,
                        &options,
                        &orientation_pattern,
                        true,
                    );
                    let before = [&context, &legacy_context].map(|context| {
                        (
                            context.shared_preparation_time,
                            context.preparations.hits,
                            context.source_preparation_builds,
                            context.generation_cache.native_generations,
                        )
                    });
                    let shared = [
                        forests.compute(&mut graph, orientation, &settings.uv, &mut context)?,
                        legacy.compute(
                            &mut graph,
                            crate::utils::vakint()?,
                            orientation,
                            &settings.uv,
                            false,
                            &mut legacy_context,
                        )?,
                    ];
                    for (index, context) in [&context, &legacy_context].into_iter().enumerate() {
                        assert_eq!(
                            context.canonical_preparation_builds,
                            canonical_builds[index]
                        );
                        assert!(context.preparations.hits > before[index].1);
                        assert_eq!(
                            shared[index],
                            context.shared_preparation_time - before[index].0
                        );
                        if representation == three_dimensional_reps::RepresentationMode::Ltd {
                            assert!(context.source_preparation_builds > before[index].2);
                            assert!(context.generation_cache.native_generations > before[index].3);
                        } else {
                            assert_eq!(context.source_preparation_builds, before[index].2);
                            assert_eq!(
                                context.generation_cache.native_generations,
                                before[index].3
                            );
                        }
                    }
                }
            }

            assert_eq!(
                forests.compute_store.entries.len(),
                forests.graph().n_nodes()
            );
            for (_, _, operation) in forests.graph().iter_nodes() {
                let computed = forests.compute_store.require(operation)?;
                assert_eq!(
                    computed.local_4d.is_some(),
                    project_local_4d,
                    "4D coefficients must be constructed exactly when the projected route needs them for {operation}"
                );
                let integrated = computed.integrated(operation)?;
                assert!(integrated.physical_pole()?.is_zero());
                assert!(integrated.physical_finite_counterterm()?.is_zero());
                assert!(
                    !computed.cuts.is_empty(),
                    "both routes must compute a cut for {operation}"
                );
                assert!(
                    computed
                        .cuts
                        .values()
                        .all(|cut| cut.final_integrands.iter().next().is_some()),
                    "every computed cut must retain a 3D final integrand for {operation}"
                );
            }
        }

        Ok(())
    }

    #[test]
    fn vacuum_subtraction_prunes_only_full_observable_spinneys() -> Result<()> {
        use crate::settings::global::MediumMode;

        test_initialise()?;
        let graphs: Vec<Graph> = Graph::unit_numerator_fixtures(
            stringify!(
            digraph sunrise {
                edge [particle="scalar_1"];
                v1 -> v2; v1 -> v2; v1 -> v2;
            }
            digraph dumbbell {
                edge [particle="scalar_1"];
                v1 -> v2;
                v2 -> v2;
                v1 -> v1; v1 -> v1;
            }
            digraph sunrise_with_bridge {
                edge [particle="scalar_1"];
                ext [style=invis];
                v1 -> v2; v1 -> v2; v1 -> v2;
                v2 -> v3;
                ext -> v3; v3 -> ext;
            }),
            &crate::utils::load_generic_model("scalars"),
        )?;
        for graph in graphs {
            // The counterterm covers the loop part, leaving any bridge untouched.
            let full = InternalSubGraph::cleaned_filter_pessimist(
                graph.full_filter().subtract(&graph.tree_edges),
                graph.as_ref(),
            );
            for subtract_uv in [false, true] {
                let mut settings = GenerationSettings::default();
                settings.uv.subtract_uv = subtract_uv;
                let ordinary =
                    Wood::new(CutStructure::empty(&graph), &graph, &graph.model, &settings)?;
                let all = ordinary
                    .graph()
                    .iter_nodes()
                    .map(|(_, _, spinney)| spinney.filter().clone())
                    .collect::<BTreeSet<_>>();
                let proper = all
                    .iter()
                    .filter(|filter| **filter != full.filter)
                    .cloned()
                    .collect::<BTreeSet<_>>();
                assert!(proper.iter().any(|filter| filter.is_empty()));
                if subtract_uv {
                    assert!(all.contains(&full.filter), "{}", graph.name);
                    assert!(proper.len() > 1, "{}", graph.name);
                }
                for mode in [
                    MediumMode::Vacuum,
                    MediumMode::ThermodynamicEquilibrium,
                    MediumMode::ZeroTemperatureEquilibrium,
                ] {
                    settings.medium.mode = mode;
                    for vacuum_subtraction in [false, true] {
                        if mode == MediumMode::Vacuum && vacuum_subtraction {
                            continue; // Generation validation rejects this combination.
                        }
                        settings.medium.vacuum_subtraction = vacuum_subtraction;
                        let expected = if vacuum_subtraction { &proper } else { &all };
                        let wood = Wood::new(
                            CutStructure::empty(&graph),
                            &graph,
                            &graph.model,
                            &settings,
                        )?;
                        assert_eq!(
                            &wood
                                .graph()
                                .iter_nodes()
                                .map(|(_, _, spinney)| spinney.filter().clone())
                                .collect::<BTreeSet<_>>(),
                            expected,
                            "{}: {mode:?}, subtract_uv={subtract_uv}, vacuum_subtraction={vacuum_subtraction}",
                            graph.name
                        );
                        // The shared filter also supplies the legacy executor.
                        let mut spinneys = graph.classified_spinneys(
                            &graph.full_filter(),
                            &graph.model,
                            &settings.uv,
                            &graph.loop_momentum_basis,
                        )?;
                        graph.remove_full_observable_spinney(&mut spinneys, &settings.medium);
                        assert_eq!(
                            &spinneys
                                .iter()
                                .map(|spinney| spinney.filter().clone())
                                .collect::<BTreeSet<_>>(),
                            expected,
                            "legacy {}: {mode:?}, subtract_uv={subtract_uv}, vacuum_subtraction={vacuum_subtraction}",
                            graph.name
                        );
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn vacuum_subtraction_commutes_with_local_and_integrated_forests() -> Result<()> {
        use crate::{
            integrands::process::param_builder::ParamBuilderGraph,
            numerator::{
                soft_energy_certificate::{EnergyZeroOutcome, SoftEnergyAlgebra},
                symbolica_ext::NumeratorAtomExt,
            },
            settings::global::MediumMode,
        };
        use symbolica::atom::AtomView;
        use three_dimensional_reps::thermal::OCCUPATION_FORM;
        test_initialise()?;
        let graphs: Vec<Graph> = Graph::unit_numerator_fixtures(
            stringify!(
            digraph nested {
                edge [particle="scalar_1"]; node [num=1];
                v1 -> v2; v1 -> v2; v1 -> v2;
                v1 -> v1;
            }
            digraph disjoint {
                edge [particle="scalar_1"]; node [num=1];
                v1 -> v2; v1 -> v2;
                v3 -> v4; v3 -> v4;
                v2 -> v3; v1 -> v4;
            }),
            &crate::utils::load_generic_model("scalars"),
        )?;
        for graph in graphs {
            for medium in [
                MediumMode::ThermodynamicEquilibrium,
                MediumMode::ZeroTemperatureEquilibrium,
            ] {
                for generate_integrated in [false, true] {
                    let mut settings = GenerationSettings::default();
                    settings.medium.mode = medium;
                    settings.medium.vacuum_subtraction = true;
                    settings.uv.softct = false;
                    settings.uv.generate_integrated = generate_integrated;
                    let mut graph = graph.clone();
                    let options = graph.production_3d_expression_options(
                        &settings,
                        three_dimensional_reps::RepresentationMode::Cff,
                    )?;
                    let contract =
                        graph.paired_edges(&graph.tree_edges.subtract(&graph.initial_state_cut));
                    let production = graph
                        .generate_3d_expression_for_integrand(&contract, &None, &options, None)?;
                    let mut results = Vec::new();
                    for early_subtraction in [false, true] {
                        let projection = OrientationProjection::exact_expression(
                            &production,
                            &options,
                            &settings.orientation_pattern,
                            true,
                        )
                        .with_vacuum_subtraction(early_subtraction);
                        // Both routes use precisely the same proper forests. Only the
                        // boundary where the full vacuum projection acts is different.
                        let mut forests = Wood::new(
                            CutStructure::empty(&graph),
                            &graph,
                            &graph.model,
                            &settings,
                        )?
                        .unfold();
                        let model = graph.model.clone();
                        let mut context = Local4dProjectionContext::default();
                        forests.prepare_four_d(
                            &graph,
                            &model,
                            crate::utils::vakint()?,
                            &settings.uv,
                            &mut context,
                        )?;
                        forests.compute(&mut graph, projection, &settings.uv, &mut context)?;
                        drop(context);
                        results.push(forests.orientation_parametric_exprs(&graph, &settings.uv)?);
                    }
                    assert_eq!(results[0].len(), results[1].len());
                    for (late, early) in results[0].iter().zip(&results[1]) {
                        let late = late.integrands.resolved()?;
                        let early = early.integrands.resolved()?;
                        assert_eq!(late.iter().count(), early.iter().count());
                        for ((late_key, raw), (early_key, actual)) in late.iter().zip(early.iter())
                        {
                            assert_eq!(late_key, early_key);
                            let vacuum = graph.make_thermal_distributions_explicit(
                                &graph.model,
                                raw,
                                MediumMode::Vacuum,
                                graph.iter_edges().map(|(_, edge, _)| edge),
                            )?;
                            assert!(
                                raw.is_zero() || raw.contains_symbol(GS.thermal_weight_wrapper)
                            );
                            assert!(
                                actual.is_zero()
                                    || actual.contains_symbol(GS.thermal_weight_wrapper)
                            );
                            // Resolve only this scalar comparison. Final energy factors can
                            // group as E^-5 or E^-1/P^2, with E^2 = P, after the two routes.
                            let difference = (actual - (raw - vacuum))
                                .replace(function!(GS.thermal_weight_wrapper, W_.a_))
                                .with(W_.a_)
                                .replace_multiple(&*OCCUPATION_FORM)
                                .expand();
                            // Certify each occupation coefficient separately so unrelated
                            // energy denominators do not multiply the proof's term budget.
                            let mut occupations = BTreeSet::new();
                            difference.visitor(&mut |view| {
                                if let AtomView::Fun(call) = view
                                    && call.get_symbol() == GS.thermal_distribution
                                {
                                    occupations.insert(view.to_owned());
                                    return false;
                                }
                                true
                            });
                            let coefficients = difference
                                .coefficient_list_exact(&occupations.into_iter().collect_vec())
                                .expect("scalar difference must be polynomial in occupations");
                            for (occupation, coefficient) in coefficients {
                                let certificate =
                                    SoftEnergyAlgebra::default().scalar_zero(&coefficient)?;
                                assert!(
                                    matches!(certificate, EnergyZeroOutcome::Zero(_)),
                                    "{}: {medium:?}, integrated={generate_integrated}, residue={late_key:?}, occupation={occupation}: {certificate:?}: {coefficient}",
                                    graph.name
                                );
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn local_leaf_operations_follow_dependency_frontiers() -> Result<()> {
        test_initialise().unwrap();
        let dumbell: Graph = finalized_runtime_dot!(
            digraph G{
                graph [projector=1]
                node [num=1]
                edge [particle="scalar_1" num=1];
                v1 -> v2 [source="{ufo_order:0}" sink="{ufo_order:0}"];
                v2 -> v2 [lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:2}"];
                v1 -> v1 [lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:2}"];
                v1 -> v1 [lmb_id=2 source="{ufo_order:3}" sink="{ufo_order:4}"];
            },"scalars"
        )?;

        let forests = Wood::new(
            CutStructure::empty(&dumbell),
            &dumbell,
            &scalar_model(),
            &GenerationSettings::default(),
        )?
        .unfold();

        let frontier_label = |operation: &OperationNode| {
            operation
                .covers()
                .map_or_else(|| "∅".to_string(), |cover| cover.string_label())
        };
        let leaf_labels = |node| {
            forests
                .local_leaf_operations(node)
                .iter()
                .map(|leaf| {
                    (
                        leaf.op.order.string_label(),
                        frontier_label(&forests.graph()[leaf.frontier]),
                    )
                })
                .collect::<Vec<_>>()
        };

        let root_disconnected = forests
            .graph()
            .iter_nodes()
            .find_map(|(node, _, _)| {
                forests
                    .normalized_node_label(node)
                    .starts_with("{36,F}:")
                    .then_some(node)
            })
            .expect("lopsided dumbbell should contain a root disconnected frontier");
        assert_eq!(
            leaf_labels(root_disconnected),
            vec![
                ("36".to_string(), "∅".to_string()),
                ("F".to_string(), "∅".to_string())
            ]
        );
        let structured_dot = forests.dot_serialize();
        let mut join_markers = Vec::new();
        for (node, _, _) in forests.graph().iter_nodes() {
            for leaf in forests.local_leaf_operations(node) {
                let frontier_depth = forests.graph()[leaf.frontier].key.op_count();
                let (current, given) = forests
                    .forests()
                    .wood()
                    .current_given_pair(leaf.op.data, frontier_depth);
                if current.subgraph() == &leaf.op.order {
                    continue;
                }
                join_markers.push(
                    UvMarker::subgraph(current.subgraph(), given.subgraph())
                        .printer(SpensoPrintSettings::typst_options())
                        .to_string(),
                );
            }
        }
        assert!(
            !join_markers.is_empty(),
            "lopsided dumbbell must contain a disconnected join"
        );
        for marker in join_markers {
            assert!(
                structured_dot.contains(&marker),
                "structured DOT must use the runtime current/given marker {marker}"
            );
        }

        let dependent_disconnected = forests
            .graph()
            .iter_nodes()
            .find_map(|(node, _, _)| {
                forests
                    .normalized_node_label(node)
                    .starts_with("{C} · {36,F}:")
                    .then_some(node)
            })
            .expect("lopsided dumbbell should contain a C-dependent disconnected frontier");
        assert_eq!(
            leaf_labels(dependent_disconnected),
            vec![
                ("F".to_string(), "C".to_string()),
                ("36".to_string(), "∅".to_string())
            ]
        );
        let replay_states = forests
            .forests()
            .union_replay_states(dependent_disconnected)?;
        assert_eq!(replay_states.len(), 6);
        assert_eq!(
            replay_states
                .iter()
                .filter(|state| !state.local_edges.is_empty())
                .count(),
            5
        );
        Ok(())
    }

    #[test]
    fn union_terms_project_factorized_typed_4d_values() -> Result<()> {
        test_initialise().unwrap();
        let graph: Graph = finalized_runtime_dot!(digraph G {
          projector = "1";
          edge	 [
        num = "1𝑖", particle = "scalar_1"]
          node	 [
        num = "UFO::SCALAR_COUPLING"]
          v1;
          v2;

          v1:0	-> v1:1	 [id=0 source="{ufo_order:0}" sink="{ufo_order:1}"  lmb_id="0"];
          v1:2	-> v1:3	 [id=1 source="{ufo_order:2}" sink="{ufo_order:3}"  lmb_id="1"];
          v1:4	-> v2:5	 [id=2 source="{ufo_order:4}" sink="{ufo_order:0}" ];
          v2:6	-> v2:7	 [id=3 source="{ufo_order:1}" sink="{ufo_order:2}"  lmb_id="2"];
        },"scalars"
                )?;
        // Disable integrated terms so each union equality isolates factorized
        // local 4D composition.
        let settings = GenerationSettings {
            uv: UVgenerationSettings {
                generate_integrated: false,
                local_uv_cts_from_expanded_4d_integrands: true,
                ..Default::default()
            },
            ..Default::default()
        };
        let cut_structure = CutStructure::empty(&graph);
        let mut forests = Wood::new(cut_structure, &graph, &graph.model, &settings)?.unfold();
        forests.integrate(&graph, &graph.model, crate::utils::vakint()?, &settings.uv)?;

        let unions = forests
            .graph()
            .iter_nodes()
            .filter_map(|(node, _, operation)| {
                (forests.graph().is_disjoint_union(node) && !operation.key.is_empty())
                    .then_some(node)
            })
            .collect::<Vec<_>>();
        assert!(!unions.is_empty());

        for union in unions {
            let components = forests.forests().disconnected_component_nodes(union)?;
            let expected = Local4dCts::from_full_product(
                components
                    .into_iter()
                    .map(|component| forests.recursion_input_4d(component))
                    .collect::<Result<Vec<_>>>()?,
            )?;
            let operation = &forests.graph()[union];
            let local = forests
                .compute_store
                .require(operation)?
                .local_4d(operation)?;

            // Compare the complete typed union expression; no per-cut parent cache
            // or component-path replay participates.
            assert!(
                (local.atom().expression().collect_factors()
                    - expected.atom().expression().collect_factors())
                .is_zero()
            );
        }
        Ok(())
    }

    #[test]
    fn union_terms_replay_component_paths_from_typed_roots() -> Result<()> {
        use crate::graph::feynman_graph::FeynmanGraph;

        test_initialise().unwrap();
        let mut graph: Graph = finalized_runtime_dot!(
            digraph G{
                graph [projector=1]
                node [num=1]
                edge [particle="scalar_1" num=1];
                v1 -> v2 [source="{ufo_order:0}" sink="{ufo_order:0}"];
                v2 -> v2 [lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:2}"];
                v1 -> v1 [lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:2}"];
                v1 -> v1 [lmb_id=2 source="{ufo_order:3}" sink="{ufo_order:4}"];
            },"scalars"
        )?;
        let settings = GenerationSettings::default();
        let cut_structure = CutStructure::empty(&graph);
        let cutset = cut_structure
            .cuts
            .first()
            .expect("empty cut structure has one cut")
            .clone();
        let mut forests = Wood::new(cut_structure, &graph, &scalar_model(), &settings)?.unfold();

        let root_disconnected = forests
            .graph()
            .iter_nodes()
            .find_map(|(node, _, _)| {
                forests
                    .normalized_node_label(node)
                    .starts_with("{36,F}:")
                    .then_some(node)
            })
            .expect("lopsided dumbbell should contain a root disconnected frontier");
        let dependent_disconnected = forests
            .graph()
            .iter_nodes()
            .find_map(|(node, _, _)| {
                forests
                    .normalized_node_label(node)
                    .starts_with("{C} · {36,F}:")
                    .then_some(node)
            })
            .expect("lopsided dumbbell should contain a C-dependent disconnected frontier");
        let non_root_frontier = forests
            .local_leaf_operations(dependent_disconnected)
            .first()
            .expect("dependent disconnected operation has a leaf")
            .frontier;
        assert_eq!(
            forests.graph()[non_root_frontier]
                .covers()
                .expect("non-root frontier has a cover")
                .string_label(),
            "C"
        );

        let operations = forests
            .graph()
            .iter_nodes()
            .map(|(_, _, operation)| operation.clone())
            .collect::<Vec<_>>();
        // Zero integrated terms isolate component-path replay in the local accumulator.
        for operation in operations {
            forests
                .compute_store
                .entry(operation)
                .or_default()
                .integrated = Some(IntegratedCts::root());
        }

        let orientation_pattern = crate::settings::global::OrientationPattern::default();
        let options = graph.denominator_only_cff_3d_expression_options();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let contract_edges = graph
            .iter_edges_of(&graph.tree_edges)
            .map(|(_, edge, _)| edge)
            .collect_vec();
        let production = graph.generate_3d_expression_for_integrand(
            &contract_edges,
            &canonization,
            &options,
            Some(&Atom::one()),
        )?;
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact_expression(
                &production,
                &options,
                &orientation_pattern,
                false,
            ),
        );
        let seed = forests.local_3d_for_node(
            forests.forests().root(),
            &mut graph,
            &cutset,
            localizer,
            &settings.uv,
            &mut Local4dProjectionContext::default(),
        )?;
        let root_store_marker = symbolica::symbol!("root_store_marker");
        let frontier_store_marker = symbolica::symbol!("frontier_store_marker");
        // Mark cached local terms so the union result proves that replay starts from typed roots.
        let marked_computation = |marker| -> Result<CutComputation> {
            Ok(CutComputation {
                local_3d: seed.local_3d.map(|_| Ok(Atom::var(marker)))?,
                final_integrands: seed.final_integrands.clone(),
            })
        };

        let root_operation = forests.graph()[forests.forests().root()].clone();
        forests
            .compute_store
            .entry(root_operation)
            .or_default()
            .cuts
            .insert(cutset.clone(), marked_computation(root_store_marker)?);
        forests
            .compute_store
            .entry(forests.graph()[non_root_frontier].clone())
            .or_default()
            .cuts
            .insert(cutset.clone(), marked_computation(frontier_store_marker)?);

        let root_result = forests.local_3d_for_node(
            root_disconnected,
            &mut graph,
            &cutset,
            localizer,
            &settings.uv,
            &mut Local4dProjectionContext::default(),
        )?;
        assert!(
            root_result
                .local_3d
                .direct()?
                .branches()?
                .iter_keys()
                .all(|(_, integrands)| integrands
                    .iter()
                    .all(|(_, term)| !term.contains_symbol(root_store_marker))),
            "an empty dependency frontier must start from the typed root"
        );

        let frontier_result = forests.local_3d_for_node(
            dependent_disconnected,
            &mut graph,
            &cutset,
            localizer,
            &settings.uv,
            &mut Local4dProjectionContext::default(),
        )?;
        let replay_states = forests
            .forests()
            .union_replay_states(dependent_disconnected)?;
        let state = replay_states
            .iter()
            .find(|state| {
                state.local_edges.len() == 2
                    && forests
                        .normalized_node_label(state.integrated)
                        .starts_with("{36}:")
            })
            .expect("the integrated 36 prefix has the two-step C-to-F suffix");
        let expected_active = state
            .local_edges
            .iter()
            .enumerate()
            .map(|(offset, edge)| {
                let step_order = forests.graph()[state.integrated].key.op_count() + offset;
                let (current, given) = forests
                    .forests()
                    .wood()
                    .current_given_pair(*edge, step_order);
                current.reduced_subgraph(&given)
            })
            .reduce(|active, reduced| active.union(&reduced))
            .expect("the selected replay state has a local suffix");
        assert!(
            expected_active.empty_intersection(forests.source_spinney(state.integrated).filter())
        );
        let sector = frontier_result
            .local_3d
            .direct()?
            .sectors()?
            .iter()
            .find(|sector| sector.active_subgraph == expected_active)
            .expect("the replay frontier contains the expected active sector");
        assert!(
            !sector.projection_paths.is_empty()
                && sector.projection_paths.iter().all(|path| {
                    path.steps.len() == state.local_edges.len()
                        && path
                            .steps
                            .windows(2)
                            .all(|pair| pair[0].topo_order <= pair[1].topo_order)
                }),
            "each disconnected sector must retain its own ordered local replay path"
        );
        assert!(
            frontier_result
                .local_3d
                .direct()?
                .branches()?
                .iter_keys()
                .all(|(_, integrands)| integrands.iter().all(|(_, term)| {
                    !term.contains_symbol(root_store_marker)
                        && !term.contains_symbol(frontier_store_marker)
                })),
            "a union must replay every component path from its typed root"
        );
        Ok(())
    }

    #[test]
    fn disconnected_soft_union_is_the_product_of_completed_component_atoms() -> Result<()> {
        test_initialise().unwrap();
        let mut graph: Graph = Graph::import_dot(
            include_str!("../../../../tests/resources/graphs/local_ir_dod2_scalar_spectacles.dot"),
            &crate::utils::load_generic_model("scalars"),
        )
        .map(|mut graphs| {
            assert_eq!(graphs.len(), 1, "one fixture graph expected");
            graphs.remove(0)
        })?;
        let settings = GenerationSettings {
            uv: UVgenerationSettings {
                generate_integrated: false,
                renormalization_prescription: RenormalizationPrescriptionSettings {
                    log_divergent: ApproximationType::IR,
                    massive_power_divergent: ApproximationType::IR,
                    massless_power_divergent: ApproximationType::IR,
                    ..Default::default()
                },
                ..Default::default()
            },
            ..Default::default()
        };
        let mut forests =
            Wood::new(CutStructure::empty(&graph), &graph, &graph.model, &settings)?.unfold();
        forests.integrate(&graph, &graph.model, crate::utils::vakint()?, &settings.uv)?;

        let union = forests
            .graph()
            .iter_nodes()
            .find_map(|(node, _, _)| forests.graph().is_disjoint_union(node).then_some(node))
            .expect("the spectacles fixture must contain a disconnected UV union");
        assert_eq!(
            forests.source_spinney(union).renormalization_scheme(),
            ApproximationType::MUV,
            "a disconnected union carries only the neutral aggregate scheme"
        );
        let components = forests.forests().disconnected_component_nodes(union)?;
        assert_eq!(components.len(), 2);
        assert!(components.iter().all(|component| {
            forests.source_spinney(*component).renormalization_scheme() == ApproximationType::IR
        }));

        let expected =
            components
                .iter()
                .try_fold(Atom::one(), |product, component| -> Result<Atom> {
                    let completed = forests.recursion_input_4d(*component)?;
                    assert!(
                        !completed.atom().is_zero(),
                        "each disconnected soft component must be completed before the join"
                    );
                    Ok(product * completed.atom().expression())
                })?;
        let operation = &forests.graph()[union];
        let actual = forests
            .compute_store
            .require(operation)?
            .local_4d(operation)?;
        assert_eq!(actual.atom().expression(), &expected);
        let provenance = forests.node_export_4d_provenance(union, actual);
        assert_eq!(provenance.node.parent_keys.len(), 2);
        let UVForestLocalProvenance::FourD {
            components: provenance_components,
            branches,
        } = &provenance.node.local
        else {
            panic!("a computed 4D union must carry 4D provenance")
        };
        assert_eq!(provenance_components.len(), 2);
        assert!(provenance_components.iter().all(|component| {
            component.scheme == ApproximationType::IR
                && !component.route_loop_edges.is_empty()
                && !component.route_signatures.is_empty()
        }));
        assert_eq!(branches.len(), 1);
        assert_eq!(branches[0].branch, Local4dBranch::Combined);
        assert!(branches[0].byte_size > 0);

        // Exercise the same factorization directly in one parent-CFF
        // orientation. With integrated generation disabled, only the replay
        // sector containing both local component operations may survive at
        // the disconnected-union node.
        let cutset = CutSet::empty(graph.n_hedges());
        let orientation_pattern = crate::settings::global::OrientationPattern::default();
        let options = graph.production_3d_expression_options(
            &Default::default(),
            three_dimensional_reps::RepresentationMode::Cff,
        )?;
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        // Match the production root: its structural bridge remains a 4D
        // denominator outside the CFF and its loop-energy contour.
        let contract_edges = graph
            .iter_edges_of(
                &graph
                    .tree_edges
                    .subtract(&graph.initial_state_cut)
                    .subtract(&graph.external_filter::<SuBitGraph>()),
            )
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        let mut production = graph.generate_3d_expression_for_integrand(
            &contract_edges,
            &canonization,
            &options,
            None,
        )?;
        assert!(
            !production.expression.orientations.is_empty(),
            "the spectacles parent CFF has at least one orientation"
        );
        production.expression.orientations.truncate(1);
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact_expression(
                &production,
                &options,
                &orientation_pattern,
                false,
            ),
        );
        let union_3d = forests.local_3d_for_node(
            union,
            &mut graph,
            &cutset,
            localizer,
            &settings.uv,
            &mut Local4dProjectionContext::default(),
        )?;
        let sectors = union_3d
            .local_3d
            .direct()?
            .sectors()
            .expect("a disconnected local 3D union retains its replay sectors");
        assert_eq!(sectors.len(), 3);
        let combined_sectors = sectors
            .iter()
            .map(DirectSector::combine)
            .collect::<Result<Vec<_>>>()?;
        let nonzero_sectors = combined_sectors
            .iter()
            .filter(|branches| {
                branches
                    .iter_keys()
                    .any(|(_, integrands)| integrands.iter().any(|(_, atom)| !atom.is_zero()))
            })
            .collect::<Vec<_>>();
        assert_eq!(
            nonzero_sectors.len(),
            1,
            "zero integrated prefixes must leave only the two-local-component cross sector"
        );
        let cross = nonzero_sectors[0];

        let assert_same_branches = |left: &DirectResidueBranches,
                                    right: &DirectResidueBranches,
                                    label: &str|
         -> Result<()> {
            let left = left.iter_keys().collect::<Vec<_>>();
            let right = right.iter_keys().collect::<Vec<_>>();
            assert_eq!(
                left.len(),
                right.len(),
                "{label} differs in residue-map count"
            );
            for ((left_key, left), (right_key, right)) in left.into_iter().zip(right) {
                assert_eq!(
                    left_key, right_key,
                    "{label} changed the complete residue map"
                );
                // Independently replayed forests have different numerator-family
                // scopes. Compare their values using the retained definitions.
                let left = left.resolved()?;
                let right = right.resolved()?;
                left.checked_zip(&right, |residue, left, right| {
                    // Compare the exact factorized forest values before selectors
                    // are materialized or any graph numerator is distributed.
                    let difference = (left - right).together().cancel();
                    assert!(
                        difference.is_zero(),
                        "{label} differs in residue {residue:?}: {difference}"
                    );
                    Ok(Atom::Zero)
                })?;
            }
            Ok(())
        };
        assert_same_branches(
            &union_3d.local_3d.direct()?.branches()?,
            cross,
            "the disconnected aggregate and its cross sector",
        )?;

        let component_edges = components
            .iter()
            .map(|component| {
                let (parent, edge) = forests
                    .graph()
                    .unique_parent(*component)
                    .expect("a primitive disconnected component has one root edge");
                assert!(
                    forests.graph()[parent].key.is_empty(),
                    "the spectacles components must both be primitive"
                );
                edge
            })
            .collect::<Vec<_>>();
        let mut replay = |edges: &[EdgeIndex]| -> Result<Direct3dCts> {
            let mut local = Direct3dCts::root(&graph, localizer)?;
            for (order, edge) in edges.iter().enumerate() {
                let (current, given) = forests.forests().wood().current_given_pair(*edge, order);
                local = Direct3dApproximation::new(localizer, &mut graph, &settings.uv)
                    .run_local(&local, &current, &given, &current, &given)?;
            }
            assert!(local.branches()?.iter_keys().all(|(_, integrands)| {
                integrands
                    .iter()
                    .all(|(_, atom)| !atom.contains_symbol(GS.rescale))
            }));
            Ok(local)
        };
        let forward = replay(&component_edges)?;
        let reverse = replay(&component_edges.iter().rev().copied().collect::<Vec<_>>())?;
        assert_same_branches(
            &forward.branches()?,
            cross,
            "forward component-local replay and the produced cross sector",
        )?;
        assert_same_branches(
            &reverse.branches()?,
            cross,
            "reverse component-local replay and the produced cross sector",
        )?;
        Ok(())
    }

    #[test]
    fn nested_soft_top_self_energy_preserves_completed_parent_poles() -> Result<()> {
        initialisation::test_initialise()?;
        let model = utils::load_generic_model("sm");
        let graph = Graph::import_dot(
            include_str!(
                "../../../../tests/resources/graphs/local_nested_soft_top_self_energy.dot"
            ),
            &model,
        )?
        .remove(0);
        let external: SuBitGraph = graph.external_filter();
        let parent = graph.full_filter().subtract(&external);
        assert_eq!(graph.n_loops(&parent), 2);
        let physical_mass = graph[EdgeIndex::from(2)].mass_atom(&model);
        assert!(!physical_mass.is_zero());
        assert_ne!(physical_mass, Atom::var(GS.m_uv_vacuum));

        let mut completed_poles = Vec::with_capacity(2);
        for scheme in [ApproximationType::MUV, ApproximationType::IR] {
            let settings = GenerationSettings {
                uv: UVgenerationSettings {
                    generate_integrated: true,
                    renormalization_prescription: RenormalizationPrescriptionSettings {
                        log_divergent: scheme,
                        massive_power_divergent: scheme,
                        massless_power_divergent: scheme,
                        ..Default::default()
                    },
                    ..Default::default()
                },
                ..Default::default()
            };
            let mut forests =
                Wood::new(CutStructure::empty(&graph), &graph, &model, &settings)?.unfold();
            forests.integrate(&graph, &model, utils::vakint()?, &settings.uv)?;
            let mut pole = Atom::Zero;
            let mut bare_parent_count = 0;
            let mut nested_parent_count = 0;
            let mut finite_prefix_retained = false;
            for (_, _, operation) in forests.graph().iter_nodes() {
                if operation.covers().as_ref() != Some(&parent) {
                    continue;
                }
                let computed = forests.compute_store.require(operation)?;
                pole += computed
                    .integrated(operation)?
                    .physical_pole()?
                    .into_expression();
                if operation.key.op_count() == 1 {
                    bare_parent_count += 1;
                } else {
                    nested_parent_count += 1;
                    finite_prefix_retained |=
                        computed.local_4d(operation)?.active_sectors().len() > 1;
                }
            }
            assert_eq!(bare_parent_count, 1);
            assert!(nested_parent_count > 0);
            assert!(
                finite_prefix_retained,
                "the nested {scheme} parent must include the integrated child's finite completion"
            );
            // Compare the complete parent forest, keeping its Feynman-rule numerator and
            // finite-prefix epsilon jets. Resolve the model's color representation
            // invariants before comparing; individual nested terms need not agree.
            let pole = SymbolicTensor::infer_indexed(
                pole.replace(GS.m_uv_expansion)
                    .with(GS.m_uv_vacuum)
                    .replace(GS.dim)
                    .with(4),
            )?
            .simplify_algebra(&AlgebraSettings {
                color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
                contract: AlgebraContraction::Minimal,
                ..Default::default()
            })?
            .to_dots()?
            .canonize()?
            .into_expression()
            .collect_factors();
            assert!(
                !pole.is_zero(),
                "the completed {scheme} parent pole must be nonzero"
            );
            completed_poles.push(pole);
        }
        let difference = (&completed_poles[0] - &completed_poles[1])
            .together()
            .cancel();
        assert!(
            difference.is_zero(),
            "MUV and H must give the same completed parent poles: {difference}"
        );
        Ok(())
    }

    #[test]
    fn heterogeneous_terminal_projects_components_before_multiplying() -> Result<()> {
        // Symbolica's exact integrated projection is intentionally recursive and
        // needs more than libtest's 2 MiB worker stack. GammaLoop's runtime main
        // thread has the platform stack, so exercise the same operation here.
        std::thread::Builder::new()
            .name("heterogeneous-terminal-projection".to_string())
            .stack_size(8 * 1024 * 1024)
            .spawn(heterogeneous_terminal_projects_components_before_multiplying_impl)
            .expect("projection test thread should start")
            .join()
            .expect("projection test thread should not panic")
    }

    fn heterogeneous_terminal_projects_components_before_multiplying_impl() -> Result<()> {
        test_initialise().unwrap();
        let graph: Graph = finalized_runtime_dot!(
            digraph G {
                num = "1";
                projector = "1";
                overall_factor = "1";
                edge [particle = "scalar_1" num=1];
                node [num = "1"];
                v1 -> v2 [source="{ufo_order:0}" sink="{ufo_order:0}"];
                v1 -> v1 [lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:2}"];
                v2 -> v3 [source="{ufo_order:1}" sink="{ufo_order:0}"];
                v2 -> v3 [lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:1}"];
            },
            "scalars"
        )?;
        let settings = GenerationSettings {
            uv: UVgenerationSettings {
                softct: false,
                renormalization_prescription: RenormalizationPrescriptionSettings {
                    log_divergent: ApproximationType::MUV,
                    massive_power_divergent: ApproximationType::PolePart,
                    massless_power_divergent: ApproximationType::PolePart,
                    ..Default::default()
                },
                ..Default::default()
            },
            ..Default::default()
        };
        let model = scalar_model();
        let mut forests =
            Wood::new(CutStructure::empty(&graph), &graph, &model, &settings)?.unfold();
        forests.integrate(&graph, &model, crate::utils::vakint()?, &settings.uv)?;

        let terminals = forests
            .graph()
            .iter_nodes()
            .filter_map(|(node, mut crown, _)| {
                (!crown.any(|hedge| forests.graph().flow(hedge).is_source())).then_some(node)
            })
            .collect::<Vec<_>>();
        assert_eq!(terminals.len(), 1);
        let terminal = terminals[0];
        assert!(forests.graph().is_disjoint_union(terminal));

        let components = forests.forests().disconnected_component_nodes(terminal)?;
        assert_eq!(components.len(), 2);
        let schemes = components
            .iter()
            .map(|component| forests.source_spinney(*component).renormalization_scheme())
            .collect::<Vec<_>>();
        assert!(schemes.contains(&ApproximationType::MUV));
        assert!(schemes.contains(&ApproximationType::PolePart));

        let expected = components.iter().try_fold(
            SymbolicTensor::infer_indexed(Atom::one())?,
            |product, &component| -> Result<_> {
                let component_key = &forests.graph()[component];
                let integrated = forests
                    .compute_store
                    .require(component_key)?
                    .integrated(component_key)?;
                let projection = match forests.source_spinney(component).renormalization_scheme() {
                    ApproximationType::MUV => integrated.physical_finite_counterterm()?,
                    ApproximationType::PolePart => integrated.physical_pole()?,
                    scheme => panic!("unexpected component scheme {scheme}"),
                };
                Ok(product.multiply(projection)?)
            },
        )?;
        let terminal_key = &forests.graph()[terminal];
        let aggregate = forests
            .compute_store
            .require(terminal_key)?
            .integrated(terminal_key)?;
        assert_ne!(
            expected.expression().expand(),
            aggregate
                .physical_finite_counterterm()?
                .expression()
                .expand()
        );
        assert_ne!(
            expected.expression().expand(),
            aggregate.physical_pole()?.expression().expand()
        );

        let wild = Atom::var(W_.x___);
        let replacements =
            graph.integrand_replacement(&graph.full_filter(), &graph.loop_momentum_basis, &[wild]);
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let expected = expected
            .simplify_algebra(&idenso::tensor::AlgebraSettings {
                color: Some(ColorSimplifySettings::default()),
                ..Default::default()
            })?
            .contract(idenso::tensor::ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })?
            .to_dots()?;
        let expected = cooking
            .uncook(expected.expression().as_view())
            .replace_multiple(&replacements)
            .replace(GS.m_uv_expansion)
            .with(GS.m_uv_vacuum);
        let actual = forests
            .renormalization_part_of_ends(&graph, &settings.uv)?
            .expression;
        assert_eq!(actual.expand(), expected.expand());

        Ok(())
    }

    #[test]
    fn triple_tadpole() -> Result<()> {
        test_initialise().unwrap();
        let dumbell: Graph = finalized_runtime_dot!(
            digraph G{
                graph [projector=1]
                node [num=1]
                edge [particle="scalar_1" num=1];
                v1 -> v2 [source="{ufo_order:0}" sink="{ufo_order:0}"];
                v2 -> v3 [source="{ufo_order:1}" sink="{ufo_order:0}"];
                v3 -> v3 [lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:2}"];
                v2 -> v2 [lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:3}"];
                v1 -> v1 [lmb_id=2 source="{ufo_order:1}" sink="{ufo_order:2}"];
            },"scalars"
        )?;

        let spinneys: Vec<_> = dumbell
            .spinneys(&dumbell.full_filter())
            .unwrap()
            .into_iter()
            .filter_map(|a| Spinney::new(a, &dumbell, &dumbell.loop_momentum_basis).transpose())
            .collect::<Result<_, _>>()?;
        let f = Wood::new(
            CutStructure::empty(&dumbell),
            &dumbell,
            &scalar_model(),
            &GenerationSettings::default(),
        )?;

        println!("{}", f);

        insta::assert_snapshot!(
            f.graph().n_nodes(),
            @"8",
            // format!("Wood does not have correct number of spinneys: \n{}",f)
        );

        let three_component_target = f
            .graph()
            .iter_nodes()
            .find_map(|(node, _, spinney)| (spinney.n_components() == 3).then_some(node))
            .expect("triple tadpole should contain its fully factorized target");
        let target_filter = f.graph()[three_component_target].filter().clone();
        let components = dumbell
            .as_ref()
            .connected_components(&target_filter)
            .into_iter()
            .collect::<BTreeSet<_>>();
        assert_eq!(components.len(), 3);
        let incoming = f
            .graph()
            .iter_crown(three_component_target)
            .filter(|hedge| f.graph().flow(*hedge) == Flow::Sink)
            .map(|hedge| f.graph()[f.graph()[&hedge]].clone())
            .collect::<BTreeSet<_>>();
        assert_eq!(incoming, components);
        assert!(
            !incoming.contains(&target_filter),
            "a fully supplied product must not retain a duplicate atomic full-target transition"
        );

        for (_, _, d) in f.graph().iter_nodes() {
            println!(
                "//Node {}: \n{}",
                d.subgraph().string_label(),
                dumbell.dot(d.subgraph())
            );
        }
        let _ff = OldWood::from_spinneys(spinneys, &dumbell); //.unfold(&g, &g.loop_momentum_basis);

        // println!("{}", ff.dot(&dumbell));

        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
            f.graph().n_nodes(),
            @"8");

        let three_component_union = f
            .graph()
            .iter_nodes()
            .find_map(|(node, _, operation)| {
                (operation.key.op_count() == 3
                    && f.forests()
                        .wood()
                        .join_factors(f.graph().source_node(node))
                        .is_some_and(|factors| factors.len() == 3))
                .then_some(node)
            })
            .expect("triple tadpole should contain a three-component union");
        assert_eq!(
            f.forests()
                .disconnected_component_nodes(three_component_union)?
                .len(),
            3
        );
        let replay_states = f.forests().union_replay_states(three_component_union)?;
        assert_eq!(replay_states.len(), 8);
        assert_eq!(
            replay_states
                .iter()
                .filter(|state| !state.local_edges.is_empty())
                .count(),
            7
        );

        Ok(())
    }

    #[test]
    fn saclay() -> Result<()> {
        test_initialise().unwrap();
        let dt: Graph = finalized_runtime_dot!(digraph GL16{

        num = "1"
                projector = "1"
                node [num=1]
                ext	 [style=invis];
                ext	-> 0  [dir=none id=0 particle="g" num=1 dod=-2 sink="{ufo_order:0}"];
                2	-> ext  [dir=none id=1 particle="g" num=1 dod=-2 source="{ufo_order:0}"];
                0 -> 1 [id=2 particle="d" num="Q(2,spenso::cind(1))" dod=-1 source="{ufo_order:1}" sink="{ufo_order:0}"];
                1 -> 2 [id=3 particle="d" num="Q(3,spenso::cind(1))" dod=-1 source="{ufo_order:1}" sink="{ufo_order:1}"];
                2 -> 3 [id=4 particle="d" num="Q(4,spenso::cind(1))" dod=-1 source="{ufo_order:2}" sink="{ufo_order:0}"];
                3 -> 0 [id=5 particle="d" num="Q(5,spenso::cind(1))" dod=-1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:2}"];
                1 -> 3 [id=6 particle="g" num=1 dod=-2 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:2}"]
                    })?;

        let f = Wood::new(
            CutStructure::empty(&dt),
            &dt,
            &standard_model(),
            &GenerationSettings::default(),
        )?;

        println!("{}", dt.dot_serialize(&DotExportSettings::default()));

        insta::assert_snapshot!(
            f.graph().n_nodes(),
            @"5",
            // format!("Wood does not have correct number of spinneys: \n{}",f)
        );

        let f = f.unfold();
        let structured_dot = f.dot_serialize();
        assert!(structured_dot.contains(r#"label="K[#expr S_44⊛0]""#));
        assert!(!structured_dot.contains("#T("));
        println!("{structured_dot}");
        insta::assert_snapshot!(
            f.graph().n_nodes(),
            @"8");

        Ok(())
    }

    #[test]
    fn dumbells() -> Result<()> {
        test_initialise().unwrap();
        let dumbell: Graph = finalized_runtime_dot!(
            digraph G{
                graph [projector=1]
                node [num=1]
                edge [particle="scalar_1" num=1];
                v1 -> v2 [source="{ufo_order:0}" sink="{ufo_order:0}"];
                v2 -> v2 [lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:2}"];
                v1 -> v1 [lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:2}"];
            },"scalars"
        )?;

        let spinneys: Vec<_> = dumbell
            .spinneys(&dumbell.full_filter())
            .unwrap()
            .into_iter()
            .filter_map(|a| Spinney::new(a, &dumbell, &dumbell.loop_momentum_basis).transpose())
            .collect::<Result<_, _>>()?;
        let f = Wood::new(
            CutStructure::empty(&dumbell),
            &dumbell,
            &scalar_model(),
            &GenerationSettings::default(),
        )?;

        println!("{}", f);

        insta::assert_snapshot!(
            f.graph().n_nodes(),
            @"4",
            // format!("Wood does not have correct number of spinneys: \n{}",f)
        );

        for (_, _, d) in f.graph().iter_nodes() {
            println!(
                "//Node {}: \n{}",
                d.subgraph().string_label(),
                dumbell.dot(d.subgraph())
            );
        }
        let _ff = OldWood::from_spinneys(spinneys, &dumbell); //.unfold(&g, &g.loop_momentum_basis);

        // println!("{}", ff.dot(&dumbell));

        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
            f.graph().n_nodes(),
            @"4");

        Ok(())
    }

    #[test]
    fn bugblatter() -> Result<()> {
        test_initialise().unwrap();

        match finalized_runtime_dot!(
            digraph G{
                graph [projector=1]
                node [num=1]
                A1 -> A2 [id=0 particle="t" num="Q(0,spenso::cind(1))" dod=-1 source="{ufo_order:0}" sink="{ufo_order:0}"];
                A2 -> A3 [id=1 particle="t" num="Q(1,spenso::cind(1))" dod=-1 source="{ufo_order:1}" sink="{ufo_order:0}"];
                A3 -> A1 [id=2 particle="t" num="Q(2,spenso::cind(1))" dod=-1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                B1 -> B2 [id=3 particle="t" num="Q(3,spenso::cind(1))" dod=-1 source="{ufo_order:0}" sink="{ufo_order:0}"];
                B2 -> B3 [id=4 particle="t" num="Q(4,spenso::cind(1))" dod=-1 source="{ufo_order:1}" sink="{ufo_order:0}"];
                B3 -> B1 [id=5 particle="t" num="Q(5,spenso::cind(1))" dod=-1 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A1 -> B1 [id=6 particle="a" num=1 dod=-2 source="{ufo_order:2}" sink="{ufo_order:2}"];
                A2 -> B2 [id=7 particle="a" num=1 dod=-2 lmb_id=2 source="{ufo_order:2}" sink="{ufo_order:2}"];
                A3 -> B3 [id=8 particle="a" num=1 dod=-2 lmb_id=3 source="{ufo_order:2}" sink="{ufo_order:2}"];
            },"sm"
        ) {
            Ok(g) => {
                let g: Graph = g;
                let spinneys: Vec<_> = g
                    .spinneys(&g.full_filter())
                    .unwrap()
                    .into_iter()
                    .filter_map(|a| Spinney::new(a, &g, &g.loop_momentum_basis).transpose())
                    .collect::<Result<_, _>>()?;
                let f = Wood::new(
                    CutStructure::empty(&g),
                    &g,
                    &standard_model(),
                    &GenerationSettings::default(),
                )?;

                println!("{}", f);

                assert_eq!(
                    20,
                    f.graph().n_nodes(),
                    "Wood does not have correct number of spinneys: \n{}",
                    f
                );

                for (_, _, d) in f.graph().iter_nodes() {
                    println!(
                        "//Node {}: \n{}",
                        d.subgraph().string_label(),
                        g.dot_lmb_of(d.subgraph(), d.lmb())
                    );
                }
                let _ff = OldWood::from_spinneys(spinneys, &g); //.unfold(&g, &g.loop_momentum_basis);

                // println!("{}", ff.dot(&g));

                let f = f.unfold();
                f.debug_walk();
                println!("{}", f);
                assert_eq!(
                    152,
                    f.graph().n_nodes(),
                    "Forest unfolds into the wrong number of terms :\n{}",
                    f
                );

                // println!("{}", f)
            }
            Err(e) => {
                eprintln!("{}", e);
            }
        }

        // let f = SpinneyWood::from_spinneys(
        //     g.spinneys(&g.full_filter()).into_iter().map(|a| a.filter),
        //     &g,
        // )
        // .unfold();

        // println!("{}", f.graph().base_dot());
        Ok(())
    }

    #[test]
    fn mercedes() -> Result<()> {
        test_initialise().unwrap();

        let mercedes: Graph = finalized_runtime_dot!(
            digraph G{
                graph [projector=1]
                node [num=1]
                edge [particle="scalar_1" num=1];
                v1 -> v2 [source="{ufo_order:0}" sink="{ufo_order:0}"];
                v2 -> v3 [source="{ufo_order:1}" sink="{ufo_order:0}"];
                v3 -> v1 [lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                v1 -> v4 [source="{ufo_order:2}" sink="{ufo_order:0}"];
                v2 -> v4 [lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:1}"];
                v3 -> v4 [lmb_id=2 source="{ufo_order:2}" sink="{ufo_order:2}"];
            },"scalars"
        )?;

        let f = Wood::new(
            CutStructure::empty(&mercedes),
            &mercedes,
            &scalar_model(),
            &GenerationSettings::default(),
        )?;
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph().n_nodes(),
        @"2",
        );
        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph().n_nodes(),
        @"2",
         );

        Ok(())
    }

    #[test]
    fn sunrise() -> Result<()> {
        test_initialise().unwrap();

        let sunrise: Graph = finalized_runtime_dot!( digraph sunrise{
            graph [projector=1]
            node [num = "1"]
            edge [particle=scalar_1 num=1]
            e        [style=invis]
            e -> A:0 [id=3 sink="{ufo_order:0}"]
            B:1 -> e [id=4 source="{ufo_order:0}"]
            A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:1}"]
            A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"]
            A -> B [id=2 lmb_id=1 source="{ufo_order:3}" sink="{ufo_order:3}"]
        },"scalars")?;
        // let spinneys = spectacles.spinneys(&spectacles.full_filter());
        let f = Wood::new(
            CutStructure::empty(&sunrise),
            &sunrise,
            &scalar_model(),
            &GenerationSettings::default(),
        )?;
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph().n_nodes(),
        @"5",
        );
        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph().n_nodes(),
        @"8",
         );

        Ok(())
    }
    #[test]
    fn dotted_sunrise() -> Result<()> {
        test_initialise().unwrap();

        let sunrise: Graph = finalized_runtime_dot!( digraph sunrise{
            graph [projector=1]
            node [num=1]
            edge [particle=scalar_1 num=1]
            e        [style=invis]
            e -> A:0 [id=3 sink="{ufo_order:0}"]
            B:1 -> e [id=4 source="{ufo_order:0}"]

            A -> C [id=0 source="{ufo_order:1}" sink="{ufo_order:0}"]
            C -> e [source="{ufo_order:1}"]
            C -> B [source="{ufo_order:2}" sink="{ufo_order:1}"]
            A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"]
            A -> B [id=2 lmb_id=1 source="{ufo_order:3}" sink="{ufo_order:3}"]
        },"scalars")?;
        // let spinneys = spectacles.spinneys(&spectacles.full_filter());
        let f = Wood::new(
            CutStructure::empty(&sunrise),
            &sunrise,
            &scalar_model(),
            &GenerationSettings::default(),
        )?;
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph().n_nodes(),
        @"3",
        );
        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph().n_nodes(),
        @"4",
         );

        Ok(())
    }

    #[test]
    fn dotted() -> Result<()> {
        test_initialise().unwrap();

        let sunrise: Graph = finalized_runtime_dot!( digraph sunrise{
            graph [projector=1]
            node [num=1]
            edge [particle=scalar_1 num=1]
            e        [style=invis]
            e -> A:0 [id=4 sink="{ufo_order:0}"]
            B:1 -> e [id=5 source="{ufo_order:0}"]
            C:2 -> e [id=6 source="{ufo_order:0}"]

            A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:1}"]
            B -> C [id=1 source="{ufo_order:2}" sink="{ufo_order:1}"]
            C -> A [id=2 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"]
            B -> C [id=3 lmb_id=1 source="{ufo_order:3}" sink="{ufo_order:3}"]

        },"scalars")?;
        // let spinneys = spectacles.spinneys(&spectacles.full_filter());
        let f = Wood::new(
            CutStructure::empty(&sunrise),
            &sunrise,
            &scalar_model(),
            &GenerationSettings::default(),
        )?;
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph().n_nodes(),
        @"3",
        );
        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph().n_nodes(),
        @"4",
         );

        Ok(())
    }

    #[test]
    fn spectacles() -> Result<()> {
        use crate::graph::feynman_graph::FeynmanGraph;

        test_initialise().unwrap();

        let mut spectacles: Graph = finalized_runtime_dot!(
            digraph G{
                graph [projector=1]
                node [num=1]
                edge [particle="scalar_1" num=1];
                v1 -> v2 [source="{ufo_order:0}" sink="{ufo_order:0}"];
                v1 -> v2 [lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];

                v3 -> v4 [source="{ufo_order:0}" sink="{ufo_order:0}"];
                v3 -> v4 [lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:1}"];

                v2 -> v3 [source="{ufo_order:2}" sink="{ufo_order:2}"];
                v1 -> v4 [lmb_id=2 source="{ufo_order:2}" sink="{ufo_order:2}"];
            },"scalars"
        )?;

        // let spinneys = spectacles.spinneys(&spectacles.full_filter());
        let settings = GenerationSettings::default();
        let cut_structure = CutStructure::empty(&spectacles);
        let cutset = cut_structure
            .cuts
            .first()
            .expect("empty cut structure has one cut")
            .clone();
        let f = Wood::new(cut_structure, &spectacles, &scalar_model(), &settings)?;
        println!("{}", f);
        let mut f = f.unfold();
        println!("{}", f);

        let (union, edge) = f
            .graph()
            .iter_nodes()
            .find_map(|(child, _, _)| {
                let (parent, edge) = f.graph().unique_parent(child)?;
                (!f.graph().is_disjoint_union(child) && f.graph().is_disjoint_union(parent))
                    .then_some((parent, edge))
            })
            .expect("spectacles has a connected child above its disconnected union");
        let operations = f
            .graph()
            .iter_nodes()
            .map(|(_, _, operation)| operation.clone())
            .collect::<Vec<_>>();
        // Zero integrated terms isolate component-path replay in the local accumulator.
        for operation in operations {
            f.compute_store.entry(operation).or_default().integrated = Some(IntegratedCts::root());
        }
        let orientation_pattern = crate::settings::global::OrientationPattern::default();
        // The replay fixture keeps its synthetic unit branch, while the
        // existing exact source supplies a valid energy map for that key.
        let options = spectacles.denominator_only_cff_3d_expression_options();
        let canonization = spectacles.get_esurface_canonization(&spectacles.loop_momentum_basis);
        let contract_edges = spectacles.paired_edges(
            &spectacles
                .tree_edges
                .subtract(&spectacles.initial_state_cut),
        );
        let production = spectacles.generate_3d_expression_for_integrand(
            &contract_edges,
            &canonization,
            &options,
            None,
        )?;
        assert!(!production.expression.orientations.is_empty());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact_expression(
                &production,
                &options,
                &orientation_pattern,
                false,
            ),
        );
        let replay_states = f
            .forests()
            .union_replay_states(union)?
            .into_iter()
            .filter(|state| !state.local_edges.is_empty())
            .collect::<Vec<_>>();
        let union_active = replay_states
            .iter()
            .map(|state| {
                state
                    .local_edges
                    .iter()
                    .enumerate()
                    .map(|(offset, edge)| {
                        let step_order = f.graph()[state.integrated].key.op_count() + offset;
                        let (current, given) =
                            f.forests().wood().current_given_pair(*edge, step_order);
                        current.reduced_subgraph(&given)
                    })
                    .reduce(|active, reduced| active.union(&reduced))
                    .expect("a proper union replay state has a local suffix")
            })
            .collect::<Vec<_>>();
        let computed = f.local_3d_for_node(
            union,
            &mut spectacles,
            &cutset,
            localizer,
            &settings.uv,
            &mut Local4dProjectionContext::default(),
        )?;
        let sectors = computed.local_3d.direct()?.sectors()?;
        assert_eq!(
            sectors.len(),
            replay_states.len(),
            "every proper integrated prefix must produce one disconnected replay sector"
        );
        assert_eq!(
            sectors
                .iter()
                .map(|sector| sector.active_subgraph.clone())
                .collect::<Vec<_>>(),
            union_active,
            "the completed union must retain the active component support of every replay path"
        );
        assert!(
            computed.final_integrands.iter().next().is_some(),
            "the successfully replayed disconnected union must reach final-integrand construction"
        );
        let union_local = Direct3dCts::from_sectors(
            union_active
                .iter()
                .cloned()
                .map(|active_subgraph| {
                    let coordinate_frames = spectacles
                        .as_ref()
                        .connected_components(&active_subgraph)
                        .into_iter()
                        .filter(|component| spectacles.n_loops(component) != 0)
                        .map(|component| DirectCoordinateFrame {
                            lmb: spectacles.lmb_of(&component),
                            active_subgraph: component,
                        })
                        .collect();
                    Ok(DirectSector {
                        active_subgraph,
                        coordinate_frames,
                        active: DirectResidueBranches::production(
                            OrientationID(0),
                            Integrands::root(),
                        )?,
                        frozen_integrands: Integrands::root(),
                        projection_paths: Vec::new(),
                    })
                })
                .collect::<Result<_>>()?,
        )?;

        let step_order = f.graph()[union].key.op_count();
        let (current, given) = f.forests().wood().current_given_pair(edge, step_order);
        let reduced = current.reduced_subgraph(&given);
        let numerator = spectacles
            .numerator(&reduced, given.subgraph())?
            .get_single_atom()
            .expect("the scalar reduced numerator is available");
        assert!(!numerator.contains_symbol(GS.emr_mom));
        // A Taylor projection leaves these momentum-independent source values
        // unchanged and supplies the subtraction sign. Compare the complete
        // result so retained and quotient coordinates can use any compatible
        // basis without fixing the number or storage of coordinate frames.
        let expected = union_local
            .branches()?
            .materialize(false)?
            .resolved()?
            .map(|atom| -(atom * &numerator));
        let given_internal =
            InternalSubGraph::try_new(given.subgraph().clone(), spectacles.as_ref())
                .expect("the connected child has an internal given prefix");
        let quotient_rank = spectacles
            .shrunken_sub_lmb(
                current.subgraph(),
                &given_internal,
                spectacles.dummy_stripped_external_flows_of(current.subgraph()),
                None,
            )?
            .loop_edges
            .len();
        let expected_frames = union_local
            .sectors()?
            .iter()
            .map(|sector| {
                let retained_region = sector.active_subgraph.intersection(given.subgraph());
                let mut retained = sector
                    .coordinate_frames
                    .iter()
                    .filter(|frame| retained_region.includes(&frame.active_subgraph))
                    .flat_map(|frame| frame.lmb.loop_edges.iter().copied())
                    .collect::<Vec<_>>();
                retained.sort();
                retained.dedup();
                let expected_rank = retained.len() + quotient_rank;
                (retained, expected_rank)
            })
            .collect::<Vec<_>>();
        let external: SuBitGraph = spectacles.external_filter();
        let full_internal = spectacles.full_filter().subtract(&external);
        let canonical_route = if current.subgraph() == &full_internal
            && spectacles.n_loops(&full_internal) == spectacles.loop_momentum_basis.loop_edges.len()
        {
            spectacles.loop_momentum_basis.clone()
        } else {
            spectacles.try_compatible_sub_lmb(
                current.subgraph(),
                spectacles.dummy_less_full_crown(current.subgraph()),
                &spectacles.loop_momentum_basis,
            )?
        };
        let expected_active = union_active
            .iter()
            .map(|active| active.union(&reduced))
            .collect::<Vec<_>>();
        let child_local = Direct3dApproximation::new(localizer, &mut spectacles, &settings.uv)
            .run(
                &union_local,
                &IntegratedCts::root(),
                &current,
                &given,
                &current,
                &given,
            )?;
        let actual = child_local.branches()?.materialize(false)?.resolved()?;
        assert_eq!(
            actual.map(|atom| atom.collect_factors()),
            expected.map(|atom| atom.collect_factors()),
            "the connected child must retain every complete source value across its enclosing Taylor operation",
        );
        let child_sectors = child_local.sectors()?;
        let child_active = child_sectors
            .iter()
            .map(|sector| sector.active_subgraph.clone())
            .collect::<Vec<_>>();
        assert_eq!(child_active.len(), 3);
        assert_eq!(child_active, expected_active);
        for (sector, (retained, expected_rank)) in child_sectors.iter().zip(expected_frames.iter())
        {
            assert_eq!(
                sector.coordinate_frames.len(),
                1,
                "one connected enclosing Taylor operation must merge every contained child frame"
            );
            let child_carriers = sector
                .coordinate_frames
                .iter()
                .flat_map(|frame| frame.lmb.loop_edges.iter())
                .copied()
                .collect::<BTreeSet<_>>();
            assert!(
                retained.iter().all(|edge| child_carriers.contains(edge)),
                "the connected child must retain every carrier from the contained union frames"
            );
            assert!(
                child_carriers.iter().all(|edge| {
                    spectacles.loop_momentum_basis.edge_signatures[*edge]
                        .external
                        .iter()
                        .all(|sign| sign.is_zero())
                }),
                "the connected child route must be homogeneous in graph-external momenta"
            );
            assert_eq!(
                sector
                    .coordinate_frames
                    .iter()
                    .map(|frame| frame.lmb.loop_edges.len())
                    .sum::<usize>(),
                *expected_rank,
                "the connected child rank is the retained prefix rank plus its quotient rank"
            );
            let canonical_carriers = canonical_route
                .loop_edges
                .iter()
                .copied()
                .collect::<BTreeSet<_>>();
            let canonical_is_admissible = retained
                .iter()
                .all(|edge| canonical_carriers.contains(edge))
                && canonical_carriers.len() == *expected_rank
                && canonical_carriers.iter().all(|edge| {
                    spectacles.loop_momentum_basis.edge_signatures[*edge]
                        .external
                        .iter()
                        .all(|sign| sign.is_zero())
                });
            if canonical_is_admissible {
                assert_eq!(
                    child_carriers, canonical_carriers,
                    "an admissible graph-canonical route must win before equivalent fallbacks"
                );
            }
        }

        Ok(())
    }

    #[test]
    fn collective_two_component_region_requires_independently_divergent_factors() -> Result<()> {
        test_initialise().unwrap();

        // The two-line component has DOD two, while the scalar triangle has DOD
        // minus two. Their disconnected union has logarithmic simultaneous
        // power counting, but the triangle is not a UV region on its own.
        // Componentwise classification therefore retains only the bubble.
        let graph: Graph = finalized_runtime_dot!(
            digraph G { graph [projector=1]
                node [num = "1"];
                edge [particle = "scalar_1" num="1𝑖"];
                a0 -> a1 [id = 0, dod = "-1", num = "gammalooprs::Q(0,spenso::mink(4,0))" source="{ufo_order:0}" sink="{ufo_order:0}" lmb_id=0];
                a0 -> a1 [id = 1, dod = "-1", num = "gammalooprs::Q(1,spenso::mink(4,0))" source="{ufo_order:1}" sink="{ufo_order:1}"];
                a1 -> b0 [id = 2 source="{ufo_order:2}" sink="{ufo_order:0}"];
                b0 -> b1 [id = 3 source="{ufo_order:1}" sink="{ufo_order:0}" lmb_id=1];
                b1 -> b2 [id = 4 source="{ufo_order:1}" sink="{ufo_order:0}"];
                b2 -> b0 [id = 5 sink="{ufo_order:2}" source="{ufo_order:1}"];
            },
            "scalars"
        )?;

        let settings = GenerationSettings {
            uv: UVgenerationSettings {
                generate_integrated: false,
                ..Default::default()
            },
            ..Default::default()
        };
        let wood = Wood::new(CutStructure::empty(&graph), &graph, &graph.model, &settings)?;
        let collective = graph
            .spinneys(&graph.full_filter())
            .unwrap()
            .into_iter()
            .find(|spinney| graph.as_ref().connected_components(spinney).len() == 2)
            .expect("the raw cycle inventory contains the collective region");
        let collective_filter = collective.filter;
        let components = graph
            .as_ref()
            .connected_components(&collective_filter)
            .into_iter()
            .collect::<BTreeSet<_>>();
        assert_eq!(components.len(), 2);
        assert_eq!(graph.compute_dod(&collective_filter), 0);
        assert_eq!(
            components
                .iter()
                .map(|component| graph.compute_dod(component))
                .sorted()
                .collect::<Vec<_>>(),
            vec![-2, 2]
        );

        let root = OwnedForestNode {
            spinney: Spinney::empty(&graph),
            topo_order: 0,
        };
        let divergent_filter = components
            .iter()
            .find(|component| graph.compute_dod(*component) >= 0)
            .expect("the bubble is the independently divergent component")
            .clone();
        let divergent = OwnedForestNode {
            spinney: Spinney::new(
                InternalSubGraph::cleaned_filter_optimist(divergent_filter, graph.as_ref()),
                &graph,
                &graph.loop_momentum_basis,
            )?
            .expect("the divergent component has a compatible source chart"),
            topo_order: 0,
        };
        let collective_source = OwnedForestNode {
            spinney: Spinney::new(
                InternalSubGraph::cleaned_filter_optimist(
                    collective_filter.clone(),
                    graph.as_ref(),
                ),
                &graph,
                &graph.loop_momentum_basis,
            )?
            .expect("the collective region has a compatible source chart"),
            topo_order: 1,
        };
        let ctx = UVCtx::new(&graph, &graph.model, &settings.uv);
        let input = Full4dCts::from_factorized_local(&Local4dCts::root());
        let prefix = local_4d::uv_limit(&input, &ctx, &divergent, &root, &divergent, &root)?;
        // The two physical forests containing U are {U} and {A,U}:
        // (-T_U) I + (-T_U)(-T_A) I. Their full Taylor scope is U;
        // its convergent complement does not supply a separate counterterm.
        let atomic = local_4d::uv_limit(
            &input,
            &ctx,
            &collective_source,
            &root,
            &collective_source,
            &root,
        )?;
        let nested = local_4d::uv_limit(
            &Full4dCts::from_factorized_local(&prefix),
            &ctx,
            &collective_source,
            &divergent,
            &collective_source,
            &divergent,
        )?;
        let expected = atomic.atom().expression() + nested.atom().expression();
        let mut forests = wood.unfold();
        forests.integrate(&graph, &graph.model, crate::utils::vakint()?, &settings.uv)?;
        let actual = forests
            .graph()
            .iter_nodes()
            .filter(|(_, _, operation)| operation.covers().as_ref() == Some(&collective_filter))
            .try_fold(Atom::Zero, |sum, (_, _, operation)| -> Result<_> {
                Ok(sum
                    + forests
                        .compute_store
                        .require(operation)?
                        .local_4d(operation)?
                        .atom()
                        .expression())
            })?;
        // These vacuum components have no nonzero crown shifts. Put both
        // complete coefficients into the same formal K coordinates, retaining
        // factorization and every denominator multiplicity before comparison.
        let replacements = graph.integrand_replacement(
            &graph.full_filter(),
            &graph.loop_momentum_basis,
            &[Atom::var(W_.x___)],
        );
        let [actual, expected] = [actual, expected].map(|expression| -> Result<Atom> {
            let expression = GS
                .erase_uv_momentum_provenance(&expression)
                .replace(GS.den(W_.a_, W_.mom_, W_.mass_, W_.prop_))
                .with(W_.prop_)
                .replace_multiple(&replacements);
            let cooking = CookSettings::indices()
                .with_mode(CookMode::ReversibleEncoding)
                .with_representation_payloads(true, true);
            let expression = SymbolicTensor::infer_indexed(
                cooking.try_cook(expression.as_view())?,
            )?
            .contract(ContractSettings {
                rank_one: false,
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })?;
            Ok(cooking
                .uncook(expression.expression().as_view())
                .collect_factors())
        });
        let (actual, expected) = (actual?, expected?);
        assert!(
            !atomic.atom().is_zero(),
            "the atomic Taylor term must be nonzero"
        );
        assert!((actual - expected).collect_factors().is_zero());

        assert!(
            forests
                .forests()
                .wood()
                .graph()
                .iter_nodes()
                .all(|(_, _, spinney)| { spinney.filter() != &collective_filter }),
            "a convergent connected factor cannot enter a disconnected UV spinney"
        );
        let convergent = components
            .iter()
            .find(|component| graph.compute_dod(*component) < 0)
            .expect("the collective region has one convergent component");
        let divergent = components
            .iter()
            .find(|component| graph.compute_dod(*component) >= 0)
            .expect("the collective region has one independently divergent component");
        assert!(
            forests
                .forests()
                .wood()
                .graph()
                .iter_nodes()
                .all(|(_, _, spinney)| spinney.filter() != convergent)
        );
        assert!(
            forests
                .forests()
                .wood()
                .graph()
                .iter_nodes()
                .any(|(_, _, spinney)| spinney.filter() == divergent)
        );

        assert!(
            forests.graph().iter_nodes().all(|(_, _, operation)| {
                operation.covers().as_ref() != Some(&collective_filter)
            }),
            "neither atomic nor nested collective counterterms may bypass factor selection"
        );
        let mut cut = CutSet::empty(graph.n_hedges());
        cut.union = convergent.clone();
        assert!(
            forests.graph().iter_nodes().any(|(_, _, operation)| {
                operation.covers().as_ref() == Some(divergent)
                    && operation.is_compatible_with(&cut.union)
            }),
            "a disjoint independently divergent prefix remains cut-compatible"
        );
        cut.union = divergent.clone();
        assert!(
            forests
                .graph()
                .iter_nodes()
                .filter(|(_, _, operation)| { operation.covers().as_ref() == Some(divergent) })
                .all(|(_, _, operation)| !operation.is_compatible_with(&cut.union)),
            "cutting the retained divergent factor must exclude its counterterms"
        );

        Ok(())
    }

    #[test]
    fn collective_three_component_region_keeps_only_its_divergent_join() -> Result<()> {
        test_initialise().unwrap();

        // The two bubbles have DOD two and the scalar triangle has DOD minus two.
        // Every component nevertheless has a divergent complement because the
        // other two components have non-negative combined DOD. The triangle itself
        // supplies no independent counterterm: this aggregate edge-label cover
        // cannot replace the collective Taylor operation by a factorized product.
        // Only the two independently divergent bubbles form a factorized join;
        // their complete Taylor cancellation permits omitting the collective
        // transition rather than classifying the triangle as a UV factor.
        let graph: Graph = finalized_runtime_dot!(
            digraph G { graph [projector=1]
                node [num = "1"];
                edge [particle = "scalar_1" num="1𝑖"];
                a0 -> a1 [id = 0, dod = "-1", num = "gammalooprs::Q(0,spenso::mink(4,0))" source="{ufo_order:0}" sink="{ufo_order:0}" lmb_id=0];
                a0 -> a1 [id = 1, dod = "-1", num = "gammalooprs::Q(1,spenso::mink(4,0))" source="{ufo_order:1}" sink="{ufo_order:1}"];
                a1 -> b0 [id = 2 source="{ufo_order:2}" sink="{ufo_order:0}"];
                b0 -> b1 [id = 3, dod = "-1", num = "gammalooprs::Q(3,spenso::mink(4,1))" source="{ufo_order:1}" sink="{ufo_order:0}" lmb_id=1];
                b0 -> b1 [id = 4, dod = "-1", num = "gammalooprs::Q(4,spenso::mink(4,1))" source="{ufo_order:2}" sink="{ufo_order:1}"];
                b1 -> c0 [id = 5 source="{ufo_order:2}" sink="{ufo_order:0}"];
                c0 -> c1 [id = 6 source="{ufo_order:1}" sink="{ufo_order:0}" lmb_id=2];
                c1 -> c2 [id = 7 source="{ufo_order:1}" sink="{ufo_order:0}"];
                c2 -> c0 [id = 8 sink="{ufo_order:2}" source="{ufo_order:1}"];
            },
            "scalars"
        )?;

        let settings = GenerationSettings {
            uv: UVgenerationSettings {
                generate_integrated: false,
                ..Default::default()
            },
            ..Default::default()
        };
        let wood = Wood::new(CutStructure::empty(&graph), &graph, &graph.model, &settings)?;
        let collective = graph
            .spinneys(&graph.full_filter())
            .unwrap()
            .into_iter()
            .find(|spinney| graph.as_ref().connected_components(spinney).len() == 3)
            .expect("the raw cycle inventory contains the collective region");
        let collective_filter = collective.filter;
        let components = graph
            .as_ref()
            .connected_components(&collective_filter)
            .into_iter()
            .collect::<BTreeSet<_>>();
        assert_eq!(components.len(), 3);
        assert_eq!(graph.compute_dod(&collective_filter), 2);
        assert_eq!(
            components
                .iter()
                .map(|component| graph.compute_dod(component))
                .sorted()
                .collect::<Vec<_>>(),
            vec![-2, 2, 2]
        );
        assert!(
            wood.graph().iter_nodes().all(|(_, _, spinney)| {
                graph
                    .as_ref()
                    .connected_components(spinney.filter())
                    .iter()
                    .all(|component| graph.compute_dod(component) >= 0)
            }),
            "a convergent connected factor cannot enter a disconnected UV spinney"
        );
        assert!(
            components.iter().all(|component| {
                graph.compute_dod(&collective_filter.subtract(component)) >= 0
            })
        );
        let retained_components = components
            .iter()
            .filter(|component| {
                wood.graph()
                    .iter_nodes()
                    .any(|(_, _, spinney)| spinney.filter() == *component)
            })
            .cloned()
            .collect::<BTreeSet<_>>();
        assert_eq!(
            retained_components.len(),
            2,
            "complete complement coverage must not stand in for independent component counterterms"
        );
        assert!(
            retained_components
                .iter()
                .all(|component| graph.compute_dod(component) >= 0)
        );
        let retained_union = retained_components
            .iter()
            .cloned()
            .reduce(|left, right| left.union(&right))
            .unwrap();
        let join = wood
            .graph()
            .iter_nodes()
            .find_map(|(node, _, spinney)| (spinney.filter() == &retained_union).then_some(node))
            .expect("the two independently divergent bubbles retain their union");
        assert_eq!(wood.wood().join_factors(join), Some(retained_components));

        let root = OwnedForestNode {
            spinney: Spinney::empty(&graph),
            topo_order: 0,
        };
        let root_local = Local4dCts::root();
        let input = Full4dCts::from_factorized_local(&root_local);
        let ctx = UVCtx::new(&graph, &graph.model, &settings.uv);
        let convergent = components
            .iter()
            .find(|component| graph.compute_dod(*component) < 0)
            .expect("the triangle is the convergent component");
        // Use one common physical chart only because these vacuum components
        // have no nonzero crown shifts. This is not a general replacement for
        // the source reconstruction certificates used by projected CFF.
        let replacements = graph.integrand_replacement(
            &graph.full_filter(),
            &graph.loop_momentum_basis,
            &[Atom::var(W_.x___)],
        );
        let mut independent = Vec::<(OwnedForestNode, Local4dCts)>::new();
        for component in components
            .iter()
            .filter(|component| graph.compute_dod(*component) >= 0)
        {
            let source = OwnedForestNode {
                spinney: Spinney::new(
                    InternalSubGraph::cleaned_filter_optimist(component.clone(), graph.as_ref()),
                    &graph,
                    &graph.loop_momentum_basis,
                )?
                .expect("each divergent bubble has a compatible source chart"),
                topo_order: 0,
            };
            let local = local_4d::uv_limit(&input, &ctx, &source, &root, &source, &root)?;
            let mixed = OwnedForestNode {
                spinney: Spinney::new(
                    InternalSubGraph::cleaned_filter_optimist(
                        component.union(convergent),
                        graph.as_ref(),
                    ),
                    &graph,
                    &graph.loop_momentum_basis,
                )?
                .expect("the bubble-triangle collective region is logarithmic"),
                topo_order: 1,
            };
            let atomic = local_4d::uv_limit(&input, &ctx, &mixed, &root, &mixed, &root)?;
            let nested = local_4d::uv_limit(
                &Full4dCts::from_factorized_local(&local),
                &ctx,
                &mixed,
                &source,
                &mixed,
                &source,
            )?;
            assert!(
                !atomic.atom().is_zero(),
                "the mixed atomic Taylor term must be nonzero"
            );
            // T_AC(1-T_A)I and T_BC(1-T_B)I vanish: the degree-two
            // bubble subtraction leaves no divergent scaling with the triangle.
            // Their complete contributions therefore cancel even if a forest
            // representation includes or removes zero-quotient outer operations.
            let mixed_remainder = GS
                .erase_uv_momentum_provenance(
                    &(atomic.atom().expression() + nested.atom().expression()),
                )
                .replace(GS.den(W_.a_, W_.mom_, W_.mass_, W_.prop_))
                .with(W_.prop_)
                .replace_multiple(&replacements);
            let cooking = CookSettings::indices()
                .with_mode(CookMode::ReversibleEncoding)
                .with_representation_payloads(true, true);
            let mixed_remainder =
                SymbolicTensor::infer_indexed(cooking.try_cook(mixed_remainder.as_view())?)?
                    .contract(ContractSettings {
                        rank_one: false,
                        collect_chains: false,
                        collect_traces: false,
                        ..Default::default()
                    })?;
            let mixed_remainder = cooking
                .uncook(mixed_remainder.expression().as_view())
                .collect_factors();
            assert!(
                mixed_remainder.is_zero(),
                "mixed collective Taylor terms must cancel"
            );
            independent.push((source, local));
        }
        let joint_scope = independent.iter().fold(
            graph.empty_subgraph::<SuBitGraph>(),
            |scope, (source, _)| scope.union(source.subgraph()),
        );
        let joint = OwnedForestNode {
            spinney: Spinney::new(
                InternalSubGraph::cleaned_filter_optimist(joint_scope, graph.as_ref()),
                &graph,
                &graph.loop_momentum_basis,
            )?
            .expect("the disjoint divergent bubbles have a compatible joint source chart"),
            topo_order: 1,
        };
        let product = Local4dCts::from_full_product(
            independent
                .iter()
                .map(|(_, local)| Full4dCts::from_factorized_local(local)),
        )?;
        let collective_source = OwnedForestNode {
            spinney: Spinney::new(
                InternalSubGraph::cleaned_filter_optimist(
                    collective_filter.clone(),
                    graph.as_ref(),
                ),
                &graph,
                &graph.loop_momentum_basis,
            )?
            .expect("the collective region has a compatible source chart"),
            topo_order: 2,
        };
        // The surviving complete value is -T_U(1-T_A)(1-T_B)I.
        // AB supplies a product of independent counterterms, while the
        // convergent triangle supplies none; no incoming-edge or join catalogue
        // determines these four physical prefix values.
        let mut expected = Atom::Zero;
        for (given, local) in std::iter::once((&root, &root_local))
            .chain(independent.iter().map(|(source, local)| (source, local)))
            .chain(std::iter::once((&joint, &product)))
        {
            let term = local_4d::uv_limit(
                &Full4dCts::from_factorized_local(local),
                &ctx,
                &collective_source,
                given,
                &collective_source,
                given,
            )?;
            debug_tags!(#uv, #forest, #trace;
                stage = "collective_three_component_oracle_term",
                given = %given.subgraph().string_label(),
                file.prefix = %local.atom().expression().to_plain_string(),
                file.term = %term.atom().expression().to_plain_string(),
                "Independent collective Taylor contribution"
            );
            if given.subgraph().is_empty() {
                assert!(
                    !term.atom().is_zero(),
                    "the full atomic Taylor term must be nonzero"
                );
            }
            expected += term.atom().expression();
        }
        let mut forests = wood.unfold();
        forests.integrate(&graph, &graph.model, crate::utils::vakint()?, &settings.uv)?;
        let actual = forests
            .graph()
            .iter_nodes()
            .filter(|(_, _, operation)| operation.covers().as_ref() == Some(&collective_filter))
            .try_fold(Atom::Zero, |sum, (_, _, operation)| -> Result<_> {
                let term = forests
                    .compute_store
                    .require(operation)?
                    .local_4d(operation)?;
                debug_tags!(#uv, #forest, #trace;
                    stage = "collective_three_component_forest_term",
                    operation = %operation,
                    file.term = %term.atom().expression().to_plain_string(),
                    "Scheduled collective Taylor contribution"
                );
                Ok(sum + term.atom().expression())
            })?;
        let [actual, expected] = [actual, expected].map(|expression| -> Result<Atom> {
            let expression = GS
                .erase_uv_momentum_provenance(&expression)
                .replace(GS.den(W_.a_, W_.mom_, W_.mass_, W_.prop_))
                .with(W_.prop_)
                .replace_multiple(&replacements);
            let cooking = CookSettings::indices()
                .with_mode(CookMode::ReversibleEncoding)
                .with_representation_payloads(true, true);
            let expression =
                SymbolicTensor::infer_indexed(cooking.try_cook(expression.as_view())?)?
                    .contract(ContractSettings {
                        rank_one: false,
                        collect_chains: false,
                        collect_traces: false,
                        ..Default::default()
                    })?
                    // Typed composition can spell a completed scalar contraction
                    // as g(p,q). Use the same dot notation before selecting the
                    // tensor basis; no indexed contraction or finite-D value is added.
                    .to_dots()?;
            Ok(cooking.uncook(expression.expression().as_view()))
        });
        let (actual, expected) = (actual?, expected?);
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let actual = SymbolicTensor::infer_indexed(cooking.try_cook(actual.as_view())?)?;
        let expected = SymbolicTensor::infer_indexed(cooking.try_cook(expected.as_view())?)?;
        // This independent diagnostic selects tensor factors whose individual
        // interfaces may be unresolved; it does not store a runtime numerator.
        let difference = SymbolicTensor::infer(actual.expression() - expected.expression())?;
        // The fixture's numerator consists of its selected Q tensor factors.
        // Compare their complete scalar denominator weights in the same formal
        // dimension; finite-component specialization would narrow the identity.
        for (tensor, coefficient) in
            difference.coefficient_list(TensorCollectFilter::<0>::Tensors)?
        {
            assert!(coefficient.is_scalar());
            assert!(
                tensor.expression().is_zero() || coefficient.expression().together().is_zero(),
                "the complete collective Taylor value must equal -T_U(1-T_A)(1-T_B)I: {}",
                coefficient.expression()
            );
        }

        assert!(
            forests.graph().iter_nodes().all(|(_, _, operation)| {
                operation.covers().is_none_or(|scope| {
                    graph
                        .as_ref()
                        .connected_components(&scope)
                        .iter()
                        .all(|component| graph.compute_dod(component) >= 0)
                })
            }),
            "no atomic or nested collective operation may include the convergent triangle"
        );
        assert!(
            forests.graph().iter_nodes().any(|(node, _, operation)| {
                operation.covers().as_ref() == Some(&retained_union)
                    && forests.graph().is_disjoint_union(node)
            }),
            "the complete independently divergent cover remains a factorized join"
        );

        Ok(())
    }

    #[test]
    fn spectacles_typed_4d_local_construction_matches_uv_limit() -> Result<()> {
        test_initialise().unwrap();

        let spectacles: Graph = finalized_runtime_dot!(digraph G {
          projector = "1";
          edge	 [
        num = "1𝑖", particle = "scalar_1"]
          node	 [
        num = "UFO::SCALAR_COUPLING"]
          v1;
          v2;
          v3;
          v4;

          v1:0	-> v2:1	 [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"  lmb_id="0"];
          v1:2	-> v2:3	 [id=1 source="{ufo_order:1}" sink="{ufo_order:1}"  lmb_id="1"];
          v1:4	-> v4:5	 [id=2 source="{ufo_order:2}" sink="{ufo_order:0}" ];
          v2:6	-> v3:7	 [id=3 source="{ufo_order:2}" sink="{ufo_order:0}" ];
          v3:8	-> v4:9	 [id=4 source="{ufo_order:1}" sink="{ufo_order:1}"  lmb_id="2"];
          v3:10	-> v4:11	 [id=5 source="{ufo_order:2}" sink="{ufo_order:2}" ];
        },"scalars"
                )?;

        // let spinneys = spectacles.spinneys(&spectacles.full_filter());
        let settings = GenerationSettings {
            uv: UVgenerationSettings {
                generate_integrated: false,
                local_uv_cts_from_expanded_4d_integrands: true,
                ..Default::default()
            },
            ..Default::default()
        };
        let cut_structure = CutStructure::empty(&spectacles);
        let f = Wood::new(cut_structure, &spectacles, &spectacles.model, &settings)?;
        println!("{}", f);
        let mut f = f.unfold();
        println!("{}", f);

        let (union, child, edge) = f
            .graph()
            .iter_nodes()
            .find_map(|(child, _, _)| {
                let (parent, edge) = f.graph().unique_parent(child)?;
                (!f.graph().is_disjoint_union(child) && f.graph().is_disjoint_union(parent))
                    .then_some((parent, child, edge))
            })
            .expect("spectacles has a connected child above its disconnected union");
        f.integrate(
            &spectacles,
            &spectacles.model,
            crate::utils::vakint()?,
            &settings.uv,
        )?;

        let step_order = f.graph()[union].key.op_count();
        let (current, given) = f.forests().wood().current_given_pair(edge, step_order);
        let expected = local_4d::uv_limit(
            &f.recursion_input_4d(union)?,
            &UVCtx::new(&spectacles, &spectacles.model, &settings.uv),
            &current,
            &given,
            &current,
            &given,
        )?;
        let operation = &f.graph()[child];
        let local = f.compute_store.require(operation)?.local_4d(operation)?;
        assert!(
            (local.atom().expression().collect_factors()
                - expected.atom().expression().collect_factors())
            .is_zero()
        );

        Ok(())
    }

    #[test]
    fn basketball() -> Result<()> {
        test_initialise().unwrap();

        let basketball: Graph = finalized_runtime_dot!(
            digraph G{
                graph [projector=1]
                node [num=1]
                edge [particle="scalar_1" num=1];
                v1 -> v2 [source="{ufo_order:0}" sink="{ufo_order:0}"];
                v1 -> v2 [lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                v1 -> v2 [lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:2}"];
                v1 -> v2 [lmb_id=2 source="{ufo_order:3}" sink="{ufo_order:3}"];
            },"scalars"
        )?;

        let f = Wood::new(
            CutStructure::empty(&basketball),
            &basketball,
            &scalar_model(),
            &GenerationSettings::default(),
        )?;
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph().n_nodes(),
        @"12",
              );
        let f = f.unfold();
        println!("{}", f.dot_serialize());
        insta::assert_snapshot!(
        f.graph().n_nodes(),
        @"46");
        Ok(())
    }

    #[test]
    fn fourloop_b() -> Result<()> {
        test_initialise().unwrap();

        let fourloop_b: Graph = finalized_runtime_dot!(
            digraph G{
                graph [projector=1]
                node [num=1]
                edge [particle="scalar_1" num=1];
                v1 -> v2 [source="{ufo_order:0}" sink="{ufo_order:0}"];
                v1 -> v2 [lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];

                v3 -> v4 [source="{ufo_order:0}" sink="{ufo_order:0}"];
                v3 -> v4 [lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:1}"];

                v2 -> v3 [source="{ufo_order:2}" sink="{ufo_order:2}"];
                v1 -> v3 [lmb_id=2 source="{ufo_order:2}" sink="{ufo_order:3}"];
                v1 -> v4 [lmb_id=3 source="{ufo_order:3}" sink="{ufo_order:2}"];
            },"scalars"
        )?;

        // let spinneys = fourloop_b.spinneys(&fourloop_b.full_filter());
        let f = Wood::new(
            CutStructure::empty(&fourloop_b),
            &fourloop_b,
            &scalar_model(),
            &GenerationSettings::default(),
        )?;
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph().n_nodes(),
        @"14",
        );
        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph().n_nodes(),
        @"80",
         );

        Ok(())
    }

    #[test]
    fn four_loop_a() -> Result<()> {
        test_initialise().unwrap();

        let four_loop_a: Graph = finalized_runtime_dot!(
            digraph G{
                graph [projector=1]
                node [num=1]
                edge [particle="scalar_1" num=1];
                v1 -> v2 [source="{ufo_order:0}" sink="{ufo_order:0}"];
                v1 -> v2 [lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                v2 -> v3 [source="{ufo_order:2}" sink="{ufo_order:0}"];
                v3 -> v1 [lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:2}"];
                v1 -> v4 [source="{ufo_order:3}" sink="{ufo_order:0}"];
                v2 -> v4 [lmb_id=2 source="{ufo_order:3}" sink="{ufo_order:1}"];
                v3 -> v4 [lmb_id=3 source="{ufo_order:2}" sink="{ufo_order:2}"];
            },"scalars"
        )?;

        let f = Wood::new(
            CutStructure::empty(&four_loop_a),
            &four_loop_a,
            &scalar_model(),
            &GenerationSettings::default(),
        )?;
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph().n_nodes(),
        @"12",
        );
        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph().n_nodes(),
        @"60",
         );

        Ok(())
    }

    #[test]
    fn triple_double_tadpole() -> Result<()> {
        test_initialise().unwrap();
        let dumbell: Graph = finalized_runtime_dot!(
            digraph G{
                graph [projector=1]
                node [num=1]
                edge [particle="scalar_1" num=1];
                v1 -> v2 [source="{ufo_order:0}" sink="{ufo_order:0}"];
                v2 -> v3 [source="{ufo_order:1}" sink="{ufo_order:0}"];
                v3 -> v3 [lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:2}"];
                v3 -> v3 [lmb_id=1 source="{ufo_order:3}" sink="{ufo_order:4}"];
                v2 -> v2 [lmb_id=2 source="{ufo_order:2}" sink="{ufo_order:3}"];
                v2 -> v2 [lmb_id=3 source="{ufo_order:4}" sink="{ufo_order:5}"];
                v1 -> v1 [lmb_id=4 source="{ufo_order:1}" sink="{ufo_order:2}"];
                v1 -> v1 [lmb_id=5 source="{ufo_order:3}" sink="{ufo_order:4}"];
            },"scalars"
        )?;

        let f = Wood::new(
            CutStructure::empty(&dumbell),
            &dumbell,
            &scalar_model(),
            &GenerationSettings::default(),
        )?;

        insta::assert_snapshot!(
            f.graph().n_nodes(),
            @"64");

        let f = f.unfold_uncached();
        insta::assert_snapshot!(
            f.graph().n_nodes(),
            @"307");

        Ok(())
    }

    mod failing {
        use super::*;

        #[test]
        fn lobsided_double_dumbell() -> Result<()> {
            test_initialise().unwrap();
            let dumbell: Graph = finalized_runtime_dot!(
                digraph G{
                    graph [projector=1]
                    node [num=1]
                    edge [particle="scalar_1" num=1];
                    v1 -> v2 [source="{ufo_order:0}" sink="{ufo_order:0}"];
                    v2 -> v2 [lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:2}"];
                    v1 -> v1 [lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:2}"];
                    v1 -> v1 [lmb_id=2 source="{ufo_order:3}" sink="{ufo_order:4}"];
                },"scalars"
            )?;

            let spinneys: Vec<_> = dumbell
                .spinneys(&dumbell.full_filter())
                .unwrap()
                .into_iter()
                .filter_map(|a| Spinney::new(a, &dumbell, &dumbell.loop_momentum_basis).transpose())
                .collect::<Result<_, _>>()?;
            let f = Wood::new(
                CutStructure::empty(&dumbell),
                &dumbell,
                &scalar_model(),
                &GenerationSettings::default(),
            )?;

            println!("{}", f);

            insta::assert_snapshot!(
                f.graph().n_nodes(),
                @"8",
                // format!("Wood does not have correct number of spinneys: \n{}",f)
            );

            for (_, _, d) in f.graph().iter_nodes() {
                println!(
                    "//Node {}: \n{}",
                    d.subgraph().string_label(),
                    dumbell.dot(d.subgraph())
                );
            }
            let _ff = OldWood::from_spinneys(spinneys, &dumbell); //.unfold(&g, &g.loop_momentum_basis);

            // println!("{}", ff.dot(&dumbell));

            let f = f.unfold();
            f.debug_walk();
            println!("{}", f);
            insta::assert_snapshot!(
                f.normalized_node_labels_with_cover("3L").join("\n"),
                @r###"
{36,F}: Approx(S_36(_))*Approx(S_F(_))
{3} · {36,F}: Approx((-1*S_3+S_F(_))*Approx(S_3(_)))*Approx(S_36(_))
{C} · {36,F}: Approx((-1*S_C+S_F(_))*Approx(S_C(_)))*Approx(S_36(_))
"###
            );
            insta::assert_snapshot!(
                f.graph().n_nodes(),
                @"12");

            let f = Wood::new(
                CutStructure::empty(&dumbell),
                &dumbell,
                &scalar_model(),
                &GenerationSettings::default(),
            )?
            .unfold_uncached();
            assert!(f.compute_store.entries.is_empty());
            insta::assert_snapshot!(
                f.normalized_node_labels_with_cover("3L").join("\n"),
                @r###"
{36,F}: Approx(S_36(_))*Approx(S_F(_))
{3} · {36,F}: Approx((-1*S_3+S_F(_))*Approx(S_3(_)))*Approx(S_36(_))
{C} · {36,F}: Approx((-1*S_C+S_F(_))*Approx(S_C(_)))*Approx(S_36(_))
"###
            );

            Ok(())
        }

        #[test]
        fn double_double_dumbell() -> Result<()> {
            test_initialise().unwrap();
            let dumbell: Graph = finalized_runtime_dot!(
                digraph G{
                    graph [projector=1]
                    node [num=1]
                    edge [particle="scalar_1" num=1];
                    v1 -> v2 [source="{ufo_order:0}" sink="{ufo_order:0}"];
                    v2 -> v2 [lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:2}"];
                    v2 -> v2 [lmb_id=1 source="{ufo_order:3}" sink="{ufo_order:4}"];
                    v1 -> v1 [lmb_id=2 source="{ufo_order:1}" sink="{ufo_order:2}"];
                    v1 -> v1 [lmb_id=3 source="{ufo_order:3}" sink="{ufo_order:4}"];
                },"scalars"
            )?;

            let _spinneys: Vec<_> = dumbell
                .spinneys(&dumbell.full_filter())
                .unwrap()
                .into_iter()
                .filter_map(|a| Spinney::new(a, &dumbell, &dumbell.loop_momentum_basis).transpose())
                .collect::<Result<_, _>>()?;
            let f = Wood::new(
                CutStructure::empty(&dumbell),
                &dumbell,
                &scalar_model(),
                &GenerationSettings::default(),
            )?;

            println!("{}", f);

            insta::assert_snapshot!(
                f.graph().n_nodes(),
                @"16",
                // format!("Wood does not have correct number of spinneys: \n{}",f)
            );

            let f = f.unfold();

            insta::assert_snapshot!(
                f.graph().n_nodes(),
                @"36");

            let f = Wood::new(
                CutStructure::empty(&dumbell),
                &dumbell,
                &scalar_model(),
                &GenerationSettings::default(),
            )?
            .unfold_uncached();

            let foata_labels = f
                .graph()
                .iter_nodes()
                .map(|(_, _, key)| key.foata_level_labels())
                .collect::<Vec<_>>();
            assert!(
                foata_labels
                    .iter()
                    .any(|label| label == "3,36;FU,F" || label == "36,3;FU,F"),
                "expected 3;F and 36;FU branch histories to combine, got:\n{}",
                foata_labels.join("\n")
            );

            Ok(())
        }
    }

    /// Known limitation of the reference enumeration, unchanged by the move:
    /// pruning grows unions only from retained ones, so a divergent union that
    /// is not the union of two cycle-space elements can be missed. The scalar
    /// Petersen graph with one incident energy factor (stored local degree one)
    /// at six vertices has a logarithmically divergent full graph, by both its
    /// pruning and its classifying degree, and only pruned proper unions. This
    /// records today's answer next to the complete one.
    #[test]
    fn known_limitation_petersen_full_graph_is_not_reached() -> Result<()> {
        test_initialise().unwrap();
        let graph: Graph = finalized_runtime_dot!(
            digraph G{
                graph [projector=1]
                node [num=1]
                edge [particle="scalar_1" num=1];
                v0 [dod=1 num="Q(0,spenso::cind(0))"]; v1 [dod=1 num="Q(1,spenso::cind(0))"];
                v2 [dod=1 num="Q(2,spenso::cind(0))"]; v3 [dod=1 num="Q(3,spenso::cind(0))"];
                v4 [dod=1 num="Q(4,spenso::cind(0))"]; v5 [dod=1 num="Q(5,spenso::cind(0))"];
                v0 -> v1 [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
                v1 -> v2 [id=1 source="{ufo_order:1}" sink="{ufo_order:0}"];
                v2 -> v3 [id=2 source="{ufo_order:1}" sink="{ufo_order:0}"];
                v3 -> v4 [id=3 source="{ufo_order:1}" sink="{ufo_order:0}"];
                v4 -> v0 [id=4 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                v0 -> v5 [id=5 source="{ufo_order:2}" sink="{ufo_order:0}"];
                v1 -> v6 [id=6 source="{ufo_order:2}" sink="{ufo_order:0}"];
                v2 -> v7 [id=7 source="{ufo_order:2}" sink="{ufo_order:0}"];
                v3 -> v8 [id=8 source="{ufo_order:2}" sink="{ufo_order:0}"];
                v4 -> v9 [id=9 source="{ufo_order:2}" sink="{ufo_order:0}"];
                v5 -> v7 [id=10 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:1}"];
                v7 -> v9 [id=11 lmb_id=2 source="{ufo_order:2}" sink="{ufo_order:1}"];
                v9 -> v6 [id=12 lmb_id=3 source="{ufo_order:2}" sink="{ufo_order:1}"];
                v6 -> v8 [id=13 lmb_id=4 source="{ufo_order:2}" sink="{ufo_order:1}"];
                v8 -> v5 [id=14 lmb_id=5 source="{ufo_order:2}" sink="{ufo_order:2}"];
            },"scalars"
        )?;
        let full = graph.full_filter();
        let unions = graph.all_cycle_unions(&full);
        assert_eq!(unions.len(), 594);
        let complete = unions
            .iter()
            .find(|union| union.filter == full)
            .expect("the full graph is a union of cycles");
        assert_eq!(graph.local_degree_of_divergence(&full)?, 0);
        assert_eq!(graph.degree_of_divergence(&full)?, 0);
        for union in unions
            .iter()
            .filter(|union| !union.is_empty() && union.filter != full)
        {
            assert!(graph.local_degree_of_divergence(&union.filter)? < 0);
        }

        // Complete answer: the full graph passes pruning. Today's answer: only
        // the empty subgraph is enumerated.
        let enumerated = graph.spinneys(&full)?;
        assert!(!enumerated.contains(complete));
        assert_eq!(enumerated.len(), 1);
        Ok(())
    }

    /// Production `Wood::new` removes the sewn initial state and each cut before
    /// classifying. A forward-scattering bubble cut through both its lines keeps
    /// only the tadpole on the initial-state line.
    #[test]
    fn production_wood_excludes_the_initial_state_and_the_cut() -> Result<()> {
        test_initialise().unwrap();
        let diagram = FeynmanDiagram::from_dot(
            FeynkitModel::phi_3_4(),
            r#"digraph sewn {
                graph [final_state="phi,phi"];
                edge [particle=phi];
                ext [style=invis];
                ext -> c [is_cut=9]; b -> ext [is_cut=9];
                c -> c; c -> a; a -> b; a -> b;
            }"#,
        )?
        .apply_feynman_rules()?;
        let graph = diagram.to_gamma_loop_graph(None, true)?;
        let edge_of = |source: &str, target: &str| {
            let vertex = |name: &str| {
                diagram
                    .vertices()
                    .find(|(_, vertex)| vertex.name == name)
                    .map(|(vertex, _)| vertex)
            };
            diagram
                .edges()
                .filter(|(_, endpoints, _)| {
                    endpoints.source == vertex(source) && endpoints.target == vertex(target)
                })
                .map(|(edge, _, _)| EdgeIndex(edge.0))
                .collect::<Vec<_>>()
        };
        let mut cut = CutSet::empty(graph.n_hedges());
        for edge in edge_of("a", "b") {
            cut.union.add(graph[&edge].1);
        }
        let wood = Wood::new(
            CutStructure { cuts: vec![cut] },
            &graph,
            &graph.model,
            &GenerationSettings::default(),
        )?;
        let inventory = wood
            .graph()
            .iter_nodes()
            .map(|(_, _, spinney)| {
                let mut edges = spinney.edges(&graph);
                edges.sort();
                edges
            })
            .collect::<BTreeSet<_>>();
        assert_eq!(inventory, BTreeSet::from([vec![], edge_of("c", "c")]));
        Ok(())
    }

    /// Two sewn initial-state lines close a bubble of their own. Production
    /// `Wood::new` must exclude the initial state as well as the cut: without the
    /// initial-state exclusion the enumeration retains the sewn bubble.
    #[test]
    fn production_wood_excludes_two_sewn_initial_state_lines() -> Result<()> {
        test_initialise().unwrap();
        let diagram = FeynmanDiagram::from_dot(
            FeynkitModel::phi4(),
            r#"digraph sewn_pair {
                graph [final_state="phi,phi"];
                edge [particle=phi];
                ext [style=invis];
                ext -> a [is_cut=1]; ext -> a [is_cut=2];
                b -> ext [is_cut=1]; b -> ext [is_cut=2];
                a -> b; a -> b;
            }"#,
        )?
        .apply_feynman_rules()?;
        let graph = diagram.to_gamma_loop_graph(None, true)?;
        let mut cut = CutSet::empty(graph.n_hedges());
        for (pair, _, _) in graph.underlying.iter_edges_of(&diagram.internal_subgraph()) {
            cut.union.add(pair);
        }

        let without_initial_state = graph.full_filter().subtract(&cut.union);
        let sewn_bubble = graph.spinneys(&without_initial_state)?;
        assert!(sewn_bubble.iter().any(|union| !union.is_empty()));

        let wood = Wood::new(
            CutStructure { cuts: vec![cut] },
            &graph,
            &graph.model,
            &GenerationSettings::default(),
        )?;
        assert_eq!(wood.graph().n_nodes(), 1);
        assert!(
            wood.graph()
                .iter_nodes()
                .all(|(_, _, spinney)| spinney.is_empty())
        );
        Ok(())
    }

    /// Structured forest DOT of base revision `055fd0c1`, before the combinatorics
    /// moved to feynkit-uv, for a graph with nested and factorised operations: a
    /// tadpole joined by a bridge to a bubble carrying a second tadpole.
    const BASE_NESTED_FACTORISED_FOREST_DOT: &str = r##"
digraph {
  node	 [shape=circle,height=0.1,label=""];
  overlap = "scale";
  layout = "neato";
  start=2;
  0	 [label="#expr" foata="" cover="0"];
  1	 [label="K[#expr S_CO⊛0]" foata="CO" cover="CO"];
  2	 [label="K[#expr S_3s⊛0]" foata="3s" cover="3s"];
  3	 [label="K[#expr S_GG⊛0]" foata="GG" cover="GG"];
  4	 [label="K[#expr S_3⊛0]" foata="3" cover="3"];
  5	 [label="K[S_GG⊛CO K[#expr S_CO⊛0]]" foata="CO;GG" cover="GG"];
  6	 [label="K[#expr S_CO⊛0] K[#expr S_CR⊛CO]" foata="CO,3" cover="CR"];
  7	 [label="K[S_GG⊛3s K[#expr S_3s⊛0]]" foata="3s;GG" cover="GG"];
  8	 [label="K[#expr S_3s⊛0] K[#expr S_3v⊛3s]" foata="3s,3" cover="3v"];
  9	 [label="K[#expr S_GG⊛0] K[#expr S_GJ⊛GG]" foata="GG,3" cover="GJ"];
  10	 [label="K[#expr S_3⊛0] K[S_GG⊛CO K[#expr S_CO⊛0]]" foata="CO;GG,3" cover="GJ"];
  11	 [label="K[#expr S_3⊛0] K[S_GG⊛3s K[#expr S_3s⊛0]]" foata="3s;GG,3" cover="GJ"];
  0:0:s	-> 1:1:s	 [id=0  color="red:blue;0.5"];
  0:2:s	-> 2:3:s	 [id=1  color="red:blue;0.5"];
  0:4:s	-> 3:5:s	 [id=2  color="red:blue;0.5"];
  0:6:s	-> 4:7:s	 [id=3  color="red:blue;0.5"];
  1:8:s	-> 5:9:s	 [id=4  color="red:blue;0.5"];
  1:10:s	-> 6:11:s	 [id=5  color="red:blue;0.5"];
  2:12:s	-> 7:13:s	 [id=6  color="red:blue;0.5"];
  2:14:s	-> 8:15:s	 [id=7  color="red:blue;0.5"];
  3:16:s	-> 9:17:s	 [id=8  color="red:blue;0.5"];
  4:18:s	-> 9:19:s	 [id=9  color="red:blue;0.5"];
  4:20:s	-> 8:21:s	 [id=10  color="red:blue;0.5"];
  4:22:s	-> 6:23:s	 [id=11  color="red:blue;0.5"];
  4:24:s	-> 10:25:s	 [id=12  color="red:blue;0.5"];
  5:26:s	-> 10:27:s	 [id=13  color="red:blue;0.5"];
  4:28:s	-> 11:29:s	 [id=14  color="red:blue;0.5"];
  7:30:s	-> 11:31:s	 [id=15  color="red:blue;0.5"];
}
"##;

    #[test]
    fn forest_dot_matches_the_base_revision_byte_for_byte() -> Result<()> {
        test_initialise().unwrap();
        let graph: Graph = finalized_runtime_dot!(
            digraph G{
                graph [projector=1]
                node [num=1]
                edge [particle="scalar_1" num=1];
                v1 -> v1 [lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:1}"];
                v1 -> v2 [source="{ufo_order:2}" sink="{ufo_order:0}"];
                v2 -> v3 [lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:0}"];
                v2 -> v3 [source="{ufo_order:2}" sink="{ufo_order:1}"];
                v3 -> v3 [lmb_id=2 source="{ufo_order:2}" sink="{ufo_order:3}"];
            },"scalars"
        )?;
        let forests = Wood::new(
            CutStructure::empty(&graph),
            &graph,
            &scalar_model(),
            &GenerationSettings::default(),
        )?
        .unfold();
        assert_eq!(
            forests.dot_serialize(),
            BASE_NESTED_FACTORISED_FOREST_DOT
                .strip_prefix('\n')
                .expect("the fixture starts on its own line")
        );
        Ok(())
    }
}

#[cfg(test)]
#[path = "overlap_poles.rs"]
pub(crate) mod overlap_poles;
