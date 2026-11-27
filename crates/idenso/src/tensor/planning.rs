//! Depth-limited planning graphs. Their leaves borrow a common evidence owner
//! until an operation needs the leaf as an independent symbolic tensor.

use super::composition::ExplicitIndexOccurrences;
use super::inference::TensorInferenceError;
use super::{SymbolicNet, SymbolicTensor, TensorProofs};
use delegate::delegate;
use linnet::half_edge::NodeIndex;
use spenso::network::graph::{NetworkLeaf, NetworkNode};
use spenso::network::library::panicing::ErroringLibrary;
use spenso::network::library::symbolic::ETS;
use spenso::network::library::{DummyKey, DummyLibrary};
use spenso::network::parsing::{
    LeafInterfaceCache, ParseSettings, StrictTensorFilter, TensorFromExpression,
};
use spenso::network::store::{NetworkStore, TensorScalarStore, TensorScalarStoreMapping};
use spenso::network::{Network, TensorNetworkError};
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::dimension::Dimension;
use spenso::structure::partial::PartialStructure;
use spenso::structure::representation::{LibraryRep, LibrarySlot, Representation};
use spenso::structure::slot::SlotMatcher;
use spenso::structure::{
    Canonicalized, HasStructure, OrderedStructure, Reindexed, SlotIndex, StructureError,
    TensorShell, TensorStructure,
};
use std::ops::Deref;
use std::sync::Arc;
use symbolica::atom::{AddView, Atom, AtomCore, AtomOrView, AtomView, Symbol};

type StoredPlan<E> = Network<NetworkStore<PlanningLeaf<E>, E>, DummyKey, Symbol>;

#[derive(Debug)]
pub(crate) struct PlanningNetwork<E> {
    pub(crate) network: StoredPlan<E>,
    pub(crate) admission: Option<Arc<TensorProofs>>,
}

impl<E> Deref for PlanningNetwork<E> {
    type Target = StoredPlan<E>;
    fn deref(&self) -> &Self::Target {
        &self.network
    }
}

impl<E> PlanningNetwork<E> {
    pub(crate) fn evidence<'a>(&'a self, leaf: &'a PlanningLeaf<E>) -> &'a TensorProofs {
        leaf.evidence(self.admission.as_deref())
    }

    pub(super) fn from_materialized(network: SymbolicNet<AbstractIndex, E>) -> Self {
        Self {
            network: network.map(|scalar| scalar, PlanningLeaf::from_tensor),
            admission: None,
        }
    }
}

/// Facts available at final publication for an exact factor of the shallow plan.
pub(super) struct EnclosingLeaf<'a> {
    pub(super) expression: AtomView<'a>,
    pub(super) boundary_published: bool,
}

impl<E: AtomCore> PlanningNetwork<E> {
    /// An exact stored leaf whose descendants and logical interface are already
    /// supplied by the final output's source. Local and graph-relabeled leaves
    /// must publish their own facts before the plan is discarded.
    pub(super) fn enclosing_leaf(
        &self,
        node: NodeIndex,
        interface: &PartialStructure,
        source: &TensorProofs,
    ) -> Option<EnclosingLeaf<'_>> {
        // The outer shallow plan stores plain scalars separately from planning
        // tensors. Their analyzed boundary is empty and any raw descendants
        // already belong to the enclosing tensor. Its admission may still lack
        // an explicit regional declaration for the scalar itself.
        let index = match self.graph.graph[node] {
            NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => {
                let expression = self.store.get_scalar_ref(index).as_atom_view();
                return Some(EnclosingLeaf {
                    expression,
                    boundary_published: source.published_interface_matches(expression, interface),
                });
            }
            NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => index,
            _ => return None,
        };
        let owner = self.admission.as_deref()?;
        let leaf = &self.store.tensors[index];
        let value = leaf.expression.as_atom_view();
        (leaf.local.is_none()
            && Arc::ptr_eq(
                owner.admitted_interfaces.get()?,
                source.admitted_interfaces.get()?,
            )
            && owner.leaf_interfaces.same_source(&source.leaf_interfaces)
            && source.published_interface_matches(value, interface))
        .then_some(EnclosingLeaf {
            expression: value,
            boundary_published: true,
        })
    }

    /// Graph layout applies to a standalone parser leaf. When emitting a
    /// larger expression, retained descendants keep their declared logical
    /// order. A genuinely new local boundary still takes precedence.
    pub(super) fn inherit_leaf_interfaces<'leaf>(
        &self,
        target: &mut TensorProofs,
        leaves: impl IntoIterator<Item = &'leaf PlanningLeaf<E>>,
    ) where
        E: 'leaf,
    {
        let mut sources = PlanningLeafSources::default();
        sources.extend(self, leaves);
        sources.publish(target);
    }
}

/// Keep precedence across leaves, including when their graphs are consumed
/// one summand at a time. Common admission is retained only once.
#[derive(Default)]
pub(super) struct PlanningLeafSources {
    local: Vec<Arc<TensorProofs>>,
    enclosing: Option<Arc<TensorProofs>>,
    layout_only: Vec<Arc<TensorProofs>>,
}

impl PlanningLeafSources {
    pub(super) fn extend<'a, E: AtomCore + 'a>(
        &mut self,
        network: &PlanningNetwork<E>,
        leaves: impl IntoIterator<Item = &'a PlanningLeaf<E>>,
    ) {
        for leaf in leaves {
            if let Some(facts) = &leaf.local {
                let value = leaf.expression.as_atom_view();
                let layout_only = network
                    .admission
                    .as_deref()
                    .and_then(|owner| owner.published_interface(value))
                    .zip(facts.published_interface(value))
                    .is_some_and(|(declared, parsed)| declared.canonical() == parsed.canonical());
                if !layout_only {
                    self.local.push(Arc::clone(facts));
                    continue;
                }
                self.layout_only.push(Arc::clone(facts));
            }
            if self.enclosing.is_none() {
                self.enclosing.clone_from(&network.admission);
            }
        }
    }

    pub(super) fn append(&mut self, other: Self) {
        self.local.extend(other.local);
        self.layout_only.extend(other.layout_only);
        if self.enclosing.is_none() {
            self.enclosing = other.enclosing;
        }
    }

    pub(super) fn publish(&self, target: &mut TensorProofs) {
        target.inherit_interfaces(
            self.local
                .iter()
                .chain(&self.enclosing)
                .chain(&self.layout_only)
                .map(Arc::as_ref),
        );
    }
}

impl PlanningNetwork<AtomOrView<'_>> {
    pub(super) fn into_owned(self) -> PlanningNetwork<Atom> {
        PlanningNetwork {
            network: self
                .network
                .map(AtomOrView::into_owned, PlanningLeaf::into_owned),
            admission: self.admission,
        }
    }
}

#[derive(Clone, Debug)]
pub(crate) struct PlanningLeaf<E> {
    pub(crate) expression: E,
    pub(crate) structure: OrderedStructure,
    pub(crate) history: Arc<ExplicitIndexOccurrences>,
    // Only a boundary whose freshly established facts differ from its enclosing
    // admission owns regional evidence at the leaf.
    pub(super) local: Option<Arc<TensorProofs>>,
    pub(super) is_metric: bool,
    pub(super) is_composite: bool,
}

impl<E> PlanningLeaf<E> {
    pub(crate) fn evidence<'a>(&'a self, admission: Option<&'a TensorProofs>) -> &'a TensorProofs {
        self.local
            .as_deref()
            .or(admission)
            .expect("plan has admission for every symbolic leaf")
    }

    fn from_tensor(tensor: SymbolicTensor<OrderedStructure, E>) -> Self {
        let history = tensor
            .proofs
            .occurrences
            .get()
            .expect("validated parser leaf history");
        let history = Arc::clone(history);
        let local = Some(Arc::new(tensor.proofs));
        Self {
            expression: tensor.expression,
            structure: tensor.structure,
            history,
            local,
            is_metric: tensor.is_metric,
            is_composite: tensor.is_composite,
        }
    }
}

/// The planner may open one checked occurrence inside a larger immutable source.
#[derive(Clone)]
pub(crate) struct PlanningSource<'src> {
    facts: Arc<TensorProofs>,
    selected: Option<(AtomView<'src>, Arc<ExplicitIndexOccurrences>)>,
}

impl PlanningSource<'_> {
    pub(crate) fn root(proofs: &TensorProofs) -> Self {
        Self {
            facts: Arc::new(proofs.clone()),
            selected: None,
        }
    }

    fn history(&self, value: AtomView<'_>) -> Option<&ExplicitIndexOccurrences> {
        match &self.selected {
            Some(_) => self.selected_history(value).map(Arc::as_ref),
            None => self.facts.parser_history(value),
        }
    }

    fn selected_history(&self, value: AtomView<'_>) -> Option<&Arc<ExplicitIndexOccurrences>> {
        self.selected
            .as_ref()
            .filter(|(selected, _)| value.get_data() == selected.get_data())
            .map(|(_, history)| history)
    }
}

impl<E: AtomCore> PlanningLeaf<E> {
    pub(crate) fn source<'src>(
        &self,
        value: AtomView<'src>,
        admission: Option<&Arc<TensorProofs>>,
    ) -> PlanningSource<'src> {
        debug_assert_eq!(value.get_data(), self.expression.as_atom_view().get_data());
        let facts = self
            .local
            .as_ref()
            .or(admission)
            .expect("plan has admission for every symbolic leaf");
        PlanningSource {
            facts: Arc::clone(facts),
            // A graph built from remapped tensors has only standalone owners.
            // Their construction history is verified on the first parser pass.
            selected: admission.map(|_| (value, Arc::clone(&self.history))),
        }
    }
}

impl PlanningLeaf<AtomOrView<'_>> {
    pub(super) fn into_owned(self) -> PlanningLeaf<Atom> {
        PlanningLeaf {
            expression: self.expression.into_owned(),
            structure: self.structure,
            history: self.history,
            local: self.local,
            is_metric: self.is_metric,
            is_composite: self.is_composite,
        }
    }
}

impl<'src> PlanningLeaf<AtomOrView<'src>> {
    fn from_parser(
        expression: AtomOrView<'src>,
        structure: Canonicalized<OrderedStructure>,
        admission: Option<&PlanningSource<'src>>,
        raw_filter: Option<StrictTensorFilter>,
    ) -> Result<Canonicalized<Self>, TensorInferenceError> {
        let value = expression.as_view();
        let storage = structure.canonical();
        let boundary = TensorProofs::parser_boundary(storage)?;
        let (history, local) = if let Some(owner) = admission
            && owner.facts.published_interface_matches(value, &boundary)
        {
            let history = if let Some(checked) = owner.selected_history(value) {
                TensorProofs::with_parser_ports(Arc::clone(checked), storage)?
            } else {
                owner
                    .facts
                    .parser_occurrences(value, storage, owner.history(value))?
            };
            (history, None)
        } else {
            let proofs = if let Some(owner) = admission {
                owner
                    .facts
                    .admitted_parser_leaf(value, storage, boundary, owner.history(value))?
            } else if let Some(filter) = raw_filter
                && !matches!(value, AtomView::Fun(_))
            {
                let raw = TensorProofs::raw_parser_regions(value, filter)?;
                raw.admitted_parser_leaf(value, storage, boundary, None)?
            } else {
                TensorProofs::parsed_leaf(value, &structure, None)?
            };
            let history = Arc::clone(proofs.occurrences.get().expect("parser leaf history"));
            (history, Some(Arc::new(proofs)))
        };
        let is_composite = !matches!(value, AtomView::Fun(_));
        let is_metric =
            matches!(value, AtomView::Fun(function) if function.get_symbol() == ETS.metric);
        Ok(structure.map_canonical(|structure| Self {
            expression,
            structure,
            history,
            local,
            is_metric,
            is_composite,
        }))
    }
}

// The parser joins storage axes. The source expression and its checked history
// identify the original occurrence, independently of graph-side axis relabelling.
impl<E> TensorStructure for PlanningLeaf<E> {
    type Slot = LibrarySlot<AbstractIndex>;
    type Indexed = OrderedStructure;
    fn reindex_storage(
        self,
        indices: &[AbstractIndex],
    ) -> Result<Reindexed<Self::Indexed>, StructureError> {
        self.structure.reindex_storage(indices)
    }
    fn dual(mut self) -> Self {
        self.structure = self.structure.dual();
        self
    }
    delegate! {
        to self.structure {
            fn is_fully_self_dual(&self) -> bool;
            fn external_structure_iter(&self) -> impl Iterator<Item=Self::Slot>;
            fn external_dims_iter(&self) -> impl Iterator<Item=Dimension>;
            fn external_reps_iter(&self) -> impl Iterator<Item=Representation<LibraryRep>>;
            fn external_indices_iter(&self) -> impl Iterator<Item=AbstractIndex>;
            fn get_aind(&self, i: impl Into<SlotIndex>) -> Option<AbstractIndex>;
            fn get_rep(&self, i: impl Into<SlotIndex>) -> Option<Representation<LibraryRep>>;
            fn get_dim(&self, i: impl Into<SlotIndex>) -> Option<Dimension>;
            fn get_slot(&self, i: impl Into<SlotIndex>) -> Option<Self::Slot>;
            fn order(&self) -> usize;
        }
    }
}

impl<E> HasStructure for PlanningLeaf<E> {
    const PRESERVE_SCALAR_METADATA: bool = true;
    type Structure = OrderedStructure;
    type Scalar = E;
    type ScalarRef<'a>
        = &'a E
    where
        Self: 'a;
    type Store<S: TensorStructure> = TensorShell<S>;
    fn structure(&self) -> &Self::Structure {
        &self.structure
    }
    fn mut_structure(&mut self) -> &mut Self::Structure {
        &mut self.structure
    }
    fn map_structure<S: TensorStructure>(
        self,
        f: impl FnOnce(Self::Structure) -> S,
    ) -> Self::Store<S> {
        f(self.structure).to_shell()
    }
    fn map_structure_result<S: TensorStructure, Er>(
        self,
        f: impl FnOnce(Self::Structure) -> Result<S, Er>,
    ) -> Result<Self::Store<S>, Er> {
        Ok(f(self.structure)?.to_shell())
    }
    fn map_same_structure(mut self, f: impl FnOnce(Self::Structure) -> Self::Structure) -> Self {
        self.structure = f(self.structure);
        self
    }
    fn scalar(self) -> Option<E> {
        self.structure.is_scalar().then_some(self.expression)
    }
    fn scalar_ref(&self) -> Option<&E> {
        self.structure.is_scalar().then_some(&self.expression)
    }
}

impl<'src>
    TensorFromExpression<
        'src,
        OrderedStructure,
        AtomOrView<'src>,
        DummyKey,
        Symbol,
        AbstractIndex,
        DummyLibrary<Self>,
        ErroringLibrary<Symbol>,
    > for PlanningLeaf<AtomOrView<'src>>
{
    type Admission = PlanningSource<'src>;
    fn admitted_structure(
        value: AtomView<'_>,
        owner: &Self::Admission,
        next_index: impl FnMut() -> AbstractIndex,
    ) -> Result<Canonicalized<OrderedStructure>, StructureError> {
        owner
            .facts
            .parser_structure(value, &mut SlotMatcher::default(), next_index)
    }
    fn admitted_structure_with_matcher(
        value: AtomView<'_>,
        owner: &Self::Admission,
        matcher: &mut SlotMatcher,
        next_index: impl FnMut() -> AbstractIndex,
    ) -> Result<Canonicalized<OrderedStructure>, StructureError> {
        owner.facts.parser_structure(value, matcher, next_index)
    }
    fn tensor_from_leaf(
        expression: AtomOrView<'src>,
        structure: Canonicalized<OrderedStructure>,
        admission: Option<&Self::Admission>,
    ) -> Result<Canonicalized<Self>, TensorNetworkError<DummyKey, Symbol>> {
        Self::from_parser(expression, structure, admission, None)
            .map_err(|error| TensorNetworkError::Other(eyre::eyre!(error)))
    }
    fn tensor_from_expression(
        expression: AtomOrView<'src>,
        structure: Canonicalized<OrderedStructure>,
        _: &DummyLibrary<Self>,
        _: &ErroringLibrary<Symbol>,
        settings: &ParseSettings,
        admission: Option<&Self::Admission>,
        _: Option<&LeafInterfaceCache>,
    ) -> Result<Canonicalized<Self>, TensorNetworkError<DummyKey, Symbol>> {
        Self::from_parser(
            expression,
            structure,
            admission,
            Some(settings.strict_tensor_filter),
        )
        .map_err(|error| TensorNetworkError::Other(eyre::eyre!(error)))
    }
}

pub(super) fn parse_summands<'src>(
    sum: AddView<'src>,
    owner: PlanningSource<'src>,
) -> impl Iterator<Item = Result<PlanningNetwork<AtomOrView<'src>>, TensorInferenceError>> {
    let admission = Arc::clone(&owner.facts);
    let cache = admission.leaf_interfaces.read_cache().unwrap_or_default();
    StoredPlan::admitted_summands::<OrderedStructure, _>(sum, DummyLibrary::new(), cache, owner)
        .map(move |network| {
            network
                .map(|network| PlanningNetwork {
                    network,
                    admission: Some(Arc::clone(&admission)),
                })
                .map_err(TensorInferenceError::from)
        })
}

pub(crate) fn parse<'src>(
    source: AtomView<'src>,
    interfaces: Option<&LeafInterfaceCache>,
    admission: Option<PlanningSource<'src>>,
    settings: &ParseSettings,
) -> Result<PlanningNetwork<AtomOrView<'src>>, TensorInferenceError> {
    type Tensor<'a> = PlanningLeaf<AtomOrView<'a>>;
    let library = DummyLibrary::<Tensor<'src>>::new();
    let network = if interfaces.is_some() || admission.is_some() {
        StoredPlan::<AtomOrView<'src>>::try_from_admitted_view_with_context::<OrderedStructure, _>(
            source,
            &library,
            settings,
            interfaces
                .cloned()
                .or_else(|| {
                    admission
                        .as_ref()
                        .map(|owner| owner.facts.leaf_interfaces.read_cache().unwrap_or_default())
                })
                .unwrap_or_default(),
            admission.as_ref(),
        )
    } else {
        StoredPlan::<AtomOrView<'src>>::try_from_view::<OrderedStructure, _>(
            source, &library, settings,
        )
    }?;
    Ok(PlanningNetwork {
        network,
        admission: admission.map(|source| source.facts),
    })
}
