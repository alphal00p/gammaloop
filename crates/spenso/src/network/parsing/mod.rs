use symbolica::atom::Symbol;

use super::*;

use crate::network::library::DummyLibrary;
use crate::network::library::FunctionLibrary;
use crate::network::library::panicing::ErroringLibrary;
use crate::network::profile::{self, Counter, Timer};
use crate::network::tags::SPENSO_TAG;

use crate::structure::abstract_index::AbstractIndex;
use crate::structure::representation::Representation;
use crate::structure::slot::{DummyAind, ParseableAind, Slot, SlotMatcher};
use crate::structure::{Canonicalized, NamedStructure, ScalarStructure, StructureError};
use crate::tensors::parametric::ParamTensor;

use std::{
    cell::{Cell, RefCell},
    collections::HashSet,
    fmt::Display,
    marker::PhantomData,
    rc::Rc,
};

use store::TensorScalarStore;
// use log::trace;

use symbolica::atom::{
    AddView, Atom, AtomOrView, AtomView, MulView, PowView, representation::FunView,
};

use crate::structure::{HasStructure, TensorStructure};

use crate::structure::representation::LibraryRep;

pub type ShadowedStructure<Aind> = NamedStructure<Symbol, Vec<Atom>, LibraryRep, Aind>;

mod construction;
use construction::Construction;
use linnet::half_edge::NodeIndex;
mod indices;
pub(crate) mod structure_inference;
pub use structure_inference::{AtomStructureExt, ChainNestingError, StructureFromAtom};
mod materialization;
mod tensor_from_expression;
pub use tensor_from_expression::{TensorFromExpression, TensorLibraryFor};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ShorthandParsing {
    /// Expand selected shorthand notation into explicit network structure.
    Expand {
        /// Controls compact Schoonschip vector materialization.
        schoonschip: SchoonschipExpansionMode,
        /// Expand `trace(...)` topology into explicit closed links.
        trace: bool,
        /// Expand `chain(...)` topology into explicit open links.
        chain: bool,
    },
    /// Keep shorthand notation as a leaf and infer its exposed structure.
    /// Scalar markers retain their literal wrapper around opaque metadata.
    Opaque,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SchoonschipExpansionMode {
    /// Expand compact vector products such as `dot(p(rep), q(rep))` or
    /// `g(p(rep), q(rep))` into explicit rank-one tensor factors.
    pub inner_products: bool,

    /// Expand compact vector arguments, for example
    /// `T(..., p(rep), ...)` into `T(..., slot, ...) * p(slot)`.
    pub expand_schoonship: bool,

    /// Expand compact Schoonschip notation inside `chain(...)` and `trace(...)`
    /// factors.
    pub expand_inside_chains: bool,
}

impl Default for SchoonschipExpansionMode {
    fn default() -> Self {
        Self::full()
    }
}

impl SchoonschipExpansionMode {
    /// Expand every compact Schoonschip shorthand recognized by the parser.
    pub const fn full() -> Self {
        Self {
            inner_products: true,
            expand_schoonship: true,
            expand_inside_chains: true,
        }
    }

    /// Keep all compact Schoonschip shorthand opaque.
    pub const fn none() -> Self {
        Self {
            inner_products: false,
            expand_schoonship: false,
            expand_inside_chains: false,
        }
    }

    /// Return this mode with chain/trace factor materialization disabled.
    pub const fn outside_chains(self) -> Self {
        Self {
            expand_inside_chains: false,
            ..self
        }
    }

    fn any(self) -> bool {
        self.inner_products || self.expand_schoonship || self.expand_inside_chains
    }

    fn for_chain_like_root(self) -> Self {
        if self.expand_inside_chains {
            self
        } else {
            Self::none()
        }
    }
}

impl Default for ShorthandParsing {
    fn default() -> Self {
        Self::expand_all()
    }
}

impl ShorthandParsing {
    /// Expand all shorthand families. This matches the historical default.
    pub const fn expand_all() -> Self {
        Self::Expand {
            schoonschip: SchoonschipExpansionMode::full(),
            trace: true,
            chain: true,
        }
    }

    /// Expand only compact Schoonschip notation while keeping chain and trace
    /// topology opaque.
    pub const fn expand_schoonschip_only() -> Self {
        Self::Expand {
            schoonschip: SchoonschipExpansionMode::full(),
            trace: false,
            chain: false,
        }
    }

    pub fn expands(self) -> bool {
        matches!(self, Self::Expand { .. })
    }

    pub(super) fn expands_chain(self) -> bool {
        matches!(self, Self::Expand { chain: true, .. })
    }

    pub(super) fn expands_trace(self) -> bool {
        matches!(self, Self::Expand { trace: true, .. })
    }

    pub(super) fn schoonschip_expansion(self) -> Option<SchoonschipExpansionMode> {
        match self {
            Self::Expand { schoonschip, .. } => Some(schoonschip),
            Self::Opaque => None,
        }
    }

    pub(super) fn with_schoonschip_expansion(self, schoonschip: SchoonschipExpansionMode) -> Self {
        match self {
            Self::Expand { trace, chain, .. } => Self::Expand {
                schoonschip,
                trace,
                chain,
            },
            Self::Opaque => Self::Opaque,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum StrictTensorFilter {
    /// Only parser syntax and heads tagged as tensors are parsed as tensorial.
    #[default]
    Tagged,
    /// Tensor-tagged heads must also expose representation syntax in their arguments.
    TaggedChecked,
    /// Any head with representation syntax in its arguments is parsed as tensorial.
    ContainsReps,
}

#[derive(Clone, Debug)]
pub struct ParseSettings {
    /// Fold factors that parse as pure scalars into one scalar factor while
    /// parsing a product.
    ///
    /// This keeps scalar-only subexpressions out of the tensor graph when they
    /// cannot affect contraction topology. Disable it when the caller needs
    /// every product factor to remain represented separately in the network.
    pub precontract_scalars: bool,

    /// Parse only the first summand of an addition.
    ///
    /// This is meant for structure discovery paths where all summands are
    /// expected to have the same external structure, for example listing
    /// dangling indices without building a full sum network.
    pub take_first_term_from_sum: bool,

    /// Stop recursive parsing once product nesting reaches this depth.
    ///
    /// At the limit, the current expression is handed to the opaque tensor
    /// expression boundary as a leaf. `None` means there is no depth limit.
    pub depth_limit: Option<usize>,

    /// Selects how shorthand notation is represented in the parsed network.
    ///
    /// `Expand` turns the selected shorthand families into explicit graph
    /// structure. Fresh dummies created by this lowering are local to the
    /// expansion.
    ///
    /// `Opaque` keeps shorthand as a leaf tensor or scalar with its exposed
    /// structure determined by syntactic observation.
    pub shorthand_parsing: ShorthandParsing,

    /// Allow parser implementations to treat composite scalar expressions as
    /// scalar-structured tensors.
    ///
    /// Depth-limited additive or multiplicative scalar composites use this to
    /// stay represented as scalar-structured tensor leaves, so recursive
    /// contraction passes can inspect their internals instead of storing them
    /// as pure scalars.
    pub parse_composite_scalars_as_tensors: bool,

    /// Controls which ordinary function heads are eligible tensor leaves.
    ///
    /// Parser-owned syntax such as slots, representations, `chain`, `trace`,
    /// `dot`, metric heads, and transparent brackets keeps its fixed meaning.
    /// This setting decides how strict the parser is for ordinary tensor heads:
    /// require tensor tags, require tensor tags plus visible representation
    /// arguments, or accept untagged heads that contain representation syntax.
    pub strict_tensor_filter: StrictTensorFilter,
}

impl Default for ParseSettings {
    fn default() -> Self {
        Self {
            precontract_scalars: true,
            take_first_term_from_sum: false,
            depth_limit: None,
            shorthand_parsing: ShorthandParsing::default(),
            parse_composite_scalars_as_tensors: false,
            strict_tensor_filter: StrictTensorFilter::Tagged,
        }
    }
}

impl ParseSettings {
    /// Return these settings with the ordinary tensor-head filter replaced.
    pub fn with_strict_tensor_filter(mut self, filter: StrictTensorFilter) -> Self {
        self.strict_tensor_filter = filter;
        self
    }

    fn with_schoonschip_expansion(mut self, mode: SchoonschipExpansionMode) -> Self {
        self.shorthand_parsing = self.shorthand_parsing.with_schoonschip_expansion(mode);
        self
    }
}

#[derive(Clone)]
pub struct ParseState<Aind = AbstractIndex, View = ()> {
    view: View,
    chain_scope_validated: bool,
    depth: usize,
    matcher: Rc<RefCell<SlotMatcher>>,
    metric: Symbol,
    next_dummy: Rc<Cell<usize>>,
    reserved_indices: Rc<RefCell<HashSet<Atom>>>,
    _aind: PhantomData<fn() -> Aind>,
}

impl<Aind, View> Debug for ParseState<Aind, View> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("ParseState")
            .field("depth", &self.depth)
            .field("next_dummy", &self.next_dummy)
            .field("reserved_indices", &self.reserved_indices)
            .finish_non_exhaustive()
    }
}

#[allow(clippy::derivable_impls)]
impl<Aind> Default for ParseState<Aind> {
    fn default() -> Self {
        Self {
            view: (),
            chain_scope_validated: false,
            depth: 0,
            matcher: Rc::new(RefCell::new(SlotMatcher::default())),
            metric: crate::network::library::symbolic::ETS.metric,
            next_dummy: Rc::new(Cell::new(1_000_000)),
            reserved_indices: Rc::default(),
            _aind: PhantomData,
        }
    }
}

impl<Aind, View> ParseState<Aind, View> {
    // Rebinding an original subtree moves the operation owners, without adding
    // per-leaf reference-count traffic. Only the dispatcher sets the view.
    fn with_view<NextView>(self, view: NextView) -> ParseState<Aind, NextView> {
        ParseState {
            view,
            chain_scope_validated: self.chain_scope_validated,
            depth: self.depth,
            matcher: self.matcher,
            metric: self.metric,
            next_dummy: self.next_dummy,
            reserved_indices: self.reserved_indices,
            _aind: self._aind,
        }
    }

    // Callback normalization and shorthand lowering produce new syntax. Its
    // leaf admission must retain the ordinary checked scope boundary.
    fn materialized(self) -> ParseState<Aind> {
        let mut state = self.with_view(());
        state.chain_scope_validated = false;
        state
    }
}

impl<'src, Aind> ParseState<Aind, AtomView<'src>> {
    /// The exact subtree currently undergoing parser leaf inference.
    pub fn current_view(&self) -> AtomView<'src> {
        self.view
    }
}

impl<Aind: DummyAind + ParseableAind, View> ParseState<Aind, View> {
    /// Reserve written index names once at an operation's input boundary.
    /// Cloned parser states share both these names and subsequent allocations.
    pub fn reserve_indices(&self, value: AtomView<'_>) {
        let mut matcher = self.matcher.borrow_mut();
        value.visitor(&mut |atom| {
            if let Ok(slot) = matcher.parse::<LibraryRep, Aind>(atom) {
                self.reserve_index(slot.aind);
            }
            true
        });
    }

    /// Reserve the serialized identity, including explicit names resembling dummies.
    pub fn reserve_index(&self, index: Aind) {
        self.reserved_indices.borrow_mut().insert(index.to_atom());
    }

    /// Allocate an operation-wide collision-free name that also remains fresh
    /// across independently constructed tensors. Parser-local deterministic
    /// materialization continues to use `next` and the same reservation set.
    pub fn fresh_index(&self) -> Aind {
        loop {
            let atom = Aind::new_dummy().to_atom();
            if self.reserved_indices.borrow_mut().insert(atom.clone()) {
                return Aind::from_view(atom.as_view()).unwrap_or_else(|_| {
                    panic!("fresh dummy atoms must remain valid abstract indices")
                });
            }
        }
    }

    fn next(&self) -> Aind {
        loop {
            let index = self.next_dummy.get();
            self.next_dummy.set(index + 1);
            let dummy = Aind::new_dummy_at(index);
            if self.reserved_indices.borrow_mut().insert(dummy.to_atom()) {
                return dummy;
            }
        }
    }

    fn slot(&self, rep: &Representation<LibraryRep>) -> Slot<LibraryRep, Aind> {
        rep.slot(self.next())
    }
}

impl<
    'src,
    Sc,
    T: HasStructure + TensorStructure,
    K: Clone + Display + Debug,
    // FK: Clone + Display + Debug,
    Str: TensorScalarStore<Tensor = T, Scalar = Sc> + Clone,
    Aind: AbsInd + DummyAind + ParseableAind,
> Network<Str, K, Symbol, Aind>
where
    Sc: TryFrom<AtomView<'src>> + TryFrom<Atom> + Clone,
    TensorNetworkError<K, Symbol>:
        From<<Sc as TryFrom<AtomView<'src>>>::Error> + From<<Sc as TryFrom<Atom>>::Error>,
{
    #[allow(clippy::result_large_err)]
    pub fn try_from_view<S, Lib: TensorLibraryFor<S, T, Key = K>>(
        value: AtomView<'src>,
        library: &Lib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, ErroringLibrary<Symbol>>,
    {
        Self::try_from_view_with_function_library(
            value,
            library,
            &ErroringLibrary::<Symbol>::new(),
            settings,
        )
    }

    /// Return the exposed ports using the same parser and materialization policy,
    /// without constructing or finalizing a network graph.
    #[allow(clippy::result_large_err)]
    pub fn try_external_slots<S, Lib>(
        value: AtomView<'src>,
        library: &Lib,
        settings: &ParseSettings,
    ) -> Result<
        Vec<crate::structure::representation::LibrarySlot<Aind>>,
        TensorNetworkError<K, Symbol>,
    >
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, ErroringLibrary<Symbol>>,
        Lib: TensorLibraryFor<S, T, Key = K>,
    {
        let (construction, root) = Self::parse_root::<S, _, _>(
            value,
            library,
            &ErroringLibrary::<Symbol>::new(),
            settings,
            false,
        )?;
        Ok(construction.slots(root))
    }

    #[allow(clippy::result_large_err)]
    pub fn try_from_view_with_function_library<S, Lib, FunLib>(
        value: AtomView<'src>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        let (construction, root) =
            Self::parse_root::<S, _, _>(value, library, function_library, settings, true)?;
        Ok(construction.finish(root))
    }

    #[allow(clippy::type_complexity)]
    fn parse_root<S, Lib, FunLib>(
        value: AtomView<'src>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        emit_graph: bool,
    ) -> Result<(Construction<Str, K, Aind>, NodeIndex), TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        value.validate_chain_like_nesting()?;
        let mut construction = if emit_graph {
            Construction::new()
        } else {
            Construction::structure_only()
        };
        let state = ParseState::<Aind> {
            chain_scope_validated: true,
            ..Default::default()
        };
        // Reserve names before any shorthand allocates a dummy.
        state.reserve_indices(value);
        let root = Self::try_from_view_impl(
            &mut construction,
            value,
            state,
            library,
            function_library,
            settings,
            AtomOrView::View,
        )?;
        Ok((construction, root))
    }

    #[allow(clippy::result_large_err)]
    fn scalar_from_expression(
        value: AtomOrView<'src>,
    ) -> Result<Sc, TensorNetworkError<K, Symbol>> {
        Ok(match value {
            AtomOrView::View(view) => view.try_into()?,
            owned => owned.into_owned().try_into()?,
        })
    }

    fn try_from_view_impl<'node, S, Lib, FunLib, View>(
        construction: &mut Construction<Str, K, Aind>,
        value: AtomView<'node>,
        state: ParseState<Aind, View>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<NodeIndex, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        let state = state.with_view(value);
        profile::bump(Counter::ParseView, 1);
        match value {
            AtomView::Mul(m) => {
                profile::bump(Counter::ParseMul, 1);
                Self::try_from_mul(
                    construction,
                    m,
                    state,
                    library,
                    function_library,
                    settings,
                    retain,
                )
            }
            AtomView::Fun(f) => {
                profile::bump(Counter::ParseFun, 1);
                Self::try_from_fun(
                    construction,
                    f,
                    state,
                    library,
                    function_library,
                    settings,
                    retain,
                )
            }
            AtomView::Add(a) => {
                profile::bump(Counter::ParseAdd, 1);
                Self::try_from_add(
                    construction,
                    a,
                    state,
                    library,
                    function_library,
                    settings,
                    retain,
                )
            }
            AtomView::Pow(p) => {
                profile::bump(Counter::ParsePow, 1);
                Self::try_from_pow(
                    construction,
                    p,
                    state,
                    library,
                    function_library,
                    settings,
                    retain,
                )
            }
            a => Ok(construction.scalar(Self::scalar_from_expression(retain(a))?)),
        }
    }

    #[allow(clippy::type_complexity, clippy::result_large_err)]
    fn as_leaf<'node, S, Lib, FunLib>(
        construction: &mut Construction<Str, K, Aind>,
        value: AtomView<'node>,
        state: &ParseState<Aind, AtomView<'node>>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<NodeIndex, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        let scalar_composite = settings.parse_composite_scalars_as_tensors
            && matches!(value, AtomView::Add(_) | AtomView::Mul(_));
        if !scalar_composite
            && !structure_inference::TensorialSyntax::is_tensorial(
                value,
                settings.strict_tensor_filter,
                &state.matcher.borrow(),
            )
        {
            return Ok(construction.scalar(Self::scalar_from_expression(retain(value))?));
        }

        Self::as_inferred_leaf::<S, Lib, FunLib>(
            construction,
            value,
            state,
            library,
            function_library,
            settings,
            retain,
        )
    }

    #[allow(
        clippy::type_complexity,
        clippy::result_large_err,
        clippy::too_many_arguments
    )]
    fn as_inferred_leaf<'node, S, Lib, FunLib>(
        construction: &mut Construction<Str, K, Aind>,
        value: AtomView<'node>,
        state: &ParseState<Aind, AtomView<'node>>,

        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<NodeIndex, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        debug_assert_eq!(state.current_view(), value);
        profile::bump(Counter::ParseStructureAttempt, 1);
        let structure = {
            let _span = profile::span(Timer::ParseStructure);
            S::structure_from_parser(state)
        };

        let structure = match structure {
            Ok(structure) => {
                profile::bump(Counter::ParseStructureOk, 1);
                structure
            }
            Err(StructureError::EmptyStructure(_)) => {
                profile::bump(Counter::ParseStructureErr, 1);
                Canonicalized::identity(S::scalar_structure())
            }
            Err(err) => {
                profile::bump(Counter::ParseStructureErr, 1);
                return Err(err.into());
            }
        };

        let layout = structure.layout().clone();
        Ok(construction.tensor(
            T::tensor_from_expression(
                retain(value),
                structure,
                library,
                function_library,
                settings,
            )?,
            layout,
        )?)
    }

    #[allow(clippy::result_large_err)]
    fn try_from_mul<'node, S, Lib, FunLib>(
        construction: &mut Construction<Str, K, Aind>,
        value: MulView<'node>,
        mut state: ParseState<Aind, AtomView<'node>>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<NodeIndex, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        let _span = profile::span(Timer::ParseMul);
        // println!("Mul");
        if let Some(a) = settings.depth_limit
            && a <= state.depth
        {
            // println!("Mul leaf");
            return Self::as_leaf::<S, Lib, FunLib>(
                construction,
                value.as_view(),
                &state,
                library,
                function_library,
                settings,
                retain,
            );
        }

        state.depth += 1;
        // println!("{} for mul {}", state.depth, value.as_view());
        let mut iter = value.iter();
        let first_atom = iter.next().unwrap();
        profile::bump(Counter::MulFactor, 1);
        let first = Self::try_from_view_impl(
            construction,
            first_atom,
            state.clone(),
            library,
            function_library,
            settings,
            retain,
        )?;

        // state

        if settings.precontract_scalars {
            let mut scalar_terms = Vec::new();

            let rest: Result<Vec<_>, _> = iter
                .filter_map(|a| {
                    match Self::try_from_view_impl(
                        construction,
                        a,
                        state.clone(),
                        library,
                        function_library,
                        settings,
                        retain,
                    ) {
                        Ok(n) => {
                            profile::bump(Counter::MulFactor, 1);
                            if let NetworkState::PureScalar = construction.state(n) {
                                scalar_terms.push(a);
                                None
                            } else {
                                Some(Ok(n))
                            }
                        }
                        Err(e) => Some(Err(e)),
                    }
                })
                .collect();

            let mut res = rest?;

            if let NetworkState::PureScalar = construction.state(first) {
                scalar_terms.push(first_atom);
            } else {
                res.push(first);
            }

            if res.is_empty() {
                Ok(construction.scalar(Self::scalar_from_expression(retain(value.as_view()))?))
            } else {
                // Fold only the coefficient that survives in a mixed network.
                // Keep the previous rest-then-first order, including inexact
                // arithmetic; an all-scalar source is retained unchanged above.
                let mut scalars = Atom::num(1);
                for term in scalar_terms {
                    let _span = profile::span(Timer::ScalarMulAccum);
                    profile::bump(Counter::ScalarMulAccum, 1);
                    scalars *= term;
                }
                let s = if scalars != Atom::num(1) {
                    construction.scalar(Self::scalar_from_expression(AtomOrView::Atom(scalars))?)
                } else {
                    res.pop().unwrap()
                };

                Ok(construction.product(std::iter::once(s).chain(res).collect()))
            }
        } else {
            let rest: Result<Vec<_>, _> = iter
                .map(|a| {
                    profile::bump(Counter::MulFactor, 1);
                    Self::try_from_view_impl(
                        construction,
                        a,
                        state.clone(),
                        library,
                        function_library,
                        settings,
                        retain,
                    )
                })
                .collect();

            Ok(construction.product(std::iter::once(first).chain(rest?).collect()))
        }
    }

    #[allow(clippy::result_large_err)]
    fn try_from_fun<'node, S, Lib, FunLib>(
        construction: &mut Construction<Str, K, Aind>,
        value: FunView<'node>,
        state: ParseState<Aind, AtomView<'node>>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<NodeIndex, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
        // <Canonicalized<S>>::Error: Debug,
    {
        let _span = profile::span(Timer::ParseFun);
        let symbol = value.get_symbol();

        if symbol == state.matcher.borrow().tags().dot && value.get_nargs() != 2 {
            return Err(TensorNetworkError::InvalidDotFunction(
                value.as_view().to_plain_string(),
            ));
        }

        if symbol == state.matcher.borrow().tags().bracket
            || symbol.has_tag(&state.matcher.borrow().tags().broadcast)
        {
            return Self::parse_expanded_function::<S, Lib, FunLib>(
                construction,
                value,
                state,
                library,
                function_library,
                settings,
                retain,
            );
        }

        if !structure_inference::TensorialSyntax::is_tensorial(
            value.as_view(),
            settings.strict_tensor_filter,
            &state.matcher.borrow(),
        ) {
            return Self::parse_scalar_function(
                construction,
                value,
                settings.shorthand_parsing,
                retain,
            );
        }

        if settings.shorthand_parsing == ShorthandParsing::Opaque
            && Self::is_shorthand_function(value, &state)
        {
            return Self::as_inferred_leaf::<S, _, _>(
                construction,
                value.as_view(),
                &state,
                library,
                function_library,
                settings,
                retain,
            );
        }

        Self::parse_expanded_function(
            construction,
            value,
            state,
            library,
            function_library,
            settings,
            retain,
        )
    }

    #[allow(clippy::result_large_err)]
    fn parse_scalar_function<'node>(
        construction: &mut Construction<Str, K, Aind>,
        value: FunView<'node>,
        shorthand_parsing: ShorthandParsing,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<NodeIndex, TensorNetworkError<K, Symbol>> {
        if value.get_symbol() == SPENSO_TAG.pure_scalar {
            if value.get_nargs() != 1 {
                return Err(TensorNetworkError::TooManyArgsFunction(
                    value.as_view().to_plain_string(),
                ));
            }

            // Opaque intake retains the literal for domain rewriting; explicit
            // network materialization consumes the scalar marker.
            if shorthand_parsing.expands() {
                return Ok(construction.scalar(Self::scalar_from_expression(retain(
                    value.iter().next().unwrap(),
                ))?));
            }
        }

        Ok(construction.scalar(Self::scalar_from_expression(retain(value.as_view()))?))
    }

    #[allow(clippy::result_large_err)]
    fn parse_broadcast_function<'node, S, Lib, FunLib>(
        construction: &mut Construction<Str, K, Aind>,
        value: FunView<'node>,
        state: ParseState<Aind, AtomView<'node>>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<NodeIndex, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        if value.get_nargs() != 1 {
            return Err(TensorNetworkError::TooManyArgsFunction(
                value.as_view().to_plain_string(),
            ));
        }

        let symbol = value.get_symbol();
        let mut state = state;
        if symbol.is_scalar() {
            state.chain_scope_validated = false;
        }
        let inner = value.iter().next().unwrap();
        let inner_tensor = Self::try_from_view_impl(
            construction,
            inner,
            state,
            library,
            function_library,
            settings,
            retain,
        )?;

        Ok(construction.function(inner_tensor, symbol))
    }

    #[allow(clippy::result_large_err)]
    fn parse_expanded_function<'node, S, Lib, FunLib>(
        construction: &mut Construction<Str, K, Aind>,
        value: FunView<'node>,
        state: ParseState<Aind, AtomView<'node>>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<NodeIndex, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        let symbol = value.get_symbol();

        if symbol == state.matcher.borrow().tags().bracket {
            let mut n_muls = value
                .iter()
                .map(|a| {
                    Self::try_from_view_impl(
                        construction,
                        a,
                        state.clone(),
                        library,
                        function_library,
                        settings,
                        retain,
                    )
                })
                .collect::<Result<Vec<_>, _>>()?;
            let Some(first) = n_muls.pop() else {
                return Err(TensorNetworkError::Other(eyre::eyre!(
                    "empty bracket expression {}",
                    value.as_view()
                )));
            };
            Ok(construction.product(std::iter::once(first).chain(n_muls).collect()))
        } else if symbol.has_tag(&state.matcher.borrow().tags().broadcast) {
            Self::parse_broadcast_function::<S, Lib, FunLib>(
                construction,
                value,
                state,
                library,
                function_library,
                settings,
                retain,
            )
        } else {
            Self::materialize_shorthand::<S, Lib, FunLib>(
                construction,
                value,
                state,
                library,
                function_library,
                settings,
                retain,
            )
        }
    }

    #[allow(clippy::result_large_err)]
    fn parse_regular_function_leaf<'node, S, Lib, FunLib>(
        construction: &mut Construction<Str, K, Aind>,
        value: FunView<'node>,
        library: &Lib,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<NodeIndex, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        profile::bump(Counter::ParseStructureAttempt, 1);
        let structure = {
            let _span = profile::span(Timer::ParseStructure);
            S::parse(value.as_view())
        };

        let structure = match structure {
            Ok(structure) => {
                profile::bump(Counter::ParseStructureOk, 1);
                structure
            }
            Err(StructureError::EmptyStructure(_)) => {
                profile::bump(Counter::ParseStructureErr, 1);
                return Ok(
                    construction.scalar(Self::scalar_from_expression(retain(value.as_view()))?)
                );
            }
            Err(err) => {
                profile::bump(Counter::ParseStructureErr, 1);
                return Err(err.into());
            }
        };

        match library.key_for_structure(&structure) {
            Ok(key) => {
                let tensor_structure = structure.canonical().clone();
                Ok(
                    construction
                        .library_tensor(&tensor_structure, structure.map_canonical(|_| key)),
                )
            }
            Err(_) if structure.canonical().is_scalar() => {
                Ok(construction.scalar(Self::scalar_from_expression(retain(value.as_view()))?))
            }
            Err(_) => {
                // Tensor dimensions may remain symbolic during analytic normalization;
                // opaque scalars and library tensors need no eager shadow here.
                // The target decides whether missing leaves need finite components.
                let layout = structure.layout().clone();
                Ok(construction.tensor(
                    T::tensor_from_leaf(retain(value.as_view()), structure)?,
                    layout,
                )?)
            }
        }
    }

    #[allow(clippy::result_large_err)]
    fn try_from_pow<'node, S, Lib, FunLib>(
        construction: &mut Construction<Str, K, Aind>,
        value: PowView<'node>,
        state: ParseState<Aind, AtomView<'node>>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> std::result::Result<NodeIndex, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        let _span = profile::span(Timer::ParsePow);
        if let Some(a) = settings.depth_limit
            && a <= state.depth
        {
            return Self::as_leaf::<S, Lib, FunLib>(
                construction,
                value.as_view(),
                &state,
                library,
                function_library,
                settings,
                retain,
            );
        }

        let (base_expression, exp) = value.get_base_exp();

        if let Ok(n) = i8::try_from(exp) {
            // println!("base:{base_expression}");
            let next_dummy = state.next_dummy.get();
            let base = Self::try_from_view_impl(
                construction,
                base_expression,
                state.clone(),
                library,
                function_library,
                settings,
                retain,
            )?;

            // println!("base state {:?}", construction.state(base));
            if settings.precontract_scalars
                && let NetworkState::PureScalar = construction.state(base)
                && !matches!(base_expression, AtomView::Fun(fun) if fun.get_symbol() == state.matcher.borrow().tags().bracket)
            {
                // println!("Pure");
                return Ok(
                    construction.scalar(Self::scalar_from_expression(retain(value.as_view()))?)
                );
            }

            if let NetworkState::Tensor = construction.state(base) {
                return Err(TensorNetworkError::NonSelfDualTensorPower(
                    value.as_view().to_plain_string(),
                ));
            }
            // An even power of a self_dual tensor, or scalar is a scalar
            if n < 0 && n % 2 != 0 && !construction.state(base).is_scalar() {
                let reason = if construction.emits_graph() {
                    let base = std::mem::replace(construction, Construction::new()).finish(base);
                    format!(
                        "Atom:{},graph of base: {}, dangling indices: {:?}",
                        value.as_view().to_plain_string(),
                        base.dot(),
                        base.graph.dangling_indices()
                    )
                } else {
                    format!(
                        "Atom:{}, dangling indices: {:?}",
                        value.as_view().to_plain_string(),
                        construction.slots(base)
                    )
                };
                return Err(TensorNetworkError::NegativeExponentNonScalar(reason));
            }

            let out = if n > 1
                && (state.next_dummy.get() != next_dummy
                    || (n == 2 && construction.state(base) == NetworkState::SelfDualTensor))
            {
                // Each lowered shorthand copy needs independent internal indices.
                // A tensor square also needs the selected product-contraction
                // strategy: scalar power execution uses the tensor's default
                // contraction, which can leave symbolic sums uncontracted.
                let rest = (1..n)
                    .map(|_| {
                        Self::try_from_view_impl(
                            construction,
                            base_expression,
                            state.clone(),
                            library,
                            function_library,
                            settings,
                            retain,
                        )
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                construction.product(std::iter::once(base).chain(rest).collect())
            } else {
                // println!("{:?}", construction.state(base));
                construction.power(base, n)
            };
            // println!("{:?}", out.state);
            Ok(out)
        } else {
            Ok(construction.scalar(Self::scalar_from_expression(retain(value.as_view()))?))
        }
    }

    #[allow(clippy::result_large_err)]
    fn try_from_add<'node, S, Lib, FunLib>(
        construction: &mut Construction<Str, K, Aind>,
        value: AddView<'node>,
        state: ParseState<Aind, AtomView<'node>>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<NodeIndex, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        let _span = profile::span(Timer::ParseAdd);
        if let Some(a) = settings.depth_limit
            && a <= state.depth
        {
            return Self::as_leaf::<S, Lib, FunLib>(
                construction,
                value.as_view(),
                &state,
                library,
                function_library,
                settings,
                retain,
            );
        }

        let mut iter = value.iter();

        let first_atom = iter.next().unwrap();
        profile::bump(Counter::AddTerm, 1);

        let first = Self::try_from_view_impl(
            construction,
            first_atom,
            state.clone(),
            library,
            function_library,
            settings,
            retain,
        )?;
        if settings.take_first_term_from_sum {
            Ok(first)
        } else if settings.precontract_scalars {
            let mut scalar_terms = Vec::new();

            let rest: Result<Vec<_>, _> = iter
                .filter_map(|a| {
                    match Self::try_from_view_impl(
                        construction,
                        a,
                        state.clone(),
                        library,
                        function_library,
                        settings,
                        retain,
                    ) {
                        Ok(n) => {
                            profile::bump(Counter::AddTerm, 1);
                            if construction
                                .state(n)
                                .is_compatible(&construction.state(first))
                            {
                                if let NetworkState::PureScalar = construction.state(n) {
                                    scalar_terms.push(a);
                                    None
                                } else {
                                    Some(Ok(n))
                                }
                            } else {
                                Some(Err(TensorNetworkError::IncompatibleSummand(format!(
                                    "{} is {:?} vs {} is {:?}",
                                    a,
                                    construction.state(n),
                                    first_atom,
                                    construction.state(first)
                                ))))
                            }
                        }
                        Err(e) => Some(Err(e)),
                    }
                })
                .collect();

            let mut res = rest?;

            if let NetworkState::PureScalar = construction.state(first) {
                scalar_terms.push(first_atom);
            } else {
                res.push(first);
            }

            if res.is_empty() {
                Ok(construction.scalar(Self::scalar_from_expression(retain(value.as_view()))?))
            } else {
                // Fold only the coefficient that survives in a mixed network.
                // Keep the previous rest-then-first order, including inexact
                // arithmetic; an all-scalar source is retained unchanged above.
                let mut scalars = Atom::Zero;
                for term in scalar_terms {
                    let _span = profile::span(Timer::ScalarAddAccum);
                    profile::bump(Counter::ScalarAddAccum, 1);
                    scalars += term;
                }
                let s = if scalars != Atom::Zero {
                    construction.scalar(Self::scalar_from_expression(AtomOrView::Atom(scalars))?)
                } else {
                    res.pop().unwrap()
                };
                construction.sum(std::iter::once(s).chain(res).collect())
            }
        } else {
            let rest: Result<Vec<_>, _> = iter
                .map(|a| {
                    match Self::try_from_view_impl(
                        construction,
                        a,
                        state.clone(),
                        library,
                        function_library,
                        settings,
                        retain,
                    ) {
                        Ok(n) => {
                            profile::bump(Counter::AddTerm, 1);
                            if construction
                                .state(n)
                                .is_compatible(&construction.state(first))
                            {
                                Ok(n)
                            } else {
                                Err(TensorNetworkError::IncompatibleSummand(format!(
                                    "{} is {:?} vs {} is {:?}",
                                    a,
                                    construction.state(n),
                                    first_atom,
                                    construction.state(first)
                                )))
                            }
                        }
                        Err(e) => Err(e),
                    }
                })
                .collect();

            construction.sum(std::iter::once(first).chain(rest?).collect())
        }
    }
}
pub type ParamNet<Aind> =
    Network<NetworkStore<ParamTensor<ShadowedStructure<Aind>>, Atom>, DummyKey, Symbol, Aind>;

impl<Aind: AbsInd + DummyAind + ParseableAind + 'static> ParamNet<Aind> {
    /// Realize independent scalar expressions with the existing component parser
    /// and executor. Each input keeps its own dummy namespace and admission error.
    ///
    /// Parsing precedes execution for a batch. Callers whose inputs can run user
    /// code must submit one input at its original evaluation position instead;
    /// the same implementation then preserves parse/execute interleaving. No
    /// unrelated expression or tensor is materialized by this operation.
    #[allow(clippy::type_complexity, clippy::result_large_err)]
    pub fn evaluate_scalar_batch(
        values: &[AtomView<'_>],
        settings: &ParseSettings,
    ) -> Result<
        Vec<Result<Atom, TensorNetworkError<DummyKey, Symbol>>>,
        TensorNetworkError<DummyKey, Symbol>,
    > {
        use crate::network::library::function_lib::SymbolLib;
        use std::{
            collections::HashMap,
            sync::atomic::{AtomicU64, Ordering},
        };
        use symbolica::{atom::SymbolBuilder, function, wrap_symbol};

        let library = DummyLibrary::<_>::new();
        let missing = ErroringLibrary::<Symbol>::new();
        let mut construction = Construction::new();
        let mut results = HashMap::new();
        let mut roots = Vec::new();
        let mut functions = SymbolLib {
            functions: HashMap::new(),
            scalar_functions: HashMap::new(),
            _missing: ErroringLibrary::<Symbol>::new(),
        };
        let mut result_heads = HashMap::new();
        // SymbolBuilder::build_group reserves all names under Symbolica's state
        // lock and rejects existing symbols, including ones with user hooks.
        // Plain symbol! lookup cannot provide that guarantee.
        static NEXT_BATCH: AtomicU64 = AtomicU64::new(0);
        let heads = if values.len() > 1 {
            loop {
                let batch = NEXT_BATCH.fetch_add(1, Ordering::Relaxed);
                let builders = (0..values.len())
                    .map(|position| {
                        SymbolBuilder::new(wrap_symbol!(format!(
                            "spenso::internal_component_result_{batch}_{position}"
                        )))
                    })
                    .collect();
                if let Ok(heads) = SymbolBuilder::build_group(builders) {
                    break heads;
                }
            }
        } else {
            Vec::new()
        };
        for (position, &value) in values.iter().enumerate() {
            let parsed = value
                .validate_chain_like_nesting()
                .map_err(Into::into)
                .and_then(|()| {
                    let state = ParseState::<Aind> {
                        chain_scope_validated: true,
                        ..Default::default()
                    };
                    state.reserve_indices(value);
                    Self::try_from_view_impl::<ShadowedStructure<Aind>, _, _, _>(
                        &mut construction,
                        value,
                        state,
                        &library,
                        &missing,
                        settings,
                        AtomOrView::View,
                    )
                });
            match parsed {
                Err(error) => {
                    results.insert(position, Err(error));
                }
                Ok(root) if values.len() == 1 => {
                    // In particular, execute a callback-produced open result
                    // before reporting NoScalar, just like ordinary execution.
                    let mut network = construction.finish(root);
                    network.execute::<Sequential, SmallestDegree, _, _, _>(&library, &missing)?;
                    return Ok(vec![network.result_scalar().map(Atom::from)]);
                }
                Ok(root) if !construction.state(root).is_scalar() => {
                    results.insert(position, Err(TensorNetworkError::NoScalar));
                }
                Ok(root) => {
                    // Private operation keys keep independent results separate;
                    // they never enter input syntax or escape this method.
                    let head = heads[position];
                    functions.insert_scalar_fallible(head, move |value| Ok(function!(head, value)));
                    result_heads.insert(head, position);
                    roots.push(construction.function(root, head));
                }
            }
        }
        if !roots.is_empty() {
            let root = construction.sum(roots)?;
            let mut network = construction.finish(root);
            network.execute::<Sequential, SmallestDegree, _, _, _>(&library, &functions)?;
            let result = Atom::from(network.result_scalar()?);
            let terms = match result.as_view() {
                AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
                value => vec![value],
            };
            for term in terms {
                let AtomView::Fun(function) = term else {
                    return Err(eyre::eyre!("component result lost its argument boundary").into());
                };
                let Some(position) = result_heads.remove(&function.get_symbol()) else {
                    return Err(eyre::eyre!("unrecognized component result boundary").into());
                };
                let mut arguments = function.iter();
                let Some(value) = arguments.next() else {
                    return Err(eyre::eyre!("empty component result boundary").into());
                };
                if arguments.next().is_some() {
                    return Err(eyre::eyre!("component result has multiple values").into());
                }
                results.insert(position, Ok(value.to_owned()));
            }
            if !result_heads.is_empty() {
                return Err(eyre::eyre!("component execution lost an independent result").into());
            }
        }
        Ok((0..values.len())
            .map(|position| {
                results
                    .remove(&position)
                    .expect("every input has an admission or component result")
            })
            .collect())
    }

    pub fn simple_execute(&mut self) {
        let lib = DummyLibrary::<_>::new();

        self.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &ErroringLibrary::new())
            .unwrap();
    }
}
pub trait NetworkParse {
    /// Parses this symbolic expression into a parameterized tensor network using `settings`.
    #[allow(clippy::result_large_err)]
    fn parse_to_atom_net<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<ParamNet<Aind>, TensorNetworkError<DummyKey, Symbol>>;
}

impl NetworkParse for Atom {
    fn parse_to_atom_net<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<ParamNet<Aind>, TensorNetworkError<DummyKey, Symbol>> {
        self.as_view().parse_to_atom_net::<Aind>(settings)
    }
}

impl NetworkParse for AtomView<'_> {
    fn parse_to_atom_net<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<ParamNet<Aind>, TensorNetworkError<DummyKey, Symbol>> {
        let lib = DummyLibrary::<ParamTensor<ShadowedStructure<Aind>>>::new();

        ParamNet::<Aind>::try_from_view::<ShadowedStructure<Aind>, _>(*self, &lib, settings)
    }
}

#[cfg(test)]
mod scalar_batch_tests {
    use super::*;
    use crate::{network::library::symbolic::ETS, structure::representation};
    use std::sync::{Arc, Mutex};
    use symbolica::{function, parse_lit};

    #[test]
    fn scalar_batch_preserves_independent_admission_and_zero_results() {
        representation::initialize();
        let p = crate::tensor_symbol!("scalar_batch_test::p");
        let q = crate::tensor_symbol!("scalar_batch_test::q");
        let dot = |rep: Atom| function!(ETS.metric, function!(p, &rep), function!(q, &rep));
        let finite = dot(parse_lit!(spenso::mink(2)));
        let symbolic = dot(parse_lit!(spenso::mink(scalar_batch_test::D)));
        let zero = Atom::Zero;
        let values = [&finite, &symbolic, &zero];
        let inputs = values
            .iter()
            .map(|value| value.as_view())
            .collect::<Vec<_>>();
        let settings = ParseSettings::default();
        let mut actual =
            ParamNet::<AbstractIndex>::evaluate_scalar_batch(&inputs, &settings).unwrap();
        let mut reference = finite
            .parse_to_atom_net::<AbstractIndex>(&settings)
            .unwrap();
        reference.simple_execute();
        assert_eq!(
            actual.remove(0).unwrap(),
            Atom::from(reference.result_scalar().unwrap())
        );
        assert!(actual.remove(0).is_err());
        assert_eq!(actual.remove(0).unwrap(), Atom::Zero);
        assert!(
            ParamNet::<AbstractIndex>::evaluate_scalar_batch(&[], &settings)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn scalar_batch_keeps_callback_created_contractions_and_fresh_input_scopes() {
        representation::initialize();
        let matrix = crate::tensor_symbol!("scalar_batch_callback::M");
        let vector = crate::tensor_symbol!("scalar_batch_callback::V");
        let other = crate::tensor_symbol!("scalar_batch_callback::B");
        let calls = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&calls);
        let p = crate::tensor_symbol!(
            "scalar_batch_callback::P",
            norm = move |view, output| {
                let AtomView::Fun(function) = view else {
                    return;
                };
                let argument = function.iter().next().unwrap();
                if SlotMatcher::default()
                    .parse::<LibraryRep, AbstractIndex>(argument)
                    .is_ok()
                {
                    observed.lock().unwrap().push(argument.to_owned());
                    let internal = parse_lit!(spenso::mink(2, 77));
                    **output =
                        function!(matrix, argument, &internal) * function!(vector, &internal);
                }
            }
        );
        let rep = parse_lit!(spenso::mink(2));
        let first = function!(ETS.metric, function!(p, &rep), function!(other, 1, &rep));
        let second = function!(ETS.metric, function!(p, &rep), function!(other, 2, &rep));
        let inputs = [first.as_view(), second.as_view()];
        let settings = ParseSettings::default();
        let reference = inputs
            .iter()
            .map(|value| {
                let mut net = value.parse_to_atom_net::<AbstractIndex>(&settings).unwrap();
                net.simple_execute();
                Atom::from(net.result_scalar().unwrap())
            })
            .collect::<Vec<_>>();
        let reference_calls = std::mem::take(&mut *calls.lock().unwrap());
        let actual = ParamNet::<AbstractIndex>::evaluate_scalar_batch(&inputs, &settings)
            .unwrap()
            .into_iter()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(actual, reference);
        assert_eq!(*calls.lock().unwrap(), reference_calls);
        assert_eq!(reference_calls.len(), 2);
        assert_eq!(reference_calls[0], reference_calls[1]);
        assert!(
            actual
                .iter()
                .all(|value| !value.to_string().contains("mink"))
        );
    }

    #[test]
    fn scalar_batch_result_keys_cannot_reuse_user_callbacks() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        static CALLS: AtomicUsize = AtomicUsize::new(0);
        let _occupied = symbolica::symbol!(
            "spenso::internal_component_result_0_0",
            norm = |_value, output| {
                CALLS.fetch_add(1, Ordering::Relaxed);
                **output = Atom::Zero;
            }
        );
        let first = Atom::num(3);
        let second = Atom::num(5);
        let actual = ParamNet::<AbstractIndex>::evaluate_scalar_batch(
            &[first.as_view(), second.as_view()],
            &ParseSettings::default(),
        )
        .unwrap()
        .into_iter()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
        assert_eq!(actual, vec![first, second]);
        assert_eq!(CALLS.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn scalar_singletons_preserve_custom_metric_execution_order() {
        static EVENTS: Mutex<Vec<String>> = Mutex::new(Vec::new());
        fn sign(index: usize) -> bool {
            EVENTS.lock().unwrap().push(format!("metric:{index}"));
            index == 1
        }
        representation::initialize();
        let rep = representation::REPS
            .write()
            .unwrap()
            .new_inline_metric("scalar_batch_custom::metric", sign)
            .unwrap();
        let p = crate::tensor_symbol!(
            "scalar_batch_custom::p",
            norm = |value, _output| {
                EVENTS.lock().unwrap().push(format!("vector:{value}"));
            }
        );
        let q = crate::tensor_symbol!("scalar_batch_custom::q");
        let compact = function!(rep.symbol(), 2);
        let inputs = [1, 2].map(|label| {
            function!(
                ETS.metric,
                function!(p, &compact),
                function!(q, label, &compact)
            )
        });
        let settings = ParseSettings::default();
        EVENTS.lock().unwrap().clear();
        let expected = inputs
            .iter()
            .map(|value| {
                let mut network = value.parse_to_atom_net::<AbstractIndex>(&settings).unwrap();
                network.simple_execute();
                Atom::from(network.result_scalar().unwrap())
            })
            .collect::<Vec<_>>();
        let expected_events = std::mem::take(&mut *EVENTS.lock().unwrap());
        let actual = inputs
            .iter()
            .map(|value| {
                ParamNet::<AbstractIndex>::evaluate_scalar_batch(&[value.as_view()], &settings)
                    .unwrap()
                    .pop()
                    .unwrap()
                    .unwrap()
            })
            .collect::<Vec<_>>();
        assert_eq!(actual, expected);
        assert_eq!(*EVENTS.lock().unwrap(), expected_events);
        assert!(
            expected_events
                .iter()
                .any(|event| event.starts_with("metric:"))
        );
        assert!(
            expected_events
                .iter()
                .any(|event| event.starts_with("vector:"))
        );
    }

    #[test]
    fn scalar_batch_does_not_enable_unknown_broadcast_execution() {
        representation::initialize();
        let p = crate::tensor_symbol!("scalar_batch_strict::p");
        let q = crate::tensor_symbol!("scalar_batch_strict::q");
        let broadcast = crate::broadcast_symbol!("scalar_batch_strict::f");
        let rep = parse_lit!(spenso::mink(2));
        let dot = function!(ETS.metric, function!(p, &rep), function!(q, &rep));
        let wrapped = function!(broadcast, &dot);
        let settings = ParseSettings::default();
        for inputs in [
            vec![wrapped.as_view()],
            vec![dot.as_view(), wrapped.as_view()],
        ] {
            assert!(ParamNet::<AbstractIndex>::evaluate_scalar_batch(&inputs, &settings).is_err());
        }
    }
}

#[cfg(test)]
mod test;
