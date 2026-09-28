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

mod indices;
pub(crate) mod structure_inference;
pub use structure_inference::{
    AtomStructureExt, ChainNestingError, StructureFromAtom, StructureInferenceMode,
};
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
    Opaque { inference: StructureInferenceMode },
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
            Self::Opaque { .. } => None,
        }
    }

    pub(super) fn with_schoonschip_expansion(self, schoonschip: SchoonschipExpansionMode) -> Self {
        match self {
            Self::Expand { trace, chain, .. } => Self::Expand {
                schoonschip,
                trace,
                chain,
            },
            Self::Opaque { inference } => Self::Opaque { inference },
        }
    }

    fn opaque_inference(self) -> Option<StructureInferenceMode> {
        match self {
            Self::Expand { .. } => None,
            Self::Opaque { inference } => Some(inference),
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

    /// Stop recursive parsing once the current parse depth reaches this value.
    ///
    /// At the limit, the current expression is handed to the opaque tensor
    /// expression boundary as a leaf. `None` means there is no depth limit.
    pub depth_limit: Option<usize>,

    /// Selects what contributes to `depth_limit`.
    ///
    /// When true, only product nesting increments parse depth. When false,
    /// additions and powers also increment depth before their children are
    /// parsed.
    pub depth_is_product_depth: bool,

    /// Selects how shorthand notation is represented in the parsed network.
    ///
    /// `Expand` turns the selected shorthand families into explicit graph
    /// structure. Fresh dummies created by this lowering are local to the
    /// expansion.
    ///
    /// `Opaque` keeps shorthand as a leaf tensor or scalar. Its inference mode
    /// controls whether the exposed structure comes from a fast syntactic walk
    /// or from expanding the shorthand and reading the resulting dangling slots.
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
            depth_is_product_depth: true,
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

#[derive(Clone, Debug)]
pub struct ParseState<Aind = AbstractIndex> {
    depth: usize,
    next_dummy: Rc<Cell<usize>>,
    reserved_indices: Rc<RefCell<HashSet<Atom>>>,
    _aind: PhantomData<fn() -> Aind>,
}

#[allow(clippy::derivable_impls)]
impl<Aind> Default for ParseState<Aind> {
    fn default() -> Self {
        Self {
            depth: 0,
            next_dummy: Rc::new(Cell::new(1_000_000)),
            reserved_indices: Rc::default(),
            _aind: PhantomData,
        }
    }
}

impl<Aind: DummyAind + ParseableAind> ParseState<Aind> {
    /// Reserve written index names once at an operation's input boundary.
    /// Cloned parser states share both these names and subsequent allocations.
    pub fn reserve_indices(&self, value: AtomView<'_>) {
        let mut matcher = SlotMatcher::default();
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
        value.validate_chain_like_nesting()?;
        let state = ParseState::<Aind>::default();
        // Reserve names before any shorthand allocates a dummy.
        state.reserve_indices(value);
        Self::try_from_view_impl(
            value,
            state,
            library,
            function_library,
            settings,
            AtomOrView::View,
        )
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

    fn try_from_view_impl<'node, S, Lib, FunLib>(
        value: AtomView<'node>,
        state: ParseState<Aind>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        profile::bump(Counter::ParseView, 1);
        match value {
            AtomView::Mul(m) => {
                profile::bump(Counter::ParseMul, 1);
                Self::try_from_mul(m, state, library, function_library, settings, retain)
            }
            AtomView::Fun(f) => {
                profile::bump(Counter::ParseFun, 1);
                Self::try_from_fun(f, state, library, function_library, settings, retain)
            }
            AtomView::Add(a) => {
                profile::bump(Counter::ParseAdd, 1);
                Self::try_from_add(a, state, library, function_library, settings, retain)
            }
            AtomView::Pow(p) => {
                profile::bump(Counter::ParsePow, 1);
                Self::try_from_pow(p, state, library, function_library, settings, retain)
            }
            a => Ok(Network::from_scalar(Self::scalar_from_expression(retain(
                a,
            ))?)),
        }
    }

    #[allow(clippy::type_complexity, clippy::result_large_err)]
    fn as_leaf<'node, S, Lib, FunLib>(
        value: AtomView<'node>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
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
        if !value.is_tensorial(settings.strict_tensor_filter) && !scalar_composite {
            return Ok(Self::from_scalar(Self::scalar_from_expression(retain(
                value,
            ))?));
        }

        Self::as_inferred_leaf::<S, Lib, FunLib>(
            value,
            StructureInferenceMode::Fast,
            library,
            function_library,
            settings,
            retain,
        )
    }

    #[allow(clippy::type_complexity, clippy::result_large_err)]
    fn as_inferred_leaf<'node, S, Lib, FunLib>(
        value: AtomView<'node>,
        mode: StructureInferenceMode,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        profile::bump(Counter::ParseStructureAttempt, 1);
        let structure = {
            let _span = profile::span(Timer::ParseStructure);
            S::structure_from_atom(value, mode)
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
        let mut network = Self::from_tensor(T::tensor_from_expression(
            retain(value),
            structure,
            library,
            function_library,
            settings,
        )?);
        let node = network.graph.graph.node_id(network.graph.head());
        network.graph.set_logical_layout(node, &layout)?;
        Ok(network)
    }

    #[allow(clippy::result_large_err)]
    fn try_from_mul<'node, S, Lib, FunLib>(
        value: MulView<'node>,
        mut state: ParseState<Aind>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
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
                value.as_view(),
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
            first_atom,
            state.clone(),
            library,
            function_library,
            settings,
            retain,
        )?;

        // state

        if settings.precontract_scalars {
            let mut scalars = Atom::num(1);

            let rest: Result<Vec<_>, _> = iter
                .filter_map(|a| {
                    match Self::try_from_view_impl(
                        a,
                        state.clone(),
                        library,
                        function_library,
                        settings,
                        retain,
                    ) {
                        Ok(n) => {
                            profile::bump(Counter::MulFactor, 1);
                            if let NetworkState::PureScalar = n.state {
                                let _span = profile::span(Timer::ScalarMulAccum);
                                profile::bump(Counter::ScalarMulAccum, 1);
                                scalars *= a;
                                None
                            } else {
                                Some(Ok((a.to_owned(), n)))
                            }
                        }
                        Err(e) => Some(Err(e)),
                    }
                })
                .collect();

            let mut res = rest?;

            if let NetworkState::PureScalar = first.state {
                let _span = profile::span(Timer::ScalarMulAccum);
                profile::bump(Counter::ScalarMulAccum, 1);
                scalars *= first_atom;
            } else {
                res.push((first_atom.to_owned(), first));
            }

            if res.is_empty() {
                Ok(Self::from_scalar(Self::scalar_from_expression(retain(
                    value.as_view(),
                ))?))
            } else {
                let s = if scalars != Atom::num(1) {
                    Self::from_scalar(Self::scalar_from_expression(AtomOrView::Atom(scalars))?)
                } else {
                    res.pop().unwrap().1
                };

                Ok(s.n_mul(res.into_iter().map(|(_, net)| net)))
            }
        } else {
            let rest: Result<Vec<_>, _> = iter
                .map(|a| {
                    profile::bump(Counter::MulFactor, 1);
                    Self::try_from_view_impl(
                        a,
                        state.clone(),
                        library,
                        function_library,
                        settings,
                        retain,
                    )
                })
                .collect();

            Ok(first.n_mul(rest?))
        }
    }

    #[allow(clippy::result_large_err)]
    fn try_from_fun<'node, S, Lib, FunLib>(
        value: FunView<'node>,
        state: ParseState<Aind>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
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

        if symbol == SPENSO_TAG.dot && value.get_nargs() != 2 {
            return Err(TensorNetworkError::InvalidDotFunction(
                value.as_view().to_plain_string(),
            ));
        }

        if symbol == SPENSO_TAG.bracket || symbol.has_tag(&SPENSO_TAG.broadcast) {
            return Self::parse_expanded_function::<S, Lib, FunLib>(
                value,
                state,
                library,
                function_library,
                settings,
                retain,
            );
        }

        if !value.as_view().is_tensorial(settings.strict_tensor_filter) {
            return Self::parse_scalar_function(value, retain);
        }

        if let Some(inference) = settings.shorthand_parsing.opaque_inference()
            && Self::is_shorthand_function(value)
        {
            return Self::as_inferred_leaf::<S, _, _>(
                value.as_view(),
                inference,
                library,
                function_library,
                settings,
                retain,
            );
        }

        Self::parse_expanded_function(value, state, library, function_library, settings, retain)
    }

    #[allow(clippy::result_large_err)]
    fn parse_scalar_function<'node>(
        value: FunView<'node>,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<Self, TensorNetworkError<K, Symbol>> {
        if value.get_symbol() == SPENSO_TAG.pure_scalar {
            if value.get_nargs() != 1 {
                return Err(TensorNetworkError::TooManyArgsFunction(
                    value.as_view().to_plain_string(),
                ));
            }

            return Ok(Self::from_scalar(Self::scalar_from_expression(retain(
                value.iter().next().unwrap(),
            ))?));
        }

        Ok(Self::from_scalar(Self::scalar_from_expression(retain(
            value.as_view(),
        ))?))
    }

    #[allow(clippy::result_large_err)]
    fn parse_broadcast_function<'node, S, Lib, FunLib>(
        value: FunView<'node>,
        state: ParseState<Aind>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
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
        let inner = value.iter().next().unwrap();
        let inner_tensor =
            Self::try_from_view_impl(inner, state, library, function_library, settings, retain)?;

        Ok(inner_tensor.fun(symbol))
    }

    #[allow(clippy::result_large_err)]
    fn parse_expanded_function<'node, S, Lib, FunLib>(
        value: FunView<'node>,
        state: ParseState<Aind>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        let symbol = value.get_symbol();

        if symbol == SPENSO_TAG.bracket {
            let mut n_muls = value
                .iter()
                .map(|a| {
                    Self::try_from_view_impl(
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
            Ok(first.n_mul(n_muls))
        } else if symbol.has_tag(&SPENSO_TAG.broadcast) {
            Self::parse_broadcast_function::<S, Lib, FunLib>(
                value,
                state,
                library,
                function_library,
                settings,
                retain,
            )
        } else {
            Self::materialize_shorthand::<S, Lib, FunLib>(
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
        value: FunView<'node>,
        library: &Lib,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
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
                return Ok(Self::from_scalar(Self::scalar_from_expression(retain(
                    value.as_view(),
                ))?));
            }
            Err(err) => {
                profile::bump(Counter::ParseStructureErr, 1);
                return Err(err.into());
            }
        };

        match library.key_for_structure(&structure) {
            Ok(key) => {
                let tensor_structure = structure.canonical().clone();
                Ok(Self::library_tensor(
                    &tensor_structure,
                    structure.map_canonical(|_| key),
                ))
            }
            Err(_) if structure.canonical().is_scalar() => Ok(Self::from_scalar(
                Self::scalar_from_expression(retain(value.as_view()))?,
            )),
            Err(_) => {
                // Tensor dimensions may remain symbolic during analytic normalization;
                // opaque scalars and library tensors need no eager shadow here.
                // The target decides whether missing leaves need finite components.
                let layout = structure.layout().clone();
                let mut network =
                    Self::from_tensor(T::tensor_from_leaf(retain(value.as_view()), structure)?);
                let node = network.graph.graph.node_id(network.graph.head());
                network.graph.set_logical_layout(node, &layout)?;
                Ok(network)
            }
        }
    }

    #[allow(clippy::result_large_err)]
    fn try_from_pow<'node, S, Lib, FunLib>(
        value: PowView<'node>,
        mut state: ParseState<Aind>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> std::result::Result<Self, TensorNetworkError<K, Symbol>>
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
                value.as_view(),
                library,
                function_library,
                settings,
                retain,
            );
        }

        if !settings.depth_is_product_depth {
            state.depth += 1;
        }
        let (base_expression, exp) = value.get_base_exp();

        if let Ok(n) = i8::try_from(exp) {
            // println!("base:{base_expression}");
            let next_dummy = state.next_dummy.get();
            let base = Self::try_from_view_impl(
                base_expression,
                state.clone(),
                library,
                function_library,
                settings,
                retain,
            )?;

            // println!("base state {:?}", base.state);
            if settings.precontract_scalars
                && let NetworkState::PureScalar = base.state
                && !matches!(base_expression, AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.bracket)
            {
                // println!("Pure");
                return Ok(Self::from_scalar(Self::scalar_from_expression(retain(
                    value.as_view(),
                ))?));
            }

            if let NetworkState::Tensor = base.state {
                return Err(TensorNetworkError::NonSelfDualTensorPower(
                    value.as_view().to_plain_string(),
                ));
            }
            // An even power of a self_dual tensor, or scalar is a scalar
            if n < 0 && n % 2 != 0 && !base.state.is_scalar() {
                return Err(TensorNetworkError::NegativeExponentNonScalar(format!(
                    "Atom:{},graph of base: {}, dangling indices: {:?}",
                    value.as_view().to_plain_string(),
                    base.dot(),
                    base.graph.dangling_indices()
                )));
            }

            let out = if n > 1
                && (state.next_dummy.get() != next_dummy
                    || (n == 2 && base.state == NetworkState::SelfDualTensor))
            {
                // Each lowered shorthand copy needs independent internal indices.
                // A tensor square also needs the selected product-contraction
                // strategy: scalar power execution uses the tensor's default
                // contraction, which can leave symbolic sums uncontracted.
                let rest = (1..n)
                    .map(|_| {
                        Self::try_from_view_impl(
                            base_expression,
                            state.clone(),
                            library,
                            function_library,
                            settings,
                            retain,
                        )
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                base.n_mul(rest)
            } else {
                // println!("{:?}", base.state);
                base.pow(n)
            };
            // println!("{:?}", out.state);
            Ok(out)
        } else {
            Ok(Self::from_scalar(Self::scalar_from_expression(retain(
                value.as_view(),
            ))?))
        }
    }

    #[allow(clippy::result_large_err)]
    fn try_from_add<'node, S, Lib, FunLib>(
        value: AddView<'node>,
        mut state: ParseState<Aind>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
        retain: fn(AtomView<'node>) -> AtomOrView<'src>,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
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
                value.as_view(),
                library,
                function_library,
                settings,
                retain,
            );
        }

        if !settings.depth_is_product_depth {
            state.depth += 1;
        }

        let mut iter = value.iter();

        let first_atom = iter.next().unwrap();
        profile::bump(Counter::AddTerm, 1);

        let first = Self::try_from_view_impl(
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
            let mut scalars = Atom::Zero;

            let rest: Result<Vec<_>, _> = iter
                .filter_map(|a| {
                    match Self::try_from_view_impl(
                        a,
                        state.clone(),
                        library,
                        function_library,
                        settings,
                        retain,
                    ) {
                        Ok(n) => {
                            profile::bump(Counter::AddTerm, 1);
                            if n.state.is_compatible(&first.state) {
                                if let NetworkState::PureScalar = n.state {
                                    let _span = profile::span(Timer::ScalarAddAccum);
                                    profile::bump(Counter::ScalarAddAccum, 1);
                                    scalars += a;
                                    None
                                } else {
                                    Some(Ok(n))
                                }
                            } else {
                                Some(Err(TensorNetworkError::IncompatibleSummand(format!(
                                    "{} is {:?} vs {} is {:?}",
                                    a, n.state, first_atom, first.state
                                ))))
                            }
                        }
                        Err(e) => Some(Err(e)),
                    }
                })
                .collect();

            let mut res = rest?;

            if let NetworkState::PureScalar = first.state {
                let _span = profile::span(Timer::ScalarAddAccum);
                profile::bump(Counter::ScalarAddAccum, 1);
                scalars += first_atom;
            } else {
                res.push(first);
            }

            if res.is_empty() {
                Ok(Self::from_scalar(Self::scalar_from_expression(retain(
                    value.as_view(),
                ))?))
            } else {
                let s = if scalars != Atom::Zero {
                    Self::from_scalar(Self::scalar_from_expression(AtomOrView::Atom(scalars))?)
                } else {
                    res.pop().unwrap()
                };
                Ok(s.n_add(res))
            }
        } else {
            let rest: Result<Vec<_>, _> = iter
                .map(|a| {
                    match Self::try_from_view_impl(
                        a,
                        state.clone(),
                        library,
                        function_library,
                        settings,
                        retain,
                    ) {
                        Ok(n) => {
                            profile::bump(Counter::AddTerm, 1);
                            if n.state.is_compatible(&first.state) {
                                Ok(n)
                            } else {
                                Err(TensorNetworkError::IncompatibleSummand(format!(
                                    "{} is {:?} vs {} is {:?}",
                                    a, n.state, first_atom, first.state
                                )))
                            }
                        }
                        Err(e) => Err(e),
                    }
                })
                .collect();

            Ok(first.n_add(rest?))
        }
    }
}
pub type ParamNet<Aind> =
    Network<NetworkStore<ParamTensor<ShadowedStructure<Aind>>, Atom>, DummyKey, Symbol, Aind>;

impl<Aind: AbsInd + DummyAind + ParseableAind + 'static> ParamNet<Aind> {
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
mod test;
