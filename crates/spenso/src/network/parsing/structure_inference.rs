//! Structure inference for symbolic parser leaves.
//!
//! The fast path is intentionally syntactic:
//! 1. dispatch known shorthand (`chain`, `trace`) to its visible-slot convention;
//! 2. otherwise infer ordinary tensor syntax from sums, products, powers, and functions;
//! 3. wrap the ordered slots in the requested structure type;
//!
//! Tests compare this syntax result to the existing fully expanded graph parser.

use symbolica::{
    atom::{
        Atom, AtomCore, AtomView, FunctionBuilder, MulView, PowView, Symbol,
        representation::FunView,
    },
    domains::rational::Rational,
};
use thiserror::Error;

use super::{ParseState, ShadowedStructure, StrictTensorFilter};
use crate::structure::{
    Canonicalized, HasName, NamedStructure, OrderedStructure, StructureContract, StructureError,
    TensorStructure,
    representation::{BaseRepName, LibraryRep, Minkowski},
    slot::{AbsInd, DummyAind, ParseableAind, Slot, SlotError, SlotMatch, SlotMatcher},
};
use crate::{
    network::tags::{SPENSO_TAG, SpensoTags},
    shadowing,
};

/// A chain or trace found inside another chain or trace.
///
/// Both shorthands use the same global `in` and `out` symbols, so treating a
/// nested shorthand as an independent placeholder scope would be ambiguous.
#[derive(Clone, Copy, Debug, Eq, Error, PartialEq)]
#[error(
    "`{inner}` cannot be nested inside `{outer}` because chain and trace share one global `in`/`out` placeholder scope"
)]
pub struct ChainNestingError {
    outer: Symbol,
    inner: Symbol,
}

impl ChainNestingError {
    /// Enter one function during an existing syntax walk. Scalar metadata is
    /// opaque to placeholder scopes; callers do not descend into it.
    pub fn enter(
        owner: Option<Symbol>,
        symbol: Symbol,
        tags: &SpensoTags,
    ) -> Result<Option<Symbol>, Self> {
        if symbol == tags.chain || symbol == tags.trace {
            if let Some(outer) = owner {
                return Err(Self {
                    outer,
                    inner: symbol,
                });
            }
            Ok(Some(symbol))
        } else {
            Ok(owner)
        }
    }
}

pub trait StructureFromAtom: Sized {
    /// Infer the permuted tensor structure exposed by `value`.
    ///
    /// This syntactic observation does not materialize dummies or build a graph.
    /// Implementations reject nested chain/trace placeholder scopes; delegating
    /// wrappers must not repeat that whole-expression validation. The caller
    /// supplies operation-local recognition state; implementations must not
    /// retain source views or a borrow of that state during materialization.
    fn structure_from_atom(
        value: AtomView<'_>,
        matcher: &mut SlotMatcher,
    ) -> Result<Canonicalized<Self>, StructureError>;

    /// Infer the parser's current subtree. Custom implementations retain the
    /// checked entry by default; the parser alone binds the source view.
    fn structure_from_parser<Aind>(
        state: &ParseState<Aind, AtomView<'_>>,
    ) -> Result<Canonicalized<Self>, StructureError> {
        Self::structure_from_atom(state.current_view(), &mut state.matcher.borrow_mut())
    }

    /// Infer structure with operation-local recognition state.
    fn parse(value: AtomView<'_>) -> Result<Canonicalized<Self>, StructureError> {
        Self::structure_from_atom(value, &mut SlotMatcher::default())
    }
}

pub trait AtomStructureExt {
    /// Whether an explicit index payload repeats in a multiplicative scope.
    ///
    /// Scans arbitrary function heads in the normalized Atom tree without
    /// constructing tensor structures; known slot payloads remain opaque.
    /// Summands are alternatives; integer powers with magnitude greater than
    /// one repeat the base's indices. Compact representations carry no explicit
    /// index. Equality ignores dimension and duality: this is a syntactic
    /// candidate check, not a certificate of a valid contraction.
    /// Stops at the first candidate.
    fn has_repeated_explicit_indices(&self) -> bool;

    /// Check repeated indices while observing the same syntactic traversal.
    ///
    /// The observer receives each visited node before its children, together
    /// with its slot classification (`Other` for non-functions). Explicit,
    /// compact and malformed slots are visited once; their dimensions and index
    /// payloads remain opaque. The root is always visited, including a root sum.
    ///
    /// At the first repetition, `on_first_repeat` is called exactly once. Return
    /// true to observe the remainder without collecting further index occurrences,
    /// or false to stop. With no repetition the decision is not called and the
    /// traversal finishes. Slot payloads remain opaque in either case.
    fn has_repeated_explicit_indices_with_observer<'a>(
        &'a self,
        observe: impl FnMut(AtomView<'a>, &SlotMatch<'a>),
        on_first_repeat: impl FnOnce() -> bool,
    ) -> bool;

    /// Replace four-dimensional Minkowski representations and slots by `dimension`.
    /// This must precede Lorentz contractions; scalar factors and other spaces are unchanged.
    fn with_lorentz_dimension(&self, dimension: AtomView<'_>) -> Atom;

    /// Convenience wrapper for `StructureFromAtom::structure_from_atom`.
    fn infer_structure<S: StructureFromAtom>(&self) -> Result<Canonicalized<S>, StructureError>;

    /// Reject more than one chain/trace placeholder consumer on any expression path.
    fn validate_chain_like_nesting(&self) -> Result<(), ChainNestingError>;

    /// Return true when this expression is valid tensor parser syntax at its root.
    ///
    /// Ordinary tensor heads expose direct structural arguments: slots,
    /// `aind(...)` bundles, or compact representation arguments. Shorthand roots
    /// (`chain`, `trace`, `dot`, and compact metrics) are tensorial because the
    /// parser gives them explicit semantics. Broadcast functions are tensorial
    /// only when they have one tensorial argument. Untagged wrappers around tensor
    /// expressions stay scalar because their head has no tensor semantics.
    fn is_tensorial(&self, filter: StrictTensorFilter) -> bool;
}

impl AtomStructureExt for Atom {
    fn has_repeated_explicit_indices(&self) -> bool {
        self.as_view().has_repeated_explicit_indices()
    }

    fn has_repeated_explicit_indices_with_observer<'a>(
        &'a self,
        observe: impl FnMut(AtomView<'a>, &SlotMatch<'a>),
        on_first_repeat: impl FnOnce() -> bool,
    ) -> bool {
        super::indices::RepeatedIndices::contains(self.as_view(), observe, on_first_repeat)
    }

    fn with_lorentz_dimension(&self, dimension: AtomView<'_>) -> Atom {
        self.as_view().with_lorentz_dimension(dimension)
    }

    fn infer_structure<S: StructureFromAtom>(&self) -> Result<Canonicalized<S>, StructureError> {
        self.as_view().infer_structure()
    }

    fn validate_chain_like_nesting(&self) -> Result<(), ChainNestingError> {
        self.as_view().validate_chain_like_nesting()
    }

    fn is_tensorial(&self, filter: StrictTensorFilter) -> bool {
        self.as_view().is_tensorial(filter)
    }
}

impl AtomStructureExt for AtomView<'_> {
    fn has_repeated_explicit_indices(&self) -> bool {
        super::indices::RepeatedIndices::contains(*self, |_, _| {}, || false)
    }

    fn has_repeated_explicit_indices_with_observer<'a>(
        &'a self,
        observe: impl FnMut(AtomView<'a>, &SlotMatch<'a>),
        on_first_repeat: impl FnOnce() -> bool,
    ) -> bool {
        super::indices::RepeatedIndices::contains(*self, observe, on_first_repeat)
    }

    fn with_lorentz_dimension(&self, dimension: AtomView<'_>) -> Atom {
        self.replace_map(|value, _, output| {
            if let AtomView::Fun(fun) = value
                && fun.get_symbol() == Minkowski::selfless_symbol()
                && matches!(fun.get_nargs(), 1 | 2)
                && fun.iter().next() == Some(Atom::num(4).as_view())
            {
                **output = fun
                    .iter()
                    .skip(1)
                    .fold(
                        FunctionBuilder::new(fun.get_symbol()).add_arg(dimension),
                        |builder, arg| builder.add_arg(arg),
                    )
                    .finish();
            }
        })
    }

    fn infer_structure<S: StructureFromAtom>(&self) -> Result<Canonicalized<S>, StructureError> {
        // StructureFromAtom owns validation, including callers that bypass this wrapper.
        S::structure_from_atom(*self, &mut SlotMatcher::default())
    }

    fn validate_chain_like_nesting(&self) -> Result<(), ChainNestingError> {
        TensorialSyntax::validate_chain_like_nesting(*self, None, &SPENSO_TAG)
    }

    fn is_tensorial(&self, filter: StrictTensorFilter) -> bool {
        TensorialSyntax::is_tensorial(*self, filter, &SlotMatcher::default())
    }
}

pub(crate) struct TensorialSyntax;

impl TensorialSyntax {
    pub(crate) fn is_tensorial(
        value: AtomView<'_>,
        filter: StrictTensorFilter,
        matcher: &SlotMatcher,
    ) -> bool {
        match value {
            AtomView::Fun(fun) => Self::function_is_tensorial(fun, filter, matcher),
            AtomView::Var(var) => var.get_symbol().has_attributes_of(matcher.tags().rep_),
            AtomView::Add(add) => add
                .iter()
                .any(|arg| Self::is_tensorial(arg, filter, matcher)),
            AtomView::Mul(mul) => mul
                .iter()
                .any(|arg| Self::is_tensorial(arg, filter, matcher)),
            AtomView::Pow(pow) => Self::is_tensorial(pow.get_base_exp().0, filter, matcher),
            _ => false,
        }
    }

    fn validate_chain_like_nesting(
        value: AtomView<'_>,
        owner: Option<Symbol>,
        tags: &SpensoTags,
    ) -> Result<(), ChainNestingError> {
        match value {
            AtomView::Add(sum) => {
                for term in sum.iter() {
                    Self::validate_chain_like_nesting(term, owner, tags)?;
                }
            }
            AtomView::Mul(product) => {
                for factor in product.iter() {
                    Self::validate_chain_like_nesting(factor, owner, tags)?;
                }
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                Self::validate_chain_like_nesting(base, owner, tags)?;
                Self::validate_chain_like_nesting(exponent, owner, tags)?;
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                if symbol.is_scalar() {
                    return Ok(());
                }
                let owner = ChainNestingError::enter(owner, symbol, tags)?;
                for argument in function.iter() {
                    Self::validate_chain_like_nesting(argument, owner, tags)?;
                }
            }
            _ => {}
        }
        Ok(())
    }

    pub(crate) fn function_is_tensorial(
        fun: FunView<'_>,
        filter: StrictTensorFilter,
        matcher: &SlotMatcher,
    ) -> bool {
        let symbol = fun.get_symbol();
        if symbol.is_scalar() {
            return false;
        }
        // Public bundle access probes Symbolica initialization. Reuse the
        // handles throughout this predicate instead of probing for every field.
        let tags = matcher.tags();
        if symbol == tags.pure_scalar || symbol == tags.scalar {
            return false;
        }

        // Registered projectors preserve tensor structure through arithmetic
        // prefactors just like brackets. Ordinary function metadata remains opaque.
        if symbol == tags.bracket || matcher.is_projector(symbol) {
            return fun
                .iter()
                .any(|arg| Self::is_tensorial(arg, filter, matcher));
        }

        if symbol.has_attributes_of(tags.rep_)
            || symbol == tags.chain
            || symbol == tags.trace
            || symbol == tags.dot
        {
            return true;
        }

        // Variance wrappers preserve representation syntax. In particular,
        // g(rep(i), dind(rep(j))) must retain its oriented tensor ports.
        if matcher.is_variance_wrapper(symbol)
            && fun.get_nargs() == 1
            && fun.iter().next().is_some_and(|arg| match arg {
                AtomView::Fun(rep) => rep.get_symbol().has_attributes_of(tags.rep_),
                AtomView::Var(rep) => rep.get_symbol().has_attributes_of(tags.rep_),
                _ => false,
            })
        {
            return true;
        }

        if symbol == matcher.metric() {
            return fun.get_nargs() == 2
                && fun
                    .iter()
                    .all(|arg| Self::is_tensorial(arg, filter, matcher));
        }

        if symbol.has_tag(&tags.broadcast) {
            let args = fun.iter().collect::<Vec<_>>();
            return matches!(args.as_slice(), [arg] if Self::is_tensorial(*arg, filter, matcher));
        }

        match filter {
            StrictTensorFilter::Tagged => symbol.has_tag(&tags.tensor),
            StrictTensorFilter::TaggedChecked => {
                symbol.has_tag(&tags.tensor)
                    && (fun.get_nargs() == 0
                        || fun
                            .iter()
                            .any(|arg| Self::contains_representation_syntax(arg, matcher)))
            }
            StrictTensorFilter::ContainsReps => fun
                .iter()
                .any(|arg| Self::contains_representation_syntax(arg, matcher)),
        }
    }

    fn contains_representation_syntax(value: AtomView<'_>, matcher: &SlotMatcher) -> bool {
        match value {
            AtomView::Fun(fun) => {
                !fun.get_symbol().is_scalar()
                    && (fun.get_symbol().has_attributes_of(matcher.tags().rep_)
                        || fun
                            .iter()
                            .any(|arg| Self::contains_representation_syntax(arg, matcher)))
            }
            AtomView::Var(var) => var.get_symbol().has_attributes_of(matcher.tags().rep_),
            AtomView::Add(add) => add
                .iter()
                .any(|arg| Self::contains_representation_syntax(arg, matcher)),
            AtomView::Mul(mul) => mul
                .iter()
                .any(|arg| Self::contains_representation_syntax(arg, matcher)),
            AtomView::Pow(pow) => {
                Self::contains_representation_syntax(pow.get_base_exp().0, matcher)
            }
            _ => false,
        }
    }
}

impl<Aind: AbsInd + DummyAind + ParseableAind> StructureFromAtom
    for OrderedStructure<LibraryRep, Aind>
{
    fn structure_from_atom(
        value: AtomView<'_>,
        matcher: &mut SlotMatcher,
    ) -> Result<Canonicalized<Self>, StructureError> {
        TensorialSyntax::validate_chain_like_nesting(value, None, matcher.tags())
            .map_err(|error| StructureError::ParsingError(error.to_string()))?;
        Self::leaf_structure_from_atom(value, matcher)
    }

    fn structure_from_parser<Index>(
        state: &ParseState<Index, AtomView<'_>>,
    ) -> Result<Canonicalized<Self>, StructureError> {
        let mut matcher = state.matcher.borrow_mut();
        if state.chain_scope_validated {
            Self::leaf_structure_from_atom(state.current_view(), &mut matcher)
        } else {
            Self::structure_from_atom(state.current_view(), &mut matcher)
        }
    }
}

impl<Aind: AbsInd + ParseableAind> OrderedStructure<LibraryRep, Aind> {
    /// Pick the fast leaf convention for the top-level atom.
    ///
    /// `chain(start, end, factors...)` and `trace(rep, factors...)` have
    /// structural arguments that are not the same as generic function slots, so
    /// they are routed to their own rules. Everything else uses ordinary
    /// syntactic tensor parsing.
    fn leaf_structure_from_atom(
        value: AtomView<'_>,
        matcher: &mut SlotMatcher,
    ) -> Result<Canonicalized<Self>, StructureError> {
        match value {
            AtomView::Fun(fun) if fun.get_symbol() == matcher.tags().chain => {
                Self::chain_structure_from_fun(fun, matcher)
            }
            AtomView::Fun(fun) if fun.get_symbol() == matcher.tags().trace => {
                Self::trace_structure_from_fun(fun, matcher)
            }
            _ => Self::from_syntactic_atom(value, matcher),
        }
    }

    /// Infer an `OrderedStructure` from ordinary tensor syntax without shorthand semantics.
    ///
    /// The dispatcher is intentionally shallow: sums, products, powers, and
    /// functions each have their own convention below. Scalar syntax is carried
    /// internally as `OrderedStructure::empty()` and converted to
    /// `EmptyStructure` only at this leaf boundary.
    fn from_syntactic_atom(
        value: AtomView<'_>,
        matcher: &mut SlotMatcher,
    ) -> Result<Canonicalized<Self>, StructureError> {
        let structure = Self::syntactic_structure_from_atom(value, matcher)?;
        if structure.canonical().is_scalar() {
            Err(StructureError::EmptyStructure(SlotError::EmptyStructure))
        } else {
            Ok(structure)
        }
    }

    pub(super) fn syntactic_structure_from_atom(
        value: AtomView<'_>,
        matcher: &mut SlotMatcher,
    ) -> Result<Canonicalized<Self>, StructureError> {
        match value {
            AtomView::Add(add) => {
                let Some(first) = add.iter().next() else {
                    return Ok(Canonicalized::identity(OrderedStructure::empty()));
                };
                Self::syntactic_structure_from_atom(first, matcher)
            }
            AtomView::Pow(pow) => Self::from_power_atom(pow, matcher),
            AtomView::Mul(mul) => Self::from_product_atom(mul, matcher),
            AtomView::Fun(fun) => Self::from_function_atom(fun, matcher),
            _ => Ok(Canonicalized::identity(OrderedStructure::empty())),
        }
    }

    /// Infer an `OrderedStructure` from a power's base structure and exponent.
    ///
    /// Scalars stay scalar. A fully self-dual tensor to an even integer power
    /// has no external structure, while an odd integer power keeps the base
    /// structure. Fractional powers and powers of non-self-dual tensors are
    /// rejected because their external structure is not well-defined here.
    fn from_power_atom(
        pow: PowView<'_>,
        matcher: &mut SlotMatcher,
    ) -> Result<Canonicalized<Self>, StructureError> {
        let (base, exp) = pow.get_base_exp();
        let base_structure = Self::syntactic_structure_from_atom(base, matcher)?;

        if base_structure.canonical().is_scalar() {
            Ok(base_structure)
        } else if base_structure.canonical().is_fully_self_dual()
            && let Ok(r) = Rational::try_from(exp)
        {
            if r.numerator() % 2 == 0 {
                Ok(Canonicalized::identity(OrderedStructure::empty()))
            } else if r.denominator() == 1 {
                Ok(base_structure)
            } else {
                Err(StructureError::ParsingError(format!(
                    "Invalid power of tensor {}",
                    pow.as_view()
                )))
            }
        } else {
            Err(StructureError::ParsingError(format!(
                "Invalid power of tensor {}",
                pow.as_view()
            )))
        }
    }

    /// Infer an `OrderedStructure` from a product by merging every factor that exposes slots.
    ///
    /// Scalar factors are empty structures, so merging them is a no-op. If no
    /// factor exposes a slot, the product remains an empty scalar structure.
    fn from_product_atom(
        product: MulView<'_>,
        matcher: &mut SlotMatcher,
    ) -> Result<Canonicalized<Self>, StructureError> {
        let mut structure = Canonicalized::identity(OrderedStructure::empty());

        for factor in product {
            structure = structure
                .merge_syntactic(&Self::syntactic_structure_from_atom(factor, matcher)?)?;
        }

        Ok(structure)
    }

    /// Infer an `OrderedStructure` from a generic function's direct structural arguments.
    ///
    /// A direct slot argument contributes one exposed slot. An `aind(...)`
    /// bundle is flattened into its slots. Other arguments are treated as
    /// metadata for the eventual named leaf and do not erase slots already seen.
    /// Recognized abstract ports with invalid dimensions or indices are errors;
    /// concrete component markers remain non-structural metadata.
    /// Chain projectors without direct structural arguments expose the combined
    /// slots of their factor sequence.
    fn from_function_atom(
        fun: FunView<'_>,
        matcher: &mut SlotMatcher,
    ) -> Result<Canonicalized<Self>, StructureError> {
        if fun.get_symbol().is_scalar() {
            return Ok(Canonicalized::identity(OrderedStructure::empty()));
        }
        if matcher.compact_inner_product_parts::<Aind>(fun).is_some() {
            let mut arguments = fun.iter();
            let left = Self::syntactic_structure_from_atom(arguments.next().unwrap(), matcher)?;
            let right = Self::syntactic_structure_from_atom(arguments.next().unwrap(), matcher)?;
            return left.merge_syntactic(&right);
        }
        // Nested chain/trace operands retain their established shorthand interface.
        if fun.get_symbol() == matcher.tags().chain || fun.get_symbol() == matcher.tags().trace {
            return match Self::leaf_structure_from_atom(fun.as_view(), matcher) {
                Err(StructureError::EmptyStructure(_)) => {
                    Ok(Canonicalized::identity(Self::empty()))
                }
                result => result,
            };
        }
        if fun.get_symbol() == matcher.index_bundle() {
            let mut slots = Vec::new();
            for arg in fun.iter() {
                slots.push(matcher.parse::<LibraryRep, Aind>(arg)?);
            }
            return Ok(OrderedStructure::new(slots));
        }

        let mut slots = Vec::new();

        for arg in fun.iter() {
            let recognized = matcher.classify(arg);
            match recognized.parse::<LibraryRep, Aind>(arg, matcher) {
                Ok(slot) => {
                    slots.push(slot);
                }
                Err(error) => {
                    // Explicit abstract ports must satisfy the caller's index
                    // grammar. Concrete component markers and non-slot metadata
                    // remain opaque; a failed abstract port cannot become scalar.
                    if matches!(recognized, SlotMatch::Explicit(slot) if !slot.is_concrete_index())
                    {
                        return Err(error.into());
                    }
                    if let AtomView::Fun(fun) = arg
                        && fun.get_symbol() == matcher.index_bundle()
                    {
                        let internal = Self::from_function_atom(fun, matcher)?;
                        slots.extend(
                            internal
                                .layout()
                                .canonical_to_logical(&internal.canonical().structure),
                        );
                    }
                }
            }
        }

        if slots.is_empty() && matcher.projectors().contains(&fun.get_symbol()) {
            let mut structure = Canonicalized::identity(OrderedStructure::empty());
            for factor in fun.iter() {
                structure = structure
                    .merge_syntactic(&Self::syntactic_structure_from_atom(factor, matcher)?)?;
            }
            return Ok(structure);
        }

        Ok(OrderedStructure::new(slots))
    }

    /// Infer an `OrderedStructure` from an opaque open chain.
    ///
    /// `args[0]` and `args[1]` are the external endpoints. Remaining factors
    /// may contain other external slots, so they are scanned recursively. The
    /// symbolic placeholders `in` and `out` are just wiring labels and are not
    /// materialized as dummies in this mode.
    fn chain_structure_from_fun(
        fun: FunView<'_>,
        matcher: &mut SlotMatcher,
    ) -> Result<Canonicalized<Self>, StructureError> {
        let args = fun.iter().collect::<Vec<_>>();
        if args.len() < 2 {
            return Err(StructureError::WrongNumberOfArguments(args.len(), 2));
        }

        let mut slots = vec![
            matcher.parse::<LibraryRep, Aind>(args[0])?,
            matcher.parse::<LibraryRep, Aind>(args[1])?,
        ];
        for factor in &args[2..] {
            Self::append_syntactic_slots(*factor, &mut slots, matcher)?;
        }

        Self::from_slots(slots)
    }

    /// Infer an `OrderedStructure` from an opaque trace.
    ///
    /// `args[0]` is the traced representation, not an exposed slot. The factors
    /// are scanned for any non-placeholder slots that remain external to the
    /// trace shorthand.
    fn trace_structure_from_fun(
        fun: FunView<'_>,
        matcher: &mut SlotMatcher,
    ) -> Result<Canonicalized<Self>, StructureError> {
        let args = fun.iter().collect::<Vec<_>>();
        if args.is_empty() {
            return Err(StructureError::WrongNumberOfArguments(0, 1));
        }

        let mut slots = Vec::new();
        for factor in shadowing::trace_factor_views(&args[1..]) {
            Self::append_syntactic_slots(factor, &mut slots, matcher)?;
        }

        Self::from_slots(slots)
    }

    /// Convert exposed slots into the canonical ordered representation.
    ///
    /// An empty slot list means the leaf is scalar, so this returns
    /// `EmptyStructure` instead of an explicit scalar structure.
    fn from_slots(
        slots: Vec<Slot<LibraryRep, Aind>>,
    ) -> Result<Canonicalized<Self>, StructureError> {
        if slots.is_empty() {
            Err(StructureError::EmptyStructure(SlotError::EmptyStructure))
        } else {
            Ok(OrderedStructure::new(slots))
        }
    }

    /// Append slots from the `OrderedStructure` inferred for a shorthand factor.
    ///
    /// This deliberately reuses ordinary syntactic inference so sums, products,
    /// powers, functions, and scalar factors follow the same conventions here.
    fn append_syntactic_slots(
        value: AtomView<'_>,
        slots: &mut Vec<Slot<LibraryRep, Aind>>,
        matcher: &mut SlotMatcher,
    ) -> Result<(), StructureError> {
        let structure = Self::syntactic_structure_from_atom(value, matcher)?;
        slots.extend(
            structure
                .layout()
                .canonical_to_logical(&structure.canonical().structure),
        );
        Ok(())
    }
}

impl<Aind: AbsInd + ParseableAind> Canonicalized<OrderedStructure<LibraryRep, Aind>> {
    /// Use the canonical merger's exact incidence decision, retaining the
    /// original left-to-right order of the surviving ports.
    fn merge_syntactic(&self, other: &Self) -> Result<Self, StructureError> {
        let (canonical, left, right, _) = self.canonical().merge(other.canonical())?;
        let surviving = |source: &Self, contracted: Vec<bool>| {
            let ports = source
                .canonical()
                .structure
                .iter()
                .copied()
                .zip(contracted)
                .collect::<Vec<_>>();
            source
                .layout()
                .canonical_to_logical(&ports)
                .into_iter()
                .filter_map(|(slot, contracted)| (!contracted).then_some(slot))
        };
        let result = OrderedStructure::new(
            surviving(self, left.iter().collect())
                .chain(surviving(other, right.iter().collect()))
                .collect(),
        );
        debug_assert_eq!(result.canonical(), &canonical);
        Ok(result)
    }
}

impl<Aind: AbsInd + DummyAind + ParseableAind> StructureFromAtom for ShadowedStructure<Aind> {
    fn structure_from_atom(
        value: AtomView<'_>,
        matcher: &mut SlotMatcher,
    ) -> Result<Canonicalized<Self>, StructureError> {
        TensorialSyntax::validate_chain_like_nesting(value, None, matcher.tags())
            .map_err(|error| StructureError::ParsingError(error.to_string()))?;
        Self::from_fast_atom(value, matcher)
    }

    fn structure_from_parser<Index>(
        state: &ParseState<Index, AtomView<'_>>,
    ) -> Result<Canonicalized<Self>, StructureError> {
        let mut matcher = state.matcher.borrow_mut();
        if state.chain_scope_validated {
            Self::from_fast_atom(state.current_view(), &mut matcher)
        } else {
            Self::structure_from_atom(state.current_view(), &mut matcher)
        }
    }
}

impl<Aind: AbsInd + ParseableAind> NamedStructure<Symbol, Vec<Atom>, LibraryRep, Aind> {
    /// Infer a named structure with the fast syntactic conventions.
    fn from_fast_atom(
        value: AtomView<'_>,
        matcher: &mut SlotMatcher,
    ) -> Result<Canonicalized<Self>, StructureError> {
        match value {
            AtomView::Fun(fun)
                if fun.get_symbol() != matcher.tags().chain
                    && fun.get_symbol() != matcher.tags().trace
                    && fun.get_symbol() != matcher.tags().dot
                    && matcher.compact_inner_product_parts::<Aind>(fun).is_none() =>
            {
                if !fun.get_symbol().has_tag(&matcher.tags().tensor) {
                    OrderedStructure::<LibraryRep, Aind>::from_syntactic_atom(value, matcher)?;
                }
                Self::from_fast_function(fun, matcher)
            }
            _ => OrderedStructure::<LibraryRep, Aind>::leaf_structure_from_atom(value, matcher)
                .map(|structure| Self::from_ordered_atom(value, structure)),
        }
    }

    /// Infer a named structure from an ordinary function leaf.
    ///
    /// Direct slot arguments define the exposed structure. Nested `aind(...)`
    /// bundles are flattened and malformed bundles return their slot parsing
    /// error, while non-structural arguments are retained as metadata on the
    /// named leaf.
    fn from_fast_function(
        value: FunView<'_>,
        matcher: &mut SlotMatcher,
    ) -> Result<Canonicalized<Self>, StructureError> {
        match value.get_symbol() {
            s if s == matcher.index_bundle() => {
                let mut structure = Vec::new();
                for arg in value.iter() {
                    structure.push(matcher.parse::<LibraryRep, Aind>(arg)?);
                }

                Ok(OrderedStructure::new(structure).map_canonical(Into::into))
            }
            name => {
                let mut args = Vec::new();
                let mut slots = Vec::new();
                let mut is_structure: Option<StructureError> =
                    Some(SlotError::EmptyStructure.into());

                for arg in value.iter() {
                    let recognized = matcher.classify(arg);
                    match recognized.parse::<LibraryRep, Aind>(arg, matcher) {
                        Ok(slot) => {
                            is_structure = None;
                            slots.push(slot);
                        }
                        Err(err) => {
                            if matches!(recognized, SlotMatch::Explicit(slot) if !slot.is_concrete_index())
                            {
                                return Err(err.into());
                            }
                            if let AtomView::Fun(fun) = arg
                                && fun.get_symbol() == matcher.index_bundle()
                            {
                                let structure = Self::from_fast_function(fun, matcher)?;
                                let internal_slots = structure.layout().canonical_to_logical(
                                    &structure.canonical().external_structure(),
                                );
                                slots.extend(internal_slots);
                                is_structure = None;
                                continue;
                            }
                            if slots.is_empty() {
                                is_structure = Some(err.into());
                            }
                            args.push(arg.to_owned());
                        }
                    }
                }

                if let Some(err) = is_structure
                    && !name.has_tag(&SPENSO_TAG.tensor)
                {
                    return Err(err);
                }

                Ok(OrderedStructure::new(slots).map_canonical(|structure| {
                    let mut structure: Self = structure.into();
                    structure.set_name(name);
                    if !args.is_empty() {
                        structure.additional_args = Some(args);
                    }
                    structure
                }))
            }
        }
    }

    /// Wrap an inferred ordered structure with the original symbolic leaf name.
    fn from_ordered_atom(
        value: AtomView<'_>,
        structure: Canonicalized<OrderedStructure<LibraryRep, Aind>>,
    ) -> Canonicalized<Self> {
        structure.map_canonical(|structure| {
            let mut named = NamedStructure::from(structure);
            if let AtomView::Fun(fun) = value {
                named.global_name = Some(fun.get_symbol());
                let args = Self::leaf_additional_args(fun);
                if !args.is_empty() {
                    named.additional_args = Some(args);
                }
            }
            named
        })
    }

    /// Keep non-structural function arguments as leaf metadata.
    ///
    /// Direct slot arguments are represented by the structure; chain endpoints
    /// are also structural and therefore not duplicated as metadata.
    fn leaf_additional_args(fun: FunView<'_>) -> Vec<Atom> {
        let mut matcher = SlotMatcher::default();
        let args = fun.iter().collect::<Vec<_>>();
        if fun.get_symbol() == SPENSO_TAG.chain {
            return args[2..].iter().map(|arg| arg.to_owned()).collect();
        }
        if fun.get_symbol() == SPENSO_TAG.trace {
            return args.iter().map(|arg| arg.to_owned()).collect();
        }

        args.into_iter()
            .filter(|arg| !Self::is_direct_structure_arg(*arg, &mut matcher))
            .map(|arg| arg.to_owned())
            .collect()
    }

    /// Return true for arguments that are represented by the inferred structure.
    fn is_direct_structure_arg(arg: AtomView<'_>, matcher: &mut SlotMatcher) -> bool {
        matcher.parse::<LibraryRep, Aind>(arg).is_ok()
            || matches!(arg, AtomView::Fun(fun) if fun.get_symbol() == matcher.index_bundle())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        bracket, chain, slot,
        structure::{
            TensorStructure,
            abstract_index::{AIND_SYMBOLS, AbstractIndex},
            representation::{Lorentz, Minkowski, RepName},
            slot::IsAbstractSlot,
        },
        tensor_symbol, trace, vector,
    };
    use symbolica::{
        atom::{Atom, AtomCore, FunctionBuilder, Symbol},
        function, symbol,
    };

    // Independent oracle: execute the existing parser's shorthand lowering,
    // then read the graph incidence instead of invoking syntax inference again.
    fn expanded_graph_structure(expression: &Atom) -> Canonicalized<OrderedStructure> {
        use crate::network::parsing::{NetworkParse, ParseSettings};
        let network = expression
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        OrderedStructure::from_slots(network.graph.dangling_indices()).unwrap()
    }

    #[test]
    fn opaque_inner_product_preserves_explicit_spectators_and_chain_boundaries() {
        let rep = Minkowski {}.new_rep(4);
        let compact = rep.to_symbolic([]);
        let left = rep.to_symbolic([Atom::num(91201)]);
        let inner = rep.to_symbolic([Atom::num(91202)]);
        let right = rep.to_symbolic([Atom::num(91203)]);
        let head = tensor_symbol!("opaque_inner_product_word");
        let direct = function!(
            SPENSO_TAG.dot,
            function!(head, &left, &inner, &compact),
            function!(head, &inner, &right, &compact)
        );
        let chained = function!(
            SPENSO_TAG.dot,
            chain!(
                &left,
                &inner,
                function!(head, SPENSO_TAG.chain_in, SPENSO_TAG.chain_out, &compact)
            ),
            chain!(
                &inner,
                &right,
                function!(head, SPENSO_TAG.chain_in, SPENSO_TAG.chain_out, &compact)
            )
        );
        let mut matcher = SlotMatcher::default();
        let expected = vec![
            matcher
                .parse::<LibraryRep, AbstractIndex>(left.as_view())
                .unwrap(),
            matcher
                .parse::<LibraryRep, AbstractIndex>(right.as_view())
                .unwrap(),
        ];
        for input in [direct, chained] {
            let structure = OrderedStructure::<LibraryRep, AbstractIndex>::structure_from_atom(
                input.as_view(),
                &mut matcher,
            )
            .unwrap();
            assert_eq!(
                structure
                    .layout()
                    .canonical_to_logical(&structure.canonical().structure),
                expected
            );
            let named = ShadowedStructure::<AbstractIndex>::structure_from_atom(
                input.as_view(),
                &mut matcher,
            )
            .unwrap();
            assert_eq!(
                named
                    .layout()
                    .canonical_to_logical(&named.canonical().external_structure()),
                expected
            );
            assert_eq!(
                structure.canonical(),
                expanded_graph_structure(&input).canonical()
            );
        }
    }

    #[test]
    fn lorentz_dimension_changes_only_four_dimensional_minkowski_metadata() {
        let input = symbolica::parse!(
            "4*F(spenso::mink(4,mu),spenso::mink(4),spenso::mink(3,nu),spenso::bis(4,i),spenso::euc(4,j))"
        );
        let d = symbolica::parse!("D");
        let expected = symbolica::parse!(
            "4*F(spenso::mink(D,mu),spenso::mink(D),spenso::mink(3,nu),spenso::bis(4,i),spenso::euc(4,j))"
        );
        let result = input.with_lorentz_dimension(d.as_view());
        assert_eq!(result, expected);
        assert_eq!(result.with_lorentz_dimension(d.as_view()), expected);
    }

    fn mink4() -> crate::structure::representation::Representation<Minkowski> {
        Minkowski {}.new_rep(4)
    }

    fn chain_factor_with_external(name: Symbol, external: Atom) -> Atom {
        FunctionBuilder::new(name)
            .add_arg(external)
            .add_arg(Atom::var(SPENSO_TAG.chain_in))
            .add_arg(Atom::var(SPENSO_TAG.chain_out))
            .finish()
    }

    #[test]
    fn tensorial_syntax_detects_representation_tags() {
        let rep = mink4();
        let compact = vector!(structure_inference_p, rep.to_symbolic([]));
        let scalar = function!(symbol!("f"), Atom::num(1));
        let scalar_with_tensor_arg = function!(symbol!("f"), rep.to_symbolic([]));
        let tagged_scalar = function!(tensor_symbol!(structure_inference_t), Atom::num(1));
        let tagged_rank_zero =
            FunctionBuilder::new(tensor_symbol!(structure_inference_scalar)).finish();
        let bracketed = bracket!(compact.clone());
        let nested = scalar.clone() + compact.clone().pow(2);

        assert!(compact.is_tensorial(StrictTensorFilter::Tagged));
        assert!(compact.as_view().is_tensorial(StrictTensorFilter::Tagged));
        assert!(bracketed.is_tensorial(StrictTensorFilter::Tagged));
        assert!(nested.is_tensorial(StrictTensorFilter::Tagged));
        assert!(tagged_scalar.is_tensorial(StrictTensorFilter::Tagged));
        assert!(!scalar.is_tensorial(StrictTensorFilter::Tagged));
        assert!(!scalar_with_tensor_arg.is_tensorial(StrictTensorFilter::Tagged));

        assert!(!tagged_scalar.is_tensorial(StrictTensorFilter::TaggedChecked));
        assert!(tagged_rank_zero.is_tensorial(StrictTensorFilter::TaggedChecked));
        assert!(compact.is_tensorial(StrictTensorFilter::TaggedChecked));

        assert!(scalar_with_tensor_arg.is_tensorial(StrictTensorFilter::ContainsReps));
        assert!(!scalar.is_tensorial(StrictTensorFilter::ContainsReps));
    }

    #[test]
    fn tensorial_syntax_recognizes_signed_projectors_without_opening_metadata() {
        let rep = mink4();
        let a = vector!(projector_syntax_a, rep.to_symbolic([]));
        let b = vector!(projector_syntax_b, rep.to_symbolic([]));
        let scalar = Atom::var(symbol!("projector_syntax_scalar"));
        let opaque = symbol!("projector_syntax_opaque");
        let metadata = symbol!("projector_syntax_metadata"; Scalar);
        for head in [*shadowing::SYM, *shadowing::ANTISYM, *shadowing::CYCLIC] {
            let projector = function!(head, &a, &b);
            let signed = -&projector;
            assert!(
                matches!(projector.as_view(), AtomView::Mul(_))
                    || matches!(signed.as_view(), AtomView::Mul(_))
            );
            for filter in [
                StrictTensorFilter::Tagged,
                StrictTensorFilter::TaggedChecked,
            ] {
                for expression in [&projector, &signed] {
                    assert!(expression.is_tensorial(filter), "{expression}");
                    assert!(!function!(opaque, expression).is_tensorial(filter));
                    assert!(!function!(metadata, expression).is_tensorial(filter));
                    assert!(!function!(SPENSO_TAG.pure_scalar, expression).is_tensorial(filter));
                }
                assert!(!function!(head, &scalar, Atom::num(2)).is_tensorial(filter));
                assert!(!FunctionBuilder::new(head).finish().is_tensorial(filter));
            }
        }
    }

    #[test]
    fn tensorial_syntax_recognizes_dual_representation_wrappers() {
        let rep = Lorentz {}.new_rep(4);
        for expression in [rep.dual().to_symbolic([]), slot!(rep.dual(), i).to_atom()] {
            for filter in [
                StrictTensorFilter::Tagged,
                StrictTensorFilter::TaggedChecked,
                StrictTensorFilter::ContainsReps,
            ] {
                assert!(expression.is_tensorial(filter), "{expression}");
            }
        }
        for wrapper in [
            AIND_SYMBOLS.dind,
            AIND_SYMBOLS.uind,
            AIND_SYMBOLS.selfdualind,
        ] {
            for argument in [rep.to_symbolic([]), slot!(rep, i).to_atom()] {
                let expression = function!(wrapper, argument);
                for filter in [
                    StrictTensorFilter::Tagged,
                    StrictTensorFilter::TaggedChecked,
                    StrictTensorFilter::ContainsReps,
                ] {
                    assert!(expression.is_tensorial(filter), "{expression}");
                }
            }
            // A variance wrapper is not a generic tensor-valued function.
            // Concrete component indices likewise do not carry tensor ports.
            for argument in [
                Atom::num(1),
                function!(AIND_SYMBOLS.cind, Atom::num(1)),
                Atom::var(symbol!("dual_scalar_argument")),
            ] {
                let expression = function!(wrapper, argument);
                assert!(
                    !expression.is_tensorial(StrictTensorFilter::Tagged),
                    "{expression}"
                );
                assert!(
                    !expression.is_tensorial(StrictTensorFilter::TaggedChecked),
                    "{expression}"
                );
            }
        }
        let vector = vector!(dual_wrapped_vector, rep.to_symbolic([]));
        let lower = function!(AIND_SYMBOLS.dind, &vector);
        assert!(!lower.is_tensorial(StrictTensorFilter::Tagged));
        assert!(!lower.is_tensorial(StrictTensorFilter::TaggedChecked));
        // The other wrappers normalize to their arguments before classification.
        for wrapper in [AIND_SYMBOLS.uind, AIND_SYMBOLS.selfdualind] {
            let expression = function!(wrapper, &vector);
            assert_eq!(expression, vector);
            assert!(expression.is_tensorial(StrictTensorFilter::Tagged));
            assert!(expression.is_tensorial(StrictTensorFilter::TaggedChecked));
        }
    }

    #[test]
    fn direct_slots_reject_unsupported_indices_but_keep_component_metadata() {
        let rep = mink4();
        let slot = |index: Atom| rep.to_symbolic([index]);
        let invalid = slot(function!(symbol!("direct_invalid_index"), 7));
        let valid = rep.slot::<AbstractIndex, _>(11).to_atom();
        let head = tensor_symbol!(direct_slot_admission);
        for arguments in [
            vec![invalid.clone()],
            vec![valid.clone(), invalid.clone()],
            vec![invalid.clone(), valid.clone()],
        ] {
            let value = arguments
                .into_iter()
                .fold(FunctionBuilder::new(head), |builder, argument| {
                    builder.add_arg(argument)
                })
                .finish();
            assert!(matches!(
                value.infer_structure::<OrderedStructure>(),
                Err(StructureError::SlotError(SlotError::AindError(_)))
            ));
            assert!(matches!(
                value.infer_structure::<ShadowedStructure<AbstractIndex>>(),
                Err(StructureError::SlotError(SlotError::AindError(_)))
            ));
        }
        let scalar = symbol!("direct_slot_metadata"; Scalar);
        let metadata = [
            function!(scalar, invalid),
            rep.to_symbolic([]),
            slot(function!(AIND_SYMBOLS.cind, 2)),
            slot(function!(AIND_SYMBOLS.find, 3)),
        ];
        let value = metadata
            .iter()
            .fold(
                FunctionBuilder::new(head).add_arg(&valid),
                |builder, argument| builder.add_arg(argument),
            )
            .finish();
        let ordered = value.infer_structure::<OrderedStructure>().unwrap();
        let named = value
            .infer_structure::<ShadowedStructure<AbstractIndex>>()
            .unwrap();
        assert_eq!(ordered.canonical().order(), 1);
        assert_eq!(
            named.canonical().external_structure(),
            ordered.canonical().external_structure()
        );
        assert_eq!(
            named.canonical().additional_args.as_deref(),
            Some(metadata.as_slice())
        );

        let index = function!(crate::index_symbol!(direct_named_index), 19, 2);
        let named_slot = slot(index.clone());
        let value = function!(head, named_slot);
        let inferred = value.infer_structure::<OrderedStructure>().unwrap();
        assert_eq!(
            inferred.canonical().external_structure()[0]
                .aind()
                .to_atom(),
            index
        );
        assert_eq!(
            value
                .infer_structure::<ShadowedStructure<AbstractIndex>>()
                .unwrap()
                .canonical()
                .order(),
            1
        );
    }

    #[test]
    fn direct_slots_use_the_requested_custom_index_grammar() {
        #[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
        struct PayloadIndex(usize);
        impl std::fmt::Display for PayloadIndex {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(f, "custom({})", self.0)
            }
        }
        impl AbsInd for PayloadIndex {}
        impl ParseableAind for PayloadIndex {
            type Error = SlotError;
            fn from_view(view: AtomView<'_>) -> Result<Self, Self::Error> {
                if let AtomView::Fun(function) = view
                    && function.get_symbol() == symbol!("direct_custom_index")
                    && function.get_nargs() == 1
                {
                    return usize::try_from(function.iter().next().unwrap())
                        .map(Self)
                        .map_err(|_| SlotError::NotNatural);
                }
                Err(SlotError::Composite)
            }
            fn to_atom(&self) -> Atom {
                function!(symbol!("direct_custom_index"), self.0)
            }
        }
        impl DummyAind for PayloadIndex {
            fn new_dummy() -> Self {
                Self(usize::MAX)
            }
            fn new_dummy_at(i: usize) -> Self {
                Self(i)
            }
            fn is_dummy(&self) -> bool {
                self.0 == usize::MAX
            }
        }
        let value = function!(
            tensor_symbol!(custom_slot_admission),
            mink4().to_symbolic([PayloadIndex(17).to_atom()])
        );
        let ordered = value
            .infer_structure::<OrderedStructure<LibraryRep, PayloadIndex>>()
            .unwrap();
        let named = value
            .infer_structure::<NamedStructure<Symbol, Vec<Atom>, LibraryRep, PayloadIndex>>()
            .unwrap();
        assert_eq!(
            ordered.canonical().external_structure()[0].aind(),
            PayloadIndex(17)
        );
        assert_eq!(
            named.canonical().external_structure(),
            ordered.canonical().external_structure()
        );
        assert!(value.infer_structure::<OrderedStructure>().is_err());
    }

    #[test]
    fn tagged_tensor_rejects_malformed_aind_bundle() {
        let malformed_aind = FunctionBuilder::new(AIND_SYMBOLS.aind)
            .add_arg(Atom::num(1))
            .finish();
        let expression = FunctionBuilder::new(tensor_symbol!(malformed_aind_tensor))
            .add_arg(malformed_aind)
            .finish();

        let error = expression
            .infer_structure::<ShadowedStructure<AbstractIndex>>()
            .unwrap_err();

        assert!(matches!(
            error,
            StructureError::SlotError(SlotError::Composite)
        ));
    }

    #[test]
    fn fast_inference_retains_logical_order_through_bundles_powers_and_products() {
        let rep = Minkowski {}.new_rep(4);
        let slots = [rep.slot(31), rep.slot(7), rep.slot(19)];
        let bundle = FunctionBuilder::new(AIND_SYMBOLS.aind)
            .add_arg(slots[0].to_atom())
            .add_arg(slots[1].to_atom())
            .finish();
        let leaf = FunctionBuilder::new(tensor_symbol!(logical_inference_leaf))
            .add_arg(bundle)
            .add_arg(slots[2].to_atom())
            .finish();
        let logical = |atom: &Atom| {
            let value = atom
                .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>()
                .unwrap();
            value
                .layout()
                .canonical_to_logical(&value.canonical().external_structure())
        };
        let expected = slots.map(|slot| slot.to_lib()).to_vec();
        assert_eq!(logical(&leaf), expected);
        assert_eq!(logical(&leaf.pow(3)), expected);
        let vector = FunctionBuilder::new(tensor_symbol!(logical_inference_partner))
            .add_arg(slots[1].to_atom())
            .finish();
        assert_eq!(
            logical(&(leaf * vector)),
            vec![slots[0].to_lib(), slots[2].to_lib()]
        );
    }

    #[test]
    fn visible_slots_use_first_sum_term_as_representative() {
        let rep = mink4();
        let mu = slot!(rep, mu).to_atom();
        let expr = FunctionBuilder::new(symbol!("A"))
            .add_arg(mu.clone())
            .finish()
            + FunctionBuilder::new(symbol!("B")).add_arg(mu).finish();
        let mut slots = Vec::new();

        OrderedStructure::<LibraryRep, AbstractIndex>::append_syntactic_slots(
            expr.as_view(),
            &mut slots,
            &mut SlotMatcher::default(),
        )
        .unwrap();

        assert_eq!(slots.len(), 1);
    }

    #[test]
    fn chain_inference_matches_expanded_graph() {
        let rep = mink4();
        let external_rep = Lorentz {}.new_rep(4);
        let expr = chain!(
            slot!(rep, i),
            slot!(rep, j),
            chain_factor_with_external(
                tensor_symbol!(structure_factor_f),
                slot!(external_rep, a).to_atom()
            ),
            chain_factor_with_external(
                tensor_symbol!(structure_factor_g),
                slot!(external_rep, b).to_atom()
            ),
        );

        let fast = expr
            .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>()
            .unwrap();
        let expanded = expanded_graph_structure(&expr);

        assert_eq!(fast.canonical().order(), expanded.canonical().order());
    }

    #[test]
    fn chain_with_schoonschipped_term_inference_matches_expanded_graph() {
        let rep = mink4();
        let compact_vector = vector!(structure_inference_p, rep.to_symbolic([]));
        let schoonschipped_term = FunctionBuilder::new(tensor_symbol!(structure_factor_f))
            .add_arg(&compact_vector)
            .add_arg(Atom::var(SPENSO_TAG.chain_in))
            .add_arg(Atom::var(SPENSO_TAG.chain_out))
            .finish();
        let expr = chain!(slot!(rep, i), slot!(rep, j), schoonschipped_term);

        let fast = expr
            .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>()
            .unwrap();
        let expanded = expanded_graph_structure(&expr);

        assert_eq!(fast.canonical().order(), expanded.canonical().order());
    }

    #[test]
    fn trace_inference_matches_expanded_graph() {
        let trace_rep = Lorentz {}.new_rep(4);
        let external_rep = mink4();
        let expr = trace!(
            &trace_rep,
            chain_factor_with_external(
                tensor_symbol!(structure_factor_f),
                slot!(external_rep, a).to_atom()
            ),
            chain_factor_with_external(
                tensor_symbol!(structure_factor_g),
                slot!(external_rep, b).to_atom()
            ),
            chain_factor_with_external(
                tensor_symbol!(structure_factor_h),
                slot!(external_rep, c).to_atom()
            ),
        );

        let fast = expr
            .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>()
            .unwrap();
        let expanded = expanded_graph_structure(&expr);

        assert_eq!(fast.canonical().order(), expanded.canonical().order());
    }

    #[test]
    fn trace_projectors_keep_all_external_slots() {
        let trace_rep = Lorentz {}.new_rep(4);
        let external_rep = mink4();
        let factors = [
            chain_factor_with_external(
                tensor_symbol!(projected_structure_f),
                slot!(external_rep, a).to_atom(),
            ),
            chain_factor_with_external(
                tensor_symbol!(projected_structure_g),
                slot!(external_rep, b).to_atom(),
            ),
            chain_factor_with_external(
                tensor_symbol!(projected_structure_h),
                slot!(external_rep, c).to_atom(),
            ),
        ];
        for expression in [
            trace!(&trace_rep; factors.clone()),
            trace!(&trace_rep, shadowing::sym(factors.clone())),
            trace!(&trace_rep, -shadowing::antisym(factors)),
        ] {
            let fast = expression
                .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>()
                .unwrap();
            let expanded = expanded_graph_structure(&expression);
            assert_eq!(fast.canonical().order(), 3);
            assert_eq!(fast.canonical(), expanded.canonical());
        }
    }
}
