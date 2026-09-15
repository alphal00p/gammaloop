//! Structure inference for symbolic parser leaves.
//!
//! The fast path is intentionally syntactic:
//! 1. dispatch known shorthand (`chain`, `trace`) to its visible-slot convention;
//! 2. otherwise infer ordinary tensor syntax from sums, products, powers, and functions;
//! 3. wrap the ordered slots in the requested structure type;
//! 4. optionally validate by expanding shorthand and comparing the graph's dangling slots.

use symbolica::{
    atom::{Atom, AtomView, MulView, PowView, Symbol, representation::FunView},
    domains::rational::Rational,
};
use thiserror::Error;

use super::{NetworkParse, ParseSettings, ShadowedStructure, ShorthandParsing, StrictTensorFilter};
use crate::network::library::symbolic::ETS;
use crate::structure::{
    Canonicalized, HasName, NamedStructure, OrderedStructure, StructureContract, StructureError,
    TensorStructure,
    abstract_index::AIND_SYMBOLS,
    representation::LibraryRep,
    slot::{AbsInd, DummyAind, ParseableAind, Slot, SlotError},
};
use crate::{network::tags::SPENSO_TAG, shadowing};

/// Chooses how tensor structure is inferred from symbolic syntax.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum StructureInferenceMode {
    /// Use a syntactic walk. This is cheap and is the production default.
    Fast,
    /// Build the expanded shorthand network and read its dangling slots.
    Expanded,
}

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

pub trait StructureFromAtom: Sized {
    /// Infer the permuted tensor structure exposed by `value`.
    ///
    /// Implementations should treat `Fast` as a non-semantic syntax pass and
    /// reserve graph construction or dummy materialization for `Expanded`.
    fn structure_from_atom(
        value: AtomView<'_>,
        mode: StructureInferenceMode,
    ) -> Result<Canonicalized<Self>, StructureError>;

    /// Infer structure with the default fast syntactic mode.
    fn parse(value: AtomView<'_>) -> Result<Canonicalized<Self>, StructureError> {
        Self::structure_from_atom(value, StructureInferenceMode::Fast)
    }
}

pub trait AtomStructureExt {
    /// Convenience wrapper for `StructureFromAtom::structure_from_atom`.
    fn infer_structure<S: StructureFromAtom>(
        &self,
        mode: StructureInferenceMode,
    ) -> Result<Canonicalized<S>, StructureError>;

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
    fn infer_structure<S: StructureFromAtom>(
        &self,
        mode: StructureInferenceMode,
    ) -> Result<Canonicalized<S>, StructureError> {
        self.as_view().infer_structure(mode)
    }

    fn validate_chain_like_nesting(&self) -> Result<(), ChainNestingError> {
        self.as_view().validate_chain_like_nesting()
    }

    fn is_tensorial(&self, filter: StrictTensorFilter) -> bool {
        self.as_view().is_tensorial(filter)
    }
}

impl AtomStructureExt for AtomView<'_> {
    fn infer_structure<S: StructureFromAtom>(
        &self,
        mode: StructureInferenceMode,
    ) -> Result<Canonicalized<S>, StructureError> {
        self.validate_chain_like_nesting()
            .map_err(|error| StructureError::ParsingError(error.to_string()))?;
        S::structure_from_atom(*self, mode)
    }

    fn validate_chain_like_nesting(&self) -> Result<(), ChainNestingError> {
        TensorialSyntax::validate_chain_like_nesting(*self, None)
    }

    fn is_tensorial(&self, filter: StrictTensorFilter) -> bool {
        match *self {
            AtomView::Fun(fun) => TensorialSyntax::function_is_tensorial(fun, filter),
            AtomView::Var(var) => var.get_symbol().has_attributes_of(SPENSO_TAG.rep_),
            AtomView::Add(add) => add.iter().any(|arg| arg.is_tensorial(filter)),
            AtomView::Mul(mul) => mul.iter().any(|arg| arg.is_tensorial(filter)),
            AtomView::Pow(pow) => pow.get_base_exp().0.is_tensorial(filter),
            _ => false,
        }
    }
}

pub(crate) struct TensorialSyntax;

impl TensorialSyntax {
    fn validate_chain_like_nesting(
        value: AtomView<'_>,
        owner: Option<Symbol>,
    ) -> Result<(), ChainNestingError> {
        match value {
            AtomView::Add(sum) => {
                for term in sum.iter() {
                    Self::validate_chain_like_nesting(term, owner)?;
                }
            }
            AtomView::Mul(product) => {
                for factor in product.iter() {
                    Self::validate_chain_like_nesting(factor, owner)?;
                }
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                Self::validate_chain_like_nesting(base, owner)?;
                Self::validate_chain_like_nesting(exponent, owner)?;
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                let chain_like = symbol == SPENSO_TAG.chain || symbol == SPENSO_TAG.trace;
                if chain_like && let Some(outer) = owner {
                    return Err(ChainNestingError {
                        outer,
                        inner: symbol,
                    });
                }
                let owner = if chain_like { Some(symbol) } else { owner };
                for argument in function.iter() {
                    Self::validate_chain_like_nesting(argument, owner)?;
                }
            }
            _ => {}
        }
        Ok(())
    }

    pub(crate) fn function_is_tensorial(fun: FunView<'_>, filter: StrictTensorFilter) -> bool {
        let symbol = fun.get_symbol();

        if symbol == SPENSO_TAG.pure_scalar || symbol == SPENSO_TAG.scalar {
            return false;
        }

        if symbol == SPENSO_TAG.bracket {
            return fun.iter().any(|arg| arg.is_tensorial(filter));
        }

        if symbol.has_attributes_of(SPENSO_TAG.rep_)
            || symbol == SPENSO_TAG.chain
            || symbol == SPENSO_TAG.trace
            || symbol == SPENSO_TAG.dot
        {
            return true;
        }

        if symbol == ETS.metric {
            return fun.get_nargs() == 2 && fun.iter().all(|arg| arg.is_tensorial(filter));
        }

        if symbol.has_tag(&SPENSO_TAG.broadcast) {
            let args = fun.iter().collect::<Vec<_>>();
            return matches!(args.as_slice(), [arg] if arg.is_tensorial(filter));
        }

        match filter {
            StrictTensorFilter::Tagged => symbol.has_tag(&SPENSO_TAG.tensor),
            StrictTensorFilter::TaggedChecked => {
                symbol.has_tag(&SPENSO_TAG.tensor)
                    && (fun.get_nargs() == 0
                        || fun.iter().any(Self::contains_representation_syntax))
            }
            StrictTensorFilter::ContainsReps => {
                fun.iter().any(Self::contains_representation_syntax)
            }
        }
    }

    fn contains_representation_syntax(value: AtomView<'_>) -> bool {
        match value {
            AtomView::Fun(fun) => {
                fun.get_symbol().has_attributes_of(SPENSO_TAG.rep_)
                    || fun.iter().any(Self::contains_representation_syntax)
            }
            AtomView::Var(var) => var.get_symbol().has_attributes_of(SPENSO_TAG.rep_),
            AtomView::Add(add) => add.iter().any(Self::contains_representation_syntax),
            AtomView::Mul(mul) => mul.iter().any(Self::contains_representation_syntax),
            AtomView::Pow(pow) => Self::contains_representation_syntax(pow.get_base_exp().0),
            _ => false,
        }
    }
}

impl<Aind: AbsInd + DummyAind + ParseableAind> StructureFromAtom
    for OrderedStructure<LibraryRep, Aind>
{
    fn structure_from_atom(
        value: AtomView<'_>,
        mode: StructureInferenceMode,
    ) -> Result<Canonicalized<Self>, StructureError> {
        TensorialSyntax::validate_chain_like_nesting(value, None)
            .map_err(|error| StructureError::ParsingError(error.to_string()))?;
        match mode {
            StructureInferenceMode::Fast => Self::leaf_structure_from_atom(value),
            StructureInferenceMode::Expanded => Self::expanded_shorthand_structure_from_atom(value),
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
    ) -> Result<Canonicalized<Self>, StructureError> {
        match value {
            AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.chain => {
                Self::chain_structure_from_fun(fun)
            }
            AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.trace => {
                Self::trace_structure_from_fun(fun)
            }
            _ => Self::from_syntactic_atom(value).map(Canonicalized::identity),
        }
    }

    /// Infer an `OrderedStructure` from ordinary tensor syntax without shorthand semantics.
    ///
    /// The dispatcher is intentionally shallow: sums, products, powers, and
    /// functions each have their own convention below. Scalar syntax is carried
    /// internally as `OrderedStructure::empty()` and converted to
    /// `EmptyStructure` only at this leaf boundary.
    fn from_syntactic_atom(value: AtomView<'_>) -> Result<Self, StructureError> {
        let structure = Self::syntactic_structure_from_atom(value)?;
        if structure.is_scalar() {
            Err(StructureError::EmptyStructure(SlotError::EmptyStructure))
        } else {
            Ok(structure)
        }
    }

    fn syntactic_structure_from_atom(value: AtomView<'_>) -> Result<Self, StructureError> {
        match value {
            AtomView::Add(add) => {
                let Some(first) = add.iter().next() else {
                    return Ok(OrderedStructure::empty());
                };
                Self::syntactic_structure_from_atom(first)
            }
            AtomView::Pow(pow) => Self::from_power_atom(pow),
            AtomView::Mul(mul) => Self::from_product_atom(mul),
            AtomView::Fun(fun) => Self::from_function_atom(fun),
            _ => Ok(OrderedStructure::empty()),
        }
    }

    /// Infer an `OrderedStructure` from a power's base structure and exponent.
    ///
    /// Scalars stay scalar. A fully self-dual tensor to an even integer power
    /// has no external structure, while an odd integer power keeps the base
    /// structure. Fractional powers and powers of non-self-dual tensors are
    /// rejected because their external structure is not well-defined here.
    fn from_power_atom(pow: PowView<'_>) -> Result<Self, StructureError> {
        let (base, exp) = pow.get_base_exp();
        let base_structure = Self::syntactic_structure_from_atom(base)?;

        if base_structure.is_scalar() {
            Ok(base_structure)
        } else if base_structure.is_fully_self_dual()
            && let Ok(r) = Rational::try_from(exp)
        {
            if r.numerator() % 2 == 0 {
                Ok(OrderedStructure::empty())
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
    fn from_product_atom(product: MulView<'_>) -> Result<Self, StructureError> {
        let mut structure = OrderedStructure::empty();

        for factor in product {
            structure = structure
                .merge(&Self::syntactic_structure_from_atom(factor)?)?
                .0;
        }

        Ok(structure)
    }

    /// Infer an `OrderedStructure` from a generic function's direct structural arguments.
    ///
    /// A direct slot argument contributes one exposed slot. An `aind(...)`
    /// bundle is flattened into its slots. Other arguments are treated as
    /// metadata for the eventual named leaf and do not erase slots already seen.
    fn from_function_atom(fun: FunView<'_>) -> Result<Self, StructureError> {
        if fun.get_symbol() == AIND_SYMBOLS.aind {
            let mut slots = Vec::new();
            for arg in fun.iter() {
                slots.push(arg.try_into()?);
            }
            return Ok(OrderedStructure::new(slots).into_canonical());
        }

        let mut slots = Vec::new();

        for arg in fun.iter() {
            match Slot::<LibraryRep, Aind>::try_from(arg) {
                Ok(slot) => {
                    slots.push(slot);
                }
                Err(_) => {
                    if let AtomView::Fun(fun) = arg
                        && fun.get_symbol() == AIND_SYMBOLS.aind
                    {
                        let internal = Self::from_function_atom(fun)?;
                        slots.extend(internal.structure);
                    }
                }
            }
        }

        Ok(OrderedStructure::new(slots).into_canonical())
    }

    /// Infer an `OrderedStructure` from expanded shorthand by reading graph dangling slots.
    ///
    /// This is an oracle/debug path: it allocates any dummies required by
    /// expansion, builds the network, and then throws away everything except the
    /// external slots.
    fn expanded_shorthand_structure_from_atom(
        value: AtomView<'_>,
    ) -> Result<Canonicalized<Self>, StructureError>
    where
        Aind: DummyAind,
    {
        let network = value
            .parse_to_atom_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::expand_all(),
                ..Default::default()
            })
            .map_err(|err| StructureError::ParsingError(err.to_string()))?;

        Self::from_slots(network.graph.dangling_indices())
    }

    /// Infer an `OrderedStructure` from an opaque open chain.
    ///
    /// `args[0]` and `args[1]` are the external endpoints. Remaining factors
    /// may contain other external slots, so they are scanned recursively. The
    /// symbolic placeholders `in` and `out` are just wiring labels and are not
    /// materialized as dummies in this mode.
    fn chain_structure_from_fun(fun: FunView<'_>) -> Result<Canonicalized<Self>, StructureError> {
        let args = fun.iter().collect::<Vec<_>>();
        if args.len() < 2 {
            return Err(StructureError::WrongNumberOfArguments(args.len(), 2));
        }

        let mut slots = vec![
            Slot::<LibraryRep, Aind>::try_from(args[0])?,
            Slot::<LibraryRep, Aind>::try_from(args[1])?,
        ];
        for factor in &args[2..] {
            Self::append_syntactic_slots(*factor, &mut slots)?;
        }

        Self::from_slots(slots)
    }

    /// Infer an `OrderedStructure` from an opaque trace.
    ///
    /// `args[0]` is the traced representation, not an exposed slot. The factors
    /// are scanned for any non-placeholder slots that remain external to the
    /// trace shorthand.
    fn trace_structure_from_fun(fun: FunView<'_>) -> Result<Canonicalized<Self>, StructureError> {
        let args = fun.iter().collect::<Vec<_>>();
        if args.is_empty() {
            return Err(StructureError::WrongNumberOfArguments(0, 1));
        }

        let mut slots = Vec::new();
        for factor in shadowing::trace_factor_views(&args[1..]) {
            Self::append_syntactic_slots(factor, &mut slots)?;
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
    ) -> Result<(), StructureError> {
        slots.extend(Self::syntactic_structure_from_atom(value)?.structure);
        Ok(())
    }
}

impl<Aind: AbsInd + DummyAind + ParseableAind> StructureFromAtom for ShadowedStructure<Aind> {
    fn structure_from_atom(
        value: AtomView<'_>,
        mode: StructureInferenceMode,
    ) -> Result<Canonicalized<Self>, StructureError> {
        TensorialSyntax::validate_chain_like_nesting(value, None)
            .map_err(|error| StructureError::ParsingError(error.to_string()))?;
        match mode {
            StructureInferenceMode::Fast => Self::from_fast_atom(value),
            StructureInferenceMode::Expanded => {
                OrderedStructure::<LibraryRep, Aind>::structure_from_atom(value, mode)
                    .map(|structure| Self::from_ordered_atom(value, structure))
            }
        }
    }
}

impl<Aind: AbsInd + ParseableAind> NamedStructure<Symbol, Vec<Atom>, LibraryRep, Aind> {
    /// Infer a named structure with the fast syntactic conventions.
    fn from_fast_atom(value: AtomView<'_>) -> Result<Canonicalized<Self>, StructureError> {
        match value {
            AtomView::Fun(fun)
                if fun.get_symbol() != SPENSO_TAG.chain && fun.get_symbol() != SPENSO_TAG.trace =>
            {
                if !fun.get_symbol().has_tag(&SPENSO_TAG.tensor) {
                    OrderedStructure::<LibraryRep, Aind>::from_syntactic_atom(value)?;
                }
                Self::from_fast_function(fun)
            }
            _ => OrderedStructure::<LibraryRep, Aind>::leaf_structure_from_atom(value)
                .map(|structure| Self::from_ordered_atom(value, structure)),
        }
    }

    /// Infer a named structure from an ordinary function leaf.
    ///
    /// Direct slot arguments define the exposed structure. Nested `aind(...)`
    /// bundles are flattened and malformed bundles return their slot parsing
    /// error, while non-structural arguments are retained as metadata on the
    /// named leaf.
    fn from_fast_function(value: FunView<'_>) -> Result<Canonicalized<Self>, StructureError> {
        match value.get_symbol() {
            s if s == AIND_SYMBOLS.aind => {
                let mut structure = Vec::new();
                for arg in value.iter() {
                    structure.push(arg.try_into()?);
                }

                Ok(OrderedStructure::new(structure).map_canonical(Into::into))
            }
            name => {
                let mut args = Vec::new();
                let mut slots = Vec::new();
                let mut is_structure: Option<StructureError> =
                    Some(SlotError::EmptyStructure.into());

                for arg in value.iter() {
                    match Slot::<LibraryRep, Aind>::try_from(arg) {
                        Ok(slot) => {
                            is_structure = None;
                            slots.push(slot);
                        }
                        Err(err) => {
                            if let AtomView::Fun(fun) = arg
                                && fun.get_symbol() == AIND_SYMBOLS.aind
                            {
                                let structure = Self::from_fast_function(fun)?;
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
        let args = fun.iter().collect::<Vec<_>>();
        if fun.get_symbol() == SPENSO_TAG.chain {
            return args[2..].iter().map(|arg| arg.to_owned()).collect();
        }
        if fun.get_symbol() == SPENSO_TAG.trace {
            return args.iter().map(|arg| arg.to_owned()).collect();
        }

        args.into_iter()
            .filter(|arg| !Self::is_direct_structure_arg(*arg))
            .map(|arg| arg.to_owned())
            .collect()
    }

    /// Return true for arguments that are represented by the inferred structure.
    fn is_direct_structure_arg(arg: AtomView<'_>) -> bool {
        Slot::<LibraryRep, Aind>::try_from(arg).is_ok()
            || matches!(arg, AtomView::Fun(fun) if fun.get_symbol() == AIND_SYMBOLS.aind)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        bracket, chain, slot,
        structure::{
            TensorStructure,
            abstract_index::AbstractIndex,
            representation::{Lorentz, Minkowski, RepName},
            slot::IsAbstractSlot,
        },
        tensor_symbol, trace, vector,
    };
    use symbolica::{
        atom::{Atom, AtomCore, FunctionBuilder, Symbol},
        function, symbol,
    };

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
    fn tagged_tensor_rejects_malformed_aind_bundle() {
        let malformed_aind = FunctionBuilder::new(AIND_SYMBOLS.aind)
            .add_arg(Atom::num(1))
            .finish();
        let expression = FunctionBuilder::new(tensor_symbol!(malformed_aind_tensor))
            .add_arg(malformed_aind)
            .finish();

        let error = expression
            .infer_structure::<ShadowedStructure<AbstractIndex>>(StructureInferenceMode::Fast)
            .unwrap_err();

        assert!(matches!(
            error,
            StructureError::SlotError(SlotError::Composite)
        ));
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
        )
        .unwrap();

        assert_eq!(slots.len(), 1);
    }

    #[test]
    fn chain_fast_and_expanded_inference_agree() {
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
            .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>(
                StructureInferenceMode::Fast,
            )
            .unwrap();
        let expanded = expr
            .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>(
                StructureInferenceMode::Expanded,
            )
            .unwrap();

        assert_eq!(fast.canonical().order(), expanded.canonical().order());
    }

    #[test]
    fn chain_with_schoonschipped_term_fast_and_expanded_inference_agree() {
        let rep = mink4();
        let compact_vector = vector!(structure_inference_p, rep.to_symbolic([]));
        let schoonschipped_term = FunctionBuilder::new(tensor_symbol!(structure_factor_f))
            .add_arg(&compact_vector)
            .add_arg(Atom::var(SPENSO_TAG.chain_in))
            .add_arg(Atom::var(SPENSO_TAG.chain_out))
            .finish();
        let expr = chain!(slot!(rep, i), slot!(rep, j), schoonschipped_term);

        let fast = expr
            .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>(
                StructureInferenceMode::Fast,
            )
            .unwrap();
        let expanded = expr
            .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>(
                StructureInferenceMode::Expanded,
            )
            .unwrap();

        assert_eq!(fast.canonical().order(), expanded.canonical().order());
    }

    #[test]
    fn trace_fast_and_expanded_inference_agree() {
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
            .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>(
                StructureInferenceMode::Fast,
            )
            .unwrap();
        let expanded = expr
            .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>(
                StructureInferenceMode::Expanded,
            )
            .unwrap();

        assert_eq!(fast.canonical().order(), expanded.canonical().order());
    }
}
