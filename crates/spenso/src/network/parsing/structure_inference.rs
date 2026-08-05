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

use super::{NetworkParse, ParseSettings, ShadowedStructure, ShorthandParsing, StrictTensorFilter};
use crate::network::library::symbolic::ETS;
use crate::structure::{
    HasName, NamedStructure, OrderedStructure, PermutedStructure, StructureContract,
    StructureError, TensorStructure,
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

pub trait StructureFromAtom: Sized {
    /// Infer the permuted tensor structure exposed by `value`.
    ///
    /// Implementations should treat `Fast` as a non-semantic syntax pass and
    /// reserve graph construction or dummy materialization for `Expanded`.
    fn structure_from_atom(
        value: AtomView<'_>,
        mode: StructureInferenceMode,
    ) -> Result<PermutedStructure<Self>, StructureError>;

    /// Infer structure with the default fast syntactic mode.
    fn parse(value: AtomView<'_>) -> Result<PermutedStructure<Self>, StructureError> {
        Self::structure_from_atom(value, StructureInferenceMode::Fast)
    }
}

pub trait AtomStructureExt {
    /// Convenience wrapper for `StructureFromAtom::structure_from_atom`.
    fn infer_structure<S: StructureFromAtom>(
        &self,
        mode: StructureInferenceMode,
    ) -> Result<PermutedStructure<S>, StructureError>;

    /// Return true when this expression is valid tensor parser syntax at its root.
    ///
    /// Ordinary tensor heads expose direct structural arguments: slots,
    /// `aind(...)` bundles, or compact representation arguments. Shorthand roots
    /// (`chain`, `trace`, `dot`, and compact metrics) are tensorial because the
    /// parser gives them explicit semantics. Broadcast functions are tensorial
    /// only when they have one tensorial argument. Untagged wrappers around tensor
    /// expressions stay scalar because their head has no tensor semantics.
    fn is_tensorial(&self, filter: StrictTensorFilter) -> bool;

    /// Return true when the parser can expose tensor topology owned by this expression.
    ///
    /// Unlike [`Self::is_tensorial`], this follows ordinary unary function wrappers so
    /// their tensor-valued argument can remain an explicit `NetworkOp::Function`
    /// operand. A bare slot is only payload outside a tensor-owning head, while
    /// brackets and non-structural projector uses deliberately stop this walk.
    fn contains_exposed_tensor_topology(&self, filter: StrictTensorFilter) -> bool;
}

impl AtomStructureExt for Atom {
    fn infer_structure<S: StructureFromAtom>(
        &self,
        mode: StructureInferenceMode,
    ) -> Result<PermutedStructure<S>, StructureError> {
        self.as_view().infer_structure(mode)
    }

    fn is_tensorial(&self, filter: StrictTensorFilter) -> bool {
        self.as_view().is_tensorial(filter)
    }

    fn contains_exposed_tensor_topology(&self, filter: StrictTensorFilter) -> bool {
        self.as_view().contains_exposed_tensor_topology(filter)
    }
}

impl AtomStructureExt for AtomView<'_> {
    fn infer_structure<S: StructureFromAtom>(
        &self,
        mode: StructureInferenceMode,
    ) -> Result<PermutedStructure<S>, StructureError> {
        S::structure_from_atom(*self, mode)
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

    fn contains_exposed_tensor_topology(&self, filter: StrictTensorFilter) -> bool {
        match *self {
            AtomView::Fun(fun) => {
                TensorialSyntax::function_contains_exposed_tensor_topology(fun, filter)
            }
            // A representation-valued variable is a direct slot only when a
            // tensor-owning parent gives it that role.
            AtomView::Var(_) => false,
            AtomView::Add(add) => add
                .iter()
                .any(|arg| arg.contains_exposed_tensor_topology(filter)),
            AtomView::Mul(mul) => mul
                .iter()
                .any(|arg| arg.contains_exposed_tensor_topology(filter)),
            AtomView::Pow(pow) => pow
                .get_base_exp()
                .0
                .contains_exposed_tensor_topology(filter),
            _ => false,
        }
    }
}

pub(crate) struct TensorialSyntax;

impl TensorialSyntax {
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

        // A projector is tensor topology only in its structural direct-slot
        // form. Projectors over chain/trace factors remain ordinary opaque
        // scalar syntax until their explicit expansion pass.
        if (symbol == *shadowing::SYM
            || symbol == *shadowing::ANTISYM
            || symbol == *shadowing::CYCLIC)
            && fun.get_nargs() > 0
            && fun.iter().all(Self::is_direct_slot_syntax)
        {
            return true;
        }

        match filter {
            StrictTensorFilter::Tagged => symbol.has_tag(&SPENSO_TAG.tensor),
            StrictTensorFilter::TaggedChecked => {
                symbol.has_tag(&SPENSO_TAG.tensor)
                    && fun.iter().any(Self::contains_representation_syntax)
            }
            StrictTensorFilter::ContainsReps => {
                fun.iter().any(Self::contains_representation_syntax)
            }
        }
    }

    fn function_contains_exposed_tensor_topology(
        fun: FunView<'_>,
        filter: StrictTensorFilter,
    ) -> bool {
        let symbol = fun.get_symbol();

        if symbol == SPENSO_TAG.pure_scalar
            || symbol == SPENSO_TAG.scalar
            || symbol == SPENSO_TAG.bracket
            || symbol.has_attributes_of(SPENSO_TAG.rep_)
        {
            return false;
        }

        if symbol == *shadowing::SYM
            || symbol == *shadowing::ANTISYM
            || symbol == *shadowing::CYCLIC
        {
            return fun.get_nargs() > 0 && fun.iter().all(Self::is_direct_slot_syntax);
        }

        if Self::function_is_tensorial(fun, filter) {
            return true;
        }

        let arguments = fun.iter().collect::<Vec<_>>();
        matches!(arguments.as_slice(), [argument]
            if argument.contains_exposed_tensor_topology(filter))
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

    fn is_direct_slot_syntax(value: AtomView<'_>) -> bool {
        match value {
            AtomView::Fun(fun) => fun.get_symbol().has_attributes_of(SPENSO_TAG.rep_),
            AtomView::Var(var) => var.get_symbol().has_attributes_of(SPENSO_TAG.rep_),
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
    ) -> Result<PermutedStructure<Self>, StructureError> {
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
    ) -> Result<PermutedStructure<Self>, StructureError> {
        match value {
            AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.chain => {
                Self::chain_structure_from_fun(fun)
            }
            AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.trace => {
                Self::trace_structure_from_fun(fun)
            }
            _ => Self::from_syntactic_atom(value).map(PermutedStructure::identity),
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
    /// bundle is flattened into its slots. On tensor-tagged heads, a structural
    /// `sym(...)`, `antisym(...)`, or `cyclic(...)` group is flattened as well.
    /// Other arguments are treated as metadata for the eventual named leaf and
    /// do not erase slots already seen.
    fn from_function_atom(fun: FunView<'_>) -> Result<Self, StructureError> {
        if fun.get_symbol() == AIND_SYMBOLS.aind {
            let mut slots = Vec::new();
            for arg in fun.iter() {
                slots.push(arg.try_into()?);
            }
            return Ok(OrderedStructure::new(slots).structure);
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
                    } else if fun.get_symbol().has_tag(&SPENSO_TAG.tensor)
                        && let Some(group_slots) = Self::structural_group_slots(arg)
                    {
                        slots.extend(group_slots?);
                    }
                }
            }
        }

        Ok(OrderedStructure::new(slots).structure)
    }

    /// Return the direct slots of a recognized partial-symmetry group.
    ///
    /// Projector uses remain opaque: a group is structural only when it is
    /// nonempty and every immediate child is a direct slot. The leading minus
    /// sign produced when Symbolica normalizes an antisymmetric group is also
    /// accepted, but other products remain opaque.
    fn structural_group_slots(
        value: AtomView<'_>,
    ) -> Option<Result<Vec<Slot<LibraryRep, Aind>>, StructureError>> {
        let fun = match value {
            AtomView::Fun(fun) => fun,
            AtomView::Mul(product) => {
                let mut factors = product.iter();
                let pair = (factors.next()?, factors.next()?);
                if factors.next().is_some() {
                    return None;
                }

                let (coefficient, fun) = match pair {
                    (coefficient @ AtomView::Num(_), AtomView::Fun(fun))
                    | (AtomView::Fun(fun), coefficient @ AtomView::Num(_)) => (coefficient, fun),
                    _ => return None,
                };
                if Rational::try_from(coefficient).ok()? != -1
                    || fun.get_symbol() != *shadowing::ANTISYM
                {
                    return None;
                }
                fun
            }
            _ => return None,
        };
        let symbol = fun.get_symbol();
        if fun.get_nargs() == 0
            || (symbol != *shadowing::SYM
                && symbol != *shadowing::ANTISYM
                && symbol != *shadowing::CYCLIC)
        {
            return None;
        }

        Some(
            fun.iter()
                .map(Slot::try_from)
                .collect::<Result<_, _>>()
                .map_err(Into::into),
        )
    }

    /// Infer an `OrderedStructure` from expanded shorthand by reading graph dangling slots.
    ///
    /// This is an oracle/debug path: it allocates any dummies required by
    /// expansion, builds the network, and then throws away everything except the
    /// external slots.
    fn expanded_shorthand_structure_from_atom(
        value: AtomView<'_>,
    ) -> Result<PermutedStructure<Self>, StructureError>
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
    fn chain_structure_from_fun(
        fun: FunView<'_>,
    ) -> Result<PermutedStructure<Self>, StructureError> {
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
    fn trace_structure_from_fun(
        fun: FunView<'_>,
    ) -> Result<PermutedStructure<Self>, StructureError> {
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
    ) -> Result<PermutedStructure<Self>, StructureError> {
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
    ) -> Result<PermutedStructure<Self>, StructureError> {
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
    fn from_fast_atom(value: AtomView<'_>) -> Result<PermutedStructure<Self>, StructureError> {
        match value {
            AtomView::Fun(fun)
                if fun.get_symbol() != SPENSO_TAG.chain && fun.get_symbol() != SPENSO_TAG.trace =>
            {
                OrderedStructure::<LibraryRep, Aind>::from_syntactic_atom(value)?;
                Self::from_fast_function(fun)
            }
            _ => OrderedStructure::<LibraryRep, Aind>::leaf_structure_from_atom(value)
                .map(|structure| Self::from_ordered_atom(value, structure)),
        }
    }

    /// Infer a named structure from an ordinary function leaf.
    ///
    /// Direct slot arguments define the exposed structure. Nested `aind(...)`
    /// bundles are flattened. Structural partial-symmetry groups are flattened
    /// only for tensor-tagged heads, while all other arguments are retained as
    /// metadata on the named leaf.
    fn from_fast_function(value: FunView<'_>) -> Result<PermutedStructure<Self>, StructureError> {
        match value.get_symbol() {
            s if s == AIND_SYMBOLS.aind => {
                let mut structure = Vec::new();
                for arg in value.iter() {
                    structure.push(arg.try_into()?);
                }

                Ok(OrderedStructure::new(structure).map_structure(Into::into))
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
                                && let Ok(structure) = Self::from_fast_function(fun)
                            {
                                let mut internal_slots = structure.structure.structure.structure;
                                structure
                                    .index_permutation
                                    .apply_slice_in_place_inv(&mut internal_slots);
                                structure
                                    .rep_permutation
                                    .apply_slice_in_place_inv(&mut internal_slots);
                                slots.extend(internal_slots);
                                is_structure = None;
                                continue;
                            }
                            if name.has_tag(&SPENSO_TAG.tensor)
                                && let Some(group_slots) =
                                    OrderedStructure::<LibraryRep, Aind>::structural_group_slots(
                                        arg,
                                    )
                            {
                                slots.extend(group_slots?);
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

                if let Some(err) = is_structure {
                    return Err(err);
                }

                let mut structure: PermutedStructure<Self> =
                    OrderedStructure::new(slots).map_structure(Into::into);
                structure.structure.set_name(name);
                if !args.is_empty() {
                    structure.structure.additional_args = Some(args);
                }
                Ok(structure)
            }
        }
    }

    /// Wrap an inferred ordered structure with the original symbolic leaf name.
    fn from_ordered_atom(
        value: AtomView<'_>,
        structure: PermutedStructure<OrderedStructure<LibraryRep, Aind>>,
    ) -> PermutedStructure<Self> {
        structure.map_structure(|structure| {
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

        let structural_groups = fun.get_symbol().has_tag(&SPENSO_TAG.tensor);
        args.into_iter()
            .filter(|arg| !Self::is_direct_structure_arg(*arg, structural_groups))
            .map(|arg| arg.to_owned())
            .collect()
    }

    /// Return true for arguments that are represented by the inferred structure.
    fn is_direct_structure_arg(arg: AtomView<'_>, structural_groups: bool) -> bool {
        Slot::<LibraryRep, Aind>::try_from(arg).is_ok()
            || matches!(arg, AtomView::Fun(fun) if fun.get_symbol() == AIND_SYMBOLS.aind)
            || (structural_groups
                && matches!(
                    OrderedStructure::<LibraryRep, Aind>::structural_group_slots(arg),
                    Some(Ok(_))
                ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        antisym, bracket, chain, cyclic, slot,
        structure::{
            TensorStructure,
            abstract_index::AbstractIndex,
            representation::{Lorentz, Minkowski, RepName},
            slot::IsAbstractSlot,
        },
        sym, tensor_symbol, trace, vector,
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
        assert!(compact.is_tensorial(StrictTensorFilter::TaggedChecked));

        assert!(scalar_with_tensor_arg.is_tensorial(StrictTensorFilter::ContainsReps));
        assert!(!scalar.is_tensorial(StrictTensorFilter::ContainsReps));
    }

    #[test]
    fn exposed_topology_follows_only_supported_unary_wrappers() {
        let rep = mink4();
        let index = slot!(rep, exposed_topology_index).to_atom();
        let tensor = function!(tensor_symbol!(exposed_topology_tensor), index.clone());
        let wrapper = symbol!("exposed_topology_wrapper");
        let wrapped = function!(wrapper, tensor.clone());
        let slot_payload = function!(wrapper, index.clone());
        let bracketed = function!(wrapper, bracket!(tensor.clone()));
        let projected = function!(wrapper, sym!(tensor.clone()));
        let structural_group = sym!(index);

        assert!(tensor.contains_exposed_tensor_topology(StrictTensorFilter::Tagged));
        assert!(wrapped.contains_exposed_tensor_topology(StrictTensorFilter::Tagged));
        assert!(!slot_payload.contains_exposed_tensor_topology(StrictTensorFilter::Tagged));
        assert!(!bracketed.contains_exposed_tensor_topology(StrictTensorFilter::Tagged));
        assert!(!projected.contains_exposed_tensor_topology(StrictTensorFilter::Tagged));
        assert!(structural_group.contains_exposed_tensor_topology(StrictTensorFilter::Tagged));
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
    fn partial_symmetry_groups_expose_all_direct_slot_children() {
        let rep = mink4();
        let parameter = Atom::var(symbol!("structure_inference_parameter"));
        let c = slot!(rep, c);
        let d = slot!(rep, d);
        let antisymmetric_group = antisym!(d, c).expand();
        assert!(
            matches!(antisymmetric_group.as_view(), AtomView::Mul(_)),
            "expected normalized signed group, got {antisymmetric_group}"
        );
        let expr = FunctionBuilder::new(tensor_symbol!(structure_inference_grouped))
            .add_arg(&parameter)
            .add_arg(sym!(slot!(rep, a), slot!(rep, b)))
            .add_arg(antisymmetric_group)
            .add_arg(cyclic!(slot!(rep, e), slot!(rep, f), slot!(rep, g)))
            .finish();

        let structure = expr
            .infer_structure::<ShadowedStructure<AbstractIndex>>(StructureInferenceMode::Fast)
            .unwrap();

        assert_eq!(structure.structure.order(), 7);
        assert_eq!(structure.structure.additional_args, Some(vec![parameter]));
    }

    #[test]
    fn tensor_tagged_heads_reject_partial_symmetry_groups_with_non_slot_children() {
        let rep = mink4();
        let invalid_group = sym!(slot!(rep, b), Atom::num(1));
        let projector_group = antisym!(
            FunctionBuilder::new(tensor_symbol!(structure_inference_projector_factor))
                .add_arg(slot!(rep, c).to_atom())
                .finish()
        );
        for group in [invalid_group, projector_group] {
            let expr = FunctionBuilder::new(tensor_symbol!(structure_inference_invalid_group))
                .add_arg(slot!(rep, a).to_atom())
                .add_arg(group)
                .finish();

            assert!(
                expr.infer_structure::<ShadowedStructure<AbstractIndex>>(
                    StructureInferenceMode::Fast
                )
                .is_err()
            );
        }
    }

    #[test]
    fn untagged_heads_keep_mixed_partial_symmetry_arguments_opaque() {
        let rep = mink4();
        let invalid_group = sym!(slot!(rep, b), Atom::num(1));
        let projector_group = antisym!(
            FunctionBuilder::new(tensor_symbol!(
                structure_inference_projector_factor_untagged
            ))
            .add_arg(slot!(rep, c).to_atom())
            .finish()
        );
        let expr = FunctionBuilder::new(symbol!("structure_inference_untagged_outer"))
            .add_arg(slot!(rep, a).to_atom())
            .add_arg(&invalid_group)
            .add_arg(&projector_group)
            .finish();

        let structure = expr
            .infer_structure::<ShadowedStructure<AbstractIndex>>(StructureInferenceMode::Fast)
            .unwrap();

        assert_eq!(structure.structure.order(), 1);
        assert_eq!(
            structure.structure.additional_args,
            Some(vec![invalid_group, projector_group])
        );
    }

    #[test]
    fn untagged_heads_keep_valid_partial_symmetry_groups_opaque() {
        let rep = mink4();
        let group = sym!(slot!(rep, b), slot!(rep, c));
        let expr = FunctionBuilder::new(symbol!("structure_inference_untagged_valid_group"))
            .add_arg(slot!(rep, a).to_atom())
            .add_arg(&group)
            .finish();

        assert!(expr.is_tensorial(StrictTensorFilter::ContainsReps));
        let ordered = expr
            .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>(
                StructureInferenceMode::Fast,
            )
            .unwrap();
        let named = expr
            .infer_structure::<ShadowedStructure<AbstractIndex>>(StructureInferenceMode::Fast)
            .unwrap();

        assert_eq!(ordered.structure.order(), 1);
        assert_eq!(named.structure.order(), 1);
        assert_eq!(named.structure.additional_args, Some(vec![group]));
    }

    #[test]
    fn tensor_tagged_heads_keep_other_nested_arguments_opaque() {
        let rep = mink4();
        let bracketed_slot = bracket!(slot!(rep, b).to_atom());
        let lookalike_group = FunctionBuilder::new(symbol!(
            "structure_inference_lookalike_sym";
            Symmetric
        ))
        .add_arg(slot!(rep, c).to_atom())
        .add_arg(slot!(rep, d).to_atom())
        .finish();
        let arbitrary_product = -sym!(slot!(rep, e), slot!(rep, f));
        let expr = FunctionBuilder::new(tensor_symbol!(structure_inference_other_nested))
            .add_arg(slot!(rep, a).to_atom())
            .add_arg(&bracketed_slot)
            .add_arg(&lookalike_group)
            .add_arg(&arbitrary_product)
            .finish();

        let structure = expr
            .infer_structure::<ShadowedStructure<AbstractIndex>>(StructureInferenceMode::Fast)
            .unwrap();

        assert_eq!(structure.structure.order(), 1);
        assert_eq!(
            structure.structure.additional_args,
            Some(vec![bracketed_slot, lookalike_group, arbitrary_product])
        );
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

        assert_eq!(fast.structure.order(), expanded.structure.order());
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

        assert_eq!(fast.structure.order(), expanded.structure.order());
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

        assert_eq!(fast.structure.order(), expanded.structure.order());
    }
}
