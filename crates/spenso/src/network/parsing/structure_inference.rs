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

use super::{NetworkParse, ParseSettings, ShadowedStructure, ShorthandParsing};
use crate::network::tags::SPENSO_TAG;
use crate::structure::{
    NamedStructure, OrderedStructure, PermutedStructure, StructureContract, StructureError,
    TensorStructure,
    abstract_index::AIND_SYMBOLS,
    representation::LibraryRep,
    slot::{AbsInd, DummyAind, ParseableAind, Slot, SlotError},
};

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
}

pub trait AtomStructureExt {
    /// Convenience wrapper for `StructureFromAtom::structure_from_atom`.
    fn infer_structure<S: StructureFromAtom>(
        &self,
        mode: StructureInferenceMode,
    ) -> Result<PermutedStructure<S>, StructureError>;
}

impl AtomStructureExt for Atom {
    fn infer_structure<S: StructureFromAtom>(
        &self,
        mode: StructureInferenceMode,
    ) -> Result<PermutedStructure<S>, StructureError> {
        self.as_view().infer_structure(mode)
    }
}

impl AtomStructureExt for AtomView<'_> {
    fn infer_structure<S: StructureFromAtom>(
        &self,
        mode: StructureInferenceMode,
    ) -> Result<PermutedStructure<S>, StructureError> {
        S::structure_from_atom(*self, mode)
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
    /// `EmptyStructure` only at this public leaf boundary.
    pub fn from_syntactic_atom(value: AtomView<'_>) -> Result<Self, StructureError> {
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
                    }
                }
            }
        }

        Ok(OrderedStructure::new(slots).structure)
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
                shorthand_parsing: ShorthandParsing::Expand,
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
        for factor in &args[1..] {
            Self::append_syntactic_slots(*factor, &mut slots)?;
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
        OrderedStructure::<LibraryRep, Aind>::structure_from_atom(value, mode)
            .map(|structure| Self::from_ordered_atom(value, structure))
    }
}

impl<Aind: AbsInd + ParseableAind> NamedStructure<Symbol, Vec<Atom>, LibraryRep, Aind> {
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
        chain, slot,
        structure::{
            TensorStructure,
            abstract_index::AbstractIndex,
            representation::{Lorentz, Minkowski, RepName},
            slot::IsAbstractSlot,
        },
        trace,
    };
    use symbolica::{
        atom::{Atom, FunctionBuilder, Symbol},
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
            chain_factor_with_external(symbol!("f"), slot!(external_rep, a).to_atom()),
            chain_factor_with_external(symbol!("g"), slot!(external_rep, b).to_atom()),
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
        let compact_vector = function!(symbol!("p"), rep.to_symbolic([]));
        let schoonschipped_term = FunctionBuilder::new(symbol!("f"))
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
            chain_factor_with_external(symbol!("f"), slot!(external_rep, a).to_atom()),
            chain_factor_with_external(symbol!("g"), slot!(external_rep, b).to_atom()),
            chain_factor_with_external(symbol!("h"), slot!(external_rep, c).to_atom()),
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
