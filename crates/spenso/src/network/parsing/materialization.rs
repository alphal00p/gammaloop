//! Materialization helpers for shorthand syntax before ordinary network parsing.
//!
//! Network-level shorthand materialization lives on the parser implementation:
//! it handles chain and trace topology and returns parsed networks. The helper
//! in this module is deliberately narrower than algebraic simplification. It
//! takes compact Schoonschip syntax that cannot be parsed as an ordinary leaf
//! yet and rewrites it into explicit slots plus ordinary tensor factors. The
//! parser can then recurse on the resulting expression and build the same graph
//! it would have built from fully expanded syntax.
//!
//! The main Schoonschip convention is:
//! 1. a tensor with one compact axis, such as `p(rep)`, used as a function
//!    argument becomes a fresh slot in that argument position; any explicit
//!    spectator slots on the tensor are preserved;
//! 2. the tensor `p(slot)` is multiplied next to the rebuilt function;
//! 3. compact scalar products `g(p(rep), q(rep.dual()))` and
//!    `dot(p(rep), q(rep.dual()))` share one fresh abstract index and become
//!    `p(rep(index)) * q(rep.dual()(index))`.
//!
//! Additional factors are accumulated beside the current atom. Compact scalar
//! factors nested inside a compact vector wrapper are recursively materialized
//! before ordinary parsing so no shorthand is left hidden inside a tensor leaf.
//!
//! Chain and trace materialization chooses the symbolic `in`/`out` replacements,
//! then lets this Schoonschip helper expand compact arguments inside each factor.

use symbolica::atom::representation::FunView;
use symbolica::atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol};
use symbolica::id::MatchSettings;

use linnet::half_edge::NodeIndex;
use std::fmt::{Debug, Display};

use super::construction::Construction;
use super::structure_inference::TensorialSyntax;
use super::{
    ParseSettings, ParseState, SchoonschipExpansionMode, StrictTensorFilter, StructureFromAtom,
    TensorFromExpression, TensorLibraryFor,
};
use crate::network::library::FunctionLibrary;
use crate::network::library::symbolic::ETS;
use crate::network::store::TensorScalarStore;
use crate::network::tags::SPENSO_TAG;
use crate::network::{Network, NetworkState, TensorNetworkError};
use crate::shadowing;
use crate::structure::representation::{LibraryRep, Representation};
use crate::structure::slot::{
    AbsInd, DualSlotTo, DummyAind, IsAbstractSlot, ParseableAind, Slot, SlotMatch, SlotMatcher,
};
use crate::structure::{HasStructure, OrderedStructure, ScalarStructure, TensorStructure};
use eyre::eyre;

struct SchoonschipMaterialization {
    /// Expression that remains at the original syntactic position.
    current: Atom,
    /// Extra factors to multiply beside `current`.
    ///
    /// These factors are explicit tensor syntax after any nested compact scalar
    /// shorthand has been recursively materialized.
    additional_factors: Vec<Atom>,
}

struct ChainEndpoint<Aind> {
    slot: Slot<LibraryRep, Aind>,
    additional_factors: Vec<Atom>,
}

impl SchoonschipMaterialization {
    /// Merge the current atom and accumulated factors into parser input.
    fn into_expression(self) -> Atom {
        self.additional_factors
            .into_iter()
            .fold(self.current, |expression, factor| expression * factor)
    }
}

/// Lowers expandable shorthand atoms into ordinary parser syntax.
///
/// The materializer owns no parse state. It borrows the parser's dummy allocator
/// so every fresh slot it creates shares the same abstract-index namespace as
/// the surrounding network parse.
pub(super) struct SchoonschipMaterializer<'a, Aind, View> {
    state: &'a ParseState<Aind, View>,
    mode: SchoonschipExpansionMode,
}

impl<'a, Aind: AbsInd + DummyAind + ParseableAind, View> SchoonschipMaterializer<'a, Aind, View> {
    /// Build a materializer with explicit Schoonschip expansion controls.
    pub(super) fn with_mode(
        state: &'a ParseState<Aind, View>,
        mode: SchoonschipExpansionMode,
    ) -> Self {
        Self { state, mode }
    }

    /// Return true when this function root contains compact Schoonschip syntax.
    ///
    /// This is the non-allocating counterpart to `materialize_shorthand`: compact
    /// scalar products are shorthand at the root, while compact vectors become
    /// shorthand only when they occur as arguments of another function.
    pub(super) fn contains_schoonschip_shorthand(&self, value: AtomView<'_>) -> bool {
        let AtomView::Fun(fun) = value else {
            return false;
        };

        if self
            .state
            .matcher
            .borrow()
            .tags()
            .is_scalar_metadata(fun.get_symbol())
        {
            return false;
        }
        self.is_compact_inner_product(fun)
            || fun
                .iter()
                .any(|value| self.contains_schoonschip_shorthand_arg(value))
    }

    /// Materialize an expandable shorthand expression, or return it unchanged.
    ///
    /// The root must be a function for rewriting to occur. When rewriting is
    /// possible, the result is one expression where the rebuilt root is
    /// multiplied by all tensor factors introduced while replacing compact
    /// vector arguments. If no shorthand is present, callers still get a valid
    /// parser input: the original atom.
    pub(super) fn materialize_shorthand(&self, value: AtomView<'_>) -> Atom {
        self.materialize_shorthand_root(value)
            .map(SchoonschipMaterialization::into_expression)
            .unwrap_or_else(|| value.to_owned())
    }

    /// Materialize shorthand at a function root.
    ///
    /// Compact scalar products get their special lowering first. Everything else
    /// is rebuilt by recursively scanning function arguments.
    fn materialize_shorthand_root(
        &self,
        value: AtomView<'_>,
    ) -> Option<SchoonschipMaterialization> {
        let AtomView::Fun(fun) = value else {
            return None;
        };

        // Scalar payloads are opaque metadata even when they contain compact
        // vectors or scalar products. Hoisting their factors changes the function.
        if self
            .state
            .matcher
            .borrow()
            .tags()
            .is_scalar_metadata(fun.get_symbol())
        {
            return None;
        }

        if self.is_chain_like_head(fun) && !self.mode.expand_inside_chains {
            return None;
        }

        if self.is_inner_product_head(fun) {
            if !self.mode.inner_products {
                return None;
            }
            if let Some(materialized) = self.compact_inner_product(fun) {
                return Some(materialized);
            }
        }

        self.materialize_shorthand_function(fun)
    }

    /// Materialize one function argument according to shorthand position rules.
    ///
    /// A compact vector has special meaning only in argument position: it
    /// contributes a fresh slot at that position and an extra rank-one tensor
    /// factor. Non-compact functions are scanned recursively.
    fn materialize_shorthand_arg(&self, value: AtomView<'_>) -> Option<SchoonschipMaterialization> {
        if self.mode.expand_schoonship
            && let Some(materialized) = self.materialize_compact_vector_arg(value)
        {
            return Some(materialized);
        }

        self.materialize_shorthand_root(value)
    }

    fn contains_schoonschip_shorthand_arg(&self, value: AtomView<'_>) -> bool {
        let compact = self
            .state
            .matcher
            .borrow_mut()
            .compact_vector_rep::<Aind>(value)
            .is_some();
        compact || self.contains_schoonschip_shorthand(value)
    }

    /// Materialize one compact vector argument as `slot` plus `vector(slot)`.
    ///
    /// The compact vector may be a single function or a sum of functions, but
    /// every visible compact representation must be the same representation.
    fn materialize_compact_vector_arg(
        &self,
        value: AtomView<'_>,
    ) -> Option<SchoonschipMaterialization> {
        let rep = self
            .state
            .matcher
            .borrow_mut()
            .compact_vector_rep::<Aind>(value)?;
        let slot = self.state.slot(&rep).to_atom();
        let factor = self.materialize_compact_vector_with_slot(value, &rep, &slot)?;
        tracing::debug!(
            target: "spenso::network::parsing",
            spenso_parser = true,
            generation = true,
            compile = true,
            inspect = true,
            stage = "schoonschip_compact_vector_slot",
            slot = %slot.to_plain_string(),
            representation = %rep,
            file.value = %value.to_plain_string(),
            file.factor = %factor.to_plain_string(),
            "Spenso parser allocated Schoonschip compact vector slot"
        );

        Some(SchoonschipMaterialization {
            current: slot,
            additional_factors: vec![factor],
        })
    }

    /// Rebuild a function after materializing any shorthand arguments.
    ///
    /// Introduced factors are accumulated next to the rebuilt function. If the
    /// rebuilt function is `dot(slot_i, slot_j)`, it is normalized to the metric
    /// spelling expected by the tensor library.
    fn materialize_shorthand_function(
        &self,
        fun: FunView<'_>,
    ) -> Option<SchoonschipMaterialization> {
        let mut changed = false;
        let mut additional_factors = Vec::new();
        let mut arguments = Vec::with_capacity(fun.get_nargs());
        let mut source_ports = Vec::new();
        let mut ports = Vec::new();
        let mut direct_arguments = true;

        for arg in fun.iter() {
            let direct_port = {
                let mut matcher = self.state.matcher.borrow_mut();
                matcher.parse::<LibraryRep, Aind>(arg).is_ok()
                    || matcher.parse_representation::<LibraryRep>(arg).is_ok()
            };
            if direct_port {
                source_ports.push(arg.to_owned());
            }
            if let Some(materialized) = self.materialize_shorthand_arg(arg) {
                changed = true;
                if self
                    .state
                    .matcher
                    .borrow_mut()
                    .parse::<LibraryRep, Aind>(materialized.current.as_view())
                    .is_ok()
                {
                    ports.push(materialized.current.clone());
                } else {
                    direct_arguments &= !TensorialSyntax::is_tensorial(
                        materialized.current.as_view(),
                        StrictTensorFilter::ContainsReps,
                        &self.state.matcher.borrow(),
                    );
                }
                arguments.push(materialized.current);
                additional_factors.extend(materialized.additional_factors);
            } else {
                arguments.push(arg.to_owned());
                if direct_port {
                    ports.push(arg.to_owned());
                } else {
                    direct_arguments &= !TensorialSyntax::is_tensorial(
                        arg,
                        StrictTensorFilter::ContainsReps,
                        &self.state.matcher.borrow(),
                    );
                }
            }
        }

        changed.then(|| {
            if let Some(cache) = &self.state.leaf_interfaces
                && let Some((true, retained)) = cache
                    .lock()
                    .expect("leaf interface cache is not poisoned")
                    .get(fun.as_view().get_data())
            {
                // A nested bundle is not a direct constructor port. Do not
                // replace its retained raw incidence with only visible slots.
                direct_arguments &= *retained == source_ports;
            }
            let replacement = FunctionBuilder::new(fun.get_symbol())
                .add_args(&arguments)
                .finish();
            if direct_arguments {
                self.state.retain_constructor_ports(
                    replacement.as_view(),
                    fun.get_symbol(),
                    &arguments,
                    (true, ports),
                );
            }
            SchoonschipMaterialization {
                current: self
                    .compact_dot_as_metric(replacement.as_view())
                    .unwrap_or(replacement),
                additional_factors,
            }
        })
    }

    /// Materialize a compact metric or dot product into two tensor factors.
    ///
    /// Both arguments must have one compact axis with matching representations.
    /// Assign one fresh abstract index with each argument's orientation. Explicit
    /// spectator slots remain unchanged, including spectators in that representation.
    fn compact_inner_product(&self, value: FunView<'_>) -> Option<SchoonschipMaterialization> {
        let (lhs, rhs, lhs_rep, rhs_rep) = self
            .state
            .matcher
            .borrow_mut()
            .compact_inner_product_parts::<Aind>(value)?;

        let lhs_slot = self.state.slot(&lhs_rep);
        let rhs_slot = rhs_rep.slot::<Aind, _>(lhs_slot.aind());
        let lhs = self.materialize_compact_vector_with_slot(lhs, &lhs_rep, &lhs_slot.to_atom())?;
        let rhs = self.materialize_compact_vector_with_slot(rhs, &rhs_rep, &rhs_slot.to_atom())?;
        tracing::debug!(
            target: "spenso::network::parsing",
            spenso_parser = true,
            generation = true,
            compile = true,
            inspect = true,
            stage = "schoonschip_compact_scalar_product_slot",
            lhs_slot = %lhs_slot.to_atom().to_plain_string(),
            rhs_slot = %rhs_slot.to_atom().to_plain_string(),
            lhs_representation = %lhs_rep,
            rhs_representation = %rhs_rep,
            file.value = %value.as_view().to_plain_string(),
            file.lhs = %lhs.to_plain_string(),
            file.rhs = %rhs.to_plain_string(),
            "Spenso parser allocated Schoonschip scalar-product slot"
        );

        Some(SchoonschipMaterialization {
            current: Atom::num(1),
            additional_factors: vec![lhs, rhs],
        })
    }

    fn is_compact_inner_product(&self, value: FunView<'_>) -> bool {
        self.state
            .matcher
            .borrow_mut()
            .compact_inner_product_parts::<Aind>(value)
            .is_some()
    }

    fn is_inner_product_head(&self, value: FunView<'_>) -> bool {
        let symbol = value.get_symbol();
        symbol == self.state.metric || symbol == self.state.matcher.borrow().tags().dot
    }

    fn is_chain_like_head(&self, value: FunView<'_>) -> bool {
        let symbol = value.get_symbol();
        symbol == self.state.matcher.borrow().tags().chain
            || symbol == self.state.matcher.borrow().tags().trace
    }

    /// Replace a compact representation with a concrete slot.
    ///
    /// For a function, this rebuilds the function with the matched compact
    /// representation argument replaced by `slot`. For a sum, every summand is
    /// rebuilt with the same slot so the expansion keeps a single dummy edge. A
    /// scalar-weighted product preserves its scalar factors around that rebuilt
    /// vector.
    fn materialize_compact_vector_with_slot(
        &self,
        value: AtomView<'_>,
        rep: &Representation<LibraryRep>,
        slot: &Atom,
    ) -> Option<Atom> {
        match value {
            AtomView::Fun(fun)
                if fun
                    .get_symbol()
                    .has_tag(&self.state.matcher.borrow().tags().broadcast) =>
            {
                let args = fun.iter().collect::<Vec<_>>();
                let [argument] = args.as_slice() else {
                    return None;
                };
                Some(
                    FunctionBuilder::new(fun.get_symbol())
                        .add_arg(self.materialize_compact_vector_with_slot(*argument, rep, slot)?)
                        .finish(),
                )
            }
            AtomView::Fun(fun) if self.state.matcher.borrow().is_projector(fun.get_symbol()) => {
                let mut changed = false;
                let mut rebuilt = FunctionBuilder::new(fun.get_symbol());
                for argument in fun.iter() {
                    if self
                        .state
                        .matcher
                        .borrow_mut()
                        .compact_vector_rep::<Aind>(argument)
                        == Some(*rep)
                    {
                        if changed {
                            return None;
                        }
                        changed = true;
                        rebuilt = rebuilt.add_arg(
                            self.materialize_compact_vector_with_slot(argument, rep, slot)?,
                        );
                    } else {
                        rebuilt = rebuilt.add_arg(argument);
                    }
                }
                changed.then(|| rebuilt.finish())
            }
            AtomView::Fun(fun) if self.is_chain_like_head(fun) => {
                let skip = if fun.get_symbol() == self.state.matcher.borrow().tags().chain {
                    2
                } else {
                    1
                };
                let (position, matched_rep) = self
                    .state
                    .matcher
                    .borrow_mut()
                    .compact_vector_sequence::<Aind>(fun.iter().enumerate().skip(skip), true)?;
                if matched_rep != *rep {
                    return None;
                }
                let mut rebuilt = FunctionBuilder::new(fun.get_symbol());
                for (i, argument) in fun.iter().enumerate() {
                    rebuilt = if i == position {
                        rebuilt.add_arg(
                            self.materialize_compact_vector_with_slot(argument, rep, slot)?,
                        )
                    } else {
                        rebuilt.add_arg(argument)
                    };
                }
                Some(rebuilt.finish())
            }
            AtomView::Fun(fun) => {
                let (matched, written) = self
                    .state
                    .matcher
                    .borrow_mut()
                    .observe_slots(|matcher| matcher.compact_tensor_rep_arg::<Aind>(fun));
                let (position, matched_rep) = matched?;
                if matched_rep != *rep {
                    return None;
                }

                let mut arguments = Vec::with_capacity(fun.get_nargs());
                let mut ports = Vec::new();
                for (arg_position, arg) in fun.iter().enumerate() {
                    if arg_position == position {
                        arguments.push(slot.clone());
                        ports.push(slot.clone());
                    } else {
                        arguments.push(arg.to_owned());
                        if written.iter().any(|port| port.as_view() == arg) {
                            ports.push(arg.to_owned());
                        }
                    }
                }
                let tensor = FunctionBuilder::new(fun.get_symbol())
                    .add_args(&arguments)
                    .finish();
                if let Some(cache) = &self.state.leaf_interfaces
                    && matches!(tensor.as_view(), AtomView::Fun(output)
                        if output.get_symbol() == fun.get_symbol()
                            && output.get_nargs() == arguments.len()
                            && output.iter().zip(&arguments).all(|(arg, supplied)| arg == supplied.as_view()))
                {
                    cache
                        .lock()
                        .expect("leaf interface cache is not poisoned")
                        .insert(tensor.clone(), (true, ports));
                }
                Some(tensor)
            }
            AtomView::Add(add) => {
                let mut terms = add
                    .iter()
                    .map(|term| self.materialize_compact_vector_with_slot(term, rep, slot));
                let first = terms.next()??;
                let rest = terms.collect::<Option<Vec<_>>>()?;
                Some(rest.into_iter().fold(first, |sum, term| sum + term))
            }
            AtomView::Mul(product) => {
                let (vector_position, matched_rep) = self
                    .state
                    .matcher
                    .borrow_mut()
                    .compact_vector_sequence::<Aind>(product.iter().enumerate(), false)?;
                if matched_rep != *rep {
                    return None;
                }

                product
                    .iter()
                    .enumerate()
                    .try_fold(Atom::num(1), |product, (position, factor)| {
                        if position == vector_position {
                            self.materialize_compact_vector_with_slot(factor, rep, slot)
                                .map(|factor| product * factor)
                        } else {
                            let scalar = self
                                .materialize_shorthand_root(factor)
                                .map(SchoonschipMaterialization::into_expression)
                                .unwrap_or_else(|| factor.to_owned());
                            Some(product * scalar)
                        }
                    })
            }
            _ => None,
        }
    }

    /// Normalize `dot(slot_i, slot_j)` to `g(slot_i, slot_j)`.
    ///
    /// This only applies after arguments have already been materialized to
    /// concrete slots. Non-slot dot products are left untouched.
    fn compact_dot_as_metric(&self, value: AtomView<'_>) -> Option<Atom> {
        let AtomView::Fun(fun) = value else {
            return None;
        };

        if fun.get_symbol() != self.state.matcher.borrow().tags().dot || fun.get_nargs() != 2 {
            return None;
        }

        let args = fun.iter().collect::<Vec<_>>();
        if args.iter().all(|arg| {
            self.state
                .matcher
                .borrow_mut()
                .parse::<LibraryRep, Aind>(*arg)
                .is_ok()
        }) {
            let arguments = args.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>();
            let metric = FunctionBuilder::new(self.state.metric)
                .add_args(&arguments)
                .finish();
            self.state.retain_constructor_ports(
                metric.as_view(),
                self.state.metric,
                &arguments,
                (true, arguments.clone()),
            );
            Some(metric)
        } else {
            None
        }
    }
}

impl SlotMatcher {
    fn is_structured_scalar<Aind: AbsInd + ParseableAind>(&mut self, value: AtomView<'_>) -> bool {
        OrderedStructure::<LibraryRep, Aind>::syntactic_structure_from_atom(value, self)
            .is_ok_and(|structure| structure.canonical().is_scalar())
    }

    /// Recognize one inner product's unique compact axis on each operand.
    /// Metrics may retain explicit spectator ports; scalar dots require rank-one
    /// operands. This observation does not allocate indices or construct components.
    #[allow(clippy::type_complexity)]
    pub fn compact_inner_product_parts<'node, Aind: AbsInd + ParseableAind>(
        &mut self,
        value: FunView<'node>,
    ) -> Option<(
        AtomView<'node>,
        AtomView<'node>,
        Representation<LibraryRep>,
        Representation<LibraryRep>,
    )> {
        if !(value.get_symbol() == self.metric() || value.get_symbol() == self.tags().dot) {
            return None;
        }

        let mut args = value.iter();
        let (lhs, rhs) = (args.next()?, args.next()?);
        if args.next().is_some() {
            return None;
        }

        let lhs_rep = self.compact_vector_rep::<Aind>(lhs)?;
        let rhs_rep = self.compact_vector_rep::<Aind>(rhs)?;
        if !lhs_rep.matches(&rhs_rep)
            || (value.get_symbol() == self.tags().dot
                && (!self.is_structured_scalar::<Aind>(lhs)
                    || !self.is_structured_scalar::<Aind>(rhs)))
        {
            return None;
        }

        Some((lhs, rhs, lhs_rep, rhs_rep))
    }

    /// Infer the unique compact axis carried by a shorthand atom.
    ///
    /// Tensor-tagged functions expose a compact representation through exactly
    /// one direct representation argument. Scalar products, unary broadcasts,
    /// projectors, and sums preserve it when they contain one compatible vector.
    /// The remaining product must have a scalar interface; its factors may
    /// carry internal contractions.
    fn compact_vector_rep<Aind: AbsInd + ParseableAind>(
        &mut self,
        value: AtomView<'_>,
    ) -> Option<Representation<LibraryRep>> {
        match value {
            AtomView::Fun(fun) if fun.get_symbol().is_scalar() => None,
            AtomView::Fun(fun) if fun.get_symbol().has_tag(&self.tags().broadcast) => {
                let args = fun.iter().collect::<Vec<_>>();
                let [argument] = args.as_slice() else {
                    return None;
                };
                self.compact_vector_rep::<Aind>(*argument)
            }
            AtomView::Fun(fun) if self.is_projector(fun.get_symbol()) => {
                let mut selected = None;
                for argument in fun.iter() {
                    if let Some(rep) = self.compact_vector_rep::<Aind>(argument)
                        && selected.replace((argument, rep)).is_some()
                    {
                        return None;
                    }
                }
                let (argument, rep) = selected?;
                for candidate in fun.iter() {
                    if candidate != argument
                        && super::structure_inference::TensorialSyntax::is_tensorial(
                            candidate,
                            StrictTensorFilter::Tagged,
                            self,
                        )
                        && !self.is_structured_scalar::<Aind>(candidate)
                    {
                        return None;
                    }
                }
                Some(rep)
            }
            AtomView::Fun(fun)
                if fun.get_symbol() == self.tags().chain
                    || fun.get_symbol() == self.tags().trace =>
            {
                let skip = if fun.get_symbol() == self.tags().chain {
                    2
                } else {
                    1
                };
                self.compact_vector_sequence::<Aind>(fun.iter().enumerate().skip(skip), true)
                    .map(|(_, rep)| rep)
            }
            AtomView::Fun(fun) => self.compact_tensor_rep_arg::<Aind>(fun).map(|(_, rep)| rep),
            AtomView::Add(add) => {
                let mut reps = add
                    .iter()
                    .map(|value| self.compact_vector_rep::<Aind>(value));
                let rep = reps.next()??;
                reps.all(|candidate| candidate == Some(rep)).then_some(rep)
            }
            AtomView::Mul(product) => self
                .compact_vector_sequence::<Aind>(product.iter().enumerate(), false)
                .map(|(_, representation)| representation),
            _ => None,
        }
    }

    /// Locate the only compact vector in a scalar-weighted product.
    ///
    /// Failure to identify a compact vector does not prove that a factor is
    /// scalar: explicit-slot tensors also have no compact representation. The
    /// ordinary syntactic structure inference is therefore the authority for
    /// the product of remaining factors, which can contain internal contractions.
    /// Chain/trace words may expose explicit spectator
    /// ports; ordinary weighted products still require scalar coefficients.
    pub fn compact_vector_sequence<'node, Aind: AbsInd + ParseableAind>(
        &mut self,
        factors: impl Iterator<Item = (usize, AtomView<'node>)>,
        allow_spectators: bool,
    ) -> Option<(usize, Representation<LibraryRep>)> {
        let mut compact_vector = None;
        let mut coefficient = Vec::new();

        for (position, factor) in factors {
            if let Some(representation) = self.compact_vector_rep::<Aind>(factor) {
                if compact_vector.is_some() {
                    return None;
                }
                compact_vector = Some((position, representation));
            } else {
                coefficient.push(factor);
            }
        }

        let structure = OrderedStructure::<LibraryRep, Aind>::syntactic_structure_from_atom(
            Atom::mul_many(coefficient).as_view(),
            self,
        )
        .ok()?;
        if !allow_spectators && !structure.canonical().is_scalar() {
            return None;
        }
        compact_vector
    }

    /// Locate the compact representation argument of one tensor function.
    ///
    /// A compact tensor function is tensor-tagged, is not itself a representation,
    /// metric or dot product, and has exactly one direct representation argument.
    /// Explicit slots are spectators: only that unique unindexed axis is replaced.
    fn compact_tensor_rep_arg<Aind: AbsInd + ParseableAind>(
        &mut self,
        value: FunView<'_>,
    ) -> Option<(usize, Representation<LibraryRep>)> {
        if value.get_symbol().is_scalar()
            || !value.get_symbol().has_tag(&self.tags().tensor)
            || value.get_symbol() == self.metric()
            || value.get_symbol() == self.tags().dot
        {
            return None;
        }

        if self.is_representation(value.as_view()) {
            return None;
        }

        let mut rep_args = value.iter().enumerate().filter_map(|(position, arg)| {
            if self.parse::<LibraryRep, Aind>(arg).is_ok() {
                None
            } else {
                self.compact_rep_pattern_match(arg)
                    .map(|rep| (position, rep))
            }
        });

        let result = rep_args.next()?;
        // Consume all arguments, including custom index readers, in the same
        // order as the previous collection into a vector.
        (rep_args.count() == 0).then_some(result)
    }

    /// Match one argument against the symbolic representation wildcard.
    ///
    /// Matching is restricted to the argument itself (`max_level = 0`) so a
    /// nested representation inside metadata does not accidentally become the
    /// tensor's compact slot.
    fn compact_rep_pattern_match(
        &mut self,
        arg: AtomView<'_>,
    ) -> Option<Representation<LibraryRep>> {
        if let Ok(representation) = self.representation_from_atom(arg) {
            return Some(representation);
        }

        let rep_pattern = Atom::var(self.tags().rep_).to_pattern();
        let settings = MatchSettings::new().max_level(0).partial(false);
        let mut matches = arg.pattern_match(&rep_pattern, None, Some(&settings));
        let matched = matches.next_detailed()?;
        let rep = rep_pattern.replace_wildcards_with_matches(matched.match_stack);
        self.representation_from_atom(rep.as_view()).ok()
    }
}

impl<Aind: AbsInd + DummyAind + ParseableAind, View> ParseState<Aind, View> {
    /// Retain ports selected by a constructor when normalization kept exactly
    /// its supplied arguments. Matching permits intrinsic argument permutations
    /// without deriving either ports or their ordering from new tensor syntax.
    fn retain_constructor_ports(
        &self,
        expression: AtomView<'_>,
        head: Symbol,
        arguments: &[Atom],
        (tensorial, mut ports): (bool, Vec<Atom>),
    ) {
        let Some(cache) = &self.leaf_interfaces else {
            return;
        };
        // Antisymmetric normalization can pull a sign out of the constructor.
        let output = match expression {
            AtomView::Fun(output) => output,
            AtomView::Mul(product)
                if product.get_nargs() == 2 && product.get_coefficient().is_some() =>
            {
                let Some(AtomView::Fun(output)) = product
                    .iter()
                    .find(|factor| matches!(factor, AtomView::Fun(_)))
                else {
                    return;
                };
                output
            }
            _ => return,
        };
        if output.get_symbol() != head || output.get_nargs() != arguments.len() {
            return;
        }
        if output
            .iter()
            .eq(arguments.iter().map(AtomCore::as_atom_view))
        {
            cache
                .lock()
                .expect("leaf interface cache is not poisoned")
                .insert(output.as_view().to_owned(), (tensorial, ports));
            return;
        }
        let mut supplied = arguments
            .iter()
            .map(AtomCore::as_atom_view)
            .collect::<Vec<_>>();
        let mut written = Vec::with_capacity(ports.len());
        for argument in output.iter() {
            let Some(position) = supplied.iter().position(|supplied| *supplied == argument) else {
                return;
            };
            supplied.swap_remove(position);
            if let Some(position) = ports.iter().position(|port| port.as_view() == argument) {
                written.push(ports.remove(position));
            }
        }
        if ports.is_empty() {
            cache
                .lock()
                .expect("leaf interface cache is not poisoned")
                .insert(output.as_view().to_owned(), (tensorial, written));
        }
    }

    /// Replace a trace by a chain with one fresh, closed channel.
    /// Factor templates remain unchanged for subsequent chain materialization.
    /// See [`Self::materialize_indexed_chain`] for opening the resulting chain.
    pub fn materialize_indexed_trace(&self, value: FunView<'_>) -> Result<Atom, String> {
        let mut arguments = value.iter();
        let representation = arguments.next().ok_or("trace requires a representation")?;
        let representation = Representation::<LibraryRep>::try_from(representation)
            .map_err(|error| error.to_string())?;
        let start = representation.slot::<Aind, _>(self.fresh_index());
        let arguments = arguments.collect::<Vec<_>>();
        Ok(FunctionBuilder::new(SPENSO_TAG.chain)
            .add_arg(start.to_atom())
            .add_arg(start.dual().to_atom())
            .add_args(shadowing::trace_factor_views(&arguments))
            .finish())
    }

    fn materialize_chain_endpoint(
        &self,
        value: AtomView<'_>,
        label: &str,
        mode: SchoonschipExpansionMode,
    ) -> Result<ChainEndpoint<Aind>, String> {
        let slot_error = match self.matcher.borrow_mut().parse::<LibraryRep, Aind>(value) {
            Ok(slot) => {
                return Ok(ChainEndpoint {
                    slot,
                    additional_factors: Vec::new(),
                });
            }
            Err(error) => error,
        };
        if let Ok(rep) = self
            .matcher
            .borrow_mut()
            .parse_representation::<LibraryRep>(value)
        {
            return Ok(ChainEndpoint {
                slot: self.slot(&rep),
                additional_factors: Vec::new(),
            });
        }
        if mode.any()
            && let Some(materialized) =
                SchoonschipMaterializer::with_mode(self, mode).materialize_shorthand_arg(value)
        {
            let slot = self
                .matcher
                .borrow_mut()
                .parse::<LibraryRep, Aind>(materialized.current.as_view())
                .map_err(|error| {
                    format!(
                        "invalid materialized chain {label} `{}` from `{value}`: {error}",
                        materialized.current
                    )
                })?;
            return Ok(ChainEndpoint {
                slot,
                additional_factors: materialized.additional_factors,
            });
        }
        Err(format!("invalid chain {label} `{value}`: {slot_error}"))
    }

    /// Unfold one indexed chain using the same placeholder substitution as the
    /// network parser, including supplied rank-one endpoints. Fresh names
    /// remain disjoint across separate conversions.
    /// Unrelated arithmetic and compact spectator ports stay opaque.
    pub fn materialize_indexed_chain(&self, value: FunView<'_>) -> Result<Atom, String> {
        self.materialize_indexed_chain_with(value, &mut |_, _| {})
    }

    /// Unfold a chain and report rebuilt functions while their sources are known.
    /// Admitted leaf ports follow the same placeholder assignments. See
    /// [`Self::materialize_indexed_chain`] for the materialization contract.
    pub fn materialize_indexed_chain_with(
        &self,
        value: FunView<'_>,
        observer: &mut impl FnMut(FunView<'_>, AtomView<'_>),
    ) -> Result<Atom, String> {
        // A projector groups ordered matrix factors, not one matrix with shared
        // endpoints. Expand its normalized permutations before wiring each chain.
        if let Some(expanded) = shadowing::expand_chain_like_projector(value.as_view()) {
            let mut error = None;
            let expression = expanded.replace_map(|atom, _, output| {
                let AtomView::Fun(chain) = atom else { return };
                if chain.get_symbol() != SPENSO_TAG.chain {
                    return;
                }
                match self.materialize_indexed_chain_with(chain, observer) {
                    Ok(indexed) => **output = indexed,
                    Err(reason) => error = Some(reason),
                }
            });
            return error.map_or(Ok(expression), Err);
        }
        if value
            .iter()
            .skip(2)
            .any(|factor| shadowing::projector_factor(factor).is_some())
        {
            return Err("cannot expand normalized chain projector".into());
        }
        let mut arguments = value.iter();
        let start = arguments.next().ok_or("chain requires two endpoints")?;
        let end = arguments.next().ok_or("chain requires two endpoints")?;
        let ChainEndpoint {
            slot: mut left,
            additional_factors: start_factors,
        } = self.materialize_chain_endpoint(start, "start", SchoonschipExpansionMode::full())?;
        let ChainEndpoint {
            slot: end_slot,
            additional_factors: end_factors,
        } = self.materialize_chain_endpoint(end, "end", SchoonschipExpansionMode::full())?;
        // Bare representations remain unresolved public ports when unfolding.
        // Their temporary slots only establish the channel and internal wiring.
        let start_atom = if self
            .matcher
            .borrow_mut()
            .compact_representation(start)
            .is_some()
        {
            start.to_owned()
        } else {
            left.to_atom()
        };
        let end_atom = if self
            .matcher
            .borrow_mut()
            .compact_representation(end)
            .is_some()
        {
            end.to_owned()
        } else {
            end_slot.to_atom()
        };
        let endpoints = start_factors.into_iter().chain(end_factors);
        if !left.rep.matches(&end_slot.rep) {
            return Err("incompatible chain endpoints".into());
        }
        if arguments.len() == 0 {
            let arguments = [start_atom, end_atom];
            let metric = FunctionBuilder::new(ETS.metric)
                .add_args(&arguments)
                .finish();
            self.retain_constructor_ports(
                metric.as_view(),
                ETS.metric,
                &arguments,
                (true, arguments.to_vec()),
            );
            return Ok(Atom::mul_many(endpoints.chain([metric])));
        }
        let count = arguments.len();
        let factors = arguments.enumerate().map(|(position, factor)| {
            let right = if position + 1 == count {
                end_slot.dual()
            } else {
                left.rep.slot(self.fresh_index())
            };
            let factor = self.replace_chain_placeholders(
                factor,
                &if position == 0 {
                    start_atom.clone()
                } else {
                    left.to_atom()
                },
                &if position + 1 == count {
                    end_atom.clone()
                } else {
                    right.dual().to_atom()
                },
                observer,
            );
            left = right;
            factor
        });
        Ok(Atom::mul_many(endpoints.chain(factors)))
    }

    /// Transpose a chain factor, carrying the raw ports of each constructor
    /// across the exchange of the two channel placeholders.
    pub fn transpose_chain_factor(&self, factor: AtomView<'_>) -> Atom {
        self.replace_chain_placeholders(
            factor,
            &Atom::var(SPENSO_TAG.chain_out),
            &Atom::var(SPENSO_TAG.chain_in),
            &mut |_, _| {},
        )
    }

    fn replace_chain_placeholders(
        &self,
        value: AtomView<'_>,
        chain_in: &Atom,
        chain_out: &Atom,
        observer: &mut impl FnMut(FunView<'_>, AtomView<'_>),
    ) -> Atom {
        ChainExpansion::replace_placeholders_with(
            value,
            chain_in,
            chain_out,
            &mut self.matcher.borrow_mut(),
            &mut |source, target, arguments| {
                if let Some(cache) = &self.leaf_interfaces {
                    let retained = cache
                        .lock()
                        .expect("leaf interface cache is not poisoned")
                        .get(source.as_view().get_data())
                        .cloned();
                    if let Some((tensorial, ports)) = retained {
                        let ports = ports
                            .into_iter()
                            .map(|port| match port.as_view() {
                                AtomView::Var(var) if var.get_symbol() == SPENSO_TAG.chain_in => {
                                    chain_in.clone()
                                }
                                AtomView::Var(var) if var.get_symbol() == SPENSO_TAG.chain_out => {
                                    chain_out.clone()
                                }
                                _ => port,
                            })
                            .collect();
                        self.retain_constructor_ports(
                            target,
                            source.get_symbol(),
                            arguments,
                            (tensorial, ports),
                        );
                    }
                }
                observer(source, target);
            },
        )
    }

    /// Open only this compact inner product using the operation's reserved indices.
    /// Callers select the algebraic domain; unrelated arguments are not traversed.
    pub fn materialize_inner_product(&self, value: FunView<'_>) -> Option<Atom> {
        SchoonschipMaterializer::with_mode(self, SchoonschipExpansionMode::none())
            .compact_inner_product(value)
            .map(SchoonschipMaterialization::into_expression)
    }
}

impl<
    'src,
    Sc,
    T: HasStructure + TensorStructure,
    K: Clone + Display + Debug,
    Str: TensorScalarStore<Tensor = T, Scalar = Sc> + Clone,
    Aind: AbsInd + DummyAind + ParseableAind,
> Network<Str, K, Symbol, Aind>
where
    Sc: TryFrom<AtomView<'src>> + TryFrom<Atom> + Clone,
    TensorNetworkError<K, Symbol>:
        From<<Sc as TryFrom<AtomView<'src>>>::Error> + From<<Sc as TryFrom<Atom>>::Error>,
{
    #[allow(clippy::result_large_err)]
    pub(super) fn is_shorthand_function(
        value: FunView<'_>,
        state: &ParseState<Aind, AtomView<'_>>,
    ) -> bool {
        let symbol = value.get_symbol();
        symbol == state.matcher.borrow().tags().chain
            || symbol == state.matcher.borrow().tags().trace
            || symbol == state.matcher.borrow().tags().dot
            || SchoonschipMaterializer::with_mode(state, SchoonschipExpansionMode::full())
                .contains_schoonschip_shorthand(value.as_view())
    }

    #[allow(clippy::result_large_err)]
    pub(super) fn materialize_shorthand<'node, S, Lib, FunLib>(
        construction: &mut Construction<Str, K, Aind, T::Admission>,
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
        let root_chain_disabled = symbol == state.matcher.borrow().tags().chain
            && !settings.shorthand_parsing.expands_chain();
        let root_trace_disabled = symbol == state.matcher.borrow().tags().trace
            && !settings.shorthand_parsing.expands_trace();

        if ((symbol == state.matcher.borrow().tags().chain && !root_chain_disabled)
            || (symbol == state.matcher.borrow().tags().trace && !root_trace_disabled))
            && let Some(expanded) = shadowing::expand_chain_like_projector(value.as_view())
        {
            return Self::try_from_view_impl(
                construction,
                expanded.as_view(),
                state.materialized(),
                library,
                function_library,
                settings,
                |value| AtomOrView::Atom(value.to_owned()),
            );
        }

        if symbol == state.matcher.borrow().tags().chain && !root_chain_disabled {
            return Self::materialize_chain_shorthand(
                construction,
                value,
                state,
                library,
                function_library,
                settings,
            );
        }

        if symbol == state.matcher.borrow().tags().trace && !root_trace_disabled {
            return Self::materialize_trace_shorthand(
                construction,
                value,
                state,
                library,
                function_library,
                settings,
            );
        }

        let has_schoonschip_shorthand =
            SchoonschipMaterializer::with_mode(&state, SchoonschipExpansionMode::full())
                .contains_schoonschip_shorthand(value.as_view());
        let schoonschip_mode = settings
            .shorthand_parsing
            .schoonschip_expansion()
            .unwrap_or_else(SchoonschipExpansionMode::none);
        let effective_schoonschip_mode = if root_chain_disabled || root_trace_disabled {
            schoonschip_mode.for_chain_like_root()
        } else {
            schoonschip_mode
        };
        let materialized = if effective_schoonschip_mode.any() {
            SchoonschipMaterializer::<Aind, _>::with_mode(&state, effective_schoonschip_mode)
                .materialize_shorthand_root(value.as_view())
                .map(SchoonschipMaterialization::into_expression)
        } else {
            None
        };

        if materialized
            .as_ref()
            .is_none_or(|atom| atom.as_view() == value.as_view())
        {
            // The atom rewriter is at a fixed point; recurse only after an actual rewrite.
            if root_chain_disabled || root_trace_disabled || has_schoonschip_shorthand {
                return Self::as_inferred_leaf::<S, Lib, FunLib>(
                    construction,
                    value.as_view(),
                    &state,
                    library,
                    function_library,
                    settings,
                    retain,
                );
            }
            return Self::parse_regular_function_leaf::<S, Lib, FunLib>(
                construction,
                value,
                &state,
                library,
                retain,
            );
        }

        let materialized = materialized.expect("changed shorthand has an owned expression");
        Self::try_from_view_impl(
            construction,
            materialized.as_view(),
            state.materialized(),
            library,
            function_library,
            settings,
            |value| AtomOrView::Atom(value.to_owned()),
        )
    }

    #[allow(clippy::result_large_err)]
    fn materialize_chain_shorthand<S, Lib, FunLib>(
        construction: &mut Construction<Str, K, Aind, T::Admission>,
        value: FunView<'_>,
        state: ParseState<Aind, AtomView<'_>>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
    ) -> Result<NodeIndex, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        let args = value.iter().collect::<Vec<_>>();
        if args.len() < 2 {
            return Err(TensorNetworkError::TooManyArgsFunction(
                value.as_view().to_plain_string(),
            ));
        }

        let schoonschip_mode = settings
            .shorthand_parsing
            .schoonschip_expansion()
            .unwrap_or_else(SchoonschipExpansionMode::none);
        let ChainEndpoint {
            slot: start,
            additional_factors: start_factors,
        } = state
            .materialize_chain_endpoint(args[0], "start", schoonschip_mode)
            .map_err(|error| eyre!(error))?;
        let ChainEndpoint {
            slot: end,
            additional_factors: end_factors,
        } = state
            .materialize_chain_endpoint(args[1], "end", schoonschip_mode)
            .map_err(|error| eyre!(error))?;
        let factors = &args[2..];

        let factor_schoonschip_mode = schoonschip_mode.for_chain_like_root();
        let factor_settings = settings
            .clone()
            .with_schoonschip_expansion(factor_schoonschip_mode);

        let mut factor_networks = Vec::new();
        for factor in start_factors.into_iter().chain(end_factors) {
            factor_networks.extend(Self::parse_chain_like_factor_networks::<S, Lib, FunLib>(
                construction,
                factor,
                state.clone(),
                library,
                function_library,
                settings,
            )?);
        }

        if factors.is_empty() {
            let arguments = [start.to_atom(), end.to_atom()];
            let metric = FunctionBuilder::new(ETS.metric)
                .add_args(&arguments)
                .finish();
            state.retain_constructor_ports(
                metric.as_view(),
                ETS.metric,
                &arguments,
                (true, arguments.to_vec()),
            );
            factor_networks.extend(Self::parse_chain_like_factor_networks::<S, Lib, FunLib>(
                construction,
                metric,
                state,
                library,
                function_library,
                settings,
            )?);
            return Ok(if factor_networks.len() == 1 {
                factor_networks.pop().unwrap()
            } else {
                construction.product(factor_networks)
            });
        }

        let mut left = start;
        for (position, factor) in factors.iter().enumerate() {
            let fresh_right = position + 1 != factors.len();
            let right = if position + 1 == factors.len() {
                end.dual()
            } else {
                state.slot(&left.rep)
            };
            if fresh_right {
                tracing::debug!(
                    target: "spenso::network::parsing",
                    spenso_parser = true,
                    generation = true,
                    compile = true,
                    inspect = true,
                    stage = "chain_shorthand_link_slot",
                    position,
                    factor_count = factors.len(),
                    left = %left,
                    right = %right,
                    right_atom = %right.to_atom().to_plain_string(),
                    start = %start,
                    end = %end,
                    file.factor = %factor.to_plain_string(),
                    "Spenso parser allocated chain shorthand link slot"
                );
            }
            let factor = state.replace_chain_placeholders(
                *factor,
                &left.to_atom(),
                &right.dual().to_atom(),
                &mut |_, _| {},
            );
            let factor = if factor_schoonschip_mode.any() {
                SchoonschipMaterializer::<Aind, _>::with_mode(&state, factor_schoonschip_mode)
                    .materialize_shorthand(factor.as_view())
            } else {
                factor
            };
            factor_networks.extend(Self::parse_chain_like_factor_networks::<S, Lib, FunLib>(
                construction,
                factor,
                state.clone(),
                library,
                function_library,
                &factor_settings,
            )?);
            left = right;
        }

        Ok(if factor_networks.len() == 1 {
            factor_networks.pop().unwrap()
        } else {
            construction.product(factor_networks)
        })
    }

    #[allow(clippy::result_large_err)]
    fn materialize_trace_shorthand<S, Lib, FunLib>(
        construction: &mut Construction<Str, K, Aind, T::Admission>,
        value: FunView<'_>,
        state: ParseState<Aind, AtomView<'_>>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
    ) -> Result<NodeIndex, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        let args = value.iter().collect::<Vec<_>>();
        let Some(rep_view) = args.first() else {
            return Err(TensorNetworkError::TooManyArgsFunction(
                value.as_view().to_plain_string(),
            ));
        };

        let rep = Representation::<LibraryRep>::try_from(*rep_view)
            .map_err(|err| eyre!("invalid trace representation `{rep_view}`: {err}"))?;
        let factors = shadowing::trace_factor_views(&args[1..]);

        if factors.is_empty() {
            return Self::try_from_view_impl(
                construction,
                rep.dim.to_symbolic().as_view(),
                state.materialized(),
                library,
                function_library,
                settings,
                |value| AtomOrView::Atom(value.to_owned()),
            );
        }

        let links = (0..factors.len())
            .map(|position| {
                let slot = state.slot(&rep);
                tracing::debug!(
                    target: "spenso::network::parsing",
                    spenso_parser = true,
                    generation = true,
                    compile = true,
                    inspect = true,
                    stage = "trace_shorthand_link_slot",
                    position,
                    factor_count = factors.len(),
                    slot = %slot,
                    slot_atom = %slot.to_atom().to_plain_string(),
                    representation = %rep,
                    file.factor = %factors[position].to_plain_string(),
                    "Spenso parser allocated trace shorthand link slot"
                );
                slot
            })
            .collect::<Vec<_>>();
        let factor_schoonschip_mode = settings
            .shorthand_parsing
            .schoonschip_expansion()
            .unwrap_or_else(SchoonschipExpansionMode::none)
            .for_chain_like_root();
        let factor_settings = settings
            .clone()
            .with_schoonschip_expansion(factor_schoonschip_mode);

        let materialized_factors = factors
            .iter()
            .enumerate()
            .map(|(position, factor)| {
                let left = links[position].to_atom();
                let right = links[(position + 1) % factors.len()].dual().to_atom();
                let factor =
                    state.replace_chain_placeholders(*factor, &left, &right, &mut |_, _| {});
                if factor_schoonschip_mode.any() {
                    SchoonschipMaterializer::<Aind, _>::with_mode(&state, factor_schoonschip_mode)
                        .materialize_shorthand(factor.as_view())
                } else {
                    factor
                }
            })
            .collect::<Vec<_>>();

        let mut factor_networks = Vec::new();
        for factor in materialized_factors {
            factor_networks.extend(Self::parse_chain_like_factor_networks::<S, Lib, FunLib>(
                construction,
                factor,
                state.clone(),
                library,
                function_library,
                &factor_settings,
            )?);
        }
        Ok(if factor_networks.len() == 1 {
            factor_networks.pop().unwrap()
        } else {
            construction.product(factor_networks)
        })
    }

    #[allow(clippy::result_large_err)]
    fn parse_chain_like_factor_networks<S, Lib, FunLib>(
        construction: &mut Construction<Str, K, Aind, T::Admission>,
        factor: Atom,
        state: ParseState<Aind, AtomView<'_>>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
    ) -> Result<Vec<NodeIndex>, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<'src, S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        let state = state.materialized();
        let AtomView::Mul(product) = factor.as_view() else {
            return Ok(vec![Self::try_from_view_impl(
                construction,
                factor.as_view(),
                state.clone(),
                library,
                function_library,
                settings,
                |value| AtomOrView::Atom(value.to_owned()),
            )?]);
        };

        let mut scalars = Vec::new();
        let mut tensors = Vec::new();
        for arg in product.iter() {
            let network = Self::try_from_view_impl(
                construction,
                arg,
                state.clone(),
                library,
                function_library,
                settings,
                |value| AtomOrView::Atom(value.to_owned()),
            )?;
            if construction.state(network) == NetworkState::PureScalar {
                scalars.push(network);
            } else {
                tensors.push(network);
            }
        }

        if scalars.is_empty() || tensors.len() != 1 {
            return Ok(vec![Self::try_from_view_impl(
                construction,
                factor.as_view(),
                state.clone(),
                library,
                function_library,
                settings,
                |value| AtomOrView::Atom(value.to_owned()),
            )?]);
        }

        scalars.extend(tensors);
        Ok(scalars)
    }
}

/// Utilities for lowering chain and trace shorthand factors.
///
/// Chain parsing chooses the concrete left and right slots for each factor.
/// These helpers rewrite the symbolic placeholders in the factor expression
/// before normal parsing resumes.
pub(super) struct ChainExpansion;

impl ChainExpansion {
    /// Replace `in` and `out` placeholders recursively.
    ///
    /// Nested chain/trace heads are rejected before lowering, so the complete
    /// factor belongs to the current placeholder consumer. This walk only
    /// substitutes its selected endpoints and never allocates dummies.
    /// Explicit slots and compact representations are opaque: their index labels
    /// and dimensions never designate contextual endpoints.
    #[cfg(test)]
    pub(super) fn replace_placeholders(
        value: AtomView<'_>,
        chain_in: &Atom,
        chain_out: &Atom,
        slots: &mut SlotMatcher,
    ) -> Atom {
        Self::replace_placeholders_with(value, chain_in, chain_out, slots, &mut |_, _, _| {})
    }

    fn replace_placeholders_with(
        value: AtomView<'_>,
        chain_in: &Atom,
        chain_out: &Atom,
        slots: &mut SlotMatcher,
        observer: &mut impl FnMut(FunView<'_>, AtomView<'_>, &[Atom]),
    ) -> Atom {
        if !matches!(slots.classify(value), SlotMatch::Other) {
            return value.to_owned();
        }
        match value {
            AtomView::Var(var) if var.get_symbol() == SPENSO_TAG.chain_in => chain_in.clone(),
            AtomView::Var(var) if var.get_symbol() == SPENSO_TAG.chain_out => chain_out.clone(),
            AtomView::Fun(fun) => {
                let arguments = fun
                    .iter()
                    .map(|arg| {
                        Self::replace_placeholders_with(arg, chain_in, chain_out, slots, observer)
                    })
                    .collect::<Vec<_>>();
                // Unaffected scalar and tensor factors retain their normalized
                // value and their raw facts without replaying custom callbacks.
                if fun.iter().eq(arguments.iter().map(AtomCore::as_atom_view)) {
                    return value.to_owned();
                }
                let rebuilt = FunctionBuilder::new(fun.get_symbol())
                    .add_args(&arguments)
                    .finish();
                observer(fun, rebuilt.as_view(), &arguments);
                rebuilt
            }
            AtomView::Add(add) => add.iter().fold(Atom::Zero, |sum, term| {
                sum + Self::replace_placeholders_with(term, chain_in, chain_out, slots, observer)
            }),
            AtomView::Mul(mul) => mul.iter().fold(Atom::num(1), |product, factor| {
                product
                    * Self::replace_placeholders_with(factor, chain_in, chain_out, slots, observer)
            }),
            AtomView::Pow(pow) => {
                let (base, exponent) = pow.get_base_exp();
                Self::replace_placeholders_with(base, chain_in, chain_out, slots, observer)
                    .pow(exponent.to_owned())
            }
            _ => value.to_owned(),
        }
    }
}

#[cfg(test)]
mod tests {
    use symbolica::{atom::Symbol, function, symbol};

    use super::*;
    use crate::{
        broadcast_symbol,
        structure::{
            abstract_index::AbstractIndex,
            representation::{Minkowski, RepName},
        },
    };

    fn mink4() -> Representation<Minkowski> {
        Minkowski {}.new_rep(4)
    }

    fn compact_vector(name: Symbol) -> Atom {
        function!(name, mink4().to_symbolic([]))
    }

    #[test]
    fn index_reservation_copies_admitted_labels_without_callbacks() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let label = symbolica::symbol!(
            "reservation_callback_label",
            tag = &SPENSO_TAG.index,
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let index = function!(label, 13, 1);
        let source = function!(
            SPENSO_TAG.tensor_symbol("reservation_callback_tensor"),
            mink4().to_symbolic([index.clone()])
        );
        let state = ParseState::<AbstractIndex>::default();
        calls.store(0, Ordering::Relaxed);
        state.reserve_indices(source.as_view());
        assert_eq!(calls.load(Ordering::Relaxed), 0);
        assert!(state.reserved_indices.borrow().contains(&index));
    }

    #[test]
    fn compact_inner_product_spectators_are_not_scalar_coefficients() {
        let rep = mink4();
        let left = rep.to_symbolic([Atom::num(91001)]);
        let right = rep.to_symbolic([Atom::num(91002)]);
        let compact = rep.to_symbolic([]);
        let head = SPENSO_TAG.tensor_symbol("compact_inner_product_matrix");
        let metric = function!(
            ETS.metric,
            function!(head, &left, &compact),
            function!(head, &right, &compact)
        );
        let vector =
            compact_vector(SPENSO_TAG.rank_one_tensor_symbol("compact_inner_product_vector"));
        let mut matcher = SlotMatcher::default();
        assert!(!matcher.is_structured_scalar::<AbstractIndex>(metric.as_view()));
        assert!(
            matcher
                .compact_vector_rep::<AbstractIndex>((&metric * &vector).as_view())
                .is_none()
        );
        let projected = function!(*shadowing::SYM, &metric, &vector);
        assert!(
            matcher
                .compact_vector_rep::<AbstractIndex>(projected.as_view())
                .is_none()
        );
        let scalar_dot = function!(SPENSO_TAG.dot, &vector, &vector);
        assert!(matcher.is_structured_scalar::<AbstractIndex>(scalar_dot.as_view()));
        assert!(
            matcher
                .compact_vector_rep::<AbstractIndex>((&scalar_dot * &vector).as_view())
                .is_some()
        );
    }

    #[test]
    fn scalar_inner_product_coefficients_normalize_before_materialization() {
        let p = compact_vector(SPENSO_TAG.rank_one_tensor_symbol("selected_inner_p"));
        let q = compact_vector(SPENSO_TAG.rank_one_tensor_symbol("selected_inner_q"));
        let r = compact_vector(SPENSO_TAG.rank_one_tensor_symbol("selected_inner_r"));
        let coefficient = function!(SPENSO_TAG.dot, &p, &q);
        let inner = function!(SPENSO_TAG.dot, &p, &r);
        let value = function!(SPENSO_TAG.dot, &coefficient * &p, &r);
        assert!(coefficient.is_scalar().is_true());
        assert_eq!(value, &coefficient * &inner);
        let state = ParseState::<AbstractIndex>::default();
        let AtomView::Fun(function) = inner.as_view() else {
            unreachable!()
        };
        let opened = &coefficient * state.materialize_inner_product(function).unwrap();
        let AtomView::Mul(product) = opened.as_view() else {
            panic!("weighted inner product remains factored")
        };
        assert_eq!(
            product
                .iter()
                .filter(|factor| *factor == coefficient.as_view())
                .count(),
            1
        );
        assert_eq!(state.next_dummy.get(), 1_000_001);
    }

    #[test]
    fn compact_inner_product_refuses_ambiguous_or_incompatible_axes() {
        let head = SPENSO_TAG.tensor_symbol("compact_inner_product_ambiguous");
        let compact = mink4().to_symbolic([]);
        let second = Minkowski {}.new_rep(3).to_symbolic([]);
        let one = function!(head, &compact);
        let two = function!(head, &compact, &second);
        let wrong = function!(head, &second);
        let mut matcher = SlotMatcher::default();
        for input in [
            function!(SPENSO_TAG.dot, &one, &two),
            function!(SPENSO_TAG.dot, &one, &wrong),
        ] {
            let AtomView::Fun(function) = input.as_view() else {
                unreachable!()
            };
            assert!(
                matcher
                    .compact_inner_product_parts::<AbstractIndex>(function)
                    .is_none()
            );
        }
    }

    #[test]
    fn compact_inner_product_observation_does_not_construct_callbacks() {
        use std::sync::{Arc, Mutex};
        let calls = Arc::new(Mutex::new(Vec::new()));
        let recorded = calls.clone();
        let head = crate::tensor_symbol!(
            "compact_inner_product_observed_callback",
            norm = move |node, _| recorded.lock().unwrap().push(node.to_owned())
        );
        let rep = mink4();
        let compact = rep.to_symbolic([]);
        let input = function!(
            ETS.metric,
            function!(head, rep.to_symbolic([Atom::num(91101)]), &compact),
            function!(head, rep.to_symbolic([Atom::num(91102)]), &compact)
        );
        calls.lock().unwrap().clear();
        let mut matcher = SlotMatcher::default();
        let structure = OrderedStructure::<LibraryRep, AbstractIndex>::structure_from_atom(
            input.as_view(),
            &mut matcher,
        )
        .unwrap();
        assert_eq!(structure.canonical().order(), 2);
        assert!(calls.lock().unwrap().is_empty());
        let state = ParseState::<AbstractIndex>::default();
        state.reserve_indices(input.as_view());
        let AtomView::Fun(function) = input.as_view() else {
            unreachable!()
        };
        let opened = state.materialize_inner_product(function).unwrap();
        assert_eq!(calls.lock().unwrap().len(), 2);
        let observed = OrderedStructure::<LibraryRep, AbstractIndex>::structure_from_atom(
            opened.as_view(),
            &mut matcher,
        )
        .unwrap();
        assert_eq!(structure, observed);
        assert_eq!(calls.lock().unwrap().len(), 2);
    }

    #[test]
    fn standalone_compact_vector_is_not_materialized() {
        let state = ParseState::<AbstractIndex>::default();
        let materializer =
            SchoonschipMaterializer::with_mode(&state, SchoonschipExpansionMode::full());

        let expression = compact_vector(SPENSO_TAG.tensor_symbol("materialized_p"));

        assert_eq!(
            materializer.materialize_shorthand(expression.as_view()),
            expression
        );
    }

    #[test]
    fn compact_vector_argument_becomes_slot_and_factor() {
        let state = ParseState::<AbstractIndex>::default();
        let materializer =
            SchoonschipMaterializer::with_mode(&state, SchoonschipExpansionMode::full());
        let visible_slot = mink4()
            .slot::<AbstractIndex, _>(AbstractIndex::from(1))
            .to_atom();
        let expression = FunctionBuilder::new(symbol!("f"))
            .add_arg(compact_vector(SPENSO_TAG.tensor_symbol("materialized_p_argument")).as_view())
            .add_arg(visible_slot.as_view())
            .finish();

        let materialized = materializer.materialize_shorthand(expression.as_view());
        let AtomView::Mul(product) = materialized.as_view() else {
            panic!("expected product");
        };
        let factors = product.iter().collect::<Vec<_>>();

        assert_eq!(factors.len(), 2);
        let f = factors
            .iter()
            .find_map(|factor| match factor {
                AtomView::Fun(fun) if fun.get_symbol() == symbol!("f") => Some(fun),
                _ => None,
            })
            .unwrap();
        let p = factors
            .iter()
            .find_map(|factor| match factor {
                AtomView::Fun(fun)
                    if fun.get_symbol() == SPENSO_TAG.tensor_symbol("materialized_p_argument") =>
                {
                    Some(fun)
                }
                _ => None,
            })
            .unwrap();

        let f_args = f.iter().collect::<Vec<_>>();
        let p_args = p.iter().collect::<Vec<_>>();
        let f_dummy = Slot::<LibraryRep, AbstractIndex>::try_from(f_args[0]).unwrap();
        let p_dummy = Slot::<LibraryRep, AbstractIndex>::try_from(p_args[0]).unwrap();

        assert_eq!(f_dummy, p_dummy);
        assert_eq!(f_args[1], visible_slot.as_view());
    }

    #[test]
    fn compact_scalar_product_still_materializes_to_two_factors() {
        let state = ParseState::<AbstractIndex>::default();
        let materializer =
            SchoonschipMaterializer::with_mode(&state, SchoonschipExpansionMode::full());
        let expression = function!(
            ETS.metric,
            compact_vector(SPENSO_TAG.tensor_symbol("materialized_metric_p")),
            compact_vector(SPENSO_TAG.tensor_symbol("materialized_metric_q"))
        );

        let materialized = materializer.materialize_shorthand(expression.as_view());
        let AtomView::Mul(product) = materialized.as_view() else {
            panic!("expected product");
        };

        assert_eq!(product.iter().count(), 2);
    }

    #[test]
    fn scalar_metadata_keeps_compact_products_sums_and_powers_opaque() {
        let state = ParseState::<AbstractIndex>::default();
        let materializer =
            SchoonschipMaterializer::with_mode(&state, SchoonschipExpansionMode::full());
        let scalar = symbol!("materialized_scalar_metadata"; Scalar);
        let wrapper = symbol!("materialized_metadata_wrapper");
        let tensor = SPENSO_TAG.tensor_symbol("materialized_metadata_owner");
        let p = compact_vector(SPENSO_TAG.rank_one_tensor_symbol("materialized_metadata_p"));
        let q = compact_vector(SPENSO_TAG.rank_one_tensor_symbol("materialized_metadata_q"));
        let dot = function!(ETS.metric, &p, &q);
        let slot = mink4()
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(31))
            .to_atom();

        for payload in [p.clone(), &p + &q, dot.clone(), dot.pow(2)] {
            let metadata = function!(scalar, payload);
            for expression in [
                metadata.clone(),
                function!(tensor, &metadata, &slot),
                function!(tensor, function!(wrapper, &metadata), &slot),
            ] {
                assert!(!materializer.contains_schoonschip_shorthand(expression.as_view()));
                assert_eq!(
                    materializer.materialize_shorthand(expression.as_view()),
                    expression
                );
            }
        }
        assert_eq!(state.next_dummy.get(), 1_000_000);
    }

    #[test]
    fn scalar_attribute_takes_priority_over_compact_tensor_tags() {
        let state = ParseState::<AbstractIndex>::default();
        let materializer =
            SchoonschipMaterializer::with_mode(&state, SchoonschipExpansionMode::full());
        let scalar_tensor = symbol!("materialized_scalar_tensor"; Scalar; tags = [
            SPENSO_TAG.tensor.clone(), SPENSO_TAG.rank1.clone()
        ]);
        let scalar_broadcast = symbol!("materialized_scalar_broadcast"; Scalar; tags = [
            SPENSO_TAG.broadcast.clone()
        ]);
        let tensor = SPENSO_TAG.tensor_symbol("materialized_conflict_owner");
        let compact = function!(scalar_tensor, mink4().to_symbolic([]));
        let broadcast = function!(
            scalar_broadcast,
            compact_vector(SPENSO_TAG.rank_one_tensor_symbol("materialized_conflict_p"))
        );

        for metadata in [compact, broadcast] {
            let AtomView::Fun(function) = metadata.as_view() else {
                unreachable!()
            };
            assert!(
                state
                    .matcher
                    .borrow_mut()
                    .compact_tensor_rep_arg::<AbstractIndex>(function)
                    .is_none()
            );
            assert!(
                state
                    .matcher
                    .borrow_mut()
                    .compact_vector_rep::<AbstractIndex>(metadata.as_view())
                    .is_none()
            );
            for expression in [metadata.clone(), function!(tensor, &metadata)] {
                assert!(!materializer.contains_schoonschip_shorthand(expression.as_view()));
                assert_eq!(
                    materializer.materialize_shorthand(expression.as_view()),
                    expression
                );
            }
        }
        assert_eq!(state.next_dummy.get(), 1_000_000);
    }

    #[test]
    fn bound_vector_lowering_preserves_adjacent_scalar_metadata() {
        let state = ParseState::<AbstractIndex>::default();
        let materializer =
            SchoonschipMaterializer::with_mode(&state, SchoonschipExpansionMode::full());
        let p = SPENSO_TAG.rank_one_tensor_symbol("materialized_bound_p");
        let q = SPENSO_TAG.rank_one_tensor_symbol("materialized_bound_q");
        let tensor = SPENSO_TAG.tensor_symbol("materialized_bound_owner");
        let scalar = symbol!("materialized_bound_metadata"; Scalar);
        let metadata = function!(
            scalar,
            function!(ETS.metric, compact_vector(p), compact_vector(q))
        );
        let visible = mink4()
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(37))
            .to_atom();
        let expression = function!(tensor, &metadata, compact_vector(p), &visible);
        assert!(materializer.contains_schoonschip_shorthand(expression.as_view()));

        let materialized = materializer.materialize_shorthand(expression.as_view());
        let AtomView::Mul(product) = materialized.as_view() else {
            panic!("expected bound-vector factors")
        };
        assert_eq!(product.iter().count(), 2);
        let owner = product
            .iter()
            .find_map(|factor| match factor {
                AtomView::Fun(function) if function.get_symbol() == tensor => Some(function),
                _ => None,
            })
            .unwrap();
        let vector = product
            .iter()
            .find_map(|factor| match factor {
                AtomView::Fun(function) if function.get_symbol() == p => Some(function),
                _ => None,
            })
            .unwrap();
        let arguments = owner.iter().collect::<Vec<_>>();
        assert_eq!(arguments[0], metadata.as_view());
        assert_eq!(arguments[2], visible.as_view());
        assert_eq!(arguments[1], vector.iter().next().unwrap());
        assert!(Slot::<LibraryRep, AbstractIndex>::try_from(arguments[1]).is_ok());
        assert_eq!(state.next_dummy.get(), 1_000_001);
    }

    #[test]
    fn untagged_representation_metadata_is_not_a_compact_vector() {
        let state = ParseState::<AbstractIndex>::default();
        let materializer =
            SchoonschipMaterializer::with_mode(&state, SchoonschipExpansionMode::full());
        let metadata = FunctionBuilder::new(symbol!("scalar_rep_metadata"))
            .add_arg(mink4().to_symbolic([]))
            .finish();
        let expression = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("metadata_tensor"))
            .add_arg(metadata)
            .add_arg(mink4().to_symbolic([]))
            .finish();

        assert_eq!(
            materializer.materialize_shorthand(expression.as_view()),
            expression
        );
    }

    #[test]
    fn placeholder_replacement_keeps_explicit_slot_payloads_opaque() {
        let input = mink4().to_symbolic([Atom::num(47)]);
        let output = mink4().to_symbolic([Atom::num(53)]);
        let head = SPENSO_TAG.tensor_symbol("placeholder_opaque_tensor");
        for index in [
            AbstractIndex::from(SPENSO_TAG.chain_in),
            AbstractIndex::Normal(59).scoped(SPENSO_TAG.chain_out),
        ] {
            let spectator = mink4().to_symbolic([index.to_atom()]);
            let factor = function!(head, SPENSO_TAG.chain_in, SPENSO_TAG.chain_out, &spectator);
            let expected = function!(head, &input, &output, &spectator);
            assert_eq!(
                ChainExpansion::replace_placeholders(
                    factor.as_view(),
                    &input,
                    &output,
                    &mut SlotMatcher::default(),
                ),
                expected
            );
        }
    }

    #[test]
    fn placeholder_replacement_descends_through_factor_wrappers() {
        let representation = mink4();
        let factor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("nested_scope_factor"))
            .add_arg(Atom::var(SPENSO_TAG.chain_in))
            .add_arg(Atom::var(SPENSO_TAG.chain_out))
            .finish();
        let wrapped = FunctionBuilder::new(broadcast_symbol!("factor_scope_broadcast"))
            .add_arg(factor)
            .finish();
        let outer_input = representation
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(47))
            .to_atom();
        let outer_output = representation
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(53))
            .to_atom();
        let expected_factor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("nested_scope_factor"))
            .add_arg(&outer_input)
            .add_arg(&outer_output)
            .finish();
        let expected = FunctionBuilder::new(broadcast_symbol!("factor_scope_broadcast"))
            .add_arg(expected_factor)
            .finish();

        assert_eq!(
            ChainExpansion::replace_placeholders(
                wrapped.as_view(),
                &outer_input,
                &outer_output,
                &mut SlotMatcher::default(),
            ),
            expected
        );
    }

    #[test]
    fn placeholder_replacement_descends_through_factor_products() {
        let representation = mink4();
        let factor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("wrapped_scope_factor"))
            .add_arg(Atom::var(SPENSO_TAG.chain_in))
            .add_arg(Atom::var(SPENSO_TAG.chain_out))
            .finish();
        let product = Atom::num(3) * factor;
        let outer_input = representation
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(47))
            .to_atom();
        let outer_output = representation
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(53))
            .to_atom();
        let expected_factor =
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol("wrapped_scope_factor"))
                .add_arg(&outer_input)
                .add_arg(&outer_output)
                .finish();
        let expected = Atom::num(3) * expected_factor;

        assert_eq!(
            ChainExpansion::replace_placeholders(
                product.as_view(),
                &outer_input,
                &outer_output,
                &mut SlotMatcher::default(),
            ),
            expected
        );
    }
}
