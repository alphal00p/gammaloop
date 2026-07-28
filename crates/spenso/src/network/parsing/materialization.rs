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
//! 1. a compact rank-one tensor `p(rep)` used as a function argument becomes a
//!    fresh slot in that argument position;
//! 2. the tensor `p(slot)` is multiplied next to the rebuilt function;
//! 3. compact scalar products `g(p(rep), q(rep))` and `dot(p(rep), q(rep))`
//!    share one fresh self-dual slot and become the product `p(slot) * q(slot)`.
//!
//! Additional factors are accumulated beside the current atom and are not
//! recursively inspected by this materialization pass.
//!
//! Chain and trace materialization chooses the symbolic `in`/`out` replacements,
//! then lets this Schoonschip helper expand compact arguments inside each factor.

use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol, representation::FunView},
    id::MatchSettings,
};

use std::fmt::{Debug, Display};

use super::{
    ParseSettings, ParseState, SchoonschipExpansionMode, StructureFromAtom, TensorFromExpression,
    TensorLibraryFor,
};
use crate::{
    network::{
        Network, NetworkState, TensorNetworkError,
        graph::NMul,
        library::{FunctionLibrary, symbolic::ETS},
        store::TensorScalarStore,
        tags::SPENSO_TAG,
    },
    shadowing,
    shadowing::Concretize,
    structure::{
        HasStructure, ScalarStructure, TensorShell, TensorStructure,
        representation::{LibraryRep, RepName, Representation},
        slot::{AbsInd, DualSlotTo, DummyAind, IsAbstractSlot, ParseableAind, Slot},
    },
};
use eyre::eyre;

struct SchoonschipMaterialization {
    /// Expression that remains at the original syntactic position.
    current: Atom,
    /// Extra factors to multiply beside `current`.
    ///
    /// These factors are already explicit tensor syntax, so the materializer
    /// does not run its shorthand pattern search on them during this pass.
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
pub(super) struct SchoonschipMaterializer<'a, Aind> {
    state: &'a ParseState<Aind>,
    mode: SchoonschipExpansionMode,
}

impl<'a, Aind: AbsInd + DummyAind + ParseableAind> SchoonschipMaterializer<'a, Aind> {
    /// Build a materializer with explicit Schoonschip expansion controls.
    pub(super) fn with_mode(state: &'a ParseState<Aind>, mode: SchoonschipExpansionMode) -> Self {
        Self { state, mode }
    }

    /// Return true when this function root contains compact Schoonschip syntax.
    ///
    /// This is the non-allocating counterpart to `materialize_shorthand`: compact
    /// scalar products are shorthand at the root, while compact vectors become
    /// shorthand only when they occur as arguments of another function.
    pub(super) fn contains_schoonschip_shorthand(value: AtomView<'_>) -> bool {
        let AtomView::Fun(fun) = value else {
            return false;
        };

        Self::is_compact_scalar_product(fun)
            || fun.iter().any(Self::contains_schoonschip_shorthand_arg)
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

        if Self::is_chain_like_head(fun) && !self.mode.expand_inside_chains {
            return None;
        }

        if Self::is_inner_product_head(fun) {
            if !self.mode.inner_products {
                return None;
            }
            if let Some(materialized) = self.compact_scalar_product(fun) {
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

    fn contains_schoonschip_shorthand_arg(value: AtomView<'_>) -> bool {
        Self::compact_vector_rep(value).is_some() || Self::contains_schoonschip_shorthand(value)
    }

    /// Materialize one compact vector argument as `slot` plus `vector(slot)`.
    ///
    /// The compact vector may be a single function or a sum of functions, but
    /// every visible compact representation must be the same representation.
    fn materialize_compact_vector_arg(
        &self,
        value: AtomView<'_>,
    ) -> Option<SchoonschipMaterialization> {
        let rep = Self::compact_vector_rep(value)?;
        let slot = self.state.slot(&rep).to_atom();
        let factor = Self::materialize_compact_vector_with_slot(value, &rep, &slot)?;
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
        let mut rebuilt = FunctionBuilder::new(fun.get_symbol());

        for arg in fun.iter() {
            if let Some(materialized) = self.materialize_shorthand_arg(arg) {
                changed = true;
                rebuilt = rebuilt.add_arg(&materialized.current);
                additional_factors.extend(materialized.additional_factors);
            } else {
                rebuilt = rebuilt.add_arg(arg);
            }
        }

        changed.then(|| {
            let replacement = rebuilt.finish();
            SchoonschipMaterialization {
                current: Self::compact_dot_as_metric(replacement.as_view()).unwrap_or(replacement),
                additional_factors,
            }
        })
    }

    /// Materialize a compact metric or dot product into two tensor factors.
    ///
    /// Both arguments must be compact vectors with the same self-dual
    /// representation. They are assigned the same fresh slot, so
    /// `g(p(rep), q(rep))` becomes `p(slot) * q(slot)`.
    fn compact_scalar_product(&self, value: FunView<'_>) -> Option<SchoonschipMaterialization> {
        let (lhs, rhs, rep) = Self::compact_scalar_product_parts(value)?;

        let slot = self.state.slot(&rep).to_atom();
        let lhs = Self::materialize_compact_vector_with_slot(lhs, &rep, &slot)?;
        let rhs = Self::materialize_compact_vector_with_slot(rhs, &rep, &slot)?;
        tracing::debug!(
            target: "spenso::network::parsing",
            spenso_parser = true,
            generation = true,
            compile = true,
            inspect = true,
            stage = "schoonschip_compact_scalar_product_slot",
            slot = %slot.to_plain_string(),
            representation = %rep,
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

    fn is_compact_scalar_product(value: FunView<'_>) -> bool {
        Self::compact_scalar_product_parts(value).is_some()
    }

    fn is_inner_product_head(value: FunView<'_>) -> bool {
        let symbol = value.get_symbol();
        symbol == ETS.metric || symbol == SPENSO_TAG.dot
    }

    fn is_chain_like_head(value: FunView<'_>) -> bool {
        let symbol = value.get_symbol();
        symbol == SPENSO_TAG.chain || symbol == SPENSO_TAG.trace
    }

    fn compact_scalar_product_parts(
        value: FunView<'_>,
    ) -> Option<(AtomView<'_>, AtomView<'_>, Representation<LibraryRep>)> {
        if !Self::is_inner_product_head(value) {
            return None;
        }

        let args = value.iter().collect::<Vec<_>>();
        let [lhs, rhs] = args.as_slice() else {
            return None;
        };

        let rep = Self::compact_vector_rep(*lhs)?;
        if rep != Self::compact_vector_rep(*rhs)? || !rep.rep.is_self_dual() {
            return None;
        }

        Some((*lhs, *rhs, rep))
    }

    /// Infer the compact representation carried by a rank-one shorthand atom.
    ///
    /// Functions expose a compact representation through exactly one direct
    /// representation argument. Sums are accepted only when every summand exposes
    /// the same representation.
    fn compact_vector_rep(value: AtomView<'_>) -> Option<Representation<LibraryRep>> {
        match value {
            AtomView::Fun(fun) => Self::compact_tensor_rep_arg(fun).map(|(_, rep)| rep),
            AtomView::Add(add) => {
                let mut reps = add.iter().map(Self::compact_vector_rep);
                let rep = reps.next()??;
                reps.all(|candidate| candidate == Some(rep)).then_some(rep)
            }
            _ => None,
        }
    }

    /// Locate the compact representation argument of one tensor function.
    ///
    /// A compact vector function is not itself a representation, is not a metric
    /// or dot product, has no explicit slot argument, and has exactly one direct
    /// argument matching the representation wildcard convention.
    fn compact_tensor_rep_arg(value: FunView<'_>) -> Option<(usize, Representation<LibraryRep>)> {
        if value.get_symbol() == ETS.metric || value.get_symbol() == SPENSO_TAG.dot {
            return None;
        }

        if Representation::<LibraryRep>::try_from(value.as_view()).is_ok() {
            return None;
        }

        let args = value.iter().collect::<Vec<_>>();
        if args
            .iter()
            .any(|arg| Slot::<LibraryRep, Aind>::try_from(*arg).is_ok())
        {
            return None;
        }

        let rep_args = args
            .iter()
            .enumerate()
            .filter_map(|(position, arg)| {
                Self::compact_rep_pattern_match(*arg).map(|rep| (position, rep))
            })
            .collect::<Vec<_>>();

        let [(position, rep)] = rep_args.as_slice() else {
            return None;
        };

        Some((*position, *rep))
    }

    /// Match one argument against the symbolic representation wildcard.
    ///
    /// Matching is restricted to the argument itself (`level_range = 0`) so a
    /// nested representation inside metadata does not accidentally become the
    /// tensor's compact slot.
    fn compact_rep_pattern_match(arg: AtomView<'_>) -> Option<Representation<LibraryRep>> {
        let rep_pattern = Atom::var(SPENSO_TAG.rep_).to_pattern();
        let settings = MatchSettings::new()
            .level_range((0, Some(0)))
            .partial(false);
        let mut matches = arg.pattern_match(&rep_pattern, None, Some(&settings));
        let matched = matches.next_detailed()?;
        let rep = rep_pattern.replace_wildcards_with_matches(matched.match_stack);
        Representation::<LibraryRep>::try_from(rep.as_view()).ok()
    }

    /// Replace a compact representation with a concrete slot.
    ///
    /// For a function, this rebuilds the function with the matched compact
    /// representation argument replaced by `slot`. For a sum, every summand is
    /// rebuilt with the same slot so the expansion keeps a single dummy edge.
    fn materialize_compact_vector_with_slot(
        value: AtomView<'_>,
        rep: &Representation<LibraryRep>,
        slot: &Atom,
    ) -> Option<Atom> {
        match value {
            AtomView::Fun(fun) => {
                let (position, matched_rep) = Self::compact_tensor_rep_arg(fun)?;
                if matched_rep != *rep {
                    return None;
                }

                let mut tensor = FunctionBuilder::new(fun.get_symbol());
                for (arg_position, arg) in fun.iter().enumerate() {
                    if arg_position == position {
                        tensor = tensor.add_arg(slot);
                    } else {
                        tensor = tensor.add_arg(arg);
                    }
                }
                Some(tensor.finish())
            }
            AtomView::Add(add) => {
                let mut terms = add
                    .iter()
                    .map(|term| Self::materialize_compact_vector_with_slot(term, rep, slot));
                let first = terms.next()??;
                let rest = terms.collect::<Option<Vec<_>>>()?;
                Some(rest.into_iter().fold(first, |sum, term| sum + term))
            }
            _ => None,
        }
    }

    /// Normalize `dot(slot_i, slot_j)` to `g(slot_i, slot_j)`.
    ///
    /// This only applies after arguments have already been materialized to
    /// concrete slots. Non-slot dot products are left untouched.
    fn compact_dot_as_metric(value: AtomView<'_>) -> Option<Atom> {
        let AtomView::Fun(fun) = value else {
            return None;
        };

        if fun.get_symbol() != SPENSO_TAG.dot || fun.get_nargs() != 2 {
            return None;
        }

        let args = fun.iter().collect::<Vec<_>>();
        if args
            .iter()
            .all(|arg| Slot::<LibraryRep, Aind>::try_from(*arg).is_ok())
        {
            Some(
                FunctionBuilder::new(ETS.metric)
                    .add_arg(args[0])
                    .add_arg(args[1])
                    .finish(),
            )
        } else {
            None
        }
    }
}

impl<
    'a,
    Sc,
    T: HasStructure + TensorStructure,
    K: Clone + Display + Debug,
    Str: TensorScalarStore<Tensor = T, Scalar = Sc> + Clone,
    Aind: AbsInd + DummyAind + ParseableAind,
> Network<Str, K, Symbol, Aind>
where
    Sc: for<'r> TryFrom<AtomView<'r>> + Clone,
    TensorNetworkError<K, Symbol>: for<'r> From<<Sc as TryFrom<AtomView<'r>>>::Error>,
{
    #[allow(clippy::result_large_err)]
    pub(super) fn is_shorthand_function(value: FunView<'a>) -> bool {
        let symbol = value.get_symbol();
        symbol == SPENSO_TAG.chain
            || symbol == SPENSO_TAG.trace
            || symbol == SPENSO_TAG.dot
            || SchoonschipMaterializer::<Aind>::contains_schoonschip_shorthand(value.as_view())
    }

    #[allow(clippy::result_large_err)]
    pub(super) fn materialize_shorthand<S, Lib, FunLib>(
        value: FunView<'a>,
        state: ParseState<Aind>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        let symbol = value.get_symbol();
        let root_chain_disabled =
            symbol == SPENSO_TAG.chain && !settings.shorthand_parsing.expands_chain();
        let root_trace_disabled =
            symbol == SPENSO_TAG.trace && !settings.shorthand_parsing.expands_trace();

        if symbol == SPENSO_TAG.chain && !root_chain_disabled {
            return Self::materialize_chain_shorthand(
                value,
                state,
                library,
                function_library,
                settings,
            );
        }

        if symbol == SPENSO_TAG.trace && !root_trace_disabled {
            return Self::materialize_trace_shorthand(
                value,
                state,
                library,
                function_library,
                settings,
            );
        }

        let has_schoonschip_shorthand =
            SchoonschipMaterializer::<Aind>::contains_schoonschip_shorthand(value.as_view());
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
            SchoonschipMaterializer::<Aind>::with_mode(&state, effective_schoonschip_mode)
                .materialize_shorthand(value.as_view())
        } else {
            value.as_view().to_owned()
        };

        if materialized == value.as_view().to_owned() {
            // The atom rewriter is at a fixed point; recurse only after an actual rewrite.
            if root_chain_disabled || root_trace_disabled || has_schoonschip_shorthand {
                return Self::as_inferred_leaf::<S, Lib, FunLib>(
                    value.as_view(),
                    super::StructureInferenceMode::Fast,
                    library,
                    function_library,
                    settings,
                );
            }
            return Self::parse_regular_function_leaf::<S, Lib>(value, library);
        }

        Self::try_from_view_impl(
            materialized.as_view(),
            state,
            library,
            function_library,
            settings,
        )
    }

    fn materialize_chain_endpoint(
        value: AtomView<'a>,
        label: &str,
        state: &ParseState<Aind>,
        schoonschip_mode: SchoonschipExpansionMode,
    ) -> Result<ChainEndpoint<Aind>, TensorNetworkError<K, Symbol>> {
        match Slot::<LibraryRep, Aind>::try_from(value) {
            Ok(slot) => Ok(ChainEndpoint {
                slot,
                additional_factors: Vec::new(),
            }),
            Err(slot_err) => {
                if schoonschip_mode.any()
                    && let Some(materialized) =
                        SchoonschipMaterializer::<Aind>::with_mode(state, schoonschip_mode)
                            .materialize_shorthand_arg(value)
                {
                    let slot =
                        match Slot::<LibraryRep, Aind>::try_from(materialized.current.as_view()) {
                            Ok(slot) => slot,
                            Err(err) => {
                                return Err(eyre!(
                                    "invalid materialized chain {label} `{}` from `{}`: {err}",
                                    materialized.current,
                                    value
                                )
                                .into());
                            }
                        };
                    return Ok(ChainEndpoint {
                        slot,
                        additional_factors: materialized.additional_factors,
                    });
                }

                Err(eyre!("invalid chain {label} `{}`: {slot_err}", value).into())
            }
        }
    }

    #[allow(clippy::result_large_err)]
    fn materialize_chain_shorthand<S, Lib, FunLib>(
        value: FunView<'a>,
        state: ParseState<Aind>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<S, Sc, K, Symbol, Aind, Lib, FunLib>,
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
        } = Self::materialize_chain_endpoint(args[0], "start", &state, schoonschip_mode)?;
        let ChainEndpoint {
            slot: end,
            additional_factors: end_factors,
        } = Self::materialize_chain_endpoint(args[1], "end", &state, schoonschip_mode)?;
        let factors = &args[2..];

        let factor_schoonschip_mode = schoonschip_mode.for_chain_like_root();
        let factor_settings = settings
            .clone()
            .with_schoonschip_expansion(factor_schoonschip_mode);

        let mut factor_networks = Vec::new();
        for factor in start_factors.into_iter().chain(end_factors) {
            factor_networks.extend(Self::parse_chain_like_factor_networks::<S, Lib, FunLib>(
                factor,
                state.clone(),
                library,
                function_library,
                settings,
            )?);
        }

        if factors.is_empty() {
            let metric = FunctionBuilder::new(ETS.metric)
                .add_arg(start.to_atom())
                .add_arg(end.to_atom())
                .finish();
            factor_networks.extend(Self::parse_chain_like_factor_networks::<S, Lib, FunLib>(
                metric,
                state,
                library,
                function_library,
                settings,
            )?);
            return Ok(Self::chain_like_network_product(factor_networks));
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
            let factor = ChainExpansion::replace_placeholders(
                *factor,
                &left.to_atom(),
                &right.dual().to_atom(),
            );
            let factor = if factor_schoonschip_mode.any() {
                SchoonschipMaterializer::<Aind>::with_mode(&state, factor_schoonschip_mode)
                    .materialize_shorthand(factor.as_view())
            } else {
                factor
            };
            factor_networks.extend(Self::parse_chain_like_factor_networks::<S, Lib, FunLib>(
                factor,
                state.clone(),
                library,
                function_library,
                &factor_settings,
            )?);
            left = right;
        }

        Ok(Self::chain_like_network_product(factor_networks))
    }

    #[allow(clippy::result_large_err)]
    fn materialize_trace_shorthand<S, Lib, FunLib>(
        value: FunView<'a>,
        state: ParseState<Aind>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<S, Sc, K, Symbol, Aind, Lib, FunLib>,
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
                rep.dim.to_symbolic().as_view(),
                state,
                library,
                function_library,
                settings,
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
                let factor = ChainExpansion::replace_placeholders(*factor, &left, &right);
                if factor_schoonschip_mode.any() {
                    SchoonschipMaterializer::<Aind>::with_mode(&state, factor_schoonschip_mode)
                        .materialize_shorthand(factor.as_view())
                } else {
                    factor
                }
            })
            .collect::<Vec<_>>();

        let mut factor_networks = Vec::new();
        for factor in materialized_factors {
            factor_networks.extend(Self::parse_chain_like_factor_networks::<S, Lib, FunLib>(
                factor,
                state.clone(),
                library,
                function_library,
                &factor_settings,
            )?);
        }
        Ok(Self::chain_like_network_product(factor_networks))
    }

    #[allow(clippy::result_large_err)]
    fn parse_chain_like_factor_networks<S, Lib, FunLib>(
        factor: Atom,
        state: ParseState<Aind>,
        library: &Lib,
        function_library: &FunLib,
        settings: &ParseSettings,
    ) -> Result<Vec<Self>, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + ScalarStructure + Clone + StructureFromAtom,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        T: TensorFromExpression<S, Sc, K, Symbol, Aind, Lib, FunLib>,
        Lib: TensorLibraryFor<S, T, Key = K>,
        FunLib: FunctionLibrary<T, Sc, Key = Symbol>,
    {
        let AtomView::Mul(product) = factor.as_view() else {
            return Ok(vec![Self::try_from_view_impl(
                factor.as_view(),
                state.clone(),
                library,
                function_library,
                settings,
            )?]);
        };

        let mut scalars = Vec::new();
        let mut tensors = Vec::new();
        for arg in product.iter() {
            let network =
                Self::try_from_view_impl(arg, state.clone(), library, function_library, settings)?;
            if network.state == NetworkState::PureScalar {
                scalars.push(network);
            } else {
                tensors.push(network);
            }
        }

        if scalars.is_empty() || tensors.len() != 1 {
            return Ok(vec![Self::try_from_view_impl(
                factor.as_view(),
                state.clone(),
                library,
                function_library,
                settings,
            )?]);
        }

        scalars.extend(tensors);
        Ok(scalars)
    }

    fn chain_like_network_product(networks: Vec<Self>) -> Self {
        let mut networks = networks.into_iter();
        let first = networks.next().unwrap();
        let rest = networks.collect::<Vec<_>>();
        if rest.is_empty() {
            first
        } else {
            first.n_mul(rest)
        }
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
    /// The recursion is purely syntactic over functions, sums, products, and
    /// powers. It does not allocate dummies or inspect tensor structure.
    pub(super) fn replace_placeholders(
        value: AtomView<'_>,
        chain_in: &Atom,
        chain_out: &Atom,
    ) -> Atom {
        match value {
            AtomView::Var(var) if var.get_symbol() == SPENSO_TAG.chain_in => chain_in.clone(),
            AtomView::Var(var) if var.get_symbol() == SPENSO_TAG.chain_out => chain_out.clone(),
            AtomView::Fun(fun) => {
                let mut rebuilt = FunctionBuilder::new(fun.get_symbol());
                for arg in fun.iter() {
                    rebuilt = rebuilt.add_arg(Self::replace_placeholders(arg, chain_in, chain_out));
                }
                rebuilt.finish()
            }
            AtomView::Add(add) => add.iter().fold(Atom::Zero, |sum, term| {
                sum + Self::replace_placeholders(term, chain_in, chain_out)
            }),
            AtomView::Mul(mul) => mul.iter().fold(Atom::num(1), |product, factor| {
                product * Self::replace_placeholders(factor, chain_in, chain_out)
            }),
            AtomView::Pow(pow) => {
                let (base, exponent) = pow.get_base_exp();
                Self::replace_placeholders(base, chain_in, chain_out).pow(exponent.to_owned())
            }
            _ => value.to_owned(),
        }
    }
}

#[cfg(test)]
mod tests {
    use symbolica::{atom::Symbol, function, symbol};

    use super::*;
    use crate::structure::{
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName},
    };

    fn mink4() -> Representation<Minkowski> {
        Minkowski {}.new_rep(4)
    }

    fn compact_vector(name: Symbol) -> Atom {
        function!(name, mink4().to_symbolic([]))
    }

    #[test]
    fn standalone_compact_vector_is_not_materialized() {
        let state = ParseState::<AbstractIndex>::default();
        let materializer =
            SchoonschipMaterializer::with_mode(&state, SchoonschipExpansionMode::full());

        let expression = compact_vector(symbol!("materialized_p"));

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
            .add_arg(compact_vector(symbol!("materialized_p")).as_view())
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
                AtomView::Fun(fun) if fun.get_symbol() == symbol!("materialized_p") => Some(fun),
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
            compact_vector(symbol!("materialized_p")),
            compact_vector(symbol!("materialized_q"))
        );

        let materialized = materializer.materialize_shorthand(expression.as_view());
        let AtomView::Mul(product) = materialized.as_view() else {
            panic!("expected product");
        };

        assert_eq!(product.iter().count(), 2);
    }
}
