//! Materialization for shorthand syntax before ordinary network parsing.
//!
//! This pass is deliberately narrower than algebraic simplification. It takes
//! compact parser syntax that cannot be parsed as an ordinary leaf yet and
//! rewrites it into explicit slots plus ordinary tensor factors. The parser can
//! then recurse on the resulting expression and build the same graph it would
//! have built from fully expanded syntax.
//!
//! The main Schoonschip convention is:
//! 1. a compact rank-one tensor `p(rep)` used as a function argument becomes a
//!    fresh slot in that argument position;
//! 2. the tensor `p(slot)` is multiplied next to the rebuilt function;
//! 3. compact scalar products `g(p(rep), q(rep))` and `dot(p(rep), q(rep))`
//!    share one fresh self-dual slot and become the product `p(slot) * q(slot)`.
//! Additional factors are accumulated beside the current atom and are not
//! recursively inspected by this materialization pass.
//!
//! Chain and trace expansion is separate. It only replaces the symbolic
//! `in`/`out` placeholders with the slots chosen by the chain parser, then lets
//! the Schoonschip materializer handle any compact arguments inside each factor.

use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, representation::FunView},
    id::MatchSettings,
};

use super::{ParseSettings, ParseState, ShorthandParsing};
use crate::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    structure::{
        representation::{LibraryRep, RepName, Representation},
        slot::{AbsInd, DummyAind, IsAbstractSlot, ParseableAind, Slot},
    },
};

struct SchoonschipMaterialization {
    /// Expression that remains at the original syntactic position.
    current: Atom,
    /// Extra factors to multiply beside `current`.
    ///
    /// These factors are already explicit tensor syntax, so the materializer
    /// does not run its shorthand pattern search on them during this pass.
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
pub(super) struct ShorthandMaterializer<'a, Aind> {
    state: &'a ParseState<Aind>,
    settings: &'a ParseSettings,
}

impl<'a, Aind: AbsInd + DummyAind + ParseableAind> ShorthandMaterializer<'a, Aind> {
    /// Build a materializer for the current parse recursion.
    pub(super) fn new(state: &'a ParseState<Aind>, settings: &'a ParseSettings) -> Self {
        Self { state, settings }
    }

    /// Materialize an expandable shorthand expression.
    ///
    /// The root must be a function. A successful materialization returns one
    /// expression where the rebuilt root is multiplied by all tensor factors
    /// introduced while replacing compact vector arguments.
    pub(super) fn materialize_shorthand(&self, value: AtomView<'_>) -> Option<Atom> {
        self.materialize_shorthand_root(value)
            .map(SchoonschipMaterialization::into_expression)
    }

    /// Materialize shorthand at a function root.
    ///
    /// Expand mode first gives compact scalar products their special lowering.
    /// Everything else is rebuilt by recursively scanning function arguments.
    fn materialize_shorthand_root(
        &self,
        value: AtomView<'_>,
    ) -> Option<SchoonschipMaterialization> {
        let AtomView::Fun(fun) = value else {
            return None;
        };

        if self.settings.shorthand_parsing.expands()
            && let Some(materialized) = self.compact_scalar_product(fun)
        {
            return Some(materialized);
        }

        self.materialize_shorthand_function(fun)
    }

    /// Materialize one function argument according to shorthand position rules.
    ///
    /// A compact vector has special meaning only in argument position under
    /// `ShorthandParsing::Expand`: it contributes a fresh slot at that position
    /// and an extra rank-one tensor factor. Non-compact functions are scanned
    /// recursively.
    fn materialize_shorthand_arg(&self, value: AtomView<'_>) -> Option<SchoonschipMaterialization> {
        if self.settings.shorthand_parsing == ShorthandParsing::Expand
            && let Some(materialized) = self.materialize_compact_vector_arg(value)
        {
            return Some(materialized);
        }

        self.materialize_shorthand_root(value)
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
        if value.get_symbol() != ETS.metric && value.get_symbol() != SPENSO_TAG.dot {
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

        let slot = self.state.slot(&rep).to_atom();
        let lhs = Self::materialize_compact_vector_with_slot(*lhs, &rep, &slot)?;
        let rhs = Self::materialize_compact_vector_with_slot(*rhs, &rep, &slot)?;

        Some(SchoonschipMaterialization {
            current: Atom::num(1),
            additional_factors: vec![lhs, rhs],
        })
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
        let settings = MatchSettings {
            level_range: (0, Some(0)),
            partial: false,
            ..Default::default()
        };
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

    /// Multiply the first actual chain/trace factor by a scalar.
    ///
    /// For open chains the first factor starts after the two endpoint slots. For
    /// traces it starts after the representation argument. Empty chains/traces
    /// have no factor to scale and are left to their dedicated parser cases.
    pub(super) fn scale_first_factor(value: AtomView<'_>, scalar: &Atom) -> Option<Atom> {
        let AtomView::Fun(fun) = value else {
            return None;
        };

        let offset = if fun.get_symbol() == SPENSO_TAG.chain {
            2
        } else if fun.get_symbol() == SPENSO_TAG.trace {
            1
        } else {
            return None;
        };

        if fun.get_nargs() <= offset {
            return None;
        }

        let mut rebuilt = FunctionBuilder::new(fun.get_symbol());
        for (position, arg) in fun.iter().enumerate() {
            if position == offset {
                rebuilt = rebuilt.add_arg(scalar.clone() * arg.to_owned());
            } else {
                rebuilt = rebuilt.add_arg(arg);
            }
        }
        Some(rebuilt.finish())
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
        let settings = ParseSettings::default();
        let state = ParseState::<AbstractIndex>::default();
        let materializer = ShorthandMaterializer::new(&state, &settings);

        assert!(
            materializer
                .materialize_shorthand(compact_vector(symbol!("p")).as_view())
                .is_none()
        );
    }

    #[test]
    fn compact_vector_argument_becomes_slot_and_factor() {
        let settings = ParseSettings::default();
        let state = ParseState::<AbstractIndex>::default();
        let materializer = ShorthandMaterializer::new(&state, &settings);
        let visible_slot = mink4()
            .slot::<AbstractIndex, _>(AbstractIndex::from(1))
            .to_atom();
        let expression = FunctionBuilder::new(symbol!("f"))
            .add_arg(compact_vector(symbol!("p")).as_view())
            .add_arg(visible_slot.as_view())
            .finish();

        let materialized = materializer
            .materialize_shorthand(expression.as_view())
            .unwrap();
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
                AtomView::Fun(fun) if fun.get_symbol() == symbol!("p") => Some(fun),
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
        let settings = ParseSettings::default();
        let state = ParseState::<AbstractIndex>::default();
        let materializer = ShorthandMaterializer::new(&state, &settings);
        let expression = function!(
            ETS.metric,
            compact_vector(symbol!("p")),
            compact_vector(symbol!("q"))
        );

        let materialized = materializer
            .materialize_shorthand(expression.as_view())
            .unwrap();
        let AtomView::Mul(product) = materialized.as_view() else {
            panic!("expected product");
        };

        assert_eq!(product.iter().count(), 2);
    }
}
