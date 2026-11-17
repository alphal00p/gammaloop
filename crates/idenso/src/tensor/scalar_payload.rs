//! Checked transformations of explicitly selected, independently closed payloads.

use spenso::network::library::function_lib::INBUILTS;
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::abstract_index::AIND_SYMBOLS;
use spenso::structure::partial::{PartialStructure, PartialStructureExt};
use spenso::structure::slot::{SlotMatch, SlotMatcher};
use std::result::Result as StdResult;
use std::sync::{Arc, OnceLock};
use symbolica::atom::representation::FunView;
use symbolica::atom::{Atom, AtomView, FunctionBuilder, PowView, Symbol};

use super::inference::{InterfaceInference, TensorInferenceError};
use super::{SymbolicTensor, composition};

type Result<T> = StdResult<T, TensorInferenceError>;
type Tensor = SymbolicTensor<PartialStructure>;

/// An original syntactic position whose value may be an independent scalar.
///
/// Selecting a position does not certify that it is scalar: the mapper checks
/// its native interface before and after the transformation. Function argument
/// numbers are zero-based. Tensor, conjugation and word operands are opaque;
/// a bracket is offered as one ordered, checked expression.
/// See the [mapping example](SymbolicTensor::map_closed_scalar_payloads).
#[derive(Clone, Copy, Debug)]
pub enum ScalarPayloadPosition<'a> {
    Root {
        expression: AtomView<'a>,
    },
    FunctionArgument {
        function: FunView<'a>,
        argument: usize,
    },
    PowerBase {
        power: PowView<'a>,
    },
    PowerExponent {
        power: PowView<'a>,
    },
    ClosedBracket {
        function: FunView<'a>,
    },
}

impl Tensor {
    /// Map caller-selected mathematical scalars without opening tensor scopes.
    ///
    /// Selection observes original positions. Selected nested payloads are
    /// transformed bottom-up, once per position; callback-generated syntax is
    /// never traversed. Unselected function arguments, power operands, brackets,
    /// native tensor/word operands and structural metadata remain opaque.
    /// Unchanged ancestors retain their original syntax and do not replay
    /// normalization callbacks. Changed expressions use native constructors.
    ///
    /// `Root` is offered once for this call, with this carrier's established
    /// interface. In particular, selecting a ranked zero fails rather than
    /// treating its expression `0` as an untyped scalar. Every selected value
    /// must be independently closed before and after `map`; the final result
    /// must retain the original outer logical interface. Scalar here describes
    /// only the interface: a closed native word retains its word algebra and
    /// is not certified to commute or to admit ordinary scalar series.
    ///
    /// # Errors
    /// Rejects unpublished inputs, selected open payloads, invalid callback
    /// results, and changes to the original outer logical interface.
    ///
    /// # Example
    /// ```
    /// use idenso::tensor::{ScalarPayloadPosition, SymbolicTensor};
    /// use idenso::tensor::inference::TensorInferenceError;
    /// use symbolica::atom::{Atom, AtomCore};
    /// use symbolica::{function, symbol};
    ///
    /// let f = symbol!("docs::scalar_payload_f");
    /// let source = SymbolicTensor::infer(function!(f, Atom::num(2)))?;
    /// let mapped = source.map_closed_scalar_payloads(
    ///     |position| matches!(position, ScalarPayloadPosition::FunctionArgument { .. }),
    ///     |_, value| value.with_checked_expression(value.expression() + Atom::one()),
    /// )?;
    /// assert_eq!(mapped.expression(), &function!(f, Atom::num(3)));
    /// # Ok::<(), TensorInferenceError>(())
    /// ```
    pub fn map_closed_scalar_payloads(
        &self,
        mut select: impl FnMut(ScalarPayloadPosition<'_>) -> bool,
        map: impl FnMut(ScalarPayloadPosition<'_>, &Self) -> Result<Self>,
    ) -> Result<Self> {
        self.explicit_index_occurrences()?;
        let position = ScalarPayloadPosition::Root {
            expression: self.expression().as_view(),
        };
        let root_selected = select(position);
        if root_selected && !self.is_scalar() {
            return Err(TensorInferenceError::invalid(
                "a selected scalar payload must have a closed scalar interface",
            ));
        }
        let mut mapper = ClosedPayloadMapper { select, map };
        let expression = mapper.expression(self.expression().as_view())?;
        let result = self.borrowed_parts().with_checked_payload(expression)?;
        if root_selected {
            mapper.apply(position, &result)
        } else {
            Ok(result)
        }
    }
}

struct ClosedPayloadMapper<S, M> {
    select: S,
    map: M,
}

impl<S, M> ClosedPayloadMapper<S, M>
where
    S: FnMut(ScalarPayloadPosition<'_>) -> bool,
    M: FnMut(ScalarPayloadPosition<'_>, &Tensor) -> Result<Tensor>,
{
    fn require_scalar(value: &Tensor) -> Result<()> {
        value.explicit_index_occurrences()?;
        if !value.is_scalar() {
            return Err(TensorInferenceError::invalid(
                "a selected scalar payload must have a closed scalar interface",
            ));
        }
        Ok(())
    }

    fn apply(&mut self, position: ScalarPayloadPosition<'_>, value: &Tensor) -> Result<Tensor> {
        Self::require_scalar(value)?;
        let mut result = (self.map)(position, value)?;
        Self::require_scalar(&result)?;
        // A callback can return a fresh carrier or cancel written indices.
        // Retain both admitted histories, without exporting nested closed
        // payloads into their enclosing expression's contraction scope.
        let source_history = value.explicit_index_occurrences()?;
        let result_history = result.explicit_index_occurrences()?;
        if !source_history.is_bounded_by(result_history) {
            let mut history = result_history.clone();
            history.branch_max(source_history);
            history.validate()?;
            result.proofs = result.proofs.algebra_result();
            result.proofs.occurrences = OnceLock::from(Arc::new(history));
        }
        Ok(result)
    }

    fn payload(
        &mut self,
        position: ScalarPayloadPosition<'_>,
        original: AtomView<'_>,
    ) -> Result<Atom> {
        // Inference also normalizes ordered products. Validate the original
        // closed syntax instead: a selected bracket is the whole wrapper, not
        // its inferred commutative product.
        let value = Tensor::from_declared_parts(
            original.to_owned(),
            PartialStructure::from_logical_slots([]),
        )?;
        Self::require_scalar(&value)?;
        let expression = match position {
            // Enter ordered factors only after the whole original bracket has
            // proved closed. Its wrapper is not offered a second time.
            ScalarPayloadPosition::ClosedBracket { function } => self
                .function(function, |mapper, _position, argument| {
                    mapper.expression(argument)
                })?,
            _ => self.expression(original)?,
        };
        let value = if expression == *value.expression() {
            value
        } else {
            value.borrowed_parts().with_checked_payload(expression)?
        };
        Ok(self.apply(position, &value)?.into_expression())
    }

    fn expression(&mut self, original: AtomView<'_>) -> Result<Atom> {
        match original {
            AtomView::Add(terms) => {
                let values = terms
                    .iter()
                    .map(|term| self.expression(term))
                    .collect::<Result<Vec<_>>>()?;
                if terms.iter().zip(&values).all(|(a, b)| a == b.as_view()) {
                    Ok(original.to_owned())
                } else {
                    Ok(Atom::add_many(values))
                }
            }
            AtomView::Mul(factors) => {
                let values = factors
                    .iter()
                    .map(|factor| self.expression(factor))
                    .collect::<Result<Vec<_>>>()?;
                if factors.iter().zip(&values).all(|(a, b)| a == b.as_view()) {
                    Ok(original.to_owned())
                } else {
                    Ok(Atom::mul_many(values))
                }
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let base_position = ScalarPayloadPosition::PowerBase { power };
                let exponent_position = ScalarPayloadPosition::PowerExponent { power };
                let map_base = (self.select)(base_position);
                let map_exponent = (self.select)(exponent_position);
                let base_result = if map_base {
                    self.payload(base_position, base)?
                } else {
                    base.to_owned()
                };
                let exponent_result = if map_exponent {
                    self.payload(exponent_position, exponent)?
                } else {
                    exponent.to_owned()
                };
                if base_result.as_view() == base && exponent_result.as_view() == exponent {
                    Ok(original.to_owned())
                } else {
                    // Typed powers protect independently closed contractions
                    // with their existing native scalar-power wrapper.
                    let base = Tensor::infer(base_result)?;
                    let exponent = Tensor::infer(exponent_result)?;
                    Ok(base.pow(&exponent)?.into_expression())
                }
            }
            AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.bracket => {
                let position = ScalarPayloadPosition::ClosedBracket { function };
                if (self.select)(position) {
                    self.payload(position, original)
                } else {
                    Ok(original.to_owned())
                }
            }
            AtomView::Fun(function) if !Self::opaque_function(function) => {
                self.function(function, |mapper, argument, value| {
                    let position = ScalarPayloadPosition::FunctionArgument { function, argument };
                    if (mapper.select)(position) {
                        mapper.payload(position, value)
                    } else {
                        Ok(value.to_owned())
                    }
                })
            }
            _ => Ok(original.to_owned()),
        }
    }

    fn function(
        &mut self,
        function: FunView<'_>,
        mut argument: impl FnMut(&mut Self, usize, AtomView<'_>) -> Result<Atom>,
    ) -> Result<Atom> {
        let values = function
            .iter()
            .enumerate()
            .map(|(position, value)| argument(self, position, value))
            .collect::<Result<Vec<_>>>()?;
        if function.iter().zip(&values).all(|(a, b)| a == b.as_view()) {
            Ok(AtomView::Fun(function).to_owned())
        } else {
            Ok(FunctionBuilder::new(function.get_symbol())
                .add_args(values)
                .finish())
        }
    }

    fn opaque_function(function: FunView<'_>) -> bool {
        let symbol = function.get_symbol();
        let indices = &*AIND_SYMBOLS;
        InterfaceInference::has_structured_syntax(AtomView::Fun(function))
            || composition::is_composite_head(symbol, &SPENSO_TAG)
            || symbol == Symbol::CONJ
            || symbol == INBUILTS.conj
            || symbol == SPENSO_TAG.scalar
            || symbol == SPENSO_TAG.pure_scalar
            || symbol.has_tag(&SPENSO_TAG.index)
            || [
                indices.aind,
                indices.uind,
                indices.dind,
                indices.selfdualind,
                indices.openind,
                indices.cind,
                indices.find,
            ]
            .contains(&symbol)
            || !matches!(
                SlotMatcher::default().classify(AtomView::Fun(function)),
                SlotMatch::Other
            )
    }
}

#[cfg(test)]
mod tests {
    use std::slice;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};

    use spenso::network::library::function_lib::INBUILTS;
    use spenso::network::tags::SPENSO_TAG;
    use spenso::structure::abstract_index::AbstractIndex;
    use spenso::structure::partial::{PartialStructure, PartialStructureExt};
    use spenso::structure::representation::{Minkowski, RepName};
    use spenso::structure::slot::IsAbstractSlot;
    use spenso::{bracket, mink, tensor_symbol};
    use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
    use symbolica::{function, symbol};

    use super::{ScalarPayloadPosition as Position, Tensor};
    use crate::tensor::composition::TensorCompositionError;
    use crate::test_support::test_initialize;

    fn vector(head: Symbol, label: usize) -> Atom {
        function!(
            head,
            Minkowski {}
                .new_rep(4)
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(label))
                .to_atom()
        )
    }

    fn expression(atom: Atom) -> Tensor {
        Tensor::infer(atom).unwrap()
    }

    #[test]
    fn root_mapping_retains_hidden_occurrences_and_unchanged_proofs() {
        test_initialize();
        let a = expression(vector(tensor_symbol!("closed_payload::history_a"), 4701));
        let b = expression(vector(tensor_symbol!("closed_payload::history_b"), 4701));
        let pair = a.multiply(&b).unwrap();
        let zero = pair.multiply(&expression(Atom::Zero)).unwrap();
        let unchanged = zero
            .map_closed_scalar_payloads(
                |position| matches!(position, Position::Root { .. }),
                |_, value| Ok(value.clone()),
            )
            .unwrap();
        assert!(unchanged.same_reduction_state(&zero));
        for source in [&pair, &zero] {
            for replacement in [Atom::one(), Atom::Zero] {
                let mapped = source
                    .map_closed_scalar_payloads(
                        |position| matches!(position, Position::Root { .. }),
                        |_, _| Ok(expression(replacement.clone())),
                    )
                    .unwrap();
                assert_eq!(mapped.expression(), &replacement);
                assert!(matches!(
                    mapped.multiply(&a),
                    Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
                ));
            }
        }

        // Callback-created hidden history also survives when the source had none.
        let mapped = expression(Atom::one())
            .map_closed_scalar_payloads(
                |position| matches!(position, Position::Root { .. }),
                |_, _| Ok(zero.clone()),
            )
            .unwrap();
        assert!(mapped.same_reduction_state(&zero));
        assert!(mapped.multiply(&a).is_err());

        let c = expression(vector(tensor_symbol!("closed_payload::history_c"), 4703));
        let d = expression(vector(tensor_symbol!("closed_payload::history_d"), 4703));
        let other_zero = c
            .multiply(&d)
            .unwrap()
            .multiply(&expression(Atom::Zero))
            .unwrap();
        let mapped = zero
            .map_closed_scalar_payloads(
                |position| matches!(position, Position::Root { .. }),
                |_, _| Ok(other_zero.clone()),
            )
            .unwrap();
        assert!(mapped.multiply(&a).is_err());
        assert!(mapped.multiply(&c).is_err());
    }

    #[test]
    fn payload_mapping_rejects_unpublished_inputs_and_callback_results() {
        test_initialize();
        let source = expression(Atom::one());
        let mut unpublished = source.clone();
        unpublished.invalidate_proofs();
        assert!(
            unpublished
                .map_closed_scalar_payloads(|_| false, |_, value| Ok(value.clone()))
                .is_err()
        );
        assert!(
            unpublished
                .scalar_algebra_input(symbol!("closed_payload::unpublished_t"; Scalar), |_| false)
                .is_err()
        );
        assert!(
            source
                .map_closed_scalar_payloads(
                    |position| matches!(position, Position::Root { .. }),
                    |_, _| Ok(unpublished.clone()),
                )
                .is_err()
        );
    }

    #[test]
    fn changed_bracket_payload_keeps_ordered_syntax_until_its_callback() {
        test_initialize();
        let f = symbol!("closed_payload::changed_bracket_f");
        let a = vector(tensor_symbol!("closed_payload::changed_bracket_a"), 4702);
        let b = vector(tensor_symbol!("closed_payload::changed_bracket_b"), 4702);
        let source = Tensor::from_declared_parts(
            bracket!(function!(f, 2), &a, &b),
            PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        let mut visited_bracket = false;
        let result = source
            .map_closed_scalar_payloads(
                |position| {
                    matches!(
                        position,
                        Position::ClosedBracket { .. } | Position::FunctionArgument { .. }
                    )
                },
                |position, value| match position {
                    Position::FunctionArgument { .. } => {
                        value.with_checked_expression(Atom::num(3))
                    }
                    Position::ClosedBracket { .. } => {
                        visited_bracket = true;
                        assert_eq!(value.expression(), &bracket!(function!(f, 3), &a, &b));
                        Ok(value.clone())
                    }
                    _ => unreachable!(),
                },
            )
            .unwrap();
        assert!(visited_bracket);
        assert_eq!(result.expression(), &bracket!(function!(f, 3), a, b));
    }

    #[test]
    fn original_positions_are_bottom_up_root_once_and_outputs_are_not_visited() {
        test_initialize();
        let (f, g, h, x) = symbol!(
            "closed_payload::f",
            "closed_payload::g",
            "closed_payload::h",
            "closed_payload::x"
        );
        let source = expression(function!(f, function!(g, x)));
        let mut positions = Vec::new();
        let result = source
            .map_closed_scalar_payloads(
                |_| true,
                |position, value| {
                    positions.push(match position {
                        Position::Root { expression } => {
                            assert_eq!(expression, source.expression().as_view());
                            "root"
                        }
                        Position::FunctionArgument { function, argument } => {
                            assert_eq!(argument, 0);
                            if function.get_symbol() == g {
                                "g"
                            } else {
                                assert_eq!(function.get_symbol(), f);
                                assert_eq!(
                                    function.iter().next().unwrap(),
                                    function!(g, x).as_view()
                                );
                                "f"
                            }
                        }
                        _ => panic!("unexpected position"),
                    });
                    value.with_checked_expression(function!(h, value.expression()))
                },
            )
            .unwrap();
        assert_eq!(positions, ["g", "f", "root"]);
        assert_eq!(
            result.expression(),
            &function!(h, function!(f, function!(h, function!(g, function!(h, x)))))
        );
    }

    #[test]
    fn unselected_arguments_and_native_tensor_word_and_structural_payloads_are_opaque() {
        test_initialize();
        let (f, g, x) = symbol!(
            "closed_payload::opaque_f",
            "closed_payload::opaque_g",
            "closed_payload::opaque_x"
        );
        let source = expression(function!(f, function!(g, x)));
        let mut calls = 0;
        let result = source
            .map_closed_scalar_payloads(
                |position| {
                    if let Position::FunctionArgument { function, .. } = position {
                        calls += 1;
                        assert_eq!(function.get_symbol(), f);
                    }
                    false
                },
                |_, _| panic!("unselected callback"),
            )
            .unwrap();
        assert_eq!(calls, 1);
        assert_eq!(result.expression(), source.expression());

        let head = tensor_symbol!("closed_payload::opaque_tensor");
        let tensor = expression(function!(head, function!(f, x), mink!(4, 4101)));
        let word_head = symbol!("closed_payload::opaque_word"; Flat, Linear;
            tag = &SPENSO_TAG.tensor_word);
        let word = Tensor::word(word_head, slice::from_ref(&tensor)).unwrap();
        for value in [
            tensor,
            word,
            expression(function!(INBUILTS.conj, vector(head, 4102))),
        ] {
            let mut callback_count = 0;
            let result = value
                .map_closed_scalar_payloads(
                    |position| matches!(position, Position::FunctionArgument { .. }),
                    |_, value| {
                        callback_count += 1;
                        Ok(value.clone())
                    },
                )
                .unwrap();
            assert_eq!(callback_count, 0);
            assert_eq!(result.expression(), value.expression());
            assert_eq!(result.structure(), value.structure());
        }
    }

    #[test]
    fn independently_closed_arguments_reuse_labels_without_outer_capture() {
        test_initialize();
        let (a, b, c) = (
            tensor_symbol!("closed_payload::a"),
            tensor_symbol!("closed_payload::b"),
            tensor_symbol!("closed_payload::c"),
        );
        let f = symbol!("closed_payload::two_args");
        let closed = vector(a, 4201) * vector(b, 4201);
        let source = expression(function!(f, &closed, &closed) * vector(c, 4201));
        let scope = symbol!("closed_payload::argument_scope");
        let mut seen = Vec::new();
        let result = source
            .map_closed_scalar_payloads(
                |position| matches!(position, Position::FunctionArgument { .. }),
                |position, value| {
                    let Position::FunctionArgument { argument, .. } = position else {
                        unreachable!()
                    };
                    seen.push(argument);
                    value.wrap_indices(scope, true)
                },
            )
            .unwrap();
        assert_eq!(seen, [0, 1]);
        assert_eq!(result.structure(), source.structure());
        let scoped = expression(closed).wrap_indices(scope, true).unwrap();
        assert_eq!(
            result.expression(),
            &(function!(f, scoped.expression(), scoped.expression()) * vector(c, 4201))
        );
    }

    #[test]
    fn selected_open_values_and_rank_changing_callbacks_fail() {
        test_initialize();
        let head = tensor_symbol!("closed_payload::open");
        let open = expression(vector(head, 4301));
        assert!(
            open.map_closed_scalar_payloads(
                |position| matches!(position, Position::Root { .. }),
                |_, value| Ok(value.clone()),
            )
            .is_err()
        );
        let zero = Tensor::from_declared_parts(Atom::Zero, open.structure().clone()).unwrap();
        assert!(
            zero.map_closed_scalar_payloads(
                |position| matches!(position, Position::Root { .. }),
                |_, value| Ok(value.clone()),
            )
            .is_err()
        );
        assert_eq!(
            zero.map_closed_scalar_payloads(|_| false, |_, _| unreachable!())
                .unwrap()
                .structure(),
            open.structure()
        );
        let scalar = expression(function!(symbol!("closed_payload::scalar"), 1));
        assert!(
            scalar
                .map_closed_scalar_payloads(
                    |position| matches!(position, Position::FunctionArgument { .. }),
                    |_, _| Ok(open.clone()),
                )
                .is_err()
        );
        let function = expression(function!(
            symbol!("closed_payload::open_argument"),
            open.expression()
        ));
        assert!(
            function
                .map_closed_scalar_payloads(
                    |position| matches!(position, Position::FunctionArgument { .. }),
                    |_, value| Ok(value.clone()),
                )
                .is_err()
        );
        let power = open.pow(&expression(Atom::num(2))).unwrap();
        assert!(
            power
                .map_closed_scalar_payloads(
                    |position| matches!(position, Position::PowerBase { .. }),
                    |_, value| Ok(value.clone()),
                )
                .is_err()
        );
    }

    #[test]
    fn scalar_metadata_is_unselected_and_selected_words_keep_their_native_algebra() {
        test_initialize();
        let (f, x) = symbol!(
            "closed_payload::metadata_body",
            "closed_payload::metadata_x"
        );
        let metadata = symbol!("closed_payload::metadata"; Scalar);
        let source = expression(function!(metadata, function!(f, x)));
        let mut selected_heads = Vec::new();
        let unchanged = source
            .map_closed_scalar_payloads(
                |position| {
                    if let Position::FunctionArgument { function, .. } = position {
                        selected_heads.push(function.get_symbol());
                    }
                    false
                },
                |_, _| unreachable!(),
            )
            .unwrap();
        assert_eq!(selected_heads, [metadata]);
        assert_eq!(unchanged.expression(), source.expression());

        let head = symbol!("closed_payload::closed_word"; Flat, Linear, Antisymmetric;
            tag = &SPENSO_TAG.tensor_word);
        let a = expression(function!(tensor_symbol!("closed_payload::scalar_odd_a")));
        let b = expression(function!(tensor_symbol!("closed_payload::scalar_odd_b")));
        for operands in [vec![a.clone()], vec![a, b]] {
            let word = Tensor::word(head, &operands).unwrap();
            assert!(word.is_scalar());
            let unchanged = word
                .map_closed_scalar_payloads(
                    |position| matches!(position, Position::Root { .. }),
                    |_, value| {
                        assert!(value.expression().contains_symbol(head));
                        Ok(value.clone())
                    },
                )
                .unwrap();
            assert_eq!(unchanged.expression(), word.expression());
            let wrapped = expression(function!(f, word.expression()));
            let unchanged = wrapped
                .map_closed_scalar_payloads(
                    |position| matches!(position, Position::FunctionArgument { .. }),
                    |_, value| {
                        assert!(value.expression().contains_symbol(head));
                        Ok(value.clone())
                    },
                )
                .unwrap();
            assert_eq!(unchanged.expression(), wrapped.expression());
        }
    }

    #[test]
    fn closed_bracket_is_one_ordered_payload_and_nested_selections_are_distinct() {
        test_initialize();
        let (a, b) = (
            tensor_symbol!("closed_payload::bracket_a"),
            tensor_symbol!("closed_payload::bracket_b"),
        );
        let f = symbol!("closed_payload::bracket_f");
        let original = bracket!(function!(f, 2), vector(a, 4401), vector(b, 4401));
        let source =
            Tensor::from_declared_parts(original.clone(), PartialStructure::from_logical_slots([]))
                .unwrap();
        let mut order = Vec::new();
        let result = source.map_closed_scalar_payloads(
            |position| matches!(position, Position::ClosedBracket { .. } | Position::FunctionArgument { .. }),
            |position, value| {
                match position {
                    Position::FunctionArgument { .. } => order.push("argument"),
                    Position::ClosedBracket { function } => {
                        assert_eq!(function.get_nargs(), 3);
                        order.push("bracket");
                        assert!(matches!(value.expression().as_view(), AtomView::Fun(f) if f.get_symbol()==SPENSO_TAG.bracket));
                    }
                    _ => unreachable!(),
                }
                Ok(value.clone())
            },
        ).unwrap();
        assert_eq!(order, ["argument", "bracket"]);
        assert_eq!(result.expression(), &original);
        let factor = expression(vector(a, 4402));
        let open =
            Tensor::from_declared_parts(bracket!(factor.expression()), factor.structure().clone())
                .unwrap();
        assert!(
            open.map_closed_scalar_payloads(
                |position| matches!(position, Position::ClosedBracket { .. }),
                |_, value| Ok(value.clone()),
            )
            .is_err()
        );
    }

    #[test]
    fn typed_power_reconstruction_preserves_closed_base_and_original_descriptors() {
        test_initialize();
        let (a, b) = (
            tensor_symbol!("closed_payload::power_a"),
            tensor_symbol!("closed_payload::power_b"),
        );
        let scalar = expression(vector(a, 4501) * vector(b, 4501));
        let source = scalar.pow(&expression(Atom::num(2))).unwrap();
        let mut callbacks = Vec::new();
        let result = source
            .map_closed_scalar_payloads(
                |position| {
                    matches!(
                        position,
                        Position::PowerBase { .. } | Position::PowerExponent { .. }
                    )
                },
                |position, value| match position {
                    Position::PowerExponent { power } => {
                        assert_eq!(power.get_base_exp().1, Atom::num(2).as_view());
                        callbacks.push("exponent");
                        Ok(expression(Atom::num(-1)))
                    }
                    Position::PowerBase { .. } => {
                        callbacks.push("base");
                        Ok(value.clone())
                    }
                    _ => unreachable!(),
                },
            )
            .unwrap();
        assert_eq!(callbacks, ["base", "exponent"]);
        assert_eq!(
            result.expression(),
            scalar.pow(&expression(Atom::num(-1))).unwrap().expression()
        );
        result.explicit_index_occurrences().unwrap();
    }

    #[test]
    fn unchanged_nodes_do_not_replay_normalizers_and_changed_nodes_use_them() {
        test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = calls.clone();
        let f = symbol!(
            "closed_payload::observed",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::SeqCst);
            }
        );
        let x = symbol!("closed_payload::observed_x");
        let source = expression(function!(f, x));
        let before = calls.load(Ordering::SeqCst);
        let result = source
            .map_closed_scalar_payloads(
                |position| matches!(position, Position::FunctionArgument { .. }),
                |_, value| Ok(value.clone()),
            )
            .unwrap();
        assert_eq!(result.expression(), source.expression());
        assert_eq!(calls.load(Ordering::SeqCst), before);
        let result = source
            .map_closed_scalar_payloads(
                |position| matches!(position, Position::FunctionArgument { .. }),
                |_, value| value.with_checked_expression(value.expression() + Atom::one()),
            )
            .unwrap();
        assert_ne!(result.expression(), source.expression());
        assert!(calls.load(Ordering::SeqCst) > before);
    }

    #[test]
    fn selected_power_base_enters_only_original_nested_arguments() {
        test_initialize();
        let (f, generated, x) = symbol!(
            "closed_payload::nested_power_f",
            "closed_payload::nested_power_generated",
            "closed_payload::nested_power_x"
        );
        let base = expression(function!(f, x));
        let source = base.pow(&expression(Atom::num(3))).unwrap();
        let mut order = Vec::new();
        let result = source
            .map_closed_scalar_payloads(
                |position| {
                    matches!(
                        position,
                        Position::PowerBase { .. } | Position::FunctionArgument { .. }
                    )
                },
                |position, value| match position {
                    Position::FunctionArgument { function, .. } => {
                        assert_eq!(function.get_symbol(), f);
                        order.push("argument");
                        value.with_checked_expression(function!(generated, value.expression()))
                    }
                    Position::PowerBase { .. } => {
                        order.push("base");
                        Ok(value.clone())
                    }
                    _ => unreachable!(),
                },
            )
            .unwrap();
        assert_eq!(order, ["argument", "base"]);
        let expected = expression(function!(f, function!(generated, x)))
            .pow(&expression(Atom::num(3)))
            .unwrap();
        assert_eq!(result.expression(), expected.expression());

        let mut nested = false;
        source
            .map_closed_scalar_payloads(
                |position| {
                    nested |= matches!(position, Position::FunctionArgument { .. });
                    false
                },
                |_, _| unreachable!(),
            )
            .unwrap();
        assert!(!nested, "an unselected power base is opaque");
    }

    #[test]
    fn changed_ancestor_normalizers_cannot_change_the_outer_interface() {
        test_initialize();
        let head = tensor_symbol!("closed_payload::normalizer_open");
        let output = vector(head, 4601);
        let source_head = symbol!(
            "closed_payload::rank_changing_normalizer",
            norm = move |value, out| {
                if let AtomView::Fun(function) = value
                    && function.iter().next().is_some_and(|arg| arg.is_one())
                {
                    **out = output.clone();
                }
            }
        );
        let source = expression(function!(source_head, 0));
        assert!(
            source
                .map_closed_scalar_payloads(
                    |position| matches!(position, Position::FunctionArgument { .. }),
                    |_, _| Ok(expression(Atom::one())),
                )
                .is_err()
        );
    }
}
