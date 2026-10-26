//! Ordered symbolic products with explicit endpoints and their scalar contractions.

use std::collections::HashMap;

use spenso::{
    network::tags::SPENSO_TAG,
    shadowing::{FactorProjector, ProjectorExpander},
    structure::{
        abstract_index::AbstractIndex,
        partial::{PartialIndex, PartialStructure, PartialStructureExt},
        representation::{LibraryRep, LibrarySlot, RepName, Representation},
        slot::IsAbstractSlot,
    },
};
use symbolica::atom::Atom;

use super::{MatrixChannel, PortPair, SymbolicTensor, TensorCompositionError, port_atom};
use crate::tensor::inference::TensorInterface;

impl SymbolicTensor<PartialStructure> {
    /// Project an ordered group of compatible matrix factors, retaining an open
    /// channel so the group can be embedded in a larger chain or trace.
    /// The returned selected channel remains usable if projection produces zero.
    /// Groups whose coincident outer endpoints close into a scalar trace are
    /// outside this open-channel operation.
    pub fn project_factors(
        factors: &[Self],
        projector: FactorProjector,
    ) -> Result<(Self, MatrixChannel), TensorCompositionError> {
        let Some((first, rest)) = factors.split_first() else {
            return Err(TensorCompositionError::NoMatrixChannel);
        };
        for factor in factors {
            let channel = factor
                .matrix_channel()
                .ok_or(TensorCompositionError::NoMatrixChannel)?;
            if factor.chain_factors(channel)?.len() != 1 {
                return Err(TensorCompositionError::InvalidResultInterface(
                    "a factor projector takes individual matrix factors or projector groups; pass a chain's factors separately".to_owned(),
                ));
            }
        }
        let first_channel = first
            .matrix_channel()
            .ok_or(TensorCompositionError::NoMatrixChannel)?;
        let slots = first.structure.logical_slots();
        let input = slots[first_channel.input];
        let output = slots[first_channel.output];
        let mut value = first.chain_form(first_channel)?;
        let mut channel = MatrixChannel {
            input: 0,
            output: 1,
        }
        .remap_output(input, output, &value.structure)?;
        for factor in rest {
            let other = factor
                .matrix_channel()
                .ok_or(TensorCompositionError::NoMatrixChannel)?;
            let input = value.structure.logical_slots()[channel.input];
            let output = factor.structure.logical_slots()[other.output];
            value = value.compose(factor, channel, other)?;
            channel = MatrixChannel {
                input: 0,
                output: 1,
            }
            .remap_output(input, output, &value.structure)?;
        }
        let slots = value.structure.logical_slots();
        let input = slots[channel.input];
        let output = slots[channel.output];
        let (coefficient, payload) = value.chain_parts(channel)?;
        let projected = projector.apply(payload);
        let expression = if value.is_zero() || coefficient.is_zero() || projected.is_zero() {
            Atom::Zero
        } else {
            coefficient * SPENSO_TAG.chain(port_atom(input), port_atom(output), [projected])
        };
        // Composition may have materialized anonymous ports into explicit
        // contractions. Their history belongs to the composed carrier, even
        // when the projected payload vanishes.
        let occurrences = value.proofs.occurrences.clone();
        let mut value = Self::new(expression, value.structure);
        value.proofs.occurrences = occurrences;
        let value =
            value.validate_rewrite("project_factors", &factors.iter().collect::<Vec<_>>())?;
        let channel = channel.remap_output(input, output, &value.structure)?;
        Ok((value, channel))
    }

    /// Expand only normalized factor projectors, preserving the tensor interface
    /// and leaving unrelated products and sums factored.
    pub fn expand_projectors(&self) -> Result<Self, TensorCompositionError> {
        Self::new(
            self.expression.expand_projectors(),
            self.structure.clone().into_algebra_result(),
        )
        .validate_rewrite("expand_projectors", &[self])
    }

    /// Choose the explicit index shared by a contracted pair, avoiding dummy capture.
    pub fn contraction_index(
        &self,
        right: &Self,
        pair: PortPair,
    ) -> Result<AbstractIndex, TensorCompositionError> {
        let left = super::validate_position(&self.structure.logical_slots(), pair.left)?;
        let other = super::validate_position(&right.structure.logical_slots(), pair.right)?;
        Ok(super::shared_index(left, other, pair)?
            .unwrap_or_else(|| Self::reserved_dummies([self, right]).fresh_index()))
    }

    /// Fill selected unresolved ports, preserving the order of surviving ports.
    pub fn indexed(
        &self,
        replacements: &HashMap<usize, AbstractIndex>,
    ) -> Result<Self, TensorCompositionError> {
        let expression = self.materialize_interface_ports(replacements)?;
        let mut logical = self.structure.logical_slots();
        for (&position, &index) in replacements {
            logical[position].set_aind(PartialIndex::Explicit(index));
        }
        let value =
            Self::from_produced_parts(expression, PartialStructure::from_logical_slots(logical))
                .map_err(|error| TensorCompositionError::InvalidResultInterface(error.to_string()))?
                .normalize_closed_root_chain()?;
        value.validate_rewrite("indexed", &[self])
    }

    /// Explicitly choose dot notation for the sole ports of two rank-one tensors.
    /// Unlike automatic multiplication, this does not require rank-one head tags.
    pub fn dot(&self, right: &Self) -> Result<Self, TensorCompositionError> {
        if self.rank() != 1 || right.rank() != 1 {
            return Err(TensorCompositionError::InvalidResultInterface(format!(
                "dot requires rank-one operands, got ranks {} and {}",
                self.rank(),
                right.rank()
            )));
        }
        super::shared_index(
            self.structure.logical_slots()[0],
            right.structure.logical_slots()[0],
            PortPair { left: 0, right: 0 },
        )?;
        let atom = if self.expression.is_zero() || right.expression.is_zero() {
            Atom::Zero
        } else {
            self.compact_dot_expression(right)?
        };
        super::validate_explicit_index_occurrences(&atom)?;
        Self::new(atom, PartialStructure::from_logical_slots([]))
            .validate_rewrite("dot", &[self, right])
    }

    /// Give one channel explicit endpoints without changing leaf arguments.
    /// Named free axes use canonical order; positional axes retain their roles.
    pub fn chain_form(&self, channel: MatrixChannel) -> Result<Self, TensorCompositionError> {
        if channel.input == channel.output {
            return Err(TensorCompositionError::DegenerateChannel {
                input: channel.input,
                output: channel.output,
            });
        }
        let slots = self.structure.logical_slots();
        let input = super::validate_position(&slots, channel.input)?;
        let output = super::validate_position(&slots, channel.output)?;
        if !input.rep().matches(&output.rep())
            || !(input.rep().rep.is_self_dual()
                || input.rep().rep.is_base() && output.rep().rep.is_dual())
        {
            return Err(TensorCompositionError::InvalidChannelOrientation {
                input: channel.input,
                output: channel.output,
            });
        }
        let expression = if self.expression.as_view().is_zero() {
            Atom::Zero
        } else {
            let (coefficient, factors) = self.chain_parts(channel)?;
            coefficient * SPENSO_TAG.chain(port_atom(input), port_atom(output), factors)
        };
        let structure = Self::composite_interface(
            [input, output].into_iter().chain(
                slots
                    .into_iter()
                    .enumerate()
                    .filter(|(position, _)| {
                        *position != channel.input && *position != channel.output
                    })
                    .map(|(_, slot)| slot),
            ),
        );
        Self::new(expression, structure).validate_rewrite("chain_form", &[self])
    }

    /// Build an ordered chain with explicit input and output slots.
    pub fn chain(
        start: LibrarySlot<AbstractIndex>,
        end: LibrarySlot<AbstractIndex>,
        factors: &[Self],
    ) -> Result<Self, TensorCompositionError> {
        let Some((first, factors)) = factors.split_first() else {
            if !start.rep().matches(&end.rep()) {
                return Err(TensorCompositionError::IncompatiblePorts { left: 0, right: 1 });
            }
            if !(start.rep().rep.is_self_dual()
                || start.rep().rep.is_base() && end.rep().rep.is_dual())
            {
                return Err(TensorCompositionError::InvalidChannelOrientation {
                    input: 0,
                    output: 1,
                });
            }
            let (expression, structure) = if start.aind() == end.aind() {
                (
                    SPENSO_TAG.trace(start.rep().to_symbolic([]), std::iter::empty::<Atom>()),
                    PartialStructure::from_logical_slots([]),
                )
            } else {
                (
                    SPENSO_TAG.chain(start.to_atom(), end.to_atom(), std::iter::empty::<Atom>()),
                    Self::composite_interface([
                        start.rep().slot(PartialIndex::Explicit(start.aind())),
                        end.rep().slot(PartialIndex::Explicit(end.aind())),
                    ]),
                )
            };
            return Self::from_produced_parts(expression, structure).map_err(|error| {
                TensorCompositionError::InvalidResultInterface(error.to_string())
            });
        };

        let channel = first
            .matrix_channel()
            .ok_or(TensorCompositionError::NoMatrixChannel)?;
        let slots = first.structure.logical_slots();
        let input = slots[channel.input];
        let output = slots[channel.output];
        let mut value = first.chain_form(channel)?;
        let mut channel = MatrixChannel {
            input: 0,
            output: 1,
        }
        .remap_output(input, output, &value.structure)?;
        for factor in factors {
            let factor_channel = factor
                .matrix_channel()
                .ok_or(TensorCompositionError::NoMatrixChannel)?;
            let input = value.structure.logical_slots()[channel.input];
            let output = factor.structure.logical_slots()[factor_channel.output];
            value = value.compose(factor, channel, factor_channel)?;
            channel = MatrixChannel {
                input: 0,
                output: 1,
            }
            .remap_output(input, output, &value.structure)?;
        }
        let slots = value.structure.logical_slots();
        if start.rep() != slots[channel.input].rep() || end.rep() != slots[channel.output].rep() {
            return Err(TensorCompositionError::InvalidResultInterface(
                "chain endpoints are incompatible with the factor channel".into(),
            ));
        }
        if start.aind() == end.aind() {
            value.trace_ports(channel)
        } else {
            let mut value = value.indexed(&HashMap::from([
                (channel.input, start.aind()),
                (channel.output, end.aind()),
            ]))?;
            value.structure = value.structure.into_algebra_result();
            Ok(value)
        }
    }

    /// Close an ordered factor sequence into a canonical cyclic trace.
    pub fn trace(
        representation: Representation<LibraryRep>,
        factors: &[Self],
    ) -> Result<Self, TensorCompositionError> {
        let Some((first, factors)) = factors.split_first() else {
            return Self::from_produced_parts(
                SPENSO_TAG.trace(representation.to_symbolic([]), std::iter::empty::<Atom>()),
                PartialStructure::from_logical_slots([]),
            )
            .map_err(|error| TensorCompositionError::InvalidResultInterface(error.to_string()));
        };
        let mut value = first.clone();
        for factor in factors {
            value = value.multiply(factor)?;
        }
        let channel = value
            .matrix_channel()
            .ok_or(TensorCompositionError::NoMatrixChannel)?;
        if representation != value.structure.logical_slots()[channel.input].rep() {
            return Err(TensorCompositionError::InvalidResultInterface(
                "trace representation does not match the factor channel".into(),
            ));
        }
        value.trace_ports(channel)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support;
    use spenso::{
        network::library::symbolic::ETS,
        structure::{dimension::Dimension, representation::ExtendibleReps},
    };
    use symbolica::atom::{AtomView, FunctionBuilder};

    fn matrix(name: &str, rep: Representation<LibraryRep>) -> SymbolicTensor<PartialStructure> {
        SymbolicTensor::new(
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_args([rep.to_symbolic([]), rep.to_symbolic([])])
                .finish(),
            PartialStructure::from_logical_slots([
                rep.slot(PartialIndex::open(0)),
                rep.slot(PartialIndex::open(1)),
            ]),
        )
    }

    #[test]
    fn zero_factor_projection_retains_materialized_contraction_history() {
        test_support::test_initialize();
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let source = matrix("projected_history_zero", rep);
        let zero = source.try_sub(&source).unwrap();
        let right = matrix("projected_history_right", rep)
            .reindex_interface_ports(&HashMap::from([(0, AbstractIndex::Normal(50))]))
            .unwrap();
        let channel = MatrixChannel {
            input: 0,
            output: 1,
        };
        let direct = zero.compose(&right, channel, channel).unwrap();
        let (projected, selected) =
            SymbolicTensor::project_factors(&[zero, right], FactorProjector::Symmetric).unwrap();
        assert!(projected.is_zero());
        assert_eq!(projected.structure, direct.structure);
        assert_eq!(selected, channel);
        let index = rep.slot(AbstractIndex::Normal(50));
        assert_eq!(direct.explicit_index_occurrences().unwrap().count(index), 2);
        assert_eq!(
            projected.explicit_index_occurrences().unwrap().count(index),
            2
        );
        let probe = SymbolicTensor::infer(
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol("projected_history_probe"))
                .add_arg(index.to_atom())
                .finish(),
        )
        .unwrap();
        for value in [direct, projected] {
            assert!(matches!(
                value.outer(&probe),
                Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
            ));
        }
    }

    #[test]
    fn factor_projection_returns_selected_zero_channel_and_rejects_closed_groups() {
        test_support::test_initialize();
        let channel = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let spectator = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let factor = |name, index| {
            SymbolicTensor::infer(
                FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                    .add_arg(channel.to_symbolic([]))
                    .add_arg(channel.to_symbolic([]))
                    .add_arg(
                        spectator
                            .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                            .to_atom(),
                    )
                    .finish(),
            )
            .unwrap()
        };
        let left = factor("selected_projection_left", 7);
        let zero = left.try_sub(&left).unwrap();
        let right = factor("selected_projection_right", 11);
        let (projected, selected) =
            SymbolicTensor::project_factors(&[zero, right], FactorProjector::Symmetric).unwrap();
        assert!(projected.is_zero());
        assert_eq!(projected.rank(), 4);
        assert_eq!(projected.matrix_channel(), None);
        assert_eq!(
            selected,
            MatrixChannel {
                input: 0,
                output: 1
            }
        );
        let slots = projected.structure.logical_slots();
        assert_eq!(slots[selected.input].rep(), channel);
        assert_eq!(slots[selected.output].rep(), channel);

        let named = |name, input, output| {
            let ports = [input, output]
                .map(|index| channel.slot(PartialIndex::Explicit(AbstractIndex::Normal(index))));
            SymbolicTensor::from_declared_parts(
                FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                    .add_args(ports.map(port_atom))
                    .finish(),
                PartialStructure::from_logical_slots(ports),
            )
            .unwrap()
        };
        let first = named("selected_projection_closed_first", 90, 50);
        let second = named("selected_projection_closed_second", 50, 90);
        assert!(
            SymbolicTensor::project_factors(&[first, second], FactorProjector::Symmetric).is_err()
        );
    }

    #[test]
    fn chains_keep_factor_order_and_close_through_the_same_trace_owner() {
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let start = rep.slot(AbstractIndex::Normal(73001));
        let end = rep.slot(AbstractIndex::Normal(73003));
        let factors = [
            matrix("shared_chain_first", rep),
            matrix("shared_chain_second", rep),
        ];
        assert!(
            factors[0]
                .chain_form(MatrixChannel {
                    input: 0,
                    output: 0
                })
                .is_err()
        );
        assert!(
            factors[0]
                .chain_form(MatrixChannel {
                    input: 0,
                    output: 2
                })
                .is_err()
        );
        let chain = SymbolicTensor::chain(start, end, &factors).unwrap();
        assert_eq!(
            chain.structure.logical_slots(),
            vec![
                rep.slot(PartialIndex::Explicit(start.aind())),
                rep.slot(PartialIndex::Explicit(end.aind())),
            ]
        );
        let AtomView::Fun(fun) = chain.expression.as_view() else {
            panic!("expected a chain")
        };
        assert_eq!(fun.get_symbol(), SPENSO_TAG.chain);
        let factor_symbols = fun
            .iter()
            .skip(2)
            .map(|factor| {
                let AtomView::Fun(fun) = factor else {
                    panic!("expected a tensor factor")
                };
                fun.get_symbol()
            })
            .collect::<Vec<_>>();
        assert_eq!(
            factor_symbols,
            factors
                .iter()
                .map(|factor| {
                    let AtomView::Fun(fun) = factor.expression.as_view() else {
                        unreachable!()
                    };
                    fun.get_symbol()
                })
                .collect::<Vec<_>>()
        );
        let closed = SymbolicTensor::chain(start, start, &factors).unwrap();
        let traced = SymbolicTensor::trace(rep, &factors).unwrap();
        assert!(closed.is_scalar());
        assert_eq!(closed.expression, traced.expression);
        assert_eq!(
            SymbolicTensor::trace(rep, &[]).unwrap().expression,
            SPENSO_TAG.trace(rep.to_symbolic([]), std::iter::empty::<Atom>())
        );
    }

    #[test]
    fn indexing_contracts_repeated_ports_and_preserves_typed_zero() {
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let index = AbstractIndex::Normal(73011);
        let value = matrix("shared_indexed_matrix", rep);
        let indexed = value
            .indexed(&HashMap::from([(0, index), (1, index)]))
            .unwrap();
        assert!(indexed.is_scalar());
        let zero = SymbolicTensor::new(Atom::Zero, value.structure.clone());
        let indexed_zero = zero.indexed(&HashMap::from([(0, index)])).unwrap();
        assert_eq!(indexed_zero.rank(), 2);
        assert!(indexed_zero.expression.is_zero());
        let incompatible = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        assert!(SymbolicTensor::chain(rep.slot(index), incompatible.slot(index), &[]).is_err());
        assert!(value.dot(&value).is_err());
        let compact = FunctionBuilder::new(ETS.metric)
            .add_args([rep.to_symbolic([]), rep.to_symbolic([])])
            .finish();
        assert_eq!(SymbolicTensor::infer(compact).unwrap().rank(), 2);
    }

    #[test]
    fn graph_reindexing_retains_explicit_storage_open_identities() {
        let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let value = matrix("shared_storage_open_matrix", rep);
        let owner = AbstractIndex::fresh_open_owner();
        let indices = [0, 1].map(|axis| AbstractIndex::Open { owner, axis });
        let indexed = value
            .reindex_interface_ports(&HashMap::from([(0, indices[0]), (1, indices[1])]))
            .unwrap();
        assert_eq!(
            indexed
                .structure
                .logical_slots()
                .iter()
                .map(|slot| slot.aind)
                .collect::<Vec<_>>(),
            indices.map(PartialIndex::Explicit),
        );
    }
}
