//! Tensor-aware interface inference shared by Rust and Python entry points.
//!
//! Cached interfaces retain raw logical ports. Only enclosing algebra contracts
//! explicit pairs; unresolved port identities are canonicalized per occurrence.

use std::collections::{HashMap, HashSet};

use spenso::{
    network::{
        library::symbolic::{ETS, ExplicitKey},
        parsing::{AtomStructureExt, StrictTensorFilter, StructureInferenceMode},
        tags::SPENSO_TAG,
    },
    shadowing,
    structure::{
        Canonicalized, OrderedStructure, StructureError, TensorStructure,
        abstract_index::AbstractIndex,
        dimension::Dimension,
        partial::{PartialIndex, PartialSlot, PartialStructure, PartialStructureExt},
        representation::{LibraryRep, Minkowski, RepName, Representation},
        slot::{IsAbstractSlot, Slot, SlotMatch, SlotMatcher},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    domains::rational::Rational,
};
use thiserror::Error;

use super::{
    SymbolicTensor,
    composition::{self, TensorCompositionError},
};
use crate::{
    color::CS,
    dirac::{AGS, DiracSimplifier},
    epsilon::EPSILON_SYMBOL,
    representations::{Bispinor, ColorAdjoint, ColorAntiFundamental, ColorFundamental},
};

#[derive(Debug, Error)]
pub enum TensorInferenceError {
    #[error("{0}")]
    Invalid(String),
    #[error("predefined tensor `{factory}` requires {signature}")]
    InvalidBuiltinSignature {
        factory: &'static str,
        signature: &'static str,
    },
    #[error(transparent)]
    Composition(#[from] TensorCompositionError),
}

impl TensorInferenceError {
    fn invalid(message: impl Into<String>) -> Self {
        Self::Invalid(message.into())
    }
}

type InferenceResult<T> = Result<T, TensorInferenceError>;

impl SymbolicTensor<PartialStructure> {
    /// Construct an unresolved tensor from a named, canonically ordered signature.
    pub fn from_signature(
        value: &Canonicalized<ExplicitKey<AbstractIndex>>,
    ) -> InferenceResult<Self> {
        let signature = value.canonical();
        let name = signature
            .global_name
            .ok_or_else(|| TensorInferenceError::invalid("tensor structure has no name"))?;
        let arguments = signature.additional_args.as_deref().unwrap_or_default();
        let canonical = signature.external_reps_iter().collect::<Vec<_>>();
        let logical = value.layout().canonical_to_logical(&canonical);
        let axes = (0..logical.len()).collect::<Vec<_>>();
        let canonical_axes = value.layout().logical_to_canonical(&axes);
        let owner = AbstractIndex::fresh_open_owner();
        let expression = FunctionBuilder::new(name)
            .add_args(arguments)
            .add_args(
                canonical
                    .iter()
                    .zip(canonical_axes)
                    .map(|(representation, axis)| {
                        representation
                            .slot::<AbstractIndex, _>(AbstractIndex::Open { owner, axis })
                            .to_atom()
                    }),
            )
            .finish();
        if !InterfaceInference::intrinsic_normalization_head(name) {
            // The constructor spells ports in canonical storage order, while
            // the declared layout retains their public logical order. Check a
            // user normalizer against that physical signature before attaching
            // the unchanged logical layout.
            let physical = PartialStructure::from_logical_slots(canonical.iter().enumerate().map(
                |(position, representation)| representation.slot(PartialIndex::open(position)),
            ));
            Self::validate_interface(&expression, &physical)?;
        }
        let structure =
            PartialStructure::from_logical_slots(logical.into_iter().enumerate().map(
                |(position, representation)| representation.slot(PartialIndex::open(position)),
            ));
        Self::checked_parts(expression, structure)
    }

    /// Infer an ordered interface from symbolic syntax, lowering unresolved powers.
    pub fn infer(atom: Atom) -> InferenceResult<Self> {
        let atom = InterfaceInference::lower_tensor_powers(atom.as_view())?.unwrap_or(atom);
        let structure = InterfaceInference::default().infer_validated(atom.as_view())?;
        Ok(Self::new(atom, structure))
    }

    /// Finish tensor-aware parts without inferring their established interface again.
    pub fn checked_parts(expression: Atom, structure: PartialStructure) -> InferenceResult<Self> {
        Self::validate_atom(&expression)?;
        let mut value = Self::new(expression, structure).normalize_closed_root_chain()?;
        value.structure =
            InterfaceInference::merge_explicit_interface_sequence(&[value.structure])?
                .canonicalize_open_ports();
        Ok(value)
    }

    /// Retain this interface after scalar algebra which treats tensor leaves as opaque.
    pub fn with_algebra_result(&self, expression: Atom) -> InferenceResult<Self> {
        if expression == self.expression {
            return Ok(Self {
                expression,
                structure: self.structure.clone(),
                is_metric: self.is_metric,
                is_composite: self.is_composite,
            });
        }
        if expression.as_view().is_zero()
            || (self.structure.open_positions().is_empty()
                && !self.expression.as_view().needs_normalization()
                && InterfaceInference::default()
                    .algebra_preserves_leaf_interfaces(self.expression.as_view()))
        {
            return Ok(Self::from_normalized_parts(
                expression,
                self.structure.clone(),
            ));
        }
        Self::validate_interface(&expression, &self.structure)?;
        Self::checked_parts(expression, self.structure.clone())
    }

    /// Retain this interface after a tensor identity, checking callback-sensitive results.
    /// Arbitrary replacements must validate or infer their new interface instead.
    pub fn with_rewritten_expression(&self, expression: Atom) -> InferenceResult<Self> {
        if expression == self.expression {
            return Ok(Self {
                expression,
                structure: self.structure.clone(),
                is_metric: self.is_metric,
                is_composite: self.is_composite,
            });
        }
        if expression.as_view().is_zero()
            || (self.structure.open_positions().is_empty()
                && !self.expression.as_view().needs_normalization()
                && {
                    let mut inference = InterfaceInference::default();
                    inference.rewrites_preserve_leaf_interfaces(self.expression.as_view())
                        // Reserve source validation for substantial trace expansion.
                        || (expression.as_view().get_byte_size()
                            > self.expression.as_view().get_byte_size().saturating_mul(2)
                            && inference.terminal_trace_preserves_interface(
                                &self.expression,
                                &self.structure,
                            ))
                })
        {
            return Ok(Self::from_normalized_parts(
                expression,
                self.structure.clone(),
            ));
        }
        Self::validate_interface(&expression, &self.structure)?;
        Self::checked_parts(expression, self.structure.clone())
    }

    /// Infer the result of an arbitrary symbolic transformation.
    pub fn with_transformed_expression(&self, expression: Atom) -> InferenceResult<Self> {
        if expression == self.expression {
            return Ok(Self {
                expression,
                structure: self.structure.clone(),
                is_metric: self.is_metric,
                is_composite: self.is_composite,
            });
        }
        let structure = if expression.as_view().is_zero() {
            self.structure.clone()
        } else if InterfaceInference::has_structured_syntax(expression.as_view()) {
            Self::infer(expression.clone())?.structure
        } else if self.is_scalar() {
            PartialStructure::from_logical_slots([])
        } else {
            return Err(TensorInferenceError::invalid(
                "Tensor transformation removed the tensor syntax of a non-scalar expression",
            ));
        };
        Self::checked_parts(expression, structure)
    }

    /// Change four-dimensional Lorentz ports in the expression and its logical interface.
    /// Newly equal explicit indices contract; callback-sensitive results must still
    /// carry the resulting interface without replaying normalization during validation.
    pub fn with_lorentz_dimension(&self, dimension: Dimension) -> InferenceResult<Self> {
        let expression = self
            .expression
            .with_lorentz_dimension(dimension.to_symbolic().as_view());
        let mapped_structure = |structure: &PartialStructure| {
            PartialStructure::from_logical_slots(structure.logical_slots().into_iter().map(
                |slot| {
                    let mut representation = slot.rep();
                    if representation.rep == (Minkowski {}).into()
                        && representation.dim == Dimension::Concrete(4)
                    {
                        representation.dim = dimension;
                    }
                    representation.slot(slot.aind)
                },
            ))
        };
        let mut value = Self::checked_parts(expression, mapped_structure(&self.structure))?;
        if !InterfaceInference::normalization_is_intrinsic(self.expression.as_view()) {
            // Named signatures may encode ports in canonical storage order while
            // exposing a different logical order. Compare the two encoded forms;
            // retain the independently mapped public layout on the result.
            let physical = Self::observed_interface(&self.expression)?;
            let expected =
                InterfaceInference::merge_explicit_interface_sequence(&[mapped_structure(
                    &physical,
                )])?;
            Self::validate_interface(&value.expression, &expected)?;
        }
        if value.expression == self.expression {
            value.is_metric = self.is_metric;
            value.is_composite = self.is_composite;
        }
        Ok(value)
    }

    /// Cook both encoded indices and their logical interface, contracting collisions.
    pub fn with_cooked_indices(&self, settings: &crate::CookSettings) -> InferenceResult<Self> {
        let cook = |atom: &Atom| {
            settings.try_cook_indices(atom.as_view()).map_err(|error| {
                TensorInferenceError::invalid(format!("cannot cook indices: {error:?}"))
            })
        };
        let expression = cook(&self.expression)?;
        if !self.structure.open_positions().is_empty() && !expression.as_view().is_zero() {
            // Unresolved interface metadata does not encode the occurrence which
            // cooking may turn into an explicit index.
            return self.with_transformed_expression(expression);
        }
        let slots = self
            .structure
            .logical_slots()
            .into_iter()
            .map(|slot| {
                if matches!(slot.aind, PartialIndex::Open(_)) {
                    return Ok(slot);
                }
                let cooked = cook(&composition::port_atom(slot))?;
                let cooked = Slot::<LibraryRep, AbstractIndex>::try_from(cooked.as_view())
                    .map_err(|error| {
                        TensorInferenceError::invalid(format!(
                            "cannot cook interface slot: {error}"
                        ))
                    })?;
                Ok(cooked.rep().slot(PartialIndex::Explicit(cooked.aind())))
            })
            .collect::<InferenceResult<Vec<_>>>()?;
        let value = Self::checked_parts(expression, PartialStructure::from_logical_slots(slots))?;
        if !InterfaceInference::default()
            .rewrites_preserve_leaf_interfaces(self.expression.as_view())
        {
            Self::validate_interface(&value.expression, &value.structure)?;
        }
        Ok(value)
    }

    /// Add tensors with compatible logical interfaces, retaining the left order.
    pub fn try_add(&self, right: &Self) -> InferenceResult<Self> {
        if !InterfaceInference::additive_interfaces_match(&self.structure, &right.structure) {
            return Err(TensorInferenceError::invalid(
                "addition requires compatible tensor interfaces",
            ));
        }
        Self::checked_parts(&self.expression + &right.expression, self.structure.clone())
    }

    /// Subtract tensors with compatible logical interfaces, retaining typed zero.
    pub fn try_sub(&self, right: &Self) -> InferenceResult<Self> {
        if !InterfaceInference::additive_interfaces_match(&self.structure, &right.structure) {
            return Err(TensorInferenceError::invalid(
                "subtraction requires compatible tensor interfaces",
            ));
        }
        Self::checked_parts(&self.expression - &right.expression, self.structure.clone())
    }

    /// Validate index multiplicity and the scopes of chain and trace placeholders.
    pub fn validate_atom(atom: &Atom) -> InferenceResult<()> {
        atom.validate_chain_like_nesting()
            .map_err(|error| TensorInferenceError::invalid(error.to_string()))?;
        InterfaceInference::validate_placeholder_scope(atom.as_view(), false)
            .map_err(TensorInferenceError::invalid)?;
        composition::validate_explicit_index_occurrences(atom)?;
        Ok(())
    }

    /// Check a rewrite against an established interface, retaining typed zeros.
    ///
    /// Index substitution can invoke a user normalizer which changes tensor rank;
    /// inspect the resulting syntax without invoking that normalizer again with
    /// synthetic indices. Constructors still materialize compact callback leaves.
    pub fn validate_interface(atom: &Atom, structure: &PartialStructure) -> InferenceResult<()> {
        if atom.as_view().is_zero() {
            return Ok(());
        }
        let inferred = Self::observed_interface(atom)?;
        if !InterfaceInference::additive_interfaces_match(structure, &inferred) {
            return Err(TensorInferenceError::invalid(
                "transformed expression does not preserve a compatible tensor interface",
            ));
        }
        Ok(())
    }

    /// Observe encoded port order without materializing compact callback leaves.
    fn observed_interface(atom: &Atom) -> InferenceResult<PartialStructure> {
        if !InterfaceInference::has_structured_syntax(atom.as_view()) {
            return Ok(PartialStructure::from_logical_slots([]));
        }
        Ok(
            InterfaceInference::merge_explicit_interface_sequence(&[InterfaceInference {
                leaf_inference: LeafInference::Observe,
                ..InterfaceInference::default()
            }
            .infer_validated(atom.as_view())?])?
            .canonicalize_open_ports(),
        )
    }
}

impl InterfaceInference {
    pub(crate) fn intrinsic_normalization_head(symbol: Symbol) -> bool {
        symbol.get_evaluation_info().is_none()
            && (symbol == ETS.metric || symbol.get_normalization_function().is_none())
    }

    /// Index substitution can retain a planned interface when reconstruction
    /// invokes only intrinsic normalization, including within scalar metadata.
    pub fn normalization_is_intrinsic(value: AtomView<'_>) -> bool {
        let mut intrinsic = true;
        value.visitor(&mut |node| {
            if let AtomView::Fun(function) = node {
                intrinsic &= Self::intrinsic_normalization_head(function.get_symbol());
            }
            intrinsic
        });
        intrinsic
    }

    /// Certify that index rewriting cannot invoke a user callback through tensor
    /// leaves or scalar metadata. The intrinsic metric normalizer is the sole
    /// allowed hook; all its vector operands are checked by the leaf proof.
    pub fn rewrites_preserve_leaf_interfaces(&mut self, value: AtomView<'_>) -> bool {
        Self::normalization_is_intrinsic(value) && self.algebra_preserves_leaf_interfaces(value)
    }

    /// Prove terminal trace identities on the short source word instead of
    /// inferring the interface of the resulting pairing polynomial.
    fn terminal_trace_preserves_interface(
        &mut self,
        expression: &Atom,
        structure: &PartialStructure,
    ) -> bool {
        let Some((cyclic, indices)) =
            DiracSimplifier::terminal_trace_interface_inputs(expression.as_view())
        else {
            return false;
        };
        // Canonical traces use this intrinsic cyclic wrapper. Its normalizer
        // only unwraps a symmetric projector, which terminal-word admission
        // excludes; it cannot invoke a component callback.
        let mut intrinsic = true;
        expression.as_view().visitor(&mut |node| {
            if let AtomView::Fun(function) = node {
                let symbol = function.get_symbol();
                intrinsic &= if Some(node) == cyclic {
                    symbol.get_evaluation_info().is_none()
                } else {
                    Self::intrinsic_normalization_head(symbol)
                };
            }
            intrinsic
        });
        if !intrinsic {
            return false;
        }
        if !indices
            .into_iter()
            .all(|index| match self.slots.classify(index) {
                SlotMatch::Explicit(_) => self
                    .slots
                    .parse::<LibraryRep, AbstractIndex>(index)
                    .is_ok_and(|slot| {
                        !matches!(slot.aind(), AbstractIndex::Open { .. })
                            && slot.rep().rep.is_self_dual()
                    }),
                SlotMatch::Other => {
                    let AtomView::Fun(vector) = index else {
                        return false;
                    };
                    self.direct_leaf_interface(vector).is_some_and(|interface| {
                        matches!(interface.logical_slots().as_slice(), [slot]
                        if matches!(slot.aind, PartialIndex::Open(_))
                            && slot.rep().rep.is_self_dual())
                    })
                }
                SlotMatch::Opaque => false,
            })
        {
            return false;
        }
        SymbolicTensor::<PartialStructure>::validate_atom(expression).is_ok()
            && SymbolicTensor::<PartialStructure>::validate_interface(expression, structure).is_ok()
    }

    /// Scalar algebra treats tensor leaves as indeterminates. Reuse their
    /// interfaces only when inference cannot run callbacks or change power semantics.
    pub fn algebra_preserves_leaf_interfaces(&mut self, value: AtomView<'_>) -> bool {
        match value {
            AtomView::Add(sum) => sum
                .iter()
                .all(|term| self.algebra_preserves_leaf_interfaces(term)),
            AtomView::Mul(product) => product
                .iter()
                .all(|factor| self.algebra_preserves_leaf_interfaces(factor)),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                // Existing occurrence validation does not multiply power counts.
                // Only opaque dots of compact vectors expose no repeated indices.
                let scalar_dot = matches!(base, AtomView::Fun(function)
                if self.direct_dot_ports(function).is_some_and(|ports| {
                    ports.iter().all(|slot| matches!(slot.aind, PartialIndex::Open(_)))
                }));
                !Self::has_structured_syntax(exponent)
                    && (!Self::has_structured_syntax(base) || scalar_dot)
                    && self.algebra_preserves_leaf_interfaces(base)
                    && self.algebra_preserves_leaf_interfaces(exponent)
            }
            AtomView::Fun(function) if Self::has_structured_syntax(value) => {
                let symbol = function.get_symbol();
                let arguments = function.iter().collect::<Vec<_>>();
                if self.direct_dot_ports(function).is_some() {
                    // Compact dots consume their operands' ports inside the
                    // opaque function. Full inference still checks compatibility.
                    return true;
                }
                if arguments.iter().all(|argument| {
                    self.slots
                        .parse::<LibraryRep, AbstractIndex>(*argument)
                        .is_ok_and(|slot| {
                            !matches!(slot.aind(), AbstractIndex::Open { .. })
                                && slot.rep().rep.is_self_dual()
                        })
                }) && self
                    .builtin_tensor_structure(symbol, &arguments)
                    .is_ok_and(|structure| structure.is_some())
                {
                    return true;
                }
                self.direct_leaf_interface(function)
                    .is_some_and(|interface| {
                        interface.logical_slots().iter().all(|slot| {
                            matches!(slot.aind, PartialIndex::Explicit(_))
                                && slot.rep().rep.is_self_dual()
                        })
                    })
            }
            // Scalar functions own their metadata, so expansion leaves it opaque.
            _ => true,
        }
    }

    fn direct_dot_ports(
        &mut self,
        function: symbolica::atom::representation::FunView<'_>,
    ) -> Option<[PartialSlot; 2]> {
        if !(function.get_symbol() == ETS.metric || function.get_symbol() == SPENSO_TAG.dot)
            || function.get_nargs() != 2
        {
            return None;
        }
        let mut ports = Vec::with_capacity(2);
        for operand in function.iter() {
            let AtomView::Fun(operand) = operand else {
                return None;
            };
            let interface = self.direct_leaf_interface(operand)?;
            let slots = interface.logical_slots();
            let [slot] = slots.as_slice() else {
                return None;
            };
            ports.push(*slot);
        }
        ports.try_into().ok()
    }

    fn builtin_tensor_structure(
        &mut self,
        symbol: Symbol,
        arguments: &[AtomView<'_>],
    ) -> InferenceResult<Option<Canonicalized<ExplicitKey<AbstractIndex>>>> {
        let builtin = symbol == *EPSILON_SYMBOL
            || symbol == ETS.metric
            || symbol == ETS.flat
            || symbol == AGS.gamma
            || symbol == AGS.gamma5
            || symbol == AGS.projm
            || symbol == AGS.projp
            || symbol == AGS.sigma
            || symbol == CS.f
            || symbol == CS.t;
        if !builtin {
            return Ok(None);
        }

        let compact_metric = symbol == ETS.metric
            && arguments.len() == 2
            && arguments.iter().all(|argument| {
                self.slots
                    .parse::<LibraryRep, AbstractIndex>(*argument)
                    .is_err()
                    && self
                        .slots
                        .parse_representation::<LibraryRep>(*argument)
                        .is_err()
                    && !Self::is_chain_placeholder(*argument)
                    && Self::has_structured_syntax(*argument)
            });
        if compact_metric {
            return Ok(None);
        }

        let ports = arguments
            .iter()
            .map(|argument| {
                if Self::is_chain_placeholder(*argument) {
                    Ok(None)
                } else {
                    self.slots
                        .parse::<LibraryRep, AbstractIndex>(*argument)
                        .map(|slot| Some(slot.rep()))
                        .or_else(|_| {
                            self.slots
                                .parse_representation::<LibraryRep>(*argument)
                                .map(Some)
                        })
                        .map_err(|_| ())
                        .or_else(|_| {
                            // A compact vector consumes this port, while its representation
                            // still participates in the predefined tensor's signature.
                            let interface = self.infer_validated(*argument).map_err(|_| ())?;
                            let slots = interface.logical_slots();
                            match slots.as_slice() {
                                [slot] if matches!(slot.aind, PartialIndex::Open(_)) => {
                                    Ok(Some(slot.rep()))
                                }
                                _ => Err(()),
                            }
                        })
                }
            })
            .collect::<Result<Vec<_>, _>>();
        let dimension = |representation: LibraryRep| {
            ports
                .as_ref()
                .ok()
                .and_then(|ports| {
                    ports
                        .iter()
                        .flatten()
                        .find(|port| port.rep == representation)
                })
                .map(|port| port.dim)
        };
        let invalid = |factory: &'static str, signature: &'static str| {
            TensorInferenceError::InvalidBuiltinSignature { factory, signature }
        };

        let (factory, signature, expected): (
            _,
            _,
            Option<Canonicalized<ExplicitKey<AbstractIndex>>>,
        ) = if symbol == ETS.metric || symbol == ETS.flat {
            let factory = if symbol == ETS.metric { "g" } else { "flat" };
            let ports = ports
                .as_ref()
                .map_err(|_| invalid(factory, "two compatible representation ports"))?;
            if arguments.len() != 2 {
                return Err(invalid(factory, "two compatible representation ports"));
            }
            let visible = ports.iter().flatten().copied().collect::<Vec<_>>();
            let compatible = visible.len() < 2
                || if symbol == ETS.metric {
                    visible[0] == visible[1] || visible[0].matches(&visible[1])
                } else {
                    visible[0] == visible[1]
                };
            if !compatible {
                return Err(invalid(factory, "two compatible representation ports"));
            }
            let expected = match visible.as_slice() {
                [left, right] => Some(ExplicitKey::from_iter([*left, *right], symbol, None)),
                _ => None,
            };
            (factory, "two compatible representation ports", expected)
        } else if symbol == *EPSILON_SYMBOL {
            let signature = "equal representation ports or contracted vectors";
            let ports = ports.as_ref().map_err(|_| invalid("epsilon", signature))?;
            let visible = ports
                .iter()
                .copied()
                .collect::<Option<Vec<_>>>()
                .ok_or_else(|| invalid("epsilon", signature))?;
            if visible.is_empty() || visible.iter().any(|port| *port != visible[0]) {
                return Err(invalid("epsilon", signature));
            }
            (
                "epsilon",
                signature,
                Some(ExplicitKey::from_iter(visible, symbol, None)),
            )
        } else if symbol == AGS.gamma {
            (
                "gamma",
                "one Minkowski and two four-dimensional bispinor ports",
                Some(
                    AGS.gamma_strct(
                        dimension(Minkowski {}.into()).unwrap_or(Dimension::Concrete(1)),
                    ),
                ),
            )
        } else if symbol == AGS.gamma5 || symbol == AGS.projm || symbol == AGS.projp {
            let dimension = dimension(Bispinor {}.into()).unwrap_or(Dimension::Concrete(1));
            let (factory, expected) = if symbol == AGS.gamma5 {
                ("gamma5", AGS.gamma5_strct(dimension))
            } else if symbol == AGS.projm {
                ("projm", AGS.projm_strct(dimension))
            } else {
                ("projp", AGS.projp_strct(dimension))
            };
            (factory, "two equal bispinor ports", Some(expected))
        } else if symbol == AGS.sigma {
            (
                "sigma",
                "two equal Minkowski and two four-dimensional bispinor ports",
                Some(
                    AGS.sigma_strct(
                        dimension(Minkowski {}.into()).unwrap_or(Dimension::Concrete(1)),
                    ),
                ),
            )
        } else if symbol == CS.f {
            (
                "f",
                "three equal adjoint ports",
                Some(
                    CS.f_strct(dimension(ColorAdjoint {}.into()).unwrap_or(Dimension::Concrete(1))),
                ),
            )
        } else if symbol == CS.t {
            let adjoint_dimension =
                dimension(ColorAdjoint {}.into()).unwrap_or(Dimension::Concrete(1));
            let fundamental_dimension = dimension(ColorFundamental {}.into())
                .or_else(|| dimension(ColorAntiFundamental {}.into()))
                .unwrap_or(Dimension::Concrete(1));
            (
                "t",
                "adjoint, fundamental, and antifundamental ports with matching color dimensions",
                Some(CS.t_strct(fundamental_dimension, adjoint_dimension)),
            )
        } else {
            unreachable!("every predefined tensor symbol was handled")
        };

        let ports = ports.map_err(|_| invalid(factory, signature))?;
        let placeholders = arguments
            .iter()
            .enumerate()
            .filter_map(|(position, argument)| {
                let AtomView::Var(variable) = *argument else {
                    return None;
                };
                let symbol = variable.get_symbol();
                Self::is_chain_placeholder(*argument).then_some((position, symbol))
            })
            .collect::<Vec<_>>();
        if !placeholders.is_empty()
            && (!matches!(placeholders.as_slice(), [(_, input), (_, output)] if *input != *output)
                || !placeholders
                    .iter()
                    .any(|(_, symbol)| *symbol == SPENSO_TAG.chain_in)
                || !placeholders
                    .iter()
                    .any(|(_, symbol)| *symbol == SPENSO_TAG.chain_out))
        {
            return Err(invalid(factory, signature));
        }
        let Some(expected) = expected else {
            return Ok(None);
        };
        let expected_ports = expected
            .canonical()
            .external_reps_iter()
            .collect::<Vec<_>>();
        if ports.len() != expected_ports.len()
            || ports
                .iter()
                .zip(&expected_ports)
                .any(|(actual, expected)| actual.is_some_and(|actual| actual != *expected))
        {
            return Err(invalid(factory, signature));
        }
        if !placeholders.is_empty() {
            let input = placeholders
                .iter()
                .find(|(_, symbol)| *symbol == SPENSO_TAG.chain_in)
                .expect("placeholder names were validated")
                .0;
            let output = placeholders
                .iter()
                .find(|(_, symbol)| *symbol == SPENSO_TAG.chain_out)
                .expect("placeholder names were validated")
                .0;
            let input = expected_ports[input];
            let output = expected_ports[output];
            let oriented =
                input.rep.is_self_dual() || (input.rep.is_base() && output.rep.is_dual());
            if !input.matches(&output) || !oriented {
                return Err(invalid(factory, signature));
            }
        }
        Ok(Some(expected))
    }

    fn validate_builtin_placeholder_channels(
        &mut self,
        value: AtomView<'_>,
        input: Representation<LibraryRep>,
        output: Representation<LibraryRep>,
    ) -> InferenceResult<()> {
        match value {
            AtomView::Add(add) => {
                for term in add.iter() {
                    self.validate_builtin_placeholder_channels(term, input, output)?;
                }
            }
            AtomView::Mul(mul) => {
                for factor in mul.iter() {
                    self.validate_builtin_placeholder_channels(factor, input, output)?;
                }
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                self.validate_builtin_placeholder_channels(base, input, output)?;
                self.validate_builtin_placeholder_channels(exponent, input, output)?;
            }
            AtomView::Fun(function) => {
                let arguments = function.iter().collect::<Vec<_>>();
                if arguments
                    .iter()
                    .any(|argument| Self::is_chain_placeholder(*argument))
                {
                    let resolved = arguments
                        .iter()
                        .map(|argument| match *argument {
                            AtomView::Var(variable)
                                if variable.get_symbol() == SPENSO_TAG.chain_in =>
                            {
                                input.to_symbolic([])
                            }
                            AtomView::Var(variable)
                                if variable.get_symbol() == SPENSO_TAG.chain_out =>
                            {
                                output.to_symbolic([])
                            }
                            argument => argument.to_owned(),
                        })
                        .collect::<Vec<_>>();
                    let resolved = resolved.iter().map(Atom::as_view).collect::<Vec<_>>();
                    self.builtin_tensor_structure(function.get_symbol(), &resolved)?;
                }
                for argument in arguments {
                    if !Self::is_chain_placeholder(argument) {
                        self.validate_builtin_placeholder_channels(argument, input, output)?;
                    }
                }
            }
            _ => {}
        }
        Ok(())
    }
}

impl InterfaceInference {
    pub fn infer_validated(&mut self, atom: AtomView<'_>) -> InferenceResult<PartialStructure> {
        #[cfg(test)]
        tests::INFERENCE_CALLS.with(|count| count.set(count.get() + 1));
        atom.validate_chain_like_nesting()
            .map_err(|error| TensorInferenceError::invalid(error.to_string()))?;
        Self::validate_placeholder_scope(atom, false).map_err(TensorInferenceError::invalid)?;

        let mut syntax_error = None;
        let mut checked_functions = HashSet::new();
        atom.visitor(&mut |value| {
        if syntax_error.is_some() {
            return false;
        }

        match value {
            AtomView::Var(variable) if variable.get_symbol().has_tag(&SPENSO_TAG.tensor) => {
                syntax_error = Some(TensorInferenceError::invalid(format!(
                    "tensor symbol `{}` must be called with at least one structural port",
                    variable.get_symbol()
                )));
            }
            AtomView::Fun(function) => {
                if value.get_byte_size() <= InterfaceInference::CACHE_KEY_BYTES
                    && checked_functions.contains(&value)
                {
                    return false;
                }
                let symbol = function.get_symbol();
                // Explicit scalar functions own their metadata, which may itself
                // contain tensor syntax without exposing any tensor ports.
                if symbol.is_scalar() {
                    return false;
                }
                let arguments = function.iter().collect::<Vec<_>>();
                // A representation carrying an index must be a valid slot. Falling back
                // to its bare representation would turn contracted indices into open ports.
                if symbol.has_tag(&SPENSO_TAG.representation)
                    && arguments.len() > 1
                    && (arguments.len() != 2
                        || self.slots.parse::<LibraryRep, AbstractIndex>(value).is_err())
                {
                    syntax_error = Some(TensorInferenceError::invalid(format!(
                        "invalid indexed representation `{value}`; cook nested indices before inferring a tensor interface"
                    )));
                    return false;
                }
                if symbol == SPENSO_TAG.bracket && arguments.is_empty() {
                    syntax_error = Some(TensorInferenceError::invalid("ordered tensor products require at least one operand"));
                    return false;
                }
                if symbol == SPENSO_TAG.dot && arguments.len() != 2 {
                    syntax_error = Some(TensorInferenceError::invalid("dot requires exactly two operands"));
                    return false;
                }
                if symbol == SPENSO_TAG.chain && arguments.len() < 2 {
                    syntax_error = Some(TensorInferenceError::invalid("chain requires explicit start and end ports"));
                    return false;
                }
                if symbol == SPENSO_TAG.trace {
                    let Some(representation) = arguments.first() else {
                        syntax_error = Some(TensorInferenceError::invalid("trace requires a representation"));
                        return false;
                    };
                    if self.slots.parse_representation::<LibraryRep>(*representation).is_err() {
                        syntax_error =
                            Some(TensorInferenceError::invalid("trace metadata is not a Spenso representation"));
                        return false;
                    }
                }
                if symbol.has_tag(&SPENSO_TAG.broadcast) && arguments.len() != 1 {
                    syntax_error = Some(TensorInferenceError::invalid(format!(
                        "broadcast function `{symbol}` requires exactly one argument"
                    )));
                    return false;
                }
                if !symbol.has_tag(&SPENSO_TAG.tensor) {
                    return true;
                }

                let compact_metric = symbol == ETS.metric
                    && arguments.len() == 2
                    && arguments.iter().all(|argument| {
                        self.slots.parse::<LibraryRep, AbstractIndex>(*argument).is_err()
                            && self.slots.parse_representation::<LibraryRep>(*argument).is_err()
                            && Self::has_structured_syntax(*argument)
                    });
                let builtin = match self.builtin_tensor_structure(symbol, &arguments) {
                    Ok(builtin) => builtin,
                    Err(error) => {
                        syntax_error = Some(error);
                        return false;
                    }
                };

                let ports = arguments
                    .iter()
                    .enumerate()
                    .filter(|(_, argument)| {
                        self.slots.parse::<LibraryRep, AbstractIndex>(**argument).is_ok()
                            || self.slots.parse_representation::<LibraryRep>(**argument).is_ok()
                            || matches!(
                                **argument,
                                AtomView::Var(variable)
                                    if variable.get_symbol() == SPENSO_TAG.chain_in
                                        || variable.get_symbol() == SPENSO_TAG.chain_out
                            )
                    })
                    .map(|(position, _)| position)
                    .collect::<Vec<_>>();
                let ports_are_final = ports
                    .iter()
                    .copied()
                    .eq(arguments.len().saturating_sub(ports.len())..arguments.len());
                if !ports_are_final && !compact_metric && builtin.is_none() {
                    syntax_error = Some(TensorInferenceError::invalid(format!(
                        "tensor function `{symbol}` requires scalar arguments before structural ports"
                    )));
                } else if symbol.has_tag(&SPENSO_TAG.rank1) && ports.len() != 1 {
                    syntax_error = Some(TensorInferenceError::invalid(format!(
                        "rank-one tensor function `{symbol}` requires exactly one final structural port"
                    )));
                }
                // Explicit built-ins and directly readable leaves have no
                // materialization effects during syntax validation. The first
                // visit still checks every child before a later occurrence is skipped.
                if syntax_error.is_none()
                    && checked_functions.len() < InterfaceInference::CACHE_ENTRIES
                    && value.get_byte_size() <= InterfaceInference::CACHE_KEY_BYTES
                    && ((builtin.is_some() && arguments.iter().all(|argument| {
                            self.slots.parse::<LibraryRep, AbstractIndex>(*argument)
                                .is_ok_and(|slot| !matches!(slot.aind(), AbstractIndex::Open { .. }))
                        }))
                        || self.direct_leaf_interface(function).is_some())
                {
                    checked_functions.insert(value);
                }
            }
            _ => {}
        }
        true
    });
        if let Some(error) = syntax_error {
            return Err(error);
        }

        self.infer_view(atom)
    }
}

#[derive(Default)]
pub struct InterfaceInference {
    reusable_interfaces: HashMap<Vec<u8>, PartialStructure>,
    slots: SlotMatcher,
    leaf_inference: LeafInference,
}

#[derive(Default, PartialEq)]
enum LeafInference {
    #[default]
    Materialize,
    Observe,
}

impl InterfaceInference {
    const CACHE_ENTRIES: usize = 256;
    const CACHE_KEY_BYTES: usize = 256;

    /// Call only for interfaces independent of materialization callbacks and
    /// fresh dummy identities. Open ports remain occurrence-local positions.
    fn cache_interface(&mut self, atom: AtomView<'_>, interface: &PartialStructure) {
        if self.reusable_interfaces.len() < Self::CACHE_ENTRIES
            && atom.get_byte_size() <= Self::CACHE_KEY_BYTES
        {
            self.reusable_interfaces
                .insert(atom.get_data().to_vec(), interface.clone());
        }
    }

    #[cfg(test)]
    fn infer(&mut self, atom: &Atom) -> InferenceResult<PartialStructure> {
        self.infer_view(atom.as_view())
    }

    /// Read ordinary leaf ports in their existing syntactic order. Compact
    /// representations denote occurrence-local open ports; no dummy atom or
    /// function reconstruction is needed to discover their representation.
    fn direct_leaf_interface(
        &mut self,
        function: symbolica::atom::representation::FunView<'_>,
    ) -> Option<PartialStructure> {
        let symbol = function.get_symbol();
        if !matches!(self.slots.classify(function.as_view()), SlotMatch::Other)
            || (self.leaf_inference == LeafInference::Materialize
                && (symbol.get_normalization_function().is_some()
                    || symbol.get_evaluation_info().is_some()
                    || symbol.is_symmetric()
                    || symbol.is_antisymmetric()
                    || symbol.is_cyclesymmetric()
                    || symbol.is_linear()))
        {
            return None;
        }
        let mut logical = Vec::new();
        let mut explicit_ports = 0;
        let mut function_metadata = false;
        for argument in function.iter() {
            if let Ok(slot) = self.slots.parse::<LibraryRep, AbstractIndex>(argument) {
                explicit_ports += 1;
                let index = match slot.aind() {
                    AbstractIndex::Open { axis, .. } => PartialIndex::open(axis),
                    index => PartialIndex::Explicit(index),
                };
                logical.push(slot.rep().slot(index));
            } else if let Ok(rep) = self.slots.parse_representation::<LibraryRep>(argument) {
                logical.push(rep.slot(PartialIndex::open(logical.len())));
            } else {
                // Nested function metadata may contain compact representations,
                // bundles, or callbacks affected by legacy materialization.
                // Keep that existing owner for those forms.
                let mut contains_function = false;
                argument.visitor(&mut |node| {
                    contains_function |= matches!(node, AtomView::Fun(_));
                    !contains_function
                });
                if contains_function && self.leaf_inference == LeafInference::Materialize {
                    return None;
                }
                function_metadata |= contains_function;
            }
        }
        if function_metadata {
            // Observation may ignore scalar metadata, but must not hide a
            // structural bundle that the existing fast leaf parser exposes.
            let inferred = match function
                .as_view()
                .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>(
                    StructureInferenceMode::Fast,
                ) {
                Ok(inferred) => inferred.canonical().order(),
                Err(StructureError::EmptyStructure(_)) => 0,
                Err(_) => return None,
            };
            if inferred != explicit_ports {
                return None;
            }
        }
        Some(PartialStructure::from_logical_slots(logical))
    }

    fn infer_view(&mut self, atom: AtomView<'_>) -> InferenceResult<PartialStructure> {
        if matches!(atom, AtomView::Fun(_))
            && atom.get_byte_size() <= Self::CACHE_KEY_BYTES
            && let Some(interface) = self.reusable_interfaces.get(atom.get_data())
        {
            return Ok(interface.clone());
        }

        if let AtomView::Add(sum) = atom {
            let mut expected = None;
            let mut has_scalar_term = false;
            for term in sum.iter() {
                if !Self::has_structured_syntax(term) {
                    has_scalar_term = true;
                    continue;
                }
                let actual = Self::merge_explicit_interface_sequence(&[self.infer_view(term)?])?;
                let Some(current) = &expected else {
                    expected = Some(actual);
                    continue;
                };
                if !Self::additive_interfaces_match(current, &actual) {
                    return Err(TensorInferenceError::invalid(
                        "tensor summands do not have compatible tensor interfaces",
                    ));
                }
            }
            let Some(expected) = expected else {
                return Err(TensorInferenceError::invalid(
                    "expression does not contain valid tagged Spenso tensor syntax",
                ));
            };
            if has_scalar_term && !expected.canonical().is_scalar() {
                return Err(TensorInferenceError::invalid(
                    "cannot add a scalar expression to a non-scalar tensor",
                ));
            }
            return Ok(expected);
        }

        if let AtomView::Mul(product) = atom {
            let mut interfaces = Vec::new();
            for factor in product.iter() {
                if !Self::has_structured_syntax(factor) {
                    continue;
                }
                interfaces.push(self.infer_view(factor)?);
            }
            if interfaces.is_empty() {
                return Err(TensorInferenceError::invalid(
                    "expression does not contain valid tagged Spenso tensor syntax",
                ));
            }
            return Self::merge_explicit_interface_sequence(&interfaces);
        }

        if let AtomView::Pow(power) = atom {
            let (base, exponent) = power.get_base_exp();
            if !Self::has_structured_syntax(base) {
                return Err(TensorInferenceError::invalid(
                    "expression does not contain valid tagged Spenso tensor syntax",
                ));
            }
            let interface = self.infer_view(base)?;
            if Self::has_structured_syntax(exponent)
                && !self.infer_view(exponent)?.canonical().is_scalar()
            {
                return Err(TensorInferenceError::invalid(
                    "a tensor exponent must be scalar",
                ));
            }
            if interface.canonical().is_scalar() {
                return Ok(interface);
            }
            if !interface
                .logical_slots()
                .iter()
                .all(|slot| slot.rep().rep.is_self_dual())
            {
                return Err(TensorInferenceError::invalid(format!(
                    "invalid power of non-self-dual tensor `{atom}`"
                )));
            }
            let exponent = Rational::try_from(exponent).map_err(|_| {
                TensorInferenceError::invalid(format!("invalid tensor power `{atom}`"))
            })?;
            if exponent.denominator() != 1 {
                return Err(TensorInferenceError::invalid(format!(
                    "fractional tensor power `{atom}` has no well-defined interface"
                )));
            }
            let repetitions = exponent.numerator().abs();
            let slots = interface.logical_slots();
            for slot in &slots {
                let PartialIndex::Explicit(index) = slot.aind else {
                    continue;
                };
                let compatible = slots
                .iter()
                .filter(|candidate| {
                    matches!(candidate.aind, PartialIndex::Explicit(candidate_index) if candidate_index == index)
                        && slot.rep().matches(&candidate.rep())
                })
                .count();
                let exceeds_einstein_multiplicity = match repetitions.to_i64() {
                    Some(0) => false,
                    Some(1) => compatible > 2,
                    Some(2) => compatible > 1,
                    _ => compatible > 0,
                };
                if exceeds_einstein_multiplicity {
                    return Err(TensorInferenceError::invalid(format!(
                        "explicit index `{index}` occurs on more than two compatible ports in `{atom}`"
                    )));
                }
            }
            if exponent.numerator() % 2 == 0 {
                return Ok(PartialStructure::from_logical_slots([]));
            }
            return Ok(interface);
        }

        if let AtomView::Fun(function) = atom {
            let symbol = function.get_symbol();
            let arguments = function.iter().collect::<Vec<_>>();
            let unmaterialized_tensor_leaf = Self::is_unmaterialized_tensor_leaf(atom);
            if symbol == *shadowing::SYM
                || symbol == *shadowing::ANTISYM
                || symbol == *shadowing::CYCLIC
            {
                let mut interfaces = Vec::new();
                for argument in arguments {
                    if !Self::has_structured_syntax(argument) {
                        continue;
                    }
                    interfaces.push(self.infer_view(argument)?);
                }
                if interfaces.is_empty() {
                    return Err(TensorInferenceError::invalid(
                        "tensor projectors require at least one tensor operand",
                    ));
                }
                return Self::merge_explicit_interface_sequence(&interfaces);
            }

            if symbol == SPENSO_TAG.bracket {
                let mut interfaces = Vec::new();
                for argument in arguments {
                    if !Self::has_structured_syntax(argument) {
                        continue;
                    }
                    interfaces.push(self.infer_view(argument)?);
                }
                if interfaces.is_empty() {
                    return Err(TensorInferenceError::invalid(
                        "ordered tensor products require at least one tensor operand",
                    ));
                }
                return Self::merge_explicit_interface_sequence(&interfaces);
            }

            let compact_metric = symbol == ETS.metric
                && arguments.len() == 2
                && arguments.iter().all(|argument| {
                    self.slots
                        .parse::<LibraryRep, AbstractIndex>(*argument)
                        .is_err()
                        && self
                            .slots
                            .parse_representation::<LibraryRep>(*argument)
                            .is_err()
                        && !Self::is_chain_placeholder(*argument)
                        && Self::has_structured_syntax(*argument)
                });
            if symbol == SPENSO_TAG.dot || compact_metric {
                let [left, right] = arguments.as_slice() else {
                    return Err(TensorInferenceError::invalid(
                        "inner products require exactly two operands",
                    ));
                };
                if !Self::has_structured_syntax(*left) || !Self::has_structured_syntax(*right) {
                    return Err(TensorInferenceError::invalid(
                        "inner products require rank-one tensor operands",
                    ));
                }
                let left = self.infer_view(*left)?;
                let right = self.infer_view(*right)?;
                let reusable = arguments
                    .iter()
                    .all(|argument| self.reusable_interfaces.contains_key(argument.get_data()));
                if left.canonical().order() != 1 || right.canonical().order() != 1 {
                    return Err(TensorInferenceError::invalid(format!(
                        "inner products require rank-one operands, got ranks {} and {}",
                        left.canonical().order(),
                        right.canonical().order()
                    )));
                }
                let left = left.logical_slots()[0];
                let right = right.logical_slots()[0];
                if !left.rep().matches(&right.rep()) {
                    return Err(TensorInferenceError::invalid(
                        "inner-product operands carry incompatible representations",
                    ));
                }
                if matches!(
                    (left.aind, right.aind),
                    (PartialIndex::Explicit(left), PartialIndex::Explicit(right)) if left != right
                ) {
                    return Err(TensorInferenceError::invalid(
                        "inner-product operands carry unequal explicit indices",
                    ));
                }
                let interface = PartialStructure::from_logical_slots([]);
                if reusable {
                    // A dot consumes both ports. Cached operands establish that
                    // repeating their inference has no materialization effects.
                    self.cache_interface(atom, &interface);
                }
                return Ok(interface);
            }

            if symbol.has_tag(&SPENSO_TAG.broadcast) {
                let argument = arguments[0];
                if !Self::has_structured_syntax(argument) {
                    return Err(TensorInferenceError::invalid(format!(
                        "broadcast function `{symbol}` does not contain a structured tensor argument"
                    )));
                }
                return self.infer_view(argument);
            }

            if let Some(structure) = self.builtin_tensor_structure(symbol, &arguments)? {
                let canonical_ports = arguments
                    .iter()
                    .enumerate()
                    .map(|(position, argument)| {
                        if let Ok(slot) = self.slots.parse::<LibraryRep, AbstractIndex>(*argument) {
                            let index = match slot.aind() {
                                AbstractIndex::Open { axis, .. } => PartialIndex::open(axis),
                                index => PartialIndex::Explicit(index),
                            };
                            Some(slot.rep().slot(index))
                        } else {
                            self.slots
                                .parse_representation::<LibraryRep>(*argument)
                                .ok()
                                .map(|representation| {
                                    representation.slot(PartialIndex::open(position))
                                })
                        }
                    })
                    .collect::<Vec<_>>();
                // Only fully explicit built-in leaves are independent of fresh
                // port materialization and normalization callbacks.
                let reusable = canonical_ports.iter().all(|port| {
                    matches!(port, Some(slot) if matches!(slot.aind, PartialIndex::Explicit(_)))
                });
                // Compact arguments are contracted vectors, and chain placeholders
                // are wiring labels. Neither exposes an external tensor port.
                let interface = PartialStructure::from_logical_slots(
                    structure
                        .layout()
                        .canonical_to_logical(&canonical_ports)
                        .into_iter()
                        .flatten(),
                );
                if reusable {
                    // Keep repeated indices and logical ordering intact: enclosing
                    // products still own contraction and multiplicity validation.
                    self.cache_interface(atom, &interface);
                }
                return Ok(interface);
            }

            if symbol == SPENSO_TAG.chain {
                let endpoints = arguments[..2]
                    .iter()
                    .map(|argument| {
                        if let Ok(slot) = self.slots.parse::<LibraryRep, AbstractIndex>(*argument) {
                            Ok(slot.rep().slot(PartialIndex::Explicit(slot.aind())))
                        } else if let Ok(representation) =
                            self.slots.parse_representation::<LibraryRep>(*argument)
                        {
                            Ok(representation.slot(PartialIndex::open(0)))
                        } else {
                            Err(TensorInferenceError::invalid(
                                "chain endpoints must be Spenso slots or representations",
                            ))
                        }
                    })
                    .collect::<InferenceResult<Vec<_>>>()?;
                let input = endpoints[0].rep();
                let output = endpoints[1].rep();
                if !input.matches(&output)
                    || !(input.rep.is_self_dual() || (input.rep.is_base() && output.rep.is_dual()))
                {
                    return Err(TensorInferenceError::invalid(
                        "chain endpoints do not form a compatible input-to-output channel",
                    ));
                }
                let mut interfaces = Vec::new();
                for factor in &arguments[2..] {
                    if !Self::has_structured_syntax(*factor) {
                        return Err(TensorInferenceError::invalid(
                            "chain factors must contain structured tensors",
                        ));
                    }
                    self.validate_builtin_placeholder_channels(*factor, input, output)?;
                    interfaces.push(self.infer_view(*factor)?);
                }
                let spectators = Self::merge_explicit_interface_sequence(&interfaces)?;
                return Ok(PartialStructure::from_logical_slots(
                    endpoints.into_iter().chain(spectators.logical_slots()),
                ));
            } else if symbol == SPENSO_TAG.trace {
                let representation = self
                    .slots
                    .parse_representation::<LibraryRep>(arguments[0])
                    .map_err(|_| {
                        TensorInferenceError::invalid(
                            "trace metadata is not a Spenso representation",
                        )
                    })?;
                let dual = representation.dual();
                let mut interfaces = Vec::new();
                for factor in shadowing::trace_factor_views(&arguments[1..]) {
                    if !Self::has_structured_syntax(factor) {
                        return Err(TensorInferenceError::invalid(
                            "trace factors must contain structured tensors",
                        ));
                    }
                    self.validate_builtin_placeholder_channels(factor, representation, dual)?;
                    interfaces.push(self.infer_view(factor)?);
                }
                return Self::merge_explicit_interface_sequence(&interfaces);
            } else if symbol.has_tag(&SPENSO_TAG.tensor) {
                for argument in arguments {
                    if self
                        .slots
                        .parse::<LibraryRep, AbstractIndex>(argument)
                        .is_ok()
                        || self
                            .slots
                            .parse_representation::<LibraryRep>(argument)
                            .is_ok()
                        || matches!(
                            argument,
                            AtomView::Var(variable)
                                if variable.get_symbol() == SPENSO_TAG.chain_in
                                    || variable.get_symbol() == SPENSO_TAG.chain_out
                        )
                        || !Self::has_structured_syntax(argument)
                    {
                        continue;
                    }
                    if !self.infer_view(argument)?.canonical().is_scalar() {
                        return Err(TensorInferenceError::invalid(format!(
                            "tensor metadata for `{symbol}` must be scalar"
                        )));
                    }
                }
                if unmaterialized_tensor_leaf {
                    return Ok(PartialStructure::from_logical_slots([]));
                }
            }
        }

        if !Self::has_structured_syntax(atom) {
            return Err(TensorInferenceError::invalid(
                "expression does not contain valid tagged Spenso tensor syntax",
            ));
        }

        if let AtomView::Fun(function) = atom
            && let Some(interface) = self.direct_leaf_interface(function)
        {
            // OpenPortIds are local positions, canonicalized when interfaces are
            // combined. Direct parsing allocates no dummy or callback-sensitive state.
            self.cache_interface(atom, &interface);
            return Ok(interface);
        }

        let atom = atom.to_owned();
        let mut open_markers = HashSet::new();
        // Constructors retain their historical compact-port materialization.
        // A result check observes the already normalized expression instead:
        // replaying a leaf callback with a fresh dummy could change its rank.
        let materialized = if self.leaf_inference == LeafInference::Observe {
            atom.clone()
        } else {
            atom.replace_map(|value, _, output| {
                if self.slots.parse::<LibraryRep, AbstractIndex>(value).is_ok() {
                    return;
                }
                if let Ok(representation) = self.slots.parse_representation::<LibraryRep>(value) {
                    let marker = loop {
                        let marker = composition::fresh_dummy_index([&atom], std::iter::empty());
                        if open_markers.insert(marker) {
                            break marker;
                        }
                    };
                    **output = representation.slot::<AbstractIndex, _>(marker).to_atom();
                }
            })
        };
        let inferred = match materialized
            .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>(
                StructureInferenceMode::Fast,
            ) {
            Ok(inferred) => inferred,
            Err(StructureError::EmptyStructure(_)) => {
                return Ok(PartialStructure::from_logical_slots([]));
            }
            Err(error) => {
                return Err(TensorInferenceError::invalid(format!(
                    "invalid Spenso expression: {error}"
                )));
            }
        };
        // `OrderedStructure` canonicalizes slots and its fast function inference
        // does not retain that permutation. Read the direct structural arguments
        // back from the atom so ordinary non-built-in leaves follow their encoded
        // syntax order.
        let logical = Self::syntactic_leaf_slots(materialized.as_view());
        if logical.len() != inferred.canonical().order() {
            return Err(TensorInferenceError::invalid(format!(
                "invalid Spenso expression: inferred {} ports but found {} direct structural arguments",
                inferred.canonical().order(),
                logical.len()
            )));
        }
        Ok(PartialStructure::from_logical_slots(
            logical.into_iter().map(|slot| {
                let index = match slot.aind() {
                    AbstractIndex::Open { axis, .. } => PartialIndex::open(axis),
                    index if open_markers.contains(&index) => PartialIndex::open(0),
                    index => PartialIndex::Explicit(index),
                };
                slot.rep().slot(index)
            }),
        ))
    }
}

impl InterfaceInference {
    pub fn additive_interfaces_match(left: &PartialStructure, right: &PartialStructure) -> bool {
        // Unresolved ports remain positional; only explicit indices identify ports across terms.
        left.logical_slots() == right.logical_slots()
            || (left.open_positions().is_empty()
                && right.open_positions().is_empty()
                && left.canonical() == right.canonical())
    }

    pub fn has_structured_syntax(value: AtomView<'_>) -> bool {
        value.is_tensorial(StrictTensorFilter::Tagged)
            || value.is_tensorial(StrictTensorFilter::ContainsReps)
            || Self::is_unmaterialized_tensor_leaf(value)
            || matches!(
                value,
                AtomView::Fun(projector)
                    if projector.get_symbol() == *shadowing::SYM
                        || projector.get_symbol() == *shadowing::ANTISYM
                        || projector.get_symbol() == *shadowing::CYCLIC
            )
    }

    fn is_unmaterialized_tensor_leaf(value: AtomView<'_>) -> bool {
        let AtomView::Fun(function) = value else {
            return false;
        };
        function.get_symbol().has_tag(&SPENSO_TAG.tensor)
            && function.iter().any(Self::is_chain_placeholder)
            && function.iter().all(|argument| {
                Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_err()
                    && Representation::<LibraryRep>::try_from(argument).is_err()
            })
    }

    fn is_chain_placeholder(value: AtomView<'_>) -> bool {
        matches!(
            value,
            AtomView::Var(variable)
                if variable.get_symbol() == SPENSO_TAG.chain_in
                    || variable.get_symbol() == SPENSO_TAG.chain_out
        )
    }

    /// Namespaced `in`/`out` placeholders are wiring labels, not free abstract indices. They may
    /// only occur as direct ports of tensor leaves inside a chain or trace factor.
    pub fn validate_placeholder_scope(
        value: AtomView<'_>,
        inside_factor: bool,
    ) -> Result<(), String> {
        match value {
            AtomView::Var(_) if Self::is_chain_placeholder(value) => Err(
                "Spenso chain placeholders are only valid inside chain or trace tensor factors"
                    .into(),
            ),
            AtomView::Add(sum) => {
                for term in sum.iter() {
                    Self::validate_placeholder_scope(term, inside_factor)?;
                }
                Ok(())
            }
            AtomView::Mul(product) => {
                for factor in product.iter() {
                    Self::validate_placeholder_scope(factor, inside_factor)?;
                }
                Ok(())
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                Self::validate_placeholder_scope(base, inside_factor)?;
                Self::validate_placeholder_scope(exponent, inside_factor)
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                if symbol.is_scalar() {
                    return Ok(());
                }
                if symbol == SPENSO_TAG.chain {
                    for (position, argument) in function.iter().enumerate() {
                        Self::validate_placeholder_scope(argument, position >= 2)?;
                    }
                    return Ok(());
                }
                if symbol == SPENSO_TAG.trace {
                    for (position, argument) in function.iter().enumerate() {
                        Self::validate_placeholder_scope(argument, position >= 1)?;
                    }
                    return Ok(());
                }

                let tensor_leaf = symbol.has_tag(&SPENSO_TAG.tensor);
                let mut inputs = 0;
                let mut outputs = 0;
                for argument in function.iter() {
                    if Self::is_chain_placeholder(argument) {
                        if !(inside_factor && tensor_leaf) {
                            return Err(
                            "Spenso chain placeholders are only valid as tensor ports inside chain or trace factors"
                                .into(),
                        );
                        }
                        let AtomView::Var(variable) = argument else {
                            unreachable!("chain placeholders are variables")
                        };
                        if variable.get_symbol() == SPENSO_TAG.chain_in {
                            inputs += 1;
                        } else {
                            outputs += 1;
                        }
                    } else {
                        Self::validate_placeholder_scope(argument, inside_factor)?;
                    }
                }
                if inputs + outputs > 0 && (inputs != 1 || outputs != 1) {
                    return Err(format!(
                        "tensor factor `{symbol}` requires exactly one direct `in` and one direct `out` placeholder"
                    ));
                }
                Ok(())
            }
            _ => Ok(()),
        }
    }

    /// Lower tensor powers, returning no replacement when the normalized atom is unchanged.
    pub fn lower_tensor_powers(value: AtomView<'_>) -> InferenceResult<Option<Atom>> {
        let lowered = match value {
            AtomView::Add(sum) => {
                let terms = sum
                    .iter()
                    .map(Self::lower_tensor_powers)
                    .collect::<InferenceResult<Vec<_>>>()?;
                if terms.iter().all(Option::is_none) {
                    return Ok(None);
                }
                sum.iter()
                    .zip(&terms)
                    .map(|(old, new)| new.as_ref().map_or(old, Atom::as_view))
                    .sum()
            }
            AtomView::Mul(product) => {
                let factors = product
                    .iter()
                    .map(Self::lower_tensor_powers)
                    .collect::<InferenceResult<Vec<_>>>()?;
                if factors.iter().all(Option::is_none) {
                    return Ok(None);
                }
                product
                    .iter()
                    .zip(&factors)
                    .map(|(old, new)| new.as_ref().map_or(old, Atom::as_view))
                    .product()
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let lowered_base = Self::lower_tensor_powers(base)?;
                let lowered_exponent = Self::lower_tensor_powers(exponent)?;
                let base = lowered_base.as_ref().map_or(base, Atom::as_view);
                let exponent = lowered_exponent.as_ref().map_or(exponent, Atom::as_view);
                if !Self::has_structured_syntax(base) {
                    if lowered_base.is_none() && lowered_exponent.is_none() {
                        return Ok(None);
                    }
                    base.pow(exponent)
                } else {
                    // Validate multiplicity before lowering through tensor-aware multiplication.
                    let mut inference = InterfaceInference::default();
                    inference.infer_view(base.pow(exponent).as_view())?;
                    let raw_interface = inference.infer_validated(base)?;
                    let interface = Self::merge_explicit_interface_sequence(&[raw_interface])?;
                    if interface.canonical().is_scalar() {
                        base.pow(exponent)
                    } else {
                        let exponent = Rational::try_from(exponent).map_err(|_| {
                            TensorInferenceError::invalid(
                                "a non-scalar tensor power requires an integer exponent",
                            )
                        })?;
                        if exponent.denominator() != 1 || exponent.numerator().is_negative() {
                            return Err(TensorInferenceError::invalid(
                                "a non-scalar tensor power requires a non-negative integer exponent",
                            ));
                        }
                        let repetitions =
                            usize::try_from(exponent.numerator().clone()).map_err(|_| {
                                TensorInferenceError::invalid(
                                    "tensor power exponent is too large to materialize",
                                )
                            })?;
                        if repetitions == 0 {
                            Atom::num(1)
                        } else if interface.open_positions().is_empty() {
                            // Explicit labels already identify the contractions; only open
                            // ports need the pair selection of tensor multiplication.
                            base.pow(Atom::num(repetitions))
                        } else {
                            let factor =
                                SymbolicTensor::<PartialStructure>::new(base.to_owned(), interface);
                            let mut result = factor.clone();
                            for _ in 1..repetitions {
                                result = result.multiply(&factor)?;
                            }
                            result.expression
                        }
                    }
                }
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                if symbol.is_scalar() {
                    return Ok(None);
                }
                let arguments = function
                    .iter()
                    .map(Self::lower_tensor_powers)
                    .collect::<InferenceResult<Vec<_>>>()?;
                if arguments.iter().all(Option::is_none)
                    && symbol.get_normalization_function().is_none()
                    && symbol.get_evaluation_info().is_none()
                {
                    return Ok(None);
                }
                // Preserve existing callback invocation when rebuilding is observable.
                FunctionBuilder::new(symbol)
                    .add_args(
                        function
                            .iter()
                            .zip(&arguments)
                            .map(|(old, new)| new.as_ref().map_or(old, Atom::as_view)),
                    )
                    .finish()
            }
            _ => return Ok(None),
        };
        Ok((lowered.as_view() != value).then_some(lowered))
    }

    fn syntactic_leaf_slots(value: AtomView<'_>) -> Vec<Slot<LibraryRep, AbstractIndex>> {
        if let Ok(slot) = Slot::<LibraryRep, AbstractIndex>::try_from(value) {
            return vec![slot];
        }
        let AtomView::Fun(function) = value else {
            return Vec::new();
        };
        function
            .iter()
            .filter_map(|argument| Slot::<LibraryRep, AbstractIndex>::try_from(argument).ok())
            .collect()
    }

    pub fn merge_explicit_interface_sequence(
        interfaces: &[PartialStructure],
    ) -> InferenceResult<PartialStructure> {
        let slots = interfaces
            .iter()
            .flat_map(PartialStructureExt::logical_slots)
            .collect::<Vec<_>>();
        let mut contracted = HashSet::new();
        for (position, slot) in slots.iter().enumerate() {
            let PartialIndex::Explicit(index) = slot.aind else {
                continue;
            };
            let compatible = slots
            .iter()
            .enumerate()
            .filter(|(candidate_position, candidate)| {
                if *candidate_position == position {
                    return true;
                }
                matches!(candidate.aind, PartialIndex::Explicit(candidate_index) if candidate_index == index)
                    && slot.rep().matches(&candidate.rep())
            })
            .map(|(candidate_position, _)| candidate_position)
            .collect::<Vec<_>>();
            if compatible.len() > 2 {
                return Err(TensorInferenceError::invalid(format!(
                    "explicit index `{index}` occurs on more than two compatible ports {compatible:?}"
                )));
            }
            if compatible.len() == 2 {
                contracted.extend(compatible);
            }
        }

        Ok(PartialStructure::from_logical_slots(
            slots
                .into_iter()
                .enumerate()
                .filter(|(position, _)| !contracted.contains(position))
                .map(|(_, slot)| slot),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shorthands::metric::MetricSimplifier;
    use spenso::structure::representation::ExtendibleReps;

    thread_local! {
        pub(super) static INFERENCE_CALLS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    }

    fn explicit_interface(
        representation: Representation<LibraryRep>,
        index: AbstractIndex,
    ) -> PartialStructure {
        PartialStructure::from_logical_slots([representation.slot(PartialIndex::Explicit(index))])
    }

    #[test]
    fn builtin_signature_errors_remain_structured_through_syntax_validation() {
        let atom = FunctionBuilder::new(AGS.gamma)
            .add_arg(Atom::one())
            .finish();
        assert!(matches!(
            InterfaceInference::default().infer_validated(atom.as_view()),
            Err(TensorInferenceError::InvalidBuiltinSignature {
                factory: "gamma",
                ..
            })
        ));
    }

    #[test]
    fn signature_construction_retains_logical_order_and_allocates_distinct_occurrences() {
        crate::representations::initialize();
        let signature = CS.t_strct(Dimension::Concrete(3), Dimension::Concrete(8));
        let first = SymbolicTensor::<PartialStructure>::from_signature(&signature).unwrap();
        let second = SymbolicTensor::<PartialStructure>::from_signature(&signature).unwrap();
        assert_ne!(first.expression, second.expression);
        assert_eq!(first.structure, second.structure);
        assert_eq!(first.rank(), 3);
        assert_eq!(first.structure.open_positions(), vec![0, 1, 2]);
        assert_eq!(
            first
                .structure
                .logical_slots()
                .iter()
                .map(IsAbstractSlot::rep)
                .collect::<Vec<_>>(),
            signature.layout().canonical_to_logical(
                &signature
                    .canonical()
                    .external_reps_iter()
                    .collect::<Vec<_>>()
            ),
        );

        let strips_ports = spenso::tensor_symbol!(
            "signature_constructor_strips_ports",
            norm = |_, output| {
                **output = Atom::one();
            }
        );
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let invalid = ExplicitKey::from_iter([representation], strips_ports, None);
        assert!(SymbolicTensor::<PartialStructure>::from_signature(&invalid).is_err());
    }

    #[test]
    fn mixed_signature_layout_survives_physical_storage_order_and_a_noop_normalizer() {
        let euclidean = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let minkowski = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let plain = spenso::tensor_symbol!("mixed_signature_order");
        let callback = spenso::tensor_symbol!("mixed_signature_order_callback", norm = |_, _| {});
        for name in [plain, callback] {
            let signature = ExplicitKey::from_iter([euclidean, minkowski], name, None);
            let value = SymbolicTensor::<PartialStructure>::from_signature(&signature).unwrap();
            assert_eq!(
                value
                    .structure
                    .logical_slots()
                    .iter()
                    .map(IsAbstractSlot::rep)
                    .collect::<Vec<_>>(),
                vec![euclidean, minkowski],
            );
            assert_eq!(value.structure.open_positions(), vec![0, 1]);
            let dimension = Dimension::Concrete(6);
            let promoted = value.with_lorentz_dimension(dimension).unwrap();
            assert_eq!(
                promoted.structure.logical_slots(),
                vec![
                    euclidean.slot(PartialIndex::open(0)),
                    ExtendibleReps::MINKOWSKI
                        .new_rep(dimension)
                        .slot(PartialIndex::open(1)),
                ],
            );
            assert_eq!(
                promoted.expression,
                value
                    .expression
                    .with_lorentz_dimension(dimension.to_symbolic().as_view()),
            );
        }
    }

    #[test]
    fn lorentz_dimension_preserves_logical_ports_typed_zero_and_dispatch_flags() {
        let minkowski = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let euclidean = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let slots = [
            euclidean.slot(PartialIndex::open(0)),
            minkowski.slot(PartialIndex::Explicit(AbstractIndex::Normal(73301))),
            ExtendibleReps::MINKOWSKI
                .new_rep(Dimension::Concrete(3))
                .slot(PartialIndex::open(1)),
        ];
        let atom = FunctionBuilder::new(spenso::tensor_symbol!("dimension_logical_ports"))
            .add_arg(4)
            .add_args(slots.into_iter().rev().map(composition::port_atom))
            .finish();
        let dimension = Dimension::from(symbolica::symbol!("shared_lorentz_dimension"));
        let expected_slots = [
            slots[0],
            ExtendibleReps::MINKOWSKI
                .new_rep(dimension)
                .slot(slots[1].aind),
            slots[2],
        ];
        for expression in [atom, Atom::Zero] {
            let value = SymbolicTensor::<PartialStructure>::new(
                expression,
                PartialStructure::from_logical_slots(slots),
            );
            INFERENCE_CALLS.with(|count| count.set(0));
            let result = value.with_lorentz_dimension(dimension).unwrap();
            assert_eq!(result.structure.logical_slots(), expected_slots);
            assert_eq!(result.structure.open_positions(), vec![0, 2]);
            assert_eq!(result.is_metric, value.is_metric);
            assert_eq!(result.is_composite, value.is_composite);
            assert_eq!(result.expression.is_zero(), value.expression.is_zero());
            assert_eq!(result.with_lorentz_dimension(dimension).unwrap(), result);
            assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        }

        let gamma = SymbolicTensor::<PartialStructure>::from_signature(
            &AGS.gamma_strct::<AbstractIndex>(Dimension::Concrete(4)),
        )
        .unwrap();
        let promoted = gamma.with_lorentz_dimension(dimension).unwrap();
        assert_eq!(
            promoted
                .structure
                .logical_slots()
                .iter()
                .map(IsAbstractSlot::dim)
                .collect::<Vec<_>>(),
            vec![4.into(), 4.into(), dimension],
        );
        assert_eq!(
            promoted.structure.open_positions(),
            gamma.structure.open_positions()
        );

        let metric =
            SymbolicTensor::<PartialStructure>::new(
                FunctionBuilder::new(ETS.metric)
                    .add_args([73303, 73305].map(|index| {
                        minkowski
                            .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                            .to_atom()
                    }))
                    .finish(),
                PartialStructure::from_logical_slots([73303, 73305].map(|index| {
                    minkowski.slot(PartialIndex::Explicit(AbstractIndex::Normal(index)))
                })),
            );
        let promoted = metric.with_lorentz_dimension(dimension).unwrap();
        assert!(promoted.is_metric);
        assert!(!promoted.is_composite);
    }

    #[test]
    fn lorentz_dimension_merges_new_index_pairs_and_rejects_excess_occurrences() {
        let index = PartialIndex::Explicit(AbstractIndex::Normal(73311));
        let minkowski4 = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let minkowski6 = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(6));
        let remaining = ExtendibleReps::EUCLIDEAN
            .new_rep(Dimension::Concrete(3))
            .slot(PartialIndex::Explicit(AbstractIndex::Normal(73313)));
        let ports = [minkowski4.slot(index), remaining, minkowski6.slot(index)];
        let head = spenso::tensor_symbol!("dimension_collision_tensor");
        let expression = FunctionBuilder::new(head)
            .add_args(ports.map(composition::port_atom))
            .finish();
        let value = SymbolicTensor::<PartialStructure>::checked_parts(
            expression,
            PartialStructure::from_logical_slots(ports.into_iter().rev()),
        )
        .unwrap();
        assert_eq!(value.rank(), 3);
        let result = value
            .with_lorentz_dimension(Dimension::Concrete(6))
            .unwrap();
        assert_eq!(result.structure.logical_slots(), vec![remaining]);
        let ports = [
            minkowski4.slot(index),
            minkowski6.slot(index),
            minkowski6.slot(index),
        ];
        let expression = FunctionBuilder::new(head)
            .add_args(ports.map(composition::port_atom))
            .finish();
        let value = SymbolicTensor::<PartialStructure>::checked_parts(
            expression,
            PartialStructure::from_logical_slots(ports),
        )
        .unwrap();
        assert_eq!(value.rank(), 1);
        assert!(
            value
                .with_lorentz_dimension(Dimension::Concrete(6))
                .is_err()
        );
    }

    #[test]
    fn lorentz_dimension_checks_actual_callback_output_without_replaying_callbacks() {
        use std::sync::{Arc, Mutex};

        let calls = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&calls);
        let head = spenso::tensor_symbol!(
            "dimension_result_normalizer",
            norm = move |node, output| {
                observed.lock().unwrap().push(node.to_owned());
                if let AtomView::Fun(function) = node
                    && let Some(AtomView::Fun(slot)) = function.iter().last()
                    && slot.iter().next() == Some(Atom::num(6).as_view())
                {
                    let mode = function.iter().next().unwrap();
                    if mode == Atom::Zero.as_view() || mode == Atom::one().as_view() {
                        **output = mode.to_owned();
                    }
                }
            }
        );
        let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        for index in [
            PartialIndex::open(0),
            PartialIndex::Explicit(AbstractIndex::Normal(73321)),
        ] {
            for mode in [0, 1, 2] {
                let port = representation.slot(index);
                let expression = FunctionBuilder::new(head)
                    .add_arg(mode)
                    .add_arg(composition::port_atom(port))
                    .finish();
                let value = SymbolicTensor::<PartialStructure>::new(
                    expression,
                    PartialStructure::from_logical_slots([port]),
                );
                calls.lock().unwrap().clear();
                let expected = value
                    .expression
                    .with_lorentz_dimension(Atom::num(6).as_view());
                let expected_calls = std::mem::take(&mut *calls.lock().unwrap());
                assert!(!expected_calls.is_empty());
                let result = value.with_lorentz_dimension(Dimension::Concrete(6));
                assert_eq!(*calls.lock().unwrap(), expected_calls);
                if mode == 1 {
                    assert_eq!(expected, Atom::one());
                    assert!(
                        result.is_err(),
                        "a scalar callback result cannot retain a tensor port"
                    );
                } else {
                    let result = result.unwrap();
                    assert_eq!(result.expression, expected);
                    assert_eq!(
                        result.structure.logical_slots(),
                        vec![
                            ExtendibleReps::MINKOWSKI
                                .new_rep(Dimension::Concrete(6))
                                .slot(index)
                        ]
                    );
                }
            }
        }

        let tensor = FunctionBuilder::new(spenso::tensor_symbol!("dimension_promoted_tensor"))
            .add_arg(ExtendibleReps::MINKOWSKI.new_rep(6).to_symbolic([]))
            .finish();
        let promoted_tensor = tensor.clone();
        let observed = Arc::clone(&calls);
        let scalar = symbolica::symbol!(
            "dimension_scalar_callback"; Scalar;
            norm = move |node, output| {
                observed.lock().unwrap().push(node.to_owned());
                if let AtomView::Fun(function) = node
                    && let Some(AtomView::Fun(slot)) = function.iter().next()
                    && slot.iter().next() == Some(Atom::num(6).as_view())
                {
                    **output = promoted_tensor.clone();
                }
            }
        );
        let expression = FunctionBuilder::new(scalar)
            .add_arg(representation.to_symbolic([]))
            .finish();
        assert!(!InterfaceInference::has_structured_syntax(
            expression.as_view()
        ));
        let value = SymbolicTensor::<PartialStructure>::new(
            expression,
            PartialStructure::from_logical_slots([]),
        );
        calls.lock().unwrap().clear();
        let expected = value
            .expression
            .with_lorentz_dimension(Atom::num(6).as_view());
        assert_eq!(expected, tensor);
        let expected_calls = std::mem::take(&mut *calls.lock().unwrap());
        assert!(
            value
                .with_lorentz_dimension(Dimension::Concrete(6))
                .is_err()
        );
        assert_eq!(*calls.lock().unwrap(), expected_calls);
    }

    #[test]
    fn cooking_preserves_logical_order_without_reinferring_and_rejects_callback_rank_loss() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let named = symbolica::symbol!("shared_cooking_index", tags = [&SPENSO_TAG.index]);
        let slots = [73101, 73103].map(|owner| {
            representation.slot(PartialIndex::Explicit(AbstractIndex::Named(
                named.into(),
                owner,
                2,
            )))
        });
        let expression = FunctionBuilder::new(spenso::tensor_symbol!("shared_cooking_tensor"))
            .add_args(slots.map(composition::port_atom))
            .finish();
        let original = SymbolicTensor::<PartialStructure>::new(
            expression.clone(),
            PartialStructure::from_logical_slots(slots.into_iter().rev()),
        );
        let settings = crate::CookSettings::indices();
        let expected = slots
            .into_iter()
            .rev()
            .map(|slot| {
                settings
                    .try_cook_indices(composition::port_atom(slot).as_view())
                    .unwrap()
            })
            .collect::<Vec<_>>();
        INFERENCE_CALLS.with(|count| count.set(0));
        let cooked = original.with_cooked_indices(&settings).unwrap();
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_ne!(cooked.expression, expression);
        assert_eq!(
            cooked
                .structure
                .logical_slots()
                .into_iter()
                .map(composition::port_atom)
                .collect::<Vec<_>>(),
            expected
        );

        let target = expected[0].clone();
        let callback = spenso::tensor_symbol!(
            "cooking_callback_strips_ports",
            norm = move |node, output| {
                if let AtomView::Fun(function) = node
                    && function.iter().any(|argument| argument == target.as_view())
                {
                    **output = Atom::one();
                }
            }
        );
        let expression = FunctionBuilder::new(callback)
            .add_arg(composition::port_atom(slots[1]))
            .finish();
        let value = SymbolicTensor::<PartialStructure>::new(
            expression,
            PartialStructure::from_logical_slots([slots[1]]),
        );
        assert!(value.with_cooked_indices(&settings).is_err());
    }

    #[test]
    fn checked_addition_preserves_left_order_and_typed_zero_without_reinference() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let slots = [73111, 73113]
            .map(|index| representation.slot(PartialIndex::Explicit(AbstractIndex::Normal(index))));
        let expression = FunctionBuilder::new(spenso::tensor_symbol!("shared_additive_tensor"))
            .add_args(slots.map(composition::port_atom))
            .finish();
        let left = SymbolicTensor::<PartialStructure>::new(
            expression.clone(),
            PartialStructure::from_logical_slots(slots.into_iter().rev()),
        );
        let right = SymbolicTensor::<PartialStructure>::new(
            expression,
            PartialStructure::from_logical_slots(slots),
        );
        INFERENCE_CALLS.with(|count| count.set(0));
        let sum = left.try_add(&right).unwrap();
        let zero = left.try_sub(&right).unwrap();
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(
            sum.structure.logical_slots(),
            left.structure.logical_slots()
        );
        assert_eq!(
            zero.structure.logical_slots(),
            left.structure.logical_slots()
        );
        assert!(zero.expression.is_zero());
        let scalar = SymbolicTensor::<PartialStructure>::new(
            Atom::one(),
            PartialStructure::from_logical_slots([]),
        );
        assert!(left.try_add(&scalar).is_err());
        assert!(left.try_sub(&scalar).is_err());
    }

    #[test]
    fn metric_result_validation_rejects_callback_rank_loss_and_retains_typed_zero() {
        let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let a = representation
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(72001))
            .to_atom();
        let b = representation
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(72003))
            .to_atom();
        let target = b.clone();
        let tensor = spenso::tensor_symbol!(
            "metric_inference_callback_rank_loss",
            norm = move |node, output| {
                if let AtomView::Fun(function) = node
                    && function.iter().any(|argument| argument == target.as_view())
                {
                    **output = Atom::one();
                }
            }
        );
        let atom = FunctionBuilder::new(ETS.metric)
            .add_arg(&a)
            .add_arg(&b)
            .finish()
            * FunctionBuilder::new(tensor).add_arg(a).finish();
        let value = SymbolicTensor::<PartialStructure>::infer(atom.clone()).unwrap();
        assert_eq!(value.rank(), 1);
        SymbolicTensor::<PartialStructure>::validate_interface(&atom, &value.structure).unwrap();

        let contracted = atom.simplify_metrics();
        assert_eq!(contracted, Atom::one());
        assert!(
            SymbolicTensor::<PartialStructure>::validate_interface(&contracted, &value.structure)
                .is_err()
        );
        assert!(value.with_rewritten_expression(contracted).is_err());
        SymbolicTensor::<PartialStructure>::validate_interface(&Atom::Zero, &value.structure)
            .unwrap();
    }

    #[test]
    fn terminal_trace_rewrites_validate_the_short_source_and_keep_logical_order() {
        use crate::dirac::GammaSimplifier;

        crate::test_support::test_initialize();
        let spin = Bispinor {}.new_rep(4).to_symbolic([]);
        let spectator = symbolica::parse_lit!((trace_interface_x + trace_interface_y) ^ 3);
        for dimension in [
            Dimension::Concrete(4),
            Dimension::from(symbolica::symbol!("trace_interface_D")),
        ] {
            let minkowski = Minkowski {}.new_rep(dimension);
            let indices = (0..8)
                .map(|position| {
                    minkowski
                        .slot::<AbstractIndex, _>(if position < 4 {
                            AbstractIndex::Normal(73501 + position)
                        } else {
                            AbstractIndex::Dummy(73501 + position)
                        })
                        .to_atom()
                })
                .collect::<Vec<_>>();
            let vector = spenso::vector_symbol!("trace_interface_vector");
            let compact = FunctionBuilder::new(vector)
                .add_arg(minkowski.to_symbolic([]))
                .finish();
            let distinct_vectors = (0..8)
                .map(|position| {
                    FunctionBuilder::new(vector)
                        .add_arg(position)
                        .add_arg(minkowski.to_symbolic([]))
                        .finish()
                })
                .collect::<Vec<_>>();
            let words = [
                indices.clone(),
                distinct_vectors,
                vec![
                    indices[0].clone(),
                    compact.clone(),
                    indices[0].clone(),
                    compact.clone(),
                ],
                vec![compact.clone(), compact.clone(), compact.clone(), compact],
            ];
            for word in words {
                let ordinary =
                    shadowing::trace(&spin, word.iter().map(|index| crate::gamma!(index)));
                let mut expressions = vec![ordinary.clone(), &spectator * ordinary];
                if dimension == Dimension::Concrete(4) && word == indices {
                    let axial = shadowing::trace(
                        &spin,
                        std::iter::once(crate::gamma5!())
                            .chain(word.iter().map(|index| crate::gamma!(index))),
                    );
                    expressions.extend([axial.clone(), &spectator * axial]);
                }
                for expression in expressions {
                    let mut value = SymbolicTensor::<PartialStructure>::infer(expression).unwrap();
                    value.structure = PartialStructure::from_logical_slots(
                        value.structure.logical_slots().into_iter().rev(),
                    );
                    assert!(
                        !InterfaceInference::default()
                            .algebra_preserves_leaf_interfaces(value.expression.as_view())
                    );
                    assert!(
                        InterfaceInference::default().terminal_trace_preserves_interface(
                            &value.expression,
                            &value.structure
                        )
                    );
                    let expected = value.expression.simplify_gamma();
                    assert_ne!(expected, value.expression);
                    if word.len() == 8 {
                        assert!(
                            expected.as_view().get_byte_size()
                                > value.expression.as_view().get_byte_size().saturating_mul(2)
                        );
                    }
                    INFERENCE_CALLS.with(|count| count.set(0));
                    SymbolicTensor::<PartialStructure>::validate_interface(
                        &value.expression,
                        &value.structure,
                    )
                    .unwrap();
                    let source_inference_calls = INFERENCE_CALLS.with(|count| count.get());
                    INFERENCE_CALLS.with(|count| count.set(0));
                    let result = value.with_rewritten_expression(expected.clone()).unwrap();
                    assert_eq!(result.expression, expected);
                    assert_eq!(result.structure, value.structure);
                    if word.len() == 8 {
                        assert_eq!(
                            INFERENCE_CALLS.with(|count| count.get()),
                            source_inference_calls,
                            "only the admitted source word needs interface inference"
                        );
                    }
                    SymbolicTensor::<PartialStructure>::validate_interface(
                        &result.expression,
                        &result.structure,
                    )
                    .unwrap();
                }
            }
        }
    }

    #[test]
    fn terminal_trace_interface_proof_rejects_unresolved_malformed_and_powered_words() {
        crate::test_support::test_initialize();
        let minkowski = Minkowski {}.new_rep(4);
        let spin = Bispinor {}.new_rep(4).to_symbolic([]);
        let [a, b] = [73521, 73523].map(|index| {
            minkowski
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let ordinary = shadowing::trace(&spin, [&a, &b].map(|index| crate::gamma!(index)));
        let structure = SymbolicTensor::<PartialStructure>::infer(ordinary.clone())
            .unwrap()
            .structure;
        let malformed_spin = shadowing::trace(
            Bispinor {}.new_rep(3).to_symbolic([]),
            [&a, &b].map(|index| crate::gamma!(index)),
        );
        let mixed_dimension = shadowing::trace(
            &spin,
            [
                crate::gamma!(&a),
                crate::gamma!(
                    Minkowski {}
                        .new_rep(6)
                        .slot::<AbstractIndex, _>(AbstractIndex::Normal(73523))
                        .to_atom()
                ),
            ],
        );
        let powered_factor = shadowing::trace(&spin, [crate::gamma!(&a).pow(2), crate::gamma!(&b)]);
        let chain = spenso::chain!(
            Bispinor {}.new_rep(4).slot::<AbstractIndex, _>(AbstractIndex::Normal(73525)).to_atom(),
            Bispinor {}.new_rep(4).slot::<AbstractIndex, _>(AbstractIndex::Normal(73527)).to_atom();
            [crate::gamma!(&a)]);
        let unresolved = shadowing::trace(
            &spin,
            [crate::gamma!(minkowski.to_symbolic([])), crate::gamma!(&b)],
        );
        let open_index = shadowing::trace(
            &spin,
            [
                crate::gamma!(
                    minkowski
                        .slot::<AbstractIndex, _>(AbstractIndex::Open {
                            owner: 73529,
                            axis: 0
                        })
                        .to_atom()
                ),
                crate::gamma!(&b),
            ],
        );
        for expression in [
            malformed_spin,
            mixed_dimension,
            powered_factor,
            chain,
            unresolved,
            open_index,
            ordinary.pow(2),
        ] {
            assert!(
                !InterfaceInference::default()
                    .terminal_trace_preserves_interface(&expression, &structure),
                "{expression}"
            );
        }
        let nested = shadowing::trace(&spin, std::iter::empty::<Atom>());
        let vector =
            FunctionBuilder::new(spenso::vector_symbol!("trace_interface_nested_metadata"))
                .add_arg(nested)
                .add_arg(minkowski.to_symbolic([]))
                .finish();
        let nested_metadata =
            shadowing::trace(&spin, [&vector, &b].map(|index| crate::gamma!(index)));
        assert!(
            DiracSimplifier::terminal_trace_interface_inputs(nested_metadata.as_view()).is_none()
        );
        assert!(
            !InterfaceInference::default()
                .terminal_trace_preserves_interface(&nested_metadata, &structure)
        );
        let cyclic_metadata =
            FunctionBuilder::new(spenso::vector_symbol!("trace_interface_cyclic_metadata"))
                .add_arg(shadowing::cyclic([Atom::var(symbolica::symbol!(
                    "trace_interface_metadata"
                ))]))
                .add_arg(minkowski.to_symbolic([]))
                .finish();
        let cyclic_metadata = shadowing::trace(
            &spin,
            [&cyclic_metadata, &b].map(|index| crate::gamma!(index)),
        );
        assert!(
            !InterfaceInference::default()
                .terminal_trace_preserves_interface(&cyclic_metadata, &structure)
        );
        let symmetric = FunctionBuilder::new(
            symbolica::symbol!("trace_interface_symmetric_vector"; Symmetric;
                tags = [&SPENSO_TAG.tensor, &SPENSO_TAG.rank1]),
        )
        .add_arg(7)
        .add_arg(minkowski.to_symbolic([]))
        .finish();
        let symmetric = shadowing::trace(&spin, [&symmetric, &b].map(|index| crate::gamma!(index)));
        assert!(
            !InterfaceInference::default()
                .terminal_trace_preserves_interface(&symmetric, &structure)
        );
        let wrong_structure = PartialStructure::from_logical_slots([]);
        assert!(
            !InterfaceInference::default().terminal_trace_preserves_interface(
                &shadowing::trace(&spin, [&a, &b].map(|index| crate::gamma!(index))),
                &wrong_structure
            )
        );
    }

    #[test]
    fn terminal_trace_callback_rank_loss_keeps_output_validation_without_replay() {
        use crate::dirac::GammaSimplifier;
        use std::sync::{Arc, Mutex};

        crate::test_support::test_initialize();
        let scalar = Atom::add_many((0..64).map(|index| {
            Atom::var(symbolica::symbol!(format!(
                "trace_interface_callback_scalar_{index}"
            )))
        }));
        let callback_scalar = scalar.clone();
        let calls = Arc::new(Mutex::new(Vec::new()));
        let recorded = Arc::clone(&calls);
        let vector = spenso::vector_symbol!(
            "trace_interface_callback_vector",
            norm = move |node, out| {
                recorded.lock().unwrap().push(node.to_owned());
                if let AtomView::Fun(function) = node
                    && let Some(AtomView::Fun(slot)) = function.iter().last()
                    && slot.get_nargs() == 2
                {
                    **out = callback_scalar.clone();
                }
            }
        );
        let minkowski = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(73531);
        let explicit = minkowski.slot::<AbstractIndex, _>(index).to_atom();
        let compact = FunctionBuilder::new(vector)
            .add_arg(minkowski.to_symbolic([]))
            .finish();
        let source = shadowing::trace(
            Bispinor {}.new_rep(4).to_symbolic([]),
            [&explicit, &compact].map(|index| crate::gamma!(index)),
        );
        let value =
            SymbolicTensor::<PartialStructure>::new(source, explicit_interface(minkowski, index));
        calls.lock().unwrap().clear();
        let rewritten = value.expression.simplify_gamma();
        assert_eq!(rewritten, Atom::num(4) * scalar);
        assert!(
            rewritten.as_view().get_byte_size()
                > value.expression.as_view().get_byte_size().saturating_mul(2)
        );
        let before_validation = calls.lock().unwrap().clone();
        assert!(!before_validation.is_empty());
        assert!(value.with_rewritten_expression(rewritten).is_err());
        assert_eq!(*calls.lock().unwrap(), before_validation);
    }

    #[test]
    fn proven_algebra_metric_and_zero_results_retain_interfaces_without_reinference() {
        let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let [a, b] = [72011, 72013].map(|index| {
            representation
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let tensor = spenso::tensor_symbol!("shared_reuse_tensor");
        let vector = |slot: &Atom| FunctionBuilder::new(tensor).add_arg(slot).finish();
        let metric = FunctionBuilder::new(ETS.metric)
            .add_arg(&a)
            .add_arg(&b)
            .finish();
        let original = SymbolicTensor::<PartialStructure>::infer(metric * vector(&a)).unwrap();
        let rewritten = original.expression.simplify_metrics();
        INFERENCE_CALLS.with(|count| count.set(0));
        let result = original
            .with_rewritten_expression(rewritten.clone())
            .unwrap();
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(result.expression, vector(&b));
        assert_eq!(
            result.structure.logical_slots(),
            original.structure.logical_slots()
        );

        let x = Atom::var(symbolica::symbol!("shared_reuse_x"));
        let y = Atom::var(symbolica::symbol!("shared_reuse_y"));
        let algebra = SymbolicTensor::<PartialStructure>::infer((x + y) * vector(&a)).unwrap();
        let expanded = algebra.expression.expand();
        assert_ne!(expanded, algebra.expression);
        INFERENCE_CALLS.with(|count| count.set(0));
        let expanded = algebra.with_algebra_result(expanded).unwrap();
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(
            expanded.structure.logical_slots(),
            algebra.structure.logical_slots()
        );

        let unresolved = SymbolicTensor::<PartialStructure>::new(
            FunctionBuilder::new(tensor)
                .add_arg(representation.to_symbolic([]))
                .finish(),
            PartialStructure::from_logical_slots([representation.slot(PartialIndex::open(0))]),
        );
        INFERENCE_CALLS.with(|count| count.set(0));
        for result in [
            unresolved.with_algebra_result(Atom::Zero).unwrap(),
            unresolved.with_rewritten_expression(Atom::Zero).unwrap(),
            unresolved
                .with_rewritten_expression(unresolved.expression.clone())
                .unwrap(),
        ] {
            assert_eq!(
                result.structure.logical_slots(),
                unresolved.structure.logical_slots()
            );
        }
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    }

    #[test]
    fn unchanged_results_preserve_layout_dispatch_flags_and_callback_state() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };

        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let head = spenso::tensor_symbol!(
            "unchanged_result_callback",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let ports = [
            ExtendibleReps::MINKOWSKI
                .new_rep(Dimension::Concrete(4))
                .slot(PartialIndex::Explicit(AbstractIndex::Normal(72017))),
            ExtendibleReps::EUCLIDEAN
                .new_rep(Dimension::Concrete(3))
                .slot(PartialIndex::open(0)),
        ];
        let expression = FunctionBuilder::new(head)
            .add_arg(17)
            .add_args(ports.map(composition::port_atom))
            .finish();
        for expression in [expression, Atom::Zero] {
            for (is_metric, is_composite) in
                [(false, false), (false, true), (true, false), (true, true)]
            {
                let value = SymbolicTensor {
                    expression: expression.clone(),
                    structure: PartialStructure::from_logical_slots(ports),
                    is_metric,
                    is_composite,
                };
                for finish in [
                    SymbolicTensor::<PartialStructure>::with_algebra_result,
                    SymbolicTensor::<PartialStructure>::with_rewritten_expression,
                    SymbolicTensor::<PartialStructure>::with_transformed_expression,
                ] {
                    calls.store(0, Ordering::Relaxed);
                    INFERENCE_CALLS.with(|count| count.set(0));
                    let result = finish(&value, value.expression.clone()).unwrap();
                    assert_eq!(result, value);
                    assert_eq!(calls.load(Ordering::Relaxed), 0);
                    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
                }
            }
        }
    }

    #[test]
    fn reused_reusable_interfaces_preserve_interface_and_multiplicity_validation() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let [a, b, c] = [71, 73, 79].map(|index| {
            representation
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let metric = |left: &Atom, right: &Atom| {
            FunctionBuilder::new(ETS.metric)
                .add_arg(left)
                .add_arg(right)
                .finish()
        };
        let ab = metric(&a, &b);
        let ac = metric(&a, &c);
        let mut inference = InterfaceInference::default();
        let expected = inference.infer(&ab).unwrap();
        assert_eq!(
            inference.infer(&ab).unwrap().logical_slots(),
            expected.logical_slots()
        );
        assert!(
            inference
                .reusable_interfaces
                .contains_key(ab.as_view().get_data())
        );

        let contracted = inference.infer(&(ab.as_ref() * ac.as_ref())).unwrap();
        assert_eq!(contracted.canonical().order(), 2);
        for invalid in [
            ab.as_ref() + ac.as_ref(),
            ab.as_ref() + Atom::num(1),
            ab.as_ref().pow(Atom::num(3)),
        ] {
            assert!(inference.infer(&invalid).is_err());
        }

        // A leaf's own repeated slots must remain visible to enclosing products.
        let diagonal = FunctionBuilder::new(ETS.flat)
            .add_arg(&a)
            .add_arg(&a)
            .finish();
        assert_eq!(inference.infer(&diagonal).unwrap().canonical().order(), 2);
        let vector = FunctionBuilder::new(spenso::vector_symbol!("inference_reuse_vector"))
            .add_arg(&a)
            .finish();
        assert!(inference.infer(&(diagonal * vector)).is_err());

        let flat = FunctionBuilder::new(ETS.flat)
            .add_arg(&a)
            .add_arg(&b)
            .finish();
        let sum = ab.as_ref() * ac.as_ref() + flat * ac.as_ref();
        let actual = InterfaceInference::default()
            .infer_validated(sum.as_view())
            .unwrap();
        assert!(InterfaceInference::additive_interfaces_match(
            &contracted,
            &actual
        ));
    }

    #[test]
    fn interface_reuse_excludes_unresolved_builtin_ports() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let compact = representation.to_symbolic([]);
        let implicit = FunctionBuilder::new(ETS.metric)
            .add_arg(&compact)
            .add_arg(&compact)
            .finish();
        let materialized = FunctionBuilder::new(ETS.metric)
            .add_args((0..2).map(|axis| {
                representation
                    .slot::<AbstractIndex, _>(AbstractIndex::Open { owner: 81, axis })
                    .to_atom()
            }))
            .finish();
        let mut inference = InterfaceInference::default();
        for (atom, rank) in [(implicit, 2), (materialized, 2)] {
            for _ in 0..2 {
                let interface = inference.infer(&atom).unwrap();
                assert_eq!(interface.open_positions(), (0..rank).collect::<Vec<_>>());
                assert!(inference.reusable_interfaces.is_empty());
            }
        }
    }

    #[test]
    fn direct_leaf_inference_preserves_open_occurrences_and_logical_order() {
        let euclidean = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let minkowski = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let leaf = FunctionBuilder::new(spenso::tensor_symbol!("direct_open_leaf"))
            .add_arg(7)
            .add_arg(minkowski.to_symbolic([]))
            .add_arg(euclidean.to_symbolic([]))
            .finish();
        let mut inference = InterfaceInference::default();
        let interface = inference.infer_validated(leaf.as_view()).unwrap();
        assert_eq!(
            interface
                .logical_slots()
                .iter()
                .map(IsAbstractSlot::rep)
                .collect::<Vec<_>>(),
            vec![minkowski, euclidean],
        );
        assert!(
            inference
                .reusable_interfaces
                .contains_key(leaf.as_view().get_data())
        );
        let repeated = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(&leaf)
            .add_arg(&leaf)
            .finish();
        let interface = inference.infer_validated(repeated.as_view()).unwrap();
        assert_eq!(interface.open_positions(), vec![0, 1, 2, 3]);
        assert_eq!(
            interface
                .logical_slots()
                .into_iter()
                .map(|slot| slot.aind)
                .collect::<Vec<_>>(),
            (0..4).map(PartialIndex::open).collect::<Vec<_>>(),
        );
        let incompatible = FunctionBuilder::new(spenso::tensor_symbol!("direct_other_leaf"))
            .add_arg(euclidean.to_symbolic([]))
            .add_arg(minkowski.to_symbolic([]))
            .finish();
        assert!(
            inference
                .infer_validated((&leaf + incompatible).as_view())
                .is_err()
        );

        let p = FunctionBuilder::new(spenso::vector_symbol!("direct_compact_p"))
            .add_arg(minkowski.to_symbolic([]))
            .finish();
        let q = FunctionBuilder::new(spenso::vector_symbol!("direct_compact_q"))
            .add_arg(minkowski.to_symbolic([]))
            .finish();
        assert!(
            inference
                .infer_validated(ETS.metric(p, q).as_view())
                .unwrap()
                .canonical()
                .is_scalar()
        );
    }

    #[test]
    fn direct_leaf_inference_retains_materialization_callback_semantics() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let callback = spenso::tensor_symbol!(
            "materialized_leaf_callback",
            norm = |node, output| {
                if let AtomView::Fun(function) = node
                    && function.iter().any(|argument| {
                        Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_ok()
                    })
                {
                    **output = Atom::Zero;
                }
            }
        );
        let leaf = FunctionBuilder::new(callback)
            .add_arg(representation.to_symbolic([]))
            .finish();
        assert!(matches!(leaf.as_view(), AtomView::Fun(_)));
        let mut inference = InterfaceInference::default();
        assert!(
            inference
                .infer_validated(leaf.as_view())
                .unwrap()
                .canonical()
                .is_scalar()
        );
        assert!(inference.reusable_interfaces.is_empty());
    }

    #[test]
    fn repeated_scalar_dot_inference_preserves_external_ports_and_errors() {
        let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let p = FunctionBuilder::new(spenso::vector_symbol!("cached_dot_p"))
            .add_arg(representation.to_symbolic([]))
            .finish();
        let q = FunctionBuilder::new(spenso::vector_symbol!("cached_dot_q"))
            .add_arg(representation.to_symbolic([]))
            .finish();
        let incompatible = FunctionBuilder::new(spenso::vector_symbol!("cached_dot_q"))
            .add_arg(ExtendibleReps::EUCLIDEAN.new_rep(4).to_symbolic([]))
            .finish();
        let mut inference = InterfaceInference::default();
        for head in [SPENSO_TAG.dot, ETS.metric] {
            let dot = FunctionBuilder::new(head).add_arg(&p).add_arg(&q).finish();
            for _ in 0..2 {
                assert!(
                    inference
                        .infer_validated(dot.as_view())
                        .unwrap()
                        .canonical()
                        .is_scalar()
                );
            }
            assert!(
                inference
                    .reusable_interfaces
                    .contains_key(dot.as_view().get_data())
            );
            let external = inference.infer_validated((&dot * &p).as_view()).unwrap();
            assert_eq!(
                external.logical_slots(),
                vec![representation.slot(PartialIndex::open(0))]
            );

            let invalid = FunctionBuilder::new(head)
                .add_arg(&p)
                .add_arg(&incompatible)
                .finish();
            assert!(inference.infer_validated(invalid.as_view()).is_err());
            assert!(
                !inference
                    .reusable_interfaces
                    .contains_key(invalid.as_view().get_data())
            );
        }
    }

    #[test]
    fn scalar_dot_inference_retains_operand_materialization_callbacks() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let callback = spenso::tensor_symbol!(
            "cached_dot_materialization_callback",
            norm = move |node, _| {
                if let AtomView::Fun(function) = node
                    && function.iter().any(|argument| {
                        Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_ok()
                    })
                {
                    observed.fetch_add(1, Ordering::Relaxed);
                }
            }
        );
        let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let p = FunctionBuilder::new(callback)
            .add_arg(representation.to_symbolic([]))
            .finish();
        let q = FunctionBuilder::new(spenso::vector_symbol!("cached_dot_plain_operand"))
            .add_arg(representation.to_symbolic([]))
            .finish();
        let dot = FunctionBuilder::new(SPENSO_TAG.dot)
            .add_arg(p)
            .add_arg(q)
            .finish();
        let mut inference = InterfaceInference::default();
        calls.store(0, Ordering::Relaxed);
        for _ in 0..2 {
            let previous = calls.load(Ordering::Relaxed);
            assert!(
                inference
                    .infer_validated(dot.as_view())
                    .unwrap()
                    .canonical()
                    .is_scalar()
            );
            assert!(calls.load(Ordering::Relaxed) > previous);
            assert!(
                !inference
                    .reusable_interfaces
                    .contains_key(dot.as_view().get_data())
            );
        }
        let previous = calls.load(Ordering::Relaxed);
        SymbolicTensor::<PartialStructure>::validate_interface(
            &dot,
            &PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        assert_eq!(calls.load(Ordering::Relaxed), previous);
    }

    #[test]
    fn result_observation_keeps_encoded_ports_without_replaying_constructor_callbacks() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let callback = spenso::tensor_symbol!(
            "observed_result_materialization_changes_rank",
            norm = move |node, output| {
                observed.fetch_add(1, Ordering::Relaxed);
                if let AtomView::Fun(function) = node
                    && function.iter().last().is_some_and(|argument| {
                        Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_ok()
                    })
                {
                    **output = Atom::Zero;
                }
            }
        );
        let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let metadata = FunctionBuilder::new(symbolica::symbol!(
            "observed_result_scalar_metadata"; Scalar
        ))
        .add_arg(representation.to_symbolic([]))
        .finish();
        let structure =
            PartialStructure::from_logical_slots([representation.slot(PartialIndex::open(0))]);
        for arguments in [vec![], vec![metadata]] {
            let atom = FunctionBuilder::new(callback)
                .add_args(arguments)
                .add_arg(representation.to_symbolic([]))
                .finish();
            calls.store(0, Ordering::Relaxed);
            let inferred = InterfaceInference::default()
                .infer_validated(atom.as_view())
                .unwrap();
            assert!(inferred.canonical().is_scalar());
            assert!(calls.load(Ordering::Relaxed) > 0);

            calls.store(0, Ordering::Relaxed);
            SymbolicTensor::<PartialStructure>::validate_interface(&atom, &structure).unwrap();
            assert_eq!(calls.load(Ordering::Relaxed), 0);

            let x = Atom::var(symbolica::symbol!("observed_result_x"));
            let y = Atom::var(symbolica::symbol!("observed_result_y"));
            let value =
                SymbolicTensor::<PartialStructure>::new((&x + &y) * &atom, structure.clone());
            let expanded = value.expression.expand();
            calls.store(0, Ordering::Relaxed);
            let result = value.with_algebra_result(expanded.clone()).unwrap();
            assert_eq!(result.expression, expanded);
            assert_eq!(result.structure, structure);
            assert_eq!(calls.load(Ordering::Relaxed), 0);
        }
        // Actual rank loss is still rejected; only zero may retain a typed interface.
        assert!(
            SymbolicTensor::<PartialStructure>::validate_interface(&Atom::one(), &structure)
                .is_err()
        );
        SymbolicTensor::<PartialStructure>::validate_interface(&Atom::Zero, &structure).unwrap();
    }

    #[test]
    fn tensor_power_lowering_preserves_normalizer_calls_on_unchanged_functions() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let head = symbolica::symbol!(
            "power_lowering_callback",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let atom = FunctionBuilder::new(head).add_arg(3).finish();
        calls.store(0, Ordering::Relaxed);
        assert!(
            InterfaceInference::lower_tensor_powers(atom.as_view())
                .unwrap()
                .is_none()
        );
        assert!(calls.load(Ordering::Relaxed) > 0);
    }

    #[test]
    fn recursive_merging_rejects_three_compatible_explicit_ports() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(11);
        let interfaces = vec![
            explicit_interface(representation, index),
            explicit_interface(representation, index),
            explicit_interface(representation, index),
        ];

        assert!(InterfaceInference::merge_explicit_interface_sequence(&interfaces).is_err());
    }

    #[test]
    fn recursive_merging_contracts_one_explicit_pair() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(13);
        let interfaces = vec![
            explicit_interface(representation, index),
            explicit_interface(representation, index),
        ];

        assert!(
            InterfaceInference::merge_explicit_interface_sequence(&interfaces)
                .unwrap()
                .canonical()
                .is_scalar()
        );
    }

    #[test]
    fn recursive_merging_normalizes_explicit_indices_within_one_interface() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(15);
        let pair = PartialStructure::from_logical_slots([
            representation.slot(PartialIndex::Explicit(index)),
            representation.slot(PartialIndex::Explicit(index)),
        ]);
        let triple = PartialStructure::from_logical_slots([
            representation.slot(PartialIndex::Explicit(index)),
            representation.slot(PartialIndex::Explicit(index)),
            representation.slot(PartialIndex::Explicit(index)),
        ]);

        assert!(
            InterfaceInference::merge_explicit_interface_sequence(&[pair])
                .unwrap()
                .canonical()
                .is_scalar()
        );
        assert!(InterfaceInference::merge_explicit_interface_sequence(&[triple]).is_err());
    }

    #[test]
    fn additive_reinference_compares_canonical_scalar_interfaces() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let diagonal = |name: &str, index| {
            let slot = representation.slot::<AbstractIndex, _>(index).to_atom();
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(&slot)
                .add_arg(&slot)
                .finish()
        };
        let left = diagonal("additive_scalar_left", AbstractIndex::Normal(21));
        let right = diagonal("additive_scalar_right", AbstractIndex::Normal(23));

        assert!(
            InterfaceInference::default()
                .infer(&(left.as_ref() + Atom::num(1)))
                .unwrap()
                .canonical()
                .is_scalar()
        );
        assert!(
            InterfaceInference::default()
                .infer(&(left + right))
                .unwrap()
                .canonical()
                .is_scalar()
        );
    }

    #[test]
    fn powers_reject_more_than_two_explicit_index_occurrences() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(17);
        let tensor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("powered_explicit_index"))
            .add_arg(representation.slot::<AbstractIndex, _>(index).to_atom())
            .finish();

        assert!(
            InterfaceInference::default()
                .infer(&tensor.as_ref().pow(Atom::num(2)))
                .unwrap()
                .canonical()
                .is_scalar()
        );
        assert!(
            InterfaceInference::default()
                .infer(&tensor.as_ref().pow(Atom::num(3)))
                .is_err()
        );
    }

    #[test]
    fn powers_count_repeated_explicit_indices_in_the_base() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(19);
        let slot = representation.slot::<AbstractIndex, _>(index).to_atom();
        let tensor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("powered_repeated_index"))
            .add_arg(&slot)
            .add_arg(&slot)
            .finish();

        assert!(
            InterfaceInference::default()
                .infer(&tensor.as_ref().pow(Atom::num(2)))
                .is_err()
        );
    }

    #[test]
    fn open_ports_retain_power_parity_semantics() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let tensor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("powered_open_index"))
            .add_arg(representation.to_symbolic([]))
            .finish();

        assert!(
            InterfaceInference::default()
                .infer(&tensor.as_ref().pow(Atom::num(2)))
                .unwrap()
                .canonical()
                .is_scalar()
        );
        assert_eq!(
            InterfaceInference::default()
                .infer(&tensor.as_ref().pow(Atom::num(3)))
                .unwrap()
                .canonical()
                .order(),
            1
        );
    }

    #[test]
    fn leaf_slots_follow_syntactic_argument_order() {
        let euclidean = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let minkowski = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let atom = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(
                euclidean
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(17))
                    .to_atom(),
            )
            .add_arg(
                minkowski
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(19))
                    .to_atom(),
            )
            .finish();

        assert_eq!(
            InterfaceInference::syntactic_leaf_slots(atom.as_view())
                .into_iter()
                .map(|slot| slot.rep())
                .collect::<Vec<_>>(),
            vec![euclidean, minkowski]
        );
    }
}
