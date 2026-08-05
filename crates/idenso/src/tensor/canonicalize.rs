//! Signed canonicalization for Spenso symbolic tensor networks.
//!
//! The network remains the semantic owner of products, sums, powers, and
//! enclosing functions. Each iteration projects the complete operation tree,
//! tensor ports, and index-line components into one graph. Graphica then
//! supplies both the canonical labeling and stabilizer generators. Hidden graph
//! data records only reconstruction origins; it does not participate in graph
//! identity.

mod driver;
mod group;
mod projection;
mod reconstruct;

use std::collections::BTreeSet;

use spenso::{
    network::{
        parsing::{AtomStructureExt, StrictTensorFilter},
        tags::SPENSO_TAG,
    },
    shadowing::{ANTISYM, CYCLIC, SYM},
    structure::{
        OrderedStructure, TensorStructure,
        abstract_index::AIND_SYMBOLS,
        representation::{LibraryRep, LibrarySlot},
        slot::{AbsInd, IsAbstractSlot, ParseableAind, Slot},
    },
};
use symbolica::{
    atom::{Atom, AtomView, FunctionBuilder, Symbol},
    domains::rational::Rational,
};
use thiserror::Error;

use crate::shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip};

use super::SymbolicTensor;

mod policy;
mod semantic;

pub(crate) use policy::CanonicalPolicyNet;
pub(crate) use semantic::semantic_atom_digest;

/// Failures detected while validating or canonicalizing a symbolic network.
#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum CanonicalizationError {
    #[error(
        "intrinsically symmetric tensor {head} requires a direct slot at argument {argument}: {expression}"
    )]
    InvalidIntrinsicArgument {
        head: Symbol,
        argument: usize,
        expression: Atom,
    },
    #[error(
        "structural symmetry group {group} requires a direct slot at member {member}: {expression}"
    )]
    InvalidPartialGroup {
        group: Symbol,
        member: usize,
        expression: Atom,
    },
    #[error(
        "ordered tensor {head} cannot hide exposed tensor topology at argument {argument}: {expression}"
    )]
    HiddenTensorTopology {
        head: Symbol,
        argument: usize,
        expression: Atom,
    },
    #[error("tensor layout does not match its network structure: {expression}")]
    StructureMismatch { expression: Atom },
    #[error("unsupported symbolic-network leaf in canonicalization")]
    UnsupportedLeaf,
    /// Stabilizer zero analysis distinguishes a global zero from local group
    /// zeros keyed by their selected canonical occurrence. This error reports
    /// that the sign cannot be assigned to either semantic class exactly.
    #[error("canonical sign cannot be assigned to one nonlinear scope")]
    AmbiguousSignScope,
    #[error("failed to parse symbolic tensor network: {0}")]
    Network(String),
    #[error(
        "failed to parse symbolic tensor network: Negative non-even power on non-scalar node:{0}"
    )]
    NegativeExponentNonScalar(String),
    #[error("failed to parse symbolic tensor network: Non self-dual tensor power{0}")]
    NonSelfDualTensorPower(String),
    #[error("failed to project the canonical tensor network: {0}")]
    Projection(String),
    #[error("cannot fold materialized Power expression {expression}: {reason}")]
    PowerReconstruction { expression: usize, reason: String },
    #[error("failed to execute canonical scope {scope}: {error}")]
    Execution { scope: String, error: String },
    #[error(
        "expanded tensor power {power} (magnitude {magnitude}) graph estimation reached {requested_vertices} vertices and {requested_edges} edges, exceeding limits {vertex_limit}/{edge_limit}: {expression}"
    )]
    GraphSizeLimit {
        expression: Atom,
        power: i8,
        magnitude: u8,
        requested_vertices: usize,
        requested_edges: usize,
        vertex_limit: usize,
        edge_limit: usize,
    },
    #[error(
        "canonical graph estimation reached {requested_vertices} vertices and {requested_edges} edges, exceeding limits {vertex_limit}/{edge_limit}: {expression}"
    )]
    WholeGraphSizeLimit {
        expression: Atom,
        requested_vertices: usize,
        requested_edges: usize,
        vertex_limit: usize,
        edge_limit: usize,
    },
    #[error(
        "canonicalization entered a cycle: iteration {repeated_iteration} repeats iteration {first_iteration} (length {cycle_length}); retry reasons: {retry_reasons:?}"
    )]
    ConvergenceCycle {
        first_iteration: usize,
        repeated_iteration: usize,
        cycle_length: usize,
        retry_reasons: Vec<String>,
    },
    #[error("canonicalization did not stabilize within {limit} iterations: {last_reason}")]
    IterationLimit { limit: usize, last_reason: String },
    #[error(
        "tensorial power requires a signed 8-bit integer exponent: {expression} (base: {base}, exponent: {exponent})"
    )]
    UnsupportedTensorPowerExponent {
        expression: Atom,
        base: Atom,
        exponent: Atom,
    },
    #[error("failed to reversibly encode a structured index: {0}")]
    StructuredIndex(String),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum SymmetryKind {
    Symmetric,
    Antisymmetric,
    Cyclic,
}

impl SymmetryKind {
    fn of(symbol: Symbol) -> Option<Self> {
        if symbol.is_antisymmetric() {
            Some(Self::Antisymmetric)
        } else if symbol.is_symmetric() {
            Some(Self::Symmetric)
        } else if symbol.is_cyclesymmetric() {
            Some(Self::Cyclic)
        } else {
            None
        }
    }

    fn structural(symbol: Symbol) -> Option<Self> {
        if symbol == *SYM {
            Some(Self::Symmetric)
        } else if symbol == *ANTISYM {
            Some(Self::Antisymmetric)
        } else if symbol == *CYCLIC {
            Some(Self::Cyclic)
        } else {
            None
        }
    }
}

#[derive(Clone, Debug)]
enum LayoutArgument {
    Opaque(Atom),
    DirectSlot(usize),
    SlotBundle {
        holes: Vec<usize>,
    },
    Group {
        symbol: Symbol,
        kind: SymmetryKind,
        prefix: Atom,
        holes: Vec<usize>,
    },
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum LayoutColorArgument {
    Opaque(semantic::SemanticAtomKey),
    DirectSlot,
    SlotBundle(usize),
    Group {
        kind: SymmetryKind,
        members: usize,
        prefix: semantic::SemanticAtomKey,
    },
}

#[derive(Clone, Debug)]
struct TensorLayout {
    head: Symbol,
    arguments: Vec<LayoutArgument>,
    color: TensorColor,
    slot_count: usize,
    structural_holes: Vec<usize>,
    intrinsic: Option<SymmetryKind>,
    outer_linear: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct TensorColor {
    head: semantic::SemanticSymbolKey,
    arguments: Vec<LayoutColorArgument>,
}

type PartialGroup<Aind> = (Symbol, SymmetryKind, Atom, Vec<LibrarySlot<Aind>>);

impl TensorLayout {
    fn scan<Aind: AbsInd + ParseableAind>(
        tensor: &SymbolicTensor<Aind>,
    ) -> Result<Self, CanonicalizationError> {
        let AtomView::Fun(function) = tensor.expression.as_view() else {
            return Err(CanonicalizationError::StructureMismatch {
                expression: tensor.expression.clone(),
            });
        };
        let head = function.get_symbol();
        let tagged = head.has_tag(&SPENSO_TAG.tensor);
        let intrinsic = tagged.then(|| SymmetryKind::of(head)).flatten();
        let root_group = SymmetryKind::structural(head);

        let mut arguments = Vec::with_capacity(function.get_nargs());
        let mut colors = Vec::with_capacity(function.get_nargs());
        let mut expression_slots = Vec::new();
        let mut slot_count = 0;

        if intrinsic.is_some() || root_group.is_some() {
            for (argument, value) in function.iter().enumerate() {
                let Ok(slot) = Slot::<LibraryRep, Aind>::try_from(value) else {
                    return Err(CanonicalizationError::InvalidIntrinsicArgument {
                        head,
                        argument,
                        expression: tensor.expression.clone(),
                    });
                };
                arguments.push(LayoutArgument::DirectSlot(slot_count));
                colors.push(LayoutColorArgument::DirectSlot);
                expression_slots.push(slot);
                slot_count += 1;
            }
        } else {
            for (argument, value) in function.iter().enumerate() {
                if let Ok(slot) = Slot::<LibraryRep, Aind>::try_from(value) {
                    arguments.push(LayoutArgument::DirectSlot(slot_count));
                    colors.push(LayoutColorArgument::DirectSlot);
                    expression_slots.push(slot);
                    slot_count += 1;
                    continue;
                }

                if let AtomView::Fun(bundle) = value
                    && bundle.get_symbol() == AIND_SYMBOLS.aind
                {
                    let start = slot_count;
                    for member in bundle.iter() {
                        let Ok(slot) = Slot::<LibraryRep, Aind>::try_from(member) else {
                            return Err(CanonicalizationError::StructureMismatch {
                                expression: tensor.expression.clone(),
                            });
                        };
                        expression_slots.push(slot);
                        slot_count += 1;
                    }
                    let holes = (start..slot_count).collect::<Vec<_>>();
                    colors.push(LayoutColorArgument::SlotBundle(holes.len()));
                    arguments.push(LayoutArgument::SlotBundle { holes });
                    continue;
                }

                if tagged
                    && let Some((symbol, kind, prefix, members)) =
                        Self::partial_group::<Aind>(value, &tensor.expression)?
                {
                    let start = slot_count;
                    slot_count += members.len();
                    expression_slots.extend(members);
                    let holes = (start..slot_count).collect::<Vec<_>>();
                    colors.push(LayoutColorArgument::Group {
                        kind,
                        members: holes.len(),
                        prefix: semantic::SemanticAtomKey::new(prefix.as_view()),
                    });
                    arguments.push(LayoutArgument::Group {
                        symbol,
                        kind,
                        prefix,
                        holes,
                    });
                    continue;
                }

                if value.contains_exposed_tensor_topology(StrictTensorFilter::Tagged) {
                    return Err(CanonicalizationError::HiddenTensorTopology {
                        head,
                        argument,
                        expression: tensor.expression.clone(),
                    });
                }

                arguments.push(LayoutArgument::Opaque(value.to_owned()));
                colors.push(LayoutColorArgument::Opaque(semantic::SemanticAtomKey::new(
                    value,
                )));
            }
        }

        if slot_count != tensor.structure.order() {
            return Err(CanonicalizationError::StructureMismatch {
                expression: tensor.expression.clone(),
            });
        }
        let structural_slots = tensor
            .structure
            .external_structure_iter()
            .collect::<Vec<_>>();
        if expression_slots.len() != structural_slots.len() {
            return Err(CanonicalizationError::StructureMismatch {
                expression: tensor.expression.clone(),
            });
        }
        let inferred = OrderedStructure::new(expression_slots);
        if inferred.structure != tensor.structure {
            return Err(CanonicalizationError::StructureMismatch {
                expression: tensor.expression.clone(),
            });
        }
        let mut structural_holes = (0..slot_count).collect::<Vec<_>>();
        inferred
            .rep_permutation
            .apply_slice_in_place(&mut structural_holes);
        inferred
            .index_permutation
            .apply_slice_in_place(&mut structural_holes);

        Ok(Self {
            head,
            arguments,
            color: TensorColor {
                head: head.into(),
                arguments: colors,
            },
            slot_count,
            structural_holes,
            intrinsic: intrinsic.or(root_group),
            outer_linear: head.is_linear(),
        })
    }

    fn partial_group<Aind: AbsInd + ParseableAind>(
        value: AtomView<'_>,
        expression: &Atom,
    ) -> Result<Option<PartialGroup<Aind>>, CanonicalizationError> {
        let (prefix, group) = match value {
            AtomView::Fun(group)
                if group.get_nargs() > 0
                    && SymmetryKind::structural(group.get_symbol()).is_some() =>
            {
                (Atom::num(1), group)
            }
            AtomView::Mul(product) => {
                let mut factors = product.iter();
                let Some(first) = factors.next() else {
                    return Ok(None);
                };
                let Some(second) = factors.next() else {
                    return Ok(None);
                };
                if factors.next().is_some() {
                    return Ok(None);
                }

                let (coefficient, group) = match (first, second) {
                    (coefficient @ AtomView::Num(_), AtomView::Fun(group))
                    | (AtomView::Fun(group), coefficient @ AtomView::Num(_)) => {
                        (coefficient, group)
                    }
                    _ => return Ok(None),
                };
                if Rational::try_from(coefficient).ok() != Some(Rational::from(-1))
                    || group.get_symbol() != *ANTISYM
                    || group.get_nargs() == 0
                {
                    return Ok(None);
                }
                (Atom::num(-1), group)
            }
            _ => return Ok(None),
        };

        let symbol = group.get_symbol();
        let kind = SymmetryKind::structural(symbol).expect("structural group was classified");
        let mut slots = Vec::with_capacity(group.get_nargs());
        for (member, value) in group.iter().enumerate() {
            let Ok(slot) = Slot::<LibraryRep, Aind>::try_from(value) else {
                return Err(CanonicalizationError::InvalidPartialGroup {
                    group: symbol,
                    member,
                    expression: expression.clone(),
                });
            };
            slots.push(slot);
        }
        Ok(Some((symbol, kind, prefix, slots)))
    }

    fn incidence_role(&self, flat_slot: usize) -> IncidenceRole {
        if let Some(kind) = self.intrinsic {
            return IncidenceRole::Group {
                key: GroupKey::Intrinsic,
                kind,
            };
        }
        for (argument, layout) in self.arguments.iter().enumerate() {
            match layout {
                LayoutArgument::DirectSlot(hole) if *hole == flat_slot => {
                    return IncidenceRole::Ordered { argument };
                }
                LayoutArgument::SlotBundle { holes } => {
                    if let Some(member) = holes.iter().position(|hole| *hole == flat_slot) {
                        return IncidenceRole::OrderedBundle { argument, member };
                    }
                }
                LayoutArgument::Group { kind, holes, .. } => {
                    if holes.contains(&flat_slot) {
                        return IncidenceRole::Group {
                            key: GroupKey::Argument(argument),
                            kind: *kind,
                        };
                    }
                }
                LayoutArgument::Opaque(_) | LayoutArgument::DirectSlot(_) => {}
            }
        }

        unreachable!("every flat slot belongs to one reversible layout hole")
    }

    fn member_position(&self, flat_slot: usize) -> usize {
        if self.intrinsic.is_some() {
            return self
                .arguments
                .iter()
                .position(
                    |layout| matches!(layout, LayoutArgument::DirectSlot(hole) if *hole == flat_slot),
                )
                .expect("every intrinsic slot is one direct argument");
        }
        for layout in &self.arguments {
            match layout {
                LayoutArgument::Group { holes, .. } | LayoutArgument::SlotBundle { holes } => {
                    if let Some(position) = holes.iter().position(|hole| *hole == flat_slot) {
                        return position;
                    }
                }
                LayoutArgument::DirectSlot(hole) if *hole == flat_slot => return 0,
                _ => {}
            }
        }
        flat_slot
    }

    fn antisymmetric_groups(&self) -> impl Iterator<Item = (GroupKey, bool)> + '_ {
        let intrinsic = (self.intrinsic == Some(SymmetryKind::Antisymmetric))
            .then_some((GroupKey::Intrinsic, true));
        intrinsic
            .into_iter()
            .chain(self.arguments.iter().enumerate().filter_map(
                |(argument, layout)| match layout {
                    LayoutArgument::Group {
                        kind: SymmetryKind::Antisymmetric,
                        ..
                    } => Some((GroupKey::Argument(argument), self.outer_linear)),
                    _ => None,
                },
            ))
    }

    fn rebuild<Aind: AbsInd + ParseableAind>(
        &self,
        slots: &[LibrarySlot<Aind>],
        zero_groups: &BTreeSet<GroupKey>,
        negative_groups: &BTreeSet<GroupKey>,
    ) -> (Atom, Vec<LibrarySlot<Aind>>) {
        let mut builder = FunctionBuilder::new(self.head);
        let mut active_slots = Vec::new();
        for (argument, layout) in self.arguments.iter().enumerate() {
            let value = match layout {
                LayoutArgument::Opaque(value) => value.clone(),
                LayoutArgument::DirectSlot(hole) => {
                    active_slots.push(slots[*hole]);
                    slots[*hole].to_atom()
                }
                LayoutArgument::SlotBundle { holes } => {
                    let mut bundle = FunctionBuilder::new(AIND_SYMBOLS.aind);
                    for hole in holes {
                        active_slots.push(slots[*hole]);
                        bundle = bundle.add_arg(slots[*hole].to_atom());
                    }
                    bundle.finish()
                }
                LayoutArgument::Group {
                    symbol,
                    prefix,
                    holes,
                    ..
                } => {
                    if zero_groups.contains(&GroupKey::Argument(argument)) {
                        Atom::Zero
                    } else {
                        active_slots.extend(holes.iter().map(|hole| slots[*hole]));
                        let mut group = FunctionBuilder::new(*symbol);
                        for hole in holes {
                            group = group.add_arg(slots[*hole].to_atom());
                        }
                        let group = group.finish();
                        match (
                            prefix == &Atom::num(1),
                            negative_groups.contains(&GroupKey::Argument(argument)),
                        ) {
                            (true, false) | (false, true) => group,
                            (true, true) => &Atom::num(-1) * group,
                            (false, false) => prefix * group,
                        }
                    }
                }
            };
            builder = builder.add_arg(value);
        }
        (builder.finish(), active_slots)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum GroupKey {
    Intrinsic,
    Argument(usize),
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum IncidenceRole {
    Ordered { argument: usize },
    OrderedBundle { argument: usize, member: usize },
    Group { key: GroupKey, kind: SymmetryKind },
}

fn validate_tensor_symmetry<Aind: AbsInd + ParseableAind>(
    value: AtomView<'_>,
) -> Result<(), CanonicalizationError> {
    // Validate parser-recognized metric shorthand after its non-allocating
    // normalization; callers still parse or retain the original expression.
    let validation = value.normalize_dots().metric_shorthand_to_dot();
    let mut stack = vec![validation.as_view()];
    while let Some(value) = stack.pop() {
        match value {
            AtomView::Fun(function) => {
                let head = function.get_symbol();
                if head.has_tag(&SPENSO_TAG.tensor) && SymmetryKind::of(head).is_some() {
                    for (argument, value) in function.iter().enumerate() {
                        if Slot::<LibraryRep, Aind>::try_from(value).is_err() {
                            return Err(CanonicalizationError::InvalidIntrinsicArgument {
                                head,
                                argument,
                                expression: function.as_view().to_owned(),
                            });
                        }
                    }
                } else if head.has_tag(&SPENSO_TAG.tensor) {
                    for (argument, value) in function.iter().enumerate() {
                        let partial_group = TensorLayout::partial_group::<Aind>(
                            value,
                            &function.as_view().to_owned(),
                        )?;
                        if function.get_nargs() != 1
                            && partial_group.is_none()
                            && value.contains_exposed_tensor_topology(StrictTensorFilter::Tagged)
                        {
                            return Err(CanonicalizationError::HiddenTensorTopology {
                                head,
                                argument,
                                expression: function.as_view().to_owned(),
                            });
                        }
                    }
                } else if function.get_nargs() > 1
                    && !function.as_view().is_tensorial(StrictTensorFilter::Tagged)
                    && head != SPENSO_TAG.pure_scalar
                    && head != SPENSO_TAG.scalar
                    && head != SPENSO_TAG.bracket
                    && head != *SYM
                    && head != *ANTISYM
                    && head != *CYCLIC
                {
                    for (argument, value) in function.iter().enumerate() {
                        if value.contains_exposed_tensor_topology(StrictTensorFilter::Tagged) {
                            return Err(CanonicalizationError::HiddenTensorTopology {
                                head,
                                argument,
                                expression: function.as_view().to_owned(),
                            });
                        }
                    }
                }
                stack.extend(function.iter());
            }
            AtomView::Add(add) => stack.extend(add.iter()),
            AtomView::Mul(product) => stack.extend(product.iter()),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                stack.push(base);
                stack.push(exponent);
            }
            AtomView::Num(_) | AtomView::Var(_) => {}
        }
    }
    Ok(())
}
