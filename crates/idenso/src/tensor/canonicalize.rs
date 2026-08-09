//! Signed canonicalization for Spenso symbolic tensor networks.
//!
//! The network remains the semantic owner of products, sums, powers, and
//! enclosing functions. Each ordinary iteration projects its complete operation
//! tree, tensor ports, and index-line components into one graph. Graphica
//! supplies both the canonical labeling and stabilizer generators. Hidden graph
//! data records only reconstruction origins; it does not participate in graph
//! identity.
//!
//! General Young tableaux first apply their exact reduced projector without
//! distributing Products. An eligible lone tensor with distinct external lines
//! sends its declared projected sum through the ordinary whole-root driver and
//! numeric normalization directly. Other eligible tensors become private linear
//! ordered-column carriers sharing one fixed head; an opaque original-head
//! payload preserves declared-head ordering. They decode with deterministic
//! numeric-content and Add-sign normalization, then fully validate and
//! canonically parse the rebuilt Atom into a graph-canonical policy, with no
//! second projector or graph pass. A carrier cycle uses exact post-projector
//! composite iteration and returns its middle graph result.
//! Young-containing Power, a normalization-bearing exposed LocalTensor head or
//! Function anywhere in a Young-containing root, repeated exposed Young heads,
//! and carrier graph- or decoration-orbit-limit failures use the staged path in
//! [`young`]. The root-wide guard includes siblings outside the Young subtree and
//! prevents carrier decode and canonical reparse from passing the graph-rebuilt
//! root back through user normalizers.
//! Strict private metadata promotes carrier slot bundles to the same graph-owned
//! signed columns as declared Young heads without granting block exchange to
//! ordinary same-sized structural groups.
//!
//! Only the successful lone-root direct route may finish as a terminal Atom,
//! which the Atom-facing entry point returns without another parse. Carrier,
//! composite, staged, and ordinary routes retain `CanonicalPolicyNet`; the
//! test-only policy entry point reparses direct terminal output. Factored Young
//! transforms and carrier decode use the same tensor-symmetry and Power-grammar
//! validation and canonical parser as other canonical-policy inputs.

mod driver;
mod group;
mod projection;
mod reconstruct;
mod young;

use std::collections::{BTreeMap, BTreeSet};

use spenso::{
    network::{
        parsing::{AtomStructureExt, StrictTensorFilter},
        tags::SPENSO_TAG,
    },
    shadowing::{ANTISYM, CYCLIC, SYM},
    structure::{
        OrderedStructure, TensorStructure, YoungTableau, YoungTableauClass,
        abstract_index::AIND_SYMBOLS,
        representation::{LibraryRep, LibrarySlot},
        slot::{AbsInd, IsAbstractSlot, ParseableAind, Slot},
    },
};
use symbolica::{
    atom::{Atom, AtomView, FunctionBuilder, Symbol, representation::FunView},
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
        "signed decoration orbit contains at least {observed_at_least} states across {sites} active sites, exceeding limit {limit}"
    )]
    SignedDecorationOrbitLimit {
        observed_at_least: usize,
        limit: usize,
        sites: usize,
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
    #[error("invalid Young-tableau metadata on tensor {head}: {reason}")]
    InvalidYoungTableauMetadata { head: Symbol, reason: String },
    #[error(
        "Young tableau on tensor {head} has arity {expected}, but the tensor has {actual} arguments: {expression}"
    )]
    InvalidYoungTableauArity {
        head: Symbol,
        expected: usize,
        actual: usize,
        expression: Atom,
    },
    #[error(
        "Young tableau on tensor {head} requires a direct slot at argument {argument}: {expression}"
    )]
    InvalidYoungTableauArgument {
        head: Symbol,
        argument: usize,
        expression: Atom,
    },
    #[error(
        "Young tableau on tensor {head} requires argument {argument} to have the same representation, dimension, and orientation as argument 0: {expression}"
    )]
    IncompatibleYoungTableauRepresentation {
        head: Symbol,
        argument: usize,
        expression: Atom,
    },
    #[error(
        "Young projector on tensor {head} requires {requested_actions} full actions, exceeding the internal limit {action_limit}: shape {shape:?}"
    )]
    YoungProjectorSizeLimit {
        head: Symbol,
        shape: Vec<usize>,
        requested_actions: usize,
        action_limit: usize,
    },
    #[error("failed to construct the Young projector on tensor {head}: {reason}")]
    YoungProjectorPlanning { head: Symbol, reason: String },
    #[error(
        "Young straightening requires {requested_terms} live product terms, exceeding the internal limit {term_limit}"
    )]
    YoungExpansionSizeLimit {
        requested_terms: usize,
        term_limit: usize,
    },
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
    young_columns: Vec<Vec<usize>>,
    outer_linear: bool,
    is_young_carrier: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct TensorColor {
    head: semantic::SemanticSymbolKey,
    arguments: Vec<LayoutColorArgument>,
}

type PartialGroup<Aind> = (Symbol, SymmetryKind, Atom, Vec<LibrarySlot<Aind>>);

struct YoungDeclaration<Aind> {
    intrinsic: Option<SymmetryKind>,
    columns: Vec<Vec<usize>>,
    slots: Vec<LibrarySlot<Aind>>,
}

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
        let carrier_declaration = young::YoungColumnCarrier::declaration(function)?;
        let is_young_carrier = carrier_declaration.is_some();
        let young_declaration = if tagged && carrier_declaration.is_none() {
            Self::young_declaration::<Aind>(head, function, &tensor.expression)?
        } else {
            None
        };
        let attribute_intrinsic = tagged.then(|| SymmetryKind::of(head)).flatten();
        let intrinsic = young_declaration
            .as_ref()
            .map(|declaration| declaration.intrinsic)
            .unwrap_or(attribute_intrinsic);
        let root_group = SymmetryKind::structural(head);

        let mut arguments = Vec::with_capacity(function.get_nargs());
        let mut colors = Vec::with_capacity(function.get_nargs());
        let mut expression_slots = Vec::<LibrarySlot<Aind>>::new();
        let mut slot_count = 0;

        if let Some(declaration) = &young_declaration {
            for &slot in &declaration.slots {
                arguments.push(LayoutArgument::DirectSlot(slot_count));
                colors.push(LayoutColorArgument::DirectSlot);
                expression_slots.push(slot);
                slot_count += 1;
            }
        } else if intrinsic.is_some() || root_group.is_some() {
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

        let carrier_columns = if let Some((original, tableau)) = &carrier_declaration {
            let manifest_columns = tableau.columns().collect::<Vec<_>>();
            let mut columns = Vec::with_capacity(manifest_columns.len());
            for (manifest, argument) in manifest_columns.iter().zip(&arguments[1..]) {
                let holes = match (manifest.as_slice(), argument) {
                    ([_], LayoutArgument::DirectSlot(hole)) => vec![*hole],
                    ([_, _, ..], LayoutArgument::SlotBundle { holes })
                        if holes.len() == manifest.len() =>
                    {
                        holes.clone()
                    }
                    _ => {
                        return Err(CanonicalizationError::Projection(format!(
                            "internal Young-column carrier for {original} does not match its declared column shape"
                        )));
                    }
                };
                columns.push(holes);
            }
            columns
        } else {
            Vec::new()
        };
        if let Some((original, _)) = &carrier_declaration
            && let Some(expected) = expression_slots.first().map(|slot| slot.rep())
            && let Some((argument, _)) = expression_slots
                .iter()
                .enumerate()
                .find(|(_, slot)| slot.rep() != expected)
        {
            return Err(
                CanonicalizationError::IncompatibleYoungTableauRepresentation {
                    head: *original,
                    argument,
                    expression: tensor.expression.clone(),
                },
            );
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
            young_columns: young_declaration
                .map(|declaration| declaration.columns)
                .unwrap_or(carrier_columns),
            outer_linear: head.is_linear(),
            is_young_carrier,
        })
    }

    fn young_declaration<Aind: AbsInd + ParseableAind>(
        head: Symbol,
        function: FunView<'_>,
        expression: &Atom,
    ) -> Result<Option<YoungDeclaration<Aind>>, CanonicalizationError> {
        let Some(tableau) = YoungTableau::from_symbol(head).map_err(|error| {
            CanonicalizationError::InvalidYoungTableauMetadata {
                head,
                reason: error.to_string(),
            }
        })?
        else {
            return Ok(None);
        };
        let class = tableau.class();
        let expected = match class {
            YoungTableauClass::Symmetric => Some(SymmetryKind::Symmetric),
            YoungTableauClass::Antisymmetric => Some(SymmetryKind::Antisymmetric),
            YoungTableauClass::General => None,
        };
        let attributes = [
            (head.is_symmetric(), SymmetryKind::Symmetric),
            (head.is_antisymmetric(), SymmetryKind::Antisymmetric),
            (head.is_cyclesymmetric(), SymmetryKind::Cyclic),
        ]
        .into_iter()
        .filter_map(|(present, kind)| present.then_some(kind))
        .collect::<Vec<_>>();
        if expected.is_none_or(|expected| {
            attributes.iter().any(|attribute| *attribute != expected) || attributes.len() > 1
        }) && !attributes.is_empty()
        {
            return Err(CanonicalizationError::InvalidYoungTableauMetadata {
                head,
                reason: format!(
                    "tableau class {class:?} conflicts with intrinsic attributes {attributes:?}"
                ),
            });
        }
        if function.get_nargs() != tableau.rank() {
            return Err(CanonicalizationError::InvalidYoungTableauArity {
                head,
                expected: tableau.rank(),
                actual: function.get_nargs(),
                expression: expression.clone(),
            });
        }

        let mut slots = Vec::<LibrarySlot<Aind>>::with_capacity(tableau.rank());
        for (argument, value) in function.iter().enumerate() {
            let slot = Slot::<LibraryRep, Aind>::try_from(value).map_err(|_| {
                CanonicalizationError::InvalidYoungTableauArgument {
                    head,
                    argument,
                    expression: expression.clone(),
                }
            })?;
            if let Some(first) = slots.first()
                && slot.rep() != first.rep()
            {
                return Err(
                    CanonicalizationError::IncompatibleYoungTableauRepresentation {
                        head,
                        argument,
                        expression: expression.clone(),
                    },
                );
            }
            slots.push(slot);
        }
        let columns = if class == YoungTableauClass::General {
            // Keep singleton columns as well so repeated height-one columns
            // can form an unsigned exchange block. A unique singleton remains
            // an ordinary ordered slot.
            tableau.columns().collect()
        } else {
            Vec::new()
        };

        Ok(Some(YoungDeclaration {
            intrinsic: expected.filter(|_| tableau.rank() > 1),
            columns,
            slots,
        }))
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
        if let Some((column, _)) = self.young_columns.iter().enumerate().find(|(_, holes)| {
            holes.contains(&flat_slot)
                && (holes.len() > 1
                    || self
                        .young_columns
                        .iter()
                        .filter(|candidate| candidate.len() == holes.len())
                        .nth(1)
                        .is_some())
        }) {
            return IncidenceRole::Group {
                key: GroupKey::YoungColumn(column),
                kind: SymmetryKind::Antisymmetric,
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
        if let Some(position) = self
            .young_columns
            .iter()
            .find_map(|column| column.iter().position(|hole| *hole == flat_slot))
        {
            return position;
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
            .chain(
                self.young_columns
                    .iter()
                    .enumerate()
                    .filter(|(_, holes)| holes.len() > 1)
                    .map(|(column, _)| (GroupKey::YoungColumn(column), true)),
            )
            .chain(self.arguments.iter().enumerate().filter_map(
                |(argument, layout)| match layout {
                    LayoutArgument::Group {
                        kind: SymmetryKind::Antisymmetric,
                        ..
                    } => {
                        let key = GroupKey::Argument(argument);
                        Some((key, self.group_lifts(key)))
                    }
                    _ => None,
                },
            ))
    }

    /// Manifest Young columns grouped into repeated-height blocks.
    ///
    /// Only these columns gain a visible owner node. Ordinary structural
    /// groups and a unique-height Young column retain their existing local
    /// incidence identity.
    fn exchangeable_young_column_blocks(&self) -> BTreeMap<usize, Vec<usize>> {
        let mut blocks = BTreeMap::<usize, Vec<usize>>::new();
        for (column, holes) in self.young_columns.iter().enumerate() {
            blocks.entry(holes.len()).or_default().push(column);
        }
        blocks.retain(|_, columns| columns.len() > 1);
        blocks
    }

    fn group_holes(&self, key: GroupKey) -> Option<Vec<usize>> {
        match key {
            GroupKey::Intrinsic => Some(
                self.arguments
                    .iter()
                    .filter_map(|argument| match argument {
                        LayoutArgument::DirectSlot(hole) => Some(*hole),
                        _ => None,
                    })
                    .collect(),
            ),
            GroupKey::Argument(argument) => match self.arguments.get(argument) {
                Some(LayoutArgument::Group { holes, .. }) => Some(holes.clone()),
                _ => None,
            },
            GroupKey::YoungColumn(column) => self.young_columns.get(column).cloned(),
        }
    }

    fn group_lifts(&self, key: GroupKey) -> bool {
        match key {
            GroupKey::Intrinsic | GroupKey::YoungColumn(_) => true,
            GroupKey::Argument(_) => self.outer_linear,
        }
    }

    fn rebuild<Aind: AbsInd + ParseableAind>(
        &self,
        slots: &[LibrarySlot<Aind>],
        zero_groups: &BTreeSet<GroupKey>,
        negative_groups: &BTreeSet<GroupKey>,
    ) -> (Atom, Vec<LibrarySlot<Aind>>) {
        if zero_groups.iter().any(|key| self.group_lifts(*key)) {
            return (Atom::Zero, Vec::new());
        }
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
    YoungColumn(usize),
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
                if head.has_tag(&SPENSO_TAG.tensor) {
                    let expression = function.as_view().to_owned();
                    if TensorLayout::young_declaration::<Aind>(head, function, &expression)?
                        .is_none()
                    {
                        if SymmetryKind::of(head).is_some() {
                            for (argument, value) in function.iter().enumerate() {
                                if Slot::<LibraryRep, Aind>::try_from(value).is_err() {
                                    return Err(CanonicalizationError::InvalidIntrinsicArgument {
                                        head,
                                        argument,
                                        expression,
                                    });
                                }
                            }
                        } else {
                            for (argument, value) in function.iter().enumerate() {
                                let partial_group =
                                    TensorLayout::partial_group::<Aind>(value, &expression)?;
                                if function.get_nargs() != 1
                                    && partial_group.is_none()
                                    && value.contains_exposed_tensor_topology(
                                        StrictTensorFilter::Tagged,
                                    )
                                {
                                    return Err(CanonicalizationError::HiddenTensorTopology {
                                        head,
                                        argument,
                                        expression,
                                    });
                                }
                            }
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
                if head == SPENSO_TAG.bracket
                    || function
                        .as_view()
                        .contains_exposed_tensor_topology(StrictTensorFilter::Tagged)
                {
                    stack.extend(function.iter());
                }
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
