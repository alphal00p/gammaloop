use std::{
    collections::{BTreeMap, BTreeSet, VecDeque},
    fmt,
};

use spenso::{
    contraction::Contract,
    network::{
        Network, NetworkState,
        graph::{NAdd, NMul},
    },
    structure::{
        OrderedStructure,
        representation::{LibraryRep, LibrarySlot, Representation},
        slot::{AbsInd, DummyAind, IsAbstractSlot, ParseableAind},
    },
};
use symbolica::atom::{Atom, AtomView};

use super::{
    CanonicalizationError, GroupKey, IncidenceRole, LayoutArgument, SymmetryKind,
    driver::execute_atom,
    group::{SignSiteFrame, SignedAction, SignedGroup, transport_site_frames},
    projection::{
        CanonicalProjection, ExpressionKind, PowerBoundaryTarget, PowerCopyDescriptor,
        UnifiedNodeColor,
    },
    semantic,
};
use crate::tensor::{SymbolicNet, SymbolicTensor};

pub(super) struct Reconstruction<Aind: AbsInd> {
    pub(super) network: SymbolicNet<Aind>,
    pub(super) retry_reason: Option<RetryReason>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum RetryReason {
    ZeroSubstitution,
    ResidualPowerPhase,
    SignedPayloadEdit,
    NormalizedTensorPayloadEdit,
    VisibleClassMergeSplit,
    ExecutionTopologyChange,
    IncompleteStabilityCertificate,
}

impl fmt::Display for RetryReason {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::ZeroSubstitution => "zero substitution",
            Self::ResidualPowerPhase => "residual Power phase",
            Self::SignedPayloadEdit => "signed payload edit",
            Self::NormalizedTensorPayloadEdit => "normalized tensor payload edit",
            Self::VisibleClassMergeSplit => "visible-class merge/split",
            Self::ExecutionTopologyChange => "execution-topology change",
            Self::IncompleteStabilityCertificate => "incomplete stability certificate",
        })
    }
}

pub(super) struct PreparedReconstruction {
    analysis: SignedAnalysis,
    pub(super) identity: SignedProblemIdentity,
}

pub(super) struct DummyAllocator<Aind: AbsInd> {
    positions: BTreeMap<(Representation<LibraryRep>, usize), usize>,
    values: Vec<Aind>,
}

impl<Aind: AbsInd> DummyAllocator<Aind> {
    pub(super) fn new() -> Self {
        Self {
            positions: BTreeMap::new(),
            values: Vec::new(),
        }
    }

    fn get<F>(
        &mut self,
        group: Representation<LibraryRep>,
        ordinal: usize,
        new_dummy: &mut F,
    ) -> Aind
    where
        F: FnMut(usize) -> Aind,
    {
        let position = if let Some(position) = self.positions.get(&(group, ordinal)) {
            *position
        } else {
            let position = self.positions.len();
            debug_assert!(position <= self.values.len());
            if position == self.values.len() {
                self.values.push(new_dummy(position));
            }
            self.positions.insert((group, ordinal), position);
            position
        };
        self.values[position]
    }
}

#[derive(Clone)]
struct SiteMeta {
    occurrence: usize,
    expression: usize,
    key: GroupKey,
    intrinsic: bool,
    lifts: bool,
}

#[derive(Clone)]
struct SignedAnalysis {
    group: SignedGroup,
    sites: Vec<SiteMeta>,
    phases: Vec<bool>,
    zero_groups: BTreeMap<usize, BTreeSet<GroupKey>>,
    zero_expressions: BTreeSet<usize>,
    zero_magnitudes: BTreeSet<usize>,
    singular_expressions: BTreeSet<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct ScopeClassKey {
    root: usize,
    kind: UnifiedNodeColor,
    boundary: Vec<(usize, semantic::SemanticAtomKey)>,
    value: semantic::SemanticAtomKey,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct SinkKey {
    kind: u8,
    vertex: usize,
    path: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct SignedProblemIdentity {
    class: ScopeClassKey,
    realization: Vec<SinkKey>,
    phases: Vec<bool>,
    zero_groups: Vec<(usize, Vec<GroupKey>)>,
    zero_expressions: Vec<usize>,
    zero_magnitudes: Vec<usize>,
    singular_expressions: Vec<usize>,
}

#[derive(Clone)]
enum SinkTarget {
    Nested { occurrence: usize, key: GroupKey },
    Leaf { expression: usize },
    Result { expression: usize },
}

#[derive(Clone, Default)]
struct NormalizedDecoration {
    sinks: BTreeMap<SinkKey, SinkTarget>,
}

impl NormalizedDecoration {
    fn signs_expression(&self, expression: usize) -> bool {
        self.sinks.values().any(|target| match target {
            SinkTarget::Leaf {
                expression: candidate,
            }
            | SinkTarget::Result {
                expression: candidate,
            } => *candidate == expression,
            SinkTarget::Nested { .. } => false,
        })
    }

    fn signs_nested(&self, occurrence: usize, key: GroupKey) -> bool {
        self.sinks.values().any(|target| {
            matches!(
                target,
                SinkTarget::Nested {
                    occurrence: candidate,
                    key: candidate_key,
                } if *candidate == occurrence && *candidate_key == key
            )
        })
    }
}

#[derive(Clone)]
struct ScopeDescriptor {
    root: usize,
    kind: UnifiedNodeColor,
    boundary: Vec<(usize, semantic::SemanticAtomKey)>,
}

struct ClassifiedScope {
    actual: ScopeClassKey,
    unsigned: ScopeClassKey,
    orientation: bool,
    realization: Vec<SinkKey>,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct PowerCopyClass {
    state: u8,
    seam: Vec<Vec<(IncidenceRole, semantic::SemanticAtomKey)>>,
    body: semantic::SemanticAtomKey,
}

struct ClassifiedPowerCopy {
    class: PowerCopyClass,
    negative: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum VisibleSyntaxClass {
    Scalar(semantic::SemanticAtomKey),
    Tensor {
        color: super::TensorColor,
        slots: Vec<(semantic::SemanticAtomKey, usize)>,
    },
    Product(Vec<usize>),
    Sum(Vec<usize>),
    Neg(usize),
    Function(semantic::SemanticSymbolKey, usize),
}

pub(super) fn prepare_reconstruction<Aind>(
    projection: &CanonicalProjection<Aind>,
) -> Result<PreparedReconstruction, CanonicalizationError>
where
    Aind: AbsInd + DummyAind + ParseableAind,
    SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
{
    let (analysis, identity) = SignedAnalysis::new(projection)?;
    Ok(PreparedReconstruction { analysis, identity })
}

pub(super) fn reconstruct<Aind, F>(
    projection: &CanonicalProjection<Aind>,
    prepared: PreparedReconstruction,
    dummy_allocator: &mut DummyAllocator<Aind>,
    new_dummy: &mut F,
) -> Result<Reconstruction<Aind>, CanonicalizationError>
where
    Aind: AbsInd + DummyAind + ParseableAind,
    F: FnMut(usize) -> Aind,
    SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
{
    let analysis = prepared.analysis;
    // Callback values are request-local, but representation/ordinal positions
    // describe only this iteration's surviving canonical line order.
    dummy_allocator.positions.clear();
    let mut rebuilder = Rebuilder::new(projection, analysis, dummy_allocator, new_dummy)?;
    let network = rebuilder.expression(projection.root_expression)?;
    let retry_reason = if !rebuilder.analysis.zero_expressions.is_empty()
        || !rebuilder.analysis.zero_groups.is_empty()
        || !rebuilder.analysis.zero_magnitudes.is_empty()
    {
        Some(RetryReason::ZeroSubstitution)
    } else if rebuilder.residual_power_phase {
        Some(RetryReason::ResidualPowerPhase)
    } else if rebuilder.analysis.phases.iter().any(|phase| *phase) {
        Some(RetryReason::SignedPayloadEdit)
    } else if rebuilder.payload_edited {
        Some(RetryReason::NormalizedTensorPayloadEdit)
    } else {
        rebuilder.certifies_selected_transport().err()
    };
    Ok(Reconstruction {
        network,
        retry_reason,
    })
}

impl SignedAnalysis {
    fn new<Aind>(
        projection: &CanonicalProjection<Aind>,
    ) -> Result<(Self, SignedProblemIdentity), CanonicalizationError>
    where
        Aind: AbsInd + DummyAind + ParseableAind,
        SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
    {
        let mut source_sites = Vec::new();
        let mut canonical_sites = Vec::<(SignSiteFrame, SiteMeta)>::new();
        for (occurrence_index, occurrence) in projection.tensors.iter().enumerate() {
            for (key, lifts) in occurrence.layout.antisymmetric_groups() {
                let mut source_members = occurrence
                    .ports
                    .iter()
                    .filter(|port| Self::port_in_group(&port.role, key))
                    .collect::<Vec<_>>();
                source_members.sort_unstable_by_key(|port| port.member);
                let source_frame = SignSiteFrame {
                    owner: occurrence.header,
                    layout_path: Self::layout_path(key),
                    members: source_members.iter().map(|port| port.vertex).collect(),
                };
                let mut canonical_members = source_members
                    .iter()
                    .map(|port| projection.vertex_map[port.vertex])
                    .collect::<Vec<_>>();
                canonical_members.sort_unstable();
                source_sites.push(source_frame);
                canonical_sites.push((
                    SignSiteFrame {
                        owner: projection.vertex_map[occurrence.header],
                        layout_path: Self::layout_path(key),
                        members: canonical_members,
                    },
                    SiteMeta {
                        occurrence: occurrence_index,
                        expression: occurrence.expression,
                        key,
                        intrinsic: key == GroupKey::Intrinsic,
                        lifts,
                    },
                ));
            }
        }
        canonical_sites.sort_by(|(left, _), (right, _)| {
            left.owner
                .cmp(&right.owner)
                .then(left.layout_path.cmp(&right.layout_path))
        });
        let (target_frames, sites): (Vec<_>, Vec<_>) = canonical_sites.into_iter().unzip();
        let images = transport_site_frames(&projection.vertex_map, &source_sites, &target_frames)
            .map_err(|error| CanonicalizationError::Projection(error.to_string()))?;
        let mut phases = vec![false; sites.len()];
        for image in images {
            phases[image.target] ^= image.odd;
        }
        let group = SignedGroup::from_graphica(
            projection.graph.nodes().len(),
            &target_frames,
            &projection.orbit_generators,
        )
        .map_err(|error| CanonicalizationError::Projection(error.to_string()))?;

        let mut analysis = Self {
            group,
            sites,
            phases,
            zero_groups: BTreeMap::new(),
            zero_expressions: BTreeSet::new(),
            zero_magnitudes: BTreeSet::new(),
            singular_expressions: BTreeSet::new(),
        };
        analysis.find_site_zeros(projection)?;
        let mut seen = BTreeSet::new();
        analysis.find_expression_zeros(projection, projection.root_expression, &mut seen)?;
        let identity = analysis.canonicalize_decoration(projection)?;
        Ok((analysis, identity))
    }

    fn normalized_decoration<Aind: AbsInd>(
        &self,
        projection: &CanonicalProjection<Aind>,
        scope: usize,
        phases: &[bool],
    ) -> Result<NormalizedDecoration, CanonicalizationError> {
        fn visit<Aind: AbsInd>(
            projection: &CanonicalProjection<Aind>,
            expression: usize,
            descendants: &mut BTreeSet<usize>,
            parents: &mut BTreeMap<usize, usize>,
        ) -> Result<(), CanonicalizationError> {
            if !descendants.insert(expression) {
                return Ok(());
            }
            let occurrence = &projection.expressions[expression];
            let children = match &occurrence.kind {
                ExpressionKind::Power { copies, .. } => copies.as_slice(),
                _ => occurrence.children.as_slice(),
            };
            for &child in children {
                if parents.insert(child, expression).is_some() {
                    return Err(CanonicalizationError::AmbiguousSignScope);
                }
                visit(projection, child, descendants, parents)?;
            }
            Ok(())
        }

        fn transparent<Aind: AbsInd>(
            projection: &CanonicalProjection<Aind>,
            expression: usize,
        ) -> bool {
            matches!(
                projection.expressions[expression].kind,
                ExpressionKind::Product | ExpressionKind::Neg
            ) || matches!(
                projection.expressions[expression].kind,
                ExpressionKind::Function(function) if function.is_linear()
            )
        }

        fn leaf_candidates<Aind: AbsInd>(
            projection: &CanonicalProjection<Aind>,
            zero_expressions: &BTreeSet<usize>,
            expression: usize,
            candidates: &mut Vec<(u8, usize, usize)>,
        ) {
            if zero_expressions.contains(&expression) {
                return;
            }
            let occurrence = &projection.expressions[expression];
            match occurrence.kind {
                ExpressionKind::Scalar(_) => {
                    candidates.push((0, projection.vertex_map[occurrence.root], expression))
                }
                ExpressionKind::Tensor(_) => {
                    candidates.push((1, projection.vertex_map[occurrence.root], expression))
                }
                ExpressionKind::Product | ExpressionKind::Neg => {
                    for &child in &occurrence.children {
                        leaf_candidates(projection, zero_expressions, child, candidates);
                    }
                }
                ExpressionKind::Function(function) if function.is_linear() => {
                    for &child in &occurrence.children {
                        leaf_candidates(projection, zero_expressions, child, candidates);
                    }
                }
                ExpressionKind::Sum
                | ExpressionKind::Function(_)
                | ExpressionKind::Power { .. } => {}
            }
        }

        fn toggle(decoration: &mut NormalizedDecoration, key: SinkKey, target: SinkTarget) {
            use std::collections::btree_map::Entry;
            match decoration.sinks.entry(key) {
                Entry::Vacant(entry) => {
                    entry.insert(target);
                }
                Entry::Occupied(entry) => {
                    entry.remove();
                }
            }
        }

        let mut descendants = BTreeSet::new();
        let mut parents = BTreeMap::new();
        visit(projection, scope, &mut descendants, &mut parents)?;

        let mut decoration = NormalizedDecoration::default();
        let mut whole_values = BTreeMap::<usize, bool>::new();
        for (site, metadata) in self.sites.iter().enumerate() {
            if !phases[site] || !descendants.contains(&metadata.expression) {
                continue;
            }
            if metadata.lifts {
                *whole_values.entry(metadata.expression).or_default() ^= true;
            } else {
                let occurrence = &projection.tensors[metadata.occurrence];
                toggle(
                    &mut decoration,
                    SinkKey {
                        kind: 0,
                        vertex: projection.vertex_map[occurrence.header],
                        path: Self::layout_path(metadata.key),
                    },
                    SinkTarget::Nested {
                        occurrence: metadata.occurrence,
                        key: metadata.key,
                    },
                );
            }
        }

        let mut transparent_regions = BTreeMap::<usize, bool>::new();
        for (mut expression, odd) in whole_values {
            if !odd {
                continue;
            }
            while let Some(parent) = parents.get(&expression).copied() {
                if !transparent(projection, parent) {
                    break;
                }
                expression = parent;
            }
            *transparent_regions.entry(expression).or_default() ^= true;
        }

        for (expression, odd) in transparent_regions {
            if !odd {
                continue;
            }
            let mut candidates = Vec::new();
            leaf_candidates(
                projection,
                &self.zero_expressions,
                expression,
                &mut candidates,
            );
            candidates.sort_unstable();
            if let Some(&(kind, vertex, leaf)) = candidates.first() {
                toggle(
                    &mut decoration,
                    SinkKey {
                        kind: kind + 1,
                        vertex,
                        path: Vec::new(),
                    },
                    SinkTarget::Leaf { expression: leaf },
                );
            } else {
                toggle(
                    &mut decoration,
                    SinkKey {
                        kind: 3,
                        vertex: projection.vertex_map[projection.expressions[expression].root],
                        path: Vec::new(),
                    },
                    SinkTarget::Result { expression },
                );
            }
        }
        Ok(decoration)
    }

    fn canonicalize_decoration<Aind>(
        &mut self,
        projection: &CanonicalProjection<Aind>,
    ) -> Result<SignedProblemIdentity, CanonicalizationError>
    where
        Aind: AbsInd + DummyAind + ParseableAind,
        SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
    {
        let descriptor = Self::scope_descriptor(projection, projection.root_expression);
        let mut active_expressions = BTreeSet::new();
        let mut pending = vec![projection.root_expression];
        while let Some(expression) = pending.pop() {
            if self.zero_expressions.contains(&expression)
                || self.zero_magnitudes.contains(&expression)
                || !active_expressions.insert(expression)
            {
                continue;
            }
            let occurrence = &projection.expressions[expression];
            match &occurrence.kind {
                ExpressionKind::Power { copies, .. } => pending.extend(copies),
                _ => pending.extend(&occurrence.children),
            }
        }
        let active_sites = self
            .sites
            .iter()
            .map(|site| {
                active_expressions.contains(&site.expression)
                    && !self
                        .zero_groups
                        .get(&site.occurrence)
                        .is_some_and(|groups| groups.contains(&site.key))
            })
            .collect::<Vec<_>>();
        let orbit = self
            .group
            .decoration_orbit(&self.phases, &active_sites)
            .map_err(|error| CanonicalizationError::Projection(error.to_string()))?;
        let mut best = None;
        for phases in orbit {
            let class =
                self.classify_scope(projection, projection.root_expression, &descriptor, &phases)?;
            let candidate = (class.actual, class.realization, phases);
            if best.as_ref().is_none_or(|current| &candidate < current) {
                best = Some(candidate);
            }
        }
        let (class, realization, phases) =
            best.expect("an affine decoration orbit contains its input");
        self.phases = phases.clone();
        let mut zero_groups = self
            .zero_groups
            .iter()
            .map(|(tensor, groups)| {
                (
                    projection.vertex_map[projection.tensors[*tensor].header],
                    groups.iter().copied().collect(),
                )
            })
            .collect::<Vec<_>>();
        zero_groups.sort();
        let canonical_expressions = |expressions: &BTreeSet<usize>| {
            let mut expressions = expressions
                .iter()
                .map(|expression| projection.vertex_map[projection.expressions[*expression].root])
                .collect::<Vec<_>>();
            expressions.sort_unstable();
            expressions
        };
        let mut zero_magnitudes = self
            .zero_magnitudes
            .iter()
            .map(|expression| {
                projection
                    .power_descriptor(*expression)
                    .expect("zero Power magnitudes have validated descriptors")
                    .magnitude
            })
            .collect::<Vec<_>>();
        zero_magnitudes.sort_unstable();
        Ok(SignedProblemIdentity {
            class,
            realization,
            phases,
            zero_groups,
            zero_expressions: canonical_expressions(&self.zero_expressions),
            zero_magnitudes,
            singular_expressions: canonical_expressions(&self.singular_expressions),
        })
    }

    fn scope_descriptor<Aind: AbsInd + ParseableAind>(
        projection: &CanonicalProjection<Aind>,
        expression: usize,
    ) -> ScopeDescriptor {
        let occurrence = &projection.expressions[expression];
        let root = projection.vertex_map[occurrence.root];
        let tensors = Self::tensor_descendants(projection, &[expression]);
        let lines = Self::lines_for_tensors(projection, &tensors);
        let mut boundary = lines
            .into_iter()
            .filter_map(|line| {
                let occurs_outside =
                    projection
                        .tensors
                        .iter()
                        .enumerate()
                        .any(|(tensor, occurrence)| {
                            !tensors.contains(&tensor)
                                && occurrence.ports.iter().any(|port| port.line == line)
                        });
                let line_data = &projection.lines[line];
                (occurs_outside || line_data.external.is_some()).then(|| {
                    let key = line_data.external.map_or_else(
                        || semantic::representation_key(line_data.group),
                        semantic::slot_key,
                    );
                    (projection.vertex_map[line_data.vertex], key)
                })
            })
            .collect::<Vec<_>>();
        boundary.sort();
        ScopeDescriptor {
            root,
            kind: projection.graph.node(root).data.data.clone(),
            boundary,
        }
    }

    fn classify_scope<Aind>(
        &self,
        projection: &CanonicalProjection<Aind>,
        expression: usize,
        descriptor: &ScopeDescriptor,
        phases: &[bool],
    ) -> Result<ClassifiedScope, CanonicalizationError>
    where
        Aind: AbsInd + DummyAind + ParseableAind,
        SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
    {
        let mut analysis = self.clone();
        analysis.phases = phases.to_vec();
        let mut dummy_allocator = DummyAllocator::new();
        let mut new_dummy = Aind::new_dummy_at;
        let mut rebuilder = Rebuilder::with_scope(
            projection,
            analysis,
            expression,
            &mut dummy_allocator,
            &mut new_dummy,
            true,
        )?;
        let network = rebuilder.expression(expression)?;
        let realization = rebuilder.realized_sinks.iter().cloned().collect();
        #[cfg(test)]
        super::driver::record_temporary_scope_execution();
        let atom = execute_atom(network).map_err(|error| match error {
            CanonicalizationError::Execution { error, .. } => CanonicalizationError::Execution {
                scope: format!("canonical scope {}", descriptor.root),
                error,
            },
            error => error,
        })?;
        let direct = semantic::SemanticAtomKey::new(atom.as_view());
        let negative = semantic::SemanticAtomKey::new((-atom).as_view());
        let actual = ScopeClassKey {
            root: descriptor.root,
            kind: descriptor.kind.clone(),
            boundary: descriptor.boundary.clone(),
            value: direct.clone(),
        };
        let (value, orientation) = if negative < direct {
            (negative, true)
        } else {
            (direct, false)
        };
        Ok(ClassifiedScope {
            actual,
            unsigned: ScopeClassKey {
                root: descriptor.root,
                kind: descriptor.kind.clone(),
                boundary: descriptor.boundary.clone(),
                value,
            },
            orientation,
            realization,
        })
    }

    fn scope_has_odd_class_stabilizer<Aind>(
        &self,
        projection: &CanonicalProjection<Aind>,
        expression: usize,
    ) -> Result<bool, CanonicalizationError>
    where
        Aind: AbsInd + DummyAind + ParseableAind,
        SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
    {
        let descriptor = Self::scope_descriptor(projection, expression);
        let mut fixed = Vec::with_capacity(descriptor.boundary.len() + 1);
        fixed.push(descriptor.root);
        fixed.extend(descriptor.boundary.iter().map(|(line, _)| *line));
        let group = self
            .group
            .pointwise_vertex_stabilizer(&fixed)
            .map_err(|error| CanonicalizationError::Projection(error.to_string()))?;
        self.group_has_odd_class_stabilizer(&group, |action| {
            let phases = action
                .apply_site_phases(&self.phases)
                .map_err(|error| CanonicalizationError::Projection(error.to_string()))?;
            let class = self.classify_scope(projection, expression, &descriptor, &phases)?;
            Ok((class.unsigned, class.orientation))
        })
    }

    fn magnitude_has_odd_class_stabilizer<Aind>(
        &self,
        projection: &CanonicalProjection<Aind>,
        expression: usize,
        copies: &[usize],
    ) -> Result<bool, CanonicalizationError>
    where
        Aind: AbsInd + DummyAind + ParseableAind,
        SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
    {
        let power = projection.power_descriptor(expression).ok_or_else(|| {
            CanonicalizationError::PowerReconstruction {
                expression,
                reason: "validated magnitude descriptor is missing".into(),
            }
        })?;
        let descriptor = ScopeDescriptor {
            root: power.magnitude,
            kind: projection.graph.node(power.magnitude).data.data.clone(),
            ..Self::scope_descriptor(projection, expression)
        };
        let mut fixed = Vec::with_capacity(descriptor.boundary.len() + 1);
        fixed.push(descriptor.root);
        fixed.extend(descriptor.boundary.iter().map(|(line, _)| *line));
        let group = self
            .group
            .pointwise_vertex_stabilizer(&fixed)
            .map_err(|error| CanonicalizationError::Projection(error.to_string()))?;
        self.group_has_odd_class_stabilizer(&group, |action| {
            let phases = action
                .apply_site_phases(&self.phases)
                .map_err(|error| CanonicalizationError::Projection(error.to_string()))?;
            let class =
                self.classify_magnitude(projection, expression, copies, &descriptor, &phases)?;
            Ok((class.unsigned, class.orientation))
        })
    }

    fn classify_magnitude<Aind>(
        &self,
        projection: &CanonicalProjection<Aind>,
        expression: usize,
        copies: &[usize],
        descriptor: &ScopeDescriptor,
        phases: &[bool],
    ) -> Result<ClassifiedScope, CanonicalizationError>
    where
        Aind: AbsInd + DummyAind + ParseableAind,
        SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
    {
        let mut analysis = self.clone();
        analysis.phases = phases.to_vec();
        let mut dummy_allocator = DummyAllocator::new();
        let mut new_dummy = Aind::new_dummy_at;
        let mut rebuilder = Rebuilder::with_scope(
            projection,
            analysis,
            expression,
            &mut dummy_allocator,
            &mut new_dummy,
            false,
        )?;
        let mut copies = copies.to_vec();
        copies
            .sort_unstable_by_key(|copy| projection.vertex_map[projection.expressions[*copy].root]);
        let mut factors = Vec::with_capacity(copies.len());
        for copy in copies {
            factors.push(rebuilder.expression(copy)?);
        }
        let realization = rebuilder.realized_sinks.iter().cloned().collect();
        let mut factors = factors.into_iter();
        let magnitude = factors
            .next()
            .map_or_else(Network::one, |first| first.n_mul(factors));
        #[cfg(test)]
        super::driver::record_temporary_scope_execution();
        let atom = execute_atom(magnitude).map_err(|error| match error {
            CanonicalizationError::Execution { error, .. } => CanonicalizationError::Execution {
                scope: format!("Power {expression} magnitude {}", descriptor.root),
                error,
            },
            error => error,
        })?;
        let direct = semantic::SemanticAtomKey::new(atom.as_view());
        let negative = semantic::SemanticAtomKey::new((-atom).as_view());
        let actual = ScopeClassKey {
            root: descriptor.root,
            kind: descriptor.kind.clone(),
            boundary: descriptor.boundary.clone(),
            value: direct.clone(),
        };
        let (value, orientation) = if negative < direct {
            (negative, true)
        } else {
            (direct, false)
        };
        Ok(ClassifiedScope {
            actual,
            unsigned: ScopeClassKey {
                root: descriptor.root,
                kind: descriptor.kind.clone(),
                boundary: descriptor.boundary.clone(),
                value,
            },
            orientation,
            realization,
        })
    }

    fn group_has_odd_class_stabilizer(
        &self,
        group: &SignedGroup,
        classify: impl Fn(&SignedAction) -> Result<(ScopeClassKey, bool), CanonicalizationError>,
    ) -> Result<bool, CanonicalizationError> {
        #[derive(Clone)]
        struct Representative {
            action: SignedAction,
            negative: bool,
        }

        let generators = group.symmetric_generators();
        if generators.is_empty() {
            return Ok(false);
        }

        let identity = group.identity();
        let (base_key, base_negative) = classify(&identity)?;
        let mut representatives = BTreeMap::from([(
            base_key.clone(),
            Representative {
                action: identity,
                negative: base_negative,
            },
        )]);
        let mut queue = VecDeque::from([base_key.clone()]);
        while let Some(key) = queue.pop_front() {
            let representative = representatives[&key].action.clone();
            for generator in &generators {
                let candidate = generator
                    .compose(&representative)
                    .map_err(|error| CanonicalizationError::Projection(error.to_string()))?;
                let (class, negative) = classify(&candidate)?;
                if let std::collections::btree_map::Entry::Vacant(entry) =
                    representatives.entry(class)
                {
                    queue.push_back(entry.key().clone());
                    entry.insert(Representative {
                        action: candidate,
                        negative,
                    });
                }
            }
        }

        for representative in representatives.values() {
            for generator in &generators {
                let candidate = generator
                    .compose(&representative.action)
                    .map_err(|error| CanonicalizationError::Projection(error.to_string()))?;
                let (target_class, target_negative) = classify(&candidate)?;
                let Some(target) = representatives.get(&target_class) else {
                    return Err(CanonicalizationError::AmbiguousSignScope);
                };
                let schreier = target
                    .action
                    .inverse()
                    .compose(&candidate)
                    .map_err(|error| CanonicalizationError::Projection(error.to_string()))?;
                let (stabilized, stabilized_negative) = classify(&schreier)?;
                if stabilized != base_key {
                    return Err(CanonicalizationError::AmbiguousSignScope);
                }
                if stabilized_negative != representatives[&base_key].negative {
                    return Ok(true);
                }
                if target_negative ^ target.negative
                    != stabilized_negative ^ representatives[&base_key].negative
                {
                    return Err(CanonicalizationError::AmbiguousSignScope);
                }
            }
        }
        Ok(false)
    }

    fn find_site_zeros<Aind: AbsInd>(
        &mut self,
        projection: &CanonicalProjection<Aind>,
    ) -> Result<(), CanonicalizationError> {
        for (site, metadata) in self.sites.iter().enumerate() {
            let occurrence = &projection.tensors[metadata.occurrence];
            let mut fixed = vec![projection.vertex_map[occurrence.header]];
            fixed.extend(
                occurrence
                    .ports
                    .iter()
                    .filter(|port| Self::port_in_group(&port.role, metadata.key))
                    .map(|port| projection.vertex_map[projection.lines[port.line].vertex]),
            );
            fixed.sort_unstable();
            fixed.dedup();
            let stabilizer = self
                .group
                .pointwise_vertex_stabilizer(&fixed)
                .map_err(|error| CanonicalizationError::Projection(error.to_string()))?;
            if stabilizer
                .has_odd_site_stabilizer(site)
                .map_err(|error| CanonicalizationError::Projection(error.to_string()))?
            {
                self.zero_groups
                    .entry(metadata.occurrence)
                    .or_default()
                    .insert(metadata.key);
            }
        }
        Ok(())
    }

    fn find_expression_zeros<Aind>(
        &mut self,
        projection: &CanonicalProjection<Aind>,
        expression: usize,
        seen: &mut BTreeSet<usize>,
    ) -> Result<bool, CanonicalizationError>
    where
        Aind: AbsInd + DummyAind + ParseableAind,
        SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
    {
        if !seen.insert(expression) {
            return Ok(self.zero_expressions.contains(&expression));
        }
        let occurrence = &projection.expressions[expression];
        let mut child_zero = Vec::new();
        for &child in &occurrence.children {
            child_zero.push(self.find_expression_zeros(projection, child, seen)?);
        }
        if let ExpressionKind::Power { copies, .. } = &occurrence.kind {
            child_zero.clear();
            for &copy in copies {
                child_zero.push(self.find_expression_zeros(projection, copy, seen)?);
            }
        }

        let child_singular = match &occurrence.kind {
            ExpressionKind::Power { copies, .. } => copies
                .iter()
                .any(|copy| self.singular_expressions.contains(copy)),
            _ => occurrence
                .children
                .iter()
                .any(|child| self.singular_expressions.contains(child)),
        };
        let mut singular = child_singular;
        let mut zero = match &occurrence.kind {
            ExpressionKind::Scalar(atom) => atom.as_view().is_zero(),
            ExpressionKind::Tensor(tensor_index) => {
                let tensor = &projection.tensors[*tensor_index];
                self.zero_groups.get(tensor_index).is_some_and(|groups| {
                    groups.contains(&GroupKey::Intrinsic)
                        || (tensor.layout.outer_linear && !groups.is_empty())
                })
            }
            ExpressionKind::Product => !child_singular && child_zero.iter().any(|zero| *zero),
            ExpressionKind::Sum => !child_zero.is_empty() && child_zero.iter().all(|zero| *zero),
            ExpressionKind::Neg => child_zero.first().copied().unwrap_or(false),
            ExpressionKind::Function(function) => {
                function.is_linear() && child_zero.first().copied().unwrap_or(false)
            }
            ExpressionKind::Power {
                exponent,
                magnitude,
                copies,
            } => {
                let canonical_magnitude = projection
                    .power_descriptor(expression)
                    .ok_or_else(|| {
                        CanonicalizationError::Projection(format!(
                            "materialized Power {expression} has no validated copy descriptor"
                        ))
                    })?
                    .magnitude;
                debug_assert_eq!(canonical_magnitude, projection.vertex_map[*magnitude]);
                let child_is_zero = child_zero.iter().any(|zero| *zero);
                let odd = !child_singular
                    && !child_is_zero
                    && self.magnitude_has_odd_class_stabilizer(projection, expression, copies)?;
                let magnitude_zero = !child_singular && (odd || child_is_zero);
                if magnitude_zero {
                    self.zero_magnitudes.insert(expression);
                    singular |= *exponent < 0;
                }
                *exponent > 0 && magnitude_zero
            }
        };

        let scalar_scope = matches!(
            occurrence.kind,
            ExpressionKind::Product | ExpressionKind::Sum | ExpressionKind::Neg
        ) || matches!(occurrence.kind, ExpressionKind::Function(function) if function.is_linear());
        if scalar_scope && !zero && !singular {
            zero = self.scope_has_odd_class_stabilizer(projection, expression)?;
        }
        if zero {
            self.zero_expressions.insert(expression);
        }
        if singular {
            self.singular_expressions.insert(expression);
        }
        Ok(zero)
    }

    fn tensor_descendants<Aind: AbsInd>(
        projection: &CanonicalProjection<Aind>,
        roots: &[usize],
    ) -> BTreeSet<usize> {
        fn visit<Aind: AbsInd>(
            projection: &CanonicalProjection<Aind>,
            expression: usize,
            tensors: &mut BTreeSet<usize>,
            seen: &mut BTreeSet<usize>,
        ) {
            if !seen.insert(expression) {
                return;
            }
            let occurrence = &projection.expressions[expression];
            match &occurrence.kind {
                ExpressionKind::Tensor(tensor) => {
                    tensors.insert(*tensor);
                }
                ExpressionKind::Power { copies, .. } => {
                    for &copy in copies {
                        visit(projection, copy, tensors, seen);
                    }
                }
                _ => {
                    for &child in &occurrence.children {
                        visit(projection, child, tensors, seen);
                    }
                }
            }
        }
        let mut tensors = BTreeSet::new();
        let mut seen = BTreeSet::new();
        for &root in roots {
            visit(projection, root, &mut tensors, &mut seen);
        }
        tensors
    }

    fn lines_for_tensors<Aind: AbsInd>(
        projection: &CanonicalProjection<Aind>,
        tensors: &BTreeSet<usize>,
    ) -> BTreeSet<usize> {
        tensors
            .iter()
            .flat_map(|tensor| {
                projection.tensors[*tensor]
                    .ports
                    .iter()
                    .map(|port| port.line)
            })
            .collect()
    }

    fn port_in_group(role: &IncidenceRole, key: GroupKey) -> bool {
        matches!(role, IncidenceRole::Group { key: candidate, .. } if *candidate == key)
    }

    fn layout_path(key: GroupKey) -> Vec<usize> {
        match key {
            GroupKey::Intrinsic => vec![0],
            GroupKey::Argument(argument) => vec![1, argument],
        }
    }
}

struct Rebuilder<'a, Aind: AbsInd, F> {
    projection: &'a CanonicalProjection<Aind>,
    analysis: SignedAnalysis,
    decoration: NormalizedDecoration,
    assigned_lines: Vec<Option<LibrarySlot<Aind>>>,
    counters: BTreeMap<Representation<LibraryRep>, usize>,
    dummy_allocator: &'a mut DummyAllocator<Aind>,
    new_dummy: &'a mut F,
    preserve_external: bool,
    realized_sinks: BTreeSet<SinkKey>,
    residual_power_phase: bool,
    payload_edited: bool,
}

impl<'a, Aind, F> Rebuilder<'a, Aind, F>
where
    Aind: AbsInd + DummyAind + ParseableAind,
    F: FnMut(usize) -> Aind,
{
    fn new(
        projection: &'a CanonicalProjection<Aind>,
        analysis: SignedAnalysis,
        dummy_allocator: &'a mut DummyAllocator<Aind>,
        new_dummy: &'a mut F,
    ) -> Result<Self, CanonicalizationError> {
        Self::with_scope(
            projection,
            analysis,
            projection.root_expression,
            dummy_allocator,
            new_dummy,
            true,
        )
    }

    fn with_scope(
        projection: &'a CanonicalProjection<Aind>,
        analysis: SignedAnalysis,
        scope: usize,
        dummy_allocator: &'a mut DummyAllocator<Aind>,
        new_dummy: &'a mut F,
        preserve_external: bool,
    ) -> Result<Self, CanonicalizationError> {
        let decoration = analysis.normalized_decoration(projection, scope, &analysis.phases)?;
        Ok(Self {
            projection,
            analysis,
            decoration,
            assigned_lines: vec![None; projection.lines.len()],
            counters: BTreeMap::new(),
            dummy_allocator,
            new_dummy,
            preserve_external,
            realized_sinks: BTreeSet::new(),
            residual_power_phase: false,
            payload_edited: false,
        })
    }

    /// Verify the visible equality partitions under Graphica's selected map.
    ///
    /// This is deliberately a conservative, one-sided certificate. It uses
    /// the already selected canonical vertex map and the exact dummy assignment
    /// made by this reconstruction; it never compares independently labeled
    /// graphs or invokes another graph search. Materialized Powers currently
    /// decline the optimization and are verified by the next whole-network
    /// iteration.
    fn certifies_selected_transport(&self) -> Result<(), RetryReason> {
        if !self.projection.powers.is_empty() {
            return Err(RetryReason::IncompleteStabilityCertificate);
        }
        let Some(assigned) = self
            .assigned_lines
            .iter()
            .copied()
            .collect::<Option<Vec<_>>>()
        else {
            return Err(RetryReason::IncompleteStabilityCertificate);
        };
        let mut indices = assigned
            .iter()
            .map(IsAbstractSlot::aind)
            .collect::<Vec<_>>();
        indices.sort_unstable();
        indices.dedup();
        let target_lines = assigned
            .iter()
            .map(|slot| indices.binary_search(&slot.aind()).unwrap())
            .collect::<Vec<_>>();

        // A merge or split of line identities can make previously distinct
        // leaves or operation children syntactically equal after execution.
        let mut reverse_lines = BTreeMap::new();
        for (source, &target) in target_lines.iter().enumerate() {
            if reverse_lines.insert(target, source).is_some() {
                return Err(RetryReason::VisibleClassMergeSplit);
            }
        }

        let Some(source) = self.visible_syntax_partition(None) else {
            return Err(RetryReason::IncompleteStabilityCertificate);
        };
        let Some(target) = self.visible_syntax_partition(Some(&target_lines)) else {
            return Err(RetryReason::IncompleteStabilityCertificate);
        };
        let mut forward = BTreeMap::new();
        let mut reverse = BTreeMap::new();
        for (source, target) in source.into_iter().zip(target) {
            if forward
                .insert(source, target)
                .is_some_and(|old| old != target)
                || reverse
                    .insert(target, source)
                    .is_some_and(|old| old != source)
            {
                return Err(RetryReason::VisibleClassMergeSplit);
            }
        }
        Ok(())
    }

    fn visible_syntax_partition(&self, target_lines: Option<&[usize]>) -> Option<Vec<usize>> {
        fn visit<Aind, F>(
            rebuilder: &Rebuilder<'_, Aind, F>,
            expression: usize,
            target_lines: Option<&[usize]>,
            assigned: &mut [Option<usize>],
            classes: &mut BTreeMap<VisibleSyntaxClass, usize>,
        ) -> Option<usize>
        where
            Aind: AbsInd + DummyAind + ParseableAind,
            F: FnMut(usize) -> Aind,
        {
            if let Some(class) = assigned[expression] {
                return Some(class);
            }
            let occurrence = &rebuilder.projection.expressions[expression];
            let mut children = occurrence
                .children
                .iter()
                .map(|child| visit(rebuilder, *child, target_lines, assigned, classes))
                .collect::<Option<Vec<_>>>()?;
            let signature = match &occurrence.kind {
                ExpressionKind::Scalar(atom) => {
                    VisibleSyntaxClass::Scalar(semantic::SemanticAtomKey::new(atom.as_view()))
                }
                ExpressionKind::Tensor(tensor) => {
                    rebuilder.tensor_syntax_class(*tensor, target_lines)?
                }
                ExpressionKind::Product => {
                    children.sort_unstable();
                    VisibleSyntaxClass::Product(children)
                }
                ExpressionKind::Sum => {
                    children.sort_unstable();
                    VisibleSyntaxClass::Sum(children)
                }
                ExpressionKind::Neg => VisibleSyntaxClass::Neg(children[0]),
                ExpressionKind::Function(function) => {
                    VisibleSyntaxClass::Function((*function).into(), children[0])
                }
                ExpressionKind::Power { .. } => return None,
            };
            let class = if let Some(class) = classes.get(&signature) {
                *class
            } else {
                let class = classes.len();
                classes.insert(signature, class);
                class
            };
            assigned[expression] = Some(class);
            Some(class)
        }

        let mut assigned = vec![None; self.projection.expressions.len()];
        let mut classes = BTreeMap::new();
        visit(
            self,
            self.projection.root_expression,
            target_lines,
            &mut assigned,
            &mut classes,
        )?;
        assigned.into_iter().collect()
    }

    fn tensor_syntax_class(
        &self,
        tensor: usize,
        target_lines: Option<&[usize]>,
    ) -> Option<VisibleSyntaxClass> {
        let occurrence = &self.projection.tensors[tensor];
        let mut slots = vec![None; occurrence.layout.slot_count];
        if let Some(target_lines) = target_lines {
            let mut groups = BTreeMap::<(GroupKey, SymmetryKind), Vec<_>>::new();
            for port in &occurrence.ports {
                let value = (
                    semantic::representation_key(port.slot.rep()),
                    target_lines[port.line],
                );
                match port.role {
                    IncidenceRole::Group { key, kind } => {
                        groups.entry((key, kind)).or_default().push((port, value));
                    }
                    _ => slots[port.flat_slot] = Some(value),
                }
            }
            for ((key, kind), members) in groups {
                let ordered = self.canonical_members(key, kind, members).ok()?;
                let holes = match key {
                    GroupKey::Intrinsic => occurrence
                        .layout
                        .arguments
                        .iter()
                        .filter_map(|argument| match argument {
                            LayoutArgument::DirectSlot(hole) => Some(*hole),
                            _ => None,
                        })
                        .collect::<Vec<_>>(),
                    GroupKey::Argument(argument) => match &occurrence.layout.arguments[argument] {
                        LayoutArgument::Group { holes, .. } => holes.clone(),
                        _ => return None,
                    },
                };
                if holes.len() != ordered.len() {
                    return None;
                }
                for (hole, (_, value)) in holes.into_iter().zip(ordered) {
                    slots[hole] = Some(value);
                }
            }
        } else {
            for port in &occurrence.ports {
                slots[port.flat_slot] =
                    Some((semantic::representation_key(port.slot.rep()), port.line));
            }
        }
        Some(VisibleSyntaxClass::Tensor {
            color: occurrence.layout.color.clone(),
            slots: slots.into_iter().collect::<Option<Vec<_>>>()?,
        })
    }

    fn expression(
        &mut self,
        expression: usize,
    ) -> Result<SymbolicNet<Aind>, CanonicalizationError> {
        if self.analysis.zero_expressions.contains(&expression) {
            return Ok(Network::zero());
        }
        let signed = self.decoration.signs_expression(expression);
        if signed {
            self.realized_sinks
                .extend(
                    self.decoration
                        .sinks
                        .iter()
                        .filter_map(|(key, target)| match target {
                            SinkTarget::Leaf {
                                expression: candidate,
                            }
                            | SinkTarget::Result {
                                expression: candidate,
                            } if *candidate == expression => Some(key.clone()),
                            _ => None,
                        }),
                );
        }
        let occurrence = &self.projection.expressions[expression];
        let mut value = match &occurrence.kind {
            ExpressionKind::Scalar(atom) => {
                let atom = if signed { -atom.clone() } else { atom.clone() };
                return Ok(Network::from_scalar(atom));
            }
            ExpressionKind::Tensor(tensor) => return self.tensor(*tensor, signed),
            ExpressionKind::Product => {
                let children = self.sorted_children(&occurrence.children);
                let mut factors = Vec::with_capacity(children.len());
                for child in children {
                    factors.push(self.expression(child)?);
                }
                Self::product(factors)
            }
            ExpressionKind::Sum => self.sum(&occurrence.children)?,
            ExpressionKind::Neg => -self.expression(occurrence.children[0])?,
            ExpressionKind::Function(function) => {
                let child = occurrence.children[0];
                let value = if self.analysis.zero_expressions.contains(&child) {
                    Network::from_scalar(Atom::Zero)
                } else {
                    self.expression(child)?
                };
                value.fun(*function)
            }
            ExpressionKind::Power {
                exponent, copies, ..
            } => self.power(expression, *exponent, copies)?,
        };
        if signed {
            value = -value;
        }
        Ok(value)
    }

    fn tensor(
        &mut self,
        tensor: usize,
        signed: bool,
    ) -> Result<SymbolicNet<Aind>, CanonicalizationError> {
        let occurrence = &self.projection.tensors[tensor];
        let mut slots = vec![None; occurrence.layout.slot_count];
        let mut groups = BTreeMap::<(GroupKey, SymmetryKind), Vec<_>>::new();
        let mut ports = occurrence.ports.iter().collect::<Vec<_>>();
        ports.sort_unstable_by_key(|port| {
            self.projection.vertex_map[self.projection.lines[port.line].vertex]
        });
        for port in ports {
            let slot = self.slot_for(port.line, port.slot.rep())?;
            match port.role {
                IncidenceRole::Group { key, kind } => {
                    groups.entry((key, kind)).or_default().push((port, slot));
                }
                _ => slots[port.flat_slot] = Some(slot),
            }
        }
        for ((key, kind), members) in groups {
            let ordered = self.canonical_members(key, kind, members)?;
            let holes = match key {
                GroupKey::Intrinsic => occurrence
                    .layout
                    .arguments
                    .iter()
                    .filter_map(|argument| match argument {
                        LayoutArgument::DirectSlot(hole) => Some(*hole),
                        _ => None,
                    })
                    .collect::<Vec<_>>(),
                GroupKey::Argument(argument) => match &occurrence.layout.arguments[argument] {
                    LayoutArgument::Group { holes, .. } => holes.clone(),
                    _ => {
                        return Err(CanonicalizationError::StructureMismatch {
                            expression: occurrence.tensor.expression.clone(),
                        });
                    }
                },
            };
            if holes.len() != ordered.len() {
                return Err(CanonicalizationError::StructureMismatch {
                    expression: occurrence.tensor.expression.clone(),
                });
            }
            for (hole, (_, slot)) in holes.into_iter().zip(ordered) {
                slots[hole] = Some(slot);
            }
        }
        let slots = slots
            .into_iter()
            .collect::<Option<Vec<_>>>()
            .ok_or_else(|| CanonicalizationError::StructureMismatch {
                expression: occurrence.tensor.expression.clone(),
            })?;
        let zero_groups = self
            .analysis
            .zero_groups
            .get(&tensor)
            .cloned()
            .unwrap_or_default();
        let mut negative_groups = BTreeSet::new();
        for metadata in &self.analysis.sites {
            debug_assert_eq!(
                metadata.expression,
                self.projection.tensors[metadata.occurrence].expression
            );
            if metadata.occurrence != tensor
                || metadata.lifts
                || zero_groups.contains(&metadata.key)
                || !self.decoration.signs_nested(tensor, metadata.key)
            {
                continue;
            }
            debug_assert!(!metadata.intrinsic);
            self.realized_sinks
                .extend(
                    self.decoration
                        .sinks
                        .iter()
                        .filter_map(|(key, target)| match target {
                            SinkTarget::Nested {
                                occurrence,
                                key: group,
                            } if *occurrence == tensor && *group == metadata.key => {
                                Some(key.clone())
                            }
                            _ => None,
                        }),
                );
            negative_groups.insert(metadata.key);
        }
        let (mut expression, active_slots) =
            occurrence
                .layout
                .rebuild(&slots, &zero_groups, &negative_groups);
        self.payload_edited |= Self::rebuilt_payload_changed(occurrence, expression.as_view());
        if signed {
            expression = -expression;
        }
        if active_slots.is_empty() {
            self.payload_edited = true;
            Ok(Network::from_scalar(expression))
        } else {
            let tensor = SymbolicTensor {
                structure: OrderedStructure::new(active_slots).structure,
                is_metric: occurrence.tensor.is_metric,
                is_composite: occurrence.tensor.is_composite,
                expression,
            };
            Ok(Network::from_tensor(tensor))
        }
    }

    fn rebuilt_payload_changed(
        occurrence: &super::projection::TensorOccurrence<Aind>,
        expression: AtomView<'_>,
    ) -> bool {
        let AtomView::Fun(function) = expression else {
            return true;
        };
        if function.get_symbol() != occurrence.layout.head
            || function.get_nargs() != occurrence.layout.arguments.len()
        {
            return true;
        }
        function.iter().zip(&occurrence.layout.arguments).any(
            |(value, argument)| match argument {
                LayoutArgument::Opaque(expected) => value != expected.as_view(),
                LayoutArgument::DirectSlot(_) => false,
                LayoutArgument::SlotBundle { holes } => {
                    !matches!(value, AtomView::Fun(bundle)
                        if bundle.get_symbol() == spenso::structure::abstract_index::AIND_SYMBOLS.aind
                            && bundle.get_nargs() == holes.len())
                }
                LayoutArgument::Group {
                    symbol,
                    kind,
                    prefix,
                    holes,
                } => !matches!(
                    super::TensorLayout::partial_group::<Aind>(value, &occurrence.tensor.expression),
                    Ok(Some((candidate_symbol, candidate_kind, candidate_prefix, members)))
                        if candidate_symbol == *symbol
                            && candidate_kind == *kind
                            && candidate_prefix == *prefix
                            && members.len() == holes.len()
                ),
            },
        )
    }

    fn canonical_members<'b, V>(
        &self,
        key: GroupKey,
        kind: SymmetryKind,
        mut members: Vec<(&'b super::projection::TensorPort<Aind>, V)>,
    ) -> Result<Vec<(&'b super::projection::TensorPort<Aind>, V)>, CanonicalizationError> {
        members.sort_by_key(|(port, _)| self.projection.vertex_map[port.vertex]);
        if kind != SymmetryKind::Cyclic || members.len() < 2 {
            return Ok(members);
        }
        let mut by_vertex = members
            .into_iter()
            .map(|member| (self.projection.vertex_map[member.0.vertex], member))
            .collect::<BTreeMap<_, _>>();
        let mut vertex = *by_vertex.keys().next().unwrap();
        let mut ordered = Vec::with_capacity(by_vertex.len());
        for _ in 0..by_vertex.len() {
            ordered.push(by_vertex.remove(&vertex).ok_or_else(|| {
                CanonicalizationError::Projection(
                    "canonical cyclic frame visits one port more than once".into(),
                )
            })?);
            vertex = self
                .projection
                .graph
                .node(vertex)
                .edges
                .iter()
                .map(|edge| self.projection.graph.edge(*edge))
                .find_map(|edge| {
                    (edge.vertices.0 == vertex
                        && matches!(edge.data.data, super::projection::UnifiedEdgeColor::Cyclic(candidate) if candidate == key))
                    .then_some(edge.vertices.1)
                })
                .ok_or_else(|| {
                    CanonicalizationError::Projection(
                        "canonical cyclic frame is not a complete directed cycle".into(),
                    )
                })?;
        }
        Ok(ordered)
    }

    fn sum(&mut self, children: &[usize]) -> Result<SymbolicNet<Aind>, CanonicalizationError> {
        let children = self.sorted_children(children);
        let child_lines = children
            .iter()
            .map(|child| {
                let tensors = SignedAnalysis::tensor_descendants(self.projection, &[*child]);
                SignedAnalysis::lines_for_tensors(self.projection, &tensors)
            })
            .collect::<Vec<_>>();
        let mut counts = BTreeMap::<usize, usize>::new();
        for lines in &child_lines {
            for line in lines {
                *counts.entry(*line).or_insert(0) += 1;
            }
        }
        let mut shared_lines = counts
            .into_iter()
            .filter_map(|(line, count)| (count > 1).then_some(line))
            .collect::<Vec<_>>();
        // Source line ids are hidden input order; shared interfaces must
        // consume dummy ordinals in canonical line order before branches fork.
        shared_lines.sort_unstable_by_key(|line| {
            self.projection.vertex_map[self.projection.lines[*line].vertex]
        });
        for line in shared_lines {
            let representation = self.projection.lines[line].group;
            let _ = self.slot_for(line, representation)?;
        }

        let start = self.counters.clone();
        let mut end = start.clone();
        let mut terms = Vec::new();
        for child in children {
            self.counters.clone_from(&start);
            if !self.analysis.zero_expressions.contains(&child) {
                terms.push(self.expression(child)?);
            }
            for (group, count) in &self.counters {
                let maximum = end.entry(*group).or_insert(0);
                *maximum = (*maximum).max(*count);
            }
        }
        self.counters = end;
        Ok(Self::sum_networks(terms))
    }

    fn power(
        &mut self,
        expression: usize,
        exponent: i8,
        copies: &[usize],
    ) -> Result<SymbolicNet<Aind>, CanonicalizationError> {
        let descriptor = self
            .projection
            .power_descriptor(expression)
            .ok_or_else(|| {
                CanonicalizationError::Projection(format!(
                    "materialized Power {expression} has no validated copy descriptor"
                ))
            })?
            .clone();
        let described_copies = descriptor
            .copies
            .iter()
            .map(|copy| copy.expression)
            .collect::<BTreeSet<_>>();
        if descriptor.exponent != exponent || described_copies != copies.iter().copied().collect() {
            return Err(CanonicalizationError::Projection(format!(
                "materialized Power {expression} changed exponent or complete-copy membership"
            )));
        }
        if self.analysis.zero_magnitudes.contains(&expression) {
            if exponent > 0 {
                return Ok(Network::zero());
            }
            return Ok(Network::from_scalar(Atom::Zero).pow(exponent));
        }
        let Some(_) = descriptor.copies.first() else {
            return Ok(Network::one());
        };

        let mut classified = BTreeMap::new();
        for copy in &descriptor.copies {
            classified.insert(copy.root, self.classify_power_copy(expression, copy)?);
        }
        let common = classified
            .values()
            .next()
            .expect("a materialized Power has at least one copy")
            .class
            .clone();
        if classified.values().any(|copy| copy.class != common) {
            return Err(CanonicalizationError::PowerReconstruction {
                expression,
                reason: "complete copies do not represent one common base up to scalar sign".into(),
            });
        }

        let interface_copies = descriptor
            .copies
            .iter()
            .filter(|copy| {
                copy.boundaries.iter().any(|boundary| {
                    matches!(boundary.target, PowerBoundaryTarget::Interface { .. })
                })
            })
            .collect::<Vec<_>>();
        if interface_copies.len() > 1 {
            return Err(CanonicalizationError::PowerReconstruction {
                expression,
                reason: "more than one complete copy owns the Power result interface".into(),
            });
        }
        let chosen = interface_copies.into_iter().next().unwrap_or_else(|| {
            descriptor
                .copies
                .iter()
                .min_by_key(|copy| copy.root)
                .expect("a materialized Power has at least one copy")
        });
        let chosen_phase = classified[&chosen.root].negative;
        let all_phase = classified
            .values()
            .fold(false, |phase, copy| phase ^ copy.negative);
        let represented_phase = chosen_phase && exponent.unsigned_abs() % 2 == 1;
        let residual = all_phase ^ represented_phase;
        let mut power = self.expression(chosen.expression)?.pow(exponent);
        if residual {
            self.residual_power_phase = true;
            let sink = SinkKey {
                kind: 3,
                vertex: descriptor.result,
                path: Vec::new(),
            };
            if !self.realized_sinks.insert(sink.clone()) {
                self.realized_sinks.remove(&sink);
            }
            power = -power;
        }
        Ok(power)
    }

    fn classify_power_copy(
        &self,
        power: usize,
        copy: &PowerCopyDescriptor,
    ) -> Result<ClassifiedPowerCopy, CanonicalizationError> {
        let mut dummy_allocator = DummyAllocator::new();
        let mut new_dummy = Aind::new_dummy_at;
        let mut rebuilder = Rebuilder::with_scope(
            self.projection,
            self.analysis.clone(),
            copy.expression,
            &mut dummy_allocator,
            &mut new_dummy,
            false,
        )?;
        let network = rebuilder.expression(copy.expression)?;
        let state = match network.state {
            NetworkState::PureScalar => 0,
            NetworkState::Scalar => 1,
            NetworkState::SelfDualTensor => 2,
            NetworkState::Tensor => 3,
        };
        #[cfg(test)]
        super::driver::record_temporary_scope_execution();
        let atom = execute_atom(network).map_err(|error| match error {
            CanonicalizationError::Execution { error, .. } => CanonicalizationError::Execution {
                scope: format!("Power {power} complete copy {}", copy.root),
                error,
            },
            error => error,
        })?;
        let direct = semantic::SemanticAtomKey::new(atom.as_view());
        let negative = semantic::SemanticAtomKey::new((-atom).as_view());
        let (body, negative) = if negative < direct {
            (negative, true)
        } else {
            (direct, false)
        };

        let mut seam = Vec::with_capacity(copy.boundaries.len());
        for boundary in &copy.boundaries {
            let mut ports = Vec::with_capacity(boundary.ports.len());
            for &canonical_port in &boundary.ports {
                let port = self
                    .projection
                    .tensors
                    .iter()
                    .flat_map(|tensor| &tensor.ports)
                    .find(|port| self.projection.vertex_map[port.vertex] == canonical_port)
                    .ok_or_else(|| CanonicalizationError::PowerReconstruction {
                        expression: power,
                        reason: format!(
                            "copy {} boundary references missing port {canonical_port}",
                            copy.root
                        ),
                    })?;
                ports.push((
                    port.role.clone(),
                    semantic::representation_key(port.slot.rep()),
                ));
            }
            ports.sort();
            seam.push(ports);
        }
        seam.sort();
        Ok(ClassifiedPowerCopy {
            class: PowerCopyClass { state, seam, body },
            negative,
        })
    }

    fn slot_for(
        &mut self,
        line: usize,
        representation: Representation<LibraryRep>,
    ) -> Result<LibrarySlot<Aind>, CanonicalizationError> {
        if let Some(slot) = self.assigned_lines[line] {
            return Ok(representation.slot(slot.aind()));
        }
        let line_data = &self.projection.lines[line];
        if self.preserve_external
            && let Some(external) = line_data.external
        {
            self.assigned_lines[line] = Some(external);
            return Ok(representation.slot(external.aind()));
        }
        let ordinal = self.counters.entry(line_data.group).or_insert(0);
        let index = self
            .dummy_allocator
            .get(line_data.group, *ordinal, self.new_dummy);
        *ordinal += 1;
        let slot = line_data.group.slot(index);
        self.assigned_lines[line] = Some(slot);
        Ok(representation.slot(slot.aind()))
    }

    fn sorted_children(&self, children: &[usize]) -> Vec<usize> {
        let mut children = children.to_vec();
        children.sort_unstable_by_key(|child| {
            self.projection.vertex_map[self.projection.expressions[*child].root]
        });
        children
    }

    fn product(factors: Vec<SymbolicNet<Aind>>) -> SymbolicNet<Aind> {
        let mut factors = factors.into_iter();
        let Some(first) = factors.next() else {
            return Network::one();
        };
        first.n_mul(factors)
    }

    fn sum_networks(terms: Vec<SymbolicNet<Aind>>) -> SymbolicNet<Aind> {
        let mut terms = terms.into_iter();
        let Some(first) = terms.next() else {
            return Network::zero();
        };
        first.n_add(terms)
    }
}

#[cfg(test)]
mod tests {
    use linnet::permutation::Permutation;
    use spenso::{
        antisym,
        network::graph::{NetworkLeaf, NetworkNode},
        slot,
        structure::{abstract_index::AbstractIndex, slot::IsAbstractSlot},
        tensor_symbol,
    };
    use symbolica::{
        atom::{AtomCore, FunctionBuilder},
        function, symbol,
    };

    use super::*;
    use crate::{tensor::canonicalize::projection, test_support::test_initialize};

    #[test]
    fn class_stabilizer_retains_a_composed_odd_word() {
        let analysis = SignedAnalysis {
            group: SignedGroup::new(0, 0, Vec::new()).unwrap(),
            sites: Vec::new(),
            phases: Vec::new(),
            zero_groups: BTreeMap::new(),
            zero_expressions: BTreeSet::new(),
            zero_magnitudes: BTreeSet::new(),
            singular_expressions: BTreeSet::new(),
        };
        let rotation =
            SignedAction::from_parts(Permutation::from_map(vec![1, 2, 3, 0]), Vec::new()).unwrap();
        let group = SignedGroup::new(4, 0, vec![rotation]).unwrap();
        let first = ScopeClassKey {
            root: 0,
            kind: UnifiedNodeColor::Root,
            boundary: Vec::new(),
            value: semantic::SemanticAtomKey::new(Atom::num(0).as_view()),
        };
        let second = ScopeClassKey {
            value: semantic::SemanticAtomKey::new(Atom::num(1).as_view()),
            ..first.clone()
        };

        let odd = analysis
            .group_has_odd_class_stabilizer(&group, |action| {
                Ok(match action.vertex_map() {
                    [0, 1, 2, 3] => (first.clone(), false),
                    [1, 2, 3, 0] => (second.clone(), false),
                    [2, 3, 0, 1] => (first.clone(), true),
                    [3, 0, 1, 2] => (second.clone(), true),
                    map => panic!("unexpected cyclic action {map:?}"),
                })
            })
            .unwrap();

        assert!(odd);
    }

    #[test]
    fn power_copy_sign_edits_are_copy_on_write() {
        let rep = test_initialize().mink4;
        let first = slot!(rep, AbstractIndex::Dummy(0));
        let second = slot!(rep, AbstractIndex::Dummy(1));
        let antisymmetric = tensor_symbol!(reconstruct_copy_on_write_tensor; Antisymmetric);
        let expression =
            function!(antisymmetric, first.to_atom(), second.to_atom()).pow(Atom::num(3));
        let policy = super::super::CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();
        assert_eq!(policy.network().store.tensors.len(), 1);
        let source_payload = policy.network().store.tensors[0].expression.clone();

        let projection = projection::project(&policy, projection::DEFAULT_GRAPH_BUDGET).unwrap();
        assert_eq!(projection.tensors.len(), 3);
        assert!(
            projection
                .tensors
                .iter()
                .all(|occurrence| occurrence.tensor.expression == source_payload)
        );
        let mut prepared = prepare_reconstruction(&projection).unwrap();
        assert_eq!(prepared.analysis.phases.len(), 3);
        // Isolate copy-on-write from which orientation Graphica selects for
        // the otherwise valid three-copy Power projection.
        prepared.analysis.phases.fill(true);
        let ExpressionKind::Power { copies, .. } =
            &projection.expressions[projection.root_expression].kind
        else {
            panic!("the test fixture projects to one Power");
        };
        let mut allocator = DummyAllocator::new();
        let mut new_dummy = AbstractIndex::Dummy;
        let mut signed = Rebuilder::with_scope(
            &projection,
            prepared.analysis.clone(),
            copies[0],
            &mut allocator,
            &mut new_dummy,
            false,
        )
        .unwrap();
        let signed = signed.expression(copies[0]).unwrap();

        assert_eq!(policy.network().store.tensors[0].expression, source_payload);
        assert!(
            projection
                .tensors
                .iter()
                .all(|occurrence| occurrence.tensor.expression == source_payload)
        );
        assert_eq!(signed.store.tensors.len(), 1);

        prepared.analysis.phases.fill(false);
        let mut sibling_allocator = DummyAllocator::new();
        let mut sibling_dummy = AbstractIndex::Dummy;
        let mut sibling = Rebuilder::with_scope(
            &projection,
            prepared.analysis,
            copies[1],
            &mut sibling_allocator,
            &mut sibling_dummy,
            false,
        )
        .unwrap();
        let sibling = sibling.expression(copies[1]).unwrap();
        assert_eq!(
            signed.store.tensors[0].expression,
            -sibling.store.tensors[0].expression.clone()
        );
    }

    #[test]
    fn power_realization_records_only_sinks_used_after_folding() {
        let rep = test_initialize().mink4;
        let first = slot!(rep, AbstractIndex::Dummy(0));
        let second = slot!(rep, AbstractIndex::Dummy(1));
        let antisymmetric = tensor_symbol!(reconstruct_power_realization_tensor; Antisymmetric);

        for exponent in [2, -2] {
            let expression = function!(antisymmetric, first.to_atom(), second.to_atom())
                .pow(Atom::num(exponent));
            let policy =
                super::super::CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();
            let projection =
                projection::project(&policy, projection::DEFAULT_GRAPH_BUDGET).unwrap();
            let (analysis, _) = SignedAnalysis::new(&projection).unwrap();
            let power = projection
                .power_descriptor(projection.root_expression)
                .unwrap();
            let chosen = power
                .copies
                .iter()
                .find(|copy| {
                    copy.boundaries.iter().any(|boundary| {
                        matches!(boundary.target, PowerBoundaryTarget::Interface { .. })
                    })
                })
                .unwrap_or_else(|| power.copies.iter().min_by_key(|copy| copy.root).unwrap());
            let unchosen = power
                .copies
                .iter()
                .find(|copy| copy.root != chosen.root)
                .unwrap();
            let unchosen_site = analysis
                .sites
                .iter()
                .position(|site| site.expression == unchosen.expression)
                .unwrap();
            let descriptor =
                SignedAnalysis::scope_descriptor(&projection, projection.root_expression);

            let mut one_phase = vec![false; analysis.phases.len()];
            one_phase[unchosen_site] = true;
            let classified = analysis
                .classify_scope(
                    &projection,
                    projection.root_expression,
                    &descriptor,
                    &one_phase,
                )
                .unwrap();
            assert_eq!(
                classified.realization,
                vec![SinkKey {
                    kind: 3,
                    vertex: power.result,
                    path: Vec::new(),
                }]
            );

            let two_phases = vec![true; analysis.phases.len()];
            let decoration = analysis
                .normalized_decoration(&projection, projection.root_expression, &two_phases)
                .unwrap();
            let chosen_sink = decoration
                .sinks
                .iter()
                .find_map(|(sink, target)| match target {
                    SinkTarget::Leaf { expression } | SinkTarget::Result { expression }
                        if *expression == chosen.expression =>
                    {
                        Some(sink.clone())
                    }
                    SinkTarget::Nested { occurrence, .. }
                        if projection.tensors[*occurrence].expression == chosen.expression =>
                    {
                        Some(sink.clone())
                    }
                    _ => None,
                })
                .unwrap();
            let classified = analysis
                .classify_scope(
                    &projection,
                    projection.root_expression,
                    &descriptor,
                    &two_phases,
                )
                .unwrap();
            assert_eq!(classified.realization, vec![chosen_sink]);
        }
    }

    #[test]
    fn residual_power_copy_phases_stay_outside_even_powers() {
        let rep = test_initialize().mink4;
        let first = slot!(rep, AbstractIndex::Dummy(0));
        let second = slot!(rep, AbstractIndex::Dummy(1));
        let antisymmetric = tensor_symbol!(reconstruct_residual_power_tensor; Antisymmetric);

        for exponent in [2, -2] {
            let expression = function!(antisymmetric, first.to_atom(), second.to_atom())
                .pow(Atom::num(exponent));
            let policy =
                super::super::CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();
            let projection =
                projection::project(&policy, projection::DEFAULT_GRAPH_BUDGET).unwrap();
            let phase_count = prepare_reconstruction(&projection)
                .unwrap()
                .analysis
                .phases
                .len();
            assert_eq!(phase_count, 2);

            let rebuild = |phases: Vec<bool>| {
                let mut prepared = prepare_reconstruction(&projection).unwrap();
                prepared.analysis.phases = phases;
                let mut allocator = DummyAllocator::new();
                let reconstructed = reconstruct(
                    &projection,
                    prepared,
                    &mut allocator,
                    &mut AbstractIndex::Dummy,
                )
                .unwrap();
                let atom = execute_atom(reconstructed.network.clone()).unwrap();
                (reconstructed.retry_reason, atom)
            };

            let (_, unsigned) = rebuild(vec![false; phase_count]);
            let mut one_phase = vec![false; phase_count];
            one_phase[0] = true;
            let (reason, one_phase) = rebuild(one_phase);
            assert_eq!(reason, Some(RetryReason::ResidualPowerPhase));
            assert_eq!(one_phase, -unsigned.clone());

            let (_, two_phases) = rebuild(vec![true; phase_count]);
            assert_eq!(two_phases, unsigned);
        }
    }

    #[test]
    fn dummy_allocator_reindexes_surviving_groups_between_iterations() {
        let reps = test_initialize();
        let original = slot!(reps.bis4, AbstractIndex::Dummy(0));
        let left = tensor_symbol!(reconstruct_retry_surviving_left);
        let right = tensor_symbol!(reconstruct_retry_surviving_right);
        let expression = function!(left, original.to_atom()) * function!(right, original.to_atom());
        let policy = super::super::CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();
        let projection = projection::project(&policy, projection::DEFAULT_GRAPH_BUDGET).unwrap();
        let prepared = prepare_reconstruction(&projection).unwrap();
        let mut allocator = DummyAllocator::new();
        let mut allocated_positions = Vec::new();
        let mut new_dummy = |position| {
            allocated_positions.push(position);
            AbstractIndex::Dummy(position + 10)
        };

        // Simulate a previous iteration in which a Minkowski line preceded
        // this Bispinor line. Its callback values remain request-local, while
        // the next reconstruction must derive positions from surviving lines.
        assert_eq!(
            allocator.get(reps.mink4.to_lib().base(), 0, &mut new_dummy),
            AbstractIndex::Dummy(10)
        );
        assert_eq!(
            allocator.get(reps.bis4.to_lib().base(), 0, &mut new_dummy),
            AbstractIndex::Dummy(11)
        );
        let reconstructed =
            reconstruct(&projection, prepared, &mut allocator, &mut new_dummy).unwrap();
        let reindexed = slot!(reps.bis4, AbstractIndex::Dummy(10));
        let expected = function!(left, reindexed.to_atom()) * function!(right, reindexed.to_atom());

        assert_eq!(allocated_positions, vec![0, 1]);
        assert_eq!(execute_atom(reconstructed.network).unwrap(), expected);
    }

    #[test]
    fn zero_partial_group_does_not_record_an_inert_nested_sign() {
        let rep = test_initialize().mink4;
        let first = slot!(rep, AbstractIndex::Dummy(0));
        let second = slot!(rep, AbstractIndex::Dummy(1));
        let outer = tensor_symbol!(reconstruct_zero_partial_group_outer);
        let expression = FunctionBuilder::new(outer)
            .add_arg(antisym!(second, first))
            .finish();
        let policy = super::super::CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();
        let projection = projection::project(&policy, projection::DEFAULT_GRAPH_BUDGET).unwrap();
        let (mut analysis, _) = SignedAnalysis::new(&projection).unwrap();
        let site = analysis
            .sites
            .iter()
            .position(|site| !site.lifts && site.key == GroupKey::Argument(0))
            .unwrap();
        let metadata = analysis.sites[site].clone();
        analysis.phases.fill(false);
        analysis.phases[site] = true;
        analysis
            .zero_groups
            .entry(metadata.occurrence)
            .or_default()
            .insert(metadata.key);
        analysis.canonicalize_decoration(&projection).unwrap();

        assert!(!analysis.phases[site]);

        let mut allocator = DummyAllocator::new();
        let mut new_dummy = AbstractIndex::Dummy;
        let mut rebuilder =
            Rebuilder::new(&projection, analysis, &mut allocator, &mut new_dummy).unwrap();
        assert!(
            !rebuilder
                .decoration
                .signs_nested(metadata.occurrence, metadata.key)
        );
        rebuilder.expression(projection.root_expression).unwrap();

        assert!(rebuilder.realized_sinks.is_empty());
    }

    #[test]
    fn product_sign_sinks_are_deterministic_without_synthesizing_minus_one() {
        let rep = test_initialize().mink4;
        let first = slot!(rep, AbstractIndex::Dummy(0));
        let second = slot!(rep, AbstractIndex::Dummy(1));
        let antisymmetric = tensor_symbol!(reconstruct_sink_antisymmetric; Antisymmetric);
        let left = tensor_symbol!(reconstruct_sink_left);
        let right = tensor_symbol!(reconstruct_sink_right);
        let tensor_product = function!(antisymmetric, first.to_atom(), second.to_atom())
            * function!(left, first.to_atom())
            * function!(right, second.to_atom());
        let scalar_policy = super::super::CanonicalPolicyNet::<AbstractIndex>::parse(
            Atom::num(2) * tensor_product.clone(),
        )
        .unwrap();
        let scalar_projection =
            projection::project(&scalar_policy, projection::DEFAULT_GRAPH_BUDGET).unwrap();
        let scalar_prepared = prepare_reconstruction(&scalar_projection).unwrap();
        assert_eq!(scalar_prepared.analysis.phases, vec![true]);
        let scalar_expression = scalar_projection
            .expressions
            .iter()
            .enumerate()
            .filter(|(_, occurrence)| matches!(occurrence.kind, ExpressionKind::Scalar(_)))
            .min_by_key(|(_, occurrence)| scalar_projection.vertex_map[occurrence.root])
            .unwrap();
        let scalar_decoration = scalar_prepared
            .analysis
            .normalized_decoration(
                &scalar_projection,
                scalar_projection.root_expression,
                &scalar_prepared.analysis.phases,
            )
            .unwrap();
        let (scalar_key, scalar_target) = scalar_decoration.sinks.iter().next().unwrap();
        assert_eq!(scalar_decoration.sinks.len(), 1);
        assert_eq!(scalar_key.kind, 1);
        assert_eq!(
            scalar_key.vertex,
            scalar_projection.vertex_map[scalar_expression.1.root]
        );
        assert!(matches!(
            scalar_target,
            SinkTarget::Leaf { expression } if *expression == scalar_expression.0
        ));
        let mut scalar_allocator = DummyAllocator::new();
        let scalar_reconstructed = reconstruct(
            &scalar_projection,
            scalar_prepared,
            &mut scalar_allocator,
            &mut AbstractIndex::Dummy,
        )
        .unwrap();
        assert_eq!(
            scalar_reconstructed.network.store.scalar,
            vec![Atom::num(-2)]
        );
        assert_eq!(
            scalar_reconstructed
                .network
                .graph
                .graph
                .iter_nodes()
                .filter(|(_, _, node)| {
                    matches!(node, NetworkNode::Leaf(NetworkLeaf::Scalar(_)))
                })
                .count(),
            1
        );

        let tensor_policy =
            super::super::CanonicalPolicyNet::<AbstractIndex>::parse(tensor_product).unwrap();
        let tensor_projection =
            projection::project(&tensor_policy, projection::DEFAULT_GRAPH_BUDGET).unwrap();
        let tensor_prepared = prepare_reconstruction(&tensor_projection).unwrap();
        assert_eq!(tensor_prepared.analysis.phases, vec![true]);
        let tensor_expression = tensor_projection
            .expressions
            .iter()
            .enumerate()
            .filter(|(_, occurrence)| matches!(occurrence.kind, ExpressionKind::Tensor(_)))
            .min_by_key(|(_, occurrence)| tensor_projection.vertex_map[occurrence.root])
            .unwrap();
        let tensor_decoration = tensor_prepared
            .analysis
            .normalized_decoration(
                &tensor_projection,
                tensor_projection.root_expression,
                &tensor_prepared.analysis.phases,
            )
            .unwrap();
        let (tensor_key, tensor_target) = tensor_decoration.sinks.iter().next().unwrap();
        assert_eq!(tensor_decoration.sinks.len(), 1);
        assert_eq!(tensor_key.kind, 2);
        assert_eq!(
            tensor_key.vertex,
            tensor_projection.vertex_map[tensor_expression.1.root]
        );
        assert!(matches!(
            tensor_target,
            SinkTarget::Leaf { expression } if *expression == tensor_expression.0
        ));
        let mut tensor_allocator = DummyAllocator::new();
        let tensor_reconstructed = reconstruct(
            &tensor_projection,
            tensor_prepared,
            &mut tensor_allocator,
            &mut AbstractIndex::Dummy,
        )
        .unwrap();
        assert!(tensor_reconstructed.network.store.scalar.is_empty());
        assert!(
            tensor_reconstructed
                .network
                .graph
                .graph
                .iter_nodes()
                .all(|(_, _, node)| !matches!(node, NetworkNode::Leaf(NetworkLeaf::Scalar(_))))
        );
    }

    #[test]
    fn branch_local_transport_that_merges_visible_syntax_classes_declines_certification() {
        let rep = test_initialize().mink4;
        let first = slot!(rep, reconstruct_certificate_class_first);
        let second = slot!(rep, reconstruct_certificate_class_second);
        let tensor = tensor_symbol!(reconstruct_certificate_class_tensor);
        let left = tensor_symbol!(reconstruct_certificate_class_left);
        let right = tensor_symbol!(reconstruct_certificate_class_right);
        let nested_left = tensor_symbol!(reconstruct_certificate_nested_left);
        let nested_right = tensor_symbol!(reconstruct_certificate_nested_right);
        let wrapper = symbol!("reconstruct_certificate_class_wrapper");
        let leaf_classes = function!(tensor, first.to_atom(), first.to_atom())
            + function!(tensor, second.to_atom(), second.to_atom());
        let product_classes = function!(left, first.to_atom()) * function!(right, first.to_atom())
            + function!(left, second.to_atom()) * function!(right, second.to_atom());
        let first_sum = function!(left, first.to_atom()) * function!(right, first.to_atom())
            + function!(nested_left, first.to_atom()) * function!(nested_right, first.to_atom());
        let second_sum = function!(left, second.to_atom()) * function!(right, second.to_atom())
            + function!(nested_left, second.to_atom()) * function!(nested_right, second.to_atom());
        let sum_classes = function!(wrapper, first_sum) + function!(wrapper, second_sum);

        for expression in [leaf_classes, product_classes, sum_classes] {
            let policy =
                super::super::CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();
            let projection =
                projection::project(&policy, projection::DEFAULT_GRAPH_BUDGET).unwrap();
            let prepared = prepare_reconstruction(&projection).unwrap();
            let mut allocator = DummyAllocator::new();
            let reconstructed = reconstruct(
                &projection,
                prepared,
                &mut allocator,
                &mut AbstractIndex::Dummy,
            )
            .unwrap();

            assert_eq!(
                reconstructed.retry_reason,
                Some(RetryReason::VisibleClassMergeSplit)
            );
        }
    }
}
