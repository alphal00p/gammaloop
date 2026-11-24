use std::borrow::Cow;
use std::collections::{BTreeSet, HashSet};
use std::sync::LazyLock;

use color_eyre::Result;
use eyre::{WrapErr, ensure, eyre};
use itertools::Itertools;
use linnet::half_edge::involution::EdgeIndex;
use linnet::half_edge::subgraph::{Inclusion, InternalSubGraph, SuBitGraph, SubSetLike, SubSetOps};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    function, symbol,
};
use three_dimensional_reps::MediumMode;

use crate::uv::thermal_zero_mode::{ThermalZeroModeRegion, ZeroModeHostKind, ZeroModeTuple};
use crate::{
    debug_tags,
    graph::{Graph, LMBext, LmbError, LoopMomentumBasis},
    integrands::process::param_builder::ParamBuilderGraph,
    numerator::symbolica_ext::NumeratorAtomExt,
    utils::{GS, W_},
    uv::{
        ApproximationType, Integrands,
        approx::{ForestNodeLike, OrientationProjection, UVCtx},
        uv_graph::UVE,
    },
};

#[cfg(test)]
use crate::uv::approx::local_3d::Local3DLoopRescaling;

use super::{
    branches::{DirectResidueBranches, DirectResidueKey},
    coefficients::numerator_family_calls,
    forest::Direct3dApproximation,
};
use thermal_soft::{SoftBranchRole, SoftBranches};
use zero_mode::{ZeroModeClass, ZeroModeProjection};

/// One independently framed connected part of a direct forest sector.
/// Disconnected replay keeps these frames separate until an enclosing Taylor
/// operation contains and combines their carriers.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct DirectCoordinateFrame {
    pub(crate) active_subgraph: SuBitGraph,
    pub(crate) lmb: LoopMomentumBasis,
}

/// Choose the loop coordinates in which one direct sector enters its next
/// Taylor operation. An enclosing operation must retain every loop carrier
/// already owned by the part of the sector contained in `given`;
/// otherwise an equivalent basis can turn a previously hard child momentum
/// into a crown-shifted momentum before the outer series is taken.
/// The result contains the canonical component route followed by the selected
/// operation route, which may differ when contraction requires a new quotient
/// carrier or when genuine prior carriers must be preserved.
pub(super) fn coordinate_lmb<S: ForestNodeLike>(
    ctx: &UVCtx<'_>,
    current: &S,
    given: &S,
    prior_active_subgraph: Option<&SuBitGraph>,
    prior_frames: &[DirectCoordinateFrame],
    active_subgraph: &SuBitGraph,
) -> Result<(LoopMomentumBasis, LoopMomentumBasis)> {
    let graph = ctx.graph;
    let active_subgraph = active_subgraph.intersection(current.subgraph());
    for frame in prior_frames {
        if given.subgraph().intersects(&frame.active_subgraph)
            && !given.subgraph().includes(&frame.active_subgraph)
        {
            return Err(eyre!(
                "a nested direct local-3D prefix partially overlaps a prior coordinate frame"
            ));
        }
    }
    let retained_region = prior_active_subgraph
        .map(|active| active.intersection(given.subgraph()))
        .unwrap_or_else(|| graph.empty_subgraph());
    let mut retained_loop_edges = prior_frames
        .iter()
        .filter(|frame| retained_region.includes(&frame.active_subgraph))
        .flat_map(|frame| frame.lmb.loop_edges.iter().copied())
        .collect::<Vec<_>>();
    retained_loop_edges.sort();
    retained_loop_edges.dedup();
    let external: SuBitGraph = graph.external_filter();
    let full_internal = graph.full_filter().subtract(&external);
    let inactive_subgraph = current.subgraph().subtract(&active_subgraph);
    // Integrated lines no longer carry momentum in the active sector. A
    // retained quotient loop must therefore be checked in that same quotient,
    // rather than against its affine routing in the uncontracted graph.
    let quotient_reference;
    let reference_lmb = if inactive_subgraph.is_empty() {
        &graph.loop_momentum_basis
    } else {
        let inactive = InternalSubGraph::try_new(inactive_subgraph, graph.as_ref())
            .ok_or_else(|| eyre!("inactive direct local-3D prefix is not an internal subgraph"))?;
        let externals = graph.dummy_stripped_external_flows_of(&full_internal);
        quotient_reference = match graph.shrunken_sub_lmb(
            &full_internal,
            &inactive,
            externals.clone(),
            Some(&graph.loop_momentum_basis),
        ) {
            Ok(lmb) => lmb,
            Err(LmbError::NoShrunkenLmb { source, .. })
                if matches!(source.as_ref(), LmbError::NoCompatibleSubLmb { .. }) =>
            {
                // Contraction can retire every original carrier and turn a
                // retained tree edge into a new quotient loop.
                graph.shrunken_sub_lmb(&full_internal, &inactive, externals, None)?
            }
            Err(error) => return Err(error.into()),
        };
        &quotient_reference
    };
    let affine_retained_loop_edges = retained_loop_edges
        .iter()
        .copied()
        .filter(|edge| {
            reference_lmb.edge_signatures[*edge]
                .external
                .iter()
                .any(|sign| !sign.is_zero())
        })
        .collect_vec();
    if !affine_retained_loop_edges.is_empty() {
        return Err(eyre!(
            "component {} relative to {} cannot preserve retained direct local-3D coordinates {:?}: they require an affine graph-external momentum shift",
            current.subgraph().string_label(),
            given.subgraph().string_label(),
            affine_retained_loop_edges,
        ));
    }
    let is_homogeneous = |candidate: &LoopMomentumBasis| {
        candidate.loop_edges.iter().all(|edge| {
            reference_lmb.edge_signatures[*edge]
                .external
                .iter()
                .all(|sign| sign.is_zero())
        })
    };
    let is_compatible = |candidate: &LoopMomentumBasis| {
        retained_loop_edges
            .iter()
            .all(|edge| candidate.loop_edges.contains(edge))
    };

    // A forest node's own LMB is enumeration metadata.  Start instead from
    // the component route induced by the graph-canonical basis, so an
    // equivalent child representation cannot select a different parent
    // Taylor chart.  Only genuine retained coordinates may force the stable
    // fallback below. A zero-mode host records the same basis as its chart.
    let canonical_component_lmb = graph
        .canonical_component_lmb(current.subgraph())
        .wrap_err_with(|| {
            format!(
                "failed to induce the graph-canonical direct local-3D route for component {}",
                current.subgraph().string_label(),
            )
        })?;

    let quotient_lmb = if given.subgraph().is_empty() {
        canonical_component_lmb.clone()
    } else {
        let given = InternalSubGraph::try_new(given.subgraph().clone(), graph.as_ref())
            .ok_or_else(|| {
                eyre!(
                    "nested direct local-3D prefix is not an internal subgraph of {}",
                    current.subgraph().string_label(),
                )
            })?;
        graph.shrunken_sub_lmb(
            current.subgraph(),
            &given,
            graph.dummy_stripped_external_flows_of(current.subgraph()),
            None,
        )?
    };
    let expected_active_rank = retained_loop_edges.len() + quotient_lmb.loop_edges.len();

    // A sector can have an already-integrated or disconnected inactive prefix.
    // Contract that prefix only after fixing the enclosing coordinates, so the
    // quotient keeps the compatible subset of the prior sector's carriers.
    let coordinate_lmb = if active_subgraph == *current.subgraph() {
        if is_compatible(&canonical_component_lmb)
            && is_homogeneous(&canonical_component_lmb)
            && canonical_component_lmb.loop_edges.len() == expected_active_rank
        {
            canonical_component_lmb.clone()
        } else {
            graph
                .generate_loop_momentum_bases_of(current.subgraph())
                .into_iter()
                .filter(|candidate| candidate.loop_edges != canonical_component_lmb.loop_edges)
                .sorted_by_key(|candidate| candidate.loop_edges.iter().copied().collect_vec())
                .find(|candidate| {
                    is_compatible(candidate)
                        && is_homogeneous(candidate)
                        && candidate.loop_edges.len() == expected_active_rank
                })
                .ok_or_else(|| {
                    eyre!(
                        "no direct local-3D enclosing frame has retained rank {} plus quotient rank {}",
                        retained_loop_edges.len(),
                        quotient_lmb.loop_edges.len(),
                    )
                })?
        }
    } else {
        let shrunken_filter = current.subgraph().subtract(&active_subgraph);
        let shrunken =
            InternalSubGraph::try_new(shrunken_filter, graph.as_ref()).ok_or_else(|| {
                eyre!(
                    "inactive direct local-3D prefix is not an internal subgraph of {}",
                    current.subgraph().string_label(),
                )
            })?;
        let externals = graph.dummy_stripped_external_flows_of(current.subgraph());
        let mut failures = Vec::new();
        let mut try_enclosing = |enclosing_lmb: &LoopMomentumBasis| {
            if !is_compatible(enclosing_lmb) || !is_homogeneous(enclosing_lmb) {
                return None;
            }
            let enclosing_loop_edges = enclosing_lmb.loop_edges.clone();
            match graph.shrunken_sub_lmb(
                current.subgraph(),
                &shrunken,
                externals.clone(),
                Some(enclosing_lmb),
            ) {
                Ok(candidate)
                    if is_compatible(&candidate)
                        && is_homogeneous(&candidate)
                        && candidate.loop_edges.len() == expected_active_rank =>
                {
                    Some(candidate)
                }
                Ok(candidate) => {
                    failures.push(format!(
                        "enclosing {:?} contracted to incompatible {:?}",
                        enclosing_loop_edges, candidate.loop_edges,
                    ));
                    None
                }
                Err(error) => {
                    failures.push(format!(
                        "enclosing {:?} failed contraction: {error}",
                        enclosing_loop_edges,
                    ));
                    None
                }
            }
        };
        let selected = try_enclosing(&canonical_component_lmb).or_else(|| {
            graph
                .generate_loop_momentum_bases_of(current.subgraph())
                .into_iter()
                .filter(|candidate| candidate.loop_edges != canonical_component_lmb.loop_edges)
                .sorted_by_key(|candidate| candidate.loop_edges.iter().copied().collect_vec())
                .find_map(|candidate| try_enclosing(&candidate))
        });
        selected.ok_or_else(|| {
            eyre!(
                "no contracted direct local-3D frame retains exactly {:?} with total rank {}; attempts: {}",
                retained_loop_edges,
                expected_active_rank,
                failures.join("; "),
            )
        })?
    };
    if !is_compatible(&coordinate_lmb) {
        return Err(eyre!(
            "contracting a direct local-3D prefix lost the prior sector frame {:?}; quotient carriers are {:?}",
            retained_loop_edges,
            coordinate_lmb.loop_edges,
        ));
    }
    if !is_homogeneous(&coordinate_lmb) {
        return Err(eyre!(
            "direct local-3D sector frame {:?} for component {} requires an affine graph-external momentum shift",
            coordinate_lmb.loop_edges,
            current.subgraph().string_label(),
        ));
    }
    if coordinate_lmb.loop_edges.len() != expected_active_rank {
        return Err(eyre!(
            "direct local-3D sector frame has {} loop carriers, expected retained rank {} plus quotient rank {} = {}",
            coordinate_lmb.loop_edges.len(),
            retained_loop_edges.len(),
            quotient_lmb.loop_edges.len(),
            expected_active_rank,
        ));
    }

    debug_tags!(#generation, #uv, #local, #direct, #trace;
        stage = "direct_3d_coordinate_lmb",
        current = %current.log_display(),
        given = %given.log_display(),
        active_subgraph = %active_subgraph.string_label(),
        retained_region = %retained_region.string_label(),
        retained_loop_edges = ?retained_loop_edges,
        quotient_loop_edges = ?quotient_lmb.loop_edges,
        expected_active_rank,
        loop_edges = ?coordinate_lmb.loop_edges,
        "Selected one coordinate frame for a direct local-3D sector"
    );
    Ok((canonical_component_lmb, coordinate_lmb))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Local3DDeformation {
    Ordinary,
    /// U of a component acting on its own completed thermal soft branch. The
    /// branch scope replaces the mass scope, so enclosing operations co-scale
    /// the result together with the soft branch it subtracts from.
    OrdinaryOnSoftBranch,
    /// U on residue rows whose outermost completed zero-mode branch has soft
    /// lines in its support: a UV counterterm of their static theory, which
    /// introduces the static-theory masses `mIR` in place of `mUV`
    /// ([`SoftBranches::expands_static_cograph`]).
    OrdinaryOnStaticCograph,
    /// The co-scaled soft chart of the vacuum and zero-temperature IR
    /// subtraction, through relative order `dod - 1`.
    Soft,
    /// The Taylor chart of a finite-temperature zero-mode host: the soft
    /// chart at fixed physical temperature through the absolute `endpoint`
    /// `J - dod` of its order tuple. It also co-scales the zero-mode
    /// temperature of every class inside its support, and writes the order
    /// tag of the host in place of its Taylor parameter.
    ZeroMode {
        endpoint: i64,
    },
}

impl Local3DDeformation {
    fn is_ordinary(self) -> bool {
        matches!(
            self,
            Self::Ordinary | Self::OrdinaryOnSoftBranch | Self::OrdinaryOnStaticCograph
        )
    }

    /// Whether the chart expands in the external momenta at fixed masses,
    /// temperature and chemical potential, co-scaling the component's mass
    /// scope with its loop momenta.
    fn is_soft(self) -> bool {
        matches!(self, Self::Soft | Self::ZeroMode { .. })
    }

    /// Last retained power of the inverted Taylor parameter.
    fn series_endpoint(self) -> i64 {
        match self {
            Self::Ordinary | Self::OrdinaryOnSoftBranch | Self::OrdinaryOnStaticCograph => 0,
            Self::Soft => -1,
            Self::ZeroMode { endpoint } => endpoint,
        }
    }

    /// The value of the Taylor parameter in the finished series of the
    /// operation on `subgraph`: one, or the order tag of the zero-mode host
    /// `subgraph`, so that its term of relative order `j` carries
    /// `lambda^(j - dod)`, also when an integrated prefix contracted part of
    /// the host.
    fn order_tag(self, subgraph: &SuBitGraph) -> Atom {
        match self {
            Self::ZeroMode { .. } => zero_mode_order(Atom::var(subgraph.symbol())),
            Self::Ordinary
            | Self::OrdinaryOnSoftBranch
            | Self::OrdinaryOnStaticCograph
            | Self::Soft => Atom::num(1),
        }
    }
}

/// Terminal masses of an ordinary on-shell-energy projection, in the scope
/// of the operation that applies it.
#[derive(Clone)]
struct OrdinaryEnergyMasses {
    vacuum_squared: Atom,
    expansion_squared: Atom,
}

impl OrdinaryEnergyMasses {
    /// The masses that the ordinary U of `deformation` introduces, in the
    /// scope `mass_scope`: the static-theory pair `mIR`, `mIRexp` on a static
    /// cograph, and the vacuum pair `mUV`, `mUVexp` otherwise.
    fn new(deformation: Local3DDeformation, mass_scope: &Atom) -> Self {
        let (vacuum, expansion) = match deformation {
            Local3DDeformation::OrdinaryOnStaticCograph => (GS.m_ir_vacuum, GS.m_ir_expansion),
            Local3DDeformation::Ordinary
            | Local3DDeformation::OrdinaryOnSoftBranch
            | Local3DDeformation::Soft
            | Local3DDeformation::ZeroMode { .. } => (GS.m_uv_vacuum, GS.m_uv_expansion),
        };
        Self {
            vacuum_squared: (Atom::var(vacuum) * mass_scope).pow(2),
            expansion_squared: (Atom::var(expansion) * mass_scope).pow(2),
        }
    }

    /// Normalize a rescaled invariant before inversion so its Laurent series
    /// never contains a square root of a negative parameter power.
    fn normalize(&self, invariant: Atom) -> Atom {
        let rescale_squared = Atom::var(GS.rescale).pow(2);
        // The terminal masses an inner U leaves in the energies it projected:
        // the vacuum mass, and the static-theory mass of a U on a static
        // cograph.
        let terminal = [GS.m_uv_vacuum, GS.m_ir_vacuum];
        if terminal.iter().any(|mass| invariant.contains_symbol(*mass)) {
            // An inner U already supplied this terminal vacuum
            // energy. Scale its actual vacuum mass with the
            // active loop before factoring the energy, rather
            // than applying the MUV rearrangement a second time.
            // This relies on mUV and mIR being reserved for
            // terminal masses inside an active OSE argument.
            terminal.iter().fold(invariant, |invariant, mass| {
                invariant.replace(*mass).with(Atom::var(*mass) * GS.rescale)
            }) / rescale_squared
        } else {
            // A fresh projection keeps the expansion mass in
            // Taylor coefficients and uses the vacuum mass in
            // the terminal on-shell energy.
            (&self.vacuum_squared * &rescale_squared + invariant - &self.expansion_squared)
                / rescale_squared
        }
    }
}

/// Unit-valued support tag retained until the last local forest operation.
///
/// Physical and UV masses are global symbols, so equal masses in disconnected
/// components would otherwise be indistinguishable after an eager Taylor
/// projection.  The argument is a representative graph edge belonging to the
/// connected component that currently owns the occurrence.  A containing
/// component rebases the tag to its own representative, whereas a disjoint
/// component leaves it untouched.  Final-integrand construction replaces the
/// tag by one before numerical compilation.
pub(crate) static LOCAL_3D_MASS_SCOPE: LazyLock<Symbol> =
    LazyLock::new(|| symbol!("gammalooprs::local_3d_mass_scope"));

/// Unit-valued mass scope of a completed finite-temperature soft branch.
///
/// The argument is the `S_<label>` symbol of the expanded support, which also
/// replaces the integer owner of every on-shell energy the branch expanded.
/// A strictly enclosing local operation treats the branch as a frozen vertex:
/// it co-scales this scope, the branch's energies and the temperature of its
/// retained occupations with its own loop momenta, so that only the reduced
/// cograph is projected and expanded. The component's own U keeps the scope as
/// its expansion scope. Final-integrand construction replaces it by one.
pub(crate) static LOCAL_3D_SOFT_SCOPE: LazyLock<Symbol> =
    LazyLock::new(|| symbol!("gammalooprs::local_3d_soft_scope"));

/// Temperature `T_zm` of a static Matsubara mode, retained until the last
/// local forest operation.
///
/// A zero-mode projection writes the temperature of the frequencies it makes
/// static as this tag, so that local operations can scale it apart from the
/// temperature of the remaining occupations. The argument is the
/// representative edge of the projected frequency class, or the `S_<label>`
/// owner of the completed zero-mode branch that took the class over.
/// Final-integrand construction replaces the tag by `1/beta`.
pub(crate) static ZERO_MODE_TEMPERATURE: LazyLock<Symbol> =
    LazyLock::new(|| symbol!("gammalooprs::zero_mode_temperature"));

/// The zero-mode temperature [`ZERO_MODE_TEMPERATURE`] of `owner`: the
/// representative edge of a frequency class as an integer, or the `S_<label>`
/// owner of the completed zero-mode branch that took the class over.
pub(super) fn zero_mode_temperature(owner: Atom) -> Atom {
    function!(*ZERO_MODE_TEMPERATURE, owner)
}

/// The zero-mode temperature of the frequency class whose representative
/// edge is `representative`.
pub(super) fn class_zero_mode_temperature(representative: EdgeIndex) -> Atom {
    zero_mode_temperature(Atom::num(usize::from(representative)))
}

/// The representative edge of `graph` that the owner argument of a
/// zero-mode temperature names, if it names a frequency class rather than a
/// completed branch.
pub(super) fn zero_mode_class_representative(
    graph: &Graph,
    owner: AtomView<'_>,
) -> Option<EdgeIndex> {
    usize::try_from(owner)
        .ok()
        .filter(|edge| *edge < graph.underlying.n_edges())
        .map(EdgeIndex)
}

/// Unit-valued Taylor-order tag `lambda` of a completed zero-mode branch.
///
/// The argument is the `S_<label>` owner of the branch. A term of relative
/// Taylor order `j` of a branch whose support has superficial degree `dod`
/// carries `lambda^(j - dod)`, so that an enclosing operation expands the
/// frozen branch at the degree of its static cograph. Final-integrand
/// construction replaces the tag by one.
pub(crate) static ZERO_MODE_ORDER: LazyLock<Symbol> =
    LazyLock::new(|| symbol!("gammalooprs::zero_mode_order"));

/// The order tag [`ZERO_MODE_ORDER`] of the zero-mode branch `owner`, the
/// `S_<label>` symbol of its support.
pub(super) fn zero_mode_order(owner: Atom) -> Atom {
    function!(*ZERO_MODE_ORDER, owner)
}

/// The final values of the tags local operations keep apart while they are
/// composed: the mass scopes, the soft scopes and the zero-mode order tags
/// are one, and the zero-mode temperatures are the temperature `1/beta`.
pub(crate) fn normalize_local_tags(atom: &Atom) -> Atom {
    atom.replace(function!(*LOCAL_3D_MASS_SCOPE, W_.a_))
        .with(Atom::one())
        .replace(function!(*LOCAL_3D_SOFT_SCOPE, W_.a_))
        .with(Atom::one())
        .replace(function!(*ZERO_MODE_ORDER, W_.a_))
        .with(Atom::one())
        .replace(function!(*ZERO_MODE_TEMPERATURE, W_.a_))
        .with(Atom::var(GS.inverse_temperature).pow(-1))
}

/// The local operation of one forest node.
#[derive(Clone, Copy)]
enum NodeOperation<'r> {
    /// The ordinary Taylor projection U.
    Ordinary,
    /// `S + U - U S` with the co-scaled soft chart of the vacuum and
    /// zero-temperature IR subtraction.
    Soft,
    /// The zero-mode operation `S0` of a finite-temperature host, followed
    /// by `U - U S0` on a connected host.
    ZeroMode(&'r ThermalZeroModeRegion),
}

impl<'r> NodeOperation<'r> {
    /// The operation of `current` in `medium_mode`.
    ///
    /// # Errors
    ///
    /// Internal errors: at finite temperature an IR node of positive degree
    /// without the zero-mode region that the direct local-3D hedge-poset
    /// policy assigns, which demotes every other IR node, and a region
    /// outside finite temperature or on a non-IR node.
    ///
    /// # Panics
    ///
    /// For the deferred local on-shell scheme and for unsubtracted nodes,
    /// which never enter the wood.
    fn of<S: ForestNodeLike>(current: &'r S, medium_mode: MediumMode) -> Result<Self> {
        let finite_temperature = medium_mode.is_finite_temperature();
        if let Some(region) = current.thermal_zero_mode() {
            ensure!(
                finite_temperature && current.renormalization_scheme() == ApproximationType::IR,
                "internal error: the zero-mode host {} reached a {medium_mode:?} operation with \
                 scheme {:?}",
                current.subgraph().string_label(),
                current.renormalization_scheme(),
            );
            return Ok(Self::ZeroMode(region));
        }
        match current.renormalization_scheme() {
            ApproximationType::MUV | ApproximationType::PolePart => Ok(Self::Ordinary),
            ApproximationType::IR if current.dod() == 0 => Ok(Self::Ordinary),
            ApproximationType::IR if finite_temperature => Err(eyre!(
                "internal error: thermal soft subtraction needs a zero-mode host from the direct \
                 local-3D hedge-poset policy: the IR component {} of degree {} hosts no \
                 zero-mode region",
                current.subgraph().string_label(),
                current.dod(),
            )),
            ApproximationType::IR => Ok(Self::Soft),
            ApproximationType::OS => unimplemented!(
                "local on-shell counterterms are deferred until local counterterms can be derived from the 4D expanded representation"
            ),
            ApproximationType::Unsubtracted => panic!("should have been kept out of the wood"),
        }
    }
}

/// Apply one local-3D Taylor operation to a complete residue-map family.
/// Selectors remain outside the atoms. Newly attached factors share one
/// coefficient-parametric body and retain each complete key's argument row.
/// This kernel does not add the current operation's subtraction minus;
/// forest composition supplies it once, retaining signs already in the branch.
///
/// The operation acts on `active_subgraph`, the part of `current` that no
/// integrated prefix of the family contracted. A zero-mode host acts there
/// as on its whole support ([`ensure_zero_mode_chart`]). Once a prefix
/// contracted the whole connected host, the family is the host's own
/// add-back, whose integrated U leaves `S0` of the contracted host alone.
pub(super) fn apply_taylor<S: ForestNodeLike>(
    ctx: &UVCtx<'_>,
    orientation: OrientationProjection<'_>,
    current: &S,
    given: &S,
    active_subgraph: Option<SuBitGraph>,
    lmb: &LoopMomentumBasis,
    integrands: &DirectResidueBranches,
) -> Result<DirectResidueBranches> {
    let active_subgraph = active_subgraph
        .as_ref()
        .map(|active| active.intersection(current.subgraph()));
    let support = active_subgraph
        .as_ref()
        .unwrap_or_else(|| current.subgraph());
    let medium_mode = orientation.cff_options()?.medium_mode;
    let operation = NodeOperation::of(current, medium_mode)?;
    match operation {
        NodeOperation::Ordinary if medium_mode.is_finite_temperature() => {
            // U first projects every unfrozen active occupation to vacuum.
            // When that annihilates the input, as for the whole-graph U on a
            // chain without a frozen soft branch, skip the numerator and
            // Taylor work.
            let soft_branches =
                SoftBranches::new(ctx.graph, current.subgraph(), Local3DDeformation::Ordinary);
            let projected = integrands.fallible_map(|_, atom| {
                Direct3dApproximation::project_unfrozen_occupations(
                    ctx,
                    atom,
                    support,
                    soft_branches,
                )
            })?;
            if projected.is_zero() {
                return Ok(projected);
            }
        }
        NodeOperation::ZeroMode(region) => {
            // Every sector meets the host in the chart of its active support,
            // a vanishing one included, so that the local and the integrated
            // routes keep one invariant.
            ensure_zero_mode_chart(ctx, current, region, support, lmb)?;
            if integrands.is_zero() {
                return Ok(integrands.clone());
            }
        }
        NodeOperation::Ordinary | NodeOperation::Soft => {}
    }
    let started = start_family(
        ctx,
        orientation,
        current,
        given,
        active_subgraph.as_ref(),
        lmb,
        integrands,
    )?;
    let ordinary = |deformation, integrands: &DirectResidueBranches| {
        Direct3dApproximation::t(
            deformation,
            ctx,
            current,
            given,
            integrands,
            active_subgraph.as_ref(),
            lmb,
        )
    };
    match operation {
        NodeOperation::Ordinary => ordinary(Local3DDeformation::Ordinary, &started.branches),
        NodeOperation::Soft => {
            let soft = ordinary(Local3DDeformation::Soft, &started.branches)?;
            // The Taylor projector U is linear, so
            // U(X) + S(X) - U(S(X)) = S(X) + U(X - S(X)). Applying U
            // once to the completed soft remainder preserves the exact
            // U/S/US forest algebra while avoiding two separately
            // materialized, nearly cancelling UV branches.
            let soft_remainder = started.branches.zip_add(&-soft.clone())?;
            let uv_remainder = ordinary(Local3DDeformation::Ordinary, &soft_remainder)?;
            soft.zip_add(&uv_remainder)
        }
        NodeOperation::ZeroMode(region) => {
            let terms = ZeroModeTerm::expand(
                ctx,
                orientation,
                current,
                given,
                region,
                support,
                lmb,
                &started,
            )?;
            let soft = ZeroModeTerm::restore(&terms, |term| Ok(Cow::Borrowed(&term.frozen)))?;
            match region.kind() {
                // The U of each factor precedes a union host in the forest.
                ZeroModeHostKind::Union => Ok(soft),
                // A connected host whose whole support an integrated prefix
                // contracted acts on its own add-back, which already holds
                // its U: the integrated form of its `U S0` is `S0` of the
                // contracted host on that add-back, the self-suffix.
                ZeroModeHostKind::Connected if support.is_empty() => {
                    // Enclosing operations tell the rows of the branch by its
                    // owner (`Direct3dApproximation::t`). Without a loop of
                    // its own the branch names it only through its order tag
                    // and soft scope.
                    let owners = SoftBranches::new(
                        ctx.graph,
                        current.subgraph(),
                        Local3DDeformation::Ordinary,
                    )
                    .row_owners(&soft)?;
                    let unowned = soft.iter_keys().zip(&owners).any(|((_, roots), owners)| {
                        roots.iter().any(|(_, root)| !root.is_zero())
                            && !owners.contains(current.subgraph())
                    });
                    ensure!(
                        !unowned,
                        "unsupported zero-mode forest: the zero mode of the connected host {} on \
                         its own integrated add-back has a term of the host's degree, which no \
                         tag ties to its branch",
                        current.subgraph().string_label()
                    );
                    Ok(soft)
                }
                // The zero-mode branch is a frozen vertex for every enclosing
                // operation, whose U therefore subtracts only the reduced
                // cograph. U(X) and U(S0(X)) must keep separate expansion
                // scopes: the ordinary scope stays fixed under enclosing
                // operations, the branch scope is co-scaled with it.
                ZeroModeHostKind::Connected => {
                    let uv = ordinary(Local3DDeformation::Ordinary, &started.branches)?;
                    let uv_on_soft = ZeroModeTerm::restore(&terms, |term| {
                        Ok(Cow::Owned(ordinary(
                            Local3DDeformation::OrdinaryOnSoftBranch,
                            &term.frozen,
                        )?))
                    })?;
                    // An enclosing U tells the rows of a static cograph by
                    // the completed branches they hold (`Direct3dApproximation::t`),
                    // so no residue row may add the U image of the unfrozen
                    // integrand to a zero-mode term: the zero-mode terms are
                    // re-keyed to static soft energies, the image keeps the
                    // keys of the integrand.
                    let zero_mode_keys = soft
                        .iter_keys()
                        .chain(uv_on_soft.iter_keys())
                        .map(|(key, _)| key)
                        .collect::<HashSet<_>>();
                    ensure!(
                        uv.iter_keys().all(|(key, _)| !zero_mode_keys.contains(key)),
                        "internal error: a residue row of the zero-mode host {} adds the U image \
                         of its unfrozen integrand to a zero-mode term",
                        support.string_label()
                    );
                    soft.zip_add(&uv)?.zip_add(&-uv_on_soft)
                }
            }
        }
    }
}

/// The input of every chart of the operation of a node on its active
/// support.
struct StartedFamily {
    /// The residue family with the numerator that the node adds to its
    /// given prefix attached and its active energies started.
    branches: DirectResidueBranches,
    /// The mass scope of the operation.
    mass_scope: Atom,
}

/// The [`StartedFamily`] of the operation of `current` on `integrands`. A
/// zero-temperature medium first takes the vacuum limit of the added
/// occupations.
fn start_family<S: ForestNodeLike>(
    ctx: &UVCtx<'_>,
    orientation: OrientationProjection<'_>,
    current: &S,
    given: &S,
    active_subgraph: Option<&SuBitGraph>,
    lmb: &LoopMomentumBasis,
    integrands: &DirectResidueBranches,
) -> Result<StartedFamily> {
    let medium_mode = orientation.cff_options()?.medium_mode;
    let reduced = current.reduced_subgraph(given);
    let mut numerator = ctx
        .graph
        .color_reduced_numerator(&reduced, given.subgraph())
        .wrap_err_with(|| {
            format!(
                "the numerator of {} relative to {} is not available",
                reduced.string_label(),
                given.subgraph().string_label(),
            )
        })?;
    let lmb_id = lmb
        .loop_edges
        .first()
        .copied()
        .unwrap_or_else(|| current.lmb_id());
    let mass_scope = function!(*LOCAL_3D_MASS_SCOPE, usize::from(lmb_id) as i64);
    let mut physical_masses = Vec::new();
    for (pair, _, edge) in ctx.graph.iter_edges_of(current.subgraph()) {
        if !pair.is_paired() {
            continue;
        }
        let mass = edge.data.mass_atom();
        if !mass.is_zero() && !physical_masses.contains(&mass) {
            physical_masses.push(mass);
        }
    }
    // Tag only the numerator introduced by this component.  Replacing the
    // same model mass in the accumulated CFF atom would also tag a disjoint
    // sibling that happens to use the same particle species.
    for mass in physical_masses {
        numerator = numerator.replace(mass.clone()).with(mass * &mass_scope);
    }
    // At finite temperature the zero-mode operation S0 reads the occupations
    // of the completed CFF, and U projects the unfrozen ones to vacuum
    // itself. A zero-temperature medium takes the vacuum limit of the
    // occupations this component adds here.
    let vacuum_integrands;
    let integrands = if medium_mode == MediumMode::Vacuum || medium_mode.is_finite_temperature() {
        integrands
    } else {
        vacuum_integrands = integrands.map_expressions(|atom| {
            ctx.graph.make_thermal_distributions_explicit(
                atom,
                MediumMode::Vacuum,
                ctx.graph.iter_edges_of(&reduced).map(|(_, edge, _)| edge),
            )
        })?;
        &vacuum_integrands
    };
    let scope = DirectResidueBranches::numerator_scope();
    let integrands = integrands.multiply_key_mapped(orientation, ctx.graph, &numerator, scope)?;
    debug_tags!(#generation, #profile, #uv, #local, #direct, #trace;
        stage = "direct_3d_taylor_family",
        current = %current.log_display(),
        given = %given.log_display(),
        reduced = %reduced.string_label(),
        loop_edges = ?lmb.loop_edges,
        branch_count = integrands.iter_keys().count(),
        "Prepared a shared numerator before its Taylor kernel"
    );
    let branches = integrands
        .map_expressions(|atom| start(ctx, current, atom, &Atom::one(), active_subgraph, lmb))?;
    Ok(StartedFamily {
        branches,
        mass_scope,
    })
}

/// Every sector of a zero-mode host expands its active `support`, the host
/// less the part that an integrated prefix contracted, possibly all of it.
/// It does so in the host's canonical chart contracted by that part, the
/// chart that the zero bins and the enclosing hosts cancelling against the
/// sector share, which is empty once the whole host is contracted. The soft
/// set stays outside the support.
///
/// The host keeps its region, its owner and its degree on a contracted
/// support. An integrated prefix inside the host leaves its soft set, and
/// with it every frequency class, untouched. A prefix in the soft set never
/// reaches the zero-mode operation: its replay keeps the U of a connected
/// host alone ([`crate::uv::hedge_poset`]).
fn ensure_zero_mode_chart<S: ForestNodeLike>(
    ctx: &UVCtx<'_>,
    current: &S,
    region: &ThermalZeroModeRegion,
    support: &SuBitGraph,
    lmb: &LoopMomentumBasis,
) -> Result<()> {
    let graph = ctx.graph;
    let host = current.subgraph().string_label();
    let contracted = current.subgraph().subtract(support);
    let canonical = region
        .coordinate_loop_edges()
        .iter()
        .copied()
        .collect::<BTreeSet<_>>();
    let expected = if contracted.is_empty() {
        canonical
    } else if support.is_empty() {
        BTreeSet::new()
    } else {
        let canonical_lmb = graph.canonical_component_lmb(current.subgraph())?;
        let internal =
            InternalSubGraph::try_new(contracted.clone(), graph.as_ref()).ok_or_else(|| {
                eyre!(
                    "internal error: the integrated part of the zero-mode host {host} is not an \
                     internal subgraph"
                )
            })?;
        graph
            .shrunken_sub_lmb(
                current.subgraph(),
                &internal,
                graph.dummy_stripped_external_flows_of(current.subgraph()),
                Some(&canonical_lmb),
            )?
            .loop_edges
            .iter()
            .copied()
            .collect()
    };
    let loop_edges = lmb.loop_edges.iter().copied().collect::<BTreeSet<_>>();
    ensure!(
        loop_edges == expected,
        "internal error: the zero-mode host {host} must expand in its canonical chart {:?}, \
         contracted by its integrated part {}, on its support {}, not in {:?}",
        expected,
        contracted.string_label(),
        support.string_label(),
        loop_edges,
    );
    let static_lines = graph
        .iter_edges_of(support)
        .filter(|(pair, edge, _)| pair.is_paired() && region.sigma().contains(edge))
        .map(|(_, edge, _)| edge)
        .collect_vec();
    ensure!(
        static_lines.is_empty(),
        "internal error: the soft lines {static_lines:?} of the zero-mode host {host} lie in its \
         rescaled support",
    );
    Ok(())
}

/// One order tuple of a zero-mode operation `S0`: its Taylor coefficient,
/// re-keyed and frozen, and the static class propagators it multiplies.
///
/// The propagators stay apart while the host's own charts act on the
/// coefficient. They do not depend on the Taylor parameter, and multiplied
/// in only afterwards they remain factors of every summand, which an
/// enclosing projection reads off without expanding a frozen branch again:
/// top-level ones, until a later operation's Taylor coefficients carry them
/// in each term of a sum, through which the projection still reads them.
struct ZeroModeTerm {
    frozen: DirectResidueBranches,
    propagators: Atom,
}

impl ZeroModeTerm {
    /// The terms of `S0` of the zero-mode host `current` with its `region` on
    /// its `started` sector: for each tuple of the region, the zero-mode
    /// coefficient of its sector, Taylor-expanded through its order `J`,
    /// re-keyed to static soft energies and frozen.
    ///
    /// The Taylor series acts on the active `support` in the chart `lmb`,
    /// which [`ensure_zero_mode_chart`] checked; an integrated prefix
    /// contracted the rest of the host, and the series co-scales its owners
    /// as hard. The branch, its order tag `lambda^(j - dod)` and its degree
    /// belong to the whole host either way.
    ///
    /// The region has one order per class in each tuple, as
    /// [`ThermalZeroModeRegion`] builds them.
    ///
    /// # Errors
    ///
    /// The projection, Taylor series, re-keying or freezing fails, in the
    /// context of the host and its tuple.
    #[expect(
        clippy::too_many_arguments,
        reason = "the zero-mode operation keeps both forest-node contexts, its region, its active \
                  support and its chart explicit"
    )]
    fn expand<S: ForestNodeLike>(
        ctx: &UVCtx<'_>,
        orientation: OrientationProjection<'_>,
        current: &S,
        given: &S,
        region: &ThermalZeroModeRegion,
        support: &SuBitGraph,
        lmb: &LoopMomentumBasis,
        started: &StartedFamily,
    ) -> Result<Vec<Self>> {
        let graph = ctx.graph;
        let host = current.subgraph().string_label();
        debug_assert!(
            region
                .tuples()
                .iter()
                .all(|tuple| tuple.orders().len() == region.classes().len()),
            "an order tuple of the zero-mode host {host} has one order per class",
        );
        let sectors = region
            .tuples()
            .iter()
            .map(ZeroModeTuple::sector)
            .collect::<BTreeSet<_>>();
        let mut terms = Vec::new();
        for sector in sectors {
            let tuples = region
                .tuples()
                .iter()
                .filter(|tuple| tuple.sector() == sector)
                .collect_vec();
            let projection = region
                .classes()
                .iter()
                .enumerate()
                .map(|(index, class)| {
                    ZeroModeClass::new(
                        graph,
                        class.edges().iter().copied(),
                        tuples.iter().map(|tuple| tuple.orders()[index]),
                    )
                })
                .collect::<Result<Vec<_>, _>>()
                .and_then(|classes| ZeroModeProjection::new(classes, sector))
                .wrap_err_with(|| {
                    format!("the {sector:?} zero-mode projection of the host {host} is invalid")
                })?;
            let mut coefficients = projection
                .stripped_coefficients(graph, &started.branches)
                .wrap_err_with(|| {
                    format!("the {sector:?} zero-mode projection of the host {host} failed")
                })?;
            for tuple in tuples {
                let context = || format!("the order {tuple} of the zero-mode host {host}");
                // The chart's series sees the coefficient's scalars resolved.
                // As shared bodies they would bound the valuation of each
                // root only through their probes, several orders below the
                // root's own: the CFF variants of a bubble have poles at its
                // coincident energies that cancel only in their sum.
                let coefficient = coefficients
                    .remove(tuple.orders())
                    .ok_or_else(|| {
                        eyre!("internal error: the zero-mode projection of {host} omits {tuple}")
                    })?
                    .resolved_scalars()
                    .wrap_err_with(context)?;
                let deformation = Local3DDeformation::ZeroMode {
                    endpoint: tuple.taylor_order() - i64::from(current.dod()),
                };
                let frozen = Direct3dApproximation::t(
                    deformation,
                    ctx,
                    current,
                    given,
                    &coefficient,
                    Some(support),
                    lmb,
                )
                .and_then(|expanded| expanded.rekey_static(orientation, region.sigma()))
                .and_then(|expanded| {
                    Direct3dApproximation::freeze_soft_branch(
                        &expanded,
                        graph,
                        current.subgraph(),
                        lmb,
                        &started.mass_scope,
                    )
                })
                .wrap_err_with(context)?;
                terms.push(Self {
                    frozen,
                    propagators: projection.class_propagators(tuple.orders()),
                });
            }
        }
        Ok(terms)
    }

    /// The sum of `map` of every term times its static propagators.
    ///
    /// # Errors
    ///
    /// `map` fails, or there is no term.
    fn restore<'t>(
        terms: &'t [Self],
        mut map: impl FnMut(&'t Self) -> Result<Cow<'t, DirectResidueBranches>>,
    ) -> Result<DirectResidueBranches> {
        let mut restored = terms
            .iter()
            .map(|term| Ok(map(term)?.map(|root| root * &term.propagators)));
        let first = restored
            .next()
            .ok_or_else(|| eyre!("internal error: a zero-mode host requests no order"))??;
        restored.try_fold(first, |sum, term: Result<DirectResidueBranches>| {
            sum.zip_add(&term?)
        })
    }
}

fn start<S: ForestNodeLike>(
    ctx: &UVCtx<'_>,
    current: &S,
    cff: &Atom,
    mapped_numerator: &Atom,
    active_subgraph: Option<&SuBitGraph>,
    lmb: &LoopMomentumBasis,
) -> Result<Atom> {
    let graph = ctx.graph;
    let rescaled_subgraph = active_subgraph.unwrap_or_else(|| current.subgraph());
    let lmb_id = lmb
        .loop_edges
        .first()
        .copied()
        .unwrap_or_else(|| current.lmb_id());
    let mass_scope = function!(*LOCAL_3D_MASS_SCOPE, usize::from(lmb_id) as i64);
    let mut numerator = mapped_numerator.clone();
    let mut physical_masses = Vec::new();
    for (pair, _, edge) in graph.iter_edges_of(current.subgraph()) {
        if !pair.is_paired() {
            continue;
        }
        let mass = edge.data.mass_atom();
        if !mass.is_zero() && !physical_masses.contains(&mass) {
            physical_masses.push(mass);
        }
    }
    // Tag only the numerator introduced by this component.  Replacing the
    // same model mass in the accumulated CFF atom would also tag a disjoint
    // sibling that happens to use the same particle species.
    for mass in physical_masses {
        numerator = numerator.replace(mass.clone()).with(mass * &mass_scope);
    }
    let mut atomarg = cff * numerator;
    debug_tags!(#generation, #profile, #uv, #local, #trace;
        stage = "local_3d_start_initial",
        byte_size = atomarg.as_view().get_byte_size(),
        file.expr = %atomarg,
        "Local 3D start expression checkpoint"
    );
    // println!("CFF: {}", cff);

    // Keep the energy opaque, with its complete routed dependence in the arguments.
    for (p, ei, e) in graph.iter_edges_of(rescaled_subgraph) {
        let eid = usize::from(ei) as i64;
        if p.is_paired() {
            // set energies from inner_t on-shell
            atomarg = atomarg.replace(function!(GS.energy, eid)).with(GS.ose(ei));

            let e_mass = e.data.mass_atom() * &mass_scope;
            atomarg = atomarg
                .replace(GS.ose(ei))
                .with(GS.ose_full(ei, lmb_id, e_mass, None));
        }
    }
    debug_tags!(#generation, #profile, #uv, #local, #trace;
        stage = "local_3d_start_after_ose_full",
        byte_size = atomarg.as_view().get_byte_size(),
        file.expr = %atomarg,
        "Local 3D start expression checkpoint"
    );

    // split numerator momenta into OSEs and spatial parts
    let mut reps = Vec::new();
    for (p, eid, e) in graph.iter_edges_of(rescaled_subgraph) {
        if p.is_paired() {
            let e_mass = e.data.mass_atom() * &mass_scope;
            let rep = GS.split_mom_pattern(eid, lmb_id, e_mass);
            debug_tags!(#uv, #local, #momentum, #trace;
                stage = "local_3d_start_split_mom_pattern",
                split_rep = %rep,
                "Local 3D start momentum split"
            );
            reps.push(rep);
        }
    }
    let atomarg = atomarg.replace_multiple(&reps);
    debug_tags!(#generation, #profile, #uv, #local, #trace;
        stage = "local_3d_start_output",
        byte_size = atomarg.as_view().get_byte_size(),
        file.expr = %atomarg,
        "Local 3D start expression checkpoint"
    );
    Ok(atomarg.collect_compact_factors())
}

#[cfg(test)]
impl Local3DLoopRescaling {
    // #[debug_instrument(
    //     current = %current.log_display(),
    //     given = %given.log_display(),
    //     reduced,
    // )]
    pub(crate) fn t<S: ForestNodeLike>(
        self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrand: &Atom,
        active_subgraph: Option<&SuBitGraph>,
        lmb: &LoopMomentumBasis,
    ) -> Result<Atom> {
        self.project(
            Local3DDeformation::Ordinary,
            ctx,
            current,
            given,
            integrand,
            active_subgraph,
            lmb,
        )
    }

    #[expect(
        clippy::too_many_arguments,
        reason = "the projection must keep both forest-node contexts, its active component, and its canonical LMB explicit"
    )]
    // The former direct-external form rescaled the external momenta in the
    // added numerator subgraph. The hard-dual chart acts on the completed CFF
    // instead, while factoring active OSEs before taking their Laurent series.
    fn project<S: ForestNodeLike>(
        self,
        deformation: Local3DDeformation,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrand: &Atom,
        active_subgraph: Option<&SuBitGraph>,
        lmb: &LoopMomentumBasis,
    ) -> Result<Atom> {
        let atomarg = Direct3dApproximation::t_rescale(
            deformation,
            ctx,
            current,
            given,
            integrand,
            active_subgraph,
            lmb,
            true,
        )?;
        let rescaled_subgraph = active_subgraph.unwrap_or_else(|| current.subgraph());
        let endpoint = deformation.series_endpoint();
        let series = atomarg
            .series(GS.rescale, Atom::Zero, endpoint)
            .map_err(|error| {
                eyre!(
                    "failed to construct the local 3D {deformation:?} series for component {} relative to {}, active support {}, route loops {:?}, route externals {:?}: {error}",
                    current.subgraph().string_label(),
                    given.subgraph().string_label(),
                    rescaled_subgraph.string_label(),
                    lmb.loop_edges,
                    lmb.ext_edges,
                )
            })?;
        let series_atom = series.to_atom();
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_after_series",
            byte_size = series_atom.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );

        debug_tags!(#uv, #local;
            expr = %series,
            ?deformation,
            "After series in local 3D hard parameter"
        );
        let order_tag = deformation.order_tag(current.subgraph()).to_pattern();
        let a = series_atom.replace(GS.rescale).with(&order_tag);

        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_output",
            byte_size = a.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );
        debug_tags!(#uv, #local;
            log.expr = a,
            ?deformation,
            "Local 3D approximation"
        );
        Ok(a)
    }
}

impl Direct3dApproximation<'_> {
    #[allow(clippy::too_many_arguments)]
    fn t_rescale<S: ForestNodeLike>(
        deformation: Local3DDeformation,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrand: &Atom,
        active_subgraph: Option<&SuBitGraph>,
        lmb: &LoopMomentumBasis,
        include_measure: bool,
    ) -> Result<Atom> {
        let graph = ctx.graph;
        let reduced = current.reduced_subgraph(given);
        let rescaled_subgraph = active_subgraph.unwrap_or_else(|| current.subgraph());
        // The roles of completed soft branches are relative to the
        // operation's subgraph: the part of it that an integrated prefix
        // contracted keeps its owners, and the branches inside it stay frozen.
        let soft_branches = SoftBranches::new(graph, current.subgraph(), deformation);
        // U takes the vacuum limit of every active occupation except those
        // retained by a frozen soft branch, whose temperature is co-scaled
        // below instead. Disjoint thermal coefficients keep their owners.
        let vacuum_integrand;
        let integrand = if deformation.is_ordinary() {
            vacuum_integrand = Self::project_unfrozen_occupations(
                ctx,
                integrand,
                rescaled_subgraph,
                soft_branches,
            )?;
            &vacuum_integrand
        } else {
            integrand
        };
        let scope_edge = lmb
            .loop_edges
            .first()
            .copied()
            .unwrap_or_else(|| current.lmb_id());
        // U on its own soft branch expands with that branch's scope, which
        // enclosing operations co-scale; every other operation uses the
        // ordinary scope.
        let mass_scope = if deformation == Local3DDeformation::OrdinaryOnSoftBranch {
            thermal_soft::soft_branch_scope(current.subgraph())
        } else {
            function!(*LOCAL_3D_MASS_SCOPE, usize::from(scope_edge) as i64)
        };

        // Only apply replacements for edges in the reduced graph. Route both
        // canonical momentum images through the same component-local LMB:
        // represented spatial vectors are Q3, while tensor execution can
        // materialize their explicit components as Q(edge,cind(i)). This is
        // independent of whether the selected deformation is U or S.
        let mut mom_reps = graph.uv_spatial_wrapped_replacement(&reduced, lmb, &[W_.x___]);
        mom_reps.extend(graph.uv_wrapped_replacement(&reduced, lmb, &[W_.x___]));
        for m in &mom_reps {
            debug_tags!(#uv,#momentum,#trace;mom_rep=%m,"Mom rep");
        }

        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_input",
            active_subgraph = %rescaled_subgraph.string_label(),
            byte_size = integrand.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );
        let mut atomarg = integrand.replace_multiple(&mom_reps);
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_after_momentum_replacements",
            byte_size = atomarg.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );

        // A containing operation takes ownership of every nested mass scale in
        // its active coefficient, including finite effective vertices whose
        // source loops were integrated already. Their frozen localizers stay
        // outside this operation. Tags belonging to a disjoint accumulated
        // sector remain different and are therefore left untouched.
        for (pair, edge, _) in graph.iter_edges_of(current.subgraph()) {
            if pair.is_paired() {
                atomarg = atomarg
                    .replace(function!(*LOCAL_3D_MASS_SCOPE, usize::from(edge) as i64))
                    .with(mass_scope.clone());
            }
        }

        // Rescale every loop momentum still active in this sector, including
        // cycles expanded by earlier local operations.  In the soft mode the
        // component-local mass scope is co-scaled as well.  By homogeneity,
        // hard (k,m) scaling is exactly the direct external-momentum soft
        // Taylor deformation, but it remains local to this CFF orientation.
        for e in &lmb.loop_edges {
            // println!("Rescale {}", e);
            for momentum_symbol in [GS.emr_vec, GS.emr_mom] {
                let momentum = function!(momentum_symbol, usize::from(*e) as i64, W_.x___);
                let soft = lmb.ext_atom(*e, momentum_symbol, &[W_.x___], true);
                let rescaled = (&momentum - &soft) * GS.rescale + soft;
                atomarg = atomarg.replace(momentum).with(rescaled);
            }
        }
        if deformation.is_soft() {
            atomarg = atomarg
                .replace(mass_scope.clone())
                .with(mass_scope.clone() * GS.rescale);
        }
        // A frozen soft branch is a vertex of this operation: its masses and
        // zero-mode temperatures scale with the loop momenta, as the soft
        // deformation scaled them. A zero-mode chart also co-scales the
        // temperature of every static class inside it, whose loop is part of
        // the hard part it expands.
        let mut tag_error = None;
        atomarg = atomarg.replace_map(|view, _, output| {
            if tag_error.is_some() {
                return;
            }
            match soft_branches.co_scales(view) {
                Ok(true) => **output = view.to_owned() * GS.rescale,
                Ok(false) => {}
                Err(error) => tag_error = Some(error),
            }
        });
        if let Some(error) = tag_error {
            return Err(error);
        }
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_after_loop_rescale",
            byte_size = atomarg.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );

        // Match the 4D forest Taylor order within the current component: vacuum
        // masses in its completed coefficients have the hard loop-momentum weight.
        // Disjoint coefficients retain their owners for a later enclosing operation.
        // Scale the complete scalar head, leaving its owner independent of Taylor t.
        // Normalized localization kernels stay outside the series, and the auxiliary
        // OSE deformation mass keeps its separate role below.
        // A zero-mode chart is the soft expansion at fixed masses. It scales
        // the algebraic vacuum masses of a contained integrated vertex, whose
        // homogeneity turns that into the soft expansion, but holds their
        // logarithms fixed, as the prefix localizer does for physical masses
        // (`Direct3dApproximation::localize_prefix`): those are logarithms
        // of the vacuum mass over the MS-bar scale, which the soft expansion
        // leaves alone. The vertex's expansion mass and physical masses carry
        // its mass scope, which the soft charts co-scaled above.
        let zero_mode_chart = matches!(deformation, Local3DDeformation::ZeroMode { .. });
        let mut invalid_mass_owner = None;
        atomarg = atomarg.replace_map(|view, _, output| {
            let AtomView::Fun(mass) = view else {
                return;
            };
            if zero_mode_chart && mass.get_symbol() == Symbol::LOG {
                // Settled as it is: its arguments are not visited.
                **output = view.to_owned();
                return;
            }
            if mass.get_symbol() != GS.m_uv_vacuum {
                return;
            }
            let owner = (mass.get_nargs() == 1)
                .then(|| thermal_soft::component_owner(mass.get(0), graph.n_hedges()))
                .flatten();
            let Some(owner) = owner else {
                invalid_mass_owner = Some("a direct vacuum mass has no valid component owner");
                return;
            };
            if current.subgraph().includes(&owner) {
                **output = view.to_owned() * GS.rescale;
            } else if current.subgraph().intersects(&owner) {
                invalid_mass_owner = Some(
                    "a direct vacuum-mass owner partially overlaps the current forest component",
                );
            }
        });
        if let Some(message) = invalid_mass_owner {
            return Err(eyre!(message));
        }

        // Normalize each active OSE invariant before inversion so its Laurent
        // series never contains a square root of a negative parameter power.
        // Ordinary U performs the established MUV rearrangement, with the
        // static-theory mass on a static cograph. Soft S factors the
        // co-scaled physical energy without introducing MUV.
        let uv_masses = OrdinaryEnergyMasses::new(deformation, &mass_scope);
        for eid in &lmb.loop_edges {
            let eid = eid.0 as i64;
            let uv_masses = uv_masses.clone();
            atomarg = atomarg
                .replace(function!(GS.on_shell_energy, eid, W_.prop_))
                .with_map(move |matched| {
                    let invariant = matched.get(W_.prop_).unwrap().to_atom();
                    let normalized_invariant = if deformation.is_ordinary() {
                        uv_masses.normalize(invariant)
                    } else {
                        invariant / Atom::var(GS.rescale).pow(2)
                    };
                    function!(GS.on_shell_energy, eid, normalized_invariant) * GS.rescale
                });
        }
        // Energies of a completed soft branch are owned by its support symbol.
        // A frozen branch scales like the soft deformation that produced it,
        // while the component's own U expands its branch like a fresh energy.
        let mut tag_error = None;
        atomarg = atomarg.replace_map(|view, _, output| {
            let AtomView::Fun(energy) = view else {
                return;
            };
            if tag_error.is_some() || energy.get_symbol() != GS.on_shell_energy {
                return;
            }
            let normalized_invariant = match soft_branches.tag_role(view) {
                Ok(Some(SoftBranchRole::Frozen)) => {
                    energy.get(1).to_owned() / Atom::var(GS.rescale).pow(2)
                }
                Ok(Some(SoftBranchRole::Own)) => uv_masses.normalize(energy.get(1).to_owned()),
                Ok(Some(SoftBranchRole::Disjoint) | None) => return,
                Err(error) => {
                    tag_error = Some(error);
                    return;
                }
            };
            **output =
                function!(GS.on_shell_energy, energy.get(0), normalized_invariant) * GS.rescale;
        });
        if let Some(error) = tag_error {
            return Err(error);
        }
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_after_ose_rescale",
            byte_size = atomarg.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );

        atomarg = Self::compensate_thermal_distributions(
            ctx,
            rescaled_subgraph,
            &atomarg,
            deformation,
            soft_branches,
        )?;

        // The supplied LMB is the integration-space authority. In particular,
        // a remainder which is a tree in the original incidence can become a
        // loop after its frozen UV prefix is contracted.
        if include_measure {
            atomarg *= Atom::var(GS.rescale).pow(3 * lmb.loop_edges.len() as i64);
        }
        atomarg = atomarg.replace(GS.rescale).with(Atom::num(1) / GS.rescale);
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_before_series",
            loop_edges = ?lmb.loop_edges,
            byte_size = atomarg.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );

        debug_tags!(#uv, #local, #before_series; log.expr = atomarg, "Before series in t");

        let started = std::time::Instant::now();
        let compact = atomarg.collect_compact_factors();
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_after_compaction",
            input_bytes = atomarg.as_view().get_byte_size(),
            byte_size = compact.as_view().get_byte_size(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Local 3D T input compaction completed"
        );
        Ok(compact)
    }

    /// The Taylor series of the chart `deformation` on `integrands`. An
    /// ordinary U is a UV counterterm of a static theory on the residue rows
    /// whose outermost completed zero-mode branch has soft lines in its
    /// support ([`SoftBranches::expands_static_cograph`]), and acts there
    /// as [`Local3DDeformation::OrdinaryOnStaticCograph`]; it stays ordinary
    /// on the other rows. A U never expands a row partly as either: a row
    /// holding a completed branch is a re-keyed static row of its `S0`, apart
    /// from every row of the unfrozen integrand, as [`apply_taylor`] checks
    /// where a connected host adds both kinds of rows.
    fn t<S: ForestNodeLike>(
        deformation: Local3DDeformation,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrands: &DirectResidueBranches,
        active_subgraph: Option<&SuBitGraph>,
        lmb: &LoopMomentumBasis,
    ) -> Result<DirectResidueBranches> {
        let chart = |deformation, integrands: &DirectResidueBranches| {
            Self::chart_series(
                deformation,
                ctx,
                current,
                given,
                integrands,
                active_subgraph,
                lmb,
            )
        };
        if deformation != Local3DDeformation::Ordinary {
            return chart(deformation, integrands);
        }
        let soft_branches = SoftBranches::new(ctx.graph, current.subgraph(), deformation);
        let mut loop_lines = None;
        let static_rows = soft_branches
            .row_owners(integrands)?
            .iter()
            .map(|owners| {
                if owners.is_empty() {
                    return Ok(false);
                }
                let loop_lines = loop_lines.get_or_insert_with(|| ctx.graph.uncut_loop_lines());
                soft_branches.expands_static_cograph(owners, loop_lines)
            })
            .collect::<Result<Vec<_>>>()?;
        if static_rows.iter().all(|is_static| !is_static) {
            return chart(deformation, integrands);
        }
        if static_rows.iter().all(|is_static| *is_static) {
            return chart(Local3DDeformation::OrdinaryOnStaticCograph, integrands);
        }
        // The two charts introduce different masses, so each part keeps its
        // own copy of the numerator and scalar bodies its rows call.
        let (static_part, ordinary_part): (Vec<_>, Vec<_>) = integrands
            .iter_keys()
            .zip(&static_rows)
            .partition(|(_, is_static)| **is_static);
        let family = |part: Vec<((&DirectResidueKey, &Integrands), &bool)>| {
            DirectResidueBranches::from_keyed(
                part.into_iter()
                    .map(|((key, roots), _)| (key.clone(), roots.clone())),
            )
        };
        chart(
            Local3DDeformation::OrdinaryOnStaticCograph,
            &family(static_part)?,
        )?
        .zip_add(&chart(deformation, &family(ordinary_part)?)?)
    }

    /// [`Self::t`] of one chart on every row of `integrands`.
    fn chart_series<S: ForestNodeLike>(
        deformation: Local3DDeformation,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrands: &DirectResidueBranches,
        active_subgraph: Option<&SuBitGraph>,
        lmb: &LoopMomentumBasis,
    ) -> Result<DirectResidueBranches> {
        let rescaled = integrands
            .fallible_map(|_, atom| {
                Self::t_rescale(
                    deformation,
                    ctx,
                    current,
                    given,
                    atom,
                    active_subgraph,
                    lmb,
                    true,
                )
            })?
            .map_numerators(|entry| {
                Self::t_rescale(
                    deformation,
                    ctx,
                    current,
                    given,
                    &entry.rhs,
                    active_subgraph,
                    lmb,
                    false,
                )
            })?;
        let endpoint = deformation.series_endpoint();
        let scope = DirectResidueBranches::numerator_scope().1;
        let center = Atom::Zero;
        // A zero-mode chart keeps a few negative orders only; its roots are
        // probed through its endpoint.
        let series = if matches!(deformation, Local3DDeformation::ZeroMode { .. }) {
            rescaled.series_preserving_numerators_to_depth(
                GS.rescale,
                center.as_view(),
                endpoint,
                scope,
            )
        } else {
            rescaled.series_preserving_numerators(GS.rescale, center.as_view(), endpoint, scope)
        }
        .map_err(|error| eyre!(
            "local 3D {deformation:?} Taylor series through order {endpoint} failed for graph `{}` at {} given {}, in loop coordinates {:?}: {error}",
            ctx.graph.name, current.subgraph().string_label(),
            given.subgraph().string_label(), lmb.loop_edges,
        ))?;
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_after_series",
            branch_count = series.iter_keys().count(),
            "Shared local 3D Taylor coefficients prepared"
        );
        let series = if deformation.is_ordinary() {
            // Recollect the known Laurent polynomial only in its bookkeeping
            // scale. Retained graph numerators remain opaque, and each residue
            // map/cut keeps its own scalar cancellation problem.
            series.fallible_map(|_key, atom| {
                let started = std::time::Instant::now();
                if atom.as_view().get_byte_size() > 1024 * 1024 {
                    return Ok(atom.clone());
                }
                let Some(coefficients) = atom.coefficient_list_exact(&[Atom::var(GS.rescale)])
                else {
                    return Ok(atom.clone());
                };
                let result =
                    Atom::add_many(coefficients.into_iter().map(|(power, coefficient)| {
                        let keys = numerator_family_calls(coefficient.as_view());
                        power * coefficient.cancel_scalar_poles(&keys)
                    }));
                let result = if result.as_view().get_byte_size() < atom.as_view().get_byte_size() {
                    #[cfg(test)]
                    pole_audit::record(
                        ctx,
                        current,
                        given,
                        _key,
                        &series,
                        atom,
                        &result,
                        active_subgraph,
                        lmb,
                    );
                    result
                } else {
                    atom.clone()
                };
                debug_tags!(#generation, #profile, #uv, #local, #summary;
                    stage = "local_3d_u_scalar_cancellation",
                    input_bytes = atom.as_view().get_byte_size(),
                    output_bytes = result.as_view().get_byte_size(),
                    elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
                    "Cancelled scalar poles in UV Laurent coefficients"
                );
                Ok(result)
            })?
        } else {
            series
        };
        // A zero-mode branch is tagged by its host, also on a contracted
        // support.
        let order_tag = deformation.order_tag(current.subgraph()).to_pattern();
        series.map_expressions(|atom| {
            let a = atom
                .replace(GS.rescale)
                .with(&order_tag)
                .collect_compact_factors();
            debug_tags!(#uv, #local; log.expr = a, "Local 3D approximation");
            Ok(a)
        })
    }
}

#[cfg(test)]
#[path = "soft_tests.rs"]
mod soft_tests;

#[cfg(test)]
#[path = "test_support.rs"]
mod test_support;

#[path = "thermal_soft.rs"]
mod thermal_soft;

#[path = "zero_mode.rs"]
mod zero_mode;

#[cfg(test)]
#[path = "pole_audit.rs"]
mod pole_audit;

#[cfg(test)]
pub(crate) use pole_audit::physical_full_h_fixture;
