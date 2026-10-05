use std::sync::LazyLock;

use color_eyre::Result;
use eyre::eyre;
use idenso::color::{ColorSimplifier, ColorSimplifySettings};
use itertools::Itertools;
use linnet::half_edge::subgraph::{Inclusion, InternalSubGraph, SuBitGraph, SubSetLike, SubSetOps};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    function, symbol,
};

use crate::{
    debug_tags,
    graph::{LMBext, LmbError, LoopMomentumBasis},
    numerator::symbolica_ext::NumeratorAtomExt,
    utils::{GS, W_},
    uv::{
        ApproximationType, UltravioletGraph,
        approx::{ForestNodeLike, OrientationProjection, UVCtx},
        uv_graph::UVE,
    },
};

#[cfg(test)]
use crate::uv::approx::local_3d::Local3DLoopRescaling;

use super::{branches::DirectResidueBranches, forest::Direct3dApproximation};

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
    // fallback below.
    let canonical_component_lmb = if current.subgraph() == &full_internal
        && graph.n_loops(&full_internal) == graph.loop_momentum_basis.loop_edges.len()
    {
        graph.loop_momentum_basis.clone()
    } else {
        graph
            .try_compatible_sub_lmb(
                current.subgraph(),
                graph.dummy_less_full_crown(current.subgraph()),
                &graph.loop_momentum_basis,
            )
            .map_err(|error| {
                eyre!(
                    "failed to induce the graph-canonical direct local-3D route for component {}: {error}",
                    current.subgraph().string_label(),
                )
            })?
    };

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

#[derive(Clone, Copy, Debug)]
enum Local3DDeformation {
    Ordinary,
    Soft,
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

/// Apply one local-3D Taylor operation to a complete residue-map family.
/// Selectors remain outside the atoms. Newly attached factors share one
/// coefficient-parametric body and retain each complete key's argument row.
/// This kernel does not add the current operation's subtraction minus;
/// forest composition supplies it once, retaining signs already in the branch.
// Keep the forest operation, coordinate frame, and residue family explicit.
#[allow(clippy::too_many_arguments)]
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
    let reduced = current.reduced_subgraph(given);
    let mut numerator = ctx
        .graph
        .numerator(&reduced, given.subgraph())
        .get_single_atom()
        .expect("graph numerator should be available")
        .simplify_color_with(ColorSimplifySettings {
            simplify_non_color: false,
            ..Default::default()
        });
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
    let started = integrands.map_expressions(|atom| {
        start(
            ctx,
            current,
            atom,
            &Atom::one(),
            active_subgraph.as_ref(),
            lmb,
        )
    })?;
    match current.renormalization_scheme() {
        ApproximationType::MUV | ApproximationType::PolePart => Direct3dApproximation::t(
            Local3DDeformation::Ordinary,
            ctx,
            current,
            given,
            &started,
            active_subgraph.as_ref(),
            lmb,
        ),
        ApproximationType::IR => {
            if current.dod() == 0 {
                return Direct3dApproximation::t(
                    Local3DDeformation::Ordinary,
                    ctx,
                    current,
                    given,
                    &started,
                    active_subgraph.as_ref(),
                    lmb,
                );
            }
            let soft = Direct3dApproximation::t(
                Local3DDeformation::Soft,
                ctx,
                current,
                given,
                &started,
                active_subgraph.as_ref(),
                lmb,
            )?;
            // The Taylor projector U is linear, so
            // U(X) + S(X) - U(S(X)) = S(X) + U(X - S(X)). Applying U
            // once to the completed soft remainder preserves the exact
            // U/S/US forest algebra while avoiding two separately
            // materialized, nearly cancelling UV branches.
            let soft_remainder = started.zip_add(&-soft.clone())?;
            let uv_remainder = Direct3dApproximation::t(
                Local3DDeformation::Ordinary,
                ctx,
                current,
                given,
                &soft_remainder,
                active_subgraph.as_ref(),
                lmb,
            )?;
            soft.zip_add(&uv_remainder)
        }
        ApproximationType::VaccuumLimit => Err(eyre!("Not yet implemented VaccuumLimit")),
        ApproximationType::OS => unimplemented!(
            "local on-shell counterterms are deferred until local counterterms can be derived from the 4D expanded representation"
        ),
        ApproximationType::Unsubtracted => panic!("should have been kept out of the wood"),
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
        let endpoint = match deformation {
            Local3DDeformation::Ordinary => 0,
            Local3DDeformation::Soft => -1,
        };
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
        let a = series_atom.replace(GS.rescale).with(Atom::num(1));

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
        let scope_edge = lmb
            .loop_edges
            .first()
            .copied()
            .unwrap_or_else(|| current.lmb_id());
        let mass_scope = function!(*LOCAL_3D_MASS_SCOPE, usize::from(scope_edge) as i64);

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
        if matches!(deformation, Local3DDeformation::Soft) {
            atomarg = atomarg
                .replace(mass_scope.clone())
                .with(mass_scope.clone() * GS.rescale);
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
        let mut invalid_mass_owner = None;
        atomarg = atomarg.replace_map(|view, _, output| {
            let AtomView::Fun(mass) = view else {
                return;
            };
            if mass.get_symbol() != GS.m_uv_vacuum {
                return;
            }
            let owner = match mass.get_nargs() {
                1 => match mass.get(0) {
                    AtomView::Var(owner) => owner
                        .get_symbol()
                        .get_stripped_name()
                        .strip_prefix("S_")
                        .and_then(|label| SuBitGraph::from_base62(label, graph.n_hedges())),
                    _ => None,
                },
                _ => None,
            };
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
        // Ordinary U performs the established MUV rearrangement. Soft S
        // factors the co-scaled physical energy without introducing MUV.
        let rescale_squared = Atom::var(GS.rescale).pow(2);
        let uv_expansion_mass_squared = (Atom::var(GS.m_uv_expansion) * &mass_scope).pow(2);
        let uv_vacuum_mass_squared = (Atom::var(GS.m_uv_vacuum) * &mass_scope).pow(2);
        for eid in &lmb.loop_edges {
            let eid = eid.0 as i64;
            let rescale_squared = rescale_squared.clone();
            let uv_expansion_mass_squared = uv_expansion_mass_squared.clone();
            let uv_vacuum_mass_squared = uv_vacuum_mass_squared.clone();
            atomarg = atomarg
                .replace(function!(GS.energy_surface, eid, W_.prop_))
                .with_map(move |matched| {
                    let invariant = matched.get(W_.prop_).unwrap().to_atom();
                    let normalized_invariant = match deformation {
                        Local3DDeformation::Ordinary
                            if invariant.contains_symbol(GS.m_uv_vacuum) =>
                        {
                            // An inner U already supplied this terminal vacuum
                            // energy. Scale its actual vacuum mass with the
                            // active loop before factoring the energy, rather
                            // than applying the MUV rearrangement a second time.
                            // This relies on mUV being reserved for terminal
                            // masses inside an active OSE argument.
                            invariant
                                .replace(GS.m_uv_vacuum)
                                .with(Atom::var(GS.m_uv_vacuum) * GS.rescale)
                                / &rescale_squared
                        }
                        Local3DDeformation::Ordinary => {
                            // A fresh projection keeps the expansion mass in
                            // Taylor coefficients and uses the vacuum mass in
                            // the terminal on-shell energy.
                            (&uv_vacuum_mass_squared * &rescale_squared + invariant
                                - &uv_expansion_mass_squared)
                                / &rescale_squared
                        }
                        Local3DDeformation::Soft => invariant / &rescale_squared,
                    };
                    function!(GS.energy_surface, eid, normalized_invariant) * GS.rescale
                });
        }
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_after_ose_rescale",
            byte_size = atomarg.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );

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

    #[allow(clippy::too_many_arguments)]
    fn t<S: ForestNodeLike>(
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
        let endpoint = match deformation {
            Local3DDeformation::Ordinary => 0,
            Local3DDeformation::Soft => -1,
        };
        let series = rescaled.series_preserving_numerators(
            GS.rescale, Atom::Zero.as_view(), endpoint, DirectResidueBranches::numerator_scope().1,
        ).map_err(|error| eyre!(
            "local 3D {deformation:?} Taylor series through order {endpoint} failed for graph `{}` at {} given {}, in loop coordinates {:?}: {error}",
            ctx.graph.name, current.subgraph().string_label(),
            given.subgraph().string_label(), lmb.loop_edges,
        ))?;
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_after_series",
            branch_count = series.iter_keys().count(),
            "Shared local 3D Taylor coefficients prepared"
        );
        let series = if matches!(deformation, Local3DDeformation::Ordinary) {
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
                        let mut keys = std::collections::BTreeSet::new();
                        let family = symbol!("gammalooprs::uv::numerator_family");
                        coefficient.visitor(&mut |part| {
                            if let AtomView::Fun(call) = part
                                && call.get_symbol() == family
                            {
                                keys.insert(part.to_owned());
                                return false;
                            }
                            true
                        });
                        let keys = keys.into_iter().collect_vec();
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
        series.map_expressions(|atom| {
            let a = atom
                .replace(GS.rescale)
                .with(Atom::num(1))
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
#[path = "pole_audit.rs"]
mod pole_audit;

#[cfg(test)]
pub(crate) use pole_audit::physical_full_h_fixture;
