use std::sync::LazyLock;

use color_eyre::Result;
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::{
    involution::HedgePair,
    subgraph::{Inclusion, InternalSubGraph, SuBitGraph, SubSetLike, SubSetOps},
};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder, Symbol},
    function,
    id::Replacement,
};

use crate::{
    debug_tags,
    graph::{LMBext, LoopMomentumBasis},
    utils::{GS, W_},
    uv::{
        ApproximationType, UltravioletGraph,
        approx::{ForestNodeLike, OrientationProjection, UVCtx},
        uv_graph::UVE,
    },
};

use super::branches::DirectResidueKey;

/// One independently framed connected part of a direct forest sector.
/// Disconnected replay keeps these frames separate until an enclosing Taylor
/// operation contains and combines their carriers.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct DirectCoordinateFrame {
    pub(crate) active_subgraph: SuBitGraph,
    pub(crate) lmb: LoopMomentumBasis,
}

#[derive(Clone, Copy, Debug)]
pub(super) enum Local3DLoopRescaling {
    FullSubgraph,
    ReducedSubgraph,
}

static OSE_FOR_LOCAL_3D_SERIES: LazyLock<Symbol> = LazyLock::new(|| {
    symbolica::symbol!(
        "gammalooprs::OSE_for_local_3d_series",
        der = |_, arg, out| {
            if arg == 2 {
                **out = Atom::num(1);
            } else {
                **out = Atom::Zero;
            }
        }
    )
});

/// Choose the loop coordinates in which one direct sector enters its next
/// Taylor operation. An enclosing operation must retain every loop carrier
/// already owned by the part of the sector contained in `given`;
/// otherwise an equivalent basis can turn a previously hard child momentum
/// into a crown-shifted momentum before the outer series is taken.
pub(super) fn coordinate_lmb<S: ForestNodeLike>(
    ctx: &UVCtx<'_>,
    current: &S,
    given: &S,
    prior_active_subgraph: Option<&SuBitGraph>,
    prior_frames: &[DirectCoordinateFrame],
    active_subgraph: &SuBitGraph,
) -> Result<LoopMomentumBasis> {
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
        .flat_map(|frame| frame.lmb.loop_edges.iter())
        .filter(|edge| retained_region.includes(&graph[*edge].1))
        .copied()
        .collect::<Vec<_>>();
    retained_loop_edges.sort();
    retained_loop_edges.dedup();
    let is_compatible = |candidate: &LoopMomentumBasis| {
        retained_loop_edges
            .iter()
            .all(|edge| candidate.loop_edges.contains(edge))
            && candidate
                .loop_edges
                .iter()
                .filter(|edge| retained_region.includes(&graph[*edge].1))
                .count()
                == retained_loop_edges.len()
    };

    let quotient_lmb = if given.subgraph().is_empty() {
        current.lmb().clone()
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
        if is_compatible(current.lmb()) && current.lmb().loop_edges.len() == expected_active_rank {
            current.lmb().clone()
        } else {
            graph
                .generate_loop_momentum_bases_of(current.subgraph())
                .into_iter()
                .filter(|candidate| candidate.loop_edges != current.lmb().loop_edges)
                .sorted_by_key(|candidate| candidate.loop_edges.iter().copied().collect_vec())
                .find(|candidate| {
                    is_compatible(candidate)
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
            if !is_compatible(enclosing_lmb) {
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
        let selected = try_enclosing(current.lmb()).or_else(|| {
            graph
                .generate_loop_momentum_bases_of(current.subgraph())
                .into_iter()
                .filter(|candidate| candidate.loop_edges != current.lmb().loop_edges)
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
    Ok(coordinate_lmb)
}

/// Apply one local-3D Taylor operation to a complete generalized residue-map
/// branch. The selector remains outside the atom, and every newly attached
/// factor uses the branch's one authoritative energy substitution.
pub(super) fn apply_taylor<S: ForestNodeLike>(
    rescaling: Local3DLoopRescaling,
    ctx: &UVCtx<'_>,
    orientation: OrientationProjection<'_>,
    current: &S,
    given: &S,
    active_subgraph: Option<SuBitGraph>,
    lmb: &LoopMomentumBasis,
    key: &DirectResidueKey,
    integrand: &Atom,
) -> Result<Atom> {
    let active_subgraph = active_subgraph
        .as_ref()
        .map(|active| active.intersection(current.subgraph()));
    let reduced = current.reduced_subgraph(given);
    debug_tags!(#generation, #profile, #uv, #local, #direct, #trace;
        stage = "direct_3d_taylor_branch",
        current = %current.log_display(),
        given = %given.log_display(),
        reduced = %reduced.string_label(),
        active_subgraph = ?active_subgraph.as_ref().map(SubSetLike::string_label),
        residue_map_key = key.selector_host.0,
        source_edge_energy_map = ?key.source_edge_energy_map(),
        loop_edges = ?lmb.loop_edges,
        file.integrand = %integrand,
        "Direct local-3D branch entering one Taylor kernel"
    );

    let numerator = ctx
        .graph
        .numerator(&reduced, given.subgraph())
        .get_single_atom()
        .expect("graph numerator should be available");
    let mapped_numerator = key.map_numerator(orientation, ctx.graph, &numerator)?;
    debug_tags!(#generation, #profile, #uv, #local, #direct, #trace;
        stage = "direct_3d_taylor_mapped_numerator",
        current = %current.log_display(),
        given = %given.log_display(),
        residue_map_key = key.selector_host.0,
        file.mapped_numerator = %mapped_numerator,
        "Mapped the direct branch-owned numerator before its Taylor kernel"
    );
    let source_map = key.source_edge_energy_map();
    match current.renormalization_scheme() {
        ApproximationType::MUV | ApproximationType::PolePart => {
            let started = start(
                ctx,
                current,
                integrand,
                &mapped_numerator,
                active_subgraph.as_ref(),
                lmb,
            )?;
            debug_tags!(#generation, #profile, #uv, #local, #direct, #trace;
                stage = "direct_3d_taylor_started",
                current = %current.log_display(),
                given = %given.log_display(),
                residue_map_key = key.selector_host.0,
                file.started = %started,
                "Direct branch after numerator attachment and before rescaling"
            );
            debug_tags!(#generation, #profile, #uv, #local, #direct, #summary;
                stage = "direct_3d_kernel_after_start",
                input_byte_size = integrand.as_view().get_byte_size(),
                output_byte_size = started.as_view().get_byte_size(),
                "Direct local-3D kernel size checkpoint"
            );
            rescaling.t(ctx, current, given, &started, active_subgraph.as_ref(), lmb)
        }
        ApproximationType::IR => {
            let t_tilde = t_tilde(
                ctx,
                orientation,
                key.selector_host,
                source_map,
                current,
                given,
                integrand,
                active_subgraph.as_ref(),
                lmb,
            )?;
            Ok(rescaling.t(
                ctx,
                current,
                given,
                &start(
                    ctx,
                    current,
                    integrand,
                    &mapped_numerator,
                    active_subgraph.as_ref(),
                    lmb,
                )?,
                active_subgraph.as_ref(),
                lmb,
            )? + &t_tilde
                - rescaling.t(ctx, current, given, &t_tilde, active_subgraph.as_ref(), lmb)?)
        }
        ApproximationType::VaccuumLimit => Err(eyre!("Not yet implemented VaccuumLimit")),
        ApproximationType::OS => Err(eyre!("Not yet implemented OS")),
        ApproximationType::Unsubtracted => panic!("should have been kept out of the wood"),
    }
}

// #[debug_instrument(
//     current = %current.log_display(),
//     given = %given.log_display(),
//     reduced,
// )]
#[allow(clippy::too_many_arguments)]
fn t_tilde<S: ForestNodeLike>(
    ctx: &UVCtx<'_>,
    orientation: OrientationProjection<'_>,
    orientation_id: crate::cff::expression::OrientationID,
    source_edge_energy_map: Option<&[crate::cff::surface::LinearEnergyExpr]>,
    current: &S,
    given: &S,
    cff: &Atom,
    active_subgraph: Option<&SuBitGraph>,
    lmb: &LoopMomentumBasis,
) -> Result<Atom> {
    let graph = ctx.graph;
    let settings = ctx.settings;
    let reduced = current.reduced_subgraph(given);
    let rescaled_subgraph = active_subgraph.unwrap_or_else(|| current.subgraph());
    let lmb_id = lmb
        .loop_edges
        .first()
        .copied()
        .unwrap_or_else(|| current.lmb_id());

    // split numerator momenta into OSEs and spatial parts
    let mut reps = Vec::new();
    for (p, eid, e) in graph.iter_edges_of(rescaled_subgraph) {
        if p.is_paired() {
            let e_mass = e.data.mass_atom();
            reps.push(GS.split_mom_pattern(eid, lmb_id, e_mass, settings.inner_products));
        }
    }

    let numerator = graph
        .numerator(&reduced, given.subgraph())
        .get_single_atom()
        .unwrap();
    let mut numerator = orientation
        .map_numerator(graph, orientation_id, source_edge_energy_map, &numerator)?
        .replace_multiple(&reps);

    // rescale the external momenta in the added numerator subgraph
    for e in &lmb.ext_edges {
        // println!("Rescale {}", e);
        numerator = numerator
            .replace(GS.emr_vec_index(*e, W_.x___))
            .with(GS.emr_vec_index(*e, W_.x___) * GS.rescale)
            .replace(GS.emr_mom(*e, W_.x___))
            .with(GS.emr_mom(*e, W_.x___) * GS.rescale);
    }

    let mut atomarg = cff * numerator;

    // add data for OSE computation and add an explicit sqrt
    for (p, ei, e) in graph.iter_edges_of(rescaled_subgraph) {
        let eid = usize::from(ei) as i64;
        if p.is_paired() {
            // set energies from inner_t on-shell
            atomarg = atomarg.replace(function!(GS.energy, eid)).with(GS.ose(ei));

            let e_mass = e.data.mass_atom();
            atomarg = atomarg.replace(GS.ose(ei)).with(GS.ose_full(
                ei,
                lmb_id,
                e_mass,
                None,
                settings.inner_products,
            ));
        }
    }

    atomarg = atomarg.replace_multiple(&reps);

    let mom_reps = graph.replacement_impl(
        |e, loops, externals| {
            Replacement::new(
                GS.emr_vec
                    .call_args([Atom::num(usize::from(e)), Atom::var(W_.x___)])
                    .to_pattern(),
                (loops
                    .replace(function!(GS.emr_vec, W_.x_))
                    .allow_new_wildcards_on_rhs(true)
                    .with(
                        FunctionBuilder::new(GS.emr_vec)
                            .add_arg(W_.x_)
                            .add_args([W_.x___])
                            .finish(),
                    )
                    + externals * GS.rescale)
                    .to_pattern(),
            )
        },
        &reduced,
        lmb,
        GS.emr_vec,
        GS.emr_vec,
        &[],
        &[W_.x___],
        HedgePair::is_paired,
        true,
    );

    atomarg = atomarg.replace_multiple(&mom_reps);
    atomarg = atomarg
        .replace(function!(GS.ose, W_.a___))
        .with(function!(*OSE_FOR_LOCAL_3D_SERIES, W_.a___));

    let a = atomarg.series(GS.rescale, Atom::Zero, -1).unwrap();

    let mut a = a
        .to_atom()
        .replace(function!(
            Symbol::DERIVATIVE,
            0,
            1,
            *OSE_FOR_LOCAL_3D_SERIES,
            W_.y___
        ))
        .with(Atom::num(1))
        .replace(function!(
            Symbol::DERIVATIVE,
            W_.x___,
            *OSE_FOR_LOCAL_3D_SERIES,
            W_.y___
        ))
        .with(Atom::num(0));
    a = a
        .replace(function!(*OSE_FOR_LOCAL_3D_SERIES, W_.a___))
        .with(function!(GS.ose, W_.a___));
    a = a.replace(GS.rescale).with(Atom::num(1));
    Ok(a)
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
    let settings = ctx.settings;
    let rescaled_subgraph = active_subgraph.unwrap_or_else(|| current.subgraph());
    let lmb_id = lmb
        .loop_edges
        .first()
        .copied()
        .unwrap_or_else(|| current.lmb_id());
    let mut atomarg = cff * mapped_numerator;
    debug_tags!(#generation, #profile, #uv, #local, #trace;
        stage = "local_3d_start_initial",
        byte_size = atomarg.as_view().get_byte_size(),
        file.expr = %atomarg,
        "Local 3D start expression checkpoint"
    );
    // println!("CFF: {}", cff);

    // add data for OSE computation and add an explicit sqrt
    for (p, ei, e) in graph.iter_edges_of(rescaled_subgraph) {
        let eid = usize::from(ei) as i64;
        if p.is_paired() {
            // set energies from inner_t on-shell
            atomarg = atomarg.replace(function!(GS.energy, eid)).with(GS.ose(ei));

            let e_mass = e.data.mass_atom();
            atomarg = atomarg.replace(GS.ose(ei)).with(GS.ose_full(
                ei,
                lmb_id,
                e_mass,
                None,
                settings.inner_products,
            ));
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
            let e_mass = e.data.mass_atom();
            let rep = GS.split_mom_pattern(eid, lmb_id, e_mass, settings.inner_products);
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
    Ok(atomarg)
}

impl Local3DLoopRescaling {
    // #[debug_instrument(
    //     current = %current.log_display(),
    //     given = %given.log_display(),
    //     reduced,
    // )]
    fn t<S: ForestNodeLike>(
        self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrand: &Atom,
        _active_subgraph: Option<&SuBitGraph>,
        lmb: &LoopMomentumBasis,
    ) -> Result<Atom> {
        let graph = ctx.graph;
        let reduced = current.reduced_subgraph(given);

        // only apply replacements for edges in the reduced graph
        let mut mom_reps = graph.uv_spatial_wrapped_replacement(&reduced, lmb, &[W_.x___]);
        // Explicit spatial components such as `Q(e, cind(1))` remain in the
        // ordinary EMR wrapper after residue-key mapping; only tensor-vector
        // occurrences use `EMRvec`. Route both wrappers through the same LMB
        // before rescaling so every hard numerator momentum is scaled exactly
        // once, without expanding or redistributing the numerator.
        mom_reps.extend(graph.uv_wrapped_replacement(&reduced, lmb, &[W_.x___]));
        for m in &mom_reps {
            debug_tags!(#uv,#momentum,#trace;mom_rep=%m,"Mom rep");
        }

        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_input",
            byte_size = integrand.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );
        let mut atomarg = integrand.replace_multiple(&mom_reps);
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_after_momentum_replacements",
            byte_size = atomarg.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );

        // Rescale every loop momentum still active in this sector, including
        // cycles expanded by earlier local operations.
        for e in &lmb.loop_edges {
            // println!("Rescale {}", e);
            atomarg = atomarg
                .replace(GS.emr_vec_index(*e, W_.x___))
                .with(GS.emr_vec_index(*e, W_.x___) * GS.rescale)
                .replace(GS.emr_mom(*e, W_.x___))
                .with(GS.emr_mom(*e, W_.x___) * GS.rescale);
        }
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_after_loop_rescale",
            byte_size = atomarg.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );

        // (re-)expand OSEs from the subgraph only
        for eid in &lmb.loop_edges {
            let eid = eid.0 as i64;
            // rescale the whole OSE so that the function itself has no poles during the expansion
            atomarg = atomarg.replace(function!(GS.ose, eid, W_.prop_)).with(
                function!(
                    GS.ose,
                    eid,
                    (GS.m_uv_expansion * GS.m_uv_expansion * GS.rescale * GS.rescale + W_.prop_
                        - GS.m_uv_expansion * GS.m_uv_expansion)
                        / GS.rescale
                        / GS.rescale
                ) * GS.rescale
                    * GS.rescale,
            )
        }
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_after_ose_rescale",
            byte_size = atomarg.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );

        atomarg = (atomarg * self.measure_scaling(lmb))
            .replace(GS.rescale)
            .with(Atom::num(1) / GS.rescale);
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_before_series",
            loop_edges = ?lmb.loop_edges,
            byte_size = atomarg.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );

        debug_tags!(#uv, #local, #before_series; log.expr = atomarg, "Before series in t");

        let series = atomarg.series(GS.rescale, Atom::Zero, 0).unwrap();
        let series_atom = series.to_atom();
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_after_series",
            byte_size = series_atom.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );

        debug_tags!(#uv, #local; expr = %series, "After series in t");
        let a = series_atom.replace(GS.rescale).with(Atom::num(1));

        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_output",
            byte_size = a.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );
        debug_tags!(#uv, #local; log.expr = a, "Local 3D approximation");
        Ok(a)
    }

    fn measure_scaling(self, lmb: &LoopMomentumBasis) -> Atom {
        // The supplied LMB is the integration-space authority. In particular,
        // a remainder which is a tree in the original incidence can become a
        // loop after its frozen UV prefix is contracted.
        Atom::var(GS.rescale).pow(3 * lmb.loop_edges.len() as i64)
    }
}
