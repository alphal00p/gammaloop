use std::{ops::Neg, sync::LazyLock};

use eyre::eyre;

use linnet::half_edge::{
    involution::{EdgeIndex, EdgeVec, HedgePair, Orientation},
    subgraph::{SuBitGraph, SubSetOps},
};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder, Symbol},
    function,
    id::Replacement,
    symbol,
};

use crate::{
    cff::{CutCFF, orientations::GraphOrientation},
    debug_tags,
    graph::{Graph, LMBext, LoopMomentumBasis, cuts::CutSet},
    utils::{GS, W_},
    uv::{
        ApproximationType, Integrands, UVgenerationSettings, UltravioletGraph,
        approx::{ForestNodeLike, OrientationProjection, UVCtx, integrated::IntegratedCts},
        marker::{UvMarker, UvOperation},
        uv_graph::UVE,
    },
};
use color_eyre::Result;

pub(crate) struct FrozenActiveCt {
    pub active: Integrands,
    pub frozen_integrands: Integrands,
}

impl FrozenActiveCt {
    pub fn combine(self) -> Result<Integrands> {
        self.active.zip_mul(&self.frozen_integrands)
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Localizer<'a> {
    cutset: &'a CutSet,
    orientation: OrientationProjection<'a>,
}

impl<'a> Localizer<'a> {
    pub fn cff(self, graph: &mut Graph, to_contract: &SuBitGraph) -> Result<CutCFF> {
        graph.cff(
            &to_contract
                .union(&graph.tree_edges)
                .subtract(&graph.initial_state_cut),
            self.cutset,
            self.orientation.orientation_pattern,
        )
    }

    pub fn new(cutset: &'a CutSet, orientation: OrientationProjection<'a>) -> Self {
        Self {
            cutset,
            orientation,
        }
    }
    fn representative_orientations(
        &self,
        reduced_orientation: &EdgeVec<Orientation>,
        internal_edges: &[EdgeIndex],
        average: bool,
    ) -> Result<Vec<&EdgeVec<Orientation>>> {
        let iterator = self
            .orientation
            .valid_orientations
            .iter()
            .filter(|orientation| orientation.is_compatible_with(reduced_orientation));

        let out: Vec<&EdgeVec<_>> = if average {
            iterator.collect()
        } else {
            let max: Option<&EdgeVec<Orientation>> =
                iterator.max_by_key(|orientation| orientation.score(internal_edges));
            max.into_iter().collect()
        };

        if out.is_empty() {
            Err(eyre!(
                "no valid global orientation matches reduced orientation {}",
                GS.orientation_delta(reduced_orientation)
            ))
        } else {
            Ok(out)
        }
    }

    /// Localize a reduced-graph orientation term onto compatible full orientations.
    ///
    /// The reduced term already contains the external-edge selector structure through
    /// `reduced_orientation.orientation_thetas()`. We only add selectors for the integrated
    /// subgraph's internal edges.
    fn localized_orientation_term(
        self,
        reduced_expression: &Atom,
        reduced_orientation: &EdgeVec<Orientation>,
        internal_edges: &[EdgeIndex],
        average_compatible_orientations: bool,
    ) -> Result<Atom> {
        let representatives = self.representative_orientations(
            reduced_orientation,
            internal_edges,
            average_compatible_orientations,
        )?;
        let reduced_selector = reduced_orientation.orientation_thetas();
        let mut localized = Atom::Zero;

        let n_reps = representatives.len();
        for representative in representatives {
            let internal_selector = representative.internal_orientation_selector(internal_edges);
            localized += reduced_expression.clone() * &reduced_selector * internal_selector;
        }
        if average_compatible_orientations {
            localized /= n_reps;
        }

        Ok(localized)
    }
    pub fn localize<S: ForestNodeLike>(
        &self,
        expr: &Atom,
        graph: &mut Graph,
        integrated_node: &S,
    ) -> Result<FrozenActiveCt> {
        if expr.is_zero() {
            let indices = self.cutset.residue_selector.generate_allowed_keys();
            return Ok(FrozenActiveCt {
                active: indices.iter().map(|index| (*index, Atom::Zero)).collect(),
                frozen_integrands: indices
                    .into_iter()
                    .map(|index| (index, Atom::one()))
                    .collect(),
            });
        }

        let to_contract = integrated_node.subgraph();
        let integrated_loop_count = graph.n_loops(to_contract);
        // Keep finite addbacks for nested multi-loop entries. Integrated expressions carry
        // their forest-composition signs where needed, so dropping the localized finite
        // representative here removes the Tint(T(...)) terms.
        let finite_ct = expr.clone();
        debug_tags!(#generation, #profile, #uv, #integrated, #local, #summary;
            stage = "localize_integrated_ct_forest_overlap",
            integrated_loop_count,
            "Applied integrated UV forest-overlap addback rule"
        );
        let cff = graph.cff(
            &to_contract
                .union(&graph.tree_edges)
                .subtract(&graph.initial_state_cut),
            self.cutset,
            self.orientation.orientation_pattern,
        )?;

        let fourddenoms = GS.wrap_tree_denoms(
            graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
        );

        let internal_edges = graph.paired_edges(to_contract);
        let localizing_integrand = GS.localizing_integrand(integrated_node.lmb());

        let mut active_integrands = Vec::new();
        let mut frozen_integrands = Vec::new();
        for (index, term) in cff.terms {
            let mut localized = Atom::Zero;
            for (cff_term, orientation) in term.expression.into_iter().zip(term.orientations) {
                localized += self.localized_orientation_term(
                    &(cff_term * &fourddenoms),
                    &orientation,
                    &internal_edges,
                    true,
                )?;
            }

            let localized_cff_byte_size = localized.as_view().get_byte_size();
            let active_ct = localized * finite_ct.clone();
            let localized_ct = &active_ct * &localizing_integrand;
            debug_tags!(#generation, #profile, #uv, #integrated, #local, #term, #summary;
                stage = "localize_integrated_ct_term",
                cut_index = ?index,
                localized_cff_byte_size,
                active_ct_byte_size = active_ct.as_view().get_byte_size(),
                localized_ct_byte_size = localized_ct.as_view().get_byte_size(),
                "Integrated UV CT localization size checkpoint"
            );

            active_integrands.push((index, active_ct));
            frozen_integrands.push((index, localizing_integrand.clone()));
        }

        Ok(FrozenActiveCt {
            active: Integrands::from_iter(active_integrands),
            frozen_integrands: Integrands::from_iter(frozen_integrands),
        })
    }
}

pub struct Local3DApproximation<'a> {
    localizer: Localizer<'a>,
    graph: &'a mut Graph,
    settings: &'a UVgenerationSettings,
}

impl<'a> Local3DApproximation<'a> {
    pub(crate) fn new(
        localizer: Localizer<'a>,
        graph: &'a mut Graph,
        settings: &'a UVgenerationSettings,
    ) -> Self {
        Self {
            localizer,
            graph,
            settings,
        }
    }

    pub(crate) fn run<S: ForestNodeLike, M: ForestNodeLike>(
        self,
        local: &Local3DCts,
        integrated: &IntegratedCts,
        current: &S,
        given: &S,
        marker_current: &M,
        marker_given: &M,
    ) -> Result<Local3DCts> {
        let integrated_t = self.localizer.localize(
            &integrated.physical_finite_counterterm_atom(),
            self.graph,
            given,
        )?;
        let ctx = UVCtx::new(self.graph, self.settings);
        let marker = UvMarker::new(ctx.settings);

        if let Some(active_sectors) = local.active_sectors() {
            let reduced_subgraph = current.reduced_subgraph(given);
            let mut next_sectors = Vec::with_capacity(active_sectors.len() + 1);

            for (active_subgraph, integrands) in active_sectors {
                let active_subgraph = active_subgraph.union(&reduced_subgraph);
                let integrands =
                    -integrands.fallible_map(Local3DLoopRescaling::FullSubgraph.map(
                        &ctx,
                        current,
                        given,
                        Some(active_subgraph.clone()),
                    ))?;
                next_sectors.push((active_subgraph, integrands));
            }

            let integrated = -(integrated_t
                .active
                .fallible_map(Local3DLoopRescaling::ReducedSubgraph.map(
                    &ctx,
                    current,
                    given,
                    Some(reduced_subgraph.clone()),
                ))?
                .zip_mul(&integrated_t.frozen_integrands)?);
            next_sectors.push((reduced_subgraph, integrated));

            return Local3DCts::from_active_sectors(next_sectors)?.map(|atom| {
                Ok(marker.apply(
                    UvOperation::Approx,
                    marker_current.subgraph(),
                    marker_given.subgraph(),
                    atom,
                ))
            });
        }

        let local = -(local.map(full(&ctx, current, given))?);
        let integrated = -(integrated_t
            .active
            .fallible_map(reduced(&ctx, current, given))?
            .zip_mul(&integrated_t.frozen_integrands)?);

        local.zip_add(&integrated)?.map(|atom| {
            Ok(marker.apply(
                UvOperation::Approx,
                marker_current.subgraph(),
                marker_given.subgraph(),
                atom,
            ))
        })
    }

    pub(crate) fn run_local<S: ForestNodeLike, M: ForestNodeLike>(
        self,
        local: &Local3DCts,
        current: &S,
        given: &S,
        marker_current: &M,
        marker_given: &M,
    ) -> Result<Local3DCts> {
        let ctx = UVCtx::new(self.graph, self.settings);
        let reduced_subgraph = current.reduced_subgraph(given);
        let active_sectors = match local.active_sectors() {
            Some(active_sectors) => active_sectors
                .iter()
                .map(|(active_subgraph, integrands)| {
                    let active_subgraph = active_subgraph.union(&reduced_subgraph);
                    // Retain the full mask for descendants, but this replay step acts only
                    // on the part of that mask covered by its component-local path.
                    let rescaled_subgraph = active_subgraph.intersection(current.subgraph());
                    Ok((
                        active_subgraph,
                        -integrands.fallible_map(Local3DLoopRescaling::FullSubgraph.map(
                            &ctx,
                            current,
                            given,
                            Some(rescaled_subgraph),
                        ))?,
                    ))
                })
                .collect::<Result<Vec<_>>>()?,
            None => vec![(
                reduced_subgraph.clone(),
                -local
                    .integrands()
                    .fallible_map(Local3DLoopRescaling::FullSubgraph.map(
                        &ctx,
                        current,
                        given,
                        Some(reduced_subgraph),
                    ))?,
            )],
        };
        Local3DCts::from_active_sectors(active_sectors)?.map(|atom| {
            Ok(UvMarker::new(ctx.settings).apply(
                UvOperation::Approx,
                marker_current.subgraph(),
                marker_given.subgraph(),
                atom,
            ))
        })
    }

    pub(crate) fn run_integrated<S: ForestNodeLike, I: ForestNodeLike, M: ForestNodeLike>(
        self,
        integrated: &IntegratedCts,
        integrated_node: &I,
        current: &S,
        given: &S,
        marker_current: &M,
        marker_given: &M,
    ) -> Result<Local3DCts> {
        let integrated = self.localizer.localize(
            &integrated.physical_finite_counterterm_atom(),
            self.graph,
            integrated_node,
        )?;
        let ctx = UVCtx::new(self.graph, self.settings);
        let active_subgraph = current.reduced_subgraph(given);
        let integrated = -(integrated
            .active
            .fallible_map(Local3DLoopRescaling::ReducedSubgraph.map(
                &ctx,
                current,
                given,
                Some(active_subgraph.clone()),
            ))?
            .zip_mul(&integrated.frozen_integrands)?);

        Local3DCts::from_active_sectors(vec![(active_subgraph, integrated)])?.map(|atom| {
            Ok(UvMarker::new(ctx.settings).apply(
                UvOperation::Approx,
                marker_current.subgraph(),
                marker_given.subgraph(),
                atom,
            ))
        })
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Local3DCts {
    integrands: Integrands,
    // A disconnected join can leave a different set of loop variables active
    // in each local/integrated cross term.
    active_sectors: Option<Vec<(SuBitGraph, Integrands)>>,
}

impl From<Integrands> for Local3DCts {
    fn from(integrands: Integrands) -> Self {
        Self {
            integrands,
            active_sectors: None,
        }
    }
}

impl Neg for Local3DCts {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            integrands: self.integrands.map(|a| a.neg()),
            active_sectors: self.active_sectors.map(|sectors| {
                sectors
                    .into_iter()
                    .map(|(active, integrands)| (active, -integrands))
                    .collect()
            }),
        }
    }
}

impl Local3DCts {
    pub fn zip_add(&self, other: &Integrands) -> Result<Self> {
        if self.active_sectors.is_some() {
            return Err(eyre!(
                "an unlabelled local term cannot be added to active UV sectors"
            ));
        }
        Ok(Self::from(self.integrands.zip_add(other)?))
    }

    pub(crate) fn integrands(&self) -> &Integrands {
        &self.integrands
    }

    pub(crate) fn active_sectors(&self) -> Option<&[(SuBitGraph, Integrands)]> {
        self.active_sectors.as_deref()
    }

    pub(crate) fn from_active_sectors(
        active_sectors: Vec<(SuBitGraph, Integrands)>,
    ) -> Result<Self> {
        let mut sectors = active_sectors.iter();
        let mut integrands = sectors
            .next()
            .ok_or_else(|| eyre!("active UV counterterm sectors cannot be empty"))?
            .1
            .clone();
        for (_, sector) in sectors {
            integrands = integrands.zip_add(sector)?;
        }

        Ok(Self {
            integrands,
            active_sectors: Some(active_sectors),
        })
    }

    pub fn map<F: FnMut(&Atom) -> Result<Atom>>(&self, mut f: F) -> Result<Self> {
        if let Some(active_sectors) = &self.active_sectors {
            let active_sectors = active_sectors
                .iter()
                .map(|(active, integrands)| Ok((active.clone(), integrands.fallible_map(&mut f)?)))
                .collect::<Result<_>>()?;
            Self::from_active_sectors(active_sectors)
        } else {
            self.integrands.fallible_map(f).map(Self::from)
        }
    }

    pub(crate) fn root(graph: &mut Graph, localizer: Localizer<'_>) -> Result<Self> {
        let cff = localizer
            .cff(graph, &graph.empty_subgraph::<SuBitGraph>())?
            .expression_with_selectors();

        let fourddenoms = GS.wrap_tree_denoms(
            graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
        );

        Ok(Local3DCts::from(cff * fourddenoms))
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum Local3DLoopRescaling {
    FullSubgraph,
    ReducedSubgraph,
}

static OSE_FOR_LOCAL_3D_SERIES: LazyLock<Symbol> = LazyLock::new(|| {
    symbol!(
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

// #[debug_instrument(
//     current = %current.log_display(),
//     given = %given.log_display(),
//     reduced,
// )]
pub(crate) fn t_tilde<S: super::ForestNodeLike>(
    ctx: &UVCtx<'_>,
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

    let mut numerator = graph
        .numerator(&reduced, given.subgraph())
        .get_single_atom()
        .unwrap()
        .replace_multiple(&reps);

    // rescale the external momenta in the added numerator subgraph
    for e in &lmb.ext_edges {
        // println!("Rescale {}", e);
        numerator = numerator
            .replace(GS.emr_vec_index(*e, W_.x___))
            .with(GS.emr_vec_index(*e, W_.x___) * GS.rescale);
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

pub(crate) fn start<S: super::ForestNodeLike>(
    ctx: &UVCtx<'_>,
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
    let mut atomarg = cff
        * graph
            .numerator(&reduced, given.subgraph())
            .get_single_atom()
            .unwrap();
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

fn full<S: super::ForestNodeLike>(
    ctx: &UVCtx<'_>,
    current: &S,
    given: &S,
) -> impl FnMut(&Atom) -> Result<Atom> {
    Local3DLoopRescaling::FullSubgraph.map(ctx, current, given, None)
}
fn reduced<S: super::ForestNodeLike>(
    ctx: &UVCtx<'_>,
    current: &S,
    given: &S,
) -> impl FnMut(&Atom) -> Result<Atom> {
    Local3DLoopRescaling::ReducedSubgraph.map(ctx, current, given, None)
}
impl Local3DLoopRescaling {
    // #[debug_instrument(
    //     current = %current.log_display(),
    //     given = %given.log_display(),
    //     reduced,
    // )]
    pub(crate) fn t<S: super::ForestNodeLike>(
        self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrand: &Atom,
        active_subgraph: Option<&SuBitGraph>,
        lmb: &LoopMomentumBasis,
    ) -> Result<Atom> {
        let graph = ctx.graph;
        let reduced = current.reduced_subgraph(given);

        // only apply replacements for edges in the reduced graph
        let mom_reps = graph.uv_spatial_wrapped_replacement(&reduced, lmb, &[W_.x___]);
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
                .with(GS.emr_vec_index(*e, W_.x___) * GS.rescale);
        }
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_after_loop_rescale",
            byte_size = atomarg.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );

        // (re-)expand OSEs from the subgraph only
        for eid in lmb.loop_edges.iter() {
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

        atomarg = (atomarg * self.measure_scaling(ctx, current, given, active_subgraph))
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

    // #[instrument(skip(self, ctx, current, given))]
    fn map<S: super::ForestNodeLike>(
        self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        active_subgraph: Option<SuBitGraph>,
    ) -> impl FnMut(&Atom) -> Result<Atom> {
        move |integrand| {
            let active_lmb = active_subgraph
                .as_ref()
                .filter(|active| *active != current.subgraph())
                .map(|active| {
                    ctx.graph.try_compatible_sub_lmb(
                        active,
                        ctx.graph.dummy_less_full_crown(active),
                        current.lmb(),
                    )
                })
                .transpose()?;
            let lmb = active_lmb.as_ref().unwrap_or_else(|| current.lmb());

            match current.renormalization_scheme() {
                ApproximationType::MUV | ApproximationType::PolePart => {
                    let started = start(
                        ctx,
                        current,
                        given,
                        integrand,
                        active_subgraph.as_ref(),
                        lmb,
                    )?;
                    crate::debug_tags!(#generation, #profile, #uv, #local, #summary;
                        stage = "local_3d_kernel_after_start",
                        input_byte_size = integrand.as_view().get_byte_size(),
                        output_byte_size = started.as_view().get_byte_size(),
                        "Local 3D kernel size checkpoint"
                    );
                    self.t(ctx, current, given, &started, active_subgraph.as_ref(), lmb)
                }
                ApproximationType::IR => {
                    let t_tilde = t_tilde(
                        ctx,
                        current,
                        given,
                        integrand,
                        active_subgraph.as_ref(),
                        lmb,
                    )?;
                    Ok(self.t(
                        ctx,
                        current,
                        given,
                        &start(
                            ctx,
                            current,
                            given,
                            integrand,
                            active_subgraph.as_ref(),
                            lmb,
                        )?,
                        active_subgraph.as_ref(),
                        lmb,
                    )? + &t_tilde
                        - self.t(ctx, current, given, &t_tilde, active_subgraph.as_ref(), lmb)?)
                }
                ApproximationType::VaccuumLimit => Err(eyre!("Not yet implemented VaccuumLimit")),
                ApproximationType::OS => Err(eyre!("Not yet implemented OS")),
                ApproximationType::Unsubtracted => {
                    panic!("should have been kept out of the wood");
                }
            }
        }
    }

    fn measure_scaling<S: super::ForestNodeLike>(
        self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        active_subgraph: Option<&SuBitGraph>,
    ) -> Atom {
        let n_rescaled_loops = match active_subgraph {
            Some(active_subgraph) => ctx.graph.n_loops(active_subgraph),
            None => match self {
                Local3DLoopRescaling::FullSubgraph => ctx.graph.n_loops(current.subgraph()),
                Local3DLoopRescaling::ReducedSubgraph => {
                    ctx.graph.n_loops(current.subgraph()) - ctx.graph.n_loops(given.subgraph())
                }
            },
        };

        Atom::var(GS.rescale).pow(3 * n_rescaled_loops as i64)
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        graph::cuts::CutSet,
        settings::global::OrientationPattern,
        utils::GS,
        uv::approx::{OrientationProjection, local_3d::Localizer},
    };
    use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
    use symbolica::function;

    fn orientation(value: i8) -> Orientation {
        match value {
            1 => Orientation::Default,
            -1 => Orientation::Reversed,
            0 => Orientation::Undirected,
            _ => panic!("invalid orientation encoding"),
        }
    }

    fn edgevec(values: impl IntoIterator<Item = i8>) -> EdgeVec<Orientation> {
        EdgeVec::from_iter(values.into_iter().map(orientation))
    }

    fn edges(values: impl IntoIterator<Item = usize>) -> Vec<EdgeIndex> {
        values.into_iter().map(EdgeIndex).collect()
    }

    #[test]
    fn orientation_term_keeps_external_selectors_and_adds_internal_ones() {
        let reduced_expression = function!(GS.ose, 0);
        let reduced_orientation = edgevec([1, 0, -1]);
        let valid = vec![edgevec([1, 1, -1]), edgevec([1, -1, -1])];
        let pat = OrientationPattern::default();
        let cutset = CutSet::empty(1);
        let loc = Localizer::new(&cutset, OrientationProjection::new(&valid, &pat));
        let localized = loc
            .localized_orientation_term(
                &reduced_expression,
                &reduced_orientation,
                &edges([1]),
                false,
            )
            .unwrap();

        let expected = reduced_expression
            * GS.sign_theta(GS.sign(EdgeIndex(0)))
            * GS.sign_theta(-GS.sign(EdgeIndex(2)))
            * GS.sign_theta(GS.sign(EdgeIndex(1)));
        assert_eq!(localized, expected);
    }
}
