use eyre::eyre;
use linnet::half_edge::{involution::HedgePair, subgraph::SuBitGraph};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder},
    function,
    id::Replacement,
    parse,
};
use tracing::{debug, instrument};

use crate::{
    graph::{Graph, LMBext, cuts::CutSet},
    utils::{GS, W_, symbolica_ext::CallSymbol},
    uv::{
        ApproximationType, UltravioletGraph,
        approx::{ApproximationKernel, UVCtx},
        uv_graph::UVE,
    },
};
use color_eyre::Result;

pub struct Local3DApproximation;

impl Local3DApproximation {
    pub(crate) fn dependent(
        graph: &mut Graph,
        to_contract: &SuBitGraph,
        cuts: &CutSet,
    ) -> Result<Vec<Atom>> {
        let cff = graph.cff(to_contract, cuts)?.expression_with_selectors();

        Ok(cff)
    }

    pub(crate) fn root(graph: &mut Graph, cuts: &CutSet) -> Result<Vec<Atom>> {
        Self::dependent(graph, &graph.empty_subgraph::<SuBitGraph>(), cuts)
    }

    pub(crate) fn t_tilde<S: super::ForestNodeLike>(
        &self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        cff: &Atom,
    ) -> Result<Atom> {
        let graph = ctx.graph;
        let settings = ctx.settings;
        let reduced = current.reduced_subgraph(given);

        // split numerator momenta into OSEs and spatial parts
        let mut reps = Vec::new();
        for (p, eid, e) in graph.iter_edges_of(current.subgraph()) {
            if p.is_paired() {
                let e_mass = e.data.mass_atom();
                reps.push(GS.split_mom_pattern(eid, e_mass, settings.inner_products));
            }
        }

        let mut numerator = graph
            .numerator(&reduced, given.subgraph())
            .get_single_atom()
            .unwrap()
            .replace_multiple(&reps);

        // rescale the external momenta in the added numerator subgraph
        for e in &current.lmb().ext_edges {
            // println!("Rescale {}", e);
            numerator = numerator
                .replace(GS.emr_vec_index(*e, W_.x___))
                .with(GS.emr_vec_index(*e, W_.x___) * GS.rescale);
        }

        let mut atomarg = cff * numerator;
        // println!("CFF: {}", cff);

        // add data for OSE computation and add an explicit sqrt
        for (p, ei, e) in graph.iter_edges_of(current.subgraph()) {
            let eid = usize::from(ei) as i64;
            if p.is_paired() {
                // set energies from inner_t on-shell
                atomarg = atomarg.replace(function!(GS.energy, eid)).with(GS.ose(ei));

                let e_mass = e.data.mass_atom();
                atomarg = atomarg.replace(GS.ose(ei)).with(GS.ose_full(
                    ei,
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
                        .f([Atom::num(usize::from(e)), Atom::var(W_.x___)])
                        .to_pattern(),
                    (loops
                        .replace(function!(GS.emr_vec, W_.x_))
                        .allow_new_wildcards_on_rhs(true)
                        .with(
                            FunctionBuilder::new(GS.emr_vec)
                                .add_arg(W_.x_)
                                .add_args(&[W_.x___])
                                .finish(),
                        )
                        + externals * GS.rescale)
                        .to_pattern(),
                )
            },
            &reduced,
            current.lmb(),
            GS.emr_vec,
            GS.emr_vec,
            &[],
            &[W_.x___],
            HedgePair::is_paired,
            true,
        );

        atomarg = atomarg.replace_multiple(&mom_reps);
        let a = atomarg
            .series(GS.rescale, Atom::Zero, (-1).into(), true)
            .unwrap();

        let mut a = a
            .to_atom()
            .replace(parse!("der(0,0,0,1, OSE(y__))"))
            .with(Atom::num(1))
            .replace(parse!("der(x__, OSE(y__))"))
            .with(Atom::num(0));
        a = a.replace(GS.rescale).with(Atom::num(1));
        Ok(a)
    }

    pub(crate) fn t<S: super::ForestNodeLike>(
        &self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrand: &Atom,
    ) -> Result<Atom> {
        let graph = ctx.graph;
        let reduced = current.reduced_subgraph(given);

        // only apply replacements for edges in the reduced graph
        let mom_reps = graph.uv_spatial_wrapped_replacement(&reduced, current.lmb(), &[W_.x___]);

        let mut atomarg = integrand.replace_multiple(&mom_reps);

        // rescale the loop momenta in the whole subgraph, including previously expanded cycles
        for e in &current.lmb().loop_edges {
            // println!("Rescale {}", e);
            atomarg = atomarg
                .replace(GS.emr_vec_index(*e, W_.x___))
                .with(GS.emr_vec_index(*e, W_.x___) * GS.rescale);
        }

        // (re-)expand OSEs from the subgraph only
        for (_, eid, _) in graph.iter_edges_of(current.subgraph()) {
            let eid = usize::from(eid) as i64;
            // rescale the whole OSE so that the function itself has no poles during the expansion
            atomarg = atomarg
                .replace(function!(GS.ose, eid, W_.mom_, W_.mass_, W_.prop_))
                .with(
                    function!(
                        GS.ose,
                        eid,
                        W_.mom_,
                        GS.m_uv * GS.m_uv,
                        (GS.m_uv * GS.m_uv * GS.rescale * GS.rescale + W_.prop_
                            - GS.m_uv * GS.m_uv)
                            / GS.rescale
                            / GS.rescale
                    ) * GS.rescale
                        * GS.rescale,
                )
                .replace(function!(GS.ose, eid, W_.mom_, W_.a___)) //rescale the momenta for the same reason
                .with_map(move |m| {
                    let mut f = FunctionBuilder::new(GS.ose);
                    f = f.add_arg(eid);
                    f = f.add_arg(
                        (m.get(W_.mom_)
                            .unwrap()
                            .to_atom()
                            .replace(GS.rescale)
                            .with(Atom::num(1) / GS.rescale)
                            * GS.rescale)
                            .expand()
                            .replace(GS.rescale)
                            .with(Atom::Zero),
                    );
                    f = f.add_arg(m.get(W_.a___).unwrap().to_atom());

                    f.finish()
                });
        }

        atomarg = (atomarg
            * Atom::var(GS.rescale).pow(3 * graph.n_loops(current.subgraph()) as i64))
        .replace(GS.rescale)
        .with(Atom::num(1) / GS.rescale);
        let a = atomarg
            .series(GS.rescale, Atom::Zero, 0.into(), true)
            .unwrap();

        let mut a = a
            .to_atom()
            .replace(parse!("der(0,0,0,1, OSE(y__))"))
            .with(Atom::num(1))
            .replace(parse!("der(x__, OSE(y__))"))
            .with(Atom::num(0));
        a = a.replace(GS.rescale).with(Atom::num(1));
        debug!("a: {}", a);
        Ok(a)
    }

    pub(crate) fn start<S: super::ForestNodeLike>(
        &self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        cff: &Atom,
    ) -> Result<Atom> {
        let graph = ctx.graph;
        let settings = ctx.settings;
        let reduced = current.reduced_subgraph(given);
        let mut atomarg = cff
            * graph
                .numerator(&reduced, given.subgraph())
                .get_single_atom()
                .unwrap();
        // println!("CFF: {}", cff);

        // add data for OSE computation and add an explicit sqrt
        for (p, ei, e) in graph.iter_edges_of(current.subgraph()) {
            let eid = usize::from(ei) as i64;
            if p.is_paired() {
                // set energies from inner_t on-shell
                atomarg = atomarg.replace(function!(GS.energy, eid)).with(GS.ose(ei));

                let e_mass = e.data.mass_atom();
                atomarg = atomarg.replace(GS.ose(ei)).with(GS.ose_full(
                    ei,
                    e_mass,
                    None,
                    settings.inner_products,
                ));
            }
        }

        // split numerator momenta into OSEs and spatial parts
        let mut reps = Vec::new();
        for (p, eid, e) in graph.iter_edges_of(current.subgraph()) {
            if p.is_paired() {
                let e_mass = e.data.mass_atom();
                reps.push(GS.split_mom_pattern(eid, e_mass, settings.inner_products));
            }
        }
        Ok(atomarg.replace_multiple(&reps))
    }
}

impl ApproximationKernel<UVCtx<'_>> for Local3DApproximation {
    #[instrument(skip(self, ctx, current, given, integrand))]
    fn kernel<S: super::ForestNodeLike>(
        &self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrand: &Atom,
    ) -> Result<Atom> {
        match current.renormalization_scheme() {
            ApproximationType::MUV => self.t(
                ctx,
                current,
                given,
                &self.start(ctx, current, given, integrand)?,
            ),
            ApproximationType::IR => Ok(self.t(
                ctx,
                current,
                given,
                &self.start(ctx, current, given, integrand)?,
            )? + self.t_tilde(ctx, current, given, integrand)?
                - self.t(
                    ctx,
                    current,
                    given,
                    &self.t_tilde(ctx, current, given, integrand)?,
                )?),
            ApproximationType::VaccuumLimit => Err(eyre!("Not yet implemented VaccuumLimit")),
            ApproximationType::OS => Err(eyre!("Not yet implemented OS")),
            ApproximationType::Unsubtracted => {
                panic!("should have been kept out of the wood");
            }
        }
    }
}
