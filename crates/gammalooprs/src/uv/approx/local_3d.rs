use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike, SubSetOps};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder},
    function, parse,
};
use tracing::{info, instrument};

use crate::{
    graph::{Graph, LMBext, cuts::CutSet},
    utils::{GS, W_},
    uv::{
        RenormalizationScheme, UltravioletGraph,
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
        let cff = graph
            .cff(
                &to_contract
                    .union(&graph.tree_edges)
                    .subtract(&graph.initial_state_cut),
                cuts,
            )?
            .expression_with_selectors();

        let fourddenoms = GS.wrap_tree_denoms(
            graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
        );

        Ok(cff.iter().map(|a| a * &fourddenoms).collect())
    }

    pub(crate) fn root(graph: &mut Graph, cuts: &CutSet) -> Result<Vec<Atom>> {
        Self::dependent(graph, &graph.empty_subgraph::<SuBitGraph>(), cuts)
    }
}

impl ApproximationKernel<UVCtx<'_>> for Local3DApproximation {
    #[instrument(skip(self, ctx, current, given, cff))]
    fn kernel<S: super::ForestNodeLike>(
        &self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        cff: &Atom,
    ) -> Result<Atom> {
        let graph = ctx.graph;
        let settings = ctx.settings;
        let reduced = current.reduced_subgraph(given);
        let mut cff = cff.clone();

        // println!("CFF: {}", cff);

        // add data for OSE computation and add an explicit sqrt
        for (p, ei, e) in graph.iter_edges_of(current.subgraph()) {
            let eid = usize::from(ei) as i64;
            if p.is_paired() {
                // set energies from inner_t on-shell
                cff = cff.replace(function!(GS.energy, eid)).with(GS.ose(ei));

                let e_mass = e.data.mass_atom();
                cff = cff.replace(GS.ose(ei)).with(GS.ose_full(
                    ei,
                    e_mass,
                    None,
                    settings.inner_products,
                ));
            }
        }

        let mut atomarg = cff
            * graph
                .numerator(&reduced, given.subgraph())
                .get_single_atom()
                .unwrap();

        // println!(
        //     "Expand-prerep {} with dod={} in {:?}",
        //     atomarg, current.dod(), current.lmb().ext_edges
        // );

        // split numerator momenta into OSEs and spatial parts
        let mut reps = Vec::new();
        for (p, eid, e) in graph.iter_edges_of(current.subgraph()) {
            if p.is_paired() {
                let e_mass = e.data.mass_atom();
                reps.push(GS.split_mom_pattern(eid, e_mass, settings.inner_products));
            }
        }

        // println!("Split reps:");
        // for r in &reps {
        //     println!("{r}");
        // }

        atomarg = atomarg.replace_multiple(&reps);

        // only apply replacements for edges in the reduced graph
        let mom_reps = graph.uv_spatial_wrapped_replacement(&reduced, current.lmb(), &[W_.x___]);

        // println!(
        //     "Mom Reps : for {}",
        //     self.simple_approx
        //         .as_ref()
        //         .unwrap()
        //         .expr(&graph.full_filter())
        // );
        // for r in &mom_reps {
        //     println!("{r}");
        // }

        atomarg = atomarg.replace_multiple(&mom_reps);

        // debug!("Before rescaling loop momenta {}",);

        // rescale the loop momenta in the whole subgraph, including previously expanded cycles
        for e in &current.lmb().loop_edges {
            // println!("Rescale {}", e);
            atomarg = atomarg
                .replace(GS.emr_vec_index(*e, W_.x___))
                .with(GS.emr_vec_index(*e, W_.x___) * GS.rescale);
        }

        // println!(
        //     "Expand {} with dod={} in {:?}",
        //     atomarg, current.dod(), current.lmb().ext_edges
        // );

        let soft_ct = current.renormalization_scheme() == RenormalizationScheme::OS
            && graph.full_crown(current.subgraph()).n_included() == 2
            && current.dod() > 0
            && settings.softct;

        if soft_ct {
            info!(
                subgraph = %current.subgraph().string_label(),
                "OS local soft counterterm path is not implemented yet; using the standard UV expansion"
            );
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

        // atomarg = atomarg
        //     .replace(function!(MS.dot, GS.rescale * W_.x_, W_.y_))
        //     .repeat()
        //     .with(function!(MS.dot, W_.x_, W_.y_) * GS.rescale);

        atomarg = (atomarg
            * Atom::var(GS.rescale).pow(3 * graph.n_loops(current.subgraph()) as i64))
        .replace(GS.rescale)
        .with(Atom::num(1) / GS.rescale);
        // .replace(Atom::var(GS.rescale).pow(2).sqrt()) //.pow((1, 2)))
        // .with(GS.rescale)
        // .replace((Atom::var(GS.rescale).pow(-2) * W_.a___).pow((1, 2))) //.pow((1, 2)))
        // .repeat()
        // .with((Atom::var(W_.a___)).sqrt() / GS.rescale)
        // .replace((Atom::var(GS.rescale).pow(2) * W_.a___).pow((-1, 2))) //.pow((1, 2)))
        // .repeat()
        // .with((Atom::var(W_.a___)).pow((-1, 2)) / GS.rescale)
        // .replace((Atom::var(GS.rescale).pow(-2) * W_.a___).pow((-1, 2))) //.pow((1, 2)))
        // .repeat()
        // .with((Atom::var(W_.a___)).pow((-1, 2)) * GS.rescale)
        // .replace((Atom::var(GS.rescale).pow(2) * W_.a___).pow((1, 2))) //.pow((1, 2)))
        // .repeat()
        // .with((Atom::var(W_.a___)).pow((1, 2)) * GS.rescale);

        // println!("atomarg:{:>}", atomarg.expand());

        // println!("Series expanding {atomarg}");
        //

        let a = atomarg
            .series(GS.rescale, Atom::Zero, 0.into(), true)
            .unwrap();

        // debug!("Series: {}", a);

        let mut a = a
            .to_atom()
            .replace(parse!("der(0,0,0,1, OSE(y__))"))
            .with(Atom::num(1))
            .replace(parse!("der(x__, OSE(y__))"))
            .with(Atom::num(0));

        if soft_ct {
            let coeffs = a.coefficient_list::<u8>(&[Atom::var(GS.rescale)]);
            let mut b = Atom::Zero;
            let dod_pow = Atom::var(GS.rescale).pow(current.dod());
            for (pow, mut i) in coeffs {
                if pow == dod_pow {
                    // FIXME: how to do this _only_ for the subgraph masses? the numerator is still
                    // only that of the subgraph
                    // set the masses in the t=dod term to 0
                    // UV rearrange the denominators
                    /*for m in &masses {
                        i = i.replace(m.clone()).with(Atom::Zero);
                    }*/

                    i = i
                        .replace(parse!("OSE(n_,q_,mass_,prop_, x___)"))
                        .with(parse!("OSE(n_,q_,mUV^2,prop_+mUV^2,x___)"));
                }

                b += i;
            }

            a = b;
        } else {
            a = a.replace(GS.rescale).with(Atom::num(1));
        }

        // println!("Expanded: {:>}", a.expand());

        Ok(a)
    }
}
