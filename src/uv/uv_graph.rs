use ahash::AHashSet;
use bitvec::vec::BitVec;
use idenso::metric::MS;
use linnet::half_edge::{
    involution::HedgePair,
    subgraph::{Cycle, InternalSubGraph, SubGraph, SubGraphOps},
    HedgeGraph, PowersetIterator,
};
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    function,
    id::Replacement,
    parse,
};

use crate::{
    momentum_sample::LoopIndex,
    new_graph::{Edge, Graph, LMBext, LoopMomentumBasis, NumHedgeData, Vertex},
    numerator::{GlobalPrefactor, Numerator},
    symbolica_ext::CallSymbol,
    utils::{GS, W_},
};

use super::{is_not_paired, spenso_lor_atom, Wood};

pub trait UltravioletGraph: LMBext {
    fn n_loops<S: SubGraph, E, V, H>(&self, subgraph: &S) -> usize
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        self.as_ref().cyclotomatic_number(subgraph)
    }
    fn numerator<S: SubGraph>(&self, subgraph: &S) -> Atom;
    fn denominator<S: SubGraph>(&self, subgraph: &S) -> Atom;
    fn all_cycle_unions<E, V, H, S: SubGraph<Base = BitVec>>(
        &self,
        subgraph: &S,
    ) -> AHashSet<InternalSubGraph>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        let ref_graph = self.as_ref();
        let init_node = ref_graph.iter_nodes_of(subgraph).next().unwrap().0;
        let all_subcycles: Vec<_> = Cycle::all_sum_powerset_filter_map(
            &ref_graph
                .paton_cycle_basis(subgraph, &init_node, None)
                .unwrap()
                .0,
            &Some,
        )
        .map(|a| a.into_iter().map(|c| c.internal_graph(ref_graph)).collect())
        .unwrap();

        // println!("{}", self.base_dot());
        let spinneys: AHashSet<_> = InternalSubGraph::all_unions_iterative(&all_subcycles);

        spinneys
    }
    fn all_limits<E, V, H, S: SubGraph>(
        &self,
        subgraph: &S,
        expr: &Atom,
        expansion: Symbol,
        lmb: &LoopMomentumBasis,
    ) -> Vec<Atom>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        let mom_reps = self.uv_spatial_wrapped_replacement(subgraph, lmb, &[W_.x___]);

        // for x in &mom_reps {
        //     println!("LMB replacement: {x}");
        // }

        let energy_reps = self.replacement_impl::<_, Atom>(
            |e, a, b| {
                Replacement::new(
                    GS.energy.f([usize::from(e) as i32]).to_pattern(),
                    (a + b).to_pattern(),
                )
            },
            subgraph,
            lmb,
            GS.energy,
            GS.energy,
            &[],
            &[],
            is_not_paired,
            true,
        );

        // for e in &energy_reps {
        //     println!("Energy replacement: {e}");
        // }

        let q3_reps = self.replacement_impl::<_, Atom>(
            |e, a, b| {
                Replacement::new(
                    GS.emr_vec.f([usize::from(e) as i32]).to_pattern(),
                    (a + b).to_pattern(),
                )
            },
            subgraph,
            lmb,
            GS.emr_vec,
            GS.emr_vec,
            &[],
            &[],
            is_not_paired,
            true,
        );

        // for e in &q3_reps {
        //     println!("Q3 replacement: {e}");
        // }

        let expr = expr
            .replace_multiple(&mom_reps)
            .replace_multiple(&energy_reps)
            .replace_multiple(&q3_reps);
        let mut loops = PowersetIterator::new(lmb.loop_edges.len() as u8).into_iter();

        let mut limits = Vec::new();

        loops.next();

        for ls in loops {
            let mut expr = expr.clone();
            for l in ls.iter_ones() {
                let e = usize::from(lmb.loop_edges[LoopIndex(l)]) as i64;
                expr = expr
                    .replace(function!(GS.emr_vec, e, W_.x___))
                    .with(function!(GS.emr_vec, e, W_.x___) / expansion);

                expr /= Atom::var(expansion).npow(3);
            }

            expr = expr
                .replace(function!(MS.dot, W_.x___))
                .with(-function!(MS.dot, W_.x___)) // make dot products positive
                .replace(function!(MS.dot, W_.x_ / expansion, W_.y_))
                .repeat()
                .with(function!(MS.dot, W_.x_, W_.y_) / expansion);

            let series = expr.series(expansion, Atom::Zero, 0.into(), true).unwrap();

            expr = series.to_atom().expand();

            let l = expr.coefficient_list::<i8>(&[Atom::var(expansion)]);

            println!(
                "LIMIT {:?}:",
                ls.iter_ones()
                    .map(|l| usize::from(lmb.loop_edges[LoopIndex(l)]) as i64)
                    .collect::<Vec<_>>(),
            );
            if l.is_empty() {
                println!("\tFull cancellation to order 1");
            }
            for (t, a) in l {
                println!("\t{}: {}", t, a);
            }

            limits.push(expr);
        }
        limits
    }

    fn wood<E, V, H, S: SubGraph<Base = BitVec>>(&self, subgraph: &S) -> Wood
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        Wood::from_spinneys(self.spinneys(subgraph), self)
    }

    fn dod<S: SubGraph>(&self, subgraph: &S) -> i32;

    fn spinneys<E, V, H, S: SubGraph<Base = BitVec>>(
        &self,
        subgraph: &S,
    ) -> AHashSet<InternalSubGraph>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        let ref_graph = self.as_ref();
        let init_node = ref_graph.iter_nodes_of(subgraph).next().unwrap().0;
        let all_subcycles: Vec<_> = Cycle::all_sum_powerset_filter_map(
            &ref_graph
                .paton_cycle_basis(subgraph, &init_node, None)
                .unwrap()
                .0,
            &Some,
        )
        .map(|a| a.into_iter().map(|c| c.internal_graph(ref_graph)).collect())
        .unwrap();

        // println!("{}", self.base_dot());
        let mut spinneys: AHashSet<_> = InternalSubGraph::all_ops_iterative_filter_map(
            &all_subcycles,
            &|a, b| a.union(b),
            &|union| {
                if self.dod(&union) >= 0 {
                    Some(union)
                } else {
                    None
                }
            },
        );

        spinneys.insert(ref_graph.empty_subgraph());

        spinneys
    }
}

impl AsRef<HedgeGraph<Edge, Vertex, NumHedgeData>> for Graph {
    fn as_ref(&self) -> &HedgeGraph<Edge, Vertex, NumHedgeData> {
        &self.underlying
    }
}

impl UltravioletGraph for Graph {
    fn denominator<S: SubGraph>(&self, subgraph: &S) -> Atom {
        let mut den = Atom::num(1);

        for (pair, eid, d) in self.underlying.iter_edges_of(subgraph) {
            if matches!(pair, HedgePair::Paired { .. }) {
                let m2 = parse!(d.data.particle.mass.name).npow(2);
                den = den
                    * function!(
                        GS.den,
                        usize::from(eid) as i64,
                        function!(GS.emr_mom, usize::from(eid) as i64),
                        m2,
                        spenso_lor_atom(usize::from(eid) as i32, usize::from(eid), GS.dim)
                            .npow(2)
                            //.to_dots()
                            - m2
                    );
            }
        }

        den.into()
    }
    fn numerator<S: SubGraph>(&self, subgraph: &S) -> Atom {
        let num = Numerator::default();

        num.from_new_graph(self, subgraph, &GlobalPrefactor::default())
            .get_single_atom()
            .unwrap()
    }

    fn dod<S: SubGraph>(&self, subgraph: &S) -> i32 {
        let mut dod: i32 = 4 * self.n_loops(subgraph) as i32;
        // println!("nloops: {}", dod / 4);

        // FIXME: does not work if subgraph has external edges that contain both nodes!
        for (p, _, e) in self.underlying.iter_edges_of(subgraph) {
            if p.is_paired() {
                dod += e.data.dod;
            }
        }

        for (_, _, n) in self.underlying.iter_nodes_of(subgraph) {
            dod += n.dod;
        }

        dod
    }
}
