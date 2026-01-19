use ahash::AHashSet;
use linnet::half_edge::{
    HedgeGraph, PowersetIterator,
    involution::{Hedge, HedgePair},
    subgraph::{
        Cycle, InternalSubGraph, ModifySubSet, SuBitGraph, SubGraphLike, SubGraphOps, SubSetLike,
        SubSetOps,
    },
};
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    function,
};

use crate::{
    gammaloop_integrand::param_builder::ParamBuilderGraph,
    graph::{Edge, FeynmanGraph, Graph, LMBext, LoopMomentumBasis, NumHedgeData, Vertex},
    momentum_sample::LoopIndex,
    numerator::{AppliedFeynmanRule, Numerator},
    utils::{GS, W_},
};

use super::{Wood, spenso_lor_atom};

pub trait UltravioletGraph: LMBext + FeynmanGraph + ParamBuilderGraph {
    fn n_loops<S: SubGraphLike, E, V, H>(&self, subgraph: &S) -> usize
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        self.as_ref().cyclotomatic_number(subgraph)
    }

    // fn concretize_spinney<S: SubGraphLike<Base = SuBitGraph>>(&self, subgraph: &S) -> Graph;

    fn dummy_less_full_crown<S: SubGraphLike>(&self, subgraph: &S) -> S::Base
    where
        S::Base: ModifySubSet<Hedge> + SubGraphOps;

    ///Get the numerator of the graph.
    /// If multiply_prefactor is true, the numerator is multiplied by the global  prefactor. (just num not projector)
    fn numerator<S: SubGraphLike>(&self, subgraph: &S) -> Numerator<AppliedFeynmanRule>;
    fn denominator<S: SubGraphLike, T: Fn(&Edge) -> isize>(
        &self,
        subgraph: &S,
        edge_powers: T,
    ) -> Atom;
    fn all_cycle_unions<E, V, H, S: SubGraphLike<Base = SuBitGraph>>(
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
    fn all_limits<E, V, H, S: SubGraphLike>(
        &self,
        subgraph: &S,
        expr: &Atom,
        expansion: Symbol,
        lmb: &LoopMomentumBasis,
    ) -> Vec<(SuBitGraph, Atom)>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        let mom_reps = self.normal_emr_replacement(subgraph, lmb, &[W_.x___], |_s| true);

        let ose_reps = self.get_ose_replacements();
        for x in &mom_reps {
            println!("LMB replacement: {x}");
        }

        for r in &ose_reps {
            println!("ose{r}");
        }

        let expr = expr
            .replace(function!(GS.broadcasting_sqrt, W_.a_))
            .with(Atom::var(W_.a_).sqrt())
            .replace_multiple(&ose_reps)
            .replace_multiple(&mom_reps);
        // .replace_multiple(&q3_reps);
        let mut loops = PowersetIterator::new(lmb.loop_edges.len() as u8);

        let mut limits = Vec::new();

        loops.next();

        for ls in loops {
            let mut expr = expr.clone();
            for l in ls.included_iter() {
                let e = usize::from(lmb.loop_edges[LoopIndex(l.0)]) as i64;
                expr = expr
                    .replace(function!(GS.emr_mom, e, W_.x___))
                    .with(function!(GS.emr_mom, e, W_.x___) / expansion);

                expr /= Atom::var(expansion).npow(3);
            }

            let series = expr.series(expansion, Atom::Zero, 0.into(), true).unwrap();

            expr = series.to_atom().expand();

            let l = expr.coefficient_list::<i8>(&[Atom::var(expansion)]);

            println!(
                "LIMIT {:?}:",
                ls.included_iter()
                    .map(|l| usize::from(lmb.loop_edges[LoopIndex(l.0)]) as i64)
                    .collect::<Vec<_>>(),
            );
            if l.is_empty() {
                println!("\tFull cancellation to order 1");
            }
            for (t, a) in l {
                println!("\t{}: {}", t, a);
            }

            limits.push((ls, expr));
        }
        limits
    }

    fn wood<E, V, H, S: SubGraphLike<Base = SuBitGraph>>(&self, subgraph: &S) -> Wood
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        Wood::from_spinneys(self.spinneys(subgraph), self)
    }

    fn dod<S: SubGraphLike>(&self, subgraph: &S) -> i32;

    fn spinneys<E, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
    ) -> AHashSet<InternalSubGraph>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        let ref_graph = self.as_ref();

        if subgraph.is_empty() {
            let mut spinneys = AHashSet::new();
            spinneys.insert(ref_graph.empty_subgraph());
            return spinneys;
        }

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

        // for s in &all_subcycles {
        //     println!("Subcycle: {}", self.as_ref().dot(s));
        // }

        let mut spinneys: AHashSet<_> = InternalSubGraph::all_ops_iterative_filter_map(
            &all_subcycles,
            &|a, b| a.union(b),
            &|union| {
                // println!("{}", self.as_ref().dot(&union));
                if self.dod(&union) >= 0 {
                    Some(union)
                } else {
                    // println!("Negative dod:{}", self.dod(&union));
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
    fn dummy_less_full_crown<S: SubGraphLike>(&self, subgraph: &S) -> S::Base
    where
        S::Base: ModifySubSet<Hedge>,
    {
        let a = self.full_crown(subgraph);
        let mut ac = a.clone();

        a.included_iter().for_each(|a| {
            if self[self[&a]].is_dummy {
                ac.sub(a);
            }
        });

        ac
    }

    // fn concretize_spinney<S: SubGraphLike<Base = SuBitGraph>>(&self, s: &S) -> Graph {
    //     let mut underlying = self.underlying.clone();
    //     let dod = self.dod(s);
    //     let cc = underlying.connected_components(s);

    //     for c in cc{

    //     let v = Vertex {

    //     }
    //     underlying.identify_nodes_without_self_edges(s, node_data_merge);
    //     };
    //     }

    //     self.clone()
    // }

    fn denominator<S: SubGraphLike, T: Fn(&Edge) -> isize>(
        &self,
        subgraph: &S,
        edge_powers: T,
    ) -> Atom {
        let mut den = Atom::num(1);

        for (pair, eid, d) in self.underlying.iter_edges_of(subgraph) {
            if matches!(pair, HedgePair::Paired { .. }) {
                let m2 = d.data.mass_atom().npow(2);
                let edge_power = edge_powers(d.data);
                let is_power_negative = edge_power < 0;
                let prop_den = function!(
                    GS.den,
                    usize::from(eid) as i64,
                    function!(GS.emr_mom, usize::from(eid) as i64),
                    m2,
                    spenso_lor_atom(usize::from(eid) as i32, usize::from(eid), GS.dim)
                            .npow(2)
                            //.to_dots()
                            - m2
                );
                for _i in 0..edge_power.abs() {
                    if is_power_negative {
                        den /= prop_den.clone();
                    } else {
                        den *= prop_den.clone();
                    }
                }
            }
        }

        den
    }
    fn numerator<S: SubGraphLike>(&self, subgraph: &S) -> Numerator<AppliedFeynmanRule> {
        let num = Numerator::default();

        num.from_new_graph(self, subgraph)
    }

    fn dod<S: SubGraphLike>(&self, subgraph: &S) -> i32 {
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

        // println!("Dod:{dod} for subgraph:{}", self.dot(subgraph));

        dod
    }
}

pub trait UVE {
    fn mass_atom(&self) -> Atom;
}
