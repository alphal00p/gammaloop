use std::{collections::BTreeSet, ops::Deref};

use ahash::AHashSet;
use idenso::metric::MetricSimplifier;
use linnet::half_edge::{
    HedgeGraph, PowersetIterator,
    involution::{Flow, Hedge, HedgePair},
    subgraph::{
        Cycle, InternalSubGraph, ModifySubSet, PairwiseSubSetOps, SuBitGraph, SubGraphLike,
        SubGraphOps, SubSetLike, SubSetOps, subset::SubSet,
    },
};
use spenso::network::library::TensorLibraryData;
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    domains::atom::AtomField,
    function,
    poly::series::Series,
};
use tracing::debug;

use crate::{
    graph::{Edge, FeynmanGraph, Graph, HedgeData, LMBext, LoopMomentumBasis, Vertex},
    integrands::process::param_builder::ParamBuilderGraph,
    momentum::sample::LoopIndex,
    numerator::{AppliedFeynmanRule, Numerator},
    utils::{
        GS, W_,
        symbolica_ext::{CallSymbol, DOD},
    },
    uv::{ApproximationType, UVgenerationSettings, settings::CTIdentifier},
};

use super::{Spinney, Wood, spenso_lor_atom};

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
    fn numerator<S: SubGraphLike + SubSetOps>(
        &self,
        subgraph: &S,
        without: &S,
    ) -> Numerator<AppliedFeynmanRule>;
    fn denominator<S: SubGraphLike, T: Fn(&Edge) -> isize>(
        &self,
        subgraph: &S,
        edge_powers: T,
    ) -> Atom;

    fn boundary_pdg_set<E: UVE, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
    ) -> BTreeSet<isize>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        let graph = self.as_ref();
        graph
            .full_crown(subgraph)
            .included_iter()
            .filter_map(|hedge| {
                let edge_id = graph[&hedge];
                graph[edge_id].particle_pdg_code().map(|pdg| {
                    if graph.flow(hedge) == Flow::Source {
                        -pdg
                    } else {
                        pdg
                    }
                })
            })
            .collect()
    }

    fn internal_pdg_set<E: UVE, V, H, S: SubGraphLike>(&self, subgraph: &S) -> BTreeSet<isize>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        self.as_ref()
            .iter_edges_of(subgraph)
            .filter_map(|(pair, edge_id, _)| {
                pair.is_paired()
                    .then(|| self.as_ref()[edge_id].particle_pdg_code())
                    .flatten()
            })
            .collect()
    }

    fn ct_identifier<E: UVE, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
    ) -> CTIdentifier
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        CTIdentifier::new(
            self.boundary_pdg_set(subgraph),
            Some(self.internal_pdg_set(subgraph)),
        )
    }

    fn approximation_scheme<E: UVE, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
        settings: &UVgenerationSettings,
    ) -> ApproximationType
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        settings.approximation_scheme_for(&self.ct_identifier(subgraph))
    }

    fn classify_spinney<E: UVE, V, H>(
        &self,
        spinney: InternalSubGraph,
        settings: &UVgenerationSettings,
        lmb: &LoopMomentumBasis,
    ) -> Option<Spinney>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        let renormalization_scheme = self.approximation_scheme(&spinney.filter, settings);

        (renormalization_scheme != ApproximationType::Unsubtracted)
            .then(|| Spinney::with_scheme(spinney, self, lmb, renormalization_scheme))
    }

    fn classified_spinneys<E: UVE, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
        settings: &UVgenerationSettings,
        lmb: &LoopMomentumBasis,
    ) -> Vec<Spinney>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        if !settings.subtract_uv {
            return vec![Spinney::empty(self)];
        }

        self.spinneys(subgraph)
            .into_iter()
            .filter_map(|spinney| self.classify_spinney(spinney, settings, lmb))
            .collect()
    }

    fn all_cycle_unions<E, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
    ) -> AHashSet<InternalSubGraph>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        let ref_graph = self.as_ref();
        let _init_node = ref_graph.iter_nodes_of(subgraph).next().unwrap().0;
        let all_subcycles: Vec<_> =
            Cycle::all_sum_powerset_filter_map(&ref_graph.cycle_basis_of(subgraph).0, &Some)
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
    ) -> Vec<(SubSet<LoopIndex>, Series<AtomField>)>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        let mom_reps = self.normal_emr_replacement(subgraph, lmb, &[W_.x___], |_s| true);

        let ose_reps = self.get_ose_replacements();
        // for x in &mom_reps {
        //     println!("LMB replacement: {x}");
        // }

        // for r in &ose_reps {
        //     println!("ose{r}");
        // }

        let expr = expr
            .replace(function!(GS.broadcasting_sqrt, W_.a_))
            .with(Atom::var(W_.a_).sqrt())
            .replace_multiple(&ose_reps)
            .replace_multiple(&mom_reps)
            .replace(GS.if_sigma.f(&[W_.a_]))
            .with(Atom::one());
        // .replace_multiple(&q3_reps);
        let mut loops = PowersetIterator::<LoopIndex>::new(lmb.loop_edges.len() as u8);

        let mut limits = Vec::new();

        loops.next();

        for ls in loops {
            let mut expr = expr.clone();
            for l in ls.included_iter() {
                let e = usize::from(lmb.loop_edges[l]) as i64;
                expr = expr
                    .replace(function!(GS.emr_mom, e, W_.x___))
                    .with(function!(GS.emr_mom, e, W_.x___) / expansion);

                expr /= Atom::var(expansion).pow(3);
            }

            let series = expr.series(expansion, Atom::Zero, 0.into(), true).unwrap();

            // expr = series.to_atom().expand();

            limits.push((ls, series));
        }
        limits
    }

    fn wood<E: UVE, V, H, S: SubGraphLike<Base = SuBitGraph>>(&self, subgraph: &S) -> Wood
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        self.wood_with_settings(
            subgraph,
            &UVgenerationSettings::default(),
            &self.as_ref().lmb_of(&self.as_ref().full_filter()),
        )
    }

    fn wood_with_settings<E: UVE, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
        settings: &UVgenerationSettings,
        lmb: &LoopMomentumBasis,
    ) -> Wood
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        Wood::from_spinneys(self.classified_spinneys(subgraph, settings, lmb), self)
    }

    fn dod<S: SubGraphLike<Base = SuBitGraph> + SubSetOps>(&self, subgraph: &S) -> i32;
    fn local_dod<S: SubGraphLike>(&self, subgraph: &S) -> i32;

    fn spinneys<E, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
    ) -> AHashSet<InternalSubGraph>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        let ref_graph: &HedgeGraph<E, V, H> = self.as_ref();
        debug!(subgraph=%ref_graph.dot(subgraph),"Spinneys of subgraph");
        let _b: SuBitGraph = ref_graph.empty_subgraph();

        if subgraph.is_empty() {
            let mut spinneys = AHashSet::new();
            spinneys.insert(ref_graph.empty_subgraph());
            return spinneys;
        }

        let cycles = ref_graph.cycle_basis_of(subgraph).0;
        for c in &cycles {
            debug!("Cycle: {}", ref_graph.dot(&c.filter));
        }

        let _init_node = ref_graph.iter_nodes_of(subgraph).next().unwrap().0;
        let all_subcycles: Vec<_> = Cycle::all_sum_powerset_filter_map(&cycles, &Some)
            .map(|a| a.into_iter().map(|c| c.internal_graph(ref_graph)).collect())
            .unwrap();

        for s in &all_subcycles {
            debug!("Subcycle: {}", self.as_ref().dot(s));
        }

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

impl AsRef<HedgeGraph<Edge, Vertex, HedgeData>> for Graph {
    fn as_ref(&self) -> &HedgeGraph<Edge, Vertex, HedgeData> {
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
                let m2 = d.data.mass_atom().pow(2);
                let edge_power = edge_powers(d.data);
                let is_power_negative = edge_power < 0;
                let prop_den = GS.den(
                    usize::from(eid),
                    function!(GS.emr_mom, usize::from(eid)),
                    &m2,
                    spenso_lor_atom(usize::from(eid) as i32, usize::from(eid), GS.dim)
                        .pow(2)
                        .to_dots()
                        - &m2,
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
    fn numerator<S: SubGraphLike + SubSetOps>(
        &self,
        subgraph: &S,
        without: &S,
    ) -> Numerator<AppliedFeynmanRule> {
        let num = Numerator::default();

        num.fill_in_reduced(self, subgraph, without)
    }

    fn dod<S: SubGraphLike<Base = SuBitGraph> + SubSetOps>(&self, subgraph: &S) -> i32 {
        let lmb = self.lmb_of(subgraph);
        let empty = self.underlying.empty_subgraph();
        let integrand = self
            .numerator(subgraph, &empty)
            .to_d_dim(GS.dim)
            .get_single_atom()
            .unwrap()
            / self.denominator(subgraph, |_| 1);
        let nloops: usize = self.n_loops(subgraph);
        self.uv_rescaled(subgraph.included(), nloops, &lmb, &integrand)
            .trailing_exponent()
    }

    fn local_dod<S: SubGraphLike>(&self, subgraph: &S) -> i32 {
        let mut dod: i32 = 4 * self.n_loops(subgraph) as i32;
        for (p, _, e) in self.underlying.iter_edges_of(subgraph) {
            if p.is_paired() {
                dod += e.data.dod.deref();
            }
        }

        for (_, _, n) in self.underlying.iter_nodes_of(subgraph) {
            dod += n.dod.deref();
        }

        dod
    }
}

pub trait UVE {
    fn mass_atom(&self) -> Atom;
    fn particle_pdg_code(&self) -> Option<isize>;
}
