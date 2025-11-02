#![allow(dead_code)]

use std::hash::Hash;

use crate::{
    cff::{
        expression::{GraphOrientation, OrientationID},
        generation::{ShiftRewrite, generate_uv_cff},
    },
    graph::{Edge, Graph, LMBext, LoopMomentumBasis, Vertex},
    momentum::Sign,
    numerator::{ParsingNet, symbolica_ext::AtomCoreExt},
    utils::{
        GS, W_,
        symbolica_ext::{CallSymbol, LOGPRINTOPTS},
    },
};
use ahash::AHashSet;
use bitvec::vec::BitVec;
use idenso::{gamma::GammaSimplifier, metric::MS};
use log::debug;

use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder, Symbol},
    function, parse, symbol,
};

use linnet::half_edge::{
    HedgeGraph,
    involution::{EdgeIndex, HedgePair},
    subgraph::{
        Inclusion, InternalSubGraph, SuBitGraph, SubGraphLike, SubGraphOps, SubSetLike, SubSetOps,
    },
};

use typed_index_collections::TiVec;
use vakint::{Vakint, VakintExpression, vakint_symbol};
// use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};

use super::{IntegrandExpr, UltravioletGraph, uv_graph::UVE};

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum ApproxOp {
    NotComputed,
    Union {
        sign: Sign,
        t_args: Vec<IntegrandExpr>,
        subgraphs: Vec<InternalSubGraph>,
    },
    Dependent {
        sign: Sign,
        t_arg: IntegrandExpr,
        subgraph: InternalSubGraph,
    },
    Root,
}

#[derive(Clone)]
pub struct SimpleApprox {
    t_args: Vec<Atom>,
    pub sign: Sign,
    graph: InternalSubGraph,
}

impl SimpleApprox {
    fn subgraph_shadow(graph: &SuBitGraph, subgraph: &InternalSubGraph) -> Symbol {
        symbol!(&format!(
            "S_{}⊛{}",
            graph.string_label(),
            subgraph.string_label()
        ))
    }

    pub(crate) fn expr(&self, bigger_graph: &SuBitGraph) -> Atom {
        let reduced = Atom::var(Self::subgraph_shadow(bigger_graph, &self.graph));
        let mut mul = Atom::num(1);
        for i in &self.t_args {
            mul = mul * i
        }
        reduced * mul
    }

    pub(crate) fn t_op(&self, bigger_graph: &SuBitGraph) -> Atom {
        function!(GS.top, self.expr(bigger_graph))
    }

    pub(crate) fn root(subgraph: InternalSubGraph) -> Self {
        if !subgraph.is_empty() {
            panic!(
                "Root approximation must be empty {} {:?}",
                subgraph.string_label(),
                subgraph
            )
        }
        SimpleApprox {
            t_args: vec![],
            sign: Sign::Positive,
            graph: subgraph,
        }
    }

    pub(crate) fn dependent(&self, bigger_graph: InternalSubGraph) -> Self {
        Self {
            t_args: vec![self.t_op(&bigger_graph.filter)],
            sign: -self.sign,
            graph: bigger_graph,
        }
    }

    pub(crate) fn union<'a>(
        subgraph: InternalSubGraph,
        union: impl IntoIterator<Item = &'a Self>,
    ) -> Self {
        let mut t_args = vec![];
        let mut sign = Sign::Positive;
        for u in union {
            if u.t_args.len() != 1 {
                panic!("Union can only be applied to dependent approximations");
            }
            t_args.push(u.t_args[0].clone());
            sign = sign * u.sign;
        }
        SimpleApprox {
            t_args,
            sign,
            graph: subgraph,
        }
    }
}

pub(crate) fn to_vakint_integrand<E: UVE, V, H, S: SubSetLike, SS: SubSetLike>(
    integrand: &Atom,
    lmb: &LoopMomentumBasis,
    graph: &HedgeGraph<E, V, H>,
    reduced: &S,
    dependent_subgraph: &SS,
    substitute_masses_to_m_uv: bool,
) -> Atom {
    let mut integrand_vakint = integrand.clone();

    //Atom::Zero

    // strip the momentum wrapper from the denominator
    integrand_vakint = integrand_vakint
        .replace(function!(
            GS.den,
            W_.prop_,
            function!(GS.emr_mom, W_.prop_, W_.mom_),
            W_.x__
        ))
        .with(function!(GS.den, W_.prop_, W_.mom_, W_.x__));

    // println!("Expanded: {:>}", integrand_vakint.expand());

    integrand_vakint = integrand_vakint.expand();

    let mut propagator_id = 1;

    let vk_prop = vakint_symbol!("prop");
    let vk_edge = vakint_symbol!("edge");
    let vk_mom = vakint_symbol!("k");
    let vk_topo = vakint_symbol!("topo");

    for (pair, index, _data) in graph.iter_edges_of(reduced) {
        if let HedgePair::Paired { source, sink } = pair {
            integrand_vakint = integrand_vakint
                .replace(function!(
                    GS.den,
                    usize::from(index) as i64,
                    W_.mom_,
                    W_.mass_,
                    W_.x___
                ))
                .with(function!(
                    vk_prop,
                    propagator_id,
                    function!(
                        vk_edge,
                        usize::from(graph.node_id(source)) as i64,
                        usize::from(graph.node_id(sink)) as i64
                    ),
                    W_.mom_,
                    if substitute_masses_to_m_uv {
                        GS.m_uv
                    } else {
                        W_.mass_
                    },
                    1
                ))
                .replace(function!(vk_prop, W_.x___, 1).pow(Atom::var(W_.e_)))
                .with(function!(vk_prop, W_.x___, -Atom::var(W_.e_)));
            propagator_id += 1;
        }
    }

    // shrink vertices of the subgraph
    for (pair, _index, _data) in graph.iter_edges_of(dependent_subgraph) {
        if let HedgePair::Paired { source, sink } = pair {
            integrand_vakint = integrand_vakint
                .replace(function!(
                    vk_prop,
                    W_.x_,
                    function!(vk_edge, usize::from(graph.node_id(source)) as i64, W_.y_),
                    W_.x___
                ))
                .with(function!(
                    vk_prop,
                    W_.x_,
                    function!(vk_edge, usize::from(graph.node_id(sink)) as i64, W_.y_),
                    W_.x___
                ))
                .replace(function!(
                    vk_prop,
                    W_.x_,
                    function!(vk_edge, W_.y_, usize::from(graph.node_id(source)) as i64),
                    W_.x___
                ))
                .with(function!(
                    vk_prop,
                    W_.x_,
                    function!(vk_edge, W_.y_, usize::from(graph.node_id(sink)) as i64),
                    W_.x___
                ));
        }
    }

    // flip edges to positive momentum
    // FIXME: how will this work for sums of momenta?
    integrand_vakint = integrand_vakint
        .replace(function!(
            vk_prop,
            W_.x_,
            function!(vk_edge, W_.a_, W_.b_),
            -Atom::var(W_.y_),
            W_.e___
        ))
        .repeat()
        .with(function!(
            vk_prop,
            W_.x_,
            function!(vk_edge, W_.b_, W_.a_),
            W_.y_,
            W_.e___
        ));

    // fuse raised edges
    integrand_vakint = integrand_vakint
        .replace(
            function!(
                vk_prop,
                W_.x_,
                function!(vk_edge, W_.a_, W_.b_),
                W_.x___,
                W_.e_
            ) * function!(
                vk_prop,
                W_.y_,
                function!(vk_edge, W_.b_, W_.c_),
                W_.x___,
                W_.f_
            ),
        )
        .repeat()
        .with(function!(
            vk_prop,
            W_.y_,
            function!(vk_edge, W_.a_, W_.c_),
            W_.x___,
            W_.e_ + W_.f_
        ));

    debug!("Integrand pre dot vakint: {:}", integrand_vakint);
    // rewrite numerator
    // linearize the numerator first
    // integrand_vakint = integrand_vakint
    //     .replace(function!(GS.emr_mom, W_.prop_, W_.a___))
    //     .with(GS.linearize.f(&[function!(GS.emr_mom, W_.a___)]));

    // // rewrite numerator
    // // linearize the numerator first
    integrand_vakint = integrand_vakint
        .replace(function!(GS.emr_mom, W_.prop_, W_.mom_, W_.x_))
        .with(function!(MS.dot, W_.mom_, W_.x_))
        .replace(function!(MS.dot, W_.mom_, W_.x_))
        .with(function!(GS.emr_mom, W_.mom_, W_.x_));

    debug!("Integrand pre vakint: {:}", integrand_vakint);
    // panic!("FUFU");
    for (i, l) in lmb.loop_edges.iter().enumerate() {
        integrand_vakint = integrand_vakint
            .replace(function!(GS.emr_mom, usize::from(*l) as i64))
            .with(function!(vk_mom, i as i64 + 1))
            .replace(function!(
                GS.emr_mom,
                function!(vk_mom, i as i64 + 1),
                W_.x___
            ))
            .with(function!(vk_mom, i as i64 + 1, W_.x___));
    }

    // collect the topology
    integrand_vakint = integrand_vakint
        .replace(function!(vk_prop, W_.x__))
        .with(function!(vk_topo, function!(vk_prop, W_.x__)))
        .replace(function!(vk_topo, W_.x_) * function!(vk_topo, W_.y_))
        .repeat()
        .with(function!(vk_topo, W_.x_ * W_.y_));

    debug!("Integrand vakint: {:#}", integrand_vakint);

    integrand_vakint
}

#[derive(Clone)]
pub struct Approximation {
    // The union of all spinneys, remaining graph is full graph minus subgraph
    pub subgraph: InternalSubGraph,
    pub dod: i32,
    pub lmb: LoopMomentumBasis,
    pub local_3d: CFFapprox, //3d denoms
    pub final_integrand: Option<ParsingNet>,
    pub integrated_4d: ApproxOp, //4d
    pub simple_approx: Option<SimpleApprox>,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum CFFapprox {
    NotComputed,
    Dependent { sign: Sign, t_arg: IntegrandExpr },
}

impl CFFapprox {
    pub(crate) fn expr(&self) -> Option<(Atom, Sign)> {
        match self {
            CFFapprox::NotComputed => None,
            CFFapprox::Dependent { sign, t_arg } => Some((t_arg.integrand.clone(), *sign)),
        }
    }

    pub(crate) fn dependent<
        E,
        V,
        H,
        S: SubGraphLike,
        SS: SubGraphLike,
        OID: OrientationID,
        O: GraphOrientation,
    >(
        graph: &HedgeGraph<E, V, H>,
        to_contract: &SS,
        amplitude_subgraph: &S,
        canonize_esurface: &Option<ShiftRewrite>,
        orientations: &TiVec<OID, O>,
        cut_edges: &[EdgeIndex],
    ) -> CFFapprox {
        let mut cff_sum = Atom::Zero;

        let mut contract_edges = Vec::new();

        for (e, eid, _) in graph.iter_edges_of(to_contract) {
            if e.is_paired() {
                contract_edges.push(eid);
            }
        }

        for o in orientations {
            cff_sum += o.orientation_thetas()
                * generate_uv_cff(
                    graph,
                    amplitude_subgraph,
                    canonize_esurface,
                    &contract_edges,
                    o.orientation(),
                    cut_edges,
                )
                .unwrap()
        }
        CFFapprox::Dependent {
            sign: Sign::Positive,
            t_arg: IntegrandExpr { integrand: cff_sum },
        }
    }

    pub(crate) fn root<E, V, H, S: SubGraphLike, OID: OrientationID, O: GraphOrientation>(
        graph: &HedgeGraph<E, V, H>,
        amplitude_subgraph: &S,
        canonize_esurface: &Option<ShiftRewrite>,
        orientations: &TiVec<OID, O>,
        cut_edges: &[EdgeIndex],
    ) -> CFFapprox {
        Self::dependent(
            graph,
            &graph.empty_subgraph::<SuBitGraph>(),
            amplitude_subgraph,
            canonize_esurface,
            orientations,
            cut_edges,
        )
    }
}
impl Approximation {
    pub(crate) fn dot(&self, graph: &Graph) -> String {
        graph.dot_lmb(&self.subgraph, &self.lmb)
    }

    pub(crate) fn root<
        H,
        G: UltravioletGraph + AsRef<HedgeGraph<Edge, Vertex, H>>,
        S: SubGraphLike<Base = SuBitGraph>,
        OID: OrientationID,
        O: GraphOrientation,
    >(
        &mut self,
        graph: &G,
        amplitude: &S,
        canonize_esurface: &Option<ShiftRewrite>,
        orientations: &TiVec<OID, O>,
        cut_edges: &[EdgeIndex],
    ) {
        self.local_3d = CFFapprox::root(
            graph.as_ref(),
            amplitude,
            canonize_esurface,
            orientations,
            cut_edges,
        );
        self.integrated_4d = ApproxOp::Root;
        self.final_integrand =
            self.final_integrand(graph, amplitude, canonize_esurface, orientations, cut_edges)
    }

    pub(crate) fn new<G, E, V, H>(
        spinney: InternalSubGraph,
        graph: &G,
        lmb: &LoopMomentumBasis,
    ) -> Approximation
    where
        G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>,
    {
        let lmb = graph.compatible_sub_lmb(&spinney, graph.dummy_less_full_crown(&spinney), lmb);
        // println!("//lmb for spinney \n{}", graph.dot_lmb(&spinney, &lmb));
        Approximation {
            dod: graph.dod(&spinney),
            subgraph: spinney,
            lmb,
            final_integrand: None,
            simple_approx: None,
            local_3d: CFFapprox::NotComputed,
            integrated_4d: ApproxOp::NotComputed,
        }
    }

    pub(crate) fn integrated_4d<
        E: UVE,
        V,
        H,
        G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>,
        S: SubGraphLike,
    >(
        &self,
        dependent: &Self,
        vakint: &Vakint,
        uv_graph: &G,
        amplitude_subgraph: &S,
    ) -> ApproxOp {
        let graph = uv_graph.as_ref();
        let reduced = self.subgraph.subtract(&dependent.subgraph);

        let Some((inner_t, sign)) = dependent.integrated_4d.expr() else {
            return ApproxOp::NotComputed;
        };
        let t_arg = uv_graph
            .numerator(&reduced)
            .to_d_dim(GS.dim)
            .get_single_atom()
            .unwrap()
            .simplify_gamma()
            / uv_graph.denominator(&reduced, |_| 1);

        let ep = vakint_symbol!("ε");

        // strip the pole part of inner_t as this is removed by MS bar
        let mut pole_stripped = inner_t
            .series(
                ep,
                Atom::Zero,
                (uv_graph.n_loops(amplitude_subgraph) as i64 + 1).into(),
                true,
            )
            .unwrap()
            .to_atom();
        for i in -(uv_graph.n_loops(&dependent.subgraph) as i64)..0 {
            pole_stripped = pole_stripped
                .replace(Atom::var(ep).npow(i))
                .with(Atom::Zero);
        }

        let mut atomarg = t_arg * pole_stripped;

        // only apply replacements for edges in the reduced graph
        let mom_reps = graph.uv_wrapped_replacement(&reduced, &self.lmb, &[W_.x___]);

        println!("Reps:");
        for r in &mom_reps {
            println!("{r}");
        }

        println!(
            "Expand-prerep {} with dod={} in {:?}",
            atomarg, self.dod, self.lmb.ext_edges
        );

        // rewrite the inner_t as well
        atomarg = atomarg.replace_multiple(&mom_reps);

        println!(
            "Expand {} with dod={} in {:?}",
            atomarg, self.dod, self.lmb.ext_edges
        );
        for e in &self.lmb.ext_edges {
            atomarg = atomarg
                .replace(function!(GS.emr_mom, usize::from(*e) as i64, W_.x___))
                .with(function!(GS.emr_mom, usize::from(*e) as i64, W_.x___) * GS.rescale);
        }
        // TODO: only enable soft CT if doing OS renormalization
        let soft_ct = graph.full_crown(&self.subgraph).n_included() == 2 && self.dod > 0;

        let mut masses = AHashSet::new();
        masses.insert(Atom::var(GS.m_uv));
        // scale all masses, including UV masses from subgraphs

        for (p, _, e) in graph.iter_edges_of(&self.subgraph) {
            if p.is_paired() {
                let e_mass = e.data.mass_atom();
                masses.insert(e_mass);
            }
        }

        if !soft_ct {
            for m in &masses {
                let rescaled = m.clone() * GS.rescale;
                atomarg = atomarg.replace(m.clone()).with(rescaled);
            }

            // expand the propagator around a propagator with a UV mass
            atomarg = atomarg
                .replace(parse!("den(n_,q_,mass_,prop_)"))
                .with(parse!(
                    "den(n_,q_,mass_ + mUV^2 - t^2*mUV^2, prop_- mUV^2 + t^2*mUV^2)"
                ));
        }

        atomarg = atomarg
            .replace(function!(MS.dot, GS.rescale * W_.x_, W_.y_))
            .repeat()
            .with(function!(MS.dot, W_.x_, W_.y_) * GS.rescale);

        // println!("atomarg:{}", atomarg);

        // den(..) tags a propagator, its first derivative is 1 and the rest is 0
        let mut a = atomarg
            .series(GS.rescale, Atom::Zero, self.dod.into(), true)
            .unwrap()
            .to_atom()
            .replace(parse!("der(0,0,0,1, den(y__))"))
            .with(Atom::num(1))
            .replace(parse!("der(x__, den(y__))"))
            .with(Atom::num(0));

        if soft_ct {
            let coeffs = a.coefficient_list::<u8>(&[Atom::var(GS.rescale)]);
            let mut b = Atom::Zero;
            let dod_pow = Atom::var(GS.rescale).npow(self.dod);
            for (pow, mut i) in coeffs {
                if pow == dod_pow {
                    // set the masses in the t=dod term to 0
                    // UV rearrange the denominators
                    for m in &masses {
                        i = i.replace(m.clone()).with(Atom::Zero);
                    }

                    i = i
                        .replace(parse!("den(n_,q_,mass_,prop_)"))
                        .with(parse!("den(n_,q_,mUV^2,prop_-mUV^2)"));
                }

                b += i;
            }

            a = b;
        } else {
            a = a.replace(GS.rescale).with(Atom::num(1));
        }

        let integrand_vakint =
            to_vakint_integrand(&a, &self.lmb, &graph, &reduced, &dependent.subgraph, true);

        let vakint_expr = VakintExpression::try_from(integrand_vakint.clone()).unwrap();
        println!("\nVakint expression:\n{:#}", vakint_expr);

        let mut res = vakint.evaluate(integrand_vakint.as_view()).unwrap();

        let vk_metric = vakint_symbol!("g");
        // apply metric
        res = res
            .replace(function!(vk_metric, W_.x_, W_.y_) * function!(GS.emr_mom, W_.x___, W_.x_))
            .with(function!(GS.emr_mom, W_.x___, W_.y_));

        res = res.replace(vakint::symbols::S.cmplx_i).with(Atom::i());

        // multiply the results with a vacuum triangle that integrates to 1
        // 1/(k^2 - m_UV^2)^3 = -i / (4 pi)^2 * 1/2 * 1/mUV^2
        // name the mUV mass mUVi as this one should not be expanded
        for l in &self.lmb.loop_edges {
            if !reduced.includes(&graph[l].1) {
                continue;
            }

            //TODO: Add orientation localisation prefactor (Sum of valid orientation thetas)/(number of valid orientations)

            res /= parse!("(-1i / (4 𝜋)^2 * 1/2 * 1/mUVI^2)");

            // multiply CFF triangle
            res *= Atom::num((3, 16))
                / ((-function!(
                    MS.dot,
                    function!(GS.emr_vec, usize::from(*l) as i64),
                    function!(GS.emr_vec, usize::from(*l) as i64)
                ) + GS.m_uv_int * GS.m_uv_int)
                    .npow((5, 2)));
        }

        println!("\nIntegrated CT:\n{}\n", res);

        ApproxOp::Dependent {
            t_arg: IntegrandExpr { integrand: res },
            sign: -sign,
            subgraph: reduced,
        }
    }

    pub(crate) fn compute_integrated<
        E: UVE,
        V,
        H,
        G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>,
        S: SubGraphLike,
    >(
        &mut self,
        graph: &G,
        vakint: &Vakint,
        amplitude_subgraph: &S,
        dependent: &Self,
    ) {
        self.integrated_4d = self.integrated_4d(dependent, vakint, graph, amplitude_subgraph);
    }

    /// Computes the 3d approximation of the UV
    pub(crate) fn compute<
        H,
        G: UltravioletGraph + AsRef<HedgeGraph<Edge, Vertex, H>>,
        S: SubGraphLike<Base = SuBitGraph>,
        OID: OrientationID,
        O: GraphOrientation,
    >(
        &mut self,
        graph: &G,
        amplitude: &S,
        canonize_esurface: &Option<ShiftRewrite>,
        orientations: &TiVec<OID, O>,
        cut_edges: &[EdgeIndex],
        dependent: &Self,
    ) {
        let Some((cff, sign)) = dependent.local_3d.expr() else {
            panic!("Should have computed the dependent cff");
        };

        let (t4, _) = if let ApproxOp::Root = dependent.integrated_4d {
            (Atom::num(0), Sign::Positive)
        } else {
            dependent
                .integrated_4d
                .expr()
                .unwrap_or((Atom::num(0), Sign::Positive))
        };

        let CFFapprox::Dependent { t_arg, .. } = CFFapprox::dependent(
            graph.as_ref(),
            &dependent.subgraph,
            amplitude,
            canonize_esurface,
            orientations,
            cut_edges,
        ) else {
            unreachable!()
        };

        let mut sum_3d = Atom::Zero;

        sum_3d += self.local_3d(dependent, graph, cff);

        let finite = t4
            .series(vakint_symbol!("ε"), Atom::Zero, 0.into(), true)
            .unwrap()
            .coefficient((0, 1).into())
            .replace(GS.m_uv_int)
            .with(GS.m_uv)
            .map_mink_dim(4);

        debug!(
            "Integrated 4d finite part: {:#}",
            finite.printer(LOGPRINTOPTS)
        );

        // TODO: multiply by the number of orientations or only apply the counterterm to
        // one orientation
        // subtract integrated CT
        sum_3d -= self.local_3d(dependent, graph, finite * t_arg.integrand);

        self.local_3d = CFFapprox::Dependent {
            sign: -sign,
            t_arg: IntegrandExpr { integrand: sum_3d },
        };

        self.final_integrand =
            self.final_integrand(graph, amplitude, canonize_esurface, orientations, cut_edges);
    }

    pub(crate) fn local_3d<E: UVE, V, H, G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>>(
        &self,
        dependent: &Self,
        uv_graph: &G,
        mut cff: Atom,
    ) -> Atom {
        let graph = uv_graph.as_ref();
        let reduced = self.subgraph.subtract(&dependent.subgraph);

        // println!("CFF: {}", cff);

        // add data for OSE computation and add an explicit sqrt
        for (p, ei, e) in graph.iter_edges_of(&self.subgraph) {
            let eid = usize::from(ei) as i64;
            if p.is_paired() {
                // set energies from inner_t on-shell
                cff = cff.replace(function!(GS.energy, eid)).with(GS.ose(ei));

                let e_mass = e.data.mass_atom();
                cff = cff.replace(GS.ose(ei)).with(GS.ose_full(ei, e_mass, None));
            }
        }

        let mut atomarg = cff * uv_graph.numerator(&reduced).get_single_atom().unwrap();

        // println!(
        //     "Expand-prerep {} with dod={} in {:?}",
        //     atomarg, self.dod, self.lmb.ext_edges
        // );

        // split numerator momenta into OSEs and spatial parts
        let mut reps = Vec::new();
        for (p, eid, e) in graph.iter_edges_of(&self.subgraph) {
            if p.is_paired() {
                let e_mass = e.data.mass_atom();
                reps.push(GS.split_mom_pattern(eid, e_mass));
            }
        }

        // println!("Split reps:");
        // for r in &reps {
        //     println!("{r}");
        // }

        atomarg = atomarg.replace_multiple(&reps);

        // only apply replacements for edges in the reduced graph
        let mom_reps = graph.uv_spatial_wrapped_replacement(&reduced, &self.lmb, &[W_.x___]);

        // println!("Reps:");
        // for r in &mom_reps {
        //     println!("{r}");
        // }

        atomarg = atomarg.replace_multiple(&mom_reps);

        // rescale the loop momenta in the whole subgraph, including previously expanded cycles
        for e in &self.lmb.loop_edges {
            // println!("Rescale {}", e);
            atomarg = atomarg
                .replace(GS.emr_vec_index(*e, W_.x___))
                .with(GS.emr_vec_index(*e, W_.x___) * GS.rescale);
        }

        // println!(
        //     "Expand {} with dod={} in {:?}",
        //     atomarg, self.dod, self.lmb.ext_edges
        // );

        let soft_ct = graph.full_crown(&self.subgraph).n_included() == 2 && self.dod > 0;

        // (re-)expand OSEs from the subgraph only
        for (_, eid, _) in graph.iter_edges_of(&self.subgraph) {
            let eid = usize::from(eid) as i64;
            if soft_ct {
                // TODO: rescale the masses in OSEs
                // TODO: also scale masses in the numerator _only_ for the subgraph
                // expand the OSEs around an OSE with a UV mass
                todo!()
            } else {
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
        }

        atomarg = atomarg
            .replace(function!(MS.dot, GS.rescale * W_.x_, W_.y_))
            .repeat()
            .with(function!(MS.dot, W_.x_, W_.y_) * GS.rescale);

        atomarg = (atomarg
            * Atom::var(GS.rescale).npow(3 * uv_graph.n_loops(&self.subgraph) as i64))
        .replace(GS.rescale)
        .with(Atom::num(1) / GS.rescale);

        // println!("atomarg:{:>}", atomarg.expand());

        // println!("Series expanding");

        let mut a = atomarg
            .series(GS.rescale, Atom::Zero, self.dod.into(), true)
            .unwrap()
            .to_atom()
            .replace(parse!("der(0,0,0,1, OSE(y__))"))
            .with(Atom::num(1))
            .replace(parse!("der(x__, OSE(y__))"))
            .with(Atom::num(0));

        if soft_ct {
            let coeffs = a.coefficient_list::<u8>(&[Atom::var(GS.rescale)]);
            let mut b = Atom::Zero;
            let dod_pow = Atom::var(GS.rescale).npow(self.dod);
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

        a
    }

    pub(crate) fn final_integrand<
        G,
        H,
        S: SubGraphLike<Base = SuBitGraph>,
        OID: OrientationID,
        O: GraphOrientation,
    >(
        &self,
        graph: &G,
        amplitude: &S,
        // orientation: &OrientationData,
        canonize_esurface: &Option<ShiftRewrite>,
        orientations: &TiVec<OID, O>,
        cut_edges: &[EdgeIndex],
    ) -> Option<ParsingNet>
    where
        G: UltravioletGraph + AsRef<HedgeGraph<Edge, Vertex, H>>,
    {
        let (t, s) = self.local_3d.expr()?;
        let (t_int, _) = if let ApproxOp::Root = self.integrated_4d {
            (Atom::num(0), Sign::Positive)
        } else {
            self.integrated_4d
                .expr()
                .unwrap_or((Atom::num(0), Sign::Positive))
        };

        let CFFapprox::Dependent { t_arg, .. } = CFFapprox::dependent(
            graph.as_ref(),
            &self.subgraph,
            amplitude,
            canonize_esurface,
            orientations,
            cut_edges,
        ) else {
            unreachable!()
        };

        let finite = t_int
            .series(vakint_symbol!("ε"), Atom::Zero, 0.into(), true)
            .unwrap()
            .coefficient((0, 1).into())
            .replace(GS.m_uv_int)
            .with(GS.m_uv)
            .map_mink_dim(4);

        debug!(
            "Integrated 4d finite part: {:#}",
            finite.printer(LOGPRINTOPTS)
        );
        let reduced = amplitude.included().subtract(&self.subgraph.included());

        let mut cff = s * t - s * finite * t_arg.integrand;

        for (p, eid, _) in graph.as_ref().iter_edges_of(&reduced) {
            let eid = usize::from(eid) as i64;
            if p.is_paired() {
                cff = cff
                    .replace(function!(GS.energy, eid))
                    .with(function!(GS.ose, eid));
            }
        }

        cff = cff.replace(function!(GS.ose, W_.a__, W_.e_)).with(W_.e_);

        let mut resnum = graph.numerator(&reduced).get_single_atom().unwrap();

        let mut reps = Vec::new();
        for (p, eid, _) in graph.as_ref().iter_edges_of(&reduced) {
            if p.is_paired() {
                reps.push(GS.add_parametric_sign(eid));
            }
        }

        resnum = resnum.replace_multiple(&reps);
        resnum *= cff;
        resnum = resnum.wrap_color(GS.color_wrap);

        debug!("Integrand before parsing:{}", resnum.printer(LOGPRINTOPTS));
        let mut res = resnum.parse_into_net().unwrap();
        // println!("Final Integrand Net:{}", res.net.dot_pretty());
        res = res.replace_multiple(&reps);

        // .contract(())
        // .unwrap()
        // .state
        // .tensor
        // .scalar()
        // .expect("Expected a scalar value when contracting integrand tensor network");

        // res = res
        //     .replace(function!(
        //         GS.emr_vec,
        //         W_.a_,
        //         function!(AIND_SYMBOLS.cind, 0)
        //     ))
        //     .with(Atom::Zero)
        //     .replace(function!(GS.ose, W_.a_, function!(AIND_SYMBOLS.cind, 1)))
        //     .with(Atom::Zero)
        //     .replace(function!(GS.ose, W_.a_, function!(AIND_SYMBOLS.cind, 2)))
        //     .with(Atom::Zero)
        //     .replace(function!(GS.ose, W_.a_, function!(AIND_SYMBOLS.cind, 3)))
        //     .with(Atom::Zero)
        //     .replace(function!(GS.ose, W_.a_, function!(AIND_SYMBOLS.cind, 0)))
        //     .with(function!(GS.ose, W_.a_))
        //     .replace(function!(GS.emr_vec, W_.a__))
        //     .with(function!(GS.emr_mom, W_.a__));

        // debug!("final_cff {res:>}");
        Some(res)
    }

    // pub(crate) fn simple_expr(
    //     &self,
    //     graph: &UVGraph,
    //     amplitude: &InternalSubGraph,
    // ) -> Option<SerializableAtom> {
    //     let simple_approx = self.simple_approx.as_ref()?;

    //     Some((simple_approx.sign * simple_approx.expr(&amplitude.filter)).into())
    // }
}

impl ApproxOp {
    pub(crate) fn sign(&self) -> Option<Sign> {
        match self {
            ApproxOp::NotComputed => None,
            ApproxOp::Union { sign, .. } => Some(*sign),
            ApproxOp::Dependent { sign, .. } => Some(*sign),
            ApproxOp::Root => Some(Sign::Positive),
        }
    }

    pub(crate) fn expr(&self) -> Option<(Atom, Sign)> {
        match self {
            ApproxOp::NotComputed => None,
            ApproxOp::Union { t_args, sign, .. } => {
                //Never gets hit now
                let mut mul = Atom::num(1);
                for t in t_args {
                    mul = mul * &t.integrand;
                }
                Some((mul, *sign))
            }
            ApproxOp::Dependent { t_arg, sign, .. } => Some((t_arg.integrand.clone(), *sign)),
            ApproxOp::Root => Some((Atom::num(1), Sign::Positive)),
        }
    }

    pub(crate) fn union(dependent: &[&Approximation]) -> Option<Self> {
        let mut t_args = vec![];
        let mut subgraphs = vec![];

        let mut final_sign = Sign::Positive;
        for d in dependent {
            match &d.integrated_4d {
                ApproxOp::Dependent {
                    t_arg,
                    sign,
                    subgraph,
                } => {
                    t_args.push(t_arg.clone());
                    final_sign = final_sign * *sign;
                    subgraphs.push(subgraph.clone())
                }
                _ => return None,
            }
        }

        Some(Self::Union {
            t_args,
            sign: final_sign,
            subgraphs,
        })
    }

    pub(crate) fn is_computed(&self) -> bool {
        match self {
            ApproxOp::NotComputed => false,
            _ => true,
        }
    }
}
