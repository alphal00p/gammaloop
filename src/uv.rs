use std::{cmp::Ordering, collections::VecDeque, hash::Hash, ops::Deref};

use crate::{
    cff::{
        expression::{GraphOrientation, OrientationData, OrientationID},
        generation::{generate_uv_cff, ShiftRewrite},
    },
    model::ArcParticle,
    momentum::Sign,
    new_graph::{Edge, LMBext, LoopMomentumBasis, Vertex},
    numerator::aind::Aind,
    utils::{sign_atom, GS, W_},
};
use ahash::AHashSet;
use bitvec::vec::BitVec;
use eyre::eyre;
use idenso::metric::MS;
use pathfinding::prelude::BfsReachable;
use serde::{Deserialize, Serialize};
use spenso::{
    network::parsing::ShadowedStructure,
    shadowing::symbolica_utils::SerializableAtom,
    structure::{
        dimension::Dimension,
        representation::{Minkowski, RepName},
        HasStructure, NamedStructure, OrderedStructure, ToSymbolic,
    },
    tensors::parametric::{atomcore::PatternReplacement, ParamTensor},
};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder, Symbol},
    function,
    id::Replacement,
    parse,
    printer::PrintOptions,
    symbol,
};

use linnet::half_edge::{
    involution::{EdgeIndex, HedgePair, SignOrZero},
    subgraph::{Inclusion, InternalSubGraph, SubGraph, SubGraphOps},
    HedgeGraph,
};

use typed_index_collections::TiVec;
use vakint::{
    vakint_symbol, EvaluationOrder, LoopNormalizationFactor, Vakint, VakintExpression,
    VakintSettings,
};
// use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};

use crate::{
    graph::{BareEdge, BareGraph, BareVertex},
    model::normalise_complex,
};

// pub static VAKINT: LazyLock<Vakint> = LazyLock::new(|| {
//     Vakint::new(Some(VakintSettings {
//         evaluation_order: EvaluationOrder::alphaloop_only(),
//         integral_normalization_factor: LoopNormalizationFactor::MSbar,
//         run_time_decimal_precision: 16,
//         ..VakintSettings::default()
//     }))
//     .unwrap()
// });

pub trait UVE {
    fn mass_atom(&self) -> Atom;
}

impl UVE for Edge {
    fn mass_atom(&self) -> Atom {
        parse!(&self.particle.mass.name)
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct UVEdge {
    pub og_edge: usize,
    pub dod: i32,
    pub particle: ArcParticle,
    // pub prop:ArcPropagator.
    pub num: Atom,
    pub den: Atom,
}

impl UVE for UVEdge {
    fn mass_atom(&self) -> Atom {
        parse!(&self.particle.mass.name)
    }
}

impl UVEdge {
    pub fn from_edge(edge: &BareEdge, id: usize, bare_graph: &BareGraph) -> Self {
        let [colorless, _] = edge.color_separated_numerator(bare_graph, id);
        UVEdge {
            particle: edge.particle.clone(),
            og_edge: id,
            dod: edge.dod() as i32,
            num: normalise_complex(&colorless).into(),
            den: Atom::num(1).into(), // edge.full_den(bare_graph, index).into(),
        }
    }

    pub fn numerator_rank(&self) -> isize {
        0
    }
}

#[derive(Clone, Debug)]
pub struct UVNode {
    dod: i32,
    num: Atom,
    color: Option<ParamTensor>,
}

impl UVNode {
    pub fn from_vertex(vertex: &BareVertex, graph: &BareGraph) -> Self {
        if let Some((colorless, color)) = vertex.contracted_colorless_vertex_rule(graph) {
            UVNode {
                dod: vertex.dod() as i32,
                num: normalise_complex(&colorless).into(),
                color: Some(ParamTensor::param(
                    color.map_structure(OrderedStructure::from),
                )),
            }
        } else {
            UVNode {
                dod: vertex.dod() as i32,
                num: Atom::num(1).into(),
                color: None,
            }
        }
    }

    pub fn numerator_rank(&self) -> isize {
        0
    }
}

#[derive(Clone, Debug)]
pub struct UVGraph {
    hedge_graph: HedgeGraph<UVEdge, UVNode>,
    pub lmb: LoopMomentumBasis,
    lmb_replacement: Vec<Replacement>,
}

impl Deref for UVGraph {
    type Target = HedgeGraph<UVEdge, UVNode>;
    fn deref(&self) -> &Self::Target {
        &self.hedge_graph
    }
}

impl AsRef<HedgeGraph<UVEdge, UVNode>> for UVGraph {
    fn as_ref(&self) -> &HedgeGraph<UVEdge, UVNode> {
        &self.hedge_graph
    }
}
pub fn spenso_lor(
    tag: i32,
    ind: impl Into<Aind>,
    dim: impl Into<Dimension>,
) -> ShadowedStructure<Aind> {
    let mink = Minkowski {}.new_slot(dim, ind);
    NamedStructure::from_iter([mink], GS.emr_mom, Some(vec![Atom::num(tag)])).structure
}

pub fn spenso_lor_atom(tag: i32, ind: impl Into<Aind>, dim: impl Into<Dimension>) -> Atom {
    spenso_lor(tag, ind, dim).to_symbolic(None).unwrap()
}

pub struct Wood {
    poset: Poset<InternalSubGraph, ()>,
    additional_unions: SecondaryMap<PosetNode, Vec<PosetNode>>,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct IntegrandExpr {
    integrand: Atom,
    // add_arg: Option<Atom>,
}

impl IntegrandExpr {
    // pub(crate) fn from_subgraph<S: SubGraph>(subgraph: &S, graph: &UVGraph) -> Self {
    //     let num = graph.numerator(subgraph);

    //     let den = graph.denominator(subgraph);

    //     IntegrandExpr {
    //         integrand: num / den,
    //     }
    // }

    // pub fn numerator_only_subgraph<S: SubGraph>(subgraph: &S, graph: &UVGraph) -> Self {
    //     let num = graph.numerator(subgraph);

    //     IntegrandExpr { integrand: num }
    // }
}

// pub fn limit(&)

pub fn is_not_paired(pair: &HedgePair) -> bool {
    !pair.is_paired()
}

pub mod uv_graph;
pub use uv_graph::UltravioletGraph;

impl Wood {
    pub fn n_spinneys(&self) -> usize {
        self.poset.n_nodes()
    }

    pub fn from_spinneys<E, V, H, I: IntoIterator<Item = InternalSubGraph>>(
        s: I,
        graph: impl AsRef<HedgeGraph<E, V, H>>,
    ) -> Self {
        let mut poset = Poset::from_iter(s.into_iter().map(|s| (s, ())));
        let ref_graph = graph.as_ref();

        poset.invert();
        poset.compute_topological_order();

        let mut unions = SecondaryMap::new();

        for (i, sg) in poset.nodes.iter() {
            let cs = ref_graph.connected_components(&sg.data);

            if cs.len() > 1 {
                // sg is a disjoint union of spinneys (at the level of half-edges) (strongly disjoint)
                let mut union = vec![];

                for &c in sg.parents.iter() {
                    let mut is_in = 0;
                    for comp in &cs {
                        let comp =
                            InternalSubGraph::cleaned_filter_optimist(comp.clone(), ref_graph);
                        if comp == poset.nodes[c].data {
                            // find the components in the wood that this union is made of
                            union.push(c);
                            is_in += 1;
                        }
                    }

                    if is_in > 1 {
                        panic!("is in too many components")
                    }
                }

                unions.insert(i, union.clone());
                // sg.children = union;
            }
        }

        // let coverset = poset.to_cover_set();
        Wood {
            poset,
            additional_unions: unions,
        }
    }

    fn unfold_bfs<E, V, H, G>(
        &self,
        graph: &G,
        lmb: &LoopMomentumBasis,
        dag: &mut DAG<Approximation, DagNode, ()>,
        unions: &mut SecondaryMap<PosetNode, Option<Vec<(PosetNode, Option<DagNode>)>>>,
        root: PosetNode,
    ) -> DagNode
    where
        G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>,
    {
        // let graph = graph.as_ref();
        let mut search_front = VecDeque::new();

        let tree_root = dag.add_node(Approximation::new(
            self.poset.data(root).clone(),
            graph,
            lmb,
        ));
        search_front.push_front((root, tree_root));

        while let Some((node, parent)) = search_front.pop_front() {
            for c in &self.poset.nodes[node].children {
                if let Some(tagged_union) = unions.get_mut(*c) {
                    // Is this node a disjoint union of spinneys
                    if let Some(mut union) = tagged_union.take() {
                        let mut all_supplied = true;
                        for (p, d) in &mut union {
                            if self.poset.nodes[*p].data == self.poset.nodes[node].data {
                                *d = Some(parent);
                            }
                            if d.is_none() {
                                all_supplied = false;
                            }
                        }
                        if all_supplied {
                            let child = dag.add_node(Approximation::new(
                                self.poset.data(*c).clone(),
                                graph,
                                lmb,
                            ));
                            for (_, d) in union {
                                dag.add_edge(d.unwrap(), child);
                            }
                            search_front.push_front((*c, child));
                        } else {
                            *tagged_union = Some(union);
                        }
                    }
                } else {
                    let child =
                        dag.add_node(Approximation::new(self.poset.data(*c).clone(), graph, lmb));
                    dag.add_edge(parent, child);
                    search_front.push_front((*c, child));
                }
            }
        }
        tree_root
    }

    pub fn unfold<E, V, H, G>(&self, graph: &G, lmb: &LoopMomentumBasis) -> Forest
    where
        G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>,
    {
        let mut dag: DAG<Approximation, DagNode, ()> = DAG::new();

        let root = self.poset.minimum().unwrap();

        let mut unions = SecondaryMap::new();

        for (p, u) in self.additional_unions.iter() {
            let union: Vec<(PosetNode, Option<DagNode>)> = u.iter().map(|i| (*i, None)).collect();
            unions.insert(p, Some(union));
        }

        let _ = self.unfold_bfs(graph, lmb, &mut dag, &mut unions, root);

        Forest { dag }
    }

    pub fn dot(&self, graph: &impl UltravioletGraph) -> String {
        let shift = self.poset.shift();
        self.poset.to_dot_impl(&|n| {
            format!(
                "label={}, dod={},topo_order = {}",
                n.dot_id(shift),
                graph.dod(&n.data),
                // graph.as_ref().count_internal_edges(&n.data),
                n.order.unwrap()
            )
        })
    }

    pub fn dot_spinneys<E, V, H, G>(&self, graph: &G)
    where
        G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>,
    {
        for s in self.poset.node_values() {
            println!(
                "found {} loop spinney with dod {}:{} ",
                graph.n_loops(s),
                graph.dod(s),
                graph.as_ref().dot(s)
            );
        }
    }
}

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
    fn subgraph_shadow(graph: &BitVec, subgraph: &InternalSubGraph) -> Symbol {
        symbol!(&format!(
            "S_{}⊛{}",
            graph.string_label(),
            subgraph.string_label()
        ))
    }

    pub fn expr(&self, bigger_graph: &BitVec) -> Atom {
        let reduced = Atom::var(Self::subgraph_shadow(bigger_graph, &self.graph));
        let mut mul = Atom::num(1);
        for i in &self.t_args {
            mul = mul * i
        }
        reduced * mul
    }

    pub fn t_op(&self, bigger_graph: &BitVec) -> Atom {
        function!(GS.top, self.expr(bigger_graph))
    }

    pub fn root(subgraph: InternalSubGraph) -> Self {
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

    pub fn dependent(&self, bigger_graph: InternalSubGraph) -> Self {
        Self {
            t_args: vec![self.t_op(&bigger_graph.filter)],
            sign: -self.sign,
            graph: bigger_graph,
        }
    }

    pub fn union<'a>(
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

#[derive(Clone)]

pub struct Approximation {
    // The union of all spinneys, remaining graph is full graph minus subgraph
    subgraph: InternalSubGraph,
    dod: i32,
    lmb: LoopMomentumBasis,
    pub local_3d: CFFapprox,     //3d denoms
    pub integrated_4d: ApproxOp, //4d
    pub simple_approx: Option<SimpleApprox>,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum CFFapprox {
    NotComputed,
    Dependent { sign: Sign, t_arg: IntegrandExpr },
}

impl CFFapprox {
    pub fn expr(&self) -> Option<(Atom, Sign)> {
        match self {
            CFFapprox::NotComputed => None,
            CFFapprox::Dependent { sign, t_arg } => Some((t_arg.integrand.clone(), *sign)),
        }
    }
    pub fn dependent<
        E,
        V,
        H,
        S: SubGraph,
        SS: SubGraph,
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

        println!("{:?}", contract_edges);

        for (oid, o) in orientations.iter_enumerated() {
            cff_sum += o.orientation_delta()
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
    pub fn root<E, V, H, S: SubGraph, OID: OrientationID, O: GraphOrientation>(
        graph: &HedgeGraph<E, V, H>,
        amplitude_subgraph: &S,
        canonize_esurface: &Option<ShiftRewrite>,
        orientations: &TiVec<OID, O>,
        cut_edges: &[EdgeIndex],
    ) -> CFFapprox {
        Self::dependent(
            graph,
            &graph.empty_subgraph::<BitVec>(),
            amplitude_subgraph,
            canonize_esurface,
            orientations,
            cut_edges,
        )
    }
}
impl Approximation {
    pub fn root<E, V, H, S: SubGraph, OID: OrientationID, O: GraphOrientation>(
        &mut self,
        graph: &HedgeGraph<E, V, H>,
        amplitude_subgraph: &S,
        canonize_esurface: &Option<ShiftRewrite>,
        orientations: &TiVec<OID, O>,
        cut_edges: &[EdgeIndex],
    ) {
        self.local_3d = CFFapprox::root(
            graph,
            amplitude_subgraph,
            canonize_esurface,
            orientations,
            cut_edges,
        );
        self.integrated_4d = ApproxOp::Root;
    }
    pub fn simplify_notation(expr: &Atom) -> Atom {
        let replacements = [(function!(GS.den, W_.a_, W_.x_), Atom::var(W_.x_))];

        let reps: Vec<_> = replacements
            .into_iter()
            .map(|(a, b)| Replacement::new(a.to_pattern(), b))
            .collect();

        expr.replace_multiple_repeat(&reps)
    }

    pub fn new<G, E, V, H>(
        spinney: InternalSubGraph,
        graph: &G,
        lmb: &LoopMomentumBasis,
    ) -> Approximation
    where
        G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>,
    {
        let lmb = graph.compatible_sub_lmb(&spinney, lmb);
        // println!("//lmb for spinney \n{}", graph.dot_lmb(&spinney, &lmb));
        Approximation {
            dod: graph.dod(&spinney),
            subgraph: spinney,
            lmb,
            simple_approx: None,
            local_3d: CFFapprox::NotComputed,
            integrated_4d: ApproxOp::NotComputed,
        }
    }

    pub fn reduced_graph(&self, subgraph: &InternalSubGraph) -> InternalSubGraph {
        self.subgraph.subtract(subgraph)
    }

    pub fn structure_approximate(&self, dependent: &Self) -> SimpleApprox {
        dependent
            .simple_approx
            .as_ref()
            .unwrap()
            .dependent(self.subgraph.clone())
    }

    pub fn integrated_4d<
        E: UVE,
        V,
        H,
        G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>,
        S: SubGraph,
    >(
        &self,
        dependent: &Self,
        uv_graph: &G,
        amplitude_subgraph: &S,
    ) -> ApproxOp {
        let graph = uv_graph.as_ref();
        let reduced = self.subgraph.subtract(&dependent.subgraph);

        let Some((inner_t, sign)) = dependent.integrated_4d.expr() else {
            return ApproxOp::NotComputed;
        };
        let t_arg = uv_graph
            .numerator(&reduced, false)
            .color_simplify()
            .gamma_simplify()
            .get_single_atom()
            .unwrap()
            / uv_graph.denominator(&reduced);

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

        let soft_ct = graph.full_crown(&self.subgraph).count_ones() == 2 && self.dod > 0;

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

        // strip the momentum wrapper from the denominator
        a = a
            .replace(function!(
                GS.den,
                W_.prop_,
                function!(GS.emr_mom, W_.prop_, W_.mom_),
                W_.x__
            ))
            .with(function!(GS.den, W_.prop_, W_.mom_, W_.x__));

        println!("Expanded: {:>}", a.expand());

        let mut integrand_vakint = a.expand();
        let mut propagator_id = 1;

        let vakint = Vakint::new(Some(VakintSettings {
            allow_unknown_integrals: false,
            evaluation_order: EvaluationOrder::alphaloop_only(),
            integral_normalization_factor: LoopNormalizationFactor::MSbar,
            run_time_decimal_precision: 32,
            number_of_terms_in_epsilon_expansion: uv_graph.n_loops(amplitude_subgraph) as i64 + 1,
            temporary_directory: Some("./form".into()),
            mu_r_sq_symbol: GS.mu_r_sq.get_name().to_string(),
            ..VakintSettings::default()
        }))
        .unwrap();

        let vk_prop = vakint_symbol!("prop");
        let vk_edge = vakint_symbol!("edge");
        let vk_mom = vakint_symbol!("k");
        let vk_topo = vakint_symbol!("topo");
        let vk_metric = vakint_symbol!("g");

        for (pair, index, _data) in graph.iter_edges_of(&reduced) {
            if let HedgePair::Paired { source, sink } = pair {
                integrand_vakint = integrand_vakint
                    .replace(function!(
                        GS.den,
                        usize::from(index) as i64,
                        W_.mom_,
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
                        GS.m_uv,
                        1
                    ))
                    .replace(function!(vk_prop, W_.x___, 1).pow(Atom::var(W_.e_)))
                    .with(function!(vk_prop, W_.x___, -Atom::var(W_.e_)));
                propagator_id += 1;
            }
        }

        // shrink vertices of the subgraph
        for (pair, _index, _data) in graph.iter_edges_of(&dependent.subgraph) {
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

        // rewrite numerator
        // linearize the numerator first
        integrand_vakint = integrand_vakint
            .replace(function!(GS.emr_mom, W_.prop_, W_.mom_, W_.x_))
            .with(function!(MS.dot, W_.mom_, W_.x_))
            .replace(function!(MS.dot, W_.mom_, W_.x_))
            .with(function!(GS.emr_mom, W_.mom_, W_.x_));

        for (i, l) in self.lmb.loop_edges.iter().enumerate() {
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

        println!("Integrand vakint: {}", integrand_vakint);

        let vakint_expr = VakintExpression::try_from(integrand_vakint.clone()).unwrap();
        println!("\nVakint expression:\n{:#}", vakint_expr);

        let mut res = vakint.evaluate(integrand_vakint.as_view()).unwrap();

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

    pub fn compute<E: UVE, V, H, G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>, S: SubGraph>(
        &mut self,
        graph: &G,
        amplitude_subgraph: &S,
        dependent: &Self,
    ) {
        self.integrated_4d = self.integrated_4d(dependent, graph, amplitude_subgraph);
    }

    pub fn compute_3d<
        E: UVE,
        V,
        H,
        G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>,
        S: SubGraph,
        OID: OrientationID,
        O: GraphOrientation,
    >(
        &mut self,
        graph: &G,
        amplitude_subgraph: &S,
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
                .expect("Should have computed the dependent integrated 4d")
        };

        let CFFapprox::Dependent { t_arg, .. } = CFFapprox::dependent(
            graph.as_ref(),
            &dependent.subgraph,
            amplitude_subgraph,
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
            .with(GS.m_uv);

        println!("Integrated 4d finite part: {}", finite);

        // TODO: multiply by the number of orientations or only apply the counterterm to
        // one orientation
        // subtract integrated CT
        sum_3d -= self.local_3d(dependent, graph, finite * t_arg.integrand);

        self.local_3d = CFFapprox::Dependent {
            sign: -sign,
            t_arg: IntegrandExpr { integrand: sum_3d },
        }
    }

    pub fn local_3d<E: UVE, V, H, G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>>(
        &self,
        dependent: &Self,
        uv_graph: &G,
        mut cff: Atom,
    ) -> Atom {
        let graph = uv_graph.as_ref();
        let reduced = self.subgraph.subtract(&dependent.subgraph);

        //println!("CFF: {}", cff);

        // add data for OSE computation and add an explicit sqrt
        for (p, eid, e) in graph.iter_edges_of(&self.subgraph) {
            let eid = usize::from(eid) as i64;
            if p.is_paired() {
                // set energies from inner_t on-shell
                cff = cff
                    .replace(function!(GS.energy, eid))
                    .with(function!(GS.ose, eid));

                let e_mass = e.data.mass_atom();
                cff = cff.replace(function!(GS.ose, eid)).with(
                    function!(
                        GS.ose,
                        eid,
                        function!(GS.emr_vec, eid),
                        &e_mass * &e_mass,
                        -function!(
                            MS.dot,
                            function!(GS.emr_vec, eid),
                            function!(GS.emr_vec, eid)
                        ) + &e_mass * &e_mass
                    )
                    .npow((1, 2)),
                );
            }
        }

        let mut atomarg = cff
            * uv_graph
                .numerator(&reduced, false)
                .color_simplify()
                .gamma_simplify()
                .get_single_atom()
                .unwrap();

        // println!(
        //     "Expand-prerep {} with dod={} in {:?}",
        //     atomarg, self.dod, self.lmb.ext_edges
        // );

        // split numerator momenta into OSEs and spatial parts
        for (p, eid, e) in graph.iter_edges_of(&self.subgraph) {
            let eidc = usize::from(eid) as i64;
            if p.is_paired() {
                let e_mass = e.data.mass_atom();

                atomarg = atomarg
                    .replace(function!(GS.emr_mom, eidc, W_.x___))
                    .with_map(move |m| {
                        let index = m.get(W_.x___).unwrap().to_atom();

                        let sign = sign_atom(eid);

                        function!(
                            GS.ose,
                            eidc,
                            function!(GS.emr_vec, eidc),
                            &e_mass * &e_mass,
                            -function!(
                                MS.dot,
                                function!(GS.emr_vec, eidc),
                                function!(GS.emr_vec, eidc)
                            ) + &e_mass * &e_mass,
                            index
                        )
                        .npow((1, 2))
                            * sign
                            + function!(GS.emr_vec, eidc, index)
                    });
            }
        }

        // only apply replacements for edges in the reduced graph
        let mom_reps = graph.uv_spatial_wrapped_replacement(&reduced, &self.lmb, &[W_.x___]);

        println!("Reps:");
        for r in &mom_reps {
            println!("{r}");
        }

        atomarg = atomarg.replace_multiple(&mom_reps);

        // rescale the loop momenta in the whole subgraph, including previously expanded cycles
        for e in &self.lmb.loop_edges {
            println!("Rescale {}", e);
            atomarg = atomarg
                .replace(function!(GS.emr_vec, usize::from(*e) as i64, W_.x___))
                .with(function!(GS.emr_vec, usize::from(*e) as i64, W_.x___) * GS.rescale);
        }

        // println!(
        //     "Expand {} with dod={} in {:?}",
        //     atomarg, self.dod, self.lmb.ext_edges
        // );

        let soft_ct = graph.full_crown(&self.subgraph).count_ones() == 2 && self.dod > 0;

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
                    .replace(function!(GS.ose, eid, W_.mom_, W_.mass_, W_.prop_, W_.a___))
                    .with(
                        function!(
                            GS.ose,
                            eid,
                            W_.mom_,
                            GS.m_uv * GS.m_uv,
                            (GS.m_uv * GS.m_uv * GS.rescale * GS.rescale + W_.prop_
                                - GS.m_uv * GS.m_uv)
                                / GS.rescale
                                / GS.rescale,
                            W_.a___
                        ) * GS.rescale
                            * GS.rescale,
                    )
                    .replace(function!(GS.ose, eid, W_.mom_, W_.a___))
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

        //println!("atomarg:{}", atomarg);

        println!("Series expanding");

        let mut a = atomarg
            .series(GS.rescale, Atom::Zero, self.dod.into(), true)
            .unwrap()
            .to_atom()
            .replace(parse!("der(0,0,0,1, OSE(y__))"))
            .with(Atom::num(1))
            .replace(parse!("der(0,0,0,1,0, OSE(y__))"))
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

        //println!("Expanded: {:>}", a.expand());

        a
    }

    pub fn final_cff<G, H, S: SubGraph<Base = BitVec>, OID: OrientationID, O: GraphOrientation>(
        &self,
        graph: &G,
        amplitude: &S,
        orientation: &OrientationData,
        canonize_esurface: &Option<ShiftRewrite>,
        orientations: &TiVec<OID, O>,
        cut_edges: &[EdgeIndex],
    ) -> Option<Atom>
    where
        G: UltravioletGraph + AsRef<HedgeGraph<Edge, Vertex, H>>,
    {
        let (t, s) = self.local_3d.expr()?;
        let (t_int, _) = if let ApproxOp::Root = self.integrated_4d {
            (Atom::num(0), Sign::Positive)
        } else {
            self.integrated_4d.expr()?
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
            .with(GS.m_uv);

        println!("Integrated 4d finite part: {}", finite);

        let reduced = amplitude.included().subtract(&self.subgraph.included());

        let mut cff = s * orientation.select(t) - s * finite * orientation.select(t_arg.integrand);

        for (p, eid, e) in graph.as_ref().iter_edges_of(&reduced) {
            let eid = usize::from(eid) as i64;
            if p.is_paired() {
                cff = cff
                    .replace(function!(GS.energy, eid))
                    .with(function!(GS.ose, eid));

                let e_mass = parse!(&e.data.particle.mass.name);
                cff = cff.replace(function!(GS.ose, eid)).with(
                    function!(
                        GS.ose,
                        eid,
                        function!(GS.emr_vec, eid),
                        &e_mass * &e_mass,
                        -function!(
                            MS.dot,
                            function!(GS.emr_vec, eid),
                            function!(GS.emr_vec, eid)
                        ) + &e_mass * &e_mass
                    )
                    .npow((1, 2)),
                );
            }
        }

        let mut res = cff
            * graph
                .numerator(&reduced, true)
                .color_simplify()
                .gamma_simplify()
                .get_single_atom()
                .unwrap();

        // set the momenta flowing through the reduced graph edges to the identity wrt the supergraph
        for (p, eid, e) in graph.as_ref().iter_edges_of(&reduced) {
            let eidc = usize::from(eid) as i64;
            if p.is_paired() {
                let e_mass = parse!(&e.data.particle.mass.name);
                let orientation = orientation.orientation.clone();
                res = res
                    .replace(function!(GS.emr_mom, eidc, W_.x___))
                    .with_map(move |m| {
                        let index = m.get(W_.x___).unwrap().to_atom();

                        let sign = SignOrZero::from(orientation[eid].clone()) * 1;

                        function!(
                            GS.ose,
                            eidc,
                            function!(GS.emr_vec, eidc),
                            &e_mass * &e_mass,
                            -function!(
                                MS.dot,
                                function!(GS.emr_vec, eidc),
                                function!(GS.emr_vec, eidc)
                            ) + &e_mass * &e_mass,
                            index
                        )
                        .npow((1, 2))
                            * sign
                            + function!(GS.emr_vec, eidc, index)
                    });
            }
        }

        Some(res)
    }

    pub fn simple_expr(
        &self,
        graph: &UVGraph,
        amplitude: &InternalSubGraph,
    ) -> Option<SerializableAtom> {
        let simple_approx = self.simple_approx.as_ref()?;

        Some((simple_approx.sign * simple_approx.expr(&amplitude.filter)).into())
    }
}

impl ApproxOp {
    pub fn sign(&self) -> Option<Sign> {
        match self {
            ApproxOp::NotComputed => None,
            ApproxOp::Union { sign, .. } => Some(*sign),
            ApproxOp::Dependent { sign, .. } => Some(*sign),
            ApproxOp::Root => Some(Sign::Positive),
        }
    }

    pub fn expr(&self) -> Option<(Atom, Sign)> {
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

    pub fn union(dependent: &[&Approximation]) -> Option<Self> {
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

    pub fn is_computed(&self) -> bool {
        match self {
            ApproxOp::NotComputed => false,
            _ => true,
        }
    }
}

pub struct Forest {
    dag: DAG<Approximation, DagNode, ()>,
}

impl Forest {
    pub fn n_terms(&self) -> usize {
        self.dag.nodes.len()
    }

    pub fn compute<
        E: UVE,
        V,
        H,
        G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>,
        S: SubGraph,
        OID: OrientationID,
        O: GraphOrientation,
    >(
        &mut self,
        graph: &G,
        amplitude_subgraph: &S,
        orientations: &TiVec<OID, O>,
        canonize_esurface: &Option<ShiftRewrite>,
        cut_edges: &[EdgeIndex],
    ) {
        let order = self.dag.compute_topological_order();

        for n in order {
            match self.dag.nodes[n].parents.len() {
                0 => {
                    self.dag.nodes[n].data.root(
                        graph.as_ref(),
                        amplitude_subgraph,
                        canonize_esurface,
                        orientations,
                        cut_edges,
                    );
                }
                1 => {
                    let [ref mut current, ref mut parent] = &mut self
                        .dag
                        .nodes
                        .get_disjoint_mut([n, self.dag.nodes[n].parents[0]])
                        .unwrap();

                    assert!(matches!(parent.data.local_3d, CFFapprox::Dependent { .. }));
                    current
                        .data
                        .compute(graph, amplitude_subgraph, &parent.data);
                    current.data.compute_3d(
                        graph,
                        amplitude_subgraph,
                        canonize_esurface,
                        orientations,
                        cut_edges,
                        &parent.data,
                    );
                }
                _ => {
                    unimplemented!("Union not implemented");
                    // let (local_4d, simple, local_3d) = {
                    //     let current = &self.dag.nodes[n];
                    //     let mut dependents =
                    //         current.parents.iter().map(|p| &self.dag.nodes[*p].data);

                    //     let mut dep_vec: Vec<&Approximation> = Vec::new();
                    //     dep_vec.push(dependents.next().unwrap());

                    //     let Some((mut approx, _)) = dep_vec[0].local_3d.expr() else {
                    //         panic!("local 3d not computed for parents")
                    //     };

                    //     let mut iterative_approx = dep_vec[0].clone();

                    //     for p in dependents {
                    //         dep_vec.push(p);

                    //         approx = p.local_3d(&iterative_approx, graph, approx);

                    //         iterative_approx = p.clone();
                    //     }
                    //     (
                    //         ApproxOp::union(&dep_vec).unwrap(),
                    //         SimpleApprox::union(
                    //             self.dag.nodes[n].data.subgraph.clone(),
                    //             dep_vec.iter().map(|s| s.simple_approx.as_ref().unwrap()),
                    //         ),
                    //         iterative_approx.local_3d,
                    //     )
                    // };

                    // self.dag.nodes[n].data.local_4d = local_4d;
                    // self.dag.nodes[n].data.simple_approx = Some(simple);
                    // self.dag.nodes[n].data.local_3d = local_3d;
                }
            }
        }
    }

    pub fn local_expr<G, H, S: SubGraph<Base = BitVec>, OID: OrientationID, O: GraphOrientation>(
        &self,
        graph: &G,
        amplitude: &S,
        orientation: &OrientationData,
        canonize_esurface: &Option<ShiftRewrite>,
        orientations: &TiVec<OID, O>,
        cut_edges: &[EdgeIndex],
    ) -> Atom
    where
        G: UltravioletGraph + AsRef<HedgeGraph<Edge, Vertex, H>>,
    {
        let mut sum = Atom::Zero;

        for (_, n) in &self.dag.nodes {
            sum += n
                .data
                .final_cff(
                    graph,
                    amplitude,
                    orientation,
                    canonize_esurface,
                    orientations,
                    cut_edges,
                )
                .unwrap();
        }

        sum
    }

    pub fn simple_expr(
        &self,
        graph: &UVGraph,
        amplitude: &InternalSubGraph,
    ) -> Option<SerializableAtom> {
        let mut sum = Atom::num(0).into();
        for (_, n) in &self.dag.nodes {
            sum = sum + n.data.simple_expr(graph, amplitude)?;
        }

        Some(sum)
    }

    pub fn graphs(&self) -> String {
        self.dag
            .to_dot_impl(&|n| format!("label=S_{}", n.data.subgraph.string_label()))
    }

    pub fn show_structure(&self, graph: &UVGraph, amplitude: &InternalSubGraph) -> Option<String> {
        let mut out = String::new();

        out.push_str(
            &self
                .simple_expr(graph, amplitude)?
                .0
                .printer(PrintOptions {
                    terms_on_new_line: true,
                    ..Default::default()
                })
                .to_string(),
        );
        Some(out)
    }
}

impl Wood {
    pub fn show_graphs<H, G>(&self, graph: &G) -> String
    where
        G: UltravioletGraph + AsRef<HedgeGraph<Edge, Vertex, H>>,
    {
        let mut out = String::new();
        out.push_str("Poset structure:\n");
        out.push_str(&self.poset.dot_structure());

        out.push_str("Graphs:\n");
        for (k, n) in self.poset.nodes.iter() {
            out.push_str(&graph.as_ref().dot_impl(
                &n.data,
                format!(
                    "dod={};nodeid ={};\n",
                    graph.dod(&n.data),
                    self.poset.dot_id(k)
                ),
                &|_h| None,
                &|e| Some(format!("dod={}", e.dod)),
                &|n| Some(format!("dod={}", n.dod)),
            ));
            out.push('\n');
        }
        out
    }
}

use slotmap::{new_key_type, Key, SecondaryMap, SlotMap};
use std::collections::HashSet;

// Define a new key type for the Poset
new_key_type! {
    pub struct PosetNode;
    pub struct DagNode;
    pub struct UnfoldedWoodNode;
    pub struct CoverSetNode;
}

/// Trait to define DOT attributes for node data.
pub trait DotAttrs {
    fn dot_attrs(&self) -> String;
}

/// A node in the poset, storing generic data, edges to child nodes, and references to parent odes.
#[derive(Debug, Eq)]
pub struct SlotNode<T, R: Key> {
    data: T,
    order: Option<u64>,
    id: R,
    parents: Vec<R>,  // References to parent nodes by key
    children: Vec<R>, // Edges to child nodes by key
}

impl<T: Hash, R: Key> Hash for SlotNode<T, R> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.data.hash(state);
        self.id.hash(state);
    }
}

impl<T: PartialEq, R: Key> PartialEq for SlotNode<T, R> {
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data && self.id == other.id
    }
}

impl<T: PartialEq, R: Key> PartialOrd for SlotNode<T, R> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.order.unwrap().partial_cmp(&other.order.unwrap())
    }
}

impl<T: PartialEq + Eq, R: Key> Ord for SlotNode<T, R> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.order.unwrap().cmp(&other.order.unwrap())
    }
}

impl<T, R: Key> SlotNode<T, R> {
    pub fn dot_id(&self, shift: u64) -> String {
        base62::encode(self.id() - shift)
    }

    pub fn to_topo_ordered(&self) -> Result<TopoOrdered<T>>
    where
        T: Clone,
    {
        Ok(TopoOrdered::new(
            self.data.clone(),
            self.order
                .ok_or_else(|| eyre!("Node has no topological order"))?,
        ))
    }

    pub fn in_degree(&self) -> usize {
        self.parents.len()
    }

    pub fn out_degree(&self) -> usize {
        self.children.len()
    }

    pub fn id(&self) -> u64 {
        self.id.data().as_ffi()
    }

    pub fn add_parent(&mut self, parent: R) {
        self.parents.push(parent);
    }

    pub fn add_child(&mut self, child: R) {
        self.children.push(child);
    }

    pub fn remove_child(&mut self, child: R) {
        self.children.retain(|&c| c != child);
    }

    pub fn remove_parent(&mut self, parent: R) {
        self.parents.retain(|&c| c != parent);
    }

    pub fn is_parent_of(&self, child: R) -> bool {
        self.children.contains(&child)
    }

    pub fn new(data: T, id: R) -> Self {
        SlotNode {
            data,
            id,
            order: None,
            parents: Vec::new(),
            children: Vec::new(),
        }
    }

    pub fn compare_data(&self, other: &Self) -> Option<std::cmp::Ordering>
    where
        T: PartialOrd,
    {
        self.data.partial_cmp(&other.data)
    }
}

/// A partially ordered set (poset) that can be built from an iterator and a slotmap.
pub struct DAG<T, R: Key, D = ()> {
    pub nodes: SlotMap<R, SlotNode<T, R>>,
    associated_data: SecondaryMap<R, D>,
}

impl<T, R: Key, D> DAG<T, R, D> {
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }
}

pub type Poset<T, D> = DAG<T, PosetNode, D>;
pub type HasseDiagram<T, D> = DAG<T, CoverSetNode, D>;

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct TopoOrdered<T> {
    pub data: T,
    order: u64,
}

impl<T> TopoOrdered<T> {
    pub fn new(data: T, order: u64) -> Self {
        TopoOrdered { data, order }
    }
}

#[allow(clippy::non_canonical_partial_ord_impl)]
impl<T> PartialOrd for TopoOrdered<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.order.partial_cmp(&other.order)
    }
}

impl<T> PartialEq for TopoOrdered<T> {
    fn eq(&self, other: &Self) -> bool {
        self.order == other.order
    }
}

impl<T> Eq for TopoOrdered<T> {}

impl<T> Ord for TopoOrdered<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.order.cmp(&other.order)
    }
}

use color_eyre::Result;

impl<T, R: Key, D> DAG<T, R, D> {
    pub fn new() -> Self {
        DAG {
            nodes: SlotMap::with_key(),
            associated_data: SecondaryMap::new(),
        }
    }

    pub fn children(&self, key: R) -> impl Iterator<Item = R> + '_ {
        self.nodes.get(key).unwrap().children.iter().copied()
    }

    pub fn dot_id(&self, key: R) -> String {
        self.nodes.get(key).unwrap().dot_id(self.shift())
    }

    pub fn node_values(&self) -> impl Iterator<Item = &T> {
        self.nodes.values().map(|node| &node.data)
    }

    pub fn add_edge(&mut self, from: R, to: R) {
        let from_node = self.nodes.get_mut(from).unwrap();
        from_node.add_child(to);
        let to_node = self.nodes.get_mut(to).unwrap();
        to_node.add_parent(from);
    }

    pub fn add_edge_if_new(&mut self, from: R, to: R) {
        if !self.nodes.get(from).unwrap().is_parent_of(to) {
            self.add_edge(from, to);
        }
    }

    pub fn remove_edge(&mut self, from: R, to: R) {
        let from_node = self.nodes.get_mut(from).unwrap();
        from_node.remove_child(to);
        let to_node = self.nodes.get_mut(to).unwrap();
        to_node.remove_parent(from);
    }

    pub fn bfs_reach<'a>(
        &'a self,
        start: &'a R,
    ) -> BfsReachable<&'a R, impl FnMut(&'a &'a R) -> &'a [R]> {
        pathfinding::directed::bfs::bfs_reach(start, |&s| self.succesors(*s))
    }

    pub fn maximum(&self) -> Option<R>
    where
        T: Eq,
    {
        Some(self.nodes.values().max()?.id)
    }

    pub fn minimum(&self) -> Option<R>
    where
        T: Eq,
    {
        Some(self.nodes.values().min()?.id)
    }

    pub fn succesors(&self, node_key: R) -> &[R] {
        &self.nodes.get(node_key).unwrap().children
    }

    /// Returns an iterator over all paths starting from the root node, traversed in BFS order.
    pub fn bfs_paths(&self) -> BfsPaths<T, R>
    where
        T: Eq,
    {
        let mut queue = VecDeque::new();
        queue.push_back(vec![self.minimum().unwrap()]);
        BfsPaths {
            queue,
            visited: HashSet::new(),
            nodes: &self.nodes,
        }
    }

    pub fn invert(&mut self) {
        self.nodes.iter_mut().for_each(|(_, node)| {
            std::mem::swap(&mut node.children, &mut node.parents);
        });
    }

    pub fn bfs_paths_inv(&self) -> BfsPaths<T, R> {
        let mut queue = VecDeque::new();
        let mut maximal_elements = Vec::new();
        let mut has_incoming = HashSet::new();

        for node in self.nodes.values() {
            for &parent in node.parents.iter() {
                has_incoming.insert(parent);
            }
        }

        for (key, _node) in self.nodes.iter() {
            if !has_incoming.contains(&key) {
                maximal_elements.push(key);
            }
        }

        for &max in &maximal_elements {
            queue.push_back(vec![max]);
        }

        BfsPaths {
            queue,
            visited: HashSet::new(),
            nodes: &self.nodes,
        }
    }

    pub fn data(&self, key: R) -> &T {
        &self.nodes.get(key).unwrap().data
    }

    fn shift(&self) -> u64 {
        self.nodes.iter().next().unwrap().1.id()
    }

    pub fn to_dot(&self, label: &impl Fn(&T) -> String) -> String {
        self.to_dot_impl(&|node| label(&node.data))
    }

    fn to_dot_impl(&self, label: &impl Fn(&SlotNode<T, R>) -> String) -> String {
        let mut dot = String::new();
        dot.push_str("digraph Poset {\n");
        dot.push_str("    node [shape=circle];\n");

        let shift = self.shift();

        for node in self.nodes.values() {
            let node_id = node.dot_id(shift);
            dot.push_str(&format!("n{} [{}];\n", node_id, label(node)));
            for &child in node.children.iter() {
                dot.push_str(&format!(
                    "n{} -> n{};\n",
                    node_id,
                    self.nodes.get(child).unwrap().dot_id(shift)
                ));
            }
        }

        dot.push_str("}\n");
        dot
    }

    pub fn dot_structure(&self) -> String {
        let shift = self.shift();
        self.to_dot_impl(&|n| format!("label={}", n.dot_id(shift)))
    }

    pub fn add_node(&mut self, data: T) -> R {
        self.nodes.insert_with_key(|key| SlotNode::new(data, key))
    }

    /// Returns all descendants of the node, used for propagating transitive relations.
    fn get_all_descendants(&self, node_key: R) -> HashSet<R> {
        let mut descendants = HashSet::new();
        let mut stack = vec![node_key];
        while let Some(current_key) = stack.pop() {
            let current_node = self.nodes.get(current_key).unwrap();
            for &child_key in current_node.children.iter() {
                if descendants.insert(child_key) {
                    stack.push(child_key);
                }
            }
        }
        descendants
    }

    /// Returns all ancestors of the node, used for propagating transitive relations.
    fn get_all_ancestors(&self, node_key: R) -> HashSet<R> {
        let mut ancestors = HashSet::new();
        let mut stack = vec![node_key];
        while let Some(current_key) = stack.pop() {
            let current_node = self.nodes.get(current_key).unwrap();
            for &parent_key in current_node.parents.iter() {
                if ancestors.insert(parent_key) {
                    stack.push(parent_key);
                }
            }
        }
        ancestors
    }

    pub fn transitive_edges(&self, a: R) -> Vec<(R, R)> {
        let mut edges = Vec::new();

        let children: AHashSet<_> = self
            .nodes
            .get(a)
            .unwrap()
            .children
            .iter()
            .cloned()
            .collect();

        for &child in &children {
            for descendant in self.get_all_descendants(child) {
                if children.contains(&descendant) {
                    edges.push((a, descendant));
                }
            }
        }
        edges
    }

    pub fn in_degree(&self, key: R) -> usize {
        self.nodes.get(key).unwrap().in_degree()
    }

    pub fn compute_topological_order(&mut self) -> Vec<R> {
        // Initialize queue with nodes having in-degree zero
        let mut queue = VecDeque::new(); //S in the wikipedia article

        let mut indegrees = SecondaryMap::new();

        for (key, node) in self.nodes.iter() {
            indegrees.insert(key, node.in_degree());
            if node.in_degree() == 0 {
                queue.push_back(key);
            }
        }

        let mut order = vec![];
        while let Some(node_key) = queue.pop_front() {
            // Assign the order number to the node
            if let Some(node) = self.nodes.get_mut(node_key) {
                node.order = Some(order.len() as u64);
                order.push(node_key);
            }

            // For each child, decrement its in-degree
            for &child_key in &self.nodes.get(node_key).unwrap().children {
                indegrees[child_key] -= 1;
                if indegrees[child_key] == 0 {
                    queue.push_back(child_key);
                }
            }
        }

        // Optional: Check if graph has cycles
        if order.len() != self.nodes.len() {
            panic!("The graph contains a cycle!");
        }
        order
    }
}

impl<T, R: Key, D> Default for DAG<T, R, D> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T, D> Poset<T, D> {
    pub fn poset_family(&self, data: &T) -> [Vec<PosetNode>; 2]
    where
        T: PartialOrd,
    {
        let mut parents = Vec::new();
        let mut children = Vec::new();
        for (key, node) in self.nodes.iter() {
            match data.partial_cmp(&node.data) {
                Some(std::cmp::Ordering::Greater) => {
                    children.push(key);
                }
                Some(std::cmp::Ordering::Less) => {
                    parents.push(key);
                }
                _ => {}
            }
        }
        [parents, children]
    }
    pub fn poset_push(&mut self, data: T, associated_data: D, flip: bool) -> PosetNode
    where
        T: PartialOrd,
    {
        let id = self.nodes.insert_with_key(|key| SlotNode::new(data, key));

        self.associated_data.insert(id, associated_data);
        let new_node = self.nodes.get(id).unwrap();

        let [mut parents, mut children] = self.poset_family(&new_node.data);

        if flip {
            std::mem::swap(&mut parents, &mut children);
        }

        for &parent_key in &parents {
            self.add_edge(parent_key, id);
        }

        for &child_key in &children {
            self.add_edge(id, child_key);
        }

        self.update_transitive_closure(id);

        id
    }

    /// Updates the transitive closure of the poset by propagating the relationships.
    fn update_transitive_closure(&mut self, new_node_key: PosetNode) {
        // Propagate relationships for all descendants of new_node
        let descendants = self.get_all_descendants(new_node_key);
        for &descendant_key in &descendants {
            self.add_edge_if_new(new_node_key, descendant_key);
        }

        // Propagate relationships for all ancestors of new_node
        let ancestors = self.get_all_ancestors(new_node_key);
        for &ancestor_key in &ancestors {
            self.add_edge_if_new(ancestor_key, new_node_key);
        }
    }

    pub fn remove_transitive_edges(mut self) -> HasseDiagram<T, D> {
        let edges_to_remove: Vec<_> = self
            .nodes
            .keys()
            .flat_map(|node_key| self.transitive_edges(node_key))
            .collect();

        // let shift = self.shift();
        for (a, b) in edges_to_remove {
            // println!(
            //     "removing edge from {} to {}",
            //     base62::encode(a.0.as_ffi() - shift),
            //     base62::encode(b.0.as_ffi() - shift)
            // );
            self.remove_edge(a, b);
        }

        let mut hasse = HasseDiagram::new();

        let mut new_map = SecondaryMap::new();

        for (key, node) in self.nodes.into_iter() {
            new_map.insert(key, hasse.add_node(node.data));
            if let Some(d) = self.associated_data.remove(key) {
                hasse.associated_data.insert(new_map[key], d);
            }
            for child in node.children {
                hasse.add_edge(new_map[key], new_map[child]);
            }
        }

        hasse
    }
}

impl<T: PartialOrd, D> FromIterator<(T, D)> for Poset<T, D> {
    fn from_iter<I: IntoIterator<Item = (T, D)>>(iter: I) -> Self {
        let mut poset = Poset::new();
        for (data, assoc) in iter {
            poset.poset_push(data, assoc, false);
        }
        poset
    }
}

/// An iterator over paths in the poset, traversed in BFS order.
pub struct BfsPaths<'a, T, R: Key> {
    queue: VecDeque<Vec<R>>,
    visited: HashSet<Vec<R>>,
    nodes: &'a SlotMap<R, SlotNode<T, R>>,
}

impl<'a, T, R: Key> Iterator for BfsPaths<'a, T, R> {
    type Item = Vec<R>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(path) = self.queue.pop_front() {
            if !self.visited.insert(path.clone()) {
                continue;
            }

            let last_node_key = path.last().unwrap();
            let last_node = self.nodes.get(*last_node_key).unwrap();

            for &child_key in last_node.children.iter() {
                if !path.contains(&child_key) {
                    let mut new_path = path.clone();
                    new_path.push(child_key);
                    self.queue.push_back(new_path);
                }
            }

            return Some(path);
        }
        None
    }
}
#[cfg(test)]
mod tests;
