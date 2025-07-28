use std::{collections::BTreeMap, ops::Deref, path::Path};

use crate::{
    model::{ArcPropagator, Model},
    momentum_sample::LoopIndex,
    numerator::{
        aind::{Aind, NewAind},
        ufo::UFO,
        GlobalPrefactor,
    },
    symbolica_ext::CallSymbol,
    utils::{GS, W_},
};
use bitvec::vec::BitVec;
use color_eyre::Result;
use eyre::eyre;
use itertools::Itertools;
use linnet::{
    half_edge::{
        involution::{EdgeVec, HedgePair, Orientation},
        nodestore::NodeStorageVec,
        subgraph::{ModifySubgraph, SubGraph},
        tree::SimpleTraversalTree,
        EdgeAccessors, HedgeGraph,
    },
    parser::{DotEdgeData, DotGraph, DotHedgeData, DotVertexData, GraphSet},
    permutation::Permutation,
};
use spenso::{
    contraction::Contract,
    iterators::IteratableTensor,
    network::library::LibraryTensor,
    structure::{
        representation::{Euclidean, RepName},
        OrderedStructure, PermutedStructure,
    },
    tensors::{data::StorageTensor, parametric::ParamTensor},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    parse,
};

use super::{
    edge::ParseEdge,
    global::ParseData,
    hedge_data::{NumIndices, ParseHedge},
    vertex::ParseVertex,
    Edge, Graph, LMBext, NumHedgeData, Vertex,
};

pub trait ToQuoted {
    fn to_quoted(&self) -> String;
}

impl<A> ToQuoted for A
where
    A: AtomCore,
{
    fn to_quoted(&self) -> String {
        format!("\"{}\"", self.to_canonical_string())
    }
}

pub trait StripParse {
    fn strip_parse(&self) -> Atom;
}

impl StripParse for String {
    fn strip_parse(&self) -> Atom {
        let a = self
            .as_str()
            .strip_prefix('"')
            .unwrap_or(&self)
            .strip_suffix('"')
            .unwrap_or(&self);
        parse!(a)
    }
}
impl StripParse for &String {
    fn strip_parse(&self) -> Atom {
        let a = self
            .as_str()
            .strip_prefix('"')
            .unwrap_or(&self)
            .strip_suffix('"')
            .unwrap_or(&self);
        parse!(a)
    }
}

#[derive(Clone, Debug)]
pub struct ParseGraph {
    pub global_data: ParseData,
    pub graph: HedgeGraph<ParseEdge, ParseVertex, ParseHedge>,
}

impl Deref for ParseGraph {
    type Target = HedgeGraph<ParseEdge, ParseVertex, ParseHedge>;

    fn deref(&self) -> &Self::Target {
        &self.graph
    }
}

impl ParseGraph {
    pub fn hedge_order(&self, model: &Model) -> Result<Vec<u8>> {
        let mut hedges = vec![None; self.n_hedges()];

        // Figure out the hedge order for each vertex to align with ufo
        for (_, neighs, v) in self.iter_nodes() {
            let mut particles = Vec::new();

            for p in &v.vertex_rule.particles {
                particles.push(Some(p.clone()));
            }

            for h in neighs {
                let eid = self[&h];

                let particle = Some(match self.orientation(h).relative_to(-self.flow(h)) {
                    Orientation::Reversed => self[eid].particle.get_anti_particle(model),
                    _ => self[eid].particle.clone(),
                });
                if let Some((pos, i)) = particles.iter_mut().find_position(|p| &&particle == p) {
                    *i = None;
                    hedges[h.0] = Some(pos as u8);
                } else {
                    return Err(eyre!("Particle not in vertex rule:{:?}", particle));
                }
            }

            if particles.iter().any(|p| p.is_some()) {
                return Err(eyre!(
                    "Particles to vertex rules no match for set: {:?}",
                    particles
                ));
            }
        }

        let hedges: Option<Vec<u8>> = hedges.into_iter().collect();
        if let Some(hedges) = hedges {
            Ok(hedges)
        } else {
            Err(eyre!("Nodes do not cover hedges"))
        }
    }

    pub fn from_parsed(
        graph: DotGraph,
        auto_detect_vertex_rule: bool,
        model: &Model,
    ) -> Result<Self> {
        let global_data = graph.global_data.into();
        let graph = graph
            .graph
            .map_data_ref_result(
                |_, _, v| Ok(v),
                ParseEdge::parse(model),
                ParseHedge::parse(),
            )?
            .map_data_ref_result(
                ParseVertex::parse(model, auto_detect_vertex_rule),
                |_, _, _, e| Ok(e.map(Clone::clone)),
                |(_, h)| Ok(h.clone()),
            )?;

        Ok(Self { graph, global_data })
    }
}

pub trait DOD {
    fn dod(&self) -> i32;
}

impl DOD for Atom {
    fn dod(&self) -> i32 {
        self.as_view().dod()
    }
}
impl DOD for AtomView<'_> {
    fn dod(&self) -> i32 {
        let rescaled = self
            .replace(GS.emr_mom.f(&[W_.a__]))
            .with(GS.emr_mom.f(&[W_.a__]) * GS.rescale);

        let series = rescaled
            .series(GS.rescale, Atom::Zero, (4, 1).into(), true)
            .unwrap();

        let dod = series.degree();

        if dod.is_integer() {
            dod.numerator().to_i64().unwrap() as i32
        } else {
            panic!("{dod} for {self}")
        }
    }
}

impl Graph {
    pub fn dot_serialize(&self) -> String {
        let mut out = String::new();
        self.dot_serialize_fmt(&mut out).unwrap();
        out
    }

    pub fn dot_serialize_io(&self, writer: &mut impl std::io::Write) -> Result<(), std::io::Error> {
        let g = DotGraph::from(self);
        g.write_io(writer)
    }

    pub fn dot_serialize_fmt(
        &self,
        writer: &mut impl std::fmt::Write,
    ) -> Result<(), std::fmt::Error> {
        let g = DotGraph::from(self);
        g.write_fmt(writer)
    }

    pub fn from_parsed(graph: ParseGraph, model: &Model) -> Result<Self> {
        let hedge_order = graph.hedge_order(model)?;
        let multiplicity = graph.global_data.multiplicity_factor.clone();
        let overall_factor = &graph.global_data.overall_factor;
        let global_prefactor = GlobalPrefactor {
            color: graph.global_data.color.clone(),
            colorless: &graph.global_data.colorless * overall_factor,
        };
        let graph = graph.graph.map_data_ref(
            |_, _, v| v.clone(),
            |_, _, _, e| e.map(|e| e.clone()),
            NumIndices::parse(&graph),
        );

        let mut full_cut: BitVec = graph.full_filter();

        let graph = graph.map(
            |_, _, v| v,
            |_, _, _, _, e| e,
            |h, num_indices| NumHedgeData {
                num_indices,
                node_order: hedge_order[h.0],
            },
        );

        let mut color_num_e: EdgeVec<_> = vec![Atom::num(1); graph.n_edges()].into();
        let mut spin_num_e: EdgeVec<_> = vec![Atom::i(); graph.n_edges()].into();

        let mut lmb_ids = BTreeMap::new();

        // Generate edge numerators
        for (p, eid, e) in graph.iter_edges() {
            match p {
                HedgePair::Paired { source, sink } => {
                    if let Some(lmb_id) = e.data.lmb_id {
                        if let Some(old_value) = lmb_ids.insert(lmb_id, eid) {
                            return Err(eyre!(
                                "lmb_id {lmb_id:?} already exists with value {old_value:?}",
                            ));
                        }
                        full_cut.sub(p);
                    }

                    let prop =
                        ArcPropagator(model.get_propagator_for_particle(&e.data.particle.name));

                    color_num_e[eid] = graph[source].color_kronekers(&graph[sink]);

                    let spin_slots = [
                        &graph[source].num_indices.spin_indices.edge_indices,
                        &graph[sink].num_indices.spin_indices.edge_indices,
                    ];

                    let momenta = [graph[&source], graph[&sink]];
                    let spin_nume =
                        UFO.reindex_spin(&spin_slots, &momenta, prop.numerator.clone(), |i| {
                            Aind::Edge(usize::from(eid) as u16, i as u16)
                        });

                    // println!("{}", spin_nume);
                    spin_num_e[eid] = spin_nume;
                }
                _ => {}
            }
        }

        let mut v_color_num: Vec<Option<ParamTensor<OrderedStructure<Euclidean, Aind>>>> =
            vec![None; graph.n_nodes()];
        let mut v_spin_num = v_color_num.clone();

        // Generate vertex numerators
        for (ni, c, v) in graph.iter_nodes() {
            let mut color_slots = vec![];
            let mut spin_slots = vec![];
            let mut order = vec![];
            let mut momenta = vec![];
            for h in c {
                color_slots.push(&graph[h].num_indices.color_indices.vertex_indices);
                spin_slots.push(&graph[h].num_indices.spin_indices.vertex_indices);
                order.push(graph[h].node_order);
                momenta.push(graph[&h]);
            }

            let perm = Permutation::sort(&order);
            perm.apply_slice_in_place(&mut color_slots);
            perm.apply_slice_in_place(&mut spin_slots);

            let [mut color_structure, mut couplings, mut spin_structure] =
                v.vertex_rule.tensors(ni.aind(1), ni.aind(0));

            spin_structure.map_data_mut(|a| {
                *a = UFO.reindex_spin(&spin_slots, &momenta, a.clone(), |u| ni.aind(u as u16));
            });

            couplings.map_data_mut(|a| *a = UFO.normalize_complex(a.clone()));

            color_structure.map_data_mut(|a| {
                *a = UFO.reindex_color(&color_slots, a.clone(), |u| ni.aind(u as u16));
            });

            v_spin_num[ni.0] = Some(spin_structure.contract(&couplings).unwrap());

            v_color_num[ni.0] = Some(color_structure);
        }

        let underlying = graph.map(
            |_, i, v| {
                let summed_slot = Euclidean {}.new_slot(1, i.aind(1));
                let num_spin = if let Some(s) = v.spin_num {
                    ParamTensor::from_dense(
                        PermutedStructure::from_iter([summed_slot]).structure,
                        vec![s],
                    )
                    .unwrap()
                } else {
                    v_spin_num[i.0].take().unwrap()
                };

                let num_color = if let Some(c) = v.color_num {
                    ParamTensor::from_dense(
                        PermutedStructure::from_iter([summed_slot]).structure,
                        vec![c],
                    )
                    .unwrap()
                } else {
                    v_color_num[i.0].take().unwrap()
                };

                let dod = if let Some(d) = v.dod {
                    d
                } else {
                    num_spin.iter_flat().map(|(i, v)| v.dod()).max().unwrap()
                };

                Vertex {
                    label: v.label.unwrap_or(i.to_string()),
                    num_color,
                    dod,
                    num_spin,
                    vertex_rule: v.vertex_rule,
                }
            },
            |_, _, _, eid, e| {
                e.map(|e| {
                    let color_num = if let Some(c) = e.color_num {
                        c
                    } else {
                        color_num_e[eid].clone()
                    };

                    let spin_num = if let Some(s) = e.spin_num {
                        s
                    } else {
                        spin_num_e[eid].clone()
                    };

                    let dod = if let Some(d) = e.dod {
                        d
                    } else {
                        spin_num.dod() - 2
                    };

                    Edge {
                        is_dummy: e.is_dummy,
                        name: e.label.unwrap_or(eid.to_string()),
                        propagator: ArcPropagator(
                            model.get_propagator_for_particle(&e.particle.name),
                        ),
                        particle: e.particle,
                        color_num,
                        dod,
                        spin_num,
                    }
                })
            },
            |e, h| h,
        );

        let mut loop_momentum_basis = if let Some(i) = full_cut.included_iter().next() {
            let tree = SimpleTraversalTree::depth_first_traverse(
                &underlying,
                &full_cut,
                &underlying.node_id(i),
                None,
            )
            .unwrap();

            let mut full = underlying.full_filter();

            for (p, _, i) in underlying.iter_edges() {
                if i.data.is_dummy {
                    full.sub(p);
                }
            }
            let covers = tree.covers(&full);
            assert_eq!(
                full,
                covers,
                "Tree does not cover all: {}, lmb specification must be wrong",
                underlying.dot(&covers)
            );
            let external = underlying.full_crown(&full);
            underlying.lmb_impl(&full, &tree.tree_subgraph, external)
        } else {
            panic!("ata")
        };

        // println!("{:#?}", loop_momentum_basis);

        let inv_lmb_ids: BTreeMap<_, _> = lmb_ids.iter().map(|(k, v)| (*v, *k)).collect();

        let mut i = 0;

        while i < loop_momentum_basis.loop_edges.len() {
            if let Some(&target_pos) =
                inv_lmb_ids.get(&loop_momentum_basis.loop_edges[LoopIndex(i)])
            {
                if target_pos.0 < loop_momentum_basis.loop_edges.len() && target_pos.0 != i {
                    loop_momentum_basis.swap_loops(LoopIndex(i), target_pos);
                    continue;
                }
            }
            i += 1;
        }
        // loop_momentum_basis.loop_edges.iter_enumerated().sorted_by_key(|(i,v)|);
        Ok(Graph {
            multiplicity,
            global_prefactor,
            name: "".to_string(),
            loop_momentum_basis,
            underlying,
        })

        // Ok(NumGraph { graph })
    }

    pub fn from_dot(graph: DotGraph, auto_detect_vertex_rule: bool, model: &Model) -> Result<Self> {
        Self::from_parsed(
            ParseGraph::from_parsed(graph, auto_detect_vertex_rule, model)?,
            model,
        )
    }
    pub fn from_file<'a, P>(p: P, model: &Model) -> Result<Vec<Self>>
    where
        P: AsRef<Path>,
    {
        let hedge_graph_set: GraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            linnet::parser::GlobalData,
            NodeStorageVec<DotVertexData>,
        > = GraphSet::from_file(p).unwrap();

        Self::from_hedge_graph_set(hedge_graph_set, model)
    }

    pub fn from_string<Str: AsRef<str>>(s: Str, model: &Model) -> Result<Vec<Self>> {
        let hedge_graph_set: GraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            linnet::parser::GlobalData,
            NodeStorageVec<DotVertexData>,
        > = GraphSet::from_string(s).unwrap();

        // println!("HOO");
        Self::from_hedge_graph_set(hedge_graph_set, model)
    }

    fn from_hedge_graph_set(
        set: GraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            linnet::parser::GlobalData,
            NodeStorageVec<DotVertexData>,
        >,
        model: &Model,
    ) -> Result<Vec<Self>> {
        let mut graphs = Vec::new();

        for (graph, global_data) in set.set.into_iter().zip(set.global_data.into_iter()) {
            graphs.push(Graph::from_parsed(
                ParseGraph::from_parsed(DotGraph { global_data, graph }, true, model).unwrap(),
                model,
            )?);
        }
        Ok(graphs)
    }
}

impl From<&Graph> for DotGraph {
    fn from(value: &Graph) -> Self {
        let global_data = value.global_data();
        let mut graph: HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData> =
            value.underlying.map_data_ref(
                |_, _, v| v.into(),
                |_, _, _, e| e.map(|e| e.into()),
                |_, d| d.into(),
            );

        for (l, i) in value.loop_momentum_basis.loop_edges.iter_enumerated() {
            graph[i].0.add_statement("lmb_id", l);
        }

        DotGraph { global_data, graph }
    }
}

#[macro_export]
macro_rules! dot {
    // ------------------ Internal Rules (Do not call directly) ------------------

    (@internal [$($code:tt)*], $model:literal) => {

        Graph::from_string(stringify!($($code)*), &load_generic_model($model))
    };

    // Internal rule: End of parsing, with an optional argument.
    // This is matched when the accumulator has collected the code block and we hit a comma.
    (@internal [$($code:tt)*], $model:expr) => {
        Graph::from_string(stringify!($($code)*), $model)
    };

    // Internal rule: End of parsing, no optional argument.
    // This is matched when the accumulator has run out of tokens to process.
    (@internal [$($code:tt)*]) => {
        Graph::from_string(stringify!($($code)*), &load_generic_model("sm"))
    };

    // Internal rule: The "accumulator".
    // It takes the next token ($next), adds it to the $code accumulator,
    // and recursively calls the macro with the rest of the tokens ($($rest)*).
    (@internal [$($code:tt)*] $next:tt $($rest:tt)*) => {
        dot!(@internal [$($code)* $next] $($rest)*)
    };

    // ------------------ Public Entry Point ------------------

    // This is the only rule users should call. It kicks off the process by
    // invoking the internal rules with an empty accumulator `[]`.
    ($($all_tokens:tt)+) => {
        dot!(@internal [] $($all_tokens)+)
    };
}

#[cfg(test)]
pub mod test {
    use idenso::metric::MetricSimplifier;
    use spenso::{
        network::{
            library::DummyLibrary, parsing::ShadowedStructure, store::NetworkStore, Network,
        },
        structure::{concrete_index::FlatIndex, HasName, HasStructure, PermutedStructure},
        tensors::{data::GetTensorData, symbolic::SymbolicTensor},
    };
    use symbolica::atom::{Atom, FunctionBuilder};

    use crate::{
        new_graph::LMBext,
        numerator::{aind::Aind, Numerator, UnInit},
        tests_from_pytest::load_generic_model,
    };

    use super::Graph;

    #[test]
    fn test_loop_momentum_basis() {
        let graphs = dot!(
            digraph triangle {
            graph [
            overall_factor = 1;
            multiplicity_factor = 1;
            ]
            edge [
            pdg=1000
            ]
            ext [style=invis]
            ext -> v4 [id=0]
            ext -> v5 [id=1]
            v6 -> ext [id=2]
            v5 -> v4 [lmb_index=0];
            v6 -> v5;
            v4 -> v6 ;
            },"scalars"
        )
        .unwrap();

        let g = &graphs[0];

        println!(
            "{}",
            g.dot_lmb(&g.underlying.full(), &g.loop_momentum_basis)
        );
    }

    #[test]
    fn parse() {
        let graphs = dot!(
            digraph G{
                ext    [style=invis]
                ext -> A   [particle=a]
                C -> A    [particle="b"]
                A -> D    [particle="b"]
                D -> B    [particle="b"]
                B -> C    [particle="b"]
                C -> D    [particle="g"]
                B -> ext   [particle=a]
            }
        )
        .unwrap();

        let g = &graphs[0];

        println!("{:#?}", g.loop_momentum_basis);

        println!(
            "{}",
            g.underlying
                .dot_lmb(&g.underlying.full(), &g.loop_momentum_basis)
        );

        println!("{}", g.dot_serialize());

        let num =
            Numerator::<UnInit>::default().from_new_graph(&g, &g.underlying.full_filter(), true);

        let expr = num.state.colorless.get_ref_linear(0.into()).unwrap();

        let lib: DummyLibrary<SymbolicTensor<Aind>> = DummyLibrary::<_>::new();
        let net =
            Network::<NetworkStore<SymbolicTensor<Aind>, Atom>, _, Aind>::try_from_view(expr, &lib)
                .unwrap();

        println!("{}", expr);
        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |_| None,
                |a| {
                    if let Ok(a) = PermutedStructure::<ShadowedStructure<Aind>>::try_from(
                        a.expression.as_view(),
                    ) {
                        a.structure
                            .name()
                            .map(|s| {
                                if let Some(a) = a.structure.args() {
                                    FunctionBuilder::new(s).add_args(&a).finish().to_string()
                                } else {
                                    s.to_string()
                                }
                            })
                            .unwrap_or("".to_string())
                    } else {
                        "".to_string()
                    }
                }
            )
        );

        let num = num.color_simplify();

        println!("{}", num.state.color);
        println!("{}", num.state.colorless);
        let num = num.gamma_simplify();

        println!("{}", num.state.colorless);
    }
    #[test]
    fn parse_local() {
        let graphs= dot!(
            digraph G{
                edge [num="1" dod=-1000 pdg=1000]
                node [num="1" dod=-1000]
                ext     [style=invis]
                A [id=1]
                B [id=0]
                C
                D [id=2]
                ext -> A  [id=1]
                ext -> C
                C -> A    [num="P(eid,spenso::mink(4,eid(0)))*Q(eid,spenso::mink(4,eid(0)))"]
                A -> D    [num="P(eid,spenso::mink(4,0))"]
                D:1 -> B:2    [id=2]
                B -> ext
                C -> D
                B -> ext
            }

            digraph G2{
                ext [style=invis]
                A [num="1"]
                B [num="1"]
                C [num="1"]
                D [num="1"]
                ext -> A  [pdg=1000]
                ext -> C  [pdg=1000]
                C -> A    [pdg=1000, num="P(eid,spenso::mink(4,eid(0)))*Q(eid,spenso::mink(4,eid(0)))"]
                A -> D    [pdg=1000, num="K(eid,spenso::mink(4,0))"]
                D -> B    [pdg=1000, num="1"]
                B -> ext   [pdg=1000]
                C -> D    [pdg=1000, num="1"]
                B -> ext   [pdg=1000]
            }
            ,
            "scalars"
        ).unwrap();

        let g = &graphs[1];

        println!(
            "{}",
            g.underlying
                .dot_lmb(&g.underlying.full(), &g.loop_momentum_basis)
        );

        let num =
            Numerator::<UnInit>::default().from_new_graph(&g, &g.underlying.full_filter(), true);

        let expr = num.state.colorless.get_ref_linear(0.into()).unwrap();

        let lib: DummyLibrary<SymbolicTensor<Aind>> = DummyLibrary::<_>::new();
        let net =
            Network::<NetworkStore<SymbolicTensor<Aind>, Atom>, _, Aind>::try_from_view(expr, &lib)
                .unwrap();

        println!("{}", expr);
        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |_| None,
                |a| {
                    if let Ok(a) = PermutedStructure::<ShadowedStructure<Aind>>::try_from(
                        a.expression.as_view(),
                    ) {
                        a.structure
                            .name()
                            .map(|s| {
                                if let Some(a) = a.structure.args() {
                                    FunctionBuilder::new(s).add_args(&a).finish().to_string()
                                } else {
                                    s.to_string()
                                }
                            })
                            .unwrap_or("".to_string())
                    } else {
                        "".to_string()
                    }
                }
            )
        );

        let num = num.color_simplify();

        println!("Hi{}", num.state.color);
        println!("{}", num.state.colorless);
        // let num = num.gamma_simplify();

        // println!(
        //     "{}",
        //     num.state
        //         .colorless
        //         .get_owned_linear(FlatIndex::from(0))
        //         .unwrap()
        //         .to_dots()
        // );
    }

    #[test]
    fn parse_lmbsetting() {
        let graphs = dot!(
            digraph G{
                graph [
                    multiplicity_factor = "1/2"
                    overall_factor = "2*x"
                    color = "spenso::t"
                    colorless = "spenso::gamma"
                ]
                ext [style=invis]
                A [num="1" color_num="a"]
                B [num="1"]
                C [num="1"]
                D [num="1"]
                E [num="1"]
                ext -> A  [pdg=1000]
                ext -> C  [pdg=1000]
                C -> A    [pdg=1000, num="P(eid,spenso::mink(4,0))"]
                A -> D    [pdg=1000, num="P(eid,spenso::mink(4,0))"]
                D -> B    [pdg=1000, num="1"]
                B -> ext   [pdg=1000]
                E -> B    [pdg=1000,lmb_id=0]
                E -> A    [pdg=1000,lmb_id=1]
                E -> C   [pdg=1000]
                C -> D    [pdg=1000, num="1"]
                B -> ext   [pdg=1000]
            },"scalars"
        )
        .unwrap();

        let g = &graphs[0];
        println!(
            "{}",
            g.underlying
                .dot_lmb(&g.underlying.full(), &g.loop_momentum_basis)
        );

        let num =
            Numerator::<UnInit>::default().from_new_graph(&g, &g.underlying.full_filter(), true);

        let expr = num.state.colorless.get_ref_linear(0.into()).unwrap();

        let lib: DummyLibrary<SymbolicTensor<Aind>> = DummyLibrary::<_>::new();
        let net =
            Network::<NetworkStore<SymbolicTensor<Aind>, Atom>, _, Aind>::try_from_view(expr, &lib)
                .unwrap();

        println!("{}", expr);
        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |_| None,
                |a| {
                    if let Ok(a) = PermutedStructure::<ShadowedStructure<Aind>>::try_from(
                        a.expression.as_view(),
                    ) {
                        a.structure
                            .name()
                            .map(|s| {
                                if let Some(a) = a.structure.args() {
                                    FunctionBuilder::new(s).add_args(&a).finish().to_string()
                                } else {
                                    s.to_string()
                                }
                            })
                            .unwrap_or("".to_string())
                    } else {
                        "".to_string()
                    }
                }
            )
        );

        let num = num.color_simplify();

        println!("{}", num.state.color);
        let num = num.gamma_simplify();

        println!("{}", num.state.colorless);
    }
}
