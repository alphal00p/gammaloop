use std::{collections::BTreeMap, ops::Deref, path::Path};

use crate::{
    model::Model,
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
use log::debug;
use spenso::{
    contraction::Contract,
    network::library::TensorLibraryData,
    structure::{representation::Euclidean, HasStructure, OrderedStructure},
    tensors::{data::StorageTensor, parametric::ParamTensor},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    parse,
    printer::PrintOptions,
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

// impl ToQuoted for SmartString<LazyCompact> {
//     fn to_quoted(&self) -> String {
//         format!("\"{}\"", self)
//     }
// }

impl<A> ToQuoted for A
where
    A: AtomCore,
{
    fn to_quoted(&self) -> String {
        let mut opts = PrintOptions::file();
        opts.hide_namespace = Some("_gammaloop");
        format!("\"{}\"", self.printer(opts))
    }
}

pub trait StripParse<T> {
    fn strip_parse(&self) -> T;
}

impl StripParse<Atom> for String {
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
impl StripParse<Atom> for &String {
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

impl StripParse<bool> for String {
    fn strip_parse(&self) -> bool {
        let a = self
            .as_str()
            .strip_prefix('"')
            .unwrap_or(&self)
            .strip_suffix('"')
            .unwrap_or(&self);
        a.parse().unwrap()
    }
}
impl StripParse<bool> for &String {
    fn strip_parse(&self) -> bool {
        let a = self
            .as_str()
            .strip_prefix('"')
            .unwrap_or(&self)
            .strip_suffix('"')
            .unwrap_or(&self);
        a.parse().unwrap()
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
    pub(crate) fn hedge_order(&self, model: &Model) -> Result<Vec<u8>> {
        let mut hedges = vec![None; self.n_hedges()];

        // Figure out the hedge order for each vertex to align with ufo
        for (_, neighs, v) in self.iter_nodes() {
            let mut particles = Vec::new();

            let name = if let Some(vertex_rule) = &v.vertex_rule {
                for p in &vertex_rule.particles {
                    particles.push(Some(p.clone()));
                }
                Some(vertex_rule.name.clone())
            } else {
                None
            };

            let mut other_order = particles.len();

            for h in neighs {
                let eid = self[&h];

                if self[eid].is_dummy {
                    hedges[h.0] = Some(other_order as u8);
                    other_order += 1;
                    continue;
                }

                let Some(particle) = self[eid].particle.particle() else {
                    hedges[h.0] = Some(other_order as u8);
                    other_order += 1;
                    continue;
                };
                let particle = Some(match self.orientation(h).relative_to(self.flow(h)) {
                    Orientation::Reversed => particle.get_anti_particle(model),
                    _ => particle.clone(),
                });

                if let Some(name) = &name {
                    if let Some((pos, i)) = particles.iter_mut().find_position(|p| &&particle == p)
                    {
                        *i = None;
                        hedges[h.0] = Some(pos as u8);
                    } else {
                        return Err(eyre!(
                            "Particle {} not in vertex rule {}",
                            particle.map(|a| a.name.clone()).unwrap_or("None".into()),
                            name
                        ));
                    }
                } else {
                    hedges[h.0] = Some(other_order as u8);
                    other_order += 1;
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

    pub(crate) fn from_parsed(graph: DotGraph, model: &Model) -> Result<Self> {
        let global_data = graph.global_data.into();
        let graph = graph
            .graph
            .map_data_ref_result(
                |_, _, v| Ok(v),
                ParseEdge::parse(model),
                ParseHedge::parse(),
            )?
            .map_data_ref_result(
                ParseVertex::parse(model, &global_data),
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
    pub(crate) fn dot_serialize(&self) -> String {
        let mut out = String::new();
        self.dot_serialize_fmt(&mut out).unwrap();
        out
    }

    pub(crate) fn dot_serialize_io(
        &self,
        writer: &mut impl std::io::Write,
    ) -> Result<(), std::io::Error> {
        let g = DotGraph::from(self);
        g.write_io(writer)
    }

    pub(crate) fn dot_serialize_fmt(
        &self,
        writer: &mut impl std::fmt::Write,
    ) -> Result<(), std::fmt::Error> {
        let g = DotGraph::from(self);
        g.write_fmt(writer)
    }

    pub(crate) fn from_parsed(graph: ParseGraph, model: &Model) -> Result<Self> {
        let hedge_order = graph.hedge_order(model)?;
        let overall_factor = graph.global_data.overall_factor.clone();
        let add_polarizations = graph.global_data.projectors.is_none();
        let projector = graph.global_data.projectors.clone().unwrap_or(Atom::one());
        let global_prefactor = GlobalPrefactor {
            num: graph.global_data.num.clone(),
            projector,
        };

        let name = graph.global_data.name.clone();
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

                    let prop = e
                        .data
                        .particle
                        .particle()
                        .map(|p| model.get_propagator_for_particle(&p.name).numerator.clone())
                        .unwrap_or(Atom::num(1));

                    color_num_e[eid] = graph[source].color_kronekers(&graph[sink]);

                    let spin_slots = [
                        &graph[source].num_indices.spin_indices.edge_indices,
                        &graph[sink].num_indices.spin_indices.edge_indices,
                    ];

                    let momenta = [graph[&source], graph[&sink]];
                    let spin_nume = UFO.reindex_spin(&spin_slots, &momenta, prop, |i| {
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

            let Some(vertex_rule) = &v.vertex_rule else {
                continue;
            };

            let [mut color_structure, mut couplings, mut spin_structure] =
                vertex_rule.tensors(ni.aind(1), ni.aind(0));

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
                let num = if let Some(num) = v.num {
                    num
                } else {
                    // println!("Num is none  for vertex {}", i);
                    v_spin_num[i.0]
                        .take()
                        .unwrap()
                        .contract(&v_color_num[i.0].take().unwrap())
                        .unwrap()
                        .scalar()
                        .unwrap()
                };

                let dod = if let Some(d) = v.dod { d } else { num.dod() };

                Vertex {
                    name: v.name.unwrap_or(i.to_string()),
                    num,
                    dod,
                    vertex_rule: v.vertex_rule,
                }
            },
            |_, _, _, eid, e| {
                e.map(|e| {
                    let num = e.num.unwrap_or(&color_num_e[eid] * &spin_num_e[eid]);

                    let dod = if let Some(d) = e.dod {
                        d
                    } else {
                        num.dod() - 2
                    };

                    Edge {
                        is_dummy: e.is_dummy,
                        name: e.label.unwrap_or(eid.to_string()),
                        particle: e.particle,
                        num,
                        dod,
                    }
                })
            },
            |_, h| h,
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

            // println!("{}", underlying.dot(&full));
            let covers = tree.covers(&full);
            assert_eq!(
                full,
                covers,
                "Tree does not cover all: {}, lmb specification must be wrong",
                underlying.dot(&covers)
            );
            let external = underlying.internal_crown(&full);
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

        let mut graph = Graph {
            overall_factor,
            global_prefactor,
            name,
            loop_momentum_basis,
            underlying,
        };

        if add_polarizations {
            let pols = graph.generate_polarizations();
            graph.global_prefactor.projector *= pols;
        }
        // loop_momentum_basis.loop_edges.iter_enumerated().sorted_by_key(|(i,v)|);
        Ok(graph)

        // Ok(NumGraph { graph })
    }

    pub(crate) fn from_dot(graph: DotGraph, model: &Model) -> Result<Self> {
        Self::from_parsed(ParseGraph::from_parsed(graph, model)?, model)
    }
    pub(crate) fn from_file<'a, P>(p: P, model: &Model) -> Result<Vec<Self>>
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

    pub(crate) fn from_string<Str: AsRef<str>>(s: Str, model: &Model) -> Result<Vec<Self>> {
        let hedge_graph_set: GraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            linnet::parser::GlobalData,
            NodeStorageVec<DotVertexData>,
        > = GraphSet::from_string(s).unwrap();

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
            let graph = DotGraph { global_data, graph };
            debug!("Parsing: \n{}", graph.debug_dot());
            graphs.push(Graph::from_parsed(
                ParseGraph::from_parsed(graph, model)?,
                model,
            )?);
        }
        Ok(graphs)
    }
}

pub trait IntoGraph<T> {
    fn into_graph(self, model: &Model) -> Result<T>;
}

impl IntoGraph<Vec<Graph>> for String {
    fn into_graph(self, model: &Model) -> Result<Vec<Graph>> {
        let hedge_graph_set: GraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            linnet::parser::GlobalData,
            NodeStorageVec<DotVertexData>,
        > = GraphSet::from_string(self).unwrap();

        Graph::from_hedge_graph_set(hedge_graph_set, model)
    }
}

impl IntoGraph<Graph> for String {
    fn into_graph(self, model: &Model) -> Result<Graph> {
        let graph: DotGraph<NodeStorageVec<DotVertexData>> = DotGraph::from_string(self).unwrap();

        Graph::from_parsed(ParseGraph::from_parsed(graph, model)?, model)
    }
}

impl IntoGraph<Vec<Graph>> for &str {
    fn into_graph(self, model: &Model) -> Result<Vec<Graph>> {
        let hedge_graph_set: GraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            linnet::parser::GlobalData,
            NodeStorageVec<DotVertexData>,
        > = GraphSet::from_string(self).unwrap();

        Graph::from_hedge_graph_set(hedge_graph_set, model)
    }
}

impl IntoGraph<Graph> for &str {
    fn into_graph(self, model: &Model) -> Result<Graph> {
        let graph: DotGraph<NodeStorageVec<DotVertexData>> = DotGraph::from_string(self).unwrap();

        Graph::from_parsed(ParseGraph::from_parsed(graph, model)?, model)
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

        // value.normal_emr_replacement(subgraph, lmb, rep_args, filter_pair)
        for (e, i, v) in value.iter_edges() {
            let loop_expr = value
                .loop_momentum_basis
                .loop_atom(i, GS.loop_mom, &[W_.a___], false);
            let external_expr =
                value
                    .loop_momentum_basis
                    .ext_atom(i, GS.external_mom, &[W_.a___], false);

            graph[i].add_statement("lmb_rep", (loop_expr + external_expr).to_quoted());
        }

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

        stringify!($($code)*).into_graph(&crate::utils::test_utils::load_generic_model($model))
    };

    // Internal rule: End of parsing, with an optional argument.
    // This is matched when the accumulator has collected the code block and we hit a comma.
    (@internal [$($code:tt)*], $model:expr) => {
       stringify!($($code)*).into_graph($model)
    };

    // Internal rule: End of parsing, no optional argument.
    // This is matched when the accumulator has run out of tokens to process.
    (@internal [$($code:tt)*]) => {
        stringify!($($code)*).into_graph( &crate::utils::test_utils::load_generic_model("sm"))
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
    use linnet::half_edge::involution::HedgePair;
    use log::info;
    use spenso::{
        network::{
            library::DummyLibrary, parsing::ShadowedStructure, store::NetworkStore, Network,
        },
        structure::{HasName, PermutedStructure},
        tensors::symbolic::SymbolicTensor,
    };
    use symbolica::atom::{Atom, FunctionBuilder};

    use crate::{
        cli::state::State,
        new_graph::{parse::IntoGraph, LMBext},
        numerator::{aind::Aind, Numerator, UnInit},
    };
    use symbolica::{parse, parse_lit};

    use super::{Graph, StripParse};

    #[test]
    fn test_load() {
        State::default();
        info!("Loading graph");
        let graph: Vec<Graph> = dot!(digraph triangle {
            graph [
                overall_factor = 1;
                multiplicity_factor = 1;
            ]
            edge [
                pdg=1000
            ]
            ext [style=invis]
            ext -> v4
            ext -> v5
            v6 -> ext
            v5 -> v4 [lmb_index=0];
            v6 -> v5;
            v4 -> v6 ;
        }

        digraph triangle2 {
            graph [
                overall_factor = 1;
                multiplicity_factor = 1;
            ]
            edge [
                pdg=1000
            ]
            ext [style=invis]
            ext -> v4
            ext -> v5
            v6 -> ext [id = 2]
            v5 -> v4 [lmb_index=0];
            v6 -> v5;
            v4 -> v6 ;
        }
        ,"scalars")
        .unwrap();

        println!("{}", graph[0].dot_serialize());

        println!(
            "{}",
            graph[0].dot_lmb(&graph[0].underlying.full(), &graph[0].loop_momentum_basis)
        );

        println!("{}", graph[0].loop_momentum_basis);

        println!(
            "{}",
            graph[1].dot_lmb(&graph[1].underlying.full(), &graph[1].loop_momentum_basis)
        );
        println!("{}", graph[1].loop_momentum_basis);
    }

    #[test]
    fn test_loop_momentum_basis() {
        let g: Graph = dot!(
            digraph triangle {
            graph [
            overall_factor = 1;
            multiplicity_factor = 1;
            ]
            edge [
            pdg=1000
            ]
            v4 [num=1]
            ext [style=invis]
            ext -> v4:0 [id=1 is_dummy=true]
            ext -> v5:1 [id=0]
            v6:2 -> ext [id=2]
            v5 -> v4 [lmb_index=0];
            v6 -> v5;
            v4 -> v6 ;
            },"scalars"
        )
        .unwrap();
        // println!("{}", g.loop_momentum_basis);

        insta::assert_ron_snapshot!(g.loop_momentum_basis);

        println!(
            "{}",
            g.dot_lmb(&g.underlying.full(), &g.loop_momentum_basis)
        );
    }

    #[test]
    fn scalar_without_model() {
        let g: Graph = dot!(
            digraph G{
                ext    [style=invis]
                node [num=1]
                ext -> A
                C -> A
                A -> D
                D -> B
                B -> C
                C -> D
                B -> ext
            }
        )
        .unwrap();

        println!("{}", g.dot_serialize())
    }

    #[test]
    fn parse() {
        let g: Graph = dot!(
            digraph G{
                ext    [style=invis]
                ext -> A   [particle=a]
                C -> A    [particle="b"]
                A -> B
                A -> D    [particle="b"]
                D -> B    [particle="b"]
                B -> C    [particle="b"]
                C -> D    [particle="g"]
                B -> ext   [particle=a]
            }
        )
        .unwrap();

        println!("{:#?}", g.loop_momentum_basis);

        println!(
            "{}",
            g.underlying
                .dot_lmb(&g.underlying.full(), &g.loop_momentum_basis)
        );

        println!("{}", g.dot_serialize());
        // return;
        let num =
            Numerator::<UnInit>::default().from_new_graph(&g, &g.underlying.full_filter(), true);

        let expr = num.state.expr.as_view();

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

        println!("{}", num.state.expr);
        // println!("{}", num.state.expr.to_dots());
        // let num = num.gamma_simplify();

        // println!("{}", num.state.expr);
    }
    #[test]
    fn parse_local() {
        let graphs:Vec<Graph>= dot!(
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

        let expr = num.state.expr.as_view();

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

        println!("{}", num.state.expr);
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
        let g: Graph = dot!(
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

        println!(
            "{}",
            g.underlying
                .dot_lmb(&g.underlying.full(), &g.loop_momentum_basis)
        );

        let num =
            Numerator::<UnInit>::default().from_new_graph(&g, &g.underlying.full_filter(), true);

        let expr = num.state.expr.as_view();

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

        println!("{}", num.state.expr);
        let num = num.gamma_simplify();

        println!("{}", num.state.expr);
    }

    #[test]
    fn parse_and_build_edgevec() {
        let graph_1: Graph = dot!(
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
            v6 -> ext
            v5 -> v4 [lmb_index=0];
            v6 -> v5;
            v4 -> v6 ;
        }, "scalars"
        )
        .unwrap();

        let graph_2: Graph = dot! (
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
        }, "scalars"
        )
        .unwrap();

        let graph_1 = &graph_1.underlying;
        let graph_2 = &graph_2.underlying;

        let test_data = vec![true, true, true];

        let graph_1_test = graph_1.new_edgevec(|_, _, p| p.is_paired());

        let graph_2_test = graph_2.new_edgevec(|_, _, p| p.is_paired());

        let mut graph_1_iter = test_data.clone().into_iter();
        let mut graph_2_iter = test_data.clone().into_iter();

        let graph_1_test_2 = graph_1
            .new_edgevec_from_iter(graph_1.iter_edges().map(|(pair, _, _)| {
                if matches!(pair, HedgePair::Paired { .. }) {
                    graph_1_iter.next().unwrap()
                } else {
                    false
                }
            }))
            .unwrap();

        let graph_2_test_2 = graph_2
            .new_edgevec_from_iter(graph_2.iter_edges().map(|(pair, _, _)| {
                if matches!(pair, HedgePair::Paired { .. }) {
                    graph_2_iter.next().unwrap()
                } else {
                    false
                }
            }))
            .unwrap();

        assert_eq!(graph_1_test, graph_1_test_2);
        assert_eq!(graph_2_test, graph_2_test_2);
    }
}
