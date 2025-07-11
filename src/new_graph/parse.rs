use std::{ops::Deref, sync::Arc};

use crate::{
    graph::{InteractionVertexInfo, VertexInfo},
    model::{ArcParticle, ArcPropagator, ArcVertexRule, Model, VertexRule},
    numerator::{aind::Aind, ufo::UFO},
    symbolica_ext::CallSymbol,
    utils::{GS, W_},
};
use color_eyre::Result;
use eyre::eyre;
use itertools::Itertools;
use linnet::{
    dot_parser::{DotEdgeData, DotGraph, DotHedgeData, DotVertexData},
    half_edge::{
        hedgevec::EdgeVec,
        involution::{EdgeData, EdgeIndex, HedgePair, Orientation},
        nodestore::{BitVecNeighborIter, NodeStorage, NodeStorageVec},
        EdgeAccessors, HedgeGraph,
    },
    permutation::Permutation,
};
use spenso::{
    contraction::Contract,
    iterators::IteratableTensor,
    structure::{
        representation::{Euclidean, RepName},
        slot::IsAbstractSlot,
        OrderedStructure, TensorStructure,
    },
    tensors::{data::StorageTensor, parametric::ParamTensor},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    domains::{atom::AtomField, integer::Z, rational::Q},
    parse,
    poly::evaluate::Variable,
};

use super::{hedge_data::NumIndices, Edge, Graph, LMBext, NumHedgeData, Vertex};

#[derive(Debug, Clone)]
pub struct ParseEdge {
    pub label: Option<String>,
    pub particle: ArcParticle,
    pub num: Option<Atom>,
}

impl ParseEdge {
    pub fn parse<'a>(
        model: &'a Model,
    ) -> impl FnMut(
        &'a DotGraph,
        EdgeIndex,
        HedgePair,
        EdgeData<&'a DotEdgeData>,
    ) -> Result<EdgeData<Self>> {
        |graph: &'a DotGraph, eid: EdgeIndex, p: HedgePair, e_data: EdgeData<&'a DotEdgeData>| {
            e_data.map_result(|e| {
                let label = e.get::<_, String>("label").transpose()?;

                let num = e.get::<_, String>("num").transpose()?.map(|a| parse!(a));

                let particle = if let Some(v) = e.get::<_, isize>("pdg") {
                    model.get_particle_from_pdg(v?)
                } else if let Some(v) = e.get::<_, String>("particle") {
                    let pname = v?;
                    let pname = pname
                        .as_str()
                        .strip_prefix('"')
                        .unwrap_or(&pname)
                        .strip_suffix('"')
                        .unwrap_or(&pname);
                    model.get_particle(pname)
                } else {
                    return Err(eyre!("no pdg or name found for edge"));
                };

                Ok(ParseEdge {
                    particle,
                    num,
                    label,
                })
            })
        }
    }
}

#[derive(Debug, Clone)]
pub struct ParseVertex {
    pub label: Option<String>,
    pub vertex_rule: ArcVertexRule,
    pub num: Option<Atom>,
}

impl From<ArcVertexRule> for ParseVertex {
    fn from(vertex_rule: ArcVertexRule) -> Self {
        ParseVertex {
            label: None,
            vertex_rule,
            num: None,
        }
    }
}

impl ParseVertex {
    pub fn parse<'a>(
        model: &'a Model,
        auto_detect_vertex_rule: bool,
    ) -> impl FnMut(
        &'a HedgeGraph<ParseEdge, &'a DotVertexData, ParseHedge>,
        BitVecNeighborIter<'a>,
        &'a &'a DotVertexData,
    ) -> Result<Self> {
        move |g, n, v| {
            let label = v.id().map(|id| id.to_string());

            let num = v.get::<_, String>("num").transpose()?.map(|a| parse!(a));

            if let Some(n) = v.get::<_, String>("int_id") {
                let vertex_rule = model.get_vertex_rule(n.unwrap());

                Ok(ParseVertex {
                    label,
                    vertex_rule,
                    num,
                })
            } else if auto_detect_vertex_rule {
                let mut particles: Vec<ArcParticle> = n
                    .map(|h| {
                        let eid = g[&h];
                        let particle = match g.orientation(h).relative_to(g.flow(h)) {
                            Orientation::Reversed => g[eid].particle.get_anti_particle(model),
                            _ => g[eid].particle.clone(),
                        };
                        particle
                    })
                    .collect();
                particles.sort();

                let res = model.particle_set_to_vertex_rules_map.get(&particles);
                if let Some(res) = res {
                    if res.len() == 1 {
                        Ok(ParseVertex {
                            label,
                            vertex_rule: res[0].clone(),
                            num,
                        })
                    } else {
                        Err(eyre!("Multiple vertex rules for {:?}", particles))
                    }
                } else {
                    let particles = particles.iter().map(|p| p.name.as_str()).collect_vec();
                    Err(eyre!(
                        "Failed to find vertex rule for particles: {:?}",
                        particles
                    ))
                }
            } else {
                Err(eyre!("Vertex rule not supplied"))
            }
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct ParseHedge {
    hedge_id: Option<usize>,
}

impl ParseHedge {
    pub fn parse<'a>() -> impl FnMut(&'a DotHedgeData) -> Result<Self> {
        |h| {
            let hedge_id = h
                .statement
                .as_ref()
                .map(|s| s.parse::<usize>().ok())
                .flatten();
            Ok(ParseHedge { hedge_id })
        }
    }
}

pub struct ParseGraph {
    pub graph: HedgeGraph<ParseEdge, ParseVertex, ParseHedge>,
}

impl Deref for ParseGraph {
    type Target = HedgeGraph<ParseEdge, ParseVertex, ParseHedge>;

    fn deref(&self) -> &Self::Target {
        &self.graph
    }
}

impl ParseGraph {
    pub fn from_parsed(
        graph: DotGraph,
        auto_detect_vertex_rule: bool,
        model: &Model,
    ) -> Result<Self> {
        let graph = graph
            .map_data_ref_result(
                |_, _, v| Ok(v),
                ParseEdge::parse(model),
                ParseHedge::parse(),
            )?
            .map_data_ref_result(
                ParseVertex::parse(model, auto_detect_vertex_rule),
                |_, _, _, e| Ok(e.map(Clone::clone)),
                |h| Ok(h.clone()),
            )?;

        Ok(Self { graph })
    }
}

pub trait DOD {
    fn dod(&self) -> i32;
}

impl DOD for Atom {
    fn dod(&self) -> i32 {
        let rescaled = self
            .replace(GS.emr_mom.f(&[W_.a__]))
            .with(GS.emr_mom.f(&[W_.a__]) * GS.rescale);

        let mut dod = 0;
        let series = rescaled
            .series(GS.rescale, Atom::Zero, (4, 1).into(), true)
            .unwrap();
        // println!("trailing: {}", series.absolute_order());
        let field = AtomField {
            statistical_zero_test: false,
            cancel_check_on_division: true,
            custom_normalization: None,
        };
        // let r = rescaled
        //     .series(GS.rescale, Atom::Zero, (4, 1).into(), true)
        //     .unwrap()
        //     .to_atom()
        //     .to_rational_polynomial::<_, _, u8>(&Q, &Z, Some(Arc::new(vec![GS.rescale.into()])));
        // dod = r.numerator.degree(0);
        // println!("{r}");
        self.as_view().dod()
    }
}
impl DOD for AtomView<'_> {
    fn dod(&self) -> i32 {
        // println!("{}has dod:", self);
        let rescaled = self
            .replace(GS.emr_mom.f(&[W_.a__]))
            .with(GS.emr_mom.f(&[W_.a__]) * GS.rescale);

        let mut dod: i32 = 0;
        let series = rescaled
            .series(GS.rescale, Atom::Zero, (4, 1).into(), true)
            .unwrap();

        for (a, b) in series
            .to_atom()
            .coefficient_list::<u8>(&[Atom::var(GS.rescale)])
        {
            // println!("a{a}\nb{b}");
            if let AtomView::Pow(a) = a.as_view() {
                if let Ok(i) = i64::try_from(a.get_exp()) {
                    if dod < i as i32 {
                        dod = i as i32
                    }
                }
            } else if let AtomView::Var(a) = a.as_view() {
                if a.get_symbol() == GS.rescale {
                    if dod < 1 {
                        dod = 1;
                    }
                }
            }
        }
        // println!("{dod}");
        dod
    }
}

impl Graph {
    pub fn from_parsed(graph: ParseGraph, model: &Model) -> Result<Self> {
        let graph = graph.map_data_ref(
            |_, _, v| v.clone(),
            |_, _, _, e| e.map(|e| e.clone()),
            NumIndices::parse(&graph),
        );
        let mut hedges = vec![None; graph.n_hedges()];

        // Figure out the hedge order for each vertex to align with ufo
        for (nid, neighs, v) in graph.iter_nodes() {
            let mut particles = Vec::new();

            for p in &v.vertex_rule.particles {
                particles.push(Some(p.clone()));
            }

            for h in neighs {
                let eid = graph[&h];
                let particle = Some(match graph.orientation(h).relative_to(graph.flow(h)) {
                    Orientation::Reversed => graph[eid].particle.get_anti_particle(model),
                    _ => graph[eid].particle.clone(),
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
            let graph = graph.map(
                |_, _, v| v,
                |_, _, _, _, e| e,
                |h, num_indices| NumHedgeData {
                    num_indices,
                    node_order: hedges[h.0],
                },
            );

            let mut color_num: EdgeVec<_> = vec![Atom::num(1); graph.n_edges()].into();
            let mut spin_num: EdgeVec<_> = vec![Atom::i(); graph.n_edges()].into();

            for (p, eid, e) in graph.iter_edges() {
                match p {
                    HedgePair::Paired { source, sink } => {
                        let prop =
                            ArcPropagator(model.get_propagator_for_particle(&e.data.particle.name));

                        for (i, j) in graph[source]
                            .num_indices
                            .color_indices
                            .edge_indices
                            .external_structure_iter()
                            .zip_eq(
                                graph[sink]
                                    .num_indices
                                    .color_indices
                                    .edge_indices
                                    .external_structure_iter(),
                            )
                        {
                            if i.rep().matches(&j.rep()) {
                                color_num[eid] *= j.rep().id(i.aind, j.aind);
                            } else {
                                panic!("Should be the same rep found:{} and {}", i, j)
                            }
                        }

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
                        spin_num[eid] = spin_nume;
                    }
                    _ => {}
                }
            }

            let mut v_color_num: Vec<Option<ParamTensor<OrderedStructure<Euclidean, Aind>>>> =
                vec![None; graph.n_nodes()];
            let mut v_spin_num = v_color_num.clone();

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

                let [mut color_structure, mut couplings, mut spin_structure] = v
                    .vertex_rule
                    .tensors(Aind::Vertex(ni.0 as u16, 1), Aind::Vertex(ni.0 as u16, 0));

                spin_structure.map_data_mut(|a| {
                    *a = UFO.reindex_spin(&spin_slots, &momenta, a.clone(), |u| {
                        Aind::Vertex(ni.0 as u16, u as u16)
                    });
                });

                couplings.map_data_mut(|a| *a = UFO.normalize_complex(a.clone()));

                color_structure.map_data_mut(|a| {
                    *a = UFO.reindex_color(&color_slots, a.clone(), |u| {
                        Aind::Vertex(ni.0 as u16, u as u16)
                    });
                });

                v_spin_num[ni.0] = Some(spin_structure.contract(&couplings).unwrap());

                v_color_num[ni.0] = Some(color_structure);
            }

            let underlying = graph.map(
                |_, i, v| {
                    let num_spin = v_spin_num[i.0].take().unwrap();

                    let dod = num_spin.iter_flat().map(|(i, v)| v.dod()).max().unwrap();

                    Vertex {
                        label: v.label.unwrap_or(i.to_string()),
                        num_color: v_color_num[i.0].take().unwrap(),
                        dod,
                        num_spin,
                        vertex_rule: v.vertex_rule,
                    }
                },
                |_, _, _, eid, e| {
                    e.map(|e| Edge {
                        name: e.label.unwrap_or(eid.to_string()),
                        propagator: ArcPropagator(
                            model.get_propagator_for_particle(&e.particle.name),
                        ),
                        particle: e.particle,
                        color_num: color_num[eid].clone(),
                        spin_num: spin_num[eid].clone(),
                        dod: spin_num[eid].dod() - 2,
                    })
                },
                |e, h| h,
            );

            Ok(Graph {
                multiplicity: Atom::num(1),
                name: "".to_string(),
                loop_momentum_basis: underlying.lmb(&underlying.full_filter()),
                underlying,
                vertex_slots: Vec::new().into(),
                external_connections: None,
            })

            // Ok(NumGraph { graph })
        } else {
            Err(eyre!("Not all hedges are set"))
        }
    }

    pub fn from_dot(graph: DotGraph, auto_detect_vertex_rule: bool, model: &Model) -> Result<Self> {
        Self::from_parsed(
            ParseGraph::from_parsed(graph, auto_detect_vertex_rule, model)?,
            model,
        )
    }
}

#[cfg(test)]
pub mod test {
    use linnet::{dot, dot_parser::DotGraph, half_edge::HedgeGraph};
    use spenso::{
        network::{
            library::DummyLibrary, parsing::ShadowedStructure, store::NetworkStore, Network,
        },
        structure::{HasName, PermutedStructure},
        tensors::{data::GetTensorData, symbolic::SymbolicTensor},
    };
    use symbolica::atom::{Atom, FunctionBuilder};

    use crate::{
        new_graph::LMBext,
        numerator::{aind::Aind, GlobalPrefactor, Numerator, UnInit},
        tests_from_pytest::load_generic_model,
    };

    use super::Graph;

    #[test]
    fn parse() {
        let model = load_generic_model("sm");
        let graph: Result<DotGraph, _> = dot!(
            digraph G{
                e1      [flow=sink]
                e2      [flow=sink]
                e3      [flow=source]
                e4      [flow=source]
                e1 -> A  [particle=a]
                e2 -> "C"  [particle="b"]
                "C" -> A    [particle="b"]
                A -> D    [particle="b"]
                D -> B    [particle="b"]
                B -> e3   [particle="b"]
                "C" -> D    [particle="g"]
                B -> e4   [particle=a]
            }
        );

        match graph {
            Ok(graph) => {
                println!("{}", graph.dot_display(&graph.full()));
                let g = Graph::from_dot(graph, true, &model);
                match g {
                    Ok(g) => {
                        println!(
                            "{}",
                            g.underlying
                                .dot_lmb(&g.underlying.full(), &g.loop_momentum_basis)
                        );

                        let num = Numerator::<UnInit>::default().from_new_graph(
                            &g,
                            &g.underlying.full_filter(),
                            &GlobalPrefactor::default(),
                        );

                        let expr = num.get_single_atom().unwrap();

                        let lib: DummyLibrary<SymbolicTensor<Aind>> = DummyLibrary::<_>::new();
                        let mut net =
                            Network::<NetworkStore<SymbolicTensor<Aind>, Atom>, _, Aind>::try_from_view(
                                expr.as_view(), &lib,
                            )
                            .unwrap();

                        println!("{}", expr);
                        println!(
                            "{}",
                            net.dot_display_impl(
                                |a| a.to_string(),
                                |_| None,
                                |a| {
                                    if let Ok(a) =
                                        PermutedStructure::<ShadowedStructure<Aind>>::try_from(
                                            a.expression.as_view(),
                                        )
                                    {
                                        a.structure
                                            .name()
                                            .map(|s| {
                                                if let Some(a) = a.structure.args() {
                                                    FunctionBuilder::new(s)
                                                        .add_args(&a)
                                                        .finish()
                                                        .to_string()
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
                    }
                    Err(e) => {
                        println!("Error parsing graph: {}", e);
                    }
                }
            }
            Err(e) => {
                println!("Error parsing graph: {:?}", e);
            }
        }
    }
}
