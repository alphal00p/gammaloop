use std::{
    collections::{BTreeMap, HashMap},
    ops::Deref,
    path::Path,
};

use crate::{
    model::{ArcParticle, ArcPropagator, ArcVertexRule, Model},
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
use dot_parser::canonical::AttrStmt;
use eyre::eyre;
use itertools::Itertools;
use linnet::{
    dot_parser::{DotEdgeData, DotGraph, DotHedgeData, DotVertexData, HedgeGraphSet},
    half_edge::{
        hedgevec::EdgeVec,
        involution::{EdgeData, EdgeIndex, Flow, HedgePair, Orientation},
        nodestore::{BitVecNeighborIter, NodeStorageVec},
        subgraph::{ModifySubgraph, SubGraph},
        tree::SimpleTraversalTree,
        EdgeAccessors, HedgeGraph,
    },
    permutation::Permutation,
};
use spenso::{
    contraction::Contract,
    iterators::IteratableTensor,
    network::library::{LibraryTensor, TensorLibraryData},
    structure::{
        representation::{Euclidean, RepName},
        slot::IsAbstractSlot,
        OrderedStructure, PermutedStructure, ScalarTensor, TensorStructure,
    },
    tensors::{data::StorageTensor, parametric::ParamTensor},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    domains::{atom::AtomField, integer::Integer},
    parse,
};

use super::{hedge_data::NumIndices, Edge, Graph, LMBext, NumHedgeData, Vertex};

#[derive(Debug, Clone)]
pub struct ParseEdge {
    pub label: Option<String>,
    pub particle: ArcParticle,
    pub lmb_id: Option<LoopIndex>,
    pub spin_num: Option<Atom>,
    pub color_num: Option<Atom>,
}

impl ParseEdge {
    pub fn new(particle: ArcParticle) -> Self {
        ParseEdge {
            label: None,
            particle,
            lmb_id: None,
            spin_num: None,
            color_num: None,
        }
    }

    pub fn with_label(mut self, label: String) -> Self {
        self.label = Some(label);
        self
    }

    pub fn with_lmb_id(mut self, lmb_id: LoopIndex) -> Self {
        self.lmb_id = Some(lmb_id);
        self
    }

    pub fn with_spin_num(mut self, spin_num: Atom) -> Self {
        self.spin_num = Some(spin_num);
        self
    }

    pub fn with_color_num(mut self, color_num: Atom) -> Self {
        self.color_num = Some(color_num);
        self
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

impl ParseEdge {
    pub fn localize_ainds(atom: impl AtomCore, eid: EdgeIndex, hedge_pair: HedgePair) -> Atom {
        let a = atom
            .replace(GS.edgeid)
            .with(Atom::num(usize::from(eid) as i64))
            .replace_map(|term, ctx, out| {
                if let AtomView::Fun(f) = term {
                    if f.get_symbol() == GS.edgeid {
                        if f.get_nargs() == 1 {
                            if let Ok(i) = i64::try_from(f.iter().next().unwrap()) {
                                if let Ok(u) = u16::try_from(i) {
                                    *out = eid.aind(u).into();
                                    return true;
                                }
                            }
                        }
                    }
                }
                false
            });

        match hedge_pair {
            HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => a
                .replace(GS.sink_id)
                .with(Atom::num(sink.0 as i64))
                .replace_map(|term, ctx, out| {
                    if let AtomView::Fun(f) = term {
                        if f.get_symbol() == GS.sink_id {
                            if f.get_nargs() == 1 {
                                if let Ok(i) = i64::try_from(f.iter().next().unwrap()) {
                                    if let Ok(u) = u16::try_from(i) {
                                        *out = sink.aind(u).into();
                                        return true;
                                    }
                                }
                            }
                        }
                    }
                    false
                })
                .replace(GS.source_id)
                .with(Atom::num(source.0 as i64))
                .replace_map(|term, ctx, out| {
                    if let AtomView::Fun(f) = term {
                        if f.get_symbol() == GS.source_id {
                            if f.get_nargs() == 1 {
                                if let Ok(i) = i64::try_from(f.iter().next().unwrap()) {
                                    if let Ok(u) = u16::try_from(i) {
                                        *out = source.aind(u).into();
                                        return true;
                                    }
                                }
                            }
                        }
                    }
                    false
                }),
            HedgePair::Unpaired { hedge, flow } => match flow {
                Flow::Source => a
                    .replace(GS.source_id)
                    .with(Atom::num(hedge.0 as i64))
                    .replace_map(|term, ctx, out| {
                        if let AtomView::Fun(f) = term {
                            if f.get_symbol() == GS.source_id {
                                if f.get_nargs() == 1 {
                                    if let Ok(i) = i64::try_from(f.iter().next().unwrap()) {
                                        if let Ok(u) = u16::try_from(i) {
                                            *out = hedge.aind(u).into();
                                            return true;
                                        }
                                    }
                                }
                            }
                        }
                        false
                    }),
                Flow::Sink => a
                    .replace(GS.sink_id)
                    .with(Atom::num(hedge.0 as i64))
                    .replace_map(|term, ctx, out| {
                        if let AtomView::Fun(f) = term {
                            if f.get_symbol() == GS.sink_id {
                                if f.get_nargs() == 1 {
                                    if let Ok(i) = i64::try_from(f.iter().next().unwrap()) {
                                        if let Ok(u) = u16::try_from(i) {
                                            *out = hedge.aind(u).into();
                                            return true;
                                        }
                                    }
                                }
                            }
                        }
                        false
                    }),
            },
        }
    }
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

                let lmb_id: Option<LoopIndex> = e
                    .get::<_, usize>("lmb_id")
                    .transpose()?
                    .map(LoopIndex::from);

                let spin_num = e
                    .get::<_, String>("num")
                    .transpose()?
                    .map(|a| Self::localize_ainds(a.strip_parse(), eid, p));
                let color_num = e
                    .get::<_, String>("color_num")
                    .transpose()?
                    .map(|a| Self::localize_ainds(a.strip_parse(), eid, p));

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
                    lmb_id,
                    particle,
                    spin_num,
                    color_num,
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
    pub spin_num: Option<Atom>,
    pub color_num: Option<Atom>,
}

impl ParseVertex {
    pub fn with_spin_num(mut self, spin_num: Atom) -> Self {
        self.spin_num = Some(spin_num);
        self
    }

    pub fn with_color_num(mut self, color_num: Atom) -> Self {
        self.color_num = Some(color_num);
        self
    }

    pub fn with_label(mut self, label: String) -> Self {
        self.label = Some(label);
        self
    }
}

impl From<ArcVertexRule> for ParseVertex {
    fn from(vertex_rule: ArcVertexRule) -> Self {
        ParseVertex {
            label: None,
            vertex_rule,
            spin_num: None,
            color_num: None,
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

            let spin_num = v.get::<_, String>("num").transpose()?.map(|a| {
                let a = a
                    .as_str()
                    .strip_prefix('"')
                    .unwrap_or(&a)
                    .strip_suffix('"')
                    .unwrap_or(&a);
                parse!(a)
            });
            let color_num = v.get::<_, String>("color_num").transpose()?.map(|a| {
                let a = a
                    .as_str()
                    .strip_prefix('"')
                    .unwrap_or(&a)
                    .strip_suffix('"')
                    .unwrap_or(&a);
                parse!(a)
            });

            if let Some(n) = v.get::<_, String>("int_id") {
                let vertex_rule = model.get_vertex_rule(n.unwrap());

                Ok(ParseVertex {
                    label,
                    vertex_rule,
                    spin_num,
                    color_num,
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
                            spin_num,
                            color_num,
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

#[derive(Clone, Debug)]
pub struct ParseGraph {
    pub global_data: ParseData,
    pub graph: HedgeGraph<ParseEdge, ParseVertex, ParseHedge>,
}

#[derive(Clone, Debug)]
pub struct ParseData {
    pub overall_factor: Atom,
    pub multiplicity_factor: Atom,
    pub color: Atom,
    pub colorless: Atom,
}

impl Default for ParseData {
    fn default() -> Self {
        ParseData {
            overall_factor: Atom::one(),
            multiplicity_factor: Atom::one(),
            color: Atom::one(),
            colorless: Atom::one(),
        }
    }
}

impl ParseData {
    pub fn with_overall_factor(self, overall_factor: Atom) -> Self {
        ParseData {
            overall_factor,
            multiplicity_factor: self.multiplicity_factor,
            color: self.color,
            colorless: self.colorless,
        }
    }

    pub fn with_multiplicity_factor(self, multiplicity_factor: Atom) -> Self {
        ParseData {
            overall_factor: self.overall_factor,
            multiplicity_factor,
            color: self.color,
            colorless: self.colorless,
        }
    }

    pub fn with_color(self, color: Atom) -> Self {
        ParseData {
            overall_factor: self.overall_factor,
            multiplicity_factor: self.multiplicity_factor,
            color,
            colorless: self.colorless,
        }
    }

    pub fn with_colorless(self, colorless: Atom) -> Self {
        ParseData {
            overall_factor: self.overall_factor,
            multiplicity_factor: self.multiplicity_factor,
            color: self.color,
            colorless,
        }
    }
}

impl<'a> TryFrom<&'a Vec<dot_parser::canonical::AttrStmt<(String, String)>>> for ParseData {
    type Error = ();
    fn try_from(
        value: &'a Vec<dot_parser::canonical::AttrStmt<(String, String)>>,
    ) -> std::result::Result<Self, Self::Error> {
        let mut parse_data = ParseData::default();
        for v in value {
            if let AttrStmt::Graph((a, b)) = v {
                // println!("Graph Attr:{a}={b}");
                match a.as_str() {
                    "overall_factor" => {
                        parse_data = parse_data.with_overall_factor(b.strip_parse());
                    }
                    "multiplicity_factor" => {
                        parse_data = parse_data.with_multiplicity_factor(b.strip_parse());
                    }
                    "color" => {
                        parse_data = parse_data.with_color(b.strip_parse());
                    }
                    "colorless" => {
                        parse_data = parse_data.with_colorless(b.strip_parse());
                    }
                    _ => {}
                }
            }
        }

        Ok(parse_data)
    }
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
        global_data: ParseData,
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

                let dod = num_spin.iter_flat().map(|(i, v)| v.dod()).max().unwrap();

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

                    Edge {
                        name: e.label.unwrap_or(eid.to_string()),
                        propagator: ArcPropagator(
                            model.get_propagator_for_particle(&e.particle.name),
                        ),
                        particle: e.particle,
                        color_num,
                        dod: spin_num.dod() - 2,
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

            let full = underlying.full_filter();
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
            vertex_slots: Vec::new().into(),
            external_connections: None,
        })

        // Ok(NumGraph { graph })
    }

    pub fn from_dot(
        graph: DotGraph,
        auto_detect_vertex_rule: bool,
        model: &Model,
        global_data: ParseData,
    ) -> Result<Self> {
        Self::from_parsed(
            ParseGraph::from_parsed(graph, auto_detect_vertex_rule, model, global_data)?,
            model,
        )
    }
    pub fn from_file<'a, P>(p: P, model: &Model) -> Result<Vec<Self>>
    where
        P: AsRef<Path>,
    {
        let hedge_graph_set: HedgeGraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            ParseData,
            NodeStorageVec<DotVertexData>,
        > = HedgeGraphSet::from_file(p).unwrap();

        Self::from_hedge_graph_set(hedge_graph_set, model)
    }

    pub fn from_string<Str: AsRef<str>>(s: Str, model: &Model) -> Result<Vec<Self>> {
        let hedge_graph_set: HedgeGraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            ParseData,
            NodeStorageVec<DotVertexData>,
        > = HedgeGraphSet::from_string(s).unwrap();

        Self::from_hedge_graph_set(hedge_graph_set, model)
    }

    fn from_hedge_graph_set(
        set: HedgeGraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            ParseData,
            NodeStorageVec<DotVertexData>,
        >,
        model: &Model,
    ) -> Result<Vec<Self>> {
        let mut graphs = Vec::new();

        for (graph, data) in set.set.into_iter().zip(set.global_data.into_iter()) {
            graphs.push(Graph::from_parsed(
                ParseGraph::from_parsed(graph, true, model, data).unwrap(),
                model,
            )?);
        }
        Ok(graphs)
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
    use linnet::{dot_parser::DotGraph, half_edge::HedgeGraph};
    use spenso::{
        network::{
            library::DummyLibrary, parsing::ShadowedStructure, store::NetworkStore, Network,
        },
        structure::{concrete_index::FlatIndex, HasName, HasStructure, PermutedStructure},
        tensors::{data::GetTensorData, symbolic::SymbolicTensor},
    };
    use symbolica::{
        atom::{Atom, FunctionBuilder},
        parse_lit,
    };

    use crate::{
        new_graph::LMBext,
        numerator::{aind::Aind, GlobalPrefactor, Numerator, UnInit},
        tests_from_pytest::load_generic_model,
    };

    use super::Graph;

    #[test]
    fn parse() {
        let graphs = dot!(
            digraph G{
                e1      [flow=sink]
                e4      [flow=source]
                e1 -> A   [particle=a]
                C -> A    [particle="b"]
                A -> D    [particle="b"]
                D -> B    [particle="b"]
                B -> C    [particle="b"]
                C -> D    [particle="g"]
                B -> e4   [particle=a]
            }
        )
        .unwrap();

        let g = &graphs[0];

        println!(
            "{}",
            g.underlying
                .dot_lmb(&g.underlying.full(), &g.loop_momentum_basis)
        );

        let num = Numerator::<UnInit>::default().from_new_graph(&g, &g.underlying.full_filter());

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
    #[test]
    fn parse_local() {
        let graphs= dot!(
            digraph G{
                e1      [flow=sink]
                e2      [flow=sink]
                e3      [flow=source]
                e4      [flow=source]
                A [num="1"]
                B [num="1"]
                C [num="1"]
                D [num="1"]
                e1 -> A  [pdg=1000]
                e2 -> C  [pdg=1000]
                C -> A    [pdg=1000, num="P(eid,spenso::mink(4,eid(0)))*Q(eid,spenso::mink(4,eid(0)))"]
                A -> D    [pdg=1000, num="P(eid,spenso::mink(4,0))"]
                D -> B    [pdg=1000, num="1"]
                B -> e3   [pdg=1000]
                C -> D    [pdg=1000, num="1"]
                B -> e4   [pdg=1000]
            }

            digraph G2{
                e1      [flow=sink]
                e2      [flow=sink]
                e3      [flow=source]
                e4      [flow=source]
                A [num="1"]
                B [num="1"]
                C [num="1"]
                D [num="1"]
                e1 -> A  [pdg=1000]
                e2 -> C  [pdg=1000]
                C -> A    [pdg=1000, num="P(eid,spenso::mink(4,eid(0)))*Q(eid,spenso::mink(4,eid(0)))"]
                A -> D    [pdg=1000, num="K(eid,spenso::mink(4,0))"]
                D -> B    [pdg=1000, num="1"]
                B -> e3   [pdg=1000]
                C -> D    [pdg=1000, num="1"]
                B -> e4   [pdg=1000]
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

        let num = Numerator::<UnInit>::default().from_new_graph(&g, &g.underlying.full_filter());

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

        println!(
            "{}",
            num.state
                .colorless
                .get_owned_linear(FlatIndex::from(0))
                .unwrap()
                .to_dots()
        );
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
                e1      [flow=sink]
                e2      [flow=sink]
                e3      [flow=source]
                e4      [flow=source]
                A [num="1" color_num="a"]
                B [num="1"]
                C [num="1"]
                D [num="1"]
                E [num="1"]
                e1 -> A  [pdg=1000]
                e2 -> C  [pdg=1000]
                C -> A    [pdg=1000, num="P(eid,spenso::mink(4,0))"]
                A -> D    [pdg=1000, num="P(eid,spenso::mink(4,0))"]
                D -> B    [pdg=1000, num="1"]
                B -> e3   [pdg=1000]
                E -> B    [pdg=1000,lmb_id=0]
                E -> A    [pdg=1000,lmb_id=1]
                E -> C   [pdg=1000]
                C -> D    [pdg=1000, num="1"]
                B -> e4   [pdg=1000]
            },"scalars"
        )
        .unwrap();

        let g = &graphs[0];
        println!(
            "{}",
            g.underlying
                .dot_lmb(&g.underlying.full(), &g.loop_momentum_basis)
        );

        let num = Numerator::<UnInit>::default().from_new_graph(&g, &g.underlying.full_filter());

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
