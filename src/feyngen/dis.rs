use std::collections::BTreeMap;

use ahash::AHashMap;
use spenso::{
    shadowing::ETS,
    structure::{
        representation::{Rep, RepName},
        slot::IsAbstractSlot,
    },
    symbolica_utils::SerializableAtom,
};
use symbolica::{
    atom::{Atom, AtomView, FunctionBuilder, Symbol},
    coefficient::Coefficient,
    fun,
    id::{MatchSettings, Pattern, PatternOrMap, Replacement},
    state::{FunctionAttribute, State},
    symb,
};

use crate::{
    graph::{
        half_edge::{
            subgraph::{Cycle, Inclusion, OrientedCut, SubGraph, SubGraphOps},
            EdgeId, Hedge, HedgeGraph, Orientation,
        },
        BareGraph, Edge, Vertex,
    },
    momentum::{Sign, SignOrZero, Signature},
    numerator::GlobalPrefactor,
};

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Embeddings {
    pub cuts: BTreeMap<Embedding, Vec<OrientedCut>>,
    pub basis: Vec<Cycle>,
}

#[derive(Debug, Clone, Eq, PartialEq, Hash, PartialOrd, Ord)]
pub struct Embedding {
    pub windings: Vec<i32>,
}

pub struct IFCuts {
    pub cuts: BTreeMap<Embedding, [Vec<OrientedCut>; 2]>,
    pub basis: Vec<Cycle>,
}

impl IFCuts {
    pub fn remove_empty(&mut self) {
        self.cuts.retain(|_, v| !v[0].is_empty() & !v[1].is_empty());
    }
}

impl Embedding {
    pub fn windings_in_field(&self, n: u32) -> Vec<u32> {
        self.windings
            .iter()
            .map(|&i| {
                if i < 0 {
                    (n as i32 + i) as u32
                } else {
                    i as u32
                }
            })
            .collect()
    }
}

impl Embeddings {
    pub fn remove_singles(&mut self) {
        self.cuts.retain(|_, v| v.len() > 1);
    }

    pub fn if_split<E, V>(self, graph: &HedgeGraph<E, V>, filter: &impl Fn(&E) -> bool) -> IFCuts {
        let cuts = self
            .cuts
            .into_iter()
            .map(|(k, v)| {
                let mut split = [Vec::new(), Vec::new()];
                for cut in v {
                    let mut is_in = false;
                    for (_, e) in cut.iter_edges(graph) {
                        if filter(e.as_ref().data.unwrap()) {
                            is_in = true;
                        }
                    }

                    if is_in {
                        split[0].push(cut);
                    } else {
                        split[1].push(cut);
                    }
                }

                (k, split)
            })
            .collect();

        IFCuts {
            cuts,
            basis: self.basis,
        }
    }

    pub fn classify(
        iter: impl IntoIterator<Item = OrientedCut>,
        basis: Vec<Cycle>,
        filter: impl Fn(&OrientedCut) -> bool,
        flip_sym: bool,
    ) -> Embeddings {
        let mut cuts = BTreeMap::new();

        for cut in iter {
            if !filter(&cut) {
                continue;
            }
            let mut windings = Vec::new();

            let mut first_non_zero = None;

            for cycle in basis.iter() {
                let mut winding_number = cut.winding_number(cycle);
                if flip_sym {
                    if let Some(sign) = first_non_zero {
                        winding_number *= sign as i32;
                    } else if winding_number > 0 {
                        first_non_zero = Some(Sign::Positive);
                    } else if winding_number < 0 {
                        first_non_zero = Some(Sign::Negative);
                        winding_number *= -1;
                    };
                }
                windings.push(winding_number);
            }
            cuts.entry(Embedding { windings })
                .or_insert_with(Vec::new)
                .push(cut);
        }

        Embeddings { cuts, basis }
    }

    // pub fn push<E, V>(&mut self, cut: SignedCut, graph: &HedgeGraph<E, V>) {
    //     let mut found = false;
    //     for (i, c) in self.cuts.iter_mut().enumerate() {
    //         for other in c.iter() {
    //             if cut.are_same_embedding(other, graph) {
    //                 found = true;
    //                 break;
    //             }
    //         }

    //         if found {
    //             c.insert(cut);
    //             return;
    //         }
    //     }

    //     let mut new_hash = HashSet::new();
    //     new_hash.insert(cut);
    //     self.cuts.push(new_hash);
    // }

    pub fn n_embeddings(&self) -> usize {
        self.cuts.len()
    }

    // pub fn from_graph<E, V>(graph: &HedgeGraph<E, V>) -> Self {
    //     let mut embeddings = Embeddings { cuts: Vec::new() };
    //     let mut ext = graph.involution.inv.iter().filter(|(_, e)| e.is_identity());

    //     let s = ext.next().unwrap().0.clone();
    //     let t = ext.next().unwrap().0.clone();

    //     let mut cuts = AHashSet::new();

    //     graph.all_s_t_cuts(&s, &t, &mut cuts);

    //     println!("n cuts {}", cuts.len());

    //     let mut n_signed = 0;

    //     for c in cuts {
    //         let node = if c.is_empty() {
    //             s.clone()
    //         } else {
    //             graph.nesting_node_from_subgraph(c)
    //         };
    //         let all = SignedCut::all_from_internal(node, graph);
    //         n_signed += all.len();

    //         for c in all {
    //             embeddings.push(c, graph);
    //         }
    //     }

    //     println!("n_signed:{n_signed}");
    //     embeddings
    // }
}

// impl SignedCut {
//     // pub fn from_vertex_set
//     //
//     //
//     pub fn signature(&self) -> Vec<Sign> {
//         let mut signature = Vec::new();
//         for h in self.cut_content.iter_ones() {
//             if self.reference[h] {
//                 // if h is outgoing
//                 signature.push(Sign::Positive);
//             } else {
//                 signature.push(Sign::Negative);
//             }
//         }
//         signature
//     }

//     pub fn bare_signature<E, V>(&self, graph: &HedgeGraph<E, V>) -> Signature {
//         self.cut_content
//             .iter()
//             .enumerate()
//             .map(|(i, h)| {
//                 if *h {
//                     graph.involution.inv[i].to_sign()
//                 } else {
//                     SignOrZero::Zero
//                 }
//             })
//             .collect()
//     }

//     pub fn all_initial_state_cuts<E, V>(graph: &HedgeGraph<E, V>) -> Vec<SignedCut> {
//         let mut all_cuts = Vec::new();

//         for c in graph.non_cut_edges() {
//             if c.count_ones() == 0 {
//                 continue;
//             }
//             let mut all_sources = graph.empty_filter();

//             for h in c.iter_ones() {
//                 match graph.involution.inv[h] {
//                     InvolutiveMapping::Identity(_) => {
//                         panic!("cut edge is identity")
//                     }
//                     InvolutiveMapping::Source(_) => {
//                         all_sources.set(h, true);
//                     }
//                     InvolutiveMapping::Sink(_) => {}
//                 }
//             }

//             let n_cut_edges: u8 = all_sources.count_ones().try_into().unwrap();

//             let pset = unsafe { PowersetIterator::new(n_cut_edges.unchecked_sub(1)) };

//             for i in pset {
//                 let mut cut_content = graph.empty_filter();
//                 for (j, h) in all_sources.iter_ones().enumerate() {
//                     if let Some(j) = j.checked_sub(1) {
//                         if i[j] {
//                             cut_content.set(graph.involution.inv(h), true);
//                         } else {
//                             cut_content.set(h, true);
//                         }
//                     } else {
//                         cut_content.set(h, true);
//                     }
//                 }
//                 all_cuts.push(SignedCut {
//                     cut_content,
//                     reference: all_sources.clone(),
//                 });
//             }
//         }

//         all_cuts
//     }

//     pub fn all_from_internal<E, V>(node: HedgeNode, graph: &HedgeGraph<E, V>) -> Vec<SignedCut> {
//         let mut all_sources = graph.empty_filter();
//         let mut cuts = Vec::new();
//         let mut all_outgoing = graph.empty_filter();

//         let init_cut = node.hairs;

//         for h in init_cut.iter_ones() {
//             match graph.involution.inv[h] {
//                 InvolutiveMapping::Identity(_) => {}
//                 InvolutiveMapping::Sink(i) => {
//                     all_outgoing.set(h, true);
//                     all_sources.set(i, true);
//                 }
//                 _ => {
//                     all_sources.set(h, true);
//                     all_outgoing.set(h, true);
//                 }
//             }
//         }

//         let n_cut_edges: u8 = all_sources.count_ones().try_into().unwrap();

//         println!("n_cut_edges: {n_cut_edges}");

//         let mut pset = PowersetIterator::new(n_cut_edges);

//         pset.next(); // skip empty set

//         for i in pset {
//             let mut cut_content = graph.empty_filter();
//             for (j, h) in all_sources.iter_ones().enumerate() {
//                 if i[j] {
//                     cut_content.set(graph.involution.inv(h), true);
//                 } else {
//                     cut_content.set(h, true);
//                 }
//             }
//             cuts.push(SignedCut {
//                 cut_content,
//                 reference: all_outgoing.clone(),
//             });
//         }

//         cuts.pop(); // remove the last one, which is the full sink cut

//         cuts
//     }

//     pub fn from_internal<E, V>(
//         internal: InternalSubGraph,
//         signs: &[Sign],
//         graph: &HedgeGraph<E, V>,
//     ) -> Option<Self> {
//         let mut cut_content = graph.empty_filter();

//         let init_cut = graph.nesting_node_from_subgraph(internal).hairs;

//         for (i, h) in init_cut.iter_ones().enumerate() {
//             match graph.involution.inv[h] {
//                 InvolutiveMapping::Identity(_) => {}
//                 InvolutiveMapping::Sink(connectedh) => {
//                     if let Sign::Positive = *signs.get(i)? {
//                         if cut_content[connectedh] {
//                             continue;
//                         }
//                         cut_content.set(h, true);
//                     } else {
//                         if cut_content[h] {
//                             continue;
//                         }
//                         cut_content.set(connectedh, true);
//                     }
//                 }
//                 _ => {}
//             }
//         }

//         Some(SignedCut {
//             cut_content,
//             reference: init_cut,
//         })
//     }

//     pub fn sign(&self, i: usize) -> SignOrZero {
//         match (self.cut_content[i], self.reference[i]) {
//             (true, true) => SignOrZero::Plus,
//             (false, false) => SignOrZero::Zero,
//             (true, false) => SignOrZero::Minus,
//             (false, true) => SignOrZero::Zero,
//         }
//     }

//     pub fn winding_number(&self, cycle: &Cycle) -> i32 {
//         let mut winding_number = 0;

//         for h in cycle.filter.iter_ones() {
//             winding_number += self.sign(h) * 1;
//         }

//         winding_number
//     }

//     pub fn are_same_embedding<E, V>(&self, other: &Self, graph: &HedgeGraph<E, V>) -> bool {
//         let sym_diff = self.cut_content.clone() ^ &other.cut_content;

//         let mut outside_sign = None;
//         for h in sym_diff.iter_ones() {
//             match graph.involution.inv[h] {
//                 InvolutiveMapping::Sink(_) => {
//                     if let Some(Sign::Positive) = outside_sign {
//                         return false;
//                     } else if outside_sign.is_none() {
//                         outside_sign = Some(Sign::Negative)
//                     }
//                 }
//                 InvolutiveMapping::Source(_) => {
//                     if let Some(Sign::Negative) = outside_sign {
//                         return false;
//                     } else if outside_sign.is_none() {
//                         outside_sign = Some(Sign::Positive)
//                     }
//                 }
//                 InvolutiveMapping::Identity(_) => {
//                     panic!("external in cut")
//                 }
//             }
//         }
//         true
//     }
// }

impl BareGraph {
    // pub fn all_cuts_contain(&self, cut_contents: &HashSet<(usize,Arc<Particle>>)>) -> bool {
    //     let h = HedgeGraph::from_bare_to_ext_half(d);

    //     // println!("{}", h.base_dot());
    //     let mut ext = h.involution.inv.iter().filter(|(d, e)| e.is_identity());

    //     let s = ext.next().unwrap().0.clone();
    //     let t = ext.next().unwrap().0.clone();

    //     let mut cuts = AHashSet::new();
    //     h.all_s_t_cuts(&s, &t, &mut cuts);

    //     for c in cuts{
    //         let cut = h.nesting_node_from_subgraph(c);

    //         let

    //     }

    //     true
    // }
}

impl BareGraph {
    pub fn dis_graph(&self) -> DisGraph {
        let mut h = self.hedge_representation.clone();

        let mut elec_node = None;

        if let Some((elec, d)) = h.iter_egdes(&h.full_filter()).find(|(e, n)| {
            if self.edges[**n.as_ref().data.unwrap()]
                .particle
                .pdg_code
                .abs()
                == 11
            {
                true
            } else {
                false
            }
        }) {
            if let EdgeId::Paired { source, sink } = elec {
                elec_node = Some(h.node_id(source).clone());
            }
        }

        let mut included_hedge = None;
        let mut hedge_filter = h.empty_filter();
        let node = if let Some(s) = elec_node {
            for i in s.hairs.included_iter() {
                if self.edges[*h.get_edge_data(i)].particle.pdg_code.abs() == 11 {
                    included_hedge = Some(i);
                    break;
                }
            }
            s
        } else {
            h.node_id(Hedge(0)).clone()
        };

        let (basis, tree) = h
            .paton_cycle_basis(&h.full_graph(), &node, included_hedge)
            .unwrap(); //TODO start basis on electron edge
        h.align_to_tree_underlying(&tree);

        println!("{}", h.base_dot());
        println!("{}", h.dot(&tree.tree));

        let mut seen_pdg22 = None;
        let mut seen_pdg11 = None;
        let lmbsymb = symb!("k");
        let graph = h.map(
            |bare_vertex_id| DisVertex {
                bare_vertex_id,
                bare_vertex: self.vertices[bare_vertex_id].clone(),
            },
            |e, d| {
                let mut mom_e = Atom::new_num(0);

                let mut first_cycle = None;
                let mut only_cycle = true;

                for (i, c) in basis.iter().enumerate() {
                    if let EdgeId::Paired { source, .. } = e {
                        if c.filter.includes(&source) {
                            if first_cycle.is_none() {
                                first_cycle = Some(i);
                            } else {
                                only_cycle = false;
                            }
                            mom_e = mom_e + fun!(lmbsymb, i as i32)
                        }
                    }
                }
                d.and_then(|bare_edge_id| {
                    let bare_edge = self.edges[bare_edge_id].clone();

                    let marked = if only_cycle {
                        if let Some(i) = first_cycle {
                            match bare_edge.particle.pdg_code.abs() {
                                11 => {
                                    if seen_pdg11.is_some() {
                                        false
                                    } else {
                                        seen_pdg11 = Some((e, i));
                                        true
                                    }
                                }
                                22 => {
                                    if seen_pdg22.is_some() {
                                        false
                                    } else {
                                        seen_pdg22 = Some((e, i));
                                        true
                                    }
                                }
                                _ => false,
                            }
                        } else {
                            false
                        }
                    } else {
                        false
                    };

                    Some(DisEdge {
                        bare_edge,
                        bare_edge_id,
                        marked,
                        momentum: mom_e,
                    })
                })
            },
        );

        let mut outer_graph = graph.empty_filter();

        for (i, e) in graph.iter_egdes(&graph.full_filter()) {
            match i {
                EdgeId::Paired { source, sink } => {
                    if e.data.as_ref().unwrap().bare_edge.particle.pdg_code.abs() == 11 {
                        outer_graph.set(source.0, true);
                        for i in graph.node_id(sink).included_iter() {
                            outer_graph.set(i.0, true);
                        }
                        outer_graph.set(sink.0, true);
                    }
                }
                _ => {}
            }
        }

        let inner_graph = outer_graph.complement(&graph);

        let mink = Rep::new_self_dual("mink").unwrap();
        let mu = mink.new_slot(4, 3).to_atom();
        let nu = mink.new_slot(4, 2).to_atom();
        let metric = fun!(ETS.metric, mu, nu);
        let p = symb!("p");
        let phat2 = Atom::new_var(symb!("phat")).pow(&Atom::new_num(2));
        let pp = fun!(p, mu) * fun!(p, nu);
        let diminv = Atom::parse("1/(2-D)").unwrap();

        let w1_proj = GlobalPrefactor {
            color: Atom::new_num(1),
            colorless: &diminv * (&metric - &pp / &phat2),
        };

        let w2_proj = GlobalPrefactor {
            color: Atom::new_num(1),
            colorless: (diminv * (metric - &pp / &phat2) + &pp / &phat2) / &phat2,
        };

        let mut w1 = crate::numerator::Numerator::default()
            .from_dis_graph(self, &graph, &inner_graph, Some(&w1_proj))
            .color_simplify()
            .gamma_simplify()
            .get_single_atom()
            .unwrap();

        let mut w2 = crate::numerator::Numerator::default()
            .from_dis_graph(self, &graph, &inner_graph, Some(&w2_proj))
            .color_simplify()
            .gamma_simplify()
            .get_single_atom()
            .unwrap();

        numerator_dis_apply(&mut w1);
        numerator_dis_apply(&mut w2);
        let mut denominator = Atom::new_num(1);
        let denomsymb = symb!("prop");
        let emr_mom_symb = symb!("Q");
        for (j, e) in graph.iter_egdes(&inner_graph) {
            let edge = &e.data.as_ref().unwrap().bare_edge;
            let i = e.data.as_ref().unwrap().bare_edge_id;
            if matches!(j, EdgeId::Paired { .. }) {
                let mass = edge.particle.mass.expression.clone().unwrap_or(Atom::Zero);
                let emr_mom = fun!(emr_mom_symb, i as i32);
                denominator = denominator * fun!(denomsymb, mass, emr_mom);
            };
        }

        DisGraph {
            graph,
            numerator: vec![w1.0.expand(), w2.0.expand()],
            denominator,
            lmb_photon: seen_pdg22.unwrap(),
            marked_electron_edge: seen_pdg11.unwrap(),
            basis,
        }
    }
}

pub fn numerator_dis_apply(num: &mut SerializableAtom) {
    let f_ = symb!("f_");
    let g_ = symb!("g_");
    let a_ = Atom::new_var(symb!("a_"));
    let b_ = Atom::new_var(symb!("b_"));
    let c_ = Atom::new_var(symb!("c_"));

    let dim = symb!("D");
    let p = symb!("p");
    let q = symb!("q");
    let emrmom = symb!("Q");
    let dot = State::get_symbol_with_attributes(
        "dot",
        &[FunctionAttribute::Symmetric, FunctionAttribute::Linear],
    )
    .unwrap();

    let reps = vec![
        (
            fun!(ETS.metric, a_, b_).pow(&Atom::new_num(2)),
            Atom::new_var(dim),
        ),
        (
            fun!(emrmom, b_, a_) * fun!(emrmom, c_, a_),
            fun!(dot, fun!(emrmom, b_), fun!(emrmom, c_)),
        ),
        (
            fun!(p, a_) * fun!(emrmom, c_, a_),
            fun!(dot, p, fun!(emrmom, c_)),
        ),
        (
            fun!(q, a_) * fun!(emrmom, c_, a_),
            fun!(dot, q, fun!(emrmom, c_)),
        ),
        (fun!(p, a_) * fun!(p, a_), fun!(dot, p, p)),
        (fun!(p, a_) * fun!(q, a_), fun!(dot, p, q)),
        (fun!(q, a_) * fun!(q, a_), fun!(dot, q, q)),
        (
            fun!(ETS.metric, a_, b_) * fun!(emrmom, c_, a_),
            fun!(emrmom, c_, b_),
        ),
        (fun!(ETS.metric, a_, b_) * fun!(p, a_), fun!(p, b_)),
        (fun!(ETS.metric, a_, b_) * fun!(q, a_), fun!(q, b_)),
    ];

    let replacements: Vec<(Pattern, PatternOrMap)> = reps
        .into_iter()
        .map(|(a, b)| (a.into_pattern(), b.into_pattern().into()))
        .collect();

    let settings = MatchSettings {
        rhs_cache_size: 0,
        ..Default::default()
    };

    let reps: Vec<Replacement> = replacements
        .iter()
        .map(|(lhs, rhs)| Replacement::new(lhs, rhs).with_settings(&settings))
        .collect();

    num.replace_repeat_multiple(&reps)
}
impl HedgeGraph<(usize, crate::graph::Edge, bool), (usize, crate::graph::Vertex)> {
    pub fn dis_graph(mut self) -> HedgeGraph<DisEdge, DisVertex> {
        let (basis, tree) = self.cycle_basis();
        self.align_to_tree_underlying(&tree);

        let lmbsymb = symb!("k");
        self.map(
            |(bare_vertex_id, bare_vertex)| DisVertex {
                bare_vertex_id,
                bare_vertex,
            },
            |e, d| {
                let mut mom_e = Atom::new_num(0);
                for (i, c) in basis.iter().enumerate() {
                    if let EdgeId::Paired { source, .. } = e {
                        if c.filter.includes(&source) {
                            mom_e = mom_e + fun!(lmbsymb, i as i32)
                        }
                    }
                }
                d.and_then(|(bare_edge_id, bare_edge, marked)| {
                    Some(DisEdge {
                        bare_edge,
                        bare_edge_id,
                        marked,
                        momentum: mom_e,
                    })
                })
            },
        )
    }

    // pub fn propagators()
}
//
pub struct DisGraph {
    graph: HedgeGraph<DisEdge, DisVertex>,
    marked_electron_edge: (EdgeId, usize),
    lmb_photon: (EdgeId, usize),
    numerator: Vec<Atom>,
    denominator: Atom,
    basis: Vec<Cycle>,
}

impl DisGraph {
    pub fn numerator(&self, cut: &OrientedCut, total: Symbol) -> Vec<Atom> {
        let emr_to_lmb_cut = self.emr_to_lmb_and_cut(cut, total);

        let reps: Vec<_> = emr_to_lmb_cut
            .iter()
            .map(|(k, v)| Replacement::new(k, v))
            .collect();

        self.numerator
            .iter()
            .map(|a| a.replace_all_multiple(&reps))
            .collect()
    }

    pub fn denominator(&self, cut: &OrientedCut, total: Symbol) -> Atom {
        let emr_to_lmb_cut = self.emr_to_lmb_and_cut(cut, total);

        let reps: Vec<_> = emr_to_lmb_cut
            .iter()
            .map(|(k, v)| Replacement::new(k, v))
            .collect();

        let mut denom = Atom::new_num(1);

        if let AtomView::Mul(m) = self.denominator.replace_all_multiple(&reps).as_view() {
            for a in m.iter() {
                if let AtomView::Fun(f) = a {
                    let mut builder = FunctionBuilder::new(f.get_symbol());
                    for arg in f.iter() {
                        builder = builder.add_arg(arg.expand().as_view())
                    }

                    denom = denom * builder.finish();
                }
            }
        }

        denom
    }

    pub fn emr_to_lmb_and_cut(
        &self,
        cut: &OrientedCut,
        total: Symbol,
    ) -> Vec<(Pattern, PatternOrMap)> {
        let (all_rest, solved_for) = self.cut_constraint(cut, total);

        let pattern = &solved_for.into_pattern();
        let rhs = &all_rest.into_pattern().into();

        let mut emr_to_lmb_cut = AHashMap::new();
        for (e, d) in self.graph.iter_egdes(&self.graph.full_graph()) {
            let data = d.data.unwrap();
            emr_to_lmb_cut.insert(
                fun!(symb!("Q"), data.bare_edge_id as i32),
                data.momentum
                    .replace_all(pattern, rhs, None, None)
                    .into_pattern(),
            );
        }

        emr_to_lmb_cut
            .into_iter()
            .map(|(k, v)| (k.into_pattern(), v.into()))
            .collect()
    }

    pub fn cut_constraint(&self, cut: &OrientedCut, total: Symbol) -> (Atom, Atom) {
        let mut sum = Atom::new_num(0);

        let mut total = Atom::new_var(total);
        let electron_momenta = fun!(symb!("k"), self.marked_electron_edge.1 as i32);
        // println!("p_e {}", electron_momenta);

        let photon_momenta = fun!(symb!("k"), self.lmb_photon.1 as i32);
        // println!("q:{}", photon_momenta);

        if let EdgeId::Paired { source, sink } = self.marked_electron_edge.0 {
            match cut.relative_orientation(source) {
                Orientation::Default => {
                    total = total + &electron_momenta;
                }
                Orientation::Reversed => {
                    total = total - &electron_momenta;
                }
                _ => {}
            }
        }

        for (o, cut_edge) in cut.iter_edges_relative(&self.graph) {
            // println!(
            //     "{}{}{}",
            //     SignOrZero::from(o),
            //     cut_edge.as_ref().data.unwrap().bare_edge_id,
            //     SignOrZero::from(o) * cut_edge.as_ref().data.unwrap().momentum.clone()
            // );
            sum = sum + SignOrZero::from(o) * cut_edge.as_ref().data.unwrap().momentum.clone();
        }

        let mut var = None;
        let mut all_rest = Atom::new_num(0);

        if let AtomView::Add(a) = sum.expand().as_view() {
            for e in a.iter() {
                if var.is_none() {
                    match e {
                        AtomView::Mul(m) => {
                            let mut iter = m.iter();
                            if let AtomView::Fun(v) = iter.next().unwrap() {
                                if photon_momenta.as_view() == v.as_view()
                                    || electron_momenta.as_view() == v.as_view()
                                {
                                    all_rest = all_rest + e;
                                } else {
                                    var = Some(e.to_owned());
                                }
                            } else {
                                panic!("{}", e)
                            }
                        }
                        AtomView::Fun(f) => {
                            if photon_momenta.as_view() == f.as_view()
                                || electron_momenta.as_view() == f.as_view()
                            {
                                all_rest = all_rest + e;
                            } else {
                                var = Some(e.to_owned());
                            }
                        }
                        _ => {
                            panic!("{}", e)
                        }
                    }
                } else {
                    all_rest = all_rest + e;
                }
            }
        }

        all_rest = total - all_rest;

        let (solved_for, coef) = match var.as_ref().unwrap().as_view() {
            AtomView::Mul(a) => {
                let mut solved = None;
                let mut coef = None;

                for i in a.iter() {
                    match i {
                        AtomView::Num(a) => match a.get_coeff_view().to_owned() {
                            symbolica::coefficient::Coefficient::Rational(a) => {
                                coef = Some(Coefficient::Rational(a.inv()));
                            }
                            _ => panic!("str"),
                        },
                        AtomView::Fun(f) => {
                            solved = Some(f.as_view().to_owned());
                        }
                        _ => panic!("str"),
                    }
                }
                (solved.unwrap(), Atom::new_num(coef.unwrap()))
            }
            AtomView::Fun(f) => (f.as_view().to_owned(), Atom::new_num(1)),
            _ => {
                panic!("should be a function or mul")
            }
        };

        // println!("coef:{coef}");

        // println!("all_rest: {}", (&all_rest * &coef).expand());
        // println!("solved_for: {}", solved_for);
        (all_rest * coef, solved_for)
    }
}

pub struct DisEdge {
    pub bare_edge_id: usize,
    pub bare_edge: Edge,
    marked: bool,
    momentum: Atom,
}

pub struct DisVertex {
    pub bare_vertex_id: usize,
    pub bare_vertex: Vertex,
}

#[cfg(test)]
mod test {
    use std::{fs, path::Path};

    use ahash::{AHashSet, HashMap, HashMapExt};
    use bitvec::vec::BitVec;
    use itertools::Itertools;
    use symbolica::symb;

    use crate::{
        feyngen::{
            diagram_generator::FeynGen, dis::Embeddings, FeynGenFilter, FeynGenOptions,
            GenerationType, SelfEnergyFilterOptions, SnailFilterOptions, TadpolesFilterOptions,
        },
        graph::{
            half_edge::{
                drawing::Decoration,
                layout::{FancySettings, LayoutParams, LayoutSettings, PositionalHedgeGraph},
                subgraph::{Cycle, OrientedCut, SubGraph, SubGraphOps},
                Flow, HedgeGraph, HedgeGraphBuilder, InvolutiveMapping, Orientation,
            },
            BareGraph,
        },
        model::Model,
        momentum::SignOrZero,
        tests_from_pytest::load_generic_model,
    };

    #[test]
    fn all_cuts_double_triangle() {
        let mut dt = HedgeGraphBuilder::new();

        let v1 = dt.add_node(());
        let v2 = dt.add_node(());
        let v3 = dt.add_node(());
        let v4 = dt.add_node(());

        dt.add_external_edge(v2, (), true, Flow::Sink);
        dt.add_external_edge(v1, (), true, Flow::Source);

        dt.add_edge(v1, v4, (), true);
        dt.add_edge(v1, v3, (), true);

        dt.add_edge(v4, v3, (), true);

        dt.add_edge(v3, v2, (), true);
        dt.add_edge(v4, v2, (), true);

        let g = dt.build();

        let mut ext = g.involution.iter().filter(|(_, (e, _))| e.is_identity());

        let s = ext.next().unwrap().1 .1.clone();
        let t = ext.next().unwrap().1 .1.clone();
        let mut cuts = AHashSet::new();

        g.all_s_t_cuts(&s, &t, &mut cuts);
        for c in cuts {
            println!("//{c:?}");
            println!("{}", g.dot(&c));
        }
    }

    #[test]
    fn all_cuts_tbt() {
        let mut dt = HedgeGraphBuilder::new();

        let v1 = dt.add_node(());
        let v2 = dt.add_node(());
        let v3 = dt.add_node(());
        let v4 = dt.add_node(());
        let v5 = dt.add_node(());
        let v6 = dt.add_node(());

        dt.add_external_edge(v1, (), true, Flow::Sink);
        dt.add_external_edge(v2, (), true, Flow::Source);

        dt.add_edge(v1, v4, (), true);
        dt.add_edge(v1, v3, (), true);

        dt.add_edge(v4, v3, (), true);

        dt.add_edge(v4, v5, (), true);
        dt.add_edge(v3, v6, (), true);

        dt.add_edge(v6, v5, (), true);

        dt.add_edge(v5, v2, (), true);
        dt.add_edge(v6, v2, (), true);

        let g = dt.build();

        let mut ext = g.involution.iter().filter(|(_, (e, _))| e.is_identity());

        let s = ext.next().unwrap().1 .1.clone();
        let t = ext.next().unwrap().1 .1.clone();

        let mut cuts = AHashSet::new();

        g.all_s_t_cuts(&s, &t, &mut cuts);
        println!("number of cuts {}", cuts.len());
        for c in cuts {
            println!("//{c:?}");
            println!("{}", g.dot(&c));
        }
    }

    #[test]
    fn all_cuts_tbbt() {
        let mut dt = HedgeGraphBuilder::new();

        let v1 = dt.add_node(());
        let v2 = dt.add_node(());
        let v3 = dt.add_node(());
        let v4 = dt.add_node(());
        let v5 = dt.add_node(());
        let v6 = dt.add_node(());
        let v7 = dt.add_node(());
        let v8 = dt.add_node(());

        dt.add_external_edge(v1, (), true, Flow::Sink);
        dt.add_external_edge(v2, (), true, Flow::Source);

        dt.add_edge(v1, v4, (), true);
        dt.add_edge(v1, v3, (), true);

        dt.add_edge(v4, v3, (), true);

        dt.add_edge(v4, v5, (), true);
        dt.add_edge(v3, v6, (), true);

        dt.add_edge(v6, v5, (), true);

        dt.add_edge(v5, v7, (), true);
        dt.add_edge(v6, v8, (), true);

        dt.add_edge(v7, v2, (), true);
        dt.add_edge(v8, v2, (), true);

        dt.add_edge(v7, v8, (), true);

        let g = dt.build();

        let mut ext = g.involution.iter().filter(|(_, (e, _))| e.is_identity());

        let s = ext.next().unwrap().1 .1.clone();
        let t = ext.next().unwrap().1 .1.clone();

        let mut cuts = AHashSet::new();

        g.all_s_t_cuts(&s, &t, &mut cuts);
        for c in &cuts {
            println!("//{c:?}");
            println!("{}", g.dot(c));
        }
        println!("number of cuts {}", cuts.len());
    }

    #[test]
    fn dis_generate() {
        let nloops = 2;
        let model = Model::from_file(String::from(
            Path::new("./src/test_resources")
                .join("gammaloop_models/sm.yaml")
                .to_str()
                .unwrap(),
        ))
        .unwrap();

        let mut coupling = HashMap::new();
        coupling.insert("QED".into(), 2);
        let options = FeynGenOptions {
            generation_type: GenerationType::CrossSection,
            initial_pdgs: vec![22],
            final_pdgs: vec![],
            loop_count_range: (nloops, nloops),
            symmetrize_final_states: true,
            symmetrize_initial_states: true,
            symmetrize_left_right_states: true,
            filters: crate::feyngen::FeynGenFilters(vec![
                FeynGenFilter::ParticleVeto(vec![
                    23, 24, 9000001, 9000002, 9000003, 9000004, 12, 14, 16, 2, 4, 6, 3, 5, 25, 250,
                    251, 11, 13, 15,
                ]),
                FeynGenFilter::SelfEnergyFilter(SelfEnergyFilterOptions::default()),
                FeynGenFilter::ZeroSnailsFilter(SnailFilterOptions::default()),
                FeynGenFilter::TadpolesFilter(TadpolesFilterOptions::default()),
                FeynGenFilter::CouplingOrders(coupling),
            ]),
        };
        let diagram_gen = FeynGen::new(options);

        let diagrams = diagram_gen
            .generate(
                &model,
                crate::feyngen::NumeratorAwareGraphGroupingOption::OnlyDetectZeroes,
                true,
                "DIS".into(),
                None,
                None,
                None,
            )
            .unwrap();

        for (j, d) in diagrams.iter().enumerate() {
            println!("{}", d.dot());
            let mut h = HedgeGraph::close_externals(d);
            println!("involution:{}", h.involution);

            let (basis, tree) = h.cycle_basis();
            h.align_to_tree_underlying(&tree);
            let settings = LayoutSettings::new(&h, LayoutParams::default(), 0, 10000, 100.);
            let init = h
                .clone()
                .layout(settings)
                .cetz_impl(&|e| e.to_string(), &|&e| d.edges[e].particle.decoration());

            println!("{init}");

            let cuts = OrientedCut::all_initial_state_cuts(&h);

            println!("N cuts{}", cuts.len());

            let basis = basis
                .into_iter()
                .map(|c| Cycle::new_circuit(c.filter, &h).unwrap())
                .collect::<Vec<_>>();

            let embeddings = Embeddings::classify(
                cuts,
                basis,
                |c| {
                    // if c.cut.count_ones() == 1 {
                    //     let e = c.cut.iter_ones().next().unwrap();
                    //     d.edges[e].particle.pdg_code.abs() != 22
                    // } else {
                    //     false
                    // }
                    true
                },
                true,
            );

            println!("N embeddings{}", embeddings.n_embeddings());

            let params = LayoutParams::default();

            let mut layouts = Vec::new();

            for (i, (e, cuts)) in embeddings.cuts.into_iter().enumerate() {
                println!(
                    "//Embedding {} {:?}",
                    i,
                    e.windings_in_field((nloops + 1) as u32)
                );

                let layouts_emb: Vec<_> = cuts
                    .into_iter()
                    .map(|c| {
                        let mut l = c.layout(&h, params, 0, 10000, 100., 20.);
                        // l.to_fancy();
                        l
                    })
                    .collect();

                layouts.push((
                    format!("embedding{}", i + 1),
                    format!("embedding {} {:?}", i + 1, e.windings),
                    layouts_emb,
                ));
            }
            fs::write(
                &format!("embeddings{j}.typ"),
                PositionalHedgeGraph::cetz_impl_collection(
                    &layouts,
                    &|(e, o)| e.to_string(),
                    &|(e, o)| d.edges[**e].particle.decoration(),
                ),
            )
            .unwrap();
            // println!("embeddings:{:#?}", embeddings);
        }

        println!("Number of diagrams {}", diagrams.len());
    }

    #[test]
    fn double_bananna() {
        let mut dt = HedgeGraphBuilder::new();

        let v1 = dt.add_node(());
        let v2 = dt.add_node(());

        dt.add_edge(v1, v2, Decoration::None, false);
        dt.add_edge(v1, v2, Decoration::Wave, false);
        dt.add_edge(v1, v2, Decoration::Arrow, false);

        let mut h = dt.build();

        let (basis, tree) = h.cycle_basis();
        h.align_to_tree_underlying(&tree);
        let cuts = OrientedCut::all_initial_state_cuts(&h);

        let basis = basis
            .into_iter()
            .map(|c| Cycle::new_circuit(c.filter, &h).unwrap())
            .collect::<Vec<_>>();

        let embeddings = Embeddings::classify(cuts, basis, |_| true, true);

        println!("N embeddings{}", embeddings.n_embeddings());

        let params = LayoutParams::default();

        let mut layouts = Vec::new();

        println!("involution {}", h.involution);

        for (i, (e, cuts)) in embeddings.cuts.into_iter().enumerate() {
            println!("Embedding {i}:{:?}", e.windings);
            // for c in &cuts {
            //     println!("{}", c)
            // }
            let layouts_emb: Vec<_> = cuts
                .into_iter()
                .map(|c| {
                    let mut l = c.layout(&h, params, 0, 10000, 100., 16.);
                    // l.to_fancy();
                    l
                })
                .collect();

            layouts.push((
                format!("embedding{}", i + 1),
                format!("embedding {} {:?}", i + 1, e.windings),
                layouts_emb,
            ));
        }
        fs::write(
            &format!("embeddingsbananna.typ"),
            PositionalHedgeGraph::cetz_impl_collection(
                &layouts,
                &|(e, o)| e.to_cetz(),
                &|(e, o)| **e,
            ),
        )
        .unwrap();
    }

    #[test]
    fn bubble() {
        let mut dt = HedgeGraphBuilder::new();

        let v1 = dt.add_node(());
        let v2 = dt.add_node(());

        dt.add_edge(v1, v2, Decoration::Arrow, true);
        dt.add_edge(v1, v2, Decoration::Arrow, Orientation::Reversed);

        let mut h = dt.build();
        let mut layouts = Vec::new();

        for i in 0..20 {
            let settings = LayoutSettings::new(&h, LayoutParams::default(), 0, 10000, 100.);
            layouts.push(h.clone().layout(settings));
        }

        println!("{}", h.base_dot());

        let (basis, tree) = h.cycle_basis();
        h.align_to_tree_underlying(&tree);

        let mut layouts_after = Vec::new();

        for i in 0..20 {
            let settings = LayoutSettings::new(&h, LayoutParams::default(), 0, 10000, 100.);

            layouts_after.push(h.clone().layout(settings));
        }

        fs::write(
            &format!("bubblebefore_and_after.typ"),
            PositionalHedgeGraph::cetz_impl_collection(
                &[
                    ("before".into(), "be".into(), layouts),
                    ("after".into(), "rts".into(), layouts_after),
                ],
                &|e| "".into(),
                &|e| *e,
            ),
        )
        .unwrap();

        let cuts = OrientedCut::all_initial_state_cuts(&h);
        println!("{}", h.base_dot());

        println!("number of cuts: {}", cuts.len());

        let basis = basis
            .into_iter()
            .map(|c| Cycle::new_circuit(c.filter, &h).unwrap())
            .collect::<Vec<_>>();

        let embeddings = Embeddings::classify(cuts, basis, |_| true, true);

        println!("number of embeddings: {}", embeddings.n_embeddings());
        for (e, c) in embeddings.cuts {
            for c in c {
                println!("n-cut {}", c.sign.count_ones());
                for (o, e) in c
                    .sign
                    .included_iter()
                    .map(|i| (c.relative_orientation(i), h.involution.edge_data(i)))
                {
                    println!("{} {}", SignOrZero::from(o), e.data.unwrap().to_cetz())
                }
            }
            println!("{:?}", e.windings);
        }

        println!("{}", h.base_dot());
    }

    #[test]
    fn self_energy_closed() {
        let mut dt = HedgeGraphBuilder::new();

        let v1 = dt.add_node(());
        let v2 = dt.add_node(());
        let v3 = dt.add_node(());
        let v4 = dt.add_node(());

        dt.add_edge(v1, v2, (), false);
        dt.add_edge(v1, v2, (), false);
        dt.add_edge(v3, v4, (), false);
        dt.add_edge(v3, v4, (), false);

        dt.add_edge(v3, v1, (), false);
        dt.add_edge(v2, v4, (), false);

        let h = dt.build();

        let mut bv: BitVec = BitVec::new();
        bv.push(true);
        bv.push(true);
        bv.push(true);
        bv.push(true);
        bv.push(true);
        bv.push(true);
        bv.push(true);
        bv.push(true);
        bv.push(false);
        bv.push(false);
        bv.push(false);
        bv.push(false);

        println!("{}", h.dot(&bv));
        // println!(":{}", h.count_connected_components(&bv));
        let cuts = OrientedCut::all_initial_state_cuts(&h);

        for c in cuts {
            let mut cut = c.sign.clone();

            for i in c.sign.included_iter() {
                cut.set(h.involution.inv(i).0, true);
            }
            let owned = c.to_owned_graph(&h);

            let complement = cut.complement(&h);

            let filter_comp_count = h.count_connected_components(&complement);
            let owned_comp_count = owned.count_connected_components(&owned.full_filter());

            let covers = complement.covers(&h) == h.full_filter();

            if filter_comp_count != owned_comp_count && covers {
                println!("//owned{owned_comp_count}:\n{}", owned.base_dot());
                println!("//cut{filter_comp_count}:\n{}", h.dot(&cut));
            }
        }
    }

    #[test]
    fn self_energy_closed_epem() {
        let mut dt = HedgeGraphBuilder::new();

        let v1 = dt.add_node(());
        let v2 = dt.add_node(());
        let v3 = dt.add_node(());
        let v4 = dt.add_node(());

        let v5 = dt.add_node(());
        let v6 = dt.add_node(());

        dt.add_edge(v1, v2, (0, Decoration::Arrow, 0), true);
        dt.add_edge(v1, v2, (0, Decoration::Coil, 1), false);
        dt.add_edge(v2, v4, (0, Decoration::Arrow, 2), true);
        dt.add_edge(v4, v3, (0, Decoration::Arrow, 3), true);
        dt.add_edge(v3, v1, (0, Decoration::Arrow, 4), true);

        dt.add_edge(v3, v5, (-1, Decoration::Wave, 5), false);
        dt.add_edge(v4, v6, (-1, Decoration::Wave, 6), false);

        dt.add_edge(v5, v6, (1, Decoration::Arrow, 7), true);
        dt.add_edge(v6, v5, (2, Decoration::Arrow, 8), true);

        let mut h = dt.build();

        let layout = h.clone().layout(LayoutSettings::new(
            &h,
            LayoutParams::default(),
            8,
            10000,
            100.,
        ));

        println!("{}", layout.cetz_impl(&|_| "".into(), &|d| d.1));

        // println!(":{}", h.count_connected_components(&bv));
        let (basis, tree) = h.cycle_basis();
        h.align_to_tree_underlying(&tree);

        println!(
            "{}",
            h.dot_impl(
                &h.full_filter(),
                "".into(),
                &|e| Some(format!("\"{}:{}:{}\"", e.0, e.1.to_cetz(), e.2)),
                &|n| None
            )
        );
        let cuts = OrientedCut::all_initial_state_cuts(&h);

        let mut electron_cut = h.empty_filter();

        for (i, (m, n)) in h.involution.iter() {
            match m {
                InvolutiveMapping::Sink { source_idx } => {}
                InvolutiveMapping::Source { data, sink_idx } => {
                    if data.data.unwrap().0 > 0 {
                        electron_cut.set(i.0, true);
                    }
                }
                _ => {}
            }
        }

        let mut filtered_cuts: Vec<_> = cuts
            .iter()
            .filter(|&c| {
                let mut has_cut_electron = false;
                let mut has_cut_photon = false;
                let mut has_cut_parton = false;

                let with_e_cut = c.reference.union(&electron_cut).complement(&h);

                let disconnecting_e = h.count_connected_components(&with_e_cut) > 1;

                let mut orientation_matches = false;
                for e in c.included_iter() {
                    if h.get_edge_data(e).0 > 0 {
                        orientation_matches =
                            c.relative_orientation(e) == h.involution.edge_data(e).orientation;
                    }
                }

                for e in h.iter_egde_data(c) {
                    match e.data.unwrap().0 {
                        -1 => has_cut_photon = true,
                        e if e > 0 => {
                            has_cut_electron = !has_cut_electron;
                        }
                        0 => {
                            has_cut_parton = true;
                        }
                        _ => {}
                    }
                }
                !has_cut_photon
                    && has_cut_electron
                    && has_cut_parton
                    && !disconnecting_e
                    && orientation_matches
            })
            .cloned()
            .collect();

        println!("n cuts {}", filtered_cuts.len());

        for (i, c) in basis.iter().enumerate() {
            println!("Cycle {i}:");

            for e in h.iter_internal_edge_data(&c.filter) {
                println!("{}", e.data.unwrap().2);
            }
        }

        filtered_cuts.push(OrientedCut::empty(h.n_hedges()));
        let mut embeddings = Embeddings::classify(filtered_cuts, basis, |_| true, true);
        embeddings.remove_singles();

        let params = LayoutParams::default();

        let mut layouts = Vec::new();

        println!("N embeddings{}", embeddings.n_embeddings());
        let mut ifsplit = embeddings.if_split(&h, &|e| e.0 == 1);

        ifsplit.remove_empty();

        for (i, (e, cuts)) in ifsplit.cuts.iter().enumerate() {
            println!("Embedding {i}:{:?}", e.windings);

            let layouts_emb_i: Vec<_> = cuts[0]
                .iter()
                .map(|c| {
                    let mut l = c.clone().layout(&h, params, 0, 10000, 100., 20.);
                    // l.to_fancy();
                    l
                })
                .collect();

            let layouts_emb_f: Vec<_> = cuts[1]
                .iter()
                .map(|c| {
                    let mut l = c.clone().layout(&h, params, 0, 10000, 100., 20.);
                    // l.to_fancy();
                    l
                })
                .collect();

            layouts.push((
                format!("embedding{}i", i + 1),
                format!("embedding {} {:?} initial", i + 1, e.windings),
                layouts_emb_i,
            ));

            layouts.push((
                format!("embedding{}f", i + 1),
                format!("embedding {} {:?} final", i + 1, e.windings),
                layouts_emb_f,
            ));
        }
        fs::write(
            &format!("embeddingselfenergyepem.typ"),
            PositionalHedgeGraph::cetz_impl_collection(
                &layouts,
                &|((_, _, n), o)| n.to_string(),
                &|&((i, d, _), o)| *d,
            ),
        )
        .unwrap();
    }

    #[test]
    fn double_triangle_dis() {
        let model = load_generic_model("sm");
        let mut symbolica_graph = symbolica::graph::Graph::new();

        let l1 = symbolica_graph.add_node((0, "V_98".into()));
        let l2 = symbolica_graph.add_node((0, "V_98".into()));
        let l3 = symbolica_graph.add_node((0, "V_71".into()));
        let l4 = symbolica_graph.add_node((0, "V_71".into()));
        let l5 = symbolica_graph.add_node((0, "V_74".into()));
        let l6 = symbolica_graph.add_node((0, "V_74".into()));
        symbolica_graph.add_edge(l1, l2, true, "e-");
        symbolica_graph.add_edge(l2, l1, true, "e-");
        symbolica_graph.add_edge(l2, l4, true, "a");
        symbolica_graph.add_edge(l1, l3, true, "a");

        symbolica_graph.add_edge(l3, l6, true, "d");
        symbolica_graph.add_edge(l6, l4, true, "d");
        symbolica_graph.add_edge(l4, l5, true, "d");
        symbolica_graph.add_edge(l5, l3, true, "d");

        symbolica_graph.add_edge(l5, l6, true, "g");
        let bare_graph = BareGraph::from_symbolica_graph(
            &model,
            "disdoubletriangle".into(),
            &symbolica_graph,
            "1".into(),
            vec![],
            None,
        )
        .unwrap();

        let dis_graph = bare_graph.dis_graph();

        let h = &dis_graph.graph;
        let cuts = OrientedCut::all_initial_state_cuts(h);

        let mut electron_cut = dis_graph.graph.empty_filter();

        for (i, (m, n)) in h.involution.iter() {
            match m {
                InvolutiveMapping::Sink { source_idx } => {}
                InvolutiveMapping::Source { data, sink_idx } => {
                    if data
                        .as_ref()
                        .data
                        .unwrap()
                        .bare_edge
                        .particle
                        .pdg_code
                        .abs()
                        == 11
                    {
                        electron_cut.set(i.0, true);
                    }
                }
                _ => {}
            }
        }

        let mut filtered_cuts: Vec<_> = cuts
            .iter()
            .filter(|&c| {
                let mut has_cut_electron = false;
                let mut has_cut_photon = false;
                let mut has_cut_parton = false;

                let with_e_cut = c.reference.union(&electron_cut).complement(&h);

                let disconnecting_e = h.count_connected_components(&with_e_cut) > 1;

                let mut orientation_matches = false;
                for e in c.included_iter() {
                    if h.get_edge_data(e).bare_edge.particle.pdg_code.abs() == 11 {
                        orientation_matches =
                            c.relative_orientation(e) == h.involution.edge_data(e).orientation;
                    }
                }

                for e in h.iter_egde_data(c) {
                    match e.data.unwrap().bare_edge.particle.pdg_code.abs() {
                        22 => has_cut_photon = true,
                        11 => {
                            has_cut_electron = !has_cut_electron;
                        }
                        _ => {
                            has_cut_parton = true;
                        }
                    }
                }
                !has_cut_photon
                    && has_cut_electron
                    && has_cut_parton
                    && !disconnecting_e
                    && orientation_matches
            })
            .cloned()
            .collect();

        println!("n cuts {}", filtered_cuts.len());

        // for (i, c) in basis.iter().enumerate() {
        //     println!("Cycle {i}:");

        //     for e in h.iter_internal_edge_data(&c.filter) {
        //         println!("{}", e.data.unwrap().2);
        //     }
        // }

        filtered_cuts.push(OrientedCut::empty(h.n_hedges()));
        let mut embeddings =
            Embeddings::classify(filtered_cuts, dis_graph.basis.clone(), |_| true, true);
        embeddings.remove_singles();

        let params = LayoutParams::default();

        let mut layouts = Vec::new();

        println!("N embeddings{}", embeddings.n_embeddings());
        let mut ifsplit = embeddings.if_split(&h, &|e| e.marked);

        ifsplit.remove_empty();

        for (i, (e, cuts)) in ifsplit.cuts.iter().enumerate() {
            println!("Embedding {}:{:?}", i + 1, e.windings);

            let first_initial = &cuts[0][0];

            let denom = dis_graph.denominator(first_initial, symb!("p"));
            let num: Vec<_> = dis_graph
                .numerator(first_initial, symb!("p"))
                .iter()
                .map(|a| a.expand())
                .collect();

            // let (all_rest, solved_for) = dis_graph.cut_constraint(first_initial, symb!("p"));

            // let pattern = &solved_for.into_pattern();
            // let rhs = &all_rest.into_pattern().into();

            println!(
                "Denominator:{} ",
                denom.printer(symbolica::printer::PrintOptions::mathematica())
            );

            println!(
                "W1:{} ",
                num[0].printer(symbolica::printer::PrintOptions::mathematica())
            );

            println!(
                "W2:{} ",
                num[1].printer(symbolica::printer::PrintOptions::mathematica())
            );
            let fancy_settings = FancySettings {
                label_shift: 0.06,
                arrow_angle_percentage: Some(0.7),
                arrow_shift: 0.06,
            };

            let layouts_emb_i: Vec<_> = cuts[0]
                .iter()
                .map(|c| {
                    let (all_rest, solved_for) = dis_graph.cut_constraint(c, symb!("p"));

                    let pattern = &solved_for.into_pattern();
                    let rhs = &all_rest.into_pattern().into();

                    let mut l = c.clone().layout(&h, params, 0, 10000, 100., 19.);

                    let mut l = l.map(
                        |e| e,
                        |e, d| {
                            d.map(|d| {
                                d.map(|(e, o)| {
                                    (
                                        e.momentum.clone(),
                                        e.momentum.clone(),
                                        e.bare_edge.particle.decoration(),
                                        o,
                                    )
                                })
                            })
                        },
                    );
                    for i in l.involution.iter_idx() {
                        let mom = &mut l.involution.edge_data_mut(i).data.as_mut().unwrap().data.0;

                        *mom = mom.replace_all(pattern, rhs, None, None).expand();
                    }

                    l.to_fancy(&fancy_settings);
                    l
                })
                .collect();

            let layouts_emb_f: Vec<_> = cuts[1]
                .iter()
                .map(|c| {
                    let (all_rest, solved_for) = dis_graph.cut_constraint(c, symb!("p"));

                    let pattern = &solved_for.into_pattern();
                    let rhs = &all_rest.into_pattern().into();

                    let mut l = c.clone().layout(&h, params, 0, 10000, 100., 19.);
                    let mut l = l.map(
                        |e| e,
                        |e, d| {
                            d.map(|d| {
                                d.map(|(e, o)| {
                                    (
                                        e.momentum.clone(),
                                        e.momentum.clone(),
                                        e.bare_edge.particle.decoration(),
                                        o,
                                    )
                                })
                            })
                        },
                    );
                    for i in l.involution.iter_idx() {
                        let mom = &mut l.involution.edge_data_mut(i).data.as_mut().unwrap().data.0;

                        *mom = mom.replace_all(pattern, rhs, None, None).expand();
                    }
                    l.to_fancy(&fancy_settings);
                    l
                })
                .collect();

            layouts.push((
                format!("embedding{}i", i + 1),
                format!("embedding {} {:?} initial", i + 1, e.windings),
                layouts_emb_i,
            ));

            layouts.push((
                format!("embedding{}f", i + 1),
                format!("embedding {} {:?} final", i + 1, e.windings),
                layouts_emb_f,
            ));
        }
        fs::write(
            &format!("doubletriangle.typ"),
            PositionalHedgeGraph::cetz_impl_collection(
                &layouts,
                &|(e, f, _, o)| {
                    if let Orientation::Reversed = o {
                        (-e.clone().expand()).to_string()
                    } else {
                        e.expand().to_string()
                    }
                },
                &|(e, _, o, d)| *o,
            ),
        )
        .unwrap();
    }

    #[test]
    fn signature_epem() {
        let mut dt = HedgeGraphBuilder::new();

        let v1 = dt.add_node(());
        let v2 = dt.add_node(());
        let v3 = dt.add_node(());
        let v4 = dt.add_node(());

        let v5 = dt.add_node(());
        let v6 = dt.add_node(());

        dt.add_edge(v1, v2, (0, Decoration::Arrow, 0), true);
        dt.add_edge(v1, v2, (0, Decoration::Coil, 1), false);
        dt.add_edge(v2, v4, (0, Decoration::Arrow, 2), true);
        dt.add_edge(v4, v3, (0, Decoration::Arrow, 3), true);
        dt.add_edge(v3, v1, (0, Decoration::Arrow, 4), true);

        dt.add_edge(v3, v5, (-1, Decoration::Wave, 5), false);
        dt.add_edge(v4, v6, (-1, Decoration::Wave, 6), false);

        dt.add_edge(v5, v6, (1, Decoration::Arrow, 7), true);
        dt.add_edge(v6, v5, (2, Decoration::Arrow, 8), true);

        let mut h = dt.build();
        let (basis, tree) = h.cycle_basis();
        h.align_to_tree_underlying(&tree);
        let mut full_cut = h.empty_filter();

        for (i, (m, n)) in h.involution.iter() {
            match h.involution.edge_data(i).data.unwrap().2 {
                3 | 8 => full_cut.set(i.0, true),
                _ => {}
            }
        }

        let mut pp = h.empty_filter();
        let mut pm = h.empty_filter();
        let mut mp = h.empty_filter();
        let mut mm = h.empty_filter();
        let mut all_sources = h.empty_filter();

        let mut first = true;
        let mut first_source = true;

        for i in full_cut.included_iter() {
            match &h.involution[i] {
                InvolutiveMapping::Identity { data, underlying } => {}
                InvolutiveMapping::Sink { source_idx } => {
                    if first {
                        mp.set(i.0, true);
                        first = false;
                    } else {
                        pm.set(i.0, true)
                    }
                    mm.set(i.0, true);
                }
                InvolutiveMapping::Source { data, sink_idx } => {
                    all_sources.set(i.0, true);
                    if first_source {
                        pm.set(i.0, true);
                        first_source = false;
                    } else {
                        mp.set(i.0, true);
                    }
                    pp.set(i.0, true);
                }
            }
        }

        let all_signed: Vec<_> = [pp, pm, mp, mm]
            .into_iter()
            .map(|c| {
                println!("Cut {:?}", c);

                let cut = OrientedCut {
                    reference: all_sources.clone(),
                    sign: c,
                };
                for i in cut.included_iter() {
                    print!(
                        "{}{}",
                        SignOrZero::from(cut.relative_orientation(i)),
                        h.involution.edge_data(i).data.unwrap().2
                    )
                }
                cut
            })
            .collect();
        let params = LayoutParams::default();
        let fancy_settings = FancySettings {
            label_shift: 0.1,
            arrow_angle_percentage: Some(0.7),
            arrow_shift: 0.1,
        };

        let layouts_emb: Vec<_> = all_signed
            .clone()
            .into_iter()
            .map(|c| {
                let mut l = c.layout(&h, params, 0, 10000, 100., 16.);
                l.to_fancy(&fancy_settings);
                l
            })
            .collect();

        fs::write(
            &format!("test.typ"),
            PositionalHedgeGraph::cetz_impl_collection(
                &[("tsrt".into(), "rt".into(), layouts_emb)],
                &|((_, _, n), o)| n.to_string(),
                &|&((i, d, _), o)| *d,
            ),
        )
        .unwrap();
        let mut layouts = Vec::new();

        let embeddings = Embeddings::classify(all_signed, basis, |_| true, true);

        for (i, (e, cuts)) in embeddings.cuts.into_iter().enumerate() {
            println!("//Embedding {i}:{:?}", e.windings);

            let layouts_emb: Vec<_> = cuts
                .into_iter()
                .map(|c| {
                    println!("");
                    print!("//");
                    for i in c.reference.included_iter() {
                        print!(
                            "{}{}",
                            SignOrZero::from(c.relative_orientation(i)),
                            h.involution.edge_data(i).data.unwrap().2
                        )
                    }
                    println!("");

                    let mut l = c.layout(&h, params, 0, 10000, 100., 16.);
                    // println!(
                    //     "{}",
                    //     l.dot_impl(
                    //         &l.full_filter(),
                    //         "srt".into(),
                    //         &|e| Some(format!("label=\"{}\"", e.data.2)),
                    //         &|_| { None }
                    //     )
                    // );
                    l.to_fancy(&fancy_settings);
                    l
                })
                .collect();

            layouts.push((
                format!("embedding{}", i + 1),
                format!("embedding {} {:?}", i + 1, e.windings),
                layouts_emb,
            ));
        }
        fs::write(
            &format!("testemb.typ"),
            PositionalHedgeGraph::cetz_impl_collection(
                &layouts,
                &|((_, _, n), o)| n.to_string(),
                &|&((i, d, _), o)| *d,
            ),
        )
        .unwrap();
    }

    #[test]
    fn self_energy_1_loop_closed_epem() {
        let mut dt = HedgeGraphBuilder::new();

        let v3 = dt.add_node(());
        let v4 = dt.add_node(());

        let v5 = dt.add_node(());
        let v6 = dt.add_node(());

        dt.add_edge(v3, v4, (0, Decoration::Arrow, 2), true);
        dt.add_edge(v4, v3, (0, Decoration::Arrow, 3), true);

        dt.add_edge(v3, v5, (-1, Decoration::Wave, 5), false);
        dt.add_edge(v4, v6, (-1, Decoration::Wave, 6), false);

        dt.add_edge(v5, v6, (1, Decoration::Arrow, 7), true);
        dt.add_edge(v6, v5, (2, Decoration::Arrow, 8), true);

        let mut h = dt.build();

        let layout = h.clone().layout(LayoutSettings::new(
            &h,
            LayoutParams::default(),
            8,
            10000,
            100.,
        ));

        println!("{}", layout.cetz_impl(&|_| "".into(), &|d| d.1));

        // println!(":{}", h.count_connected_components(&bv));
        let (basis, tree) = h.cycle_basis();
        h.align_to_tree_underlying(&tree);

        println!(
            "{}",
            h.dot_impl(
                &h.full_filter(),
                "".into(),
                &|e| Some(format!("label=\"{}\"", e.2)),
                &|n| None
            )
        );

        println!(
            "{}",
            h.dot_impl(
                &tree.tree,
                "".into(),
                &|e| Some(format!("label=\"{}\"", e.2)),
                &|n| None
            )
        );
        let cuts = OrientedCut::all_initial_state_cuts(&h);

        let mut electron_cut = h.empty_filter();

        for (i, (m, n)) in h.involution.iter() {
            match m {
                InvolutiveMapping::Sink { source_idx } => {}
                InvolutiveMapping::Source { data, sink_idx } => {
                    if data.data.unwrap().0 > 0 {
                        electron_cut.set(i.0, true);
                    }
                }
                _ => {}
            }
        }

        let mut filtered_cuts: Vec<_> = cuts
            .iter()
            .filter(|&c| {
                let mut has_cut_electron = false;
                let mut has_cut_photon = false;
                let mut has_cut_parton = false;

                let with_e_cut = c.reference.union(&electron_cut).complement(&h);

                let disconnecting_e = h.count_connected_components(&with_e_cut) > 1;

                let mut orientation_matches = false;
                for e in c.included_iter() {
                    if h.get_edge_data(e).0 > 0 {
                        orientation_matches =
                            c.relative_orientation(e) == h.involution.edge_data(e).orientation;
                    }
                }

                for e in h.iter_egde_data(c) {
                    match e.data.unwrap().0 {
                        -1 => has_cut_photon = true,
                        e if e > 0 => {
                            has_cut_electron = !has_cut_electron;
                        }
                        0 => {
                            has_cut_parton = true;
                        }
                        _ => {}
                    }
                }
                !has_cut_photon
                    && has_cut_electron
                    && has_cut_parton
                    && !disconnecting_e
                    && orientation_matches
            })
            .cloned()
            .collect();

        println!("n cuts {}", filtered_cuts.len());

        for (i, c) in basis.iter().enumerate() {
            println!("Cycle {i}:");

            for e in h.iter_internal_edge_data(&c.filter) {
                println!("{}", e.data.unwrap().2);
            }
        }

        filtered_cuts.push(OrientedCut::empty(h.n_hedges()));
        let embeddings = Embeddings::classify(filtered_cuts, basis, |_| true, true);

        println!("N embeddings{}", embeddings.n_embeddings());

        let params = LayoutParams::default();

        let mut layouts = Vec::new();

        let embeddings = embeddings
            .cuts
            .into_iter()
            .sorted()
            .filter(|(e, l)| l.len() > 1)
            .collect::<Vec<_>>();

        for (i, (e, cuts)) in embeddings.into_iter().enumerate() {
            println!("//Embedding {i}:{:?}", e.windings);

            let layouts_emb: Vec<_> = cuts
                .into_iter()
                .map(|c| {
                    println!("");
                    print!("//");
                    for i in c.included_iter() {
                        print!(
                            "{}{}",
                            SignOrZero::from(c.relative_orientation(i)),
                            h.involution.edge_data(i).data.unwrap().2
                        )
                    }
                    println!("");

                    let mut l = c.layout(&h, params, 0, 10000, 100., 16.);
                    // println!(
                    //     "{}",
                    //     l.dot_impl(
                    //         &l.full_filter(),
                    //         "srt".into(),
                    //         &|e| Some(format!("label=\"{}\"", e.data.2)),
                    //         &|_| { None }
                    //     )
                    // );
                    // l.to_fancy();
                    l
                })
                .collect();

            layouts.push((
                format!("embedding{}", i + 1),
                format!("embedding {} {:?}", i + 1, e.windings),
                layouts_emb,
            ));
        }
        fs::write(
            &format!("embeddingselfenergyepem1l.typ"),
            PositionalHedgeGraph::cetz_impl_collection(
                &layouts,
                &|((_, _, n), o)| n.to_string(),
                &|&((i, d, _), o)| *d,
            ),
        )
        .unwrap();
    }
}
