use idenso::IndexTooling;
use idenso::color::{CS, ColorSimplifier, SelectiveExpand};
use idenso::gamma::{AGS, GammaSimplifier};
use indicatif::ProgressBar;
use indicatif::{ParallelProgressIterator, ProgressStyle};

use linnet::half_edge::tree::SimpleTraversalTree;
use linnet::permutation::Permutation;
use rayon::ThreadPoolBuilder;
use rayon::iter::ParallelIterator;
use rayon::prelude::*;

use spenso::network::library::LibraryTensor;
use spenso::network::library::symbolic::{ExplicitKey, TensorLibrary};
use spenso::network::parsing::ParseSettings;
use spenso::network::{Sequential, SmallestDegree};

// use spenso::network::Network;

// use spenso::shadowing::symbolica_utils::AtomCoreExt;
use spenso::structure::representation::{LibraryRep, Minkowski, RepName};
use spenso::structure::{PermutedStructure, TensorStructure};
use spenso::tensors::data::DataTensor;
use spenso::tensors::parametric::ParamTensor;
use spenso_hep_lib::{gamma_data_weyl, gamma5_weyl_data, proj_m_data_weyl, proj_p_data_weyl};
use std::collections::hash_map::Entry;
use std::collections::{HashSet, VecDeque};

use std::ops::{Deref, RangeInclusive};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;
use symbolica::atom::AtomView;
use symbolica::coefficient::Coefficient;
use symbolica::domains::algebraic_number::AlgebraicExtension;
use symbolica::domains::finite_field::PrimeIteratorU64;
use symbolica::domains::float::Complex as SymbolicaComplex;
use symbolica::function;
use symbolica::graph::{GenerationSettings, HalfEdge};
use symbolica::id::Replacement;
use tracing::{error, event_enabled, info, instrument};

use ahash::AHashMap;
use ahash::AHashSet;
use ahash::HashMap;
use colored::Colorize;
use smartstring::{LazyCompact, SmartString};
use symbolica::atom::AtomCore;
use symbolica::{parse, symbol};
use tracing::debug;
use tracing::warn;

use super::FeynGenError;
use super::NumeratorAwareGraphGroupingOption;
use super::SelfEnergyFilterOptions;
use super::SnailFilterOptions;
use super::TadpolesFilterOptions;
use crate::graph::ext::HedgeGraphExt;
use crate::graph::parse::ParseGraph;
use crate::graph::parse::string_utils::ToOrderedSimple;
use crate::graph::{FeynmanGraph, Graph, LMBext};
use crate::model::ArcVertexRule;
use crate::model::VertexRule;
use crate::model::{ArcParticle, ColorStructure};
use crate::momentum::{Pow, Sign, SignOrZero};
use crate::momentum_sample::LoopIndex;
use crate::numerator::ParamParsingNet;
use crate::numerator::aind::Aind;
use crate::numerator::graph::ReversibleEdge;
use crate::numerator::symbolica_ext::AtomCoreExt;
use crate::processes::ProcessDefinition;
use crate::settings::GlobalSettings;
use crate::utils::symbolica_ext::{COMPLEXRATPOLYFIELD, LOGPRINTOPTS, PrimeGenerate, Q_I};
use crate::utils::{self, GS, PARAM_FUN_LIB, W_};
use crate::uv::UltravioletGraph;
use crate::{INTERRUPTED, is_interrupted, set_interrupted};
use crate::{
    feyngen::{FeynGenFilter, GenerationType},
    model::Model,
};
use eyre::eyre;

type NumeratorSample = (
    Vec<Replacement>,
    TensorLibrary<ParamTensor<ExplicitKey<Aind>>, Aind>,
);
type RationalPoly = symbolica::poly::polynomial::MultivariatePolynomial<
    AlgebraicExtension<
        symbolica::domains::rational::FractionField<symbolica::domains::integer::IntegerRing>,
    >,
>;
type SampleEvaluationsAsPolynomial = (Vec<RationalPoly>, Vec<bool>);

use color_eyre::Result;
use itertools::Itertools;
use linnet::half_edge::HedgeGraph;
use linnet::half_edge::NodeIndex;
use linnet::half_edge::involution::{EdgeData, EdgeIndex, Flow, Orientation};
use linnet::half_edge::subgraph::{
    InternalSubGraph, ModifySubSet, OrientedCut, SuBitGraph, SubSetLike, SubSetOps,
};
use symbolica::{atom::Atom, graph::Graph as SymbolicaGraph};

const CANONIZE_GRAPH_FLOWS: bool = true;
const ANALYZE_RATIO_AS_RATIONAL_POLYNOMIAL: bool = true;
const EXPAND_NUMERICAL_SAMPLES_BEFORE_COMPARISON: bool = false;
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeColorWithVertexRule {
    pub external_tag: i32,
    pub vertex_rule: ArcVertexRule,
}

impl NodeColorWithVertexRule {
    pub fn from_particles<B: AsRef<str>, A: IntoIterator<Item = B>>(
        iter: A,
        model: &Model,
    ) -> Self {
        let mut particles: Vec<ArcParticle> = iter
            .into_iter()
            .map(|p| model.get_particle(p.as_ref()))
            .collect();
        particles.sort_by_key(|p| p.0.pdg_code);
        let vertex_rule = model.particle_set_to_vertex_rules_map[&particles][0].clone();
        Self {
            external_tag: 0,
            vertex_rule,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Default, Copy)]
pub struct NodeColorWithoutVertexRule {
    pub external_tag: i32,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct EdgeColor {
    pub pdg: isize,
}

impl std::fmt::Display for EdgeColor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.pdg)
    }
}

impl EdgeColor {
    pub(crate) fn from_particle(particle: ArcParticle) -> Self {
        Self {
            pdg: particle.0.pdg_code,
        }
    }
}

fn polyrat_to_atom(
    polyrat: &symbolica::domains::rational::Fraction<
        symbolica::poly::polynomial::PolynomialRing<
            AlgebraicExtension<
                symbolica::domains::rational::FractionField<
                    symbolica::domains::integer::IntegerRing,
                >,
            >,
            u16,
        >,
    >,
) -> Atom {
    let polyrat_numerator = polyrat.numerator();
    let num = polyrat_numerator.to_expression_with_coeff_map(|_x, a, b| {
        let a_poly = a.clone().into_poly();
        *b = Atom::num(match a_poly.coefficients.len() {
            0 => Coefficient::Complex(SymbolicaComplex::new(0.into(), 0.into())),
            1 => {
                if a_poly.exponents[0] == 1 {
                    Coefficient::Complex(SymbolicaComplex::new(
                        0.into(),
                        a_poly.coefficients[0].clone(),
                    ))
                } else {
                    Coefficient::Complex(SymbolicaComplex::new(
                        a_poly.coefficients[0].clone(),
                        0.into(),
                    ))
                }
            }
            2 => Coefficient::Complex(SymbolicaComplex::new(
                a_poly.coefficients[0].clone(),
                a_poly.coefficients[1].clone(),
            )),
            _ => unreachable!(),
        });
    });
    let polyrat_denominator = polyrat.denominator();
    let den = polyrat_denominator.to_expression_with_coeff_map(|_x, a, b| {
        let a_poly = a.clone().into_poly();

        *b = Atom::num(match a_poly.coefficients.len() {
            0 => Coefficient::Complex(SymbolicaComplex::new(0.into(), 0.into())),
            1 => {
                if a_poly.exponents[0] == 1 {
                    Coefficient::Complex(SymbolicaComplex::new(
                        0.into(),
                        a_poly.coefficients[0].clone(),
                    ))
                } else {
                    Coefficient::Complex(SymbolicaComplex::new(
                        a_poly.coefficients[0].clone(),
                        0.into(),
                    ))
                }
            }
            2 => Coefficient::Complex(SymbolicaComplex::new(
                a_poly.coefficients[0].clone(),
                a_poly.coefficients[1].clone(),
            )),
            _ => unreachable!(),
        });
    });
    num / den
}

#[derive(Clone)]
pub(crate) struct CanonizedGraphInfo {
    pub canonized_graph: SymbolicaGraph<NodeColorWithoutVertexRule, String>,
    #[allow(dead_code)]
    pub graph: SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
    pub graph_with_canonized_flow: SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
    pub gammaloop_graph: ParseGraph,
}

pub trait NodeColorFunctions: Sized + std::fmt::Display {
    fn get_external_tag(&self) -> i32;
    fn set_external_tag(&mut self, external_tag: i32);
    fn is_external(&self) -> bool {
        self.get_external_tag() > 0
    }

    fn is_incoming(&self, n_initial_states: usize) -> bool {
        self.is_external() && (self.get_external_tag() <= (n_initial_states as i32))
    }

    fn is_outgoing(&self, n_initial_states: usize) -> bool {
        self.is_external() && (self.get_external_tag() > (n_initial_states as i32))
    }

    /// Only applicable for XS
    fn pairing_tag(&self, n_initial_states: usize) -> i32 {
        let tag = self.get_external_tag();

        if tag > 0 {
            if tag > (n_initial_states as i32) {
                tag - (n_initial_states as i32)
            } else {
                tag
            }
        } else {
            tag
        }
    }

    fn coupling_orders(&self, _model: &Model) -> AHashMap<SmartString<LazyCompact>, usize> {
        AHashMap::default()
    }

    fn get_sign(&self, n_initial_states: usize) -> SignOrZero {
        if self.get_external_tag() > 0 {
            if self.get_external_tag() <= n_initial_states as i32 {
                SignOrZero::Plus
            } else {
                SignOrZero::Minus
            }
        } else {
            SignOrZero::Zero
        }
    }

    #[allow(clippy::type_complexity)]
    fn passes_amplitude_filter(
        _model: &Model,
        _amplitude_subgraph: &SuBitGraph,
        _graph: &HedgeGraph<ArcParticle, Self>,
        _amp_couplings: Option<
            &std::collections::HashMap<String, (usize, Option<usize>), ahash::RandomState>,
        >,
    ) -> bool;
}

impl NodeColorFunctions for NodeColorWithVertexRule {
    fn get_external_tag(&self) -> i32 {
        self.external_tag
    }

    fn set_external_tag(&mut self, external_tag: i32) {
        self.external_tag = external_tag;
    }

    fn coupling_orders(&self, model: &Model) -> AHashMap<SmartString<LazyCompact>, usize> {
        // info!("looking at :{}", self.vertex_rule.name);
        let mut coupling_orders = AHashMap::default();
        let vr = self.vertex_rule.clone();
        if vr.0.name == "external" {
            return coupling_orders;
        }
        for (k, v) in vr.0.coupling_orders(model) {
            *coupling_orders.entry(k).or_insert(0) += v;
        }
        coupling_orders
    }

    fn passes_amplitude_filter(
        model: &Model,
        amplitude_subgraph: &SuBitGraph,
        graph: &HedgeGraph<ArcParticle, Self>,
        amp_couplings: Option<
            &std::collections::HashMap<String, (usize, Option<usize>), ahash::RandomState>,
        >,
    ) -> bool {
        // info!(
        //     "//looking at \n{}\n",
        //     graph.dot_impl(
        //         amplitude_subgraph,
        //         "",
        //         &|a| Some(format!("label=\"{}\"", a.name)),
        //         &|b| { None }
        //     )
        // );
        if let Some(amp_couplings) = amp_couplings {
            let mut coupling_orders = AHashMap::default();
            for (_, _, s) in graph.iter_nodes_of(amplitude_subgraph) {
                // println!("node {}:{}", s.vertex_rule.name, s.get_external_tag());
                if !s.is_external() {
                    for (k, v) in s.coupling_orders(model) {
                        *coupling_orders.entry(k).or_insert(0) += v;
                    }
                }
            }

            // info!("Coupling orders: {:?}", coupling_orders);
            // info!("Amplitude couplings: {:?}", amp_couplings);

            // if ans {
            //     info!("Passes amplitude filter");
            // }
            amp_couplings.iter().all(|(k, (lower_bound, upper_bound))| {
                coupling_orders
                    .get(&SmartString::from(k))
                    .map_or(*lower_bound == 0, |o| {
                        lower_bound <= o && upper_bound.map(|a| *o <= a).unwrap_or(true)
                    })
            })
        } else {
            true
        }
    }
}

impl NodeColorFunctions for NodeColorWithoutVertexRule {
    fn get_external_tag(&self) -> i32 {
        self.external_tag
    }
    fn set_external_tag(&mut self, external_tag: i32) {
        self.external_tag = external_tag;
    }

    fn passes_amplitude_filter(
        _model: &Model,
        _amplitude_subgraph: &SuBitGraph,
        _graph: &HedgeGraph<ArcParticle, Self>,
        _amp_couplings: Option<
            &std::collections::HashMap<String, (usize, Option<usize>), ahash::RandomState>,
        >,
    ) -> bool {
        panic!("Cannot apply amplitude filters without vertex information.");
    }
}

impl std::fmt::Display for NodeColorWithoutVertexRule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.external_tag)
    }
}

impl std::fmt::Display for NodeColorWithVertexRule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "({}|{})",
            self.external_tag,
            self.vertex_rule.0.name.clone()
        )
    }
}

#[allow(clippy::type_complexity)]
pub(crate) fn group_isomorphic_graphs<
    NodeColor: Clone + PartialEq + Eq + PartialOrd + Ord + std::hash::Hash,
>(
    graphs: &[SymbolicaGraph<NodeColor, EdgeColor>],
) -> Vec<(SymbolicaGraph<NodeColor, EdgeColor>, usize)> {
    if graphs.len() == 1 {
        return vec![(graphs[0].clone(), 1)];
    }
    let mut iso_buckets: HashMap<SymbolicaGraph<NodeColor, EdgeColor>, usize> = HashMap::default();
    for g in graphs.iter() {
        let canonized_g = g.canonize();
        *iso_buckets.entry(canonized_g.graph).or_insert_with(|| 1) += 1;
    }

    iso_buckets
        .iter()
        .map(|(g, count)| (g.clone(), *count))
        .collect()
}

pub(crate) fn contains_particles(
    graph: &SymbolicaGraph<NodeColorWithoutVertexRule, EdgeColor>,
    particles: &[isize],
) -> bool {
    let mut particles_stack = particles.to_vec();
    for edge in graph.edges().iter() {
        if let Some(pos) = particles_stack.iter().position(|&x| x == edge.data.pdg) {
            particles_stack.swap_remove(pos);
        }
        if particles_stack.is_empty() {
            return true;
        }
    }
    false
}

pub(crate) fn find_edge_position(
    graph: &SymbolicaGraph<NodeColorWithoutVertexRule, EdgeColor>,
    edge_vertices: (usize, usize),
    oriented: bool,
) -> Option<usize> {
    for (i_e, edge) in graph.edges().iter().enumerate() {
        if edge.vertices == edge_vertices
            || (!oriented && edge.vertices == (edge_vertices.1, edge_vertices.0))
        {
            return Some(i_e);
        }
    }
    None
}

pub(crate) fn follow_chain(
    current_node: usize,
    vetos: &mut Vec<bool>,
    adj_map: &HashMap<usize, Vec<(usize, usize)>>,
    one_step_only: bool,
) -> Result<(Option<usize>, usize), FeynGenError> {
    if let Some(connected_edges) = adj_map.get(&current_node) {
        let targets = connected_edges
            .iter()
            .filter_map(|(i_e, next_node)| {
                if !vetos[*i_e] {
                    Some((i_e, next_node))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        if targets.is_empty() {
            Ok((None, current_node))
        } else if targets.len() == 1 {
            let (next_edge, next_node) = targets.first().unwrap();
            vetos[**next_edge] = true;
            if one_step_only {
                Ok((Some(**next_edge), **next_node))
            } else {
                follow_chain(**next_node, vetos, adj_map, one_step_only)
            }
        } else {
            Ok((None, *targets.first().unwrap().1))
            // return Err(FeynGenError::GenericError(
            //     "GammaLoop does not support four-fermion vertices yet".to_string(),
            // ));
        }
    } else {
        Ok((None, current_node))
    }
}

#[allow(clippy::type_complexity)]
pub(crate) fn assign_node_colors(
    model: &Model,
    graph: &SymbolicaGraph<NodeColorWithoutVertexRule, EdgeColor>,
    node_colors: &HashMap<
        Vec<(Option<bool>, SmartString<LazyCompact>)>,
        Vec<SmartString<LazyCompact>>,
    >,
) -> Result<Vec<(SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>, usize)>, FeynGenError> {
    // println!("graph = {}", graph.to_dot());
    let mut colored_nodes: Vec<Vec<NodeColorWithVertexRule>> = vec![vec![]];
    let edges = graph.edges();
    for (i_n, node) in graph.nodes().iter().enumerate() {
        let mut node_edges = vec![];
        for e in node.edges.iter() {
            let orientation: Option<bool> = if edges[*e].directed {
                Some(edges[*e].vertices.0 == i_n)
            } else {
                None
            };
            let particle = model.get_particle_from_pdg(edges[*e].data.pdg);
            node_edges.push((orientation, particle.0.name.clone()));
            // self-loop edge must be counted twice
            if edges[*e].vertices.0 == edges[*e].vertices.1 {
                node_edges.push((orientation.map(|o| !o), particle.0.name.clone()));
            }
        }
        node_edges.sort();
        let dummy_external_vertex_rule = ArcVertexRule(Arc::new(VertexRule {
            name: "external".into(),
            couplings: vec![],
            lorentz_structures: vec![],
            particles: vec![],
            color_structures: ColorStructure {
                color_structure: vec![],
            },
        }));
        let colors = if node_edges.len() == 1 {
            &vec![dummy_external_vertex_rule]
        } else if let Some(cs) = node_colors.get(&node_edges) {
            &cs.iter()
                .map(|c| model.get_vertex_rule(c))
                .collect::<Vec<_>>()
        } else {
            return Err(FeynGenError::GenericError(format!(
                "Could not find node colors for node edges {:?}",
                node_edges
            )));
        };
        let mut new_colored_nodes: Vec<Vec<NodeColorWithVertexRule>> = vec![];
        for current_colors in colored_nodes.iter_mut() {
            for color in colors {
                current_colors.push(NodeColorWithVertexRule {
                    external_tag: node.data.external_tag,
                    vertex_rule: color.clone(),
                });
                new_colored_nodes.push(current_colors.clone());
            }
        }

        colored_nodes = new_colored_nodes;
    }

    let mut colored_graphs = vec![];
    for color in colored_nodes.clone() {
        let mut g: SymbolicaGraph<_, _> = SymbolicaGraph::new();
        for colored_node in color {
            g.add_node(colored_node);
        }
        for edge in graph.edges() {
            _ = g.add_edge(edge.vertices.0, edge.vertices.1, edge.directed, edge.data);
        }
        colored_graphs.push(g);
    }
    Ok(group_isomorphic_graphs(&colored_graphs))
}

fn group_isomorphic_graphs_after_node_color_change<NC>(
    graphs: &HashMap<SymbolicaGraph<NC, EdgeColor>, Atom>,
    node_colors_to_change: &HashMap<i32, i32>,
    pool: &rayon::ThreadPool,
    progress_bar_style: &ProgressStyle,
) -> HashMap<SymbolicaGraph<NC, EdgeColor>, Atom>
where
    NC: NodeColorFunctions + Clone + PartialOrd + Ord + Eq + std::hash::Hash + Send + Sync,
{
    #[allow(clippy::type_complexity)]
    let iso_buckets: Arc<
        Mutex<
            HashMap<
                SymbolicaGraph<NC, EdgeColor>,
                (usize, (SymbolicaGraph<NC, EdgeColor>, Vec<usize>, Atom)),
            >,
        >,
    > = Arc::new(Mutex::new(HashMap::default()));

    let iso_buckets_clone = iso_buckets.clone();
    let bar = ProgressBar::new(graphs.len() as u64);
    bar.set_style(progress_bar_style.clone());
    bar.set_message("Grouping isomorphic topologies including external leg symmetrization...");
    pool.install(|| {
        graphs
            .par_iter()
            .progress_with(bar.clone())
            .for_each(|(g, symmetry_factor)| {
                let mut g_node_color_modified = g.clone();
                let mut modifications: Vec<(usize, i32)> = vec![];
                for (i_n, node) in g.nodes().iter().enumerate() {
                    for (src_node_color, trgt_node_color) in node_colors_to_change {
                        if node.data.get_external_tag() == *src_node_color {
                            modifications.push((i_n, *trgt_node_color));
                        }
                    }
                }
                for (i_n, new_color) in modifications {
                    let mut node_data = g.nodes()[i_n].data.clone();
                    node_data.set_external_tag(new_color);
                    g_node_color_modified.set_node_data(i_n, node_data);
                }
                let canonized_g = g_node_color_modified.canonize();
                let ext_ordering = g
                    .nodes()
                    .iter()
                    .filter_map(|n| {
                        if n.data.get_external_tag() > 0 {
                            Some(*n.edges.first().unwrap())
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>();

                {
                    let mut iso_buckets_guard = iso_buckets_clone.lock().unwrap();

                    if let Some((count, (other_g, external_ordering, sym_fact))) =
                        iso_buckets_guard.get_mut(&canonized_g.graph)
                    {
                        // Make sure to pick a representative with a canonical external ordering so that additional flavour grouping afterwards is possible
                        if ext_ordering < *external_ordering {
                            *external_ordering = ext_ordering;
                            *other_g = g.clone();
                        }
                        *count += 1;

                        assert!(sym_fact == symmetry_factor);
                    } else {
                        iso_buckets_guard.insert(
                            canonized_g.graph,
                            (1, (g.clone(), ext_ordering, symmetry_factor.clone())),
                        );
                    }
                }
            });
    });
    bar.finish_and_clear();
    let iso_buckets_guard = iso_buckets.lock().unwrap();
    iso_buckets_guard
        .iter()
        .map(
            |(_canonized_g, (count, (g, _external_ordering, symmetry_factor)))| {
                let new_symmetry_factor = if *count != 1 {
                    symmetry_factor
                        * function!(
                            symbol!("NumeratorIndependentSymmetryGrouping"),
                            Atom::num(*count as i64)
                        )
                } else {
                    symmetry_factor.clone()
                };
                (g.clone(), new_symmetry_factor)
            },
        )
        .collect()
}

#[instrument(skip_all)]
pub(crate) fn veto_special_topologies(
    model: &Model,
    graph: &SymbolicaGraph<NodeColorWithoutVertexRule, EdgeColor>,
    veto_self_energy: Option<&SelfEnergyFilterOptions>,
    veto_tadpole: Option<&TadpolesFilterOptions>,
    veto_snails: Option<&SnailFilterOptions>,
    factorized_loop_topologies_count_range: Option<&(usize, usize)>,
) -> bool {
    let graph_nodes: &[symbolica::graph::Node<NodeColorWithoutVertexRule>] = graph.nodes();
    let graph_edges: &[symbolica::graph::Edge<EdgeColor>] = graph.edges();

    let max_external = graph_nodes
        .iter()
        .filter(|n| n.data.external_tag > 0)
        .map(|n| n.data.external_tag)
        .max()
        .unwrap_or(0) as usize;

    let mut external_partices: Vec<ArcParticle> = vec![model.particles[0].clone(); max_external];
    for e in graph_edges {
        if graph_nodes[e.vertices.0].data.external_tag != 0 {
            external_partices[(graph_nodes[e.vertices.0].data.external_tag - 1) as usize] =
                model.get_particle_from_pdg(e.data.pdg);
        } else if graph_nodes[e.vertices.1].data.external_tag != 0 {
            external_partices[(graph_nodes[e.vertices.1].data.external_tag - 1) as usize] =
                model.get_particle_from_pdg(e.data.pdg);
        }
    }

    // Test vetoing of from all external spanning tree root positions to test that there are no issues from spanning tree directions
    // TODO rewrite and improve the vetoing logic of special topologies
    (0..=((max_external as isize) - 1).max(0)).all(|shift| {
        // (0..=0).all(|shift: usize| {
        let spanning_tree_root_node_position = graph_nodes
            .iter()
            .position(|n| n.data.external_tag == ((max_external - shift as usize) as i32))
            .unwrap();
        // println!(
        //     "Spanning tree root position: external_tag={},node_position={}",
        //     max_external - shift as usize,
        //     spanning_tree_root_node_position
        // );

        veto_special_topologies_with_spanning_tree_root(
            model,
            graph,
            veto_self_energy,
            veto_tadpole,
            veto_snails,
            factorized_loop_topologies_count_range,
            &external_partices,
            spanning_tree_root_node_position,
        )
    })
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn veto_special_topologies_with_spanning_tree_root(
    model: &Model,
    graph: &SymbolicaGraph<NodeColorWithoutVertexRule, EdgeColor>,
    veto_self_energy: Option<&SelfEnergyFilterOptions>,
    veto_tadpole: Option<&TadpolesFilterOptions>,
    veto_snails: Option<&SnailFilterOptions>,
    factorized_loop_topologies_count_range: Option<&(usize, usize)>,
    external_particles: &[ArcParticle],
    spanning_tree_root: usize,
) -> bool {
    let max_external = external_particles.len();

    if graph.nodes().iter().any(|n| n.data.external_tag < 0) {
        panic!(
            "External tag must be positive, but found negative as obtained when performing external state symmetrization"
        );
    }
    const DEBUG_VETO: bool = false;
    if DEBUG_VETO {
        debug!(
            "\n\n>> Vetoing special topologies for the following {}-loop graph:\n{}",
            graph.num_loops(),
            graph.to_dot()
        );
        debug!(
            "Graph edges raw representation:\n{}",
            graph
                .edges()
                .iter()
                .enumerate()
                .map(|(i_e, e)| format!(
                    "[ #{} vertices: {:?}, e.directed: {}, e.data: {:?} ]",
                    i_e, e.vertices, e.directed, e.data
                ))
                .join("\n")
        );
        debug!(
            "Graph nodes raw representation:\n{}",
            graph
                .nodes()
                .iter()
                .enumerate()
                .map(|(i_n, n)| format!("[ #{} edges: {:?}, n.data: {:?} ]", i_n, n.edges, n.data))
                .join("\n")
        );
    }
    if veto_self_energy.is_none()
        && veto_tadpole.is_none()
        && veto_snails.is_none()
        && factorized_loop_topologies_count_range.is_none()
    {
        return false;
    }
    let graph_edges: &[symbolica::graph::Edge<EdgeColor>] = graph.edges();
    let spanning_tree_node_external_tag =
        graph.nodes()[spanning_tree_root].data.external_tag as usize;
    let mut spanning_tree = graph.get_spanning_tree(spanning_tree_root);
    spanning_tree.chain_decomposition();

    if DEBUG_VETO {
        debug!(
            "Spanning tree: order={:?}\n{}",
            spanning_tree.order,
            spanning_tree
                .nodes
                .iter()
                .enumerate()
                .map(|(i_n, n)| format!(
                    "[ #{} position: {:?}, parent: {:?}, back_edges: {:?}, external: {}, chain_id: {:?} ]",
                    i_n, n.position, n.parent, n.back_edges, n.external, n.chain_id
                ))
                .join("\n")
        );
    }
    let mut node_children: Vec<Vec<usize>> = vec![vec![]; spanning_tree.nodes.len()];
    for (i_n, node) in spanning_tree.nodes.iter().enumerate() {
        node_children[node.parent].push(i_n);
    }
    if DEBUG_VETO {
        debug!("node_children={:?}", node_children);
    }

    let mut external_momenta_routing: Vec<Vec<usize>> = vec![vec![]; spanning_tree.nodes.len()];
    for (i_n, node) in graph.nodes().iter().enumerate() {
        if (node.edges.len() != 1)
            || node.data.external_tag == (spanning_tree_node_external_tag as i32)
        {
            continue;
        }
        let external_index = node.data.external_tag as usize;
        external_momenta_routing[i_n].push(external_index);
        let mut next_node = spanning_tree.nodes[i_n].parent;
        // println!("max_external={},max_external_node_position={},sink_node_position_in_spanning_tree={}", max_external, max_external_node_position, sink_node_position_in_spanning_tree);
        // println!("external_momenta_routing={:?}", external_momenta_routing);
        // println!("spanning_tree={:?}", spanning_tree);
        // println!("Starting from node #{} = {:?}", i_n, node);
        while next_node != spanning_tree_root {
            external_momenta_routing[next_node].push(external_index);
            next_node = spanning_tree.nodes[next_node].parent;
        }
    }
    for route in external_momenta_routing.iter_mut() {
        if max_external > 0 && route.len() == max_external - 1 {
            *route = vec![spanning_tree_node_external_tag];
        }
    }
    if DEBUG_VETO {
        debug!("external_momenta_routing={:?}", external_momenta_routing);
    }

    // See https://arxiv.org/pdf/1209.0700 for information on the logic of this algorithm

    // Tuple format: (external_leg_id, back_edge_start_node_position, back_edge_position_in_list, chain_id)
    let mut self_energy_attachments: HashSet<(usize, usize, usize, usize)> = HashSet::default();
    // Tuple format: (back_edge_start_node_position, back_edge_position_in_list, chain_id)
    let mut vacuum_attachments: HashSet<(usize, usize, usize)> = HashSet::default();
    // Tuple format: (back_edge_start_node_position, back_edge_position_in_list, chain_id)
    let mut self_loops: HashSet<(usize, usize, usize)> = HashSet::default();
    let mut n_factorizable_loops = 0;

    let mut visited_nodes = vec![false; spanning_tree.nodes.len()];
    for &i_n in &spanning_tree.order {
        let node = &spanning_tree.nodes[i_n];
        for (i_back_edge, &back_edge) in node.back_edges.iter().enumerate() {
            let i_chain = i_n;
            if back_edge.target == i_n {
                n_factorizable_loops += 1;

                self_loops.insert((i_n, i_back_edge, i_chain));
                continue;
            }
            let mut self_energy_external_leg_id: Option<usize> = None;
            let mut curr_chain_node = back_edge.target;
            let mut is_valid_chain = true;
            'follow_chain: loop {
                if curr_chain_node == i_n {
                    // This conditional is wrong!
                    // if i_back_edge > 0 {
                    n_factorizable_loops += 1;
                    // }
                    break 'follow_chain;
                }
                let moms = &external_momenta_routing[curr_chain_node];
                if moms.len() == 1 {
                    if let Some(se_leg) = self_energy_external_leg_id
                        && se_leg != moms[0]
                    {
                        is_valid_chain = false;
                        break 'follow_chain;
                    }
                    self_energy_external_leg_id = Some(moms[0]);
                    for child in node_children[curr_chain_node].iter() {
                        if (!external_momenta_routing[*child].is_empty())
                            && external_momenta_routing[*child][0] != moms[0]
                        {
                            is_valid_chain = false;
                            break 'follow_chain;
                        }
                    }
                } else if !moms.is_empty() {
                    is_valid_chain = false;
                    break 'follow_chain;
                }
                // if let Some(chain_id) = spanning_tree.nodes[curr_chain_node].chain_id {
                //     if chain_id != i_chain {
                //         is_valid_chain = false;
                //         break 'follow_chain;
                //     }
                // } else {
                //     is_valid_chain = false;
                //     break 'follow_chain;
                // }
                if spanning_tree.nodes[curr_chain_node].chain_id.is_none() {
                    is_valid_chain = false;
                    break 'follow_chain;
                }

                if visited_nodes[curr_chain_node] {
                    is_valid_chain = false;
                    break 'follow_chain;
                } else {
                    visited_nodes[curr_chain_node] = true;
                    curr_chain_node = spanning_tree.nodes[curr_chain_node].parent;
                }
            }

            if is_valid_chain {
                if let Some(leg_id) = self_energy_external_leg_id {
                    // Make sure the attachment point of the self-energy does not receive any other external momenta
                    if external_momenta_routing[i_n].len() == 1
                        && external_momenta_routing[i_n][0] == leg_id
                    {
                        // Also For 1 -> 1 processes, we must also verify that it is not the whole graph
                        // Also For 1 -> 1 processes, we must also verify that it is not the whole graph
                        if max_external > 2
                            || spanning_tree
                                .nodes
                                .iter()
                                .filter(|n| (!n.external) && n.chain_id.is_none())
                                .count()
                                > 1
                        {
                            self_energy_attachments.insert((leg_id, i_n, i_back_edge, i_chain));
                        }
                    }
                } else {
                    vacuum_attachments.insert((i_n, i_back_edge, i_chain));
                }
            }
        }
    }
    if DEBUG_VETO {
        debug!("self_energy_attachments={:?}", self_energy_attachments);
        debug!("vacuum_attachments={:?}", vacuum_attachments);
        debug!("self_loops={:?}", self_loops);
        debug!("n_factorizable_loops={:?}", n_factorizable_loops);
    }

    if let Some((min_n_fact_loops, max_n_fact_loops)) =
        factorized_loop_topologies_count_range.as_ref()
        && (n_factorizable_loops < *min_n_fact_loops || n_factorizable_loops > *max_n_fact_loops)
    {
        if DEBUG_VETO {
            debug!(
                "Vetoing graph due to having a number of factorizable loops ({}) outside the range [{}, {}]",
                n_factorizable_loops, min_n_fact_loops, max_n_fact_loops
            );
        }
        return true;
    }

    let mut tree_bridge_node_indices: HashSet<usize> = HashSet::default();
    for (i_n, node) in spanning_tree.nodes.iter().enumerate() {
        if node.chain_id.is_none()
            && !node.external
            && !external_momenta_routing[i_n].is_empty()
            && !spanning_tree.nodes[node.parent]
                .back_edges
                .iter()
                .any(|&end| i_n == end.target)
        {
            tree_bridge_node_indices.insert(i_n);
        }
    }
    if DEBUG_VETO {
        debug!("bridge_node_positions={:?}", tree_bridge_node_indices);
    }

    // For self-energies we must confirm that they are self-energies by checking if the back edge start node is a bridge
    for (leg_id, back_edge_start_node_index, back_edge_position_in_list, _chain_id) in
        self_energy_attachments.iter()
    {
        if tree_bridge_node_indices.contains(back_edge_start_node_index)
            && let Some(veto_self_energy_options) = veto_self_energy
        {
            if DEBUG_VETO {
                debug!(
                    "Vetoing self-energy for leg_id={}, back_edge_start_node_index={}, back_edge_position_in_list={}, chain_id={}, with options:\n{:?}",
                    leg_id,
                    back_edge_start_node_index,
                    back_edge_position_in_list,
                    _chain_id,
                    veto_self_energy_options
                );
            }
            if veto_self_energy_options.veto_only_scaleless_self_energy {
                panic!("Option to only remove scaleless self-energies is not implemented yet");
            } else {
                #[allow(clippy::unnecessary_unwrap)]
                if external_particles[leg_id - 1].0.is_massive() {
                    if veto_self_energy_options.veto_self_energy_of_massive_lines {
                        return true;
                    }
                } else {
                    #[allow(clippy::unnecessary_unwrap)]
                    if veto_self_energy_options.veto_self_energy_of_massless_lines {
                        return true;
                    }
                }
            }
        }
    }
    // For vaccuum attachments, we must differentiate if there are snails (start node is a tree bridge) or tadpoles (start node *is not* a tree bridge)
    for (back_edge_start_node_index, back_edge_position_in_list, chain_id) in
        vacuum_attachments.iter().chain(self_loops.iter())
    {
        let attachment_particle_is_massive = if max_external > 0 {
            let mut first_tree_attachment_node_index = *back_edge_start_node_index;
            while external_momenta_routing[first_tree_attachment_node_index].is_empty()
                && spanning_tree.nodes[first_tree_attachment_node_index]
                    .chain_id
                    .is_none()
            {
                first_tree_attachment_node_index =
                    spanning_tree.nodes[first_tree_attachment_node_index].parent;
            }
            let attachment_edge = &graph_edges[find_edge_position(
                graph,
                (
                    first_tree_attachment_node_index,
                    spanning_tree.nodes[first_tree_attachment_node_index].parent,
                ),
                false,
            )
            .unwrap_or_else(|| {
                panic!("Could not find edge between bridge node parent and grandparent")
            })];

            model
                .get_particle_from_pdg(attachment_edge.data.pdg)
                .0
                .is_massive()
        } else {
            // Always consider the attachment particle as massive for vaccuum graphs as it does not matter in that case
            // given that there is no support for differentiating massive and massless attachments in that case.
            true
        };

        if !tree_bridge_node_indices.contains(back_edge_start_node_index)
            && spanning_tree.nodes[*back_edge_start_node_index]
                .chain_id
                .is_none()
        {
            // Tadpole
            if let Some(veto_tadpole_options) = veto_tadpole {
                if DEBUG_VETO {
                    debug!(
                        "Vetoing tadpole for back_edge_start_node_index={}, back_edge_position_in_list={}, chain_id={}, with options:\n{:?}",
                        back_edge_start_node_index,
                        back_edge_position_in_list,
                        chain_id,
                        veto_tadpole_options
                    );
                }
                if veto_tadpole_options.veto_only_scaleless_tadpoles {
                    panic!("Option to only remove scaleless self-energies is not implemented yet");
                } else {
                    #[allow(clippy::unnecessary_unwrap)]
                    if attachment_particle_is_massive {
                        if veto_tadpole_options.veto_tadpoles_attached_to_massive_lines {
                            return true;
                        }
                    } else {
                        #[allow(clippy::unnecessary_unwrap)]
                        if veto_tadpole_options.veto_tadpoles_attached_to_massless_lines {
                            return true;
                        }
                    }
                }
            }
        } else {
            #[allow(clippy::unnecessary_unwrap)]
            // Snail
            if let Some(veto_snails_options) = veto_snails {
                if DEBUG_VETO {
                    debug!(
                        "Vetoing snail for back_edge_start_node_index={}, back_edge_position_in_list={}, chain_id={}, with options:\n{:?}",
                        back_edge_start_node_index,
                        back_edge_position_in_list,
                        chain_id,
                        veto_snails_options
                    );
                }
                #[allow(clippy::unnecessary_unwrap)]
                if veto_snails_options.veto_only_scaleless_snails {
                    panic!("Option to only remove scaleless self-energies is not implemented yet");
                } else {
                    #[allow(clippy::unnecessary_unwrap)]
                    if attachment_particle_is_massive {
                        if veto_snails_options.veto_snails_attached_to_massive_lines {
                            return true;
                        }
                    } else {
                        #[allow(clippy::unnecessary_unwrap)]
                        if veto_snails_options.veto_snails_attached_to_massless_lines {
                            return true;
                        }
                    }
                }
            }
        }
    }

    if DEBUG_VETO {
        debug!(">> No special topology veto applied to this graph");
    }

    false
}

pub fn evaluate_overall_factor(factor: AtomView) -> Atom {
    let mut res = factor.to_owned();
    for header in [
        "AutG",
        "CouplingsMultiplicity",
        "InternalFermionLoopSign",
        "ExternalFermionOrderingSign",
        "AntiFermionSpinSumSign",
        "NumeratorIndependentSymmetryGrouping",
    ] {
        res = res
            .replace(function!(symbol!(header), Atom::var(symbol!("x_"))).to_pattern())
            .with(Atom::var(symbol!("x_")).to_pattern());
    }
    res = res
        .replace(
            function!(
                symbol!("NumeratorDependentGrouping"),
                Atom::var(symbol!("GraphId_")),
                Atom::var(symbol!("ratio_")),
                Atom::var(symbol!("GraphSymmetryFactor_"))
            )
            .to_pattern(),
        )
        .with(
            (Atom::var(symbol!("ratio_")) * Atom::var(symbol!("GraphSymmetryFactor_")))
                .to_pattern(),
        );
    res.expand()
}

pub fn evaluate_sign_origin(factor: AtomView) -> Atom {
    let mut res = factor.to_owned();
    for header in [
        "AutG",
        "CouplingsMultiplicity",
        "InternalFermionLoopSign",
        "ExternalFermionOrderingSign",
        "AntiFermionSpinSumSign",
        "NumeratorIndependentSymmetryGrouping",
    ] {
        res = res
            .replace(function!(symbol!(header), Atom::var(symbol!("x_"))).to_pattern())
            .with(Atom::var(symbol!("x_")).to_pattern());
    }
    res = res
        .replace(function!(
            symbol!("NumeratorDependentGrouping"),
            Atom::var(symbol!("GraphId_")),
            Atom::var(symbol!("ratio_")),
            Atom::var(symbol!("GraphSymmetryFactor_"))
        ))
        .with(function!(
            symbol!("Group"),
            Atom::var(symbol!("GraphId_")),
            Atom::var(symbol!("ratio_")),
            Atom::var(symbol!("GraphSymmetryFactor_"))
        ));
    res.expand()
}

impl ProcessDefinition {
    pub fn sample_lib(
        &self,
        sample_iterator: &mut PrimeIteratorU64,
        add_model_params: bool,
        symmetric_polarizations: bool,
        model: &Model,
    ) -> (
        Vec<Replacement>,
        TensorLibrary<ParamTensor<ExplicitKey<Aind>>, Aind>,
    ) {
        // let ff = Zp64::new(PrimeIteratorU64::new(100000000).next().unwrap());

        let mut weyl = TensorLibrary::new();
        weyl.update_ids();

        let gamma_key = PermutedStructure::identity(ParamTensor::composite(DataTensor::Sparse(
            gamma_data_weyl(AGS.gamma_strct::<Aind>(4), Atom::num(1), Atom::num(0))
                .map_data(|a| a.re + Atom::i() * a.im),
        )));
        // println!("permutation{}", gamma_key.rep_permutation);
        weyl.insert_explicit(gamma_key);

        let gamma5_key = PermutedStructure::identity(ParamTensor::composite(DataTensor::Sparse(
            gamma5_weyl_data(AGS.gamma5_strct::<Aind>(4), Atom::num(1), Atom::num(0))
                .map_data(|a| a.re + Atom::i() * a.im),
        )));
        weyl.insert_explicit(gamma5_key);

        let projm_key = PermutedStructure::identity(ParamTensor::composite(DataTensor::Sparse(
            proj_m_data_weyl(AGS.projm_strct::<Aind>(4), Atom::num(1), Atom::num(0))
                .map_data(|a| a.re + Atom::i() * a.im),
        )));
        weyl.insert_explicit(projm_key);

        let projp_key = PermutedStructure::identity(ParamTensor::composite(DataTensor::Sparse(
            proj_p_data_weyl(AGS.projp_strct::<Aind>(4), Atom::num(1), Atom::num(0))
                .map_data(|a| a.re + Atom::i() * a.im),
        )));
        weyl.insert_explicit(projp_key);
        let mut lib = weyl;

        for i in 0..self.loop_count_range.1 {
            let key = ExplicitKey::from_iter(
                [Minkowski {}.new_rep(4)],
                GS.loop_mom,
                Some(vec![Atom::num(i)]),
            );

            //debug!("lib_loop:{}", key.clone().permute_with_metric());
            let key = ParamTensor::from_dense(
                key.structure,
                (0..4)
                    .map(|_| Atom::prime_generate_rat(sample_iterator))
                    .collect(),
            )
            .unwrap();

            lib.insert_explicit(PermutedStructure::identity(key));
        }

        let mut pol_vals = vec![];
        for (i, pdg) in self.initial_pdgs.iter().enumerate() {
            let additional_args = Some(vec![Atom::num(i)]);
            let key = ExplicitKey::from_iter(
                [Minkowski {}.new_rep(4)],
                GS.external_mom,
                additional_args.clone(),
            );

            //debug!("lib_ext:{}", key.clone().permute_with_metric());

            let key = ParamTensor::from_dense(
                key.structure,
                (0..4)
                    .map(|_| Atom::prime_generate_rat(sample_iterator))
                    .collect(),
            )
            .unwrap();

            lib.insert_explicit(PermutedStructure::identity(key));

            let p = model.get_particle_from_pdg(*pdg as isize);

            let structure = p.spin_reps();
            let global_name = EdgeData::new(p, Orientation::Default).pol_symbol(Flow::Sink);

            let len = structure.size().unwrap();

            let key = PermutedStructure::identity(ExplicitKey {
                structure,
                global_name,
                additional_args,
            });

            pol_vals.push(
                (0..len)
                    .map(|_| Atom::prime_generate_rat_complex(sample_iterator))
                    .collect::<Vec<_>>(),
            );

            //debug!("lib_pol:{}", key.clone().permute_with_metric());
            let key =
                ParamTensor::from_dense(key.structure, pol_vals.last().unwrap().clone()).unwrap();

            lib.insert_explicit(PermutedStructure::identity(key));
        }
        match self.generation_type {
            GenerationType::Amplitude => {
                let ext_shift = self.initial_pdgs.len();
                for (i, pdg) in self.final_pdgs_lists[0].iter().enumerate() {
                    let additional_args = Some(vec![Atom::num(i + ext_shift)]);
                    let key = ExplicitKey::from_iter(
                        [Minkowski {}.new_rep(4)],
                        GS.external_mom,
                        additional_args.clone(),
                    );

                    let key = ParamTensor::from_dense(
                        key.structure,
                        (0..4)
                            .map(|_| Atom::prime_generate_rat(sample_iterator))
                            .collect(),
                    )
                    .unwrap();

                    lib.insert_explicit(PermutedStructure::identity(key));

                    let p = model.get_particle_from_pdg(*pdg as isize);

                    let structure = p.spin_reps();
                    let global_name = EdgeData::new(p, Orientation::Default).pol_symbol(Flow::Sink);

                    let len = structure.size().unwrap();

                    let key = ExplicitKey {
                        structure,
                        global_name,
                        additional_args,
                    };

                    let key = ParamTensor::from_dense(
                        key,
                        (0..len)
                            .map(|_| Atom::prime_generate_rat_complex(sample_iterator))
                            .collect(),
                    )
                    .unwrap();

                    lib.insert_explicit(PermutedStructure::identity(key));
                }
            }
            GenerationType::CrossSection => {
                for (i, pdg) in self.initial_pdgs.iter().enumerate() {
                    let additional_args = Some(vec![Atom::num(i)]);

                    let p = model.get_particle_from_pdg(*pdg as isize);

                    let structure = p.spin_reps();
                    let global_name =
                        EdgeData::new(p, Orientation::Default).pol_symbol(Flow::Source);

                    // let len = structure.size().unwrap();

                    let key = PermutedStructure::identity(ExplicitKey {
                        structure,
                        global_name,
                        additional_args,
                    });

                    //debug!("lib_pol:{}", key.clone().permute_with_metric());
                    let key = ParamTensor::from_dense(
                        key.structure,
                        pol_vals[i]
                            .iter()
                            .map(|a| {
                                if symmetric_polarizations {
                                    a.clone()
                                } else {
                                    Atom::prime_generate_rat_complex(sample_iterator)
                                }
                            })
                            .collect(),
                    )
                    .unwrap();

                    lib.insert_explicit(PermutedStructure::identity(key));
                }
            }
        }
        let mut reps = vec![];

        if add_model_params {
            for param in model.parameters.values().filter(|p| p.value.is_some()) {
                if param.value.is_some() {
                    let name: Atom = param.name.into();
                    reps.push(Replacement::new(
                        name.to_pattern(),
                        Atom::prime_generate_rat_complex(sample_iterator),
                    ));
                }
            }
        }
        // for r in &reps {
        //     info!("Model replacements: {}", r);
        // }
        (reps, lib)
    }

    #[instrument(skip_all)]

    pub(crate) fn unresolved_cut_content(&self, model: &Model) -> (usize, AHashSet<ArcParticle>) {
        if let Some(p) = self.cross_section_filters.get_perturbative_orders() {
            let mut unresolved = AHashSet::new();
            for k in p.keys() {
                let k: SmartString<_> = k.clone().into();
                if let Some(p) = model.unresolved_particles.get(&k) {
                    unresolved = unresolved.union(p).cloned().collect();
                }
            }
            (p.values().sum(), unresolved)
        } else {
            (0, AHashSet::new())
        }
    }

    pub(crate) fn half_edge_filters<NodeColor>(
        &self,
        model: &Model,
        graph: &SymbolicaGraph<NodeColor, EdgeColor>,
        external_connections: &[(Option<usize>, Option<usize>)],
        n_unresolved: usize,
        unresolved_type: &AHashSet<ArcParticle>,
    ) -> bool
    where
        NodeColor: NodeColorFunctions + Clone,
    {
        #[allow(clippy::too_many_arguments)]
        #[allow(clippy::type_complexity)]
        fn is_valid_cut<NodeColor: NodeColorFunctions>(
            cut: &(SuBitGraph, OrientedCut, SuBitGraph),
            blob_range: &RangeInclusive<usize>,
            spectator_range: &RangeInclusive<usize>,
            model: &Model,
            n_unresolved: usize,
            unresolved_type: &AHashSet<ArcParticle>,
            particle_content_options: &[Vec<ArcParticle>],
            amp_couplings: Option<
                &std::collections::HashMap<String, (usize, Option<usize>), ahash::RandomState>,
            >,
            amp_loop_count: Option<(usize, usize)>,
            graph: &HedgeGraph<ArcParticle, NodeColor>,
        ) -> bool {
            if validate_connectivity(&cut.0, blob_range, spectator_range, graph)
                && validate_connectivity(&cut.2, blob_range, spectator_range, graph)
            {
                let cut_content: Vec<_> = cut
                    .1
                    .iter_left_hedges()
                    .map(|h| {
                        let o = graph.flow(h);
                        if matches!(o, Flow::Sink) {
                            graph[[&h]].0.as_ref().get_anti_particle(model)
                        } else {
                            graph[[&h]].clone()
                        }
                    })
                    .collect();

                // info!(
                //     "//left\n{}\n",
                //     graph.dot_impl(
                //         &cut.0,
                //         "",
                //         &|a| Some(format!("label=\"{}\"", a.name)),
                //         &|b| None
                //     )
                // // );

                // info!(
                //     "//cut\n{}\n",
                //     graph.dot_impl(
                //         &cut.1.reference,
                //         "",
                //         &|a| Some(format!("label=\"{}\"", a.name)),
                //         &|b| None
                //     )
                // );
                // for p in cut_content.iter() {
                //     info!("Particle {} in cut", p.name);
                // }
                // info!(
                //     "//right\n{}\n",
                //     graph.dot_impl(
                //         &cut.2,
                //         "",
                //         &|a| Some(format!("label=\"{}\"", a.name)),
                //         &|b| None
                //     )
                // );
                if !particle_content_options.iter().any(|particle_content| {
                    let mut cut_content_clone = cut_content.clone();
                    for p in particle_content.iter() {
                        if let Some(pos) = cut_content_clone.iter().position(|c| c == p) {
                            cut_content_clone.swap_remove(pos);
                        } else {
                            // info!("Particle {} not found in cut content", p.name);
                            return false;
                        }
                    }
                    if cut_content_clone.len() > n_unresolved {
                        // info!(
                        //     "Cut content has {} more particles than {} allowed unresolved particles",
                        //     cut_content.len() - n_unresolved,
                        //     n_unresolved
                        // );
                        return false;
                    }

                    for p in cut_content_clone.iter() {
                        if !unresolved_type.contains(p) {
                            // info!("Particle {} not found in unresolved type", p.name);
                            return false;
                        }
                    }

                    true
                }) {
                    // info!("Cut content not right");
                    return false;
                }

                if let Some((min_loop, max_loop)) = amp_loop_count {
                    let left_internal_subgraph =
                        InternalSubGraph::cleaned_filter_pessimist(cut.0.clone(), graph);

                    let right_internal_subgraph =
                        InternalSubGraph::cleaned_filter_pessimist(cut.2.clone(), graph);

                    let left_loop = graph.cyclotomatic_number(&left_internal_subgraph);
                    let right_loop = graph.cyclotomatic_number(&right_internal_subgraph);

                    let sum = left_loop + right_loop;
                    if sum < min_loop || sum > max_loop {
                        // info!(
                        //     "Loop count sum {} not within range [{}, {}]",
                        //     sum, min_loop, max_loop
                        // );
                        return false;
                    }
                }
                let left = NodeColor::passes_amplitude_filter(model, &cut.0, graph, amp_couplings);

                let right = NodeColor::passes_amplitude_filter(model, &cut.2, graph, amp_couplings);

                if !left {
                    // info!("Left node color not valid");
                    return false;
                }
                if !right {
                    // info!("Right node color not valid");
                    return false;
                }

                left && right
            } else {
                // info!("isn't s-channel");
                false
            }
        }

        fn validate_connectivity<NodeColor>(
            subgraph: &SuBitGraph,
            blob_range: &RangeInclusive<usize>,
            spectator_range: &RangeInclusive<usize>,
            graph: &HedgeGraph<ArcParticle, NodeColor>,
        ) -> bool {
            let components = graph.connected_components(subgraph);

            let mut n_blobs = 0;
            let mut n_spectators = 0;

            for component in components {
                if component.n_included() > 1 {
                    n_blobs += 1;
                } else {
                    n_spectators += 1;
                }
            }

            blob_range.contains(&n_blobs) && spectator_range.contains(&n_spectators)
        }

        let amp_couplings = self.amplitude_filters.get_coupling_orders();
        let amp_loop_count = self.amplitude_filters.get_loop_count_range();
        let blob_range = self.cross_section_filters.get_blob_range().unwrap();
        let spectator_range = self.cross_section_filters.get_spectator_range().unwrap();

        let n_particles = self.initial_pdgs.len();

        let mut he_graph = HedgeGraph::<EdgeColor, NodeColor>::from_sym(graph.clone()).map(
            |_, _, node_color| node_color,
            |_, _, _, _, d| d.map(|d| model.get_particle_from_pdg(d.pdg)),
            |_, n| n,
        );
        // info!(
        //     "Looking at\n{}",
        //     he_graph.dot_impl(
        //         &he_graph.full_filter(),
        //         "",
        //         &|a| Some(format!("label=\"{}\"", a.name)),
        //         &|b| None
        //     )
        // );
        let mut s_set = Vec::new();
        let mut t_set = Vec::new();

        for (id, _, f) in he_graph.iter_nodes() {
            match f.get_sign(n_particles) {
                SignOrZero::Plus => {
                    s_set.push(id);
                }
                SignOrZero::Minus => {
                    t_set.push(id);
                }
                _ => {}
            }
        }

        // info!("HI");
        if !s_set.is_empty() && !t_set.is_empty() {
            // info!("{}-{}", s, t);
            let cuts = he_graph.all_cuts_from_ids(&s_set, &t_set);

            // info!("found {} cuts", cuts.len());
            let particle_content_options = self
                .final_pdgs_lists
                .iter()
                .map(|pdg_list| {
                    pdg_list
                        .iter()
                        .map(|pdg| model.get_particle_from_pdg((*pdg) as isize))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let pass_cut_filter = cuts.iter().any(|c| {
                is_valid_cut(
                    c,
                    blob_range,
                    spectator_range,
                    model,
                    n_unresolved,
                    unresolved_type,
                    &particle_content_options,
                    amp_couplings,
                    amp_loop_count,
                    &he_graph,
                )
            });

            if !pass_cut_filter {
                return false;
            }

            if self.cross_section_filters.filter_cross_section_tadpoles() {
                let externals: Vec<_> = he_graph
                    .iter_nodes()
                    .filter_map(|(id, _, n)| if n.is_external() { Some(id) } else { None })
                    .collect();

                let connected_components_before = he_graph.tadpoles(&externals).len() + 1;
                for (i, f) in external_connections {
                    if let (Some(i), Some(f)) = (i, f) {
                        let i = he_graph
                            .iter_nodes()
                            .find_position(|a| a.2.get_external_tag() == *i as i32);
                        let f = he_graph
                            .iter_nodes()
                            .find_position(|a| a.2.get_external_tag() == *f as i32);

                        if let (Some((id_i, (_, _, c_i))), Some((id_f, (_, _, _)))) = (i, f) {
                            he_graph
                                .identify_nodes(&[NodeIndex(id_i), NodeIndex(id_f)], c_i.clone());
                        }
                    }
                }

                // info!("Identified nodes");

                let non_bridges = he_graph.non_bridges();

                // info!("\\he_graph:\n{}", he_graph.dot(&non_bridges));

                let connected_components = he_graph.count_connected_components(&non_bridges);
                if connected_components == connected_components_before {
                    return true;
                } else {
                    // info!(
                    //     "Vetoed: \n{}\n because it had {} connected components",
                    //     he_graph.dot_impl(
                    //         &non_bridges,
                    //         "",
                    //         &|a| Some(format!("label=\"{}\"", a.name)),
                    //         &|b| Some(format!("label=\"{}\"", b.get_external_tag()))
                    //     ),
                    //     connected_components
                    // );
                    // info!(
                    //     "OG:\n {}",
                    //     og_he.dot_impl(
                    //         &non_bridges,
                    //         "",
                    //         &|a| Some(format!("label=\"{}\"", a.name)),
                    //         &|b| Some(format!("label=\"{}\"", b.get_external_tag()))
                    //     ),
                    // );
                    // info!("{}", graph.to_dot());

                    return false;
                }
            }

            true
        } else {
            warn!("No external particles found");
            true //TODO still check the amplitude level filters in the case where there is no initial-state specified
        }
    }

    // This fast cut checker does not enumerate all cuts, but rather checks if the graph can contain a cut with the right particles
    // It also does not consider amplitude-level filters
    pub(crate) fn contains_cut_fast<NodeColor: NodeColorFunctions>(
        &self,
        model: &Model,
        graph: &SymbolicaGraph<NodeColor, EdgeColor>,
        n_unresolved: usize,
        unresolved_type: &AHashSet<ArcParticle>,
        particles: &[isize],
    ) -> bool {
        let n_initial_states = self.initial_pdgs.len();
        // Quick filter to check if the diagram can possibly contain a cut with the right particles and splitting the graph in two parts,
        // as expected for a final-state cut of a forward scattering graph.

        // TODO: replace this with s_and_t_cut from hedge, so as to have proper handling of particle vs anti-particles.
        //let he_graph = HedgeGraph::from(graph.clone());
        //he_graph.all_s_t_cuts(s, t, regions)

        // for now normalize it all to particles
        let mut cut_map: std::collections::HashMap<
            isize,
            std::vec::Vec<usize>,
            ahash::RandomState,
        > = HashMap::default();
        for &p in particles.iter().collect::<HashSet<_>>() {
            let particle = model.get_particle_from_pdg(p);
            if particle.0.is_antiparticle() {
                cut_map.insert(particle.0.get_anti_particle(model).0.pdg_code, vec![]);
            } else {
                cut_map.insert(p, vec![]);
            }
        }
        let mut unresolved_set: Vec<usize> = vec![];
        for (i_e, edge) in graph.edges().iter().enumerate() {
            // filter out external edges
            if graph.nodes()[edge.vertices.0].data.get_external_tag() != 0
                || graph.nodes()[edge.vertices.1].data.get_external_tag() != 0
            {
                continue;
            }
            let particle = model.get_particle_from_pdg(edge.data.pdg);
            let e = if particle.0.is_antiparticle() {
                cut_map.get_mut(&particle.0.get_anti_particle(model).0.pdg_code)
            } else {
                cut_map.get_mut(&particle.0.pdg_code)
            };
            if let Some(cut_entry) = e {
                cut_entry.push(i_e);
            }
            if unresolved_type.contains(&particle) {
                unresolved_set.push(i_e);
            }
        }
        // Generate unique, unordered combinations, for each possible unresolved multiplicity
        let mut unique_combinations: Vec<
            std::collections::HashSet<Vec<usize>, ahash::RandomState>,
        > = vec![];
        for unresolved_multiplicity in 0..=n_unresolved {
            let mut unique_combinations_for_this_multiplicity: std::collections::HashSet<
                Vec<usize>,
                ahash::RandomState,
            > = HashSet::default();
            let mut particles_for_this_multiplicity =
                particles.iter().map(Some).collect::<Vec<_>>();
            for _ in 0..unresolved_multiplicity {
                particles_for_this_multiplicity.push(None);
            }
            for combination in particles_for_this_multiplicity
                .iter()
                .map(|p| {
                    if let Some(particle) = p {
                        cut_map.get(particle).unwrap().iter()
                    } else {
                        unresolved_set.iter()
                    }
                })
                .multi_cartesian_product()
            {
                // Sort the combination and insert into HashSet
                let mut sorted_combination: Vec<_> = combination.into_iter().cloned().collect();
                if sorted_combination.iter().collect::<HashSet<_>>().len()
                    != sorted_combination.len()
                {
                    continue;
                }
                sorted_combination.sort_unstable();
                unique_combinations_for_this_multiplicity.insert(sorted_combination);
            }
            unique_combinations.push(unique_combinations_for_this_multiplicity);
        }

        let mut adj_list: HashMap<usize, Vec<(usize, usize)>> = HashMap::default();
        for (i_e, e) in graph.edges().iter().enumerate() {
            adj_list
                .entry(e.vertices.0)
                .or_default()
                .push((i_e, e.vertices.1));
            if e.vertices.0 != e.vertices.1 {
                adj_list
                    .entry(e.vertices.1)
                    .or_default()
                    .push((i_e, e.vertices.0));
            }
        }

        let left_initial_state_positions = (1..=n_initial_states)
            .map(|leg_id| {
                graph
                    .nodes()
                    .iter()
                    .position(|n| n.data.get_external_tag() == (leg_id as i32))
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let right_initial_state_positions = (n_initial_states + 1..=2 * n_initial_states)
            .map(|leg_id| {
                graph
                    .nodes()
                    .iter()
                    .position(|n| n.data.get_external_tag() == (leg_id as i32))
                    .unwrap()
            })
            .collect::<Vec<_>>();

        fn are_incoming_connected_to_outgoing<NodeColor>(
            cut: &[usize],
            graph: &SymbolicaGraph<NodeColor, EdgeColor>,
            adj_list: &HashMap<usize, Vec<(usize, usize)>>,
            left_nodes: &[usize],
            right_nodes: &[usize],
        ) -> bool {
            for left_node in left_nodes.iter() {
                for right_node in right_nodes.iter() {
                    let mut visited: Vec<bool> = vec![false; graph.nodes().len()];
                    let mut stack: Vec<usize> = vec![*left_node];
                    while let Some(node) = stack.pop() {
                        if node == *right_node {
                            return true;
                        }
                        visited[node] = true;
                        for (i_e, neighbor) in adj_list[&node].iter() {
                            if !cut.contains(i_e) && !visited[*neighbor] {
                                stack.push(*neighbor);
                            }
                        }
                    }
                }
            }
            false
        }

        fn are_nodes_connected<NodeColor>(
            cut: &[usize],
            graph: &SymbolicaGraph<NodeColor, EdgeColor>,
            adj_list: &HashMap<usize, Vec<(usize, usize)>>,
            nodes: &[usize],
        ) -> bool {
            for (node_a, node_b) in nodes.iter().tuple_combinations() {
                let mut visited: Vec<bool> = vec![false; graph.nodes().len()];
                let mut stack: Vec<usize> = vec![*node_a];
                let mut found_connecting_path = false;
                'dfs_search: while let Some(node) = stack.pop() {
                    if node == *node_b {
                        found_connecting_path = true;
                        break 'dfs_search;
                    }
                    visited[node] = true;
                    for (i_e, neighbor) in adj_list[&node].iter() {
                        if !cut.contains(i_e) && !visited[*neighbor] {
                            stack.push(*neighbor);
                        }
                    }
                }
                if !found_connecting_path {
                    return false;
                }
            }
            true
        }

        fn is_valid_cut<NodeColor>(
            cut: &[usize],
            graph: &SymbolicaGraph<NodeColor, EdgeColor>,
            adj_list: &HashMap<usize, Vec<(usize, usize)>>,
            left_nodes: &[usize],
            right_nodes: &[usize],
        ) -> bool {
            if are_incoming_connected_to_outgoing(cut, graph, adj_list, left_nodes, right_nodes) {
                return false;
            }
            if !are_nodes_connected(cut, graph, adj_list, left_nodes) {
                return false;
            }
            if !are_nodes_connected(cut, graph, adj_list, right_nodes) {
                return false;
            }
            true
        }

        for unique_combination_with_fixed_multiplicity in unique_combinations {
            'cutloop: for cut in unique_combination_with_fixed_multiplicity {
                if !is_valid_cut(
                    &cut,
                    graph,
                    &adj_list,
                    &left_initial_state_positions,
                    &right_initial_state_positions,
                ) {
                    continue 'cutloop;
                }
                // Make sure the cut is minimal and that no subcut is a valid cut.
                for subcut in (1..=(cut.len() - 1))
                    .flat_map(|offset| cut.iter().combinations(cut.len() - offset))
                {
                    if is_valid_cut(
                        subcut.iter().map(|&id| *id).collect::<Vec<_>>().as_slice(),
                        graph,
                        &adj_list,
                        &left_initial_state_positions,
                        &right_initial_state_positions,
                    ) {
                        continue 'cutloop;
                    }
                }
                return true;
            }
        }

        false
    }

    pub(crate) fn canonize_external_momenta_assignment(
        &self,
        model: &Model,
        node_colors_for_external_symmetrization: &HashMap<i32, i32>,
        graph: &mut SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
    ) {
        let mut initial_pdgs = self
            .initial_pdgs
            .iter()
            .enumerate()
            .map(|(i_n, pdg)| (pdg, i_n + 1))
            .collect::<Vec<_>>();
        let mut final_pdgs = if matches!(self.generation_type, GenerationType::CrossSection) {
            self.initial_pdgs
                .iter()
                .enumerate()
                .map(|(i_n, pdg)| (pdg, i_n + 1 + initial_pdgs.len()))
                .collect::<Vec<_>>()
        } else {
            self.final_pdgs_lists[0]
                .iter()
                .enumerate()
                .map(|(i_n, pdg)| (pdg, i_n + 1 + initial_pdgs.len()))
                .collect::<Vec<_>>()
        };

        let mut all_pdgs = initial_pdgs.clone();
        all_pdgs.extend(final_pdgs.clone());

        let mut new_node_data = vec![];
        let mut new_edge_data = vec![];

        // We do two passes, first assigning "specified externals" only (i.e. with tags > 0) and finally the remaining non specified ones (tags < 0)
        // The three pass steps are:
        // 0: distribute all externals with forced assignment, i.e. external_tag > 0
        // 1: distribute all externals with forced assignment up to left-right symmetry, i.e. external_tag < -2000
        // 2: distribute symmetrized externals, i.e. external_tag > -2000 && external_tag < 0

        for pass_steps in [0, 1, 2] {
            for (i_e, e) in graph.edges().iter().enumerate() {
                // All edges supposed to be incoming at this stage
                assert!(graph.nodes()[e.vertices.1].data.external_tag == 0);
                assert!(graph.nodes()[e.vertices.0].data.external_tag >= 0);
                let p = model.get_particle_from_pdg(e.data.pdg);
                let external_tag = graph.nodes()[e.vertices.0].data.external_tag;
                let symmetrized_external_tag = node_colors_for_external_symmetrization
                    .get(&external_tag)
                    .copied()
                    .unwrap_or(external_tag);
                let is_initial_state = external_tag <= self.initial_pdgs.len() as i32;
                let container = if is_initial_state {
                    &mut initial_pdgs
                } else {
                    &mut final_pdgs
                };
                if pass_steps == 0 && symmetrized_external_tag > 0 {
                    if self.symmetrize_left_right_states {
                        let matched_external_pos = all_pdgs
                            .iter()
                            .position(|(_pdg, i_ext)| (*i_ext as i32) == symmetrized_external_tag)
                            .unwrap();
                        all_pdgs.remove(matched_external_pos);
                    } else {
                        let matched_external_pos = container
                            .iter()
                            .position(|(_pdg, i_ext)| (*i_ext as i32) == symmetrized_external_tag)
                            .unwrap();
                        container.remove(matched_external_pos);
                    };
                } else if pass_steps == 1 && symmetrized_external_tag < -2000 {
                    // Node colors below -2000 indicate a left-right symmetrization of the external legs for a forward-scattering diagrams
                    // without symmetrizing the initial-states.
                    let external_leg_position = (-symmetrized_external_tag) % 1000;
                    // Try and find this either in initial or final states
                    let (matched_position, is_initial_match) = if let Some(matched_initial_pos) =
                        all_pdgs
                            .iter()
                            .position(|(_pdg, i_ext)| (*i_ext as i32) == external_leg_position)
                    {
                        (matched_initial_pos, true)
                    } else if let Some(matched_final_pos) =
                        all_pdgs.iter().position(|(_pdg, i_ext)| {
                            (*i_ext as i32)
                                == external_leg_position + (self.initial_pdgs.len() as i32)
                        })
                    {
                        (matched_final_pos, false)
                    } else {
                        unreachable!(
                            "Logical mistake in feyngen: external legs in canonicalized graphs should always be matchable."
                        )
                    };

                    let mut new_data = graph.nodes()[e.vertices.0].data.clone();
                    let new_external_tag = all_pdgs[matched_position].1 as i32;
                    // If we swapped initial and final state assignment, then we must also flip the pdg code of the corresponding half-edges
                    if is_initial_state != is_initial_match {
                        let mut e_data = e.data;
                        e_data.pdg = model
                            .get_particle_from_pdg(e_data.pdg)
                            .0
                            .get_anti_particle(model)
                            .0
                            .pdg_code;
                        new_edge_data.push((i_e, e_data));
                    }
                    new_data.set_external_tag(new_external_tag);
                    new_node_data.push((e.vertices.0, new_data));
                    all_pdgs.remove(matched_position);
                } else if pass_steps == 2 && (-2000..0).contains(&symmetrized_external_tag) {
                    let pdg_code = if is_initial_state {
                        p.0.pdg_code
                    } else {
                        p.0.get_anti_particle(model).0.pdg_code
                    };
                    if self.symmetrize_left_right_states {
                        let matched_external_pos: usize = all_pdgs
                            .iter()
                            .position(|(pdg, _i_ext)| **pdg == pdg_code as i64)
                            .unwrap();
                        let mut new_data = graph.nodes()[e.vertices.0].data.clone();
                        let new_external_tag = all_pdgs[matched_external_pos].1 as i32;
                        // If we swapped initial and final state assignment, then we must also flip the pdg code of the corresponding half-edges
                        if (new_external_tag > initial_pdgs.len() as i32
                            && new_data.external_tag <= initial_pdgs.len() as i32)
                            || (new_external_tag <= initial_pdgs.len() as i32
                                && new_data.external_tag > initial_pdgs.len() as i32)
                        {
                            let mut e_data = e.data;
                            e_data.pdg = model
                                .get_particle_from_pdg(e_data.pdg)
                                .0
                                .get_anti_particle(model)
                                .0
                                .pdg_code;
                            new_edge_data.push((i_e, e_data));
                        }
                        new_data.set_external_tag(new_external_tag);
                        new_node_data.push((e.vertices.0, new_data));
                        all_pdgs.remove(matched_external_pos);
                    } else {
                        let matched_external_pos = container
                            .iter()
                            .position(|(pdg, _i_ext)| **pdg == pdg_code as i64)
                            .unwrap();
                        let mut new_data = graph.nodes()[e.vertices.0].data.clone();
                        new_data.set_external_tag(container[matched_external_pos].1 as i32);
                        new_node_data.push((e.vertices.0, new_data));
                        container.remove(matched_external_pos);
                    }
                }
            }
        }

        for (node_pos, new_data) in new_node_data {
            graph.set_node_data(node_pos, new_data);
        }
        for (edge_pos, new_data) in new_edge_data {
            graph.set_edge_data(edge_pos, new_data);
        }
    }

    /// This function canonizes the edge and vertex ordering of a graph based on the skeletton graph with only propagator mass as edge color.
    /// This is useful to then allow for further grouping of isomorphic graphs, incl numerator.
    #[allow(clippy::type_complexity)]
    pub(crate) fn canonicalize_edge_and_vertex_ordering(
        &self,
        model: &Model,
        input_graph: &SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
        node_colors_for_external_symmetrization: &HashMap<i32, i32>,
        numerator_aware_isomorphism_grouping: &NumeratorAwareGraphGroupingOption,
        // This option contains: (self.symmetrize_initial_states, self.symmetrize_left_right_states) if any is true, else None.
        manually_canonalize_initial_states_cross_section_ordering: Option<(bool, bool)>,
    ) -> Result<
        (
            SymbolicaGraph<NodeColorWithoutVertexRule, std::string::String>,
            SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
        ),
        FeynGenError,
    > {
        // println!("INPUT GRAPH:\n{}", input_graph.to_dot());
        let mut canonized_skelettons = vec![];
        let external_node_positions = input_graph
            .nodes()
            .iter()
            .enumerate()
            .filter_map(|(i_n, n)| {
                if n.data.external_tag > 0 {
                    Some((n.data.external_tag, i_n))
                } else {
                    None
                }
            })
            .collect::<HashMap<_, _>>();
        // Contains the two tuple (are_left_and_right_swapped, remapping)
        let mut external_remappings = vec![];
        if let Some((symmetrize_initial_states, symmetrize_left_right_states)) =
            manually_canonalize_initial_states_cross_section_ordering
        {
            let external_ids = (1..=self.initial_pdgs.len()).collect::<Vec<_>>();
            let mut initial_states_orders = vec![];
            if symmetrize_initial_states {
                'permutation_loop: for permutation in external_ids
                    .iter()
                    .cloned()
                    .permutations(external_ids.len())
                {
                    if permutation
                        .iter()
                        .enumerate()
                        .any(|(src, trgt)| self.initial_pdgs[src] != self.initial_pdgs[*trgt - 1])
                    {
                        continue 'permutation_loop;
                    }
                    initial_states_orders.push(permutation);
                }
            } else {
                initial_states_orders.push(external_ids);
            };
            for initial_state_order in initial_states_orders {
                let mut remapping = AHashMap::<usize, usize>::default();
                for ext_id in 1..=self.initial_pdgs.len() {
                    remapping.insert(
                        *external_node_positions.get(&(ext_id as i32)).unwrap(),
                        *external_node_positions
                            .get(&((initial_state_order[ext_id - 1]) as i32))
                            .unwrap(),
                    );
                    remapping.insert(
                        *external_node_positions
                            .get(&((ext_id + self.initial_pdgs.len()) as i32))
                            .unwrap(),
                        *external_node_positions
                            .get(
                                &((initial_state_order[ext_id - 1] + self.initial_pdgs.len())
                                    as i32),
                            )
                            .unwrap(),
                    );
                }
                external_remappings.push((false, remapping));
                if symmetrize_left_right_states {
                    let mut remapping = AHashMap::<usize, usize>::default();
                    for ext_id in 1..=self.initial_pdgs.len() {
                        remapping.insert(
                            *external_node_positions.get(&(ext_id as i32)).unwrap(),
                            *external_node_positions
                                .get(
                                    &((initial_state_order[ext_id - 1] + self.initial_pdgs.len())
                                        as i32),
                                )
                                .unwrap(),
                        );
                        remapping.insert(
                            *external_node_positions
                                .get(&((ext_id + self.initial_pdgs.len()) as i32))
                                .unwrap(),
                            *external_node_positions
                                .get(&(initial_state_order[ext_id - 1] as i32))
                                .unwrap(),
                        );
                    }
                    external_remappings.push((true, remapping));
                }
            }
        } else {
            external_remappings.push((false, AHashMap::<usize, usize>::default()));
        }

        for (are_left_and_right_swapped, remapping) in external_remappings {
            // Make sure to canonize the edge ordering based on the skeletton graph with only propagator mass as edge color
            let mut skeletton_graph = SymbolicaGraph::new();
            for node in input_graph.nodes() {
                skeletton_graph.add_node(NodeColorWithoutVertexRule {
                    external_tag: *node_colors_for_external_symmetrization
                        .get(&node.data.external_tag)
                        .unwrap_or(&node.data.external_tag),
                });
            }
            let color_according_to_mass: bool = if let Some(grouping_options) =
                numerator_aware_isomorphism_grouping.get_options()
            {
                grouping_options.differentiate_particle_masses_only
            } else {
                true
            };
            for edge in input_graph.edges() {
                // We must maintain the directionality of the external edges to make sure that antifermions are not mapped to fermions
                // when grouping is done based on the skeletton graph and only the mass of the propagators
                // Also external edges must always retain their PDG since exchanging the momentum of e.g. a massless electron and muon
                // would not be a valid isomorphism w.r.t the process definition
                let is_edge_external = [edge.vertices.0, edge.vertices.1]
                    .iter()
                    .any(|e| input_graph.nodes()[*e].data.external_tag != 0);
                let mut remapped_edge_vertices: (usize, usize) = (
                    *remapping.get(&edge.vertices.0).unwrap_or(&edge.vertices.0),
                    *remapping.get(&edge.vertices.1).unwrap_or(&edge.vertices.1),
                );
                // Force all externals incoming in the canonicalization
                let particle = if input_graph.nodes()[edge.vertices.1].data.external_tag > 0 {
                    remapped_edge_vertices = (remapped_edge_vertices.1, remapped_edge_vertices.0);
                    model
                        .get_particle_from_pdg(edge.data.pdg)
                        .0
                        .get_anti_particle(model)
                } else {
                    model.get_particle_from_pdg(edge.data.pdg)
                };

                let keep_direction = is_edge_external;

                // // We want to forget the orientation of fermions (except for externals) to capture furry, but we need the canonization to be sensitive to charged vector bosons when it is connected to fermion lines
                // if !keep_direction && edge.directed && particle.is_vector() {
                //     let src_ferm = input_graph.node(edge.vertices.0).edges.iter().any(|a| {
                //         model
                //             .get_particle_from_pdg(input_graph.edge(*a).data.pdg)
                //             .is_fermion()
                //     });
                //     let sink_ferm = input_graph.node(edge.vertices.1).edges.iter().any(|a| {
                //         model
                //             .get_particle_from_pdg(input_graph.edge(*a).data.pdg)
                //             .is_fermion()
                //     });

                //     if src_ferm || sink_ferm {
                //         keep_direction = true;
                //     }
                // }

                // WARNING: It is important to note that the colouring choice below dictates what representative diagram (i.e. "sorted_g") will be used
                // for performing numerical comparisons for grouping isomorphic graphs. If we do not include the spin in the colouring, we may incorrectly
                // sort two isomorphic graphs with and interchange of massless quarks and gluons which will prevent their grouping (final result still correct, but more diagrams).
                // This is avoided by including the spin in the color, which will prevent massless quark and gluons from ever being interchanged.
                // Of course, even then, it could be that we pick two sorted representatives that are not isomorphic, even though there exist a different isomorphic skeletton graph
                // which would have given a match when doing the numerical comparison. In the SM, this should never happen, and for BSM we can afford to lose some grouping.
                // The only way to avoid this would be to numerically test all isomorphic permutations of the skeletton graph, which is prohibitively slow.
                skeletton_graph
                    .add_edge(
                        remapped_edge_vertices.0,
                        remapped_edge_vertices.1,
                        keep_direction,
                        if color_according_to_mass && !is_edge_external {
                            format!("{} | {}", particle.0.mass.0, particle.0.spin)
                        } else {
                            particle.0.name.to_string()
                        },
                    )
                    .unwrap();
            }
            canonized_skelettons.push((
                (are_left_and_right_swapped, remapping),
                skeletton_graph.canonize(),
            ));
        }
        canonized_skelettons.sort_by(|a, b| {
            (a.1.graph.nodes(), a.1.graph.edges()).cmp(&(b.1.graph.nodes(), b.1.graph.edges()))
        });
        let ((was_left_and_right_swapped, selected_external_remapping), canonized_skeletton) =
            canonized_skelettons.first().unwrap();

        let mut can_graph_node_pos_to_input_graph_node_pos: Vec<usize> =
            vec![0; input_graph.nodes().len()];
        for (input_graph_node_position, node_order) in
            canonized_skeletton.vertex_map.iter().enumerate()
        {
            can_graph_node_pos_to_input_graph_node_pos[*node_order] = *selected_external_remapping
                .get(&input_graph_node_position)
                .unwrap_or(&input_graph_node_position);
        }
        let mut input_graph_node_pos_to_can_graph_node_pos: Vec<usize> =
            vec![0; input_graph.nodes().len()];
        for (can_graph_node_pos, input_graph_node_pos) in can_graph_node_pos_to_input_graph_node_pos
            .iter()
            .enumerate()
        {
            input_graph_node_pos_to_can_graph_node_pos[*input_graph_node_pos] = can_graph_node_pos;
        }

        // Sort nodes according to the canonized skeleton graph
        // This will also ensure that EMR variables line up
        let mut sorted_g: SymbolicaGraph<NodeColorWithVertexRule, EdgeColor> =
            SymbolicaGraph::new();
        for &node_order in can_graph_node_pos_to_input_graph_node_pos.iter() {
            let node_data = input_graph.nodes()[*selected_external_remapping
                .get(&node_order)
                .unwrap_or(&node_order)]
            .data
            .clone();
            sorted_g.add_node(node_data);
        }

        let input_graph_nodes = input_graph.nodes();
        let mut reordered_edges = input_graph
            .edges()
            .iter()
            .map(|e| {
                // Set all external edges incoming and internal edges going from lower to higher node order
                let is_external = input_graph_nodes[e.vertices.0].data.external_tag > 0
                    || input_graph_nodes[e.vertices.1].data.external_tag > 0;
                let is_external_outgoing = input_graph_nodes[e.vertices.1].data.external_tag > 0;
                let is_flipped = is_external_outgoing
                    || (!is_external
                        && input_graph_node_pos_to_can_graph_node_pos[e.vertices.0]
                            > input_graph_node_pos_to_can_graph_node_pos[e.vertices.1]);
                let mut particle = if is_flipped {
                    model
                        .get_particle_from_pdg(e.data.pdg)
                        .0
                        .get_anti_particle(model)
                } else {
                    model.get_particle_from_pdg(e.data.pdg)
                };
                // Apply the switch from particle to anti-particle (CP symmetry) if the canonicalization swapped initial and final states.
                // IMPORTANT: we must here to apply CP not just to the external particles since the assignment of fermion flow matters, for
                // internal edges connecting these external fermions.
                // If one would fix those edges (necessary e.g. for e- d > e- d DIS process) then one could add '&& is_external' below, which
                // would allow to capture additional groupings involving the charged. W bosons.
                if *was_left_and_right_swapped {
                    particle = particle.0.get_anti_particle(model);
                }
                (
                    (e, is_flipped),
                    // This key will serve to give a unique ordering of the edges
                    (
                        if !is_flipped {
                            input_graph_node_pos_to_can_graph_node_pos[e.vertices.0]
                        } else {
                            input_graph_node_pos_to_can_graph_node_pos[e.vertices.1]
                        },
                        if !is_flipped {
                            input_graph_node_pos_to_can_graph_node_pos[e.vertices.1]
                        } else {
                            input_graph_node_pos_to_can_graph_node_pos[e.vertices.0]
                        },
                        particle.0.pdg_code,
                    ),
                )
            })
            .collect::<Vec<_>>();
        reordered_edges.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        for &((_e, _is_flipped), (v_in, v_out, pdg)) in reordered_edges.iter() {
            // We must canonize the fermion flow as well, so we force the fermion flow to always go from the lower to the higher node order,
            // and adjust the edge colour (particle vs anti-particle) accordingly.
            // This is necessary to ensure that the meaning of the edge momentum representation (which always aligns the monentum with the fermion flow)
            // is preserved
            sorted_g
                .add_edge(
                    v_in,
                    v_out,
                    // Now that the edge orientation is canonized we set all edges as directed because
                    // it is semantic even for bosons as it dictates the flow of the momentum
                    true,
                    EdgeColor::from_particle(model.get_particle_from_pdg(pdg)),
                )
                .unwrap();
        }

        // In order for the external assignment momenta to line up between the canonized versions of these graphs
        self.canonize_external_momenta_assignment(
            model,
            node_colors_for_external_symmetrization,
            &mut sorted_g,
        );

        Ok((canonized_skeletton.graph.to_owned(), sorted_g))
    }

    // This function normalize fermion and ghost flows and makes sure that all virtual charged particles without flows (like the W boson) are particles
    pub(crate) fn normalize_flows(
        &self,
        graph: &SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
        model: &Model,
    ) -> Result<(SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>, bool), FeynGenError> {
        let mut adj_map: HashMap<usize, Vec<(usize, usize)>> = HashMap::default();
        for (i_e, e) in graph.edges().iter().enumerate() {
            // Build an adjacency list including only fermions
            let p = model.get_particle_from_pdg(e.data.pdg);
            if !(p.0.is_fermion() || p.0.is_ghost()) {
                continue;
            }
            adj_map
                .entry(e.vertices.0)
                .or_default()
                .push((i_e, e.vertices.1));
            if e.vertices.0 != e.vertices.1 {
                adj_map
                    .entry(e.vertices.1)
                    .or_default()
                    .push((i_e, e.vertices.0));
            }
        }

        let mut vetoed_edges: Vec<bool> = graph
            .edges()
            .iter()
            .map(|e| {
                let p = model.get_particle_from_pdg(e.data.pdg);
                !(p.0.is_fermion() || p.0.is_ghost())
            })
            .collect();
        let mut new_edges: AHashMap<usize, (usize, usize, bool, EdgeColor)> = AHashMap::default();
        let mut normalized_graph: SymbolicaGraph<NodeColorWithVertexRule, EdgeColor> =
            SymbolicaGraph::new();
        for n in graph.nodes().iter() {
            normalized_graph.add_node(n.data.clone());
        }

        let graph_edges = graph.edges();
        // Pairing of the external fermion flows. The keys are (sorted) two-tuple of the PDGs of the external fermions (assuming all incoming)
        // and the values are the external leg ids of the fermions connected and in the order of the sorted key.
        #[allow(clippy::type_complexity)]
        let mut external_fermion_flow_pairings: AHashMap<
            (ArcParticle, ArcParticle),
            Vec<(usize, usize)>,
        > = AHashMap::default();
        // First fix flows connected to external only and after that fix all internal fermion/ghost flows
        let mut external_tags_ordering = (1..=self.initial_pdgs.len()).collect::<Vec<_>>();
        if self.generation_type == GenerationType::CrossSection {
            external_tags_ordering
                .extend((1..=self.initial_pdgs.len()).map(|i| i + self.initial_pdgs.len()));
        } else {
            external_tags_ordering
                .extend((1..=self.final_pdgs_lists[0].len()).map(|i| i + self.initial_pdgs.len()));
        }
        external_tags_ordering.push(0);
        for external_tag_to_consider in external_tags_ordering {
            for (i_e, e) in graph_edges.iter().enumerate() {
                let edge_tag = if graph.nodes()[e.vertices.0].data.external_tag != 0 {
                    graph.nodes()[e.vertices.0].data.external_tag
                } else {
                    graph.nodes()[e.vertices.1].data.external_tag
                };
                assert!(edge_tag >= 0);

                let is_a_virtual_edge = edge_tag == 0;
                if external_tag_to_consider != (edge_tag as usize) {
                    continue;
                }
                if vetoed_edges[i_e] {
                    if !new_edges.contains_key(&i_e) {
                        new_edges.insert(i_e, (e.vertices.0, e.vertices.1, e.directed, e.data));
                    }
                    continue;
                }
                let starting_particle = model.get_particle_from_pdg(e.data.pdg);
                let mut is_starting_antiparticle = starting_particle.0.is_antiparticle();
                let starting_vertices = if is_a_virtual_edge && is_starting_antiparticle {
                    // Force all virtual closed fermion loops to be particles
                    is_starting_antiparticle = false;
                    if !new_edges.contains_key(&i_e) {
                        new_edges.insert(
                            i_e,
                            (
                                e.vertices.1,
                                e.vertices.0,
                                e.directed,
                                EdgeColor::from_particle(
                                    model
                                        .get_particle_from_pdg(e.data.pdg)
                                        .0
                                        .get_anti_particle(model),
                                ),
                            ),
                        );
                        (e.vertices.1, e.vertices.0)
                    } else {
                        continue;
                    }
                } else if !new_edges.contains_key(&i_e) {
                    new_edges.insert(i_e, (e.vertices.0, e.vertices.1, e.directed, e.data));
                    (e.vertices.0, e.vertices.1)
                } else {
                    continue;
                };

                vetoed_edges[i_e] = true;
                let mut connected_leg_ids: AHashSet<usize> = AHashSet::new();
                for read_to_the_right in [true, false] {
                    let mut previous_node_position = if read_to_the_right {
                        starting_vertices.1
                    } else {
                        starting_vertices.0
                    };
                    // println!(
                    //     "Starting reading chain from edge {}->{}, from {}",
                    //     starting_vertices.0, starting_vertices.1, previous_node_position
                    // );
                    if graph.nodes()[previous_node_position].data.external_tag > 0 {
                        connected_leg_ids.insert(
                            graph.nodes()[previous_node_position].data.external_tag as usize,
                        );
                    }

                    'fermion_chain: loop {
                        let (next_fermion_edge_position, next_fermion_node_position) =
                            follow_chain(
                                previous_node_position,
                                &mut vetoed_edges,
                                &adj_map,
                                true,
                            )?;
                        if graph.nodes()[next_fermion_node_position].data.external_tag > 0 {
                            connected_leg_ids.insert(
                                graph.nodes()[next_fermion_node_position].data.external_tag
                                    as usize,
                            );
                        }
                        // println!(
                        //     "Next edge: {}",
                        //     if let Some(nfep) = next_fermion_edge_position {
                        //         format!(
                        //             "{}, {} -> {}",
                        //             graph_edges[nfep].vertices.0, graph_edges[nfep].vertices.1, nfep
                        //         )
                        //     } else {
                        //         "None".to_string()
                        //     }
                        // );
                        if let Some(nfep) = next_fermion_edge_position {
                            let this_has_same_orientation_starting = (read_to_the_right
                                && graph_edges[nfep].vertices.1 == next_fermion_node_position)
                                || (!read_to_the_right
                                    && graph_edges[nfep].vertices.0 == next_fermion_node_position);
                            if !new_edges.contains_key(&nfep) {
                                if this_has_same_orientation_starting {
                                    new_edges.insert(
                                        nfep,
                                        (
                                            graph_edges[nfep].vertices.0,
                                            graph_edges[nfep].vertices.1,
                                            graph_edges[nfep].directed,
                                            graph_edges[nfep].data,
                                        ),
                                    );
                                } else {
                                    let this_edge_particle =
                                        model.get_particle_from_pdg(graph_edges[nfep].data.pdg);
                                    new_edges.insert(
                                        nfep,
                                        (
                                            graph_edges[nfep].vertices.1,
                                            graph_edges[nfep].vertices.0,
                                            graph_edges[nfep].directed,
                                            // Force fermion to be the same species as the one starting the chain
                                            if this_edge_particle.0.is_antiparticle()
                                                != is_starting_antiparticle
                                            {
                                                EdgeColor::from_particle(
                                                    this_edge_particle.0.get_anti_particle(model),
                                                )
                                            } else {
                                                EdgeColor::from_particle(this_edge_particle)
                                            },
                                        ),
                                    );
                                }
                            }
                            previous_node_position = next_fermion_node_position;
                        } else {
                            break 'fermion_chain;
                        }
                    }
                }
                if external_tag_to_consider > 0 && starting_particle.0.is_fermion() {
                    if connected_leg_ids.len() != 2 {
                        return Err(FeynGenError::GenericError(
                            "External fermion flow must have exactly two legs".to_string(),
                        ));
                    }
                    let connected_leg_ids_vec =
                        connected_leg_ids.iter().copied().collect::<Vec<_>>();
                    let connected_leg_pdgs = (0..=1)
                        .map(|i| {
                            if connected_leg_ids_vec[i] <= self.initial_pdgs.len() {
                                model.get_particle_from_pdg(
                                    self.initial_pdgs[connected_leg_ids_vec[i] - 1] as isize,
                                )
                            } else {
                                // Assign line types assuming all incoming particles
                                let right_side_pdgs =
                                    if self.generation_type == GenerationType::CrossSection {
                                        &self.initial_pdgs
                                    } else {
                                        &self.final_pdgs_lists[0]
                                    };
                                model
                                    .get_particle_from_pdg(
                                        right_side_pdgs
                                            [connected_leg_ids_vec[i] - self.initial_pdgs.len() - 1]
                                            as isize,
                                    )
                                    .0
                                    .get_anti_particle(model)
                            }
                        })
                        .collect::<Vec<_>>();
                    if connected_leg_pdgs
                        .iter()
                        .any(|particle| !particle.0.is_fermion())
                    {
                        return Err(FeynGenError::GenericError(
                            "External fermion flow must connect two fermions".to_string(),
                        ));
                    }
                    if connected_leg_pdgs[0].0.is_antiparticle()
                        && !connected_leg_pdgs[1].0.is_antiparticle()
                    {
                        external_fermion_flow_pairings
                            .entry((connected_leg_pdgs[0].clone(), connected_leg_pdgs[1].clone()))
                            .or_default()
                            .push((connected_leg_ids_vec[0], connected_leg_ids_vec[1]));
                    } else if connected_leg_pdgs[1].0.is_antiparticle()
                        && !connected_leg_pdgs[0].0.is_antiparticle()
                    {
                        external_fermion_flow_pairings
                            .entry((connected_leg_pdgs[1].clone(), connected_leg_pdgs[0].clone()))
                            .or_default()
                            .push((connected_leg_ids_vec[1], connected_leg_ids_vec[0]));
                    } else {
                        return Err(FeynGenError::GenericError(
                            "External fermion flow must connect a fermion and an anti-fermion. GammaLoop has no support for Majorana particles yet.".to_string(),
                        ));
                    }
                }
            }
        }

        let mut external_fermion_flow_sign = 1;
        let mut concatenated_lines = vec![];
        for (_line_type, lines) in external_fermion_flow_pairings.iter().sorted() {
            for (a, b) in lines.iter() {
                concatenated_lines.push(*a);
                concatenated_lines.push(*b);
            }
        }
        if Permutation::sort(&concatenated_lines)
            .transpositions()
            .len()
            % 2
            == 1
        {
            external_fermion_flow_sign *= -1;
        }

        // Finally make sure that all virtual non-self-antiparticles that are not fermions or ghost (like the W-boson) are set to particles
        for (i_e, e) in graph_edges.iter().enumerate() {
            let is_a_virtual_edge = graph.nodes()[e.vertices.0].data.external_tag == 0
                && graph.nodes()[e.vertices.1].data.external_tag == 0;
            let this_edge_particle = model.get_particle_from_pdg(e.data.pdg);
            if is_a_virtual_edge
                && this_edge_particle.0.is_antiparticle()
                && !(this_edge_particle.0.is_fermion() || this_edge_particle.0.is_ghost())
            {
                new_edges.insert(
                    i_e,
                    (
                        e.vertices.1,
                        e.vertices.0,
                        e.directed,
                        EdgeColor::from_particle(this_edge_particle.0.get_anti_particle(model)),
                    ),
                );
            }
        }

        let mut sorted_new_edges = new_edges.iter().collect::<Vec<_>>();
        sorted_new_edges.sort_by(|(i_e_0, _), (i_e_1, _)| i_e_0.cmp(i_e_1));

        for (v0, v1, directed, edge_data) in sorted_new_edges.iter().map(|(_, v)| v) {
            normalized_graph
                .add_edge(*v0, *v1, *directed, *edge_data)
                .map_err(|e| FeynGenError::GenericError(e.to_string()))?;
        }
        assert!(normalized_graph.edges().len() == graph.edges().len());

        Ok((normalized_graph, external_fermion_flow_sign == -1))
    }

    // Note, this function will not work as intended with four-fermion vertices, and only aggregated self-loops or fermion-loops not involving four-femion vertices
    pub(crate) fn count_closed_fermion_loops(
        &self,
        graph: &SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
        model: &Model,
    ) -> Result<usize, FeynGenError> {
        let mut adj_map: HashMap<usize, Vec<(usize, usize)>> = HashMap::default();
        for (i_e, e) in graph.edges().iter().enumerate() {
            // Build an adjacency list including only fermions
            if !model.get_particle_from_pdg(e.data.pdg).0.is_fermion() {
                continue;
            }
            adj_map
                .entry(e.vertices.0)
                .or_default()
                .push((i_e, e.vertices.1));
            if e.vertices.0 != e.vertices.1 {
                adj_map
                    .entry(e.vertices.1)
                    .or_default()
                    .push((i_e, e.vertices.0));
            }
        }

        let mut vetoed_edges: Vec<bool> = graph
            .edges()
            .iter()
            .map(|e| !model.get_particle_from_pdg(e.data.pdg).0.is_fermion())
            .collect();
        let mut n_fermion_loops = 0;
        for (i_e, e) in graph.edges().iter().enumerate() {
            if vetoed_edges[i_e] {
                continue;
            }
            vetoed_edges[i_e] = true;
            let (_, left_trail_end) =
                follow_chain(e.vertices.0, &mut vetoed_edges, &adj_map, false)?;
            let (_, right_trail_end) =
                follow_chain(e.vertices.1, &mut vetoed_edges, &adj_map, false)?;
            if left_trail_end == right_trail_end {
                n_fermion_loops += 1;
            }
        }
        Ok(n_fermion_loops)
    }

    #[instrument(skip_all)]
    pub fn generate(
        &self,
        model: &Model,
        settings: &GlobalSettings,
    ) -> Result<Vec<Graph>, FeynGenError> {
        let num_threads = Some(settings.n_cores.feyngen);
        let progress_bar_style = ProgressStyle::with_template(
            "[{elapsed_precise} | ETA: {eta_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} ({percent}%) {msg}",
        )
        .unwrap();
        let mut cpl_reps: Vec<Replacement> = vec![];
        for cpl in model.couplings.values() {
            let [lhs, rhs] = cpl.rep_rule();
            cpl_reps.push(Replacement::new(lhs.to_pattern(), rhs));
        }

        if self.final_pdgs_lists.is_empty() {
            return Err(FeynGenError::GenericError(
                "No specification of final states is provided".to_string(),
            ));
        }

        let representative_final_pdgs = if self.generation_type == GenerationType::CrossSection {
            if self.final_pdgs_lists.len() > 1
                && self.final_pdgs_lists.iter().any(|fs| fs.is_empty())
            {
                return Err(FeynGenError::GenericError(
                        "When specifying multiple possible set of final states in cross-section generation, each must contain at least one particle".to_string(),
                    ));
            }
            self.final_pdgs_lists[0].clone()
        } else {
            if self.final_pdgs_lists.len() > 1 {
                return Err(FeynGenError::GenericError(
                    "Multiple selection of possible final states is not supported for amplitude generation".to_string(),
                ));
            }
            self.final_pdgs_lists[0].clone()
        };

        // Build a custom ThreadPool. If `Some(threads)` is given, use that;
        // otherwise, default to the number of logical CPUs on the system.
        let pool = match num_threads {
            Some(threads) => ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("Failed to build custom Rayon ThreadPool"),
            None => ThreadPoolBuilder::new()
                .build()
                .expect("Failed to build default Rayon ThreadPool"),
        };

        debug!(
            "Generating Feynman diagrams over {} cores for model {} and process:\n{}",
            pool.current_num_threads(),
            model.name,
            self
        );

        let filters = match self.generation_type {
            GenerationType::Amplitude => &self.amplitude_filters,
            GenerationType::CrossSection => &self.cross_section_filters,
        };

        let vertex_vetoes_filter = filters
            .0
            .iter()
            .filter_map(|f| {
                if let FeynGenFilter::VertexVeto(filter) = f {
                    Some(filter)
                } else {
                    None
                }
            })
            .next();

        const SB_INCOMING: bool = true;
        const SB_OUTGOING: bool = false;
        // const SB_INCOMING: bool = false;
        // const SB_OUTGOING: bool = true;

        #[allow(clippy::type_complexity)]
        let mut vertex_signatures: HashMap<
            Vec<(Option<bool>, SmartString<LazyCompact>)>,
            Vec<SmartString<LazyCompact>>,
        > = HashMap::default();
        'add_vertex_rules: for vertex_rule in model.vertex_rules.iter() {
            if let Some(veto) = vertex_vetoes_filter {
                if veto.contains(&vertex_rule.0.name.clone().into()) {
                    continue 'add_vertex_rules;
                }
            }
            let mut oriented_particles = vec![];
            for p in vertex_rule.0.particles.iter() {
                if p.0.is_self_antiparticle() {
                    oriented_particles.push((None, p.0.name.clone()));
                } else if p.0.is_antiparticle() {
                    oriented_particles.push((
                        Some(SB_INCOMING),
                        // Some(SB_OUTGOING),
                        p.0.get_anti_particle(model).0.name.clone(),
                    ));
                } else {
                    oriented_particles.push((
                        Some(SB_OUTGOING),
                        // Some(SB_INCOMING),
                        p.0.name.clone(),
                    ));
                }
                if let Some(vetoed_particles) = filters.get_particle_vetos()
                    && (vetoed_particles.contains(&(p.0.pdg_code as i64))
                        || vetoed_particles
                            .contains(&(p.0.get_anti_particle(model).0.pdg_code as i64)))
                {
                    continue 'add_vertex_rules;
                }
            }
            vertex_signatures
                .entry(oriented_particles)
                .or_default()
                .push(vertex_rule.0.name.clone());
        }

        // In Symbolica generation *all* particles have to be *particle* and never *anti-particles*.
        // The differentiation in symbolica graphs happens through the edge directions.
        // For that reason we must correctly set the directions of the external edges according to the nature of
        // the particle specified in the process definition.

        let mut external_edges = self
            .initial_pdgs
            .iter()
            .enumerate()
            .map(|(i_initial, pdg)| {
                let p = model.get_particle_from_pdg(*pdg as isize);
                (
                    NodeColorWithoutVertexRule {
                        external_tag: (i_initial + 1) as i32,
                    },
                    if p.0.is_self_antiparticle() {
                        HalfEdge {
                            direction: None,
                            data: EdgeColor::from_particle(p),
                        }
                    } else if p.0.is_antiparticle() {
                        HalfEdge {
                            direction: Some(SB_OUTGOING),
                            data: EdgeColor::from_particle(p.0.get_anti_particle(model)),
                        }
                    } else {
                        HalfEdge {
                            direction: Some(SB_INCOMING),
                            data: EdgeColor::from_particle(p),
                        }
                    },
                )
            })
            .collect::<Vec<_>>();

        let mut external_connections = vec![];
        match self.generation_type {
            GenerationType::Amplitude => {
                for i_initial in 1..=self.initial_pdgs.len() {
                    external_connections.push((Some(i_initial), None));
                }
                for i_final in (self.initial_pdgs.len() + 1)
                    ..=(self.initial_pdgs.len() + representative_final_pdgs.len())
                {
                    external_connections.push((None, Some(i_final)));
                }
                external_edges.extend(
                    representative_final_pdgs
                        .iter()
                        .enumerate()
                        .map(|(i_final, pdg)| {
                            let p = model.get_particle_from_pdg(*pdg as isize);
                            (
                                NodeColorWithoutVertexRule {
                                    external_tag: (self.initial_pdgs.len() + i_final + 1) as i32,
                                },
                                if p.0.is_self_antiparticle() {
                                    HalfEdge {
                                        direction: None,
                                        data: EdgeColor::from_particle(p),
                                    }
                                } else if p.0.is_antiparticle() {
                                    HalfEdge {
                                        direction: Some(SB_INCOMING),
                                        data: EdgeColor::from_particle(
                                            p.0.get_anti_particle(model),
                                        ),
                                    }
                                } else {
                                    HalfEdge {
                                        direction: Some(SB_OUTGOING),
                                        data: EdgeColor::from_particle(p),
                                    }
                                },
                            )
                        })
                        .collect::<Vec<_>>(),
                );
            }
            GenerationType::CrossSection => {
                let mut i_final = self.initial_pdgs.len();
                for (i_initial, pdg) in self.initial_pdgs.iter().enumerate() {
                    i_final += 1;
                    let p = model.get_particle_from_pdg(*pdg as isize);
                    external_connections.push((Some(i_initial + 1), Some(i_final)));
                    external_edges.push((
                        NodeColorWithoutVertexRule {
                            external_tag: i_final as i32,
                        },
                        if p.0.is_self_antiparticle() {
                            HalfEdge {
                                direction: None,
                                data: EdgeColor::from_particle(p),
                            }
                        } else if p.0.is_antiparticle() {
                            HalfEdge {
                                direction: Some(SB_INCOMING),
                                data: EdgeColor::from_particle(p.0.get_anti_particle(model)),
                            }
                        } else {
                            HalfEdge {
                                direction: Some(SB_OUTGOING),
                                data: EdgeColor::from_particle(p),
                            }
                        },
                    ));
                }
            }
        }

        // debug!("external_edges = {:?}", external_edges);
        // debug!("vertex_signatures = {:?}", vertex_signatures);

        // let external_edges_for_generation = external_edges
        //     .iter()
        //     .map(|(i, (orientation, name))| (i.clone(), (*orientation, name)))
        //     .collect::<Vec<_>>();
        let external_edges_for_generation = external_edges.clone();
        let mut vertex_signatures_for_generation = vertex_signatures
            .keys()
            .map(|v| {
                v.iter()
                    .map(|(orientation, p)| HalfEdge {
                        direction: *orientation,
                        data: EdgeColor::from_particle(model.get_particle(p)),
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        // Sort vertices to improve reproducibility.
        vertex_signatures_for_generation.sort();
        debug!(
            "generation_external_edges = {:?}",
            external_edges_for_generation
        );
        debug!(
            "generation_vertex_signatures = {:?}",
            vertex_signatures_for_generation
        );

        debug!("feygen options: {:#?}", self);

        // Record the start time
        let start = Instant::now();
        let mut last_step = start;

        info!(
            "{} | Δ={} | {:<95}",
            format!("{:<6}", utils::format_wdhms(0)).blue().bold(),
            format!("{:<6}", utils::format_wdhms(0)).blue(),
            "Starting Feynman graph generation with Symbolica..."
        );

        let graph_gen_bar = Arc::new(ProgressBar::new(0_u64));
        graph_gen_bar.set_style(progress_bar_style.clone());
        graph_gen_bar.set_message("Starting Feynman graphs generation with Symbolica...");
        graph_gen_bar.tick();
        let graph_gen_bar_arc_clone = graph_gen_bar.clone();
        let symbolica_generation_settings = if let Some(max_bridges) = filters.get_max_bridge() {
            GenerationSettings::new().max_bridges(max_bridges)
        } else {
            GenerationSettings::new()
        }
        .progress_fn(Box::new(move |_| {
            graph_gen_bar_arc_clone.inc_length(1);
            graph_gen_bar_arc_clone.inc(1);
            !INTERRUPTED.swap(false, Ordering::SeqCst)
        }))
        .max_loops(self.loop_count_range.1)
        .allow_self_loops(!self.filter_self_loop)
        .allow_zero_flow_edges(!self.filter_zero_flow_edges)
        .abort_check(Box::new(|| {
            // Check for where keyboard interrupt has been triggered
            if INTERRUPTED.swap(false, Ordering::SeqCst) {
                // status_warn!(
                //     "Generation aborted by the user via keyboard interrupt (Ctrl+C). Returning graphs generated thus far, BUT THE RESULT WILL BE INCOMPLETE"
                // );
                true
            } else {
                false
            }
        }));

        // println!("max_bridges = {:#?}", filters.get_max_bridge());
        // println!("loop_count = {:#?}", self.loop_count_range.1);
        // println!("self.filter_self_loop = {:#?}", self.filter_self_loop);
        // println!(
        //     "len external_edges_for_generation = {}",
        //     external_edges_for_generation.len()
        // );
        // println!(
        //     "len external_edges_for_generation = {}",
        //     vertex_signatures_for_generation.len()
        // );

        // println!(
        //     "external_edges_for_generation={:?}",
        //     external_edges_for_generation
        // );
        // println!(
        //     "vertex_signatures_for_generation=\n{:?}",
        //     vertex_signatures_for_generation
        // );

        // println!(
        //     "vertex_signatures_for_generation bis=\nvec![{:?}]",
        //     vertex_signatures_for_generation
        //         .iter()
        //         .map(|v| {
        //             format!(
        //                 "vec![{}]",
        //                 v.iter()
        //                     .map(|e| format!(
        //                         "{:?}",
        //                         HalfEdge {
        //                             direction: e.direction,
        //                             data: e.data.pdg
        //                         }
        //                     ),)
        //                     .collect::<Vec<_>>()
        //                     .join(", ")
        //             )
        //         })
        //         .collect::<Vec<_>>()
        //         .join(", ")
        // );

        // println!(
        //     "vertex_signatures_for_generation bis=\n{:?}",
        //     vertex_signatures_for_generation
        //         .iter()
        //         .map(|v| {
        //             v.iter()
        //                 .map(|he| {
        //                     (
        //                         he.direction
        //                             .map(|d| format!("{}", d))
        //                             .unwrap_or("None".to_string()),
        //                         format!(
        //                             "{}|{}",
        //                             he.data.pdg,
        //                             model.get_particle_from_pdg(he.data.pdg).0.name.to_string()
        //                         ),
        //                     )
        //                 })
        //                 .collect::<Vec<_>>()
        //         })
        //         .collect::<Vec<_>>()
        // );
        // println!(
        //     "Symbolica generation settings:\n{:?}",
        //     symbolica_generation_settings
        // );
        let mut graphs = match SymbolicaGraph::generate(
            external_edges_for_generation.as_slice(),
            vertex_signatures_for_generation.as_slice(),
            symbolica_generation_settings,
        ) {
            Ok(gs) => {
                if gs.is_empty() {
                    error!(
                        "\nSymbolica graph generation for this process did not yield any graph. This likely indicates a problem with the model or process specification. Please check the model and process and try again.\n"
                    );
                    return Err(FeynGenError::GenericError(
                        "Symbolica graph generation failed to generate any graphs".to_string(),
                    ));
                }
                gs
            }
            Err(gs_thus_far) => {
                if gs_thus_far.is_empty() {
                    error!(
                        "\nSymbolica graph generation for this process did not yield any graph before being interrupted by the user.\n"
                    );
                    return Err(FeynGenError::GenericError(
                        "Symbolica graph generation failed to generate any graphs".to_string(),
                    ));
                }
                error!(
                    "\nSymbolica graph generation was aborted by the user after generating {} graphs. Generation will continue with these only, {}\n",
                    gs_thus_far.len(),
                    "BUT THE RESULT WILL BE INCOMPLETE".red().bold()
                );
                gs_thus_far
            }
        };
        graph_gen_bar.finish_and_clear();

        // Immediately drop lower loop count contributions
        graphs.retain(|g, _| g.num_loops() >= self.loop_count_range.0);
        let mut step = Instant::now();
        info!(
            "{} | Δ={} | {:<95}{}",
            format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                .blue()
                .bold(),
            format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
            "Number of graphs resulting from Symbolica generation:",
            format!("{}", graphs.len()).green().bold()
        );
        last_step = step;

        let unoriented_final_state_particles_lists = if self.generation_type
            == GenerationType::CrossSection
            && !representative_final_pdgs.is_empty()
        {
            self.final_pdgs_lists
                .iter()
                .map(|final_pdgs| {
                    final_pdgs
                        .iter()
                        .map(|pdg| {
                            let p = model.get_particle_from_pdg(*pdg as isize);
                            if p.0.is_antiparticle() {
                                p.0.get_anti_particle(model).0.pdg_code
                            } else {
                                p.0.pdg_code
                            }
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        } else {
            vec![]
        };

        if !unoriented_final_state_particles_lists.is_empty() {
            let mut unoriented_final_state_particles_always_present: HashSet<isize> =
                HashSet::from_iter(unoriented_final_state_particles_lists[0].iter().cloned());
            for final_pdgs in unoriented_final_state_particles_lists.iter().skip(1) {
                unoriented_final_state_particles_always_present =
                    unoriented_final_state_particles_always_present
                        .intersection(&HashSet::from_iter(final_pdgs.iter().cloned()))
                        .cloned()
                        .collect();
            }
            let unoriented_final_state_particles_always_present_vec =
                unoriented_final_state_particles_always_present
                    .iter()
                    .cloned()
                    .collect::<Vec<_>>();
            if !unoriented_final_state_particles_always_present_vec.is_empty() {
                let bar = ProgressBar::new(graphs.len() as u64);
                bar.set_style(progress_bar_style.clone());
                bar.set_message("Enforcing particle content...");
                pool.install(|| {
                    graphs = graphs
                        .par_iter()
                        .progress_with(bar.clone())
                        .filter(|(g, _symmetry_factor)| {
                            contains_particles(
                                g,
                                unoriented_final_state_particles_always_present_vec.as_slice(),
                            )
                        })
                        .map(|(g, sf)| (g.clone(), sf.clone()))
                        .collect::<HashMap<_, _>>()
                });
                bar.finish_and_clear();

                step = Instant::now();
                info!(
                    "{} | Δ={} | {:<95}{}",
                    format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                        .blue()
                        .bold(),
                    format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
                    "Number of graphs retained after enforcing supergraph particle content:",
                    format!("{}", graphs.len()).green()
                );
                last_step = step;
            }
        }

        let tadpole_filter = filters
            .0
            .iter()
            .filter_map(|f| {
                if let FeynGenFilter::TadpolesFilter(opts) = f {
                    Some(opts)
                } else {
                    None
                }
            })
            .next();
        let external_self_energy_filter = filters
            .0
            .iter()
            .filter_map(|f| {
                if let FeynGenFilter::SelfEnergyFilter(opts) = f {
                    Some(opts)
                } else {
                    None
                }
            })
            .next();
        let zero_snails_filter: Option<&SnailFilterOptions> = filters
            .0
            .iter()
            .filter_map(|f| {
                if let FeynGenFilter::ZeroSnailsFilter(opts) = f {
                    Some(opts)
                } else {
                    None
                }
            })
            .next();
        let factorized_loop_topologies_count_range = filters
            .0
            .iter()
            .filter_map(|f| {
                if let FeynGenFilter::FactorizedLoopTopologiesCountRange(range) = f {
                    Some(range)
                } else {
                    None
                }
            })
            .next();

        if tadpole_filter.is_some()
            || external_self_energy_filter.is_some()
            || zero_snails_filter.is_some()
            || factorized_loop_topologies_count_range.is_some()
        {
            let bar = ProgressBar::new(graphs.len() as u64);
            bar.set_style(progress_bar_style.clone());
            bar.set_message("Vetoing special topologies...");
            pool.install(|| {
                graphs = graphs
                    .par_iter()
                    .progress_with(bar.clone())
                    .filter(|(g, _symmetry_factor)| {
                        !veto_special_topologies(
                            model,
                            g,
                            external_self_energy_filter,
                            tadpole_filter,
                            zero_snails_filter,
                            factorized_loop_topologies_count_range,
                        )
                    })
                    .map(|(g, sf)| (g.clone(), sf.clone()))
                    .collect::<HashMap<_, _>>()
            });
            bar.finish_and_clear();
        }

        // Re-interpret the symmetry factor as a multiplier where the symbolica factor from symbolica appears in the denominator
        let graphs = graphs
            .iter()
            .map(|(g, symmetry_factor)| {
                (
                    g.clone(),
                    Atom::num(1) / function!(symbol!("AutG"), Atom::num(symmetry_factor.clone())),
                )
            })
            .collect::<HashMap<_, _>>();

        step = Instant::now();
        info!(
            "{} | Δ={} | {:<95}{}",
            format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                .blue()
                .bold(),
            format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
            "Number of graphs retained after removal of vetoed topologies:",
            format!("{}", graphs.len()).green()
        );
        last_step = step;

        #[allow(clippy::type_complexity)]
        let mut node_colors: HashMap<
            Vec<(Option<bool>, SmartString<LazyCompact>)>,
            Vec<SmartString<LazyCompact>>,
        > = HashMap::default();
        for (v_legs, v_colors) in vertex_signatures.iter() {
            let mut sorted_ps = v_legs.clone();
            sorted_ps.sort();
            node_colors.insert(sorted_ps, v_colors.clone());
        }

        let mut processed_graphs = vec![];
        for (g, symmetry_factor) in graphs.iter() {
            for (colored_g, multiplicity) in assign_node_colors(model, g, &node_colors)? {
                processed_graphs.push((
                    colored_g.canonize().graph,
                    (Atom::num(multiplicity as i64) * symmetry_factor),
                ));
            }
        }

        let bar = ProgressBar::new(graphs.len() as u64);
        bar.set_style(progress_bar_style.clone());
        bar.set_message("Assigning interactions to vertices...");
        let mut processed_graphs = vec![];
        let lists_of_entries = pool.install(|| {
            graphs
                .par_iter()
                .progress_with(bar.clone())
                .map(|(g, symmetry_factor)| {
                    assign_node_colors(model, g, &node_colors).map(|colored_graphs| {
                        colored_graphs
                            .iter()
                            .map(|(colored_g, multiplicity)| {
                                (
                                    colored_g.canonize().graph,
                                    if *multiplicity != 1 {
                                        function!(
                                            symbol!("CouplingsMultiplicity"),
                                            Atom::num(*multiplicity as i64)
                                        ) * symmetry_factor
                                    } else {
                                        symmetry_factor.clone()
                                    },
                                )
                            })
                            .collect::<Vec<_>>()
                    })
                })
                .collect::<Result<Vec<_>, FeynGenError>>()
        })?;
        bar.finish_and_clear();

        for entries in lists_of_entries {
            processed_graphs.extend(entries);
        }
        processed_graphs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        step = Instant::now();
        info!(
            "{} | Δ={} | {:<95}{}",
            format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                .blue()
                .bold(),
            format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
            "Number of graphs after vertex info assignment:",
            format!("{}", processed_graphs.len()).green()
        );
        last_step = step;

        filters.apply_filters(&mut processed_graphs, model, &pool, &progress_bar_style)?;

        step = Instant::now();
        info!(
            "{} | Δ={} | {:<95}{}",
            format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                .blue()
                .bold(),
            format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
            "Number of graphs after all complete graphs filters are applied:",
            format!("{}", processed_graphs.len()).green()
        );
        last_step = step;

        let fermion_loop_count_range_filter = filters.get_fermion_loop_count_range();
        let bar = ProgressBar::new(processed_graphs.len() as u64);
        bar.set_style(progress_bar_style.clone());
        bar.set_message(
            "Analyzing closed fermion chains to capture antisymmetry and apply fermion filters...",
        );
        processed_graphs = pool.install(|| {
            processed_graphs
                .iter()
                .filter_map(|(g, symmetry_factor)| {
                    match self.count_closed_fermion_loops(g, model) {
                        Ok(n_closed_fermion_loops) => {
                            let new_symmetry_factor = if n_closed_fermion_loops % 2 == 1 {
                                function!(symbol!("InternalFermionLoopSign"), -1) * symmetry_factor
                            } else {
                                symmetry_factor.clone()
                            };
                            if let Some((min_n_fermion_loops, max_n_fermion_loops)) =
                                fermion_loop_count_range_filter
                            {
                                if n_closed_fermion_loops >= min_n_fermion_loops
                                    && n_closed_fermion_loops <= max_n_fermion_loops
                                {
                                    Some(Ok((g.clone(), new_symmetry_factor)))
                                } else {
                                    None
                                }
                            } else {
                                Some(Ok((g.clone(), new_symmetry_factor)))
                            }
                        }
                        Err(e) => Some(Err(e)),
                    }
                })
                .collect::<Result<Vec<_>, _>>()
        })?;
        bar.finish_and_clear();
        step = Instant::now();
        info!(
            "{} | Δ={} | {:<95}{}",
            format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                .blue()
                .bold(),
            format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
            "Number of graphs after closed fermion chains analysis:",
            format!("{}", processed_graphs.len()).green()
        );
        last_step = step;

        let (n_unresolved, unresolved_type) = self.unresolved_cut_content(model);

        // The fast cutksoky filter is only fast for up to ~ 6 particles to check
        let mut applied_fast_cutksosky_cut_filter = false;
        if self.generation_type == GenerationType::CrossSection
            && !representative_final_pdgs.is_empty()
            && self
                .final_pdgs_lists
                .iter()
                .map(|pdgs| pdgs.len())
                .max()
                .unwrap_or(0)
                + n_unresolved
                <= self.max_multiplicity_for_fast_cut_filter
        {
            applied_fast_cutksosky_cut_filter = true;

            let bar = ProgressBar::new(processed_graphs.len() as u64);
            bar.set_style(progress_bar_style.clone());
            bar.set_message("Applying fast Cutkosky cut filter...");
            pool.install(|| {
                processed_graphs = processed_graphs
                    .par_iter()
                    .progress_with(bar.clone())
                    .filter(|(g, _symmetry_factor)| {
                        unoriented_final_state_particles_lists.iter().any(
                            |unoriented_final_state_particles| {
                                self.contains_cut_fast(
                                    model,
                                    g,
                                    n_unresolved,
                                    &unresolved_type,
                                    unoriented_final_state_particles.as_slice(),
                                )
                            },
                        )
                    })
                    .map(|(g, sf)| (g.clone(), sf.clone()))
                    .collect::<Vec<_>>()
            });
            bar.finish_and_clear();

            step = Instant::now();
            info!(
                "{} | Δ={} | {:<95}{}",
                format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                    .blue()
                    .bold(),
                format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
                "Number of graphs after fast Cutkosky cut filter is applied:",
                format!("{}", processed_graphs.len()).green()
            );
            last_step = step;
        }

        // This secondary cutkosky cut filter is always necessary as the fast one can have false positives
        // even in the absence of constraints on amplitudes on either side of the cut
        if self.generation_type == GenerationType::CrossSection
            && !representative_final_pdgs.is_empty()
            && (self.amplitude_filters.get_coupling_orders().is_some()
                || self.amplitude_filters.get_loop_count_range().is_some()
                || !applied_fast_cutksosky_cut_filter)
        {
            // if self.generation_type == GenerationType::CrossSection
            //     && !self.final_pdgs.is_empty()
            // {
            let bar = ProgressBar::new(processed_graphs.len() as u64);
            bar.set_style(progress_bar_style.clone());
            bar.set_message("Applying secondary exact Cutkosky cut filter...");

            pool.install(|| {
                processed_graphs = processed_graphs
                    .par_iter()
                    .progress_with(bar.clone())
                    .filter(|(g, _)| {
                        self.half_edge_filters(
                            model,
                            g,
                            &external_connections,
                            n_unresolved,
                            &unresolved_type,
                        )
                    })
                    .map(|(g, sf)| (g.clone(), sf.clone()))
                    .collect::<Vec<_>>()
            });
            bar.finish_and_clear();

            step = Instant::now();
            info!(
                "{} | Δ={} | {:<95}{}",
                format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                    .blue()
                    .bold(),
                format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
                "Number of graphs after exact Cutkosky cut filter is applied:",
                format!("{}", processed_graphs.len()).green()
            );
            last_step = step;
        }

        // Because of the interplay with the cutkosky cut filter and left-right canonization when using symmetrize_left_right_states
        // we must do to canonicalizations here and we will select the "smallest one"
        let mut node_colors_for_canonicalization: HashMap<i32, i32> = HashMap::default();
        let mut perform_graph_pregrouping_without_numerator_and_left_right_symmetry = true;
        match self.generation_type {
            GenerationType::CrossSection => {
                match (
                    self.symmetrize_initial_states,
                    self.symmetrize_left_right_states,
                ) {
                    // (true, true) | (true, false) => {
                    //     for initial_color in 1..=self.initial_pdgs.len() {
                    //         node_colors_for_canonicalization.insert(initial_color as i32, -2);
                    //     }
                    //     for final_color in self.initial_pdgs.len() + 1
                    //         ..=2 * self.initial_pdgs.len()
                    //     {
                    //         node_colors_for_canonicalization.insert(final_color as i32, -3);
                    //     }
                    // }
                    (false, true) | (true, true) | (true, false) => {
                        // In this case we only care about left-righ  the pre-grouping does nothing so we can skip it
                        perform_graph_pregrouping_without_numerator_and_left_right_symmetry = false;
                        for initial_color in 1..=self.initial_pdgs.len() {
                            node_colors_for_canonicalization
                                .insert(initial_color as i32, -2000 - (initial_color as i32));
                        }
                        for final_color in self.initial_pdgs.len() + 1..=2 * self.initial_pdgs.len()
                        {
                            node_colors_for_canonicalization.insert(
                                final_color as i32,
                                -3000 - ((final_color - self.initial_pdgs.len()) as i32),
                            );
                        }
                    }
                    (false, false) => {}
                }
            }
            GenerationType::Amplitude => {
                match (
                    self.symmetrize_initial_states,
                    self.symmetrize_final_states,
                    self.symmetrize_left_right_states,
                ) {
                    (true, true, true) => {
                        for initial_color in 1..=self.initial_pdgs.len() {
                            if !model
                                .get_particle_from_pdg(
                                    self.initial_pdgs[initial_color - 1] as isize,
                                )
                                .0
                                .is_fermion()
                                || self.allow_symmetrization_of_external_fermions_in_amplitudes
                            {
                                node_colors_for_canonicalization.insert(initial_color as i32, -1);
                            }
                        }
                        for final_color in self.initial_pdgs.len() + 1
                            ..=self.initial_pdgs.len() + representative_final_pdgs.len()
                        {
                            if !model
                                .get_particle_from_pdg(
                                    representative_final_pdgs
                                        [final_color - self.initial_pdgs.len() - 1]
                                        as isize,
                                )
                                .0
                                .is_fermion()
                                || self.allow_symmetrization_of_external_fermions_in_amplitudes
                            {
                                node_colors_for_canonicalization.insert(final_color as i32, -1);
                            }
                        }
                    }
                    (true, false, false) => {
                        for initial_color in 1..=self.initial_pdgs.len() {
                            if !model
                                .get_particle_from_pdg(
                                    self.initial_pdgs[initial_color - 1] as isize,
                                )
                                .0
                                .is_fermion()
                                || self.allow_symmetrization_of_external_fermions_in_amplitudes
                            {
                                node_colors_for_canonicalization.insert(initial_color as i32, -2);
                            }
                        }
                    }
                    (false, true, false) => {
                        for final_color in self.initial_pdgs.len() + 1
                            ..=self.initial_pdgs.len() + representative_final_pdgs.len()
                        {
                            if !model
                                .get_particle_from_pdg(
                                    representative_final_pdgs
                                        [final_color - self.initial_pdgs.len() - 1]
                                        as isize,
                                )
                                .0
                                .is_fermion()
                                || self.allow_symmetrization_of_external_fermions_in_amplitudes
                            {
                                node_colors_for_canonicalization.insert(final_color as i32, -3);
                            }
                        }
                    }
                    (true, true, false) => {
                        for initial_color in 1..=self.initial_pdgs.len() {
                            if !model
                                .get_particle_from_pdg(
                                    self.initial_pdgs[initial_color - 1] as isize,
                                )
                                .0
                                .is_fermion()
                                || self.allow_symmetrization_of_external_fermions_in_amplitudes
                            {
                                node_colors_for_canonicalization.insert(initial_color as i32, -2);
                            }
                        }
                        for final_color in self.initial_pdgs.len() + 1
                            ..=self.initial_pdgs.len() + representative_final_pdgs.len()
                        {
                            if !model
                                .get_particle_from_pdg(
                                    representative_final_pdgs
                                        [final_color - self.initial_pdgs.len() - 1]
                                        as isize,
                                )
                                .0
                                .is_fermion()
                                || self.allow_symmetrization_of_external_fermions_in_amplitudes
                            {
                                node_colors_for_canonicalization.insert(final_color as i32, -3);
                            }
                        }
                    }
                    (false, false, false) => {
                        // No external symmetrization needed
                        perform_graph_pregrouping_without_numerator_and_left_right_symmetry = false;
                    }
                    _ => {
                        return Err(FeynGenError::GenericError(format!(
                            "Option symmetrize_initial_states={}, symmetrize_final_states={} and symmetrize_left_right_states={} not valid for amplitude generation.",
                            self.symmetrize_initial_states,
                            self.symmetrize_final_states,
                            self.symmetrize_left_right_states
                        )));
                    }
                }
            }
        }

        if perform_graph_pregrouping_without_numerator_and_left_right_symmetry
            && !node_colors_for_canonicalization.is_empty()
        {
            processed_graphs = group_isomorphic_graphs_after_node_color_change(
                &processed_graphs.iter().cloned().collect::<HashMap<_, _>>(),
                &node_colors_for_canonicalization,
                &pool,
                &progress_bar_style,
            )
            .iter()
            .map(|(g, m)| (g.clone(), m.clone()))
            .collect::<Vec<_>>();
        }

        step = Instant::now();
        info!(
            "{} | Δ={} | {:<95}{}",
            format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                .blue()
                .bold(),
            format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
            "Number of graphs after symmetrization of external states:",
            format!("{}", processed_graphs.len()).green()
        );
        last_step = step;

        let bar = ProgressBar::new(graphs.len() as u64);
        bar.set_style(progress_bar_style.clone());
        bar.set_message("Canonizing edge and nodes ordering of selected graphs, including leg symmetrization...");
        // println!(
        //     "Before edge and vertex canonization, first graph:\n{}",
        //     processed_graphs.first().unwrap().0.to_dot()
        // );
        let mut canonized_processed_graphs = pool.install(|| {
            processed_graphs
                .par_iter()
                .progress_with(bar.clone())
                .map(|(g, symmetry_factor)| {
                    let manually_canonalize_initial_states_cross_section_ordering =
                        if self.generation_type == GenerationType::CrossSection
                            && (self.symmetrize_initial_states || self.symmetrize_left_right_states)
                        {
                            Some((
                                self.symmetrize_initial_states,
                                self.symmetrize_left_right_states,
                            ))
                        } else {
                            None
                        };
                    let (mut canonical_repr, mut sorted_g) = self
                        .canonicalize_edge_and_vertex_ordering(
                            model,
                            g,
                            &node_colors_for_canonicalization,
                            &self.numerator_grouping,
                            manually_canonalize_initial_states_cross_section_ordering,
                        )
                        .unwrap();
                    // If we are symmetrizing the left-right states in the context of a cross-section, the canonaliztion above
                    // has canonized the choice of which nodes to assign to the initial and final state.
                    // We now do a second pass to canonalize the vertex ordering for that particular choice.
                    if manually_canonalize_initial_states_cross_section_ordering.is_some() {
                        (canonical_repr, sorted_g) = self
                            .canonicalize_edge_and_vertex_ordering(
                                model,
                                &sorted_g,
                                &node_colors_for_canonicalization,
                                &self.numerator_grouping,
                                None,
                            )
                            .unwrap();
                    }
                    // println!("NON CANONIZED FLOW: {}", sorted_g.to_dot());
                    // NOTE: The normalization of the fermion flow cannot be performed at this stage yet because
                    // it involves flipping edge orientation which modifies the sign of the momenta for the local
                    // numerator comparison during the grouping of graphs. This is done later in the process then.
                    let (g_with_canonical_flows, is_external_fermion_flow_sign_negative) = self
                        .normalize_flows(&sorted_g, model)
                        .expect("Failed to normalize fermion flow");
                    // println!("CANONIZED FLOW: {}", g_with_canonical_flows.to_dot());

                    let mut bare_graph = ParseGraph::from_symbolica_graph(
                        model,
                        "",
                        &sorted_g,
                        symmetry_factor.clone(),
                        &external_connections,
                    )?;

                    let fermion_sign = if self.generation_type == GenerationType::Amplitude {
                        if (!self.allow_symmetrization_of_external_fermions_in_amplitudes)
                            || (!self.symmetrize_initial_states && !self.symmetrize_final_states)
                        {
                            function!(
                                symbol!("ExternalFermionOrderingSign"),
                                Atom::num(if is_external_fermion_flow_sign_negative {
                                    -1
                                } else {
                                    1
                                })
                            )
                        } else {
                            Atom::num(1)
                        }
                    } else {
                        self.cross_section_external_fermion_ordering_sign(&mut bare_graph, model)?
                    };

                    bare_graph.global_data.overall_factor =
                        &bare_graph.global_data.overall_factor * fermion_sign;

                    Ok(CanonizedGraphInfo {
                        canonized_graph: canonical_repr,
                        graph: sorted_g,
                        graph_with_canonized_flow: g_with_canonical_flows,
                        gammaloop_graph: bare_graph,
                    })
                })
                .collect::<Result<Vec<_>>>()
        })?;
        bar.finish_and_clear();

        // Combine duplicates
        let mut combined_canonized_processed_graphs = HashMap::default();
        let numerator_independent_symmetry_pattern = function!(
            symbol!("NumeratorIndependentSymmetryGrouping"),
            Atom::var(symbol!("x_"))
        )
        .to_pattern();
        for canonized_graph in canonized_processed_graphs {
            let g_with_canonical_flows_clone = canonized_graph.graph_with_canonized_flow.clone();
            combined_canonized_processed_graphs
                .entry(g_with_canonical_flows_clone)
                .and_modify(|entry: &mut CanonizedGraphInfo| {
                    // NumeratorIndependentSymmetryGrouping
                    let ratio = (evaluate_overall_factor(
                        canonized_graph
                            .gammaloop_graph
                            .global_data
                            .overall_factor
                            .as_view(),
                    ) / evaluate_overall_factor(
                        entry
                            .gammaloop_graph
                            .global_data
                            .overall_factor
                            .replace(&numerator_independent_symmetry_pattern)
                            .with(Atom::num(1).to_pattern())
                            .as_view(),
                    ))
                    .expand();
                    if entry
                        .gammaloop_graph
                        .global_data
                        .overall_factor
                        .pattern_match(&numerator_independent_symmetry_pattern, None, None)
                        .next()
                        .is_some()
                    {
                        entry.gammaloop_graph.global_data.overall_factor = entry
                            .gammaloop_graph
                            .global_data
                            .overall_factor
                            .replace(&numerator_independent_symmetry_pattern)
                            .with(
                                function!(
                                    symbol!("NumeratorIndependentSymmetryGrouping"),
                                    (Atom::var(symbol!("x_")) + ratio).expand()
                                )
                                .to_pattern(),
                            );
                    } else {
                        entry.gammaloop_graph.global_data.overall_factor =
                            &entry.gammaloop_graph.global_data.overall_factor
                                * function!(
                                    symbol!("NumeratorIndependentSymmetryGrouping"),
                                    (Atom::num(1) + ratio).expand()
                                );
                    }
                })
                .or_insert(canonized_graph);
        }
        canonized_processed_graphs = combined_canonized_processed_graphs
            .values()
            .map(|v| v.to_owned())
            .collect::<Vec<_>>();

        // println!(
        //     "After edge and vertex canonization, first graph:\n{}",
        //     canonized_processed_graphs.first().unwrap().1.to_dot()
        // );

        // Make the graph order deterministic and already filter out the identical ones using the canonical representation with normalized fermion flow
        canonized_processed_graphs.sort_by(|a, b| {
            let ordering = (a.graph_with_canonized_flow)
                .partial_cmp(&b.graph_with_canonized_flow)
                .unwrap();
            if ordering.is_eq() {
                panic!(
                    "Two graphs are identical after canonization with normalized fermion flow. This should never happen."
                );
            }
            ordering
        });

        step = Instant::now();
        info!(
            "{} | Δ={} | {:<95}{}",
            format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                .blue()
                .bold(),
            format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
            "Number of graphs after canonization of vertices, edges and flows (graph labels assigned here):",
            format!("{}", canonized_processed_graphs.len()).green()
        );
        last_step = step;

        let n_zeroes_color = Arc::new(Mutex::new(0));
        let n_zeroes_lorentz = Arc::new(Mutex::new(0));
        let n_groupings = Arc::new(Mutex::new(0));
        // the pooled bare graphs have keys being the skeletton graphs identifying the topology
        #[allow(clippy::type_complexity)]
        let pooled_bare_graphs: Arc<
            Mutex<
                HashMap<
                    SymbolicaGraph<NodeColorWithoutVertexRule, std::string::String>,
                    Vec<Vec<PooledGraphData>>,
                >,
            >,
        > = Arc::new(Mutex::new(HashMap::default()));
        let pooled_bare_graphs_clone = pooled_bare_graphs.clone();

        let bar = ProgressBar::new(canonized_processed_graphs.len() as u64);
        bar.set_style(progress_bar_style.clone());
        bar.set_message(format!(
            "Final numerator-aware processing of remaining graphs ({} found: {} | {} found: {})...",
            "#zeros".green(),
            format!("{}", 0).green().bold(),
            "#groupings".green(),
            format!("{}", 0).green().bold(),
        ));

        // The quantity 'cuts_per_graph' below is only for debugging, remove its filling if too slow
        #[allow(unused_mut)]
        let mut cuts_per_graph: AHashMap<String, Vec<Vec<isize>>> = AHashMap::default();
        // This function 'half_edge_filters' no longer returns the number of active cuts (it used to returned all of them as Vec<Vec<isize>>)
        // A modified version without early termination would need to be used instead
        // canonized_processed_graphs
        //     .iter()
        //     .enumerate()
        //     .for_each(|(i_g, canonical_graph)| {
        //         let graph_name = format!("{}{}", self.graph_prefix, i_g);
        //         cuts_per_graph.insert(
        //             graph_name.clone(),
        //             self.half_edge_filters(
        //                 model,
        //                 &canonical_graph.graph,
        //                 &external_connections,
        //                 n_unresolved,
        //                 &unresolved_type,
        //             ),
        //         );
        //     });

        let samples: Vec<_> = if let Some(opts) = self.numerator_grouping.get_options() {
            let mut sample_iterator = PrimeIteratorU64::new(1);
            sample_iterator.nth(opts.numerical_sample_seed as usize);
            (0..opts.number_of_numerical_samples)
                .map(|_| {
                    self.sample_lib(
                        &mut sample_iterator,
                        opts.fully_numerical_substitution_when_comparing_numerators,
                        opts.symmetric_polarizations,
                        model,
                    )
                })
                .collect()
        } else {
            vec![]
        };

        let was_interrupted = Arc::new(AtomicBool::new(false));
        let graph_count = canonized_processed_graphs.len();
        let padding_width = if graph_count <= 1 {
            1
        } else {
            (graph_count - 1).to_string().len()
        };
        pool.install(|| {
            canonized_processed_graphs
                .into_par_iter()
                .progress_with(bar.clone())
                .enumerate()
                .map({
                    |(i_g, canonical_graph)| {
                        let was_interrupted = Arc::clone(&was_interrupted);
                        let graph_name: String = format!("{}{:0width$}", self.graph_prefix, i_g, width = padding_width);
                        if let Some(selected_graphs) = &self.selected_graphs
                            && !selected_graphs.contains(&graph_name) {
                                return Ok(())
                            }
                        if let Some(vetoed_graphs) = &self.vetoed_graphs
                            && vetoed_graphs.contains(&graph_name) {
                                return Ok(())
                            }
                        let mut bare_graph = Graph::from_parsed(canonical_graph.gammaloop_graph,model)?;
                        bare_graph.name = graph_name.clone();

                        // The step below is optional, but it is nice to have all internal fermion edges canonized as particles.
                        // Notice that we cannot do this on the bare graph used for numerator local comparisons and diagram grouping because
                        // it induces a misalignment of the LMB (w.r.t to their sign/orientation) due to the flip of the edges.
                        // In principle this could be fixed by forcing to pick an LMB for edges that have not been flipped and allowing closed loops
                        // of antiparticles as well, but I prefer to leave this a post-re-processing option instead.
                        let canonized_fermion_flow_bare_graph = if CANONIZE_GRAPH_FLOWS {
                            Graph::from_symbolica_graph(
                                model,
                                graph_name,
                                &canonical_graph.graph_with_canonized_flow,
                                bare_graph.overall_factor.clone(),
                                &external_connections,
                            )?
                        } else {
                            bare_graph.clone()
                        };

                        // println!(
                        //     "bare_graph:\n{}",
                        //     bare_graph.dot()
                        // );
                        // When disabling numerator-aware graph isomorphism, each graph is added separately

                        let pooled_graph = PooledGraphData {
                            graph_id: i_g,
                            numerator_data: None,
                            ratio: Atom::num(1),
                            bare_graph:canonized_fermion_flow_bare_graph.clone(),
                        };
                        if was_interrupted.load(Ordering::Relaxed) && is_interrupted(){
                            was_interrupted.store(true, Ordering::Relaxed);
                            eprintln!("Numerator-aware comparison of graphs interrupted by user, finishing current operations {}...","WITHOUT NUMERATOR COMPARISONS".red().bold());
                        }
                        if was_interrupted.load(Ordering::Relaxed) || matches!(
                            self.numerator_grouping,
                            NumeratorAwareGraphGroupingOption::NoGrouping
                        ) {
                            {
                                let mut pooled_bare_graphs_lock = pooled_bare_graphs_clone.lock().unwrap();
                                match pooled_bare_graphs_lock.entry(canonical_graph.canonized_graph.clone()) {
                                    Entry::Vacant(entry) => {
                                        entry.insert(vec![vec![pooled_graph]]);
                                    }
                                    Entry::Occupied(mut entry) => {
                                        entry.get_mut().push(vec![pooled_graph]);
                                    }
                                }
                            }
                        } else {
                            // println!("Processing graph #{}...", i_g);
                            // println!("Bare graph: {}",bare_graph.dot_serialize());
                            let mut numerator = bare_graph.numerator(&bare_graph.no_dummy(),&bare_graph.empty_subgraph());


                            // TODO Check if we include overall factor in main
                            numerator.state.expr *=&bare_graph.global_prefactor.num * &bare_graph.global_prefactor.projector;// * &bare_graph.overall_factor;
                            numerator.state.expr = numerator.state.expr.replace_multiple(&cpl_reps);

                            debug!(num = %numerator.state.expr.to_ordered_simple(),"Initial numerator",);

                            // println!("HEEEEy");
                            let numerator_color_simplified =
                                numerator.clone().get_single_atom().unwrap().to_param_color().simplify_color();

                            // println!("numerator_color_simplified=\n{}",numerator_color_simplified.to_plain_string());

                            // Color part need to be expanded to identify zeroes properly
                            if numerator_color_simplified.expand_color().iter().fold(Atom::Zero,|acc,(a, b)| a*b+acc).is_zero()
                            {
                                {
                                    let mut n_zeroes_color_value = n_zeroes_color.lock().unwrap();
                                    *n_zeroes_color_value += 1;
                                    let n_groupings_value = n_groupings.lock().unwrap();
                                    let n_zeroes_lorentz_value = n_zeroes_lorentz.lock().unwrap();
                                    bar.set_message(format!("Final numerator-aware processing of remaining graphs ({} found: {} | {} found: {})...",
                                        "#zeros".green(),
                                        format!("{}",*n_zeroes_color_value + *n_zeroes_lorentz_value).green().bold(),
                                        "#groupings".green(),
                                        format!("{}",n_groupings_value).green().bold(),
                                    ));
                                }
                                return Ok(())
                            }
                            if matches!(
                                self.numerator_grouping,
                                NumeratorAwareGraphGroupingOption::OnlyDetectZeroes
                            ) {
                                {
                                    let mut pooled_bare_graphs_lock = pooled_bare_graphs_clone.lock().unwrap();
                                    match pooled_bare_graphs_lock.entry(canonical_graph.canonized_graph.clone()) {
                                        Entry::Vacant(entry) => {
                                            entry.insert(vec![vec![pooled_graph]]);
                                        }
                                        Entry::Occupied(mut entry) => {
                                            entry.get_mut().push(vec![pooled_graph]);
                                        }
                                    }
                                }
                            } else {



                                // println!(
                                //     "I have:\n{}",
                                //     numerator_color_simplified
                                //         .clone()
                                //         .parse()
                                //         .contract(ContractionSettings::<Rational>::Normal)
                                //         .unwrap()
                                //         .state
                                //         .tensor
                                //         .iter_flat()
                                //         .map(|(id, d)| format!("{}: {}", id, d))
                                //         .collect::<Vec<_>>()
                                //         .join("\n")
                                // );
                                // Important: The numerator-aware grouping is done with the simplified color structure and *not* with the fermion flow canonized bare graph
                                let numerator_data = Some(
                                    ProcessedNumeratorForComparison::from_numerator_symbolic_expression(
                                        i_g,
                                        &bare_graph,
                                       numerator_color_simplified,
                                        &samples,
                                        settings,
                                        &self.numerator_grouping,
                                    )?,
                                );

                                // Test if Lorentz evaluations are zero
                                if !numerator_data.as_ref().unwrap().sample_evaluations.is_empty()
                                    && numerator_data.as_ref().unwrap().sample_evaluations_are_zero.iter().all(|&b| b) {
                                        {
                                            let n_zeroes_color_value = n_zeroes_color.lock().unwrap();
                                            let mut n_zeroes_lorentz_value = n_zeroes_lorentz.lock().unwrap();
                                            *n_zeroes_lorentz_value += 1;
                                            let n_groupings_value = n_groupings.lock().unwrap();
                                            bar.set_message(format!("Final numerator-aware processing of remaining graphs ({} found: {} | {} found: {})...",
                                                "#zeros".green(),
                                                format!("{}",*n_zeroes_color_value + *n_zeroes_lorentz_value).green().bold(),
                                                "#groupings".green(),
                                                format!("{}",n_groupings_value).green().bold(),
                                            ));
                                        }
                                        return Ok(())
                                    }

                                // println!("Skeletton G#{}:\n{}", i_g, canonical_repr.to_dot());
                                {
                                    let mut pooled_bare_graphs_lock = pooled_bare_graphs_clone.lock().unwrap();

                                    match pooled_bare_graphs_lock.entry(canonical_graph.canonized_graph.clone()) {
                                        Entry::Vacant(entry) => {
                                            entry.insert(vec![vec![
                                                PooledGraphData {
                                                    graph_id: i_g,
                                                    numerator_data,
                                                    ratio: Atom::num(1),
                                                    bare_graph: canonized_fermion_flow_bare_graph,
                                                }
                                                ]]);
                                        }
                                        Entry::Occupied(mut entry) => {

                                            let match_found = entry.get().iter().enumerate().find_map(|(i_entry, pooled_graphs_lists_for_this_topology)| {
                                                let reference_pooled_graph_data = &pooled_graphs_lists_for_this_topology[0];
                                                Self::compare_numerator_tensors(
                                                    &self.numerator_grouping,
                                                    numerator_data.as_ref().unwrap(),
                                                    reference_pooled_graph_data.numerator_data.as_ref().unwrap(),
                                                ).map(|ratio| {
                                                    (i_entry, PooledGraphData {
                                                        graph_id: i_g,
                                                        numerator_data: None,
                                                        ratio,
                                                        bare_graph: canonized_fermion_flow_bare_graph.clone(),
                                                    })
                                                })
                                            });
                                            if let Some((i_entry, new_entry)) = match_found {
                                                {
                                                    let n_zeroes_color_value = n_zeroes_color.lock().unwrap();
                                                    let n_zeroes_lorentz_value = n_zeroes_lorentz.lock().unwrap();
                                                    let mut n_groupings_value =
                                                        n_groupings.lock().unwrap();
                                                    *n_groupings_value += 1;
                                                    bar.set_message(format!("Final numerator-aware processing of remaining graphs ({} found: {} | {} found: {})...",
                                                        "#zeros".green(),
                                                        format!("{}",*n_zeroes_color_value+ *n_zeroes_lorentz_value).green().bold(),
                                                        "#groupings".green(),
                                                        format!("{}",n_groupings_value).green().bold(),
                                                    ));
                                                }
                                                entry.get_mut()[i_entry].push(new_entry);
                                            } else {
                                                entry.get_mut().push(vec![
                                                    PooledGraphData {
                                                        graph_id: i_g,
                                                        numerator_data,
                                                        ratio: Atom::num(1),
                                                        bare_graph: canonized_fermion_flow_bare_graph,
                                                    }
                                                ]);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        Ok(())
                    }
                }).collect::<Result<Vec<()>, FeynGenError>>()
        })?;
        // Reset the interrupt flag
        set_interrupted(false);
        bar.finish_and_clear();

        // Now combine the pooled graphs identified to be combined.
        let mut bare_graphs: Vec<(usize, Graph)> = Vec::default();
        let mut pooled_bare_graphs_len = 0;
        let mut n_cancellations: i32 = 0;
        for (_canonical_repr, pooled_graphs_lists_for_this_topology) in
            pooled_bare_graphs.lock().unwrap().iter()
        {
            pooled_bare_graphs_len += 1;
            for pooled_graphs_list in pooled_graphs_lists_for_this_topology {
                let sorted_graphs_to_combine = pooled_graphs_list
                    .iter()
                    .sorted_by_key(|pooled_graph| pooled_graph.graph_id)
                    .collect::<Vec<_>>();
                let new_reference_ratio = sorted_graphs_to_combine[0].ratio.clone();
                let mut combined_overall_factor = Atom::num(0);
                let mut bare_graph_representative = sorted_graphs_to_combine[0].bare_graph.clone();
                if sorted_graphs_to_combine.len() > 1 {
                    for graph_to_combine in sorted_graphs_to_combine {
                        // The computation of &graph_to_combine.ratio / &previous_reference_ratio is necessary so that if we started with this order of graphs to combine
                        //   (A, B, C), with ratios (r_1 = 1, r_2 = B/A, r_3 = C/A)
                        // And when forcing the reference diagram to be e.g. C, (so that we get a predictive ref graph in the multi-thread case), we need to pick C as the
                        // reference and have A and B be multiplied by the ratios A/C and B/C respectively. respectively.
                        // These will here be computed as A/C = r_1 / r_3 and B/C = r_2 / r_3. In general `r_i / r_ref` always yield the desired ratio.
                        combined_overall_factor += function!(
                            symbol!("NumeratorDependentGrouping"),
                            Atom::num(graph_to_combine.graph_id as i64),
                            (&graph_to_combine.ratio / &new_reference_ratio).expand(),
                            graph_to_combine.bare_graph.overall_factor.clone()
                        );
                    }
                } else {
                    combined_overall_factor = bare_graph_representative.overall_factor.clone();
                }
                if evaluate_overall_factor(combined_overall_factor.as_view())
                    .expand()
                    .is_zero()
                {
                    println!("{combined_overall_factor}");
                    n_cancellations += 1;
                } else {
                    bare_graph_representative.overall_factor = combined_overall_factor;
                    bare_graphs.push((pooled_graphs_list[0].graph_id, bare_graph_representative));
                }
            }
        }
        bare_graphs.sort_by(|a: &(usize, Graph), b| (a.0).cmp(&b.0));

        let (n_zeroes_color_value, n_zeroes_lorentz_value, n_groupings_value) = {
            (
                *n_zeroes_color.lock().unwrap(),
                *n_zeroes_lorentz.lock().unwrap(),
                *n_groupings.lock().unwrap(),
            )
        };
        step = Instant::now();
        info!(
            "{} | Δ={} | {:<95}{} ({} isomorphically unique graph{}, {} color zero{}, {} lorentz zero{}, {} grouped and {} cancellation{})",
            format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                .blue()
                .bold(),
            format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
            format!(
                "Number of graphs after numerator-aware grouping with strategy '{}':",
                self.numerator_grouping
            ),
            format!("{}", bare_graphs.len()).green().bold(),
            format!("{}", pooled_bare_graphs_len).blue().bold(),
            if pooled_bare_graphs_len > 1 { "s" } else { "" },
            format!("{}", n_zeroes_color_value).blue().bold(),
            if n_zeroes_color_value > 1 { "s" } else { "" },
            format!("{}", n_zeroes_lorentz_value).blue().bold(),
            if n_zeroes_lorentz_value > 1 { "s" } else { "" },
            format!("{}", n_groupings_value).blue().bold(),
            format!("{}", n_cancellations).blue().bold(),
            if n_cancellations > 1 { "s" } else { "" },
        );

        if let Some(lmbs) = self.loop_momentum_bases.as_ref() {
            for (_graph_id, graph) in bare_graphs.iter_mut() {
                if let Some(lmb) = lmbs.get(&graph.name) {
                    let full_filter = graph.full_filter();
                    let mut cut_graph = full_filter.subtract(&graph.initial_state_cut.right);

                    for e in lmb.iter() {
                        cut_graph.sub(graph[&EdgeIndex(*e)].1);
                    }

                    let mut loop_momentum_basis = if let Some(i) = cut_graph.included_iter().next()
                    {
                        let tree = SimpleTraversalTree::depth_first_traverse(
                            graph,
                            &cut_graph,
                            &graph.node_id(i),
                            None,
                        )
                        .unwrap();

                        let external = graph.internal_crown(&full_filter);
                        graph.lmb_impl(&full_filter, &tree.tree_subgraph, external)
                    } else {
                        return Err(FeynGenError::Eyre(eyre!(
                            "No included edges found in full_cut for loop momentum basis setup"
                        )));
                    };

                    for e in graph.initial_state_cut.left.included_iter() {
                        let e = graph[&e];
                        let (l, _) = loop_momentum_basis
                            .loop_edges
                            .iter()
                            .find_position(|a| *a == &e)
                            .unwrap();

                        loop_momentum_basis.put_loop_to_ext(LoopIndex(l));
                    }

                    let lmb_init_loop_ids = lmb
                        .iter()
                        .map(|e| {
                            LoopIndex(
                                loop_momentum_basis
                                    .loop_edges
                                    .iter()
                                    .find_position(|a| a.0 == *e)
                                    .unwrap()
                                    .0,
                            )
                        })
                        .collect::<Vec<_>>();

                    let p = Permutation::sort(&lmb_init_loop_ids)
                        .inverse()
                        .transpositions();

                    for (a, b) in p {
                        loop_momentum_basis.swap_loops(LoopIndex(a), LoopIndex(b));
                    }
                    graph.loop_momentum_basis = loop_momentum_basis;
                }
            }
        }
        debug!(
            "Graphs: [\n{}\n]",
            bare_graphs
                .iter()
                .map(|(_graph_id, graph)| format!(
                    "{:-6} @ {} = {{{}}}",
                    graph.name.clone(),
                    evaluate_overall_factor(graph.overall_factor.as_view())
                        .expand()
                        .to_canonical_string(),
                    graph.overall_factor
                ))
                .collect::<Vec<_>>()
                .join("\n"),
        );
        let mut total_sym_factor = Atom::num(0);
        for (_i_g, g) in bare_graphs.iter() {
            total_sym_factor += evaluate_overall_factor(g.overall_factor.as_view());
        }
        debug!(
            "Graphs: [\n{}\n]",
            bare_graphs
                .iter()
                .map(|(_graph_id, graph)| {
                    let n_cuts_str = if let Some(cuts) = cuts_per_graph.get(&graph.name.to_string())
                    {
                        format!(" ({} cuts)", cuts.len())
                    } else {
                        String::from("")
                    };
                    format!(
                        "{:-6} @ {} = {{{}}}{}",
                        graph.name.clone(),
                        evaluate_overall_factor(graph.overall_factor.as_view())
                            .expand()
                            .to_canonical_string(),
                        graph.overall_factor,
                        n_cuts_str
                    )
                })
                .collect::<Vec<_>>()
                .join("\n"),
        );
        info!(
            "( Sum of the symmetry factors from each graph generated = {} ) ",
            format!("{}", total_sym_factor).green()
        );

        Ok(bare_graphs
            .iter()
            .map(|(_graph_id, graph)| graph.clone())
            .collect::<Vec<_>>())
    }

    #[instrument(skip_all)]
    fn compare_numerator_tensors(
        numerator_aware_isomorphism_grouping: &NumeratorAwareGraphGroupingOption,
        numerator_a: &ProcessedNumeratorForComparison,
        numerator_b: &ProcessedNumeratorForComparison,
    ) -> Option<Atom> {
        if numerator_a.canonized_numerator.is_none() && numerator_b.canonized_numerator.is_some()
            || numerator_a.canonized_numerator.is_some()
                && numerator_b.canonized_numerator.is_none()
        {
            panic!(
                "Inconsistent state: one sample has canonalized numerator while the other does not."
            );
        }

        if numerator_a.sample_evaluations.len() != numerator_b.sample_evaluations.len() {
            panic!(
                "Inconsistent state: the two samples have different number of numerical evaluations."
            );
        }

        match numerator_aware_isomorphism_grouping {
            NumeratorAwareGraphGroupingOption::NoGrouping => None,
            NumeratorAwareGraphGroupingOption::OnlyDetectZeroes => None,
            NumeratorAwareGraphGroupingOption::GroupIdenticalGraphUpToSign(_) => {
                numerator_a.compare_with_sign_only(numerator_b)
            }
            NumeratorAwareGraphGroupingOption::GroupIdenticalGraphUpToScalarRescaling(_) => {
                numerator_a.compare_with_scalar_rescaling(numerator_b)
            }
        }
    }

    pub(crate) fn substitute_color_factors(expr: AtomView) -> Atom {
        let replacements = vec![
            (Atom::var(CS.nc), Atom::num(3)),
            (Atom::var(CS.tr), parse!("1/2")),
            // (parse!("CA"), parse!("3")),
            // (parse!("CF"), parse!("4/3")),
        ];
        let mut res = expr.to_owned();
        for (src, trgt) in replacements {
            res = res.replace(src.to_pattern()).with(trgt.to_pattern());
        }
        res
    }

    fn cross_section_external_fermion_ordering_sign(
        &self,
        graph: &mut ParseGraph,
        model: &Model,
    ) -> Result<Atom> {
        let n_external_fermion_loops = graph.n_external_fermion_loops()?;

        let number_of_initial_antifermions = self
            .initial_pdgs
            .iter()
            .filter(|&pdg| {
                let p = model.get_particle_from_pdg(*pdg as isize);

                p.0.is_antiparticle() && p.0.is_fermion()
            })
            .count();

        let sign = Sign::Negative.pow(n_external_fermion_loops);
        let antifermion_spinsum_sign = Sign::Negative.pow(number_of_initial_antifermions);

        Ok(function!(
            symbol!("ExternalFermionOrderingSign"),
            Atom::num(match sign {
                Sign::Positive => 1,
                Sign::Negative => -1,
            })
        ) * function!(
            symbol!("AntiFermionSpinSumSign"),
            Atom::num(match antifermion_spinsum_sign {
                Sign::Positive => 1,
                Sign::Negative => -1,
            })
        ))
    }
}

struct PooledGraphData {
    graph_id: usize,
    numerator_data: Option<ProcessedNumeratorForComparison>,
    ratio: Atom,
    bare_graph: Graph,
}

#[derive(Clone)]
pub(crate) struct ProcessedNumeratorForComparison {
    diagram_id: usize,
    /// Canonized numerator expression for symbolic comparison.
    ///
    /// This is `None` when:
    /// - `--compare-canonized-numerator` flag is disabled (test_canonized_numerator = false)
    /// - Numerator grouping is disabled (NoGrouping or OnlyDetectZeroes variants)
    /// - No grouping options are available (numerator_aware_isomorphism_grouping.get_options() returns None)
    ///
    /// This is `Some(canonized_atom)` when:
    /// - `--compare-canonized-numerator` flag is enabled (test_canonized_numerator = true)
    /// - The numerator has been processed through loop momentum replacements
    /// - Optionally processed through color factor substitution (if --fully-numerical-substitution-when-comparing-numerators)
    /// - Optionally processed through factor collection (for GroupIdenticalGraphUpToScalarRescaling mode)
    canonized_numerator: Option<Atom>,

    /// Numerical evaluations of the numerator at sample points for numerical comparison.
    ///
    /// This is empty (`vec![]`) when:
    /// - `--number-of-samples-for-numerator-comparisons 0` is specified
    /// - Numerator grouping is disabled (NoGrouping or OnlyDetectZeroes variants)
    /// - No grouping options are available (numerator_aware_isomorphism_grouping.get_options() returns None)
    /// - Neither test_canonized_numerator nor number_of_numerical_samples > 0
    ///
    /// This contains evaluated `Atom`s when:
    /// - `--number-of-samples-for-numerator-comparisons N` where N > 0
    /// - Each element represents the numerator evaluated at a different sample point
    /// - Sample points are generated using a prime-based iterator seeded by numerical_sample_seed
    /// - Evaluations include color expansion, tensor network execution, and color factor substitution
    /// - May have common factors collected depending on the grouping mode and ANALYZE_RATIO_AS_RATIONAL_POLYNOMIAL setting
    sample_evaluations: Vec<Atom>,
    sample_evaluations_as_polynomial: Vec<
        symbolica::poly::polynomial::MultivariatePolynomial<
            AlgebraicExtension<
                symbolica::domains::rational::FractionField<
                    symbolica::domains::integer::IntegerRing,
                >,
            >,
        >,
    >,
    sample_evaluations_are_zero: Vec<bool>,
}

impl ProcessedNumeratorForComparison {
    /// Compare two numerators allowing only sign changes (±1 ratios).
    #[instrument(skip_all, fields(self_id = %self.diagram_id, other_id = %other.diagram_id))]
    fn compare_with_sign_only(&self, other: &ProcessedNumeratorForComparison) -> Option<Atom> {
        debug!(
            self_diagram_id = %self.diagram_id,
            other_diagram_id = %other.diagram_id,
            "Starting sign-only comparison between diagrams"
        );
        fn analyze_diff_and_sum(a: AtomView, b: AtomView) -> Option<Atom> {
            if (a - b).expand().is_zero() {
                return Some(Atom::num(1));
            }
            if (a + b).expand().is_zero() {
                return Some(Atom::num(-1));
            }
            debug!(a = %a.floatify(13).to_canonical_string(),b=%b.floatify(13).to_canonical_string(),"compared but no luck");
            None
        }

        // Try canonized numerator comparison first
        if let (Some(canonized_num_a), Some(canonized_num_b)) = (
            self.canonized_numerator.as_ref(),
            other.canonized_numerator.as_ref(),
        ) {
            debug!("Attempting symbolic comparison using canonized numerators");
            debug!(
                numerator = %canonized_num_a.to_canonical_string(),
                numerator_diagram_id = %self.diagram_id,
                "Canonized numerator for comparison"
            );
            debug!(
                denominator = %canonized_num_b.to_canonical_string(),
                denominator_diagram_id = %other.diagram_id,
                "Canonized denominator for comparison"
            );
            if let Some(ratio) =
                analyze_diff_and_sum(canonized_num_a.as_view(), canonized_num_b.as_view())
            {
                debug!(
                    ratio = %ratio.to_canonical_string(),
                    method = "canonized_numerators",
                    "Successfully matched diagrams using canonized numerators"
                );
                return Some(ratio);
            } else {
                debug!(
                    comparison_type = "canonized_numerator",
                    result = "failed",
                    reason = "expressions_not_identical_up_to_sign",
                    "Canonized numerator comparison failed"
                );
            }
        } else {
            debug!(
                comparison_type = "canonized_numerator",
                result = "skipped",
                reason = "missing_numerators",
                "Skipping canonized numerator comparison"
            );
        }

        // Fall back to sample evaluations
        if !self.sample_evaluations.is_empty() {
            debug!(
                comparison_type = "numerical_samples",
                sample_count = %self.sample_evaluations.len(),
                "Attempting numerical comparison using sample evaluations"
            );
            let ratios = self
                .sample_evaluations
                .iter()
                .zip(other.sample_evaluations.iter())
                .enumerate()
                .map(|(idx, (a, b))| {
                    debug!(
                        sample_idx = %idx,
                        numerator = %a.to_canonical_string(),
                        numerator_diagram_id = %self.diagram_id,
                        "Sample numerator evaluation"
                    );
                    debug!(
                        sample_idx = %idx,
                        denominator = %b.to_canonical_string(),
                        denominator_diagram_id = %other.diagram_id,
                        "Sample denominator evaluation"
                    );
                    let ratio = analyze_diff_and_sum(a.as_view(), b.as_view());
                    if let Some(a) = &ratio {
                        debug!(
                            sample_idx = %idx,
                            ratio = %a.to_canonical_string(),
                            "Sign-only comparison result for sample"
                        );
                    } else {
                        debug!(
                            sample_idx = %idx,
                            result = "none_ratio",
                            reason = "expressions_not_identical_up_to_sign",
                            "Sign-only comparison result for sample is None"
                        );
                    }

                    ratio
                })
                .collect::<HashSet<_>>();

            debug!(
                unique_ratios_found = %ratios.len(),
                "Found unique ratios from sample evaluations"
            );

            if ratios.len() == 1 {
                if let Some(ratio) = ratios.iter().next().unwrap().to_owned() {
                    debug!(
                        ratio = %ratio.to_canonical_string(),
                        method = "numerical_evaluation",
                        "Successfully matched diagrams using numerical evaluation"
                    );
                    return Some(ratio);
                } else {
                    debug!(
                        result = "none_ratio",
                        reason = "likely_zeros",
                        "Sample evaluations yielded None ratio"
                    );
                }
            } else {
                debug!(
                    result = "inconsistent_ratios",
                    unique_ratios_count = %ratios.len(),
                    "Sample evaluations yielded inconsistent ratios - cannot group diagrams"
                );
            }
        } else {
            debug!(
                comparison_type = "numerical_samples",
                result = "skipped",
                reason = "no_sample_evaluations",
                "Skipping numerical comparison"
            );
        }

        debug!(
            self_diagram_id = %self.diagram_id,
            other_diagram_id = %other.diagram_id,
            final_result = "no_grouping_possible",
            "No grouping possible between diagrams"
        );
        None
    }

    /// Compare two numerators allowing arbitrary scalar rescaling.
    #[instrument(skip_all, fields(self_id = %self.diagram_id, other_id = %other.diagram_id))]
    pub(crate) fn compare_with_scalar_rescaling(
        &self,
        other: &ProcessedNumeratorForComparison,
    ) -> Option<Atom> {
        debug!(
            self_diagram_id = %self.diagram_id,
            other_diagram_id = %other.diagram_id,
            "Starting scalar rescaling comparison between diagrams"
        );
        fn analyze_ratios(ratios: &HashSet<Option<Atom>>) -> Option<Atom> {
            if ratios.len() > 1 {
                None
            } else {
                let ratio = ratios.iter().next().unwrap().to_owned();
                if let Some(r) = ratio {
                    for head in LibraryRep::all_self_duals()
                        .chain(LibraryRep::all_inline_metrics())
                        .chain(LibraryRep::all_dualizables())
                        .map(|a| a.to_symbolic([W_.a__]).to_pattern())
                    {
                        if r.pattern_match(&head, None, None).next().is_some() {
                            return None;
                        }
                    }
                    Some(r)
                } else {
                    None
                }
            }
        }

        // Try canonized numerator comparison first
        if let (Some(canonized_num_a), Some(canonized_num_b)) = (
            self.canonized_numerator.as_ref(),
            other.canonized_numerator.as_ref(),
        ) {
            debug!(
                numerator = %canonized_num_a.to_canonical_string(),
                numerator_diagram_id = %self.diagram_id,
                denominator = %canonized_num_b.to_canonical_string(),
                denominator_diagram_id = %other.diagram_id,
                "Attempting symbolic comparison using canonized numerators",
            );
            let mut ratios = HashSet::<Option<Atom>>::default();
            let r = if canonized_num_a == canonized_num_b {
                debug!("Canonized numerators are identical");
                Some(Atom::num(1))
            } else if *canonized_num_a == canonized_num_b * Atom::num(-1) {
                debug!("Canonized numerators differ by sign only");
                Some(Atom::num(-1))
            } else if canonized_num_b.is_zero() {
                debug!("Reference canonized numerator is zero - cannot compute ratio");
                None
            } else {
                let ratio = canonized_num_a / canonized_num_b;
                debug!(
                    computed_ratio = %ratio.to_canonical_string(),
                    "Computed canonized numerator ratio"
                );
                Some(ratio)
            };
            ratios.insert(r);
            if let Some(ratio) = analyze_ratios(&ratios) {
                debug!(
                    ratio = %ratio.to_canonical_string(),
                    method = "canonized_numerators",
                    "Successfully matched diagrams using canonized numerators"
                );
                return Some(ratio);
            } else {
                debug!(
                    comparison_type = "canonized_numerator",
                    result = "rejected",
                    reason = "problematic_patterns",
                    rejection_stage = "analyze_ratios",
                    "Canonized numerator ratio was rejected"
                );
            }
        } else {
            debug!(
                comparison_type = "canonized_numerator",
                result = "skipped",
                reason = "missing_numerators",
                "Skipping canonized numerator comparison"
            );
        }

        // Fall back to sample evaluations
        if !self.sample_evaluations.is_empty() {
            debug!(
                comparison_type = "numerical_samples",
                sample_count = %self.sample_evaluations.len(),
                "Attempting numerical comparison using sample evaluations"
            );
            let evaluations_a = &self.sample_evaluations;
            let evaluations_b = &other.sample_evaluations;
            if evaluations_a.is_empty() {
                debug!(
                    result = "cannot_proceed",
                    reason = "empty_self_evaluations",
                    "Self sample evaluations are empty"
                );
                return None;
            }

            let ratios = evaluations_a
                .iter()
                .zip(evaluations_b.iter())
                .enumerate()
                .map(|(idx, (a, b))| {
                    debug!(
                        sample_id= %idx,
                        numerator = %a.to_canonical_string(),
                        numerator_diagram_id = %self.diagram_id,
                        denominator = %b.to_canonical_string(),
                        denominator_diagram_id = %other.diagram_id,
                        "Sample evaluation A"
                    );
                    if a == b {
                        Some(Atom::num(1))
                    } else if *a == b * Atom::num(-1) {
                        Some(Atom::num(-1))
                    } else if b.is_zero() || a.is_zero() {
                        debug!(
                            sample_idx = %idx,
                            "Skipping sample due to zero value"
                        );
                        None
                    } else {
                        let ratio = if ANALYZE_RATIO_AS_RATIONAL_POLYNOMIAL {
                            // let a_poly = a.to_polynomial(&Q_I.clone(), None);
                            // let b_poly = b.to_polynomial(&Q_I.clone(), None);
                            let a_poly = &self.sample_evaluations_as_polynomial[idx];
                            let b_poly = &other.sample_evaluations_as_polynomial[idx];
                            if a_poly.is_zero() || b_poly.is_zero() {
                                debug!(
                                    sample_idx = %idx,
                                    "Skipping sample due to zero value after expansion"
                                );
                                None
                            } else {
                                let element = COMPLEXRATPOLYFIELD.to_element(a_poly.clone(), b_poly.clone(), true);
                                Some(polyrat_to_atom(&element))
                            }
                            // let element = COMPLEXRATPOLYFIELD.to_element(a.to_polynomial(&Q_I, None), b.to_polynomial(&Q_I, None), true);
                            // Some(polyrat_to_atom(&element))
                        } else {
                            Some((a / b).cancel())
                        };
                        debug!(
                            sample_idx = %idx,
                            computed_ratio = %ratio.clone().map(|r| r.to_ordered_simple().to_string()).unwrap_or("None".into()),
                            "Computed ratio for sample"
                        );
                        ratio
                    }
                })
                .collect::<HashSet<_>>();

            debug!(
                unique_ratios_found = %ratios.len(),
                "Found unique ratios from sample evaluations"
            );
            if event_enabled!(tracing::Level::DEBUG, parent: &Span::current()) {
                for ((rat, a), b) in ratios.iter().zip(evaluations_a).zip(evaluations_b) {
                    debug!(
                        self_diagram_id = %self.diagram_id,
                        other_diagram_id = %other.diagram_id,
                        ratio_value = ?rat.as_ref().map(|ra| ra.floatify(13).to_ordered_simple()).unwrap_or("None".into()),
                        numerator_value = %a.floatify(13).to_ordered_simple(),
                        denominator_value = %b.floatify(13).to_ordered_simple(),
                        "Detailed sample evaluation ratio information"
                    );
                }
            }
            if let Some(ratio) = analyze_ratios(&ratios) {
                debug!(
                    ratio = %ratio.to_canonical_string(),
                    method = "numerical_evaluation",
                    "Successfully matched diagrams using numerical evaluation"
                );
                return Some(ratio);
            } else {
                debug!(
                    comparison_type = "numerical_samples",
                    result = "rejected",
                    reason = "inconsistent_or_problematic_patterns",
                    rejection_stage = "analyze_ratios",
                    "Sample evaluation ratios were rejected"
                );
            }
        } else {
            debug!(
                comparison_type = "numerical_samples",
                result = "skipped",
                reason = "no_sample_evaluations",
                "Skipping numerical comparison"
            );
        }

        debug!(
            self_diagram_id = %self.diagram_id,
            other_diagram_id = %other.diagram_id,
            final_result = "no_grouping_possible",
            "No grouping possible between diagrams"
        );
        None
    }
    #[instrument(skip_all,fields(graph_name = %graph.name))]
    pub fn from_numerator_symbolic_expression(
        diagram_id: usize,
        graph: &Graph,
        mut numerator: Atom,
        samples: &[NumeratorSample],
        settings: &GlobalSettings,
        numerator_aware_isomorphism_grouping: &NumeratorAwareGraphGroupingOption,
    ) -> Result<Self, FeynGenError> {
        let default_processed_data = ProcessedNumeratorForComparison {
            diagram_id,
            canonized_numerator: None,
            sample_evaluations: vec![],
            sample_evaluations_as_polynomial: vec![],
            sample_evaluations_are_zero: vec![],
        };
        let res = if let Some(group_options) = numerator_aware_isomorphism_grouping.get_options() {
            if group_options.test_canonized_numerator
                || group_options.number_of_numerical_samples > 0
            {
                let lmb_reps = graph.integrand_replacement(
                    &graph.full_filter(),
                    &graph.loop_momentum_basis,
                    &[W_.x___],
                );

                numerator = numerator.replace_multiple(&lmb_reps);

                debug!(numerator=%numerator.to_ordered_simple(),diagram_id=%diagram_id,debug_dot=%graph.debug_dot(),"Initial Numerator");

                let canonized_numerator = if group_options.test_canonized_numerator {
                    let mut canonized_numerator_to_consider = numerator.canonize(Aind::Dummy);
                    if group_options.fully_numerical_substitution_when_comparing_numerators {
                        canonized_numerator_to_consider =
                            ProcessDefinition::substitute_color_factors(
                                canonized_numerator_to_consider.as_atom_view(),
                            )
                    };
                    // IMPORTANT: we must make sure to collect all common coefficients first
                    // with `collect_factors()` to ensure that common factors get simplified when looking at ratios.
                    if matches!(
                        numerator_aware_isomorphism_grouping,
                        NumeratorAwareGraphGroupingOption::GroupIdenticalGraphUpToScalarRescaling(
                            _
                        )
                    ) {
                        canonized_numerator_to_consider =
                            canonized_numerator_to_consider.collect_factors();
                    }
                    // numerator = canonized_numerator_to_consider.clone();
                    Some(canonized_numerator_to_consider)
                } else {
                    None
                };

                let mut numerators = vec![numerator];

                if settings
                    .generation
                    .feyngen
                    .gamma_simplification_closure_check
                {
                    debug!(numerator = %numerators[0].to_ordered_simple(),"Gamma Simplifying");
                    numerators.push(numerators[0].simplify_gamma());
                    debug!("Done Simplifying");
                }

                let mut sample_evals = VecDeque::new();
                for numerator in numerators {
                    let expanded = numerator.expand_color();

                    let sample_evaluations = samples
                    .iter()
                    .map(|(reps, lib)| {
                        let mut sample_evaluation = expanded
                            .iter()
                            .map(|(c, l)| {
                                debug!("Sample evaluation inputs c:{c},l:{l}");
                                let mut net = ParamParsingNet::try_from_view(l.as_view(), lib,&ParseSettings::default()).unwrap();
                                net.store
                                    .scalar
                                    .iter_mut()
                                    .for_each(|a| *a = a.replace_multiple(reps));

                                // debug!(net=?net.dot_pretty());
                                net.execute::<Sequential, SmallestDegree, _, _,_>(lib,PARAM_FUN_LIB.deref())
                                    .unwrap_or_else(|_| panic!("failed to execute net:{}", net.dot_pretty()));

                                // let c = ProcessDefinition::substitute_color_factors(c.as_view())
                                //     .expand();

                                let scalar = net
                                    .result_scalar();

                                if scalar.is_err(){
                                    debug!("Failed scalar for {} for graph:{}",net.dot_pretty(),graph.debug_dot());
                                }

                                let scalar:Atom = scalar
                                    .unwrap_or_else(|_| panic!("Expected scalar for c:{c} l:{l} that yields net:{} for graph {}",
                                        net.dot_pretty()
                                        ,graph.debug_dot()))
                                    .into();

                                // println!("Trying to canonize:{c}");
                                let canonized_color = c.canonize::<Aind>(Aind::Dummy);
                                debug!("canonizing \n{c}\n gives\n{canonized_color}");
                                let a = ProcessDefinition::substitute_color_factors(
                                    (canonized_color * scalar).as_view(),
                                );

                                debug!(evaluated=%a.printer(LOGPRINTOPTS),"evaluated{}",a.floatify(13).printer(LOGPRINTOPTS));
                                a
                            })
                            .fold(Atom::Zero, |acc, l| acc + l);

                        if EXPAND_NUMERICAL_SAMPLES_BEFORE_COMPARISON {
                            sample_evaluation = sample_evaluation.expand();
                        }
                        // TODO: optimize the above by instead directly storing "a.to_polynomial(&Q_I.clone(), None)" in the sample record, and not the expanded atom.
                        // Only do that for the if-branch below though.
                        if ANALYZE_RATIO_AS_RATIONAL_POLYNOMIAL || !matches!(numerator_aware_isomorphism_grouping,NumeratorAwareGraphGroupingOption::GroupIdenticalGraphUpToScalarRescaling(_))
                        {
                            sample_evaluation
                        } else {
                            // When not looking at the ratio of samples as a rational polynomial, we must make sure to collect all common coefficients first
                            // with `collect_factors()` to ensure that common factors get simplified when looking at ratios.
                            sample_evaluation.collect_factors()
                        }
                    })
                    .collect();

                    sample_evals.push_back(sample_evaluations);
                }

                let sample_evaluations = sample_evals.pop_front().unwrap();
                let sample_gamma_simplified_evaluations =
                    sample_evals.pop_front().unwrap_or(vec![]);

                if settings
                    .generation
                    .feyngen
                    .gamma_simplification_closure_check
                {
                    for (i, (a, b)) in sample_evaluations
                        .iter()
                        .zip(sample_gamma_simplified_evaluations.iter())
                        .enumerate()
                    {
                        if !(a - b).expand().is_zero() {
                            return Err(FeynGenError::Eyre(eyre!(
                                "Gamma simplification closure check failed for diagram ID {} on sample evaluation #{}. Numerator evaluation before gamma simplification: {}. After gamma simplification: {}.",
                                diagram_id,
                                i,
                                a.to_ordered_simple(),
                                b.to_ordered_simple()
                            )));
                        }
                    }
                }
                let samples_evals_as_polynomial: SampleEvaluationsAsPolynomial =
                    if ANALYZE_RATIO_AS_RATIONAL_POLYNOMIAL
                        || !matches!(
                        numerator_aware_isomorphism_grouping,
                        NumeratorAwareGraphGroupingOption::GroupIdenticalGraphUpToScalarRescaling(
                            _
                        )
                    ) {
                        let se_as_poly = sample_evaluations
                            .iter()
                            .map(|a| a.to_polynomial(&Q_I.clone(), None))
                            .collect::<Vec<_>>();
                        let se_are_zero =
                            se_as_poly.iter().map(|p| p.is_zero()).collect::<Vec<_>>();
                        (se_as_poly, se_are_zero)
                    } else {
                        (
                            vec![],
                            sample_evaluations
                                .iter()
                                .map(|a| a.expand().is_zero())
                                .collect::<Vec<_>>(),
                        )
                    };
                ProcessedNumeratorForComparison {
                    diagram_id,
                    canonized_numerator,
                    sample_evaluations,
                    sample_evaluations_as_polynomial: samples_evals_as_polynomial.0,
                    sample_evaluations_are_zero: samples_evals_as_polynomial.1,
                }
            } else {
                default_processed_data
            }
        } else {
            default_processed_data
        };
        Ok(res)
    }
}

// #[test]
// pub(crate) fn symbolica_symm_factors_bug() {
//     let external_edges_for_generation: Vec<(usize, (Option<bool>, usize))> = vec![];

//     let vertex_signatures_for_generation_a = vec![
//         vec![(None, 21), (None, 21), (None, 21)],
//         vec![(None, 21), (None, 21), (None, 21), (None, 21)],
//         vec![(Some(true), 1), (Some(false), 1), (None, 21)],
//         vec![(Some(true), 6), (Some(false), 6), (None, 21)],
//         vec![(Some(false), 9000005), (Some(true), 9000005), (None, 21)],
//         vec![(Some(true), 2), (Some(false), 2), (None, 21)],
//     ];

//     let settings = GenerationSettings::new()
//         .max_loops(5)
//         .max_bridges(0)
//         .allow_self_loops(true);

//     //let mut graphs_a = SymbolicaGraph::generate(
//     //    &external_edges_for_generation,
//     //    vertex_signatures_for_generation_a.as_slice(),
//     //    None,
//     //    Some(5),
//     //    Some(0),
//     //    true,
//     //);

//     let mut graphs_a = SymbolicaGraph::generate(
//         &external_edges_for_generation,
//         vertex_signatures_for_generation_a.as_slice(),
//         &settings,
//     )
//     .unwrap();

//     graphs_a.retain(|g, _| g.num_loops() >= 5);

//     let mut tot_symm_fact_graphs_a = Atom::num(0);
//     for (_g, s) in graphs_a.iter() {
//         tot_symm_fact_graphs_a = tot_symm_fact_graphs_a + Atom::num(1) / Atom::num(s.clone());
//     }
//     println!("tot_symm_fact_graphs_A = {}", tot_symm_fact_graphs_a);

//     let vertex_signatures_for_generation_b = vec![
//         vec![(Some(true), 1), (Some(false), 1), (None, 21)],
//         vec![(None, 21), (None, 21), (None, 21), (None, 21)],
//         vec![(Some(true), 6), (Some(false), 6), (None, 21)],
//         vec![(None, 21), (None, 21), (None, 21)],
//         vec![(Some(false), 9000005), (Some(true), 9000005), (None, 21)],
//         vec![(Some(true), 2), (Some(false), 2), (None, 21)],
//     ];
//     let mut graphs_b = SymbolicaGraph::generate(
//         &external_edges_for_generation,
//         vertex_signatures_for_generation_b.as_slice(),
//         &settings,
//     )
//     .unwrap();

//     graphs_b.retain(|g, _| g.num_loops() >= 5);

//     let mut tot_symm_fact_graphs_b = Atom::num(0);
//     for (_g, s) in graphs_b.iter() {
//         tot_symm_fact_graphs_b = tot_symm_fact_graphs_b + Atom::num(1) / Atom::num(s.clone());
//     }
//     println!("tot_symm_fact_graphs_B = {}", tot_symm_fact_graphs_b);

//     assert!(tot_symm_fact_graphs_a == tot_symm_fact_graphs_b);
// }
