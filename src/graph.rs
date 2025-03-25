use crate::{
    cff::{
        esurface::{
            generate_esurface_data, get_existing_esurfaces, EsurfaceDerivedData, ExistingEsurfaces,
            ExternalShift,
        },
        expression::CFFExpression,
        generation::generate_cff_expression,
    },
    feyngen::{
        diagram_generator::{EdgeColor, NodeColorWithVertexRule},
        FeynGenError,
    },
    gammaloop_integrand::{BareSample, DefaultSample},
    ltd::{generate_ltd_expression, LTDExpression},
    model::{self, ColorStructure, EdgeSlots, Model, Particle, VertexSlots},
    momentum::{FourMomentum, Polarization, Rotation, SignOrZero, Signature, ThreeMomentum},
    numerator::{
        ufo::{preprocess_ufo_color_wrapped, preprocess_ufo_spin_wrapped, UFO},
        AppliedFeynmanRule, ContractionSettings, Evaluate, Evaluators, ExtraInfo, GammaAlgebraMode,
        Numerator, NumeratorParseMode, NumeratorState, NumeratorStateError, PythonState,
        RepeatingIteratorTensorOrScalar, TypedNumeratorState, UnInit,
    },
    subtraction::{
        overlap::{find_maximal_overlap, OverlapStructure},
        static_counterterm::{self, CounterTerm},
    },
    utils::{self, sorted_vectorize, FloatLike, F, GS},
    ExportSettings, Settings, TropicalSubgraphTableSettings,
};

use linnet::half_edge::{
    builder::HedgeGraphBuilder,
    involution::{HedgePair, Orientation},
    nodestorage::NodeStorageOps,
    subgraph::SubGraphOps,
    HedgeGraph, HedgeGraphError, NodeIndex,
};

use ahash::{AHashMap, HashSet, RandomState};

use bincode::{Decode, Encode};
use color_eyre::Result;
use color_eyre::{Help, Report};
use enum_dispatch::enum_dispatch;
use eyre::eyre;
use itertools::Itertools;
use log::{debug, warn};
use momtrop::SampleGenerator;
use nalgebra::DMatrix;

use gat_lending_iterator::LendingIterator;
#[allow(unused_imports)]
use spenso::contraction::Contract;
use spenso::{
    arithmetic::ScalarMul,
    complex::Complex,
    contraction::{IsZero, RefZero},
    data::{DataTensor, DenseTensor, GetTensorData, SetTensorData, SparseTensor, StorageTensor},
    permutation::Permutation,
    scalar::Scalar,
    shadowing::{Shadowable, ETS},
    structure::{
        abstract_index::AbstractIndex,
        representation::{
            BaseRepName, ColorAdjoint, ColorFundamental, ColorSextet, Euclidean, Minkowski,
        },
        slot::{DualSlotTo, IsAbstractSlot},
        CastStructure, HasStructure, NamedStructure, ScalarTensor, ToSymbolic, VecStructure,
    },
};
use uuid::Uuid;

use core::panic;
use serde::{Deserialize, Serialize};
use smallvec::{smallvec, SmallVec};
use smartstring::{LazyCompact, SmartString};
use std::{
    collections::{HashMap, VecDeque},
    fmt::{Display, Formatter},
    fs,
    ops::{AddAssign, Neg, Not},
    path::{Path, PathBuf},
    sync::Arc,
};

use symbolica::{
    atom::{Atom, Symbol},
    coefficient::CoefficientView,
    domains::{float::NumericalFloatLike, rational::Rational},
    function,
    id::{Pattern, Replacement},
    state::State,
    symb,
};
use symbolica::{
    atom::{AtomCore, AtomView},
    graph::Graph as SymbolicaGraph,
    printer::{AtomPrinter, PrintOptions},
};

pub trait HedgeGraphExt<N, E> {
    type Error;
    fn from_sym(graph: SymbolicaGraph<N, E>) -> Self;

    fn to_sym(graph: &Self) -> Result<SymbolicaGraph<&N, &E>, Self::Error>;
}

impl<N: Clone, E: Clone, S: NodeStorageOps<NodeData = N>> HedgeGraphExt<N, E>
    for HedgeGraph<E, N, S>
{
    fn from_sym(graph: SymbolicaGraph<N, E>) -> Self {
        let mut builder = HedgeGraphBuilder::new();
        let mut map = AHashMap::new();

        for (i, node) in graph.nodes().iter().enumerate() {
            map.insert(i, builder.add_node(node.data.clone()));
        }

        // let mut edges = graph.edges().to_vec();
        // edges.sort_by(|a, b| a.vertices.cmp(&b.vertices));

        for edge in graph
            .edges()
            .iter()
            .sorted_by(|a, b| a.vertices.cmp(&b.vertices))
        {
            let vertices = edge.vertices;
            let source = map[&vertices.0];
            let sink = map[&vertices.1];
            builder.add_edge(source, sink, edge.data.clone(), edge.directed);
        }

        builder.into()
    }

    type Error = HedgeGraphError;

    fn to_sym(value: &HedgeGraph<E, N, S>) -> Result<SymbolicaGraph<&N, &E>, Self::Error> {
        let mut graph = SymbolicaGraph::new();
        let mut map = AHashMap::new();

        for (n, (_, node)) in value.iter_nodes().enumerate() {
            map.insert(NodeIndex(n), graph.add_node(node));
        }

        for (i, _, d) in value.iter_all_edges() {
            if let HedgePair::Paired { source, sink } = i {
                let source = map[&value.node_id(source)];
                let sink = map[&value.node_id(sink)];

                let data = d.data;
                let orientation = d.orientation;

                match orientation {
                    Orientation::Default => {
                        graph
                            .add_edge(source, sink, true, data)
                            .map_err(HedgeGraphError::SymbolicaError)?;
                    }
                    Orientation::Reversed => {
                        graph
                            .add_edge(sink, source, true, data)
                            .map_err(HedgeGraphError::SymbolicaError)?;
                    }
                    Orientation::Undirected => {
                        graph
                            .add_edge(source, sink, false, data)
                            .map_err(HedgeGraphError::SymbolicaError)?;
                    }
                }
            } else {
                return Err(HedgeGraphError::HasIdentityHedge);
            }
        }

        Ok(graph)
    }
}

//use symbolica::{atom::Symbol,state::State};
//pub mod half_edge;
#[derive(Debug, Default, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum EdgeType {
    #[serde(rename = "in")]
    Incoming,
    #[serde(rename = "out")]
    Outgoing,
    #[default]
    #[serde(rename = "virtual")]
    Virtual,
}

impl Display for EdgeType {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            EdgeType::Incoming => write!(f, "Incoming"),
            EdgeType::Outgoing => write!(f, "Outgoing"),
            EdgeType::Virtual => write!(f, "Virtual"),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[allow(non_snake_case)]
#[serde(tag = "type")]
pub enum SerializableVertexInfo {
    #[serde(rename = "external_vertex_info")]
    ExternalVertexInfo(SerializableExternalVertexInfo),
    #[serde(rename = "interacton_vertex_info")]
    InteractonVertexInfo(SerializableInteractionVertexInfo),
}

#[enum_dispatch]
pub trait HasVertexInfo {
    fn get_type(&self) -> SmartString<LazyCompact>;

    fn apply_vertex_rule(
        &self,
        edges: &[isize],
        vertex_pos: usize,
        vertex_slots: &VertexSlots,
    ) -> Option<[DataTensor<Atom>; 3]>;
}

#[derive(Debug, Clone)]
// #[enum_dispatch(HasVertexInfo)]
pub enum VertexInfo {
    ExternalVertexInfo(ExternalVertexInfo),
    InteractonVertexInfo(InteractionVertexInfo),
}

impl HasVertexInfo for VertexInfo {
    fn get_type(&self) -> SmartString<LazyCompact> {
        match self {
            VertexInfo::ExternalVertexInfo(e) => e.get_type(),
            VertexInfo::InteractonVertexInfo(i) => i.get_type(),
        }
    }

    fn apply_vertex_rule(
        &self,
        edges: &[isize],
        vertex_pos: usize,
        vertex_slots: &VertexSlots,
    ) -> Option<[DataTensor<Atom>; 3]> {
        match self {
            VertexInfo::ExternalVertexInfo(e) => {
                e.apply_vertex_rule(edges, vertex_pos, vertex_slots)
            }
            VertexInfo::InteractonVertexInfo(i) => {
                i.apply_vertex_rule(edges, vertex_pos, vertex_slots)
            }
        }
    }
}

impl VertexInfo {
    pub fn generate_vertex_slots(&self, shifts: Shifts, model: &Model) -> (VertexSlots, Shifts) {
        match self {
            VertexInfo::ExternalVertexInfo(e) => {
                let (e, mut updated_shifts) = match e.direction {
                    EdgeType::Outgoing => e.particle.get_anti_particle(model).slots(shifts),
                    EdgeType::Incoming => e.particle.slots(shifts),
                    EdgeType::Virtual => panic!("Virtual external vertex not supported"),
                };

                if updated_shifts.color == shifts.color {
                    updated_shifts.color += 1;
                }
                if updated_shifts.spin == shifts.spin {
                    updated_shifts.spin += 1;
                }
                if updated_shifts.lorentz == shifts.lorentz {
                    updated_shifts.lorentz += 1;
                }
                (e.into(), updated_shifts)
            }
            VertexInfo::InteractonVertexInfo(i) => i.vertex_rule.generate_vertex_slots(shifts),
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct SerializableInteractionVertexInfo {
    vertex_rule: SmartString<LazyCompact>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct SerializableExternalVertexInfo {
    direction: EdgeType,
    particle: SmartString<LazyCompact>,
}

#[derive(Debug, Clone)]
pub struct ExternalVertexInfo {
    pub direction: EdgeType,
    pub particle: Arc<model::Particle>,
}

impl ExternalVertexInfo {
    pub fn get_concrete_polarization_atom(
        &self,
        vertex_pos: usize,
        vertex_slots: &VertexSlots,
    ) -> Vec<Atom> {
        match self.direction {
            EdgeType::Incoming => self
                .particle
                .incoming_polarization_atom_concrete(&vertex_slots[0].dual(), vertex_pos),
            EdgeType::Outgoing => self
                .particle
                .outgoing_polarization_atom_concrete(&vertex_slots[0].dual(), vertex_pos),
            EdgeType::Virtual => panic!("Virtual external vertex not supported"),
        }
    }
}

impl HasVertexInfo for ExternalVertexInfo {
    fn get_type(&self) -> SmartString<LazyCompact> {
        SmartString::<LazyCompact>::from("external_vertex_info")
    }

    fn apply_vertex_rule(
        &self,
        _edges: &[isize],
        vertex_pos: usize,
        vertex_slots: &VertexSlots,
    ) -> Option<[DataTensor<Atom>; 3]> {
        let polarization = match self.direction {
            EdgeType::Incoming => self
                .particle
                .incoming_polarization_atom(&vertex_slots[0].dual(), vertex_pos),
            EdgeType::Outgoing => self
                .particle
                .outgoing_polarization_atom(&vertex_slots[0].dual(), vertex_pos),
            EdgeType::Virtual => panic!("Virtual external vertex not supported"),
        };

        Some([
            DataTensor::new_scalar(polarization),
            DataTensor::new_scalar(Atom::one()),
            DataTensor::new_scalar(Atom::one()),
        ])
    }
}

#[derive(Debug, Clone)]
pub struct InteractionVertexInfo {
    #[allow(unused)]
    pub vertex_rule: Arc<model::VertexRule>,
}

impl HasVertexInfo for InteractionVertexInfo {
    fn get_type(&self) -> SmartString<LazyCompact> {
        SmartString::<LazyCompact>::from("interacton_vertex_info")
    }

    fn apply_vertex_rule(
        &self,
        edges: &[isize],
        _vertex_pos: usize,
        vertex_slots: &VertexSlots,
    ) -> Option<[DataTensor<Atom>; 3]> {
        let spin_structure = self
            .vertex_rule
            .lorentz_structures
            .iter()
            .map(|ls| {
                let mut atom = ls.structure.clone();

                for (i, e) in edges.iter().enumerate() {
                    let momentum_in_pattern = Pattern::parse(&format!("P(x_,{})", i + 1)).unwrap();

                    let momentum_out_pattern = if e < &0 {
                        Pattern::parse(&format!("-Q({},mink(4,indexid(x_)))", -e)).unwrap()
                        //TODO flip based on flow
                    } else {
                        Pattern::parse(&format!("Q({},mink(4,indexid(x_)))", e)).unwrap()
                        //TODO flip based on flow
                    };

                    atom =
                        atom.replace_all(&momentum_in_pattern, &momentum_out_pattern, None, None);
                }

                atom = preprocess_ufo_spin_wrapped(atom);

                for (i, _) in edges.iter().enumerate() {
                    let replacements = vertex_slots[i].replacements(i + 1);

                    atom = atom.replace_all_multiple(&replacements);
                }

                let n_dummies = ls.number_of_dummies();
                for i in 0..n_dummies {
                    let pat: Pattern = Atom::parse(&format!("indexid({})", -1 - i as i64))
                        .unwrap()
                        .to_pattern();

                    atom = atom.replace_all(
                        &pat,
                        Atom::new_num(
                            usize::from(vertex_slots.internal_dummy.lorentz_and_spin[i]) as i64
                        )
                        .to_pattern(),
                        None,
                        None,
                    );
                }
                atom
            })
            .collect_vec();

        let mut color_dummy_shift = 0;
        let color_structure: Vec<Atom> = self
            .vertex_rule
            .color_structures
            .iter()
            .map(|cs| {
                let n_dummies = ColorStructure::number_of_dummies_in_atom(cs.as_view());
                let mut atom = cs.clone();

                atom = preprocess_ufo_color_wrapped(atom);

                let spins: Vec<isize> =
                    self.vertex_rule.particles.iter().map(|s| s.color).collect();

                for (i, s) in spins.iter().enumerate() {
                    let id1 = function!(UFO.identity, Atom::new_num((i + 1) as i32), symb!("x_"))
                        .to_pattern();
                    let id2 =
                        function!(ETS.id, symb!("x_"), Atom::new_num((i + 1) as i32)).to_pattern();

                    let ind = match s {
                        1 => Euclidean::slot(1, i + 1).to_symbolic_wrapped(),
                        3 => ColorFundamental::slot(3, i + 1).to_symbolic_wrapped(),
                        -3 => ColorFundamental::slot(3, i + 1)
                            .dual()
                            .to_symbolic_wrapped(),
                        6 => ColorSextet::slot(6, i + 1).to_symbolic_wrapped(),
                        -6 => ColorSextet::slot(6, i + 1).dual().to_symbolic_wrapped(),
                        8 => ColorAdjoint::slot(8, i + 1).to_symbolic_wrapped(),
                        i => panic!("Color {i} not supported "),
                    };

                    atom = atom.replace_all(
                        &id1,
                        function!(ETS.id, ind, symb!("x_")).to_pattern(),
                        None,
                        None,
                    );

                    atom = atom.replace_all(
                        &id2,
                        function!(ETS.id, symb!("x_"), ind).to_pattern(),
                        None,
                        None,
                    );
                }

                for (i, _) in edges.iter().enumerate() {
                    let replacements = vertex_slots[i].replacements(i + 1);

                    atom = atom.replace_all_multiple(&replacements);
                }

                for i in 0..n_dummies {
                    let pat: Pattern = Atom::parse(&format!("indexid({})", -1 - i as i64))
                        .unwrap()
                        .to_pattern();

                    atom = atom.replace_all(
                        &pat,
                        Atom::new_num(usize::from(vertex_slots.internal_dummy.color[i]) as i64)
                            .to_pattern(),
                        None,
                        None,
                    );
                }

                color_dummy_shift += n_dummies;

                atom
            })
            .collect();

        let [i, j] = vertex_slots.coupling_indices.unwrap();

        let color_structure = DataTensor::Dense(
            DenseTensor::from_data(color_structure, VecStructure::from(vec![i.cast()])).unwrap(),
        );

        let spin_structure = DataTensor::Dense(
            DenseTensor::from_data(spin_structure, VecStructure::from(vec![j.cast()])).unwrap(),
        );

        let mut couplings: DataTensor<Atom> =
            DataTensor::Sparse(SparseTensor::empty(VecStructure::from(vec![
                i.cast(),
                j.cast(),
            ])));

        for (i, row) in self.vertex_rule.couplings.iter().enumerate() {
            for (j, col) in row.iter().enumerate() {
                if let Some(atom) = col {
                    couplings.set(&[i, j], atom.expression.clone()).unwrap();
                }
            }
        }

        Some([spin_structure, color_structure, couplings])
    }
}

fn serializable_vertex_info(vertex_info: &VertexInfo) -> SerializableVertexInfo {
    match vertex_info {
        VertexInfo::ExternalVertexInfo(external_vertex_info) => {
            SerializableVertexInfo::ExternalVertexInfo(SerializableExternalVertexInfo {
                direction: external_vertex_info.direction,
                particle: external_vertex_info.particle.name.clone(),
            })
        }
        VertexInfo::InteractonVertexInfo(interaction_vertex_info) => {
            SerializableVertexInfo::InteractonVertexInfo(SerializableInteractionVertexInfo {
                vertex_rule: interaction_vertex_info.vertex_rule.name.clone(),
            })
        }
    }
}

fn deserialize_vertex_info(
    model: &model::Model,
    serialized_vertex_info: &SerializableVertexInfo,
) -> VertexInfo {
    match serialized_vertex_info {
        SerializableVertexInfo::ExternalVertexInfo(external_vertex_info) => {
            VertexInfo::ExternalVertexInfo(
                ExternalVertexInfo {
                    direction: external_vertex_info.direction,
                    particle: model.get_particle(&external_vertex_info.particle),
                }
                .clone(),
            )
        }
        SerializableVertexInfo::InteractonVertexInfo(interaction_vertex_info) => {
            VertexInfo::InteractonVertexInfo(InteractionVertexInfo {
                vertex_rule: model
                    .get_vertex_rule(&interaction_vertex_info.vertex_rule)
                    .clone(),
            })
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct SerializableEdge {
    name: SmartString<LazyCompact>,
    edge_type: EdgeType,
    particle: SmartString<LazyCompact>,
    propagator: SmartString<LazyCompact>,
    vertices: [SmartString<LazyCompact>; 2],
}

impl SerializableEdge {
    pub fn from_edge(graph: &BareGraph, edge: &Edge) -> SerializableEdge {
        SerializableEdge {
            name: edge.name.clone(),
            edge_type: edge.edge_type,
            particle: edge.particle.name.clone(),
            propagator: edge.propagator.name.clone(),
            vertices: [
                graph.vertices[edge.vertices[0]].name.clone(),
                graph.vertices[edge.vertices[1]].name.clone(),
            ],
        }
    }
}

#[derive(Debug, Clone)]
pub struct Edge {
    pub name: SmartString<LazyCompact>,
    pub edge_type: EdgeType,
    pub propagator: Arc<model::Propagator>,
    pub particle: Arc<model::Particle>,
    pub vertices: [usize; 2],
    pub internal_index: Vec<AbstractIndex>,
}

impl Edge {
    pub fn n_dummies(&self) -> usize {
        5
    }
    pub fn from_serializable_edge(
        model: &model::Model,
        graph: &BareGraph,
        serializable_edge: &SerializableEdge,
    ) -> Edge {
        Edge {
            name: serializable_edge.name.clone(),
            edge_type: serializable_edge.edge_type,
            particle: model.get_particle(&serializable_edge.particle),
            propagator: model.get_propagator(&serializable_edge.propagator),
            vertices: [
                graph
                    .get_vertex_position(&serializable_edge.vertices[0])
                    .unwrap(),
                graph
                    .get_vertex_position(&serializable_edge.vertices[1])
                    .unwrap(),
            ],
            internal_index: vec![],
        }
    }

    pub fn is_incoming_to(&self, vertex: usize) -> bool {
        self.vertices[1] == vertex
    }

    pub fn denominator(&self, graph: &BareGraph) -> (Atom, Atom) {
        let num = *graph.edge_name_to_position.get(&self.name).unwrap();
        let mom = Atom::parse(&format!("Q({num})")).unwrap();
        let mass = self
            .particle
            .mass
            .expression
            .clone()
            .unwrap_or(Atom::new_num(0));

        (mom, mass)
    }

    pub fn substitute_lmb(&self, atom: Atom, graph: &BareGraph, lmb: &LoopMomentumBasis) -> Atom {
        let num = *graph.edge_name_to_position.get(&self.name).unwrap();
        let mom = Pattern::parse(&format!("Q({num},x_)")).unwrap();
        let mom_rep = lmb.pattern(num);
        atom.replace_all(&mom, &mom_rep, None, None)
    }

    // pub fn edge_momentum_symbol(&self, graph: &Graph) -> Symbol {
    //     let num = *graph.edge_name_to_position.get(&self.name).unwrap();
    //     State::get_symbol(format!("Q{num}"))
    // }

    pub fn in_slot(&self, graph: &BareGraph) -> EdgeSlots<Minkowski> {
        let local_pos_in_sink_vertex =
            graph.vertices[self.vertices[0]].get_local_edge_position(self, graph, false);

        graph.vertex_slots[self.vertices[0]][local_pos_in_sink_vertex].dual()
    }

    pub fn out_slot(&self, graph: &BareGraph) -> EdgeSlots<Minkowski> {
        let local_pos_in_sink_vertex = graph.vertices[self.vertices[1]].get_local_edge_position(
            self,
            graph,
            self.vertices[0] == self.vertices[1],
        );

        graph.vertex_slots[self.vertices[1]][local_pos_in_sink_vertex].dual()
    }

    pub fn numerator(&self, graph: &BareGraph, num: usize) -> Atom {
        let [colorless, color] = self.color_separated_numerator(graph, num);

        colorless * color
    }

    pub fn color_separated_numerator(&self, graph: &BareGraph, num: usize) -> [Atom; 2] {
        // let num = *graph.edge_name_to_position.get(&self.name).unwrap();
        let in_slots = self.in_slot(graph);
        let out_slots = self.out_slot(graph);

        match self.edge_type {
            EdgeType::Incoming => {
                let [lorentz, spin, color] = in_slots.dual().kroneker(&out_slots);

                [lorentz * spin, color]
            }
            EdgeType::Outgoing => {
                let [lorentz, spin, color] = out_slots.dual().kroneker(&in_slots);

                [lorentz * spin, color]
            }
            EdgeType::Virtual => {
                let mut atom = self.propagator.numerator.clone();

                let pfun = Pattern::parse("P(x_)").unwrap();
                if self.particle.is_antiparticle() {
                    atom = atom.replace_all(
                        &pfun,
                        Pattern::parse(&format!("-Q({},mink(4,indexid(x_)))", num)).unwrap(),
                        None,
                        None,
                    );
                } else {
                    atom = atom.replace_all(
                        &pfun,
                        Pattern::parse(&format!("Q({},mink(4,indexid(x_)))", num)).unwrap(),
                        None,
                        None,
                    );
                }

                let pslashfun = Pattern::parse("PSlash(i_,j_)").unwrap();
                let pindex_num: usize = self.internal_index[0].into();
                if self.particle.is_antiparticle() {
                    atom = atom.replace_all(
                        &pslashfun,
                        Pattern::parse(&format!(
                            "-Q({},mink(4,{}))*Gamma({},i_,j_)",
                            num, pindex_num, pindex_num
                        ))
                        .unwrap(),
                        None,
                        None,
                    );
                } else {
                    atom = atom.replace_all(
                        &pslashfun,
                        Pattern::parse(&format!(
                            "Q({},mink(4,{}))*Gamma({},i_,j_)",
                            num, pindex_num, pindex_num
                        ))
                        .unwrap(),
                        None,
                        None,
                    );
                }

                atom = preprocess_ufo_spin_wrapped(atom);
                let indexidpat = Pattern::parse("indexid(x_)").unwrap();

                let dummies: HashSet<_> = atom
                    .pattern_match(&indexidpat, None, None)
                    .filter_map(|a| {
                        if let AtomView::Num(n) = a[&GS.x_].as_view() {
                            let e = if let CoefficientView::Natural(a, b) = n.get_coeff_view() {
                                if b == 1 {
                                    a
                                } else {
                                    0
                                }
                            } else {
                                0
                            };
                            if e < 0 {
                                Some(e)
                            } else {
                                None
                            }
                        } else {
                            None
                        }
                    })
                    .collect();

                let (replacements_in, replacements_out) = if self.particle.is_antiparticle() {
                    (in_slots.replacements(2), out_slots.replacements(1))
                } else {
                    (in_slots.replacements(1), out_slots.replacements(2))
                };

                // replacements_out.push(Replacement::new(
                //     Atom::parse("indexid(x_)").unwrap().to_pattern(),
                //     Atom::parse("x_").unwrap().to_pattern(),
                // ));

                let mut color_atom = Atom::new_num(1);
                for (&cin, &cout) in in_slots.color.iter().zip(out_slots.color.iter()) {
                    let id: NamedStructure<String, ()> =
                        NamedStructure::from_iter([cin, cout], "id".into(), None);
                    color_atom = color_atom * &id.to_symbolic().unwrap();
                }

                let reps: Vec<Replacement> = replacements_in
                    .into_iter()
                    .chain(replacements_out)
                    .collect();

                let atom = atom.replace_all_multiple(&reps);
                let color_atom = color_atom.replace_all_multiple(&reps);

                let indexid_reps: Vec<_> = dummies
                    .into_iter()
                    .enumerate()
                    .sorted()
                    .map(|(i, d)| {
                        Replacement::new(
                            Atom::parse(&format!("indexid({})", d))
                                .unwrap()
                                .to_pattern(),
                            Atom::parse(&format!("{}", self.internal_index[i + 1]))
                                .unwrap()
                                .to_pattern(),
                        )
                    })
                    .collect();

                let atom = atom.replace_all_multiple(&indexid_reps);
                let color_atom = color_atom.replace_all_multiple(&indexid_reps);

                [
                    atom.replace_all(
                        &Pattern::parse("indexid(x_)").unwrap(),
                        Atom::new_var(GS.x_).to_pattern(),
                        None,
                        None,
                    ),
                    color_atom.replace_all(
                        &Pattern::parse("indexid(x_)").unwrap(),
                        Atom::new_var(GS.x_).to_pattern(),
                        None,
                        None,
                    ),
                ]
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SerializableVertex {
    name: SmartString<LazyCompact>,
    vertex_info: SerializableVertexInfo,
    edges: Vec<SmartString<LazyCompact>>,
}

impl SerializableVertex {
    pub fn from_vertex(graph: &BareGraph, vertex: &Vertex) -> SerializableVertex {
        SerializableVertex {
            name: vertex.name.clone(),
            vertex_info: serializable_vertex_info(&vertex.vertex_info),
            edges: vertex
                .edges
                .iter()
                .map(|e| graph.edges[*e].name.clone())
                .collect(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Vertex {
    pub name: SmartString<LazyCompact>,
    pub vertex_info: VertexInfo,
    pub edges: Vec<usize>,
}

impl Vertex {
    pub fn get_local_edge_position(&self, edge: &Edge, graph: &BareGraph, skip_one: bool) -> usize {
        let global_id = graph.edge_name_to_position[&edge.name];
        let skip_n = if skip_one { 1 } else { 0 };

        self.edges
            .iter()
            .enumerate()
            .filter(|(_, &e)| e == global_id)
            .nth(skip_n)
            .unwrap()
            .0
    }

    pub fn generate_vertex_slots(&self, shifts: Shifts, model: &Model) -> (VertexSlots, Shifts) {
        self.vertex_info.generate_vertex_slots(shifts, model)
    }

    pub fn order_edges_following_interaction(
        &mut self,
        edge_id_and_pdgs_of_current_order: Vec<(usize, Arc<Particle>)>,
    ) -> Result<(), FeynGenError> {
        let mut new_edges_order = vec![];
        match self.vertex_info {
            VertexInfo::InteractonVertexInfo(ref i) => {
                let mut pdgs_to_position_map =
                    edge_id_and_pdgs_of_current_order.iter().collect::<Vec<_>>();
                for p in i.vertex_rule.particles.iter() {
                    let matched_pos = if let Some(pos) =
                        pdgs_to_position_map.iter().position(|(_, x)| *x == *p)
                    {
                        pos
                    } else {
                        return Err(FeynGenError::GenericError(
                                format!("Could not match some particles vertex ({}). It is incompatible with the ones in the interaction info ({})",
                                    edge_id_and_pdgs_of_current_order.iter().map(|(_,x)| x.name.clone()).join(","),
                                    i.vertex_rule.particles.iter().map(|x| x.name.clone()).join(","),
                                ),
                            ));
                    };
                    new_edges_order.push(pdgs_to_position_map[matched_pos].0);
                    pdgs_to_position_map.remove(matched_pos);
                }
                if !pdgs_to_position_map.is_empty() {
                    return Err(FeynGenError::GenericError(
                        format!("Not all particle of vertex ({}) were matched with the ones in the interaction info ({})",
                            edge_id_and_pdgs_of_current_order.iter().map(|(_,x)| x.name.clone()).join(","),
                            i.vertex_rule.particles.iter().map(|x| x.name.clone()).join(","),
                        ),
                    ));
                }
                self.edges = new_edges_order;
                Ok(())
            }

            VertexInfo::ExternalVertexInfo(_) => Ok(()),
        }
    }

    pub fn from_serializable_vertex(model: &model::Model, vertex: &SerializableVertex) -> Vertex {
        Vertex {
            name: vertex.name.clone(),
            vertex_info: deserialize_vertex_info(model, &vertex.vertex_info),
            // This will be filled in later during deserialization of the complete graph
            edges: vec![],
        }
    }

    pub fn is_edge_incoming(&self, edge: usize, graph: &BareGraph) -> bool {
        graph.is_edge_incoming(edge, graph.get_vertex_position(&self.name).unwrap())
    }

    fn add_signs_to_edges(&self, graph: &BareGraph) -> Vec<isize> {
        self.edges
            .iter()
            .map(|&e| {
                if !self.is_edge_incoming(e, graph) {
                    -(e as isize)
                } else {
                    e as isize
                }
            })
            .collect()
    }

    pub fn apply_vertex_rule(&self, graph: &BareGraph) -> Option<[DataTensor<Atom>; 3]> {
        let pos = graph.get_vertex_position(&self.name).unwrap();
        self.vertex_info.apply_vertex_rule(
            &self.add_signs_to_edges(graph),
            pos,
            &graph.vertex_slots[pos],
        )
    }

    pub fn contracted_vertex_rule(&self, graph: &BareGraph) -> Option<Atom> {
        let all = self.apply_vertex_rule(graph)?;
        let scalar = all
            .into_iter()
            .reduce(|acc, tensor| acc.contract(&tensor).unwrap())
            .unwrap()
            .get_owned([])
            .unwrap()
            .clone();

        Some(scalar)
    }

    pub fn colorless_vertex_rule(&self, graph: &BareGraph) -> Option<[DataTensor<Atom>; 2]> {
        let [spin, color, couplings] = self.apply_vertex_rule(graph)?;

        let colorless = spin.contract(&couplings).unwrap();

        Some([colorless, color])
    }

    #[allow(clippy::type_complexity)]
    pub fn contracted_colorless_vertex_rule(
        &self,
        graph: &BareGraph,
    ) -> Option<(Atom, DataTensor<Atom, NamedStructure<Symbol, usize>>)> {
        let [spin, color, couplings] = self.apply_vertex_rule(graph)?;

        let v = graph.get_vertex_position(&self.name).unwrap();

        let color = color.map_structure(|s| s.to_named(Symbol::new("Col"), Some(v)));

        let color_shadow = color.expanded_shadow().unwrap().cast_structure().into();

        let colorless = spin
            .contract(&couplings)
            .unwrap()
            .contract(&color_shadow)
            .unwrap()
            .scalar()
            .unwrap();

        Some((colorless, color))
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SerializableGraph {
    name: SmartString<LazyCompact>,
    vertices: Vec<SerializableVertex>,
    edges: Vec<SerializableEdge>,
    overall_factor: String,
    #[allow(clippy::type_complexity)]
    external_connections: Vec<(
        Option<SmartString<LazyCompact>>,
        Option<SmartString<LazyCompact>>,
    )>,
    loop_momentum_basis: Vec<SmartString<LazyCompact>>,
    #[serde(with = "sorted_vectorize")]
    edge_signatures: HashMap<SmartString<LazyCompact>, (Vec<isize>, Vec<isize>)>,
}

impl SerializableGraph {
    pub fn from_graph(graph: &BareGraph) -> SerializableGraph {
        SerializableGraph {
            name: graph.name.clone(),
            vertices: graph
                .vertices
                .iter()
                .map(|v| SerializableVertex::from_vertex(graph, v))
                .collect(),
            edges: graph
                .edges
                .iter()
                .map(|e| SerializableEdge::from_edge(graph, e))
                .collect(),

            overall_factor: graph.overall_factor.to_canonical_string(),
            external_connections: graph
                .external_connections
                .iter()
                .map(|(v1, v2)| {
                    (
                        v1.as_ref().map(|&v| graph.vertices[v].name.clone()),
                        v2.as_ref().map(|&v| graph.vertices[v].name.clone()),
                    )
                })
                .collect(),
            loop_momentum_basis: graph
                .loop_momentum_basis
                .basis
                .iter()
                .map(|&e| graph.edges[e].name.clone())
                .collect(),
            edge_signatures: graph
                .loop_momentum_basis
                .edge_signatures
                .iter()
                .enumerate()
                .map(|(i_e, sig)| (graph.edges[i_e].name.clone(), sig.to_momtrop_format()))
                .collect(),
        }
    }

    pub fn from_yaml_str(yaml_str: String) -> Result<SerializableGraph, Report> {
        serde_yaml::from_str(yaml_str.as_str())
            .map_err(|e| eyre!(format!("Error parsing graph yaml: {}", e)))
            .suggestion("Is it a correct yaml file")
    }
}

#[derive(Debug, Clone)]
pub struct Graph<S: NumeratorState = Evaluators> {
    pub bare_graph: BareGraph,
    pub derived_data: Option<DerivedGraphData<S>>,
}

impl Graph<PythonState> {
    pub fn load_derived_data<S: NumeratorState>(&mut self, path: &Path) -> Result<()> {
        let derived_data = DerivedGraphData::load_python_from_path::<S>(path)?;

        self.derived_data = Some(derived_data);

        Ok(())
    }

    pub fn numerator_apply<F, S: TypedNumeratorState, T: TypedNumeratorState>(
        &mut self,
        f: F,
    ) -> Result<()>
    where
        F: FnMut(Numerator<S>) -> Numerator<T> + Copy,
    {
        self.statefull_apply::<_, S, T>(|d, _| d.map_numerator(f))
    }

    pub fn forgetfull_apply<F, S: TypedNumeratorState>(&mut self, mut f: F) -> Result<()>
    where
        F: FnMut(DerivedGraphData<S>, &mut BareGraph) -> DerivedGraphData<PythonState>,
    {
        if let Some(d) = self.derived_data.take() {
            self.derived_data = Some(f(
                DerivedGraphData::<S>::try_from_python(d)?,
                &mut self.bare_graph,
            ));
            Ok(())
        } else {
            Err(eyre!("Derived data is None"))
        }
    }
    pub fn statefull_apply<F, S: TypedNumeratorState, T: TypedNumeratorState>(
        &mut self,
        mut f: F,
    ) -> Result<()>
    where
        F: FnMut(DerivedGraphData<S>, &mut BareGraph) -> DerivedGraphData<T>,
    {
        if let Some(d) = self.derived_data.take() {
            self.derived_data = Some(
                f(
                    DerivedGraphData::<S>::try_from_python(d)?,
                    &mut self.bare_graph,
                )
                .forget_type(),
            );
            Ok(())
        } else {
            Err(eyre!("Derived data is None"))
        }
    }

    pub fn statefull_apply_res<F, S: TypedNumeratorState, T: TypedNumeratorState>(
        &mut self,
        mut f: F,
    ) -> Result<()>
    where
        F: FnMut(
            DerivedGraphData<S>,
            &mut BareGraph,
        ) -> Result<DerivedGraphData<T>, DerivedGraphData<PythonState>>,
    {
        if let Some(d) = self.derived_data.take() {
            let res = f(
                DerivedGraphData::<S>::try_from_python(d)?,
                &mut self.bare_graph,
            );
            self.derived_data = Some(match res {
                Ok(r) => r.forget_type(),
                Err(e) => e,
            });
            Ok(())
        } else {
            Err(eyre!("Derived data is None"))
        }
    }
}

#[derive(Debug, Clone)]
pub struct BareGraph {
    pub name: SmartString<LazyCompact>,
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub external_edges: Vec<usize>,
    pub overall_factor: Atom,
    pub external_connections: Vec<(Option<usize>, Option<usize>)>,
    pub loop_momentum_basis: LoopMomentumBasis,
    pub vertex_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub edge_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub vertex_slots: Vec<VertexSlots>,
    pub shifts: Shifts,
    pub hedge_representation: HedgeGraph<usize, usize>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default)]
pub struct Shifts {
    pub lorentz: usize,
    pub lorentzdummy: usize,
    pub color: usize,
    pub colordummy: usize,
    pub spin: usize,
    pub coupling: usize,
}

impl BareGraph {
    pub fn rep_rules_print(&self, printer_ops: PrintOptions) -> Vec<(String, String)> {
        self.generate_lmb_replacement_rules("Q(<i>,x___)", "K(<i>,x___)", "P(<i>,x___)")
            .iter()
            .map(|(lhs, rhs)| {
                (
                    format!(
                        "{}",
                        AtomPrinter::new_with_options(lhs.as_view(), printer_ops)
                    ),
                    format!(
                        "{}",
                        AtomPrinter::new_with_options(rhs.as_view(), printer_ops)
                    ),
                )
            })
            .collect()
    }
    pub fn denominator_print(&self, printer_ops: PrintOptions) -> Vec<(String, String)> {
        self.edges
            .iter()
            .filter_map(|e| {
                if let EdgeType::Virtual = e.edge_type {
                    let (mom, mass) = e.denominator(self);
                    Some((
                        format!(
                            "{}",
                            AtomPrinter::new_with_options(mom.as_view(), printer_ops)
                        ),
                        format!(
                            "{}",
                            AtomPrinter::new_with_options(mass.as_view(), printer_ops)
                        ),
                    ))
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn external_slots(&self) -> Vec<EdgeSlots<Minkowski>> {
        self.vertices
            .iter()
            .enumerate()
            .filter_map(|(i, v)| {
                if matches!(v.vertex_info, VertexInfo::ExternalVertexInfo(_)) {
                    Some(self.vertex_slots[i][0].dual())
                } else {
                    None
                }
            })
            .collect()
        // self.external_edges()
        //     .iter()
        //     .map(|e| match e.edge_type {
        //         EdgeType::Incoming => e.out_slot(self),
        //         EdgeType::Outgoing => e.in_slot(self),
        //         _ => panic!("External edge is not incoming or outgoing"),
        //     })
        //     .collect()
    }

    pub fn external_edges(&self) -> Vec<&Edge> {
        self.external_edges
            .iter()
            .map(|&i| &self.edges[i])
            .collect()
    }
    pub fn external_particle_spin(&self) -> Vec<isize> {
        self.external_edges
            .iter()
            .map(|&i| self.edges[i].particle.spin)
            .collect()
    }

    pub fn external_particles(&self) -> Vec<Arc<Particle>> {
        self.external_edges
            .iter()
            .map(|&i| self.edges[i].particle.clone())
            .collect()
    }

    pub fn external_particle_spin_and_masslessness(&self) -> Vec<(isize, bool)> {
        self.external_edges
            .iter()
            .map(|&i| {
                (
                    self.edges[i].particle.spin,
                    self.edges[i].particle.mass.value.unwrap() == Complex::new_zero(),
                )
            })
            .collect()
    }
    pub fn is_tree(&self) -> bool {
        self.loop_momentum_basis.basis.is_empty()
    }

    pub fn external_in_or_out_signature(&self) -> Signature {
        self.external_edges
            .iter()
            .map(|&i| match self.edges[i].edge_type {
                EdgeType::Incoming => -1,
                EdgeType::Outgoing => 1i8,
                _ => panic!("External edge is not incoming or outgoing"),
            })
            .collect()
    }

    pub fn get_dependent_externals<T: FloatLike>(
        &self,
        indep_externals: &[FourMomentum<F<T>>],
    ) -> Vec<FourMomentum<F<T>>> {
        self.external_edges
            .iter()
            .map(|&i| {
                self.loop_momentum_basis.edge_signatures[i]
                    .external
                    .apply(indep_externals)
            })
            .collect()
    }

    pub fn dot_preamble(&self) -> String {
        [
            "digraph G {",
            format!("label=\"{}\";", self.name).as_str(),
            "noverlap=\"scale\"; layout=\"neato\";",
            "graph [ fontsize=10 ];",
            "node [ fontsize=7,shape=circle,margin=0,height=0.01 ];",
            "edge [ fontsize=5,arrowsize=0.4 ];",
        ]
        .join("\n")
    }

    pub fn dot(&self) -> String {
        let mut dot = String::new();
        dot.push_str(format!("{}\n", self.dot_preamble()).as_str());
        for (i, edge) in self.edges.iter().enumerate() {
            let from = self.vertices[edge.vertices[0]].name.clone();
            let to = self.vertices[edge.vertices[1]].name.clone();
            dot.push_str(&format!(
                "\"{}\" -> \"{}\" [label=\"{} | {} | Q({}) \"];\n",
                from, to, edge.name, edge.particle.name, i
            ));
        }
        dot.push_str("}\n");
        dot
    }

    pub fn dot_lmb(&self) -> String {
        let mut dot = String::new();
        dot.push_str(format!("{}\n", self.dot_preamble()).as_str());
        for (i, edge) in self.edges.iter().enumerate() {
            let from = self.vertices[edge.vertices[0]].name.clone();
            let to = self.vertices[edge.vertices[1]].name.clone();
            dot.push_str(&format!(
                "\"{}\" -> \"{}\" [label=\"{}\"];\n",
                from,
                to,
                self.loop_momentum_basis.edge_signatures[i].format_momentum()
            ));
        }
        dot.push_str("}\n");
        dot
    }

    pub fn dot_internal(&self) -> String {
        let mut dot = String::new();

        dot.push_str(format!("{}\n", self.dot_preamble()).as_str());
        for edge in &self.edges {
            let from = self.vertices[edge.vertices[0]].name.clone();
            let to = self.vertices[edge.vertices[1]].name.clone();
            dot.push_str(&format!(
                "\"{}\" -> \"{}\" [label=\"{}\"];\n",
                self.vertex_name_to_position[&from],
                self.vertex_name_to_position[&to],
                self.edge_name_to_position[&edge.name]
            ));
        }
        dot.push_str("}\n");
        dot
    }

    pub fn dot_internal_vertices(&self) -> String {
        let mut dot = String::new();
        // let mut seen_edges = HashSet::new();
        dot.push_str(format!("{}\n", self.dot_preamble()).as_str());
        for (vi, v) in self.vertices.iter().enumerate() {
            dot.push_str(&format!("\"{}\" [label=\"{}: {:?}\"] \n", vi, vi, v.edges));
            for e in &v.edges {
                let [from, to] = self.edges[*e].vertices;
                // if seen_edges.insert(e) {
                dot.push_str(&format!("\"{}\" -> \"{}\" [label=\"{}\"];\n", from, to, e));
                // }
            }
        }
        dot.push_str("}\n");
        dot
    }

    pub fn from_serializable_graph(model: &model::Model, graph: &SerializableGraph) -> BareGraph {
        // First build vertices
        let mut vertices: Vec<Vertex> = vec![];
        let mut vertex_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState> =
            HashMap::default();
        for vertex in &graph.vertices {
            let vertex = Vertex::from_serializable_vertex(model, vertex);
            vertex_name_to_position.insert(vertex.name.clone(), vertices.len());
            vertices.push(vertex);
        }

        let mut g = BareGraph {
            name: graph.name.clone(),
            vertices,
            edges: vec![],
            external_edges: vec![],
            overall_factor: Atom::parse(&graph.overall_factor).unwrap(),
            external_connections: vec![],
            loop_momentum_basis: LoopMomentumBasis {
                basis: vec![],
                edge_signatures: vec![],
            },
            vertex_name_to_position,
            edge_name_to_position: HashMap::default(),
            vertex_slots: vec![],
            shifts: Shifts::default(),
            hedge_representation: HedgeGraphBuilder::new().build(),
        };

        // let mut edges: Vec<Edge> = vec![];
        let mut edge_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState> =
            HashMap::default();
        for edge in &graph.edges {
            let edge = Edge::from_serializable_edge(model, &g, edge);
            edge_name_to_position.insert(edge.name.clone(), g.edges.len());
            g.edges.push(edge);
        }

        debug!("Loaded graph: {}", g.dot());
        g.external_edges = g
            .edges
            .iter()
            .enumerate()
            .filter_map(|(i, e)| {
                if e.edge_type == EdgeType::Incoming || e.edge_type == EdgeType::Outgoing {
                    Some(i)
                } else {
                    None
                }
            })
            .collect();
        for (i_e, e) in g.edges.iter().enumerate() {
            edge_name_to_position.insert(e.name.clone(), i_e);
        }
        for (vertex, serializable_vertex) in g.vertices.iter_mut().zip(graph.vertices.iter()) {
            vertex.edges = serializable_vertex
                .edges
                .iter()
                .map(|e| *edge_name_to_position.get(e).unwrap())
                .collect();
        }
        g.edge_name_to_position = edge_name_to_position;

        // set the half-edge graph representation
        g.hedge_representation = HedgeGraph::<usize, usize>::from(&g);

        let mut edge_signatures: Vec<_> = vec![None; graph.edges.len()];
        for (e_name, sig) in graph.edge_signatures.iter() {
            edge_signatures[g.get_edge_position(e_name).unwrap()] = Some(sig.clone().into());
        }
        g.loop_momentum_basis.edge_signatures =
            edge_signatures.into_iter().collect::<Option<_>>().unwrap();

        g.external_connections = graph
            .external_connections
            .iter()
            .map(|(v1, v2)| {
                (
                    v1.as_ref().map(|v| {
                        g.get_vertex_position(v).unwrap_or_else(|| {
                            panic!("Cannot find vertex named {} in this graph.", v)
                        })
                    }),
                    v2.as_ref().map(|v| {
                        g.get_vertex_position(v).unwrap_or_else(|| {
                            panic!("Cannot find vertex named {} in this graph.", v)
                        })
                    }),
                )
            })
            .collect();

        g.loop_momentum_basis.basis = graph
            .loop_momentum_basis
            .iter()
            .map(|e| g.get_edge_position(e).unwrap())
            .collect();

        g.generate_vertex_slots(model);

        // panic!("{:?}", g.edge_name_to_position);
        g.generate_internal_indices_for_edges();

        g
    }

    /*
    pub fn generate_lmb_replacements(&self) {
        let loop_basis = (0..self.loop_momentum_basis.basis.len())
            .map(|i_k| Atom::parse(&format!("K{}", i_k)).unwrap())
            .collect::<Vec<_>>();
        for (i_e, e) in self.edges.iter().enumerate() {
            i_e = self.edge_name_to_position.get(e.name).unwrap();

            self.loop_momentum_basis.edge_signatures[i_e]
                .internal
                .apply(basis);
            let mom = Atom::parse(&format!("Q{}", i)).unwrap();
            let mass = e
                .particle
                .mass
                .expression
                .clone()
                .unwrap_or(Atom::new_num(0));
            self.loop_momentum_basis.edge_signatures[i].external = vec![mom, mass];
        }
    }
    */

    pub fn select_cp_variant(
        model: &Model,
        vertex_info: &mut VertexInfo,
        edge_particles: &[Arc<Particle>],
    ) -> Result<(), FeynGenError> {
        let mut vertex_particles = edge_particles.to_vec();
        vertex_particles.sort();
        let vertex_rule = match vertex_info {
            VertexInfo::InteractonVertexInfo(ref mut i) => i.vertex_rule.clone(),
            VertexInfo::ExternalVertexInfo(_) => return Ok(()),
        };
        let mut this_interaction_particles = vertex_rule.particles.clone();
        this_interaction_particles.sort();
        if this_interaction_particles == vertex_particles {
            Ok(())
        } else {
            let mut this_interaction_particles_cp_conjugate = this_interaction_particles
                .iter()
                .map(|p| p.get_anti_particle(model))
                .collect::<Vec<_>>();
            this_interaction_particles_cp_conjugate.sort();
            if this_interaction_particles_cp_conjugate != vertex_particles {
                Err(FeynGenError::GenericError(
                    format!("Neither the vertex rule currently assigned nor its CP conjugate connect particles that match the PDG of the edges connected to this node: ({}) != ({})",
                        this_interaction_particles.iter().map(|p| model.get_particle_from_pdg(p.pdg_code).name.clone()).collect::<Vec<_>>().join(", "),
                        vertex_particles.iter().map(|p| model.get_particle_from_pdg(p.pdg_code).name.clone()).collect::<Vec<_>>().join(", "),
                    ),
                ))
            } else if let Some(vertex_rules) = model
                .particle_set_to_vertex_rules_map
                .get(&this_interaction_particles_cp_conjugate)
            {
                let this_rule_coupling_orders = vertex_rule.get_coupling_orders();
                let candidates = vertex_rules
                    .iter()
                    .filter_map(|vr| {
                        let vr_coupling_orders = vr.get_coupling_orders();
                        if vr_coupling_orders == this_rule_coupling_orders {
                            Some(vr.clone())
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>();
                if candidates.is_empty() {
                    return Err(FeynGenError::GenericError(
                        format!("Could not find the CP conjugate of this vertex in the Feynman rules of the model: ({}). None have matching couplings orders. Consider generating without the option '--symmetrize_left_right_states'.",
                        vertex_particles.iter().map(|p| model.get_particle_from_pdg(p.pdg_code).name.clone()).collect::<Vec<_>>().join(", "),
                        ),
                    ));
                } else if candidates.len() > 1 {
                    return Err(FeynGenError::GenericError(
                        format!("Could not find the CP conjugate of this vertex in the Feynman rules of the model: ({}). Multiple candidate interactions have matching coupling orders: {}. Consider generating without the option '--symmetrize_left_right_states'.",
                        vertex_particles.iter().map(|p| model.get_particle_from_pdg(p.pdg_code).name.clone()).collect::<Vec<_>>().join(", "),
                            candidates.iter().map(|vr| vr.name.clone()).collect::<Vec<_>>().join(", ")
                        ),
                    ));
                } else {
                    *vertex_info = VertexInfo::InteractonVertexInfo(InteractionVertexInfo {
                        vertex_rule: candidates[0].clone(),
                    });
                    Ok(())
                }
            } else {
                return Err(FeynGenError::GenericError(
                        format!("Could not find the CP conjugate of this vertex in the Feynman rules of the model: ({}). Consider generating without the option '--symmetrize_left_right_states'.",
                        vertex_particles.iter().map(|p| model.get_particle_from_pdg(p.pdg_code).name.clone()).collect::<Vec<_>>().join(", "),
                        ),
                    ));
            }
        }
    }

    pub fn from_symbolica_graph(
        model: &model::Model,
        name: String,
        graph: &SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
        symmetry_factor: Atom,
        external_connections: Vec<(Option<usize>, Option<usize>)>,
        forced_lmb: Option<Vec<SmartString<LazyCompact>>>,
    ) -> Result<BareGraph, FeynGenError> {
        let graph_nodes = graph.nodes().to_owned();
        let mut graph_edges = graph.edges().to_owned();

        // Useful adjacency matrix matrix
        let mut graph_adj_matrix: HashMap<usize, Vec<usize>> = HashMap::default();
        for e in graph_edges.iter() {
            graph_adj_matrix
                .entry(e.vertices.0)
                .or_default()
                .push(e.vertices.1);
            graph_adj_matrix
                .entry(e.vertices.1)
                .or_default()
                .push(e.vertices.0);
        }
        // Fix external edge directions, as for now *particle* fermion flow is observed, but we instead want:
        // > incoming antiparticles to be actually incoming (*particle* fermion flow would be outgoing)
        // > outgoing antiparticles to be actually outgoing (*particle* fermion flow would be incoming)
        let mut external_nodes: HashMap<usize, (usize, EdgeType, SmartString<LazyCompact>)> =
            HashMap::default();
        let mut external_edges: HashMap<usize, (EdgeType, SmartString<LazyCompact>)> =
            HashMap::default();
        let mut external_directions: HashMap<usize, EdgeType> = HashMap::default();
        for connection in external_connections.iter() {
            if let Some(leg_id) = connection.0 {
                external_directions.insert(leg_id, EdgeType::Incoming);
            }
            if let Some(leg_id) = connection.1 {
                external_directions.insert(leg_id, EdgeType::Outgoing);
            }
        }
        for (i_n, n) in graph_nodes.iter().enumerate() {
            if n.edges.len() == 1 {
                let external_edge = &mut graph_edges[n.edges[0]];
                let mut particle = model.get_particle_from_pdg(external_edge.data.pdg);
                let edge_type_from_symbolica = if graph_adj_matrix
                    .get(&external_edge.vertices.0)
                    .unwrap()
                    .len()
                    == 1
                {
                    EdgeType::Incoming
                } else if graph_adj_matrix
                    .get(&external_edge.vertices.1)
                    .unwrap()
                    .len()
                    == 1
                {
                    EdgeType::Outgoing
                } else {
                    panic!("Graph inconsistency in Feyngen.")
                };
                let physical_edge_type = *external_directions
                    .get(&(n.data.external_tag as usize))
                    .unwrap_or_else(|| {
                        panic!(
                            "External edge direction not specified for external leg {} in Feyngen.",
                            n.data.external_tag
                        )
                    });
                if edge_type_from_symbolica != physical_edge_type {
                    external_edge.vertices = (external_edge.vertices.1, external_edge.vertices.0);
                    particle = particle.get_anti_particle(model);
                }
                assert!(
                    n.data.external_tag > 0,
                    "External tag is not positive in Feyngen.",
                );
                external_nodes.insert(
                    n.data.external_tag as usize,
                    (i_n, physical_edge_type, particle.name.clone()),
                );
                external_edges.insert(n.edges[0], (physical_edge_type, particle.name.clone()));
            }
        }
        // First build vertices
        let mut vertices: Vec<Vertex> = vec![];
        let mut vertex_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState> =
            HashMap::default();
        // println!(
        //     "Symbolica graph edges = {:?}",
        //     graph
        //         .edges()
        //         .iter()
        //         .map(|e| format!("{:?} | {}", e.vertices, e.data))
        //         .collect::<Vec<_>>()
        // );
        // println!(
        //     "Symbolica graph nodes = {:?}",
        //     graph
        //         .nodes()
        //         .iter()
        //         .map(|n| format!("{:?} | {} | {}", n.edges.clone(), n.data.0, n.data.1))
        //         .collect::<Vec<_>>()
        // );
        for (i_n, node) in graph_nodes.iter().enumerate() {
            let vertex_info = if node.data.vertex_rule.name == "external" {
                let (_external_node_position, external_direction, external_particle_name) =
                    external_nodes
                        .get(&(node.data.external_tag as usize))
                        .unwrap();
                SerializableVertexInfo::ExternalVertexInfo(SerializableExternalVertexInfo {
                    direction: *external_direction,
                    particle: external_particle_name.clone(),
                })
            } else {
                if node.data.external_tag != 0 {
                    panic!("Internal vertex with non-zero integer colour data found in Feyngen.")
                }

                SerializableVertexInfo::InteractonVertexInfo(SerializableInteractionVertexInfo {
                    vertex_rule: node.data.vertex_rule.name.clone(), //node.data.vertex_rule.name.clone(),
                })
            };
            let vertex = Vertex::from_serializable_vertex(
                model,
                &SerializableVertex {
                    name: format!("v{}", i_n).into(),
                    vertex_info,
                    edges: vec![],
                },
            );
            vertex_name_to_position.insert(vertex.name.clone(), vertices.len());
            vertices.push(vertex);
        }

        let mut g = BareGraph {
            name: name.into(),
            vertices,
            edges: vec![],
            external_edges: vec![],
            overall_factor: symmetry_factor,
            external_connections: vec![],
            loop_momentum_basis: LoopMomentumBasis {
                basis: vec![],
                edge_signatures: vec![],
            },
            vertex_name_to_position: vertex_name_to_position.clone(),
            edge_name_to_position: HashMap::default(),
            vertex_slots: vec![],
            shifts: Shifts::default(),
            hedge_representation: HedgeGraphBuilder::new().build(),
        };

        let mut i_edge_internal = 0;
        let mut symbolica_edge_position_to_edge_name: HashMap<usize, SmartString<LazyCompact>> =
            HashMap::default();
        let mut edges_sorting_priority: HashMap<SmartString<LazyCompact>, usize> =
            HashMap::default();
        for (i_e, edge) in graph_edges.iter().enumerate() {
            let (edge_type, particle) =
                if let Some((edge_direction, particle_name)) = external_edges.get(&i_e) {
                    (*edge_direction, model.get_particle(particle_name))
                } else {
                    (
                        EdgeType::Virtual,
                        model.get_particle_from_pdg(edge.data.pdg),
                    )
                };

            let propagator: Arc<model::Propagator> =
                model.get_propagator_for_particle(&particle.name);

            let (start_vertex, end_vertex) = (edge.vertices.0, edge.vertices.1);

            let name = match edge_type {
                EdgeType::Incoming => {
                    edges_sorting_priority.insert(
                        format!("p{}", graph_nodes[edge.vertices.0].data.external_tag).into(),
                        graph_nodes[edge.vertices.0].data.external_tag as usize,
                    );
                    format!("p{}", graph_nodes[edge.vertices.0].data.external_tag)
                }
                EdgeType::Outgoing => {
                    edges_sorting_priority.insert(
                        format!("p{}", graph_nodes[edge.vertices.1].data.external_tag).into(),
                        graph_nodes[edge.vertices.1].data.external_tag as usize,
                    );
                    format!("p{}", graph_nodes[edge.vertices.1].data.external_tag)
                }
                _ => {
                    i_edge_internal += 1;
                    edges_sorting_priority.insert(
                        format!("q{}", i_edge_internal).into(),
                        i_edge_internal + external_edges.len() + 1,
                    );
                    format!("q{}", i_edge_internal)
                }
            };
            let edge = Edge {
                name: name.clone().into(),
                edge_type,
                particle,
                propagator,
                vertices: [start_vertex, end_vertex],
                internal_index: vec![],
            };
            symbolica_edge_position_to_edge_name.insert(i_e, name.into());
            g.edges.push(edge);
        }

        // Sort edges with external first according to their definition in the process
        g.edges.sort_by(|a, b| {
            let a_priority = edges_sorting_priority.get(&a.name).unwrap();
            let b_priority = edges_sorting_priority.get(&b.name).unwrap();
            a_priority.cmp(b_priority)
        });

        g.edge_name_to_position = HashMap::default();
        for (i_e, e) in g.edges.iter().enumerate() {
            g.edge_name_to_position.insert(e.name.clone(), i_e);
        }

        //debug!("Loaded graph: {}", g.dot());
        g.external_edges = g
            .edges
            .iter()
            .enumerate()
            .filter_map(|(i, e)| {
                if e.edge_type == EdgeType::Incoming || e.edge_type == EdgeType::Outgoing {
                    Some(i)
                } else {
                    None
                }
            })
            .collect();

        for (i_v, (vertex, node)) in g.vertices.iter_mut().zip(graph_nodes.iter()).enumerate() {
            vertex.edges = node
                .edges
                .iter()
                .map(|i_e| {
                    g.edge_name_to_position[symbolica_edge_position_to_edge_name.get(i_e).unwrap()]
                })
                .collect::<Vec<_>>();
            // We must honour the ordering given by the UFO interaction vertex info here
            let mut current_vertex_edge_order = vec![];
            for e_pos in &vertex.edges {
                let edge = &g.edges[*e_pos];
                // For a self-loop, the particle must be counted twice.
                if edge.vertices[0] == edge.vertices[1] {
                    if !edge.particle.is_self_antiparticle() {
                        // return Err(FeynGenError::GenericError(format!(
                        //     "Self-loop of edge {} *must* be a self-antiparticle",
                        //     edge.name
                        // )));
                        current_vertex_edge_order.push((*e_pos, edge.particle.clone()));
                        current_vertex_edge_order
                            .push((*e_pos, edge.particle.get_anti_particle(model).clone()));
                    } else {
                        current_vertex_edge_order.push((*e_pos, edge.particle.clone()));
                        current_vertex_edge_order.push((*e_pos, edge.particle.clone()));
                    }
                } else if edge.vertices[1] == i_v {
                    current_vertex_edge_order.push((*e_pos, edge.particle.clone()));
                } else if edge.vertices[0] == i_v {
                    current_vertex_edge_order
                        .push((*e_pos, edge.particle.get_anti_particle(model).clone()));
                } else {
                    return Err(FeynGenError::GenericError(format!(
                        "Edge {} is not connected to vertex {}",
                        edge.name, vertex.name
                    )));
                }
            }
            // When the graphs have been generated with --symmetrize_left_right_states, it may be that
            // the representative graphs contains the CP conjugate of the vertex rule.
            // In this case, we must select the correct vertex rule.
            // You can test this by generating the following process:
            //   generate a > d d~ [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling
            // where you will see that the left-right symmetrization will turn the original vertex (d~,u,W-) into (u~,d,W+)
            // so that we must select the correct CP correspondant vertex rule.
            // Carefully keeping track of the CP variant is important for capturing the right coupling phase in case of an complex CKM matrix for instance.
            BareGraph::select_cp_variant(
                model,
                &mut vertex.vertex_info,
                &current_vertex_edge_order
                    .iter()
                    .map(|(_, p)| p.clone())
                    .collect::<Vec<_>>(),
            )?;
            vertex.order_edges_following_interaction(current_vertex_edge_order)?;
        }

        let bare_graph_external_connections = external_connections
            .iter()
            .map(|(v1, v2)| {
                (
                    v1.map(|leg_id| {
                        assert!(
                            g.edges[g.external_edges[leg_id - 1]].edge_type == EdgeType::Incoming
                        );
                        g.edges[g.external_edges[leg_id - 1]].vertices[0]
                    }),
                    v2.map(|leg_id| {
                        assert!(
                            g.edges[g.external_edges[leg_id - 1]].edge_type == EdgeType::Outgoing
                        );
                        g.edges[g.external_edges[leg_id - 1]].vertices[1]
                    }),
                )
            })
            .collect::<Vec<_>>();
        g.external_connections = bare_graph_external_connections;

        // Set the half-edge graph representation
        g.hedge_representation = HedgeGraph::from(&g);

        g.set_loop_momentum_basis(&forced_lmb)?;

        g.generate_vertex_slots(model);

        // panic!("{:?}", g.edge_name_to_position);
        g.generate_internal_indices_for_edges();

        Ok(g)
    }

    pub fn set_loop_momentum_basis(
        &mut self,
        forced_lmb: &Option<Vec<SmartString<LazyCompact>>>,
    ) -> Result<(), FeynGenError> {
        let lmb_basis = if let Some(user_selected_lmb) = forced_lmb {
            let mut user_basis = vec![];
            for e_lmb in user_selected_lmb {
                if let Some(e_pos) = self.edge_name_to_position.get(e_lmb) {
                    user_basis.push(*e_pos);
                } else {
                    return Err(FeynGenError::LoopMomentumBasisError(format!(
                        "Specified loop momentum basis edge named '{}' not found in the graph.",
                        e_lmb
                    )));
                }
            }
            user_basis
        } else {
            // !g.hedge_representation
            //     .cycle_basis().1.
            //     .iter()
            //     .map(|cycle| {
            //         *g.hedge_representation
            //             .iter_egde_data(&cycle.internal_graph)
            //             .next()
            //             .unwrap()
            //     })
            //     .collect::<Vec<_>>()
            let spanning_tree = self.hedge_representation.cycle_basis().1;
            // println!("hedge graph: \n{}", g.hedge_representation.base_dot());
            // let spanning_tree_half_edge_node = sefl
            //     .hedge_representation
            //     .nesting_node_from_subgraph(spanning_tree.clone());
            // println!(
            //     "spanning tree: \n{}",
            //     self.hedge_representation.dot(&spanning_tree_half_edge_node)
            // );
            self.hedge_representation
                .iter_internal_edge_data(
                    &spanning_tree
                        .tree_subgraph
                        .complement(&self.hedge_representation),
                )
                .map(|e| *e.data)
                .collect::<Vec<_>>()
        };
        debug!(
            "Loop momentum basis edge positions selected: {:?}",
            lmb_basis
        );
        debug!(
            "Loop momentum basis selected: {}",
            lmb_basis
                .iter()
                .map(|i_e| self.edges[*i_e].clone().name)
                .collect::<Vec<_>>()
                .join(", ")
        );
        let mut lmb: LoopMomentumBasis = LoopMomentumBasis {
            basis: lmb_basis,
            edge_signatures: vec![],
        };
        lmb.set_edge_signatures(self).map_err(|e| {
            FeynGenError::LoopMomentumBasisError(format!(
                "{} | Error: {}",
                lmb.basis
                    .iter()
                    .map(|i_e| format!("{}", self.edges[*i_e].name))
                    .collect::<Vec<_>>()
                    .join(", "),
                e
            ))
        })?;
        self.loop_momentum_basis = lmb;
        Ok(())
    }

    pub fn verify_external_edge_order(&self) -> Result<Vec<usize>> {
        if self.external_edges.is_empty() {
            return Ok(vec![]);
        }
        let last = self.external_edges.len() - 1;
        let mut external_vertices_in_external_edge_order = vec![];
        for (i, ext) in self.external_edges.iter().enumerate() {
            let edge = &self.edges[*ext];

            match edge.edge_type {
                EdgeType::Incoming => match self.vertices[edge.vertices[0]].vertex_info {
                    VertexInfo::ExternalVertexInfo(_) => {
                        external_vertices_in_external_edge_order.push(edge.vertices[0]);
                    }
                    _ => {
                        return Err(eyre!(
                                "Incoming edge {} at position {i} is not connected to an external vertex with the correct direction",
                                edge.name
                            ));
                    }
                },
                EdgeType::Outgoing => match self.vertices[edge.vertices[1]].vertex_info {
                    VertexInfo::ExternalVertexInfo(_) => {
                        external_vertices_in_external_edge_order.push(edge.vertices[1]);
                    }
                    _ => {
                        return Err(eyre!(
                                "Outgoing edge {} at position {i} is not connected to an external vertex with the correct direction",
                                edge.name
                            ));
                    }
                },

                _ => {
                    return Err(eyre!(
                        "External edge {} at position {i} is not incoming or outgoing",
                        edge.name
                    ));
                }
            }

            // Skip verifications if signatures have not yet been assigned
            // if self.loop_momentum_basis.edge_signatures.is_empty() {
            //     continue;
            // }

            for s in &self.loop_momentum_basis.edge_signatures[*ext].internal {
                if !s.is_zero() {
                    return Err(eyre!(
                        "External edge {} at position {i} has a non-zero internal momentum signature",
                        edge.name
                    ));
                }
            }

            let edge_name = &edge.name;
            let signatures = &self.loop_momentum_basis.edge_signatures[*ext].external;

            for (j, &s) in signatures.iter().enumerate() {
                if j != i && s.is_sign() && i != last {
                    return Err(eyre!(
                        "External edge {} at position {i} has non-zero sign at position {} in its signature, expected zero",
                        edge_name,
                        j
                    ));
                } else if j != i && s.is_zero() && i == last {
                    return Err(eyre!(
                        "Last external edge {} at position {} in its signature has zero sign, expected non-zero, as it should be a sum of all other external edges",
                        edge_name,
                        j
                    ));
                }

                if j == i {
                    match (i, s) {
                        (idx, SignOrZero::Plus) if idx != last => {}
                        (idx, SignOrZero::Zero) if idx == last => {}
                        (idx, SignOrZero::Plus) if idx == last => {
                            return Err(eyre!(
                        "Last external edge {} at position {} should have zero sign, found positive",
                        edge_name, idx
                    ));
                        }
                        (idx, SignOrZero::Minus) => {
                            return Err(eyre!(
                                "External edge {} at position {} has a negative sign",
                                edge_name,
                                idx
                            ));
                        }
                        _ => {
                            return Err(eyre!(
                                "External edge {} at position {} has an unexpected sign",
                                edge_name,
                                i
                            ));
                        }
                    }
                }
            }
        }
        let mut sorted = external_vertices_in_external_edge_order.clone();

        // sorted.as_slice().is_sorted(); wait for 1.82
        // warn!(
        //     "External vertices are not in the same order as the external edges. This may cause issues."
        // }
        sorted.sort();

        let validate_external_vertices = self
            .vertices
            .iter()
            .enumerate()
            .filter(|(_, v)| matches!(v.vertex_info, VertexInfo::ExternalVertexInfo(_)))
            .map(|(i, _)| i)
            .collect::<Vec<_>>();

        if sorted != validate_external_vertices {
            return Err(eyre!("External vertices do not match the external edges."));
        }

        Ok(external_vertices_in_external_edge_order)
    }

    /// Generates the index ids needed based on the vertex rules.
    /// A first pass, passing around the shifts, is done to generate the vertex slots, that are seen by the edges.
    /// A second pass is done to generate the internal indices.
    fn generate_vertex_slots(&mut self, model: &Model) {
        let initial_shifts = Shifts {
            spin: 0,
            lorentzdummy: 0,
            colordummy: 0,
            color: 0,
            lorentz: 0,
            coupling: 0,
        };

        let ext_vertices = self.verify_external_edge_order().unwrap();

        let perm = Permutation::sort(&ext_vertices);

        let (mut external, s) = ext_vertices.into_iter().map(|i| &self.vertices[i]).fold(
            (Vec::new(), initial_shifts),
            |(mut acc, shifts), v| {
                let (e, new_shifts) = v.generate_vertex_slots(shifts, model);
                acc.push(e);
                (acc, new_shifts)
            },
        );

        perm.apply_slice_in_place_inv(&mut external);
        let mut external = VecDeque::from(external);

        let (mut v, s) = self
            .vertices
            .iter()
            .fold((vec![], s), |(mut acc, shifts), v| {
                if matches!(v.vertex_info, VertexInfo::ExternalVertexInfo(_)) {
                    acc.push(external.pop_front().unwrap());
                    return (acc, shifts);
                }
                let (e, new_shifts) = v.generate_vertex_slots(shifts, model);
                acc.push(e);
                (acc, new_shifts)
            });
        self.shifts = s;
        for slot in &mut v {
            slot.shift_internals(&self.shifts);
        }

        self.vertex_slots = v;
    }

    fn generate_internal_indices_for_edges(&mut self) {
        self.shifts.lorentzdummy =
            self.shifts.spin + self.shifts.lorentzdummy + self.shifts.lorentz;
        for edge in self.edges.iter_mut() {
            self.shifts.lorentzdummy += 1;
            edge.internal_index = vec![(self.shifts.lorentzdummy + 1).into()];
        }
    }

    #[inline]
    pub fn get_vertex(&self, name: &SmartString<LazyCompact>) -> Option<&Vertex> {
        match self.vertex_name_to_position.get(name) {
            Some(position) => Some(&self.vertices[*position]),
            None => None,
        }
    }

    #[inline]
    pub fn get_vertex_position(&self, name: &SmartString<LazyCompact>) -> Option<usize> {
        self.vertex_name_to_position.get(name).copied()
    }

    #[inline]
    pub fn get_edge(&self, name: &SmartString<LazyCompact>) -> Option<&Edge> {
        match self.edge_name_to_position.get(name) {
            Some(position) => Some(&self.edges[*position]),
            None => None,
        }
    }

    #[inline]
    pub fn get_edge_position(&self, name: &SmartString<LazyCompact>) -> Option<usize> {
        self.edge_name_to_position.get(name).copied()
    }

    #[inline]
    pub fn compute_emr<T: FloatLike>(
        &self,
        loop_moms: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
    ) -> Vec<ThreeMomentum<F<T>>> {
        self.compute_emr_in_lmb(loop_moms, external_moms, &self.loop_momentum_basis)
    }

    #[inline]
    pub fn compute_emr_in_lmb<T: FloatLike>(
        &self,
        loop_moms: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
        lmb: &LoopMomentumBasis,
    ) -> Vec<ThreeMomentum<F<T>>> {
        lmb.edge_signatures
            .iter()
            .map(|sig| sig.compute_three_momentum_from_four(loop_moms, external_moms))
            .collect()
    }

    #[inline]
    pub fn compute_onshell_energies<T: FloatLike>(
        &self,
        loop_moms: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
    ) -> Vec<F<T>> {
        self.compute_onshell_energies_in_lmb(loop_moms, external_moms, &self.loop_momentum_basis)
    }

    #[inline]
    pub fn get_mass_vector(&self) -> Vec<Option<Complex<F<f64>>>> {
        self.edges.iter().map(|e| e.particle.mass.value).collect()
    }

    #[inline]
    pub fn get_real_mass_vector(&self) -> Vec<F<f64>> {
        self.edges
            .iter()
            .map(|e| match e.particle.mass.value {
                Some(mass) => mass.re,
                None => F::from_f64(0.0),
            })
            .collect()
    }

    #[inline]
    pub fn compute_onshell_energies_in_lmb<T: FloatLike>(
        &self,
        loop_moms: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
        lmb: &LoopMomentumBasis,
    ) -> Vec<F<T>> {
        lmb.edge_signatures
            .iter()
            .map(|sig| sig.compute_four_momentum_from_three(loop_moms, external_moms))
            .zip(self.edges.iter())
            .map(|(emr_mom, edge)| match edge.edge_type {
                EdgeType::Virtual => {
                    emr_mom
                        .spatial
                        .on_shell_energy(edge.particle.mass.value.map(|m| {
                            if m.im.is_non_zero() {
                                panic!("Complex masses not yet supported in gammaLoop")
                            }
                            F::<T>::from_ff64(m.re)
                        }))
                        .value
                }
                _ => emr_mom.temporal.value, // a wierd way of just obtaining the energy of the external particles
            })
            .collect()
    }

    #[inline]
    pub fn compute_energy_product<T: FloatLike>(
        &self,
        loop_moms: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
    ) -> F<T> {
        self.compute_energy_product_in_lmb(loop_moms, external_moms, &self.loop_momentum_basis)
    }

    #[inline]
    pub fn compute_energy_product_in_lmb<T: FloatLike>(
        &self,
        loop_moms: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
        lmb: &LoopMomentumBasis,
    ) -> F<T> {
        let all_energies = self.compute_onshell_energies_in_lmb(loop_moms, external_moms, lmb);

        self.edges
            .iter()
            .zip(all_energies.iter())
            .filter(|(e, _)| e.edge_type == EdgeType::Virtual)
            .map(|(_, val)| F::<T>::from_f64(2.) * val) // why times 2? -lucien
            .fold(F::<T>::from_f64(1.), |acc, x| acc * x)
    }

    #[inline]
    pub fn get_edge_type_list(&self) -> Vec<EdgeType> {
        self.edges.iter().map(|e| e.edge_type).collect()
    }

    pub fn generate_lmb_replacement_rules(
        &self,
        q_format: &str,
        k_format: &str,
        p_format: &str,
    ) -> Vec<(Atom, Atom)> {
        self.loop_momentum_basis_replacement_rule(
            &self.loop_momentum_basis,
            q_format,
            k_format,
            p_format,
        )
    }

    fn loop_momentum_basis_replacement_rule(
        &self,
        lmb: &LoopMomentumBasis,
        q_format: &str,
        k_format: &str,
        p_format: &str,
    ) -> Vec<(Atom, Atom)> {
        let mut rule = vec![];

        for (i, signature) in lmb.edge_signatures.iter().enumerate() {
            rule.push((
                Atom::parse(
                    &q_format
                        .replace("<i>", &format!("{}", i))
                        .replace("<j>", &format!("{}", i)),
                )
                .unwrap(),
                self.replacement_rule_from_signature(i, signature, k_format, p_format),
            ));
        }

        rule
    }

    fn replacement_rule_from_signature(
        &self,
        index: usize,
        signature: &LoopExtSignature,
        k_format: &str,
        p_format: &str,
    ) -> Atom {
        let mut acc = Atom::new_num(0);
        for (i_l, &sign) in signature.internal.iter().enumerate() {
            let k = sign
                * Atom::parse(
                    &k_format
                        .replace("<i>", &format!("{}", i_l))
                        .replace("<j>", &format!("{}", index)),
                )
                .unwrap();
            acc = &acc + &k;
        }

        for (i_e, &sign) in signature.external.iter().enumerate() {
            let p = sign
                * Atom::parse(
                    &p_format
                        .replace("<i>", &format!("{}", i_e))
                        .replace("<j>", &format!("{}", index)),
                )
                .unwrap();
            acc = &acc + &p;
        }
        acc
    }

    pub fn denominator(self) -> Vec<(Atom, Atom)> {
        self.edges.iter().map(|e| e.denominator(&self)).collect()
    }

    pub fn cff_emr_from_lmb<T: FloatLike>(
        &self,
        sample: &BareSample<T>,
        lmb: &LoopMomentumBasis,
    ) -> Vec<FourMomentum<F<T>>> {
        let massless = lmb.to_massless_emr(sample);
        self.edges
            .iter()
            .zip(massless)
            .map(|(edge, emr)| match edge.edge_type {
                EdgeType::Virtual => {
                    emr.spatial
                        .into_on_shell_four_momentum(edge.particle.mass.value.map(|m| {
                            if m.im.is_non_zero() {
                                panic!("Complex masses not yet supported in gammaLoop")
                            }
                            F::<T>::from_ff64(m.re)
                        }))
                }
                _ => emr,
            })
            .collect()
    }

    pub fn evaluate_model_params(&mut self, _model: &Model) {}

    #[inline]
    pub fn get_virtual_edges_iterator(&self) -> impl Iterator<Item = (usize, &Edge)> {
        self.edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.edge_type == EdgeType::Virtual)
    }

    /// iterate over all edges which are part of a loop, (tree-level attachments removed)
    #[inline]
    pub fn get_loop_edges_iterator(&self) -> impl Iterator<Item = (usize, &Edge)> {
        self.edges.iter().enumerate().filter(|(index, e)| {
            e.edge_type == EdgeType::Virtual && {
                self.loop_momentum_basis.edge_signatures[*index]
                    .internal
                    .iter()
                    .any(|x| x.is_sign())
            }
        })
    }

    /// iterate over all edges which are virtual but do not carry a loop momentum.
    #[inline]
    pub fn get_tree_level_edges_iterator(&self) -> impl Iterator<Item = (usize, &Edge)> {
        self.edges.iter().enumerate().filter(|(index, e)| {
            e.edge_type == EdgeType::Virtual && {
                self.loop_momentum_basis.edge_signatures[*index]
                    .internal
                    .iter()
                    .all(|x| x.is_zero())
            }
        })
    }

    /// For a given edge, return the indices of the edges which have the same signature
    #[inline]
    pub fn is_edge_raised(&self, edge_index: usize) -> SmallVec<[usize; 2]> {
        let edge_signature = &self.loop_momentum_basis.edge_signatures[edge_index];
        let virtual_edges = self.get_virtual_edges_iterator();

        virtual_edges
            .filter(|(index, _)| {
                self.loop_momentum_basis.edge_signatures[*index] == *edge_signature
                    && *index != edge_index
            })
            .map(|(index, _)| index)
            .collect()
    }

    /// Returns groups of edges which all have the same signature
    pub fn group_edges_by_signature(&self) -> Vec<SmallVec<[usize; 3]>> {
        let mut edges: Vec<usize> = self
            .get_virtual_edges_iterator()
            .map(|(index, _)| index)
            .collect();

        let mut grouped_edges = Vec::with_capacity(edges.len());

        while !edges.is_empty() {
            let current_edge = edges.remove(0);

            let mut group = smallvec![current_edge];
            let mut index = 0;

            while index < edges.len() {
                if self.loop_momentum_basis.edge_signatures[current_edge]
                    == self.loop_momentum_basis.edge_signatures[edges[index]]
                {
                    group.push(edges.remove(index));
                } else {
                    index += 1;
                }
            }

            grouped_edges.push(group);
        }

        grouped_edges
    }

    pub fn is_edge_incoming(&self, edge: usize, vertex: usize) -> bool {
        self.edges[edge].is_incoming_to(vertex)
    }

    pub fn build_params_for_cff(&self) -> Vec<Atom> {
        self.edges
            .iter()
            .enumerate()
            .map(|(id, edge)| match edge.edge_type {
                EdgeType::Virtual => Atom::parse(&format!("E{}", id)).unwrap(),
                _ => Atom::parse(&format!("p{}", id)).unwrap(),
            })
            .collect()
    }

    pub fn get_dep_mom_expr(&self) -> (usize, ExternalShift) {
        let external_edges = self
            .edges
            .iter()
            .enumerate()
            .filter(|(_index, edge)| edge.edge_type != EdgeType::Virtual)
            .collect_vec();

        // find the external leg which does not appear in it's own signature
        if let Some((_, (dep_mom, _))) =
            external_edges
                .iter()
                .enumerate()
                .find(|(external_index, (index, _edge))| {
                    self.loop_momentum_basis.edge_signatures[*index].external[*external_index]
                        .is_positive()
                        .not()
                })
        {
            let dep_mom_signature = &self.loop_momentum_basis.edge_signatures[*dep_mom].external;

            let external_shift = external_edges
                .iter()
                .zip(dep_mom_signature.iter())
                .filter(|(_external_edge, dep_mom_sign)| dep_mom_sign.is_sign())
                .map(|((external_edge, _), dep_mom_sign)| (*external_edge, *dep_mom_sign as i64))
                .collect();

            (*dep_mom, external_shift)
        } else {
            // hack for vacuum diagrams
            (usize::MAX, ExternalShift::new())
        }
    }

    pub fn generate_loop_momentum_bases(&self) -> Vec<LoopMomentumBasis> {
        let loop_number = self.loop_momentum_basis.basis.len();
        let num_edges = self.edges.len();
        let external_signature_length = self.loop_momentum_basis.edge_signatures[0].external.len();

        // the full virtual signature matrix in the form of a flattened vector
        let signature_matrix_flattened = self
            .loop_momentum_basis
            .edge_signatures
            .iter()
            .flat_map(|sig| sig.internal.iter().map(|s| (*s as i8) as f64).collect_vec())
            .collect_vec();

        // convert to dmatrix
        let signature_matrix =
            DMatrix::from_row_slice(num_edges, loop_number, &signature_matrix_flattened);

        // the full external signature matrix in the form of a flattened vector
        let external_signature_matrix_flattened = self
            .loop_momentum_basis
            .edge_signatures
            .iter()
            .flat_map(|sig| sig.external.iter().map(|&s| (s as i8) as f64).collect_vec())
            .collect_vec();

        // convert to dmatrix
        let external_signature_matrix = DMatrix::from_row_slice(
            num_edges,
            external_signature_length,
            &external_signature_matrix_flattened,
        );

        let possible_lmbs = self
            .edges
            .iter()
            .filter(|e| e.edge_type == EdgeType::Virtual)
            .map(|e| self.get_edge_position(&e.name).unwrap())
            .combinations(loop_number);

        let valid_lmbs = possible_lmbs
            .map(|basis| {
                let reduced_signature_matrix_flattened = basis
                    .iter()
                    .flat_map(|e| {
                        self.loop_momentum_basis.edge_signatures[*e]
                            .internal
                            .iter()
                            .map(|s| (*s as i8) as f64)
                    })
                    .collect_vec();

                (
                    basis,
                    DMatrix::from_row_slice(
                        loop_number,
                        loop_number,
                        &reduced_signature_matrix_flattened,
                    ),
                )
            })
            .filter(|(_basis, reduced_signature_matrix)| {
                reduced_signature_matrix.determinant() != 0. // nonzero determinant means valid lmb
            })
            .map(|(basis, reduced_signature_matrix)| {
                let mut sorted_basis = basis;
                sorted_basis.sort();
                (sorted_basis, reduced_signature_matrix)
            });

        let lmbs = valid_lmbs
            .map(|(basis, reduced_signature_matrix)| {
                let mut new_signatures = self.loop_momentum_basis.edge_signatures.clone();

                // construct the signatures

                // for this we need the reduced external signatures
                let reduced_external_signatures_vec = basis
                    .iter()
                    .flat_map(|&e| {
                        self.loop_momentum_basis.edge_signatures[e]
                            .external
                            .iter()
                            .map(|s| (*s as i8) as f64)
                    })
                    .collect_vec();

                let reduced_external_signatures = DMatrix::from_row_slice(
                    loop_number,
                    external_signature_length,
                    &reduced_external_signatures_vec,
                );

                // also the inverse of the reduced_signature matrix
                let reduced_signature_matrix_inverse =
                    reduced_signature_matrix.clone().try_inverse().unwrap();

                let new_virtual_signatures =
                    signature_matrix.clone() * reduced_signature_matrix_inverse.clone();
                let new_external_signatures = external_signature_matrix.clone()
                    - signature_matrix.clone()
                        * reduced_signature_matrix_inverse.clone()
                        * reduced_external_signatures;

                // convert back to isize vecs
                for (edge_id, new_signature) in new_signatures.iter_mut().enumerate() {
                    let new_virtual_signature = new_virtual_signatures
                        .row(edge_id)
                        .iter()
                        .map(|s| s.round() as i8)
                        .collect();
                    let new_external_signature = new_external_signatures
                        .row(edge_id)
                        .iter()
                        .map(|s| s.round() as i8)
                        .collect();

                    *new_signature = LoopExtSignature {
                        internal: new_virtual_signature,
                        external: new_external_signature,
                    };
                }
                LoopMomentumBasis {
                    basis,
                    edge_signatures: new_signatures,
                }
            })
            .collect_vec();

        lmbs
    }

    pub fn generate_tropical_subgraph_table(
        &mut self,
        settings: &TropicalSubgraphTableSettings,
    ) -> Result<SampleGenerator<3>> {
        let num_virtual_loop_edges = self.get_loop_edges_iterator().count();
        let num_loops = self.loop_momentum_basis.basis.len();
        let target_omega = settings.target_omega;

        let weight = (target_omega + (3 * num_loops) as f64 / 2.) / num_virtual_loop_edges as f64;

        debug!(
            "Building tropical subgraph table with all edge weights set to: {}",
            weight
        );

        let tropical_edges = self
            .get_loop_edges_iterator()
            .map(|(_edge_id, edge)| {
                let is_massive = match edge.particle.mass.value {
                    Some(complex_mass) => !complex_mass.is_zero(),
                    None => false,
                };

                momtrop::Edge {
                    is_massive,
                    weight,
                    vertices: (edge.vertices[0] as u8, edge.vertices[1] as u8),
                }
            })
            .collect_vec();

        // Candidates to be external vertices of the tree-stripped graph
        let mut external_vertices_pool = HashSet::default();

        for edge in self
            .edges
            .iter()
            .filter(|e| e.edge_type == EdgeType::Incoming)
            .chain(
                self.edges
                    .iter()
                    .filter(|e| e.edge_type == EdgeType::Outgoing),
            )
        {
            external_vertices_pool.insert(edge.vertices[0] as u8);
            external_vertices_pool.insert(edge.vertices[1] as u8);
        }

        for (_, edge) in self.get_tree_level_edges_iterator() {
            external_vertices_pool.insert(edge.vertices[0] as u8);
            external_vertices_pool.insert(edge.vertices[1] as u8);
        }

        let mut external_vertices = vec![];

        for tropical_edge in &tropical_edges {
            if external_vertices_pool.contains(&tropical_edge.vertices.0) {
                external_vertices.push(tropical_edge.vertices.0);
            }

            if external_vertices_pool.contains(&tropical_edge.vertices.1) {
                external_vertices.push(tropical_edge.vertices.1);
            }
        }

        let tropical_graph = momtrop::Graph {
            edges: tropical_edges,
            externals: external_vertices,
        };

        let loop_part = self
            .get_loop_edges_iterator()
            .map(|(edge_id, _edge)| {
                self.loop_momentum_basis.edge_signatures[edge_id]
                    .internal
                    .clone()
                    .to_momtrop_format()
            })
            .collect_vec();

        tropical_graph
            .build_sampler(loop_part)
            .map_err(|e| eyre!(e))
    }

    pub fn load_derived_data<NumState: NumeratorState>(
        self,
        model: &Model,
        path: &Path,
        settings: &Settings,
    ) -> Result<Graph<NumState>, Report> {
        let derived_data_path = path.join(format!("derived_data_{}.bin", self.name.as_str()));
        let state_path = path.join("state.bin");
        debug!("Loading derived data from {:?}", derived_data_path);
        let mut derived_data = DerivedGraphData::load_from_path(&derived_data_path, &state_path)?;
        debug!("updating model in numerator");

        derived_data.numerator.update_model(model)?;

        derived_data
            .cff_expression
            .as_mut()
            .unwrap()
            .load_compiled(path.into(), self.name.clone(), settings)?;

        // if the user has edited the lmb in amplitude.yaml, this will set the right signature.
        let lmb_indices = self.loop_momentum_basis.basis.clone();

        let mut graph = Graph {
            bare_graph: self,
            derived_data: Some(derived_data),
        };
        graph.set_lmb(&lmb_indices)?;

        Ok(graph)
    }
}

impl Graph<UnInit> {
    pub fn from_serializable_graph(model: &model::Model, graph: &SerializableGraph) -> Self {
        Graph {
            derived_data: Some(DerivedGraphData::new_empty()),
            bare_graph: BareGraph::from_serializable_graph(model, graph),
        }
    }

    pub fn process_numerator(
        mut self,
        model: &Model,
        contraction_settings: ContractionSettings<Rational>,
        export_path: PathBuf,
        export_settings: &ExportSettings,
    ) -> Graph {
        let processed_data = self.derived_data.map(|d| {
            d.process_numerator(
                &mut self.bare_graph,
                model,
                contraction_settings,
                export_path,
                export_settings,
            )
            .unwrap()
        });
        Graph {
            bare_graph: self.bare_graph,
            derived_data: processed_data,
        }
    }

    pub fn apply_feynman_rules(
        mut self,
        export_settings: &ExportSettings,
    ) -> Graph<AppliedFeynmanRule> {
        let processed_data = self
            .derived_data
            .map(|d| d.apply_feynman_rules(&mut self.bare_graph, export_settings));
        Graph {
            bare_graph: self.bare_graph,
            derived_data: processed_data,
        }
    }
}

impl<S: TypedNumeratorState> Graph<S> {
    pub fn try_from_python(g: Graph<PythonState>) -> Result<Self> {
        let derived_data = if let Some(d) = g.derived_data {
            Some(DerivedGraphData::<S>::try_from_python(d)?)
        } else {
            None
        };
        Ok(Graph {
            bare_graph: g.bare_graph,
            derived_data,
        })
    }
}
impl<S: NumeratorState> Graph<S> {
    pub fn forget_type(self) -> Graph<PythonState> {
        Graph {
            bare_graph: self.bare_graph,
            derived_data: self.derived_data.map(|d| d.forget_type()),
        }
    }
    pub fn generate_cff(&mut self) {
        self.derived_data.as_mut().unwrap().cff_expression =
            Some(generate_cff_expression(&self.bare_graph).unwrap());
    }

    pub fn generate_ltd(&mut self) {
        self.derived_data.as_mut().unwrap().ltd_expression = Some(generate_ltd_expression(self));
    }
    pub fn generate_edge_groups(&mut self) {
        self.derived_data.as_mut().unwrap().edge_groups =
            Some(self.bare_graph.group_edges_by_signature());
    }

    pub fn generate_esurface_data(&mut self) -> Result<(), Report> {
        let data = generate_esurface_data(self, &self.get_cff().esurfaces)?;
        self.derived_data.as_mut().unwrap().esurface_derived_data = Some(data);

        Ok(())
    }

    // helper function
    #[inline]
    pub fn get_cff(&self) -> &CFFExpression {
        self.derived_data
            .as_ref()
            .unwrap()
            .cff_expression
            .as_ref()
            .unwrap()
    }

    #[inline]
    pub fn get_tropical_subgraph_table(&self) -> &SampleGenerator<3> {
        self.derived_data
            .as_ref()
            .unwrap()
            .tropical_subgraph_table
            .as_ref()
            .unwrap()
    }

    #[inline]
    pub fn get_esurface_derived_data(&self) -> &EsurfaceDerivedData {
        self.derived_data
            .as_ref()
            .unwrap()
            .esurface_derived_data
            .as_ref()
            .unwrap()
    }

    #[inline]
    pub fn get_existing_esurfaces<T: FloatLike>(
        &self,
        externals: &[FourMomentum<F<T>>],
        e_cm: F<f64>,
        settings: &Settings,
    ) -> ExistingEsurfaces {
        get_existing_esurfaces(
            &self.get_cff().esurfaces,
            self.get_esurface_derived_data(),
            externals,
            &self.bare_graph.loop_momentum_basis,
            settings.general.debug,
            e_cm,
        )
    }

    #[inline]
    pub fn get_maximal_overlap(
        &self,
        externals: &[FourMomentum<F<f64>>],
        e_cm: F<f64>,
        settings: &Settings,
    ) -> OverlapStructure {
        let existing_esurfaces = self.get_existing_esurfaces(externals, e_cm, settings);
        find_maximal_overlap(
            &self.bare_graph.loop_momentum_basis,
            &existing_esurfaces,
            &self.get_cff().esurfaces,
            &self.bare_graph.get_mass_vector(),
            externals,
            settings,
        )
    }

    pub fn build_compiled_expression(
        &mut self,
        export_path: PathBuf,
        export_settings: &ExportSettings,
    ) -> Result<(), Report> {
        let params = self.bare_graph.build_params_for_cff();
        match self.derived_data.as_mut().unwrap().cff_expression.as_mut() {
            Some(cff) => cff.build_compiled_expression::<f64>(
                &params,
                export_path,
                self.bare_graph.name.clone(),
                export_settings,
            ),
            None => {
                self.generate_cff();
                self.build_compiled_expression(export_path, export_settings)
            }
        }
    }
    // == Generation fns

    pub fn generate_loop_momentum_bases(&mut self) {
        self.derived_data.as_mut().unwrap().loop_momentum_bases =
            Some(self.bare_graph.generate_loop_momentum_bases());
    }
    pub fn generate_loop_momentum_bases_if_not_exists(&mut self) {
        if self
            .derived_data
            .as_ref()
            .unwrap()
            .loop_momentum_bases
            .is_none()
        {
            self.generate_loop_momentum_bases();
        }
    }

    // attempt to set a new loop momentum basis
    pub fn set_lmb(&mut self, lmb: &[usize]) -> Result<(), Report> {
        match &self.derived_data.as_ref().unwrap().loop_momentum_bases {
            None => Err(eyre!("lmbs not yet generated")),
            Some(lmbs) => {
                for (position, lmb_from_list) in lmbs.iter().enumerate() {
                    // search a matching lmb
                    if let Some(permutation_map) = utils::is_permutation(lmb, &lmb_from_list.basis)
                    {
                        // obtain the edge signatures
                        let mut new_edge_signatures = self
                            .derived_data
                            .as_ref()
                            .unwrap()
                            .loop_momentum_bases
                            .as_ref()
                            .unwrap()[position]
                            .edge_signatures
                            .clone();

                        // permutate the elemements of the loop part to match the ordering in the basis
                        for edge in new_edge_signatures.iter_mut() {
                            let new_loop_signature = edge
                                .internal
                                .iter()
                                .enumerate()
                                .map(|(ind, _)| edge.internal[permutation_map.right_to_left(ind)])
                                .collect();
                            edge.internal = new_loop_signature;
                        }

                        return Ok(());
                    }
                }

                Err(eyre!("lmb does not exist"))
            }
        }
    }

    #[inline]
    pub fn generate_params<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        _settings: &Settings,
    ) -> Vec<Complex<F<T>>> {
        self.derived_data
            .as_mut()
            .unwrap()
            .generate_params(&self.bare_graph, sample.possibly_rotated_sample())
    }

    pub fn generate_tropical_subgraph_table(&mut self, settings: &TropicalSubgraphTableSettings) {
        let table = self.bare_graph.generate_tropical_subgraph_table(settings);

        match table {
            Ok(table) => {
                debug!("min dod: {}", table.get_smallest_dod());
                if let Some(d) = &mut self.derived_data {
                    d.tropical_subgraph_table = Some(table);
                } else {
                    panic!("Derived data not initialized")
                }
            }
            Err(error) => {
                if settings.panic_on_fail {
                    panic!(
                        "Tropical subgraph table generation failed for graph: {} 🥥\n error: {:?}",
                        self.bare_graph.name, error
                    );
                } else {
                    warn!(
                        "Tropical subgraph table generation failed for graph: {} 🥥\n error: {:?}",
                        self.bare_graph.name, error
                    );
                }
            }
        }
    }

    pub fn map_numerator<F, T: NumeratorState>(mut self, f: F) -> Graph<T>
    where
        F: FnOnce(Numerator<S>, &mut BareGraph) -> Numerator<T>,
    {
        let mut_bare_graph = &mut self.bare_graph;
        let derived_data = self
            .derived_data
            .map(|d| d.map_numerator(|n| f(n, mut_bare_graph)));

        Graph {
            bare_graph: self.bare_graph,
            derived_data,
        }
    }

    pub fn map_numerator_res<E, F, T: NumeratorState>(mut self, f: F) -> Result<Graph<T>, E>
    where
        F: FnOnce(Numerator<S>, &mut BareGraph) -> Result<Numerator<T>, E>,
    {
        let mut_bare_graph = &mut self.bare_graph;
        let derived_data = self
            .derived_data
            .map(|d| d.map_numerator_res(|n| f(n, mut_bare_graph)))
            .transpose()?;

        Ok(Graph {
            bare_graph: self.bare_graph,
            derived_data,
        })
    }
}
impl Graph<Evaluators> {
    pub fn evaluate_fourd_expr<T: FloatLike>(
        &mut self,
        loop_mom: &[FourMomentum<F<T>>],
        external_mom: &[FourMomentum<F<T>>],
        polarizations: &[Polarization<Complex<F<T>>>],
        settings: &Settings,
    ) -> DataTensor<Complex<F<T>>> {
        self.derived_data.as_mut().unwrap().evaluate_fourd_expr(
            loop_mom,
            external_mom,
            polarizations,
            settings,
            &self.bare_graph,
        )
    }

    #[inline]
    pub fn evaluate_cff_expression<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        settings: &Settings,
    ) -> Complex<F<T>> {
        self.derived_data
            .as_mut()
            .unwrap()
            .evaluate_cff_expression(&self.bare_graph, sample, settings)
            .scalar()
            .unwrap()
    }

    #[inline]
    pub fn evaluate_cff_expression_in_lmb<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        lmb_specification: &LoopMomentumBasisSpecification,
        settings: &Settings,
    ) -> Complex<F<T>> {
        self.derived_data
            .as_mut()
            .unwrap()
            .evaluate_cff_expression_in_lmb(&self.bare_graph, sample, lmb_specification, settings)
            .scalar()
            .unwrap()
    }
    #[inline]
    pub fn evaluate_cff_all_orientations<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        settings: &Settings,
    ) -> Complex<F<T>> {
        self.derived_data
            .as_mut()
            .unwrap()
            .evaluate_cff_all_orientations(
                &self.bare_graph,
                sample.possibly_rotated_sample(),
                settings,
            )
    }

    #[inline]
    pub fn evaluate_numerator_all_orientations<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        settings: &Settings,
    ) -> DataTensor<Complex<F<T>>> {
        self.derived_data
            .as_mut()
            .unwrap()
            .evaluate_numerator_all_orientations(
                &self.bare_graph,
                sample.possibly_rotated_sample(),
                sample.uuid(),
                settings,
            )
    }

    #[inline]
    pub fn evaluate_ltd_expression<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        settings: &Settings,
    ) -> Complex<F<T>> {
        self.derived_data
            .as_mut()
            .unwrap()
            .evaluate_ltd_expression(sample, settings, &self.bare_graph)
            .scalar()
            .unwrap()
    }

    #[inline]
    pub fn evaluate_ltd_expression_in_lmb<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        lmb: &LoopMomentumBasis,
        settings: &Settings,
    ) -> Complex<F<T>> {
        self.derived_data
            .as_mut()
            .unwrap()
            .evaluate_ltd_expression_in_lmb(sample, settings, lmb, &self.bare_graph)
            .scalar()
            .unwrap()
    }

    #[inline]
    /// evaluates the cff expression at the given loop momenta and external momenta. The loop momenta are assumed to be in the loop momentum basis specified, and have irrelevant energy components. The output is a vector of complex numbers, one for each orientation.
    pub fn evaluate_cff_orientations<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        lmb_specification: &LoopMomentumBasisSpecification,
        settings: &Settings,
    ) -> Vec<F<T>> {
        let lmb = lmb_specification.basis(self);
        let energy_cache = self.bare_graph.compute_onshell_energies_in_lmb(
            sample.loop_moms(),
            sample.external_moms(),
            lmb,
        );

        self.derived_data
            .as_ref()
            .unwrap()
            .cff_expression
            .as_ref()
            .unwrap()
            .evaluate_orientations(&energy_cache, settings)
    }

    pub fn evaluate_threshold_counterterm<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        rotation_for_overlap: &Rotation,
        settings: &Settings,
    ) -> Complex<F<T>> {
        self.derived_data
            .as_mut()
            .unwrap()
            .evaluate_threshold_counterterm(
                &self.bare_graph,
                sample,
                rotation_for_overlap,
                settings,
            )
    }
}

#[allow(dead_code)]
#[derive(Debug, Clone, Deserialize, Serialize, Encode, Decode)]
pub struct DerivedGraphData<NumState> {
    pub loop_momentum_bases: Option<Vec<LoopMomentumBasis>>,
    pub cff_expression: Option<CFFExpression>,
    pub ltd_expression: Option<LTDExpression>,
    #[bincode(with_serde)]
    pub tropical_subgraph_table: Option<SampleGenerator<3>>,
    #[bincode(with_serde)]
    pub edge_groups: Option<Vec<SmallVec<[usize; 3]>>>,
    pub esurface_derived_data: Option<EsurfaceDerivedData>,
    pub static_counterterm: Option<static_counterterm::CounterTerm>,
    pub numerator: Numerator<NumState>,
}

impl DerivedGraphData<PythonState> {
    pub fn load_python_from_path<S: NumeratorState>(path: &Path) -> Result<Self, Report> {
        let mut source = std::fs::File::open(path)?;
        let mut statemap = State::import(&mut source, None)?;
        match std::fs::read(path) {
            Ok(derived_data_bytes) => {
                let derived_data: DerivedGraphData<S> = bincode::decode_from_slice_with_context(
                    &derived_data_bytes,
                    bincode::config::standard(),
                    &mut statemap,
                )?
                .0;
                Ok(derived_data.forget_type())
            }
            Err(_) => {
                Err(eyre!(
                    "Could not read derived data from path: {}",
                    path.display()
                ))
                // Ok(Self::new_empty())
            }
        }
    }

    pub fn apply<F, S: TypedNumeratorState, T: TypedNumeratorState>(
        &mut self,
        f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Numerator<S>) -> Numerator<T>,
    {
        self.numerator.apply(f)
    }
}

// impl<NumState: NumeratorState> Serialize for DerivedGraphData<NumState> {
//     fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
//     where
//         S: serde::Serializer,
//     {
//         let mut state = serializer.serialize_struct("DerivedGraphData", 7)?;
//         state.serialize_field("loop_momentum_bases", &self.loop_momentum_bases)?;
//         state.serialize_field("cff_expression", &self.cff_expression)?;
//         state.serialize_field("ltd_expression", &self.ltd_expression)?;
//         state.serialize_field("tropical_subgraph_table", &self.tropical_subgraph_table)?;
//         state.serialize_field("edge_groups", &self.edge_groups)?;
//         state.serialize_field("esurface_derived_data", &self.esurface_derived_data)?;
//         state.serialize_field("static_counterterm", &self.static_counterterm)?;
//         state.serialize_field("numerator", &self.numerator)?;
//         state.end()
//     }
// }

impl<NumState: NumeratorState> DerivedGraphData<NumState> {
    pub fn map_numerator<F, T: NumeratorState>(self, f: F) -> DerivedGraphData<T>
    where
        F: FnOnce(Numerator<NumState>) -> Numerator<T>,
    {
        DerivedGraphData {
            loop_momentum_bases: self.loop_momentum_bases,
            cff_expression: self.cff_expression,
            ltd_expression: self.ltd_expression,
            tropical_subgraph_table: self.tropical_subgraph_table,
            edge_groups: self.edge_groups,
            esurface_derived_data: self.esurface_derived_data,
            static_counterterm: self.static_counterterm,
            numerator: f(self.numerator),
        }
    }

    pub fn map_numerator_res<E, F, T: NumeratorState>(self, f: F) -> Result<DerivedGraphData<T>, E>
    where
        F: FnOnce(Numerator<NumState>) -> Result<Numerator<T>, E>,
    {
        Ok(DerivedGraphData {
            loop_momentum_bases: self.loop_momentum_bases,
            cff_expression: self.cff_expression,
            ltd_expression: self.ltd_expression,
            tropical_subgraph_table: self.tropical_subgraph_table,
            edge_groups: self.edge_groups,
            esurface_derived_data: self.esurface_derived_data,
            static_counterterm: self.static_counterterm,
            numerator: f(self.numerator)?,
        })
    }
}

impl<NumState: TypedNumeratorState> DerivedGraphData<NumState> {
    fn try_from_python(value: DerivedGraphData<PythonState>) -> Result<Self> {
        value.map_numerator_res(|n| n.try_from())
    }
}

impl DerivedGraphData<Evaluators> {
    pub fn evaluate_fourd_expr<T: FloatLike>(
        &mut self,
        loop_moms: &[FourMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
        polarizations: &[Polarization<Complex<F<T>>>],
        settings: &Settings,
        bare_graph: &BareGraph,
    ) -> DataTensor<Complex<F<T>>> {
        let emr = bare_graph
            .loop_momentum_basis
            .edge_signatures
            .iter()
            .map(|sig| sig.compute_momentum(loop_moms, external_moms))
            .collect_vec();

        let mut den = Complex::new_re(F::from_f64(1.));
        for (e, q) in bare_graph.edges.iter().zip(emr.iter()) {
            if e.edge_type == EdgeType::Virtual {
                // println!("q: {}", q);
                if let Some(mass) = e.particle.mass.value {
                    let m2 = mass.norm_squared();
                    let m2: F<T> = F::from_ff64(m2);
                    den *= &q.square() - &m2;
                } else {
                    den *= q.square();
                }
            }
        }
        // println!("den: {}", den);
        let den = den.inv();

        let num = self
            .numerator
            .evaluate_single(&emr, polarizations, None, settings);
        num.scalar_mul(&den).unwrap()
    }

    #[inline]
    pub fn evaluate_ltd_expression<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        settings: &Settings,
        bare_graph: &BareGraph,
    ) -> DataTensor<Complex<F<T>>> {
        let one = sample.one();
        let zero = one.zero();
        let i = Complex::new(zero, one);
        let loop_number = bare_graph.loop_momentum_basis.basis.len();
        // Unexplained overall (-1)^(L+1) sign to match with cFF which we know is correct
        let prefactor = (Complex::new(-sample.one(), sample.zero())).pow((loop_number + 1) as u64)
            * i.pow(loop_number as u64);

        self.ltd_expression
            .as_ref()
            .unwrap()
            .evaluate(sample, bare_graph, &mut self.numerator, settings)
            .scalar_mul(&prefactor)
            .unwrap()
    }

    #[inline]
    pub fn evaluate_ltd_expression_in_lmb<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        settings: &Settings,
        lmb: &LoopMomentumBasis,
        bare_graph: &BareGraph,
    ) -> DataTensor<Complex<F<T>>> {
        let one = sample.one();
        let zero = one.zero();
        let i = Complex::new(zero, one);
        let loop_number = bare_graph.loop_momentum_basis.basis.len();
        // Unexplained overall (-1)^(L+1) sign to match with cFF which we know is correct
        let prefactor = (Complex::new(-sample.one(), sample.zero())).pow((loop_number + 1) as u64)
            * i.pow(loop_number as u64);

        self.ltd_expression
            .as_ref()
            .unwrap()
            .evaluate_in_lmb(sample, bare_graph, lmb, &mut self.numerator, settings)
            .scalar_mul(&prefactor)
            .unwrap()
    }

    #[inline]
    /// evaluates the cff expression at the given loop momenta and external momenta. The loop momenta are assumed to be in the loop momentum basis specified, and have irrelevant energy components. The output is a vector of complex numbers, one for each orientation.
    pub fn evaluate_cff_expression_in_lmb<T: FloatLike>(
        &mut self,
        graph: &BareGraph,
        sample: &DefaultSample<T>,
        lmb_specification: &LoopMomentumBasisSpecification,
        settings: &Settings,
    ) -> DataTensor<Complex<F<T>>> {
        let one = sample.one();
        let zero = one.zero();
        let ni = Complex::new(zero.clone(), -one.clone());

        let loop_number = graph.loop_momentum_basis.basis.len();

        let prefactor = ni.pow((loop_number) as u64); // (ni).pow(loop_number as u64) * (-(ni.ref_one())).pow(internal_vertex_number as u64 - 1);

        let mut cff = self
            .evaluate_cff_orientations(
                graph,
                sample.possibly_rotated_sample(),
                lmb_specification,
                settings,
            )
            .into_iter();

        let (num_sample, uuid) = sample.numerator_sample(settings);
        let num_iter = self.evaluate_numerator_orientations(
            graph,
            num_sample,
            uuid,
            lmb_specification,
            settings,
        );

        match num_iter {
            RepeatingIteratorTensorOrScalar::Scalars(mut s) => {
                if let Some(i) = s.next() {
                    // println!("num: {}", i);
                    let c = Complex::new_re(cff.next().unwrap());
                    // println!("cff: {}", c);
                    let mut sum = i * &c;

                    for j in cff.by_ref() {
                        let c = Complex::new_re(j);
                        // println!("cff: {}", c);
                        let num = s.next().unwrap();
                        // println!("num: {}", num);
                        let summand = &c * num;
                        sum += summand;
                    }
                    sum *= prefactor;
                    DataTensor::new_scalar(sum)
                } else {
                    panic!("Empty iterator in sum");
                }
            }
            RepeatingIteratorTensorOrScalar::Tensors(mut s) => {
                if let Some(i) = s.next() {
                    let c = Complex::new_re(cff.next().unwrap());
                    let mut sum = i.map_data_ref(|n| n * &c);

                    for j in cff {
                        let c = Complex::new_re(j);
                        sum += s.next().unwrap().map_data_ref(|n| n * &c);
                    }
                    sum.scalar_mul(&prefactor).unwrap()
                } else {
                    panic!("Empty iterator in sum");
                }
            }
        }
    }

    pub fn evaluate_numerator_all_orientations<T: FloatLike>(
        &mut self,
        graph: &BareGraph,
        sample: &BareSample<T>,
        tag: Option<Uuid>,
        settings: &Settings,
    ) -> DataTensor<Complex<F<T>>> {
        let emr = graph.cff_emr_from_lmb(sample, &graph.loop_momentum_basis);

        let rep = self
            .numerator
            .evaluate_all_orientations(&emr, &sample.polarizations, tag, settings)
            .unwrap();

        match rep {
            RepeatingIteratorTensorOrScalar::Tensors(mut t) => {
                if let Some(i) = t.next() {
                    let mut sum = i.clone();

                    while let Some(j) = t.next() {
                        sum += j;
                    }
                    sum
                } else {
                    panic!("Empty iterator in sum");
                }
            }
            RepeatingIteratorTensorOrScalar::Scalars(mut s) => {
                if let Some(i) = s.next() {
                    let mut sum = i.clone();

                    while let Some(j) = s.next() {
                        sum += j;
                    }
                    DataTensor::new_scalar(sum)
                } else {
                    panic!("Empty iterator in sum");
                }
            }
        }
    }

    pub fn evaluate_threshold_counterterm<T: FloatLike>(
        &mut self,
        graph: &BareGraph,
        sample: &DefaultSample<T>,
        rotation_for_overlap: &Rotation,
        settings: &Settings,
    ) -> Complex<F<T>> {
        match self.static_counterterm.as_ref() {
            Some(ct) => CounterTerm::evaluate(
                sample,
                graph,
                &self.cff_expression.as_ref().unwrap().esurfaces,
                ct,
                &mut self.numerator,
                rotation_for_overlap,
                settings,
            ),
            None => Complex::new(sample.zero(), sample.zero()),
        }
    }

    #[inline]
    /// evaluates the numerator at the given loop momenta and external momenta. The loop momenta are assumed to be in the loop momentum basis specified, and have irrelevant energy components. The output is a vector of complex numbers, one for each orientation.
    pub fn evaluate_numerator_orientations<T: FloatLike>(
        &mut self,
        graph: &BareGraph,
        sample: &BareSample<T>,
        tag: Option<Uuid>,
        lmb_specification: &LoopMomentumBasisSpecification,
        settings: &Settings,
    ) -> RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<T>>>> {
        let lmb = lmb_specification.basis_from_derived(self);
        let emr = graph.cff_emr_from_lmb(sample, lmb);

        self.numerator
            .evaluate_all_orientations(&emr, &sample.polarizations, tag, settings)
            .unwrap()
    }

    #[inline]
    pub fn evaluate_cff_expression<T: FloatLike>(
        &mut self,
        graph: &BareGraph,
        sample: &DefaultSample<T>,
        settings: &Settings,
    ) -> DataTensor<Complex<F<T>>> {
        let lmb_specification = LoopMomentumBasisSpecification::Literal(&graph.loop_momentum_basis);
        self.evaluate_cff_expression_in_lmb(graph, sample, &lmb_specification, settings)
    }
}

impl DerivedGraphData<UnInit> {
    fn new_empty() -> Self {
        Self {
            loop_momentum_bases: None,
            cff_expression: None,
            ltd_expression: None,
            tropical_subgraph_table: None,
            edge_groups: None,
            esurface_derived_data: None,
            numerator: Numerator::default(),
            static_counterterm: None,
        }
    }

    pub fn process_numerator_no_eval(
        self,
        base_graph: &mut BareGraph,
        contraction_settings: ContractionSettings<Rational>,
        export_path: PathBuf,
        export_settings: &ExportSettings,
    ) -> Result<DerivedGraphData<PythonState>> {
        let expr_path = export_path.join("expressions");
        std::fs::create_dir_all(&expr_path)?;
        let extra_info = self.generate_extra_info(export_path);

        let dump_numerator = export_settings.numerator_settings.dump_expression;

        if let Some(n) = dump_numerator {
            let opts = PrintOptions::from(n);
            let dens: Vec<(String, String)> = base_graph.denominator_print(opts);
            let rep_rules: Vec<(String, String)> = base_graph.rep_rules_print(opts);

            fs::write(
                expr_path.join(format!("{}_model_dens.json", base_graph.name)),
                serde_json::to_string_pretty(&(rep_rules, dens)).unwrap(),
            )?;
        }
        let color_simplified =
            if let Some(global) = &export_settings.numerator_settings.global_numerator {
                debug!("Using global numerator: {}", global);
                let global = Atom::parse(global).unwrap();
                self.map_numerator(|n| {
                    let n = n.from_global(
                        global,
                        // base_graph,
                        &export_settings.numerator_settings.global_prefactor,
                    );

                    if let Some(numerator_format) = dump_numerator {
                        n.write(&extra_info, base_graph, numerator_format).unwrap();
                    }

                    n.color_simplify()
                    // .color_project()
                })
            } else {
                self.map_numerator(|n| {
                    let n = n.from_graph(
                        base_graph,
                        &export_settings.numerator_settings.global_prefactor,
                    );

                    if let Some(numerator_format) = dump_numerator {
                        n.write(&extra_info, base_graph, numerator_format).unwrap();
                    }
                    n.color_simplify()
                })
                //.color_project())
            };
        Ok(match &export_settings.numerator_settings.gamma_algebra {
            GammaAlgebraMode::Symbolic => color_simplified
                .map_numerator_res(|n| {
                    if let Some(numerator_format) = dump_numerator {
                        n.write(&extra_info, base_graph, numerator_format).unwrap();
                    }
                    let n = n.gamma_simplify();
                    if let Some(numerator_format) = dump_numerator {
                        n.write(&extra_info, base_graph, numerator_format).unwrap();
                    }
                    let n = n.parse().contract(contraction_settings)?;
                    if let Some(numerator_format) = dump_numerator {
                        n.write(&extra_info, base_graph, numerator_format).unwrap();
                    }
                    Result::<_, Report>::Ok(n)
                })?
                .forget_type(),
            GammaAlgebraMode::Concrete => match &export_settings.numerator_settings.parse_mode {
                NumeratorParseMode::Polynomial => color_simplified
                    .map_numerator_res(|n| {
                        if let Some(numerator_format) = dump_numerator {
                            n.write(&extra_info, base_graph, numerator_format).unwrap();
                        }
                        let n = n.parse_poly(base_graph).contract()?;
                        if let Some(numerator_format) = dump_numerator {
                            n.write(&extra_info, base_graph, numerator_format).unwrap();
                        }
                        Result::<_, Report>::Ok(n)
                    })?
                    .forget_type(),
                NumeratorParseMode::Direct => color_simplified
                    .map_numerator_res(|n| {
                        if let Some(numerator_format) = dump_numerator {
                            n.write(&extra_info, base_graph, numerator_format).unwrap();
                        }
                        let n = n.parse().contract(contraction_settings)?;
                        if let Some(numerator_format) = dump_numerator {
                            n.write(&extra_info, base_graph, numerator_format).unwrap();
                        }
                        Result::<_, Report>::Ok(n)
                    })?
                    .forget_type(),
            },
        })
    }

    pub fn process_numerator(
        self,
        base_graph: &mut BareGraph,
        model: &Model,
        contraction_settings: ContractionSettings<Rational>,
        export_path: PathBuf,
        export_settings: &ExportSettings,
    ) -> Result<DerivedGraphData<Evaluators>> {
        let expr_path = export_path.join("expressions");
        std::fs::create_dir_all(&expr_path)?;
        let extra_info = self.generate_extra_info(export_path);

        let dump_numerator = export_settings.numerator_settings.dump_expression;

        if let Some(n) = dump_numerator {
            let opts = PrintOptions::from(n);
            let dens: Vec<(String, String)> = base_graph.denominator_print(opts);
            let rep_rules: Vec<(String, String)> = base_graph.rep_rules_print(opts);

            fs::write(
                expr_path.join(format!("{}_model_dens.json", base_graph.name)),
                serde_json::to_string_pretty(&(rep_rules, dens)).unwrap(),
            )?;
        }

        let color_simplified =
            if let Some(global) = &export_settings.numerator_settings.global_numerator {
                debug!("Using global numerator: {}", global);
                let global = Atom::parse(global).unwrap();
                self.map_numerator(|n| {
                    let n = n.from_global(
                        global,
                        // base_graph,
                        &export_settings.numerator_settings.global_prefactor,
                    );

                    if let Some(numerator_format) = dump_numerator {
                        n.write(&extra_info, base_graph, numerator_format).unwrap();
                    }
                    n.color_simplify()
                    // .color_project()
                })
            } else {
                self.map_numerator(|n| {
                    let n = n.from_graph(
                        base_graph,
                        &export_settings.numerator_settings.global_prefactor,
                    );

                    if let Some(numerator_format) = dump_numerator {
                        n.write(&extra_info, base_graph, numerator_format).unwrap();
                    }
                    n.color_simplify()
                })
                //.color_project())
            };
        Ok(match &export_settings.numerator_settings.gamma_algebra {
            GammaAlgebraMode::Symbolic => color_simplified
                .map_numerator_res(|n| {
                    if let Some(numerator_format) = dump_numerator {
                        n.write(&extra_info, base_graph, numerator_format).unwrap();
                    }
                    let n = n.gamma_simplify();
                    if let Some(numerator_format) = dump_numerator {
                        n.write(&extra_info, base_graph, numerator_format).unwrap();
                    }
                    n.parse().contract(contraction_settings)
                })?
                .map_numerator(|n| {
                    if let Some(numerator_format) = dump_numerator {
                        n.write(&extra_info, base_graph, numerator_format).unwrap();
                    }
                    n.generate_evaluators(model, base_graph, &extra_info, export_settings)
                }),
            GammaAlgebraMode::Concrete => match &export_settings.numerator_settings.parse_mode {
                NumeratorParseMode::Polynomial => color_simplified
                    .map_numerator_res(|n| {
                        if let Some(numerator_format) = dump_numerator {
                            n.write(&extra_info, base_graph, numerator_format).unwrap();
                        }
                        n.parse_poly(base_graph).contract()
                    })?
                    .map_numerator(|n| {
                        if let Some(numerator_format) = dump_numerator {
                            n.write(&extra_info, base_graph, numerator_format).unwrap();
                        }
                        n.generate_evaluators(model, base_graph, &extra_info, export_settings)
                    }),
                NumeratorParseMode::Direct => color_simplified
                    .map_numerator_res(|n| {
                        if let Some(numerator_format) = dump_numerator {
                            n.write(&extra_info, base_graph, numerator_format).unwrap();
                        }
                        n.parse().contract(contraction_settings)
                    })?
                    .map_numerator(|n| {
                        if let Some(numerator_format) = dump_numerator {
                            n.write(&extra_info, base_graph, numerator_format).unwrap();
                        }
                        n.generate_evaluators(model, base_graph, &extra_info, export_settings)
                    }),
            },
        })
    }

    fn apply_feynman_rules(
        self,
        base_graph: &mut BareGraph,
        export_settings: &ExportSettings,
    ) -> DerivedGraphData<AppliedFeynmanRule> {
        self.map_numerator(|n| {
            n.from_graph(
                base_graph,
                &export_settings.numerator_settings.global_prefactor,
            )
        })
    }
}

impl<NumState: NumeratorState> DerivedGraphData<NumState> {
    pub fn forget_type(self) -> DerivedGraphData<PythonState> {
        self.map_numerator(|n| n.forget_type())
    }
    pub fn generate_extra_info(&self, export_path: PathBuf) -> ExtraInfo {
        ExtraInfo {
            orientations: self
                .cff_expression
                .as_ref()
                .unwrap()
                .orientations
                .iter()
                .map(|a| a.orientation.clone())
                .collect(),
            path: export_path,
        }
    }
    pub fn evaluate_cff_all_orientations<T: FloatLike>(
        &mut self,
        graph: &BareGraph,
        sample: &BareSample<T>,
        settings: &Settings,
    ) -> Complex<F<T>> {
        let lmb_specification = LoopMomentumBasisSpecification::Literal(&graph.loop_momentum_basis);
        Complex {
            re: self
                .evaluate_cff_orientations(graph, sample, &lmb_specification, settings)
                .into_iter()
                .reduce(|acc, e| acc + &e)
                .unwrap_or_else(|| panic!("no orientations to evaluate")),
            im: F::new_zero(),
        }
    }

    #[inline]
    /// evaluates the cff expression at the given loop momenta and external momenta. The loop momenta are assumed to be in the loop momentum basis specified, and have irrelevant energy components. The output is a vector of complex numbers, one for each orientation.
    pub fn evaluate_cff_orientations<T: FloatLike>(
        &self,
        graph: &BareGraph,
        sample: &BareSample<T>,
        lmb_specification: &LoopMomentumBasisSpecification,
        settings: &Settings,
    ) -> Vec<F<T>> {
        let lmb = lmb_specification.basis_from_derived(self);
        let energy_cache =
            graph.compute_onshell_energies_in_lmb(&sample.loop_moms, &sample.external_moms, lmb);

        self.cff_expression
            .as_ref()
            .unwrap()
            .evaluate_orientations(&energy_cache, settings)
    }

    pub fn generate_params<T: FloatLike>(
        &mut self,
        graph: &BareGraph,
        sample: &BareSample<T>,
    ) -> Vec<Complex<F<T>>> {
        let lmb_specification = LoopMomentumBasisSpecification::Literal(&graph.loop_momentum_basis);
        let lmb = lmb_specification.basis_from_derived(self);
        let emr = graph.cff_emr_from_lmb(sample, lmb);

        let mut params: Vec<Complex<F<T>>> = emr
            .into_iter()
            .flat_map(|p| p.into_iter().map(Complex::new_re))
            .collect_vec();

        for p in &sample.polarizations {
            for pi in p {
                params.push(pi.clone());
            }
        }
        params
    }

    pub fn load_from_path(derived_data_path: &Path, state_path: &Path) -> Result<Self, Report> {
        let mut source = std::fs::File::open(state_path)?;
        let mut statemap = State::import(&mut source, None)?;
        match std::fs::read(derived_data_path) {
            Ok(derived_data_bytes) => {
                let derived_data: Self = bincode::decode_from_slice_with_context(
                    &derived_data_bytes,
                    bincode::config::standard(),
                    &mut statemap,
                )?
                .0;
                Ok(derived_data)
            }
            Err(_) => {
                Err(eyre!(
                    "Could not read derived data from path: {}",
                    derived_data_path.display()
                ))
                // Ok(Self::new_empty())
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct LoopMomentumBasis {
    pub basis: Vec<usize>,
    pub edge_signatures: Vec<LoopExtSignature>,
}

#[derive(
    Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, PartialOrd, Eq, Ord, Hash,
)]
pub struct LoopExtSignature {
    pub internal: Signature,
    pub external: Signature,
}

impl From<(Signature, Signature)> for LoopExtSignature {
    fn from((loops, external_signature): (Signature, Signature)) -> Self {
        Self {
            internal: loops,
            external: external_signature,
        }
    }
}

impl From<(Vec<isize>, Vec<isize>)> for LoopExtSignature {
    fn from((loops, external_signature): (Vec<isize>, Vec<isize>)) -> Self {
        Self {
            internal: Signature::from_iter(loops),
            external: Signature::from_iter(external_signature),
        }
    }
}

impl LoopExtSignature {
    pub fn compute_momentum<'a, 'b: 'a, T>(&self, loop_moms: &'a [T], external_moms: &'b [T]) -> T
    where
        T: RefZero + Clone + Neg<Output = T> + AddAssign<T>,
    {
        if loop_moms.is_empty() {
            return self.external.apply(external_moms);
        }
        if external_moms.is_empty() {
            return self.internal.apply(loop_moms);
        }
        let mut res = self.internal.apply(loop_moms);
        res += self.external.apply(external_moms);
        res
    }

    pub fn to_momtrop_format(&self) -> (Vec<isize>, Vec<isize>) {
        (
            self.internal.to_momtrop_format(),
            self.external.to_momtrop_format(),
        )
    }

    /// Usefull for debugging
    pub fn format_momentum(&self) -> String {
        let mut res = String::new();
        let mut first = true;

        for (i, sign) in (&self.internal).into_iter().enumerate() {
            if !first {
                res.push_str(&sign.to_string());
            } else {
                first = false;
            }
            if sign.is_sign() {
                res.push_str(&format!("k_{}", i));
            }
        }

        for (i, sign) in (&self.external).into_iter().enumerate() {
            if !first {
                res.push_str(&sign.to_string());
            } else {
                first = false;
            }
            if sign.is_sign() {
                res.push_str(&format!("l_{}", i));
            }
        }

        res
    }

    #[allow(unused)]
    pub fn compute_four_momentum_from_three<T: FloatLike>(
        &self,
        loop_moms: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
    ) -> FourMomentum<F<T>> {
        let loop_moms = loop_moms
            .iter()
            .map(|m| m.clone().into_on_shell_four_momentum(None))
            .collect_vec();
        self.compute_momentum(&loop_moms, external_moms)
    }

    pub fn compute_three_momentum_from_four<T: FloatLike>(
        &self,
        loop_moms: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
    ) -> ThreeMomentum<F<T>> {
        let external_moms = external_moms
            .iter()
            .map(|m| m.spatial.clone())
            .collect_vec();
        self.compute_momentum(loop_moms, &external_moms)
    }
}

impl LoopMomentumBasis {
    pub fn spatial_emr<T: FloatLike>(&self, sample: &BareSample<T>) -> Vec<ThreeMomentum<F<T>>> {
        let three_externals = sample
            .external_moms
            .iter()
            .map(|m| m.spatial.clone())
            .collect_vec();
        self.edge_signatures
            .iter()
            .map(|sig| sig.compute_momentum(&sample.loop_moms, &three_externals))
            .collect()
    }

    pub fn to_massless_emr<T: FloatLike>(&self, sample: &BareSample<T>) -> Vec<FourMomentum<F<T>>> {
        self.edge_signatures
            .iter()
            .map(|sig| {
                sig.compute_four_momentum_from_three(&sample.loop_moms, &sample.external_moms)
            })
            .collect()
    }

    pub fn pattern(&self, edge_id: usize) -> Pattern {
        let signature = self.edge_signatures[edge_id].clone();

        let mut atom = Atom::new_num(0);

        for (i, sign) in signature.internal.into_iter().enumerate() {
            let k = sign * Atom::parse(&format!("K({},x_)", i)).unwrap();

            atom = &atom + &k;
        }

        for (i, sign) in signature.external.into_iter().enumerate() {
            let p = sign * Atom::parse(&format!("P({},x_)", i)).unwrap();
            atom = &atom + &p;
        }

        atom.to_pattern()
    }

    pub fn set_edge_signatures(&mut self, graph: &BareGraph) -> Result<(), Report> {
        // Initialize signature
        self.edge_signatures = vec![
            LoopExtSignature {
                internal: Signature(vec![SignOrZero::Zero; self.basis.len()]),
                external: Signature(vec![SignOrZero::Zero; graph.external_edges.len()])
            };
            graph.edges.len()
        ];

        // Build the adjacency list excluding vetoed edges
        let mut adj_list: HashMap<usize, Vec<(usize, usize, bool)>> = HashMap::new();
        for (edge_index, edge) in graph.edges.iter().enumerate() {
            if self.basis.contains(&edge_index) {
                continue;
            }
            let (u, v) = (edge.vertices[0], edge.vertices[1]);

            // Original orientation
            adj_list.entry(u).or_default().push((v, edge_index, false));
            // Flipped orientation
            adj_list.entry(v).or_default().push((u, edge_index, true));
        }

        // Route internal LMB momenta
        for (i_lmb, lmb_edge_id) in self.basis.iter().enumerate() {
            let edge = &graph.edges[*lmb_edge_id];
            let (u, v) = (edge.vertices[0], edge.vertices[1]);

            self.edge_signatures[*lmb_edge_id].internal.0[i_lmb] = SignOrZero::Plus;
            if let Some(path) = self.find_shortest_path(&adj_list, v, u) {
                for (edge_index, is_flipped) in path {
                    if self.edge_signatures[edge_index].internal.0[i_lmb] != SignOrZero::Zero {
                        return Err(eyre!(
                            "Inconsitency in edge momentum lmb signature assignment."
                        ));
                    }
                    self.edge_signatures[edge_index].internal.0[i_lmb] = if is_flipped {
                        SignOrZero::Minus
                    } else {
                        SignOrZero::Plus
                    };
                }
            } else {
                return Err(eyre!(
                    "No path found between vertices {} and {} for LMB: {:?}",
                    u,
                    v,
                    self.basis
                ));
            }
        }

        let sink_node = if let Some(last_external) = graph.external_edges.last() {
            match graph.edges[*last_external].edge_type {
                EdgeType::Outgoing => graph.edges[*last_external].vertices[1],
                EdgeType::Incoming => graph.edges[*last_external].vertices[0],
                _ => {
                    return Err(eyre!(
                        "External edge {} is not incoming or outgoing.",
                        graph.edges[*last_external].name
                    ))
                }
            }
        } else {
            0
        };

        // Route external momenta
        if graph.external_edges.len() >= 2 {
            for i_ext in 0..=(graph.external_edges.len() - 2) {
                let external_edge_index = graph.external_edges[i_ext];
                let external_edge = &graph.edges[external_edge_index];
                let (u, v) = match external_edge.edge_type {
                    EdgeType::Outgoing => (sink_node, external_edge.vertices[1]),
                    EdgeType::Incoming => (external_edge.vertices[0], sink_node),
                    _ => {
                        return Err(eyre!(
                            "External edge {} is not incoming or outgoing.",
                            external_edge.name
                        ))
                    }
                };

                if let Some(path) = self.find_shortest_path(&adj_list, u, v) {
                    //println!("External path from {}->{}: {} {:?}", u, v, i_ext, path);
                    for (edge_index, is_flipped) in path {
                        if self.edge_signatures[edge_index].external.0[i_ext] != SignOrZero::Zero {
                            return Err(eyre!(
                                "Inconsitency in edge momentum signature assignment."
                            ));
                        }
                        self.edge_signatures[edge_index].external.0[i_ext] = if is_flipped {
                            SignOrZero::Minus
                        } else {
                            SignOrZero::Plus
                        };
                    }
                } else {
                    return Err(eyre!(
                        "No path found between vertices {} and {} for LMB: {:?}",
                        u,
                        v,
                        self.basis
                    ));
                }
                if self.edge_signatures[external_edge_index].external.0[i_ext] != SignOrZero::Plus {
                    return Err(eyre!(
                        "Inconsitency in edge momentum external signature assignment."
                    ));
                }
            }
        }
        Ok(())
    }

    fn find_shortest_path(
        &self,
        adjacency_list: &HashMap<usize, Vec<(usize, usize, bool)>>,
        start: usize,
        end: usize,
    ) -> Option<Vec<(usize, bool)>> {
        if start == end {
            return Some(vec![]);
        }

        // Initialize BFS
        let mut queue = VecDeque::new();
        let mut visited: HashMap<usize, Option<(usize, usize, bool)>> = HashMap::new();

        queue.push_back(start);
        visited.insert(start, None);

        // Perform BFS
        while let Some(u) = queue.pop_front() {
            if u == end {
                break;
            }
            if let Some(neighbors) = adjacency_list.get(&u) {
                for &(v, edge_index, is_flipped) in neighbors {
                    #[allow(clippy::map_entry)]
                    if !visited.contains_key(&v) {
                        visited.insert(v, Some((u, edge_index, is_flipped)));
                        queue.push_back(v);
                    }
                }
            }
        }

        // Reconstruct the path if end is reached
        if !visited.contains_key(&end) {
            return None;
        }

        let mut path = Vec::new();
        let mut current = end;

        while let Some(Some((prev, edge_index, is_flipped))) = visited.get(&current) {
            path.push((*edge_index, *is_flipped));
            current = *prev;
        }

        path.reverse();
        Some(path)
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SerializableLoopMomentumBasis {
    pub basis: Vec<usize>,
    pub edge_signatures: Vec<(Vec<isize>, Vec<isize>)>,
}

// helper enum for different ways of specifying an lmb to compute stuff in.
pub enum LoopMomentumBasisSpecification<'a> {
    Literal(&'a LoopMomentumBasis),
    FromList(usize),
}

impl<'a> LoopMomentumBasisSpecification<'a> {
    pub fn basis(&self, graph: &'a Graph) -> &'a LoopMomentumBasis {
        match self {
            LoopMomentumBasisSpecification::Literal(basis) => basis,
            LoopMomentumBasisSpecification::FromList(idx) => &graph
                .derived_data
                .as_ref()
                .unwrap()
                .loop_momentum_bases
                .as_ref()
                .unwrap_or_else(|| panic!("Loop momentum bases not yet generated"))[*idx],
        }
    }

    pub fn basis_from_derived<S: NumeratorState>(
        &self,
        derived: &'a DerivedGraphData<S>,
    ) -> &'a LoopMomentumBasis {
        match self {
            LoopMomentumBasisSpecification::Literal(basis) => basis,
            LoopMomentumBasisSpecification::FromList(idx) => &derived
                .loop_momentum_bases
                .as_ref()
                .unwrap_or_else(|| panic!("Loop momentum bases not yet generated"))[*idx],
        }
    }
}

impl From<&BareGraph> for HedgeGraph<usize, usize> {
    fn from(value: &BareGraph) -> Self {
        let mut builder = HedgeGraphBuilder::new();
        let mut map = HashMap::new();

        for (i, _) in value.vertices.iter().enumerate() {
            map.insert(i, builder.add_node(i));
        }

        for (i, edge) in value.edges.iter().enumerate() {
            let source = map[&edge.vertices[0]];
            let sink = map[&edge.vertices[1]];
            builder.add_edge(source, sink, i, true);
        }

        builder.into()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn model_sm() -> Model {
        Model::from_file(String::from("src/test_resources/gammaloop_models/sm.yaml")).unwrap()
    }

    #[test]
    fn vertex_rule() {
        let model = model_sm();

        let v = &model.vertex_rules[model.vertex_rule_name_to_position["V_44"]];

        let _vertex = Vertex {
            name: "v".into(),
            vertex_info: VertexInfo::InteractonVertexInfo(InteractionVertexInfo {
                vertex_rule: v.clone(),
            }),
            edges: vec![2, 3],
        };

        let _lorentz = v.lorentz_structures[0].structure.clone();
        // println!("{}");
    }
}
