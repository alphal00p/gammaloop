use crate::{
    cff::{
        esurface::{
            generate_esurface_data, get_existing_esurfaces, EsurfaceDerivedData, ExistingEsurfaces,
            ExternalShift,
        },
        expression::CFFExpression,
        generation::generate_cff_expression,
    },
    gammaloop_integrand::{BareSample, DefaultSample},
    ltd::{generate_ltd_expression, LTDExpression},
    model::{self, EdgeSlots, Model, Particle, VertexSlots},
    momentum::{FourMomentum, Polarization, Rotation, Signature, ThreeMomentum},
    numerator::{
        AppliedFeynmanRule, AtomStructure, ContractionSettings, Evaluate, Evaluators, ExtraInfo,
        GammaAlgebraMode, Numerator, NumeratorState, NumeratorStateError, PythonState,
        RepeatingIteratorTensorOrScalar, TypedNumeratorState, UnInit,
    },
    subtraction::{
        overlap::{find_maximal_overlap, OverlapStructure},
        static_counterterm::{self, CounterTerm},
    },
    utils::{self, sorted_vectorize, FloatLike, F},
    ExportSettings, Settings, TropicalSubgraphTableSettings,
};

use ahash::{HashSet, RandomState};

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
    data::{DataTensor, DenseTensor, GetTensorData, SetTensorData, SparseTensor},
    structure::{
        AbstractIndex, BaseRepName, Euclidean, HasStructure, Lorentz, NamedStructure, PhysReps,
        Representation, ScalarTensor, ToSymbolic, VecStructure, COLORADJ, COLORANTIFUND,
        COLORANTISEXT, COLORFUND, COLORSEXT, EUCLIDEAN,
    },
    ufo::{preprocess_ufo_color_wrapped, preprocess_ufo_spin_wrapped},
};
use uuid::Uuid;

use core::panic;
use serde::{de::DeserializeOwned, Deserialize, Serialize};
use smallvec::{smallvec, SmallVec};
use smartstring::{LazyCompact, SmartString};
use std::{
    collections::HashMap,
    fmt::{Display, Formatter},
    ops::{AddAssign, Neg, Not},
    path::{Path, PathBuf},
    str::FromStr,
    sync::Arc,
};

use symbolica::{
    atom::Atom,
    domains::{float::NumericalFloatLike, rational::Rational},
    id::{Pattern, Replacement},
};
//use symbolica::{atom::Symbol,state::State};

use constcat::concat;

const MAX_COLOR_INNER_CONTRACTIONS: usize = 3;

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
}

#[derive(Debug, Clone)]
#[enum_dispatch(HasVertexInfo)]
pub enum VertexInfo {
    ExternalVertexInfo(ExternalVertexInfo),
    InteractonVertexInfo(InteractionVertexInfo),
}

impl VertexInfo {
    pub fn generate_vertex_slots(
        &self,
        shifts: (usize, usize, usize),
    ) -> (VertexSlots, (usize, usize, usize)) {
        match self {
            VertexInfo::ExternalVertexInfo(e) => {
                let (e, shifts) = e.particle.slots(shifts);
                (e.into(), shifts)
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
    direction: EdgeType,
    particle: Arc<model::Particle>,
}

impl HasVertexInfo for ExternalVertexInfo {
    fn get_type(&self) -> SmartString<LazyCompact> {
        SmartString::<LazyCompact>::from("external_vertex_info")
    }
}

#[derive(Debug, Clone)]
pub struct InteractionVertexInfo {
    #[allow(unused)]
    pub vertex_rule: Arc<model::VertexRule>,
}

impl InteractionVertexInfo {
    pub fn apply_vertex_rule(
        &self,
        edges: &[isize],
        vertex_pos: usize,
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
                        Pattern::parse(&format!("-Q({},aind(lord(4,indexid(x_))))", -e))
                            .unwrap()
                            .into() //TODO flip based on flow
                    } else {
                        Pattern::parse(&format!("Q({},aind(loru(4,indexid(x_))))", e))
                            .unwrap()
                            .into() //TODO flip based on flow
                    };

                    atom = momentum_in_pattern.replace_all(
                        atom.as_view(),
                        &momentum_out_pattern,
                        None,
                        None,
                    );
                }

                atom = preprocess_ufo_spin_wrapped(atom, true);

                for (i, _) in edges.iter().enumerate() {
                    let replacements = vertex_slots[i].replacements(i + 1);

                    let reps: Vec<Replacement> = replacements
                        .iter()
                        .map(|(pat, rhs)| Replacement::new(pat, rhs))
                        .collect();

                    atom = atom.replace_all_multiple(&reps);
                }

                atom
            })
            .collect_vec();

        let color_structure: Vec<Atom> = self
            .vertex_rule
            .color_structures
            .iter()
            .map(|cs| {
                let mut atom = cs.clone();
                //a is adjoint index, i is fundamental index, ia is antifundamental index, s is sextet ,sa is antisextet

                atom = preprocess_ufo_color_wrapped(atom);
                //T(1,2,3) Fundamental representation matrix (T a1 )ı  ̄3 i2
                // f(1,2,3) Antisymmetric structure constant f a1a2a3
                // d(1,2,3) Symmetric structure constant da1 a2 a3
                // Epsilon(1,2,3) Fundamental Levi-Civita tensor εi1 i2 i3 EpsilonBar(1,2,3) Antifundamental Levi-Civita tensor εı  ̄1 ı  ̄2 ı  ̄3
                // T6(1,2,3) Sextet representation matrix (T a1 6 ) ̄ α3 α2
                // K6(1,2,3) Sextet Clebsch-Gordan coefficient (K6)ı  ̄2 ı  ̄3 α1 K6Bar(1,2,3) Antisextet Clebsch-Gordan coefficient (K6) ̄ α1 i2i3

                // First process kronkers, with respect to spin:

                let spins: Vec<isize> =
                    self.vertex_rule.particles.iter().map(|s| s.color).collect();

                for (i, s) in spins.iter().enumerate() {
                    let id1 = Pattern::parse(&format!("Identity({},x_)", i + 1)).unwrap();

                    let id2 = Pattern::parse(&format!("id(aind(x_,{}))", i + 1)).unwrap();

                    let ind = match s {
                        1 => concat!(EUCLIDEAN, "(1,"),
                        3 => concat!(COLORFUND, "(3,"),
                        -3 => concat!(COLORANTIFUND, "(3,"),
                        6 => concat!(COLORSEXT, "(6,"),
                        -6 => concat!(COLORANTISEXT, "(6,"),
                        8 => concat!(COLORADJ, "(8,"),
                        i => panic!("Color {i}not supported "),
                    };

                    atom = id1.replace_all(
                        atom.as_view(),
                        &Pattern::parse(&format!("id(aind({}indexid({})),x_))", ind, i + 1))
                            .unwrap()
                            .into(),
                        None,
                        None,
                    );

                    atom = id2.replace_all(
                        atom.as_view(),
                        &Pattern::parse(&format!("id(aind(x_,{}indexid({}))))", ind, i + 1))
                            .unwrap()
                            .into(),
                        None,
                        None,
                    );
                }

                for (i, _) in edges.iter().enumerate() {
                    let replacements = vertex_slots[i].replacements(i + 1);

                    let reps: Vec<Replacement> = replacements
                        .iter()
                        .map(|(pat, rhs)| Replacement::new(pat, rhs))
                        .collect();

                    atom = atom.replace_all_multiple(&reps);
                }

                for i in 0..MAX_COLOR_INNER_CONTRACTIONS {
                    let pat: Pattern = Atom::parse(&format!("indexid({})", -1 - i as i64))
                        .unwrap()
                        .into_pattern();

                    atom = pat.replace_all(
                        atom.as_view(),
                        &Atom::new_num(i as i64).into_pattern().into(),
                        None,
                        None,
                    );
                }

                atom
            })
            .collect();

        let irep: Representation<PhysReps> =
            Euclidean::new_dimed_rep_selfless(color_structure.len()).cast();
        let i = irep.new_slot(AbstractIndex::try_from(format!("i{}", vertex_pos)).unwrap());

        let color_structure = DataTensor::Dense(
            DenseTensor::from_data(color_structure, VecStructure::from(vec![i])).unwrap(),
        );

        let jrep: Representation<PhysReps> =
            Euclidean::new_dimed_rep_selfless(spin_structure.len()).cast();

        let j = jrep.new_slot(AbstractIndex::try_from(format!("j{}", vertex_pos)).unwrap());

        let spin_structure = DataTensor::Dense(
            DenseTensor::from_data(spin_structure, VecStructure::from(vec![j])).unwrap(),
        );

        let mut couplings: DataTensor<Atom> =
            DataTensor::Sparse(SparseTensor::empty(VecStructure::from(vec![i, j])));

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

impl HasVertexInfo for InteractionVertexInfo {
    fn get_type(&self) -> SmartString<LazyCompact> {
        SmartString::<LazyCompact>::from("interacton_vertex_info")
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
}

impl Edge {
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
        }
    }

    pub fn is_incoming_to(&self, vertex: usize) -> bool {
        self.vertices[1] == vertex
    }

    pub fn denominator(&self, graph: &BareGraph) -> (Atom, Atom) {
        let num = *graph.edge_name_to_position.get(&self.name).unwrap();
        let mom = Atom::parse(&format!("Q{num}")).unwrap();
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
        let mom_rep = lmb.pattern(num).into();
        atom.replace_all(&mom, &mom_rep, None, None)
    }

    // pub fn edge_momentum_symbol(&self, graph: &Graph) -> Symbol {
    //     let num = *graph.edge_name_to_position.get(&self.name).unwrap();
    //     State::get_symbol(format!("Q{num}"))
    // }

    pub fn in_slot(&self, graph: &BareGraph) -> EdgeSlots<Lorentz> {
        let local_pos_in_sink_vertex =
            graph.vertices[self.vertices[0]].get_local_edge_position(self, graph);

        graph.vertex_slots[self.vertices[0]][local_pos_in_sink_vertex].dual()
    }

    pub fn out_slot(&self, graph: &BareGraph) -> EdgeSlots<Lorentz> {
        let local_pos_in_sink_vertex =
            graph.vertices[self.vertices[1]].get_local_edge_position(self, graph);

        graph.vertex_slots[self.vertices[1]][local_pos_in_sink_vertex].dual()
    }

    pub fn numerator(&self, graph: &BareGraph) -> (Atom, usize) {
        let num = *graph.edge_name_to_position.get(&self.name).unwrap();
        let in_slots = self.in_slot(graph);
        let out_slots = self.out_slot(graph);

        match self.edge_type {
            EdgeType::Incoming => (self.particle.incoming_polarization_atom(&out_slots, num), 0),
            EdgeType::Outgoing => (self.particle.outgoing_polarization_atom(&in_slots, num), 0),
            EdgeType::Virtual => {
                let mut atom = self.propagator.numerator.clone();

                let pfun = Pattern::parse("P(x_)").unwrap();
                if self.particle.is_antiparticle() {
                    atom = pfun.replace_all(
                        atom.as_view(),
                        &Pattern::parse(&format!("-Q({},aind(loru(4,x_)))", num))
                            .unwrap()
                            .into(),
                        None,
                        None,
                    );
                } else {
                    atom = pfun.replace_all(
                        atom.as_view(),
                        &Pattern::parse(&format!("Q({},aind(loru(4,x_)))", num))
                            .unwrap()
                            .into(),
                        None,
                        None,
                    );
                }

                let pslashfun = Pattern::parse("PSlash(i_,j_)").unwrap();
                let pindex_num = graph.shifts.0 + 1;
                if self.particle.is_antiparticle() {
                    atom = pslashfun.replace_all(
                        atom.as_view(),
                        &Pattern::parse(&format!(
                            "-Q({},aind(lord(4,{})))*Gamma({},i_,j_)",
                            num, pindex_num, pindex_num
                        ))
                        .unwrap()
                        .into(),
                        None,
                        None,
                    );
                } else {
                    atom = pslashfun.replace_all(
                        atom.as_view(),
                        &Pattern::parse(&format!(
                            "Q({},aind(lord(4,{})))*Gamma({},i_,j_)",
                            num, pindex_num, pindex_num
                        ))
                        .unwrap()
                        .into(),
                        None,
                        None,
                    );
                }

                atom = preprocess_ufo_spin_wrapped(atom, false);

                let (replacements_in, mut replacements_out) = if self.particle.is_antiparticle() {
                    (in_slots.replacements(2), out_slots.replacements(1))
                } else {
                    (in_slots.replacements(1), out_slots.replacements(2))
                };

                replacements_out.push((
                    Atom::parse("indexid(x_)").unwrap().into_pattern(),
                    Atom::parse("x_").unwrap().into_pattern().into(),
                ));

                for (&cin, &cout) in in_slots.color.iter().zip(out_slots.color.iter()) {
                    let id: NamedStructure<String, ()> =
                        NamedStructure::from_iter([cin, cout], "id".into(), None);
                    atom = atom * &id.to_symbolic().unwrap();
                }

                let reps: Vec<Replacement> = replacements_in
                    .iter()
                    .chain(replacements_out.iter())
                    .map(|(pat, rhs)| Replacement::new(pat, rhs))
                    .collect();

                (atom.replace_all_multiple(&reps), pindex_num + 1)
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
    pub fn get_local_edge_position(&self, edge: &Edge, graph: &BareGraph) -> usize {
        let global_id: usize = graph.edge_name_to_position[&edge.name];
        self.edges
            .iter()
            .enumerate()
            .find(|(_, &e)| e == global_id)
            .unwrap()
            .0
    }

    pub fn generate_vertex_slots(
        &self,
        shifts: (usize, usize, usize),
    ) -> (VertexSlots, (usize, usize, usize)) {
        self.vertex_info.generate_vertex_slots(shifts)
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
        match &self.vertex_info {
            VertexInfo::ExternalVertexInfo(_) => None,
            VertexInfo::InteractonVertexInfo(interaction_vertex_info) => interaction_vertex_info
                .apply_vertex_rule(
                    &self.add_signs_to_edges(graph),
                    pos,
                    &graph.vertex_slots[pos],
                ),
        }
    }

    pub fn contracted_vertex_rule(&self, graph: &BareGraph) -> Option<Atom> {
        let all = self.apply_vertex_rule(graph)?;
        let scalar = all
            .into_iter()
            .reduce(|acc, tensor| acc.contract(&tensor).unwrap())
            .unwrap()
            .get(&[])
            .unwrap()
            .clone();

        Some(scalar)
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

            overall_factor: graph.overall_factor.clone(),
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
}

#[derive(Debug, Clone)]
pub struct BareGraph {
    pub name: SmartString<LazyCompact>,
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub external_edges: Vec<usize>,
    pub overall_factor: String,
    pub external_connections: Vec<(Option<usize>, Option<usize>)>,
    pub loop_momentum_basis: LoopMomentumBasis,
    pub vertex_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub edge_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub vertex_slots: Vec<VertexSlots>,
    pub shifts: (usize, usize, usize),
}

impl BareGraph {
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

    pub fn dot(&self) -> String {
        let mut dot = String::new();
        dot.push_str("digraph G {\n");
        for (i, edge) in self.edges.iter().enumerate() {
            let from = self.vertices[edge.vertices[0]].name.clone();
            let to = self.vertices[edge.vertices[1]].name.clone();
            dot.push_str(&format!(
                "\"{}\" -> \"{}\" [label=\"name: {} particle:{}  Q({}) {} \"];\n",
                from, to, edge.name, edge.particle.name, i, edge.edge_type
            ));
        }
        dot.push_str("}\n");
        dot
    }

    pub fn dot_lmb(&self) -> String {
        let mut dot = String::new();
        dot.push_str("digraph G {\n");
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
        dot.push_str("digraph G {\n");
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
        dot.push_str("digraph G {\n");
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
            overall_factor: graph.overall_factor.clone(),
            external_connections: vec![],
            loop_momentum_basis: LoopMomentumBasis {
                basis: vec![],
                edge_signatures: vec![],
            },
            vertex_name_to_position,
            edge_name_to_position: HashMap::default(),
            vertex_slots: vec![],
            shifts: (0, 0, MAX_COLOR_INNER_CONTRACTIONS),
        };

        // let mut edges: Vec<Edge> = vec![];
        let mut edge_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState> =
            HashMap::default();
        for edge in &graph.edges {
            let edge = Edge::from_serializable_edge(model, &g, edge);
            edge_name_to_position.insert(edge.name.clone(), g.edges.len());
            g.edges.push(edge);
        }

        debug!("Loaded {} edges", g.edges.len());
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

        g.generate_vertex_slots();

        // panic!("{:?}", g.edge_name_to_position);

        g
    }

    fn generate_vertex_slots(&mut self) {
        let (v, s) = self
            .vertices
            .iter()
            .fold((vec![], self.shifts), |(mut acc, shifts), v| {
                let (e, new_shifts) = v.generate_vertex_slots(shifts);
                acc.push(e);
                (acc, new_shifts)
            });
        self.shifts = s;
        self.vertex_slots = v;
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

    pub fn generate_lmb_replacement_rules(&self) -> Vec<(Atom, Atom)> {
        self.loop_momentum_basis_replacement_rule(&self.loop_momentum_basis)
    }

    fn loop_momentum_basis_replacement_rule(&self, lmb: &LoopMomentumBasis) -> Vec<(Atom, Atom)> {
        let mut rule = vec![];

        for (i, signature) in lmb.edge_signatures.iter().enumerate() {
            rule.push((
                Atom::parse(&format!("Q{}(x{}__)", i, i)).unwrap(),
                self.replacement_rule_from_signature(i, signature),
            ));
        }

        rule
    }

    fn replacement_rule_from_signature(&self, index: usize, signature: &LoopExtSignature) -> Atom {
        let mut acc = Atom::new_num(0);
        for (i_l, &sign) in signature.internal.iter().enumerate() {
            let k = sign * Atom::parse(&format!("K{}(x{}__)", i_l, index)).unwrap();
            acc = &acc + &k;
        }

        for (i_e, &sign) in signature.external.iter().enumerate() {
            let p = sign * Atom::parse(&format!("P{}(x{}__)", i_e, index)).unwrap();
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
        let (_, (dep_mom, _)) = external_edges
            .iter()
            .enumerate()
            .find(|(external_index, (index, _edge))| {
                self.loop_momentum_basis.edge_signatures[*index].external[*external_index]
                    .is_positive()
                    .not()
            })
            .unwrap_or_else(|| panic!("could not determine dependent momenta"));

        let dep_mom_signature = &self.loop_momentum_basis.edge_signatures[*dep_mom].external;

        let external_shift = external_edges
            .iter()
            .zip(dep_mom_signature.iter())
            .filter(|(_external_edge, dep_mom_sign)| dep_mom_sign.is_sign())
            .map(|((external_edge, _), dep_mom_sign)| (*external_edge, *dep_mom_sign as i64))
            .collect();

        (*dep_mom, external_shift)
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
        let dimension = 3;
        let num_virtual_loop_edges = self.get_loop_edges_iterator().count();
        let num_loops = self.loop_momentum_basis.basis.len();
        let target_omega = settings.target_omega;

        let weight =
            (target_omega + (dimension * num_loops) as f64 / 2.) / num_virtual_loop_edges as f64;

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
            .build_sampler(loop_part, dimension)
            .map_err(|e| eyre!(e))
    }

    pub fn load_derived_data<NumState: NumeratorState + DeserializeOwned>(
        self,
        model: &Model,
        path: &Path,
        settings: &Settings,
    ) -> Result<Graph<NumState>, Report> {
        let derived_data_path = path.join(format!("derived_data_{}.bin", self.name.as_str()));
        debug!("Loading derived data from {:?}", derived_data_path);
        let mut derived_data = DerivedGraphData::load_from_path(&derived_data_path)?;
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
        });
        Graph {
            bare_graph: self.bare_graph,
            derived_data: processed_data,
        }
    }

    pub fn apply_feynman_rules(mut self) -> Graph<AppliedFeynmanRule> {
        let processed_data = self
            .derived_data
            .map(|d| d.apply_feynman_rules(&mut self.bare_graph));
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
            settings.general.debug,
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

        if let Ok(table) = table {
            debug!("min dod: {}", table.get_smallest_dod());
            if let Some(d) = &mut self.derived_data {
                d.tropical_subgraph_table = Some(table);
            }
        } else if settings.panic_on_fail {
            panic!("Tropical subgraph table generation failed 🥥");
        } else {
            warn!("Tropical subgraph table generation failed 🥥");
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
    ) -> DataTensor<Complex<F<T>>, AtomStructure> {
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
    ) -> DataTensor<Complex<F<T>>, AtomStructure> {
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
        match std::fs::read(path) {
            Ok(derived_data_bytes) => {
                let derived_data: DerivedGraphData<S> =
                    bincode::decode_from_slice(&derived_data_bytes, bincode::config::standard())?.0;
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
    ) -> DataTensor<Complex<F<T>>, AtomStructure> {
        let emr = bare_graph
            .loop_momentum_basis
            .edge_signatures
            .iter()
            .map(|sig| sig.compute_momentum(loop_moms, external_moms))
            .collect_vec();

        let mut den = Complex::new_re(F::from_f64(1.));
        for (e, q) in bare_graph.edges.iter().zip(emr.iter()) {
            if e.edge_type == EdgeType::Virtual {
                println!("q: {}", q);
                if let Some(mass) = e.particle.mass.value {
                    let m2 = mass.norm_squared();
                    let m2: F<T> = F::from_ff64(m2);
                    den *= &q.square() - &m2;
                } else {
                    den *= q.square();
                }
            }
        }
        println!("den: {}", den);
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
    ) -> DataTensor<Complex<F<T>>, AtomStructure> {
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
    ) -> DataTensor<Complex<F<T>>, AtomStructure> {
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
    ) -> DataTensor<Complex<F<T>>, AtomStructure> {
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
            RepeatingIteratorTensorOrScalar::Tensors(mut num_iter) => {
                if let Some(i) = num_iter.next() {
                    let c = Complex::new_re(cff.next().unwrap());
                    let mut sum = i.map_data_ref(|n| n * &c);

                    for j in cff {
                        let c = Complex::new_re(j);
                        sum += num_iter.next().unwrap().map_data_ref(|n| n * &c);
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
    ) -> DataTensor<Complex<F<T>>, AtomStructure> {
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
    ) -> RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<T>>, AtomStructure>> {
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
    ) -> DataTensor<Complex<F<T>>, AtomStructure> {
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

    pub fn process_numerator(
        self,
        base_graph: &mut BareGraph,
        model: &Model,
        contraction_settings: ContractionSettings<Rational>,
        export_path: PathBuf,
        export_settings: &ExportSettings,
    ) -> DerivedGraphData<Evaluators> {
        let extra_info = self.generate_extra_info(export_path);

        let color_simplified =
            if let Some(global) = &export_settings.numerator_settings.global_numerator {
                debug!("Using global numerator: {}", global);
                let global = Atom::parse(global).unwrap();
                self.map_numerator(|n| {
                    n.from_global(global, base_graph)
                        .color_symplify()
                        .color_project()
                })
            } else {
                self.map_numerator(|n| n.from_graph(base_graph).color_symplify().color_project())
            };

        let parsed = match &export_settings.numerator_settings.gamma_algebra {
            GammaAlgebraMode::Symbolic => {
                color_simplified.map_numerator(|n| n.gamma_symplify().parse())
            }
            GammaAlgebraMode::Concrete => color_simplified.map_numerator(|n| n.parse()),
        };

        parsed
            .map_numerator_res(|n| {
                Result::<_, Report>::Ok(n.contract(contraction_settings)?.generate_evaluators(
                    model,
                    base_graph,
                    &extra_info,
                    export_settings,
                ))
            })
            .unwrap()
    }

    fn apply_feynman_rules(
        self,
        base_graph: &mut BareGraph,
    ) -> DerivedGraphData<AppliedFeynmanRule> {
        self.map_numerator(|n| n.from_graph(base_graph))
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

    pub fn load_from_path(path: &Path) -> Result<Self, Report> {
        match std::fs::read(path) {
            Ok(derived_data_bytes) => {
                let derived_data: Self =
                    bincode::decode_from_slice(&derived_data_bytes, bincode::config::standard())?.0;
                Ok(derived_data)
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

        atom.into_pattern()
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

#[cfg(test)]
mod tests {
    use std::fs::File;

    use crate::cross_section::OutputMetaData;

    use super::*;

    fn model_sm() -> Model {
        let path = Path::new("./src/test_resources/lbl/");
        let output_meta_data: OutputMetaData =
            serde_yaml::from_reader(File::open(path.join("output_metadata.yaml")).unwrap())
                .unwrap();
        Model::from_file(String::from(
            path.join(format!(
                "sources/model/{}.yaml",
                output_meta_data.model_name
            ))
            .to_str()
            .unwrap(),
        ))
        .unwrap()
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
