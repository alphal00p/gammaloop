use bitvec::vec::BitVec;
use symbolica::atom::{representation::InlineNum, Atom};

use crate::{
    graph::{
        half_edge::{subgraph::OrientedCut, HedgeGraph},
        BareGraph, DerivedGraphData, Edge, LoopMomentumBasis, Vertex,
    },
    numerator::NumeratorState,
};

pub struct Graph<S: NumeratorState> {
    multiplicity: Atom,
    underlying: HedgeGraph<Edge, Vertex>,
    loop_momentum_basis: LoopMomentumBasis,
    derived_data: DerivedGraphData<S>,
}

struct Process<S: NumeratorState> {
    pub initial_pdgs: Vec<i64>,
    pub final_pdgs: Vec<i64>,
    pub collection: ProcessCollection<S>,
}

enum ProcessCollection<S: NumeratorState> {
    Amplitude(Vec<Amplitude<S>>),
    CrossSection(Vec<CrossSection<S>>),
}

struct Amplitude<S: NumeratorState> {
    graphs: Vec<Graph<S>>,
}
struct CrossSection<S: NumeratorState> {
    supergraphs: Vec<CrossSectionGraph<S>>,
}

pub struct CrossSectionGraph<S: NumeratorState> {
    graph: Graph<S>,
    cuts: Vec<CrossSectionCut>,
}

pub struct CrossSectionCut {
    cut: OrientedCut,
    left: BitVec,
    right: BitVec,
}
