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
use log::debug;
use pathfinding::prelude::BfsReachable;
use serde::{Deserialize, Serialize};
use spenso::{
    algebra::ScalarMul,
    network::parsing::ShadowedStructure,
    shadowing::symbolica_utils::SerializableAtom,
    structure::{
        dimension::Dimension,
        representation::{Minkowski, RepName},
        HasStructure, NamedStructure, OrderedStructure, ToSymbolic,
    },
    tensors::parametric::{
        atomcore::{PatternReplacement, TensorAtomMaps},
        ParamTensor,
    },
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

pub(crate) fn spenso_lor(
    tag: i32,
    ind: impl Into<Aind>,
    dim: impl Into<Dimension>,
) -> ShadowedStructure<Aind> {
    let mink = Minkowski {}.new_slot(dim, ind);
    NamedStructure::from_iter([mink], GS.emr_mom, Some(vec![Atom::num(tag)])).structure
}

pub(crate) fn spenso_lor_atom(tag: i32, ind: impl Into<Aind>, dim: impl Into<Dimension>) -> Atom {
    spenso_lor(tag, ind, dim).to_symbolic(None).unwrap()
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct IntegrandExpr {
    integrand: Atom,
    // add_arg: Option<Atom>,
}

pub(crate) fn is_not_paired(pair: &HedgePair) -> bool {
    !pair.is_paired()
}

pub mod uv_graph;
pub use uv_graph::UltravioletGraph;

pub mod poset;
pub use poset::Poset;

pub mod wood;
pub use wood::Wood;

pub mod approx;
pub use approx::ApproxOp;

pub mod forest;
pub use forest::Forest;

#[cfg(test)]
mod tests;
