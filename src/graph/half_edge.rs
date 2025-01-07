use std::num::TryFromIntError;
use std::ops::{Index, IndexMut, Neg};
use std::{collections::VecDeque, fmt::Display, hash::Hash};

use ahash::{AHashMap, AHashSet};
use bitvec::{slice::IterOnes, vec::BitVec};
use itertools::Itertools;
use rand::{rngs::SmallRng, Rng, SeedableRng};
use serde::{Deserialize, Serialize};

use bitvec::prelude::*;
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub struct NodeIndex(pub usize);

impl std::fmt::Display for NodeIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

pub struct PowersetIterator {
    size: usize,
    current: usize,
}

impl PowersetIterator {
    pub fn new(n_elements: u8) -> Self {
        PowersetIterator {
            size: 1 << n_elements,
            current: 0,
        }
    }
}

impl Iterator for PowersetIterator {
    type Item = BitVec;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current < self.size {
            let out = BitVec::<_, Lsb0>::from_element(self.current);

            self.current += 1;
            Some(out)
        } else {
            None
        }
    }
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Hedge(pub usize);

impl Hedge {
    pub fn to_edge_id<N, E>(self, involution: &Involution<N, E>) -> EdgeId {
        EdgeId::from_half_edge(self, involution)
    }

    pub fn to_edge_id_with_subgraph<N, E, S: SubGraph>(
        self,
        involution: &Involution<N, E>,
        subgraph: &S,
    ) -> Option<EdgeId> {
        EdgeId::from_half_edge_with_subgraph(self, involution, subgraph)
    }
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub enum EdgeId {
    Unpaired {
        hedge: Hedge,
        flow: Flow,
    },
    Paired {
        source: Hedge,
        sink: Hedge,
    },
    Split {
        source: Hedge,
        sink: Hedge,
        split: Flow,
    },
}

impl EdgeId {
    pub fn from_half_edge<N, E>(hedge: Hedge, involution: &Involution<N, E>) -> Self {
        match involution[hedge] {
            InvolutiveMapping::Source { sink_idx, .. } => EdgeId::Paired {
                source: hedge,
                sink: sink_idx,
            },
            InvolutiveMapping::Sink { source_idx } => EdgeId::Paired {
                source: source_idx,
                sink: hedge,
            },
            InvolutiveMapping::Identity { underlying, .. } => EdgeId::Unpaired {
                hedge,
                flow: underlying,
            },
        }
    }

    pub fn from_source<N, E>(hedge: Hedge, involution: &Involution<N, E>) -> Option<Self> {
        match involution[hedge] {
            InvolutiveMapping::Source { sink_idx, .. } => Some(EdgeId::Paired {
                source: hedge,
                sink: sink_idx,
            }),
            InvolutiveMapping::Identity { underlying, .. } => Some(EdgeId::Unpaired {
                hedge,
                flow: underlying,
            }),
            _ => None,
        }
    }

    pub fn from_source_with_subgraph<N, E, S: SubGraph>(
        hedge: Hedge,
        involution: &Involution<N, E>,
        subgraph: &S,
    ) -> Option<Self> {
        if subgraph.includes(&hedge) {
            match involution[hedge] {
                InvolutiveMapping::Source { sink_idx, .. } => {
                    if subgraph.includes(&sink_idx) {
                        Some(EdgeId::Paired {
                            source: hedge,
                            sink: sink_idx,
                        })
                    } else {
                        Some(EdgeId::Split {
                            source: hedge,
                            sink: sink_idx,
                            split: Flow::Source,
                        })
                    }
                }
                InvolutiveMapping::Sink { source_idx } => {
                    if subgraph.includes(&source_idx) {
                        None
                    } else {
                        Some(EdgeId::Split {
                            source: source_idx,
                            sink: hedge,
                            split: Flow::Sink,
                        })
                    }
                }
                InvolutiveMapping::Identity { underlying, .. } => Some(EdgeId::Unpaired {
                    hedge,
                    flow: underlying,
                }),
            }
        } else {
            None
        }
    }

    pub fn from_half_edge_with_subgraph<N, E, S: SubGraph>(
        hedge: Hedge,
        involution: &Involution<N, E>,
        subgraph: &S,
    ) -> Option<Self> {
        if subgraph.includes(&hedge) {
            Some(match involution[hedge] {
                InvolutiveMapping::Source { sink_idx, .. } => {
                    if subgraph.includes(&sink_idx) {
                        EdgeId::Paired {
                            source: hedge,
                            sink: sink_idx,
                        }
                    } else {
                        EdgeId::Split {
                            source: hedge,
                            sink: sink_idx,
                            split: Flow::Source,
                        }
                    }
                }
                InvolutiveMapping::Sink { source_idx } => {
                    if subgraph.includes(&source_idx) {
                        EdgeId::Paired {
                            source: source_idx,
                            sink: hedge,
                        }
                    } else {
                        EdgeId::Split {
                            source: source_idx,
                            sink: hedge,
                            split: Flow::Sink,
                        }
                    }
                }
                InvolutiveMapping::Identity { underlying, .. } => EdgeId::Unpaired {
                    hedge,
                    flow: underlying,
                },
            })
        } else {
            None
        }
    }
}

impl Display for Hedge {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum InvolutiveMapping<E> {
    Identity { data: EdgeData<E>, underlying: Flow },
    Source { data: EdgeData<E>, sink_idx: Hedge },
    Sink { source_idx: Hedge },
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EdgeData<E> {
    pub orientation: Orientation,
    pub data: Option<E>,
}

impl<E> From<Orientation> for EdgeData<E> {
    fn from(orientation: Orientation) -> Self {
        EdgeData {
            orientation,
            data: None,
        }
    }
}

impl<E> EdgeData<E> {
    pub fn un_orient(&mut self) {
        self.orientation = Orientation::Undirected;
    }

    pub fn orient(&mut self) {
        if let Orientation::Undirected = self.orientation {
            self.orientation = Orientation::Default;
        }
    }
    pub fn flip_orientation(&mut self) {
        self.orientation = self.orientation.reverse();
    }
    pub fn new(data: E, orientation: Orientation) -> Self {
        EdgeData {
            data: Some(data),
            orientation,
        }
    }
    pub fn take(&mut self) -> Self {
        EdgeData {
            orientation: self.orientation,
            data: self.data.take(),
        }
    }
    pub const fn none() -> Self {
        Self {
            orientation: Orientation::Undirected,
            data: None,
        }
    }
    pub fn map<E2>(self, f: impl FnOnce(E) -> E2) -> EdgeData<E2> {
        EdgeData {
            orientation: self.orientation,
            data: self.data.map(f),
        }
    }

    pub fn and_then<U, F>(self, f: F) -> EdgeData<U>
    where
        F: FnOnce(E) -> Option<U>,
    {
        match self.data {
            Some(x) => EdgeData {
                data: f(x),
                orientation: self.orientation,
            },
            None => EdgeData {
                data: None,
                orientation: self.orientation,
            },
        }
    }

    pub const fn as_ref(&self) -> EdgeData<&E> {
        EdgeData {
            orientation: self.orientation,
            data: self.data.as_ref(),
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, Copy, PartialEq, Eq, Hash)]
pub enum Orientation {
    Default,  //incoming for externals
    Reversed, //outgoing for externals
    Undirected,
}

#[derive(Clone, Debug, Serialize, Deserialize, Copy, PartialEq, Eq, Hash)]
pub enum Flow {
    Source,
    Sink,
}

impl Neg for Flow {
    type Output = Flow;

    fn neg(self) -> Self::Output {
        match self {
            Flow::Source => Flow::Sink,
            Flow::Sink => Flow::Source,
        }
    }
}

impl From<bool> for Flow {
    fn from(value: bool) -> Self {
        if value {
            Flow::Source
        } else {
            Flow::Sink
        }
    }
}

impl From<Flow> for Orientation {
    fn from(value: Flow) -> Self {
        match value {
            Flow::Source => Orientation::Default,
            Flow::Sink => Orientation::Reversed,
        }
    }
}

impl From<Flow> for SignOrZero {
    fn from(value: Flow) -> Self {
        match value {
            Flow::Source => SignOrZero::Plus,
            Flow::Sink => SignOrZero::Minus,
        }
    }
}

impl TryFrom<Orientation> for Flow {
    type Error = &'static str;

    fn try_from(value: Orientation) -> Result<Self, Self::Error> {
        match value {
            Orientation::Default => Ok(Flow::Source),
            Orientation::Reversed => Ok(Flow::Sink),
            Orientation::Undirected => Err("Cannot convert undirected orientation to flow"),
        }
    }
}

impl Orientation {
    pub fn reverse(self) -> Orientation {
        match self {
            Orientation::Default => Orientation::Reversed,
            Orientation::Reversed => Orientation::Default,
            Orientation::Undirected => Orientation::Undirected,
        }
    }

    pub fn relative_to(&self, other: Flow) -> Orientation {
        match (self, other) {
            (Orientation::Default, Flow::Source) => Orientation::Default,
            (Orientation::Default, Flow::Sink) => Orientation::Reversed,
            (Orientation::Reversed, Flow::Source) => Orientation::Reversed,
            (Orientation::Reversed, Flow::Sink) => Orientation::Default,
            (Orientation::Undirected, _) => Orientation::Undirected,
        }
    }
}

impl From<bool> for Orientation {
    fn from(value: bool) -> Self {
        if value {
            Orientation::Default
        } else {
            Orientation::Undirected
        }
    }
}

impl From<SignOrZero> for Orientation {
    fn from(value: SignOrZero) -> Self {
        match value {
            SignOrZero::Zero => Orientation::Undirected,
            SignOrZero::Plus => Orientation::Default,
            SignOrZero::Minus => Orientation::Reversed,
        }
    }
}

impl From<Orientation> for SignOrZero {
    fn from(value: Orientation) -> Self {
        match value {
            Orientation::Default => SignOrZero::Plus,
            Orientation::Reversed => SignOrZero::Minus,
            Orientation::Undirected => SignOrZero::Zero,
        }
    }
}

impl<E> InvolutiveMapping<E> {
    pub fn is_identity(&self) -> bool {
        matches!(self, InvolutiveMapping::Identity { .. })
    }

    pub fn to_sign(&self) -> SignOrZero {
        match self {
            InvolutiveMapping::Identity { .. } => SignOrZero::Zero,
            InvolutiveMapping::Sink { .. } => SignOrZero::Minus,
            InvolutiveMapping::Source { .. } => SignOrZero::Plus,
        }
    }

    pub fn map_data_ref<E2>(&self, f: impl FnOnce(&E) -> E2) -> InvolutiveMapping<E2> {
        match self {
            InvolutiveMapping::Identity { data, underlying } => InvolutiveMapping::Identity {
                data: data.as_ref().map(f),
                underlying: *underlying,
            },
            InvolutiveMapping::Source { data, sink_idx } => InvolutiveMapping::Source {
                data: data.as_ref().map(f),
                sink_idx: *sink_idx,
            },
            InvolutiveMapping::Sink { source_idx } => InvolutiveMapping::Sink {
                source_idx: *source_idx,
            },
        }
    }

    pub fn map_data<E2>(self, f: impl FnOnce(E) -> E2) -> InvolutiveMapping<E2> {
        match self {
            InvolutiveMapping::Identity { data, underlying } => InvolutiveMapping::Identity {
                data: data.map(f),
                underlying,
            },
            InvolutiveMapping::Source { data, sink_idx } => InvolutiveMapping::Source {
                data: data.map(f),
                sink_idx,
            },
            InvolutiveMapping::Sink { source_idx } => InvolutiveMapping::Sink { source_idx },
        }
    }

    // pub fn wrap_fn_mut<E2>(
    //     f: impl FnMut(&E) -> E2,
    // ) -> impl FnMut(&InvolutiveMapping<E>) -> InvolutiveMapping<E2> {
    //     move |mapping| match mapping {
    //         InvolutiveMapping::Identity { data, underlying } => InvolutiveMapping::Identity {
    //             data: data.as_ref().map(f),
    //             underlying: *underlying,
    //         },
    //         InvolutiveMapping::Source { data, sink_idx } => InvolutiveMapping::Source {
    //             data: data.as_ref().map(f),
    //             sink_idx: *sink_idx,
    //         },
    //         InvolutiveMapping::Sink { source_idx } => InvolutiveMapping::Sink {
    //             source_idx: *source_idx,
    //         },
    //     }
    // }

    pub fn map_data_none<E2>(&self) -> InvolutiveMapping<E2> {
        match self {
            InvolutiveMapping::Identity { underlying, .. } => InvolutiveMapping::Identity {
                data: EdgeData::none(),
                underlying: *underlying,
            },
            InvolutiveMapping::Source { sink_idx, .. } => InvolutiveMapping::Source {
                data: EdgeData::none(),
                sink_idx: *sink_idx,
            },
            InvolutiveMapping::Sink { source_idx, .. } => InvolutiveMapping::Sink {
                source_idx: *source_idx,
            },
        }
    }

    pub fn is_internal(&self) -> bool {
        !matches!(self, InvolutiveMapping::Identity { .. })
    }

    pub fn is_source(&self) -> bool {
        matches!(self, InvolutiveMapping::Source { .. })
    }

    pub fn is_sink(&self) -> bool {
        matches!(self, InvolutiveMapping::Sink { .. })
    }

    pub fn make_source(&mut self, sink_idx: Hedge) -> Option<InvolutiveMapping<E>> {
        let data = match self {
            InvolutiveMapping::Identity { data, .. } => data.take(),
            _ => EdgeData::none(),
        };
        Some(InvolutiveMapping::Source { data, sink_idx })
    }

    pub fn new_identity(data: E, orientation: impl Into<Orientation>, underlying: Flow) -> Self {
        let o = orientation.into();
        InvolutiveMapping::Identity {
            data: EdgeData::new(data, o),
            underlying,
        }
    }

    pub fn identity_dot(
        edge_id: Hedge,
        source: usize,
        attr: Option<&GVEdgeAttrs>,
        orientation: Orientation,
        flow: Flow,
    ) -> String {
        let mut out = "".to_string();
        out.push_str(&format!("ext{} [shape=none, label=\"\"];\n ", edge_id));

        out.push_str(&format!("ext{} -> {}[", edge_id, source));
        match (orientation, flow) {
            (Orientation::Default, Flow::Source) => {
                out.push_str("dir=back ");
            }
            (Orientation::Default, Flow::Sink) => {
                out.push_str("dir=forward ");
            }
            (Orientation::Reversed, Flow::Sink) => {
                out.push_str("dir=back ");
            }
            (Orientation::Reversed, Flow::Source) => {
                out.push_str("dir=forward ");
            }
            (Orientation::Undirected, _) => {
                out.push_str("dir=none ");
            }
        }
        if let Some(attr) = attr {
            out.push_str(&format!("{}", attr));
        }
        out.push_str("];\n");
        out
    }

    pub fn pair_dot(
        source: usize,
        sink: usize,
        attr: Option<&GVEdgeAttrs>,
        orientation: Orientation,
    ) -> String {
        let mut out = "".to_string();

        out.push_str(&format!("{} -> {}[", source, sink));
        match orientation {
            Orientation::Default => {
                out.push_str(" dir=forward ");
            }
            Orientation::Reversed => {
                out.push_str(" dir=back ");
            }
            Orientation::Undirected => {
                out.push_str(" dir=none ");
            }
        }
        if let Some(attr) = attr {
            out.push_str(&format!("{}];\n", attr));
        } else {
            out.push_str(" color=\"red:blue;0.5 \" ];\n");
        }
        out
    }

    pub fn default_dot(
        &self,
        edge_id: Hedge,
        source: Option<usize>,
        sink: Option<usize>,
        attr: Option<&GVEdgeAttrs>,
    ) -> String {
        let mut out = "".to_string();
        match self {
            InvolutiveMapping::Identity { data, underlying } => {
                out.push_str(&Self::identity_dot(
                    edge_id,
                    source.unwrap(),
                    attr,
                    data.orientation,
                    *underlying,
                ));
            }
            InvolutiveMapping::Source { data, .. } => {
                out.push_str(&Self::pair_dot(
                    source.unwrap(),
                    sink.unwrap(),
                    attr,
                    data.orientation,
                ));
            }
            InvolutiveMapping::Sink { .. } => {}
        }
        out
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Involution<N, E> {
    inv: Vec<InvolutiveMapping<E>>,
    hedge_data: Vec<N>,
    // pub inv: Vec<(N, InvolutiveMapping<E>)>,
}

impl<N, E> IntoIterator for Involution<N, E> {
    type Item = (InvolutiveMapping<E>, N);
    type IntoIter = std::iter::Zip<std::vec::IntoIter<InvolutiveMapping<E>>, std::vec::IntoIter<N>>;
    fn into_iter(self) -> Self::IntoIter {
        self.inv.into_iter().zip(self.hedge_data)
    }
}

impl<'a, N, E> IntoIterator for &'a Involution<N, E> {
    type Item = (&'a InvolutiveMapping<E>, &'a N);
    type IntoIter =
        std::iter::Zip<std::slice::Iter<'a, InvolutiveMapping<E>>, std::slice::Iter<'a, N>>;
    fn into_iter(self) -> Self::IntoIter {
        self.inv.iter().zip(self.hedge_data.iter())
    }
}

impl<'a, N, E> IntoIterator for &'a mut Involution<N, E> {
    type Item = (&'a mut InvolutiveMapping<E>, &'a mut N);
    type IntoIter =
        std::iter::Zip<std::slice::IterMut<'a, InvolutiveMapping<E>>, std::slice::IterMut<'a, N>>;
    fn into_iter(self) -> Self::IntoIter {
        self.inv.iter_mut().zip(self.hedge_data.iter_mut())
    }
}

impl<N, E> IndexMut<Hedge> for Involution<N, E> {
    fn index_mut(&mut self, index: Hedge) -> &mut Self::Output {
        &mut self.inv[index.0]
    }
}
impl<N, E> Index<Hedge> for Involution<N, E> {
    type Output = InvolutiveMapping<E>;
    fn index(&self, index: Hedge) -> &Self::Output {
        &self.inv[index.0]
    }
}

pub trait Get<H> {
    type Output<'a>
    where
        Self: 'a;

    type MutOutput<'a>
    where
        Self: 'a;
    fn get(&self, h: H) -> Self::Output<'_>;

    fn get_mut(&mut self, h: H) -> Self::MutOutput<'_>;
}

impl<N, E> Get<Hedge> for Involution<N, E> {
    type Output<'a>
        = &'a InvolutiveMapping<E>
    where
        Self: 'a;
    type MutOutput<'a>
        = &'a mut InvolutiveMapping<E>
    where
        Self: 'a;
    fn get(&self, h: Hedge) -> Self::Output<'_> {
        &self.inv[h.0]
    }

    fn get_mut(&mut self, h: Hedge) -> Self::MutOutput<'_> {
        &mut self.inv[h.0]
    }
}

impl<N, E> Get<EdgeId> for Involution<N, E> {
    type Output<'a>
        = Option<(Flow, &'a EdgeData<E>)>
    where
        Self: 'a;

    type MutOutput<'a>
        = Option<(Flow, &'a mut EdgeData<E>)>
    where
        Self: 'a;

    fn get(&self, h: EdgeId) -> Self::Output<'_> {
        let (source, flow) = match h {
            EdgeId::Paired { source, .. } => (source, Flow::Source),
            EdgeId::Unpaired { hedge, flow } => (hedge, flow),
            EdgeId::Split { source, split, .. } => (source, split),
        };

        match &self[source] {
            InvolutiveMapping::Source { data, .. } => Some((flow, data)),
            InvolutiveMapping::Identity { data, underlying } => Some((*underlying, data)),
            _ => None,
        }
    }

    fn get_mut(&mut self, h: EdgeId) -> Self::MutOutput<'_> {
        let (source, flow) = match h {
            EdgeId::Paired { source, .. } => (source, Flow::Source),
            EdgeId::Unpaired { hedge, flow } => (hedge, flow),
            EdgeId::Split { source, split, .. } => (source, split),
        };

        match &mut self[source] {
            InvolutiveMapping::Source { data, .. } => Some((flow, data)),
            InvolutiveMapping::Identity { data, underlying } => Some((*underlying, data)),
            _ => None,
        }
    }
}

impl<N, E> Display for Involution<N, E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut out = "".to_string();
        for (i, e) in self.inv.iter().enumerate() {
            match e {
                InvolutiveMapping::Identity { .. } => {
                    out.push_str(&format!("{}\n", i));
                }
                InvolutiveMapping::Source { sink_idx, .. } => {
                    out.push_str(&format!("{}->{}\n", i, sink_idx));
                }
                InvolutiveMapping::Sink { source_idx } => {
                    out.push_str(&format!("{}<-{}\n", i, source_idx));
                }
            }
        }
        write!(f, "{}", out)
    }
}

#[allow(dead_code)]
impl<N, E> Involution<N, E> {
    pub fn iter_idx(&self) -> impl Iterator<Item = Hedge> {
        (0..self.inv.len()).map(Hedge)
    }

    fn new() -> Self {
        Involution {
            inv: Vec::new(),
            hedge_data: Vec::new(),
        }
    }

    pub fn set_as_source(&mut self, hedge: Hedge) {
        let is_sink = self.inv[hedge.0].is_sink();
        if is_sink {
            self.flip_underlying(hedge)
        }
    }

    pub fn set_as_sink(&mut self, hedge: Hedge) {
        let is_source = self.inv[hedge.0].is_source();
        if is_source {
            self.flip_underlying(hedge)
        }
    }

    pub fn flip_orientation(&mut self, hedge: Hedge) {
        match &mut self[hedge] {
            InvolutiveMapping::Sink { source_idx } => {
                let source = *source_idx;
                self.flip_orientation(source);
            }
            InvolutiveMapping::Source { data, .. } => data.orientation = data.orientation.reverse(),
            InvolutiveMapping::Identity { data, .. } => {
                data.orientation = data.orientation.reverse()
            }
        }
    }

    pub fn flip_underlying(&mut self, hedge: Hedge) {
        let pair = self.inv(hedge);

        self.inv.swap(hedge.0, pair.0);

        match &mut self[hedge] {
            InvolutiveMapping::Identity { underlying, .. } => {
                *underlying = -*underlying;
            }
            InvolutiveMapping::Source { data, sink_idx } => {
                *sink_idx = pair;
                data.orientation = data.orientation.reverse();
            }
            InvolutiveMapping::Sink { source_idx } => {
                *source_idx = pair;
            }
        }

        match &mut self[pair] {
            InvolutiveMapping::Identity { underlying, .. } => {
                *underlying = -*underlying;
            }
            InvolutiveMapping::Source { data, sink_idx } => {
                *sink_idx = hedge;
                data.orientation = data.orientation.reverse();
            }
            InvolutiveMapping::Sink { source_idx } => {
                *source_idx = hedge;
            }
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = (Hedge, (&InvolutiveMapping<E>, &N))> {
        self.into_iter().enumerate().map(|(i, a)| (Hedge(i), a))
    }

    fn print<S: SubGraph>(
        &self,
        subgraph: &S,
        n_label: &impl Fn(&N) -> Option<String>,
        h_label: &impl Fn(&E) -> Option<String>,
    ) -> String {
        let mut out = "".to_string();
        for (i, (e, n)) in self.iter() {
            if !subgraph.includes(&i) {
                continue;
            }
            if let Some(l) = n_label(n) {
                out.push_str(&format!("{}:", l));
            }
            match e {
                InvolutiveMapping::Identity { .. } => {
                    out.push_str(&format!("{}\n", i));
                }
                InvolutiveMapping::Source { data, sink_idx } => {
                    if let Some(d) = &data.data {
                        if let Some(l) = h_label(d) {
                            out.push_str(&format!("{}-{}->{}\n", i, sink_idx, l));
                        } else {
                            out.push_str(&format!("{}->{}\n", i, sink_idx));
                        }
                    } else {
                        out.push_str(&format!("{}->{}\n", i, sink_idx));
                    }
                }
                InvolutiveMapping::Sink { source_idx } => {
                    out.push_str(&format!("{}<-{}\n", i, source_idx));
                }
            }
        }
        out
    }

    fn map_data_ref<F, G, N2, E2>(&self, f: F, g: &G) -> Involution<N2, E2>
    where
        F: FnMut(&N) -> N2,
        G: FnMut(&E) -> E2 + Clone,
    {
        let inv = self
            .inv
            .iter()
            .map(|e| (e.map_data_ref(g.clone())))
            .collect_vec();

        let hedge_data = self.hedge_data.iter().map(f).collect_vec();

        Involution { inv, hedge_data }
    }

    // fn map_edge_data_ref<'a, E2>(&'a self, g: impl FnMut(&'a E) -> E2) -> Involution<N, E2>
    // where
    //     N: Clone,
    // {
    //     let inv = self
    //         .inv
    //         .iter()
    //         .map(|(e)| (e.map_data_ref(&g)))
    //         .collect_vec();

    //     Involution {
    //         inv,
    //         hedge_data: self.hedge_data.clone(),
    //     }
    // }

    fn map_edge_data<E2>(
        self,
        mut g: impl FnMut(EdgeId, EdgeData<E>) -> EdgeData<E2>,
    ) -> Involution<N, E2>
    where
        N: Clone,
    {
        let inv = self
            .inv
            .into_iter()
            .enumerate()
            .map(|(i, e)| match e {
                InvolutiveMapping::Identity { data, underlying } => InvolutiveMapping::Identity {
                    data: g(
                        EdgeId::Unpaired {
                            hedge: Hedge(i),
                            flow: underlying,
                        },
                        data,
                    ),
                    underlying,
                },
                InvolutiveMapping::Source { data, sink_idx } => InvolutiveMapping::Source {
                    data: g(
                        EdgeId::Paired {
                            source: Hedge(i),
                            sink: sink_idx,
                        },
                        data,
                    ),

                    sink_idx,
                },
                InvolutiveMapping::Sink { source_idx } => InvolutiveMapping::Sink { source_idx },
            })
            .collect_vec();

        Involution {
            inv,
            hedge_data: self.hedge_data.clone(),
        }
    }

    fn forgetful_map_node_data_ref<N2, E2>(&self, f: impl FnMut(&N) -> N2) -> Involution<N2, E2> {
        let inv = self.inv.iter().map(|e| e.map_data_none()).collect_vec();

        let hedge_data = self.hedge_data.iter().map(f).collect_vec();

        Involution { inv, hedge_data }
    }

    fn n_internals<S: SubGraph>(&self, subgraph: &S) -> usize {
        subgraph
            .included_iter()
            .filter(|i| self.is_internal(*i, subgraph))
            .count()
            / 2
    }

    fn first_internal(&self) -> Option<usize> {
        self.inv.iter().position(|e| e.is_internal())
    }

    fn random(len: usize, seed: u64) -> Involution<Option<N>, ()> {
        let mut rng = SmallRng::seed_from_u64(seed);

        let mut inv = Involution::new();

        for _ in 0..len {
            let r = rng.gen_bool(0.1);
            if r {
                inv.add_identity((), None, Orientation::Undirected, Flow::Sink);
            } else {
                inv.add_pair((), None, None, false);
            }
        }

        inv
    }

    fn add_identity(
        &mut self,
        data: E,
        id: N,
        orientation: impl Into<Orientation>,
        underlying: Flow,
    ) -> Hedge {
        let index = self.inv.len();
        self.inv.push(InvolutiveMapping::new_identity(
            data,
            orientation,
            underlying,
        ));
        self.hedge_data.push(id);
        Hedge(index)
    }

    pub fn iter_edge_data(&self) -> impl Iterator<Item = &EdgeData<E>> {
        self.inv.iter().filter_map(|e| match e {
            InvolutiveMapping::Source { data, .. } => Some(data),
            InvolutiveMapping::Identity { data, .. } => Some(data),
            _ => None,
        })
    }

    fn is_internal<S: SubGraph>(&self, index: Hedge, subgraph: &S) -> bool {
        if !subgraph.includes(&index) {
            return false;
        }
        match &self.inv[index.0] {
            InvolutiveMapping::Identity { .. } => false,
            InvolutiveMapping::Source { sink_idx, .. } => subgraph.includes(sink_idx),
            InvolutiveMapping::Sink { source_idx } => subgraph.includes(source_idx),
        }
    }

    pub fn find_from_data(&self, data: &E) -> Option<usize>
    where
        E: PartialEq,
    {
        (0..self.inv.len()).find(|&i| self.edge_data(Hedge(i)).as_ref().data == Some(data))
    }

    fn connect_to_identity(&mut self, source_idx: Hedge, hedge_data: N) -> Hedge {
        let sink = self.inv.len();
        self.inv.push(InvolutiveMapping::Sink { source_idx });
        self.hedge_data.push(hedge_data);
        self.inv[source_idx.0] = self.inv[source_idx.0]
            .make_source(Hedge(sink))
            .unwrap_or_else(|| panic!("Source must be an identity mapping"));
        Hedge(sink)
    }

    fn connect_identities(
        &mut self,
        source_idx: Hedge,
        sink_idx: Hedge,
        merge_fn: &impl Fn(Flow, EdgeData<E>, Flow, EdgeData<E>) -> (Flow, EdgeData<E>),
    ) {
        let (source, sink) = if source_idx.0 < sink_idx.0 {
            let (a, b) = self.inv.split_at_mut(sink_idx.0);

            (&mut a[source_idx.0], &mut b[0])
        } else {
            let (a, b) = self.inv.split_at_mut(source_idx.0);

            (&mut b[0], &mut a[sink_idx.0])
        };
        let (flow, new_data) = if let (
            InvolutiveMapping::Identity { data, underlying },
            InvolutiveMapping::Identity {
                data: sink_data,
                underlying: sink_underlying,
            },
        ) = (source, sink)
        {
            merge_fn(*underlying, data.take(), *sink_underlying, sink_data.take())
        } else {
            return;
        };

        let (source, sink) = if source_idx.0 < sink_idx.0 {
            let (a, b) = self.inv.split_at_mut(sink_idx.0);

            (&mut a[source_idx.0], &mut b[0])
        } else {
            let (a, b) = self.inv.split_at_mut(source_idx.0);

            (&mut b[0], &mut a[sink_idx.0])
        };

        match flow {
            Flow::Sink => {
                *source = InvolutiveMapping::Sink {
                    source_idx: sink_idx,
                };
                *sink = InvolutiveMapping::Source {
                    data: new_data,
                    sink_idx: source_idx,
                };
            }
            Flow::Source => {
                *source = InvolutiveMapping::Source {
                    data: new_data,
                    sink_idx,
                };
                *sink = InvolutiveMapping::Sink { source_idx };
            }
        }
    }

    fn connect(&mut self, source_idx: Hedge, sink_idx: Hedge)
    where
        E: Clone,
    {
        self.inv[source_idx.0] = self.inv[source_idx.0].make_source(sink_idx).unwrap();

        if let InvolutiveMapping::Identity { .. } = &self.inv[sink_idx.0] {
            self.inv[sink_idx.0] = InvolutiveMapping::Sink { source_idx };
        } else {
            panic!("Sink must be an identity mapping")
        }
    }

    fn add_pair(
        &mut self,
        data: E,
        source_id: N,
        sink_id: N,
        directed: impl Into<Orientation>,
    ) -> (Hedge, Hedge) {
        let orientation = directed.into();
        let source = self.add_identity(data, source_id, orientation, Flow::Sink);
        let sink = self.connect_to_identity(source, sink_id);
        (source, sink)
    }

    fn len(&self) -> usize {
        self.inv.len()
    }

    pub fn inv(&self, index: Hedge) -> Hedge {
        match &self.inv[index.0] {
            InvolutiveMapping::Source { sink_idx, .. } => *sink_idx,
            InvolutiveMapping::Sink { source_idx } => *source_idx,
            _ => index,
        }
    }

    pub fn superficial_orientation(&self, index: Hedge) -> Orientation {
        self.edge_data(index).orientation
    }

    pub fn edge_data(&self, index: Hedge) -> &EdgeData<E> {
        match &self.inv[index.0] {
            InvolutiveMapping::Source { data, .. } => data,
            InvolutiveMapping::Identity { data, .. } => data,
            InvolutiveMapping::Sink { source_idx } => self.edge_data(*source_idx),
        }
    }

    #[allow(clippy::needless_return)]
    pub fn edge_data_mut(&mut self, index: Hedge) -> &mut EdgeData<E> {
        if let InvolutiveMapping::Sink { source_idx } = self.inv[index.0] {
            return self.edge_data_mut(source_idx);
        }
        match &mut self.inv[index.0] {
            InvolutiveMapping::Source { data, .. } => return data,
            InvolutiveMapping::Identity { data, .. } => return data,
            _ => unreachable!(),
        };
    }

    fn smart_data<S: SubGraph>(&self, index: Hedge, subgraph: &S) -> Option<&EdgeData<E>> {
        if subgraph.includes(&index) {
            match &self[index] {
                InvolutiveMapping::Identity { .. } => Some(self.edge_data(index)),
                InvolutiveMapping::Source { .. } => Some(self.edge_data(index)),
                InvolutiveMapping::Sink { source_idx } => {
                    if subgraph.includes(source_idx) {
                        None
                    } else {
                        Some(self.edge_data(*source_idx))
                    }
                }
            }
        } else {
            None
        }
    }

    #[allow(clippy::needless_return)]
    fn smart_data_mut<S: SubGraph>(
        &mut self,
        index: Hedge,
        subgraph: &S,
    ) -> Option<&mut EdgeData<E>> {
        if subgraph.includes(&index) {
            if let InvolutiveMapping::Sink { source_idx } = &self.inv[index.0] {
                if subgraph.includes(source_idx) {
                    return None;
                } else {
                    return Some(self.edge_data_mut(*source_idx));
                }
            }
            match &mut self.inv[index.0] {
                InvolutiveMapping::Source { data, .. } => return Some(data),
                InvolutiveMapping::Identity { data, .. } => return Some(data),
                _ => unreachable!(),
            };
        } else {
            None
        }
    }

    // pub fn edge_data_mut(&mut self, index: InvolutionIdx) -> &mut EdgeData<E> {
    //     let source = match &mut self.inv[index.0] {
    //         InvolutiveMapping::Source { data, .. } => return data,
    //         InvolutiveMapping::Identity { data, .. } => return data,
    //         InvolutiveMapping::Sink { source_idx } => source_idx.0,
    //     };

    //     self.edge_data_mut(InvolutionIdx(source))
    // }

    fn set_edge_data(&mut self, index: Hedge, new_data: E) -> EdgeData<E> {
        match &mut self[index] {
            InvolutiveMapping::Source { data, .. } => {
                let old = data.take();
                *data = EdgeData::new(new_data, old.orientation);
                old
            }
            InvolutiveMapping::Identity { data, .. } => {
                let old = data.take();
                *data = EdgeData::new(new_data, old.orientation);
                old
            }
            InvolutiveMapping::Sink { source_idx } => {
                let source = *source_idx;
                self.set_edge_data(source, new_data)
            }
        }
    }

    pub fn hedge_data(&self, index: Hedge) -> &N {
        &self.hedge_data[index.0]
    }

    fn set_hedge_data(&mut self, index: Hedge, id: N) {
        self.hedge_data[index.0] = id;
    }

    pub fn involved_hedge_data(&self, index: Hedge) -> Option<&N> {
        match &self.inv[index.0] {
            InvolutiveMapping::Source { sink_idx, .. } => Some(self.hedge_data(*sink_idx)),
            InvolutiveMapping::Sink { source_idx } => Some(self.hedge_data(*source_idx)),
            InvolutiveMapping::Identity { .. } => None,
        }
    }

    pub fn is_identity(&self, index: Hedge) -> bool {
        self.inv[index.0].is_identity()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
pub struct GVEdgeAttrs {
    pub label: Option<String>,
    pub color: Option<String>,
    pub other: Option<String>,
}

impl std::fmt::Display for GVEdgeAttrs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let out = [
            ("label=", self.label.as_ref()),
            ("color=", self.color.as_ref()),
            ("", self.other.as_ref()),
        ]
        .iter()
        .filter_map(|(prefix, x)| x.map(|s| format!("{}{}", prefix, s)))
        .join(",")
        .to_string();
        write!(f, "{}", out)
    }
}
pub mod subgraph;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HedgeGraph<E, V> {
    node_data: Vec<V>,
    nodes: Vec<HedgeNode>,
    base_nodes: usize,
    pub involution: Involution<NodeIndex, E>, // Involution of half-edges
}

impl<E, V> Index<&HedgeNode> for HedgeGraph<E, V> {
    type Output = V;
    fn index(&self, index: &HedgeNode) -> &Self::Output {
        let id = self.id_from_hairs(index).unwrap();
        &self[id]
    }
}

impl<E, V> IndexMut<&HedgeNode> for HedgeGraph<E, V> {
    fn index_mut(&mut self, index: &HedgeNode) -> &mut Self::Output {
        let id = self.id_from_hairs(index).unwrap();
        &mut self[id]
    }
}

impl<E, V> Index<NodeIndex> for HedgeGraph<E, V> {
    type Output = V;
    fn index(&self, index: NodeIndex) -> &Self::Output {
        &self.node_data[index.0]
    }
}

impl<E, V> IndexMut<NodeIndex> for HedgeGraph<E, V> {
    fn index_mut(&mut self, index: NodeIndex) -> &mut Self::Output {
        &mut self.node_data[index.0]
    }
}

#[derive(Clone, Debug)]
pub struct HedgeNodeBuilder<V> {
    data: V,
    hedges: Vec<Hedge>,
}

#[derive(Clone, Debug)]
pub struct HedgeGraphBuilder<E, V> {
    nodes: Vec<HedgeNodeBuilder<V>>,
    involution: Involution<NodeIndex, E>,
}

pub struct NodeIterator<'a, E, V, I = IterOnes<'a, usize, Lsb0>> {
    graph: &'a HedgeGraph<E, V>,
    edges: I,
    seen: BitVec,
}

impl<'a, E, V, I: Iterator<Item = Hedge>> Iterator for NodeIterator<'a, E, V, I> {
    type Item = (&'a HedgeNode, &'a V);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(next) = self.edges.next() {
            let node = self.graph.node_hairs(next);
            let node_pos = self.graph.id_from_hairs(node).unwrap().0;

            if self.seen[node_pos] {
                self.next()
            } else {
                self.seen.set(node_pos, true);
                Some((node, &self.graph[node]))
            }
        } else {
            None
        }
    }
}

impl HedgeGraph<usize, usize> {
    pub fn close_externals(value: &BareGraph) -> Self {
        let mut builder = HedgeGraphBuilder::new();
        let mut map = AHashMap::new();

        let mut ext_map = AHashMap::new();

        for (i, v) in value.vertices.iter().enumerate() {
            if matches!(v.vertex_info, VertexInfo::InteractonVertexInfo(_)) {
                map.insert(i, builder.add_node(i));
            }
        }

        for (i, edge) in value.edges.iter().enumerate() {
            match edge.edge_type {
                EdgeType::Incoming => {
                    let key = edge.particle.pdg_code;
                    let sink = map[&edge.vertices[1]];
                    ext_map
                        .entry(key)
                        .or_insert_with(|| (i, [None, Some(sink)]))
                        .1[1] = Some(sink);
                }
                EdgeType::Outgoing => {
                    let key = if edge.particle.is_self_antiparticle() {
                        edge.particle.pdg_code
                    } else {
                        -edge.particle.pdg_code
                    };
                    let source = map[&edge.vertices[0]];
                    ext_map
                        .entry(key)
                        .or_insert_with(|| (i, [Some(source), None]))
                        .1[0] = Some(source);
                }
                EdgeType::Virtual => {
                    let sink = map[&edge.vertices[1]];
                    let source = map[&edge.vertices[0]];
                    builder.add_edge(source, sink, i, true);
                }
            }
        }

        for (_, (i, [source, sink])) in ext_map {
            println!("{source:?},{sink:?}");
            if let (Some(source), Some(sink)) = (source, sink) {
                builder.add_edge(source, sink, i, true);
            }
        }

        builder.into()
    }
}

impl<E, V> HedgeGraphBuilder<E, V> {
    pub fn new() -> Self {
        HedgeGraphBuilder {
            nodes: Vec::new(),
            involution: Involution::new(),
        }
    }

    pub fn build(self) -> HedgeGraph<E, V> {
        self.into()
    }

    pub fn add_node(&mut self, data: V) -> NodeIndex {
        let index = self.nodes.len();
        self.nodes.push(HedgeNodeBuilder {
            data,
            hedges: Vec::new(),
        });
        NodeIndex(index)
    }

    pub fn add_edge(
        &mut self,
        source: NodeIndex,
        sink: NodeIndex,
        data: E,
        directed: impl Into<Orientation>,
    ) {
        let id = self.involution.add_pair(data, source, sink, directed);
        self.nodes[source.0].hedges.push(id.0);

        self.nodes[sink.0].hedges.push(id.1);
    }

    pub fn add_external_edge(
        &mut self,
        source: NodeIndex,
        data: E,
        orientation: impl Into<Orientation>,
        underlying: Flow,
    ) {
        let id = self
            .involution
            .add_identity(data, source, orientation, underlying);
        self.nodes[source.0].hedges.push(id);
    }
}

impl<E, V> Default for HedgeGraphBuilder<E, V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<E, V> From<HedgeGraphBuilder<E, V>> for HedgeGraph<E, V> {
    fn from(builder: HedgeGraphBuilder<E, V>) -> Self {
        let len = builder.involution.len();
        let nodes: Vec<HedgeNode> = builder
            .nodes
            .iter()
            .map(|x| (HedgeNode::from_builder(x, len)))
            .collect();

        let node_data: Vec<V> = builder.nodes.into_iter().map(|x| x.data).collect();
        let involution = Involution {
            inv: builder.involution.inv,
            hedge_data: builder.involution.hedge_data,
        };
        HedgeGraph {
            base_nodes: nodes.len(),
            nodes,
            node_data,
            involution,
        }
    }
}

impl<N: Clone, E: Clone> From<symbolica::graph::Graph<N, E>> for HedgeGraph<E, N> {
    fn from(graph: symbolica::graph::Graph<N, E>) -> Self {
        let mut builder = HedgeGraphBuilder::new();
        let mut map = AHashMap::new();

        for (i, node) in graph.nodes().iter().enumerate() {
            map.insert(i, builder.add_node(node.data.clone()));
        }

        for edge in graph.edges() {
            let vertices = edge.vertices;
            let source = map[&vertices.0];
            let sink = map[&vertices.1];
            builder.add_edge(source, sink, edge.data.clone(), edge.directed);
        }

        builder.into()
    }
}

impl From<&BareGraph> for HedgeGraph<usize, usize> {
    fn from(value: &BareGraph) -> Self {
        let mut builder = HedgeGraphBuilder::new();
        let mut map = AHashMap::new();

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

use subgraph::{
    Cycle, HedgeNode, Inclusion, InternalSubGraph, SubGraph, SubGraphHedgeIter, SubGraphOps,
};

use thiserror::Error;

use crate::momentum::SignOrZero;

use super::{BareGraph, EdgeType, VertexInfo};

#[derive(Error, Debug)]

pub enum HedgeError {
    #[error("Invalid start node")]
    InvalidStart,
}

pub struct EdgeIdIter<'a, E, V, S> {
    graph: &'a HedgeGraph<E, V>,
    included_iter: SubGraphHedgeIter<'a>,
    subgraph: &'a S,
}
impl<'a, E, V, S> EdgeIdIter<'a, E, V, S>
where
    S: SubGraph,
{
    pub fn new(graph: &'a HedgeGraph<E, V>, subgraph: &'a S) -> Self {
        EdgeIdIter {
            graph,
            subgraph,
            included_iter: subgraph.included_iter(),
        }
    }
}

impl<E, V, S> Iterator for EdgeIdIter<'_, E, V, S>
where
    S: SubGraph,
{
    type Item = EdgeId;

    fn next(&mut self) -> Option<Self::Item> {
        let i = self.included_iter.next()?;
        if let Some(e) = EdgeId::from_source_with_subgraph(i, &self.graph.involution, self.subgraph)
        {
            Some(e)
        } else {
            self.next()
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, Error)]
pub enum HedgeGraphError {
    #[error("Nodes do not partition")]
    NodesDoNotPartition,
    #[error("Invalid node")]
    NoNode,
    #[error("Invalid edge")]
    NoEdge,
}

impl<E, V> HedgeGraph<E, V> {
    pub fn iter_nodes(&self) -> impl Iterator<Item = (&HedgeNode, &V)> {
        self.nodes.iter().zip(self.node_data.iter())
    }

    pub fn base_nodes_iter(&self) -> impl Iterator<Item = NodeIndex> {
        (0..self.base_nodes).map(NodeIndex)
    }

    /// Splits the edge that hedge is a part of into two dangling hedges, adding the data to the side given by hedge. The underlying orientation of the new edges is the same as the original edge, i.e. the source will now have `Flow::Source` and the sink will have `Flow::Sink`. The superficial orientation has to be given knowing this.
    pub fn split_edge(&mut self, hedge: Hedge, data: EdgeData<E>) -> Result<(), HedgeGraphError> {
        let (data_hedge, other) = match &mut self.involution[hedge] {
            InvolutiveMapping::Source { data, sink_idx } => {
                let old_data = data.take();
                (Some(old_data), *sink_idx)
            }
            InvolutiveMapping::Sink { source_idx } => (None, *source_idx),
            _ => return Err(HedgeGraphError::NoEdge),
        };

        if let Some(d) = data_hedge {
            self.involution.inv[hedge.0] = InvolutiveMapping::Identity {
                data,
                underlying: Flow::Source,
            };
            self.involution.inv[other.0] = InvolutiveMapping::Identity {
                data: d,
                underlying: Flow::Sink,
            };
        } else {
            self.involution.inv[hedge.0] = InvolutiveMapping::Identity {
                data,
                underlying: Flow::Sink,
            };

            let data = self.involution.edge_data_mut(other).take();

            self.involution.inv[other.0] = InvolutiveMapping::Identity {
                data,
                underlying: Flow::Source,
            };
        }

        Ok(())
    }

    /// Joins two graphs together, matching edges with the given function and merging them with the given function.
    /// The function `matching_fn` should return true if the two dangling half edges should be matched.
    /// The function `merge_fn` should return the new data for the merged edge, given the data of the two edges being merged.
    ///
    pub fn join(
        self,
        other: Self,
        matching_fn: impl Fn(Flow, EdgeData<&E>, Flow, EdgeData<&E>) -> bool,
        merge_fn: impl Fn(Flow, EdgeData<E>, Flow, EdgeData<E>) -> (Flow, EdgeData<E>),
    ) -> Result<Self, HedgeGraphError> {
        let self_inv_len = self.empty_filter();
        let mut full_self = self.full_filter();
        let other_inv_len = other.empty_filter();
        let mut full_other = self_inv_len.clone();
        full_self.extend(other_inv_len.clone());
        full_other.extend(other.full_filter());

        let mut node_data = self.node_data;
        node_data.extend(other.node_data);

        let nodes: Vec<_> = self
            .nodes
            .into_iter()
            .map(|mut k| {
                k.hairs.extend(other_inv_len.clone());
                k.internal_graph.filter.extend(other_inv_len.clone());
                k
            })
            .chain(other.nodes.into_iter().map(|mut k| {
                let mut new_hairs = self_inv_len.clone();
                new_hairs.extend(k.hairs.clone());
                k.hairs = new_hairs;

                let mut internal = self_inv_len.clone();
                internal.extend(k.internal_graph.filter.clone());
                k.internal_graph.filter = internal;

                k
            }))
            .collect();

        let involution = Involution {
            inv: self
                .involution
                .inv
                .into_iter()
                .chain(other.involution.inv)
                .collect(),
            hedge_data: self
                .involution
                .hedge_data
                .into_iter()
                .chain(other.involution.hedge_data)
                .collect(),
        };

        let mut found_match = true;

        let mut g = HedgeGraph {
            base_nodes: self.base_nodes + other.base_nodes, // need to fix
            nodes,
            node_data,
            involution,
        };

        while found_match {
            let mut matching_ids = None;

            for i in full_self.included_iter() {
                if let InvolutiveMapping::Identity {
                    data: datas,
                    underlying: underlyings,
                } = &g.involution.inv[i.0]
                {
                    for j in full_other.included_iter() {
                        if let InvolutiveMapping::Identity { data, underlying } =
                            &g.involution.inv[j.0]
                        {
                            if matching_fn(*underlyings, datas.as_ref(), *underlying, data.as_ref())
                            {
                                matching_ids = Some((i, j));
                                break;
                            }
                        }
                    }
                }
            }

            if let Some((source, sink)) = matching_ids {
                g.involution.connect_identities(source, sink, &merge_fn);
            } else {
                found_match = false;
            }
        }

        g.check_and_set_nodes()?;

        Ok(g)
    }

    fn check_and_set_nodes(&mut self) -> Result<(), HedgeGraphError> {
        let mut cover = self.empty_filter();
        for i in 0..self.base_nodes {
            let node = self.nodes.get(i).unwrap();
            for h in node.hairs.included_iter() {
                if cover.includes(&h) {
                    return Err(HedgeGraphError::NodesDoNotPartition);
                } else {
                    cover.set(h.0, true);
                    self.involution.set_hedge_data(h, NodeIndex(i));
                }
            }
        }

        if cover.sym_diff(&self.full_filter()).count_ones() > 0 {
            return Err(HedgeGraphError::NodesDoNotPartition);
        }

        Ok(())
    }

    /// Adds a dangling edge to the specified node with the given data and superficial orientation. /// The underlying orientation is always `Flow::Source`.
    pub fn add_dangling_edge(
        self,
        source: NodeIndex,
        data: E,
        orientation: impl Into<Orientation>,
    ) -> Result<(Hedge, Self), HedgeGraphError> {
        let mut involution = self.involution;
        let o = orientation.into();
        let mut data = Some(data);
        let mut hedge = None;

        let nodes: Vec<_> = self
            .nodes
            .into_iter()
            .enumerate()
            .map(|(i, mut k)| {
                k.internal_graph.filter.push(false);
                if NodeIndex(i) == source {
                    if let Some(d) = data.take() {
                        k.hairs.push(true);
                        hedge = Some(involution.add_identity(d, source, o, Flow::Source));
                    } else {
                        k.hairs.push(false);
                    }
                } else {
                    k.hairs.push(false);
                }

                k
            })
            .collect();

        let mut g = HedgeGraph {
            base_nodes: self.base_nodes,
            node_data: self.node_data,
            nodes,
            involution,
        };

        g.check_and_set_nodes()?;

        if let Some(hedge) = hedge {
            Ok((hedge, g))
        } else {
            Err(HedgeGraphError::NoNode)
        }
    }

    pub fn iter_edge_id<'a, S: SubGraph>(&'a self, subgraph: &'a S) -> EdgeIdIter<'a, E, V, S> {
        EdgeIdIter::new(self, subgraph)
    }

    // pub fn map_edges_ref<E2>(&self, f: &impl Fn(&E) -> E2) -> HedgeGraph<E2, V>
    // where
    //     V: Clone,
    // {
    //     let involution = self.involution.map_edge_data_ref(f);
    //     let nodes = self.nodes.clone();
    //     HedgeGraph {
    //         base_nodes: self.base_nodes,
    //         nodes,
    //         involution,
    //     }
    // }

    pub fn map_nodes_ref<V2>(&self, f: &impl Fn(&Self, &V, &HedgeNode) -> V2) -> HedgeGraph<E, V2>
    where
        E: Clone,
    {
        let involution = self.involution.clone();
        let node_data = self
            .node_data
            .iter()
            .zip(self.nodes.iter())
            .map(|(v, h)| f(self, v, h))
            .collect();
        HedgeGraph {
            base_nodes: self.base_nodes,
            node_data,
            nodes: self.nodes.clone(),
            involution,
        }
    }

    pub fn map<E2, V2>(
        self,
        f: impl FnMut(V) -> V2,
        g: impl FnMut(EdgeId, EdgeData<E>) -> EdgeData<E2>,
    ) -> HedgeGraph<E2, V2> {
        let involution = self.involution.map_edge_data(g);
        let node_data = self.node_data.into_iter().map(f).collect();
        HedgeGraph {
            base_nodes: self.base_nodes,
            nodes: self.nodes,
            node_data,
            involution,
        }
    }

    pub fn iter_egde_data<'a, S: SubGraph>(
        &'a self,
        subgraph: &'a S,
    ) -> impl Iterator<Item = EdgeData<&'a E>> + 'a {
        subgraph
            .included_iter()
            .map(|i| self.involution.edge_data(i).as_ref())
    }

    pub fn iter_egdes<'a, S: SubGraph>(
        &'a self,
        subgraph: &'a S,
    ) -> impl Iterator<Item = (EdgeId, EdgeData<&'a E>)> + 'a {
        subgraph.included_iter().flat_map(|i| {
            self.involution.smart_data(i, subgraph).map(|d| {
                (
                    EdgeId::from_half_edge_with_subgraph(i, &self.involution, subgraph).unwrap(),
                    d.as_ref(),
                )
            })
        })
    }

    pub fn iter_internal_edge_data<'a, S: SubGraph>(
        &'a self,
        subgraph: &'a S,
    ) -> impl Iterator<Item = EdgeData<&'a E>> + 'a {
        subgraph
            .included_iter()
            .flat_map(|i| self.involution.smart_data(i, subgraph).map(|d| d.as_ref()))
    }

    pub fn is_connected<S: SubGraph>(&self, subgraph: &S) -> bool {
        let n_edges = subgraph.nedges(self);
        if let Some(start) = subgraph.included_iter().next() {
            TraversalTree::dfs(self, subgraph, self.node_hairs(start), None)
                .covers()
                .nedges(self)
                == n_edges
        } else {
            true
        }
    }

    pub fn count_connected_components<S: SubGraph>(&self, subgraph: &S) -> usize {
        self.connected_components(subgraph).len()
    }

    pub fn node_hairs(&self, hedge: Hedge) -> &HedgeNode {
        self.hairs_from_id(*self.involution.hedge_data(hedge))
    }

    pub fn hairs_from_id(&self, id: NodeIndex) -> &HedgeNode {
        self.nodes.get(id.0).unwrap()
    }

    pub fn id_from_hairs(&self, id: &HedgeNode) -> Option<NodeIndex> {
        let e = id.hairs.included_iter().next()?;
        let nodeid = self.involution.hedge_data(e);

        Some(*nodeid)
    }

    pub fn involved_node_hairs(&self, hedge: Hedge) -> Option<&HedgeNode> {
        self.involution
            .involved_hedge_data(hedge)
            .map(|i| self.hairs_from_id(*i))
    }

    pub fn involved_node_id(&self, hedge: Hedge) -> Option<NodeIndex> {
        self.involution.involved_hedge_data(hedge).copied()
    }

    pub fn connected_components<S: SubGraph>(&self, subgraph: &S) -> Vec<BitVec> {
        let mut visited_edges = self.empty_filter();

        let mut components = vec![];

        // Iterate over all edges in the subgraph
        for hedge_index in subgraph.included_iter() {
            if !visited_edges.includes(&hedge_index) {
                // Perform DFS to find all reachable edges from this edge

                //
                let root_node = self.node_hairs(hedge_index);
                let reachable_edges = TraversalTree::dfs(self, subgraph, root_node, None).covers();

                visited_edges.union_with(&reachable_edges);

                components.push(reachable_edges);
            }
        }
        components
    }

    pub fn cyclotomatic_number(&self, subgraph: &InternalSubGraph) -> usize {
        let n_hedges = self.count_internal_edges(subgraph);
        // println!("n_hedges: {}", n_hedges);
        let n_nodes = self.number_of_nodes_in_subgraph(subgraph);
        // println!("n_nodes: {}", n_nodes);
        let n_components = self.count_connected_components(subgraph);
        // println!("n_components: {}", n_components);

        n_hedges - n_nodes + n_components
    }

    pub fn count_internal_edges(&self, subgraph: &InternalSubGraph) -> usize {
        let mut internal_edge_count = 0;
        // Iterate over all half-edges in the subgraph
        for hedge_index in subgraph.included_iter() {
            let inv_hedge_index = self.involution.inv(hedge_index);

            // Check if the involuted half-edge is also in the subgraph
            if subgraph.includes(&inv_hedge_index) {
                // To avoid double-counting, only count when hedge_index < inv_hedge_index
                if hedge_index < inv_hedge_index {
                    internal_edge_count += 1;
                }
            }
        }
        internal_edge_count
    }

    fn non_cut_edges_impl(
        &self,
        connected_components: usize,
        cyclotomatic_number: usize,
        start: Hedge,
        current: &mut BitVec,
        set: &mut AHashSet<BitVec>,
    ) {
        if current.count_ones() > 2 * cyclotomatic_number {
            return;
        }

        let complement = current.complement(self);

        if current.count_ones() > 0
            && self.count_connected_components(&complement) == connected_components
            && complement.covers(self) == self.full_filter()
        {
            // println!("//inserted with {con_comp}");
            set.insert(current.clone());
        }

        for i in (start.0..self.involution.len()).map(Hedge) {
            let j = self.involution.inv(i);
            if i > j {
                current.set(i.0, true);
                current.set(j.0, true);
                self.non_cut_edges_impl(
                    connected_components,
                    cyclotomatic_number,
                    Hedge(i.0 + 1),
                    current,
                    set,
                );
                current.set(i.0, false);
                current.set(j.0, false);
            }
        }
    }

    pub fn non_cut_edges(&self) -> AHashSet<BitVec> {
        let connected_components = self.count_connected_components(&self.full_filter());

        let cyclotomatic_number = self.cyclotomatic_number(&self.full_node().internal_graph);

        let mut current = self.empty_filter();
        let mut set = AHashSet::new();

        self.non_cut_edges_impl(
            connected_components,
            cyclotomatic_number,
            Hedge(0),
            &mut current,
            &mut set,
        );

        set
    }

    /// including pos
    pub fn neighbors<S: SubGraph>(&self, subgraph: &S, pos: Hedge) -> BitVec {
        subgraph.hairs(self.node_hairs(pos))
    }

    pub fn connected_neighbors<S: SubGraph>(&self, subgraph: &S, pos: Hedge) -> Option<BitVec> {
        Some(subgraph.hairs(self.involved_node_hairs(pos)?))
    }

    pub fn iter_egde_node<'a, S: SubGraph>(
        &'a self,
        subgraph: &'a S,
    ) -> impl Iterator<Item = &'a HedgeNode> + 'a {
        subgraph.included_iter().map(|i| self.node_hairs(i))
    }

    pub fn iter_node_data<'a, S: SubGraph>(
        &'a self,
        subgraph: &'a S,
    ) -> impl Iterator<Item = (&'a HedgeNode, &'a V)> {
        NodeIterator {
            graph: self,
            edges: subgraph.included_iter(),
            seen: bitvec![usize, Lsb0; 0; self.base_nodes],
        }
    }

    pub fn n_hedges(&self) -> usize {
        self.involution.len()
    }

    pub fn n_nodes(&self) -> usize {
        self.base_nodes
    }

    pub fn n_externals(&self) -> usize {
        self.involution
            .inv
            .iter()
            .filter(|e| e.is_identity())
            .count()
    }

    pub fn n_internals(&self) -> usize {
        self.involution
            .inv
            .iter()
            .filter(|e| e.is_internal())
            .count()
            / 2
    }

    pub fn n_base_nodes(&self) -> usize {
        self.nodes.iter().filter(|n| n.is_node()).count()
    }

    pub fn superficial_hedge_orientation(&self, hedge: Hedge) -> Option<Flow> {
        match &self.involution[hedge] {
            InvolutiveMapping::Identity { data, underlying } => {
                data.orientation.relative_to(*underlying).try_into().ok()
            }
            InvolutiveMapping::Sink { source_idx } => self
                .superficial_hedge_orientation(*source_idx)
                .map(Neg::neg),
            InvolutiveMapping::Source { data, .. } => data.orientation.try_into().ok(),
        }
    }

    pub fn underlying_hedge_orientation(&self, hedge: Hedge) -> Flow {
        match &self.involution[hedge] {
            InvolutiveMapping::Identity { underlying, .. } => *underlying,
            InvolutiveMapping::Sink { .. } => Flow::Sink,
            InvolutiveMapping::Source { .. } => Flow::Source,
        }
    }

    pub fn random(nodes: usize, edges: usize, seed: u64) -> HedgeGraph<(), ()> {
        let mut inv: Involution<Option<NodeIndex>, ()> =
            Involution::<NodeIndex, ()>::random(edges, seed);

        let mut rng = SmallRng::seed_from_u64(seed);

        let mut externals = Vec::new();
        let mut sources = Vec::new();
        let mut sinks = Vec::new();

        for (i, e) in inv.inv.iter().enumerate() {
            let nodeid = HedgeNode::node_from_pos(&[i], inv.inv.len());
            match e {
                InvolutiveMapping::Identity { .. } => externals.push(nodeid),
                InvolutiveMapping::Source { .. } => sources.push(nodeid),
                InvolutiveMapping::Sink { .. } => sinks.push(nodeid),
            }
        }

        while !externals.is_empty() {
            if rng.gen_bool(0.5) {
                let source_i = rng.gen_range(0..sources.len());

                sources[source_i].union_with(&externals.pop().unwrap());
            } else {
                let sink_i = rng.gen_range(0..sinks.len());

                sinks[sink_i].union_with(&externals.pop().unwrap());
            }
        }

        let mut lengthone = false;

        while sources.len() + sinks.len() > nodes {
            if rng.gen_bool(0.5) {
                if sources.len() <= 1 {
                    if lengthone {
                        break;
                    }
                    lengthone = true;
                    continue;
                }

                let idx1 = rng.gen_range(0..sources.len());
                let idx2 = rng.gen_range(0..sources.len() - 1);

                let n_i = sources.swap_remove(idx1);
                sources[idx2].union_with(&n_i);
            } else {
                if sinks.len() <= 1 {
                    if lengthone {
                        break;
                    }
                    lengthone = true;
                    continue;
                }
                let idx1 = rng.gen_range(0..sinks.len());
                let idx2 = rng.gen_range(0..sinks.len() - 1);
                let n_i = sinks.swap_remove(idx1);
                sinks[idx2].union_with(&n_i);
            }
        }

        let mut nodes = Vec::new();
        let mut node_data = Vec::new();

        for (nid, n) in sources.into_iter().enumerate() {
            nodes.push(n.clone());
            node_data.push(());
            for i in n.hairs.included_iter() {
                inv.set_hedge_data(i, Some(NodeIndex(nid)));
            }
        }

        let len = nodes.len();

        for (nid, n) in sinks.into_iter().enumerate() {
            nodes.push(n.clone());
            node_data.push(());

            for i in n.hairs.included_iter() {
                inv.set_hedge_data(i, Some(NodeIndex(nid + len)));
            }
        }

        let new_inv = Involution {
            inv: inv.inv,
            hedge_data: inv.hedge_data.into_iter().map(|n| n.unwrap()).collect(),
        };

        HedgeGraph {
            base_nodes: nodes.len(),
            nodes,
            node_data,
            involution: new_inv,
        }
    }

    pub fn base_dot(&self) -> String {
        let mut out = "digraph {\n ".to_string();
        out.push_str("  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\";\n ");
        for (i, (e, n)) in self.involution.iter() {
            out.push_str(
                &e.default_dot(
                    i,
                    Some(n.0),
                    self.involved_node_hairs(i)
                        .map(|x| self.id_from_hairs(x).unwrap().0),
                    None,
                ),
            );
        }
        out += "}";
        out
    }

    pub fn base_dot_underlying(&self) -> String {
        let mut out = "digraph {\n ".to_string();
        out.push_str("  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\";\n ");
        for (i, (e, n)) in self.involution.iter() {
            let m = match e {
                InvolutiveMapping::Identity { data, underlying } => {
                    let mut d = data.as_ref();
                    d.orientation = (*underlying).into();
                    InvolutiveMapping::Identity {
                        data: d,
                        underlying: *underlying,
                    }
                }
                InvolutiveMapping::Source { sink_idx, data } => {
                    let mut d = data.as_ref();
                    d.orientation = Orientation::Default;
                    InvolutiveMapping::Source {
                        sink_idx: *sink_idx,
                        data: d,
                    }
                }
                InvolutiveMapping::Sink { source_idx } => InvolutiveMapping::Sink {
                    source_idx: *source_idx,
                },
            };
            out.push_str(
                &m.default_dot(
                    i,
                    Some(n.0),
                    self.involved_node_hairs(i)
                        .map(|x| self.id_from_hairs(x).unwrap().0),
                    None,
                ),
            );
        }
        out += "}";
        out
    }

    pub fn nesting_node_from_subgraph(&self, internal_graph: InternalSubGraph) -> HedgeNode {
        let mut hairs = bitvec![usize, Lsb0; 0; self.involution.len()];

        if !internal_graph.valid::<E, V>(self) {
            panic!("Invalid subgraph")
        }

        for i in internal_graph.included_iter() {
            hairs |= &self.node_hairs(i).hairs;
        }

        HedgeNode {
            hairs: !(!hairs | &internal_graph.filter),
            internal_graph,
        }
    }

    fn nesting_node_fix(&self, node: &mut HedgeNode) {
        let mut externalhedges = bitvec![usize, Lsb0; 0; self.involution.len()];

        for i in node.internal_graph.filter.included_iter() {
            externalhedges |= &self.node_hairs(i).hairs;
        }

        node.hairs = !(!externalhedges | &node.internal_graph.filter);
    }

    pub fn paired_filter_from_pos(&self, pos: &[Hedge]) -> BitVec {
        let mut filter = bitvec![usize, Lsb0; 0; self.involution.len()];

        for &i in pos {
            filter.set(i.0, true);
            filter.set(self.involution.inv(i).0, true);
        }

        filter
    }

    // pub fn filter_from_pos(&self, pos: &[usize]) -> BitVec {
    //     Nested<T>::filter_from_pos(pos, self.involution.len())
    // }

    // pub fn nesting_node_from_pos(&self, pos: &[usize]) -> Nested<T> {
    //     self.nesting_node_from_subgraph(SubGraph::from(self.filter_from_pos(pos)))
    // }

    fn remove_externals(&self, subgraph: &mut HedgeNode) {
        let externals = self.external_filter();

        subgraph.internal_graph.filter &= !externals;
    }

    pub fn external_filter(&self) -> BitVec {
        let mut filter = bitvec![usize, Lsb0; 0; self.involution.len()];

        for (i, edge) in self.involution.inv.iter().enumerate() {
            if edge.is_identity() {
                filter.set(i, true);
            }
        }

        filter
    }

    pub fn full_filter(&self) -> BitVec {
        bitvec![usize, Lsb0; 1; self.involution.len()]
    }

    pub fn empty_filter(&self) -> BitVec {
        bitvec![usize, Lsb0; 0; self.involution.len()]
    }

    pub fn clean_subgraph(&self, filter: BitVec) -> InternalSubGraph {
        InternalSubGraph::cleaned_filter_optimist(filter, self)
    }

    pub fn full_node(&self) -> HedgeNode {
        self.nesting_node_from_subgraph(self.full_graph())
    }

    pub fn cycle_basis(&self) -> (Vec<Cycle>, TraversalTree) {
        self.paton_cycle_basis(&self.full_graph(), self.node_hairs(Hedge(0)), None)
            .unwrap()
    }

    pub fn align_superficial_to_underlying(&mut self) {
        for i in self.involution.iter_idx() {
            let orientation = self.involution.edge_data(i).orientation;
            if let Orientation::Reversed = orientation {
                self.involution.flip_underlying(i);
            }
        }
    }

    pub fn align_to_tree_underlying(&mut self, tree: &TraversalTree) {
        for (i, (_, p)) in tree.inv.iter() {
            match p {
                Parent::Root => {
                    if tree.tree.includes(&i) {
                        self.involution.set_as_sink(i)
                    } else {
                        self.involution.set_as_source(i)
                    }
                }
                Parent::Hedge {
                    hedge_to_root,
                    traversal_order,
                } => {
                    if tree.tree.includes(&i) {
                        if *hedge_to_root == i {
                            self.involution.set_as_source(i)
                        } else {
                            self.involution.set_as_sink(i)
                        }
                    } else {
                        let tord = *traversal_order;
                        if let Parent::Hedge {
                            traversal_order, ..
                        } = tree.inv.hedge_data(tree.inv.inv(i))
                        {
                            if tord > *traversal_order {
                                self.involution.set_as_sink(i);
                            }
                        }
                    }
                }
                Parent::Unset => {
                    println!("unset{i}");
                }
            }
        }
    }

    pub fn align_to_tree_superficial(&mut self, tree: &TraversalTree) {
        for (i, (_, p)) in tree.inv.iter() {
            match &mut self.involution[i] {
                InvolutiveMapping::Identity { data, .. } => match p {
                    Parent::Root => {}
                    Parent::Hedge { hedge_to_root, .. } => {
                        if *hedge_to_root == i {
                            data.orientation = Orientation::Default;
                        } else {
                            data.orientation = Orientation::Reversed;
                        }
                    }
                    Parent::Unset => {}
                },
                InvolutiveMapping::Source { data, .. } => match p {
                    Parent::Root => {}
                    Parent::Hedge { hedge_to_root, .. } => {
                        if *hedge_to_root == i {
                            data.orientation = Orientation::Default;
                        } else {
                            data.orientation = Orientation::Reversed;
                        }
                    }
                    Parent::Unset => {}
                },
                _ => {}
            }
        }
    }

    pub fn order_basis(&self, basis: &[HedgeNode]) -> Vec<Vec<InternalSubGraph>> {
        let mut seen = vec![basis[0].internal_graph.clone()];
        let mut partitions = vec![seen.clone()];

        for cycle in basis.iter() {
            if seen
                .iter()
                .any(|p| !p.empty_intersection(&cycle.internal_graph))
            {
                partitions
                    .last_mut()
                    .unwrap()
                    .push(cycle.internal_graph.clone());
            } else {
                for p in partitions.last().unwrap() {
                    seen.push(p.clone());
                }
                partitions.push(vec![cycle.internal_graph.clone()]);
            }
        }

        partitions
    }

    pub fn all_cycles(&self) -> Vec<Cycle> {
        Cycle::all_sum_powerset_filter_map(&self.cycle_basis().0, &|mut c| {
            if c.is_circuit(self) {
                c.loop_count = Some(1);
                Some(c)
            } else {
                None
            }
        })
        .unwrap()
        .into_iter()
        .collect()
    }

    pub fn all_cycle_sym_diffs(&self) -> Result<Vec<InternalSubGraph>, TryFromIntError> {
        Cycle::all_sum_powerset_filter_map(&self.cycle_basis().0, &Some)
            .map(|a| a.into_iter().map(|c| c.internal_graph(self)).collect())
    }

    pub fn all_cycle_unions(&self) -> AHashSet<InternalSubGraph> {
        InternalSubGraph::all_unions_iterative(&self.all_cycle_sym_diffs().unwrap())
    }

    pub fn full_graph(&self) -> InternalSubGraph {
        InternalSubGraph::cleaned_filter_optimist(self.full_filter(), self)
    }

    pub fn empty_subgraph<S: SubGraph>(&self) -> S {
        S::empty(self.n_hedges())
    }

    pub fn paton_count_loops(
        &self,
        subgraph: &InternalSubGraph,
        start: &HedgeNode,
    ) -> Result<usize, HedgeError> {
        let tree = TraversalTree::dfs(self, subgraph, start, None);

        let cuts = subgraph.subtract(&tree.tree);
        Ok(self.involution.n_internals(&cuts))
    }

    pub fn number_of_nodes_in_subgraph(&self, subgraph: &InternalSubGraph) -> usize {
        self.iter_node_data(subgraph).count()
    }

    pub fn node_degrees_in_subgraph(&self, subgraph: &InternalSubGraph) -> AHashMap<usize, usize> {
        let mut degrees = AHashMap::new();

        for (node, _) in self.iter_node_data(subgraph) {
            let node_pos = self.id_from_hairs(node).unwrap().0;

            // Count the number of edges in the subgraph incident to this node
            let incident_edges = node.hairs.clone() & &subgraph.filter;
            let degree = incident_edges.count_ones();

            degrees.insert(node_pos, degree);
        }

        degrees
    }

    pub fn hairy_from_filter(&self, filter: BitVec) -> HedgeNode {
        self.nesting_node_from_subgraph(InternalSubGraph::cleaned_filter_pessimist(filter, self))
    }

    pub fn paton_cycle_basis(
        &self,
        subgraph: &InternalSubGraph,
        start: &HedgeNode,
        included_hedge: Option<Hedge>,
    ) -> Result<(Vec<Cycle>, TraversalTree), HedgeError> {
        if !subgraph.intersects(&start.hairs) {
            return Err(HedgeError::InvalidStart);
        }

        let tree = TraversalTree::dfs(self, subgraph, start, included_hedge);

        let cuts = subgraph.subtract(&tree.tree);

        let mut cycle_basis = Vec::new();

        for c in cuts.included_iter() {
            if c > self.involution.inv(c) {
                cycle_basis.push(tree.cycle(c).unwrap());
            }
        }

        Ok((cycle_basis, tree))
    }

    pub fn dot_impl<S: SubGraph>(
        &self,
        node_as_graph: &S,
        graph_info: String,
        edge_attr: &impl Fn(&E) -> Option<String>,
        node_attr: &impl Fn(&V) -> Option<String>,
    ) -> String {
        node_as_graph.dot(self, graph_info, edge_attr, node_attr)
    }

    pub fn dot<S: SubGraph>(&self, node_as_graph: &S) -> String {
        self.dot_impl(node_as_graph, "".to_string(), &|_| None, &|_| None)
    }

    pub fn cut_branches(&self, subgraph: &mut HedgeNode) {
        let nodes = AHashSet::<&HedgeNode>::from_iter(
            subgraph
                .internal_graph
                .included_iter()
                .map(|i| self.node_hairs(i)),
        );
        self.remove_externals(subgraph);

        let mut has_branch = true;
        while has_branch {
            has_branch = false;

            for n in &nodes {
                let int = n.hairs.clone() & &subgraph.internal_graph.filter;

                if int.count_ones() == 1 {
                    subgraph
                        .internal_graph
                        .filter
                        .set(int.first_one().unwrap(), false);
                    subgraph.internal_graph.filter.set(
                        self.involution.inv(Hedge(int.first_one().unwrap())).0,
                        false,
                    );
                    has_branch = true;
                }
            }
        }

        self.nesting_node_fix(subgraph);
    }

    pub fn get_edge_data(&self, edge: Hedge) -> &E {
        self.involution.edge_data(edge).data.as_ref().unwrap()
    }

    pub fn all_spinneys_with_basis(&self, basis: &[&InternalSubGraph]) -> AHashSet<HedgeNode> {
        let mut spinneys = AHashSet::new();
        let mut base_cycle: InternalSubGraph = self.empty_subgraph();

        for cycle in basis {
            base_cycle.sym_diff_with(cycle);
        }

        spinneys.insert(self.nesting_node_from_subgraph(base_cycle.clone()));

        if basis.len() == 1 {
            return spinneys;
        }

        for i in 0..basis.len() {
            for s in self.all_spinneys_with_basis(
                &basis
                    .iter()
                    .enumerate()
                    .filter_map(|(j, s)| if j != i { Some(*s) } else { None })
                    .collect_vec(),
            ) {
                spinneys
                    .insert(self.nesting_node_from_subgraph(s.internal_graph.union(&base_cycle)));
                spinneys.insert(s);
            }
        }

        spinneys
    }

    pub fn all_spinneys_rec(&self, spinneys: &mut AHashSet<HedgeNode>, cycle_sums: Vec<HedgeNode>) {
        let _len = spinneys.len();

        let mut pset = PowersetIterator::new(cycle_sums.len() as u8);

        pset.next(); //Skip empty set

        for (ci, cj) in cycle_sums.iter().tuple_combinations() {
            let _union = ci.internal_graph.union(&cj.internal_graph);

            // spinneys.insert(union);
        }
    }

    pub fn all_spinneys(
        &self,
    ) -> AHashMap<InternalSubGraph, Vec<(InternalSubGraph, Option<InternalSubGraph>)>> {
        let basis_cycles = self.cycle_basis().0;

        let mut all_combinations = PowersetIterator::new(basis_cycles.len() as u8);
        all_combinations.next(); //Skip empty set

        let mut spinneys: AHashMap<
            InternalSubGraph,
            Vec<(InternalSubGraph, Option<InternalSubGraph>)>,
        > = AHashMap::new();

        let mut cycles: Vec<InternalSubGraph> = Vec::new();
        for p in all_combinations {
            let mut base_cycle: InternalSubGraph = self.empty_subgraph();

            for i in p.iter_ones() {
                base_cycle.sym_diff_with(&basis_cycles[i].clone().internal_graph(self));
            }

            cycles.push(base_cycle);
        }

        for (ci, cj) in cycles.iter().tuple_combinations() {
            let union = ci.union(cj);

            if let Some(v) = spinneys.get_mut(&union) {
                v.push((ci.clone(), Some(cj.clone())));
            } else {
                spinneys.insert(union, vec![(ci.clone(), Some(cj.clone()))]);
            }
        }

        for c in cycles {
            spinneys.insert(c.clone(), vec![(c.clone(), None)]);
        }
        spinneys
    }

    pub fn all_spinneys_alt(&self) -> AHashSet<InternalSubGraph> {
        let mut spinneys = AHashSet::new();
        let cycles = self.all_cycles();

        let mut pset = PowersetIterator::new(cycles.len() as u8);
        pset.next(); //Skip empty set

        for p in pset {
            let mut union: InternalSubGraph = self.empty_subgraph();

            for i in p.iter_ones() {
                union.union_with(&cycles[i].clone().internal_graph(self));
            }

            spinneys.insert(union);
        }

        for c in cycles {
            spinneys.insert(c.internal_graph(self));
        }
        spinneys
    }

    pub fn all_s_t_cuts(
        &self,
        s: &HedgeNode,
        t: &HedgeNode,
        regions: &mut AHashSet<InternalSubGraph>,
    ) {
        let mut new_internals = vec![];
        for h in s.hairs.included_iter() {
            let invh = self.involution.inv(h);

            if h > invh && s.hairs.includes(&self.involution.inv(h)) {
                new_internals.push(h);
            }
        }

        let mut new_node = s.clone();

        for h in new_internals {
            new_node.hairs.set(h.0, false);
            new_node.hairs.set(self.involution.inv(h).0, false);
            new_node.internal_graph.filter.set(h.0, true);
            new_node
                .internal_graph
                .filter
                .set(self.involution.inv(h).0, true);
        }

        let complement = new_node.complement(self).hairs;

        let count = self.count_connected_components(&complement);

        if count == 1 && !regions.insert(new_node.internal_graph.clone()) {
            return;
        }

        for h in new_node.included_iter() {
            let invh = self.involution.inv(h);

            if invh != h && !t.hairs.includes(&invh) {
                let mut new_node = s.clone();
                new_node.hairs.union_with(&self.node_hairs(invh).hairs);

                new_node.hairs.set(h.0, false);
                new_node.hairs.set(invh.0, false);
                new_node.internal_graph.filter.set(h.0, true);
                new_node.internal_graph.filter.set(invh.0, true);
                self.all_s_t_cuts(&new_node, t, regions);
            }
        }
    }
}

pub struct TraversalTree {
    pub traversal: Vec<Hedge>,
    inv: Involution<Parent, ()>,
    pub tree: InternalSubGraph,
    pub covers: BitVec,
}

pub enum Parent {
    Unset,
    Root,
    Hedge {
        hedge_to_root: Hedge,
        traversal_order: usize,
    },
}

impl TraversalTree {
    fn covers(&self) -> BitVec {
        let mut covers = <BitVec as SubGraph>::empty(self.inv.inv.len());
        for (i, m) in self.inv.hedge_data.iter().enumerate() {
            match m {
                Parent::Unset => {}
                _ => {
                    covers.set(i, true);
                }
            }
        }
        covers
    }

    fn path_to_root(&self, start: Hedge) -> BitVec {
        let mut path = <BitVec as SubGraph>::empty(self.inv.inv.len());
        let mut current = start;
        path.set(current.0, true);

        while let Parent::Hedge { hedge_to_root, .. } = self.inv.hedge_data(current) {
            path.set(hedge_to_root.0, true);
            current = self.inv.inv(*hedge_to_root);
            path.set(current.0, true);
        }
        path
    }

    pub fn cycle(&self, cut: Hedge) -> Option<Cycle> {
        match self.inv.hedge_data(cut) {
            Parent::Hedge { hedge_to_root, .. } => {
                if *hedge_to_root == cut {
                    //if cut is in the tree, no cycle can be formed
                    return None;
                }
            }
            Parent::Root => {}
            _ => return None,
        }

        let cut_pair = self.inv.inv(cut);
        match self.inv.hedge_data(cut_pair) {
            Parent::Hedge { hedge_to_root, .. } => {
                if *hedge_to_root == cut {
                    //if cut is in the tree,no cycle can be formed
                    return None;
                }
            }
            Parent::Root => {}
            _ => return None,
        }

        let mut cycle = self.path_to_root(cut);
        cycle.sym_diff_with(&self.path_to_root(cut_pair));
        let mut cycle = Cycle::new_unchecked(cycle);
        cycle.loop_count = Some(1);

        Some(cycle)
    }

    pub fn bfs<E, V, S: SubGraph>(
        graph: &HedgeGraph<E, V>,
        subgraph: &S,
        root_node: &HedgeNode,
        // target: Option<&HedgeNode>,
    ) -> Self {
        let mut queue = VecDeque::new();
        let mut seen = subgraph.hairs(root_node);

        let mut traversal: Vec<Hedge> = Vec::new();
        let mut involution: Involution<Parent, ()> = graph
            .involution
            .forgetful_map_node_data_ref(|_| Parent::Unset);

        // add all hedges from root node that are not self loops
        // to the queue
        // They are all potential branches
        for i in seen.included_iter() {
            involution.set_hedge_data(i, Parent::Root);
            if !seen.includes(&graph.involution.inv(i)) {
                // if not self loop
                queue.push_back(i)
            }
        }
        while let Some(hedge) = queue.pop_front() {
            // if the hedge is not external get the neighbors of the paired hedge
            if let Some(cn) = graph.connected_neighbors(subgraph, hedge) {
                let connected = involution.inv(hedge);

                if !seen.includes(&connected) && subgraph.includes(&connected) {
                    // if this new hedge hasn't been seen before, it means the node it belongs to
                    //  a new node in the traversal
                    traversal.push(connected);
                } else {
                    continue;
                }
                // mark the new node as seen
                seen.union_with(&cn);

                // for all hedges in this new node, they have a parent, the initial hedge
                for i in cn.included_iter() {
                    if let Parent::Unset = involution.hedge_data(i) {
                        involution.set_hedge_data(
                            i,
                            Parent::Hedge {
                                hedge_to_root: connected,
                                traversal_order: traversal.len(),
                            },
                        );
                    }
                    // if they lead to a new node, they are potential branches, add them to the queue
                    if !seen.includes(&involution.inv(i)) && subgraph.includes(&i) {
                        queue.push_back(i);
                    }
                }
            }
        }

        TraversalTree::new(graph, traversal, seen, involution)
    }

    pub fn new<E, V>(
        graph: &HedgeGraph<E, V>,
        traversal: Vec<Hedge>,
        covers: BitVec,
        inv: Involution<Parent, ()>,
    ) -> Self {
        let mut tree = graph.empty_filter();

        for (i, j) in traversal.iter().map(|x| (*x, inv.inv(*x))) {
            tree.set(i.0, true);
            tree.set(j.0, true);
        }

        TraversalTree {
            traversal,
            covers,
            inv,
            tree: InternalSubGraph::cleaned_filter_optimist(tree, graph),
        }
    }

    pub fn dfs<E, V, S: SubGraph>(
        graph: &HedgeGraph<E, V>,
        subgraph: &S,
        root_node: &HedgeNode,
        include_hegde: Option<Hedge>,
        // target: Option<&HedgeNode>,
    ) -> Self {
        let mut stack = Vec::new();
        let mut seen = subgraph.hairs(root_node);

        let mut traversal: Vec<Hedge> = Vec::new();
        let mut involution: Involution<Parent, ()> = graph
            .involution
            .forgetful_map_node_data_ref(|_| Parent::Unset);

        let mut included_hedge_is_possible = false;

        // add all hedges from root node that are not self loops
        // to the stack
        // They are all potential branches

        for i in seen.included_iter() {
            involution.set_hedge_data(i, Parent::Root);
            if !seen.includes(&graph.involution.inv(i)) {
                // if not self loop
                if let Some(hedge) = include_hegde {
                    if hedge != i {
                        stack.push(i);
                    } else {
                        println!("skipping{i}");
                        included_hedge_is_possible = true;
                    }
                } else {
                    stack.push(i);
                }
            }
        }

        if included_hedge_is_possible {
            stack.push(include_hegde.unwrap());
        }
        while let Some(hedge) = stack.pop() {
            // println!("looking at {hedge}");
            // if the hedge is not external get the neighbors of the paired hedge
            if let Some(cn) = graph.connected_neighbors(subgraph, hedge) {
                let connected = involution.inv(hedge);

                if !seen.includes(&connected) && subgraph.includes(&connected) {
                    // if this new hedge hasn't been seen before, it means the node it belongs to
                    // is a new node in the traversal
                    traversal.push(connected);
                } else {
                    continue;
                }

                // mark the new node as seen
                seen.union_with(&cn);

                for i in cn.included_iter() {
                    if let Parent::Unset = involution.hedge_data(i) {
                        involution.set_hedge_data(
                            i,
                            Parent::Hedge {
                                hedge_to_root: connected,
                                traversal_order: traversal.len(),
                            },
                        );
                    }

                    if !seen.includes(&involution.inv(i)) && subgraph.includes(&i) {
                        stack.push(i);
                    }
                }
            }
        }

        TraversalTree::new(graph, traversal, seen, involution)
    }
}

pub mod drawing;
pub mod layout;
#[cfg(test)]
mod tests;
