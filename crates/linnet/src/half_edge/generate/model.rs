//! Dense form of the generation input: stub types, the vertex automaton, external-leg classes and
//! edge labels.
//!
//! A stub type is `3 * field + direction` (undirected, incoming, outgoing), with fields ranked in
//! their `Ord` order. A partial vertex is a state of the [`Automaton`] built from the signatures:
//! the multiset of the stubs placed at it so far, while some signature of its degree contains it.
//! Adding a stub is a lookup in the state's transitions, and the edge labels whose stubs a state
//! can still take are precomputed per state as bit masks, so the dressing search visits only the
//! labels that fit both ends of an edge.
//!
//! External legs are packed by label class (interchangeable legs: equal label and stub) on complete
//! topologies, and by slot class while the topology search places them, where distinct labels with
//! a common stub can share an anonymous class.

use ahash::AHashMap;

use super::{
    budget::{Budgets, Request},
    GenerationError,
};

/// Largest vertex degree: a [`StubSet`] holds 15 stubs and its degree.
pub(crate) const MAX_DEGREE: usize = 15;
/// Largest number of fields: stub types are bytes and `StubSet` stores them shifted by one.
pub(crate) const MAX_FIELDS: usize = 84;
pub(crate) const UNDIRECTED: u8 = 0;
pub(crate) const INCOMING: u8 = 1;
pub(crate) const OUTGOING: u8 = 2;

/// Sorted multiset of stub types at a vertex of known degree; byte `15` is the degree and bytes
/// `0..len` the sorted types plus one, so the packed key is exact. It names the states of the
/// [`Automaton`] while they are interned.
#[derive(Clone, Copy, Debug)]
struct StubSet {
    bytes: [u8; 16],
}

impl StubSet {
    fn new(degree: u8) -> Self {
        let mut bytes = [0; 16];
        bytes[15] = degree;
        StubSet { bytes }
    }

    fn insert(&mut self, stub: u8) {
        let b = stub + 1;
        let i = self.bytes[..15]
            .iter()
            .position(|&x| x == 0 || x > b)
            .expect("a vertex holds at most its degree in stubs");
        self.bytes.copy_within(i..14, i + 1);
        self.bytes[i] = b;
    }

    fn key(&self) -> u128 {
        u128::from_le_bytes(self.bytes)
    }
}

/// A state of the [`Automaton`]: a partial vertex.
pub(crate) type VertexState = u32;

/// No state: the stubs fit no signature.
const NONE: VertexState = VertexState::MAX;

/// No stub type: they are below `3 * MAX_FIELDS`.
const EMPTY: u8 = u8::MAX;

/// The partial vertices of every signature degree as the states of one automaton, a component per
/// degree. A state is a multiset of stubs that some signature of its degree contains, and adding a
/// stub leads to the state with that stub added while a signature of the degree still contains it.
/// The roots are the empty multisets; the full states (as many stubs as their degree) are the
/// signatures themselves, and their vertex rules are the indices of the equal signatures.
///
/// Most states take few of the stub types, so each state keeps its transitions in a small hash
/// table of its own, keyed by stub type with linear probing, at most half full: memory grows with
/// the transitions rather than with the states times the stub types.
///
/// Signatures may carry integer weights (coupling orders, counterterm orders, any tag); every state
/// then carries, per weight, the least and the most weight of the signatures containing it, which
/// bound the weight of every vertex completing it.
#[derive(Debug)]
struct Automaton {
    /// Per degree `0..=MAX_DEGREE`: its root, `NONE` without a signature of that degree.
    roots: [VertexState; MAX_DEGREE + 1],
    /// Per state: the first slot of its table and its size less one, a power of two less one.
    tables: Vec<(u32, u32)>,
    /// `(stub, target)`: a transition by `stub` to `target`, or `(EMPTY, NONE)`.
    slots: Vec<(u8, VertexState)>,
    /// The signatures equal to state `s` are `rules[starts[s]..starts[s + 1]]`, ascending.
    starts: Vec<u32>,
    rules: Vec<u32>,
    /// Per state and weight (`n_weights` per state): its least and most weight.
    weights: Vec<[i32; 2]>,
    n_weights: usize,
}

/// A transition of the [`Automaton`]: `(state, stub, target)`.
type Transition = (VertexState, u8, VertexState);

impl Automaton {
    /// The automaton of `signatures`, given as stub types below `stubs` and of degree
    /// `1..=MAX_DEGREE`, with every transition once. `weights` holds `n_weights` weights per
    /// signature.
    fn new(
        signatures: &[Vec<u8>],
        stubs: usize,
        weights: &[i32],
        n_weights: usize,
    ) -> (Self, Vec<Transition>) {
        let width = stubs.div_ceil(64);
        let mut ids: AHashMap<u128, VertexState> = AHashMap::new();
        // per state, `width` words of a bit per stub type it takes, to list each transition once
        let mut takes: Vec<u64> = vec![];
        let mut state_weights: Vec<[i32; 2]> = vec![];
        let mut transitions = vec![];
        let mut roots = [NONE; MAX_DEGREE + 1];
        let mut full = vec![];
        // the sub-multisets of a signature, by index
        let mut states = vec![];
        for (i, signature) in signatures.iter().enumerate() {
            let degree = signature.len() as u8;
            let mut sorted = signature.clone();
            sorted.sort_unstable();
            // every sub-multiset: choose how many of each run of equal stubs to keep, `keep[r]` of
            // run `r`, visited in the order of the index `sum keep[r] * stride[r]`
            let runs: Vec<(u8, usize)> = sorted
                .chunk_by(|a, b| a == b)
                .map(|run| (run[0], run.len()))
                .collect();
            let stride: Vec<usize> = runs
                .iter()
                .scan(1, |n, &(_, len)| {
                    let stride = *n;
                    *n *= len + 1;
                    Some(stride)
                })
                .collect();
            let mut keep = vec![0usize; runs.len()];
            states.clear();
            loop {
                let mut set = StubSet::new(degree);
                for (&(stub, _), &k) in runs.iter().zip(&keep) {
                    for _ in 0..k {
                        set.insert(stub);
                    }
                }
                let id = ids.len() as VertexState;
                let state = *ids.entry(set.key()).or_insert(id);
                if state == id {
                    takes.resize(takes.len() + width, 0);
                    state_weights.resize(state_weights.len() + n_weights, [i32::MAX, i32::MIN]);
                }
                let own = &weights[i * n_weights..(i + 1) * n_weights];
                for (bounds, &w) in state_weights[state as usize * n_weights..]
                    .iter_mut()
                    .zip(own)
                {
                    *bounds = [bounds[0].min(w), bounds[1].max(w)];
                }
                // the sub-multiset with one stub of run `r` fewer takes that stub to this one
                for (r, &(stub, _)) in runs.iter().enumerate() {
                    if keep[r] > 0 {
                        let from = states[states.len() - stride[r]];
                        let word = &mut takes[from as usize * width + stub as usize / 64];
                        let bit = 1 << (stub % 64);
                        if *word & bit == 0 {
                            *word |= bit;
                            transitions.push((from, stub, state));
                        }
                    }
                }
                states.push(state);
                let Some(r) = (0..runs.len()).find(|&r| keep[r] < runs[r].1) else {
                    full.push((state, i as u32));
                    break;
                };
                keep[r] += 1;
                keep[..r].fill(0);
            }
            roots[degree as usize] = states[0];
        }
        let n_states = ids.len();
        // a table per state, of the least power of two that is at least twice its transitions, so
        // a probe always ends at an empty slot
        let mut tables = vec![(0, 0u32); n_states];
        for &(state, ..) in &transitions {
            tables[state as usize].1 += 2;
        }
        let mut n_slots = 0;
        for (first, size) in &mut tables {
            let len = size.next_power_of_two();
            (*first, *size) = (n_slots, len - 1);
            n_slots += len;
        }
        let mut slots = vec![(EMPTY, NONE); n_slots as usize];
        for &(state, stub, target) in &transitions {
            let (first, mask) = tables[state as usize];
            let mut i = stub as u32 & mask;
            while slots[(first + i) as usize].0 != EMPTY {
                i = (i + 1) & mask;
            }
            slots[(first + i) as usize] = (stub, target);
        }
        // signatures are visited in order, so the rules of a state stay ascending
        full.sort_by_key(|&(state, _)| state);
        let mut starts = vec![0u32; n_states + 1];
        for &(state, _) in &full {
            starts[state as usize + 1] += 1;
        }
        for s in 0..n_states {
            starts[s + 1] += starts[s];
        }
        let automaton = Automaton {
            roots,
            tables,
            slots,
            starts,
            rules: full.into_iter().map(|(_, rule)| rule).collect(),
            weights: state_weights,
            n_weights,
        };
        (automaton, transitions)
    }

    fn len(&self) -> usize {
        self.starts.len() - 1
    }

    /// The empty vertex of this degree, if a signature has the degree.
    fn root(&self, degree: u8) -> Option<VertexState> {
        let root = self.roots[degree as usize];
        (root != NONE).then_some(root)
    }

    /// The state with `stub` added, if a signature of its degree contains it.
    #[inline(always)]
    fn step(&self, state: VertexState, stub: u8) -> Option<VertexState> {
        let (first, mask) = self.tables[state as usize];
        let mut i = stub as u32 & mask;
        loop {
            match self.slots[(first + i) as usize] {
                (key, target) if key == stub => return Some(target),
                (EMPTY, _) => return None,
                _ => i = (i + 1) & mask,
            }
        }
    }

    /// The signatures equal to `state`, ascending: empty unless the state is full.
    fn rules(&self, state: VertexState) -> &[u32] {
        let s = state as usize;
        &self.rules[self.starts[s] as usize..self.starts[s + 1] as usize]
    }

    /// The least and the most weight `k` of the signatures containing `state`.
    #[inline(always)]
    fn weight(&self, state: VertexState, k: usize) -> [i32; 2] {
        self.weights[state as usize * self.n_weights + k]
    }

    /// `label` fits an edge from a vertex in state `a` to a higher one in state `b`.
    fn label_fits(&self, label: &Label, [a, b]: [VertexState; 2]) -> bool {
        self.step(a, label.at_u).is_some() && self.step(b, label.at_v).is_some()
    }

    /// `label` fits a self-loop at a vertex in state `a`.
    fn loop_fits(&self, label: &LoopLabel, a: VertexState) -> bool {
        let [s, t] = label.stubs;
        self.step(a, s).and_then(|a| self.step(a, t)).is_some()
    }
}

/// The labels a state can take, a bit per label: paired labels by the stub they put at the lower
/// vertex of their pair (`at_u`) and at the higher one (`at_v`), and self-loops by both stubs.
#[derive(Clone, Copy, Debug, Default)]
struct Masks {
    at_u: u128,
    at_v: u128,
    loops: u128,
}

impl Masks {
    /// The masks of every state of `vertices`, whose transitions are `transitions`, if the labels
    /// fit into them; a state takes the labels of its transitions.
    fn of(
        vertices: &Automaton,
        transitions: &[Transition],
        labels: &[Label],
        loop_labels: &[LoopLabel],
    ) -> Option<Vec<Self>> {
        let width = u128::BITS as usize;
        if labels.len() > width || loop_labels.len() > width {
            return None;
        }
        // per stub type (a byte): the labels putting it at either vertex, and the self-loops
        // starting with it
        let mut by_stub = vec![Masks::default(); 1 << u8::BITS];
        for (l, label) in labels.iter().enumerate() {
            by_stub[label.at_u as usize].at_u |= 1 << l;
            by_stub[label.at_v as usize].at_v |= 1 << l;
        }
        for (l, label) in loop_labels.iter().enumerate() {
            by_stub[label.stubs[0] as usize].loops |= 1 << l;
        }
        let mut masks = vec![Masks::default(); vertices.len()];
        for &(state, stub, target) in transitions {
            let labels = by_stub[stub as usize];
            let mask = &mut masks[state as usize];
            mask.at_u |= labels.at_u;
            mask.at_v |= labels.at_v;
            for l in Fitting::Mask(labels.loops) {
                if vertices.step(target, loop_labels[l].stubs[1]).is_some() {
                    mask.loops |= 1 << l;
                }
            }
        }
        Some(masks)
    }
}

/// The labels from some label on that fit an edge, ascending: the set bits of the intersection of
/// the masks of its ends, or, for models without masks, a scan of the labels.
pub(crate) enum Fitting<'m> {
    Mask(u128),
    Scan(Scan<'m>),
}

/// The labels from `next` on whose stubs both ends of an edge can take.
pub(crate) struct Scan<'m> {
    model: &'m Model,
    ends: [VertexState; 2],
    self_loop: bool,
    next: usize,
}

impl Iterator for Fitting<'_> {
    type Item = usize;

    #[inline(always)]
    fn next(&mut self) -> Option<usize> {
        match self {
            Fitting::Mask(bits) => {
                if *bits == 0 {
                    return None;
                }
                let l = bits.trailing_zeros() as usize;
                *bits &= *bits - 1;
                Some(l)
            }
            Fitting::Scan(scan) => scan.next(),
        }
    }
}

impl Iterator for Scan<'_> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        let Scan {
            model,
            ends,
            self_loop,
            next,
        } = *self;
        let vertices = &model.vertices;
        let fits = |&l: &usize| {
            if self_loop {
                vertices.loop_fits(&model.loop_labels[l], ends[0])
            } else {
                vertices.label_fits(&model.labels[l], ends)
            }
        };
        let end = if self_loop {
            model.loop_labels.len()
        } else {
            model.labels.len()
        };
        let l = (next..end).find(fits);
        self.next = l.map_or(end, |l| l + 1);
        l
    }
}

/// A paired edge between two distinct vertices `u < v` as seen from the pair: the stub it puts at
/// `u` and at `v`, its image when `u` and `v` swap, and how to materialize it.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Label {
    pub(crate) at_u: u8,
    pub(crate) at_v: u8,
    pub(crate) flip: u8,
    pub(crate) field: u8,
    pub(crate) direction: Direction,
}

/// A self-loop: the two stubs it puts at its vertex.
#[derive(Clone, Copy, Debug)]
pub(crate) struct LoopLabel {
    pub(crate) stubs: [u8; 2],
    pub(crate) field: u8,
    pub(crate) directed: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Direction {
    Undirected,
    /// From the lower to the higher vertex of the pair.
    Forward,
    Backward,
}

/// External legs grouped into classes of interchangeable legs (equal label and stub). A multiset
/// of legs is packed into a `u64`: class `c` owns `width[c]` bits at `offset[c]`, enough for its
/// multiplicity, so sums and differences of packed multisets never carry.
#[derive(Clone, Debug, Default)]
pub(crate) struct Legs {
    pub(crate) stub: Vec<u8>,
    pub(crate) mult: Vec<u8>,
    offset: Vec<u8>,
    width: Vec<u8>,
    /// All legs, packed.
    pub(crate) all: u64,
    pub(crate) total: u32,
}

impl Legs {
    pub(crate) fn new(classes: &[(u8, u8)]) -> Result<Self, GenerationError> {
        let mut legs = Legs::default();
        let mut offset = 0u32;
        for &(stub, mult) in classes {
            let width = u8::BITS - mult.leading_zeros();
            if offset + width > u64::BITS {
                return Err(GenerationError::TooManyLegs);
            }
            legs.stub.push(stub);
            legs.mult.push(mult);
            legs.offset.push(offset as u8);
            legs.width.push(width as u8);
            legs.all |= (mult as u64) << offset;
            legs.total += mult as u32;
            offset += width;
        }
        Ok(legs)
    }

    pub(crate) fn n_classes(&self) -> usize {
        self.stub.len()
    }

    pub(crate) fn count(&self, packed: u64, class: usize) -> u8 {
        ((packed >> self.offset[class]) & ((1 << self.width[class]) - 1)) as u8
    }

    pub(crate) fn one(&self, class: usize) -> u64 {
        1 << self.offset[class]
    }

    /// Adding the legs `add` to `used` stays within `total`, class by class.
    pub(crate) fn fits(&self, used: u64, add: u64, total: u64) -> bool {
        (0..self.n_classes()).all(|c| {
            self.count(used, c) as u32 + self.count(add, c) as u32 <= self.count(total, c) as u32
        })
    }

    /// Number of legs in `packed`.
    pub(crate) fn size(&self, packed: u64) -> u32 {
        (0..self.n_classes())
            .map(|c| self.count(packed, c) as u32)
            .sum()
    }

    /// Every sub-multiset of `packed` with at most `max_len` legs, with its size.
    pub(crate) fn subsets(&self, packed: u64, max_len: u32, out: &mut Vec<(u64, u32)>) {
        out.clear();
        out.push((0, 0));
        for c in 0..self.n_classes() {
            let available = self.count(packed, c) as u32;
            for i in 0..out.len() {
                let (base, len) = out[i];
                for k in 1..=available.min(max_len.saturating_sub(len)) {
                    out.push((base + k as u64 * self.one(c), len + k));
                }
            }
        }
    }
}

/// The dense model shared by both stages.
#[derive(Debug)]
pub(crate) struct Model {
    /// Distinct signature degrees, ascending.
    pub(crate) degrees: Vec<u8>,
    vertices: Automaton,
    pub(crate) labels: Vec<Label>,
    pub(crate) loop_labels: Vec<LoopLabel>,
    /// Per state of `vertices`: the labels it can take. `None` with more than 128 labels or
    /// self-loop labels, where [`Fitting`] scans the transitions instead.
    masks: Option<Vec<Masks>>,
    /// The external legs by label class, as complete topologies carry them.
    pub(crate) legs: Legs,
    /// The external legs as the topology search places them: three or more label classes of
    /// multiplicity one with a common stub form one anonymous slot class, whose labels are assigned
    /// once the topology is shaped; every other label class is its own slot class. With two such
    /// legs, labelling costs about what the symmetric topologies save.
    pub(crate) slots: Legs,
    /// Per slot class: its label class, `None` for an anonymous one.
    slot_labels: Vec<Option<usize>>,
    /// The labels of anonymous slots in assignment order: `(slot class, label class)`.
    pub(crate) labelling: Vec<(usize, usize)>,
    /// Every vertex has at most one rule for its stubs.
    pub(crate) unique_rules: bool,
    /// The budgets of the request.
    pub(crate) budgets: Budgets,
    /// Vertices can balance the charges and parities of the legs (see [`Model::balances`]);
    /// otherwise no graph exists.
    pub(crate) balanced: bool,
}

impl Model {
    /// `signatures` and `legs` are given as stub types; legs as classes `(stub, multiplicity)`.
    pub(crate) fn new(
        n_fields: usize,
        signatures: &[Vec<u8>],
        legs: &[(u8, u8)],
        request: &Request,
    ) -> Result<Self, GenerationError> {
        if n_fields > MAX_FIELDS {
            return Err(GenerationError::TooManyFields);
        }
        let mut degrees = vec![];
        let mut present = vec![false; 3 * n_fields];
        for signature in signatures {
            let degree = signature.len();
            if degree == 0 || degree > MAX_DEGREE {
                return Err(GenerationError::UnsupportedDegree(degree));
            }
            degrees.push(degree as u8);
            for &s in signature {
                present[s as usize] = true;
            }
        }
        degrees.sort_unstable();
        degrees.dedup();
        let weights = request.signature_weights(signatures.len());
        let (vertices, transitions) =
            Automaton::new(signatures, 3 * n_fields, &weights, request.weights.len());

        let mut labels = vec![];
        let mut loop_labels = vec![];
        for field in 0..n_fields as u8 {
            let stub = |direction: u8| 3 * field + direction;
            if present[stub(UNDIRECTED) as usize] {
                let s = stub(UNDIRECTED);
                labels.push(Label {
                    at_u: s,
                    at_v: s,
                    flip: labels.len() as u8,
                    field,
                    direction: Direction::Undirected,
                });
                loop_labels.push(LoopLabel {
                    stubs: [s, s],
                    field,
                    directed: false,
                });
            }
            if present[stub(INCOMING) as usize] && present[stub(OUTGOING) as usize] {
                let (i, o) = (stub(INCOMING), stub(OUTGOING));
                let forward = labels.len() as u8;
                labels.push(Label {
                    at_u: o,
                    at_v: i,
                    flip: forward + 1,
                    field,
                    direction: Direction::Forward,
                });
                labels.push(Label {
                    at_u: i,
                    at_v: o,
                    flip: forward,
                    field,
                    direction: Direction::Backward,
                });
                loop_labels.push(LoopLabel {
                    stubs: [i, o],
                    field,
                    directed: true,
                });
            }
        }
        let unique_rules = (0..vertices.len() as VertexState).all(|s| vertices.rules(s).len() <= 1);
        // slot classes in the order of their first label class
        let anonymous = |&(stub, mult): &(u8, u8)| {
            mult == 1 && legs.iter().filter(|&&other| other == (stub, 1)).count() > 2
        };
        let mut slots: Vec<(u8, u8)> = vec![];
        let mut slot_labels = vec![];
        let mut labelling = vec![];
        for (c, class) in legs.iter().enumerate() {
            if !anonymous(class) {
                slots.push(*class);
                slot_labels.push(Some(c));
                continue;
            }
            let s = match (0..slots.len())
                .find(|&s| slot_labels[s].is_none() && slots[s].0 == class.0)
            {
                Some(s) => {
                    slots[s].1 += 1;
                    s
                }
                None => {
                    slots.push(*class);
                    slot_labels.push(None);
                    slots.len() - 1
                }
            };
            labelling.push((s, c));
        }
        Ok(Model {
            degrees,
            masks: Masks::of(&vertices, &transitions, &labels, &loop_labels),
            budgets: Budgets::new(request, weights, signatures, &labels, &loop_labels, legs)?,
            balanced: Self::balances(n_fields, signatures, legs),
            vertices,
            labels,
            loop_labels,
            legs: Legs::new(legs)?,
            slots: Legs::new(&slots)?,
            slot_labels,
            labelling,
            unique_rules,
        })
    }

    /// Whether vertices can balance the legs. The two stubs of an edge are an incoming and an
    /// outgoing stub of a charged field or two stubs of a self-conjugate one, so summed over the
    /// vertices, the flow of every charged field (incoming less outgoing stubs) equals that of the
    /// legs, and the stubs of every field have the parity of the legs'. The legs' flows must then
    /// lie in the span of the signatures' flows, and their parities in the span of the signatures'
    /// parities, checked modulo a large prime and modulo two: the legs carry none of the charges
    /// and parities that all signatures conserve.
    fn balances(n_fields: usize, signatures: &[Vec<u8>], legs: &[(u8, u8)]) -> bool {
        // a column per charged field: those with a directed stub
        let mut columns = vec![usize::MAX; n_fields];
        let mut n_charged = 0;
        let leg_stubs = legs.iter().map(|(stub, _)| stub);
        for &stub in signatures.iter().flatten().chain(leg_stubs) {
            let column = &mut columns[stub as usize / 3];
            if stub % 3 != UNDIRECTED && *column == usize::MAX {
                *column = n_charged;
                n_charged += 1;
            }
        }
        // the flows over the charged fields and the parities of all fields, a bit each
        fn vectors(
            stubs: impl Iterator<Item = (u8, u8)>,
            columns: &[usize],
            flows: &mut [i64],
        ) -> u128 {
            flows.fill(0);
            let mut parities = 0;
            for (stub, mult) in stubs {
                parities ^= (mult as u128 & 1) << (stub / 3);
                let column = columns[stub as usize / 3];
                match stub % 3 {
                    INCOMING => flows[column] += mult as i64,
                    OUTGOING => flows[column] -= mult as i64,
                    _ => {}
                }
            }
            parities
        }
        let mut flows = vec![0; n_charged];
        let mut signature_flows = vec![];
        let mut parity_basis = [0u128; 128];
        for signature in signatures {
            let stubs = signature.iter().map(|&stub| (stub, 1));
            let parities = reduce_bits(&parity_basis, vectors(stubs, &columns, &mut flows));
            if parities != 0 {
                parity_basis[127 - parities.leading_zeros() as usize] = parities;
            }
            if flows.iter().any(|&f| f != 0) {
                signature_flows.push(flows.clone());
            }
        }
        signature_flows.sort_unstable();
        signature_flows.dedup();
        let mut span = Span::new(n_charged);
        for flows in &signature_flows {
            span.insert(flows);
        }
        let parities = vectors(legs.iter().copied(), &columns, &mut flows);
        span.contains(&flows) && reduce_bits(&parity_basis, parities) == 0
    }

    /// The legs of the labelled slot classes in `slots` (packed by slot class), packed by label
    /// class.
    pub(crate) fn labelled(&self, slots: u64) -> u64 {
        self.slot_labels
            .iter()
            .enumerate()
            .filter_map(|(s, c)| Some(self.slots.count(slots, s) as u64 * self.legs.one((*c)?)))
            .sum()
    }

    /// The legs of the anonymous slot classes in `slots`.
    pub(crate) fn anonymous(&self, slots: u64) -> u64 {
        self.slot_labels
            .iter()
            .enumerate()
            .filter(|(_, c)| c.is_none())
            .map(|(s, _)| self.slots.count(slots, s) as u64 * self.slots.one(s))
            .sum()
    }

    /// The signatures equal to the stubs of a vertex in `state`; empty until all its stubs are
    /// placed.
    pub(crate) fn rules(&self, state: VertexState) -> &[u32] {
        self.vertices.rules(state)
    }

    /// The least and the most weight `k` of a vertex completing `state`.
    #[inline(always)]
    pub(crate) fn weight(&self, state: VertexState, k: usize) -> [i32; 2] {
        self.vertices.weight(state, k)
    }

    /// The state of a vertex of this degree with the legs in `slots` (packed by slot class), if
    /// a signature of the degree contains them.
    pub(crate) fn leg_state(&self, degree: u8, slots: u64) -> Option<VertexState> {
        let mut state = self.vertices.root(degree)?;
        let (mut rest, mut c) = (slots, 0);
        while rest != 0 {
            let count = self.slots.count(rest, c);
            for _ in 0..count {
                state = self.vertices.step(state, self.slots.stub[c])?;
            }
            rest -= count as u64 * self.slots.one(c);
            c += 1;
        }
        Some(state)
    }

    /// The paired labels from `from` on that fit an edge from a vertex in state `a` to a higher
    /// one in state `b`.
    #[inline(always)]
    pub(crate) fn fitting(&self, [a, b]: [VertexState; 2], from: usize) -> Fitting<'_> {
        match &self.masks {
            Some(masks) => {
                Fitting::Mask(masks[a as usize].at_u & masks[b as usize].at_v & at_least(from))
            }
            None => Fitting::Scan(self.scan([a, b], false, from)),
        }
    }

    /// The self-loop labels from `from` on that fit a vertex in state `a`.
    #[inline(always)]
    pub(crate) fn fitting_loops(&self, a: VertexState, from: usize) -> Fitting<'_> {
        match &self.masks {
            Some(masks) => Fitting::Mask(masks[a as usize].loops & at_least(from)),
            None => Fitting::Scan(self.scan([a, a], true, from)),
        }
    }

    /// The labels from `from` on that fit, without masks.
    fn scan(&self, ends: [VertexState; 2], self_loop: bool, from: usize) -> Scan<'_> {
        Scan {
            model: self,
            ends,
            self_loop,
            next: from,
        }
    }

    /// The states of the lower and the higher vertex of a pair after adding paired label `l`,
    /// which fits them.
    #[inline(always)]
    pub(crate) fn add_label(&self, [a, b]: [VertexState; 2], l: usize) -> [VertexState; 2] {
        let label = self.labels[l];
        let step = |s, stub| self.vertices.step(s, stub).expect("the label fits");
        [step(a, label.at_u), step(b, label.at_v)]
    }

    /// The state after adding self-loop label `l`, which fits it.
    #[inline(always)]
    pub(crate) fn add_loop(&self, a: VertexState, l: usize) -> VertexState {
        let [s, t] = self.loop_labels[l].stubs;
        let step = |s, stub| self.vertices.step(s, stub).expect("the label fits");
        step(step(a, s), t)
    }

    /// Edges can take only one label and vertices one rule: a topology has one diagram.
    pub(crate) fn single_diagram(&self) -> bool {
        self.labels.len() <= 1 && self.loop_labels.len() <= 1 && self.unique_rules
    }
}

/// The mask of the labels from `from` on.
fn at_least(from: usize) -> u128 {
    u128::MAX.checked_shl(from as u32).unwrap_or(0)
}

/// The prime modulus of [`Span`].
const PRIME: u64 = (1 << 31) - 1;

/// A basis of the span of some vectors modulo [`PRIME`] in row echelon form: per pivot column, the
/// row that is zero before it and one at it.
#[derive(Debug)]
struct Span(Vec<Option<Vec<u64>>>);

impl Span {
    /// The span of no vector of length `n`.
    fn new(n: usize) -> Self {
        Span(vec![None; n])
    }

    /// `x` modulo [`PRIME`] less the rows that make it zero at every pivot: zero if `x` is in the
    /// span.
    fn reduced(&self, x: &[i64]) -> Vec<u64> {
        let mut x: Vec<u64> = x
            .iter()
            .map(|&a| a.rem_euclid(PRIME as i64) as u64)
            .collect();
        // a row changes the columns from its pivot on only
        for (pivot, row) in self.0.iter().enumerate() {
            let Some(row) = row else { continue };
            if x[pivot] != 0 {
                let factor = PRIME - x[pivot];
                for (a, &b) in x[pivot..].iter_mut().zip(&row[pivot..]) {
                    *a = (*a + factor * b) % PRIME;
                }
            }
        }
        x
    }

    fn contains(&self, x: &[i64]) -> bool {
        self.reduced(x).iter().all(|&a| a == 0)
    }

    fn insert(&mut self, x: &[i64]) {
        let mut x = self.reduced(x);
        let Some(pivot) = x.iter().position(|&a| a != 0) else {
            return;
        };
        // the inverse by Fermat's little theorem
        let (mut inverse, mut base, mut exponent) = (1, x[pivot], PRIME - 2);
        while exponent > 0 {
            if exponent & 1 == 1 {
                inverse = inverse * base % PRIME;
            }
            base = base * base % PRIME;
            exponent >>= 1;
        }
        x[pivot..]
            .iter_mut()
            .for_each(|a| *a = *a * inverse % PRIME);
        self.0[pivot] = Some(x);
    }
}

/// `x` reduced by a basis modulo two, stored by leading bit: zero if `x` is in the span.
fn reduce_bits(basis: &[u128; 128], mut x: u128) -> u128 {
    while x != 0 {
        let top = basis[127 - x.leading_zeros() as usize];
        if top == 0 {
            break;
        }
        x ^= top;
    }
    x
}

#[cfg(test)]
mod tests {
    use ahash::AHashMap;

    use super::{Automaton, Direction, Legs, Model, Request, StubSet, VertexState};

    #[test]
    fn stub_sets_are_sorted_multisets() {
        let mut a = StubSet::new(4);
        for s in [5, 1, 5, 0] {
            a.insert(s);
        }
        let mut b = StubSet::new(4);
        for s in [0, 5, 1, 5] {
            b.insert(s);
        }
        assert_eq!(a.key(), b.key());
        assert_ne!(a.key(), StubSet::new(3).key());
    }

    #[test]
    fn automaton_states_are_the_sub_multisets_of_the_signatures() {
        // stub types 0 = a (undirected), 4 = e in, 5 = e out: e e~ a, a a a, a a a a twice
        let model = Model::new(
            2,
            &[vec![4, 5, 0], vec![0; 3], vec![0; 4], vec![0; 4]],
            &[],
            &Request::default(),
        )
        .unwrap();
        let v = &model.vertices;
        // degree 3: {}, {a}, {e}, {e~}, {a e}, {a e~}, {e e~}, {a e e~}, {a a}, {a a a};
        // degree 4: {}, {a}, {a a}, {a a a}, {a a a a}
        assert_eq!(v.len(), 15);
        let walk = |degree, stubs: &[u8]| {
            stubs
                .iter()
                .try_fold(v.root(degree).unwrap(), |s, &stub| v.step(s, stub))
        };
        // states are multisets: the order of the stubs does not matter
        assert_eq!(walk(3, &[4, 0, 5]), walk(3, &[5, 4, 0]));
        assert_eq!(walk(3, &[0, 4]).map(|s| v.rules(s)), Some(&[][..]));
        assert_eq!(walk(3, &[0, 5, 4]).map(|s| v.rules(s)), Some(&[0][..]));
        assert_eq!(walk(3, &[0, 0, 0]).map(|s| v.rules(s)), Some(&[1][..]));
        assert_eq!(walk(4, &[0; 4]).map(|s| v.rules(s)), Some(&[2, 3][..]));
        assert_eq!(walk(3, &[0, 0, 4]), None);
        assert_eq!(walk(3, &[4, 4]), None);
        assert_eq!(walk(4, &[4]), None);
        assert_eq!(v.root(5), None);
        assert!(!model.unique_rules);
    }

    #[test]
    fn states_bound_the_weights_of_the_signatures_containing_them() {
        // e e~ a of order one, a a a of order zero, a a a a twice with orders zero and two
        let request = Request {
            weights: vec![(vec![1, 0, 0, 2], [0, 4])],
            ..Request::default()
        };
        let signatures = [vec![4, 5, 0], vec![0; 3], vec![0; 4], vec![0; 4]];
        let model = Model::new(2, &signatures, &[], &request).unwrap();
        let walk = |degree, stubs: &[u8]| {
            let v = &model.vertices;
            let state = stubs
                .iter()
                .try_fold(v.root(degree).unwrap(), |s, &stub| v.step(s, stub));
            model.weight(state.unwrap(), 0)
        };
        assert_eq!(walk(3, &[]), [0, 1]);
        assert_eq!(walk(3, &[0]), [0, 1]);
        assert_eq!(walk(3, &[4]), [1, 1]);
        assert_eq!(walk(3, &[0, 0]), [0, 0]);
        // a full state with two rules
        assert_eq!(walk(4, &[0, 0, 0, 0]), [0, 2]);
    }

    #[test]
    fn legs_carry_only_what_vertices_balance() {
        // QED, e e~ a: stub types 0 = a, 4 = e in, 5 = e out
        let qed = |legs: &[(u8, u8)]| {
            Model::new(2, &[vec![4, 5, 0]], legs, &Request::default())
                .unwrap()
                .balanced
        };
        assert!(qed(&[(4, 1), (5, 1)]));
        assert!(qed(&[(0, 1)]));
        assert!(qed(&[]));
        // the electron charge
        assert!(!qed(&[(4, 2)]));
        assert!(!qed(&[(4, 1), (0, 1)]));
        // phi^3 and phi^2 chi^2 (stub types 0 = phi, 3 = chi): chi comes in pairs
        let scalars = |legs: &[(u8, u8)]| {
            Model::new(
                2,
                &[vec![0; 3], vec![3, 3, 0, 0]],
                legs,
                &Request::default(),
            )
            .unwrap()
            .balanced
        };
        assert!(scalars(&[(3, 2)]));
        assert!(scalars(&[(0, 1)]));
        assert!(!scalars(&[(3, 1), (0, 1)]));
    }

    #[test]
    fn state_tables_hold_exactly_the_transitions() {
        // 65 charged fields with a photon: 198 stub types, of which each state takes few
        let signatures: Vec<Vec<u8>> = (1..=65).map(|f| vec![3 * f + 1, 3 * f + 2, 0]).collect();
        let (v, transitions) = Automaton::new(&signatures, 198, &[], 0);
        let listed: AHashMap<_, _> = transitions
            .iter()
            .map(|&(state, stub, target)| ((state, stub), target))
            .collect();
        assert_eq!(listed.len(), transitions.len());
        for s in 0..v.len() as VertexState {
            for stub in 0..198 {
                assert_eq!(v.step(s, stub), listed.get(&(s, stub)).copied());
            }
        }
        // fewer than four slots per transition and one per state without any, far fewer than a
        // slot per state and stub type
        assert!(v.slots.len() < 4 * transitions.len() + v.len());
        assert!(4 * transitions.len() + v.len() < v.len() * 198 / 10);
    }

    #[test]
    fn masks_select_the_labels_that_the_transitions_admit() {
        // QCD-like: g g g, g g g g, q~ q g, c~ c g with q = field 1, c = field 2
        let signatures = [vec![0; 3], vec![0; 4], vec![4, 5, 0], vec![7, 8, 0]];
        let model = Model::new(3, &signatures, &[(0, 2)], &Request::default()).unwrap();
        assert!(model.masks.is_some());
        let scan = |ends, self_loop, from| model.scan(ends, self_loop, from).collect::<Vec<_>>();
        let n = model.vertices.len() as VertexState;
        let mut fitting = 0;
        for a in 0..n {
            for from in 0..=model.loop_labels.len() {
                let loops: Vec<_> = model.fitting_loops(a, from).collect();
                assert_eq!(loops, scan([a, a], true, from));
            }
            for b in 0..n {
                for from in 0..=model.labels.len() {
                    let labels: Vec<_> = model.fitting([a, b], from).collect();
                    assert_eq!(labels, scan([a, b], false, from));
                    fitting += labels.len();
                }
            }
        }
        assert!(fitting > 0);
    }

    #[test]
    fn many_labels_scan_the_transitions() {
        // 65 charged fields with a photon: 131 paired labels
        let signatures: Vec<Vec<u8>> = (1..=65).map(|f| vec![3 * f + 1, 3 * f + 2, 0]).collect();
        let model = Model::new(66, &signatures, &[(0, 2)], &Request::default()).unwrap();
        assert_eq!(model.labels.len(), 131);
        assert!(model.masks.is_none());
        let photon = model.leg_state(3, model.slots.one(0)).unwrap();
        // from a vertex with a photon and an incoming fermion of field 65 to a vertex with a
        // photon, only that fermion fits, flowing out of the first vertex
        let e = model.vertices.step(photon, 3 * 65 + 1).unwrap();
        let fits: Vec<_> = model.fitting([e, photon], 0).collect();
        assert_eq!(fits.len(), 1);
        let label = model.labels[fits[0]];
        assert_eq!((label.field, label.direction), (65, Direction::Forward));
        assert_eq!(model.fitting([e, photon], fits[0] + 1).count(), 0);
        let [full, partial] = model.add_label([e, photon], fits[0]);
        assert_eq!(model.rules(full), &[64]);
        assert!(model.rules(partial).is_empty());
    }

    #[test]
    fn packed_legs_count_and_enumerate_sub_multisets() {
        // classes of multiplicity 2, 1 and 3
        let legs = Legs::new(&[(0, 2), (1, 1), (3, 3)]).unwrap();
        assert_eq!(legs.total, 6);
        assert_eq!(
            (0..3).map(|c| legs.count(legs.all, c)).collect::<Vec<_>>(),
            [2, 1, 3]
        );
        let mut subsets = vec![];
        legs.subsets(legs.all, u32::MAX, &mut subsets);
        assert_eq!(subsets.len(), 3 * 2 * 4);
        legs.subsets(legs.all, 1, &mut subsets);
        assert_eq!(subsets.len(), 4);
        assert!(subsets
            .iter()
            .all(|&(s, n)| (0..3).map(|c| legs.count(s, c) as u32).sum::<u32>() == n));
    }
}
