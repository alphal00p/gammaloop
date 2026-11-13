//! Memoized reduction of closed colour networks to group invariants.

use std::{cell::RefCell, collections::HashMap};

use super::graph::{Edges, Ix, Net, NodeGraph, Obj, Rep, SmallCut, canonical, rebuild};
use super::names;
use super::poly::{C, NA, Poly};
use super::rules;

/// Why a network could not be reduced.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Failure {
    /// A fully reduced network of symmetric tensors without a name.
    UnnamedInvariant(Vec<u32>),
    /// No rule applies.
    Stuck,
    /// A cycle of rewrites returns its own network with coefficient one.
    Cycle,
    /// Nonlinear dependence on networks being reduced.
    Nonlinear,
    Depth,
    /// A monomial outgrew its fixed size (see `poly::take_overflow`).
    Overflow,
    /// An invariance system with a row 0 = nonzero: a relation is wrong.
    Inconsistent,
}

type Key = u128;

/// A value, possibly depending linearly on networks still being reduced.
#[derive(Clone, Debug, Default)]
struct Value {
    known: Poly,
    unknown: Vec<(Key, Poly)>,
}

impl Value {
    fn known(p: Poly) -> Value {
        Value { known: p, unknown: Vec::new() }
    }
    fn zero() -> Value {
        Value::default()
    }
    fn is_zero(&self) -> bool {
        self.known.is_zero() && self.unknown.is_empty()
    }
    fn add_scaled(&mut self, other: &Value, c: &Poly) {
        self.known += &other.known.mul(c);
        for (k, p) in &other.unknown {
            let p = p.mul(c);
            match self.unknown.iter_mut().find(|(j, _)| j == k) {
                Some((_, q)) => *q += &p,
                None => self.unknown.push((*k, p)),
            }
        }
        self.unknown.retain(|(_, p)| !p.is_zero());
    }
    fn scale(&self, c: &Poly) -> Value {
        let mut out = Value::zero();
        out.add_scaled(self, c);
        out
    }
    fn mul(&self, other: &Value) -> Result<Value, Failure> {
        match (self.unknown.is_empty(), other.unknown.is_empty()) {
            (true, _) => Ok(other.scale(&self.known)),
            (_, true) => Ok(self.scale(&other.known)),
            _ => Err(Failure::Nonlinear),
        }
    }
}

/// Per-thread memo budget (estimated bytes). Checked before each top-level
/// network, so no entry is dropped while a reduction may still read it: a
/// long-running caller holds at most about this much per thread plus the
/// working set of one network.
pub(crate) const MEMO_BUDGET: usize = 32 << 20;

/// IDENSO_CNR_TRACE: print the networks of an unsolvable cycle.
static TRACE: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| std::env::var_os("IDENSO_CNR_TRACE").is_some());
const DEPTH_LIMIT: usize = 400;

thread_local! {
    static MEMO: RefCell<HashMap<Key, Poly>> = RefCell::new(HashMap::new());
    static MEMO_BYTES: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    static BUDGET: std::cell::Cell<usize> = const { std::cell::Cell::new(MEMO_BUDGET) };
    static VISITS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// Drop the memo and release its table (`clear` would keep the capacity).
pub(crate) fn clear_memo() {
    MEMO.with(|m| *m.borrow_mut() = HashMap::new());
    MEMO_BYTES.with(|b| b.set(0));
}

/// Estimated bytes of one memo entry: the table slot (with hashbrown's
/// control byte and 7/8 load factor) and the terms of the polynomial.
fn entry_bytes(value: &Poly) -> usize {
    use std::mem::size_of;
    (size_of::<(Key, Poly)>() + 1) * 8 / 7 + value.0.capacity() * size_of::<(super::poly::Mono, C)>()
}

/// (entries, estimated bytes) of this thread's memo.
pub(crate) fn memo_stats() -> (usize, usize) {
    (MEMO.with(|m| m.borrow().len()), MEMO_BYTES.with(|b| b.get()))
}

/// Set this thread's memo budget (tests).
#[cfg(test)]
pub(crate) fn set_memo_budget(bytes: usize) {
    BUDGET.with(|b| b.set(bytes));
}

/// The estimated bytes of the memo, recounted entry by entry (tests).
#[cfg(test)]
pub(crate) fn recount_memo_bytes() -> usize {
    MEMO.with(|m| m.borrow().values().map(entry_bytes).sum())
}

/// Before a top-level network: drop the memo if it is over budget.
pub(crate) fn trim_memo() {
    if MEMO_BYTES.with(|b| b.get()) > BUDGET.with(|b| b.get()) {
        clear_memo();
    }
}

pub(crate) fn take_visits() -> usize {
    VISITS.with(|v| v.replace(0))
}

fn key_of(certificate: &[u32]) -> Key {
    let mut hasher = blake3::Hasher::new();
    for c in certificate {
        hasher.update(&c.to_le_bytes());
    }
    let bytes = hasher.finalize();
    u128::from_le_bytes(bytes.as_bytes()[..16].try_into().unwrap())
}

pub(crate) struct Reducer {
    stack: Vec<Key>,
    /// Equations of networks in a cycle through a network on the stack.
    pending: HashMap<Key, Value>,
    pub(crate) visited: usize,
    /// With IDENSO_CNR_TRACE: the networks by key.
    traced: HashMap<Key, Vec<Obj>>,
}

/// Replace the unknowns of `value` that have been solved (memoized) since it
/// was set up: a pending equation can outlive the cycle that solved them.
fn substitute_solved(value: &mut Value) {
    let mut i = 0;
    while i < value.unknown.len() {
        match MEMO.with(|m| m.borrow().get(&value.unknown[i].0).cloned()) {
            Some(solved) => {
                let (_, c) = value.unknown.swap_remove(i);
                value.known += &solved.mul(&c);
            }
            None => i += 1,
        }
    }
}

fn memoize(key: Key, value: &Poly) {
    // Hard bound inside one reduction (a runaway guard).
    if MEMO_BYTES.with(|b| b.get()) >= 8 * BUDGET.with(|b| b.get()) {
        clear_memo();
    }
    let value = value.clone();
    let added = entry_bytes(&value);
    let removed = MEMO.with(|m| m.borrow_mut().insert(key, value)).map_or(0, |old| entry_bytes(&old));
    MEMO_BYTES.with(|b| b.set(b.get() + added - removed));
}

/// A row 0 = rhs with a nonzero rhs that depends on no outer unknown: the
/// relations that produced it are wrong.
fn contradicts(rhs: &Value) -> bool {
    rhs.unknown.is_empty() && !rhs.known.is_zero()
}

/// Gauss-Jordan elimination of sum_j a_j x_j = rhs rows; the values of the
/// determined unknowns. A row that reduces to 0 = nonzero is an error.
fn eliminate(m: usize, rows: Vec<(Vec<(usize, C)>, Value)>) -> Result<Vec<Option<Value>>, Failure> {
    let mut rows: Vec<(HashMap<usize, C>, Value)> = rows.into_iter().map(|(r, v)| (r.into_iter().collect(), v)).collect();
    let mut pivots: Vec<Option<usize>> = vec![None; m];
    let mut done = vec![false; rows.len()];
    for col in 0..m {
        let Some(r) = (0..rows.len()).find(|&r| !done[r] && rows[r].0.get(&col).is_some_and(|x| !x.is_zero())) else {
            continue;
        };
        done[r] = true;
        let inv = rows[r].0[&col].inv();
        let pivot_row: HashMap<usize, C> = rows[r].0.iter().map(|(k, x)| (*k, x.mul(&inv))).collect();
        let pivot_rhs = rows[r].1.scale(&Poly::constant(inv));
        rows[r] = (pivot_row.clone(), pivot_rhs.clone());
        for t in 0..rows.len() {
            if t == r {
                continue;
            }
            let Some(f) = rows[t].0.get(&col).cloned() else { continue };
            if f.is_zero() {
                continue;
            }
            for (k, x) in &pivot_row {
                let e = rows[t].0.entry(*k).or_insert_with(C::zero);
                e.add_assign(&x.mul(&f).neg());
            }
            rows[t].0.retain(|_, x| !x.is_zero());
            rows[t].1.add_scaled(&pivot_rhs, &Poly::constant(f.neg()));
        }
        pivots[col] = Some(r);
    }
    // Rows without a pivot have no unknown left.
    if rows.iter().zip(&done).any(|((row, rhs), &d)| !d && row.is_empty() && contradicts(rhs)) {
        return Err(Failure::Inconsistent);
    }
    Ok((0..m)
        .map(|col| {
            let r = pivots[col]?;
            // Determined if the pivot row has no other unknown.
            (rows[r].0.len() == 1).then(|| rows[r].1.clone())
        })
        .collect())
}

/// Solve x_i = known_i + sum_j c_ij x_j for the members (numeric c).
fn solve(members: &[Key], equations: &[Value]) -> Result<Vec<Poly>, Failure> {
    let m = members.len();
    let index: HashMap<Key, usize> = members.iter().enumerate().map(|(i, k)| (*k, i)).collect();
    let mut a = vec![vec![C::zero(); m]; m];
    let mut b: Vec<Poly> = Vec::with_capacity(m);
    for (i, eq) in equations.iter().enumerate() {
        a[i][i] = C::one();
        for (k, c) in &eq.unknown {
            let j = *index.get(k).ok_or(Failure::Nonlinear)?;
            let c = c.as_constant().ok_or(Failure::Nonlinear)?;
            a[i][j].add_assign(&c.neg());
        }
        b.push(eq.known.clone());
    }
    for col in 0..m {
        let pivot = (col..m).find(|&r| !a[r][col].is_zero()).ok_or(Failure::Cycle)?;
        a.swap(col, pivot);
        b.swap(col, pivot);
        let inv = a[col][col].inv();
        for x in a[col].iter_mut() {
            *x = x.mul(&inv);
        }
        b[col] = b[col].scale(&inv);
        for r in 0..m {
            if r == col || a[r][col].is_zero() {
                continue;
            }
            let f = a[r][col].clone();
            for c in col..m {
                let t = a[col][c].mul(&f);
                a[r][c].add_assign(&t.neg());
            }
            let t = b[col].scale(&f.neg());
            b[r] += &t;
        }
    }
    Ok(b)
}

impl Reducer {
    pub(crate) fn new() -> Reducer {
        Reducer { stack: Vec::new(), pending: HashMap::new(), visited: 0, traced: HashMap::new() }
    }

    pub(crate) fn evaluate(&mut self, net: Net) -> Result<Poly, Failure> {
        let value = self.reduce(net, 0)?;
        if !value.unknown.is_empty() {
            if *TRACE {
                eprintln!("unresolved unknowns: {}", value.unknown.len());
            }
            return Err(Failure::Cycle);
        }
        Ok(value.known)
    }

    fn reduce(&mut self, net: Net, depth: usize) -> Result<Value, Failure> {
        let Some((factor, net)) = net.normalize() else {
            return Ok(Value::zero());
        };
        if net.objs.is_empty() {
            return Ok(Value::known(factor));
        }
        let graph = NodeGraph::of(&net);
        let edges = Edges::of(&graph);
        let (count, comp) = edges.components(&[]);
        if count > 1 {
            let mut value = Value::known(factor);
            for part in split_components(&net, &graph, &comp, count) {
                let v = self.reduce(part, depth + 1)?;
                value = value.mul(&v)?;
                if value.is_zero() {
                    break;
                }
            }
            return Ok(value);
        }
        let canon = canonical(&graph);
        if canon.sign == 0 {
            return Ok(Value::zero());
        }
        let key = key_of(&canon.certificate);
        let factor = factor.scale(&C::rat(canon.sign as i64, 1));
        if let Some(v) = MEMO.with(|m| m.borrow().get(&key).cloned()) {
            return Ok(Value::known(v.mul(&factor)));
        }
        if self.stack.contains(&key) || self.pending.contains_key(&key) {
            return Ok(Value { known: Poly::zero(), unknown: vec![(key, factor)] });
        }
        if depth > DEPTH_LIMIT {
            return Err(Failure::Depth);
        }
        self.visited += 1;
        VISITS.with(|v| v.set(v.get() + 1));
        let cnet = rebuild(&net, &graph, &canon.labels);
        if *TRACE {
            self.traced.insert(key, cnet.objs.clone());
        }
        self.stack.push(key);
        let result = self.expand(cnet, &canon.certificate, depth);
        self.stack.pop();
        let mut value = result?;
        substitute_solved(&mut value);
        if value.unknown.is_empty() {
            memoize(key, &value.known);
            return Ok(value.scale(&factor));
        }
        // The networks this value depends on, through unsolved equations.
        let mut members = vec![key];
        let mut seen = std::collections::HashSet::from([key]);
        let mut queue: Vec<Key> = value.unknown.iter().map(|(k, _)| *k).collect();
        let mut on_stack = false;
        while let Some(k) = queue.pop() {
            if !seen.insert(k) {
                continue;
            }
            if self.stack.contains(&k) {
                on_stack = true;
                break;
            }
            // An unknown is on the stack or pending, unless its network was
            // solved and then dropped with a full memo.
            let Some(equation) = self.pending.get_mut(&k) else {
                return Err(Failure::Cycle);
            };
            substitute_solved(equation);
            members.push(k);
            queue.extend(equation.unknown.iter().map(|(j, _)| *j));
        }
        if on_stack {
            // Part of a cycle through a network still being reduced.
            self.pending.insert(key, value);
            return Ok(Value { known: Poly::zero(), unknown: vec![(key, factor)] });
        }
        let mut equations = Vec::with_capacity(members.len());
        equations.push(value);
        for k in &members[1..] {
            equations.push(self.pending.remove(k).unwrap());
        }
        let solution = solve(&members, &equations).inspect_err(|_| {
            if *TRACE {
                eprintln!("singular cycle of {} networks", members.len());
                for (k, eq) in members.iter().zip(&equations) {
                    eprintln!("  {k:032x} {:?}", self.traced.get(k));
                    eprintln!("     = {:?} + {:?}", eq.known.0.len(), eq.unknown.iter().map(|(j, c)| (format!("{:08x}", j >> 96), c.as_constant())).collect::<Vec<_>>());
                }
            }
        })?;
        for (k, v) in members.iter().zip(&solution) {
            memoize(*k, v);
        }
        Ok(Value::known(solution[0].mul(&factor)))
    }

    fn sum(&mut self, terms: rules::Terms, depth: usize) -> Result<Value, Failure> {
        let mut value = Value::zero();
        for (c, net) in terms {
            if c.is_zero() {
                continue;
            }
            let v = self.reduce(net, depth + 1)?;
            value.add_scaled(&v, &c);
        }
        Ok(value)
    }

    fn expand(&mut self, net: Net, certificate: &[u32], depth: usize) -> Result<Value, Failure> {
        let graph = NodeGraph::of(&net);
        let edges = Edges::of(&graph);
        match edges.small_cut() {
            Some(SmallCut::Bridge) => return Ok(Value::zero()),
            Some(SmallCut::Two(e1, e2)) => return self.split_two_cut(&net, &graph, &edges, e1, e2, depth),
            None => {}
        }
        let rules: [fn(&Net) -> Option<rules::Terms>; 4] =
            [rules::symmetrize_step, rules::trace_local, rules::f_cycle, rules::trace_pairs];
        for rule in rules {
            if let Some(terms) = rule(&net) {
                return self.sum(terms, depth);
            }
        }
        if let Some(terms) = rules::symmetrize_start(&net) {
            return self.sum(terms, depth);
        }
        if let Some(terms) = rules::symmetric_local(&net) {
            return self.sum(terms, depth);
        }
        if rules::odd_structure_constants(&net) {
            return Ok(Value::zero());
        }
        if let Some(value) = rules::dddff(&net) {
            return Ok(Value::known(value));
        }
        if net.count_f() > 0 && rules::only_symmetric_or_f(&net) {
            return self.solve_invariance(net, certificate, depth);
        }
        if rules::only_symmetric(&net) {
            let entry = names::lookup(certificate).ok_or_else(|| Failure::UnnamedInvariant(certificate.to_vec()))?;
            return Ok(Value::known(names::value(entry)));
        }
        Err(Failure::Stuck)
    }

    /// Whether no rule other than an invariance move applies.
    fn is_core(net: &Net) -> bool {
        if !rules::only_symmetric_or_f(net) || net.count_f() == 0 {
            return false;
        }
        let graph = NodeGraph::of(net);
        let edges = Edges::of(&graph);
        if edges.components(&[]).0 != 1 || edges.small_cut().is_some() {
            return false;
        }
        rules::f_cycle(net).is_none()
            && rules::symmetric_local(net).is_none()
            && !rules::odd_structure_constants(net)
            && rules::dddff(net).is_none()
    }

    /// A network of symmetric tensors and structure constants to which only
    /// invariance moves apply: every invariance relation among the networks
    /// they connect (all with the same tensors), with the other networks
    /// reduced, solved as one linear system (a local Laporta reduction; it
    /// replaces color.h's choices of moves and its special identities for
    /// holding patterns).
    fn solve_invariance(&mut self, net: Net, certificate: &[u32], depth: usize) -> Result<Value, Failure> {
        const CONFIG_LIMIT: usize = 4000;
        let root = key_of(certificate);
        let mut index: HashMap<Key, usize> = HashMap::from([(root, 0)]);
        let mut configs: Vec<Net> = vec![net];
        let mut rows: Vec<(Vec<(usize, C)>, Value)> = Vec::new();
        let mut next = 0;
        while next < configs.len() {
            let config = configs[next].clone();
            next += 1;
            for relation in rules::invariance_relations(&config) {
                let mut row: Vec<(usize, C)> = Vec::new();
                let mut rhs = Value::zero();
                for (c, term) in relation {
                    let Some((factor, term)) = term.normalize() else { continue };
                    let c = c.mul(&factor);
                    if Self::is_core(&term) {
                        let (cert, sign, canon) = term.canonicalize();
                        if sign == 0 {
                            continue;
                        }
                        let key = key_of(&cert);
                        let coefficient = c.as_constant().ok_or(Failure::Nonlinear)?.mul(&C::rat(sign as i64, 1));
                        if let Some(v) = MEMO.with(|m| m.borrow().get(&key).cloned()) {
                            rhs.add_scaled(&Value::known(v), &Poly::constant(coefficient.neg()));
                            continue;
                        }
                        let j = match index.get(&key) {
                            Some(&j) => j,
                            None => {
                                if configs.len() >= CONFIG_LIMIT {
                                    return Err(Failure::Depth);
                                }
                                index.insert(key, configs.len());
                                configs.push(canon);
                                configs.len() - 1
                            }
                        };
                        match row.iter_mut().find(|(i, _)| *i == j) {
                            Some((_, x)) => x.add_assign(&coefficient),
                            None => row.push((j, coefficient)),
                        }
                    } else {
                        let v = self.reduce(term, depth + 1)?;
                        rhs.add_scaled(&v, &c.scale(&C::rat(-1, 1)));
                    }
                }
                row.retain(|(_, x)| !x.is_zero());
                if !row.is_empty() {
                    rows.push((row, rhs));
                } else if contradicts(&rhs) {
                    // Every network of the relation is known: it must hold.
                    return Err(Failure::Inconsistent);
                }
            }
        }
        let solution = eliminate(configs.len(), rows)?;
        let value = solution[0].clone().ok_or(Failure::Cycle)?;
        for (key, &j) in &index {
            if let Some(v) = &solution[j]
                && v.unknown.is_empty()
            {
                memoize(*key, &v.known);
            }
        }
        Ok(value)
    }

    /// Schur's lemma: a subnetwork with two adjoint legs is proportional to
    /// delta^{ab}, one with an R line through it to the identity.
    fn split_two_cut(
        &mut self,
        net: &Net,
        graph: &NodeGraph,
        edges: &Edges,
        e1: u32,
        e2: u32,
        depth: usize,
    ) -> Result<Value, Failure> {
        let (_, comp) = edges.components(&[e1, e2]);
        let side = |node: u32| comp[node as usize];
        let (fundamental, rep) = if edges.fundamental[e1 as usize] {
            let (o, _) = graph.origin[edges.ends[e1 as usize][0].0 as usize];
            let rep = match &net.objs[o as usize] {
                Obj::Tr(rep, _) | Obj::St(rep, _, _) => *rep,
                _ => unreachable!(),
            };
            (true, rep)
        } else {
            (false, Rep::A)
        };
        let mut parts = [Vec::new(), Vec::new()];
        if !fundamental {
            let label = |e: u32| -> (Ix, Ix) {
                // (index, side of end 0)
                let (node, slot) = edges.ends[e as usize][0];
                let (o, p) = graph.origin[node as usize];
                let idx = match &net.objs[o as usize] {
                    Obj::Tr(_, s) => s[p as usize],
                    Obj::St(_, _, s) => {
                        if graph.kind[node as usize] >= 128 { s[slot as usize - 2] } else { s[p as usize] }
                    }
                    obj => obj.indices()[slot as usize],
                };
                (idx, 0)
            };
            let (i, _) = label(e1);
            let (j, _) = label(e2);
            let mut seen = vec![false; net.objs.len()];
            for (node, &(o, _)) in graph.origin.iter().enumerate() {
                if !seen[o as usize] {
                    seen[o as usize] = true;
                    parts[side(node as u32) as usize].push(net.objs[o as usize].clone());
                }
            }
            let mut value = Value::known(Poly::mono(super::poly::Mono::var(NA, -1), C::one()));
            for objs in parts {
                let mut part = Net { objs, next: net.next };
                part.rename(j, i);
                let v = self.reduce(part, depth + 1)?;
                value = value.mul(&v)?;
                if value.is_zero() {
                    break;
                }
            }
            return Ok(value);
        }
        // Both edges on one trace: split it into its two arcs. Sorted
        // positions of an arc are a rotation of it; a partially symmetrized
        // trace keeps its block, followed by its remaining generators.
        let (o, _) = graph.origin[edges.ends[e1 as usize][0].0 as usize];
        let mut seen = vec![false; net.objs.len()];
        seen[o as usize] = true;
        for (node, &(obj, _)) in graph.origin.iter().enumerate() {
            if !seen[obj as usize] {
                seen[obj as usize] = true;
                parts[side(node as u32) as usize].push(net.objs[obj as usize].clone());
            }
        }
        let mut arcs: [Vec<(u16, Ix)>; 2] = [Vec::new(), Vec::new()];
        let mut block = None;
        let (indices, split) = match &net.objs[o as usize] {
            Obj::Tr(_, seq) => (seq.clone(), 0),
            Obj::St(_, split, seq) => (seq.clone(), *split),
            _ => unreachable!(),
        };
        for (node, &(obj, p)) in graph.origin.iter().enumerate() {
            if obj != o {
                continue;
            }
            if split > 0 && p == 0 {
                block = Some(side(node as u32) as usize);
            } else {
                arcs[side(node as u32) as usize].push((p, indices[p as usize]));
            }
        }
        let mut value = Value::known(Poly::mono(super::poly::Mono::var(rep.dimension(), -1), C::one()));
        for (k, mut objs) in parts.into_iter().enumerate() {
            arcs[k].sort_unstable();
            let ord: Vec<Ix> = arcs[k].iter().map(|x| x.1).collect();
            if block == Some(k) {
                let mut seq = indices[..split].to_vec();
                seq.extend(ord);
                objs.push(Obj::St(rep, split, seq));
            } else {
                objs.push(Obj::Tr(rep, ord));
            }
            let v = self.reduce(Net { objs, next: net.next }, depth + 1)?;
            value = value.mul(&v)?;
            if value.is_zero() {
                break;
            }
        }
        Ok(value)
    }
}

fn split_components(net: &Net, graph: &NodeGraph, comp: &[u32], count: usize) -> Vec<Net> {
    let mut parts = vec![Vec::new(); count];
    let mut seen = vec![false; net.objs.len()];
    for (node, &(o, _)) in graph.origin.iter().enumerate() {
        if !seen[o as usize] {
            seen[o as usize] = true;
            parts[comp[node] as usize].push(net.objs[o as usize].clone());
        }
    }
    parts.into_iter().map(|objs| Net { objs, next: net.next }).collect()
}

#[cfg(test)]
mod tests {
    use super::{C, Failure, Poly, Value, eliminate};

    fn row(entries: &[(usize, i64)], rhs: i64) -> (Vec<(usize, C)>, Value) {
        (entries.iter().map(|&(j, c)| (j, C::rat(c, 1))).collect(), Value::known(Poly::rat(rhs, 1)))
    }

    #[test]
    fn elimination_detects_contradictions() {
        // x0 + x1 = 3, x0 - x1 = 1, 2 x0 = 4: consistent and redundant.
        let solution = eliminate(2, vec![row(&[(0, 1), (1, 1)], 3), row(&[(0, 1), (1, -1)], 1), row(&[(0, 2)], 4)]).unwrap();
        assert_eq!(solution[0].as_ref().map(|v| v.known.clone()), Some(Poly::rat(2, 1)));
        assert_eq!(solution[1].as_ref().map(|v| v.known.clone()), Some(Poly::rat(1, 1)));
        // ... and 2 x0 = 5 instead: 0 = 1 after elimination.
        let contradiction = eliminate(2, vec![row(&[(0, 1), (1, 1)], 3), row(&[(0, 1), (1, -1)], 1), row(&[(0, 2)], 5)]);
        assert_eq!(contradiction.err(), Some(Failure::Inconsistent));
        // An underdetermined system leaves x0, x1 open and is no error.
        let open = eliminate(2, vec![row(&[(0, 1), (1, 1)], 3), row(&[(0, 2), (1, 2)], 6)]).unwrap();
        assert!(open.iter().all(Option::is_none));
    }
}
