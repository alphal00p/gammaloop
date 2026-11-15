//! Conversion between Symbolica colour atoms and closed colour networks.

use std::collections::HashMap;

use spenso::{network::library::symbolic::ETS, network::tags::SPENSO_TAG as T, shadowing};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    coefficient::CoefficientView,
};

use spenso::structure::{abstract_index::AIND_SYMBOLS, representation::RepName};

use super::graph::{Ix, Net, Obj, Rep};
use crate::representations::{ColorAdjoint, ColorFundamental};
use super::names;
use super::poly::{C, CA, CR, FIRST_NAMED, I2R, NA, NR, Poly};
use crate::color::CS;

/// The representations of a network.
#[derive(Clone, Debug, Default)]
pub(crate) struct Context {
    /// cof(N), without index.
    pub(crate) fundamental: Option<Atom>,
    /// coad(N_A), without index.
    pub(crate) adjoint: Option<Atom>,
}

#[derive(Default)]
struct Builder {
    ctx: Context,
    coefficient: C,
    adjoint: HashMap<Atom, Ix>,
    adjoint_count: Vec<u8>,
    fundamental: HashMap<Atom, Ix>,
    /// (adjoint label or None for identity, row, col)
    lines: Vec<(Option<Ix>, Ix, Ix)>,
    objs: Vec<Obj>,
    metric_pairs: Vec<(Ix, Ix)>,
    fresh_fundamental: usize,
}

fn fun(view: AtomView<'_>) -> Option<(symbolica::atom::Symbol, Vec<AtomView<'_>>)> {
    let AtomView::Fun(f) = view else { return None };
    Some((f.get_symbol(), f.iter().collect()))
}

fn number(view: AtomView<'_>) -> Option<C> {
    let AtomView::Num(n) = view else { return None };
    match n.get_coeff_view() {
        CoefficientView::Natural(num, den, inum, iden) => Some(C {
            re: symbolica::domains::rational::Rational::new(num, den),
            im: symbolica::domains::rational::Rational::new(inum, iden),
        }),
        CoefficientView::Large(r, i) => Some(C { re: r.to_rat(), im: i.to_rat() }),
        _ => None,
    }
}

fn is_marker(view: AtomView<'_>, marker: symbolica::atom::Symbol) -> bool {
    matches!(view, AtomView::Var(v) if v.get_symbol() == marker)
}

/// A fixed component cind(..) or find(..): equal components are not summed.
fn is_concrete_index(index: AtomView<'_>) -> bool {
    fun(index).is_some_and(|(s, _)| s == AIND_SYMBOLS.cind || s == AIND_SYMBOLS.find)
}

impl Builder {
    fn rep_of(&mut self, rep: AtomView<'_>, adjoint: bool) -> Option<()> {
        let slot = if adjoint { &mut self.ctx.adjoint } else { &mut self.ctx.fundamental };
        match slot {
            Some(existing) => (existing.as_view() == rep).then_some(()),
            None => {
                *slot = Some(rep.to_owned());
                Some(())
            }
        }
    }

    /// An adjoint slot coad(N, i) with an abstract (summable) index i.
    fn adjoint_slot(&mut self, slot: AtomView<'_>) -> Option<Ix> {
        let (s, args) = fun(slot)?;
        if s != CS.adjoint_rep || args.len() != 2 || is_concrete_index(args[1]) {
            return None;
        }
        let rep = ColorAdjoint {}.to_symbolic([args[0].to_owned()]);
        self.rep_of(rep.as_view(), true)?;
        let label = match self.adjoint.get(&args[1].to_owned()) {
            Some(&label) => label,
            None => {
                let label = self.adjoint_count.len() as Ix;
                self.adjoint_count.push(0);
                self.adjoint.insert(args[1].to_owned(), label);
                label
            }
        };
        self.adjoint_count[label as usize] += 1;
        Some(label)
    }

    /// A summed adjoint index of an expanded line or trace.
    fn fresh_adjoint_label(&mut self) -> Ix {
        let label = self.adjoint_count.len() as Ix;
        self.adjoint_count.push(2);
        label
    }

    fn fresh_fundamental_label(&mut self) -> Ix {
        self.fresh_fundamental += 1;
        (u16::MAX as usize - self.fresh_fundamental) as Ix
    }

    /// A fundamental slot cof(N, i) (row) or dind(cof(N, i)) (column).
    fn fundamental_slot(&mut self, slot: AtomView<'_>) -> Option<(Ix, bool)> {
        let (s, args) = fun(slot)?;
        let (inner, dual) = if s == CS.fundamental_rep {
            (slot, false)
        } else if s == AIND_SYMBOLS.dind && args.len() == 1 && fun(args[0]).is_some_and(|(s, _)| s == CS.fundamental_rep) {
            (args[0], true)
        } else {
            return None;
        };
        let (_, args) = fun(inner)?;
        if args.len() != 2 || is_concrete_index(args[1]) {
            return None;
        }
        let rep = ColorFundamental {}.to_symbolic([args[0].to_owned()]);
        self.rep_of(rep.as_view(), false)?;
        let next = self.fundamental.len() as Ix;
        Some((*self.fundamental.entry(args[1].to_owned()).or_insert(next), dual))
    }

    /// A factor with an optional numeric prefactor.
    fn split_number<'a>(&mut self, factor: AtomView<'a>) -> Option<AtomView<'a>> {
        match factor {
            AtomView::Mul(m) => {
                let mut body = None;
                for x in m.iter() {
                    if let Some(c) = number(x) {
                        self.coefficient = self.coefficient.mul(&c);
                    } else if body.replace(x).is_some() {
                        return None;
                    }
                }
                body
            }
            v => Some(v),
        }
    }

    /// An adjoint matrix factor f(.., in, .., out, ..): (slot, stored args with in/out positions).
    fn adjoint_matrix(&mut self, factor: AtomView<'_>) -> Option<(Ix, [Option<Ix>; 3], usize, usize)> {
        let factor = self.split_number(factor)?;
        let (s, args) = fun(factor)?;
        if s != CS.f || args.len() != 3 {
            return None;
        }
        let pin = args.iter().position(|a| is_marker(*a, T.chain_in))?;
        let pout = args.iter().position(|a| is_marker(*a, T.chain_out))?;
        let pa = (0..3).find(|&p| p != pin && p != pout)?;
        let a = self.adjoint_slot(args[pa])?;
        let mut stored = [None; 3];
        stored[pa] = Some(a);
        Some((a, stored, pin, pout))
    }

    fn fundamental_matrix(&mut self, factor: AtomView<'_>) -> Option<Ix> {
        let factor = self.split_number(factor)?;
        let (s, args) = fun(factor)?;
        if s != CS.t || args.len() != 3 || !is_marker(args[1], T.chain_in) || !is_marker(args[2], T.chain_out) {
            return None;
        }
        self.adjoint_slot(args[0])
    }

    fn factor(&mut self, factor: AtomView<'_>) -> Option<()> {
        if let Some(c) = number(factor) {
            self.coefficient = self.coefficient.mul(&c);
            return Some(());
        }
        if let AtomView::Pow(p) = factor {
            // An open base squared: its two copies are contracted.
            let (base, exp) = p.get_base_exp();
            if number(exp) != Some(C::rat(2, 1)) {
                return None;
            }
            self.factor(base)?;
            return self.factor(base);
        }
        let (s, args) = fun(factor)?;
        if s == CS.f && args.len() == 3 {
            let a = self.adjoint_slot(args[0])?;
            let b = self.adjoint_slot(args[1])?;
            let c = self.adjoint_slot(args[2])?;
            self.objs.push(Obj::F([a, b, c]));
            return Some(());
        }
        if s == CS.t && args.len() == 3 {
            let a = self.adjoint_slot(args[0])?;
            let (r, rd) = self.fundamental_slot(args[1])?;
            let (c, cd) = self.fundamental_slot(args[2])?;
            if rd || !cd {
                return None;
            }
            self.lines.push((Some(a), r, c));
            return Some(());
        }
        if s == ETS.metric && args.len() == 2 {
            if let (Some(x), Some(y)) = (self.adjoint_slot_peek(args[0]), self.adjoint_slot_peek(args[1])) {
                let _ = (x, y);
                let x = self.adjoint_slot(args[0])?;
                let y = self.adjoint_slot(args[1])?;
                self.metric_pairs.push((x, y));
                return Some(());
            }
            let (x, xd) = self.fundamental_slot(args[0])?;
            let (y, yd) = self.fundamental_slot(args[1])?;
            match (xd, yd) {
                (false, true) => self.lines.push((None, x, y)),
                (true, false) => self.lines.push((None, y, x)),
                _ => return None,
            }
            return Some(());
        }
        if s == T.chain && args.len() >= 2 {
            let (start, end, factors) = (args[0], args[1], &args[2..]);
            if self.adjoint_slot_peek(start).is_some() {
                let x = self.adjoint_slot(start)?;
                let y = self.adjoint_slot(end)?;
                if factors.is_empty() {
                    self.metric_pairs.push((x, y));
                    return Some(());
                }
                let mut prev = x;
                for (t, factor) in factors.iter().enumerate() {
                    let (_, mut stored, pin, pout) = self.adjoint_matrix(*factor)?;
                    let next = if t + 1 == factors.len() { y } else { self.fresh_adjoint_label() };
                    stored[pin] = Some(prev);
                    stored[pout] = Some(next);
                    self.objs.push(Obj::F(stored.map(|s| s.unwrap())));
                    prev = next;
                }
                return Some(());
            }
            let (r, rd) = self.fundamental_slot(start)?;
            let (c, cd) = self.fundamental_slot(end)?;
            if rd || !cd {
                return None;
            }
            if factors.is_empty() {
                self.lines.push((None, r, c));
                return Some(());
            }
            let mut prev = r;
            for (t, factor) in factors.iter().enumerate() {
                let a = self.fundamental_matrix(*factor)?;
                let next = if t + 1 == factors.len() { c } else { self.fresh_fundamental_label() };
                self.lines.push((Some(a), prev, next));
                prev = next;
            }
            return Some(());
        }
        if s == T.trace {
            let AtomView::Fun(f) = factor else { return None };
            let (rep, factors) = shadowing::trace_parts(f)?;
            let (rs, rargs) = fun(rep)?;
            let adjoint = if rs == CS.adjoint_rep {
                true
            } else if rs == CS.fundamental_rep {
                false
            } else {
                return None;
            };
            if rargs.len() != 1 {
                return None;
            }
            let rep_atom = rep.to_owned();
            self.rep_of(rep_atom.as_view(), adjoint)?;
            let symmetric = matches!(factors.as_slice(), [AtomView::Fun(p)] if p.get_symbol() == *shadowing::SYM);
            let factors: Vec<AtomView<'_>> = if symmetric {
                let AtomView::Fun(p) = factors[0] else { unreachable!() };
                p.iter().collect()
            } else {
                factors
            };
            if factors.is_empty() {
                self.objs.push(Obj::Tr(if adjoint { Rep::A } else { Rep::R }, Vec::new()));
                return Some(());
            }
            if adjoint {
                if symmetric {
                    // Str of f matrices: prod(sigma) i^n D_A.
                    let mut labels = Vec::new();
                    let mut phase = C::one();
                    for factor in &factors {
                        let (a, _, pin, pout) = self.adjoint_matrix(*factor)?;
                        let pa = 3 - pin - pout;
                        // f(a,in,out) orientation
                        let even = (pa + 1) % 3 == pin;
                        phase = phase.mul(&if even { C::i() } else { C::imag(-1, 1) });
                        labels.push(a);
                    }
                    self.coefficient = self.coefficient.mul(&phase);
                    self.objs.push(Obj::D(Rep::A, labels));
                    return Some(());
                }
                let n = factors.len();
                let ks: Vec<Ix> = (0..n).map(|_| self.fresh_adjoint_label()).collect();
                for (t, factor) in factors.iter().enumerate() {
                    let (_, mut stored, pin, pout) = self.adjoint_matrix(*factor)?;
                    stored[pin] = Some(ks[(t + n - 1) % n]);
                    stored[pout] = Some(ks[t]);
                    self.objs.push(Obj::F(stored.map(|s| s.unwrap())));
                }
                return Some(());
            }
            let labels = factors.iter().map(|x| self.fundamental_matrix(*x)).collect::<Option<Vec<_>>>()?;
            self.objs.push(if symmetric { Obj::D(Rep::R, labels) } else { Obj::Tr(Rep::R, labels) });
            return Some(());
        }
        None
    }

    fn adjoint_slot_peek(&self, slot: AtomView<'_>) -> Option<()> {
        let (s, args) = fun(slot)?;
        (s == CS.adjoint_rep && args.len() == 2).then_some(())
    }

    fn finish(mut self) -> Option<(Net, C, Context)> {
        // Every adjoint label twice, after metric joins.
        let n = self.adjoint_count.len();
        let mut parent: Vec<Ix> = (0..n as Ix).collect();
        fn find(p: &mut [Ix], x: Ix) -> Ix {
            let mut r = x;
            while p[r as usize] != r {
                r = p[r as usize];
            }
            let mut y = x;
            while p[y as usize] != r {
                let next = p[y as usize];
                p[y as usize] = r;
                y = next;
            }
            r
        }
        // Every label occurs twice, counting metric slots; metrics then join
        // labels, and a class held by metrics only is a loop, delta^{aa} = N_A.
        if self.adjoint_count.iter().any(|&c| c != 2) {
            return None;
        }
        let mut held: Vec<i32> = self.adjoint_count.iter().map(|&c| c as i32).collect();
        for &(x, y) in &self.metric_pairs {
            held[x as usize] -= 1;
            held[y as usize] -= 1;
        }
        for &(x, y) in &self.metric_pairs {
            let (rx, ry) = (find(&mut parent, x), find(&mut parent, y));
            if rx != ry {
                parent[rx as usize] = ry;
            }
        }
        let mut class_held = vec![0i32; n];
        for x in 0..n as Ix {
            let r = find(&mut parent, x);
            class_held[r as usize] += held[x as usize];
        }
        for x in 0..n as Ix {
            if find(&mut parent, x) == x && class_held[x as usize] == 0 {
                self.objs.push(Obj::Tr(Rep::A, Vec::new()));
            }
        }
        // Fundamental lines into traces.
        let mut by_row: HashMap<Ix, usize> = HashMap::new();
        for (k, &(_, r, _)) in self.lines.iter().enumerate() {
            if by_row.insert(r, k).is_some() {
                return None;
            }
        }
        let mut used = vec![false; self.lines.len()];
        for start in 0..self.lines.len() {
            if used[start] {
                continue;
            }
            let mut seq = Vec::new();
            let mut k = start;
            loop {
                if used[k] {
                    if k != start {
                        return None;
                    }
                    break;
                }
                used[k] = true;
                let (a, _, c) = self.lines[k];
                if let Some(a) = a {
                    seq.push(a);
                }
                k = *by_row.get(&c)?;
            }
            self.objs.push(Obj::Tr(Rep::R, seq));
        }
        let mut objs = self.objs;
        for obj in objs.iter_mut() {
            match obj {
                Obj::F(i) => i.iter_mut().for_each(|x| *x = find(&mut parent, *x)),
                Obj::Tr(_, i) | Obj::D(_, i) | Obj::St(_, _, i) => i.iter_mut().for_each(|x| *x = find(&mut parent, *x)),
            }
        }
        // Dense labels.
        let mut dense: HashMap<Ix, Ix> = HashMap::new();
        for obj in objs.iter_mut() {
            let idx: &mut [Ix] = match obj {
                Obj::F(i) => i,
                Obj::Tr(_, i) | Obj::D(_, i) | Obj::St(_, _, i) => i,
            };
            for x in idx.iter_mut() {
                let next = dense.len() as Ix;
                *x = *dense.entry(*x).or_insert(next);
            }
        }
        Some((Net::new(objs), self.coefficient, self.ctx))
    }
}

/// Parse a product of colour factors into a closed network.
pub(crate) fn parse(factors: &[AtomView<'_>]) -> Option<(Net, C, Context)> {
    let mut b = Builder { coefficient: C::one(), ..Default::default() };
    for f in factors {
        b.factor(*f)?;
    }
    b.finish()
}

fn c_atom(c: &C) -> Atom {
    let re = Atom::num(c.re.clone());
    if c.im.is_zero() {
        return re;
    }
    re + Atom::num(c.im.clone()) * Atom::i()
}

fn representation(rep: Rep, ctx: &Context) -> Atom {
    match rep {
        Rep::R => ctx.fundamental.clone().unwrap_or_else(|| ColorFundamental {}.to_symbolic([Atom::var(CS.nc)])),
        Rep::A => ctx.adjoint.clone().unwrap_or_else(|| ColorAdjoint {}.to_symbolic([Atom::var(symbolica::symbol!("spenso::dA"))])),
    }
}

fn rep_dimension(rep: &Atom) -> Atom {
    let AtomView::Fun(f) = rep.as_view() else { unreachable!() };
    f.iter().next().unwrap().to_owned()
}

/// The atom of a named invariant.
pub(crate) fn invariant_atom(inv: &names::Invariant, ctx: &Context) -> Atom {
    let reps: Vec<Atom> = inv.reps.iter().map(|r| representation(*r, ctx)).collect();
    if inv.reps.len() == 2 && inv.name.len() == 3 {
        let degree = (inv.name.as_bytes()[1] - b'0') as i64;
        // gram(k, A, R): adjoint first.
        let (l, r) = if inv.reps[0] == Rep::R && inv.reps[1] == Rep::A { (reps[1].clone(), reps[0].clone()) } else { (reps[0].clone(), reps[1].clone()) };
        return CS.gram(Atom::num(degree), l, r);
    }
    let symbol = super::invariant_symbol(inv.name);
    let mut f = symbolica::atom::FunctionBuilder::new(symbol);
    for r in &reps {
        f = f.add_arg(r);
    }
    f.finish()
}

fn variable_atom(var: u16, ctx: &Context) -> Atom {
    let fundamental = representation(Rep::R, ctx);
    let adjoint = representation(Rep::A, ctx);
    match var {
        NA => rep_dimension(&adjoint),
        NR => rep_dimension(&fundamental),
        I2R => CS.idx(Atom::num(2), fundamental),
        CR => CS.cas(Atom::num(2), fundamental),
        CA => CS.cas(Atom::num(2), adjoint),
        v if v >= FIRST_NAMED => invariant_atom(names::variable_of(v).expect("registered invariant"), ctx),
        _ => unreachable!(),
    }
}

pub(crate) fn to_atom(poly: &Poly, ctx: &Context) -> Atom {
    let mut cache: HashMap<u16, Atom> = HashMap::new();
    let terms: Vec<Atom> = poly
        .0
        .iter()
        .map(|(m, c)| {
            let mut t = c_atom(c);
            for (v, e) in m.iter() {
                let base = cache.entry(v).or_insert_with(|| variable_atom(v, ctx)).clone();
                t = t * base.pow(Atom::num(e as i64));
            }
            t
        })
        .collect();
    Atom::add_many(terms)
}
