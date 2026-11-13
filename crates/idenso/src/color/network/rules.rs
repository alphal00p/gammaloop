//! Exact rewrite rules of closed colour networks (after color.h, van Ritbergen,
//! Schellekens and Vermaseren, hep-ph/9802376). Each rule returns a linear
//! combination of networks equal to its input, or `None` if it does not apply.

use super::graph::{Ix, Net, Obj, Rep};
use super::poly::{C, CA, CR, I2R, Mono, NA, Poly};

pub(crate) type Terms = Vec<(Poly, Net)>;

fn rep_casimir(rep: Rep) -> Poly {
    Poly::var(rep.casimir())
}

/// Sign of the permutation taking the stored order `f` to `(x, y, z)`.
fn orientation(f: &[Ix; 3], x: Ix, y: Ix) -> i64 {
    let px = f.iter().position(|&i| i == x).unwrap();
    let py = f.iter().position(|&i| i == y).unwrap();
    if (px + 1) % 3 == py { 1 } else { -1 }
}

fn third(f: &[Ix; 3], x: Ix, y: Ix) -> Ix {
    *f.iter().find(|&&i| i != x && i != y).unwrap()
}

fn with_objects(net: &Net, remove: &[usize], add: Vec<Obj>) -> Net {
    let mut objs: Vec<Obj> = net
        .objs
        .iter()
        .enumerate()
        .filter(|(k, _)| !remove.contains(k))
        .map(|(_, o)| o.clone())
        .collect();
    objs.extend(add);
    Net { objs, next: net.next }
}

fn rotate(seq: &[Ix], start: usize) -> Vec<Ix> {
    let n = seq.len();
    (0..n).map(|k| seq[(start + k) % n]).collect()
}

/// Index -> F objects containing it.
fn f_holder(net: &Net) -> Vec<Option<usize>> {
    let mut h = vec![None; net.next as usize];
    for (k, o) in net.objs.iter().enumerate() {
        if let Obj::F(f) = o {
            for &i in f {
                h[i as usize] = Some(k);
            }
        }
    }
    h
}

/// Generator pairs of one trace at cyclic distance one or two.
pub(crate) fn trace_local(net: &Net) -> Option<Terms> {
    let holder = f_holder(net);
    for (k, obj) in net.objs.iter().enumerate() {
        let Obj::Tr(rep, seq) = obj else { continue };
        let n = seq.len();
        // T^a T^a = c_X.
        for p in 0..n {
            if seq[p] == seq[(p + 1) % n] {
                let rest: Vec<Ix> = rotate(seq, p + 2).into_iter().take(n - 2).collect();
                return Some(vec![(rep_casimir(*rep), with_objects(net, &[k], vec![Obj::Tr(*rep, rest)]))]);
            }
        }
        // T^a T^b f^{abc} = i C_A / 2 T^c.
        for p in 0..n {
            let (a, b) = (seq[p], seq[(p + 1) % n]);
            if let (Some(fa), Some(fb)) = (holder[a as usize], holder[b as usize])
                && fa == fb
            {
                let Obj::F(f) = &net.objs[fa] else { unreachable!() };
                let c = third(f, a, b);
                let sign = orientation(f, a, b);
                let mut rest = vec![c];
                rest.extend(rotate(seq, p + 2).into_iter().take(n - 2));
                let coef = Poly::mono(Mono::var(CA, 1), C::imag(sign, 2));
                return Some(vec![(coef, with_objects(net, &[k, fa], vec![Obj::Tr(*rep, rest)]))]);
            }
        }
        // T^a T^b T^a = (c_X - C_A/2) T^b.
        if n >= 3 {
            for p in 0..n {
                if seq[p] == seq[(p + 2) % n] && seq[p] != seq[(p + 1) % n] {
                    let mut rest = vec![seq[(p + 1) % n]];
                    rest.extend(rotate(seq, p + 3).into_iter().take(n - 3));
                    let coef = rep_casimir(*rep) + &Poly::mono(Mono::var(CA, 1), C::rat(-1, 2));
                    return Some(vec![(coef, with_objects(net, &[k], vec![Obj::Tr(*rep, rest)]))]);
                }
            }
        }
    }
    None
}

/// Contracted generator pairs of a trace at larger distance: a repeated
/// index or both legs of one structure constant. The pair with the
/// smallest distance moves together by commutation.
pub(crate) fn trace_pairs(net: &Net) -> Option<Terms> {
    let holder = f_holder(net);
    // (distance, kind (0 delta, 1 f), trace, start position, trace length)
    let mut best: Option<(usize, u8, usize, usize, usize)> = None;
    for (k, obj) in net.objs.iter().enumerate() {
        let Obj::Tr(_, seq) = obj else { continue };
        let n = seq.len();
        for p in 0..n {
            for q in p + 1..n {
                let kind = if seq[p] == seq[q] {
                    0
                } else if let (Some(fa), Some(fb)) = (holder[seq[p] as usize], holder[seq[q] as usize])
                    && fa == fb
                {
                    1
                } else {
                    continue;
                };
                let forward = q - p;
                let (d, start) = if forward <= n - forward { (forward, p) } else { (n - forward, q) };
                let candidate = (d, kind, k, start, n);
                if best.is_none_or(|b| (candidate.0, candidate.1, candidate.4, candidate.2, candidate.3) < (b.0, b.1, b.4, b.2, b.3)) {
                    best = Some(candidate);
                }
            }
        }
    }
    let (d, kind, k, start, _) = best?;
    let Obj::Tr(rep, seq) = &net.objs[k] else { unreachable!() };
    let rep = *rep;
    let s = rotate(seq, start);
    let n = s.len();
    let m = d - 1;
    let w: Vec<Ix> = s[1..=m].to_vec();
    let z: Vec<Ix> = s[d + 1..].to_vec();
    let mut terms = Vec::new();
    let mut base = net.clone();
    if kind == 0 {
        let a = s[0];
        // c_X Tr(W Z)
        let mut rest = w.clone();
        rest.extend(z.iter().copied());
        terms.push((rep_casimir(rep), with_objects(&base, &[k], vec![Obj::Tr(rep, rest)])));
        // i f(w_j, a, k) Tr(a, W[w_j -> k], Z)
        for j in 0..m {
            let kk = base.fresh();
            let mut seq2 = vec![a];
            seq2.extend(w.iter().enumerate().map(|(t, &x)| if t == j { kk } else { x }));
            seq2.extend(z.iter().copied());
            let net2 = with_objects(&base, &[k], vec![Obj::Tr(rep, seq2), Obj::F([w[j], a, kk])]);
            terms.push((Poly::constant(C::i()), Net { next: base.next, ..net2 }));
        }
    } else {
        let (x, y) = (s[0], s[d]);
        let fk = holder[x as usize].unwrap();
        let Obj::F(f) = &net.objs[fk] else { unreachable!() };
        let c = third(f, x, y);
        let sigma = orientation(f, x, y);
        // sigma i C_A/2 Tr(c, W, Z)
        let mut rest = vec![c];
        rest.extend(w.iter().copied());
        rest.extend(z.iter().copied());
        terms.push((
            Poly::mono(Mono::var(CA, 1), C::imag(sigma, 2)),
            with_objects(&base, &[k, fk], vec![Obj::Tr(rep, rest)]),
        ));
        // sigma i f(x,y,c) f(w_j, y, k) Tr(x, W[w_j -> k], Z)
        for j in 0..m {
            let kk = base.fresh();
            let mut seq2 = vec![x];
            seq2.extend(w.iter().enumerate().map(|(t, &v)| if t == j { kk } else { v }));
            seq2.extend(z.iter().copied());
            let net2 = with_objects(
                &base,
                &[k, fk],
                vec![Obj::Tr(rep, seq2), Obj::F([x, y, c]), Obj::F([w[j], y, kk])],
            );
            terms.push((Poly::constant(C::imag(sigma, 1)), Net { next: base.next, ..net2 }));
        }
    }
    let next = base.next;
    for t in terms.iter_mut() {
        t.1.next = next;
    }
    debug_assert!(n >= 2);
    Some(terms)
}

/// The shortest cycle of structure constants becomes an adjoint trace:
/// prod f(i_m, k_m, k_{m+1}) = i^n Tr(F^{i_1} ... F^{i_n}).
pub(crate) fn f_cycle(net: &Net) -> Option<Terms> {
    let fs: Vec<usize> = net.objs.iter().enumerate().filter(|(_, o)| matches!(o, Obj::F(_))).map(|(k, _)| k).collect();
    if fs.len() < 2 {
        return None;
    }
    // F-graph edges: (index, a, b)
    let holder = f_holder(net);
    let occ = net.occurrences();
    let mut edges = Vec::new();
    for (i, pair) in occ.iter().enumerate() {
        if pair[0].0 == u16::MAX {
            continue;
        }
        let (a, b) = (pair[0].0 as usize, pair[1].0 as usize);
        if holder[i] .is_some() && matches!(net.objs[a], Obj::F(_)) && matches!(net.objs[b], Obj::F(_)) && a != b {
            edges.push((i as Ix, a, b));
        }
    }
    let mut best: Option<Vec<(usize, Ix)>> = None; // cycle: (F object, index to next)
    for (e, &(idx, a, b)) in edges.iter().enumerate() {
        // shortest path b -> a avoiding edge e
        let mut prev: std::collections::HashMap<usize, (usize, Ix)> = Default::default();
        let mut queue = std::collections::VecDeque::from([b]);
        let mut seen = std::collections::HashSet::from([b]);
        let mut found = false;
        while let Some(v) = queue.pop_front() {
            if v == a {
                found = true;
                break;
            }
            for (e2, &(idx2, x, y)) in edges.iter().enumerate() {
                if e2 == e {
                    continue;
                }
                let w = if x == v { y } else if y == v { x } else { continue };
                if seen.insert(w) {
                    prev.insert(w, (v, idx2));
                    queue.push_back(w);
                }
            }
        }
        if !found {
            continue;
        }
        // path a <- ... <- b, then edge a -> b via idx
        let mut cycle = vec![(a, idx)];
        let mut v = a;
        while v != b {
            let (u, via) = prev[&v];
            cycle.push((u, via));
            v = u;
        }
        // cycle: a -(idx)-> b? We built [(a, idx), (p1, via a<-p1), ...]; reorder so each entry
        // holds the index leading to the next entry.
        // Entries: a links to b by idx; prev chain links b..a reversed.
        // Build forward order: a, b, ..., (via prev from b towards a reversed)
        let len = cycle.len();
        if best.as_ref().is_none_or(|bc| len < bc.len()) {
            // Construct forward cycle explicitly.
            let mut path = vec![a];
            let mut v = a;
            let mut links = Vec::new();
            while v != b {
                let (u, via) = prev[&v];
                links.push(via);
                path.push(u);
                v = u;
            }
            // path: a, ..., b with links between consecutive; closing link b -> a is idx
            let mut forward = Vec::new();
            for (t, &node) in path.iter().enumerate() {
                let to_next = if t + 1 < path.len() { links[t] } else { idx };
                forward.push((node, to_next));
            }
            best = Some(forward);
        }
        if best.as_ref().is_some_and(|c| c.len() == 2) {
            break;
        }
    }
    let cycle = best?;
    let len = cycle.len();
    let mut legs = Vec::with_capacity(len);
    let mut sign = 1i64;
    for t in 0..len {
        let (node, out) = cycle[t];
        let incoming = cycle[(t + len - 1) % len].1;
        let Obj::F(f) = &net.objs[node] else { unreachable!() };
        let leg = third(f, incoming, out);
        // want f(leg, in, out)
        sign *= orientation(f, leg, incoming);
        legs.push(leg);
    }
    // i^len
    let phase = match len % 4 {
        0 => C::rat(sign, 1),
        1 => C::imag(sign, 1),
        2 => C::rat(-sign, 1),
        _ => C::imag(-sign, 1),
    };
    let remove: Vec<usize> = cycle.iter().map(|c| c.0).collect();
    Some(vec![(Poly::constant(phase), with_objects(net, &remove, vec![Obj::Tr(Rep::A, legs)]))])
}

/// B_x / x! with B_1 = +1/2: the weight of one ordered selection of x
/// symmetrized generators (checked against explicit SU(3) generators).
/// `None` exactly when the weight is zero (odd x > 1). Exact for every x,
/// from sum_{k=0}^{n} (B_k/k!) / (n+1-k)! = 0 (n >= 1, B_1 = -1/2).
fn bernoulli(x: usize) -> Option<C> {
    use symbolica::domains::rational::Rational;
    if x == 1 {
        return Some(C::rat(1, 2));
    }
    if x % 2 == 1 {
        return None;
    }
    let mut inv_factorial = vec![Rational::one()];
    for j in 1..=x + 1 {
        let next = &inv_factorial[j - 1] * &Rational::new(1, j as i64);
        inv_factorial.push(next);
    }
    let mut a = vec![Rational::one()];
    for n in 1..=x {
        let mut sum = Rational::zero();
        for (k, ak) in a.iter().enumerate() {
            sum += &(ak * &inv_factorial[n + 1 - k]);
        }
        a.push(-sum);
    }
    Some(C { re: a[x].clone(), im: Rational::zero() })
}

/// A trace without contracted pairs starts its symmetrization (color.h):
/// Tr(T^{a1} T^{a2} ...) = Tr(S(T^{a1}) T^{a2} ...). The shortest goes first.
pub(crate) fn symmetrize_start(net: &Net) -> Option<Terms> {
    let (k, rep, seq) = net
        .objs
        .iter()
        .enumerate()
        .filter_map(|(k, o)| match o {
            Obj::Tr(rep, seq) => Some((k, *rep, seq.clone())),
            _ => None,
        })
        .min_by_key(|(_, _, s)| s.len())?;
    Some(vec![(Poly::one(), with_objects(net, &[k], vec![Obj::St(rep, 1, seq)]))])
}

/// One step of the symmetrization: S(A) T^b = S(A, b) + sum_x B_x E(x, A, b),
/// E(x, A, b) = sum_{a in A} i f(a, b, k) E(x-1, A - a, k), E(0, A, k) = S(A, k);
/// under the trace S(A) T^b = S(A, b) exactly.
pub(crate) fn symmetrize_step(net: &Net) -> Option<Terms> {
    let (k, rep, split, idx) = net.objs.iter().enumerate().find_map(|(k, o)| match o {
        Obj::St(rep, split, idx) => Some((k, *rep, *split, idx.clone())),
        _ => None,
    })?;
    let set = &idx[..split];
    let b = idx[split];
    let ord = &idx[split + 1..];
    let mut grown = set.to_vec();
    grown.push(b);
    grown.extend_from_slice(ord);
    let mut terms = vec![(Poly::one(), with_objects(net, &[k], vec![Obj::St(rep, split + 1, grown)]))];
    if ord.is_empty() {
        return Some(terms);
    }
    let mut base = net.clone();
    let m = set.len();
    for x in 1..=m {
        let Some(bernoulli_weight) = bernoulli(x) else { continue };
        let weight = bernoulli_weight.mul(&match x % 4 {
            0 => C::one(),
            1 => C::i(),
            2 => C::rat(-1, 1),
            _ => C::imag(-1, 1),
        });
        let mut selections = Vec::new();
        let mut chosen = Vec::with_capacity(x);
        let mut used = vec![false; m];
        fn rec(depth: usize, x: usize, m: usize, used: &mut Vec<bool>, chosen: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
            if depth == x {
                out.push(chosen.clone());
                return;
            }
            for p in 0..m {
                if !used[p] {
                    used[p] = true;
                    chosen.push(p);
                    rec(depth + 1, x, m, used, chosen, out);
                    chosen.pop();
                    used[p] = false;
                }
            }
        }
        rec(0, x, m, &mut used, &mut chosen, &mut selections);
        for sel in selections {
            let mut add = Vec::with_capacity(x + 1);
            let mut prev = b;
            for &p in &sel {
                let kk = base.fresh();
                add.push(Obj::F([set[p], prev, kk]));
                prev = kk;
            }
            let mut rest: Vec<Ix> = (0..m).filter(|p| !sel.contains(p)).map(|p| set[p]).collect();
            rest.push(prev);
            let split2 = rest.len();
            rest.extend_from_slice(ord);
            add.push(Obj::St(rep, split2, rest));
            terms.push((Poly::constant(weight.clone()), with_objects(net, &[k], add)));
        }
    }
    for t in terms.iter_mut() {
        t.1.next = base.next;
    }
    Some(terms)
}

/// Rules for networks of symmetric tensors and structure constants.
pub(crate) fn symmetric_local(net: &Net) -> Option<Terms> {
    let occ = net.occurrences();
    for (k, obj) in net.objs.iter().enumerate() {
        match obj {
            Obj::F(f) => {
                // Two legs on one symmetric tensor.
                for (x, y) in [(f[0], f[1]), (f[1], f[2]), (f[0], f[2])] {
                    let ox = other(&occ, x, k as u16);
                    let oy = other(&occ, y, k as u16);
                    if ox == oy && matches!(net.objs[ox as usize], Obj::D(..)) {
                        return Some(Vec::new());
                    }
                }
            }
            Obj::D(rep, idx) => {
                // Self contraction.
                for p in 0..idx.len() {
                    for q in p + 1..idx.len() {
                        if idx[p] == idx[q] {
                            let rest: Vec<Ix> = idx.iter().enumerate().filter(|(t, _)| *t != p && *t != q).map(|(_, &i)| i).collect();
                            return contract(net, k, *rep, rest);
                        }
                    }
                }
            }
            _ => {}
        }
    }
    // d_R^{abx} d_A^{abyz} = C_A^2/6 d_R^{xyz}.
    for (k, obj) in net.objs.iter().enumerate() {
        let Obj::D(Rep::R, small) = obj else { continue };
        if small.len() != 3 {
            continue;
        }
        for (l, other_obj) in net.objs.iter().enumerate() {
            let Obj::D(Rep::A, big) = other_obj else { continue };
            if big.len() != 4 {
                continue;
            }
            let common: Vec<Ix> = small.iter().copied().filter(|i| big.contains(i)).collect();
            if common.len() != 2 {
                continue;
            }
            let x = *small.iter().find(|i| !common.contains(i)).unwrap();
            let yz: Vec<Ix> = big.iter().copied().filter(|i| !common.contains(i)).collect();
            let coef = Poly::mono(Mono::var(CA, 2), C::rat(1, 6));
            return Some(vec![(coef, with_objects(net, &[k, l], vec![Obj::D(Rep::R, vec![x, yz[0], yz[1]])]))]);
        }
    }
    // d(x,y,z) f(x,a,k) f(y,b,k) = C_A/2 d(a,b,z).
    for (k, obj) in net.objs.iter().enumerate() {
        let Obj::D(rep, idx) = obj else { continue };
        if idx.len() != 3 {
            continue;
        }
        for p in 0..3 {
            for q in 0..3 {
                if p == q {
                    continue;
                }
                let (x, y) = (idx[p], idx[q]);
                let z = idx[3 - p - q];
                let ox = other(&occ, x, k as u16) as usize;
                let oy = other(&occ, y, k as u16) as usize;
                if ox == oy {
                    continue;
                }
                let (Obj::F(f1), Obj::F(f2)) = (&net.objs[ox], &net.objs[oy]) else { continue };
                let shared: Vec<Ix> = f1.iter().copied().filter(|i| *i != x && f2.contains(i)).collect();
                if shared.len() != 1 {
                    continue;
                }
                let kk = shared[0];
                let a = third(f1, x, kk);
                let b = third(f2, y, kk);
                let s = orientation(f1, x, a) * orientation(f2, y, b);
                let coef = Poly::mono(Mono::var(CA, 1), C::rat(s, 2));
                return Some(vec![(coef, with_objects(net, &[k, ox, oy], vec![Obj::D(*rep, vec![a, b, z])]))]);
            }
        }
    }
    None
}

fn other(occ: &[[(u16, u16); 2]], index: Ix, holder: u16) -> u16 {
    let pair = occ[index as usize];
    if pair[0].0 == holder && pair[1].0 != holder { pair[1].0 } else if pair[1].0 == holder && pair[0].0 != holder { pair[0].0 } else { holder }
}

/// d_X with one self-contracted pair (color.h `contract`).
fn contract(net: &Net, k: usize, rep: Rep, rest: Vec<Ix>) -> Option<Terms> {
    let n = rest.len();
    let base = with_objects(net, &[k], Vec::new());
    let cr = Poly::var(CR);
    let ca = Poly::var(CA);
    let lin = |a: i64, b: i64| -> Poly {
        // cR - a/b cA  (rep R) or the adjoint numbers below
        cr.clone() + &ca.scale(&C::rat(-a, b))
    };
    let add = |objs: Vec<Obj>| -> Net {
        let mut o = base.objs.clone();
        o.extend(objs);
        Net { objs: o, next: base.next }
    };
    // delta^{ab}: join the two legs, or N_A for a closed pair.
    let delta = |coef: Poly| -> Terms {
        let mut net = base.clone();
        if rest[0] == rest[1] {
            vec![(coef.mul(&Poly::var(NA)), net)]
        } else {
            net.rename(rest[1], rest[0]);
            vec![(coef, net)]
        }
    };
    match (rep, n) {
        (_, 0) => Some(vec![(Poly::var(rep.casimir()).mul(&Poly::var(rep.dimension())), base.clone())]),
        (_, 1) => Some(Vec::new()),
        (Rep::R, 2) => Some(delta(lin(1, 6).mul(&Poly::var(I2R)))),
        (Rep::R, 3) => Some(vec![(lin(1, 4), add(vec![Obj::D(Rep::R, rest)]))]),
        (Rep::R, 4) => Some(vec![
            (lin(1, 3), add(vec![Obj::D(Rep::R, rest.clone())])),
            (Poly::var(I2R).scale(&C::rat(1, 30)), add(vec![Obj::D(Rep::A, rest)])),
        ]),
        (Rep::R, 5) => {
            let mut terms = vec![(lin(5, 12), add(vec![Obj::D(Rep::R, rest.clone())]))];
            let mut b = base.clone();
            for x in 0..5 {
                for y in x + 1..5 {
                    let kk = b.fresh();
                    let others: Vec<Ix> = (0..5).filter(|t| *t != x && *t != y).map(|t| rest[t]).collect();
                    let mut o = base.objs.clone();
                    o.push(Obj::D(Rep::R, vec![rest[x], rest[y], kk]));
                    let mut da = vec![kk];
                    da.extend(others);
                    o.push(Obj::D(Rep::A, da));
                    // A pure number: both sides are homogeneous of degree 7
                    // in the generators (color.h's contract has I2R/120).
                    terms.push((Poly::rat(1, 120), Net { objs: o, next: 0 }));
                }
            }
            for t in terms.iter_mut() {
                t.1.next = b.next;
            }
            Some(terms)
        }
        (Rep::A, 2) => Some(delta(ca.mul(&ca).scale(&C::rat(5, 6)))),
        (Rep::A, 4) => Some(vec![(ca.scale(&C::rat(7, 10)), add(vec![Obj::D(Rep::A, rest)]))]),
        (Rep::A, n) if n % 2 == 1 => Some(Vec::new()),
        _ => None,
    }
}

/// The odd number of structure constants in a network of symmetric tensors
/// (color.h `simpli`).
pub(crate) fn odd_structure_constants(net: &Net) -> bool {
    net.objs.iter().all(|o| matches!(o, Obj::F(_) | Obj::D(..))) && net.count_f() % 2 == 1
}

/// Every invariance relation of the network: for each structure constant f
/// with a leg k on a symmetric tensor D, and each choice of its leg p that
/// moves, sum over the slots s of D of the network in which D's slot of k
/// holds p, slot s holds k and f takes D's former slot-s index:
/// sum_s f(.., k, .., e_s, ..) D(.., p, .., k at s, ..) = 0, the original
/// network being the term s = slot of k.
pub(crate) fn invariance_relations(net: &Net) -> Vec<Terms> {
    let occ = net.occurrences();
    let mut relations = Vec::new();
    for (fk, obj) in net.objs.iter().enumerate() {
        let Obj::F(f) = obj else { continue };
        for &k in f {
            let d = other(&occ, k, fk as u16) as usize;
            let Obj::D(rep, dp) = &net.objs[d] else { continue };
            let u = dp.iter().position(|&i| i == k).unwrap();
            for &p in f {
                if p == k {
                    continue;
                }
                let mut terms = vec![(Poly::one(), net.clone())];
                for (s, &e) in dp.iter().enumerate() {
                    if s == u {
                        continue;
                    }
                    let mut dp2 = dp.clone();
                    dp2[u] = p;
                    dp2[s] = k;
                    let f2 = [f[0], f[1], f[2]].map(|i| if i == p { e } else { i });
                    let mut objs: Vec<Obj> = net.objs.clone();
                    objs[fk] = Obj::F(f2);
                    objs[d] = Obj::D(*rep, dp2);
                    terms.push((Poly::one(), Net { objs, next: net.next }));
                }
                relations.push(terms);
            }
        }
    }
    relations
}

pub(crate) fn only_symmetric_or_f(net: &Net) -> bool {
    net.objs.iter().all(|o| matches!(o, Obj::D(..) | Obj::F(_)))
}

pub(crate) fn only_symmetric(net: &Net) -> bool {
    net.objs.iter().all(|o| matches!(o, Obj::D(..)))
}

#[allow(dead_code)]
pub(crate) fn fresh_mono_one() -> Mono {
    Mono::ONE
}

/// color.h's `dddff` pattern: three rank-four tensors D1, D2, D3 and two
/// structure constants f(a1, b1, i) f(a2, b2, i) with a's on D1, b's on D2,
/// and two legs of each of D1, D2 on D3.
pub(crate) fn dddff(net: &Net) -> Option<Poly> {
    if net.objs.len() != 5 {
        return None;
    }
    let fs: Vec<&[Ix; 3]> = net.objs.iter().filter_map(|o| if let Obj::F(f) = o { Some(f) } else { None }).collect();
    let ds: Vec<(usize, Rep, &Vec<Ix>)> = net
        .objs
        .iter()
        .enumerate()
        .filter_map(|(k, o)| if let Obj::D(r, i) = o { Some((k, *r, i)) } else { None })
        .collect();
    if fs.len() != 2 || ds.len() != 3 || ds.iter().any(|d| d.2.len() != 4) {
        return None;
    }
    let shared: Vec<Ix> = fs[0].iter().copied().filter(|i| fs[1].contains(i)).collect();
    if shared.len() != 1 {
        return None;
    }
    let i = shared[0];
    let holder = |x: Ix| ds.iter().position(|d| d.2.contains(&x));
    let legs = |f: &[Ix; 3]| -> Option<(usize, usize, Ix, Ix)> {
        let others: Vec<Ix> = f.iter().copied().filter(|&x| x != i).collect();
        Some((holder(others[0])?, holder(others[1])?, others[0], others[1]))
    };
    let (p0, q0, x0, y0) = legs(fs[0])?;
    let (p1, q1, x1, y1) = legs(fs[1])?;
    if p0 == q0 {
        return None;
    }
    // Orient both as f(D1 leg, D2 leg, i).
    let (d1, d2) = (p0, q0);
    let (a0, b0) = (x0, y0);
    let (a1, b1) = if p1 == d1 && q1 == d2 { (x1, y1) } else if p1 == d2 && q1 == d1 { (y1, x1) } else { return None };
    let d3 = 3 - d1 - d2;
    let common = |u: usize, v: usize| ds[u].2.iter().filter(|x| ds[v].2.contains(x)).count();
    if common(d1, d3) != 2 || common(d2, d3) != 2 || common(d1, d2) != 0 {
        return None;
    }
    let s = orientation(fs[0], a0, b0) * orientation(fs[1], a1, b1);
    let reps: Vec<Rep> = ds.iter().map(|d| d.1).collect();
    Some(super::names::dddff_value(&reps).scale(&C::rat(s, 1)))
}

#[cfg(test)]
mod tests {
    use super::{C, bernoulli};

    /// B_x / x! for every x, past the old table (which ended at x = 16 and
    /// silently dropped the terms of larger even x).
    #[test]
    fn bernoulli_weights_are_exact() {
        let factorial = |x: i64| (1..=x).map(|k| k as i128).product::<i128>();
        assert_eq!(bernoulli(1), Some(C::rat(1, 2)));
        for x in [3, 5, 7, 17, 19, 21] {
            assert_eq!(bernoulli(x), None, "B_{x}");
        }
        // B_x numerators and denominators (OEIS A027641 / A027642).
        let table: [(usize, i128, i128); 12] = [
            (2, 1, 6),
            (4, -1, 30),
            (6, 1, 42),
            (8, -1, 30),
            (10, 5, 66),
            (12, -691, 2730),
            (14, 7, 6),
            (16, -3617, 510),
            (18, 43867, 798),
            (20, -174611, 330),
            (22, 854513, 138),
            (24, -236364091, 2730),
        ];
        for (x, n, d) in table {
            let weight = bernoulli(x).unwrap_or_else(|| panic!("B_{x} is nonzero"));
            let expected = symbolica::domains::rational::Rational::new(n, d * factorial(x as i64));
            assert_eq!(weight.re, expected, "B_{x}/{x}!");
            assert!(weight.im.is_zero());
        }
    }
}
