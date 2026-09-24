//! Factor trace recipes before attaching user indices. The cached graph stores
//! integer factor identifiers and coefficients only; atoms live in one evaluation.

use std::collections::HashMap;

use symbolica::atom::{Atom, AtomView};

type Polynomial = Vec<(Vec<usize>, i32)>;

enum Node {
    Constant(i32),
    Product {
        factors: Vec<usize>,
        coefficient: i32,
        child: usize,
    },
    Sum(Vec<usize>),
}

pub(super) struct FactoredTrace {
    nodes: Vec<Node>,
    root: usize,
}

struct Compiler {
    nodes: Vec<Node>,
    memo: HashMap<Polynomial, usize>,
    factor_count: usize,
}

impl Compiler {
    fn push(&mut self, node: Node) -> usize {
        let id = self.nodes.len();
        self.nodes.push(node);
        id
    }

    fn product(
        &mut self,
        mut factors: Vec<usize>,
        mut coefficient: i32,
        mut child: usize,
    ) -> usize {
        if let Node::Product {
            factors: inner,
            coefficient: inner_coefficient,
            child: next,
        } = &self.nodes[child]
        {
            factors.extend(inner);
            coefficient *= inner_coefficient;
            child = *next;
        }
        if factors.is_empty() && coefficient == 1 {
            return child;
        }
        factors.sort_unstable();
        self.push(Node::Product {
            factors,
            coefficient,
            child,
        })
    }

    fn build(&mut self, mut terms: Polynomial) -> usize {
        for (factors, _) in &mut terms {
            factors.sort_unstable();
        }
        terms.sort_unstable();
        if let Some(&id) = self.memo.get(&terms) {
            return id;
        }
        let key = terms.clone();
        // Equivalent subpolynomials can differ by an overall sign.
        let negative = terms[0].1 < 0;
        if negative {
            for (_, coefficient) in &mut terms {
                *coefficient = -*coefficient;
            }
        }
        let result = if terms.len() == 1 && terms[0].0.is_empty() {
            self.push(Node::Constant(terms[0].1))
        } else {
            // Trace recipes are square-free in their positional factor IDs,
            // even when user indices later make two metric atoms equal.
            let common: Vec<_> = terms[0]
                .0
                .iter()
                .copied()
                .filter(|factor| terms.iter().all(|(ids, _)| ids.contains(factor)))
                .collect();
            if !common.is_empty() {
                for (ids, _) in &mut terms {
                    ids.retain(|id| !common.contains(id));
                }
                let child = self.build(terms);
                self.product(common, 1, child)
            } else {
                let mut children = Vec::new();
                loop {
                    let mut counts = vec![0; self.factor_count];
                    for (ids, _) in &terms {
                        for &id in ids {
                            counts[id] += 1;
                        }
                    }
                    let selected = counts
                        .iter()
                        .enumerate()
                        .max_by_key(|&(id, count)| (*count, std::cmp::Reverse(id)));
                    let Some((factor, &count)) = selected.filter(|(_, count)| **count >= 2) else {
                        for (ids, coefficient) in terms {
                            let child = self.push(Node::Constant(coefficient));
                            children.push(self.product(ids, 1, child));
                        }
                        break;
                    };
                    let mut inside = Vec::with_capacity(count);
                    let mut outside = Vec::with_capacity(terms.len() - count);
                    for (mut ids, coefficient) in terms {
                        if let Some(position) = ids.iter().position(|&id| id == factor) {
                            ids.remove(position);
                            inside.push((ids, coefficient));
                        } else {
                            outside.push((ids, coefficient));
                        }
                    }
                    let child = self.build(inside);
                    children.push(self.product(vec![factor], 1, child));
                    if outside.is_empty() {
                        break;
                    }
                    terms = outside;
                }
                if children.len() == 1 {
                    children[0]
                } else {
                    self.push(Node::Sum(children))
                }
            }
        };
        let result = if negative {
            self.product(Vec::new(), -1, result)
        } else {
            result
        };
        self.memo.insert(key, result);
        result
    }

    fn finish(self, root: usize) -> FactoredTrace {
        // Product-chain compression leaves some intermediate nodes unused.
        // Remove them once, so every later evaluation builds only live atoms.
        let mut live = vec![false; self.nodes.len()];
        let mut pending = vec![root];
        while let Some(id) = pending.pop() {
            if std::mem::replace(&mut live[id], true) {
                continue;
            }
            match &self.nodes[id] {
                Node::Product { child, .. } => pending.push(*child),
                Node::Sum(children) => pending.extend(children),
                Node::Constant(_) => {}
            }
        }
        let mut remap = vec![0; self.nodes.len()];
        let mut nodes = Vec::new();
        for (old, mut node) in self.nodes.into_iter().enumerate() {
            if !live[old] {
                continue;
            }
            match &mut node {
                Node::Product { child, .. } => *child = remap[*child],
                Node::Sum(children) => {
                    for child in children {
                        *child = remap[*child];
                    }
                }
                Node::Constant(_) => {}
            }
            remap[old] = nodes.len();
            nodes.push(node);
        }
        FactoredTrace {
            nodes,
            root: remap[root],
        }
    }
}

impl FactoredTrace {
    /// Build a runtime recipe without the square-free positional compiler.
    pub(super) fn unit_recipe() -> Self {
        Self {
            nodes: vec![Node::Constant(1)],
            root: 0,
        }
    }

    pub(super) fn push_product(
        &mut self,
        factors: Vec<usize>,
        coefficient: i32,
        child: usize,
    ) -> usize {
        if factors.is_empty() && coefficient == 1 {
            return child;
        }
        // Do not fuse integer coefficients here: generic words may exceed i32.
        let id = self.nodes.len();
        self.nodes.push(Node::Product {
            factors,
            coefficient,
            child,
        });
        id
    }

    pub(super) fn push_sum(&mut self, children: Vec<usize>) -> usize {
        if children.len() == 1 {
            return children[0];
        }
        let id = self.nodes.len();
        self.nodes.push(Node::Sum(children));
        id
    }

    pub(super) fn recipe_leaf_count(&self, root: usize) -> usize {
        let mut counts = Vec::with_capacity(self.nodes.len());
        for node in &self.nodes {
            counts.push(match node {
                Node::Constant(c) => usize::from(*c != 0),
                Node::Product {
                    coefficient, child, ..
                } => {
                    if *coefficient == 0 {
                        0
                    } else {
                        counts[*child]
                    }
                }
                Node::Sum(children) => children.iter().map(|&child| counts[child]).sum(),
            });
        }
        counts[root]
    }

    /// Visit each final monomial with one reusable factor stack. Repeated IDs
    /// are retained, and exact coefficients multiply in the arbitrary-size domain.
    pub(super) fn visit_recipe_leaves(
        &self,
        root: usize,
        unit: &symbolica::domains::integer::Integer,
        emit: &mut impl FnMut(&[usize], &symbolica::domains::integer::Integer),
    ) {
        use symbolica::domains::{
            RingOps,
            integer::{Integer, Z},
        };
        fn visit(
            nodes: &[Node],
            id: usize,
            factors: &mut Vec<usize>,
            coefficient: &Integer,
            emit: &mut impl FnMut(&[usize], &Integer),
        ) {
            match &nodes[id] {
                Node::Constant(c) => {
                    if *c == 1 {
                        emit(factors, coefficient);
                    } else if *c != 0 {
                        emit(factors, &Z.mul(coefficient, &Integer::from(*c)));
                    }
                }
                Node::Product {
                    factors: ids,
                    coefficient: c,
                    child,
                } => {
                    if *c == 0 {
                        return;
                    }
                    let old_len = factors.len();
                    factors.extend(ids);
                    if *c == 1 {
                        visit(nodes, *child, factors, coefficient, emit);
                    } else {
                        visit(
                            nodes,
                            *child,
                            factors,
                            &Z.mul(coefficient, &Integer::from(*c)),
                            emit,
                        );
                    }
                    factors.truncate(old_len);
                }
                Node::Sum(children) => {
                    for &child in children {
                        visit(nodes, child, factors, coefficient, emit);
                    }
                }
            }
        }
        visit(&self.nodes, root, &mut Vec::new(), unit, emit);
    }
    /// Compile square-free positional factor recipes, with IDs below `factor_count`.
    pub(super) fn new(terms: Polynomial, factor_count: usize) -> Self {
        if terms.is_empty() {
            return Self {
                nodes: vec![Node::Constant(0)],
                root: 0,
            };
        }
        let mut compiler = Compiler {
            nodes: Vec::new(),
            memo: HashMap::new(),
            factor_count,
        };
        let root = compiler.build(terms);
        compiler.finish(root)
    }

    pub(super) fn evaluate(&self, factors: &[Atom], unit: AtomView<'_>) -> Atom {
        let mut values: Vec<Atom> = Vec::with_capacity(self.nodes.len());
        for node in &self.nodes {
            let value = match node {
                Node::Constant(coefficient) => Atom::num(*coefficient),
                Node::Product {
                    factors: ids,
                    coefficient,
                    child,
                } => {
                    let coefficient_atom = Atom::num(*coefficient);
                    Atom::mul_many(
                        std::iter::once(values[*child].as_view())
                            .chain((*coefficient != 1).then_some(coefficient_atom.as_view()))
                            .chain(ids.iter().map(|&id| factors[id].as_view())),
                    )
                }
                Node::Sum(children) => {
                    Atom::add_many(children.iter().map(|&child| values[child].as_view()))
                }
            };
            values.push(value);
        }
        Atom::mul_many([unit, values[self.root].as_view()])
    }

    /// Reconstruct formal integer coefficients for proof tests without expanding
    /// or otherwise rewriting a user's tensor expression.
    #[cfg(test)]
    pub(super) fn coefficient_map(&self) -> std::collections::BTreeMap<Vec<usize>, i32> {
        use std::collections::BTreeMap;

        let mut values: Vec<BTreeMap<Vec<usize>, i32>> = Vec::new();
        for node in &self.nodes {
            let polynomial = match node {
                Node::Constant(coefficient) => BTreeMap::from([(Vec::new(), *coefficient)]),
                Node::Product {
                    factors,
                    coefficient,
                    child,
                } => values[*child]
                    .iter()
                    .map(|(ids, value)| {
                        let mut ids = ids.clone();
                        ids.extend(factors);
                        ids.sort_unstable();
                        (ids, value * coefficient)
                    })
                    .collect(),
                Node::Sum(children) => {
                    let mut out = BTreeMap::new();
                    for &child in children {
                        for (ids, coefficient) in &values[child] {
                            *out.entry(ids.clone()).or_default() += coefficient;
                        }
                    }
                    out
                }
            };
            values.push(polynomial);
        }
        let mut result = values.swap_remove(self.root);
        result.retain(|_, coefficient| *coefficient != 0);
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    #[test]
    fn shared_integer_subpolynomials_preserve_every_coefficient() {
        // Dense square-free polynomials exercise shared factors, signs, constant
        // terms and magnitudes beyond the current short-trace recipe tables.
        for count in 1..=8 {
            let terms: Polynomial = (0usize..1 << count)
                .filter_map(|mask| {
                    let coefficient = (mask * 7 % 11) as i32 - 5;
                    (coefficient != 0).then(|| {
                        (
                            (0..count).filter(|i| mask & (1 << i) != 0).collect(),
                            coefficient,
                        )
                    })
                })
                .collect();
            let expected: BTreeMap<_, _> = terms.iter().cloned().collect();
            assert_eq!(FactoredTrace::new(terms, count).coefficient_map(), expected);
        }
        assert!(
            FactoredTrace::new(Vec::new(), 0)
                .coefficient_map()
                .is_empty()
        );
    }

    #[test]
    fn evaluation_respects_trace_unit_and_normalizes_numeric_values() {
        let terms = vec![(vec![0, 1], 2), (vec![0, 2], -3), (vec![1, 2], 1)];
        let factored = FactoredTrace::new(terms.clone(), 3);
        for sample in -2..=2 {
            let factors = [
                Atom::num(sample),
                Atom::num(sample + 1),
                Atom::num(2 - sample),
            ];
            for unit in [Atom::num(4), Atom::var(symbolica::symbol!("trace_unit"))] {
                let expected = &unit
                    * Atom::add_many(terms.iter().map(|(ids, coefficient)| {
                        Atom::mul_many(
                            std::iter::once(Atom::num(*coefficient))
                                .chain(ids.iter().map(|&id| factors[id].clone())),
                        )
                    }));
                assert_eq!(factored.evaluate(&factors, unit.as_view()), expected);
            }
        }
    }
}
