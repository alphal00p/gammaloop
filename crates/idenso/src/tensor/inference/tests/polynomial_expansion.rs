use super::*;
use spenso::network::tags::SPENSO_TAG;
use symbolica::symbol;

/// Whether `value` divides by an expression.
fn divides(value: &Atom) -> bool {
    let mut divides = false;
    value.visitor(&mut |node| {
        if let AtomView::Pow(power) = node {
            let (base, exponent) = power.get_base_exp();
            divides |= !matches!(base, AtomView::Num(_))
                && Rational::try_from(exponent)
                    .is_ok_and(|exponent| exponent.numerator().is_negative());
        }
        !divides
    });
    divides
}

/// The selected polynomial expansion divides by nothing, keeps the interface,
/// and has the value of `raw` at three exact component assignments.
#[track_caller]
fn assert_polynomial_expansion(raw: &Atom, selector: &Atom, heads: &[Symbol], scalars: &[Atom]) {
    let source = SymbolicTensor::infer(raw.clone()).unwrap();
    let result = source.expanded(Some(selector.as_view()), true).unwrap();
    assert!(
        !divides(&result.expression),
        "{raw} in {selector}: {}",
        result.expression
    );
    assert_eq!(result.structure, source.structure);
    let admitted = SymbolicTensor::infer(result.expression.clone()).unwrap();
    assert_eq!(admitted.structure, source.structure);
    let symbols = scalars
        .iter()
        .map(|scalar| scalar.get_symbol().unwrap())
        .collect::<Vec<_>>();
    for shift in 0..3 {
        let mut heads = heads.to_vec();
        let mut scalars = symbols.clone();
        let rotation = shift % heads.len().max(1);
        heads.rotate_left(rotation);
        scalars.rotate_left(shift % symbols.len());
        let components = VectorComponents {
            heads,
            matrices: Vec::new(),
            scalars,
        };
        let bound = HashMap::new();
        assert_eq!(
            components.value(result.expression.as_view(), &bound),
            components.value(raw.as_view(), &bound),
            "{raw} in {selector}"
        );
    }
}

#[test]
fn selected_polynomial_expansion_adds_no_denominator() {
    test_support::test_initialize();
    let [a, b, m, x, y, z] =
        ["a", "b", "m", "x", "y", "z"].map(|name| Atom::var(symbol!(&format!("poly_den_{name}"))));
    // Scalar inputs whose polynomial powers Symbolica divides by the leading
    // coefficient of the raised base.
    let scalar = [
        ((&x * (&z + &m).pow(2) + &y).pow(2), &m),
        ((&x * (&x + &z) * &y + &z).pow(2), &x),
        ((&x * (&m * &a + &z).pow(2) + &y).pow(2), &m),
        ((&x * (&m * &a * &b + &z).pow(2) + &y).pow(2), &m),
        ((&y * m.pow(2) + &z * &m + &x).pow(2), &m),
        (((&x + &y) * m.pow(2) + &z * &m + 1).pow(2), &m),
        ((&x * m.pow(2) + &z * &m + &y).pow(3), &m),
    ];
    let scalars = [a, b, m.clone(), x.clone(), y.clone(), z.clone()];
    for (raw, selector) in &scalar {
        assert_polynomial_expansion(raw, selector, &[], &scalars);
    }

    let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
    let [e, i, j] = [98_281, 98_282, 98_283].map(|index| {
        rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
            .to_atom()
    });
    let heads = ["v0", "v1", "v2", "v3", "v4"]
        .map(|name| SPENSO_TAG.rank_one_tensor_symbol(&format!("poly_den_{name}")));
    let [v0, v1, v2, v3, v4] =
        heads.map(|head| move |index: &Atom| FunctionBuilder::new(head).add_arg(index).finish());
    // Without the fallback these divide by `v3(e)`, which is not admissible,
    // or by closed products, which vanish at lightlike points.
    let tensor = [
        ((&x * (&x + &z) * v3(&e) + v1(&e)).pow(2), x.clone()),
        (
            ((Atom::num(2) * &x + &z) * 3 * &x * v3(&e) + v1(&e)).pow(2),
            x.clone(),
        ),
        (
            ((&x + &z) * &x * v3(&e) + v1(&e)).pow(2) * v2(&i) * v4(&i),
            x.clone(),
        ),
        (
            (&x * (&z + v0(&j) * v1(&j)).pow(2) * v2(&i) * v3(&i) + &y).pow(2),
            v0(&j),
        ),
    ];
    for (raw, selector) in &tensor {
        assert_polynomial_expansion(raw, selector, &heads, &scalars);
    }

    // An input denominator in a changed result also takes the plain route.
    let rational = SymbolicTensor::infer((&x * &m + &y).pow(2) / (&x + &z)).unwrap();
    assert_eq!(
        rational.expanded(Some(m.as_view()), true).unwrap(),
        rational.expanded(Some(m.as_view()), false).unwrap()
    );
}

#[test]
fn selected_polynomial_expansion_stays_finite_where_the_input_is() {
    test_support::test_initialize();
    let [m, w, x, y, z] =
        ["m", "w", "x", "y", "z"].map(|name| Atom::var(symbol!(&format!("poly_fin_{name}"))));
    let opaque = FunctionBuilder::new(symbol!("poly_fin_f"))
        .add_arg(y.pow(-1))
        .finish();
    let square = (&y * m.pow(2) + &z * &m + &x).pow(2);
    // An existing `1/y` does not license the backend's spurious `1/y`: at
    // `w = 0` or `f(1/y) = 0`, and then `y = 0`, the input is `(z m + x)^2`.
    for (raw, discarded) in [(&square + &w / &y, &w), (&square + &opaque, &opaque)] {
        let source = SymbolicTensor::infer(raw).unwrap();
        let result = source.expanded(Some(m.as_view()), true).unwrap();
        assert_eq!(result, source.expanded(Some(m.as_view()), false).unwrap());
        let finite = result
            .expression
            .replace(discarded.clone())
            .with(Atom::num(0))
            .replace(y.clone())
            .with(Atom::num(0));
        assert_eq!(finite, (&z * &m + &x).pow(2).expand_in(m.as_view()));
    }
}

#[test]
fn selected_polynomial_expansion_falls_back_above_the_tenth_power() {
    test_support::test_initialize();
    let [m, x, z] = ["m", "x", "z"].map(|name| Atom::var(symbol!(&format!("poly_high_{name}"))));
    let source = SymbolicTensor::infer((&x * m.pow(3) + &z * &m + 1).pow(11)).unwrap();
    let result = source.expanded(Some(m.as_view()), true).unwrap();
    assert!(!divides(&result.expression), "{}", result.expression);
    assert_eq!(result, source.expanded(Some(m.as_view()), false).unwrap());
}

#[test]
fn selected_polynomial_fallback_uses_the_plain_grouping() {
    test_support::test_initialize();
    let [u, v, x, y, z] =
        ["u", "v", "x", "y", "z"].map(|name| Atom::var(symbol!(&format!("poly_group_{name}"))));
    let m = Atom::var(symbol!("poly_group_m"));
    let call = |head: &str, argument: &Atom| {
        FunctionBuilder::new(symbol!(&format!("poly_group_{head}")))
            .add_arg(argument)
            .finish()
    };
    let selected = call("f", &m);
    // Only the first summand divides on the polynomial route, but the
    // fallback also regroups the second one.
    let raw = (&x * (&z + &selected).pow(2) + &y).pow(2)
        + (&u + &v) * (Atom::num(1) + &selected + call("h", &selected));
    let source = SymbolicTensor::infer(raw.clone()).unwrap();
    let result = source.expanded(Some(selected.as_view()), true).unwrap();
    assert_eq!(
        result,
        source.expanded(Some(selected.as_view()), false).unwrap()
    );
    assert_ne!(
        result.expression,
        raw.expand_via_poly::<i32, Atom>(Some(selected.clone()))
    );
}

/// Seeded inputs shaped like the randomized campaigns' polynomial-route
/// failures: powers of sums with up to three terms in the selector, whose
/// coefficients are scalars, scalar sums or closed tensor pairs.
struct PolynomialCases {
    state: u64,
    next_index: usize,
    scalars: [Atom; 3],
    selector: Atom,
    heads: [Symbol; 3],
    open: [Symbol; 2],
}

impl PolynomialCases {
    fn new(seed: u64) -> Self {
        Self {
            state: seed,
            next_index: 0,
            scalars: ["y", "z", "w"].map(|name| Atom::var(symbol!(&format!("poly_prop_{name}")))),
            selector: Atom::var(symbol!("poly_prop_x")),
            heads: ["p", "q", "r"]
                .map(|name| SPENSO_TAG.rank_one_tensor_symbol(&format!("poly_prop_{name}"))),
            open: ["u", "v"]
                .map(|name| SPENSO_TAG.rank_one_tensor_symbol(&format!("poly_prop_{name}"))),
        }
    }

    /// The splitmix64 stream.
    fn next(&mut self, bound: usize) -> usize {
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut value = self.state;
        value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        value ^= value >> 31;
        (value % bound as u64) as usize
    }

    fn index(&mut self) -> Atom {
        self.next_index += 1;
        ExtendibleReps::EUCLIDEAN
            .new_rep(Dimension::Concrete(2))
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(98_300 + self.next_index))
            .to_atom()
    }

    fn vector(head: Symbol, index: &Atom) -> Atom {
        FunctionBuilder::new(head).add_arg(index).finish()
    }

    fn pair(&mut self) -> Atom {
        let index = self.index();
        let [left, right] = [self.next(3), self.next(3)].map(|choice| self.heads[choice]);
        Self::vector(left, &index) * Self::vector(right, &index)
    }

    fn scalar(&mut self) -> Atom {
        self.scalars[self.next(3)].clone()
    }

    /// A coefficient free of the selector.
    fn coefficient(&mut self) -> Atom {
        match self.next(6) {
            0 => Atom::num(self.next(3) as i64 + 1),
            1 => self.scalar(),
            2 => self.scalar() + self.scalar() + 1,
            3 => self.pair(),
            4 => self.scalar() * (self.pair() + self.scalar()),
            _ => self.scalar() * self.scalar() + self.pair(),
        }
    }

    /// A sum of up to three powers of the selector, with an optional open
    /// vector on each term, and the power to raise it to.
    fn case(&mut self) -> (Atom, Atom) {
        self.next_index = 0;
        let open = (self.next(3) == 0).then(|| self.index());
        let degree = self.next(2) + 1;
        let mut base = Atom::num(0);
        for power in 0..=degree {
            let mut term = self.coefficient() * self.selector.pow(power as i64);
            if let Some(index) = &open {
                term *= Self::vector(self.open[self.next(2)], index);
            }
            base += term;
        }
        let exponent = if open.is_some() {
            2
        } else {
            self.next(2) as i64 + 2
        };
        let mut source = base.pow(exponent);
        if self.next(2) == 0 {
            source = source * self.coefficient() + self.coefficient() * &self.selector;
        }
        (source, self.selector.clone())
    }
}

#[test]
fn selected_polynomial_expansion_adds_no_denominator_on_seeded_inputs() {
    test_support::test_initialize();
    let mut cases = PolynomialCases::new(0x0281);
    let heads = [cases.heads.as_slice(), cases.open.as_slice()].concat();
    let scalars = [cases.scalars.as_slice(), &[cases.selector.clone()]].concat();
    for _ in 0..120 {
        let (raw, selector) = cases.case();
        assert_polynomial_expansion(&raw, &selector, &heads, &scalars);
    }
}
