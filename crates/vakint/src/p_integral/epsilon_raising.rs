//! Exact integer-family representation of analytic line raisings. A unit
//! massless G-scheme bubble at momentum q equals ep(q,1)/ε. The auxiliary
//! loops exist only at the terminal integration boundary.
use super::{PGraph, PLine, invalid};
use crate::routing::linear_momentum;
use crate::{
    CallerProgress, CandidateRoutedApplicationObserver, EvaluationMethod, EvaluationOrder,
    LoopNormalizationFactor, RustRedEvaluationOptions, RustRedMasterMode, Vakint, VakintError,
    VakintExpression, VakintSettings, get_integer_from_atom, symbols::S,
    utils::vakint_macros::vk_symbol,
};
use std::collections::BTreeMap;
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};

pub(crate) fn has_marker(input: AtomView<'_>) -> bool {
    if !input.contains_symbol(S.ep) {
        return false;
    }
    let mut found = false;
    input.visitor(&mut |part| {
        found |= matches!(part, AtomView::Fun(f) if f.get_symbol() == S.ep);
        !found
    });
    found
}

fn factors(view: AtomView<'_>) -> Vec<AtomView<'_>> {
    match view {
        AtomView::Mul(product) => product.iter().collect(),
        _ => vec![view],
    }
}

/// Collect multiplicative ep factors into each topology before matching its
/// integer family. Coefficients are grouped at fixed powers of all markers.
fn collect_markers(numerator: &Atom, integral: &Atom) -> Result<Atom, VakintError> {
    Ok(Atom::add_many(
        marker_terms(numerator, integral)?
            .into_iter()
            .map(|(topology, coefficient)| topology * coefficient),
    ))
}

/// One topology and its numerator per monomial of analytic weights.
pub(super) fn marker_terms(
    numerator: &Atom,
    integral: &Atom,
) -> Result<Vec<(Atom, Atom)>, VakintError> {
    let props = super::input::topology_argument(integral.as_view())
        .ok_or_else(|| invalid("ep requires one topo argument containing a propagator product"))?;
    if !has_marker(numerator.as_view()) {
        return Ok(vec![(integral.clone(), numerator.clone())]);
    }
    let terms = marker_polynomial(numerator.as_view())?;
    let mut result = Vec::with_capacity(terms.len());
    for (powers, coefficient) in terms {
        let mut tagged = props.to_owned();
        for power in factors(powers.as_view()) {
            if power.is_one() {
                continue;
            }
            let (q, m) = marker_power(power)?;
            tagged *= S.ep.call((q, m));
        }
        result.push((S.topo.call(tagged), coefficient));
    }
    Ok(result)
}

type MarkerPolynomial = BTreeMap<Atom, Atom>;

fn sum_coefficients(groups: BTreeMap<Atom, Vec<Atom>>) -> MarkerPolynomial {
    groups
        .into_iter()
        .filter_map(|(power, terms)| {
            let coefficient = Atom::add_many(terms);
            (!coefficient.is_zero()).then_some((power, coefficient))
        })
        .collect()
}

fn multiply_markers(lhs: &MarkerPolynomial, rhs: &MarkerPolynomial) -> MarkerPolynomial {
    let mut groups = BTreeMap::<Atom, Vec<Atom>>::new();
    for (x, a) in lhs {
        for (y, b) in rhs {
            groups.entry(x * y).or_default().push(a * b);
        }
    }
    sum_coefficients(groups)
}

/// Distribute only factors that contain ep. Scalar coefficients of a
/// closed tensor can themselves be large factored polynomials. Even when a
/// marked sum is nested inside a product, each coefficient is summed once.
fn marker_polynomial(value: AtomView<'_>) -> Result<MarkerPolynomial, VakintError> {
    if !has_marker(value) {
        return Ok(BTreeMap::from([(Atom::one(), value.to_owned())]));
    }
    match value {
        AtomView::Fun(f) if f.get_symbol() == S.ep => {
            let _ = marker_power(value)?;
            Ok(BTreeMap::from([(value.to_owned(), Atom::one())]))
        }
        AtomView::Add(sum) => {
            let mut by_weight = BTreeMap::<Atom, Vec<Atom>>::new();
            for term in sum.iter() {
                let (weight, scalar): (Vec<_>, Vec<_>) = factors(term)
                    .into_iter()
                    .partition(|factor| has_marker(*factor));
                by_weight
                    .entry(Atom::mul_many(weight))
                    .or_default()
                    .push(Atom::mul_many(scalar));
            }
            let mut groups = BTreeMap::<Atom, Vec<Atom>>::new();
            for (weight, terms) in by_weight {
                let scalar = Atom::add_many(terms);
                for (power, coefficient) in marker_polynomial(weight.as_view())? {
                    groups.entry(power).or_default().push(coefficient * &scalar);
                }
            }
            Ok(sum_coefficients(groups))
        }
        AtomView::Mul(product) => {
            let (weights, scalars): (Vec<_>, Vec<_>) =
                product.iter().partition(|factor| has_marker(*factor));
            let mut poly = BTreeMap::from([(Atom::one(), Atom::one())]);
            for weight in weights {
                poly = multiply_markers(&poly, &marker_polynomial(weight)?);
            }
            let scalar = Atom::mul_many(scalars);
            for coefficient in poly.values_mut() {
                *coefficient *= &scalar;
            }
            Ok(poly)
        }
        AtomView::Pow(power) => {
            let exp = get_integer_from_atom(power.get_exp())
                .filter(|n| *n >= 0)
                .ok_or_else(|| invalid("the numerator must be polynomial in ep markers"))?;
            if matches!(power.get_base(), AtomView::Fun(f) if f.get_symbol() == S.ep) {
                let _ = marker_power(value)?;
                return Ok(BTreeMap::from([(value.to_owned(), Atom::one())]));
            }
            let mut base = marker_polynomial(power.get_base())?;
            let mut n = exp;
            let mut poly = BTreeMap::from([(Atom::one(), Atom::one())]);
            while n > 0 {
                if n % 2 != 0 {
                    poly = multiply_markers(&poly, &base);
                }
                n /= 2;
                if n > 0 {
                    base = multiply_markers(&base, &base);
                }
            }
            Ok(poly)
        }
        _ => Err(invalid("the numerator must be polynomial in ep markers")),
    }
}

fn marker_power(power: AtomView<'_>) -> Result<(Atom, i64), VakintError> {
    let (base, n) = if let AtomView::Pow(p) = power {
        (
            p.get_base(),
            get_integer_from_atom(p.get_exp())
                .filter(|n| *n >= 0)
                .ok_or_else(|| invalid("ep multiplicity must be a nonnegative integer"))?,
        )
    } else {
        (power, 1)
    };
    let AtomView::Fun(f) = base else {
        return Err(invalid("invalid ep monomial"));
    };
    if f.get_symbol() != S.ep || f.get_nargs() != 2 {
        return Err(invalid("expected ep(momentum,count)"));
    }
    let mut args = f.iter();
    let q = args.next().unwrap();
    let m = get_integer_from_atom(args.next().unwrap())
        .filter(|m| *m >= 0)
        .and_then(|m| m.checked_mul(n))
        .ok_or_else(|| invalid("invalid ep count"))?;
    Ok((q.expand(), m))
}

/// Assemble the routed, massless p-integral at a tensor evaluator's terminal
/// boundary. Recognize a dimensionless power only when its base equals an
/// actual propagator square, including its expanded bilinear form.
pub fn with_analytic_line_powers(
    numerator: Atom,
    topology: Atom,
    epsilon: Symbol,
    mu_squared: Symbol,
) -> Result<Atom, VakintError> {
    Ok(mark_internal_powers(numerator, &topology, epsilon, mu_squared)? * topology)
}

pub(super) fn mark_internal_powers(
    numerator: Atom,
    topology: &Atom,
    epsilon: Symbol,
    mu_squared: Symbol,
) -> Result<Atom, VakintError> {
    let mut candidates = std::collections::BTreeSet::new();
    numerator.visitor(&mut |part| {
        if let AtomView::Pow(p) = part
            && p.get_exp().contains_symbol(epsilon)
            && p.get_base().contains_symbol(S.k)
        {
            candidates.insert(part.to_owned());
            return false;
        }
        true
    });
    if candidates.is_empty() {
        return Ok(numerator);
    }
    let (graph, _) = marked_graph(topology)?;
    let mu = Atom::var(mu_squared);
    let mut squares = std::collections::BTreeMap::new();
    for line in graph.lines {
        let parts = linear_momentum(line.momentum.as_view())?;
        let mut expanded = Atom::Zero;
        for (a, ca) in &parts {
            for (b, cb) in &parts {
                expanded += Atom::num(ca * cb) * S.dot.call((a, b));
            }
        }
        for square in [expanded, S.dot.call((&line.momentum, &line.momentum))] {
            squares.insert((-square / &mu).expand(), line.momentum.clone());
        }
    }
    let mut replacements = Vec::new();
    for candidate in candidates {
        let AtomView::Pow(power) = candidate.as_view() else {
            unreachable!()
        };
        let slope = -power.get_exp().to_owned() / Atom::var(epsilon);
        let m = get_integer_from_atom(slope.as_view())
            .filter(|m| *m >= 0)
            .ok_or_else(|| {
                invalid(format!(
                    "unsupported momentum-dependent exponent {}",
                    power.get_exp()
                ))
            })?;
        let q = squares.get(&power.get_base().expand()).ok_or_else(|| {
            invalid(format!(
                "analytic power does not match a routed propagator square: {candidate}"
            ))
        })?;
        replacements.push(symbolica::id::Replacement::new(
            candidate.to_pattern(),
            S.ep.call((q, m)).to_pattern(),
        ));
    }
    Ok(numerator.replace_multiple(&replacements))
}

pub(super) fn lifted(graph: &mut PGraph, momentum: &Atom, count: usize) -> Result<(), VakintError> {
    if graph.loops.saturating_add(count) > 5 {
        return Err(invalid(format!(
            "ep({momentum},{count}) on a {}-loop graph exceeds the five-loop RustRed family library",
            graph.loops
        )));
    }
    let opposite = -momentum;
    let original = graph
        .lines
        .iter_mut()
        .find(|line| line.momentum.expand() == *momentum || line.momentum.expand() == opposite)
        .ok_or_else(|| {
            invalid(format!(
                "ep momentum {momentum} matches no propagator up to sign"
            ))
        })?;
    let q = original.momentum.clone();
    let original_id = original.id;
    let mut line_id = graph.lines.iter().map(|l| l.id).max().unwrap();
    let mut vertex = graph
        .lines
        .iter()
        .map(|l| l.ends.0.max(l.ends.1))
        .max()
        .unwrap();
    let mut heads = std::collections::BTreeSet::new();
    for line in &graph.lines {
        heads.extend(linear_momentum(line.momentum.as_view())?.into_keys());
    }
    let mut candidate = 0i64;
    for _ in 0..count {
        let fresh = loop {
            candidate += 1;
            let k = S.k.call(candidate);
            if heads.insert(k.clone()) {
                break k;
            }
        };
        vertex = vertex
            .checked_add(1)
            .ok_or_else(|| invalid("no unused vertex id for ep"))?;
        let line = graph
            .lines
            .iter_mut()
            .find(|l| l.id == original_id)
            .unwrap();
        let sink = std::mem::replace(&mut line.ends.1, vertex);
        for momentum in [fresh.clone(), &q - fresh] {
            line_id = line_id
                .checked_add(1)
                .ok_or_else(|| invalid("no unused line id for ep"))?;
            graph.lines.push(PLine {
                id: line_id,
                ends: (vertex, sink),
                momentum,
                power: 1,
            });
        }
    }
    graph.loops += count;
    Ok(())
}

/// A two-line scalar bubble is known for any integer multiple of ε on either
/// line. This also covers a one-loop skeleton carrying five deleted loops,
/// beyond the ordinary five-loop family library.
fn analytic_bubble(
    graph: &PGraph,
    markers: &[(Atom, usize)],
    numerator: &Atom,
    settings: &VakintSettings,
) -> Result<Option<Atom>, VakintError> {
    if graph.loops != 1
        || graph.lines.len() != 2
        || graph
            .lines
            .iter()
            .any(|l| l.power != 1 || l.ends.0 == l.ends.1)
        || numerator.contains_symbol(S.k)
    {
        return Ok(None);
    }
    let mut shifts = [0i64; 2];
    for (q, m) in markers {
        let i = graph
            .lines
            .iter()
            .position(|l| {
                let p = l.momentum.expand();
                &p == q || p == -q
            })
            .ok_or_else(|| invalid(format!("ep momentum {q} matches no bubble line")))?;
        shifts[i] = shifts[i]
            .checked_add(i64::try_from(*m).map_err(|_| invalid("ep count overflow"))?)
            .ok_or_else(|| invalid("ep count overflow"))?;
    }
    let total = shifts[0]
        .checked_add(shifts[1])
        .filter(|n| *n <= i64::MAX - 2)
        .ok_or_else(|| invalid("ep count overflow"))?;
    let mut exact_settings = settings.clone();
    exact_settings.use_dot_product_notation = true;
    let reduced = super::ReducedPropagatorValue {
        scale_power: 0,
        reduced: super::pg_symbol().call((shifts[0], shifts[1])) / Atom::num(1 + total),
        spectator: None,
    };
    let value = super::finalize::finalize(
        &[reduced],
        1,
        &graph.external,
        &exact_settings,
        RustRedMasterMode::Retain,
        "analytic ep bubble",
    )?;
    let (p, _) = super::external_head(&graph.external)?;
    let log =
        Symbol::LOG.call(-S.dot.call((&p, &p)) / Atom::var(vk_symbol!(&settings.mu_r_sq_symbol)));
    Ok(Some(
        value
            * numerator
            * Symbol::EXP
                .call(-Atom::num(total) * Atom::var(vk_symbol!(&settings.epsilon_symbol)) * log),
    ))
}

pub(super) fn marked_graph(topology: &Atom) -> Result<(PGraph, Vec<(Atom, usize)>), VakintError> {
    let product = super::input::topology_argument(topology.as_view())
        .ok_or_else(|| invalid("ep requires one topo argument containing a propagator product"))?;
    let mut props = Atom::one();
    let mut markers = Vec::new();
    for factor in factors(product) {
        if !has_marker(factor) {
            props *= factor.to_owned();
            continue;
        }
        let (q, m) = marker_power(factor)?;
        let count = usize::try_from(m).map_err(|_| invalid("ep count overflow"))?;
        markers.push((q, count));
    }
    Ok((PGraph::parse(S.topo.call(props).as_view())?, markers))
}

impl Vakint {
    pub(crate) fn evaluate_epsilon_raised(
        &self,
        settings: &VakintSettings,
        input: AtomView,
        mut observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
        progress: Option<&CallerProgress>,
    ) -> Result<Atom, VakintError> {
        self.validate_settings(settings)?;
        let mut result = Atom::Zero;
        let eps = Atom::var(vk_symbol!(&settings.epsilon_symbol));
        for term in VakintExpression::split_integrals(input)? {
            if has_marker(term.numerator.as_view()) {
                let marked = collect_markers(&term.numerator, &term.integral)?;
                result += self.evaluate_observed(
                    settings,
                    marked.as_view(),
                    observation.as_deref_mut(),
                    progress,
                )?;
                continue;
            }
            if !has_marker(term.integral.as_view()) {
                result += self.evaluate_observed(
                    settings,
                    Atom::from(term).as_view(),
                    observation.as_deref_mut(),
                    progress,
                )?;
                continue;
            }
            let Some(EvaluationMethod::RustRed(options)) = settings.evaluation_order.0.first()
            else {
                return Err(invalid(
                    "ep requires RustRed as the first evaluation method",
                ));
            };
            let (mut graph, markers) = marked_graph(&term.integral)?;
            let physical = graph.loops;
            let extra = markers
                .iter()
                .try_fold(0usize, |n, (_, m)| n.checked_add(*m))
                .ok_or_else(|| invalid("ep count overflow"))?;
            let analytic = analytic_bubble(&graph, &markers, &term.numerator, settings)?;
            let retained = if let Some(value) = analytic {
                value
            } else {
                for (q, m) in markers {
                    lifted(&mut graph, &q, m)?;
                }

                let mut inner = settings.clone();
                inner.evaluation_order =
                    EvaluationOrder(vec![EvaluationMethod::RustRed(RustRedEvaluationOptions {
                        master_mode: RustRedMasterMode::Retain,
                    })]);
                inner.integral_normalization_factor = LoopNormalizationFactor::GScheme;
                inner.use_dot_product_notation = true;
                let retained = self.evaluate_observed(
                    &inner,
                    (graph.to_atom() * term.numerator).as_view(),
                    observation.as_deref_mut(),
                    progress,
                )?;
                let norm = (settings.get_staged_integral_normalization_factor_atom()?
                    / inner.get_staged_integral_normalization_factor_atom()?)
                .replace(S.n_loops.to_pattern())
                .with(Atom::num(physical).to_pattern())
                .replace(Atom::var(vk_symbol!("log_mu_sq")).to_pattern())
                .with(
                    Symbol::LOG
                        .call(Atom::var(vk_symbol!(&settings.mu_r_sq_symbol)))
                        .to_pattern(),
                );
                retained * norm * eps.pow(extra as u64)
            };
            if options.master_mode == RustRedMasterMode::Retain {
                result += if settings.use_dot_product_notation {
                    retained
                } else {
                    Self::convert_from_dot_notation(retained.as_view())
                };
            } else {
                // refinalise uses target = terms - catalogue_loops - 1. The
                // caller's requested top order belongs to the physical graph.
                let mut final_settings = settings.clone();
                final_settings.number_of_terms_in_epsilon_expansion = final_settings
                    .number_of_terms_in_epsilon_expansion
                    .checked_add((graph.loops - physical) as i64)
                    .ok_or_else(|| invalid("epsilon depth overflow"))?;
                result +=
                    self.refinalise_p_retained(&final_settings, retained.as_view(), graph.loops)?;
            }
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::vakint_macros::vk_parse;

    #[test]
    fn malformed_marker_topologies_are_rejected_before_matching() {
        let marker = S.ep.call((S.k.call(1), 1));
        for topology in [
            vk_parse!("topo()").unwrap(),
            vk_parse!("topo(1,2)").unwrap(),
            vk_parse!("not_topo(1)").unwrap(),
        ] {
            assert!(marked_graph(&topology).is_err());
            assert!(marker_terms(&marker, &topology).is_err());
        }
    }

    #[test]
    fn batched_marker_collection_preserves_nested_polynomials_and_cancellations() {
        let q = S.k.call(1);
        let r = S.k.call(2);
        let x = S.ep.call((&q, 1));
        let y = S.ep.call((&r, 1));
        let a = vk_parse!("a").unwrap();
        let b = vk_parse!("b").unwrap();
        let props = vk_parse!("prop(1,edge(0,1),k(1),0,1)").unwrap();
        let topo = S.topo.call(&props);
        let numerator = (&x + &y).pow(2) * &a - &x * &x * &a + &x * &y * &b + &a + &b;
        let mut expected = BTreeMap::new();
        expected.insert(topo.clone(), &a + &b);
        expected.insert(S.topo.call(&props * &x * &y), Atom::num(2) * &a + &b);
        expected.insert(S.topo.call(&props * S.ep.call((&r, 2))), a);
        for input in [numerator.clone(), numerator.expand()] {
            let got = marker_terms(&input, &topo)
                .unwrap()
                .into_iter()
                .map(|(topology, coefficient)| (topology, coefficient.expand()))
                .collect::<BTreeMap<_, _>>();
            assert_eq!(got, expected);
        }
        assert!(marker_terms(&Symbol::LOG.call(x), &topo).is_err());
    }

    #[test]
    fn nested_weighted_sum_keeps_scalar_epsilon_poles_and_opaque_functions() {
        let marker = S.ep.call((S.k.call(1), 1));
        let epsilon = Atom::var(S.ep);
        let opaque = vk_symbol!("test_coefficient").call(&epsilon);
        let a = &opaque / epsilon.pow(4);
        let inner = &marker * &a + &opaque;
        let value = &a * inner.pow(2);
        let expected = value
            .coefficient_list::<i16>(&[&marker])
            .into_iter()
            .map(|(m, c)| (m, c.expand()))
            .collect::<BTreeMap<_, _>>();
        let actual = marker_polynomial(value.as_view())
            .unwrap()
            .into_iter()
            .map(|(m, c)| (m, c.expand()))
            .collect::<BTreeMap<_, _>>();
        assert_eq!(actual, expected);
    }

    #[test]
    #[ignore = "timing comparison on a large closed numerator"]
    fn time_many_terms_sharing_an_analytic_weight() {
        let q = S.k.call(1);
        let x = S.ep.call((&q, 1));
        let tagged = S.topo.call(&x);
        let scalar = vk_symbol!("bench_scalar");
        let terms = 16_384;
        // Scalar epsilon and the ep function share a symbol in this convention.
        // An outer product keeps the marked sum nested, as in the cut numerator.
        let eps = Atom::var(S.ep);
        let numerator = (Atom::one() + &eps)
            * Atom::add_many((0..terms).map(|i| &x * scalar.call((i, &eps)) / eps.pow(i % 6)));
        let start = std::time::Instant::now();
        let reference = numerator.coefficient_list::<i16>(&[&x]);
        let unbatched = start.elapsed();
        let start = std::time::Instant::now();
        let grouped = marker_terms(&numerator, &S.topo.call(Atom::one())).unwrap();
        let batched = start.elapsed();
        assert_eq!(grouped.len(), 1);
        assert_eq!(grouped[0].0, tagged);
        assert!((&grouped[0].1 - &reference[0].1).expand().is_zero());
        println!("{terms} terms: unbatched {unbatched:?}, batched {batched:?}");
    }

    #[test]
    fn raised_bubble_supplies_the_six_loop_simple_pole() {
        let engine = Vakint::new().unwrap();
        let mut settings = VakintSettings {
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
                RustRedEvaluationOptions {
                    master_mode: RustRedMasterMode::Symbolic,
                },
            )]),
            integral_normalization_factor: LoopNormalizationFactor::GScheme,
            number_of_terms_in_epsilon_expansion: 1,
            use_dot_product_notation: true,
            ..Default::default()
        };
        let eps = vk_symbol!(&settings.epsilon_symbol);
        let e = Atom::var(eps);
        let input = vk_parse!(
            "topo(prop(1,edge(0,1),k(7),0,1)*prop(2,edge(0,1),p(3)-k(7),0,1)*ep(-k(7),1))"
        )
        .unwrap()
            / e.pow(5);
        let unit = |a: Atom| {
            a.replace(vk_parse!("dot(p(3),p(3))").unwrap().to_pattern())
                .with(Atom::num(-1).to_pattern())
                .replace(Atom::var(vk_symbol!(&settings.mu_r_sq_symbol)).to_pattern())
                .with(Atom::one().to_pattern())
                .expand()
        };
        let pole = unit(engine.evaluate(&settings, input.as_view()).unwrap());
        // Expansion of Γ(2ε)Γ(1-2ε)Γ(1-ε)/(Γ(1+ε)Γ(2-3ε)),
        // divided by εG(1,1), at Q²=μ²=1, with five spectator poles.
        let z3 = symbolica::transcendental::zeta().call(3);
        let z5 = symbolica::transcendental::zeta().call(5);
        let pi4 = Atom::var(Symbol::PI).pow(4);
        let expected = Atom::num((1, 2)) / e.pow(6)
            + Atom::num((1, 2)) / e.pow(5)
            + Atom::num((3, 2)) / e.pow(4)
            + (Atom::num((9, 2)) - Atom::num(3) * &z3) / e.pow(3)
            + (Atom::num((27, 2)) - Atom::num(3) * &z3 - &pi4 / Atom::num(20)) / e.pow(2)
            + (Atom::num((81, 2)) - Atom::num(9) * z3 - Atom::num(21) * z5 - pi4 / Atom::num(20))
                / &e;
        assert!((&pole - expected).expand().is_zero(), "{pole}");
        settings.number_of_terms_in_epsilon_expansion += 2;
        let deep = unit(
            engine
                .evaluate_integral(&settings, input.as_view())
                .unwrap(),
        );
        assert!(
            (deep.series(eps, 0, -1).unwrap().to_atom() - pole)
                .expand()
                .is_zero()
        );
    }
    #[test]
    fn two_raised_lines_keep_the_scale_and_the_loop_budget() {
        let engine = Vakint::new().unwrap();
        let eps = vk_symbol!("ep");
        let e = Atom::var(eps);
        let mut settings = VakintSettings {
            epsilon_symbol: eps.get_name().into(),
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
                RustRedEvaluationOptions {
                    master_mode: RustRedMasterMode::Symbolic,
                },
            )]),
            integral_normalization_factor: LoopNormalizationFactor::GScheme,
            number_of_terms_in_epsilon_expansion: 2,
            use_dot_product_notation: true,
            ..Default::default()
        };
        let topo =
            vk_parse!("topo(prop(1,edge(0,1),k(8),0,1)*prop(2,edge(0,1),p(3)-k(8),0,1))").unwrap();
        let k = S.k.call(8);
        let p = S.p.call(3);
        let mu = vk_symbol!(&settings.mu_r_sq_symbol);
        let q = &p - &k;
        let qm = -S.dot.call((&q, &q)) / Atom::var(mu);
        let physical =
            with_analytic_line_powers(qm.pow(-&e) * S.ep.call((&k, 1)), topo.clone(), eps, mu)
                .unwrap();
        let v = engine.evaluate(&settings, physical.as_view()).unwrap();
        let log = Symbol::LOG.call(-S.dot.call((&p, &p)) / Atom::var(mu));
        // 1/(3ε) × (1-2ε)/(1-4ε) × exp(-3ε log(Q²/μ²)).
        assert!(
            (v - (Atom::one() / (Atom::num(3) * &e) + Atom::num((2, 3)) - &log))
                .expand()
                .is_zero()
        );
        // Five excised bubbles on one line: the gamma integral supplies
        // this scalar skeleton beyond the five-loop integer-family library.
        let six = topo * S.ep.call((&k, 5));
        let v = engine.evaluate(&settings, six.as_view()).unwrap();
        assert!(
            (v - (Atom::one() / (Atom::num(6) * &e) + Atom::num((5, 6)) - &log))
                .expand()
                .is_zero()
        );
        settings.number_of_terms_in_epsilon_expansion = 1;
        let v = engine
            .evaluate(&settings, (six / e.pow(5)).as_view())
            .unwrap();
        assert!(!v.is_zero());
    }
    #[test]
    fn five_renormalized_scalar_bubbles_and_a_two_loop_rustred_skeleton() {
        let engine = Vakint::new().unwrap();
        let mut settings = VakintSettings {
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
                RustRedEvaluationOptions {
                    master_mode: RustRedMasterMode::Symbolic,
                },
            )]),
            integral_normalization_factor: LoopNormalizationFactor::GScheme,
            use_dot_product_notation: true,
            number_of_terms_in_epsilon_expansion: 1,
            ..Default::default()
        };
        let eps = vk_symbol!(&settings.epsilon_symbol);
        let e = Atom::var(eps);
        let marker = vk_parse!("ep(k(1),1)").unwrap();
        let bubble =
            vk_parse!("prop(1,edge(0,1),k(1),0,1)*prop(2,edge(0,1),p(1)-k(1),0,1)").unwrap();
        let mu = Atom::var(vk_symbol!(&settings.mu_r_sq_symbol));
        let unit = |v: Atom| {
            v.replace(vk_parse!("dot(p(1),p(1))").unwrap().to_pattern())
                .with(Atom::num(-1).to_pattern())
                .replace(mu.to_pattern())
                .with(Atom::one().to_pattern())
                .expand()
        };
        let kernel = S.topo.call(&bubble) * ((&marker - Atom::one()) / &e).pow(5);
        let actual = unit(engine.evaluate(&settings, kernel.as_view()).unwrap());
        let z3 = symbolica::transcendental::zeta().call(3) / Atom::num(3);
        let z4 = Atom::var(Symbol::PI).pow(4) / Atom::num(180);
        let z5 = symbolica::transcendental::zeta().call(5);
        let a = Atom::num((1, 6));
        // Fifth binomial difference of the one-loop gamma integral G(1+jε,1).
        let expected = -&a / e.pow(6)
            + &a / e.pow(5)
            + &a / e.pow(4)
            + (&a - &z3) / e.pow(3)
            + (&a + &z3 - &z4) / e.pow(2)
            + (&a + z3 + z4 - z5) / &e;
        assert!((&actual - expected).expand().is_zero(), "{actual}");
        // Two successive massless two-point loops. The second unweighted
        // bubble is 1/ε at unit scale. This case uses the integer-family lift.
        let other =
            vk_parse!("prop(3,edge(1,2),k(2),0,1)*prop(4,edge(1,2),p(1)-k(2),0,1)").unwrap();
        let two = S.topo.call(bubble * other * marker);
        settings.number_of_terms_in_epsilon_expansion = 2;
        let actual = unit(engine.evaluate(&settings, two.as_view()).unwrap());
        assert!(
            (actual - Atom::num((1, 2)) / e.pow(2) - Atom::num((1, 2)) / e)
                .expand()
                .is_zero()
        );
    }
}
