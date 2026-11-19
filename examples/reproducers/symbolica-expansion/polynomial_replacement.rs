//! Compare monomial-at-a-time replacement with grouped sparse Horner evaluation.
//! Usage: polynomial_replacement check | monomial|horner [N=256] [degree=4]
//! Only polynomial substitution is timed; the dependency is unmodified Symbolica.
use std::{hint::black_box, sync::Arc, time::Instant};
use symbolica::{
    domains::{Ring, finite_field::Zp, integer::Z, rational::Q},
    poly::{PolyVariable, polynomial::MultivariatePolynomial},
    symbol,
};

type Poly<R> = MultivariatePolynomial<R, u16>;

fn grouped_horner<R: Ring>(source: &Poly<R>, variable: usize, rhs: &Poly<R>) -> Poly<R> {
    assert_eq!(source.get_vars_ref(), rhs.get_vars_ref());
    if rhs.is_constant() {
        return source.replace_with_poly(variable, rhs);
    }
    // Reuse the public sparse coefficient-list operation. Its keys have the
    // original variable-map width; coefficient polynomials have x exponent 0.
    // Unlike a degree-indexed Vec, storage does not grow with exponent gaps.
    let mut coefficients: Vec<_> = source
        .to_multivariate_polynomial_list(&[variable], true)
        .into_iter()
        .map(|(powers, coefficient)| (powers[variable], coefficient))
        .collect();
    coefficients.sort_unstable_by_key(|(degree, _)| std::cmp::Reverse(*degree));
    let mut coefficients = coefficients.into_iter();
    let Some((mut degree, mut result)) = coefficients.next() else {
        return source.zero();
    };
    for (next_degree, coefficient) in coefficients {
        result = &result * &rhs.pow(usize::from(degree - next_degree)) + coefficient;
        degree = next_degree;
    }
    if degree != 0 {
        result = &result * &rhs.pow(usize::from(degree));
    }
    result
}

fn variables() -> Arc<Vec<PolyVariable>> {
    vec![
        symbol!("polynomial_replacement::x").into(),
        symbol!("polynomial_replacement::y").into(),
        symbol!("polynomial_replacement::z").into(),
    ]
    .into()
}

fn polynomial<R: Ring>(ring: &R, terms: &[(i64, [u16; 3])]) -> Poly<R> {
    let mut result = Poly::new(ring, Some(terms.len()), variables());
    for (coefficient, powers) in terms {
        result.append_monomial(ring.nth((*coefficient).into()), powers);
    }
    result
}

fn check_domain<R: Ring>(ring: &R) -> usize {
    let zero = polynomial(ring, &[]);
    let one = polynomial(ring, &[(1, [0, 0, 0])]);
    let x = polynomial(ring, &[(1, [1, 0, 0])]);
    let y = polynomial(ring, &[(1, [0, 1, 0])]);
    let z = polynomial(ring, &[(1, [0, 0, 1])]);
    let constant = polynomial(ring, &[(2, [0, 0, 0])]);
    let source = &(&x.pow(3) * &y) + &(&x.pow(2) * &z) + one.clone();
    let right_sides = [
        zero.clone(),
        one.clone(),
        constant.clone(),
        &y + &z,
        &x + &y,
    ];
    let sources = [zero.clone(), one.clone(), y.clone(), source];
    let mut checks = 0;
    for source in &sources {
        for rhs in &right_sides {
            let monomial = source.replace_with_poly(0, rhs);
            let horner = grouped_horner(source, 0, rhs);
            assert_eq!(monomial, horner);
            // Independent exact evaluations also cover self-reference: only
            // the source's x coordinate is replaced, not the RHS recursively.
            for point in [[1, 2, 3], [2, -1, 4], [0, 3, -2]] {
                let original = point.map(|v| ring.nth(v.into()));
                let mut substituted = original.clone();
                substituted[0] = rhs.replace_all(&original);
                assert_eq!(
                    horner.replace_all(&original),
                    source.replace_all(&substituted)
                );
            }
            checks += 1;
        }
    }
    // Explicit independent identities; these do not use the old replacement
    // route as their sole oracle. Keep a positive minimum and a large gap.
    let sparse = &(&x.pow(257) * &y) + &(&x.pow(3) * &z);
    assert_eq!(
        grouped_horner(&sparse, 0, &y),
        &y.pow(258) + &(&y.pow(3) * &z)
    );
    assert_eq!(grouped_horner(&x.pow(2), 0, &(&x + &y)), (&x + &y).pow(2));
    let cancels = &(&x.pow(2) - &(&x * &y)) + &z;
    assert_eq!(grouped_horner(&cancels, 0, &y), z);
    checks + 3
}

fn checks() {
    let count = check_domain(&Q) + check_domain(&Z) + check_domain(&Zp::new(17));
    let x = polynomial(&Q, &[(1, [1, 0, 0])]);
    let y = polynomial(&Q, &[(1, [0, 1, 0])]);
    let rational = x.pow(3).mul_coeff((2, 3).into()) + y.clone().mul_coeff((-5, 7).into());
    let expected = (&x + &y).pow(3).mul_coeff((2, 3).into()) + y.clone().mul_coeff((-5, 7).into());
    assert_eq!(grouped_horner(&rational, 0, &(&x + &y)), expected);
    assert_eq!(rational.replace_with_poly(0, &(&x + &y)), expected);
    let field = Zp::new(17);
    let x = polynomial(&field, &[(1, [1, 0, 0])]);
    let y = polynomial(&field, &[(1, [0, 1, 0])]);
    // In characteristic 17 all intermediate binomial coefficients vanish.
    assert_eq!(
        grouped_horner(&x.pow(17), 0, &(&x + &y)),
        &x.pow(17) + &y.pow(17)
    );
    println!(
        "{{\"correctness_cases\":{},\"domains\":[\"Q\",\"Z\",\"F17\"],\"exact\":true}}",
        count + 2
    );
}

fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.first().is_some_and(|arg| arg == "check") {
        assert_eq!(args.len(), 1);
        checks();
        return;
    }
    assert!(
        (1..=3).contains(&args.len()),
        "usage: polynomial_replacement check | monomial|horner [N=256] [degree=4]"
    );
    let route = args[0].as_str();
    assert!(matches!(route, "monomial" | "horner"));
    let n = args
        .get(1)
        .map_or(256, |v| v.parse::<u16>().expect("integer N"));
    let degree = args
        .get(2)
        .map_or(4, |v| v.parse::<u16>().expect("integer degree"));
    assert!((1..=4096).contains(&n) && (1..=16).contains(&degree));
    // P(x,y,z) = (sum_i (i mod 7 + 1)y^i)(x + ... + x^degree).
    // Three variables at every size; replacement x -> y+z shares one map.
    let mut coefficient = polynomial(&Q, &[]);
    for i in 0..n {
        coefficient.append_monomial((i32::from(i % 7 + 1), 1).into(), &[0, i, 0]);
    }
    let x = polynomial(&Q, &[(1, [1, 0, 0])]);
    let rhs = polynomial(&Q, &[(1, [0, 1, 0]), (1, [0, 0, 1])]);
    let mut powers = x.zero();
    for d in 1..=degree {
        powers = powers + x.pow(usize::from(d));
    }
    let source = &coefficient * &powers;
    let source_before = source.clone();
    let rhs_before = rhs.clone();
    let start = Instant::now();
    let result = black_box(match route {
        "monomial" => black_box(&source).replace_with_poly(0, black_box(&rhs)),
        "horner" => grouped_horner(black_box(&source), 0, black_box(&rhs)),
        _ => unreachable!(),
    });
    let elapsed_ns = start.elapsed().as_nanos();
    // An independently factored closed form is built outside the clock. No
    // Atom conversion, callback or tensor initialization enters either route.
    let mut expected_powers = source.zero();
    for d in 1..=degree {
        expected_powers = expected_powers + rhs.pow(usize::from(d));
    }
    assert_eq!(result, &coefficient * &expected_powers);
    assert_eq!(source, source_before);
    assert_eq!(rhs, rhs_before);
    println!(
        "{{\"route\":\"{route}\",\"coefficient_terms\":{n},\"degree\":{degree},\"variables\":3,\"input_terms\":{},\"output_terms\":{},\"elapsed_ns\":{elapsed_ns},\"exact\":true,\"inputs_unchanged\":true}}",
        source.nterms(),
        result.nterms()
    );
}
