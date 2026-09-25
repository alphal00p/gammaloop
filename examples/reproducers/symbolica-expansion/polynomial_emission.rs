//! Isolate Symbolica's polynomial-to-expression conversion, without tensor libraries.
//! Run: cargo run --release --bin polynomial_emission -- [N=12] [tensor|shallow|variables] [samples=5] [roundtrip]
//! Only `profile_run` performs the measured conversion; setup and validation are outside it.
//! Inverse polynomial roundtrips run for N<=10; `roundtrip` forces this slow check for larger N.
use std::{hint::black_box, time::Instant};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder},
    domains::{
        Ring,
        integer::{Integer, Z},
    },
    poly::{CoefficientToExpression, PolyVariable, polynomial::MultivariatePolynomial},
    symbol,
};

// Enumerate signed perfect matchings: pair the first remaining index with each
// possible partner, alternating signs. Each leaf contributes one squarefree term.
fn pairings(
    remaining: u16,
    negative: bool,
    variables: &[Vec<usize>],
    powers: &mut [u8],
    coefficients: &mut Vec<Integer>,
    exponents: &mut Vec<u8>,
) {
    if remaining == 0 {
        coefficients.push(Integer::from(if negative { -4 } else { 4 }));
        exponents.extend_from_slice(powers);
        return;
    }
    let first = remaining.trailing_zeros() as usize;
    let rest = remaining & !(1 << first);
    let mut candidates = rest;
    let mut negative = negative;
    while candidates != 0 {
        let partner = candidates.trailing_zeros() as usize;
        candidates &= candidates - 1;
        let variable = variables[first][partner];
        powers[variable] = 1;
        pairings(
            rest & !(1 << partner),
            negative,
            variables,
            powers,
            coefficients,
            exponents,
        );
        powers[variable] = 0;
        negative = !negative;
    }
}

#[inline(never)]
fn profile_run<R: Ring>(polynomial: &MultivariatePolynomial<R, u8>) -> Atom
where
    R::Element: CoefficientToExpression<R>,
{
    polynomial.to_expression()
}

fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert!(
        args.len() <= 4,
        "usage: polynomial_emission [N=12] [tensor|shallow|variables] [samples=5] [roundtrip]"
    );
    let n: usize = args
        .first()
        .map_or(12, |s| s.parse().expect("N must be an integer"));
    let payload = args.get(1).map(String::as_str).unwrap_or("tensor");
    let samples: usize = args
        .get(2)
        .map_or(5, |s| s.parse().expect("samples must be an integer"));
    assert!(
        n.is_multiple_of(2) && (2..=14).contains(&n),
        "N must be even, from 2 to 14"
    );
    assert!(
        (1..=1000).contains(&samples),
        "samples must be from 1 to 1000"
    );
    assert!(matches!(payload, "tensor" | "shallow" | "variables"));
    assert!(args.get(3).is_none_or(|arg| arg == "roundtrip"));
    let roundtrip_checked = n <= 10 || args.get(3).is_some();

    // These functions are opaque commuting polynomial variables, with no custom
    // normalization callbacks. Their payload changes; the polynomial stays fixed.
    let g = symbol!("emission_mre::g"; Symmetric);
    let mink = symbol!("emission_mre::mink");
    let idx = symbol!("emission_mre::idx");
    let endpoints: Vec<_> = (0..n)
        .map(|i| Atom::var(symbol!(&format!("emission_mre::mu{i}"))))
        .collect();
    let slots: Vec<_> = endpoints
        .iter()
        .map(|endpoint| {
            let index = FunctionBuilder::new(idx)
                .add_arg(endpoint)
                .add_arg(0)
                .finish();
            FunctionBuilder::new(mink)
                .add_arg(4)
                .add_arg(index)
                .finish()
        })
        .collect();
    let variable_count = n * (n - 1) / 2;
    let mut variables: Vec<PolyVariable> = Vec::with_capacity(variable_count);
    let mut variable_ids = vec![vec![0; n]; n];
    for i in 0..n {
        for j in i + 1..n {
            let variable = match payload {
                "tensor" => FunctionBuilder::new(g)
                    .add_arg(&slots[i])
                    .add_arg(&slots[j])
                    .finish(),
                "shallow" => FunctionBuilder::new(g)
                    .add_arg(&endpoints[i])
                    .add_arg(&endpoints[j])
                    .finish(),
                "variables" => Atom::var(symbol!(&format!("emission_mre::m_{i}_{j}"))),
                _ => unreachable!(),
            };
            variable_ids[i][j] = variables.len();
            variables.push(variable.try_into().unwrap());
        }
    }
    let terms = (1..n).step_by(2).product::<usize>();
    let mut coefficients = Vec::with_capacity(terms);
    let mut exponents = Vec::with_capacity(terms * variable_count);
    pairings(
        (1 << n) - 1,
        false,
        &variable_ids,
        &mut vec![0; variable_count],
        &mut coefficients,
        &mut exponents,
    );
    let polynomial = MultivariatePolynomial::<_, u8>::from_coefficient_list(
        coefficients,
        exponents,
        variables.into(),
        &Z,
    );
    assert_eq!(polynomial.nterms(), terms);

    // Warm caches before any measured sample. The inverse conversion is much
    // slower at large N, so make that semantic check opt-in above N=10.
    let expected = polynomial.to_expression();
    assert_eq!(expected.nterms(), terms);
    if roundtrip_checked {
        let roundtrip = expected.to_polynomial::<_, u8>(&Z, polynomial.get_vars());
        assert_eq!(roundtrip, polynomial, "exact polynomial roundtrip failed");
    }
    let roundtrip_equal = if roundtrip_checked { "true" } else { "null" };
    for sample in 0..samples {
        let start = Instant::now();
        let output = black_box(profile_run(black_box(&polynomial)));
        let emission_ns = start.elapsed().as_nanos();
        assert_eq!(output, expected); // Compare with warmup outside both clocks.
        let start = Instant::now();
        drop(output);
        let drop_ns = start.elapsed().as_nanos();
        let emission_and_drop_ns = emission_ns + drop_ns;
        println!(
            "{{\"n\":{n},\"payload\":\"{payload}\",\"sample\":{sample},\"variables\":{variable_count},\"terms\":{terms},\"emission_ns\":{emission_ns},\"drop_ns\":{drop_ns},\"emission_and_drop_ns\":{emission_and_drop_ns},\"exact_equal\":true,\"roundtrip_checked\":{roundtrip_checked},\"roundtrip_equal\":{roundtrip_equal}}}"
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::domains::float::{Float, FloatField};

    #[test]
    fn compound_variables_preserve_rounded_coefficient_multiplication_order() {
        let field = FloatField::from_rep(Float::with_val(53, 0));
        let x = Atom::var(symbol!("emission_boundary::x"));
        let y = Atom::var(symbol!("emission_boundary::y"));
        for [a, b, c] in [[0.1, 0.2, 0.3], [1.0e16, 1.0e-16, 0.3], [1.1, 1.3, 1.7]] {
            let mut variables = [
                Atom::num(Float::with_val(53, a)) * &x,
                Atom::num(Float::with_val(53, b)) * &y,
            ];
            variables.sort();
            let coefficient = Float::with_val(53, c);
            // The general emitter appends its coefficient after variable factors.
            // Moving it first changes rounded arithmetic, even for one monomial.
            let expected = (&variables[0] * &variables[1]) * Atom::num(coefficient.clone());
            let polynomial = MultivariatePolynomial::<_, u8>::from_coefficient_list(
                vec![coefficient],
                vec![1, 1],
                vec![
                    PolyVariable::Function(
                        symbol!("emission_boundary::map_a"),
                        variables[0].clone(),
                    ),
                    PolyVariable::Function(
                        symbol!("emission_boundary::map_b"),
                        variables[1].clone(),
                    ),
                ]
                .into(),
                &field,
            );
            assert_eq!(polynomial.to_expression(), expected, "{a} * {b} * {c}");
        }
    }
}
