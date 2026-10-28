//! Compare Symbolica scalar expansion APIs without tensor libraries.
//! Usage: scalar_expansion_routes ROUTE --synthetic N | ROUTE --input FILE
//! Routes: ordinary, via_poly, whole_poly, balanced. One timed call per process.
use std::{hint::black_box, time::Instant};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    domains::rational::{Q, RationalField},
    parser::ParseSettings,
    poly::polynomial::MultivariatePolynomial,
};

type Poly = MultivariatePolynomial<RationalField, i16>;

fn convert(value: AtomView<'_>) -> Poly {
    value.try_to_polynomial(&Q, None).unwrap()
}

fn join(mut left: Poly, mut right: Poly) -> Poly {
    left.unify_variables(&mut right);
    &left + &right
}

fn balanced(value: AtomView<'_>) -> Poly {
    let AtomView::Add(sum) = value else {
        return convert(value);
    };
    // Binary reduction of the existing root-Add sequence. Nested conversion
    // and polynomial addition still use the ordinary public APIs.
    let mut levels: Vec<Option<Poly>> = Vec::new();
    for term in sum {
        let mut carry = convert(term);
        let mut level = 0;
        loop {
            if level == levels.len() {
                levels.push(Some(carry));
                break;
            }
            if let Some(left) = levels[level].take() {
                carry = join(left, carry);
                level += 1;
            } else {
                levels[level] = Some(carry);
                break;
            }
        }
    }
    let mut result = None;
    for polynomial in levels.into_iter().rev().flatten() {
        result = Some(match result {
            None => polynomial,
            Some(left) => join(left, polynomial),
        });
    }
    result.unwrap_or_else(|| convert(Atom::Zero.as_view()))
}

fn synthetic(n: usize) -> String {
    const VARIABLES: usize = 20;
    // There are binomial(20 + 5 - 1, 5) = 42,504 degree-five monomials.
    assert!((VARIABLES..=42_504).contains(&n), "N must be 20..=42504");
    // Seed every variable, so even the smallest admitted N has the same basis.
    let mut monomials: Vec<_> = (0..VARIABLES).map(|i| format!("x{i}^5")).collect();
    'enumerate: for a in 0..VARIABLES {
        for b in a..VARIABLES {
            for c in b..VARIABLES {
                for d in c..VARIABLES {
                    for e in d..VARIABLES {
                        if monomials.len() == n {
                            break 'enumerate;
                        }
                        if a == e {
                            continue; // Already seeded pure fifth power.
                        }
                        monomials.push([a, b, c, d, e].map(|i| format!("x{i}")).join("*"));
                    }
                }
            }
        }
    }
    assert_eq!(monomials.len(), n);
    monomials
        .iter()
        .enumerate()
        .map(|(i, monomial)| format!("{}/{}*({monomial})*(x0+x1)*(x2-x3)", i % 7 + 1, i % 5 + 1))
        .collect::<Vec<_>>()
        .join("+")
}

fn terms(value: &Atom) -> usize {
    match value.as_view() {
        AtomView::Add(sum) => sum.get_nargs(),
        _ => 1,
    }
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(
        args.len(),
        4,
        "usage: scalar_expansion_routes ordinary|via_poly|whole_poly|balanced --synthetic N|--input FILE"
    );
    let route = args[1].as_str();
    assert!(matches!(
        route,
        "ordinary" | "via_poly" | "whole_poly" | "balanced"
    ));
    let (text, requested_terms) = match args[2].as_str() {
        "--synthetic" => {
            let n = args[3].parse::<usize>().expect("N must be an integer");
            (synthetic(n), Some(n))
        }
        "--input" => (
            std::fs::read_to_string(&args[3]).expect("read input file"),
            None,
        ),
        _ => panic!("use --synthetic N or --input FILE"),
    };
    // Parsing is outside all operation clocks. There are deliberately no
    // tensor tags, symmetry declarations, callbacks or metadata registrations.
    let source = Atom::parse(&text, "scalar_expansion_mre", ParseSettings::symbolica()).unwrap();
    let start = Instant::now();
    let (result, conversion_ns, emission_ns, polynomial_drop_ns) = match route {
        "ordinary" => (black_box(&source).expand(), None, None, None),
        "via_poly" => (
            black_box(&source).expand_via_poly::<i16, Atom>(None),
            None,
            None,
            None,
        ),
        "whole_poly" | "balanced" => {
            let polynomial = if route == "whole_poly" {
                convert(black_box(source.as_view()))
            } else {
                balanced(black_box(source.as_view()))
            };
            let converted = start.elapsed().as_nanos();
            let emission_start = Instant::now();
            let expression = polynomial.to_expression();
            let emitted = emission_start.elapsed().as_nanos();
            let drop_start = Instant::now();
            drop(polynomial);
            (
                expression,
                Some(converted),
                Some(emitted),
                Some(drop_start.elapsed().as_nanos()),
            )
        }
        _ => unreachable!(),
    };
    let result = black_box(result);
    let total_ns = start.elapsed().as_nanos();
    // This is a pure-parser oracle, not a tensor/FORM certificate. All checks,
    // term counting and output destruction are outside the operation clock.
    let expected = source.expand();
    assert_eq!(result, expected, "ordinary expansion mismatch: {route}");
    assert_eq!(result.expand(), result, "scalar expansion fixed point");
    let ns = |n: Option<u128>| n.map_or_else(|| "null".to_string(), |n| n.to_string());
    println!(
        "{{\"route\":\"{route}\",\"input_kind\":\"{}\",\"synthetic_terms\":{},\"synthetic_variables\":{},\"input_top_terms\":{},\"input_bytes\":{},\"output_terms\":{},\"exact_ordinary\":true,\"fixedpoint\":true,\"total_ns\":{total_ns},\"conversion_ns\":{},\"emission_ns\":{},\"polynomial_drop_ns\":{}}}",
        if requested_terms.is_some() {
            "synthetic"
        } else {
            "file"
        },
        ns(requested_terms.map(|n| n as u128)),
        ns(requested_terms.map(|_| 20)),
        terms(&source),
        source.as_view().get_byte_size(),
        terms(&result),
        ns(conversion_ns),
        ns(emission_ns),
        ns(polynomial_drop_ns),
    );
}
