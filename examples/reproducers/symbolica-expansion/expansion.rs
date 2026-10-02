//! Expand an ordinary Clifford pairing polynomial using only Symbolica.
//! Metrics are independent commuting indeterminates; only their payload changes.
use std::{collections::HashMap, hint::black_box, time::Instant};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder},
    symbol,
};

fn pairing(
    word: &[usize],
    metrics: &[Atom],
    n: usize,
    memo: &mut HashMap<Vec<usize>, Atom>,
) -> Atom {
    if word.is_empty() {
        return Atom::num(4);
    }
    if let Some(result) = memo.get(word) {
        return result.clone();
    }
    let mut terms = Vec::with_capacity(word.len() - 1);
    for j in 1..word.len() {
        let rest: Vec<_> = word[1..j].iter().chain(&word[j + 1..]).copied().collect();
        let mut term = &metrics[word[0] * n + word[j]] * pairing(&rest, metrics, n, memo);
        if j % 2 == 0 {
            term = -term;
        }
        terms.push(term);
    }
    let result = Atom::add_many(terms);
    memo.insert(word.to_vec(), result.clone());
    result
}

#[inline(never)]
fn profile_run(input: &Atom, polynomial: bool) -> Atom {
    if polynomial {
        input.expand_via_poly::<u8, Atom>(None)
    } else {
        input.expand()
    }
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let n = args
        .get(1)
        .expect("usage: expansion N [tensor|shallow|variables] [expand|poly]")
        .parse::<usize>()
        .expect("N must be an even integer from 2 to 14");
    assert!(n.is_multiple_of(2) && (2..=14).contains(&n));
    let payload = args.get(2).map(String::as_str).unwrap_or("tensor");
    let method = args.get(3).map(String::as_str).unwrap_or("expand");
    assert!(matches!(method, "expand" | "poly"));
    let g = symbol!("expansion_mre::g"; Symmetric);
    let mink = symbol!("expansion_mre::mink");
    let dimension = Atom::var(symbol!("expansion_mre::D"));
    let endpoints: Vec<_> = (0..n)
        .map(|i| Atom::var(symbol!(&format!("expansion_mre::mu{i}"))))
        .collect();
    let slots: Vec<_> = endpoints
        .iter()
        .map(|a| {
            FunctionBuilder::new(mink)
                .add_arg(&dimension)
                .add_arg(a)
                .finish()
        })
        .collect();
    let metrics: Vec<_> = (0..n * n)
        .map(|k| {
            let (i, j) = (k / n, k % n);
            match payload {
                "tensor" => FunctionBuilder::new(g)
                    .add_arg(&slots[i])
                    .add_arg(&slots[j])
                    .finish(),
                "shallow" => FunctionBuilder::new(g)
                    .add_arg(&endpoints[i])
                    .add_arg(&endpoints[j])
                    .finish(),
                "variables" => Atom::var(symbol!(&format!(
                    "expansion_mre::m_{}_{}",
                    i.min(j),
                    i.max(j)
                ))),
                _ => panic!("use tensor, shallow or variables"),
            }
        })
        .collect();
    let input = pairing(
        &(0..n).collect::<Vec<_>>(),
        &metrics,
        n,
        &mut HashMap::new(),
    );
    let expected_terms = (1..n).step_by(2).product::<usize>();
    let expected = input.expand();
    assert_eq!(expected.nterms(), expected_terms);
    for sample in 0..5 {
        let start = Instant::now();
        let output = black_box(profile_run(black_box(&input), method == "poly"));
        let expansion_ns = start.elapsed().as_nanos();
        assert_eq!(output, expected);
        let start = Instant::now();
        drop(output);
        let drop_ns = start.elapsed().as_nanos();
        println!(
            "{{\"n\":{n},\"payload\":\"{payload}\",\"method\":\"{method}\",\"sample\":{sample},\"expansion_ns\":{expansion_ns},\"drop_ns\":{drop_ns},\"terms\":{expected_terms},\"exact_equal\":true}}"
        );
    }
}
