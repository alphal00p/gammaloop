//! Public-API MRE: both operations return exactly (function head ID, arity).
//! rustc --edition=2024 -C opt-level=2 header.rs --extern symbolica=... -L dependency=...
//! Pin the produced executable to a CPU when comparing timings.
use std::{hint::black_box, time::Instant};
use symbolica::atom::{Atom, AtomView, representation::FunView};

#[inline(always)]
fn through_iterator(fun: FunView<'_>) -> (u32, usize) {
    (fun.get_symbol_id(), fun.iter().len())
}
#[inline(always)]
fn through_accessors(fun: FunView<'_>) -> (u32, usize) {
    (fun.get_symbol_id(), fun.get_nargs())
}
fn batch<const ITER: bool>(functions: &[FunView<'_>], repetitions: usize) -> f64 {
    let start = Instant::now();
    let mut checksum = 0u64;
    for _ in 0..repetitions {
        for &function in functions {
            let function = black_box(function);
            let (id, n) = if ITER {
                through_iterator(function)
            } else {
                through_accessors(function)
            };
            checksum = checksum.wrapping_add(id as u64 + n as u64);
        }
    }
    black_box(checksum);
    start.elapsed().as_secs_f64() * 1e9 / (repetitions * functions.len()) as f64
}
fn main() {
    let sources = [
        ("single_short_header", vec!["f(4,mu)".to_string()]),
        (
            "mixed_heads_and_arities",
            (0..128)
                .map(|i| {
                    format!(
                        "f{i}({})",
                        (0..i % 9 + 1)
                            .map(|j| format!("x{j}"))
                            .collect::<Vec<_>>()
                            .join(",")
                    )
                })
                .collect(),
        ),
        (
            "wide_header",
            vec![format!(
                "fwide({})",
                (0..300)
                    .map(|i| format!("z{i}"))
                    .collect::<Vec<_>>()
                    .join(",")
            )],
        ),
        ("wide_head_id", vec!["last(4,mu)".to_string()]),
    ];
    for (name, sources) in sources {
        let atoms: Vec<_> = sources
            .iter()
            .map(|s| {
                Atom::parse(
                    s,
                    "header_mre",
                    symbolica::parser::ParseSettings::symbolica(),
                )
                .unwrap()
            })
            .collect();
        let functions: Vec<_> = atoms
            .iter()
            .map(|a| match a.as_view() {
                AtomView::Fun(f) => f,
                _ => unreachable!(),
            })
            .collect();
        for &f in &functions {
            assert_eq!(through_iterator(f), through_accessors(f));
        }
        let repetitions = 10_000_000 / functions.len();
        batch::<true>(&functions, repetitions / 10);
        batch::<false>(&functions, repetitions / 10);
        let mut iterator = Vec::new();
        let mut accessors = Vec::new();
        for round in 0..9 {
            if round % 2 == 0 {
                iterator.push(batch::<true>(&functions, repetitions));
                accessors.push(batch::<false>(&functions, repetitions));
            } else {
                accessors.push(batch::<false>(&functions, repetitions));
                iterator.push(batch::<true>(&functions, repetitions));
            }
        }
        let mut sorted_iterator = iterator.clone();
        sorted_iterator.sort_by(f64::total_cmp);
        let mut sorted_accessors = accessors.clone();
        sorted_accessors.sort_by(f64::total_cmp);
        println!(
            "{{\"case\":\"{name}\",\"functions\":{},\"first_id\":{},\"first_arity\":{},\"repetitions\":{repetitions},\"iterator_samples_ns\":{iterator:?},\"accessor_samples_ns\":{accessors:?},\"iterator_median_ns\":{},\"accessor_median_ns\":{}}}",
            functions.len(),
            functions[0].get_symbol_id(),
            functions[0].get_nargs(),
            sorted_iterator[4],
            sorted_accessors[4]
        );
    }
}
