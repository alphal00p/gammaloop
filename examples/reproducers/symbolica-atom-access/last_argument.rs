//! MRE: locating the end of a function's final argument can be unnecessary.
//!
//! Dependency: symbolica = "=3.0.0". No Spenso, private encoding constants, or
//! cached argument boundaries. Both routines check arity and decode argument 0.
//! The span variant borrows the rest of the exact function bytes as argument 1.
//! AtomView::from requires valid complete atom bytes; that follows here from a
//! valid FunView with exactly two arguments and the first argument's exact span.
//!
//! Run optimized on one CPU. Output is CSV: seven alternating samples/case.
//! Construction, validation, and iteration calibration are outside samples.
//! Optional argument: target batch duration in milliseconds (default 30).
use std::{
    hint::black_box,
    time::{Duration, Instant},
};
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, representation::FunView};

#[inline(never)]
fn second_iterator(fun: FunView<'_>) -> Option<AtomView<'_>> {
    let mut args = fun.iter();
    if args.len() != 2 {
        return None;
    }
    args.next()?;
    args.next()
}

#[inline(never)]
fn second_remaining_span(fun: FunView<'_>) -> Option<AtomView<'_>> {
    let mut args = fun.iter();
    if args.len() != 2 {
        return None;
    }
    let first = args.next()?.get_data();
    let whole = fun.as_view().get_data();
    // Both slices belong to this same valid function's immutable allocation.
    let offset = first.as_ptr() as usize - whole.as_ptr() as usize + first.len();
    Some(AtomView::from(&whole[offset..]))
}

fn parse(source: &str) -> Atom {
    Atom::parse(
        source,
        "last_argument_mre",
        symbolica::parser::ParseSettings::symbolica(),
    )
    .unwrap()
}
fn as_function(atom: &Atom) -> FunView<'_> {
    let AtomView::Fun(fun) = atom.as_view() else {
        panic!("expected function")
    };
    fun
}
fn batch(fun: FunView<'_>, repetitions: usize, method: bool) -> Duration {
    let start = Instant::now();
    if method {
        for _ in 0..repetitions {
            black_box(second_remaining_span(black_box(fun)));
        }
    } else {
        for _ in 0..repetitions {
            black_box(second_iterator(black_box(fun)));
        }
    }
    start.elapsed()
}
fn main() {
    let target_ms = std::env::args()
        .nth(1)
        .map(|s| s.parse::<u64>().unwrap())
        .unwrap_or(30);
    let target = Duration::from_millis(target_ms.max(1));
    let holder = parse("f(4,mu)");
    let head = as_function(&holder).get_symbol();
    let mut cases = vec![("symbol".to_string(), holder)];
    for depth in [1usize, 4, 16, 64] {
        let mut power = parse("z");
        for _ in 0..depth {
            power = parse("x").pow(power);
        }
        let atom = FunctionBuilder::new(head)
            .add_arg(4)
            .add_arg(power)
            .finish();
        cases.push((format!("power_depth_{depth}"), atom));
    }
    for source in [
        "f()",
        "f(4)",
        "f(4,mu,extra)",
        "f(4,1/256)",
        "f(4,1+2*𝑖)",
        "f(4,g(a,b,c))",
        "f(4,x+y)",
        "f(4,x*y)",
        "f(g(a,b),x^(y^z))",
        "f(4,184467440737095516160000000001)",
    ] {
        let atom = parse(source);
        assert_eq!(
            second_iterator(as_function(&atom)),
            second_remaining_span(as_function(&atom)),
            "{source}"
        );
    }
    for (_, atom) in &cases {
        let fun = as_function(atom);
        assert_eq!(second_iterator(fun), second_remaining_span(fun));
    }
    println!("case,method,sample,repetitions,ns_per_call,last_arg_bytes");
    for (name, atom) in &cases {
        let fun = as_function(atom);
        let mut repetitions = 1usize;
        while batch(fun, repetitions, false) < target {
            repetitions *= 2;
        }
        let bytes = second_iterator(fun).unwrap().get_data().len();
        for sample in 0..7 {
            for method in [sample % 2 == 1, sample % 2 == 0] {
                let elapsed = batch(fun, repetitions, method);
                println!(
                    "{name},{},{sample},{repetitions},{:.3},{bytes}",
                    if method { "remaining_span" } else { "iterator" },
                    elapsed.as_secs_f64() * 1e9 / repetitions as f64
                );
            }
        }
    }
}
