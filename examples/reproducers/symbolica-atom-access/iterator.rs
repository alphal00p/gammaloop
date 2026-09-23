//! Reproduce second-argument access costs with Symbolica 3.0.0 and std only.
//! Parsing and equality checks are outside timing. Both readers receive the
//! same opaque, copied iterator and their returned AtomViews are consumed equally.
use std::{hint::black_box, time::Instant};
use symbolica::atom::{Atom, AtomView, representation::ListIterator};

#[inline]
fn nth(mut arguments: ListIterator<'_>) -> AtomView<'_> {
    arguments.nth(1).unwrap()
}

#[inline]
fn next_twice(mut arguments: ListIterator<'_>) -> AtomView<'_> {
    arguments.next();
    arguments.next().unwrap()
}

#[inline(never)]
fn nth_boundary(arguments: ListIterator<'_>) -> AtomView<'_> {
    nth(arguments)
}

#[inline(never)]
fn next_boundary(arguments: ListIterator<'_>) -> AtomView<'_> {
    next_twice(arguments)
}

#[derive(Clone, Copy)]
enum Method {
    Nth,
    Next,
    NthBoundary,
    NextBoundary,
}
impl Method {
    const ALL: [Self; 4] = [Self::Nth, Self::Next, Self::NthBoundary, Self::NextBoundary];
    fn name(self) -> &'static str {
        match self {
            Self::Nth => "nth",
            Self::Next => "next",
            Self::NthBoundary => "nth_boundary",
            Self::NextBoundary => "next_boundary",
        }
    }
    fn measure(self, arguments: ListIterator<'_>, iterations: usize) -> f64 {
        let start = Instant::now();
        match self {
            Self::Nth => {
                for _ in 0..iterations {
                    black_box(nth(black_box(arguments)));
                }
            }
            Self::Next => {
                for _ in 0..iterations {
                    black_box(next_twice(black_box(arguments)));
                }
            }
            Self::NthBoundary => {
                for _ in 0..iterations {
                    black_box(nth_boundary(black_box(arguments)));
                }
            }
            Self::NextBoundary => {
                for _ in 0..iterations {
                    black_box(next_boundary(black_box(arguments)));
                }
            }
        }
        start.elapsed().as_secs_f64() * 1e9 / iterations as f64
    }
}

fn main() {
    let rounds: usize = std::env::args().nth(1).map_or(9, |v| v.parse().unwrap());
    for source in [
        "f(4,mu)",
        "f(D,mu)",
        "f(1/3,mu)",
        "f(g(a,b,c),mu)",
        "f(a+b+c,mu)",
    ] {
        let input =
            Atom::parse(source, "mre", symbolica::parser::ParseSettings::symbolica()).unwrap();
        let AtomView::Fun(function) = input.as_view() else {
            panic!("expected a function")
        };
        let arguments = function.iter();
        let expected = function.get(1);
        assert_eq!(nth(arguments), expected);
        assert_eq!(next_twice(arguments), expected);
        assert_eq!(nth_boundary(arguments), expected);
        assert_eq!(next_boundary(arguments), expected);
        let mut iterations = [0usize; 4];
        for (m, method) in Method::ALL.into_iter().enumerate() {
            method.measure(arguments, 20_000);
            let ns = method.measure(arguments, 100_000);
            iterations[m] = (20_000_000.0 / ns).clamp(10_000.0, 5_000_000.0) as usize;
        }
        for round in 0..rounds {
            for step in 0..4 {
                let m = (round + step) % 4;
                let method = Method::ALL[m];
                let ns = method.measure(arguments, iterations[m]);
                println!(
                    "{{\"input\":\"{source}\",\"method\":\"{}\",\"round\":{round},\"iterations\":{},\"ns\":{ns}}}",
                    method.name(),
                    iterations[m]
                );
            }
        }
    }
}
