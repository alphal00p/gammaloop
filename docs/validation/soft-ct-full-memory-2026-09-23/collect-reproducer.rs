use std::time::Instant;
use symbolica::{
    atom::{Atom, AtomCore},
    function, symbol,
};

fn main() {
    let n: usize = std::env::args().nth(1).unwrap().parse().unwrap();
    let a = symbol!("memory_probe::a");
    let b = symbol!("memory_probe::b");
    let c = symbol!("memory_probe::c");
    let expression =
        Atom::add_many((0..n).map(|i| (function!(a, i) + function!(b, i)) * function!(c, i)));
    let input_bytes = expression.as_view().get_byte_size();
    let started = Instant::now();
    let result = expression.collect_factors();
    let elapsed = started.elapsed().as_secs_f64();
    let status = std::fs::read_to_string("/proc/self/status").unwrap();
    let peak_kib = status
        .lines()
        .find(|line| line.starts_with("VmHWM:"))
        .unwrap()
        .split_whitespace()
        .nth(1)
        .unwrap();
    println!(
        "{{\"terms\":{n},\"input_bytes\":{input_bytes},\"result_bytes\":{},\"unchanged\":{},\"elapsed_seconds\":{elapsed},\"peak_rss_kib\":{peak_kib}}}",
        result.as_view().get_byte_size(),
        result == expression
    );
}
