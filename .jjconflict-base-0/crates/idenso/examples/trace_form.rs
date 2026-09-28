//! Reproducible native trace comparison: cargo run -p idenso --profile dev-optim
//! --example trace_form -- /absolute/path/to/form [maximum_even_length]
//! Idenso times include typed admission, gamma simplification, alias resolution,
//! and result destruction; FORM times include process startup, parsing,
//! sorting and the exact scalar check. Neither includes input construction.
//! Term counts describe the output representation, not algebraic correctness.

use std::{hint::black_box, process::Command, time::Instant};

use idenso::{gamma, representations::Bispinor};
use spenso::{
    network::library::symbolic::ETS,
    p, q,
    structure::representation::{Minkowski, RepName},
    trace,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    symbol,
};
use symbolica_utils::AtomPrintExt;

fn median(mut samples: Vec<f64>) -> f64 {
    samples.sort_by(f64::total_cmp);
    samples[samples.len() / 2]
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let form = args
        .get(1)
        .expect("usage: trace_form /path/to/form [maximum_even_length]");
    let maximum: usize = args.get(2).map_or(14, |n| n.parse().unwrap());
    assert!((2..=14).contains(&maximum));
    idenso::representations::initialize();
    let mink = Minkowski {}.new_rep(4);
    let spin = Bispinor {}.new_rep(4).to_symbolic([]);
    let p = p!(mink.to_symbolic([]));
    let q = q!(mink.to_symbolic([]));
    let [pp, qq, pq] =
        ["pp", "qq", "pq"].map(|name| Atom::var(symbol!(&format!("trace_benchmark::{name}"))));
    let version = Command::new(form).arg("-v").output().expect("launch FORM");
    assert!(version.status.success());
    eprintln!("{}", String::from_utf8_lossy(&version.stdout).trim());
    eprintln!(
        "Idenso: profile selected by cargo; cold includes lazy kernel generation; warm median of 5. FORM: median of 5 processes after one warm-up, including scalar verification."
    );
    let directory = std::env::temp_dir().join(format!("idenso-form-{}", std::process::id()));
    std::fs::create_dir(&directory).unwrap();
    eprintln!("FORM programs and output: {}", directory.display());
    println!(
        "case,length,terms,idenso_cold_ms,idenso_warm_ms,form_trace4_process_ms,form_tracen_process_ms,form_trace4_terms,form_tracen_terms"
    );
    for n in (2..=maximum).step_by(2) {
        for case in ["free", "paired", "alternating"] {
            let is_p = |i: usize| {
                if case == "alternating" {
                    i.is_multiple_of(2)
                } else {
                    i < 2
                }
            };
            let indices: Vec<_> = (0..n)
                .map(|i| {
                    if case == "free" {
                        mink.pattern(symbol!(&format!("trace_benchmark::mu{i}")))
                    } else if is_p(i) {
                        p.clone()
                    } else {
                        q.clone()
                    }
                })
                .collect();
            let input = trace!(&spin; indices.iter().map(|mu| gamma!(mu)));
            let start = Instant::now();
            let result = idenso::tensor::SymbolicTensor::infer(
                (black_box(&input)).as_atom_view().to_owned(),
            )
            .unwrap()
            .simplify_gamma(idenso::dirac::GammaSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression();
            let cold = start.elapsed().as_secs_f64();
            let mut times = Vec::new();
            for _ in 0..5 {
                let start = Instant::now();
                let _ = black_box(
                    idenso::tensor::SymbolicTensor::infer(
                        (black_box(&input)).as_atom_view().to_owned(),
                    )
                    .unwrap()
                    .simplify_gamma(idenso::dirac::GammaSimplifySettings::default())
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression(),
                );
                times.push(start.elapsed().as_secs_f64());
            }
            // Verify the full scalar polynomial in p^2, q^2 and p.q,
            // outside the Idenso timer. No on-shell or numeric specialization.
            let expected = if case == "alternating" {
                let (mut previous, mut current) = (Atom::num(4), Atom::num(4) * &pq);
                for _ in 2..=n / 2 {
                    (previous, current) = (
                        current.clone(),
                        (Atom::num(2) * &pq * &current - &pp * &qq * &previous).expand(),
                    );
                }
                current
            } else {
                Atom::num(4) * &pp * Atom::mul_many(std::iter::repeat_n(&qq, n / 2 - 1))
            };
            let projected = result.replace_map(|a, _, out| {
                if let AtomView::Fun(f) = a
                    && f.get_symbol() == ETS.metric
                {
                    let pair: Vec<_> = f
                        .iter()
                        .map(|mu| {
                            indices
                                .iter()
                                .position(|index| index.as_view() == mu)
                                .expect("known trace argument")
                        })
                        .collect();
                    let value = match (is_p(pair[0]), is_p(pair[1])) {
                        (true, true) => &pp,
                        (false, false) => &qq,
                        _ => &pq,
                    };
                    **out = value.clone();
                }
            });
            assert!((&projected - &expected).expand().is_zero(), "{case}/{n}");
            let expected = expected.to_bare_ordered_string();
            let word = (0..n)
                .map(|i| {
                    if case == "free" {
                        format!("mu{i}")
                    } else if is_p(i) {
                        "p".into()
                    } else {
                        "q".into()
                    }
                })
                .collect::<Vec<_>>()
                .join(",");
            let declarations = (0..n)
                .map(|i| format!("mu{i}"))
                .collect::<Vec<_>>()
                .join(",");
            let projection = if case == "free" {
                format!(
                    "Multiply {};\n.sort\n",
                    (0..n)
                        .map(|i| format!("{}(mu{i})", if is_p(i) { "p" } else { "q" }))
                        .collect::<Vec<_>>()
                        .join("*")
                )
            } else {
                String::new()
            };
            let mut form_times = Vec::new();
            let mut form_counts = Vec::new();
            for mode in ["trace4", "tracen"] {
                let source = format!(
                    "Off Statistics;\nVectors p,q;\nSymbols pp,qq,pq;\nIndices {declarations};\nLocal F=g_(1,{word});\n{mode},1;\n.sort\n#$nterms=termsin_(F);\n#write \"NTERMS=%$\",$nterms\n{projection}id p.p=pp;\nid q.q=qq;\nid p.q=pq;\n.sort\nLocal Check=F-({expected});\n.sort\n#write \"CHECK=%E\",Check\n.end\n"
                );
                let path = directory.join(format!("{case}-{n}-{mode}.frm"));
                std::fs::write(&path, source).unwrap();
                let mut samples = Vec::new();
                for round in 0..6 {
                    let start = Instant::now();
                    let output = Command::new(form)
                        .arg("-q")
                        .arg(&path)
                        .current_dir(&directory)
                        .output()
                        .unwrap();
                    let seconds = start.elapsed().as_secs_f64();
                    assert!(
                        output.status.success(),
                        "{}",
                        String::from_utf8_lossy(&output.stderr)
                    );
                    let stdout = String::from_utf8(output.stdout).unwrap();
                    assert!(
                        stdout.lines().any(|line| line
                            .strip_prefix("CHECK=")
                            .is_some_and(|value| value.trim().trim_end_matches(';') == "0")),
                        "{stdout}"
                    );
                    if round == 0 {
                        let count: usize = stdout
                            .lines()
                            .find_map(|line| line.strip_prefix("NTERMS="))
                            .unwrap()
                            .trim()
                            .parse()
                            .unwrap();
                        form_counts.push(count);
                    }
                    if round > 0 {
                        samples.push(seconds);
                    }
                    std::fs::write(path.with_extension("log"), stdout).unwrap();
                }
                form_times.push(median(samples));
            }
            println!(
                "{case},{n},{},{:.6},{:.6},{:.6},{:.6},{},{}",
                result.expand().nterms(),
                cold * 1000.,
                median(times) * 1000.,
                form_times[0] * 1000.,
                form_times[1] * 1000.,
                form_counts[0],
                form_counts[1]
            );
        }
    }
}
