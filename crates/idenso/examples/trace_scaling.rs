//! Amortized native comparison, with raw per-call timings in CSV:
//! cargo run -p idenso --profile dev-optim --example trace_scaling -- /path/to/form
//! [maximum_even_length] > trace_scaling.csv
//!
//! Idenso measures typed admission, gamma simplification, alias resolution, and
//! result destruction in each warm loop.
//! FORM wall times divide complete processes by all traces, including amortized
//! startup, parsing, sorting and disposal. Its separate CPU timer covers only
//! trace4 and sorting. Input construction and correctness checks are excluded
//! from Idenso timing; this diagnostic performs no independent algebra check.
//! Validate with the HEP-library trace tests and the trace_form scalar checks.

use std::{hint::black_box, process::Command, time::Instant};

use idenso::{gamma, representations::Bispinor};
use spenso::{
    structure::representation::{Minkowski, RepName},
    trace,
};
use symbolica::{atom::AtomCore, symbol};

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let form = args
        .get(1)
        .expect("usage: trace_scaling /path/to/form [maximum_even_length]");
    let maximum: usize = args.get(2).map_or(14, |n| n.parse().unwrap());
    assert!((8..=14).contains(&maximum) && maximum.is_multiple_of(2));
    let version = Command::new(form).arg("-v").output().expect("launch FORM");
    assert!(version.status.success());
    eprintln!("{}", String::from_utf8_lossy(&version.stdout).trim());
    eprintln!(
        "Idenso: five warm loop samples. FORM: three processes, six batches each; wall includes all six, CPU is median of the last five. Counts are diagnostics."
    );
    let directory =
        std::env::temp_dir().join(format!("idenso-trace-scaling-{}", std::process::id()));
    std::fs::create_dir(&directory).unwrap();
    eprintln!("FORM programs and raw output: {}", directory.display());

    idenso::representations::initialize();
    let mink = Minkowski {}.new_rep(4);
    let spin = Bispinor {}.new_rep(4).to_symbolic([]);
    println!("engine,length,terms,sample,repetitions,wall_ms_per_trace,cpu_ms_per_trace");
    for (n, form_repeats) in [(8, 2000), (10, 500), (12, 100), (14, 30)] {
        if n > maximum {
            break;
        }
        let indices: Vec<_> = (0..n)
            .map(|i| mink.pattern(symbol!(&format!("trace_benchmark::mu{i}"))))
            .collect();
        let input = trace!(&spin; indices.iter().map(|mu| gamma!(mu)));
        let terms = idenso::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
            .unwrap()
            .simplify_gamma(idenso::dirac::GammaSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
            .nterms();
        let repeats = if n < 12 { 100 } else { 20 };
        for sample in 0..5 {
            let start = Instant::now();
            for _ in 0..repeats {
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
            }
            let per_trace = start.elapsed().as_secs_f64() * 1000. / f64::from(repeats);
            println!("idenso,{n},{terms},{sample},{repeats},{per_trace:.9},");
        }

        let word = (0..n)
            .map(|i| format!("mu{i}"))
            .collect::<Vec<_>>()
            .join(",");
        let source = format!(
            "Off Statistics;\nIndices {word};\n#do batch=1,6\n#do j=1,{form_repeats}\nLocal F`j'=g_(1,{word});\n#enddo\n.sort\n#reset timer\ntrace4,1;\n.sort\n#write \"BATCH_CPU_MS=`timer_'\"\n#$n=termsin_(F1);\n#write \"TERMS=%$\",$n\nDrop;\n.sort\n#enddo\n.end\n"
        );
        let path = directory.join(format!("free-{n}.frm"));
        std::fs::write(&path, source).unwrap();
        for sample in 0..3 {
            let start = Instant::now();
            let output = Command::new(form)
                .arg("-q")
                .arg(&path)
                .current_dir(&directory)
                .output()
                .unwrap();
            let per_trace = start.elapsed().as_secs_f64() * 1000. / f64::from(6 * form_repeats);
            assert!(
                output.status.success(),
                "{}\n{}",
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr)
            );
            let stdout = String::from_utf8(output.stdout).unwrap();
            let mut cpu: Vec<f64> = stdout
                .lines()
                .filter_map(|line| line.strip_prefix("BATCH_CPU_MS="))
                .map(|value| value.trim().parse().unwrap())
                .skip(1)
                .collect();
            assert_eq!(cpu.len(), 5, "expected six FORM timing batches");
            cpu.sort_by(f64::total_cmp);
            let per_trace_cpu = cpu[2] / f64::from(form_repeats);
            let form_terms: usize = stdout
                .lines()
                .find_map(|line| line.strip_prefix("TERMS="))
                .unwrap()
                .trim()
                .parse()
                .unwrap();
            println!(
                "form,{n},{form_terms},{sample},{},{per_trace:.9},{per_trace_cpu:.9}",
                6 * form_repeats
            );
            std::fs::write(path.with_extension(format!("{sample}.log")), stdout).unwrap();
        }
    }
}
