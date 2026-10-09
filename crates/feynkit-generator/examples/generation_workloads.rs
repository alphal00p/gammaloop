//! Diagram-generation workloads for instruction-count gates; `bin/generation_ir.py` runs them
//! under callgrind and attributes the instructions to the generation stages.
//!
//! ```text
//! cargo run -p feynkit-generator --release --example generation_workloads -- MODE CASE
//! ```
//!
//! - `full CASE` runs the whole pipeline: enumeration, interaction filters, numerators,
//!   selection and grouping.
//! - `generation CASE` cancels at the first progress report after enumeration, so only the
//!   enumeration and the certificate sort run.
//! - `list [reference]` prints the case names.
//!
//! Both runs use one generation thread and print one `RESULT` line with the report counts. The
//! cases are those of the generation oracle tests (`tests/generation_oracle/cases.rs`).

#[path = "../tests/generation_oracle/cases.rs"]
mod cases;

use cases::{CASES, Case, REFERENCE};
use feynkit_generator::GenerationControl;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let args: Vec<&str> = args.iter().map(String::as_str).collect();
    match args.as_slice() {
        ["list"] => CASES.iter().for_each(|case| println!("{}", case.name)),
        ["list", "reference"] => REFERENCE.iter().for_each(|name| println!("{name}")),
        [mode @ ("full" | "generation"), name] => run(Case::named(name), *mode == "generation"),
        _ => {
            eprintln!("usage: generation_workloads (full|generation) CASE | list [reference]");
            std::process::exit(2);
        }
    }
}

fn run(case: &Case, generation_only: bool) {
    let model = case.model();
    let mut options = case.options(&model, 1);
    if generation_only {
        options = options.progress(|progress| {
            if progress.stage == "topologies" {
                GenerationControl::Continue
            } else {
                GenerationControl::Cancel
            }
        });
    }
    let report = case
        .generate_with(model, &options)
        .unwrap_or_else(|error| panic!("{}: {error}", case.name))
        .report;
    println!(
        "RESULT case={} mode={} topologies={} assignments={} retained={} completed={}",
        case.name,
        if generation_only {
            "generation"
        } else {
            "full"
        },
        report.topology_count,
        report.interaction_assignment_count,
        report.retained_count,
        report.completed,
    );
}
