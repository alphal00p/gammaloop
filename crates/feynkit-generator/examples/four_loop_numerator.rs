//! Run with `cargo run --profile dev-optim -p feynkit-generator --example
//! four_loop_numerator -- [threads=1] [diagrams.jsonl|-] [filter_zero_color=false]`.
use num_bigint::BigUint;
use serde_json::Value;
use std::env;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Mutex;
use std::time::Instant;

use feynkit_generator::{
    Bridges, GenerationControl, Loops, Naming, NumeratorOptions, Process, Restrict, Run,
    TopologyRules, Veto,
};
use feynkit_model::Model;

fn main() -> Result<(), Box<dyn Error>> {
    let mut arguments = env::args().skip(1);
    let threads = arguments
        .next()
        .map(|value| value.parse())
        .transpose()?
        .unwrap_or(1);
    let output = arguments.next().filter(|path| path != "-");
    let filter_zero_color = arguments
        .next()
        .map(|value| value.parse())
        .transpose()?
        .unwrap_or(false);
    let model = Model::qcd();
    let request = Process::new(["g"], ["g"])
        .amplitude(Loops::total(4..=4))
        .with_restrict(Restrict::new().exclude_particles(["u", "c", "s", "t", "b"]))
        .with_topology(
            // The original options: no self-loops and every snail vetoed.
            TopologyRules::new()
                .with_self_loops(false)
                .with_max_vertices(8)
                .with_bridges(Bridges::AtMost(0))
                .with_snails(Veto::all()),
        );
    let stage = Mutex::new((None, Instant::now()));
    let run = Run::default().threads(threads).progress(move |progress| {
        let mut stage = stage.lock().unwrap();
        if stage.0 != Some(progress.stage) {
            if let Some(previous) = stage.0 {
                eprintln!("{previous:?}: {:.3}s", stage.1.elapsed().as_secs_f64());
            }
            *stage = (Some(progress.stage), Instant::now());
        }
        GenerationControl::Continue
    });
    let start = Instant::now();
    let result = request.generate(model.into(), &run)?.build(
        &NumeratorOptions::default()
            .with_naming(Naming::prefix("D"))
            .with_zero_color_filter(filter_zero_color),
        &run,
    )?;
    let report = result.report();
    eprintln!(
        "{} diagrams, {} zero numerators in {:.3}s",
        result.len(),
        report.zero_numerator_count,
        start.elapsed().as_secs_f64()
    );
    assert!(report.completed);
    assert_eq!(
        report.zero_numerator_count.clone() + result.len(),
        4970_u32.into()
    );
    if filter_zero_color {
        assert_eq!(result.len(), 4664);
        assert_eq!(report.zero_numerator_count, BigUint::from(306_u32));
    }
    if let Some(output) = output {
        let mut output = BufWriter::new(File::create(output)?);
        for diagram in result.diagrams() {
            let diagram: Value = serde_json::from_str(&diagram.to_json()?)?;
            serde_json::to_writer(&mut output, &diagram)?;
            writeln!(output)?;
        }
    }
    Ok(())
}
