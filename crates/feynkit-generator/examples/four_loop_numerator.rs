//! Run with `cargo run --profile dev-optim -p feynkit-generator --example
//! four_loop_numerator -- [threads=1] [diagrams.jsonl|-] [filter_zero_color=false]`.
use std::{io::Write, sync::Mutex, time::Instant};

use feynkit_generator::{
    GenerationControl, GenerationFilter, GenerationOptions, Process, SnailFilterOptions,
};
use feynkit_model::Model;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut arguments = std::env::args().skip(1);
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
    let process = Process::new(["g"], ["g"]).with_filters(
        ["u", "c", "s", "t", "b"]
            .into_iter()
            .map(Into::into)
            .collect(),
        None,
        vec![],
    );
    let stage = Mutex::new(("start", Instant::now()));
    let options = GenerationOptions::default()
        .threads(threads)
        .filter_zero_color(filter_zero_color)
        .with_loop_count(4, 4)?
        .max_vertices(8)
        .with_graph_filter(GenerationFilter::MaxNumberOfBridges(0))
        .with_graph_filter(GenerationFilter::ZeroSnails(SnailFilterOptions {
            veto_attached_to_massive: true,
            ..Default::default()
        }))
        .progress(move |progress| {
            let mut stage = stage.lock().unwrap();
            if stage.0 != progress.stage {
                eprintln!("{}: {:.3}s", stage.0, stage.1.elapsed().as_secs_f64());
                *stage = (progress.stage, Instant::now());
            }
            GenerationControl::Continue
        });
    let start = Instant::now();
    let result = process.generate_diagrams(model, &options)?;
    eprintln!(
        "{} diagrams, {} zero numerators in {:.3}s",
        result.diagrams.len(),
        result.report.zero_numerator_count,
        start.elapsed().as_secs_f64()
    );
    assert!(result.report.completed);
    assert_eq!(
        result.diagrams.len() + result.report.zero_numerator_count,
        4970
    );
    if let Some(output) = output {
        let mut output = std::io::BufWriter::new(std::fs::File::create(output)?);
        for diagram in result.diagrams {
            let diagram: serde_json::Value = serde_json::from_str(&diagram.to_json()?)?;
            serde_json::to_writer(&mut output, &diagram)?;
            writeln!(output)?;
        }
    }
    Ok(())
}
