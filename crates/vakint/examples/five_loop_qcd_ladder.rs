//! Slow native preparation only; saved-union census lives in the full slow test.
#[path = "../tests/fixtures/qcd_ladder/mod.rs"]
mod fixture;
use ::vakint;
use symbolica::atom::AtomCore;

fn main() -> fixture::Result<()> {
    let arguments = std::env::args_os().skip(1).collect::<Vec<_>>();
    let evidence = arguments
        .first()
        .ok_or("usage: five_loop_qcd_ladder EVIDENCE_DIRECTORY [--dress-only]")?;
    if arguments
        .get(1)
        .is_some_and(|value| value == "--dress-only")
    {
        let started = std::time::Instant::now();
        let generated = fixture::dress(std::path::Path::new(evidence))?;
        println!(
            "native_qcd_ladder_dressed seconds={:.6} diagrams={} loops={} vertices={} edges={}",
            started.elapsed().as_secs_f64(),
            generated.diagrams.len(),
            generated.diagrams[0].loop_count(),
            generated.diagrams[0].vertices().count(),
            generated.diagrams[0].edges().count()
        );
        return Ok(());
    }
    let prepared = fixture::prepare(std::path::Path::new(&evidence))?;
    println!(
        "native_qcd_ladder_prepared seconds={:.6} scalar_terms={} input_bytes={} contributions={} guards={}",
        prepared.preparation_seconds,
        prepared.scalar_terms,
        prepared.input.to_canonical_string().len(),
        prepared.contributions,
        prepared.nonzero_conditions.len()
    );
    Ok(())
}
