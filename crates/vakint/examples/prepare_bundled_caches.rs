//! Explicit offline preparation; ordinary Vakint evaluation never writes these packs.
//! cargo run --release -p vakint --example prepare_bundled_caches -- LOOP OUTPUT
use rustred_app::CachePreparationStatus;
use std::{path::PathBuf, sync::atomic::AtomicBool};
use vakint::{RustRedCachePreparationOptions, Vakint};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args_os().skip(1);
    let loops: usize = args
        .next()
        .ok_or("missing loop count")?
        .to_str()
        .ok_or("invalid loop count")?
        .parse()?;
    let output = PathBuf::from(args.next().ok_or("missing output directory")?);
    if args.next().is_some() || !(1..=4).contains(&loops) {
        return Err("usage: prepare_bundled_caches LOOP(1..4) NEW_OUTPUT_DIRECTORY".into());
    }
    if output.exists() {
        return Err("output directory must be new; preserve existing packs".into());
    }
    let engine = Vakint::new()?;
    let options = RustRedCachePreparationOptions::default();
    let report =
        engine.prepare_cache(
            loops,
            &output,
            &options,
            &AtomicBool::new(false),
            |_| Ok(()),
        )?;
    println!("{}", report.to_document());
    if report.status != CachePreparationStatus::Complete {
        return Err("preparation incomplete; preserve and inspect its native report".into());
    }
    let expected = [5, 77, 2320, 164352][loops - 1];
    if report.total_roots != expected || report.completed_roots != expected {
        return Err("completed primer has an unexpected selected-root count".into());
    }
    engine.load_cache(loops, &output, &Default::default())?;
    if !engine.cache_statistics(loops)?.attached {
        return Err("fresh generated bundle did not attach".into());
    }
    Ok(())
}
