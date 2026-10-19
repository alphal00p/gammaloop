//! Small opt-in workflow. Both output directories must be new.
//! cargo run -p vakint --example native_cache -- PREPARED_DIR SNAPSHOT_DIR
use rustred_app::CachePreparationStatus;
use std::{path::PathBuf, sync::atomic::AtomicBool};
use symbolica::atom::{Atom, AtomCore};
use vakint::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedCacheSaveMode,
    RustRedEvaluationOptions, Vakint, VakintSettings, vakint_parse,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let paths = std::env::args_os()
        .skip(1)
        .map(PathBuf::from)
        .collect::<Vec<_>>();
    if paths.len() != 2 {
        return Err("usage: native_cache PREPARED_DIR SNAPSHOT_DIR (new directories)".into());
    }
    if paths.iter().any(|path| path.exists()) {
        return Err(
            "both output directories must be new; existing caches are never overwritten".into(),
        );
    }
    let engine = Vakint::new()?;
    let report = engine.prepare_cache(
        1,
        &paths[0],
        &Default::default(),
        &AtomicBool::new(false),
        |_| Ok(()),
    )?;
    println!("preparation={}", report.to_document());
    if report.status != CachePreparationStatus::Complete {
        return Err("preparation is partial; inspect its report and resume explicitly".into());
    }
    engine.load_cache(1, &paths[0], &Default::default())?;
    let settings = VakintSettings {
        form_exe_path: "/native-cache-example-does-not-use-form".into(),
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: vakint::RustRedMasterMode::Retain,
            },
        )]),
        use_dot_product_notation: true,
        ..Default::default()
    };
    let prepared = vakint_parse!("topo(I1L(4,3))")?;
    let first = engine.evaluate_integral(&settings, prepared.as_view())?;
    let warm = engine.evaluate_integral(&settings, prepared.as_view())?;
    require_same(&first, &warm)?;
    // This root lies outside the small primer. It performs native reduction,
    // then becomes part of the resident overlay of the explicit snapshot.
    let held_out = vakint_parse!("topo(I1L(4,5))")?;
    let expected = engine.evaluate_integral(&settings, held_out.as_view())?;
    println!("cache={:?}", engine.cache_statistics(1)?);
    let saved = engine.save_cache(
        1,
        &paths[1],
        RustRedCacheSaveMode::FullSnapshot,
        &Default::default(),
    )?;
    println!("snapshot={}", saved.to_document());
    if !saved.complete {
        return Err(
            "snapshot is partial; the report distinguishes copied base and resident overlay".into(),
        );
    }
    let restored = Vakint::new()?;
    restored.load_cache(1, &paths[1], &Default::default())?;
    require_same(
        &expected,
        &restored.evaluate_integral(&settings, held_out.as_view())?,
    )?;
    assert_eq!(
        restored
            .cache_statistics(1)?
            .native
            .work
            .rule_applications(),
        0
    );
    println!("exact_parity=passed snapshot_hit_without_native_rules=true");
    Ok(())
}

fn require_same(left: &Atom, right: &Atom) -> Result<(), Box<dyn std::error::Error>> {
    if !(left.clone() - right).expand().together().is_zero() {
        return Err("cached value differs from its exact reference".into());
    }
    Ok(())
}
