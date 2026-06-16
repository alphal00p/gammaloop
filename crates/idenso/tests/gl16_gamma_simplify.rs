use std::time::Instant;

use idenso::{
    Cookable, IndexTooling,
    color::ColorSimplifier,
    dirac::GammaSimplifier,
    representations::initialize,
    shorthands::schoonschip::{Schoonschip, SchoonschipSettings},
};
use itertools::Itertools;
use spenso::shadowing::{TensorCollectExt, symbolica_utils::SpensoPrintSettings};
use symbolica::prelude::*;

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/aa_aa_2l_gl16_integrated_uv_start_after_simplify_metrics.sym"
);

#[test]
#[ignore = "diagnostic timing for GL16 integrated-UV gamma simplification"]
fn aa_aa_2l_gl16_integrated_uv_simplify_gamma_fixture() {
    initialize();

    let start = Instant::now();
    let input = std::fs::read_to_string(FIXTURE)
        .unwrap_or_else(|error| panic!("failed to read {FIXTURE}: {error}"));
    eprintln!(
        "read fixture: source_bytes={} elapsed={:.3?}",
        input.len(),
        start.elapsed()
    );

    let start = Instant::now();
    let expr = parse!(input)
        .replace(parse_lit!(gammalooprs::hedge(x_)))
        .with(parse_lit!(h(x_)))
        .cook_indices()
        .collect_gamma_chains()
        .schoonschip_with_settings(&SchoonschipSettings {
            simplify_chain_like_functions: true,
            schoonschip_rank1_tensors: true,
            ..Default::default()
        })
        .collect_chains_and_traces()
        .collect_color()
        .collect_factors();
    eprintln!(
        "pre_gamma: terms={} atom_bytes={} elapsed={:.3?}, printer={}",
        expr.nterms(),
        expr.as_view().get_byte_size(),
        start.elapsed(),
        expr.printer(SpensoPrintSettings::compact().nice_symbolica()),
    );

    let start = Instant::now();
    let after_gamma = expr.simplify_gamma();

    eprintln!(
        "simplify_gamma: terms={} atom_bytes={} count={} elapsed={:.3?}",
        after_gamma.nterms(),
        after_gamma.as_view().get_byte_size(),
        after_gamma
            .alias_subexpressions(|_a, _count, _i| None)
            .count_operations(),
        start.elapsed(),
    );
    let after_aliasing = after_gamma.alias_subtensors("T");
    eprintln!(
        "after_aliasing: terms={} atom_bytes={} nop={}, \n{},\n{}",
        after_aliasing.get_root().nterms(),
        after_aliasing.get_root().as_view().get_byte_size(),
        after_aliasing.count_operations(),
        after_aliasing
            .get_root()
            .printer(SpensoPrintSettings::compact().nice_symbolica()),
        after_aliasing
            .get_aliases()
            .iter()
            .map(|(k, v)| format!(
                "{}: {}",
                k.printer(SpensoPrintSettings::compact().nice_symbolica()),
                v.printer(SpensoPrintSettings::compact().nice_symbolica())
            ))
            .join("\n")
    );

    let after_horner = after_aliasing.collect_horner::<Symbol>(None);
    eprintln!(
        "after_horner: terms={} atom_bytes={} nop={}",
        after_horner.get_root().nterms(),
        after_horner.get_root().as_view().get_byte_size(),
        after_horner.count_operations(),
    );

    let after_expand = after_aliasing.get_root().collect_metrics();
    eprintln!(
        "after_expand: terms={} atom_bytes={}\n {}",
        after_expand.nterms(),
        after_expand.as_view().get_byte_size(),
        after_expand.printer(SpensoPrintSettings::compact().nice_symbolica()),
    );

    let start = Instant::now();
    let after_collect_chains = after_gamma.collect_chains_and_traces();
    eprintln!(
        "collect_chains_and_traces: terms={} atom_bytes={} elapsed={:.3?}",
        after_collect_chains.nterms(),
        after_collect_chains.as_view().get_byte_size(),
        start.elapsed()
    );
    assert_eq!(after_collect_chains.nterms(), after_gamma.nterms());
    assert_eq!(
        after_collect_chains.as_view().get_byte_size(),
        after_gamma.as_view().get_byte_size()
    );

    let start = Instant::now();
    let after_collect_factors = after_gamma.collect_factors();
    eprintln!(
        "collect_factors: terms={} atom_bytes={} elapsed={:.3?}",
        after_collect_factors.nterms(),
        after_collect_factors.as_view().get_byte_size(),
        start.elapsed()
    );
    assert_eq!(after_collect_factors.nterms(), 1);
    assert!(
        after_collect_factors.as_view().get_byte_size() < after_gamma.as_view().get_byte_size()
    );
}
