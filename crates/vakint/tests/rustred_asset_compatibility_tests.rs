//! Native codec/replay regression only: no rule generation or numerical oracle.
use std::{collections::BTreeSet, fs, io::Read, path::PathBuf};

use rustred::foundry::artifact::ClosedArtifact;
use rustred::persistence::{BinaryIoLimits, ExactTerminalCatalog, TerminalCatalogCoverage};
use rustred::reduction::{ReductionLimits, terminal_normalization::TerminalNormalizationPlan};
use rustred_app::{CandidateBundleLimits, load_generated_candidate_bundle};

fn asset(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("data/rustred")
        .join(name)
}

fn compressed(name: &str, limit: usize) -> Vec<u8> {
    let bytes = fs::read(asset(name)).unwrap();
    assert!(bytes.len() <= limit);
    let mut decoder = flate2::bufread::GzDecoder::new(bytes.as_slice());
    let mut output = Vec::new();
    decoder
        .by_ref()
        .take(limit as u64 + 1)
        .read_to_end(&mut output)
        .unwrap();
    assert!(output.len() <= limit);
    assert!(decoder.into_inner().is_empty(), "single gzip member only");
    output
}

#[test]
fn sealed_assets_cold_replay_and_roundtrip_with_current_native_codec() {
    let _ = symbolica::parse!(
        "asset_codec_noise::f(a,b)",
        default_namespace = "asset_codec_noise"
    );
    for arity in [1, 3, 6] {
        let bytes = fs::read(asset(&format!("unit_mass_vacuum_k{arity}.rrbin"))).unwrap();
        let loaded = ClosedArtifact::decode_durable(&bytes).unwrap();
        assert_eq!(loaded.arity(), arity);
        let encoded = loaded.encode_durable().unwrap();
        let repeated = ClosedArtifact::decode_durable(&encoded).unwrap();
        assert_eq!(repeated.family_fingerprint(), loaded.family_fingerprint());
        assert!(
            rustred::persistence::equivalent_generated_programs(
                &bytes,
                &encoded,
                Default::default()
            )
            .unwrap()
        );
    }
}

#[test]
fn four_loop_assets_replay_normalizers_and_exact_catalog_output_sets() {
    let _ = symbolica::parse!(
        "asset_codec_noise::g(c,d)",
        default_namespace = "asset_codec_noise"
    );
    for (name, raw_count, output_count) in [
        ("h", 386, 22),
        ("x", 445, 19),
        ("bmw", 179, 17),
        ("fg", 145, 16),
    ] {
        let bytes = compressed(
            &format!("four_loop/{name}.candidates.rrbin.gz"),
            512 * 1024 * 1024,
        );
        let (family, mut reducer) = load_generated_candidate_bundle::<10>(
            &bytes,
            CandidateBundleLimits {
                max_collection_entries: 8_000_000,
                max_total_coefficient_bytes: 512 * 1024 * 1024,
                ..Default::default()
            },
            ReductionLimits::default(),
        )
        .unwrap();
        let plan = TerminalNormalizationPlan::decode_generated(
            &fs::read(asset(&format!("four_loop/{name}.rrnorm.bin"))).unwrap(),
            &family,
            reducer.terminals(),
            reducer.ordering(),
            Default::default(),
            BinaryIoLimits::default(),
        )
        .unwrap();
        assert_eq!(plan.raw_terminals(), reducer.terminals());
        assert_eq!(plan.raw_terminals().len(), raw_count);
        assert_eq!(plan.canonical_terminals().len(), output_count);
        let catalog = ExactTerminalCatalog::decode_generated(
            &compressed(&format!("four_loop/{name}.rrcat.bin.gz"), 64 * 1024),
            family.fingerprint(),
            10,
            Default::default(),
        )
        .unwrap();
        assert_eq!(catalog.coverage(), TerminalCatalogCoverage::Complete);
        assert_eq!(
            catalog.terms().keys().cloned().collect::<BTreeSet<_>>(),
            *plan.canonical_terminals()
        );
        reducer.install_terminal_normalization(plan).unwrap();
    }
}
