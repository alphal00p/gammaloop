//! Offline oracle revalidation of already-declared terminals; never solve IBPs.

use super::{OfflineTerminalCatalog, acceptance_inputs::ParentDescriptor, input::ParentInput};
use std::{collections::BTreeMap, path::PathBuf};
use symbolica::atom::{Atom, AtomCore};
use vakint::{
    EvaluationMethod, EvaluationOrder, FMFTOptions, LoopNormalizationFactor, TensorReductionMethod,
    Vakint, VakintSettings, vakint_parse,
};

#[test]
#[ignore = "offline FMFT catalog audit; requires explicit oracle executable"]
fn shipped_terminal_projections_match_fmft() {
    super::test_utils::run_multi_lane_acceptance(|| {
        let settings = VakintSettings {
            form_exe_path: std::env::var("VAKINT_4L_CANDIDATE_ORACLE_FORM_PATH").unwrap(),
            epsilon_symbol: "ep".into(),
            tensor_reduction_method: TensorReductionMethod::FeynKit,
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {
                expand_masters: false,
                susbstitute_masters: false,
            })]),
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
            use_dot_product_notation: true,
            ..VakintSettings::default()
        };
        let vakint = Vakint::new().unwrap();
        let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/rustred/four_loop");
        let output = std::env::var_os("VAKINT_4L_CATALOG_AUDIT_OUTPUT").map(PathBuf::from);
        let filter = std::env::var("VAKINT_4L_CANDIDATE_FAMILY_FILTER").ok();
        let mut checked = 0;
        for (name, descriptor) in [
            ("H", ParentDescriptor::H),
            ("FG", ParentDescriptor::Fg),
            ("BMW", ParentDescriptor::Bmw),
            ("X", ParentDescriptor::X),
        ] {
            if filter.as_deref().is_some_and(|f| f != name) {
                continue;
            }
            let parent = ParentInput::from_csv(descriptor.csv());
            let filename = format!("{}.rrcat.bin", name.to_ascii_lowercase());
            let old = OfflineTerminalCatalog::decode_generated(
                &std::fs::read(directory.join(&filename)).unwrap(),
                parent.family.fingerprint(),
                10,
            )
            .unwrap();
            old.require_complete().unwrap();
            let mut updated = BTreeMap::new();
            let mut changes = 0;
            for (key, previous) in old.terms() {
                let input = parent.integral(key.powers());
                let value = vakint
                    .evaluate(&settings, input.as_view())
                    .unwrap()
                    .replace(vakint_parse!("mursq").unwrap().to_pattern())
                    .with(Atom::num(1).to_pattern())
                    .replace(vakint_parse!("muvsq").unwrap().to_pattern())
                    .with(Atom::num(1).to_pattern())
                    .together();
                let difference = (&value - previous).together().expand();
                if !difference.is_zero() {
                    println!(
                        "catalog CHANGE {name} {:?}: {previous} -> {value}",
                        key.powers()
                    );
                    changes += 1;
                }
                updated.insert(key.clone(), value);
            }
            println!(
                "catalog {name}: {} checked; {changes} differences",
                old.terms().len()
            );
            if let Some(output) = &output {
                let refreshed = OfflineTerminalCatalog::from_keys(
                    parent.family.fingerprint(),
                    10,
                    old.terms().keys().cloned(),
                    |key| updated.get(key).cloned().ok_or("missing oracle result"),
                )
                .unwrap();
                std::fs::create_dir_all(output).unwrap();
                std::fs::write(output.join(filename), refreshed.encode().unwrap()).unwrap();
            } else {
                assert_eq!(changes, 0, "shipped {name} catalog differs from FMFT");
            }
            checked += 1;
        }
        assert_eq!(checked, if filter.is_some() { 1 } else { 4 });
    });
}
