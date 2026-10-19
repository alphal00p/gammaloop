//! The shipped packs must bind to the current embedded programs and serve exact results.
use super::four_loop::cache as four_loop;
use crate::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions,
    VakintSettings, vakint_parse,
};
use symbolica::atom::AtomCore;

#[test]
fn bundled_prepared_caches_match_uncached_reductions() {
    let settings = VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: crate::RustRedMasterMode::Retain,
            },
        )]),
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        form_exe_path: "/bundled-cache-must-not-use-form".into(),
        use_dot_product_notation: true,
        ..Default::default()
    };
    for loops in 1..=4 {
        let cached = super::cache_tests::external_engine();
        let ordinary = super::cache_tests::external_engine();
        let path = super::cache_tests::external_path().join(format!("{loops}L/cache"));
        cached
            .load_cache(loops, &path, &Default::default())
            .unwrap();
        let mut inputs = Vec::new();
        for mass in [1, 4] {
            match loops {
                1 => inputs.push(vakint_parse!(&format!("topo(I1L({mass},2))")).unwrap()),
                2 => {
                    inputs.push(vakint_parse!(&format!("topo(I2L({mass},2,1,1))")).unwrap());
                    inputs
                        .push(vakint_parse!(&format!("topo(I2L_pinch_3({mass},2,1,0))")).unwrap());
                }
                3 => {
                    inputs.push(vakint_parse!(&format!("topo(I3L({mass},2,1,1,1,1,1))")).unwrap());
                    inputs.push(
                        vakint_parse!(&format!("topo(I3L_pinch_1_6({mass},0,1,1,1,1,0))")).unwrap(),
                    );
                }
                4 => {
                    for (namespace, mask) in four_loop::parent_masks(&cached).unwrap() {
                        let mut powers = mask.into_iter().map(i64::from).collect::<Vec<_>>();
                        powers[0] = 2;
                        inputs.push(four_loop::test_input(namespace, &powers, mass));
                    }
                }
                _ => unreachable!(),
            }
        }
        for input in inputs {
            let actual = cached
                .evaluate_integral(&settings, input.as_view())
                .unwrap();
            let expected = ordinary
                .evaluate_integral(&settings, input.as_view())
                .unwrap();
            assert!(
                (actual - expected).expand().together().is_zero(),
                "loops={loops}, input={input}"
            );
        }
        let stats = cached.cache_statistics(loops).unwrap();
        assert_eq!(
            stats.native.completed_reductions, 0,
            "prepared input missed at loops={loops}"
        );
        assert_eq!(
            stats.native.work.rule_applications(),
            0,
            "cached input ran IBPs at loops={loops}"
        );
        assert!(
            stats.native.backing_hits > 0,
            "bundle was not consumed at loops={loops}"
        );
    }
}

#[test]
fn bundled_five_loop_application_loads_from_its_portable_directory() {
    use sha2::{Digest, Sha256};
    let directory = super::cache_tests::external_path().join("5L/ibps");
    let bytes = std::fs::read(directory.join("application.json")).unwrap();
    assert_eq!(
        format!("{:x}", Sha256::digest(&bytes)),
        "01425e1cecc71893d2e3a592a162600100326f285914926df61e43dd0c870e4b"
    );
    let document: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(document["metadata"]["remaining_terminals"], 905);
    assert_eq!(document["metadata"]["source_generation"], 54);
    let mut settings = VakintSettings::default();
    assert!(settings.rustred_artifact_path.is_none());
    let engine = super::cache_tests::external_engine()
        .with_rustred_application(&mut settings, &directory)
        .unwrap();
    assert_eq!(
        settings.rustred_artifact_path,
        Some(directory.canonicalize().unwrap())
    );
    // Repeated explicit attachment reuses the existing authenticated application.
    let _ = engine
        .with_rustred_application(&mut settings, &directory)
        .unwrap();
    assert!(VakintSettings::default().rustred_artifact_path.is_none());
}
