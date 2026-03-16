use color_eyre::Result;
use gammaloop_integration_tests::{clean_test, get_tests_workspace_path};
use gammalooprs::{
    integrands::evaluation::GraphEvaluationResult,
    momentum::{Energy, FourMomentum, ThreeMomentum},
    observables::{
        AdditionalWeightKey, CutInfo, EntrySelection, Event, EventGroupList,
        EventProcessingRuntime, FilterQuantity, GenericEventGroupList, HistogramSettings,
        ObservablePhase, ObservableSettings, ObservableSnapshotBundle, ObservableValueTransform,
        ParticleScalarQuantitySettings, QuantitiesSettings, QuantitySettings,
        SelectorDefinitionSettings, SelectorReduction, SelectorSettings,
        ValueRangeSelectorSettings,
    },
    settings::RuntimeSettings,
    utils::{F, f128},
};
use spenso::algebra::complex::Complex;
use std::{collections::BTreeMap, fs, path::PathBuf, time::Duration};

fn new_artifact_dir(name: &str) -> Result<PathBuf> {
    let path = get_tests_workspace_path().join(name);
    if path.exists() {
        clean_test(&path);
    }
    fs::create_dir_all(&path)?;
    Ok(path)
}

fn momentum(px: f64, py: f64, pz: f64, energy: f64) -> FourMomentum<F<f64>> {
    FourMomentum {
        temporal: Energy { value: F(energy) },
        spatial: ThreeMomentum {
            px: F(px),
            py: F(py),
            pz: F(pz),
        },
    }
}

fn make_event(
    pt_x: f64,
    pt_y: f64,
    pz: f64,
    energy: f64,
    pdg: isize,
    weight: (f64, f64),
    additional_weights: impl IntoIterator<Item = (AdditionalWeightKey, (f64, f64))>,
) -> Event {
    let mut weights = BTreeMap::new();
    for (key, (re, im)) in additional_weights {
        weights.insert(key, Complex::new(F(re), F(im)));
    }

    Event {
        kinematic_configuration: (
            vec![momentum(0.0, 0.0, 1.0, 1.0), momentum(0.0, 0.0, -1.0, 1.0)].into(),
            vec![momentum(pt_x, pt_y, pz, energy)].into(),
        ),
        cut_info: CutInfo {
            particle_pdgs: (Default::default(), vec![pdg].into()),
            cut_id: 0,
            graph_id: 0,
        },
        weight: Complex::new(F(weight.0), F(weight.1)),
        additional_weights: gammalooprs::observables::GenericAdditionalWeightInfo { weights },
    }
}

fn singleton_groups(event: Event) -> EventGroupList {
    let mut groups = EventGroupList::default();
    groups.push_singleton(event);
    groups
}

fn runtime_settings() -> RuntimeSettings {
    let mut settings = RuntimeSettings::default();
    settings.quantities = QuantitiesSettings::default();
    settings.observables.clear();
    settings.selectors.clear();

    settings.quantities.insert(
        "pt".to_string(),
        QuantitySettings::ParticleScalar(ParticleScalarQuantitySettings {
            pdgs: vec![1],
            quantity: FilterQuantity::PT,
        }),
    );

    settings.observables.insert(
        "pt_real".to_string(),
        ObservableSettings {
            quantity: "pt".to_string(),
            entry_selection: EntrySelection::All,
            entry_index: 0,
            value_transform: ObservableValueTransform::Identity,
            phase: ObservablePhase::Real,
            misbinning_max_normalized_distance: None,
            histogram: HistogramSettings {
                x_min: 0.0,
                x_max: 1.0,
                n_bins: 2,
                log_x_axis: false,
                log_y_axis: true,
            },
        },
    );

    settings.observables.insert(
        "pt_imag".to_string(),
        ObservableSettings {
            quantity: "pt".to_string(),
            entry_selection: EntrySelection::All,
            entry_index: 0,
            value_transform: ObservableValueTransform::Identity,
            phase: ObservablePhase::Imag,
            misbinning_max_normalized_distance: None,
            histogram: HistogramSettings {
                x_min: 0.0,
                x_max: 1.0,
                n_bins: 2,
                log_x_axis: false,
                log_y_axis: true,
            },
        },
    );

    settings.selectors.insert(
        "pt_window".to_string(),
        SelectorSettings {
            quantity: "pt".to_string(),
            entry_selection: EntrySelection::All,
            entry_index: 0,
            selector: SelectorDefinitionSettings::ValueRange(ValueRangeSelectorSettings {
                min: 0.2,
                max: Some(0.3),
                reduction: SelectorReduction::AnyInRange,
            }),
        },
    );

    settings
}

fn assert_close(lhs: f64, rhs: f64) {
    assert!(
        (lhs - rhs).abs() < 1.0e-12,
        "expected {rhs:.16e}, got {lhs:.16e}"
    );
}

#[test]
fn graph_evaluation_result_merges_groups_and_downcasts() {
    let first_event = make_event(
        0.25,
        0.0,
        0.0,
        0.25,
        1,
        (1.5, 0.5),
        [
            (AdditionalWeightKey::Original, (1.5, 0.5)),
            (
                AdditionalWeightKey::ThresholdCounterterm { subset_index: 3 },
                (0.25, -0.75),
            ),
        ],
    );
    let second_event = make_event(
        0.75,
        0.0,
        0.0,
        0.75,
        1,
        (2.0, -1.0),
        [(AdditionalWeightKey::Original, (2.0, -1.0))],
    );

    let mut lhs = GraphEvaluationResult {
        integrand_result: Complex::new(F::<f128>::from_f64(1.0), F::<f128>::from_f64(2.0)),
        event_groups: GenericEventGroupList::<f128>::from_f64(&singleton_groups(first_event)),
        event_processing_time: Duration::from_millis(5),
        generated_event_count: 1,
        accepted_event_count: 1,
    };
    let rhs = GraphEvaluationResult {
        integrand_result: Complex::new(F::<f128>::from_f64(3.0), F::<f128>::from_f64(-4.0)),
        event_groups: GenericEventGroupList::<f128>::from_f64(&singleton_groups(second_event)),
        event_processing_time: Duration::from_millis(7),
        generated_event_count: 2,
        accepted_event_count: 1,
    };

    lhs.merge_in_place(rhs);
    let merged = lhs.into_f64();

    assert_close(merged.integrand_result.re.0, 4.0);
    assert_close(merged.integrand_result.im.0, -2.0);
    assert_eq!(merged.event_groups.len(), 2);
    assert_eq!(merged.generated_event_count, 3);
    assert_eq!(merged.accepted_event_count, 2);
    assert_eq!(merged.event_processing_time, Duration::from_millis(12));
    assert!(
        merged.event_groups[0][0]
            .additional_weights
            .weights
            .contains_key(&AdditionalWeightKey::Original)
    );
    assert!(
        merged.event_groups[0][0]
            .additional_weights
            .weights
            .contains_key(&AdditionalWeightKey::ThresholdCounterterm { subset_index: 3 })
    );
}

#[test]
fn event_processing_runtime_filters_events_and_builds_snapshots() -> Result<()> {
    let settings = runtime_settings();
    let mut selector_runtime = EventProcessingRuntime::from_settings(&settings)?;
    let mut observable_runtime = EventProcessingRuntime::from_settings(&settings)?;

    let accepted_event = make_event(
        0.25,
        0.0,
        0.0,
        0.25,
        1,
        (2.0, 4.0),
        [(AdditionalWeightKey::Original, (2.0, 4.0))],
    );
    let rejected_event = make_event(
        0.75,
        0.0,
        0.0,
        0.75,
        1,
        (2.0, 4.0),
        [(AdditionalWeightKey::Original, (2.0, 4.0))],
    );

    assert!(selector_runtime.process_event(&accepted_event));
    assert!(!selector_runtime.process_event(&rejected_event));

    observable_runtime.process_event_groups(&singleton_groups(accepted_event.clone()));
    observable_runtime.process_event_groups(&singleton_groups(accepted_event));
    observable_runtime.update_results(0);

    let snapshot = observable_runtime.snapshot_bundle();
    let real_histogram = snapshot
        .histograms
        .get("pt_real")
        .expect("missing real-valued histogram");
    let imag_histogram = snapshot
        .histograms
        .get("pt_imag")
        .expect("missing imag-valued histogram");

    assert_eq!(real_histogram.bins.len(), 2);
    assert_close(
        real_histogram.bins[0].average(real_histogram.sample_count),
        2.0,
    );
    assert_close(
        real_histogram.bins[1].average(real_histogram.sample_count),
        0.0,
    );
    assert_close(
        imag_histogram.bins[0].average(imag_histogram.sample_count),
        4.0,
    );
    assert_close(
        imag_histogram.bins[1].average(imag_histogram.sample_count),
        0.0,
    );
    assert_eq!(real_histogram.underflow_bin.entry_count, 0);
    assert_eq!(real_histogram.overflow_bin.entry_count, 0);

    let artifact_dir = new_artifact_dir("test_differential_runtime_snapshots")?;
    let json_path = artifact_dir.join("observables.json");
    let hwu_path = artifact_dir.join("observables.hwu");

    snapshot.to_json_file(&json_path)?;
    let round_tripped = ObservableSnapshotBundle::from_json_file(&json_path)?;
    assert_eq!(round_tripped, snapshot);

    snapshot.write_hwu_file(&hwu_path)?;
    let hwu_contents = fs::read_to_string(&hwu_path)?;
    assert!(hwu_contents.contains("pt_real"));
    assert!(hwu_contents.contains("pt_imag"));

    clean_test(&artifact_dir);
    Ok(())
}

#[test]
fn event_processing_runtime_merges_worker_results_and_batch_histograms() -> Result<()> {
    let settings = runtime_settings();

    let left_event = make_event(
        0.25,
        0.0,
        0.0,
        0.25,
        1,
        (2.0, 0.0),
        [(AdditionalWeightKey::Original, (2.0, 0.0))],
    );
    let right_event = make_event(
        0.75,
        0.0,
        0.0,
        0.75,
        1,
        (4.0, 0.0),
        [(AdditionalWeightKey::Original, (4.0, 0.0))],
    );

    let mut master = EventProcessingRuntime::from_settings(&settings)?;
    master.process_event_groups(&singleton_groups(left_event.clone()));

    let mut worker = master.cleared_observable_clone();
    worker.process_event_groups(&singleton_groups(right_event.clone()));
    master.merge_samples(&mut worker)?;
    master.update_results(0);

    let merged_snapshot = master.snapshot_bundle();
    let merged_real_histogram = merged_snapshot
        .histograms
        .get("pt_real")
        .expect("missing merged real histogram");
    assert_close(
        merged_real_histogram.bins[0].average(merged_real_histogram.sample_count),
        1.0,
    );
    assert_close(
        merged_real_histogram.bins[1].average(merged_real_histogram.sample_count),
        2.0,
    );

    let mut batch_master = EventProcessingRuntime::from_settings(&settings)?;
    batch_master.process_event_groups(&singleton_groups(left_event));
    let mut batch_worker = batch_master.cleared_observable_clone();
    batch_worker.process_event_groups(&singleton_groups(right_event));
    let mut bundle = batch_worker.accumulator_bundle();
    batch_master.merge_accumulator_bundle(&mut bundle)?;
    batch_master.update_results(0);

    let batch_snapshot = batch_master.snapshot_bundle();
    let batch_real_histogram = batch_snapshot
        .histograms
        .get("pt_real")
        .expect("missing batch-merged real histogram");
    assert_close(
        batch_real_histogram.bins[0].average(batch_real_histogram.sample_count),
        1.0,
    );
    assert_close(
        batch_real_histogram.bins[1].average(batch_real_histogram.sample_count),
        2.0,
    );

    let rebinned = batch_real_histogram.rebin(2)?;
    assert_eq!(rebinned.bins.len(), 1);
    assert_eq!(rebinned.bins[0].entry_count, 2);
    assert_close(rebinned.bins[0].average(rebinned.sample_count), 3.0);

    let merged_again = merged_real_histogram.merge(batch_real_histogram)?;
    assert_eq!(merged_again.sample_count, 4);
    assert_eq!(merged_again.bins[0].entry_count, 2);
    assert_eq!(merged_again.bins[1].entry_count, 2);

    Ok(())
}
