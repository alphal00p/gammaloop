use color_eyre::Result;
use gammaloop_integration_tests::{clean_test, get_tests_workspace_path};
use gammalooprs::{
    integrands::evaluation::GraphEvaluationResult,
    momentum::{Energy, FourMomentum, ThreeMomentum},
    observables::{
        AdditionalWeightKey, CutInfo, EntrySelection, Event, EventGroupList,
        EventProcessingRuntime, FilterQuantity, GenericEventGroupList, HistogramSettings,
        JetClusteringSettings, JetQuantitySettings, ObservablePhase, ObservableSettings,
        ObservableSnapshotBundle, ObservableValueTransform, PairQuantity, ParticleQuantitySettings,
        QuantitiesSettings, QuantityComputation, QuantityComputationSettings, QuantityOrder,
        QuantityOrdering, QuantitySettings, SelectorDefinitionSettings, SelectorReduction,
        SelectorSettings, ValueRangeSelectorSettings,
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
    make_multi_event(&[(pt_x, pt_y, pz, energy, pdg)], weight, additional_weights)
}

fn make_multi_event(
    outgoing: &[(f64, f64, f64, f64, isize)],
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
            outgoing
                .iter()
                .map(|(px, py, pz, energy, _)| momentum(*px, *py, *pz, *energy))
                .collect::<Vec<_>>()
                .into(),
        ),
        cut_info: CutInfo {
            particle_pdgs: (
                Default::default(),
                outgoing
                    .iter()
                    .map(|(_, _, _, _, pdg)| *pdg)
                    .collect::<Vec<_>>()
                    .into(),
            ),
            cut_id: 0,
            graph_id: 0,
        },
        weight: Complex::new(F(weight.0), F(weight.1)),
        additional_weights: gammalooprs::observables::GenericAdditionalWeightInfo { weights },
        derived_observable_data: Default::default(),
    }
}

fn singleton_groups(event: Event) -> EventGroupList {
    let mut groups = EventGroupList::default();
    groups.push_singleton(event);
    groups
}

fn single_group(events: impl IntoIterator<Item = Event>) -> EventGroupList {
    let mut groups = EventGroupList::default();
    groups.push(gammalooprs::observables::GenericEventGroup(
        events.into_iter().collect(),
    ));
    groups
}

fn process_all_events(runtime: &mut EventProcessingRuntime, groups: &mut EventGroupList) {
    for event_group in groups.iter_mut() {
        for event in event_group.iter_mut() {
            assert!(runtime.process_event(event));
        }
    }
}

fn runtime_settings() -> RuntimeSettings {
    let mut settings = RuntimeSettings {
        quantities: QuantitiesSettings::default(),
        observables: BTreeMap::new(),
        selectors: BTreeMap::new(),
        ..Default::default()
    };

    settings.quantities.insert(
        "pt".to_string(),
        QuantitySettings::Particle(ParticleQuantitySettings {
            pdgs: vec![1],
            computation: scalar_computation(FilterQuantity::PT),
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
                title: None,
                type_description: "AL".to_string(),
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
                title: None,
                type_description: "AL".to_string(),
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
                min: Some(0.2),
                max: Some(0.3),
                reduction: SelectorReduction::AnyInRange,
            }),
        },
    );

    settings
}

fn jet_count_runtime_settings(with_misbinning: bool) -> RuntimeSettings {
    let mut settings = RuntimeSettings::default();
    settings.quantities.insert(
        "jet_count".to_string(),
        QuantitySettings::Jet(JetQuantitySettings {
            clustering: JetClusteringSettings {
                clustered_pdgs: Some(vec![21]),
                ..JetClusteringSettings::default()
            },
            computation: count_computation(),
        }),
    );
    settings.observables.insert(
        "jet_count_hist".to_string(),
        ObservableSettings {
            quantity: "jet_count".to_string(),
            entry_selection: EntrySelection::All,
            entry_index: 0,
            value_transform: ObservableValueTransform::Identity,
            phase: ObservablePhase::Real,
            misbinning_max_normalized_distance: with_misbinning.then_some(0.1),
            histogram: HistogramSettings {
                x_min: 0.0,
                x_max: 6.0,
                n_bins: 6,
                log_x_axis: false,
                log_y_axis: true,
                title: None,
                type_description: "AL".to_string(),
            },
        },
    );
    settings
}

fn insert_histogram_observable(
    settings: &mut RuntimeSettings,
    name: &str,
    quantity: &str,
    x_min: f64,
    x_max: f64,
    n_bins: usize,
) {
    settings.observables.insert(
        name.to_string(),
        ObservableSettings {
            quantity: quantity.to_string(),
            entry_selection: EntrySelection::All,
            entry_index: 0,
            value_transform: ObservableValueTransform::Identity,
            phase: ObservablePhase::Real,
            misbinning_max_normalized_distance: None,
            histogram: HistogramSettings {
                x_min,
                x_max,
                n_bins,
                log_x_axis: false,
                log_y_axis: true,
                title: None,
                type_description: "AL".to_string(),
            },
        },
    );
}

fn scalar_computation(quantity: FilterQuantity) -> QuantityComputationSettings {
    QuantityComputationSettings::scalar(quantity)
}

fn scalar_computation_with_ordering(
    quantity: FilterQuantity,
    ordering: QuantityOrdering,
    order: QuantityOrder,
) -> QuantityComputationSettings {
    QuantityComputationSettings {
        quantity: Some(quantity),
        ordering: Some(ordering),
        order,
        ..QuantityComputationSettings::scalar(quantity)
    }
}

fn count_computation() -> QuantityComputationSettings {
    QuantityComputationSettings {
        computation: QuantityComputation::Count,
        ..Default::default()
    }
}

fn pair_computation(quantity: PairQuantity) -> QuantityComputationSettings {
    QuantityComputationSettings {
        computation: QuantityComputation::Pair,
        pair_quantity: Some(quantity),
        ..Default::default()
    }
}

fn pair_computation_with_order(order: QuantityOrder) -> QuantityComputationSettings {
    QuantityComputationSettings {
        order,
        ..pair_computation(PairQuantity::DeltaR)
    }
}

fn scalar_projection_runtime_settings() -> RuntimeSettings {
    let mut settings = RuntimeSettings::default();

    for (name, quantity) in [
        ("particle_px", FilterQuantity::Px),
        ("particle_py", FilterQuantity::Py),
        ("particle_pz", FilterQuantity::Pz),
        ("particle_mass", FilterQuantity::Mass),
    ] {
        settings.quantities.insert(
            name.to_string(),
            QuantitySettings::Particle(ParticleQuantitySettings {
                pdgs: vec![1],
                computation: scalar_computation(quantity),
            }),
        );
    }

    for (name, quantity) in [
        ("jet_px", FilterQuantity::Px),
        ("jet_mass", FilterQuantity::Mass),
    ] {
        settings.quantities.insert(
            name.to_string(),
            QuantitySettings::Jet(JetQuantitySettings {
                clustering: JetClusteringSettings {
                    clustered_pdgs: Some(vec![21]),
                    ..JetClusteringSettings::default()
                },
                computation: scalar_computation(quantity),
            }),
        );
    }

    insert_histogram_observable(
        &mut settings,
        "particle_px_hist",
        "particle_px",
        -1.0,
        1.0,
        4,
    );
    insert_histogram_observable(
        &mut settings,
        "particle_py_hist",
        "particle_py",
        -1.0,
        1.0,
        4,
    );
    insert_histogram_observable(
        &mut settings,
        "particle_pz_hist",
        "particle_pz",
        0.0,
        1.0,
        4,
    );
    insert_histogram_observable(
        &mut settings,
        "particle_mass_hist",
        "particle_mass",
        0.0,
        1.0,
        4,
    );
    insert_histogram_observable(&mut settings, "jet_px_hist", "jet_px", -1.0, 1.0, 4);
    insert_histogram_observable(&mut settings, "jet_mass_hist", "jet_mass", 0.0, 1.0, 4);

    settings
}

fn count_and_pair_runtime_settings() -> RuntimeSettings {
    let mut settings = RuntimeSettings::default();

    settings.quantities.insert(
        "particle_count".to_string(),
        QuantitySettings::Particle(ParticleQuantitySettings {
            pdgs: vec![1, -1],
            computation: count_computation(),
        }),
    );
    settings.quantities.insert(
        "particle_delta_r".to_string(),
        QuantitySettings::Particle(ParticleQuantitySettings {
            pdgs: vec![1, -1],
            computation: pair_computation(PairQuantity::DeltaR),
        }),
    );
    settings.quantities.insert(
        "jet_count".to_string(),
        QuantitySettings::Jet(JetQuantitySettings {
            clustering: JetClusteringSettings {
                dR: 0.4,
                clustered_pdgs: Some(vec![21]),
                ..JetClusteringSettings::default()
            },
            computation: count_computation(),
        }),
    );
    settings.quantities.insert(
        "jet_delta_r".to_string(),
        QuantitySettings::Jet(JetQuantitySettings {
            clustering: JetClusteringSettings {
                dR: 0.4,
                clustered_pdgs: Some(vec![21]),
                ..JetClusteringSettings::default()
            },
            computation: pair_computation(PairQuantity::DeltaR),
        }),
    );

    insert_histogram_observable(
        &mut settings,
        "particle_count_hist",
        "particle_count",
        0.0,
        4.0,
        4,
    );
    insert_histogram_observable(
        &mut settings,
        "particle_delta_r_hist",
        "particle_delta_r",
        0.0,
        4.0,
        4,
    );
    insert_histogram_observable(&mut settings, "jet_count_hist", "jet_count", 0.0, 4.0, 4);
    insert_histogram_observable(
        &mut settings,
        "jet_delta_r_hist",
        "jet_delta_r",
        0.0,
        4.0,
        4,
    );

    settings
}

fn leading_particle_quantity_runtime_settings(
    computation: QuantityComputationSettings,
    observable_name: &str,
    x_min: f64,
    x_max: f64,
    n_bins: usize,
) -> RuntimeSettings {
    let mut settings = RuntimeSettings::default();
    settings.quantities.insert(
        "ordered_particles".to_string(),
        QuantitySettings::Particle(ParticleQuantitySettings {
            pdgs: vec![1],
            computation,
        }),
    );
    settings.observables.insert(
        observable_name.to_string(),
        ObservableSettings {
            quantity: "ordered_particles".to_string(),
            entry_selection: EntrySelection::LeadingOnly,
            entry_index: 0,
            value_transform: ObservableValueTransform::Identity,
            phase: ObservablePhase::Real,
            misbinning_max_normalized_distance: None,
            histogram: HistogramSettings {
                x_min,
                x_max,
                n_bins,
                log_x_axis: false,
                log_y_axis: true,
                title: None,
                type_description: "AL".to_string(),
            },
        },
    );
    settings
}

fn leading_pair_quantity_runtime_settings(order: QuantityOrder) -> RuntimeSettings {
    let mut settings = RuntimeSettings::default();
    settings.quantities.insert(
        "particle_delta_r".to_string(),
        QuantitySettings::Particle(ParticleQuantitySettings {
            pdgs: vec![1],
            computation: pair_computation_with_order(order),
        }),
    );
    settings.observables.insert(
        "particle_delta_r_leading".to_string(),
        ObservableSettings {
            quantity: "particle_delta_r".to_string(),
            entry_selection: EntrySelection::LeadingOnly,
            entry_index: 0,
            value_transform: ObservableValueTransform::Identity,
            phase: ObservablePhase::Real,
            misbinning_max_normalized_distance: None,
            histogram: HistogramSettings {
                x_min: 0.0,
                x_max: 6.0,
                n_bins: 6,
                log_x_axis: false,
                log_y_axis: true,
                title: None,
                type_description: "AL".to_string(),
            },
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

    let mut accepted_event = accepted_event;
    let mut rejected_event = rejected_event;
    assert!(selector_runtime.process_event(&mut accepted_event));
    assert!(!selector_runtime.process_event(&mut rejected_event));

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

#[test]
fn event_group_contributions_fill_multiple_bins_with_one_sample_per_bin() -> Result<()> {
    let mut settings = runtime_settings();
    settings.selectors.clear();

    let mut runtime = EventProcessingRuntime::from_settings(&settings)?;
    let mut groups = single_group([
        make_event(0.25, 0.0, 0.0, 0.25, 1, (2.0, 0.0), []),
        make_event(0.75, 0.0, 0.0, 0.75, 1, (3.0, 0.0), []),
    ]);
    process_all_events(&mut runtime, &mut groups);
    runtime.process_event_groups(&groups);
    runtime.update_results(0);

    let histogram = &runtime.snapshot_bundle().histograms["pt_real"];
    assert_eq!(histogram.sample_count, 1);
    assert_eq!(histogram.bins[0].entry_count, 1);
    assert_eq!(histogram.bins[1].entry_count, 1);
    assert_close(histogram.bins[0].average(histogram.sample_count), 2.0);
    assert_close(histogram.bins[1].average(histogram.sample_count), 3.0);
    assert_eq!(
        histogram
            .bins
            .iter()
            .map(|bin| bin.entry_count)
            .sum::<usize>(),
        2
    );

    Ok(())
}

#[test]
fn grouped_histogram_updates_keep_same_bin_entries_correlated() -> Result<()> {
    let mut settings = runtime_settings();
    settings.selectors.clear();

    let mut runtime = EventProcessingRuntime::from_settings(&settings)?;

    let mut first_groups = single_group([
        make_event(0.25, 0.0, 0.0, 0.25, 1, (1.0, 0.0), []),
        make_event(0.35, 0.0, 0.0, 0.35, 1, (2.0, 0.0), []),
        make_event(0.75, 0.0, 0.0, 0.75, 1, (4.0, 0.0), []),
    ]);
    process_all_events(&mut runtime, &mut first_groups);
    runtime.process_event_groups(&first_groups);

    let mut second_groups = single_group([make_event(0.25, 0.0, 0.0, 0.25, 1, (5.0, 0.0), [])]);
    process_all_events(&mut runtime, &mut second_groups);
    runtime.process_event_groups(&second_groups);
    runtime.update_results(0);

    let histogram = &runtime.snapshot_bundle().histograms["pt_real"];
    assert_eq!(histogram.sample_count, 2);

    assert_eq!(histogram.bins[0].entry_count, 3);
    assert_close(histogram.bins[0].sum_weights, 8.0);
    assert_close(histogram.bins[0].sum_weights_squared, 34.0);
    assert_close(histogram.bins[0].average(histogram.sample_count), 4.0);
    assert_close(histogram.bins[0].error(histogram.sample_count), 1.0);

    assert_eq!(histogram.bins[1].entry_count, 1);
    assert_close(histogram.bins[1].sum_weights, 4.0);
    assert_close(histogram.bins[1].sum_weights_squared, 16.0);
    assert_close(histogram.bins[1].average(histogram.sample_count), 2.0);
    assert_close(histogram.bins[1].error(histogram.sample_count), 2.0);

    Ok(())
}

#[test]
fn jet_count_histograms_disable_misbinning_and_reject_misbinning_settings() -> Result<()> {
    let runtime = EventProcessingRuntime::from_settings(&jet_count_runtime_settings(false))?;
    let snapshot = runtime.snapshot_bundle();
    let histogram = snapshot
        .histograms
        .get("jet_count_hist")
        .expect("missing jet-count histogram");
    assert!(!histogram.supports_misbinning_mitigation);

    let error = match EventProcessingRuntime::from_settings(&jet_count_runtime_settings(true)) {
        Ok(_) => panic!("jet-count observable must reject misbinning mitigation"),
        Err(error) => error,
    };
    assert!(
        error
            .to_string()
            .contains("does not support misbinning mitigation")
    );
    Ok(())
}

#[test]
fn scalar_quantities_project_particle_and_jet_momenta() -> Result<()> {
    let mut runtime = EventProcessingRuntime::from_settings(&scalar_projection_runtime_settings())?;

    let mut particle_groups = singleton_groups(make_event(0.3, -0.4, 0.6, 0.5, 1, (1.0, 0.0), []));
    process_all_events(&mut runtime, &mut particle_groups);
    runtime.process_event_groups(&particle_groups);

    let mut jet_groups = singleton_groups(make_event(0.3, -0.4, 0.6, 0.5, 21, (1.0, 0.0), []));
    process_all_events(&mut runtime, &mut jet_groups);
    runtime.process_event_groups(&jet_groups);
    runtime.update_results(0);

    let snapshot = runtime.snapshot_bundle();

    assert_eq!(
        snapshot.histograms["particle_px_hist"].bins[2].entry_count,
        1
    );
    assert_eq!(
        snapshot.histograms["particle_py_hist"].bins[1].entry_count,
        1
    );
    assert_eq!(
        snapshot.histograms["particle_pz_hist"].bins[2].entry_count,
        1
    );
    assert_eq!(
        snapshot.histograms["particle_mass_hist"].bins[2].entry_count,
        1
    );
    assert_eq!(snapshot.histograms["jet_px_hist"].bins[2].entry_count, 1);
    assert_eq!(snapshot.histograms["jet_mass_hist"].bins[2].entry_count, 1);

    Ok(())
}

#[test]
fn count_and_pair_quantities_handle_exact_counts_and_delta_r() -> Result<()> {
    let mut runtime = EventProcessingRuntime::from_settings(&count_and_pair_runtime_settings())?;

    let outgoing = [
        (0.6, 0.0, 0.8, 1.0, 1),
        (0.0, 0.8, 0.6, 1.0, -1),
        (0.7, 0.0, 0.7, 1.0, 21),
        (0.0, -0.7, 0.7, 1.0, 21),
    ];
    let mut groups = singleton_groups(make_multi_event(&outgoing, (1.0, 0.0), []));
    process_all_events(&mut runtime, &mut groups);
    runtime.process_event_groups(&groups);
    runtime.update_results(0);

    let snapshot = runtime.snapshot_bundle();

    assert_eq!(
        snapshot.histograms["particle_count_hist"].bins[2].entry_count,
        1
    );
    assert_eq!(snapshot.histograms["jet_count_hist"].bins[2].entry_count, 1);
    assert_eq!(
        snapshot.histograms["particle_delta_r_hist"].bins[1].entry_count,
        1
    );
    assert_eq!(
        snapshot.histograms["jet_delta_r_hist"].bins[1].entry_count,
        1
    );

    Ok(())
}

#[test]
fn leading_only_uses_configured_particle_ordering() -> Result<()> {
    let outgoing = [
        (0.9, 0.0, 0.1, 1.2, 1),
        (0.1, 0.0, 2.0, 2.1, 1),
        (0.5, 0.0, 0.3, 0.8, 1),
    ];

    let scenarios = [
        (
            "particle_energy_default",
            scalar_computation(FilterQuantity::Energy),
            4usize,
        ),
        (
            "particle_energy_pt_desc",
            scalar_computation_with_ordering(
                FilterQuantity::Energy,
                QuantityOrdering::PT,
                QuantityOrder::Descending,
            ),
            2usize,
        ),
        (
            "particle_energy_abs_rapidity_desc",
            scalar_computation_with_ordering(
                FilterQuantity::Energy,
                QuantityOrdering::AbsRapidity,
                QuantityOrder::Descending,
            ),
            4usize,
        ),
        (
            "particle_energy_quantity_asc",
            scalar_computation_with_ordering(
                FilterQuantity::Energy,
                QuantityOrdering::Quantity,
                QuantityOrder::Ascending,
            ),
            1usize,
        ),
    ];

    for (observable_name, computation, expected_bin) in scenarios {
        let mut runtime = EventProcessingRuntime::from_settings(
            &leading_particle_quantity_runtime_settings(computation, observable_name, 0.0, 3.0, 6),
        )?;
        let mut groups = singleton_groups(make_multi_event(&outgoing, (1.0, 0.0), []));
        process_all_events(&mut runtime, &mut groups);
        runtime.process_event_groups(&groups);
        runtime.update_results(0);

        let histogram = &runtime.snapshot_bundle().histograms[observable_name];
        assert_eq!(histogram.sample_count, 1);
        assert_eq!(
            histogram.bins[expected_bin].entry_count, 1,
            "{observable_name}"
        );
        assert_eq!(
            histogram
                .bins
                .iter()
                .map(|bin| bin.entry_count)
                .sum::<usize>(),
            1,
            "{observable_name}"
        );
    }

    Ok(())
}

#[test]
fn leading_only_uses_pt_for_default_jet_ordering() -> Result<()> {
    let mut settings = RuntimeSettings::default();
    settings.quantities.insert(
        "jet_energy".to_string(),
        QuantitySettings::Jet(JetQuantitySettings {
            clustering: JetClusteringSettings {
                dR: 0.4,
                min_jpt: 0.0,
                clustered_pdgs: Some(vec![21]),
                ..JetClusteringSettings::default()
            },
            computation: scalar_computation(FilterQuantity::Energy),
        }),
    );
    settings.observables.insert(
        "jet_energy_leading".to_string(),
        ObservableSettings {
            quantity: "jet_energy".to_string(),
            entry_selection: EntrySelection::LeadingOnly,
            entry_index: 0,
            value_transform: ObservableValueTransform::Identity,
            phase: ObservablePhase::Real,
            misbinning_max_normalized_distance: None,
            histogram: HistogramSettings {
                x_min: 0.0,
                x_max: 3.0,
                n_bins: 6,
                log_x_axis: false,
                log_y_axis: true,
                title: None,
                type_description: "AL".to_string(),
            },
        },
    );

    let mut runtime = EventProcessingRuntime::from_settings(&settings)?;
    let outgoing = [(0.9, 0.0, 0.1, 1.2, 21), (0.1, 0.0, 2.0, 2.1, 21)];
    let mut groups = singleton_groups(make_multi_event(&outgoing, (1.0, 0.0), []));
    process_all_events(&mut runtime, &mut groups);
    runtime.process_event_groups(&groups);
    runtime.update_results(0);

    let histogram = &runtime.snapshot_bundle().histograms["jet_energy_leading"];
    assert_eq!(histogram.sample_count, 1);
    assert_eq!(histogram.bins[2].entry_count, 1);
    assert_eq!(
        histogram
            .bins
            .iter()
            .map(|bin| bin.entry_count)
            .sum::<usize>(),
        1
    );

    Ok(())
}

#[test]
fn pair_quantities_respect_order_direction_for_leading_only() -> Result<()> {
    let outgoing = [
        (1.0, 0.0, 0.0, 1.0, 1),
        (0.9950041652780258, 0.09983341664682815, 0.0, 1.0, 1),
        (-1.0, 0.0, 0.0, 1.0, 1),
    ];

    let mut descending_runtime = EventProcessingRuntime::from_settings(
        &leading_pair_quantity_runtime_settings(QuantityOrder::Descending),
    )?;
    let mut descending_groups = singleton_groups(make_multi_event(&outgoing, (1.0, 0.0), []));
    process_all_events(&mut descending_runtime, &mut descending_groups);
    descending_runtime.process_event_groups(&descending_groups);
    descending_runtime.update_results(0);
    let descending_histogram =
        &descending_runtime.snapshot_bundle().histograms["particle_delta_r_leading"];
    assert_eq!(descending_histogram.bins[3].entry_count, 1);

    let mut ascending_runtime = EventProcessingRuntime::from_settings(
        &leading_pair_quantity_runtime_settings(QuantityOrder::Ascending),
    )?;
    let mut ascending_groups = singleton_groups(make_multi_event(&outgoing, (1.0, 0.0), []));
    process_all_events(&mut ascending_runtime, &mut ascending_groups);
    ascending_runtime.process_event_groups(&ascending_groups);
    ascending_runtime.update_results(0);
    let ascending_histogram =
        &ascending_runtime.snapshot_bundle().histograms["particle_delta_r_leading"];
    assert_eq!(ascending_histogram.bins[0].entry_count, 1);

    Ok(())
}
