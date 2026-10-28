//! Test fault injection for process payloads and exports.

use std::cell::Cell;
use std::fs;
use std::path::Path;
use std::thread_local;

use std::io::Result as IoResult;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SaveStep {
    Integrand,
    Metadata,
}

thread_local! {
    static SAVE_FAULT: Cell<Option<SaveStep>> = const { Cell::new(None) };
}

pub(super) fn save_step(step: SaveStep, path: &Path) -> IoResult<()> {
    if SAVE_FAULT.get() == Some(step) {
        match step {
            SaveStep::Integrand => fs::write(path, b"test-only integrand write obstruction")?,
            SaveStep::Metadata => fs::create_dir(path)?,
        }
    }
    Ok(())
}

pub(super) struct SaveFaultGuard;

impl SaveFaultGuard {
    pub(super) fn new(step: SaveStep) -> Self {
        assert!(SAVE_FAULT.replace(Some(step)).is_none());
        Self
    }
}

impl Drop for SaveFaultGuard {
    fn drop(&mut self) {
        SAVE_FAULT.set(None);
    }
}

mod tests {
    use super::*;

    use alphal00p_persistence::{OutputPath, OutputPathError};
    use std::env;
    use std::error::Error as StdError;
    use std::io::{Cursor, Error as IoError, ErrorKind};
    #[cfg(unix)]
    use std::os::unix::fs as unix_fs;
    use std::path::PathBuf;
    use std::process;
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::time::{SystemTime, UNIX_EPOCH};

    use crate::GammaLoopContextContainer;
    use crate::graph::Graph;
    use crate::initialisation::test_initialise;
    use crate::integrands::process::ProcessIntegrand;
    use crate::model::Model;
    use crate::processes::{
        DotExportSettings, INTEGRAND_BINARY_FILE, INTEGRAND_SETTINGS_FILE, IntegrandPayload,
        IntegrandSnapshotEntry, PROCESSES_DIRECTORY, Process, ProcessCollection, ProcessDefinition,
        ProcessList, ProcessLoadIntegrandSelector, ProcessLoadSelection, ProcessPayloadError,
        ProcessSnapshotError,
    };
    use crate::settings::{GlobalSettings, RuntimeSettings};
    use crate::utils;
    use crate::utils::serde_utils::{ShowDefaultsGuard, SmartSerde};
    use bincode::config;
    use bincode::error::DecodeError;
    use feynkit_generator::GenerationType;
    use rayon::ThreadPoolBuilder;
    use symbolica::evaluate::FunctionMap;
    use symbolica::state::State as SymbolicaState;
    use toml::Value;
    use toml::de::Error as SettingsDecodeError;

    static NEXT_ROOT: AtomicU64 = AtomicU64::new(0);

    struct TestRoot(PathBuf);

    impl TestRoot {
        fn new() -> Self {
            let nonce = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos();
            let sequence = NEXT_ROOT.fetch_add(1, Ordering::Relaxed);
            let root = env::temp_dir().join(format!(
                "process-output-{}-{nonce}-{sequence}",
                process::id()
            ));
            fs::create_dir(&root).unwrap();
            Self(root)
        }
    }

    impl Drop for TestRoot {
        fn drop(&mut self) {
            fs::remove_dir_all(&self.0).unwrap();
        }
    }

    fn fixture(kind: GenerationType) -> (Model, ProcessList) {
        test_initialise().unwrap();
        let (model, source) = match kind {
            GenerationType::Amplitude => (
                utils::load_generic_model("scalars"),
                include_str!("../../../../tests/resources/graphs/scalar_bubble.dot"),
            ),
            GenerationType::CrossSection => (
                Model::from_json(include_str!(
                    "../../../../assets/models/json/scalars/scalars_2p_3p.json"
                ))
                .unwrap(),
                include_str!(
                    "../../../../tests/resources/graphs/mass_approach_scalar_self_energy.dot"
                ),
            ),
        };
        let mut graphs = Graph::import_dot(source, &model).unwrap();
        graphs.truncate(1);
        let definition = ProcessDefinition::from_graph_list(&graphs, kind, &model).unwrap();
        let process = Process::from_graph_list(
            "process".into(),
            "child".into(),
            graphs,
            kind,
            Some(definition),
            None,
            &model,
        )
        .unwrap();
        let mut list = ProcessList::new();
        list.add_process(process);
        (model, list)
    }

    // ParamBuilder::fn_map is reconstructed on decode; its recorded `reps` remain
    // part of this comparison. Do not normalize any other persisted field.
    fn persisted_data_without_rebuilt_maps(integrand: &mut ProcessIntegrand) -> Vec<u8> {
        match integrand {
            ProcessIntegrand::Amplitude(integrand) => {
                for term in &mut integrand.data.graph_terms {
                    for builder in [
                        &mut term.param_builder,
                        &mut term.graph.param_builder,
                        &mut term.multi_channeling_setup.graph.param_builder,
                    ] {
                        builder.fn_map = FunctionMap::new();
                    }
                }
                bincode::encode_to_vec(&integrand.data, config::standard()).unwrap()
            }
            ProcessIntegrand::CrossSection(integrand) => {
                for term in &mut integrand.data.graph_terms {
                    for builder in [
                        &mut term.param_builder,
                        &mut term.graph.param_builder,
                        &mut term.multi_channeling_setup.graph.param_builder,
                    ] {
                        builder.fn_map = FunctionMap::new();
                    }
                }
                bincode::encode_to_vec(&integrand.data, config::standard()).unwrap()
            }
        }
    }

    #[test]
    fn process_snapshot_write_failures_retain_generated_integrands() {
        for kind in [GenerationType::Amplitude, GenerationType::CrossSection] {
            let (model, mut list) = fixture(kind);
            let pool = ThreadPoolBuilder::new().num_threads(1).build().unwrap();
            let mut settings = GlobalSettings::default();
            settings.generation.evaluator.compile = false;
            settings.generation.uv.subtract_uv = false;
            settings.generation.uv.generate_integrated = false;
            settings.generation.uv.softct = false;
            settings.generation.threshold_subtraction.enable_thresholds = false;
            let runtime = RuntimeSettings::default();
            list.preprocess(&model, &settings, &(&runtime).into(), &pool)
                .unwrap();
            list.generate_integrands(&model, &settings, (&runtime).into(), &pool)
                .unwrap();
            let before = list.processes[0]
                .get_integrand("child")
                .unwrap()
                .require_generated()
                .unwrap()
                .resume_fingerprint()
                .unwrap();
            let before_settings = toml::to_string(
                list.processes[0]
                    .get_integrand("child")
                    .unwrap()
                    .require_generated()
                    .unwrap()
                    .get_settings(),
            )
            .unwrap();
            for step in [SaveStep::Integrand, SaveStep::Metadata] {
                let root = TestRoot::new();
                let output = OutputPath::new(&root.0).unwrap();
                let fault = SaveFaultGuard::new(step);
                let error = list.write_snapshot(&output).err().unwrap();
                assert!(
                    matches!(error, ProcessSnapshotError::Payload { .. }),
                    "{kind}/{step:?}: {error:?}"
                );
                let mut source = error.source();
                let mut original_io = false;
                while let Some(error) = source {
                    original_io |= error.is::<IoError>();
                    source = error.source();
                }
                assert!(original_io, "{kind}/{step:?}: {error:?}");
                drop(fault);
                let integrand = list.processes[0].get_integrand("child").unwrap();
                let integrand = integrand.require_generated().unwrap();
                assert_eq!(
                    integrand.resume_fingerprint().unwrap(),
                    before,
                    "{kind}/{step:?}"
                );
                assert_eq!(
                    toml::to_string(integrand.get_settings()).unwrap(),
                    before_settings
                );
            }
            let root = TestRoot::new();
            let output = OutputPath::new(&root.0).unwrap();
            let index = list.write_snapshot(&output).unwrap();
            assert_eq!(
                index.processes[0].children[0].payload,
                IntegrandPayload::Generated
            );
            let definition_path = root
                .0
                .join(PROCESSES_DIRECTORY)
                .join(index.processes[0].kind_directory())
                .join("process/def.bin");
            let old_definition = fs::read(&definition_path).unwrap();
            assert!(list.write_snapshot(&output).is_err());
            assert_eq!(fs::read(definition_path).unwrap(), old_definition);

            let child = output
                .join(PROCESSES_DIRECTORY)
                .unwrap()
                .join(index.processes[0].kind_directory())
                .unwrap()
                .join("process")
                .unwrap()
                .join("child")
                .unwrap();
            let settings_path = child
                .join_relative(Path::new(INTEGRAND_SETTINGS_FILE))
                .unwrap();
            let settings_text = fs::read_to_string(settings_path.as_path()).unwrap();
            let schema = runtime.schema_path().unwrap();
            assert_eq!(
                settings_text.lines().next().unwrap(),
                format!("#:schema {}", schema.display())
            );
            let integrand = list.get_integrand(0, "child").unwrap();
            let settings = integrand.require_generated().unwrap().get_settings();
            let expected_settings = {
                let _guard = ShowDefaultsGuard::new(true);
                Value::try_from(settings).unwrap()
            };
            let persisted_settings: Value = toml::from_str(&settings_text).unwrap();
            assert_eq!(persisted_settings, expected_settings);
            let required_files = index.processes[0].required_files();
            assert!(!required_files.is_empty());
            for relative in required_files {
                assert!(root.0.join(relative).is_file());
            }

            let mut symbols = Vec::new();
            SymbolicaState::export(&mut symbols).unwrap();
            let state_map = SymbolicaState::import(&mut Cursor::new(symbols), None).unwrap();
            let context = GammaLoopContextContainer {
                model: &model,
                state_map: &state_map,
            };
            let mut loaded = ProcessList::load(&output, &index, context).unwrap();
            assert_eq!(loaded.processes.len(), list.processes.len());
            for (actual, expected) in loaded.processes.iter().zip(&list.processes) {
                assert_eq!(actual.definition, expected.definition);
                assert_eq!(
                    actual.snapshot_entry().unwrap(),
                    expected.snapshot_entry().unwrap(),
                );
            }
            let original = list.get_integrand(0, "child").unwrap();
            let original = original.require_generated().unwrap();
            let restored = loaded.get_integrand(0, "child").unwrap();
            let restored = restored.require_generated().unwrap();
            let restored_settings = {
                let _guard = ShowDefaultsGuard::new(true);
                Value::try_from(restored.get_settings()).unwrap()
            };
            assert_eq!(restored_settings, expected_settings);
            assert_eq!(restored.graph_count(), original.graph_count());
            assert_eq!(restored.graph_group_count(), original.graph_group_count());
            assert_eq!(
                restored.graph_group_master_names(),
                original.graph_group_master_names(),
            );

            // These test-owned integrands are no longer used after this comparison.
            // Compare the complete wire payload, including private and future fields,
            // while excluding only the hash-ordered maps rebuilt from recorded entries.
            let original = list.get_integrand_mut(0, "child").unwrap();
            let restored = loaded.get_integrand_mut(0, "child").unwrap();
            assert_eq!(
                restored.kind_name(),
                original.kind_name(),
                "saved {kind} integrand changed its process kind",
            );
            assert_eq!(
                persisted_data_without_rebuilt_maps(restored),
                persisted_data_without_rebuilt_maps(original),
                "saved {kind} integrand changed its persisted fields",
            );

            fs::remove_file(settings_path.as_path()).unwrap();
            let error = ProcessList::load(&output, &index, context).err().unwrap();
            let ProcessSnapshotError::Payload(failure) = error else {
                panic!("expected a missing-settings payload failure: {error:?}");
            };
            assert!(matches!(
                &failure.source,
                ProcessPayloadError::Read { path, source }
                    if path == settings_path.as_path() && source.kind() == ErrorKind::NotFound
            ));
            assert!(failure.source.source().unwrap().is::<IoError>());

            fs::write(settings_path.as_path(), "[invalid").unwrap();
            let error = ProcessList::load(&output, &index, context).err().unwrap();
            let ProcessSnapshotError::Payload(failure) = error else {
                panic!("expected malformed-settings payload failure: {error:?}");
            };
            assert!(matches!(
                &failure.source,
                ProcessPayloadError::SettingsDecode { path, .. } if path == settings_path.as_path()
            ));
            assert!(failure.source.source().unwrap().is::<SettingsDecodeError>());

            let binary_path = child
                .join_relative(Path::new(INTEGRAND_BINARY_FILE))
                .unwrap();
            fs::write(binary_path.as_path(), []).unwrap();
            let error = ProcessList::load(&output, &index, context).err().unwrap();
            let ProcessSnapshotError::Payload(failure) = error else {
                panic!("expected malformed-binary payload failure: {error:?}");
            };
            assert!(matches!(
                &failure.source,
                ProcessPayloadError::Decode { path, .. } if path == binary_path.as_path()
            ));
            assert!(failure.source.source().unwrap().is::<DecodeError>());
        }
    }

    #[test]
    fn process_output_rejects_escaping_names_and_preserves_directory_targets() {
        let temporary = TestRoot::new();
        let outside = temporary.0.join("outside");
        fs::create_dir(&outside).unwrap();
        fs::write(outside.join("sentinel"), b"preserve").unwrap();
        let root = temporary.0.join("output");
        let output = OutputPath::new(&root).unwrap();
        fs::create_dir(&root).unwrap();
        for bad in [
            "",
            ".",
            "..",
            "../outside",
            "a/b",
            "a\\b",
            outside.to_str().unwrap(),
        ] {
            let (_, mut list) = fixture(GenerationType::Amplitude);
            list.processes[0].definition.folder_name = bad.into();
            assert!(
                list.export_dot(&output, &DotExportSettings::default())
                    .is_err(),
                "process {bad:?}"
            );
            let (_, mut list) = fixture(GenerationType::Amplitude);
            if let ProcessCollection::Amplitudes(amplitudes) = &mut list.processes[0].collection {
                amplitudes.get_mut("child").unwrap().name = bad.into();
            }
            assert!(
                list.export_dot(&output, &DotExportSettings::default())
                    .is_err(),
                "child {bad:?}"
            );
            let (_, mut list) = fixture(GenerationType::Amplitude);
            if let ProcessCollection::Amplitudes(amplitudes) = &mut list.processes[0].collection {
                amplitudes.get_mut("child").unwrap().graphs[0].graph.name = bad.into();
            }
            assert!(
                list.export_dot(&output, &DotExportSettings::default())
                    .is_err(),
                "graph {bad:?}"
            );
            assert_eq!(fs::read(outside.join("sentinel")).unwrap(), b"preserve");
        }
        let (_, mut list) = fixture(GenerationType::Amplitude);
        if let ProcessCollection::Amplitudes(amplitudes) = &mut list.processes[0].collection {
            amplitudes.get_mut("child").unwrap().graphs[0].graph.name = "victim".into();
        }
        let target = root
            .join(PROCESSES_DIRECTORY)
            .join("amplitudes/process/child/victim.dot");
        fs::create_dir_all(&target).unwrap();
        fs::write(target.join("sentinel"), b"preserve directory").unwrap();
        assert!(
            list.export_dot(&output, &DotExportSettings::default())
                .is_err()
        );
        assert_eq!(
            fs::read(target.join("sentinel")).unwrap(),
            b"preserve directory"
        );
        assert_eq!(fs::read(outside.join("sentinel")).unwrap(), b"preserve");
        let outside_target = outside.join("victim.dot");
        fs::create_dir(&outside_target).unwrap();
        fs::write(outside_target.join("sentinel"), b"outside directory").unwrap();
        for escaped_graph in [
            outside.join("victim").to_string_lossy().into_owned(),
            "../../../../../outside/victim".to_owned(),
        ] {
            if let ProcessCollection::Amplitudes(amplitudes) = &mut list.processes[0].collection {
                amplitudes.get_mut("child").unwrap().graphs[0].graph.name = escaped_graph;
            }
            assert!(
                list.export_dot(&output, &DotExportSettings::default())
                    .is_err()
            );
            assert_eq!(
                fs::read(outside_target.join("sentinel")).unwrap(),
                b"outside directory"
            );
        }
    }

    #[cfg(unix)]
    #[test]
    fn process_output_refuses_symlinked_parents_and_files() {
        let temporary = TestRoot::new();
        let outside = temporary.0.join("outside");
        fs::create_dir(&outside).unwrap();
        fs::write(outside.join("sentinel"), b"preserve").unwrap();
        let root = temporary.0.join("output");
        let output = OutputPath::new(&root).unwrap();
        fs::create_dir(&root).unwrap();
        unix_fs::symlink(&outside, root.join(PROCESSES_DIRECTORY)).unwrap();
        let (_, list) = fixture(GenerationType::Amplitude);
        assert!(
            list.export_dot(&output, &DotExportSettings::default())
                .is_err()
        );
        fs::remove_file(root.join(PROCESSES_DIRECTORY)).unwrap();
        let link = root.join("file.dot");
        unix_fs::symlink(outside.join("sentinel"), &link).unwrap();
        assert!(matches!(
            output.join("file.dot"),
            Err(OutputPathError::Symlink(_))
        ));
        assert_eq!(fs::read(outside.join("sentinel")).unwrap(), b"preserve");
        assert_eq!(fs::read_dir(&outside).unwrap().count(), 1);
    }

    #[test]
    fn process_snapshot_rejects_inconsistent_ids_and_names_before_writing() {
        let root = TestRoot::new();
        let output = OutputPath::new(&root.0).unwrap();
        let (_, mut list) = fixture(GenerationType::Amplitude);
        list.processes[0].definition.process_id = 7;
        assert!(matches!(
            list.write_snapshot(&output),
            Err(ProcessSnapshotError::InvalidProcessId {
                expected: 0,
                actual: 7,
                ..
            })
        ));
        assert_eq!(fs::read_dir(&root.0).unwrap().count(), 0);
        list.processes[0].definition.process_id = 0;
        if let ProcessCollection::Amplitudes(amplitudes) = &mut list.processes[0].collection {
            amplitudes.get_mut("child").unwrap().name = "other".into();
        }
        let error = list.write_snapshot(&output).unwrap_err();
        assert!(matches!(
            error,
            ProcessSnapshotError::ChildNameMismatch { expected, actual, .. }
                if expected == "child" && actual == "other"
        ));
        assert_eq!(fs::read_dir(&root.0).unwrap().count(), 0);
    }

    #[test]
    fn process_snapshot_errors_distinguish_membership_and_preserve_decode_sources() {
        let root = TestRoot::new();
        let output = OutputPath::new(&root.0).unwrap();
        let (model, mut list) = fixture(GenerationType::Amplitude);
        let index = list.write_snapshot(&output).unwrap();
        let mut invalid = index.clone();
        invalid.processes.push(invalid.processes[0].clone());
        let error = invalid.validate().unwrap_err();
        assert!(matches!(
            error,
            ProcessSnapshotError::DuplicateProcess { process } if process == "process"
        ));
        invalid = index.clone();
        invalid.processes[0].children.push(IntegrandSnapshotEntry {
            name: "child".into(),
            payload: IntegrandPayload::MetadataOnly,
        });
        let error = invalid.validate().unwrap_err();
        assert!(matches!(
            error,
            ProcessSnapshotError::InvalidChildOrder { previous, child, .. } if previous == child
        ));
        invalid = index.clone();
        invalid.processes[0].folder_name = "../outside".into();
        assert!(matches!(
            invalid.validate(),
            Err(ProcessSnapshotError::Output(
                OutputPathError::InvalidComponent(_)
            ))
        ));

        let mut symbols = Vec::new();
        SymbolicaState::export(&mut symbols).unwrap();
        let state_map = SymbolicaState::import(&mut Cursor::new(symbols), None).unwrap();
        let context = GammaLoopContextContainer {
            model: &model,
            state_map: &state_map,
        };
        let mut selection = ProcessLoadSelection::default();
        selection.process_names.insert("missing".into());
        let error = ProcessList::load_with_selection(&output, &index, context, Some(&selection))
            .err()
            .unwrap();
        assert!(matches!(
            error,
            ProcessSnapshotError::UnknownProcess { process, available }
                if process == "missing" && available == ["process"]
        ));
        selection.process_names.clear();
        selection
            .integrand_selectors
            .push(ProcessLoadIntegrandSelector::Name("missing".into()));
        assert!(matches!(
            ProcessList::load_with_selection(&output, &index, context, Some(&selection)),
            Err(ProcessSnapshotError::UnknownIntegrand { .. })
        ));

        let definition = root.0.join("processes/amplitudes/process/def.bin");
        fs::write(&definition, []).unwrap();
        let error = ProcessList::load(&output, &index, context).err().unwrap();
        assert!(matches!(
            &error,
            ProcessSnapshotError::DefinitionDecode { path, .. } if path == &definition
        ));
        assert!(error.source().unwrap().is::<DecodeError>());
        fs::remove_file(&definition).unwrap();
        let error = ProcessList::load(&output, &index, context).err().unwrap();
        assert!(matches!(
            &error,
            ProcessSnapshotError::Read { path, source }
                if path == &definition && source.kind() == ErrorKind::NotFound
        ));
        assert!(error.source().unwrap().is::<IoError>());
    }

    #[test]
    fn graph_list_constructor_preserves_the_collection_kind_in_its_definition() {
        let model = Model::empty("fixture");
        let process = Process::from_graph_list(
            "process".into(),
            "child".into(),
            Vec::new(),
            GenerationType::CrossSection,
            None,
            None,
            &model,
        )
        .unwrap();
        assert_eq!(
            process.definition.generation_type,
            GenerationType::CrossSection
        );
        assert!(
            Process::from_graph_list(
                "process".into(),
                "child".into(),
                Vec::new(),
                GenerationType::CrossSection,
                Some(ProcessDefinition::default()),
                None,
                &model
            )
            .is_err()
        );
    }
}
