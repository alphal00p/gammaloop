//! Filesystem boundaries for process payloads and exports.

use std::fs::{self, File};
use std::io::{Error as IoError, ErrorKind, Write};
use std::path::{Component, Path, PathBuf};

#[cfg(test)]
use color_eyre::Result as ReportResult;
#[cfg(test)]
use std::cell::Cell;
#[cfg(test)]
use std::thread_local;
use thiserror::Error;

/// A rejected process output path. See [`OutputPath`] for the filesystem contract.
#[derive(Debug, Error)]
pub enum OutputPathError {
    #[error("Invalid output component {0:?}: expected one nonempty filename")]
    InvalidComponent(String),
    #[error("Refusing symlink in process output path {0}")]
    Symlink(PathBuf),
    #[error("Expected a regular file at {0}")]
    NotFile(PathBuf),
    #[error("Process output I/O at {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: IoError,
    },
}

/// An output path whose existing ancestors contain no symlinks.
///
/// Names are preserved verbatim and must be single components; this never applies
/// a lossy slug conversion. Existing symlinks, including ancestors of the supplied
/// root, are refused. The caller must exclude concurrent hostile filesystem changes;
/// these checks do not provide descriptor-relative protection against symlink races.
///
/// ```no_run
/// # use gammalooprs::processes::OutputPath;
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let root = OutputPath::new("export".as_ref())?;
/// root.ensure_directory()?;
/// let child = root.join("process")?;
/// child.create_directory()?;
/// child.join("graph.dot")?.create_file()?;
/// child.join("template.typ")?.write(b"complete template")?;
/// # Ok(()) }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OutputPath(PathBuf);

impl OutputPath {
    /// Check a root without creating it. See the example on [`Self`].
    pub fn new(path: &Path) -> Result<Self, OutputPathError> {
        let output = Self(path.to_path_buf());
        output.check()?;
        Ok(output)
    }

    /// Reject separators, absolute paths, dot components and empty names.
    /// See the example on [`Self`].
    pub fn validate_component(name: &str) -> Result<(), OutputPathError> {
        if name.is_empty()
            || name.contains(['/', '\\', '\0'])
            || !matches!(
                Path::new(name).components().next(),
                Some(Component::Normal(_))
            )
            || Path::new(name).components().count() != 1
        {
            return Err(OutputPathError::InvalidComponent(name.to_owned()));
        }
        Ok(())
    }

    /// Append one validated, unchanged component. See the example on [`Self`].
    pub fn join(&self, name: &str) -> Result<Self, OutputPathError> {
        Self::validate_component(name)?;
        Self::new(&self.0.join(name))
    }

    /// Borrow the checked path. See the example on [`Self`].
    pub fn as_path(&self) -> &Path {
        &self.0
    }

    /// Recheck all existing ancestors, without writing. See the example on [`Self`].
    pub fn check(&self) -> Result<(), OutputPathError> {
        for path in self.0.ancestors() {
            if path.as_os_str().is_empty() {
                continue;
            }
            match fs::symlink_metadata(path) {
                Ok(metadata) if metadata.file_type().is_symlink() => {
                    return Err(OutputPathError::Symlink(path.to_owned()));
                }
                Ok(_) => {}
                Err(error) if error.kind() == ErrorKind::NotFound => {}
                Err(source) => return Err(Self::io(path, source)),
            }
        }
        Ok(())
    }

    // Lower compiler/export serializers own their filenames. Refuse symlinks in
    // their existing subtree before handing them a mutable output directory.
    pub(super) fn check_tree(&self) -> Result<(), OutputPathError> {
        self.check()?;
        let mut pending = vec![self.0.clone()];
        while let Some(path) = pending.pop() {
            let metadata = match fs::symlink_metadata(&path) {
                Ok(metadata) => metadata,
                Err(error) if error.kind() == ErrorKind::NotFound => continue,
                Err(source) => return Err(Self::io(&path, source)),
            };
            if metadata.file_type().is_symlink() {
                return Err(OutputPathError::Symlink(path));
            }
            if metadata.is_dir() {
                for entry in fs::read_dir(&path).map_err(|source| Self::io(&path, source))? {
                    pending.push(entry.map_err(|source| Self::io(&path, source))?.path());
                }
            }
        }
        Ok(())
    }

    /// Create directories, refusing existing symlinked parents. See [`Self`].
    pub fn ensure_directory(&self) -> Result<(), OutputPathError> {
        self.check()?;
        fs::create_dir_all(&self.0).map_err(|source| Self::io(&self.0, source))
    }

    /// Exclusively create a staging directory. See the example on [`Self`].
    pub fn create_directory(&self) -> Result<(), OutputPathError> {
        self.check()?;
        fs::create_dir(&self.0).map_err(|source| Self::io(&self.0, source))
    }

    /// Create a fresh file, replacing a regular-file entry; never remove a directory.
    /// Existing hard links keep their contents. See [`Self`].
    pub fn create_file(&self) -> Result<File, OutputPathError> {
        self.check()?;
        match fs::symlink_metadata(&self.0) {
            Ok(metadata) if !metadata.is_file() => {
                return Err(OutputPathError::NotFile(self.0.clone()));
            }
            Ok(_) => {
                fs::remove_file(&self.0).map_err(|source| Self::io(&self.0, source))?;
            }
            Err(error) if error.kind() == ErrorKind::NotFound => {}
            Err(source) => return Err(Self::io(&self.0, source)),
        }
        File::create_new(&self.0).map_err(|source| Self::io(&self.0, source))
    }

    /// Write and flush a checked regular-file target. See the example on [`Self`].
    pub fn write(&self, contents: &[u8]) -> Result<(), OutputPathError> {
        let mut file = self.create_file()?;
        file.write_all(contents)
            .map_err(|source| Self::io(&self.0, source))?;
        file.flush().map_err(|source| Self::io(&self.0, source))
    }

    fn io(path: &Path, source: IoError) -> OutputPathError {
        OutputPathError::Io {
            path: path.to_owned(),
            source,
        }
    }
}

#[cfg(test)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SaveStep {
    Integrand,
    Metadata,
}

#[cfg(test)]
thread_local! {
    static SAVE_FAULT: Cell<Option<SaveStep>> = const { Cell::new(None) };
}

#[cfg(test)]
pub(super) fn save_step(step: SaveStep, path: &Path) -> ReportResult<()> {
    if SAVE_FAULT.get() == Some(step) {
        match step {
            SaveStep::Integrand => fs::write(path, b"test-only integrand write obstruction")?,
            SaveStep::Metadata => fs::create_dir(path)?,
        }
    }
    Ok(())
}

#[cfg(test)]
pub(super) struct SaveFaultGuard;

#[cfg(test)]
impl SaveFaultGuard {
    pub(super) fn new(step: SaveStep) -> Self {
        assert!(SAVE_FAULT.replace(Some(step)).is_none());
        Self
    }
}

#[cfg(test)]
impl Drop for SaveFaultGuard {
    fn drop(&mut self) {
        SAVE_FAULT.set(None);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::env;
    use std::error::Error as StdError;
    use std::io::Cursor;
    #[cfg(unix)]
    use std::os::unix::fs as unix_fs;
    use std::process;
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::time::{SystemTime, UNIX_EPOCH};

    use crate::GammaLoopContextContainer;
    use crate::graph::Graph;
    use crate::initialisation::test_initialise;
    use crate::model::Model;
    use crate::processes::{
        DotExportSettings, IntegrandSnapshotEntry, Process, ProcessCollection, ProcessDefinition,
        ProcessList, ProcessLoadIntegrandSelector, ProcessLoadSelection, ProcessSnapshotError,
    };
    use crate::settings::{GlobalSettings, RuntimeSettings};
    use crate::utils;
    use bincode::error::DecodeError;
    use feynkit_generator::GenerationType;
    use rayon::ThreadPoolBuilder;
    use symbolica::state::State as SymbolicaState;

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
                let fault = SaveFaultGuard::new(step);
                let error = list.write_snapshot(&root.0).err().unwrap();
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
            let index = list.write_snapshot(&root.0).unwrap();
            assert!(index.processes[0].children[0].generated);
            let definition_path = root
                .0
                .join("processes")
                .join(index.processes[0].kind_directory())
                .join("process/def.bin");
            let old_definition = fs::read(&definition_path).unwrap();
            assert!(list.write_snapshot(&root.0).is_err());
            assert_eq!(fs::read(definition_path).unwrap(), old_definition);
        }
    }

    #[test]
    fn process_output_rejects_escaping_names_and_preserves_directory_targets() {
        let temporary = TestRoot::new();
        let outside = temporary.0.join("outside");
        fs::create_dir(&outside).unwrap();
        fs::write(outside.join("sentinel"), b"preserve").unwrap();
        let root = temporary.0.join("output");
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
                list.export_dot(&root, &DotExportSettings::default())
                    .is_err(),
                "process {bad:?}"
            );
            let (_, mut list) = fixture(GenerationType::Amplitude);
            if let ProcessCollection::Amplitudes(amplitudes) = &mut list.processes[0].collection {
                amplitudes.get_mut("child").unwrap().name = bad.into();
            }
            assert!(
                list.export_dot(&root, &DotExportSettings::default())
                    .is_err(),
                "child {bad:?}"
            );
            let (_, mut list) = fixture(GenerationType::Amplitude);
            if let ProcessCollection::Amplitudes(amplitudes) = &mut list.processes[0].collection {
                amplitudes.get_mut("child").unwrap().graphs[0].graph.name = bad.into();
            }
            assert!(
                list.export_dot(&root, &DotExportSettings::default())
                    .is_err(),
                "graph {bad:?}"
            );
            assert_eq!(fs::read(outside.join("sentinel")).unwrap(), b"preserve");
        }
        let (_, mut list) = fixture(GenerationType::Amplitude);
        if let ProcessCollection::Amplitudes(amplitudes) = &mut list.processes[0].collection {
            amplitudes.get_mut("child").unwrap().graphs[0].graph.name = "victim".into();
        }
        let target = root.join("processes/amplitudes/process/child/victim.dot");
        fs::create_dir_all(&target).unwrap();
        fs::write(target.join("sentinel"), b"preserve directory").unwrap();
        assert!(
            list.export_dot(&root, &DotExportSettings::default())
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
                list.export_dot(&root, &DotExportSettings::default())
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
        fs::create_dir(&root).unwrap();
        unix_fs::symlink(&outside, root.join("processes")).unwrap();
        let (_, list) = fixture(GenerationType::Amplitude);
        assert!(
            list.export_dot(&root, &DotExportSettings::default())
                .is_err()
        );
        fs::remove_file(root.join("processes")).unwrap();
        let link = root.join("file.dot");
        unix_fs::symlink(outside.join("sentinel"), &link).unwrap();
        assert!(OutputPath::new(&link).is_err());
        assert_eq!(fs::read(outside.join("sentinel")).unwrap(), b"preserve");
        assert_eq!(fs::read_dir(&outside).unwrap().count(), 1);
    }

    #[test]
    fn process_snapshot_rejects_inconsistent_ids_and_names_before_writing() {
        let root = TestRoot::new();
        let (_, mut list) = fixture(GenerationType::Amplitude);
        list.processes[0].definition.process_id = 7;
        assert!(matches!(
            list.write_snapshot(&root.0),
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
        assert!(
            matches!(list.write_snapshot(&root.0), Err(ProcessSnapshotError::ChildNameMismatch { expected, actual, .. }) if expected == "child" && actual == "other")
        );
        assert_eq!(fs::read_dir(&root.0).unwrap().count(), 0);
    }

    #[test]
    fn process_snapshot_errors_distinguish_membership_and_preserve_decode_sources() {
        let root = TestRoot::new();
        let (model, mut list) = fixture(GenerationType::Amplitude);
        let index = list.write_snapshot(&root.0).unwrap();
        let mut invalid = index.clone();
        invalid.processes.push(invalid.processes[0].clone());
        assert!(
            matches!(invalid.validate(), Err(ProcessSnapshotError::DuplicateProcess { process }) if process == "process")
        );
        invalid = index.clone();
        invalid.processes[0].children.push(IntegrandSnapshotEntry {
            name: "child".into(),
            generated: false,
        });
        assert!(
            matches!(invalid.validate(), Err(ProcessSnapshotError::InvalidChildOrder { previous, child, .. }) if previous == child)
        );
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
        assert!(
            matches!(ProcessList::load_with_selection(&root.0, &index, context, Some(&selection)), Err(ProcessSnapshotError::UnknownProcess { process, available }) if process == "missing" && available == ["process"])
        );
        selection.process_names.clear();
        selection
            .integrand_selectors
            .push(ProcessLoadIntegrandSelector::Name("missing".into()));
        assert!(matches!(
            ProcessList::load_with_selection(&root.0, &index, context, Some(&selection)),
            Err(ProcessSnapshotError::UnknownIntegrand { .. })
        ));

        let definition = root.0.join("processes/amplitudes/process/def.bin");
        fs::write(&definition, []).unwrap();
        let error = ProcessList::load(&root.0, &index, context).err().unwrap();
        assert!(
            matches!(&error, ProcessSnapshotError::DefinitionDecode { path, .. } if path == &definition)
        );
        assert!(error.source().unwrap().is::<DecodeError>());
        fs::remove_file(&definition).unwrap();
        let error = ProcessList::load(&root.0, &index, context).err().unwrap();
        assert!(
            matches!(&error, ProcessSnapshotError::Read { path, source } if path == &definition && source.kind() == ErrorKind::NotFound)
        );
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
