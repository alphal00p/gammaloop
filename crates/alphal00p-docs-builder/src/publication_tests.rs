use std::fs;
use std::path::{Path, PathBuf};

use alphal00p_persistence::{DirectoryPublication, DirectoryPublicationError};
use tempfile::TempDir;

use super::{
    BuildChannel, CliTypstRenderer, GENERATED_OUTPUT_MARKER, GENERATED_OUTPUT_SIGNATURE,
    GeneratedTreePublication, SiteBuilder, clear_stale_staging, copy_tree, directory_files,
};

const SNAPSHOT: &str = "products/linnet/snapshots/v1/index.html";

pub(super) fn initialize_publication_lock(destination: &Path) -> PathBuf {
    let owner = DirectoryPublication::new(destination).unwrap();
    let error = DirectoryPublication::new(destination).unwrap_err();
    let DirectoryPublicationError::Busy { path } = error else {
        panic!("expected contention on the destination's persistent lock: {error}");
    };
    drop(owner);
    assert!(path.is_file());
    path
}

struct Fixture {
    root: TempDir,
    destination: PathBuf,
}

impl Fixture {
    fn new() -> Self {
        let root = TempDir::new().unwrap();
        let destination = root.path().canonicalize().unwrap().join("site");
        Self { root, destination }
    }

    fn write_site(path: &Path, index: &[u8], snapshot: &[u8]) {
        let snapshot_path = path.join(SNAPSHOT);
        fs::create_dir_all(snapshot_path.parent().unwrap()).unwrap();
        fs::write(
            path.join(GENERATED_OUTPUT_MARKER),
            GENERATED_OUTPUT_SIGNATURE,
        )
        .unwrap();
        fs::write(path.join("index.html"), index).unwrap();
        fs::write(snapshot_path, snapshot).unwrap();
    }

    fn journal(&self, destination: &Path) -> PathBuf {
        // Scope depends on the destination filename. Obtain its real journal name in an isolated
        // parent so constructing a second fixture cannot recover the first fixture's backup.
        let template = TempDir::new().unwrap();
        let publication =
            DirectoryPublication::new(&template.path().join(destination.file_name().unwrap()))
                .unwrap();
        let name = publication
            .staging_path()
            .parent()
            .unwrap()
            .file_name()
            .unwrap();
        let journal = self.root.path().join(name);
        fs::create_dir(&journal).unwrap();
        fs::write(
            journal.join("destination"),
            destination.as_os_str().as_encoded_bytes(),
        )
        .unwrap();
        journal
    }
}

#[test]
fn interrupted_publication_recovers_snapshots_before_candidate_copy() {
    let fixture = Fixture::new();
    Fixture::write_site(
        &fixture.destination,
        b"previous index",
        b"immutable snapshot",
    );
    let previous = directory_files(&fixture.destination).unwrap();
    let journal = fixture.journal(&fixture.destination);
    fs::rename(&fixture.destination, journal.join("backup")).unwrap();
    Fixture::write_site(&journal.join("candidate"), b"interrupted", b"incomplete");

    let publication =
        GeneratedTreePublication::new(fixture.root.path(), &fixture.destination).unwrap();

    assert_eq!(
        directory_files(publication.destination()).unwrap(),
        previous
    );
    assert!(!journal.exists());
    let candidate = fixture.root.path().join("next");
    copy_tree(publication.destination(), &candidate).unwrap();
    fs::write(candidate.join("index.html"), b"next index").unwrap();
    let error =
        GeneratedTreePublication::new(fixture.root.path(), &fixture.destination).unwrap_err();
    assert!(matches!(
        error.downcast_ref::<DirectoryPublicationError>(),
        Some(DirectoryPublicationError::Busy { .. })
    ));
    publication.publish(&candidate).unwrap();

    assert_eq!(
        fs::read(fixture.destination.join("index.html")).unwrap(),
        b"next index"
    );
    assert_eq!(
        fs::read(fixture.destination.join(SNAPSHOT)).unwrap(),
        b"immutable snapshot"
    );
    // Publication releases exclusion even though the persistent lock file remains.
    let next = GeneratedTreePublication::new(fixture.root.path(), &fixture.destination).unwrap();
    drop(next);
}

#[test]
fn live_publisher_refuses_another_builder_without_touching_the_site() {
    let fixture = Fixture::new();
    Fixture::write_site(
        &fixture.destination,
        b"previous index",
        b"immutable snapshot",
    );
    let before = directory_files(&fixture.destination).unwrap();
    let owner = DirectoryPublication::new(&fixture.destination).unwrap();
    let journal = owner.staging_path().parent().unwrap().to_owned();

    let error =
        GeneratedTreePublication::new(fixture.root.path(), &fixture.destination).unwrap_err();

    assert!(matches!(
        error.downcast_ref::<DirectoryPublicationError>(),
        Some(DirectoryPublicationError::Busy { .. })
    ));
    assert_eq!(directory_files(&fixture.destination).unwrap(), before);
    assert!(journal.is_dir());
    drop(owner);
    assert!(!journal.exists());
    let retry = GeneratedTreePublication::new(fixture.root.path(), &fixture.destination).unwrap();
    drop(retry);
    assert_eq!(directory_files(&fixture.destination).unwrap(), before);
}

#[test]
fn ambiguous_previous_trees_refuse_publication_and_retain_every_snapshot() {
    for destination_exists in [false, true] {
        let fixture = Fixture::new();
        if destination_exists {
            Fixture::write_site(&fixture.destination, b"current index", b"current snapshot");
        }
        // Initialize the destination's real persistent lock before creating stale journals.
        let lock = initialize_publication_lock(&fixture.destination);
        let first = fixture.journal(&fixture.destination);
        Fixture::write_site(&first.join("backup"), b"first index", b"first snapshot");
        if !destination_exists {
            let second = fixture.journal(&fixture.destination);
            Fixture::write_site(&second.join("backup"), b"second index", b"second snapshot");
        }
        let before = directory_files(fixture.root.path()).unwrap();

        let error =
            GeneratedTreePublication::new(fixture.root.path(), &fixture.destination).unwrap_err();

        let Some(DirectoryPublicationError::RecoveryRequired { destination, .. }) =
            error.downcast_ref::<DirectoryPublicationError>()
        else {
            panic!("expected ambiguous recovery refusal: {error:#}");
        };
        assert_eq!(destination, &fixture.destination);
        assert!(lock.is_file());
        assert_eq!(directory_files(fixture.root.path()).unwrap(), before);
    }
}

#[test]
fn unrelated_destination_journal_is_preserved_during_recovery() {
    let fixture = Fixture::new();
    let other = fixture
        .root
        .path()
        .canonicalize()
        .unwrap()
        .join("other-site");
    let foreign = fixture.journal(&other);
    Fixture::write_site(
        &foreign.join("backup"),
        b"foreign index",
        b"foreign snapshot",
    );
    let before = directory_files(&foreign).unwrap();
    let matching = fixture.journal(&fixture.destination);
    Fixture::write_site(&matching.join("backup"), b"old index", b"retained snapshot");

    let publication =
        GeneratedTreePublication::new(fixture.root.path(), &fixture.destination).unwrap();

    assert_eq!(
        fs::read(publication.destination().join(SNAPSHOT)).unwrap(),
        b"retained snapshot"
    );
    assert_eq!(directory_files(&foreign).unwrap(), before);
    assert!(!other.exists());
    assert!(!matching.exists());
    drop(publication);
    assert_eq!(directory_files(&foreign).unwrap(), before);
}

#[test]
fn recovered_destination_still_requires_documentation_ownership() {
    let fixture = Fixture::new();
    let journal = fixture.journal(&fixture.destination);
    let backup = journal.join("backup");
    Fixture::write_site(&backup, b"unowned index", b"keep snapshot");
    fs::remove_file(backup.join(GENERATED_OUTPUT_MARKER)).unwrap();
    let before = directory_files(&backup).unwrap();

    let error =
        GeneratedTreePublication::new(fixture.root.path(), &fixture.destination).unwrap_err();

    assert!(format!("{error:#}").contains("without generated ownership"));
    assert_eq!(directory_files(&fixture.destination).unwrap(), before);
}

#[test]
fn copied_publication_metadata_refuses_replacement_and_preserves_snapshots() {
    for relative in [
        ".publication-lock-old",
        "products/linnet/.publication-old/backup",
    ] {
        let fixture = Fixture::new();
        Fixture::write_site(&fixture.destination, b"old index", b"immutable snapshot");
        let stale = fixture.destination.join(relative);
        fs::create_dir_all(stale.parent().unwrap()).unwrap();
        fs::write(&stale, b"stale coordination data").unwrap();
        let before = directory_files(&fixture.destination).unwrap();
        let source = fixture.root.path().join("candidate");
        copy_tree(&fixture.destination, &source).unwrap();
        fs::write(source.join("index.html"), b"new index").unwrap();

        let error = GeneratedTreePublication::new(fixture.root.path(), &fixture.destination)
            .unwrap()
            .publish(&source)
            .unwrap_err();

        assert!(format!("{error:#}").contains("internal publication metadata"));
        assert_eq!(directory_files(&fixture.destination).unwrap(), before);
        assert_eq!(
            fs::read(source.join(SNAPSHOT)).unwrap(),
            b"immutable snapshot"
        );
    }
}

#[test]
fn published_site_contains_no_internal_publication_entries() {
    let fixture = Fixture::new();
    let source = fixture.root.path().join("work/site");
    Fixture::write_site(&source, b"site index", b"immutable snapshot");
    let builder = SiteBuilder::discover().unwrap();
    let mut renderer = CliTypstRenderer::new(&builder.root);
    builder
        .write_developer_docs(
            &source,
            BuildChannel::Latest,
            None,
            false,
            None,
            &mut renderer,
        )
        .unwrap();

    clear_stale_staging(&source).unwrap();
    GeneratedTreePublication::new(fixture.root.path(), &fixture.destination)
        .unwrap()
        .publish(&source)
        .unwrap();

    assert!(fixture.destination.join("developers/index.html").is_file());
    assert_eq!(
        fs::read(fixture.destination.join(SNAPSHOT)).unwrap(),
        b"immutable snapshot"
    );
    for entry in fs::read_dir(&fixture.destination).unwrap() {
        let entry = entry.unwrap();
        assert!(
            !entry
                .file_name()
                .to_string_lossy()
                .starts_with(".publication"),
            "internal publication metadata leaked into the site: {}",
            entry.path().display()
        );
    }
}
