use std::collections::BTreeSet;
use std::ffi::OsString;
use std::fs;
#[cfg(unix)]
use std::os::unix::ffi::OsStringExt;
use std::path::{Component, Path, PathBuf};
use std::process::{Command, Output};

use eyre::{Context, ContextCompat, Result, bail, ensure};
use serde::Deserialize;
use serde_json::Value;
use sha2::{Digest, Sha256};

use crate::normalize_repository_path;

const SCHEMA_VERSION: u32 = 1;
const EVIDENCE_ROOT: &str = "docs/architecture/evidence/";

#[derive(Clone, Debug, Deserialize)]
pub(crate) struct GeneratedEvidenceRegistry {
    schema: u32,
    artifact: Vec<GeneratedEvidence>,
}

#[derive(Clone, Debug, Deserialize)]
struct GeneratedEvidence {
    path: PathBuf,
    sha256: String,
    format: String,
    role: String,
    source: PathBuf,
    source_sha256: String,
    artifact_owner: PathBuf,
    artifact_owner_sha256: String,
    artifact_binding: ArtifactBinding,
    source_owner: PathBuf,
    source_owner_sha256: String,
    source_binding: SourceBinding,
    producer_binding: ProducerBinding,
    producer_revision: String,
    tool_owner: PathBuf,
    tool_owner_sha256: String,
    tool: String,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
enum ArtifactBinding {
    ReceiptEvidenceSha256 { key: String },
    InventoryEntry { index: usize },
}

#[derive(Clone, Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
enum SourceBinding {
    ImplementationFile { index: usize },
    ChangedCodeFile { index: usize },
    RenderBinding,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[serde(rename_all = "kebab-case")]
enum ProducerBinding {
    SourceMappingCodeCommit,
    CodeCommit,
    Code,
    PreservedCandidate,
}

impl GeneratedEvidenceRegistry {
    #[cfg(test)]
    pub(crate) fn empty_for_tests() -> Self {
        Self {
            schema: SCHEMA_VERSION,
            artifact: Vec::new(),
        }
    }

    pub(crate) fn validate(
        &self,
        root: &Path,
        canonical_sources: &BTreeSet<PathBuf>,
    ) -> Result<BTreeSet<PathBuf>> {
        ensure!(
            self.schema == SCHEMA_VERSION,
            "generated rendered evidence schema {} does not match checker schema {SCHEMA_VERSION}",
            self.schema
        );
        let mut paths = BTreeSet::new();
        for artifact in &self.artifact {
            artifact.validate(root, canonical_sources)?;
            ensure!(
                paths.insert(artifact.path.clone()),
                "duplicate generated rendered evidence path {}",
                artifact.path.display()
            );
        }
        Ok(paths)
    }
}

impl GeneratedEvidence {
    fn validate(&self, root: &Path, canonical_sources: &BTreeSet<PathBuf>) -> Result<()> {
        for (kind, path) in [
            ("artifact", &self.path),
            ("canonical source", &self.source),
            ("artifact owner", &self.artifact_owner),
            ("source owner", &self.source_owner),
            ("tool owner", &self.tool_owner),
        ] {
            ensure_normalized(path, kind)?;
            ensure_regular_file_without_symlinks(root, path, kind)?;
        }
        ensure!(
            self.path.starts_with(EVIDENCE_ROOT),
            "generated rendered evidence {} must be under {EVIDENCE_ROOT}",
            self.path.display()
        );
        ensure!(
            self.path.extension().and_then(|value| value.to_str()) == Some("html")
                && self.format == "wrapped-html",
            "generated rendered evidence {} has unsupported format {}",
            self.path.display(),
            self.format
        );
        ensure!(
            canonical_sources.contains(&self.source),
            "generated rendered evidence {} binds unregistered developer source {}",
            self.path.display(),
            self.source.display()
        );
        ensure!(
            self.role == "historical-render-proof",
            "generated rendered evidence {} has unsupported role {}",
            self.path.display(),
            self.role
        );
        ensure_hash(root, &self.path, &self.sha256, "artifact")?;
        ensure_sha256(&self.source_sha256, "historical source", &self.source)?;
        ensure_revision(&self.producer_revision, &self.path)?;
        ensure_hash(
            root,
            &self.artifact_owner,
            &self.artifact_owner_sha256,
            "artifact owner",
        )?;
        ensure_hash(
            root,
            &self.source_owner,
            &self.source_owner_sha256,
            "source owner",
        )?;
        ensure_hash(
            root,
            &self.tool_owner,
            &self.tool_owner_sha256,
            "tool owner",
        )?;

        let artifact_owner: Value = serde_json::from_slice(&fs::read(
            root.join(&self.artifact_owner),
        )?)
        .wrap_err_with(|| {
            format!(
                "invalid generated evidence owner {}",
                self.artifact_owner.display()
            )
        })?;
        let source_owner: Value = serde_json::from_slice(&fs::read(root.join(&self.source_owner))?)
            .wrap_err_with(|| {
                format!(
                    "invalid generated evidence owner {}",
                    self.source_owner.display()
                )
            })?;
        self.artifact_binding.validate(
            &artifact_owner,
            &self.artifact_owner,
            &self.path,
            &self.sha256,
        )?;
        self.source_binding.validate(
            &source_owner,
            &self.source_owner,
            &self.source,
            &self.source_sha256,
        )?;
        self.producer_binding.validate(
            &source_owner,
            &self.source_owner,
            &self.producer_revision,
        )?;
        let tool = fs::read_to_string(root.join(&self.tool_owner))?;
        ensure!(
            !self.tool.trim().is_empty() && tool.trim() == self.tool,
            "generated evidence tool owner {} does not bind declared tool {}",
            self.tool_owner.display(),
            self.tool
        );
        Ok(())
    }
}

fn ensure_normalized(path: &Path, kind: &str) -> Result<()> {
    let normalized = normalize_repository_path(path);
    ensure!(
        !path.as_os_str().is_empty()
            && !path.is_absolute()
            && normalized
                .as_ref()
                .is_some_and(|normalized| normalized.as_os_str() == path.as_os_str())
            && path
                .components()
                .all(|component| matches!(component, Component::Normal(_))),
        "generated evidence {kind} path {} must be normalized and repository-relative",
        path.display()
    );
    Ok(())
}

fn ensure_regular_file_without_symlinks(root: &Path, path: &Path, kind: &str) -> Result<()> {
    let mut current = root.to_path_buf();
    for component in path.components() {
        let Component::Normal(part) = component else {
            bail!(
                "generated evidence {kind} path {} is not normalized",
                path.display()
            );
        };
        current.push(part);
        let metadata = fs::symlink_metadata(&current)
            .wrap_err_with(|| format!("missing generated evidence {kind} {}", path.display()))?;
        ensure!(
            !metadata.file_type().is_symlink(),
            "generated evidence {kind} {} contains symlink component {}",
            path.display(),
            current.display()
        );
    }
    ensure!(
        current.is_file(),
        "generated evidence {kind} {} must be a regular file",
        path.display()
    );
    Ok(())
}

fn ensure_hash(root: &Path, path: &Path, expected: &str, kind: &str) -> Result<()> {
    ensure_sha256(expected, kind, path)?;
    let actual = format!("{:x}", Sha256::digest(fs::read(root.join(path))?));
    ensure!(
        actual == expected,
        "generated evidence {kind} {} changed (expected {expected}, found {actual})",
        path.display()
    );
    Ok(())
}

fn ensure_sha256(value: &str, kind: &str, path: &Path) -> Result<()> {
    ensure!(
        value.len() == 64
            && value
                .bytes()
                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte)),
        "generated evidence {kind} {} has an invalid lowercase SHA-256",
        path.display()
    );
    Ok(())
}

fn ensure_revision(value: &str, path: &Path) -> Result<()> {
    ensure!(
        value.len() == 40
            && value
                .bytes()
                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte)),
        "generated rendered evidence {} has an invalid lowercase producer revision",
        path.display()
    );
    Ok(())
}

impl ArtifactBinding {
    fn validate(
        &self,
        owner: &Value,
        owner_path: &Path,
        artifact_path: &Path,
        expected_hash: &str,
    ) -> Result<()> {
        let owner_parent = owner_path
            .parent()
            .context("generated evidence owner has no parent")?;
        let (relative, actual_hash) = match self {
            Self::ReceiptEvidenceSha256 { key } => {
                ensure!(
                    normalize_repository_path(Path::new(key)).is_some_and(|normalized| normalized
                        .as_os_str()
                        == Path::new(key).as_os_str()),
                    "generated evidence receipt key {key} is not normalized"
                );
                let hash = owner
                    .get("evidence_sha256")
                    .and_then(|value| value.get(key))
                    .and_then(Value::as_str)
                    .context("generated evidence receipt hash binding is missing")?;
                (PathBuf::from(key), hash)
            }
            Self::InventoryEntry { index } => {
                let entry = owner
                    .get(*index)
                    .context("generated evidence inventory entry is missing")?;
                let path = entry
                    .get("path")
                    .and_then(Value::as_str)
                    .context("generated evidence inventory path binding is missing")?;
                ensure!(
                    normalize_repository_path(Path::new(path)).is_some_and(|normalized| normalized
                        .as_os_str()
                        == Path::new(path).as_os_str()),
                    "generated evidence inventory path {path} is not normalized"
                );
                let hash = entry
                    .get("sha256")
                    .and_then(Value::as_str)
                    .context("generated evidence inventory hash binding is missing")?;
                (PathBuf::from(path), hash)
            }
        };
        let bound = normalize_repository_path(&owner_parent.join(relative))
            .context("generated evidence artifact binding escapes the repository")?;
        ensure!(
            bound == artifact_path && actual_hash == expected_hash,
            "generated evidence owner {} has a stale artifact path or hash binding",
            owner_path.display()
        );
        Ok(())
    }
}

impl SourceBinding {
    fn validate(
        &self,
        owner: &Value,
        owner_path: &Path,
        source_path: &Path,
        expected_hash: &str,
    ) -> Result<()> {
        let binding = match self {
            Self::ImplementationFile { index } => owner
                .get("implementation_files")
                .and_then(|value| value.get(*index)),
            Self::ChangedCodeFile { index } => owner
                .get("changed_code_files")
                .and_then(|value| value.get(*index)),
            Self::RenderBinding => owner.get("render_binding"),
        }
        .context("generated evidence historical source binding is missing")?;
        let (path, hash) = match self {
            Self::ImplementationFile { .. } | Self::ChangedCodeFile { .. } => (
                binding.get("path").and_then(Value::as_str),
                binding.get("sha256").and_then(Value::as_str),
            ),
            Self::RenderBinding => (
                binding.get("source").and_then(Value::as_str),
                binding.get("source_sha256").and_then(Value::as_str),
            ),
        };
        ensure!(
            path == source_path.to_str() && hash == Some(expected_hash),
            "generated evidence owner {} has a stale historical source path or hash binding",
            owner_path.display()
        );
        Ok(())
    }
}

impl ProducerBinding {
    fn validate(self, owner: &Value, owner_path: &Path, expected: &str) -> Result<()> {
        let actual = match self {
            Self::SourceMappingCodeCommit => owner
                .get("source_mapping")
                .and_then(|value| value.get("code_commit")),
            Self::CodeCommit => owner.get("code_commit"),
            Self::Code => owner.get("code"),
            Self::PreservedCandidate => owner.get("preserved_candidate"),
        }
        .and_then(Value::as_str);
        ensure!(
            actual == Some(expected),
            "generated evidence owner {} has a stale producer revision binding",
            owner_path.display()
        );
        Ok(())
    }
}

pub(crate) enum RepositorySources {
    MetadataFree,
    Listed(BTreeSet<PathBuf>),
}

impl RepositorySources {
    pub(crate) fn discover(root: &Path) -> Result<Self> {
        if root.join(".jj").exists() {
            return Self::from_jj(root);
        }
        if root.join(".git").exists() {
            return Self::from_git(root);
        }
        Ok(Self::MetadataFree)
    }

    fn from_git(root: &Path) -> Result<Self> {
        let git_dir = selected_git_dir(root)?;
        ensure_no_index_gitlinks(root, &git_dir)?;
        let output = controlled_git_command(root)
            .arg("--git-dir")
            .arg(&git_dir)
            .arg("--work-tree")
            .arg(root)
            .args([
                "ls-files",
                "-z",
                "--cached",
                "--others",
                "--exclude-standard",
            ])
            .output()
            .context("failed to enumerate Git repository sources")?;
        ensure_command_success(&output, "Git repository source enumeration")?;
        ensure_empty_stderr(&output, "Git repository source enumeration")?;
        Ok(Self::Listed(parse_complete_file_paths(
            &output.stdout,
            "Git repository source enumeration",
        )?))
    }

    fn from_jj(root: &Path) -> Result<Self> {
        let tracked = controlled_jj_command(root)
            .args([
                "--ignore-working-copy",
                "--quiet",
                "file",
                "list",
                "-r",
                "@",
                "-T",
                "path ++ \"\\0\" ++ file_type ++ \"\\0\"",
            ])
            .output()
            .context("failed to enumerate the JJ tracked tree")?;
        ensure_command_success(&tracked, "JJ tracked-tree enumeration")?;
        ensure_empty_stderr(&tracked, "JJ tracked-tree enumeration")?;
        let mut listed = parse_jj_tracked_paths(&tracked.stdout)?;

        let git_root = controlled_jj_command(root)
            .args(["--ignore-working-copy", "--quiet", "git", "root"])
            .output()
            .context("failed to locate the JJ Git backend")?;
        ensure_command_success(&git_root, "JJ Git-backend discovery")?;
        ensure_empty_stderr(&git_root, "JJ Git-backend discovery")?;
        let git_dir = PathBuf::from(
            String::from_utf8(git_root.stdout)
                .context("JJ Git backend path is not UTF-8")?
                .trim(),
        );
        ensure!(git_dir.is_dir(), "JJ Git backend does not exist");

        let empty_index_directory = tempfile::tempdir()?;
        let empty_index = empty_index_directory.path().join("index");
        let working_files = controlled_git_command(root)
            .arg("--git-dir")
            .arg(&git_dir)
            .arg("--work-tree")
            .arg(root)
            .args(["ls-files", "-z", "--others", "--exclude-standard"])
            .env("GIT_INDEX_FILE", empty_index)
            .output()
            .context("failed to enumerate nonignored JJ working-tree files")?;
        ensure_command_success(&working_files, "JJ working-tree source enumeration")?;
        ensure_empty_stderr(&working_files, "JJ working-tree source enumeration")?;
        listed.extend(parse_complete_file_paths(
            &working_files.stdout,
            "JJ working-tree source enumeration",
        )?);
        Ok(Self::Listed(listed))
    }

    pub(crate) fn contains(&self, path: &Path) -> bool {
        match self {
            Self::MetadataFree => true,
            Self::Listed(paths) => paths.contains(path),
        }
    }

    pub(crate) fn may_contain(&self, directory: &Path) -> bool {
        match self {
            Self::MetadataFree => true,
            Self::Listed(paths) => paths.iter().any(|path| path.starts_with(directory)),
        }
    }
}

const GIT_SELECTOR_ENVIRONMENT: [&str; 12] = [
    "GIT_DIR",
    "GIT_WORK_TREE",
    "GIT_INDEX_FILE",
    "GIT_COMMON_DIR",
    "GIT_OBJECT_DIRECTORY",
    "GIT_ALTERNATE_OBJECT_DIRECTORIES",
    "GIT_QUARANTINE_PATH",
    "GIT_PREFIX",
    "GIT_CEILING_DIRECTORIES",
    "GIT_DISCOVERY_ACROSS_FILESYSTEM",
    "GIT_CONFIG_PARAMETERS",
    "GIT_CONFIG_COUNT",
];

fn controlled_git_command(root: &Path) -> Command {
    controlled_repository_command("git", root)
}

fn controlled_jj_command(root: &Path) -> Command {
    controlled_repository_command("jj", root)
}

fn controlled_repository_command(program: &str, root: &Path) -> Command {
    let mut command = Command::new(program);
    command.current_dir(root);
    for variable in GIT_SELECTOR_ENVIRONMENT {
        command.env_remove(variable);
    }
    command
}

fn selected_git_dir(root: &Path) -> Result<PathBuf> {
    let output = controlled_git_command(root)
        .args(["rev-parse", "--absolute-git-dir"])
        .output()
        .context("failed to resolve selected Git checkout metadata")?;
    ensure_command_success(&output, "Git checkout metadata discovery")?;
    ensure_empty_stderr(&output, "Git checkout metadata discovery")?;
    let git_dir = PathBuf::from(
        String::from_utf8(output.stdout)
            .context("selected Git checkout metadata path is not UTF-8")?
            .trim(),
    );
    ensure!(
        git_dir.is_dir(),
        "selected Git checkout metadata does not exist"
    );
    Ok(git_dir)
}

fn ensure_no_index_gitlinks(root: &Path, git_dir: &Path) -> Result<()> {
    let output = controlled_git_command(root)
        .arg("--git-dir")
        .arg(git_dir)
        .arg("--work-tree")
        .arg(root)
        .args(["ls-files", "--stage", "-z"])
        .output()
        .context("failed to inspect Git index entry kinds")?;
    ensure_command_success(&output, "Git index entry-kind inspection")?;
    ensure_empty_stderr(&output, "Git index entry-kind inspection")?;
    ensure_no_gitlinks(&output.stdout, "Git index")
}

fn ensure_no_gitlinks(bytes: &[u8], inventory: &str) -> Result<()> {
    for entry in bytes
        .split(|byte| *byte == 0)
        .filter(|entry| !entry.is_empty())
    {
        if entry.starts_with(b"160000 ") {
            let path = entry
                .iter()
                .position(|byte| *byte == b'\t')
                .map_or(entry, |separator| &entry[separator + 1..]);
            bail!(
                "{inventory} contains gitlink {}; nested tracked repositories cannot be treated as complete source inventory",
                String::from_utf8_lossy(path)
            );
        }
    }
    Ok(())
}

fn ensure_command_success(output: &Output, operation: &str) -> Result<()> {
    ensure!(
        output.status.success(),
        "{operation} failed: {}",
        String::from_utf8_lossy(&output.stderr).trim()
    );
    Ok(())
}

fn ensure_empty_stderr(output: &Output, operation: &str) -> Result<()> {
    ensure!(
        output.stderr.is_empty(),
        "{operation} produced a warning and cannot be trusted: {}",
        String::from_utf8_lossy(&output.stderr).trim()
    );
    Ok(())
}

fn parse_complete_file_paths(bytes: &[u8], operation: &str) -> Result<BTreeSet<PathBuf>> {
    for path in bytes
        .split(|byte| *byte == 0)
        .filter(|path| !path.is_empty())
    {
        ensure!(
            !path.ends_with(b"/"),
            "{operation} returned opaque directory {}; nested repositories cannot be treated as complete source inventory",
            String::from_utf8_lossy(path)
        );
    }
    parse_nul_paths(bytes)
}

fn parse_jj_tracked_paths(bytes: &[u8]) -> Result<BTreeSet<PathBuf>> {
    let fields = bytes
        .split(|byte| *byte == 0)
        .filter(|field| !field.is_empty())
        .collect::<Vec<_>>();
    let (pairs, remainder) = fields.as_chunks::<2>();
    ensure!(
        remainder.is_empty(),
        "JJ tracked-tree enumeration returned malformed path/type framing"
    );
    let mut paths = BTreeSet::new();
    for pair in pairs {
        ensure!(
            matches!(pair[1], b"file" | b"symlink"),
            "JJ tracked tree contains unsupported {} entry {}; nested tracked repositories cannot be treated as complete source inventory",
            String::from_utf8_lossy(pair[1]),
            String::from_utf8_lossy(pair[0])
        );
        paths.extend(parse_nul_paths(&[pair[0], b"\0"].concat())?);
    }
    Ok(paths)
}

#[cfg(unix)]
fn parse_nul_paths(bytes: &[u8]) -> Result<BTreeSet<PathBuf>> {
    Ok(bytes
        .split(|byte| *byte == 0)
        .filter(|path| !path.is_empty())
        .map(|path| PathBuf::from(OsString::from_vec(path.to_vec())))
        .collect())
}

#[cfg(not(unix))]
fn parse_nul_paths(bytes: &[u8]) -> Result<BTreeSet<PathBuf>> {
    Ok(String::from_utf8(bytes.to_vec())?
        .split('\0')
        .filter(|path| !path.is_empty())
        .map(PathBuf::from)
        .collect())
}

#[cfg(test)]
mod tests {
    use std::env;
    use std::fs;
    #[cfg(unix)]
    use std::os::unix::ffi::OsStringExt;
    #[cfg(unix)]
    use std::os::unix::fs::{PermissionsExt, symlink};
    use std::process::Command;

    use serde_json::json;
    use tempfile::TempDir;

    use super::*;

    fn write(root: &Path, path: &str, bytes: &[u8]) {
        let path = root.join(path);
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(path, bytes).unwrap();
    }

    fn hash(root: &Path, path: &str) -> String {
        format!("{:x}", Sha256::digest(fs::read(root.join(path)).unwrap()))
    }

    fn fixture() -> (TempDir, GeneratedEvidenceRegistry) {
        let temporary = tempfile::tempdir().unwrap();
        let root = temporary.path();
        let artifact = "docs/architecture/evidence/review/note.html";
        let source = "docs/architecture/note.typ";
        let owner = "docs/architecture/evidence/review/receipt.json";
        let tool_owner = "docs/architecture/evidence/review/typst-version.log";
        write(root, artifact, b"<html>rendered</html>");
        write(root, source, b"= Current note");
        write(root, tool_owner, b"typst 0.15.1\n");
        let artifact_sha256 = hash(root, artifact);
        let source_sha256 = "1234567890abcdef1234567890abcdef1234567890abcdef1234567890abcdef";
        let producer_revision = "1234567890abcdef1234567890abcdef12345678";
        write(
            root,
            owner,
            serde_json::to_vec(&json!({
                "evidence_sha256": {"note.html": artifact_sha256},
                "implementation_files": [{"path": source, "sha256": source_sha256}],
                "source_mapping": {"code_commit": producer_revision},
            }))
            .unwrap()
            .as_slice(),
        );
        let entry = GeneratedEvidence {
            path: artifact.into(),
            sha256: artifact_sha256,
            format: "wrapped-html".into(),
            role: "historical-render-proof".into(),
            source: source.into(),
            source_sha256: source_sha256.into(),
            artifact_owner: owner.into(),
            artifact_owner_sha256: hash(root, owner),
            artifact_binding: ArtifactBinding::ReceiptEvidenceSha256 {
                key: "note.html".into(),
            },
            source_owner: owner.into(),
            source_owner_sha256: hash(root, owner),
            source_binding: SourceBinding::ImplementationFile { index: 0 },
            producer_binding: ProducerBinding::SourceMappingCodeCommit,
            producer_revision: producer_revision.into(),
            tool_owner: tool_owner.into(),
            tool_owner_sha256: hash(root, tool_owner),
            tool: "typst 0.15.1".into(),
        };
        (
            temporary,
            GeneratedEvidenceRegistry {
                schema: SCHEMA_VERSION,
                artifact: vec![entry],
            },
        )
    }

    #[test]
    fn registry_validates_complete_historical_bindings_and_all_entries() {
        let (temporary, registry) = fixture();
        let canonical_sources = BTreeSet::from([registry.artifact[0].source.clone()]);
        let paths = registry
            .validate(temporary.path(), &canonical_sources)
            .unwrap();
        assert_eq!(paths, BTreeSet::from([registry.artifact[0].path.clone()]));

        let mut duplicate = registry.clone();
        duplicate.artifact.push(duplicate.artifact[0].clone());
        assert!(
            duplicate
                .validate(temporary.path(), &canonical_sources)
                .is_err()
        );

        for mutate in ["artifact", "source", "producer", "role", "location"] {
            let mut invalid = registry.clone();
            match mutate {
                "artifact" => invalid.artifact[0].sha256.replace_range(0..1, "g"),
                "source" => invalid.artifact[0].source_sha256.replace_range(0..1, "g"),
                "producer" => invalid.artifact[0]
                    .producer_revision
                    .replace_range(0..1, "g"),
                "role" => invalid.artifact[0].role = "authored-prose".into(),
                "location" => invalid.artifact[0].path = "docs/rendered.html".into(),
                _ => unreachable!(),
            }
            assert!(
                invalid
                    .validate(temporary.path(), &canonical_sources)
                    .is_err(),
                "{mutate}"
            );
        }

        let mut traversal = registry.clone();
        traversal.artifact[0].artifact_owner = PathBuf::from("docs/../receipt.json");
        assert!(
            traversal
                .validate(temporary.path(), &canonical_sources)
                .is_err()
        );
        fs::remove_file(temporary.path().join(&registry.artifact[0].artifact_owner)).unwrap();
        assert!(
            registry
                .validate(temporary.path(), &canonical_sources)
                .is_err()
        );
    }

    #[test]
    fn registry_rejects_synchronized_malformed_and_unregistered_bindings() {
        let (temporary, mut registry) = fixture();
        let root = temporary.path();
        let canonical_sources = BTreeSet::from([registry.artifact[0].source.clone()]);

        registry.artifact[0].source_sha256 = "x".repeat(64);
        let owner = root.join(&registry.artifact[0].source_owner);
        let mut receipt: Value = serde_json::from_slice(&fs::read(&owner).unwrap()).unwrap();
        receipt["implementation_files"][0]["sha256"] = Value::String("x".repeat(64));
        fs::write(&owner, serde_json::to_vec(&receipt).unwrap()).unwrap();
        registry.artifact[0].source_owner_sha256 =
            hash(root, registry.artifact[0].source_owner.to_str().unwrap());
        registry.artifact[0].artifact_owner_sha256 =
            registry.artifact[0].source_owner_sha256.clone();
        assert!(registry.validate(root, &canonical_sources).is_err());

        let (temporary, mut registry) = fixture();
        let root = temporary.path();
        registry.artifact[0].producer_revision = "2026-10-28".into();
        let owner = root.join(&registry.artifact[0].source_owner);
        let mut receipt: Value = serde_json::from_slice(&fs::read(&owner).unwrap()).unwrap();
        receipt["source_mapping"]["code_commit"] = Value::String("2026-10-28".into());
        fs::write(&owner, serde_json::to_vec(&receipt).unwrap()).unwrap();
        registry.artifact[0].source_owner_sha256 =
            hash(root, registry.artifact[0].source_owner.to_str().unwrap());
        registry.artifact[0].artifact_owner_sha256 =
            registry.artifact[0].source_owner_sha256.clone();
        let canonical_sources = BTreeSet::from([registry.artifact[0].source.clone()]);
        assert!(registry.validate(root, &canonical_sources).is_err());

        let (temporary, registry) = fixture();
        assert!(
            registry
                .validate(temporary.path(), &BTreeSet::new())
                .is_err()
        );
    }

    #[test]
    fn typed_source_bindings_ignore_cross_role_aliases() {
        let temporary = tempfile::tempdir().unwrap();
        let root = temporary.path();
        let artifact = "docs/architecture/evidence/review/note.html";
        let source = "docs/architecture/note.typ";
        let distractor = "docs/architecture/distractor.typ";
        let owner = "docs/architecture/evidence/review/receipt.json";
        let tool_owner = "docs/architecture/evidence/review/typst-version.log";
        write(root, artifact, b"<html>rendered</html>");
        write(root, source, b"= Historical source");
        write(root, distractor, b"= Distractor");
        write(root, tool_owner, b"typst 0.15.1\n");
        let artifact_sha256 = hash(root, artifact);
        let source_sha256 = hash(root, source);
        let distractor_sha256 = hash(root, distractor);
        let producer_revision = "1234567890abcdef1234567890abcdef12345678";
        write(
            root,
            owner,
            serde_json::to_vec(&json!({
                "evidence_sha256": {"note.html": artifact_sha256},
                "render_binding": {
                    "source": source,
                    "source_sha256": source_sha256,
                    "path": distractor,
                    "sha256": distractor_sha256,
                },
                "preserved_candidate": producer_revision,
            }))
            .unwrap()
            .as_slice(),
        );
        let entry = GeneratedEvidence {
            path: artifact.into(),
            sha256: artifact_sha256,
            format: "wrapped-html".into(),
            role: "historical-render-proof".into(),
            source: source.into(),
            source_sha256,
            artifact_owner: owner.into(),
            artifact_owner_sha256: hash(root, owner),
            artifact_binding: ArtifactBinding::ReceiptEvidenceSha256 {
                key: "note.html".into(),
            },
            source_owner: owner.into(),
            source_owner_sha256: hash(root, owner),
            source_binding: SourceBinding::RenderBinding,
            producer_binding: ProducerBinding::PreservedCandidate,
            producer_revision: producer_revision.into(),
            tool_owner: tool_owner.into(),
            tool_owner_sha256: hash(root, tool_owner),
            tool: "typst 0.15.1".into(),
        };
        let canonical_sources = BTreeSet::from([source.into(), distractor.into()]);
        let registry = GeneratedEvidenceRegistry {
            schema: SCHEMA_VERSION,
            artifact: vec![entry],
        };
        registry.validate(root, &canonical_sources).unwrap();

        let mut wrong = registry.clone();
        wrong.artifact[0].source = distractor.into();
        wrong.artifact[0].source_sha256 = distractor_sha256;
        assert!(wrong.validate(root, &canonical_sources).is_err());

        let inventory_owner = json!({
            "implementation_files": [{
                "path": source,
                "sha256": hash(root, source),
                "source": distractor,
                "source_sha256": hash(root, distractor),
            }]
        });
        SourceBinding::ImplementationFile { index: 0 }
            .validate(
                &inventory_owner,
                Path::new(owner),
                Path::new(source),
                &hash(root, source),
            )
            .unwrap();
        assert!(
            SourceBinding::ImplementationFile { index: 0 }
                .validate(
                    &inventory_owner,
                    Path::new(owner),
                    Path::new(distractor),
                    &hash(root, distractor),
                )
                .is_err()
        );
    }

    #[cfg(unix)]
    #[test]
    fn registry_rejects_symlinked_ancestor_inside_repository() {
        let (temporary, registry) = fixture();
        let root = temporary.path();
        let artifact = root.join(&registry.artifact[0].path);
        let directory = artifact.parent().unwrap();
        let retained = directory.with_file_name("retained-review");
        fs::rename(directory, &retained).unwrap();
        symlink(&retained, directory).unwrap();
        let canonical_sources = BTreeSet::from([registry.artifact[0].source.clone()]);
        assert!(registry.validate(root, &canonical_sources).is_err());
    }

    #[test]
    fn git_enumeration_keeps_tracked_ignored_and_nonignored_untracked_sources() {
        let temporary = tempfile::tempdir().unwrap();
        let root = temporary.path();
        assert!(
            Command::new("git")
                .arg("init")
                .current_dir(root)
                .status()
                .unwrap()
                .success()
        );
        write(
            root,
            ".gitignore",
            b".local/\nignored/*\n!ignored/keep.html\n",
        );
        for path in [
            ".local/scratch.html",
            "ordinary.html",
            "ignored/drop.html",
            "ignored/keep.html",
            "ignored/tracked.html",
        ] {
            write(root, path, b"candidate");
        }
        assert!(
            Command::new("git")
                .args(["add", ".gitignore", "-f", "ignored/tracked.html"])
                .current_dir(root)
                .status()
                .unwrap()
                .success()
        );
        let sources = RepositorySources::discover(root).unwrap();
        assert!(sources.contains(Path::new("ordinary.html")));
        assert!(sources.contains(Path::new("ignored/keep.html")));
        assert!(sources.contains(Path::new("ignored/tracked.html")));
        assert!(!sources.contains(Path::new("ignored/drop.html")));
        assert!(!sources.contains(Path::new(".local/scratch.html")));

        for path in [
            "café.md",
            "quoted\"name.md",
            "tab\tname.md",
            "back\\slash.md",
            "new\nline.md",
        ] {
            write(root, path, b"candidate");
        }
        #[cfg(unix)]
        let non_utf8 = { PathBuf::from(OsString::from_vec(b"non-utf8-\xff.md".to_vec())) };
        #[cfg(unix)]
        fs::write(root.join(&non_utf8), b"candidate").unwrap();
        let sources = RepositorySources::discover(root).unwrap();
        for path in [
            "café.md",
            "quoted\"name.md",
            "tab\tname.md",
            "back\\slash.md",
            "new\nline.md",
        ] {
            assert!(sources.contains(Path::new(path)), "{path:?}");
        }
        #[cfg(unix)]
        assert!(sources.contains(&non_utf8));
    }

    #[test]
    fn jj_enumeration_is_read_only_and_independent_of_auto_tracking() {
        let temporary = tempfile::tempdir().unwrap();
        let root = temporary.path();
        let init = Command::new("jj")
            .args(["git", "init", "--no-colocate", "."])
            .current_dir(root)
            .output()
            .unwrap();
        assert!(init.status.success(), "{init:?}");
        let configure = Command::new("jj")
            .args(["config", "set", "--repo", "snapshot.auto-track", "none()"])
            .current_dir(root)
            .output()
            .unwrap();
        assert!(configure.status.success(), "{configure:?}");
        write(
            root,
            ".gitignore",
            b".local/\nignored/*\n!ignored/keep.html\n",
        );
        write(root, "ordinary.html", &vec![b'x'; 2_100_000]);
        write(root, "accepted.typ", b"= Accepted");
        write(root, ".local/scratch.html", b"ignored");
        write(root, "ignored/keep.html", b"not ignored");
        let bookmark = Command::new("jj")
            .args(["bookmark", "create", "frozen", "-r", "@"])
            .current_dir(root)
            .output()
            .unwrap();
        assert!(bookmark.status.success(), "{bookmark:?}");
        let before = jj_state(root);

        let sources = RepositorySources::discover(root).unwrap();

        assert!(sources.contains(Path::new("ordinary.html")));
        assert!(sources.contains(Path::new("accepted.typ")));
        assert!(sources.contains(Path::new("ignored/keep.html")));
        assert!(!sources.contains(Path::new(".local/scratch.html")));
        assert_eq!(jj_state(root), before);
    }

    #[test]
    fn repository_enumeration_ignores_inherited_git_selectors() {
        const CASE: &str = "ALPHAL00P_DOCS_GIT_SELECTOR_TEST";
        if env::var_os(CASE).is_none() {
            let temporary = tempfile::tempdir().unwrap();
            let alternate = temporary.path().join("alternate.git");
            let initialized = Command::new("git")
                .args(["init", "--bare"])
                .arg(&alternate)
                .output()
                .unwrap();
            assert!(initialized.status.success(), "{initialized:?}");
            for (name, variables) in [
                (
                    "index",
                    vec![("GIT_INDEX_FILE", temporary.path().join("missing-index"))],
                ),
                ("dir", vec![("GIT_DIR", alternate.clone())]),
                (
                    "combined",
                    vec![
                        ("GIT_DIR", alternate.clone()),
                        ("GIT_WORK_TREE", temporary.path().join("other-worktree")),
                        ("GIT_COMMON_DIR", alternate.clone()),
                        ("GIT_INDEX_FILE", temporary.path().join("missing-index")),
                    ],
                ),
            ] {
                let mut command = Command::new(env::current_exe().unwrap());
                command
                    .args([
                        "--exact",
                        "generated_evidence::tests::repository_enumeration_ignores_inherited_git_selectors",
                        "--nocapture",
                    ])
                    .env(CASE, name);
                for (variable, value) in variables {
                    command.env(variable, value);
                }
                assert!(command.status().unwrap().success(), "{name}");
            }
            return;
        }

        let temporary = tempfile::tempdir().unwrap();
        let root = temporary.path();
        assert!(
            controlled_git_command(root)
                .arg("init")
                .status()
                .unwrap()
                .success()
        );
        write(root, "tracked.html", b"tracked prose");
        assert!(
            controlled_git_command(root)
                .args(["add", "tracked.html"])
                .status()
                .unwrap()
                .success()
        );
        write(root, ".gitignore", b"tracked.html\n");
        let sources = RepositorySources::discover(root).unwrap();
        assert!(sources.contains(Path::new("tracked.html")));
    }

    #[test]
    fn jj_working_files_are_root_relative_from_any_caller_directory() {
        const CASE: &str = "ALPHAL00P_DOCS_JJ_CALLER_DIRECTORY_TEST";
        let Ok(case) = env::var(CASE) else {
            for case in ["subdirectory", "external"] {
                let status = Command::new(env::current_exe().unwrap())
                    .args([
                        "--exact",
                        "generated_evidence::tests::jj_working_files_are_root_relative_from_any_caller_directory",
                        "--nocapture",
                    ])
                    .env(CASE, case)
                    .status()
                    .unwrap();
                assert!(status.success(), "{case}");
            }
            return;
        };

        let temporary = tempfile::tempdir().unwrap();
        let root = temporary.path().join("repository");
        fs::create_dir(&root).unwrap();
        let init = controlled_jj_command(&root)
            .args(["git", "init", "--no-colocate", "."])
            .output()
            .unwrap();
        assert!(init.status.success(), "{init:?}");
        write(&root, "root.html", b"root prose");
        write(&root, "sub/nested.html", b"nested prose");
        let caller = if case == "subdirectory" {
            root.join("sub")
        } else {
            temporary.path().join("external")
        };
        fs::create_dir_all(&caller).unwrap();
        env::set_current_dir(&caller).unwrap();

        let sources = RepositorySources::discover(&root).unwrap();
        assert!(sources.contains(Path::new("root.html")));
        assert!(sources.contains(Path::new("sub/nested.html")));
    }

    #[cfg(unix)]
    #[test]
    fn git_warning_from_unreadable_directory_fails_closed() {
        let temporary = tempfile::tempdir().unwrap();
        let root = temporary.path();
        assert!(
            controlled_git_command(root)
                .arg("init")
                .status()
                .unwrap()
                .success()
        );
        write(root, "hidden/rogue.html", b"prose");
        let directory = root.join("hidden");
        let original = fs::metadata(&directory).unwrap().permissions();
        let mut unreadable = original.clone();
        unreadable.set_mode(0o000);
        fs::set_permissions(&directory, unreadable).unwrap();

        let result = RepositorySources::discover(root);

        fs::set_permissions(&directory, original).unwrap();
        assert!(result.is_err());
        let sources = RepositorySources::discover(root).unwrap();
        assert!(sources.contains(Path::new("hidden/rogue.html")));
    }

    #[test]
    fn nested_repositories_fail_closed_but_ignored_nested_scratch_stays_excluded() {
        for vcs in ["git", "jj"] {
            let temporary = tempfile::tempdir().unwrap();
            let root = temporary.path();
            initialize_vcs(root, vcs);
            write(root, "nested/unauthorized.html", b"authored prose");
            assert!(
                RepositorySources::discover(root)
                    .unwrap()
                    .contains(Path::new("nested/unauthorized.html")),
                "{vcs} pre-marker control"
            );
            assert!(
                controlled_git_command(&root.join("nested"))
                    .arg("init")
                    .status()
                    .unwrap()
                    .success()
            );
            assert!(RepositorySources::discover(root).is_err(), "{vcs}");

            let ignored = tempfile::tempdir().unwrap();
            let ignored_root = ignored.path();
            initialize_vcs(ignored_root, vcs);
            write(ignored_root, ".gitignore", b"nested/\n");
            write(
                ignored_root,
                "nested/ignored.html",
                b"ignored scratch prose",
            );
            assert!(
                controlled_git_command(&ignored_root.join("nested"))
                    .arg("init")
                    .status()
                    .unwrap()
                    .success()
            );
            let sources = RepositorySources::discover(ignored_root).unwrap();
            assert!(!sources.contains(Path::new("nested/ignored.html")));
        }
    }

    #[test]
    fn linked_worktree_directory_sentinel_fails_closed() {
        let temporary = tempfile::tempdir().unwrap();
        let source = temporary.path().join("source");
        let root = temporary.path().join("outer");
        fs::create_dir(&source).unwrap();
        fs::create_dir(&root).unwrap();
        initialize_vcs(&source, "git");
        assert!(
            controlled_git_command(&source)
                .args([
                    "-c",
                    "user.name=Docs Test",
                    "-c",
                    "user.email=docs@example.invalid"
                ])
                .args(["commit", "--allow-empty", "-m", "initial"])
                .status()
                .unwrap()
                .success()
        );
        initialize_vcs(&root, "git");
        let nested = root.join("nested-linked");
        assert!(
            controlled_git_command(&source)
                .args(["worktree", "add", "--detach"])
                .arg(&nested)
                .arg("HEAD")
                .status()
                .unwrap()
                .success()
        );
        write(&nested, "unauthorized.html", b"authored prose");
        assert!(nested.join(".git").is_file());
        assert!(RepositorySources::discover(&root).is_err());
    }

    #[test]
    fn tracked_gitlinks_fail_closed_in_git_and_jj() {
        for vcs in ["git", "jj"] {
            let temporary = tempfile::tempdir().unwrap();
            let root = temporary.path();
            initialize_vcs(root, "git");
            assert!(
                controlled_git_command(root)
                    .args([
                        "-c",
                        "user.name=Docs Test",
                        "-c",
                        "user.email=docs@example.invalid",
                    ])
                    .args(["commit", "--allow-empty", "-m", "initial"])
                    .status()
                    .unwrap()
                    .success()
            );
            let commit = controlled_git_command(root)
                .args(["rev-parse", "HEAD"])
                .output()
                .unwrap();
            assert!(commit.status.success(), "{commit:?}");
            let commit = String::from_utf8(commit.stdout).unwrap();
            assert!(
                controlled_git_command(root)
                    .args(["update-index", "--add", "--cacheinfo"])
                    .arg(format!("160000,{},nested-gitlink", commit.trim()))
                    .status()
                    .unwrap()
                    .success()
            );
            assert!(
                controlled_git_command(root)
                    .args([
                        "-c",
                        "user.name=Docs Test",
                        "-c",
                        "user.email=docs@example.invalid",
                    ])
                    .args(["commit", "-m", "add gitlink"])
                    .status()
                    .unwrap()
                    .success()
            );
            if vcs == "jj" {
                let initialized = controlled_jj_command(root)
                    .args(["git", "init", "--colocate", "."])
                    .output()
                    .unwrap();
                assert!(initialized.status.success(), "{initialized:?}");
            }
            assert!(RepositorySources::discover(root).is_err(), "{vcs}");
        }
    }

    fn initialize_vcs(root: &Path, vcs: &str) {
        let mut command = if vcs == "git" {
            let mut command = controlled_git_command(root);
            command.arg("init");
            command
        } else {
            let mut command = controlled_jj_command(root);
            command.args(["git", "init", "--no-colocate", "."]);
            command
        };
        let output = command.output().unwrap();
        assert!(output.status.success(), "{vcs}: {output:?}");
    }

    fn jj_state(root: &Path) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
        let run = |arguments: &[&str]| {
            let output = Command::new("jj")
                .args(arguments)
                .current_dir(root)
                .output()
                .unwrap();
            assert!(output.status.success(), "{output:?}");
            output.stdout
        };
        (
            run(&[
                "--ignore-working-copy",
                "log",
                "-r",
                "@",
                "--no-graph",
                "-T",
                "commit_id",
            ]),
            run(&["--ignore-working-copy", "op", "log", "-n", "1"]),
            run(&["--ignore-working-copy", "bookmark", "list"]),
        )
    }

    #[test]
    fn metadata_free_exports_conservatively_check_every_file() {
        let temporary = tempfile::tempdir().unwrap();
        let sources = RepositorySources::discover(temporary.path()).unwrap();
        assert!(sources.contains(Path::new("docs/unexpected.html")));
    }
}
