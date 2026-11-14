//! Engines on a RustRed p library of any scope, and the regeneration check of
//! a library directory (M3-V8).
//!
//! A library is a `rustred p-package` directory: the built-in one
//! (`data/p/<library>`, [`built_in_p_library`]), M3's scope 1-4 package
//! (`VAKINT_TEST_P_LIBRARY_1_4`, [`p_library_1_4`]) or any other. An engine
//! reaches it injected ([`PLibraryRoute::Injected`]) or listed in the
//! `p_integral/` area of a temporary external-data root, the way
//! vakint-external-data ships it ([`PLibraryRoute::Listed`]).
//!
//! Include it with `#[path = "support/p_library4.rs"] mod p_library4;` (it
//! needs no other support file).

#![allow(dead_code)]

use std::fs;
use std::ops::{Deref, DerefMut};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use rustred::propagator::custom::inspect_custom_payload;
use rustred::propagator::families::{P_FAMILIES_V1_JSON, PFamiliesFile};
use rustred::propagator::library::PLibrary;
use rustred_app::propagator::generate::custom_owner_lp_count;
use rustred_app::propagator::loader::read_custom_payloads;
use rustred_app::propagator::package::{assemble_p_package, p_library_differences};
use rustred_app::propagator::{
    AssembleOptions, CustomPayloads, PLibraryOptions, PLibrarySource, equivalent_p_libraries,
    load_p_library,
};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use vakint::Vakint;

/// The scope 1-4 package of M3 (set in `ops/env.sh` from the PP-M3b pin).
pub const P_LIBRARY_1_4_VARIABLE: &str = "VAKINT_TEST_P_LIBRARY_1_4";
/// `1`: recount the LP of every Custom owner (an 11-line owner takes about
/// 24 min); otherwise only owners of at most [`RECOUNT_MAX_LINES`] lines.
pub const RECOUNT_ALL_VARIABLE: &str = "P_LP_RECOUNT_ALL";
pub const RECOUNT_MAX_LINES: u32 = 9;

fn data() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data")
}

fn json_file(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).unwrap())
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// The built-in library directory `data/p/<library>` that
/// `data/p/manifest.json` names.
pub fn built_in_p_library() -> PathBuf {
    let manifest = json_file(&data().join("p/manifest.json"));
    data().join("p").join(manifest["library"].as_str().unwrap())
}

/// The scope 1-4 package, if [`P_LIBRARY_1_4_VARIABLE`] is set.
pub fn p_library_1_4() -> Option<PathBuf> {
    std::env::var_os(P_LIBRARY_1_4_VARIABLE).map(PathBuf::from)
}

/// [`p_library_1_4`] for a test that cannot run without it (such tests are
/// ignored and run explicitly, so a missing package is a failure).
pub fn require_p_library_1_4(test: &str) -> PathBuf {
    p_library_1_4().unwrap_or_else(|| panic!("{test} needs {P_LIBRARY_1_4_VARIABLE}"))
}

/// `[a, b]` of the library in `directory` (its `library.json`).
pub fn p_library_scope(directory: &Path) -> [u8; 2] {
    serde_json::from_value(json_file(&directory.join("library.json"))["scope"].clone()).unwrap()
}

/// `p-library-<a>-<b>` of the library in `directory`.
pub fn p_library_name(directory: &Path) -> String {
    let [start, end] = p_library_scope(directory);
    format!("p-library-{start}-{end}")
}

/// The library in `directory`, loaded without the descent certificate.
pub fn load_p_library_directory(directory: &Path) -> Arc<PLibrary> {
    load_p_library(
        PLibrarySource::Directory(directory.to_owned()),
        &PLibraryOptions {
            verify_descent: false,
        },
    )
    .unwrap_or_else(|error| panic!("{}: {error}", directory.display()))
}

/// How [`engine_with_p_library`] gives the engine its library.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PLibraryRoute {
    /// `Vakint::inject_p_library` (its own step-cache directory, no shipped cache).
    Injected,
    /// Listed in the `p_integral/` area of a temporary external-data root
    /// ([`listed_p_library_root`]), loaded by `load_external_data`.
    Listed,
}

/// An engine and, on the listed route, the external-data root it reads.
pub struct PLibraryEngine {
    pub engine: Vakint,
    pub root: Option<tempfile::TempDir>,
}

impl Deref for PLibraryEngine {
    type Target = Vakint;

    fn deref(&self) -> &Vakint {
        &self.engine
    }
}

impl DerefMut for PLibraryEngine {
    fn deref_mut(&mut self) -> &mut Vakint {
        &mut self.engine
    }
}

/// An engine whose p-integrals within the scope of the library in
/// `directory` go to that library (any p-integral of that scope under the
/// default order).
pub fn engine_with_p_library(directory: &Path, route: PLibraryRoute) -> PLibraryEngine {
    let mut engine = Vakint::new().unwrap();
    match route {
        PLibraryRoute::Injected => {
            engine
                .inject_p_library(load_p_library_directory(directory))
                .unwrap();
            PLibraryEngine { engine, root: None }
        }
        PLibraryRoute::Listed => {
            let root = listed_p_library_root(directory);
            engine.load_external_data(root.path(), 3).unwrap();
            PLibraryEngine {
                engine,
                root: Some(root),
            }
        }
    }
}

fn copy_tree(from: &Path, to: &Path, skip: &dyn Fn(&Path) -> bool) {
    for entry in fs::read_dir(from).unwrap() {
        let path = entry.unwrap().path();
        let target = to.join(path.file_name().unwrap());
        if skip(&path) {
            continue;
        }
        if path.is_dir() {
            fs::create_dir_all(&target).unwrap();
            copy_tree(&path, &target, skip);
        } else {
            fs::copy(&path, &target).unwrap();
        }
    }
}

/// A temporary external-data root: vakint's `data/` without its `p/` area
/// (the vacuum 1L-3L packs, so `load_external_data(root, 3)` finds the root),
/// and `p_integral/p-library-<a>-<b>/`, a copy of `directory`, listed in
/// `p_integral/manifest.json` as `tools/p_integral_area.py stage-library` of
/// vakint-external-data lists it, with the root key.
pub fn listed_p_library_root(directory: &Path) -> tempfile::TempDir {
    let root = tempfile::tempdir().unwrap();
    let p_area = data().join("p");
    copy_tree(&data(), root.path(), &|path| path == p_area);
    let name = p_library_name(directory);
    let target = root.path().join("p_integral").join(&name);
    fs::create_dir_all(&target).unwrap();
    copy_tree(directory, &target, &|_| false);
    let mut files = serde_json::Map::new();
    let mut stack = vec![target.clone()];
    while let Some(next) = stack.pop() {
        for entry in fs::read_dir(next).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                stack.push(path);
                continue;
            }
            let bytes = fs::read(&path).unwrap();
            let relative = path.strip_prefix(&target).unwrap().to_str().unwrap();
            files.insert(
                relative.replace('\\', "/"),
                json!({"sha256": format!("{:x}", Sha256::digest(&bytes)), "bytes": bytes.len()}),
            );
        }
    }
    let description = json_file(&directory.join("library.json"));
    let area = json!({
        "schema": "vakint.external-data.p-integral.v1",
        "libraries": [{
            "directory": name,
            "schema": description["schema"],
            "scope": description["scope"],
            "identity": description["identity"],
            "files": files,
        }],
        "files": {},
    });
    let bytes = serde_json::to_vec_pretty(&area).unwrap();
    fs::write(root.path().join("p_integral/manifest.json"), &bytes).unwrap();
    let release_path = root.path().join("manifest.json");
    let mut release = json_file(&release_path);
    release["p_integral"] = json!({
        "manifest": "p_integral/manifest.json",
        "manifest_sha256": format!("{:x}", Sha256::digest(&bytes)),
    });
    fs::write(&release_path, serde_json::to_vec_pretty(&release).unwrap()).unwrap();
    root
}

/// The LP count of one Custom owner of a library: the one its record stores
/// and, for an owner the recount takes, the one counted again.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OwnerLp {
    pub owner: String,
    pub lines: u32,
    pub stored: u32,
    pub recounted: Option<u32>,
}

/// What [`regenerate_p_library`] checked, with its wall seconds.
#[derive(Clone, Debug)]
pub struct Regeneration {
    pub scope: [u8; 2],
    pub certificate_seconds: f64,
    pub assembly_seconds: f64,
    pub owners: Vec<OwnerLp>,
}

/// The regeneration check of the library in `directory` (plan M3 D-M3-3,
/// I-8, D-V4-2):
/// 1. it loads WITH the descent certificate (the once-per-library run: the
///    certificate, the census's strong completeness and the basis layer
///    derived again);
/// 2. the LP count each Custom payload records is counted again for owners of
///    at most [`RECOUNT_MAX_LINES`] lines (every owner when
///    [`RECOUNT_ALL_VARIABLE`] is `1`), as generation counts it
///    (`custom_owner_lp_count` on the owner presentation and the payload's
///    support; rustred's `RecountProvidedLp` filter once the pin has it);
/// 3. the current RustRed assembles, from RustRed's families and the
///    library's own payloads (`Provided`, stored LP), a library equivalent to
///    it (`equivalent_p_libraries`), on `threads` census threads.
pub fn regenerate_p_library(directory: &Path, threads: usize) -> Regeneration {
    let source = PLibrarySource::Directory(directory.to_owned());
    let started = Instant::now();
    let library = load_p_library(
        source.clone(),
        &PLibraryOptions {
            verify_descent: true,
        },
    )
    .unwrap_or_else(|error| panic!("{} with the certificate: {error}", directory.display()));
    let certificate_seconds = started.elapsed().as_secs_f64();

    let recount_all = std::env::var(RECOUNT_ALL_VARIABLE).as_deref() == Ok("1");
    let mut owners = Vec::new();
    for owner in library.owners() {
        let (Some(meta), Some(payload)) = (&owner.spec.custom, &owner.spec.payload) else {
            continue;
        };
        let name = format!("{:?}", owner.spec.class);
        let support = inspect_custom_payload(&payload.bytes().unwrap())
            .unwrap()
            .support;
        let lines = support.count_ones();
        let recounted = (recount_all || lines <= RECOUNT_MAX_LINES).then(|| {
            custom_owner_lp_count(library.presentation(owner.spec.presentation), support)
                .unwrap_or_else(|error| panic!("{name}: {error}"))
        });
        assert!(
            recounted.is_none_or(|count| count == meta.lp_count),
            "{name} ({lines} lines) records LP {}, counted {recounted:?}",
            meta.lp_count
        );
        owners.push(OwnerLp {
            owner: name,
            lines,
            stored: meta.lp_count,
            recounted,
        });
    }

    let families: PFamiliesFile = serde_json::from_str(P_FAMILIES_V1_JSON).unwrap();
    let payloads = read_custom_payloads(&source).unwrap();
    assert_eq!(payloads.len(), owners.len(), "Custom payloads and owners");
    let started = Instant::now();
    let assembled = assemble_p_package(
        &families,
        *library.scope().end(),
        CustomPayloads::Provided(payloads),
        &AssembleOptions {
            threads,
            verify_descent: false,
        },
    )
    .unwrap();
    let assembly_seconds = started.elapsed().as_secs_f64();
    assert_equivalent_up_to_revision(&library, &assembled.library);
    let scope = library.scope();
    Regeneration {
        scope: [*scope.start(), *scope.end()],
        certificate_seconds,
        assembly_seconds,
        owners,
    }
}

/// `shipped` is the program of `assembled` (`equivalent_p_libraries`). Batch 6:
/// a shipped package may carry a content revision (`,P_PACKAGE_REVISION=<name>`,
/// assembled under RUSTRED_P_PACKAGE_REVISION) that an in-process assembly does
/// not have; then the identities must agree without the revision (and the
/// ordering tag, which p_library_differences ignores) and nothing else may differ.
pub fn assert_equivalent_up_to_revision(shipped: &PLibrary, assembled: &PLibrary) {
    if shipped.revision().is_none() {
        assert!(
            equivalent_p_libraries(shipped, assembled).unwrap(),
            "{:?}",
            p_library_differences(shipped, assembled).unwrap()
        );
        return;
    }
    let bare = |l: &PLibrary| {
        let mut identity = l.identity();
        if let Some(tag) = l.ordering() {
            identity = identity.replacen(&format!(",P_ORDERING_VERSION={tag}"), "", 1);
        }
        if let Some(name) = l.revision() {
            identity = identity.replacen(&format!(",P_PACKAGE_REVISION={name}"), "", 1);
        }
        identity
    };
    assert_eq!(bare(shipped), bare(assembled));
    let differences: Vec<String> = p_library_differences(shipped, assembled)
        .unwrap()
        .into_iter()
        .filter(|difference| !difference.starts_with("identity: "))
        .collect();
    assert!(differences.is_empty(), "{differences:?}");
}
