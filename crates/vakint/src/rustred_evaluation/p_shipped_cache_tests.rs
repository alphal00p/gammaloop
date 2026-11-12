//! The shipped p cache (M2-V6) on a synthetic scope-1 area: RustRed's one-loop
//! reference library plays a non-injected library (as one from external data
//! or the built-in one), a cold read-write run records its steps, the cache is
//! exported and staged into a writable copy of the vendored release the way
//! `tools/p_integral_area.py stage` of vakint-external-data does. The area
//! round-trips in place and installed; every rejected item runs cold with the
//! cache-off values; the none cases ship nothing.
use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
    sync::Arc,
};

use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use symbolica::atom::AtomCore;
use tempfile::TempDir;

use super::{
    p_data::{PLibraryHandle, PLibraryOrigin},
    p_shipped_cache::{AREA_SCHEMA, METADATA_SCHEMA},
};
use crate::{
    EvaluationMethod, EvaluationOrder, OrderedCacheMode, OrderedCacheProgramStatistics,
    OrderedCacheSettings, RustRedEvaluationOptions, RustRedMasterMode, ShippedOrderedCacheStatus,
    Vakint, VakintSettings, vakint_parse,
};

const INPUTS: [&str; 2] = [
    "(dot(k(1),p(1))^3)*topo(P1L(p(1),2,1))",
    "(1)*topo(P1L(p(1),2,2))",
];

fn settings(ordered_cache: OrderedCacheSettings) -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Retain,
            },
        )]),
        form_exe_path: "/p-shipped-cache-tests-must-not-call-form".into(),
        use_dot_product_notation: true,
        n_cores: 1,
        ordered_cache,
        ..Default::default()
    }
}

fn read_only_in_memory() -> OrderedCacheSettings {
    OrderedCacheSettings {
        mode: OrderedCacheMode::ReadOnly,
        ..Default::default()
    }
}

/// A writable copy of the vendored 1L-3L release.
fn data_root() -> TempDir {
    let root = tempfile::tempdir().unwrap();
    for (relative, bytes) in super::data::vendored_files() {
        let path = root.path().join(relative);
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(path, bytes).unwrap();
    }
    root
}

/// An engine reading its one-loop data from `root` whose p library is the
/// reference library, injected or as non-injected data.
fn engine(root: &Path, injected: bool) -> Vakint {
    let mut engine = Vakint::new().unwrap();
    engine.load_external_data(root, 1).unwrap();
    let library = rustred::propagator::testing::one_loop_reference_library().unwrap();
    if injected {
        engine.inject_p_library(library).unwrap();
    } else {
        let handle = PLibraryHandle::new(library, PLibraryOrigin::ExternalData).unwrap();
        engine.p_library = Some(Arc::new(handle));
    }
    engine
}

fn values(engine: &Vakint, settings: &VakintSettings) -> Vec<String> {
    INPUTS
        .iter()
        .map(|input| {
            engine
                .evaluate(settings, vakint_parse!(input).unwrap().as_view())
                .unwrap()
                .to_canonical_string()
        })
        .collect()
}

/// The statistics row of the p program.
fn p_row(engine: &Vakint) -> OrderedCacheProgramStatistics {
    engine
        .ordered_cache_statistics()
        .programs
        .into_iter()
        .find(|program| program.program == "p-library-1-1")
        .expect("the p program is attached")
}

fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn files_below(directory: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    for entry in fs::read_dir(directory).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() {
            out.extend(files_below(&path));
        } else {
            out.push(path);
        }
    }
    out.sort();
    out
}

/// The cache-off values of [`INPUTS`] and the export of the p cache a cold
/// read-write run recorded.
struct Exported {
    cold: Vec<String>,
    export: TempDir,
}

fn exported() -> Exported {
    let root = data_root();
    let engine = engine(root.path(), false);
    let cold = values(&engine, &settings(OrderedCacheSettings::off()));
    let cache = tempfile::tempdir().unwrap();
    let recording = OrderedCacheSettings {
        shipped: false,
        ..OrderedCacheSettings::read_write(cache.path())
    };
    assert_eq!(values(&engine, &settings(recording)), cold);
    let export = tempfile::tempdir().unwrap();
    let report = engine
        .export_ordered_caches(
            cache.path(),
            &export.path().join("out"),
            json!({"command": "p_shipped_cache_tests"}),
            vec![json!({"name": "bubbles", "loop_count": 1})],
        )
        .unwrap();
    let entry = report
        .p_integral
        .as_ref()
        .expect("the cold run recorded the p cache");
    assert!(entry.records > 0 && entry.shards > 0, "{entry:?}");
    assert!(report.packs.is_empty(), "{report:?}");
    assert!(
        export
            .path()
            .join("out/p_integral/ordered-cache.json")
            .is_file()
    );
    Exported { cold, export }
}

/// Stage the export into a fresh data root like `tools/p_integral_area.py
/// stage`: `p_integral/`, its manifest listing every file with sha256 and
/// bytes, and the root key. `metadata` and `area` edit `ordered-cache.json`
/// and the area manifest before they are hashed.
fn staged(
    exported: &Exported,
    metadata: impl FnOnce(&mut Value),
    area: impl FnOnce(&mut Value),
) -> TempDir {
    let root = data_root();
    let source = exported.export.path().join("out/p_integral");
    let target = root.path().join("p_integral");
    for file in files_below(&source) {
        let to = target.join(file.strip_prefix(&source).unwrap());
        fs::create_dir_all(to.parent().unwrap()).unwrap();
        fs::copy(&file, &to).unwrap();
    }
    let metadata_path = target.join("ordered-cache.json");
    let mut json: Value = serde_json::from_slice(&fs::read(&metadata_path).unwrap()).unwrap();
    assert_eq!(json["schema"], METADATA_SCHEMA);
    metadata(&mut json);
    fs::write(&metadata_path, serde_json::to_vec_pretty(&json).unwrap()).unwrap();
    let files: BTreeMap<String, Value> = files_below(&target)
        .into_iter()
        .map(|file| {
            let bytes = fs::read(&file).unwrap();
            let relative = file.strip_prefix(&target).unwrap().to_str().unwrap();
            let record = json!({"sha256": sha256(&bytes), "bytes": bytes.len()});
            (relative.to_owned(), record)
        })
        .collect();
    let mut manifest = json!({
        "schema": AREA_SCHEMA,
        "libraries": [],
        "ordered_cache": {"path": "ordered-cache", "metadata": "ordered-cache.json"},
        "files": files,
    });
    area(&mut manifest);
    let manifest = serde_json::to_vec_pretty(&manifest).unwrap();
    fs::write(target.join("manifest.json"), &manifest).unwrap();
    let root_manifest = root.path().join("manifest.json");
    let mut release: Value = serde_json::from_slice(&fs::read(&root_manifest).unwrap()).unwrap();
    release["p_integral"] = json!({
        "manifest": "p_integral/manifest.json",
        "manifest_sha256": sha256(&manifest),
    });
    fs::write(&root_manifest, serde_json::to_vec_pretty(&release).unwrap()).unwrap();
    root
}

/// Step-record shards and tail entries the p cache under `settings` read
/// from its files. These one-loop inputs replay from RustRed's tail table
/// (complete reductions of tail identities) without any step lookup.
fn loaded_entries(engine: &Vakint, settings: &VakintSettings) -> usize {
    let handle = engine.p_library().unwrap();
    let cache = engine
        .ordered_step_cache(settings, &handle.cache)
        .unwrap()
        .unwrap();
    let tail = cache
        .tail_statistics()
        .map_or(0, |tail| tail.loaded_entries);
    cache.statistics().loaded_shards + tail
}

/// L3: the staged area replays the cold values in place (read-only, no cache
/// directory, nothing written into the data) and installed (read-write, an
/// empty cache directory receives a copy), reading the shipped files and
/// without a single miss.
#[test]
fn the_p_cache_round_trips_in_place_and_installed() {
    let exported = exported();
    let root = staged(&exported, |_| {}, |_| {});
    let area = root.path().join("p_integral");
    let digests = || {
        files_below(&area)
            .iter()
            .map(|file| sha256(&fs::read(file).unwrap()))
            .collect::<Vec<_>>()
    };
    let before = digests();

    let in_place = engine(root.path(), false);
    let read_only = settings(read_only_in_memory());
    assert_eq!(values(&in_place, &read_only), exported.cold);
    let row = p_row(&in_place);
    let shipped = row.shipped.as_ref().unwrap();
    assert_eq!(
        shipped.status,
        ShippedOrderedCacheStatus::InPlace,
        "{row:?}"
    );
    assert!(shipped.source.starts_with(&area));
    assert_eq!((row.native.misses, row.native.inserts), (0, 0), "{row:?}");
    assert!(loaded_entries(&in_place, &read_only) > 0, "{row:?}");
    assert_eq!(digests(), before);

    let directory = tempfile::tempdir().unwrap();
    let installed = engine(root.path(), false);
    let read_write = settings(OrderedCacheSettings::read_write(directory.path()));
    assert_eq!(values(&installed, &read_write), exported.cold);
    let row = p_row(&installed);
    let shipped = row.shipped.as_ref().unwrap();
    assert_eq!(
        shipped.status,
        ShippedOrderedCacheStatus::Installed,
        "{row:?}"
    );
    assert_eq!((row.native.misses, row.native.inserts), (0, 0), "{row:?}");
    assert!(loaded_entries(&installed, &read_write) > 0, "{row:?}");
    let program = row.directory.unwrap();
    assert!(program.starts_with(directory.path()));
    assert!(
        program
            .join(super::shipped_cache::RUSTRED_MANIFEST)
            .is_file()
    );
    assert_eq!(digests(), before);
}

/// L2: a cache of another stamp or library identity, a flipped byte and an
/// unknown area schema are rejected (the reason is reported, the stderr line
/// printed) and the program runs cold with the cache-off values.
#[test]
fn a_mismatching_or_corrupt_p_cache_is_rejected_and_runs_cold() {
    let exported = exported();
    let flip = |root: &Path| {
        let shards = root.join("p_integral/ordered-cache");
        let shard = files_below(&shards)
            .into_iter()
            .find(|file| {
                file.extension()
                    .is_some_and(|extension| extension == "rrstep")
            })
            .unwrap();
        let mut bytes = fs::read(&shard).unwrap();
        let last = bytes.len() - 1;
        bytes[last] ^= 1;
        fs::write(shard, bytes).unwrap();
    };
    let cases: [(&str, TempDir, &str); 4] = [
        (
            "stamp",
            staged(
                &exported,
                |m| m["programs"][0]["stamp"] = "0".repeat(64).into(),
                |_| {},
            ),
            "another cache version",
        ),
        (
            "identity",
            staged(
                &exported,
                |m| m["programs"][0]["identity"] = "p-library:scope=1-1:P_OTHER=1".into(),
                |_| {},
            ),
            "prepared for the library",
        ),
        (
            "flipped byte",
            {
                let root = staged(&exported, |_| {}, |_| {});
                flip(root.path());
                root
            },
            "does not match its sha256",
        ),
        (
            "unknown schema",
            staged(
                &exported,
                |_| {},
                |a| a["schema"] = "vakint.external-data.p-integral.v2".into(),
            ),
            "has schema",
        ),
    ];
    for (name, root, reason) in cases {
        let engine = engine(root.path(), false);
        assert_eq!(
            values(&engine, &settings(OrderedCacheSettings::default())),
            exported.cold,
            "{name}"
        );
        let row = p_row(&engine);
        let shipped = row.shipped.as_ref().unwrap();
        assert_eq!(
            shipped.status,
            ShippedOrderedCacheStatus::Rejected,
            "{name}"
        );
        let detail = shipped.detail.as_deref().unwrap();
        assert!(detail.contains(reason), "{name}: {detail}");
        assert!(row.native.misses > 0, "{name}: {row:?}");
    }
}

/// L2: no shipped p cache without an area, with a root key whose area is
/// absent (as packed data keep it), or for an injected library.
#[test]
fn without_an_area_or_for_an_injected_library_no_p_cache_is_shipped() {
    let exported = exported();
    let key_without_area = staged(&exported, |_| {}, |_| {});
    fs::remove_dir_all(key_without_area.path().join("p_integral")).unwrap();
    let cases = [
        ("no area", data_root(), false),
        ("key without area", key_without_area, false),
        ("injected", staged(&exported, |_| {}, |_| {}), true),
    ];
    for (name, root, injected) in cases {
        let engine = engine(root.path(), injected);
        assert_eq!(
            values(&engine, &settings(OrderedCacheSettings::default())),
            exported.cold,
            "{name}"
        );
        let row = p_row(&engine);
        assert_eq!(row.shipped, None, "{name}");
        assert!(row.native.misses > 0, "{name}: {row:?}");
    }
}

/// L1: once the p library is loaded, `shipped_ordered_caches` and
/// `load_ordered_cache` list the p program with its shipped cache in place
/// (`export_ordered_caches` lists it in every [`exported`]).
#[test]
fn the_public_entry_points_list_the_loaded_p_program() {
    let exported = exported();
    let root = staged(&exported, |_| {}, |_| {});
    let engine = engine(root.path(), false);
    let checks = engine.shipped_ordered_caches().unwrap();
    let check = checks
        .iter()
        .find(|check| check.program == "p-library-1-1")
        .expect("shipped_ordered_caches lists the p program");
    assert_eq!(check.loop_count, 1);
    assert_eq!(
        check.shipped.as_ref().map(|shipped| shipped.status),
        Some(ShippedOrderedCacheStatus::InPlace),
        "{check:?}"
    );
    let statistics = engine
        .load_ordered_cache(&settings(read_only_in_memory()))
        .unwrap();
    let row = statistics
        .programs
        .iter()
        .find(|row| row.program == "p-library-1-1")
        .expect("load_ordered_cache attaches the p program");
    assert_eq!(
        row.shipped.as_ref().map(|shipped| shipped.status),
        Some(ShippedOrderedCacheStatus::InPlace),
        "{row:?}"
    );
}
