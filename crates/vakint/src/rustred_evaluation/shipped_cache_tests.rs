//! Shipped Ordered step caches: caches exported from a cold run and staged
//! into a data directory replay the cold results exactly in the default
//! in-memory mode (opened in place, nothing written into the data) and from a
//! cache directory (installed copy); caches of another version or rule tables
//! and corrupt files are rejected and the evaluation runs cold.
use std::{
    fs,
    path::{Path, PathBuf},
};

use sha2::{Digest, Sha256};
use symbolica::atom::Atom;

use super::*;
use crate::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, OrderedCacheMode,
    OrderedCacheProgramStatistics, OrderedCacheSettings, RustRedEvaluationOptions,
    RustRedMasterMode, Vakint, VakintSettings, vakint_parse,
};

const INPUTS: [&str; 3] = [
    "topo(I3L(muvsq,2,2,2,1,1,1))",
    "dot(k(1),k(2))*dot(k(2),k(3))*topo(I3L(muvsq,2,1,1,1,1,1))",
    "topo(I2L(muvsq,2,1,2))",
];

fn settings(cache: OrderedCacheSettings) -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Retain,
            },
        )]),
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        use_dot_product_notation: true,
        form_exe_path: "/shipped-cache-tests-must-not-call-form".into(),
        ordered_cache: cache,
        ..Default::default()
    }
}

fn evaluate_all(engine: &Vakint, settings: &VakintSettings) -> Vec<Atom> {
    INPUTS
        .iter()
        .map(|input| {
            engine
                .evaluate(settings, vakint_parse!(input).unwrap().as_view())
                .unwrap()
        })
        .collect()
}

/// A writable copy of the vendored 1L-3L release.
fn data_directory() -> tempfile::TempDir {
    let directory = tempfile::tempdir().unwrap();
    for (relative, bytes) in super::super::data::vendored_files() {
        let path = directory.path().join(relative);
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(path, bytes).unwrap();
    }
    directory
}

fn load_engine(data: &Path) -> Vakint {
    let mut engine = Vakint::new().unwrap();
    engine.load_external_data(data, 3).unwrap();
    engine
}

fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn walkdir(directory: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    for entry in fs::read_dir(directory).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() {
            out.extend(walkdir(&path));
        } else {
            out.push(path);
        }
    }
    out.sort();
    out
}

/// Every file under `root` with its bytes.
fn listing(root: &Path) -> Vec<(PathBuf, Vec<u8>)> {
    walkdir(root)
        .into_iter()
        .map(|path| {
            (
                path.strip_prefix(root).unwrap().to_owned(),
                fs::read(&path).unwrap(),
            )
        })
        .collect()
}

/// Refresh every `files` record of the pack's `ordered-cache*` files and the
/// release manifest's digest of the loop manifest, as the data repository's
/// tooling does after staging.
fn rebind(data: &Path, loops: usize) {
    let pack = data.join(format!("{loops}L"));
    let manifest_path = pack.join("manifest.json");
    let mut manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(&manifest_path).unwrap()).unwrap();
    let files = manifest["files"].as_object_mut().unwrap();
    files.retain(|name, _| !name.starts_with("ordered-cache"));
    for path in walkdir(&pack) {
        let relative = path
            .strip_prefix(&pack)
            .unwrap()
            .to_str()
            .unwrap()
            .to_owned();
        if relative.starts_with("ordered-cache") {
            let bytes = fs::read(&path).unwrap();
            files.insert(
                relative,
                serde_json::json!({"sha256": sha256(&bytes), "bytes": bytes.len()}),
            );
        }
    }
    manifest["ordered_cache"] =
        serde_json::json!({"path": "ordered-cache", "metadata": "ordered-cache.json"});
    let bytes = serde_json::to_vec_pretty(&manifest).unwrap();
    fs::write(&manifest_path, &bytes).unwrap();
    let root_path = data.join("manifest.json");
    let mut root: serde_json::Value =
        serde_json::from_slice(&fs::read(&root_path).unwrap()).unwrap();
    for entry in root["loop_packs"].as_array_mut().unwrap() {
        if entry["loop_count"] == loops {
            entry["manifest_sha256"] = sha256(&bytes).into();
        }
    }
    fs::write(&root_path, serde_json::to_vec_pretty(&root).unwrap()).unwrap();
}

/// Copy an export into the data packs and rebind their manifests.
fn stage(data: &Path, export: &Path) {
    for loops in 1..=3 {
        let source = export.join(format!("{loops}L"));
        if !source.is_dir() {
            continue;
        }
        for (relative, bytes) in listing(&source) {
            let target = data.join(format!("{loops}L")).join(relative);
            fs::create_dir_all(target.parent().unwrap()).unwrap();
            fs::write(target, bytes).unwrap();
        }
        rebind(data, loops);
    }
}

/// Cold results, and a staged data directory shipping the caches they
/// recorded (2L and 3L programs).
fn shipping_data(work: &Path) -> (tempfile::TempDir, Vec<Atom>) {
    let data = data_directory();
    let engine = load_engine(data.path());
    let root = work.join("cold-root");
    let cold = evaluate_all(
        &engine,
        &settings(OrderedCacheSettings {
            shipped: false,
            ..OrderedCacheSettings::read_write(&root)
        }),
    );
    engine.save_ordered_cache().unwrap();
    let export = |output: &Path| {
        engine
            .export_ordered_caches(
                &root,
                output,
                ordered_cache_producer("test", "test", "shipped_cache_tests"),
                vec![serde_json::json!({"inputs": INPUTS})],
            )
            .unwrap()
    };
    let report = export(&work.join("export"));
    // The 2L and 3L programs were used (the 1L one may be, through
    // factorized products).
    assert!(
        report.packs.contains_key(&2) && report.packs.contains_key(&3),
        "{report:?}"
    );
    assert_eq!(
        report.packs.values().map(Vec::len).sum::<usize>() + report.missing.len(),
        3
    );
    // Exports are deterministic.
    export(&work.join("export-again"));
    assert_eq!(
        listing(&work.join("export")),
        listing(&work.join("export-again"))
    );
    let program = &report.packs[&3][0];
    assert!(program.complete && program.records > 0 && program.tail_bytes == 0);
    assert_eq!(
        program.directory,
        format!(
            "{PACK_DIRECTORY}/{}.{}.{}",
            program.label,
            &program.rule_tables[..16],
            &program.stamp[..16]
        )
    );
    stage(data.path(), &work.join("export"));
    (data, cold)
}

/// The number of programs of the vendored packs that ship no cache.
fn report_missing(data: &Path) -> usize {
    (1..=3)
        .filter(|loops| {
            !data
                .join(format!("{loops}L/{METADATA_FILE}"))
                .is_file()
        })
        .count()
}

fn program(engine: &Vakint, loops: usize) -> OrderedCacheProgramStatistics {
    engine
        .ordered_cache_statistics()
        .programs
        .into_iter()
        .find(|program| program.loop_count == loops)
        .unwrap()
}

fn status(program: &OrderedCacheProgramStatistics) -> Option<ShippedOrderedCacheStatus> {
    program.shipped.as_ref().map(|shipped| shipped.status)
}

#[test]
fn shipped_caches_replay_in_memory_and_from_a_cache_directory() {
    let work = tempfile::tempdir().unwrap();
    let (data, cold) = shipping_data(work.path());
    let before = listing(data.path());

    // Default settings: in memory, shipped caches opened in place.
    let engine = load_engine(data.path());
    assert_eq!(
        engine
            .shipped_ordered_caches()
            .unwrap()
            .iter()
            .filter(|check| check
                .shipped
                .as_ref()
                .is_some_and(|shipped| shipped.status == ShippedOrderedCacheStatus::InPlace))
            .count(),
        3 - report_missing(data.path())
    );
    let warm = evaluate_all(&engine, &settings(OrderedCacheSettings::default()));
    assert_eq!(warm, cold);
    for loops in [2, 3] {
        let program = program(&engine, loops);
        assert_eq!(
            status(&program),
            Some(ShippedOrderedCacheStatus::InPlace),
            "{program:?}"
        );
        assert_eq!(program.directory, None);
        assert!(program.native.hits > 0, "{program:?}");
        assert_eq!(program.native.inserts, 0, "{program:?}");
    }
    // New records stay in memory; the data are never written.
    let other = engine
        .evaluate(
            &settings(OrderedCacheSettings::default()),
            vakint_parse!("topo(I3L(muvsq,3,2,2,1,1,1))")
                .unwrap()
                .as_view(),
        )
        .unwrap();
    assert!(program(&engine, 3).native.inserts > 0);
    engine.save_ordered_cache().unwrap();
    drop(engine);
    assert_eq!(listing(data.path()), before);

    // Read-only without a directory replays the shipped records only.
    let engine = load_engine(data.path());
    let read_only = OrderedCacheSettings {
        mode: OrderedCacheMode::ReadOnly,
        ..OrderedCacheSettings::default()
    };
    assert_eq!(evaluate_all(&engine, &settings(read_only)), cold);
    assert!(program(&engine, 3).native.hits > 0);
    drop(engine);

    // A cache directory gets an installed copy, byte for byte.
    let root = work.path().join("warm-root");
    let engine = load_engine(data.path());
    let directory = settings(OrderedCacheSettings::read_write(&root));
    assert_eq!(evaluate_all(&engine, &directory), cold);
    let installed = program(&engine, 3);
    assert_eq!(
        status(&installed),
        Some(ShippedOrderedCacheStatus::Installed)
    );
    assert!(installed.native.hits > 0 && installed.native.inserts == 0);
    let name = installed
        .directory
        .as_ref()
        .unwrap()
        .file_name()
        .unwrap()
        .to_owned();
    assert_eq!(
        listing(&root.join(&name)),
        listing(&data.path().join("3L/ordered-cache").join(&name))
    );
    // Later records go to the directory, never to the data.
    assert_eq!(
        engine
            .evaluate(
                &directory,
                vakint_parse!("topo(I3L(muvsq,3,2,2,1,1,1))")
                    .unwrap()
                    .as_view()
            )
            .unwrap(),
        other
    );
    engine.save_ordered_cache().unwrap();
    drop(engine);
    assert_eq!(listing(data.path()), before);
    // An existing cache is kept.
    let engine = load_engine(data.path());
    assert_eq!(evaluate_all(&engine, &directory), cold);
    assert_eq!(
        status(&program(&engine, 3)),
        Some(ShippedOrderedCacheStatus::AlreadyPresent)
    );

    // `shipped: false` ignores them.
    let engine = load_engine(data.path());
    let ignored = settings(OrderedCacheSettings {
        shipped: false,
        ..OrderedCacheSettings::default()
    });
    assert_eq!(evaluate_all(&engine, &ignored), cold);
    assert_eq!(status(&program(&engine, 3)), None);
    assert!(program(&engine, 3).native.inserts > 0);
}

/// Edit the staged 3L cache with `edit(pack directory, cache directory)`,
/// optionally rebinding the manifests, and check that the cache is rejected
/// with a reason containing `reason` while the results stay the cold ones.
fn assert_rejected(edit: impl Fn(&Path, &Path), rebound: bool, reason: &str) {
    let work = tempfile::tempdir().unwrap();
    let (data, cold) = shipping_data(work.path());
    let pack = data.path().join("3L");
    let metadata: ShippedOrderedCacheMetadata =
        serde_json::from_slice(&fs::read(pack.join(METADATA_FILE)).unwrap()).unwrap();
    edit(&pack, &pack.join(&metadata.programs[0].directory));
    if rebound {
        rebind(data.path(), 3);
    }
    for cache in [
        OrderedCacheSettings::default(),
        OrderedCacheSettings::read_write(work.path().join("root")),
    ] {
        let engine = load_engine(data.path());
        let directory = cache.directory.clone();
        assert_eq!(evaluate_all(&engine, &settings(cache)), cold);
        let program = program(&engine, 3);
        let shipped = program.shipped.as_ref().unwrap();
        assert_eq!(shipped.status, ShippedOrderedCacheStatus::Rejected);
        assert!(
            shipped.detail.as_deref().unwrap().contains(reason),
            "{shipped:?}"
        );
        // Cold: every record is recomputed.
        assert!(program.native.inserts > 0, "{program:?}");
        // The 2L cache is unaffected.
        assert_eq!(
            status(&self::program(&engine, 2)),
            Some(if directory.is_some() {
                ShippedOrderedCacheStatus::Installed
            } else {
                ShippedOrderedCacheStatus::InPlace
            })
        );
    }
}

/// Replace `from` by `to` in a text file.
fn replace(path: &Path, from: &str, to: &str) {
    let text = fs::read_to_string(path).unwrap();
    assert!(text.contains(from), "{} lacks {from:?}", path.display());
    fs::write(path, text.replace(from, to)).unwrap();
}

#[test]
fn stale_or_corrupt_shipped_caches_are_rejected_and_run_cold() {
    // Another cache version: the stamp of manifest and metadata edited
    // consistently (sha256 records refreshed).
    assert_rejected(
        |pack, cache| {
            let manifest = cache.join("ordered-step-cache.txt");
            let text = fs::read_to_string(&manifest).unwrap();
            let stamp = text
                .lines()
                .find_map(|line| line.strip_prefix("stamp "))
                .unwrap()
                .to_owned();
            let edited = format!("{}{}", "0".repeat(16), &stamp[16..]);
            replace(&manifest, &stamp, &edited);
            replace(&pack.join(METADATA_FILE), &stamp, &edited);
        },
        true,
        "another cache version",
    );
    // An older step-semantics version in the RustRed manifest.
    assert_rejected(
        |_, cache| {
            replace(
                &cache.join("ordered-step-cache.txt"),
                &format!("semantics {}", rustred::reduction::ORDERED_STEP_SEMANTICS),
                "semantics rustred-ordered-steps-2000-01-01",
            )
        },
        true,
        "this build reads",
    );
    // Other rule tables.
    assert_rejected(
        |pack, _| {
            let path = pack.join(METADATA_FILE);
            let mut metadata: serde_json::Value =
                serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
            metadata["programs"][0]["rule_tables"] = "f".repeat(64).into();
            fs::write(&path, serde_json::to_vec_pretty(&metadata).unwrap()).unwrap();
        },
        true,
        "other rule tables",
    );
    // A corrupt shard: its sha256 no longer matches the loop manifest.
    assert_rejected(
        |_, cache| {
            let shard = walkdir(&cache.join("shards")).pop().unwrap();
            let mut bytes = fs::read(&shard).unwrap();
            let last = bytes.len() - 1;
            bytes[last] ^= 1;
            fs::write(&shard, bytes).unwrap();
        },
        false,
        "SHA256 mismatch",
    );
    // A shard file the RustRed manifest does not list.
    assert_rejected(
        |_, cache| fs::write(cache.join("shards/extra.rrstep"), b"extra").unwrap(),
        true,
        "differ from the shards",
    );
}

#[test]
fn tail_table_headers_are_checked_by_version() {
    let identity = format!(
        "schema {}\nsemantics {}\nadapter SavedApplication\n",
        rustred::reduction::ORDERED_TAIL_SCHEMA,
        rustred::reduction::ORDERED_TAIL_SEMANTICS
    );
    let file = |header: &str| {
        let mut bytes = b"RR-ORDERED-TAIL2".to_vec();
        bytes.extend((header.len() as u64).to_le_bytes());
        bytes.extend(header.as_bytes());
        bytes.extend([0; 40]);
        bytes
    };
    check_tail_header(&file(&identity)).unwrap();
    let older = identity.replace(
        rustred::reduction::ORDERED_TAIL_SCHEMA,
        "rustred.ordered-tail.v2",
    );
    assert!(
        check_tail_header(&file(&older))
            .unwrap_err()
            .contains("this build reads")
    );
    assert!(check_tail_header(&file(&identity)[..30]).is_err());
}

#[test]
fn rustred_manifests_parse_their_stamp_and_shards() {
    let manifest = parse_manifest(
        "schema s\nsemantics t\nadapter closed-artifact\nfamily f%20g\nvariables 1 d\n\
         unresolved none\nstamp abc\nshard integral 3 0000000000000000 shards/a.rrstep 5 17 \
         00\nend\n",
    )
    .unwrap();
    assert_eq!(
        (
            manifest.schema.as_str(),
            manifest.semantics.as_str(),
            manifest.adapter.as_str(),
            manifest.stamp.as_str()
        ),
        ("s", "t", "closed-artifact", "abc")
    );
    assert_eq!(manifest.shards, [("shards/a.rrstep".to_owned(), 5, 17)]);
    assert!(parse_manifest("schema s\nstamp abc\n").is_err());
}
