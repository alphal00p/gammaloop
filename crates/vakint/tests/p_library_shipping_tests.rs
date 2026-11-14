//! The built-in 1L-3L p library (M2-V4; `docs/p_integrals.md`, section 6,
//! "Built-in library"): the shipped `data/p/` against the current RustRed,
//! the default path against Forcer's frozen values, fail-closed data and the
//! step cache of the built-in program. The crate-internal checks (embedded
//! files against the manifest, vakint's families against RustRed's and the
//! library's, the listed-library resolver) are unit tests of
//! `rustred_evaluation::p_data`. The checks of a library take its scope from
//! the library directory (M3-V8); the scope 1-4 package of M3 runs them in
//! the ignored tests below (`VAKINT_TEST_P_LIBRARY_1_4`).
//!
//! `VAKINT_P_DATA=write` rewrites `data/p/manifest.json` from the shipped
//! library directory (after `rustred p-package`, see the docs).

#![cfg(feature = "feynkit-ingress")]

#[path = "support/p_backends.rs"]
mod p_backends;
#[path = "support/p_integrals.rs"]
mod p_integrals;
#[path = "support/p_library4.rs"]
mod p_library4;
#[path = "support/p_rustred.rs"]
mod p_rustred;

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use p_backends::values_agree;
use p_library4::{
    P_LIBRARY_1_4_VARIABLE, PLibraryRoute, built_in_p_library, engine_with_p_library,
    load_p_library_directory, p_library_1_4, p_library_name, p_library_scope, regenerate_p_library,
    require_p_library_1_4,
};
use rustred::propagator::families::{P_FAMILIES_V1_JSON, PFamiliesFile};
use rustred::propagator::library::{PBasis, PLibrary};
use rustred::propagator::present::top_topology;
use rustred::propagator::tag::GammaTag;
use rustred::propagator::topology::{ChildMode, reduction_graph, topologies_json};
use rustred::propagator::{PTerminalId, p_versions};
use rustred_app::propagator::loader::read_custom_payloads;
use rustred_app::propagator::package::assemble_p_package;
use rustred_app::propagator::{
    AssembleOptions, CustomOwnerOptions, CustomPayloads, PLibrarySource,
};
use sha2::{Digest, Sha256};
use symbolica::atom::{Atom, AtomCore};
use vakint::p_integral::forcer_fixtures::{ForcerCase, ForcerFixture};
use vakint::{
    EvaluationMethod, OrderedCacheProgramStatistics, OrderedCacheSettings, RustRedMasterMode,
    Vakint, VakintSettings,
};

fn data() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data")
}

/// The built-in library directory, `data/p/p-library-<a>-<b>`.
fn shipped_directory() -> PathBuf {
    built_in_p_library()
}

fn shipped_source() -> PLibrarySource {
    PLibrarySource::Directory(shipped_directory())
}

fn shipped() -> Arc<PLibrary> {
    load_p_library_directory(&shipped_directory())
}

fn json(path: &Path) -> serde_json::Value {
    serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
}

/// Every file below `directory`, by `/`-separated relative path.
fn files(directory: &Path) -> BTreeMap<String, Vec<u8>> {
    fn visit(root: &Path, directory: &Path, files: &mut BTreeMap<String, Vec<u8>>) {
        for entry in fs::read_dir(directory).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                visit(root, &path, files);
            } else {
                let name = path.strip_prefix(root).unwrap().to_str().unwrap();
                files.insert(name.replace('\\', "/"), fs::read(&path).unwrap());
            }
        }
    }
    let mut files = BTreeMap::new();
    visit(directory, directory, &mut files);
    files
}

/// A writable external-data root: vakint's `data/` without (`with_p` false)
/// or with its `p/` area.
fn external_data_root(with_p: bool) -> tempfile::TempDir {
    let root = tempfile::tempdir().unwrap();
    for (name, bytes) in files(&data()) {
        if with_p || !name.starts_with("p/") {
            let path = root.path().join(&name);
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            fs::write(path, bytes).unwrap();
        }
    }
    root
}

fn battery(loops: usize) -> ForcerFixture {
    let text = match loops {
        2 => include_str!("fixtures/forcer/p2l_battery.json"),
        3 => include_str!("fixtures/forcer/p3l_battery.json"),
        _ => unreachable!(),
    };
    ForcerFixture::parse(text).unwrap()
}

fn case(loops: usize, name: &str) -> ForcerCase {
    let mut battery = battery(loops);
    let index = battery
        .cases
        .iter()
        .position(|case| case.name == name)
        .unwrap();
    battery.cases.swap_remove(index)
}

/// A dotted 2L kite, a 3L input whose value carries `PMaster(no)` and one
/// carrying `PMaster(t105)` (the M2-V2 batteries).
fn default_path_cases() -> [ForcerCase; 3] {
    [case(2, "kite_dot1"), case(3, "no_top"), case(3, "be_top")]
}

/// The p program of `engine`: the one cache program `p-library-<a>-<b>`.
fn p_program(engine: &Vakint) -> OrderedCacheProgramStatistics {
    let mut programs = (engine.ordered_cache_statistics().programs.into_iter())
        .filter(|program| program.program.starts_with("p-library-"))
        .collect::<Vec<_>>();
    assert_eq!(programs.len(), 1, "{programs:?}");
    programs.pop().unwrap()
}

/// The record, hit, miss and insert counts and the file bytes of a program.
fn counts(program: &OrderedCacheProgramStatistics) -> String {
    let native = &program.native;
    format!(
        "records {}, hits {}, misses {}, inserts {}, file bytes {}",
        native.records, native.hits, native.misses, native.inserts, native.file_bytes
    )
}

/// `settings` in `mode` for both p backends, through `ε^(L+2)` for `loops`.
fn in_mode(settings: &VakintSettings, mode: RustRedMasterMode, loops: usize) -> VakintSettings {
    let mut settings = VakintSettings {
        number_of_terms_in_epsilon_expansion: 2 * loops as i64 + 3,
        run_time_decimal_precision: 40,
        ..settings.clone()
    };
    for method in settings.evaluation_order.0.iter_mut() {
        match method {
            EvaluationMethod::RustRed(options) => options.master_mode = mode,
            EvaluationMethod::Forcer(options) => options.master_mode = mode,
            _ => {}
        }
    }
    settings
}

fn evaluate(engine: &Vakint, settings: &VakintSettings, case: &ForcerCase) -> Atom {
    engine
        .evaluate(settings, case.expression().unwrap().as_view())
        .unwrap_or_else(|error| panic!("{}: {error}", case.name))
}

/// `data/p/manifest.json` describes the shipped library directory: its name,
/// the identity and versions of its `library.json`, every file with its
/// sha256 and size. `VAKINT_P_DATA=write` rewrites it.
#[test]
fn the_p_data_manifest_lists_the_shipped_library_files() {
    let library = json(&shipped_directory().join("library.json"));
    let expected = serde_json::json!({
        "schema": "vakint.p-data.v1",
        "library": p_library_name(&shipped_directory()),
        "identity": library["identity"],
        "versions": library["versions"],
        "files": files(&shipped_directory()).iter().map(|(path, bytes)| serde_json::json!({
            "path": path,
            "sha256": format!("{:x}", Sha256::digest(bytes)),
            "bytes": bytes.len(),
        })).collect::<Vec<_>>(),
    });
    let manifest = data().join("p/manifest.json");
    if std::env::var("VAKINT_P_DATA").as_deref() == Ok("write") {
        let text = serde_json::to_string_pretty(&expected).unwrap() + "\n";
        fs::write(&manifest, text).unwrap();
    }
    assert_eq!(json(&manifest), expected);
}

/// The shipped identity and versions are this RustRed's (`p_versions()`),
/// in the manifest, in `library.json` and in the loaded library.
#[test]
fn p_library_versions_match_rustred() {
    let versions = p_versions();
    let [start, end] = p_library_scope(&shipped_directory());
    let identity = format!(
        "p-library:scope={start}-{end}:{}",
        (versions.iter())
            .map(|(name, value)| format!("{name}={value}"))
            .collect::<Vec<_>>()
            .join(",")
    );
    // A shipped rule ordering (p artefact v2, `ordering.json`) is part of the
    // identity, after the versions: `,P_ORDERING_VERSION=<version>:<heuristic>`.
    let identity = match shipped().ordering() {
        Some(tag) => format!("{identity},P_ORDERING_VERSION={tag}"),
        None => identity,
    };
    // Batch 6: a package content revision ends the identity, after the ordering.
    let identity = match shipped().revision() {
        Some(name) => format!("{identity},P_PACKAGE_REVISION={name}"),
        None => identity,
    };
    let versions: serde_json::Map<_, _> = (versions.into_iter())
        .map(|(name, value)| (name.to_owned(), value.into()))
        .collect();
    for file in [
        data().join("p/manifest.json"),
        shipped_directory().join("library.json"),
    ] {
        let content = json(&file);
        assert_eq!(content["identity"], identity, "{}", file.display());
        assert_eq!(
            content["versions"],
            serde_json::Value::Object(versions.clone()),
            "{}",
            file.display()
        );
    }
    assert_eq!(shipped().identity(), identity);
}

/// Every terminal of the shipped library maps onto `1` or a Forcer master
/// through its basis layer, never a residual: the 1-3 library onto `1`,
/// `PMaster(no)` and `PMaster(t105)`, each of them reached. vakint's own
/// load check (masters it tabulates) accepts it.
#[test]
fn every_terminal_maps_onto_one_or_a_forcer_master() {
    let library = shipped();
    engine_with_p_library(&shipped_directory(), PLibraryRoute::Injected);
    let map = library
        .basis_map()
        .expect("the shipped library has a basis layer");
    let mut bases = BTreeSet::new();
    for id in 0..library.terminals().len() {
        let row = map.map(PTerminalId(id as u32), &GammaTag::one()).unwrap();
        assert!(!row.is_empty(), "terminal {id}");
        for term in row {
            bases.insert(match term.basis {
                PBasis::Unit => "1".to_owned(),
                PBasis::Master(name) => name.to_string(),
                PBasis::Residual(residual) => panic!("terminal {id} -> residual {}", residual.0),
            });
        }
    }
    if p_library_scope(&shipped_directory()) == [1, 3] {
        assert_eq!(bases, BTreeSet::from(["1", "no", "t105"].map(String::from)));
    }
}

/// The shipped library is the program the current RustRed assembles
/// (`assemble_p_package`: assembly and basis layer), compared structurally
/// (`equivalent_p_libraries`), never by bytes: the regeneration check
/// (`regenerate_p_library`: the descent certificate once, the stored LP of
/// its Custom owners NO and t1star05 counted again, the assembly from its
/// payloads), and with NO and t1star05 generated again.
#[test]
fn embedded_p_library_matches_current_rustred_assembly() {
    let regeneration = regenerate_p_library(&shipped_directory(), 1);
    println!("{regeneration:?}");
    assert_eq!(
        regeneration.owners.len(),
        2,
        "the Custom owners NO and t1star05"
    );
    assert!(
        regeneration
            .owners
            .iter()
            .all(|owner| owner.recounted.is_some())
    );
    let families: PFamiliesFile = serde_json::from_str(P_FAMILIES_V1_JSON).unwrap();
    let options = AssembleOptions {
        threads: 1,
        verify_descent: false,
    };
    let generated = CustomPayloads::Generate(CustomOwnerOptions::default());
    let [_, end] = p_library_scope(&shipped_directory());
    let assembled = assemble_p_package(&families, end, generated, &options).unwrap();
    let shipped = shipped();
    p_library4::assert_equivalent_up_to_revision(&shipped, &assembled.library);
    assert_eq!(read_custom_payloads(&shipped_source()).unwrap().len(), 2);
}

/// The regeneration check of M3's scope 1-4 package (`VAKINT_TEST_P_LIBRARY_1_4`,
/// from the PP-M3b pin on): the descent certificate once, the LP of its
/// Custom owners of at most 9 lines counted again (every owner with
/// `P_LP_RECOUNT_ALL=1`), the assembly from its payloads equal to it.
#[test]
#[ignore = "slow: the scope 1-4 certificate and assembly (about 10 min)"]
fn the_one_to_four_library_regenerates_from_its_payloads() {
    let directory = require_p_library_1_4("the_one_to_four_library_regenerates_from_its_payloads");
    let regeneration = regenerate_p_library(&directory, 4);
    println!("{regeneration:#?}");
    assert_eq!(regeneration.scope, [1, 4]);
}

/// Graphica skew: `rustred_topology` run in vakint's build (gammaloop's
/// graphica) on the shipped tops reproduces the shipped DAG node by node,
/// keys and actions included.
#[test]
fn p_topology_keys_match_shipped_library() {
    let families: PFamiliesFile = serde_json::from_slice(
        &fs::read(shipped_directory().join("families/p-families.json")).unwrap(),
    )
    .unwrap();
    let roots = families
        .families
        .iter()
        .map(top_topology)
        .collect::<Vec<_>>();
    let library = shipped();
    let dag = library.topology().expect("the shipped library has its DAG");
    let counts = json(&shipped_directory().join("library.json"))["counts"].clone();
    assert_eq!(serde_json::json!(dag.nodes.len()), counts["dag_nodes"]);
    assert_eq!(
        topologies_json(&reduction_graph(&roots, ChildMode::ForcerTree)),
        topologies_json(dag)
    );
}

/// A library directory under an external-data root whose `library.json`
/// declares another version fails closed with the regeneration hint: its
/// p-integrals neither use the built-in library nor go to Forcer.
#[test]
fn a_library_of_another_version_fails_closed() {
    let root = external_data_root(true);
    let library = root
        .path()
        .join("p")
        .join(p_library_name(&shipped_directory()))
        .join("library.json");
    let text = fs::read_to_string(&library).unwrap();
    let tampered = text.replacen(
        "\"P_CLASS_KEY_VERSION\": \"1\"",
        "\"P_CLASS_KEY_VERSION\": \"2\"",
        1,
    );
    assert_ne!(tampered, text);
    fs::write(&library, tampered).unwrap();
    let mut engine = Vakint::new().unwrap();
    engine.load_external_data(root.path(), 3).unwrap();
    let kite = case(2, "kite");
    let error = engine
        .evaluate(
            &VakintSettings::default(),
            kite.expression().unwrap().as_view(),
        )
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("P_CLASS_KEY_VERSION")
            && error.contains("regenerate with rustred p-package"),
        "{error}"
    );
}

/// D-V3-5 (review R2-S2): a `p/` area that cannot be resolved counts as data
/// for every loop count, so a 1L and a 4L p-integral both fail with its error
/// under the default order instead of going to the built-in library or to
/// Forcer.
#[test]
fn a_broken_p_area_fails_closed_for_every_loop_count() {
    let root = external_data_root(false);
    fs::create_dir(root.path().join("p")).unwrap();
    let mut engine = Vakint::new().unwrap();
    engine.load_external_data(root.path(), 3).unwrap();
    for input in [
        "topo(P1L(p(1),1,1))",
        "topo(P4L_LA4(p(1),1,1,1,1,1,1,1,1,1,1,1))",
    ] {
        let error = engine
            .evaluate(
                &VakintSettings::default(),
                p_integrals::parse(input).as_view(),
            )
            .unwrap_err()
            .to_string();
        assert!(error.contains("holds no p-library-"), "{input}: {error}");
    }
}

/// The default path without external data: `Vakint::new()` and
/// `VakintSettings::default()` (default order, normalization and arithmetic)
/// reduce a dotted 2L kite and 3L inputs carrying `PMaster(no)` and
/// `PMaster(t105)` with the built-in 1L-3L library, equal to Forcer's frozen
/// values finalized alike: exactly in `Symbolic` mode through `ε^(L+2)`,
/// within 1e-30 in the default `Numerical` mode. No FORM is called. Four
/// loops need the scope 1-4 library of vakint-external-data (next test);
/// without it they stay with Forcer.
#[test]
fn the_default_path_matches_forcer_through_three_loops() {
    let engine = Vakint::new().unwrap();
    let default = VakintSettings::default();
    assert!(matches!(
        default.evaluation_order.0[0],
        EvaluationMethod::RustRed(_)
    ));
    for case in default_path_cases() {
        let loops = case.loops().unwrap();
        for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
            let settings = in_mode(&default, mode, loops);
            let actual = evaluate(&engine, &settings, &case);
            let expected = case.finalized(&settings, mode).unwrap();
            values_agree(&actual, &expected, mode)
                .unwrap_or_else(|difference| panic!("{} ({mode}): {difference}", case.name));
        }
    }
    let program = p_program(&engine);
    assert_eq!(program.data_version, "built-in");
    assert!(program.native.inserts > 0);
}

/// The default path at four loops (plan M3 D-M3-7, the default switch):
/// `Vakint::new()`, `load_external_data` of a root whose `p_integral/` area
/// lists M3's scope 1-4 library as vakint-external-data's `p_integral` branch
/// ships it (`VAKINT_TEST_P_LIBRARY_1_4`; a skip on stderr without it) and
/// `VakintSettings::default()`: the four-loop `bubu` top, whose value carries
/// the masters `bebe`, `cross` and `t145`, is reduced by RustRed with that
/// library and equals Forcer's frozen value finalized alike, exactly in
/// `Symbolic` mode and within 1e-30 in `Numerical` mode, through `ε^(L+2)`
/// or the deepest order its masters reach.
#[test]
fn the_default_path_with_external_data_matches_forcer_at_four_loops() {
    let Some(package) = p_library_1_4() else {
        eprintln!("{P_LIBRARY_1_4_VARIABLE} is unset: skipped the four-loop default path");
        return;
    };
    let engine = engine_with_p_library(&package, PLibraryRoute::Listed);
    let case = ForcerFixture::parse(include_str!("fixtures/forcer/p4l_anchors.json"))
        .unwrap()
        .cases
        .into_iter()
        .find(|case| case.name == "bubu_top")
        .unwrap();
    assert_eq!(case.loops().unwrap(), 4);
    let default = VakintSettings::default();
    let terms = case.deepest().unwrap().terms.min(2 * 4 + 3);
    assert!(terms > 4, "{} does not finalize through ε^0", case.name);
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        let settings = VakintSettings {
            number_of_terms_in_epsilon_expansion: terms,
            ..in_mode(&default, mode, 4)
        };
        let actual = evaluate(&engine, &settings, &case);
        let expected = case.finalized(&settings, mode).unwrap();
        values_agree(&actual, &expected, mode)
            .unwrap_or_else(|difference| panic!("{} ({mode}): {difference}", case.name));
    }
    let program = p_program(&engine);
    assert_eq!(
        (program.program.as_str(), program.data_version.as_str()),
        ("p-library-1-4", "external-data")
    );
    assert!(program.native.inserts > 0, "{program:?}");
}

/// A p-integral whose numerator keeps a tensor spectator with native `dim`
/// metric slots (as gammaloop's R* hands it over, `uv::rstar`) evaluates alike
/// with RustRed (the built-in library) and with Forcer: both p backends take
/// one numerator path. Needs FORM 5 for the Forcer side.
#[test]
fn a_tensor_spectator_takes_one_numerator_path_with_both_p_backends() {
    let Some(form) = p_backends::live_form("the tensor spectator through Forcer") else {
        return;
    };
    let input = p_integrals::parse(
        "g(spenso::mink(dim,mu),spenso::mink(dim,nu))*dot(k(1),p(1))*topo(P2L(p(1),1,1,1,1,1))",
    );
    let value = |method| {
        let settings = VakintSettings {
            evaluation_order: vakint::EvaluationOrder(vec![method]),
            form_exe_path: form.clone(),
            use_dot_product_notation: true,
            ..in_mode(&VakintSettings::default(), RustRedMasterMode::Symbolic, 2)
        };
        Vakint::new()
            .unwrap()
            .evaluate_integral(&settings, input.as_view())
            .unwrap_or_else(|error| panic!("{error}"))
    };
    let mode = RustRedMasterMode::Symbolic;
    let rustred = value(EvaluationMethod::RustRed(
        vakint::RustRedEvaluationOptions { master_mode: mode },
    ));
    let forcer = value(EvaluationMethod::Forcer(vakint::ForcerOptions {
        master_mode: mode,
    }));
    assert!(!rustred.is_zero());
    assert!(
        (&rustred - &forcer).expand().is_zero(),
        "{rustred}\n{forcer}"
    );
}

/// The built-in program's step cache on disk: directory `p-library-1-3_...`
/// (no `_injected`, at most 200 bytes); a read-write evaluation of a 2L and a
/// 3L p-integral writes it, a fresh read-only engine replays both without a
/// miss or an insert, and so do fresh engines whose library of the same
/// identity comes from an external-data `p/` area or is listed in a
/// `p_integral/` area (M3-V8; same directory).
#[test]
fn the_built_in_step_cache_is_replayed_by_fresh_engines() {
    let directory = tempfile::tempdir().unwrap();
    let inputs = [case(2, "kite_dot1"), case(3, "no_top")];
    let symbolic = |cache| VakintSettings {
        ordered_cache: cache,
        ..in_mode(&VakintSettings::default(), RustRedMasterMode::Symbolic, 3)
    };
    let (read_write, read_only) = (
        symbolic(OrderedCacheSettings::read_write(directory.path())),
        symbolic(OrderedCacheSettings::read_only(directory.path())),
    );
    let run = |engine: &Vakint, settings: &VakintSettings| {
        let values = inputs
            .iter()
            .map(|case| evaluate(engine, settings, case))
            .collect::<Vec<_>>();
        (values, p_program(engine))
    };

    let (values, cold) = run(&Vakint::new().unwrap(), &read_write);
    let cache = cold.directory.clone().unwrap();
    let name = cache.file_name().unwrap().to_str().unwrap();
    assert!(
        name.starts_with(&format!("{}_", p_library_name(&shipped_directory())))
            && !name.ends_with("_injected"),
        "{name}"
    );
    assert!(name.len() <= 200, "{} bytes", name.len());
    assert_eq!(cache.parent(), Some(directory.path()));
    assert!(cold.native.inserts > 0 && cache.is_dir(), "{cold:?}");
    println!("cold read-write: {}", counts(&cold));

    let root = external_data_root(true);
    let mut external = Vakint::new().unwrap();
    external.load_external_data(root.path(), 3).unwrap();
    let listed = engine_with_p_library(&shipped_directory(), PLibraryRoute::Listed);
    for (engine, origin) in [
        (Vakint::new().unwrap(), "built-in"),
        (external, "external-data"),
        (listed.engine, "external-data"),
    ] {
        let (replayed, warm) = run(&engine, &read_only);
        assert_eq!(replayed, values, "{origin}");
        assert_eq!(
            (warm.data_version.as_str(), &warm.directory),
            (origin, &cold.directory)
        );
        assert_eq!(
            (warm.native.misses, warm.native.inserts),
            (0, 0),
            "{origin}: {warm:?}"
        );
        println!("warm read-only ({origin}): {}", counts(&warm));
    }
}

/// Grant U2 (review R2-S3): the p program the engine resolves to keeps its
/// in-memory records across an explicit external-data (re-)activation, for
/// the built-in, an external-data and an injected library; the next
/// evaluation replays them.
#[test]
fn the_resolved_p_program_keeps_its_records_across_external_data_activation() {
    let (without_p, with_p) = (external_data_root(false), external_data_root(true));
    let kite = case(2, "kite_dot1");
    let bubble = ForcerFixture::parse(include_str!("fixtures/forcer/p1l_battery.json"))
        .unwrap()
        .cases
        .swap_remove(0);
    let settings = in_mode(&VakintSettings::default(), RustRedMasterMode::Symbolic, 2);

    let mut built_in = Vakint::new().unwrap();
    let mut external = Vakint::new().unwrap();
    external.load_external_data(with_p.path(), 3).unwrap();
    let mut injected = Vakint::new().unwrap();
    injected
        .inject_p_library(rustred::propagator::testing::one_loop_reference_library().unwrap())
        .unwrap();
    for (engine, root, input, origin) in [
        (&mut built_in, without_p.path(), &kite, "built-in"),
        (&mut external, with_p.path(), &kite, "external-data"),
        (&mut injected, without_p.path(), &bubble, "injected"),
    ] {
        let value = evaluate(engine, &settings, input);
        let before = p_program(engine);
        assert!(before.native.records > 0, "{origin}: {before:?}");
        engine.load_external_data(root, 3).unwrap();
        let after = p_program(engine);
        assert_eq!(after.data_version, origin);
        assert_eq!(after.native.records, before.native.records, "{origin}");
        assert_eq!(evaluate(engine, &settings, input), value, "{origin}");
        assert_eq!(
            p_program(engine).native.misses,
            before.native.misses,
            "{origin}"
        );
    }
}
