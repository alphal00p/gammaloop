//! The adapter `evaluate_factorized` (factorised value -> idenso contraction x
//! prefactor -> vakint's existing p evaluation) against the reference path
//! (lane refpath, `vakint::refpath::run_graph`) on every background-field
//! `B -> B` graph of one and two loops of the beta thread's prepared DOT files,
//! every colour branch through `FactorizedNumerator::from_feynkit`, three
//! projections: (1) the value, exactly, against refpath run live in the same
//! process; (2) the value against lane refpath's ARCHIVED values
//! (B/num-refpath/reference, `<set>.values.tar.zst`) where an archive exists;
//! (3) what RustRed receives (`p_gate` views, zero-sector keys dropped) summed
//! over the branches against refpath's view of the summed numerator.
//!
//! Data: the DOT root (`VAKINT_TEST_FACTORIZED_DOT_ROOT`, default
//! `/common/ben/beta/results/prepared`), the archive root
//! (`VAKINT_TEST_REFPATH_REFERENCE`, default lane refpath's), the p library of
//! the archive (`VAKINT_TEST_P_LIBRARY_1_4`, else the built-in one). Without
//! the DOT data the test prints SKIP and passes; without an archive it checks
//! the live comparison only.
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use feynkit_graph::FeynmanDiagram;
use feynkit_model::Model;
use symbolica::atom::{Atom, AtomCore};
use vakint::Vakint;
use vakint::p_integral::FactorizedNumerator;
use vakint::p_integral::factorized::{
    FromFeynkitOptions, Projector, ProjectorChoice, evaluate_factorized,
};
use vakint::refpath::{self, Projection, RunOptions};
use vakint::rustred_evaluation::p_gate::{self, GateView};

const DEFAULT_ROOT: &str = "/common/ben/beta/results/prepared";
const DEFAULT_REFERENCE: &str =
    "/common/build-cache/vhirschi/p_integral_helper_VH/num-refpath/reference";

fn root() -> Option<PathBuf> {
    let root = std::env::var("VAKINT_TEST_FACTORIZED_DOT_ROOT")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from(DEFAULT_ROOT));
    if root.join("L1").join("model.json").is_file() {
        Some(root)
    } else {
        eprintln!("SKIP: no DOT data under {}", root.display());
        None
    }
}

fn dots(dir: &Path) -> Vec<PathBuf> {
    let mut files = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with("diagram-") && n.ends_with(".dot"))
        })
        .collect::<Vec<_>>();
    files.sort();
    files
}

/// The archived values of one set (`graph number -> projection -> atom`), or
/// `None` when the archive is not readable.
fn archived(set: &str) -> Option<BTreeMap<String, BTreeMap<String, Atom>>> {
    let reference = std::env::var("VAKINT_TEST_REFPATH_REFERENCE")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from(DEFAULT_REFERENCE));
    let archive = reference.join(format!("{set}.values.tar.zst"));
    if !archive.is_file() {
        eprintln!("NOTE: no archive {}", archive.display());
        return None;
    }
    let directory =
        std::env::temp_dir().join(format!("factorized-evaluate-{set}-{}", std::process::id()));
    std::fs::create_dir_all(&directory).ok()?;
    let status = std::process::Command::new("bash")
        .arg("-c")
        .arg(format!(
            "zstd -qdc '{}' | tar -x -C '{}'",
            archive.display(),
            directory.display()
        ))
        .status()
        .ok()?;
    if !status.success() {
        return None;
    }
    let mut out = BTreeMap::new();
    for entry in std::fs::read_dir(directory.join("values")).ok()? {
        let path = entry.ok()?.path();
        let graph = path.file_stem()?.to_str()?.to_string();
        let mut values = BTreeMap::new();
        for line in std::fs::read_to_string(&path).ok()?.lines() {
            let Some((key, atom)) = line.split_once('\t') else {
                continue;
            };
            if key.contains('.') {
                continue; // diagnostics (`numerator.metric.bytes`)
            }
            Vakint::initialize_vakint_symbols();
            let atom = symbolica::try_parse!(atom, default_namespace = "vakint").unwrap();
            values.insert(key.to_string(), atom);
        }
        out.insert(graph, values);
    }
    let _ = std::fs::remove_dir_all(&directory);
    Some(out)
}

fn engine() -> Vakint {
    let mut vakint = Vakint::new().unwrap();
    if let Some(library) = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4") {
        vakint
            .load_p_library_directory(Path::new(&library), false)
            .unwrap();
    }
    vakint
}

/// The refpath convention of a value's prefactor: couplings expanded by the
/// model at `gs = 1`, `refpath::nf`, colour invariants in refpath's symbols.
fn refpath_prefactor(model: &Model, value: &FactorizedNumerator) -> Atom {
    use idenso::color::CS;
    let total = model
        .expand_couplings(&value.prefactor.total())
        .replace(Atom::var(symbolica::symbol!("UFO::gs")))
        .with(Atom::num(1))
        .replace(Atom::var(symbolica::symbol!("vakint::nf")))
        .with(Atom::var(refpath::SYMBOLS.nf));
    refpath::map_invariants(
        &total,
        &Atom::var(symbolica::symbol!("spenso::dA")),
        &Atom::var(CS.nc),
    )
    .expand()
    .together()
}

/// The value with refpath's projection (`p2_power`, tensor) and the extra
/// normalisation scalar (`1/(i (D-1))` for metric / transverse, `1/i` for the
/// longitudinal projection).
fn projected(value: &FactorizedNumerator, projection: Projection) -> (FactorizedNumerator, Atom) {
    let mut value = value.clone();
    let d = Atom::var(value.dimension);
    let i = Atom::i();
    match projection {
        Projection::Metric => {
            value.projector.as_mut().unwrap().p2_power = -1;
            (value, Atom::num(1) / (i * (&d - Atom::num(1))))
        }
        Projection::Transverse => {
            value.projector.as_mut().unwrap().p2_power = -2;
            (value, Atom::num(1) / (i * (&d - Atom::num(1))))
        }
        Projection::Longitudinal => {
            let x = value.external_edges[0].unwrap() as i64;
            let q = feynkit_graph::momentum_symbol();
            let tensor = symbolica::function!(q, x, &value.ports[0])
                * symbolica::function!(q, x, &value.ports[1]);
            value.projector = Some(Projector {
                tensor,
                p2_power: -2,
            });
            (value, Atom::num(1) / i)
        }
        Projection::Both => unreachable!(),
    }
}

fn add_views(total: &mut GateView, view: GateView) {
    for (key, value) in view {
        let entry = total.entry(key).or_insert(Atom::Zero);
        *entry = (&*entry + value).expand();
    }
}

/// One graph and projection: the adapter's value (sum over the colour
/// branches) and refpath's, compared exactly, and their gate views; returns
/// (graph number, adapter value, refpath loops).
fn compare_graph(
    vakint: &Vakint,
    model: &Arc<Model>,
    path: &Path,
    projection: Projection,
    choice: ProjectorChoice,
) -> (String, Atom, usize) {
    let graph = path
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .trim_start_matches("diagram-")
        .to_string();
    let text = std::fs::read_to_string(path).unwrap();
    let options = RunOptions {
        projection,
        extra_terms: 1,
        count_terms: false,
        contract_only: false,
    };
    let reference = refpath::run_graph(vakint, model, &text, &options).unwrap();
    let settings = refpath::vakint_settings(reference.loops as i64 + 1);
    let (expected, expected_view) = match reference.results.first() {
        Some(result) => (result.value.clone(), refpath::gate_view(result).0),
        None => (Atom::Zero, GateView::new()),
    };
    let diagram = FeynmanDiagram::from_dot(Arc::clone(model), &text).unwrap();
    let values = FactorizedNumerator::from_feynkit(
        &diagram,
        model,
        &FromFeynkitOptions {
            projector: choice,
            ..FromFeynkitOptions::default()
        },
    )
    .unwrap();
    let mut total = Atom::Zero;
    let mut view = GateView::new();
    for value in &values {
        let (value, normalisation) = projected(value, projection);
        let prefactor = refpath_prefactor(model, &value) * normalisation;
        let (result, observation) =
            p_gate::observe(|| evaluate_factorized(vakint, &settings, &value, &prefactor));
        total += result.unwrap().value;
        add_views(&mut view, observation.view_without_zero_sectors().0);
    }
    view.retain(|_, value| !value.is_zero());
    let total = total.expand();
    assert!(
        (&total - &expected).expand().is_zero(),
        "{graph} {}: adapter {total} vs refpath {expected}",
        projection.name()
    );
    let difference = p_gate::compare_gate_views(&expected_view, &view);
    assert!(
        difference.is_empty(),
        "{graph} {}: gate views differ: {difference:?}",
        projection.name()
    );
    (graph, total, reference.loops)
}

const PROJECTIONS: [(Projection, ProjectorChoice); 3] = [
    (Projection::Metric, ProjectorChoice::Metric),
    (Projection::Transverse, ProjectorChoice::Transverse),
    (Projection::Longitudinal, ProjectorChoice::Metric),
];

#[test]
fn the_adapter_equals_refpath_on_every_l1_and_l2_graph() {
    let Some(root) = root() else { return };
    let vakint = engine();
    let mut checked = 0;
    let mut archived_checked = 0;
    for level in ["L1", "L2"] {
        let directory = root.join(level);
        let model = refpath::load_model(&directory).unwrap();
        for (projection, choice) in PROJECTIONS {
            let set = format!("{level}-{}", projection.name());
            let archive = archived(&set);
            for path in dots(&directory) {
                let (graph, total, _) = compare_graph(&vakint, &model, &path, projection, choice);
                if let Some(archive) = &archive {
                    let stored = archive
                        .get(&graph)
                        .and_then(|values| values.get(projection.name()))
                        .cloned()
                        .unwrap_or(Atom::Zero);
                    assert!(
                        (&total - &stored).expand().is_zero(),
                        "{level} {graph} {}: adapter {total} vs archived {stored}",
                        projection.name()
                    );
                    archived_checked += 1;
                }
                checked += 1;
            }
            eprintln!(
                "factorized_evaluate {set}: {} graphs equal to refpath, archive {}",
                dots(&directory).len(),
                if archive.is_some() {
                    "compared"
                } else {
                    "absent"
                }
            );
        }
    }
    eprintln!(
        "factorized_evaluate: {checked} graph projections equal to refpath, {archived_checked} also to the archive"
    );
    assert_eq!(checked, 3 * (3 + 29));
}

/// The bare part of the D5 acceptance fixture (N/notes/D5-RSTAR-INTERFACE.md
/// section 11): L3 diagram 185 (a one-loop quark bubble on a quantum gluon
/// line and a closed ghost loop: odd loop count, colour I2R cA^2 nf, overall
/// sign +1 = (-1)_quark (-1)_ghost) through the adapter equals refpath in
/// all three projections; its weights are counted on the ORIGINAL graph.
#[test]
fn the_acceptance_fixture_bare_graph_equals_refpath() {
    let Some(root) = root() else { return };
    let directory = root.join("L3");
    let path = directory.join("diagram-000185.dot");
    if !path.is_file() {
        eprintln!("SKIP: no {}", path.display());
        return;
    }
    let vakint = engine();
    let model = refpath::load_model(&directory).unwrap();
    let diagram =
        FeynmanDiagram::from_dot(Arc::clone(&model), &std::fs::read_to_string(&path).unwrap())
            .unwrap();
    let values =
        FactorizedNumerator::from_feynkit(&diagram, &model, &FromFeynkitOptions::default())
            .unwrap();
    for value in &values {
        let loops = value.closed_loops().unwrap();
        assert_eq!(loops.len(), 2, "a quark and a ghost loop");
        assert_eq!(value.prefactor.loop_sign, 1);
        assert_eq!(value.prefactor.nf_power, 1);
    }
    for (projection, choice) in PROJECTIONS {
        let (_, total, loops) = compare_graph(&vakint, &model, &path, projection, choice);
        assert_eq!(loops, 3);
        eprintln!("FIXTURE 185 {}: {total}", projection.name());
    }
}
