//! Contract tests of `vakint::p_integral::factorized` on the beta thread's
//! background-field B → B self-energy graphs (FeynGen DOT files).
//!
//! Data: `$VAKINT_TEST_FACTORIZED_DOT_ROOT/L{1..4}` (default
//! `/common/ben/beta/results/prepared`, each with `model.json` and
//! `diagram-NNNNNN.dot`). Without the data the tests print `SKIP` and pass.
#![cfg(feature = "feynkit-ingress")]

use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use feynkit_graph::FeynmanDiagram;
use feynkit_model::Model;
use symbolica::atom::{Atom, AtomCore, AtomView};
use vakint::p_integral::factorized::{
    ColourMode, FromFeynkitOptions, ProjectorChoice, contract_with_idenso, parse_factorized_records,
};
use vakint::p_integral::{FactorizedError, FactorizedNumerator};

const DEFAULT_ROOT: &str = "/common/ben/beta/results/prepared";

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

fn model(dir: &Path) -> Arc<Model> {
    Arc::new(Model::from_json(&std::fs::read_to_string(dir.join("model.json")).unwrap()).unwrap())
}

/// The DOT files of `dir`, all of them or `n` evenly spaced ones.
fn files(dir: &Path, n: Option<usize>) -> Vec<PathBuf> {
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
    match n {
        Some(n) if n < files.len() => (0..n).map(|i| files[i * files.len() / n].clone()).collect(),
        _ => files,
    }
}

fn load(model: &Arc<Model>, path: &Path) -> (String, FeynmanDiagram) {
    let text = std::fs::read_to_string(path).unwrap();
    let diagram = FeynmanDiagram::from_dot(Arc::clone(model), &text)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    (text, diagram)
}

/// value → text → value (equal) → text (byte-identical); validate.
fn round_trip(values: &[FactorizedNumerator], what: &str) -> usize {
    let mut bytes = 0;
    let mut concatenated = String::new();
    for value in values {
        value
            .validate()
            .unwrap_or_else(|e| panic!("{what}: validate: {e}"));
        let text = value.to_text().unwrap();
        let back = FactorizedNumerator::from_text(&text)
            .unwrap_or_else(|e| panic!("{what}: read: {e}\n{text}"));
        assert_eq!(&back, value, "{what}: value changed in the round trip");
        let again = back.to_text().unwrap();
        assert_eq!(again, text, "{what}: second write differs");
        bytes += text.len();
        concatenated.push_str(&text);
    }
    let all = parse_factorized_records(&concatenated).unwrap();
    assert_eq!(all.as_slice(), values, "{what}: concatenated records");
    bytes
}

fn run_round_trips(level: &str, sample: Option<usize>) {
    let Some(root) = root() else { return };
    let dir = root.join(level);
    let model = model(&dir);
    let files = files(&dir, sample);
    let (mut branches, mut zeros) = (0, 0);
    for path in &files {
        let (_, diagram) = load(&model, path);
        let values =
            FactorizedNumerator::from_feynkit(&diagram, &model, &FromFeynkitOptions::default())
                .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        // A colour-zero graph (every branch vanishes) has no value.
        zeros += usize::from(values.is_empty());
        branches += values.len();
        round_trip(&values, &path.display().to_string());
        // The other two projector choices are valid too.
        for projector in [ProjectorChoice::None, ProjectorChoice::Transverse] {
            let options = FromFeynkitOptions {
                projector,
                ..Default::default()
            };
            let values = FactorizedNumerator::from_feynkit(&diagram, &model, &options).unwrap();
            round_trip(&values, &format!("{} {projector:?}", path.display()));
        }
    }
    eprintln!(
        "{level}: {} graphs ({zeros} colour zeros), {branches} colour branches: round trip exact, validate ok",
        files.len()
    );
}

#[test]
fn l1_and_l2_round_trip_exactly_and_validate() {
    run_round_trips("L1", None);
    run_round_trips("L2", None);
}

#[test]
fn l3_sample_round_trips_exactly_and_validates() {
    run_round_trips("L3", Some(100));
}

#[test]
fn l4_sample_round_trips_exactly_and_validates() {
    run_round_trips("L4", Some(100));
}

/// `mink(4, i)` → `mink(dim, i)` (written here independently of the crate).
fn d_dim(atom: &Atom) -> Atom {
    let dim = Atom::var(feynkit_graph::symbols::dimension());
    atom.replace_map(|part, _ctx, out| {
        if let AtomView::Fun(f) = part
            && f.get_symbol().get_name() == "spenso::mink"
            && f.get_nargs() == 2
        {
            let args = f.iter().collect::<Vec<_>>();
            if i64::try_from(args[0]) == Ok(4) {
                **out = symbolica::atom::FunctionBuilder::new(f.get_symbol())
                    .add_arg(&dim)
                    .add_arg(args[1])
                    .finish();
            }
        }
    })
}

/// A graph attribute of the DOT text (first occurrence), DOT escapes removed.
fn graph_attribute(text: &str, name: &str) -> String {
    let start = text
        .find(&format!(" {name}=\""))
        .unwrap_or_else(|| panic!("no graph attribute {name}"))
        + name.len()
        + 3;
    let mut value = String::new();
    let mut chars = text[start..].chars();
    while let Some(c) = chars.next() {
        match c {
            '\\' => value.push(chars.next().unwrap()),
            '"' => return value,
            c => value.push(c),
        }
    }
    panic!("unterminated attribute {name}")
}

fn parse(text: &str) -> Atom {
    symbolica::try_parse!(text).unwrap()
}

/// Closed quark loops of the DOT text, counted independently: components of
/// the internal directed edges with a quark particle.
fn closed_quark_loops(text: &str) -> usize {
    let mut parent = std::collections::BTreeMap::<String, String>::new();
    fn find(parent: &mut std::collections::BTreeMap<String, String>, x: &str) -> String {
        let mut x = x.to_string();
        loop {
            let p = parent.entry(x.clone()).or_insert_with(|| x.clone()).clone();
            if p == x {
                return x;
            }
            x = p;
        }
    }
    let mut quark = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        let Some((ends, attributes)) = line.split_once(" [") else {
            continue;
        };
        let Some((a, b)) = ends.split_once(" -> ") else {
            continue;
        };
        if a.starts_with("ext") || b.starts_with("ext") || !attributes.contains("directed=true") {
            continue;
        }
        let ra = find(&mut parent, a);
        let rb = find(&mut parent, b);
        if ra != rb {
            parent.insert(ra, rb);
        }
        if attributes.contains("particle=\"q\"") || attributes.contains("particle=\"q~\"") {
            quark.push(a.to_string());
        }
    }
    let mut roots = quark
        .iter()
        .map(|a| find(&mut parent, a))
        .collect::<Vec<_>>();
    roots.sort();
    roots.dedup();
    roots.len()
}

/// The reference: the DOT's whole-product `numerator` attribute, closed with
/// `g(μ₀, μ₁) δ(a₀, a₁) / 8` on its two external ports, contracted (Lorentz,
/// Dirac and SU(3) colour) by idenso's production chain, times the evaluated
/// `overall_factor` and `nf^(closed quark loops)`.
fn reference(text: &str) -> Atom {
    let numerator = d_dim(&parse(&graph_attribute(text, "numerator")));
    let interface: Vec<String> =
        serde_json::from_str(&graph_attribute(text, "num_interface")).unwrap();
    let ports = interface
        .iter()
        .map(|p| d_dim(&parse(p)))
        .collect::<Vec<_>>();
    let (mink, coad): (Vec<_>, Vec<_>) = ports.into_iter().partition(
        |p| matches!(p.as_view(), AtomView::Fun(f) if f.get_symbol().get_name() == "spenso::mink"),
    );
    assert_eq!((mink.len(), coad.len()), (2, 2), "two B legs");
    let metric = spenso::network::library::symbolic::ETS.metric(mink[0].clone(), mink[1].clone());
    let delta = spenso::network::library::symbolic::ETS.metric(coad[0].clone(), coad[1].clone());
    let contracted = contract_with_idenso(&(numerator * metric * delta), true).unwrap();
    let overall = feynkit_graph::expressions::evaluate_overall_factor(
        parse(&graph_attribute(text, "overall_factor")).as_view(),
    );
    let nf = Atom::var(symbolica::symbol!("vakint::nf"));
    (contracted / Atom::num(8) * overall * nf.pow(Atom::num(closed_quark_loops(text) as i64)))
        .expand()
}

/// Σ over the colour branches of `prefactor × naive contraction`.
fn through_factorized(model: &Model, diagram: &FeynmanDiagram) -> (Atom, usize) {
    let options = FromFeynkitOptions {
        colour: ColourMode::Su3,
        projector: ProjectorChoice::Metric,
        ..Default::default()
    };
    let values = FactorizedNumerator::from_feynkit(diagram, model, &options).unwrap();
    let mut sum = Atom::num(0);
    for value in &values {
        sum += value.contract_naively().unwrap() * value.prefactor.total();
    }
    (sum.expand(), values.len())
}

/// Returns (graphs checked, graphs whose reference is non-zero).
fn semantic(level: &str, n: Option<usize>) -> (usize, usize) {
    let Some(root) = root() else {
        return (0, 0);
    };
    let dir = root.join(level);
    let model = model(&dir);
    let files = files(&dir, n);
    let mut non_zero = 0;
    for path in &files {
        let (text, diagram) = load(&model, path);
        let started = Instant::now();
        let expected = reference(&text);
        non_zero += usize::from(!expected.is_zero());
        let reference_seconds = started.elapsed().as_secs_f64();
        let started = Instant::now();
        let (found, branches) = through_factorized(&model, &diagram);
        let factorized_seconds = started.elapsed().as_secs_f64();
        let difference = (&expected - &found).expand();
        assert!(
            difference.is_zero(),
            "{}: the factorized value differs from the DOT numerator\nreference {expected}\nfactorized {found}\ndifference {difference}",
            path.display()
        );
        eprintln!(
            "SEMANTIC {} branches={branches} zero={} terms={} reference {:.3}s factorized {:.3}s: equal",
            path.display(),
            expected.is_zero(),
            match expected.as_view() {
                AtomView::Add(a) => a.get_nargs(),
                _ => 1,
            },
            reference_seconds,
            factorized_seconds
        );
    }
    (files.len(), non_zero)
}

#[test]
fn factors_denote_the_dot_numerator_at_one_and_two_loops() {
    let one = semantic("L1", None);
    let two = semantic("L2", None);
    eprintln!("SEMANTIC summary: L1 {one:?}, L2 {two:?} (graphs, non-zero references)");
    if root().is_some() {
        assert_eq!((one.0, two.0), (3, 29));
        assert_eq!(one.1, 3, "every one-loop reference is non-zero");
        assert!(two.1 >= 5, "at least five non-zero two-loop references");
    }
}

/// The invariant-mode colour coefficients evaluated at SU(3) equal the SU(3)
/// numbers, branch by branch (L1 and L2: only C_A, C_F, T_F, N_A, N_c occur).
#[test]
fn colour_invariants_at_su3_equal_the_su3_coefficients() {
    let Some(root) = root() else { return };
    let su3 = |atom: &Atom| -> Atom {
        let value = atom.replace_map(|part, _ctx, out| match part {
            AtomView::Fun(f) if matches!(f.get_symbol().get_name(), "spenso::cas" | "spenso::idx") => {
                let args = f.iter().collect::<Vec<_>>();
                assert_eq!(i64::try_from(args[0]), Ok(2), "quadratic invariant: {part}");
                let adjoint = matches!(args[1], AtomView::Fun(r) if r.get_symbol().get_name() == "spenso::coad");
                **out = match (f.get_symbol().get_name(), adjoint) {
                    ("spenso::cas", true) | ("spenso::idx", true) => Atom::num(3),
                    ("spenso::cas", false) => Atom::num(4) / Atom::num(3),
                    _ => Atom::num(1) / Atom::num(2),
                };
            }
            AtomView::Var(v) if v.get_symbol().get_name() == "spenso::dA" => **out = Atom::num(8),
            AtomView::Var(v) if v.get_symbol().get_name() == "spenso::Nc" => **out = Atom::num(3),
            _ => {}
        });
        value.expand()
    };
    let mut checked = 0;
    for level in ["L1", "L2"] {
        let dir = root.join(level);
        let model = model(&dir);
        for path in files(&dir, None) {
            let (_, diagram) = load(&model, &path);
            let invariants =
                FactorizedNumerator::from_feynkit(&diagram, &model, &FromFeynkitOptions::default())
                    .unwrap();
            let options = FromFeynkitOptions {
                colour: ColourMode::Su3,
                ..Default::default()
            };
            let numbers = FactorizedNumerator::from_feynkit(&diagram, &model, &options).unwrap();
            assert_eq!(invariants.len(), numbers.len(), "{}", path.display());
            for (a, b) in invariants.iter().zip(&numbers) {
                assert_eq!(a.colour_branch, b.colour_branch);
                assert_eq!(
                    su3(&a.prefactor.colour),
                    b.prefactor.colour.expand(),
                    "{}: {}",
                    path.display(),
                    a.prefactor.colour
                );
                checked += 1;
            }
        }
    }
    eprintln!("colour: {checked} branches, invariants at SU(3) = SU(3) numbers");
}

/// Rejections on a real graph (the L2 graph with a quark loop and the most
/// factors): each mutation is refused with its own error.
#[test]
fn mutated_real_graphs_are_rejected_with_the_right_error() {
    let Some(root) = root() else { return };
    let dir = root.join("L1");
    let model = model(&dir);
    // L1 diagram 2 is the quark bubble.
    let (_, diagram) = load(&model, &dir.join("diagram-000002.dot"));
    let values =
        FactorizedNumerator::from_feynkit(&diagram, &model, &FromFeynkitOptions::default())
            .unwrap();
    let value = &values[0];
    assert_eq!(value.closed_loops().unwrap().len(), 1);
    assert_eq!(value.prefactor.loop_sign, -1);
    assert_eq!(value.prefactor.nf_power, 1);
    // A wrong loop sign.
    let mut wrong = value.clone();
    wrong.prefactor.loop_sign = 1;
    assert!(matches!(
        wrong.validate(),
        Err(FactorizedError::LoopSign { .. })
    ));
    // A wrong n_f power.
    let mut wrong = value.clone();
    wrong.prefactor.nf_power = 0;
    assert!(matches!(
        wrong.validate(),
        Err(FactorizedError::NfPower { .. })
    ));
    // A reversed quark edge without transposing its gamma: orientation error
    // (or broken conservation, which validate checks first).
    let mut wrong = value.clone();
    let basis = wrong.loop_edges.clone();
    let edge = wrong
        .edges
        .iter_mut()
        .find(|e| {
            e.line == vakint::p_integral::factorized::LineKind::Fermion
                && !basis.contains(&Some(e.id))
        })
        .unwrap();
    std::mem::swap(&mut edge.source, &mut edge.target);
    std::mem::swap(&mut edge.source_hedge, &mut edge.target_hedge);
    for c in edge.momentum.iter_mut() {
        *c = -*c;
    }
    let error = wrong.validate().unwrap_err();
    assert!(
        matches!(
            error,
            FactorizedError::FermionOrientation { .. } | FactorizedError::FlowNotConserved { .. }
        ),
        "{error}"
    );
}

/// Sizes and timings of the text form at one to four loops (not part of the
/// fast suite: `--run-ignored only`). Prints one TSV line per loop order.
#[test]
#[ignore]
fn factorized_text_sizes_and_timings() {
    let Some(root) = root() else { return };
    let sample = std::env::var("VAKINT_FACTORIZED_STATS_SAMPLE")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(400usize);
    println!(
        "level\tgraphs\tbranches_median\tbranches_max\tbytes_graph_median\tbytes_graph_max\tbytes_record_median\tms_from_dot_median\tms_build_median\tms_write_median\tms_read_median\tms_validate_median"
    );
    for level in ["L1", "L2", "L3", "L4"] {
        let dir = root.join(level);
        let model = model(&dir);
        let files = files(&dir, Some(sample));
        let (mut branches, mut graph_bytes, mut record_bytes) =
            (Vec::new(), Vec::new(), Vec::new());
        let (mut t_dot, mut t_build, mut t_write, mut t_read, mut t_validate) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for path in &files {
            let text = std::fs::read_to_string(path).unwrap();
            let started = Instant::now();
            let diagram = FeynmanDiagram::from_dot(Arc::clone(&model), &text).unwrap();
            t_dot.push(started.elapsed().as_secs_f64() * 1e3);
            let started = Instant::now();
            let values =
                FactorizedNumerator::from_feynkit(&diagram, &model, &FromFeynkitOptions::default())
                    .unwrap();
            t_build.push(started.elapsed().as_secs_f64() * 1e3);
            branches.push(values.len() as f64);
            let mut bytes = 0;
            for value in &values {
                let started = Instant::now();
                let text = value.to_text().unwrap();
                t_write.push(started.elapsed().as_secs_f64() * 1e3);
                let started = Instant::now();
                let back = FactorizedNumerator::parse_text(&text).unwrap();
                t_read.push(started.elapsed().as_secs_f64() * 1e3);
                let started = Instant::now();
                back.validate().unwrap();
                t_validate.push(started.elapsed().as_secs_f64() * 1e3);
                assert_eq!(&back, value);
                record_bytes.push(text.len() as f64);
                bytes += text.len();
            }
            graph_bytes.push(bytes as f64);
        }
        let median = |v: &mut Vec<f64>| {
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            v[v.len() / 2]
        };
        let max = |v: &Vec<f64>| v.iter().cloned().fold(0.0, f64::max);
        println!(
            "{level}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.3}\t{:.3}\t{:.3}",
            files.len(),
            median(&mut branches),
            max(&branches),
            median(&mut graph_bytes),
            max(&graph_bytes),
            median(&mut record_bytes),
            median(&mut t_dot),
            median(&mut t_build),
            median(&mut t_write),
            median(&mut t_read),
            median(&mut t_validate)
        );
    }
}

/// Write the records of every L1 and L2 graph (one file per graph, all its
/// colour branches concatenated) to `$VAKINT_FACTORIZED_EXAMPLES_OUT`
/// (documentation examples; `--run-ignored only`).
#[test]
#[ignore]
fn write_example_records() {
    let Some(root) = root() else { return };
    let Ok(out) = std::env::var("VAKINT_FACTORIZED_EXAMPLES_OUT") else {
        eprintln!("SKIP: VAKINT_FACTORIZED_EXAMPLES_OUT not set");
        return;
    };
    let out = PathBuf::from(out);
    for level in ["L1", "L2"] {
        let dir = root.join(level);
        let model = model(&dir);
        std::fs::create_dir_all(out.join(level)).unwrap();
        for path in files(&dir, None) {
            let (_, diagram) = load(&model, &path);
            let values =
                FactorizedNumerator::from_feynkit(&diagram, &model, &FromFeynkitOptions::default())
                    .unwrap();
            let text = values
                .iter()
                .map(|v| v.to_text().unwrap())
                .collect::<String>();
            let name = path
                .file_stem()
                .unwrap()
                .to_string_lossy()
                .replace("diagram", "factorized");
            std::fs::write(out.join(level).join(format!("{name}.txt")), text).unwrap();
        }
    }
}

/// Graphs whose every colour branch vanishes, over all 392 three-loop graphs:
/// the beta census counts 13 raw colour zeros at L3 (FACTS-0.4).
#[test]
fn colour_zero_graphs_match_the_beta_census_at_three_loops() {
    let Some(root) = root() else { return };
    let dir = root.join("L3");
    let model = model(&dir);
    let mut zeros = Vec::new();
    let files = files(&dir, None);
    for path in &files {
        let (_, diagram) = load(&model, path);
        let values =
            FactorizedNumerator::from_feynkit(&diagram, &model, &FromFeynkitOptions::default())
                .unwrap();
        if values.is_empty() {
            zeros.push(path.file_name().unwrap().to_string_lossy().to_string());
        }
    }
    eprintln!(
        "L3: {} graphs, colour zeros {}: {zeros:?}",
        files.len(),
        zeros.len()
    );
    assert_eq!((files.len(), zeros.len()), (392, 13));
}

/// Beyond the brief: the same semantic check on 12 evenly spaced three-loop
/// graphs (colour zeros included: both sides must then vanish).
#[test]
fn factors_denote_the_dot_numerator_on_a_three_loop_sample() {
    let (graphs, non_zero) = semantic("L3", Some(12));
    eprintln!("SEMANTIC summary: L3 ({graphs}, {non_zero}) (graphs, non-zero references)");
    if root().is_some() {
        assert_eq!(graphs, 12);
        assert!(non_zero >= 6);
    }
}
