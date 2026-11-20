//! The sub-step timing counters of `vakint::timing` (contraction, lowering,
//! reduction, finalisation) on real graphs through the graph entry
//! `evaluate_graph`, single core: the four numbers are present, add up to at
//! most the wall time, and switching the resource counters on changes no
//! value. Data: the beta's prepared DOT files
//! (`VAKINT_TEST_FACTORIZED_DOT_ROOT`, default
//! /common/ben/beta/results/prepared) and, for four loops, the 1-4 loop p
//! library (`VAKINT_TEST_P_LIBRARY_1_4`); SKIP without them.
//! `VAKINT_TEST_TIMING_GRAPHS=L4/2519,L3/354` overrides the graphs.
#![cfg(feature = "feynkit-ingress")]

use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use feynkit_graph::FeynmanDiagram;
use vakint::p_integral::factorized::{FromFeynkitOptions, evaluate_graph};
use vakint::{Vakint, refpath, timing};

fn root() -> Option<PathBuf> {
    let root = std::env::var("VAKINT_TEST_FACTORIZED_DOT_ROOT")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("/common/ben/beta/results/prepared"));
    root.join("L2")
        .join("model.json")
        .is_file()
        .then_some(root)
        .or_else(|| {
            eprintln!("SKIP: no DOT data");
            None
        })
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

/// One graph, metric projection, through the graph entry: the value and the
/// interval's counters.
fn run(vakint: &Vakint, root: &Path, level: &str, number: usize) -> (symbolica::atom::Atom, timing::Delta, f64) {
    let directory = root.join(level);
    let model = refpath::load_model(&directory).unwrap();
    let text =
        std::fs::read_to_string(directory.join(format!("diagram-{number:06}.dot"))).unwrap();
    let diagram = FeynmanDiagram::from_dot(Arc::clone(&model), &text).unwrap();
    let loops = level[1..].parse::<i64>().unwrap();
    let settings = refpath::vakint_settings(loops + 1);
    let start = timing::Snapshot::now();
    let wall = Instant::now();
    let value = evaluate_graph(vakint, &settings, &diagram, &model, &FromFeynkitOptions::default())
        .unwrap();
    let wall = wall.elapsed().as_secs_f64();
    (value.value, start.elapsed(), wall)
}

fn graphs() -> Vec<(String, usize)> {
    let list = std::env::var("VAKINT_TEST_TIMING_GRAPHS").unwrap_or_else(|_| {
        if std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4").is_some() {
            "L4/100".into()
        } else {
            "L3/354".into()
        }
    });
    list.split(',')
        .map(|item| {
            let (level, number) = item.split_once('/').unwrap();
            (level.to_string(), number.parse().unwrap())
        })
        .collect()
}

/// The worked example: per graph the four exclusive wall numbers, their sum
/// against the wall time, and (resources on) the thread CPU and the increase
/// of the process peak RSS per sub-step; the value with resources on equals
/// the value with them off.
#[test]
fn the_four_sub_steps_of_a_graph_through_the_graph_entry() {
    let Some(root) = root() else { return };
    let vakint = engine();
    for (level, number) in graphs() {
        timing::set_resources(false);
        let (off, delta_off, wall_off) = run(&vakint, &root, &level, number);
        timing::set_resources(true);
        let (on, delta, wall) = run(&vakint, &root, &level, number);
        timing::set_resources(false);
        assert_eq!(off, on, "{level}/{number}: resources changed the value");
        for (d, w) in [(&delta_off, wall_off), (&delta, wall)] {
            let s = |stage| d.stage(stage).wall_ns as f64 * 1e-9;
            let sum: f64 = timing::Stage::ALL.iter().map(|x| s(*x)).sum();
            assert!(sum <= w + 1e-3, "{level}/{number}: {d:?} > wall {w}");
            for stage in timing::Stage::ALL {
                assert!(s(stage) > 0.0, "{level}/{number}: {stage:?} is 0 ({d:?})");
            }
            // P12's two counters keep their meaning: lower = lowering + the
            // sub-steps inside it, reduce = reduction + the finalisation inside it
            assert!(d.lower_seconds + 1e-9 >= s(timing::Stage::Lowering));
            assert!(d.reduce_seconds + 1e-9 >= s(timing::Stage::Reduction) + s(timing::Stage::Finalization) - 1e-3);
        }
        let mut line = format!("TIMING_SPLIT graph={level}/{number} wall_s={wall:.3}");
        let mut sum = 0.0;
        for stage in timing::Stage::ALL {
            let t = delta.stage(stage);
            sum += t.wall_ns as f64 * 1e-9;
            line.push_str(&format!(
                " {stage:?}=wall:{:.3},cpu:{:.3},peak_rss_gib:{:.2}",
                t.wall_ns as f64 * 1e-9,
                t.cpu_ns.unwrap_or(0) as f64 * 1e-9,
                t.peak_rss_bytes.unwrap_or(0) as f64 / (1u64 << 30) as f64
            ));
        }
        line.push_str(&format!(
            " sum_s={sum:.3} cpu_s={:.3} off_wall_s={wall_off:.3}{}",
            delta.cpu_seconds().unwrap_or(0.0),
            delta.fields()
        ));
        eprintln!("{line}");
    }
}
