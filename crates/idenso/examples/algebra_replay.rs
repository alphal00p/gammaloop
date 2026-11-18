//! Replay saved `simplify_algebra` inputs (helper VH, lane idenso-audit).
//!
//! Inputs: `<name>.bin` = the atom handed to `SymbolicTensor::infer_indexed`
//! (Symbolica binary export with its symbol state, written by vakint's
//! `examples/idenso_prof.rs --dump-dir`), optionally `<name>.out` = the
//! expected `to_dots` result as a canonical string (byte-identity check).
//! Settings: the reference chain's, `gamma: canonical(), epsilon: true`, then
//! `to_dots`.
//!
//! ```text
//! cargo run -p idenso --profile ci-optim --features reference-cases --example algebra_replay -- \
//!     <dir or list file of .bin paths> <out.jsonl> [reps=3] [phases]
//! ```
//! One JSON line per input: wall / thread-CPU medians of `infer_indexed`,
//! `simplify_algebra`, `to_dots`, sizes, whether the output equals `.out`, and
//! (with `phases`, feature `reference-cases`) idenso's phase clocks.
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

use idenso::dirac::GammaSimplifySettings;
use idenso::tensor::{AlgebraSettings, SymbolicTensor};
use symbolica::atom::{Atom, AtomCore};

fn cpu_seconds() -> f64 {
    let mut spec = Timespec::default();
    // SAFETY: clock_gettime writes one timespec into a valid pointer.
    unsafe { clock_gettime(3, &mut spec) }; // CLOCK_THREAD_CPUTIME_ID = 3 on Linux
    spec.tv_sec as f64 + spec.tv_nsec as f64 * 1e-9
}

#[repr(C)]
#[derive(Default)]
struct Timespec {
    tv_sec: i64,
    tv_nsec: i64,
}

unsafe extern "C" {
    fn clock_gettime(clock: i32, spec: *mut Timespec) -> i32;
}

fn median(mut values: Vec<f64>) -> f64 {
    values.sort_by(f64::total_cmp);
    let n = values.len();
    if n == 0 {
        0.0
    } else if n % 2 == 1 {
        values[n / 2]
    } else {
        (values[n / 2 - 1] + values[n / 2]) / 2.0
    }
}

fn inputs(source: &Path) -> Vec<PathBuf> {
    let mut paths = if source.is_dir() {
        std::fs::read_dir(source)
            .unwrap()
            .filter_map(|entry| entry.ok().map(|entry| entry.path()))
            .filter(|path| path.extension().is_some_and(|ext| ext == "bin"))
            .collect::<Vec<_>>()
    } else {
        std::fs::read_to_string(source)
            .unwrap()
            .split_whitespace()
            .map(PathBuf::from)
            .collect()
    };
    paths.sort();
    paths
}

fn main() {
    let arguments = std::env::args().skip(1).collect::<Vec<_>>();
    let source = PathBuf::from(&arguments[0]);
    let out_path = PathBuf::from(&arguments[1]);
    let reps: usize = arguments.get(2).map_or(3, |s| s.parse().unwrap());
    let phases = arguments.get(3).is_some_and(|s| s == "phases");
    #[cfg(not(feature = "reference-cases"))]
    assert!(!phases, "phases needs the feature reference-cases");
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::canonical()),
        epsilon: true,
        ..Default::default()
    };
    let mut out = std::io::BufWriter::new(std::fs::File::create(&out_path).unwrap());
    for path in inputs(&source) {
        let atom = Atom::import(&mut std::fs::File::open(&path).unwrap(), None).unwrap();
        let expected = std::fs::read_to_string(path.with_extension("out")).ok();
        let mut clocks = [(); 3].map(|_| (Vec::new(), Vec::new()));
        let mut output = None;
        let mut status = String::new();
        for _ in 0..reps {
            let input = atom.clone();
            let (w, c) = (Instant::now(), cpu_seconds());
            let inferred = SymbolicTensor::infer_indexed(input).unwrap();
            clocks[0].0.push(w.elapsed().as_secs_f64());
            clocks[0].1.push(cpu_seconds() - c);
            let (w, c) = (Instant::now(), cpu_seconds());
            let simplified = inferred.simplify_algebra(&settings).unwrap();
            clocks[1].0.push(w.elapsed().as_secs_f64());
            clocks[1].1.push(cpu_seconds() - c);
            status = format!("{:?}", simplified.reduction_status());
            let (w, c) = (Instant::now(), cpu_seconds());
            let reduced = simplified.to_dots().unwrap();
            clocks[2].0.push(w.elapsed().as_secs_f64());
            clocks[2].1.push(cpu_seconds() - c);
            output = Some(reduced.expression().clone());
        }
        let output = output.unwrap_or_default();
        let canonical = output.to_canonical_string();
        let identical = expected.as_ref().map(|text| *text == canonical);
        #[allow(unused_mut)]
        let mut phase_json = String::from("null");
        #[cfg(feature = "reference-cases")]
        if phases {
            let input = atom.clone();
            let ((), measurement) =
                idenso::reference_cases::timing::measure("replay", || {
                    let inferred = SymbolicTensor::infer_indexed(input).unwrap();
                    let simplified = inferred.simplify_algebra(&settings).unwrap();
                    let _ = simplified.to_dots().unwrap();
                });
            phase_json = measurement.to_string();
        }
        let stage = |index: usize| {
            format!(
                "{{\"wall\":{},\"cpu\":{}}}",
                median(clocks[index].0.clone()),
                median(clocks[index].1.clone())
            )
        };
        writeln!(
            out,
            "{{\"input\":\"{}\",\"input_bytes\":{},\"output_bytes\":{},\"output_terms\":{},\"status\":\"{}\",\"identical\":{},\"infer_indexed\":{},\"simplify_algebra\":{},\"to_dots\":{},\"phases\":{}}}",
            path.file_name().unwrap().to_string_lossy(),
            atom.as_view().get_byte_size(),
            output.as_view().get_byte_size(),
            match output.as_view() {
                symbolica::atom::AtomView::Add(add) => add.get_nargs(),
                _ => 1,
            },
            status,
            identical.map_or("null".to_string(), |b| b.to_string()),
            stage(0),
            stage(1),
            stage(2),
            phase_json
        )
        .unwrap();
        out.flush().unwrap();
        if identical == Some(false) {
            eprintln!("[algebra_replay] {}: output differs from .out", path.display());
        }
    }
}
