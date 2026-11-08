//! Benchmark v1 of the p-integrals: RustRed vs Forcer. Appends raw measurements
//! to `docs/benchmark/data/p_results.json`, from which `docs/benchmark/benchmark.typ`
//! builds the report (M2 plan, section 9). See `--help`.
//!
//! Every measurement is one fresh child process of this executable, one input
//! at a time, under the thread count and CPU list of its table. A measurement
//! already recorded for the same gammaloop and rustred shas, build profile and
//! tree state is skipped, so an interrupted run resumes by running the same
//! command (results committed in between do not count as a source change).
//! SIGTERM and SIGINT kill the measuring process group before the harness
//! exits: the group sits outside the harness's own, so it would survive.
//!
//! Build with `just vakint-p-benchmark-build`, which stamps the binary with the
//! source it is built from; a release measurement refuses any other binary.

#[cfg(all(feature = "experimental-rustred", target_os = "linux"))]
#[path = "benchmark_p_integrals/results.rs"]
mod results;

#[cfg(all(feature = "experimental-rustred", target_os = "linux"))]
fn main() -> eyre::Result<()> {
    harness::main()
}

#[cfg(not(all(feature = "experimental-rustred", target_os = "linux")))]
fn main() {
    eprintln!("benchmark_p_integrals needs Linux and the feature vakint/experimental-rustred");
    std::process::exit(2);
}

#[cfg(all(feature = "experimental-rustred", target_os = "linux"))]
mod harness {
    use std::collections::{BTreeMap, BTreeSet};
    use std::fs::{self, File};
    use std::os::unix::process::CommandExt;
    use std::path::{Path, PathBuf};
    use std::process::{Child, Command, ExitStatus, Stdio};
    use std::sync::atomic::{AtomicI32, Ordering};
    use std::time::{Duration, Instant};

    use eyre::{Context, Result, bail, ensure};
    use serde::Deserialize;
    use serde_json::{Value, json};
    use sha2::{Digest, Sha256};
    use symbolica::atom::{AtomCore, AtomView};
    use vakint::{
        EvaluationMethod, EvaluationOrder, ForcerOptions, OrderedCacheMode, OrderedCacheSettings,
        OrderedCacheStatistics, OrderedWeightArithmetic, RustRedEvaluationOptions,
        RustRedMasterMode, Vakint, VakintSettings,
        p_integral::{FORCER_UPSTREAM_COMMIT, forcer_benchmark},
        rustred_evaluation::experimental::capture_reduction_timing,
        vakint_parse,
    };

    use super::results::{self, Key, NO_ARITHMETIC, Source};

    const GRID: &str = include_str!("../tests/fixtures/benchmark_p_integrals.json");
    const PHYSICS: &str = include_str!("../tests/fixtures/benchmark_p_physics.json");
    const BUBBLE: &str = "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,1))";
    const SCOPE: &str = "Wall time of Vakint::evaluate in Retain mode, in a fresh process, input \
        parsed, engine (and RustRed library) loaded and the one-loop bubble evaluated once \
        beforehand by either backend, on directories of its own: same frontend and finalizer for \
        both backends. RustRed cold: empty ordered-cache directory, shipped caches off. \
        Warm-same: fresh process replaying a read-only copy of the cache its cold run filled. \
        Warm-shared: per L one chain over the cells in table order sharing one cache from empty. \
        Forcer cold: private empty materialisation directory; warm: materialised (Forcer has no \
        result cache). Raw: the generated FORM programs run directly, one FORM process each. OS \
        pages are not flushed.";

    /// The variable that swaps RustRed's Ordered adjoint tail; recorded per
    /// RustRed child, never set by the harness.
    const ADJOINT_TAIL: &str = "RUSTRED_ORDERED_ADJOINT_TAIL";

    /// The source this binary was built from; see `results::check_build`.
    const BUILD_STAMP: Option<&str> = option_env!("VAKINT_P_BENCHMARK_SOURCE");

    /// The signal that asked this harness to stop, 0 before any.
    static STOP_SIGNAL: AtomicI32 = AtomicI32::new(0);

    extern "C" fn note_stop_signal(signal: libc::c_int) {
        STOP_SIGNAL.store(signal, Ordering::SeqCst);
    }

    /// Note SIGTERM and SIGINT, for `run_capped` to act on.
    fn catch_stop_signals() {
        for signal in [libc::SIGTERM, libc::SIGINT] {
            // SAFETY: the handler only stores into an atomic.
            unsafe { libc::signal(signal, note_stop_signal as *const () as libc::sighandler_t) };
        }
    }

    /// SIGKILL the process group `child` leads, then reap it.
    fn kill_group(child: &mut Child) -> Result<()> {
        // SAFETY: a plain system call; the group is the one of the unreaped
        // `child`, so its id cannot have been reused.
        unsafe { libc::kill(-(child.id() as libc::pid_t), libc::SIGKILL) };
        child.wait()?;
        Ok(())
    }

    #[derive(Clone, Debug, Deserialize)]
    struct Cell {
        name: String,
        loops: usize,
        #[serde(default)]
        dots: Option<usize>,
        #[serde(default)]
        rank: Option<usize>,
        input: String,
    }

    #[derive(Deserialize)]
    struct Cells {
        cells: Vec<Cell>,
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum Backend {
        RustRed,
        Forcer,
    }

    impl Backend {
        const fn name(self) -> &'static str {
            match self {
                Self::RustRed => "rustred",
                Self::Forcer => "forcer",
            }
        }

        fn parse(name: &str) -> Result<Self> {
            match name {
                "rustred" => Ok(Self::RustRed),
                "forcer" => Ok(Self::Forcer),
                _ => bail!("unknown backend {name}"),
            }
        }

        const fn states(self) -> &'static [&'static str] {
            match self {
                Self::RustRed => &["cold", "warm-same", "warm-shared"],
                Self::Forcer => &["cold", "warm", "raw"],
            }
        }
    }

    /// The directory of one child: the Forcer materialisation base, and for
    /// RustRed the ordered cache and how it is opened.
    struct Placement {
        base: PathBuf,
        cache: Option<PathBuf>,
        read_only: bool,
    }

    struct Options {
        output: PathBuf,
        exact: PathBuf,
        work_dir: PathBuf,
        form: PathBuf,
        child_program: Option<PathBuf>,
        build_profile: Option<String>,
        repetitions: usize,
        threads: Vec<usize>,
        cpus: BTreeMap<usize, String>,
        numa_node: Option<usize>,
        numactl: PathBuf,
        loops: Vec<usize>,
        dots: Vec<usize>,
        ranks: Vec<usize>,
        physics: Vec<usize>,
        backends: Vec<Backend>,
        states: Vec<String>,
        arithmetics: Vec<String>,
        cell_names: Vec<String>,
        cap_seconds: Option<f64>,
        stop_file: Option<PathBuf>,
        select_current: bool,
        child: Option<String>,
        backend: Backend,
        child_threads: usize,
        base: Option<PathBuf>,
        cache: Option<PathBuf>,
        read_only: bool,
        report: Option<PathBuf>,
        programs: Option<PathBuf>,
    }

    fn list<T: std::str::FromStr>(text: &str) -> Result<Vec<T>>
    where
        T::Err: std::error::Error + Send + Sync + 'static,
    {
        text.split(',')
            .filter(|item| !item.is_empty())
            .map(|item| Ok(item.parse()?))
            .collect()
    }

    const USAGE: &str = "benchmark_p_integrals --build-profile release|ci-optim [OPTIONS]\n\
        --output PATH        results (default: crate docs/benchmark/data/p_results.json)\n\
        --exact PATH         exact marks (default: crate docs/benchmark/data/p_exact.json)\n\
        --work-dir PATH      scratch directory, kept for resuming (default: target/p_benchmark of the repository)\n\
        --form PATH          FORM executable; tform is its sibling (default: form on PATH)\n\
        --repetitions N      fresh processes per measurement (default: 3)\n\
        --threads LIST       thread counts (default: 1,16)\n\
        --cpus T=CPUS,...    taskset list of a thread count, e.g. 1=128,16=128-143\n\
        --numa-node N        numactl --cpunodebind=N --membind=N around every child\n\
        --numactl PATH       numactl executable (default: numactl on PATH)\n\
        --loops LIST         grid loop counts (default: 1,2,3)\n\
        --dots LIST          grid dots (default: 0,1,2)\n\
        --ranks LIST         grid ranks (default: 0,2,4,6,8)\n\
        --physics LIST       physics cells by L (default: 2,3; empty for none)\n\
        --backends LIST      rustred,forcer (default: both)\n\
        --states LIST        rustred: cold,warm-same,warm-shared; forcer: cold,warm,raw (default: all)\n\
        --arithmetic LIST    RustRed Ordered weight arithmetics: exact,modular_polynomial,modular_sampled\n\
                             (default: modular_polynomial); Forcer is measured once, whatever the list\n\
        --cells LIST         only these cells (names as in the fixtures) of the selected loop counts\n\
        --cap-seconds S      stop a measurement after S seconds and record it as capped\n\
        --stop-file PATH     stop cleanly between measurements when PATH exists (exit status 3)\n\
        --child-program PATH run PATH instead of this executable as the measuring child (harness tests,\n\
                             refused with --build-profile release)\n\
        --select-current     select this run for the report; recorded runs stay\n\
        A release measurement needs a clean tree (crates/vakint/docs/benchmark excepted) and a binary\n\
        of that tree, built by `just vakint-p-benchmark-build`.\n\
        Never removes anything outside its work directory; every Forcer run materialises under it.";

    impl Options {
        fn parse() -> Result<Option<Self>> {
            let crate_path = Path::new(env!("CARGO_MANIFEST_DIR"));
            let mut result = Self {
                output: crate_path.join("docs/benchmark/data/p_results.json"),
                exact: crate_path.join("docs/benchmark/data/p_exact.json"),
                work_dir: crate_path.join("../../target/p_benchmark"),
                form: "form".into(),
                child_program: None,
                build_profile: None,
                repetitions: 3,
                threads: vec![1, 16],
                cpus: BTreeMap::new(),
                numa_node: None,
                numactl: "numactl".into(),
                loops: vec![1, 2, 3],
                dots: vec![0, 1, 2],
                ranks: vec![0, 2, 4, 6, 8],
                physics: vec![2, 3],
                backends: vec![Backend::RustRed, Backend::Forcer],
                states: Vec::new(),
                arithmetics: vec!["modular_polynomial".to_owned()],
                cell_names: Vec::new(),
                cap_seconds: None,
                stop_file: None,
                select_current: false,
                child: None,
                backend: Backend::Forcer,
                child_threads: 1,
                base: None,
                cache: None,
                read_only: false,
                report: None,
                programs: None,
            };
            let mut args = std::env::args().skip(1);
            while let Some(arg) = args.next() {
                match arg.as_str() {
                    "--help" | "-h" => {
                        println!("{USAGE}");
                        return Ok(None);
                    }
                    "--select-current" => {
                        result.select_current = true;
                        continue;
                    }
                    "--read-only" => {
                        result.read_only = true;
                        continue;
                    }
                    _ => {}
                }
                let value = args
                    .next()
                    .ok_or_else(|| eyre::eyre!("missing value for {arg}"))?;
                match arg.as_str() {
                    "--output" => result.output = value.into(),
                    "--exact" => result.exact = value.into(),
                    "--work-dir" => result.work_dir = value.into(),
                    "--form" => result.form = value.into(),
                    "--child-program" => result.child_program = Some(value.into()),
                    "--build-profile" => {
                        ensure!(
                            matches!(value.as_str(), "release" | "ci-optim"),
                            "build profile must be release or ci-optim"
                        );
                        result.build_profile = Some(value);
                    }
                    "--repetitions" => result.repetitions = value.parse()?,
                    "--threads" => result.threads = list(&value)?,
                    "--cpus" => {
                        for item in value.split(',').filter(|item| !item.is_empty()) {
                            let (threads, cpus) = item
                                .split_once('=')
                                .ok_or_else(|| eyre::eyre!("--cpus takes THREADS=CPULIST items"))?;
                            result.cpus.insert(threads.parse()?, cpus.into());
                        }
                    }
                    "--numa-node" => result.numa_node = Some(value.parse()?),
                    "--numactl" => result.numactl = value.into(),
                    "--loops" => result.loops = list(&value)?,
                    "--dots" => result.dots = list(&value)?,
                    "--ranks" => result.ranks = list(&value)?,
                    "--physics" => result.physics = list(&value)?,
                    "--backends" => {
                        result.backends = value
                            .split(',')
                            .map(Backend::parse)
                            .collect::<Result<_>>()?;
                    }
                    "--states" => result.states = list(&value)?,
                    "--arithmetic" => result.arithmetics = list(&value)?,
                    "--cells" => result.cell_names = list(&value)?,
                    "--cap-seconds" => result.cap_seconds = Some(value.parse()?),
                    "--stop-file" => result.stop_file = Some(value.into()),
                    "--child" => result.child = Some(value),
                    "--backend" => result.backend = Backend::parse(&value)?,
                    "--child-threads" => result.child_threads = value.parse()?,
                    "--base" => result.base = Some(value.into()),
                    "--cache" => result.cache = Some(value.into()),
                    "--report" => result.report = Some(value.into()),
                    "--programs" => result.programs = Some(value.into()),
                    _ => bail!("unknown option {arg}; see --help"),
                }
            }
            ensure!(result.repetitions > 0, "repetitions must be positive");
            ensure!(!result.arithmetics.is_empty(), "--arithmetic needs a mode");
            for name in &result.arithmetics {
                let selection: OrderedWeightArithmetic = name.parse()?;
                ensure!(
                    selection.name() == name,
                    "--arithmetic takes mode names (exact, modular_polynomial, modular_sampled), not {name}"
                );
            }
            ensure!(
                result.child.is_some() || result.build_profile.is_some(),
                "declare the actual Cargo profile with --build-profile release or ci-optim"
            );
            for state in &result.states {
                ensure!(
                    result
                        .backends
                        .iter()
                        .any(|backend| backend.states().contains(&state.as_str())),
                    "no selected backend has the state {state}"
                );
            }
            Ok(Some(result))
        }

        fn selected(&self, backend: Backend, state: &str) -> bool {
            self.backends.contains(&backend)
                && backend.states().contains(&state)
                && (self.states.is_empty() || self.states.iter().any(|name| name == state))
        }

        /// The selected cells of one L in table order: the grid, then physics;
        /// with `--cells`, exactly the named cells of that L.
        fn cells(&self, cells: &[Cell], loops: usize) -> Vec<Cell> {
            if !self.cell_names.is_empty() {
                return cells
                    .iter()
                    .filter(|cell| {
                        cell.loops == loops
                            && self.loops.contains(&loops)
                            && self.cell_names.contains(&cell.name)
                    })
                    .cloned()
                    .collect();
            }
            let grid = cells.iter().filter(|cell| {
                cell.loops == loops
                    && self.loops.contains(&loops)
                    && cell.dots.is_some_and(|dots| self.dots.contains(&dots))
                    && cell.rank.is_some_and(|rank| self.ranks.contains(&rank))
            });
            let physics = cells.iter().filter(|cell| {
                cell.loops == loops && cell.dots.is_none() && self.physics.contains(&loops)
            });
            grid.chain(physics).cloned().collect()
        }
    }

    fn all_cells() -> Result<Vec<Cell>> {
        let mut cells = serde_json::from_str::<Cells>(GRID)?.cells;
        cells.extend(serde_json::from_str::<Cells>(PHYSICS)?.cells);
        Ok(cells)
    }

    /// The settings of one measurement; both backends retain the masters and
    /// keep every other default but RustRed's arithmetic.
    fn settings(options: &Options) -> Result<VakintSettings> {
        let retain = RustRedMasterMode::Retain;
        let (method, form) = match options.backend {
            Backend::RustRed => (
                EvaluationMethod::RustRed(RustRedEvaluationOptions {
                    master_mode: retain,
                }),
                "/rustred-benchmark-must-not-use-form".to_owned(),
            ),
            Backend::Forcer => (
                EvaluationMethod::Forcer(ForcerOptions {
                    master_mode: retain,
                }),
                options.form.to_string_lossy().into_owned(),
            ),
        };
        let mode = match (&options.cache, options.read_only) {
            (None, _) => OrderedCacheMode::Off,
            (Some(_), true) => OrderedCacheMode::ReadOnly,
            (Some(_), false) => OrderedCacheMode::ReadWrite,
        };
        let ordered_weight_arithmetic = match options.backend {
            Backend::RustRed => options.arithmetics[0].parse()?,
            Backend::Forcer => OrderedWeightArithmetic::default(),
        };
        Ok(VakintSettings {
            evaluation_order: EvaluationOrder(vec![method]),
            ordered_weight_arithmetic,
            form_exe_path: form,
            temporary_directory: options
                .base
                .as_ref()
                .map(|base| base.to_string_lossy().into_owned()),
            n_cores: options.child_threads,
            ordered_cache: OrderedCacheSettings {
                mode,
                directory: options.cache.clone(),
                shipped: false,
                ..OrderedCacheSettings::off()
            },
            ..Default::default()
        })
    }

    /// The cache counters of a RustRed child after its measured evaluate: the
    /// totals, and per attached program its own with how a shipped cache was
    /// used (`null`: none, as the harness turns shipped caches off).
    fn cache_record(statistics: &OrderedCacheStatistics) -> Value {
        let programs = statistics
            .programs
            .iter()
            .map(|program| {
                let native = &program.native;
                json!({
                    "loop_count": program.loop_count,
                    "program": program.program,
                    "mode": format!("{:?}", program.mode),
                    "shipped": program.shipped.as_ref().map(|shipped| format!("{shipped:?}")),
                    "hits": native.hits,
                    "misses": native.misses,
                    "inserts": native.inserts,
                    "records": native.records,
                    "file_bytes": native.file_bytes,
                    "loaded_shards": native.loaded_shards,
                    "load_seconds": native.load_seconds,
                })
            })
            .collect::<Vec<_>>();
        json!({
            "hits": statistics.hits,
            "misses": statistics.misses,
            "inserts": statistics.inserts,
            "records": statistics.records,
            "file_bytes": statistics.file_bytes,
            "resident_bytes": statistics.resident_bytes,
            "programs": programs,
        })
    }

    /// One measurement in this fresh process.
    fn child(options: &Options, name: &str) -> Result<()> {
        let cells = all_cells()?;
        let cell = cells
            .iter()
            .find(|cell| cell.name == name)
            .ok_or_else(|| eyre::eyre!("no cell {name}"))?;
        let report = options
            .report
            .as_ref()
            .ok_or_else(|| eyre::eyre!("a child needs --report"))?;
        let started = Instant::now();
        let engine = Vakint::new()?;
        let input = vakint_parse!(&cell.input).map_err(eyre::Report::msg)?;
        let settings = settings(options)?;
        // One evaluate before the clock, on either backend: what a process does
        // once (the RustRed library load, the FORM version probe, tables built
        // on first use) stays out of the measurement. It touches no directory
        // of the measurement: no cache, and Forcer materialises apart.
        let base = options
            .base
            .as_ref()
            .ok_or_else(|| eyre::eyre!("a child needs --base"))?;
        let apart = tempfile::tempdir_in(base)?;
        let load = VakintSettings {
            ordered_cache: OrderedCacheSettings::off(),
            temporary_directory: Some(apart.path().to_string_lossy().into_owned()),
            ..settings.clone()
        };
        let bubble = vakint_parse!(BUBBLE).map_err(eyre::Report::msg)?;
        drop(engine.evaluate(&load, bubble.as_view())?);
        apart.close()?;
        let load_seconds = started.elapsed().as_secs_f64();
        forcer_benchmark::enable();
        let started = Instant::now();
        let capture = capture_reduction_timing(|| engine.evaluate(&settings, input.as_view()));
        let evaluate_seconds = started.elapsed().as_secs_f64();
        let value = capture.result?;
        let (form, programs) = forcer_benchmark::take();
        if let Some(directory) = &options.programs {
            write_programs(directory, &programs)?;
        }
        let timing = capture.timing;
        let reduction_seconds = match options.backend {
            Backend::RustRed => results::reduction_seconds(
                timing.complete,
                timing.interval_count,
                timing.reduction_math.map(|math| math.as_secs_f64()),
            ),
            Backend::Forcer => Some(form.as_secs_f64()),
        };
        let output = value.to_canonical_string();
        let terms = match value.as_view() {
            AtomView::Add(sum) => sum.get_nargs(),
            _ => 1,
        };
        let mut record = json!({
            "evaluate_seconds": evaluate_seconds,
            "reduction_seconds": reduction_seconds,
            "load_seconds": load_seconds,
            "form_runs": programs.len(),
            "output_sha256": format!("{:x}", Sha256::digest(output.as_bytes())),
            "output_terms": terms,
        });
        match options.backend {
            Backend::RustRed => {
                // The arithmetic as this child resolved it, which must be the
                // requested mode.
                let resolved = settings.rustred_ordered_weight_arithmetic();
                ensure!(
                    resolved.name() == options.arithmetics[0],
                    "asked for {}, resolved {resolved}",
                    options.arithmetics[0]
                );
                record["arithmetic"] = json!(resolved.name());
                record["arithmetic_resolved"] = json!(resolved.to_string());
                record["cache"] = cache_record(&engine.ordered_cache_statistics());
                record["adjoint_tail_env"] = json!(std::env::var(ADJOINT_TAIL).ok());
            }
            Backend::Forcer => record["arithmetic"] = json!(NO_ARITHMETIC),
        }
        if options.backend == Backend::RustRed {
            record["reduction_gross_seconds"] = json!(timing.gross.as_secs_f64());
            record["reduction_interval_count"] = json!(timing.interval_count);
            record["reduction_timing_complete"] = json!(timing.complete);
            record["reduction_incompleteness_reason"] = json!(timing.incompleteness_reason);
        }
        serde_json::to_writer(File::create(report)?, &record)?;
        Ok(())
    }

    /// Replace the program directory at once: a killed child leaves either the
    /// former programs or none, never a part of them.
    fn write_programs(directory: &Path, programs: &[String]) -> Result<()> {
        let partial = directory.with_extension(format!("partial-{}", std::process::id()));
        fs::create_dir_all(&partial)?;
        for (index, program) in programs.iter().enumerate() {
            fs::write(partial.join(format!("{index:04}.frm")), program)?;
        }
        if directory.exists() {
            fs::remove_dir_all(directory)?;
        }
        fs::rename(partial, directory)?;
        Ok(())
    }

    fn command_output(program: &str, args: &[&str]) -> Result<String> {
        let output = Command::new(program).args(args).output()?;
        ensure!(output.status.success(), "{program} {args:?} failed");
        Ok(String::from_utf8(output.stdout)?.trim().to_owned())
    }

    fn utc_now() -> Result<String> {
        command_output("date", &["-u", "+%Y-%m-%dT%H:%M:%SZ"])
    }

    fn sha256_of_file(path: &Path) -> Result<String> {
        Ok(format!("{:x}", Sha256::digest(fs::read(path)?)))
    }

    fn executable(path: &Path) -> Result<PathBuf> {
        if path.is_file() {
            return Ok(path.canonicalize()?);
        }
        if path.components().count() == 1 {
            for directory in std::env::split_paths(&std::env::var_os("PATH").unwrap_or_default()) {
                if directory.join(path).is_file() {
                    return Ok(directory.join(path).canonicalize()?);
                }
            }
        }
        bail!("FORM executable not found: {}", path.display())
    }

    /// The rustred commit pinned in the workspace `Cargo.lock`.
    fn rustred_sha(crate_path: &Path) -> Result<String> {
        let lock = fs::read_to_string(crate_path.join("../../Cargo.lock"))?;
        lock.lines()
            .find_map(|line| {
                let source = line.strip_prefix("source = \"git+")?;
                source.contains("/rustred?").then_some(source)?;
                Some(source.rsplit_once('#')?.1.trim_end_matches('"').to_owned())
            })
            .ok_or_else(|| eyre::eyre!("Cargo.lock pins no rustred commit"))
    }

    /// The thread variables of every child, explicit whatever the caller's
    /// environment holds.
    fn thread_environment(threads: usize) -> [(&'static str, String); 3] {
        ["RAYON_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"]
            .map(|name| (name, threads.to_string()))
    }

    /// The digest of the exact-agreement file, which must mark `source` (or be
    /// the stub, or absent).
    fn check_exact(path: &Path, source: Source) -> Result<Option<String>> {
        if !path.exists() {
            return Ok(None);
        }
        let exact: Value = serde_json::from_reader(File::open(path)?)?;
        results::check_exact(&exact, source).with_context(|| path.display().to_string())?;
        Ok(Some(sha256_of_file(path)?))
    }

    /// `form_command` falls back to plain `form` without a sibling `tform`,
    /// which would measure one worker under a thread count.
    fn check_tform(form: &Path) -> Result<()> {
        let tform = form.with_file_name("tform");
        let usable = Command::new(&tform)
            .arg("-version")
            .output()
            .is_ok_and(|output| output.status.success());
        ensure!(usable, "no usable tform at {}", tform.display());
        Ok(())
    }

    /// Forcer, which has no arithmetic to choose.
    const FORCER: (Backend, &str) = (Backend::Forcer, NO_ARITHMETIC);

    struct Harness<'a> {
        options: &'a Options,
        form: PathBuf,
        run: String,
        output: PathBuf,
        executable: PathBuf,
        recorded: BTreeSet<Key>,
        outputs: BTreeMap<(String, String), String>,
        integrals: BTreeMap<String, usize>,
        work: PathBuf,
    }

    /// Remove a directory of this harness, read-only copies included.
    fn remove_tree(path: &Path) -> Result<()> {
        if path.exists() {
            let status = Command::new("chmod")
                .args(["-R", "u+w"])
                .arg(path)
                .status()?;
            ensure!(status.success(), "chmod failed on {}", path.display());
            fs::remove_dir_all(path)?;
        }
        Ok(())
    }

    impl Harness<'_> {
        fn stop_requested(&self) -> bool {
            self.options
                .stop_file
                .as_ref()
                .is_some_and(|path| path.exists())
        }

        /// `arithmetic` of RustRed; Forcer's is always [`NO_ARITHMETIC`].
        fn key(
            &self,
            cell: &Cell,
            (backend, arithmetic): (Backend, &str),
            threads: usize,
            state: &str,
            rep: usize,
        ) -> Key {
            let arithmetic = match backend {
                Backend::RustRed => arithmetic,
                Backend::Forcer => NO_ARITHMETIC,
            };
            Key {
                cell: cell.name.clone(),
                backend: backend.name().into(),
                arithmetic: arithmetic.into(),
                threads,
                state: state.into(),
                repetition: rep,
            }
        }

        /// `numactl` and `taskset` words for the CPU list of a thread count,
        /// then `program`.
        fn pinned(&self, threads: usize, program: &Path) -> Command {
            let mut words = Vec::new();
            if let Some(node) = self.options.numa_node {
                words.extend([
                    self.options.numactl.to_string_lossy().into_owned(),
                    format!("--cpunodebind={node}"),
                    format!("--membind={node}"),
                ]);
            }
            if let Some(cpus) = self.options.cpus.get(&threads) {
                words.extend(["taskset".to_owned(), "-c".to_owned(), cpus.clone()]);
            }
            words.push(program.to_string_lossy().into_owned());
            let mut command = Command::new(&words[0]);
            command.args(&words[1..]);
            command
        }

        /// Run `command` in its own process group, so that the cap also stops
        /// its FORM processes. `None` when it hit the cap after `spent` seconds
        /// of the same measurement. A stop signal kills the group and ends the
        /// harness: the measurement in flight is not recorded. The wait polls
        /// every millisecond, which is also the resolution of the raw seconds.
        fn run_capped(
            &self,
            command: &mut Command,
            spent: f64,
        ) -> Result<Option<(ExitStatus, f64)>> {
            let stopped = |signal| eyre::eyre!("stopped by signal {signal}");
            match STOP_SIGNAL.load(Ordering::SeqCst) {
                0 => {}
                signal => return Err(stopped(signal)),
            }
            let started = Instant::now();
            let mut child = command.process_group(0).spawn()?;
            let status = loop {
                if let Some(status) = child.try_wait()? {
                    break status;
                }
                let signal = STOP_SIGNAL.load(Ordering::SeqCst);
                if signal != 0 {
                    kill_group(&mut child)?;
                    return Err(stopped(signal));
                }
                let elapsed = started.elapsed().as_secs_f64();
                if self
                    .options
                    .cap_seconds
                    .is_some_and(|cap| spent + elapsed > cap)
                {
                    kill_group(&mut child)?;
                    return Ok(None);
                }
                std::thread::sleep(Duration::from_millis(1));
            };
            Ok(Some((status, started.elapsed().as_secs_f64())))
        }

        /// One fresh measuring process. `None` when it hit the cap.
        fn run_child(
            &self,
            cell: &Cell,
            (backend, arithmetic): (Backend, &str),
            threads: usize,
            placement: &Placement,
            label: &str,
        ) -> Result<Option<(Value, f64)>> {
            let logs = self.work.join("logs");
            fs::create_dir_all(&logs)?;
            let report = logs.join(format!("{label}.json"));
            let log = File::create(logs.join(format!("{label}.log")))?;
            let mut command = self.pinned(threads, &self.executable);
            command
                .arg("--child")
                .arg(&cell.name)
                .arg("--backend")
                .arg(backend.name())
                .arg("--child-threads")
                .arg(threads.to_string())
                .arg("--base")
                .arg(&placement.base)
                .arg("--report")
                .arg(&report)
                .arg("--form")
                .arg(&self.form)
                .envs(thread_environment(threads))
                .stdin(Stdio::null())
                .stdout(Stdio::from(log.try_clone()?))
                .stderr(Stdio::from(log));
            if let Some(cache) = &placement.cache {
                command.arg("--cache").arg(cache);
            }
            if placement.read_only {
                command.arg("--read-only");
            }
            if backend == Backend::RustRed {
                command.arg("--arithmetic").arg(arithmetic);
            }
            if backend == Backend::Forcer {
                command
                    .arg("--programs")
                    .arg(self.work.join("programs").join(&cell.name));
            }
            let Some((status, process_seconds)) = self.run_capped(&mut command, 0.0)? else {
                return Ok(None);
            };
            ensure!(
                status.success(),
                "{label} failed ({status}); see {}",
                logs.join(format!("{label}.log")).display()
            );
            let report = serde_json::from_reader(File::open(report)?)?;
            Ok(Some((report, process_seconds)))
        }

        /// The fields every record of `key` has.
        fn measurement(&self, key: &Key, cell: &Cell) -> Result<Value> {
            Ok(json!({
                "cell": key.cell,
                "loops": cell.loops,
                "dots": cell.dots,
                "rank": cell.rank,
                "n": self.integrals[&cell.name],
                "backend": key.backend,
                "arithmetic": key.arithmetic,
                "threads": key.threads,
                "state": key.state,
                "repetition": key.repetition,
                "finished_utc": utc_now()?,
            }))
        }

        /// Record `key` as not measurable after a capped run, as `reason` says.
        fn record_skipped(&mut self, key: Key, cell: &Cell, reason: String) -> Result<()> {
            let mut measurement = self.measurement(&key, cell)?;
            measurement["status"] = json!("skipped");
            measurement["reason"] = json!(reason);
            measurement["cap_seconds"] = json!(self.options.cap_seconds);
            eprintln!("{} {} {}: {reason}", key.cell, key.backend, key.state);
            results::append(&self.output, &self.run, measurement)?;
            self.recorded.insert(key);
            Ok(())
        }

        fn record(&mut self, key: Key, cell: &Cell, measured: Option<(Value, f64)>) -> Result<()> {
            let mut measurement = self.measurement(&key, cell)?;
            match measured {
                Some((report, process_seconds)) => {
                    measurement["status"] = json!("ok");
                    measurement["process_seconds"] = json!(process_seconds);
                    let digest = report["output_sha256"].as_str().map(str::to_owned);
                    if let Some(digest) = digest {
                        let result = (key.cell.clone(), key.backend.clone());
                        let first = self.outputs.entry(result).or_insert_with(|| digest.clone());
                        ensure!(
                            *first == digest,
                            "{key:?} returned another result than the earlier measurements"
                        );
                    }
                    for (name, value) in report.as_object().expect("a child report is an object") {
                        measurement[name] = value.clone();
                    }
                    eprintln!(
                        "{} {} {} {} threads {} repetition {}: {} s",
                        key.cell,
                        key.backend,
                        key.arithmetic,
                        key.state,
                        key.threads,
                        key.repetition,
                        measurement["evaluate_seconds"]
                    );
                }
                None => {
                    measurement["status"] = json!("capped");
                    measurement["cap_seconds"] = json!(self.options.cap_seconds);
                    eprintln!("{} {} {}: capped", key.cell, key.backend, key.state);
                }
            }
            results::append(&self.output, &self.run, measurement)?;
            self.recorded.insert(key);
            Ok(())
        }

        /// The materialised Forcer base of the warm state, filled on first use
        /// by an unmeasured run of the smallest cell.
        fn warm_base(&self, threads: usize) -> Result<PathBuf> {
            let base = self.work.join("forcer-warm");
            let materialised = base.exists()
                && fs::read_dir(&base)?.flatten().any(|entry| {
                    entry
                        .file_name()
                        .to_string_lossy()
                        .starts_with("vakint-forcer-")
                });
            if !materialised {
                fs::create_dir_all(&base)?;
                let smallest = all_cells()?.swap_remove(0);
                let placement = Placement {
                    base: base.clone(),
                    cache: None,
                    read_only: false,
                };
                let label = "forcer-warm-fill";
                let filled = self.run_child(&smallest, FORCER, threads, &placement, label)?;
                ensure!(filled.is_some(), "the warm-state fill hit the cap");
            }
            Ok(base)
        }

        fn forcer(&mut self, cell: &Cell, threads: usize, rep: usize) -> Result<()> {
            let scratch = self.work.join("forcer-cold");
            if self.options.selected(Backend::Forcer, "cold") {
                let key = self.key(cell, FORCER, threads, "cold", rep);
                if !self.recorded.contains(&key) {
                    let base = scratch.join(format!("{}-{threads}-{rep}", cell.name));
                    remove_tree(&base)?;
                    fs::create_dir_all(&base)?;
                    let placement = Placement {
                        base: base.clone(),
                        cache: None,
                        read_only: false,
                    };
                    let label = format!("{}-forcer-cold-{threads}-{rep}", cell.name);
                    let measured = self.run_child(cell, FORCER, threads, &placement, &label)?;
                    self.record(key, cell, measured)?;
                    remove_tree(&base)?;
                }
            }
            if self.options.selected(Backend::Forcer, "warm") {
                let key = self.key(cell, FORCER, threads, "warm", rep);
                if !self.recorded.contains(&key) {
                    let placement = Placement {
                        base: self.warm_base(threads)?,
                        cache: None,
                        read_only: false,
                    };
                    let label = format!("{}-forcer-warm-{threads}-{rep}", cell.name);
                    let measured = self.run_child(cell, FORCER, threads, &placement, &label)?;
                    self.record(key, cell, measured)?;
                }
            }
            if self.options.selected(Backend::Forcer, "raw") {
                let key = self.key(cell, FORCER, threads, "raw", rep);
                if !self.recorded.contains(&key) {
                    let measured = self.raw(cell, threads, rep)?;
                    self.record(key, cell, measured)?;
                }
            }
            Ok(())
        }

        /// The generated FORM programs of `cell`, each run directly by one
        /// `form` or `tform -wN` process. `None` when the sum hit the cap.
        fn raw(&self, cell: &Cell, threads: usize, rep: usize) -> Result<Option<(Value, f64)>> {
            let base = self.warm_base(threads)?;
            let programs = self.work.join("programs").join(&cell.name);
            if !programs.exists() {
                let placement = Placement {
                    base: base.clone(),
                    cache: None,
                    read_only: false,
                };
                let label = format!("{}-forcer-programs-{threads}-{rep}", cell.name);
                if self
                    .run_child(cell, FORCER, threads, &placement, &label)?
                    .is_none()
                {
                    return Ok(None);
                }
            }
            let sources = fs::read_dir(&base)?
                .flatten()
                .find(|entry| {
                    entry
                        .file_name()
                        .to_string_lossy()
                        .starts_with("vakint-forcer-")
                })
                .ok_or_else(|| eyre::eyre!("no materialised Forcer sources"))?
                .path();
            let mut files = fs::read_dir(&programs)?
                .flatten()
                .map(|entry| entry.path())
                .collect::<Vec<_>>();
            files.sort();
            let tform = self.form.with_file_name("tform");
            let mut seconds = 0.0;
            for file in &files {
                let run = self.work.join("raw-run");
                remove_tree(&run)?;
                fs::create_dir_all(&run)?;
                for name in ["forcer.h", "forcer"] {
                    std::os::unix::fs::symlink(sources.join(name), run.join(name))?;
                }
                fs::copy(file, run.join("run_forcer.frm"))?;
                let mut command = if threads > 1 {
                    let mut command = self.pinned(threads, &tform);
                    command.arg(format!("-w{threads}"));
                    command
                } else {
                    self.pinned(threads, &self.form)
                };
                command
                    .arg("run_forcer.frm")
                    .current_dir(&run)
                    .stdin(Stdio::null())
                    .stdout(Stdio::null())
                    .stderr(Stdio::null());
                let Some((status, elapsed)) = self.run_capped(&mut command, seconds)? else {
                    return Ok(None);
                };
                seconds += elapsed;
                ensure!(status.success(), "raw FORM run failed: {}", file.display());
                remove_tree(&run)?;
            }
            Ok(Some((
                json!({
                    "evaluate_seconds": seconds,
                    "reduction_seconds": seconds,
                    "load_seconds": 0.0,
                    "form_runs": files.len(),
                }),
                seconds,
            )))
        }

        fn rustred(
            &mut self,
            cell: &Cell,
            arithmetic: &str,
            threads: usize,
            rep: usize,
        ) -> Result<()> {
            let backend = (Backend::RustRed, arithmetic);
            let cold_key = self.key(cell, backend, threads, "cold", rep);
            let warm_key = self.key(cell, backend, threads, "warm-same", rep);
            let want_cold = self.options.selected(Backend::RustRed, "cold")
                && !self.recorded.contains(&cold_key);
            let want_warm = self.options.selected(Backend::RustRed, "warm-same")
                && !self.recorded.contains(&warm_key);
            if !(want_cold || want_warm) {
                return Ok(());
            }
            let directory = self
                .work
                .join("rustred")
                .join(format!("{}-{arithmetic}-{threads}-{rep}", cell.name));
            let cache = directory.join("cache");
            remove_tree(&directory)?;
            fs::create_dir_all(&cache)?;
            let placement = Placement {
                base: directory.clone(),
                cache: Some(cache.clone()),
                read_only: false,
            };
            let label = format!("{}-rustred-{arithmetic}-cold-{threads}-{rep}", cell.name);
            let measured = self.run_child(cell, backend, threads, &placement, &label)?;
            let capped = measured.is_none();
            if want_cold {
                self.record(cold_key, cell, measured)?;
            }
            if want_warm && capped {
                let reason = "cold run capped, its cache is partial".to_owned();
                self.record_skipped(warm_key, cell, reason)?;
            } else if want_warm {
                let copy = directory.join("cache-read-only");
                let status = Command::new("cp")
                    .arg("-a")
                    .arg(&cache)
                    .arg(&copy)
                    .status()?;
                ensure!(status.success(), "copying {} failed", cache.display());
                let status = Command::new("chmod")
                    .args(["-R", "a-w"])
                    .arg(&copy)
                    .status()?;
                ensure!(status.success(), "chmod failed on {}", copy.display());
                let placement = Placement {
                    base: directory.clone(),
                    cache: Some(copy),
                    read_only: true,
                };
                let label = format!(
                    "{}-rustred-{arithmetic}-warm-same-{threads}-{rep}",
                    cell.name
                );
                let measured = self.run_child(cell, backend, threads, &placement, &label)?;
                self.record(warm_key, cell, measured)?;
            }
            remove_tree(&directory)
        }

        /// The warm-shared chain of one L: one cache from empty over the cells
        /// in table order. The ledger names the cells the cache holds and the
        /// cell in flight; a chain that cannot continue is rebuilt from empty,
        /// its recorded cells unmeasured. A capped cell leaves a partial cache,
        /// so it ends the chain: the cells after it are recorded as skipped,
        /// and it stays in flight for a later run under a longer cap.
        fn chain(
            &mut self,
            cells: &[Cell],
            arithmetic: &str,
            threads: usize,
            rep: usize,
        ) -> Result<()> {
            let backend = (Backend::RustRed, arithmetic);
            let measured = |harness: &Self, cell: &Cell| {
                harness
                    .recorded
                    .contains(&harness.key(cell, backend, threads, "warm-shared", rep))
            };
            if cells.iter().all(|cell| measured(self, cell)) {
                return Ok(());
            }
            let loops = cells[0].loops;
            let directory = self
                .work
                .join("chains")
                .join(format!("L{loops}-{arithmetic}-{threads}-{rep}"));
            let cache = directory.join("cache");
            let ledger_path = directory.join("chain.json");
            let mut ledger: Value = if ledger_path.exists() {
                serde_json::from_reader(File::open(&ledger_path)?)?
            } else {
                json!({"done": [], "in_flight": null})
            };
            let holds = |ledger: &Value, cell: &Cell| {
                ledger["done"]
                    .as_array()
                    .is_some_and(|done| done.iter().any(|name| name == &json!(cell.name)))
            };
            let continuable = ledger["in_flight"].is_null()
                && cache.exists()
                && cells
                    .iter()
                    .all(|cell| !measured(self, cell) || holds(&ledger, cell));
            if !continuable {
                remove_tree(&directory)?;
                fs::create_dir_all(&cache)?;
                ledger = json!({"done": [], "in_flight": null});
            }
            for (position, cell) in cells.iter().enumerate() {
                if holds(&ledger, cell) {
                    continue;
                }
                if self.stop_requested() {
                    return Ok(());
                }
                ledger["in_flight"] = json!(cell.name);
                serde_json::to_writer(File::create(&ledger_path)?, &ledger)?;
                let placement = Placement {
                    base: directory.clone(),
                    cache: Some(cache.clone()),
                    read_only: false,
                };
                let label = format!("{}-rustred-{arithmetic}-chain-{threads}-{rep}", cell.name);
                let result = self.run_child(cell, backend, threads, &placement, &label)?;
                let capped = result.is_none();
                if !measured(self, cell) {
                    let key = self.key(cell, backend, threads, "warm-shared", rep);
                    self.record(key, cell, result)?;
                }
                if capped {
                    for later in &cells[position + 1..] {
                        if !measured(self, later) {
                            let key = self.key(later, backend, threads, "warm-shared", rep);
                            let reason = format!("chain capped at {}", cell.name);
                            self.record_skipped(key, later, reason)?;
                        }
                    }
                    return Ok(());
                }
                ledger["in_flight"] = Value::Null;
                ledger["done"]
                    .as_array_mut()
                    .expect("the ledger holds a list")
                    .push(json!(cell.name));
                serde_json::to_writer(File::create(&ledger_path)?, &ledger)?;
            }
            Ok(())
        }

        fn run(&mut self, by_loops: &BTreeMap<usize, Vec<Cell>>) -> Result<bool> {
            // Forcer's keys carry no arithmetic: it is measured under the first.
            let arithmetics = self.options.arithmetics.clone();
            for arithmetic in &arithmetics {
                for &threads in &self.options.threads {
                    for rep in 0..self.options.repetitions {
                        for cells in by_loops.values() {
                            for cell in cells {
                                if self.stop_requested() {
                                    return Ok(false);
                                }
                                // Alternate the order of the backends between repetitions.
                                if rep % 2 == 0 {
                                    self.forcer(cell, threads, rep)?;
                                    self.rustred(cell, arithmetic, threads, rep)?;
                                } else {
                                    self.rustred(cell, arithmetic, threads, rep)?;
                                    self.forcer(cell, threads, rep)?;
                                }
                            }
                            if self.options.selected(Backend::RustRed, "warm-shared") {
                                self.chain(cells, arithmetic, threads, rep)?;
                            }
                        }
                    }
                }
            }
            Ok(!self.stop_requested())
        }
    }

    fn parent(options: &Options) -> Result<bool> {
        let build_profile = options.build_profile.clone().expect("checked by parse");
        if build_profile == "release" {
            ensure!(
                !cfg!(debug_assertions),
                "the release profile conflicts with enabled debug assertions"
            );
        }
        let crate_path = Path::new(env!("CARGO_MANIFEST_DIR"));
        let root = crate_path.join("../..").to_string_lossy().into_owned();
        // The last commit of the sources: committing results between two
        // sessions of a measurement does not make it another measurement.
        let sources = [".", ":(exclude)crates/vakint/docs/benchmark"];
        let gammaloop_sha = command_output(
            "git",
            &[
                &["-C", &root, "log", "-1", "--format=%H", "--"],
                &sources[..],
            ]
            .concat(),
        )?;
        let rustred = rustred_sha(crate_path)?;
        let dirty = !command_output(
            "git",
            &[&["-C", &root, "status", "--porcelain", "--"], &sources[..]].concat(),
        )?
        .is_empty();
        let source = Source {
            gammaloop_sha: &gammaloop_sha,
            rustred_sha: &rustred,
            build_profile: &build_profile,
            source_dirty: dirty,
        };
        results::check_build(source, BUILD_STAMP, options.child_program.is_some())?;
        let form = executable(&options.form)?;
        if options.backends.contains(&Backend::Forcer) && options.threads.iter().any(|&n| n > 1) {
            check_tform(&form)?;
        }
        let exe = match &options.child_program {
            Some(program) => program.clone(),
            None => std::env::current_exe()?,
        };
        let exact_sha256 = check_exact(&options.exact, source)?;
        let cpuinfo = fs::read_to_string("/proc/cpuinfo").unwrap_or_default();
        let status = fs::read_to_string("/proc/self/status").unwrap_or_default();
        let status_field = |name: &str| {
            status
                .lines()
                .find_map(|line| line.strip_prefix(name))
                .map(str::trim)
                .map(str::to_owned)
        };
        let binary_sha256 = sha256_of_file(&exe)?;
        let started_utc = utc_now()?;
        let work = options.work_dir.join(format!(
            "{build_profile}-{}-{}",
            &gammaloop_sha[..8],
            &rustred[..8]
        ));
        fs::create_dir_all(&work)?;
        let work = work.canonicalize()?;
        let protocol = json!({
            "gammaloop_sha": gammaloop_sha,
            "rustred_sha": rustred,
            "source_dirty": dirty,
            "build_profile": build_profile,
            "debug_assertions": cfg!(debug_assertions),
            "forcer_commit": FORCER_UPSTREAM_COMMIT,
            "form_version": command_output(&form.to_string_lossy(), &["-v"])?,
            "form_sha256": sha256_of_file(&form)?,
            "binary_sha256": binary_sha256,
            "build_stamp": BUILD_STAMP,
            "grid_sha256": format!("{:x}", Sha256::digest(GRID.as_bytes())),
            "physics_sha256": format!("{:x}", Sha256::digest(PHYSICS.as_bytes())),
            "exact_sha256": exact_sha256,
            "rustc": command_output("rustc", &["--version"]).ok(),
            "host": command_output("hostname", &[])?,
            "cpu_model": cpuinfo.lines().find_map(|line| {
                line.strip_prefix("model name")?.split_once(':').map(|(_, model)| model.trim())
            }),
            "cpus_allowed": status_field("Cpus_allowed_list:"),
            "mems_allowed": status_field("Mems_allowed_list:"),
            "cpus": options.cpus,
            "numactl": options.numa_node.map(|node| format!("--cpunodebind={node} --membind={node}")),
            "threads": options.threads,
            "thread_environment": options.threads.iter().map(|&threads| {
                (threads.to_string(), thread_environment(threads).into_iter()
                    .map(|(name, value)| (name.to_owned(), value))
                    .collect::<BTreeMap<_, _>>())
            }).collect::<BTreeMap<_, _>>(),
            "arithmetics": options.arithmetics,
            "ordered_weight_arithmetic": options.arithmetics.iter().map(|name| {
                let by_threads = options.threads.iter().map(|&threads| {
                    let settings = VakintSettings {
                        n_cores: threads,
                        ordered_weight_arithmetic: name.parse().expect("checked by parse"),
                        ..Default::default()
                    };
                    (threads.to_string(), settings.rustred_ordered_weight_arithmetic().to_string())
                }).collect::<BTreeMap<_, _>>();
                (name.clone(), by_threads)
            }).collect::<BTreeMap<_, _>>(),
            "adjoint_tail_env": std::env::var(ADJOINT_TAIL).ok(),
            "cells": options.cell_names,
            "repetitions": options.repetitions,
            "cap_seconds": options.cap_seconds,
            "work_filesystem": command_output(
                "stat", &["-f", "-c", "%T", &work.to_string_lossy()]).ok(),
            "invocation": std::env::args().collect::<Vec<_>>(),
            "scope": SCOPE,
            "started_utc": started_utc,
        });
        let by_loops = {
            let cells = all_cells()?;
            (1..=4)
                .map(|loops| (loops, options.cells(&cells, loops)))
                .filter(|(_, cells)| !cells.is_empty())
                .collect::<BTreeMap<_, _>>()
        };
        Vakint::initialize_vakint_symbols();
        let mut integrals = BTreeMap::new();
        for cell in by_loops.values().flatten() {
            let input = vakint_parse!(&cell.input).map_err(eyre::Report::msg)?;
            integrals.insert(cell.name.clone(), results::integral_count(&input));
        }
        let document = results::read(&options.output)?;
        let recorded = results::recorded(&document, source, options.cap_seconds)?;
        let outputs = results::outputs(&document, source);
        let run = format!("{started_utc}-{}", &binary_sha256[..12]);
        results::begin_run(&options.output, &run, &protocol, options.select_current)?;
        let output = options.output.canonicalize()?;
        let mut harness = Harness {
            options,
            form,
            run: run.clone(),
            output,
            executable: exe,
            recorded,
            outputs,
            integrals,
            work,
        };
        catch_stop_signals();
        let finished = harness.run(&by_loops)?;
        eprintln!(
            "{} run {run} in {}",
            if finished { "finished" } else { "stopped" },
            harness.output.display()
        );
        Ok(finished)
    }

    pub fn main() -> Result<()> {
        let Some(options) = Options::parse()? else {
            return Ok(());
        };
        if let Some(cell) = &options.child {
            return child(&options, cell);
        }
        if !parent(&options)? {
            std::process::exit(3);
        }
        Ok(())
    }
}
