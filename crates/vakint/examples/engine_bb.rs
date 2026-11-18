//! The numerator engine (`vakint::rustred_evaluation::p_contract`) on
//! background-field `B -> B` DOT graphs, gated and timed against the
//! reference path R2 (`vakint::refpath`), helper VH numerator lane engine.
//!
//! ```text
//! engine_bb <prepared dir with model.json + diagram-NNNNNN.dot> --out <dir>
//!     [--graphs 1,2,7-9 | --graph-list <file>]
//!     [--projection metric|transverse|longitudinal|both]   (default metric)
//!     [--reps <n>]           paired repetitions, each on fresh engines (default 1)
//!     [--ref]                run the reference R2 live, paired (gate + timing + values)
//!     [--ref-dir <dir>]      lane refpath's stored run (gate/<NNNNNN>.<projection>.gate, values/)
//!     [--per-branch]         also gate every colour branch against R2 on that branch
//!     [--ref-dump]           with --ref: write R2's gate views and values in
//!                            direct_bb's format (`<out>/ref/gate/<graph>.<projection>.gate`,
//!                            `<out>/ref/values/<graph>.txt`): reference dumps for later gates
//!     [--no-reduce]          engine: stop after the folded inputs (no value)
//!     [--order greedy|natural] [--no-prune] [--trace idenso|idenso-basis|naive]
//!     [--extra-terms <k>] [--p-library <dir>] [--warmup] [--steps]
//! ```
//! Writes `<out>/results.jsonl` (one object per graph), the engine's gate views
//! `<out>/gate/<graph>.<projection>.gate` and, for every gate mismatch,
//! `<out>/mismatch/<graph>.<projection>.txt`.
use std::collections::BTreeMap;
use std::io::Write as _;
use std::path::PathBuf;
use std::time::Instant;

use serde_json::json;
use symbolica::atom::{Atom, AtomCore};
use vakint::Vakint;
use vakint::engine_gate::{self, EngineRunOptions};
use vakint::refpath::{self, Projection, RunOptions};
use vakint::rustred_evaluation::p_contract::ContractionOrder;
use vakint::rustred_evaluation::p_contract::trace::TraceBackend;
use vakint::rustred_evaluation::p_gate::{self, GateDifference, GateView};

fn parse_graphs(text: &str) -> Vec<usize> {
    let mut out = Vec::new();
    for part in text.split(',').filter(|part| !part.is_empty()) {
        match part.split_once('-') {
            Some((a, b)) => out.extend(a.parse::<usize>().unwrap()..=b.parse::<usize>().unwrap()),
            None => out.push(part.parse().unwrap()),
        }
    }
    out
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

/// This thread's CPU seconds (user + system, nanosecond clock).
fn cpu_seconds() -> f64 {
    let mut spec = libc::timespec {
        tv_sec: 0,
        tv_nsec: 0,
    };
    // SAFETY: clock_gettime writes one timespec into a valid pointer.
    unsafe { libc::clock_gettime(libc::CLOCK_THREAD_CPUTIME_ID, &mut spec) };
    spec.tv_sec as f64 + spec.tv_nsec as f64 * 1e-9
}

/// User-space instructions retired by this thread (Linux `perf_event_open`,
/// `PERF_COUNT_HW_INSTRUCTIONS`, kernel and hypervisor excluded; lane prune's
/// counter); `None` when the counter cannot be opened.
struct InstructionCounter(i32);

impl InstructionCounter {
    fn open() -> Option<Self> {
        #[repr(C)]
        #[derive(Default)]
        struct Attributes {
            kind: u32,
            size: u32,
            config: u64,
            sample_period: u64,
            sample_type: u64,
            read_format: u64,
            flags: u64,
            wakeup_events: u32,
            bp_type: u32,
            config1: u64,
            config2: u64,
            branch_sample_type: u64,
            sample_regs_user: u64,
            sample_stack_user: u32,
            clockid: i32,
            sample_regs_intr: u64,
            aux_watermark: u32,
            sample_max_stack: u16,
            reserved2: u16,
            aux_sample_size: u32,
            reserved3: u32,
            sig_data: u64,
            config3: u64,
        }
        let attributes = Attributes {
            kind: 0,
            size: std::mem::size_of::<Attributes>() as u32,
            config: 1,
            flags: (1 << 5) | (1 << 6),
            ..Default::default()
        };
        // SAFETY: a valid attribute struct of the stated size; this thread,
        // any CPU, no group, no flags.
        let fd = unsafe {
            libc::syscall(
                libc::SYS_perf_event_open,
                &attributes as *const Attributes,
                0,
                -1,
                -1,
                0,
            )
        };
        (fd >= 0).then_some(Self(fd as i32))
    }

    fn read(&self) -> u64 {
        let mut value = 0_u64;
        // SAFETY: reads 8 bytes into a u64 from the counter's descriptor.
        let read = unsafe { libc::read(self.0, (&mut value as *mut u64).cast(), 8) };
        if read == 8 { value } else { 0 }
    }
}

impl Drop for InstructionCounter {
    fn drop(&mut self) {
        // SAFETY: closes the descriptor this value owns.
        unsafe {
            libc::close(self.0);
        }
    }
}

fn difference_json(difference: &GateDifference, keys: usize) -> serde_json::Value {
    json!({
        "pass": difference.is_empty(),
        "reference_keys": keys,
        "identical": difference.identical,
        "equal": difference.equal,
        "only_reference": difference.only_left.len(),
        "only_engine": difference.only_right.len(),
        "differ": difference.differ.len(),
    })
}

fn write_mismatch(
    out: &std::path::Path,
    graph: &str,
    tag: &str,
    difference: &GateDifference,
    reference: &GateView,
    engine: &GateView,
) -> eyre::Result<()> {
    let mut text = String::new();
    for key in &difference.only_left {
        text.push_str(&format!("only_reference {key:?} {}\n", reference[key]));
    }
    for key in &difference.only_right {
        text.push_str(&format!("only_engine {key:?} {}\n", engine[key]));
    }
    for (key, delta) in &difference.differ {
        text.push_str(&format!(
            "differ {key:?}\n  reference {}\n  engine    {}\n  delta     {delta}\n",
            reference[key], engine[key]
        ));
    }
    std::fs::create_dir_all(out.join("mismatch"))?;
    std::fs::write(
        out.join("mismatch").join(format!("{graph}.{tag}.txt")),
        text,
    )?;
    Ok(())
}

fn main() -> eyre::Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let mut directory = None;
    let mut out = None;
    let mut graphs = None;
    let mut options = EngineRunOptions::default();
    let mut reps = 1usize;
    let mut library = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4").map(PathBuf::from);
    let mut live = false;
    let mut ref_dir: Option<PathBuf> = None;
    let mut per_branch = false;
    let mut ref_dump = false;
    let mut warmup = false;
    let mut steps = false;
    let mut index = 0;
    while index < args.len() {
        let value = || args.get(index + 1).cloned().expect("missing value");
        match args[index].as_str() {
            "--out" => {
                out = Some(PathBuf::from(value()));
                index += 1;
            }
            "--graphs" => {
                graphs = Some(parse_graphs(&value()));
                index += 1;
            }
            "--graph-list" => {
                let text = std::fs::read_to_string(value())?;
                graphs = Some(
                    text.lines()
                        .filter(|l| !l.trim().is_empty())
                        .map(|l| l.trim().parse::<usize>())
                        .collect::<Result<Vec<_>, _>>()?,
                );
                index += 1;
            }
            "--projection" => {
                options.projection = value().parse::<Projection>().map_err(|e| eyre::eyre!(e))?;
                index += 1;
            }
            "--reps" => {
                reps = value().parse()?;
                index += 1;
            }
            "--extra-terms" => {
                options.extra_terms = value().parse()?;
                index += 1;
            }
            "--p-library" => {
                library = Some(PathBuf::from(value()));
                index += 1;
            }
            "--ref-dir" => {
                ref_dir = Some(PathBuf::from(value()));
                index += 1;
            }
            "--order" => {
                options.contract.order = match value().as_str() {
                    "greedy" => ContractionOrder::Greedy,
                    "natural" => ContractionOrder::Natural,
                    other => eyre::bail!("unknown order {other}"),
                };
                index += 1;
            }
            "--trace" => {
                options.contract.trace = match value().as_str() {
                    "idenso" => TraceBackend::Idenso,
                    "naive" => TraceBackend::Naive,
                    "idenso-basis" => TraceBackend::IdensoBasis,
                    other => eyre::bail!("unknown trace engine {other}"),
                };
                index += 1;
            }
            "--ref" => live = true,
            "--per-branch" => per_branch = true,
            "--ref-dump" => ref_dump = true,
            "--no-reduce" => options.reduce = false,
            "--no-prune" => options.contract.prune = false,
            "--warmup" => warmup = true,
            "--steps" => steps = true,
            other if directory.is_none() && !other.starts_with("--") => {
                directory = Some(PathBuf::from(other))
            }
            other => eyre::bail!("unknown argument {other}"),
        }
        index += 1;
    }
    let directory = directory.ok_or_else(|| eyre::eyre!("no prepared directory"))?;
    let out = out.ok_or_else(|| eyre::eyre!("no --out"))?;
    std::fs::create_dir_all(&out)?;
    let model = refpath::load_model(&directory)?;
    let mut files = std::fs::read_dir(&directory)?
        .filter_map(|entry| entry.ok().map(|entry| entry.path()))
        .filter(|path| {
            path.file_name()
                .and_then(|name| name.to_str())
                .is_some_and(|name| name.starts_with("diagram-") && name.ends_with(".dot"))
        })
        .collect::<Vec<_>>();
    files.sort();
    // Graph ids: the file stem without `diagram-` (`000003`, or
    // `000000.xid1.t000` for the xi-derivative files); `--graphs` selects
    // numeric ids, else positions in the sorted list.
    let id_of = |path: &PathBuf| -> String {
        let name = path.file_stem().unwrap().to_str().unwrap();
        name.trim_start_matches("diagram-").to_string()
    };
    let selected = match &graphs {
        Some(list) => {
            let by_number = files
                .iter()
                .filter_map(|p| id_of(p).parse::<usize>().ok().map(|n| (n, p.clone())))
                .collect::<BTreeMap<_, _>>();
            list.iter()
                .map(|n| {
                    by_number
                        .get(n)
                        .or_else(|| (by_number.is_empty()).then(|| files.get(*n)).flatten())
                        .cloned()
                        .ok_or_else(|| eyre::eyre!("no diagram {n}"))
                })
                .collect::<eyre::Result<Vec<_>>>()?
        }
        None => files.clone(),
    };
    let mut results = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(out.join("results.jsonl"))?;
    let host = std::fs::read_to_string("/proc/sys/kernel/hostname").unwrap_or_default();
    let mut template = Vakint::new()?;
    if let Some(library) = &library {
        template.load_p_library_directory(library, false)?;
    }
    let ref_options = RunOptions {
        projection: options.projection,
        extra_terms: options.extra_terms,
        count_terms: false,
        contract_only: false,
    };
    eprintln!(
        "[engine_bb] {} graphs of {} projection {} reps {reps} live-ref {live} ref-dir {:?} options {:?}",
        selected.len(),
        directory.display(),
        options.projection.name(),
        ref_dir,
        options.contract
    );
    if warmup && let Some(path) = selected.first() {
        let dot = std::fs::read_to_string(path)?;
        let engine = refpath::fresh_engine(&template)?;
        let _ = engine_gate::run_graph(&engine, &model, &dot, &options);
        if live {
            let engine = refpath::fresh_engine(&template)?;
            let _ = refpath::run_graph(&engine, &model, &dot, &ref_options);
        }
    }
    let (mut passed, mut failed, mut errors) = (0usize, 0usize, 0usize);
    let counter = InstructionCounter::open();
    eprintln!(
        "[engine_bb] instruction counter: {}",
        if counter.is_some() {
            "open"
        } else {
            "unavailable"
        }
    );
    for path in &selected {
        let dot = std::fs::read_to_string(path)?;
        let graph = id_of(path);
        let mut engine_runs = Vec::new();
        let mut ref_runs = Vec::new();
        let mut engine_wall = Vec::new();
        let mut engine_cpu = Vec::new();
        let mut ref_wall = Vec::new();
        let mut ref_cpu = Vec::new();
        let mut engine_rss = Vec::new();
        let mut engine_ins = Vec::new();
        let mut ref_ins = Vec::new();
        let mut ref_rss = Vec::new();
        let mut error = None;
        for rep in 0..reps {
            let order: &[bool] = if rep % 2 == 0 {
                &[true, false]
            } else {
                &[false, true]
            };
            for &engine_turn in order {
                if engine_turn {
                    let engine = refpath::fresh_engine(&template)?;
                    refpath::start_peak_rss();
                    let start_ins = counter.as_ref().map_or(0, InstructionCounter::read);
                    let mut marks = BTreeMap::<&'static str, u64>::new();
                    let (wall, cpu) = (Instant::now(), cpu_seconds());
                    let run = engine_gate::run_graph_marked(
                        &engine,
                        &model,
                        &dot,
                        &options,
                        &mut |label| {
                            if let Some(counter) = &counter {
                                let now = counter.read();
                                let entry = marks.entry(label).or_insert(0);
                                // Opening marks subtract, closing marks add: the
                                // sums are the instructions inside the stages.
                                if label.ends_with('d') {
                                    *entry += now
                                } else {
                                    *entry = entry.wrapping_sub(now)
                                }
                            }
                        },
                    );
                    let total_ins =
                        counter.as_ref().map_or(0, InstructionCounter::read) - start_ins;
                    let lower_ins = marks
                        .get("lower")
                        .copied()
                        .unwrap_or(0)
                        .wrapping_add(marks.get("lowered").copied().unwrap_or(0));
                    let reduce_ins = marks
                        .get("reduce")
                        .copied()
                        .unwrap_or(0)
                        .wrapping_add(marks.get("reduced").copied().unwrap_or(0));
                    engine_ins.push([total_ins, lower_ins, reduce_ins]);
                    engine_wall.push(wall.elapsed().as_secs_f64());
                    engine_cpu.push(cpu_seconds() - cpu);
                    engine_rss.push(refpath::peak_rss_bytes().unwrap_or(0));
                    match run {
                        Ok(run) => engine_runs.push(run),
                        Err(e) => error = Some(format!("engine: {e:?}")),
                    }
                } else if live {
                    let engine = refpath::fresh_engine(&template)?;
                    refpath::start_peak_rss();
                    let start_ins = counter.as_ref().map_or(0, InstructionCounter::read);
                    let (wall, cpu) = (Instant::now(), cpu_seconds());
                    let run = refpath::run_graph(&engine, &model, &dot, &ref_options);
                    ref_ins.push(counter.as_ref().map_or(0, InstructionCounter::read) - start_ins);
                    ref_wall.push(wall.elapsed().as_secs_f64());
                    ref_cpu.push(cpu_seconds() - cpu);
                    ref_rss.push(refpath::peak_rss_bytes().unwrap_or(0));
                    match run {
                        Ok(run) => ref_runs.push(run),
                        Err(e) => error = Some(format!("reference: {e:?}")),
                    }
                }
            }
            if error.is_some() {
                break;
            }
        }
        // The reference's instructions up to its contraction (parse, colour,
        // idenso contraction, normalisation), once: with its whole-run count
        // and the engine's reduce + finalise count this splits its run.
        let mut ref_contract_ins = None;
        if live
            && error.is_none()
            && let Some(counter) = &counter
        {
            let engine = refpath::fresh_engine(&template)?;
            let options = RunOptions {
                contract_only: true,
                ..ref_options.clone()
            };
            let start = counter.read();
            if refpath::run_graph(&engine, &model, &dot, &options).is_ok() {
                ref_contract_ins = Some(counter.read() - start);
            }
        }
        if let Some(error) = error {
            errors += 1;
            eprintln!("[engine_bb] diagram {graph}: ERROR {error}");
            writeln!(
                results,
                "{}",
                json!({"graph": graph, "file": path.display().to_string(), "error": error})
            )?;
            continue;
        }
        let run = &engine_runs[0];
        // The reference's gate views (zero-sector keys dropped), once per
        // projection: compared below and, with --ref-dump, written as dumps.
        let ref_views = ref_runs
            .first()
            .map(|reference| {
                reference
                    .results
                    .iter()
                    .map(refpath::gate_view)
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        if ref_dump && let Some(reference) = ref_runs.first() {
            let dump = out.join("ref");
            std::fs::create_dir_all(dump.join("gate"))?;
            std::fs::create_dir_all(dump.join("values"))?;
            let mut values_text = String::new();
            for (result, (view, count)) in reference.results.iter().zip(&ref_views) {
                // direct_bb's header and text form, so that these files are
                // byte-identical to what direct_bb writes for the graph.
                let header = BTreeMap::from([
                    ("graph".to_string(), reference.name.clone()),
                    ("file".to_string(), path.display().to_string()),
                    ("projection".to_string(), result.projection.to_string()),
                    ("loops".to_string(), reference.loops.to_string()),
                    ("zero_sectors_dropped".to_string(), count.zero.to_string()),
                    ("keys_before_drop".to_string(), count.total.to_string()),
                    ("census_unknown_kept".to_string(), count.unknown.to_string()),
                    (
                        "producer".to_string(),
                        "vakint refpath R2 (helper VH lane refpath)".to_string(),
                    ),
                ]);
                std::fs::write(
                    dump.join("gate")
                        .join(format!("{graph}.{}.gate", result.projection)),
                    p_gate::write_gate_view(view, &header),
                )?;
                values_text.push_str(&format!(
                    "{}\t{}\n",
                    result.projection,
                    result.value.to_canonical_string()
                ));
            }
            for (name, numerator) in &reference.numerators {
                values_text.push_str(&format!(
                    "numerator.{name}.bytes\t{}\n",
                    numerator.to_canonical_string().len()
                ));
            }
            std::fs::write(
                dump.join("values").join(format!("{graph}.txt")),
                values_text,
            )?;
        }
        let mut gates = Vec::new();
        let mut graph_pass = true;
        std::fs::create_dir_all(out.join("gate"))?;
        for (position, result) in run.results.iter().enumerate() {
            let header = BTreeMap::from([
                ("graph".to_string(), run.name.clone()),
                ("file".to_string(), path.display().to_string()),
                ("projection".to_string(), result.projection.to_string()),
                ("loops".to_string(), run.loops.to_string()),
                (
                    "zero_sectors_dropped".to_string(),
                    result.zero.zero.to_string(),
                ),
                (
                    "producer".to_string(),
                    "vakint p_contract engine (helper VH lane engine)".to_string(),
                ),
            ]);
            std::fs::write(
                out.join("gate")
                    .join(format!("{graph}.{}.gate", result.projection)),
                p_gate::write_gate_view(&result.view, &header),
            )?;
            let mut entry = json!({
                "projection": result.projection,
                "engine_keys": result.view.len(),
                "engine_zero_dropped": result.zero.zero,
                "channels": result.inputs.channels(),
                "entries": result.inputs.entries(),
                "reductions": result.reductions,
            });
            if let Some(reference) = ref_runs.first() {
                if let Some(ref_result) = reference.results.get(position) {
                    let (view, count) = &ref_views[position];
                    let difference = p_gate::compare_gate_views(view, &result.view);
                    if !difference.is_empty() {
                        graph_pass = false;
                        write_mismatch(
                            &out,
                            &graph,
                            &format!("{}.live", result.projection),
                            &difference,
                            view,
                            &result.view,
                        )?;
                    }
                    entry["live"] = difference_json(&difference, view.len());
                    entry["ref_keys_all"] = json!(count.total);
                    entry["ref_zero_dropped"] = json!(count.zero);
                    entry["ref_reductions"] = json!(ref_result.observation.reductions);
                    if let Some(value) = &result.value {
                        let equal = (value - &ref_result.value).expand().is_zero();
                        graph_pass &= equal;
                        entry["value_equal_live"] = json!(equal);
                    }
                } else if reference.zero_reason.is_none() {
                    // The reference skips a numerator that sums to zero over
                    // the colour branches: the engine's inputs must be empty.
                    let empty = result.view.is_empty();
                    graph_pass &= empty;
                    entry["live"] = json!({"pass": empty, "reference": "zero numerator (not evaluated)", "engine_keys": result.view.len()});
                }
            }
            if let Some(dir) = &ref_dir {
                let file = dir
                    .join("gate")
                    .join(format!("{graph}.{}.gate", result.projection));
                match std::fs::read_to_string(&file) {
                    Ok(text) => {
                        let (view, _) =
                            p_gate::read_gate_view(&text).map_err(|e| eyre::eyre!(e))?;
                        let difference = p_gate::compare_gate_views(&view, &result.view);
                        if !difference.is_empty() {
                            graph_pass = false;
                            write_mismatch(
                                &out,
                                &graph,
                                &format!("{}.stored", result.projection),
                                &difference,
                                &view,
                                &result.view,
                            )?;
                        }
                        entry["stored"] = difference_json(&difference, view.len());
                    }
                    Err(_) => {
                        // lane refpath writes no dump for a zero numerator.
                        let empty = result.view.is_empty();
                        graph_pass &= empty;
                        entry["stored"] = json!({"pass": empty, "reference": "no stored dump", "engine_keys": result.view.len()});
                    }
                }
                let values = dir.join("values").join(format!("{graph}.txt"));
                if let (Some(value), Ok(text)) = (&result.value, std::fs::read_to_string(&values)) {
                    for line in text.lines() {
                        if let Some((name, stored)) = line.split_once('\t')
                            && name == result.projection
                        {
                            let stored =
                                symbolica::try_parse!(stored).map_err(|e| eyre::eyre!("{e}"))?;
                            let equal = (value - &stored).expand().is_zero();
                            graph_pass &= equal;
                            entry["value_equal_stored"] = json!(equal);
                        }
                    }
                }
            }
            let stats = &result.inputs.stats;
            entry["engine_stats"] = json!({
                "factors": stats.factors, "renames": stats.renames, "trace_lengths": stats.trace_lengths,
                "peak_terms": stats.peak_terms, "final_terms": stats.final_terms, "formed": stats.formed,
                "dropped": stats.dropped, "emitted": stats.emitted, "oracle": stats.oracle, "ring": stats.ring,
                "values": result.inputs.values, "zero_values": result.inputs.zeros,
                "seconds": {"match": stats.match_seconds, "prepare": stats.prepare_seconds, "trace": stats.trace_seconds,
                            "contract": stats.contract_seconds, "fold": stats.fold_seconds},
            });
            if steps {
                entry["steps"] = json!(
                    stats
                        .steps
                        .iter()
                        .map(|s| json!([
                            s.factor,
                            s.factor_terms,
                            s.terms,
                            s.open,
                            s.formed,
                            s.dropped
                        ]))
                        .collect::<Vec<_>>()
                );
            }
            gates.push(entry);
        }
        if let Some(reference) = ref_runs.first()
            && reference.zero_reason.is_some() != run.zero_reason.is_some()
        {
            graph_pass = false;
        }
        let mut branch_gates = Vec::new();
        if per_branch && run.zero_reason.is_none() {
            let engine = refpath::fresh_engine(&template)?;
            for (projection, branch, difference, keys) in engine_gate::branch_gate(
                &engine,
                &model,
                &dot,
                options.projection,
                &options.contract,
                options.extra_terms,
            )? {
                if !difference.is_empty() {
                    graph_pass = false;
                    let tag = format!("{projection}.branch{branch}");
                    std::fs::create_dir_all(out.join("mismatch"))?;
                    std::fs::write(
                        out.join("mismatch").join(format!("{graph}.{tag}.txt")),
                        format!("{difference:?}\n"),
                    )?;
                }
                branch_gates.push(json!({"projection": projection, "branch": branch, "gate": difference_json(&difference, keys)}));
            }
        }
        if graph_pass {
            passed += 1;
        } else {
            failed += 1;
        }
        let e = |get: fn(&engine_gate::EngineSeconds) -> f64| {
            median(engine_runs.iter().map(|r| get(&r.seconds)).collect())
        };
        let r = |get: fn(&refpath::StageSeconds) -> f64| {
            median(ref_runs.iter().map(|r| get(&r.seconds)).collect())
        };
        let record = json!({
            "graph": graph,
            "file": path.display().to_string(),
            "name": run.name,
            "loops": run.loops,
            "branches": run.branches,
            "zero_reason": run.zero_reason,
            "pass": graph_pass,
            "gates": gates,
            "branch_gates": branch_gates,
            "reps": reps,
            "engine_seconds": {
                "parse": e(|s| s.parse), "colour": e(|s| s.colour), "match": e(|s| s.matching),
                "lower": e(|s| s.lower), "trace": e(|s| s.trace), "reduce": e(|s| s.reduce),
                "prepare": e(|s| s.prepare), "library": e(|s| s.library), "finalise": e(|s| s.finalise),
                "wall": median(engine_wall.clone()), "cpu": median(engine_cpu.clone()),
            },
            "engine_reps": engine_runs.iter().map(|r| json!([r.seconds.parse, r.seconds.colour, r.seconds.matching, r.seconds.lower, r.seconds.trace, r.seconds.reduce, r.seconds.finalise])).collect::<Vec<_>>(),
            "engine_wall_reps": engine_wall, "engine_cpu_reps": engine_cpu,
            "ref_seconds": if ref_runs.is_empty() { json!(null) } else { json!({
                "parse": r(|s| s.parse), "colour": r(|s| s.colour), "contract": r(|s| s.contract),
                "lower": r(|s| s.lower), "reduce": r(|s| s.reduce), "prepare": r(|s| s.prepare),
                "library": r(|s| s.library), "finalise": r(|s| s.finalise),
                "wall": median(ref_wall.clone()), "cpu": median(ref_cpu.clone()),
            }) },
            "ref_reps": ref_runs.iter().map(|r| json!([r.seconds.parse, r.seconds.colour, r.seconds.contract, r.seconds.lower, r.seconds.reduce, r.seconds.finalise])).collect::<Vec<_>>(),
            "ref_wall_reps": ref_wall, "ref_cpu_reps": ref_cpu,
            "engine_peak_rss": engine_rss.iter().max(), "ref_peak_rss": ref_rss.iter().max(),
            "engine_instructions_reps": engine_ins, "ref_instructions_reps": ref_ins,
            "ref_contract_only_instructions": ref_contract_ins,
            "host": host.trim(),
        });
        writeln!(results, "{record}")?;
        results.flush()?;
        eprintln!(
            "[engine_bb] diagram {graph} L{} branches {} {} | engine match {:.4} lower {:.4} (trace {:.4}) reduce {:.4} | ref contract {:.4} lower {:.4} reduce {:.4} | {}",
            run.loops,
            run.branches,
            if graph_pass { "PASS" } else { "FAIL" },
            e(|s| s.matching),
            e(|s| s.lower),
            e(|s| s.trace),
            e(|s| s.reduce),
            r(|s| s.contract),
            r(|s| s.lower),
            r(|s| s.reduce),
            run.zero_reason.as_deref().unwrap_or("")
        );
        let _ = Atom::Zero;
    }
    eprintln!("[engine_bb] done: {passed} pass, {failed} fail, {errors} errors");
    Ok(())
}
