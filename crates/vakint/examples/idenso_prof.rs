//! Per-stage profile of idenso's work on background-field `B -> B` DOT graphs
//! (helper VH, numerator assignment, lane idenso-audit).
//!
//! For every graph: the reference chain R2's whole-numerator contraction
//! (`vakint::refpath::contract`, replayed stage by stage: product of the
//! residuals, routing, `mink(4)` -> `mink(D)`, cook, projector,
//! `infer_indexed`, `simplify_algebra{gamma: canonical(), epsilon: true}`,
//! `to_dots`, uncook, couplings), checked equal to `refpath::contract` on the
//! same branches; R2's colour split (`split_local_colour` and the closed
//! network reductions `colour_coefficient` separately); the engine's colour
//! stage (`FactorizedNumerator::from_feynkit`, per-channel and merged).
//!
//! ```text
//! idenso_prof <prepared dir with model.json + diagram-NNNNNN.dot> --out <file.jsonl>
//!     [--graphs 1,2,7-9 | --graph-list <file>] [--reps <n>] (default 3)
//!     [--projection metric|both] (default metric) [--phases] [--no-engine-colour]
//!     [--evaluate [--p-library <dir>]]  also refpath::run_graph end to end (fresh engine per rep;
//!                          set VAKINT_P_PHASES=1 for vakint's front-end phases on stderr, between
//!                          the markers `[idenso_prof] EVAL graph G rep r` and `EVAL-END`)
//!     [--direct]           with --evaluate: also R3 (`refpath::direct`: edge-basis contraction, denominator-first
//!                          match, direct read without expansion), paired with R2, gated against it
//!     [--r3-routed]        R3 contracts the ROUTED product exactly as R2 (loop basis); only the hand-over differs
//!     [--r3-compact]       R3 routes with the spanning tree of the fewest loop + root terms (exact enumeration)
//!     [--r3-closing]       R3 with idenso's frontier-level closing of quark loops (DirectOptions::closing)
//!     [--r3-order <o>]     R3 with idenso's frontier order (DirectOptions::frontier_order): boundary | close | merge | merge:<cap>
//!     [--r3-dots]          R3 asks idenso for dot notation directly (AlgebraContraction::Dots, no to_dots)
//!     [--r3-prune]         R3's reader prunes zero sectors early (DirectOptions::prune_early; the library default)
//!     [--edge-contract]    idenso alone: R3's edge-basis contraction timed next to R2's (no vakint; L5)
//!     [--ref-dir <dir>]    with --direct: also gate R3 against lane refpath's stored dumps (gate/<NNNNNN>.<projection>.gate)
//!     [--no-staged]        skip the stage-by-stage replay of R2's contraction and the colour split timing
//!     [--stem-list <file>]  inputs by file stem (`diagram-<stem>.dot`, dumps `<stem>.<projection>.gate`), e.g. the
//!                          xi-derivative sets `000000.xid1.t000`; the record's `graph` is the stem's leading number
//!     [--merge-colour]     R3 with merged colour terms (refpath::direct::colour_terms)
//!     [--dump-dir <dir>]   per graph and branch: `<graph>.<b>.bin` = the `infer_indexed` input
//!                          (Symbolica binary export with its symbol state), `<graph>.<b>.out` = the
//!                          `to_dots` result (canonical string): inputs of idenso's own replay driver
//! ```
//! `--phases` (needs the vakint feature `idenso-phase-timing`): one extra
//! repetition with idenso's own phase clocks (`idenso::reference_cases::timing`)
//! around `simplify_algebra` and `to_dots`.
use std::io::Write as _;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use feynkit_graph::FeynmanDiagram;
use feynkit_kinematics::Kinematics;
use idenso::dirac::GammaSimplifySettings;
use idenso::tensor::{AlgebraSettings, SymbolicTensor};
use idenso::{CookMode, CookSettings, IndexTooling, W_};
use serde_json::json;
use spenso::network::library::symbolic::ETS;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::representation::{LibrarySlot, Minkowski, RepName};
use spenso::structure::slot::IsAbstractSlot;
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder};
use symbolica::symbol;
use vakint::p_integral::factorized::{
    ColourMode, FactorizedNumerator, FromFeynkitOptions, colour_coefficient, split_local_colour,
};
use vakint::refpath::{self, Projection};

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

fn cpu_seconds() -> f64 {
    let mut spec = libc::timespec {
        tv_sec: 0,
        tv_nsec: 0,
    };
    // SAFETY: clock_gettime writes one timespec into a valid pointer.
    unsafe { libc::clock_gettime(libc::CLOCK_THREAD_CPUTIME_ID, &mut spec) };
    spec.tv_sec as f64 + spec.tv_nsec as f64 * 1e-9
}

/// User-space instructions retired by this thread (as engine_bb).
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
        // SAFETY: valid attribute struct; this thread, any CPU, no group.
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

/// Stage clocks: wall seconds, thread CPU seconds, user instructions.
const STAGES: [&str; 13] = [
    "product",
    "route",
    "to_d",
    "cook",
    "projector",
    "infer_indexed",
    "simplify_algebra",
    "to_dots",
    "uncook",
    "check_tags",
    "couplings",
    "colour_split",
    "colour_reduce",
];

#[derive(Clone, Default)]
struct Clocks {
    wall: [f64; STAGES.len()],
    cpu: [f64; STAGES.len()],
    instructions: [u64; STAGES.len()],
}

struct Timer<'a> {
    counter: Option<&'a InstructionCounter>,
    wall: Instant,
    cpu: f64,
    instructions: u64,
}

impl<'a> Timer<'a> {
    fn start(counter: Option<&'a InstructionCounter>) -> Self {
        Self {
            counter,
            instructions: counter.map_or(0, InstructionCounter::read),
            cpu: cpu_seconds(),
            wall: Instant::now(),
        }
    }

    fn stop(self, clocks: &mut Clocks, stage: usize) {
        clocks.wall[stage] += self.wall.elapsed().as_secs_f64();
        clocks.cpu[stage] += cpu_seconds() - self.cpu;
        clocks.instructions[stage] += self
            .counter
            .map_or(0, |c| c.read().saturating_sub(self.instructions));
    }
}

fn indexed(momentum: &Atom, slot: &Atom) -> Atom {
    match momentum.as_view() {
        AtomView::Var(variable) => variable.get_symbol().call(slot),
        AtomView::Fun(fun) => FunctionBuilder::new(fun.get_symbol())
            .add_args(fun.iter())
            .add_arg(slot)
            .finish(),
        _ => panic!("expected a momentum head, got {momentum}"),
    }
}

fn terms_of(atom: &Atom) -> usize {
    terms(atom)
}

fn byte_size(atom: &Atom) -> usize {
    atom.as_view().get_byte_size()
}

fn terms(atom: &Atom) -> usize {
    match atom.as_view() {
        AtomView::Add(add) => add.get_nargs(),
        view if view.is_zero() => 0,
        _ => 1,
    }
}

#[derive(Default)]
struct Sizes {
    cooked_bytes: usize,
    simplified_bytes: usize,
    scalar_terms: usize,
}

#[cfg(feature = "idenso-phase-timing")]
type PhaseLog = Vec<String>;
#[cfg(not(feature = "idenso-phase-timing"))]
type PhaseLog = ();

/// R2's `contract` (refpath.rs), replayed with stage clocks. Returns the
/// projections' scalar numerators (before the overall weight).
#[allow(clippy::too_many_arguments)]
fn contract_staged(
    diagram: &FeynmanDiagram,
    branches: &[refpath::ColourBranch],
    momentum: &Atom,
    q_squared: &Atom,
    metric_only: bool,
    clocks: &mut Clocks,
    counter: Option<&InstructionCounter>,
    sizes: &mut Sizes,
    dump: Option<&mut Vec<(Vec<u8>, String)>>,
    phases: Option<&mut PhaseLog>,
) -> eyre::Result<Vec<(&'static str, Atom)>> {
    let s = &*refpath::SYMBOLS;
    let d = Atom::var(s.d);
    let mut g_sum = Atom::Zero;
    let mut qq_sum = Atom::Zero;
    let basis = diagram.loop_momentum_basis();
    let mut dump = dump;
    #[allow(unused_mut, unused_variables)]
    let mut phases = phases;
    for branch in branches {
        let t = Timer::start(counter);
        let product = branch
            .residuals
            .iter()
            .fold(Atom::num(1), |product, residual| product * residual);
        t.stop(clocks, 0);
        let t = Timer::start(counter);
        let routed = basis.route_expression(&product);
        t.stop(clocks, 1);
        let t = Timer::start(counter);
        let continued = routed
            .replace(Minkowski {}.to_symbolic([Atom::var(W_.d_), Atom::var(W_.a___)]))
            .with(Minkowski {}.to_symbolic([d.clone(), Atom::var(W_.a___)]));
        t.stop(clocks, 2);

        let t = Timer::start(counter);
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let cooked = cooking.try_cook(continued.as_view())?;
        t.stop(clocks, 3);
        let t = Timer::start(counter);
        let ports = cooked
            .list_dangling::<AbstractIndex>()?
            .iter()
            .map(|port| LibrarySlot::<AbstractIndex>::try_from(port.as_view()))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| eyre::eyre!("{e:?}"))?;
        eyre::ensure!(
            ports.len() == 2
                && ports
                    .iter()
                    .all(|slot| slot.rep_name() == (Minkowski {}).into()),
            "two dangling Lorentz ports expected"
        );
        let (left, right) = (ports[0].to_atom(), ports[1].to_atom());
        let metric = ETS.metric(&left, &right);
        let projector = if metric_only {
            metric
        } else {
            Atom::var(s.metric_tag) * metric
                + Atom::var(s.qq_tag) * indexed(momentum, &left) * indexed(momentum, &right)
        };
        let input = &cooked * &projector;
        t.stop(clocks, 4);
        sizes.cooked_bytes += byte_size(&input);
        let mut exported = Vec::new();
        if dump.is_some() {
            input.as_view().export(&mut exported)?;
        }

        let t = Timer::start(counter);
        let inferred = SymbolicTensor::infer_indexed(input)?;
        t.stop(clocks, 5);
        let settings = AlgebraSettings {
            gamma: Some(GammaSimplifySettings::canonical()),
            epsilon: true,
            ..Default::default()
        };
        let t = Timer::start(counter);
        #[cfg(feature = "idenso-phase-timing")]
        let simplified = if let Some(log) = phases.as_deref_mut() {
            let (result, measurement) =
                idenso::reference_cases::timing::measure("simplify_algebra", || {
                    inferred.simplify_algebra(&settings)
                });
            log.push(measurement.to_string());
            result?
        } else {
            inferred.simplify_algebra(&settings)?
        };
        #[cfg(not(feature = "idenso-phase-timing"))]
        let simplified = inferred.simplify_algebra(&settings)?;
        t.stop(clocks, 6);
        sizes.simplified_bytes += byte_size(simplified.expression());
        let t = Timer::start(counter);
        #[cfg(feature = "idenso-phase-timing")]
        let reduced = if let Some(log) = phases.as_deref_mut() {
            let (result, measurement) =
                idenso::reference_cases::timing::measure("to_dots", || simplified.to_dots());
            log.push(measurement.to_string());
            result?
        } else {
            simplified.to_dots()?
        };
        #[cfg(not(feature = "idenso-phase-timing"))]
        let reduced = simplified.to_dots()?;
        t.stop(clocks, 7);
        let t = Timer::start(counter);
        let scalar = cooking.uncook(reduced.expression().as_view());
        t.stop(clocks, 8);
        if let Some(dump) = dump.as_deref_mut() {
            dump.push((exported, reduced.expression().to_canonical_string()));
        }
        sizes.scalar_terms += terms(&scalar);
        let t = Timer::start(counter);
        eyre::ensure!(
            scalar
                .list_dangling::<AbstractIndex>()
                .is_ok_and(|ports| ports.is_empty()),
            "the contraction kept open indices"
        );
        let (g, qq) = if metric_only {
            (scalar, Atom::Zero)
        } else {
            (
                scalar
                    .replace(Atom::var(s.qq_tag))
                    .with(Atom::Zero)
                    .replace(Atom::var(s.metric_tag))
                    .with(Atom::num(1)),
                scalar
                    .replace(Atom::var(s.metric_tag))
                    .with(Atom::Zero)
                    .replace(Atom::var(s.qq_tag))
                    .with(Atom::num(1)),
            )
        };
        g_sum += &branch.coefficient * g;
        qq_sum += &branch.coefficient * qq;
        t.stop(clocks, 9);
    }
    let t = Timer::start(counter);
    let couplings = |value: Atom| {
        diagram
            .model()
            .expand_couplings(&value)
            .replace(Atom::var(symbol!("UFO::gs")))
            .with(Atom::num(1))
    };
    let i_q2 = Atom::i() * q_squared;
    let mut out = Vec::new();
    if metric_only {
        out.push(("metric", couplings(&g_sum / (&d - Atom::num(1)) / &i_q2)));
    } else {
        out.push((
            "transverse",
            couplings((&g_sum - &qq_sum / q_squared) / (&d - Atom::num(1)) / &i_q2),
        ));
        out.push(("longitudinal", couplings(&qq_sum / q_squared / &i_q2)));
    }
    t.stop(clocks, 10);
    Ok(out)
}

/// R2's colour split with the network reductions timed apart.
fn colour_staged(
    diagram: &FeynmanDiagram,
    clocks: &mut Clocks,
    counter: Option<&InstructionCounter>,
) -> eyre::Result<usize> {
    use feynkit_graph::NumeratorInterface;
    let t = Timer::start(counter);
    let mut locals = Vec::new();
    for (_, vertex) in diagram.vertices() {
        locals.push(split_local_colour(
            &NumeratorInterface::export(&vertex.numerator).0,
        )?);
    }
    for (_, _, edge) in diagram.edges() {
        locals.push(split_local_colour(
            &NumeratorInterface::export(&edge.numerator).0,
        )?);
    }
    let mut partial: Vec<Atom> = vec![Atom::num(1)];
    for rows in &locals {
        let mut next = Vec::with_capacity(partial.len() * rows.len());
        for colour in &partial {
            for row in rows {
                next.push(colour * &row.colour);
            }
        }
        partial = next;
    }
    t.stop(clocks, 11);
    let t = Timer::start(counter);
    let mut nonzero = 0;
    for colour in &partial {
        let coefficient = colour_coefficient(colour, ColourMode::Invariants)?;
        if !coefficient.is_zero() {
            nonzero += 1;
        }
    }
    t.stop(clocks, 12);
    Ok(nonzero)
}

fn main() -> eyre::Result<()> {
    let args = std::env::args().collect::<Vec<_>>();
    let directory = PathBuf::from(&args[1]);
    let mut out_path = None;
    let mut graphs = None;
    let mut reps = 3usize;
    let mut projection = Projection::Metric;
    let mut want_phases = false;
    let mut engine_colour = true;
    let mut dump_dir = None;
    let mut evaluate = false;
    let mut direct = false;
    let mut merge_colour = false;
    let mut r3_dots = false;
    let mut r3_closing = false;
    let mut r3_order = idenso::tensor::FrontierOrder::Boundary;
    let mut r3_routed = false;
    let mut r3_compact = false;
    let mut r3_prune = false;
    let mut edge_contract = false;
    let mut ref_dir: Option<PathBuf> = None;
    let mut no_staged = false;
    let mut stems: Option<Vec<String>> = None;
    let mut library = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4").map(PathBuf::from);
    let mut i = 2;
    while i < args.len() {
        match args[i].as_str() {
            "--out" => {
                out_path = Some(PathBuf::from(&args[i + 1]));
                i += 1;
            }
            "--graphs" => {
                graphs = Some(parse_graphs(&args[i + 1]));
                i += 1;
            }
            "--graph-list" => {
                let text = std::fs::read_to_string(&args[i + 1])?;
                graphs = Some(
                    text.split_whitespace()
                        .map(|t| t.parse().unwrap())
                        .collect(),
                );
                i += 1;
            }
            "--reps" => {
                reps = args[i + 1].parse()?;
                i += 1;
            }
            "--projection" => {
                projection = args[i + 1].parse().map_err(|e: String| eyre::eyre!(e))?;
                i += 1;
            }
            "--phases" => want_phases = true,
            "--no-engine-colour" => engine_colour = false,
            "--evaluate" => evaluate = true,
            "--direct" => direct = true,
            "--merge-colour" => merge_colour = true,
            "--r3-dots" => r3_dots = true,
            "--r3-closing" => r3_closing = true,
            "--r3-order" => {
                r3_order = match args[i + 1].as_str() {
                    "boundary" => idenso::tensor::FrontierOrder::Boundary,
                    "close" => idenso::tensor::FrontierOrder::Close,
                    "merge" => idenso::tensor::FrontierOrder::Merge { cap: None },
                    other => idenso::tensor::FrontierOrder::Merge {
                        cap: Some(
                            other
                                .strip_prefix("merge:")
                                .and_then(|cap| cap.parse().ok())
                                .unwrap_or_else(|| panic!("unknown --r3-order {other}")),
                        ),
                    },
                };
                i += 1;
            }
            "--r3-routed" => r3_routed = true,
            "--r3-compact" => r3_compact = true,
            "--r3-prune" => r3_prune = true,
            "--edge-contract" => edge_contract = true,
            "--no-staged" => no_staged = true,
            "--stem-list" => {
                let text = std::fs::read_to_string(&args[i + 1])?;
                stems = Some(text.split_whitespace().map(str::to_owned).collect());
                i += 1;
            }
            "--ref-dir" => {
                ref_dir = Some(PathBuf::from(&args[i + 1]));
                i += 1;
            }
            "--p-library" => {
                library = Some(PathBuf::from(&args[i + 1]));
                i += 1;
            }
            "--dump-dir" => {
                dump_dir = Some(PathBuf::from(&args[i + 1]));
                i += 1;
            }
            other => eyre::bail!("unknown option {other}"),
        }
        i += 1;
    }
    #[cfg(not(feature = "idenso-phase-timing"))]
    if want_phases {
        eyre::bail!("--phases needs the vakint feature idenso-phase-timing");
    }
    let out_path = out_path.ok_or_else(|| eyre::eyre!("--out <file> is required"))?;
    if let Some(dir) = &dump_dir {
        std::fs::create_dir_all(dir)?;
    }
    vakint::Vakint::initialize_vakint_symbols();
    let model = refpath::load_model(&directory)?;
    let graphs = match graphs {
        Some(list) => list,
        None => {
            let mut list = std::fs::read_dir(&directory)?
                .filter_map(|e| e.ok())
                .filter_map(|e| {
                    let name = e.file_name().to_string_lossy().to_string();
                    name.strip_prefix("diagram-")?
                        .strip_suffix(".dot")?
                        .parse::<usize>()
                        .ok()
                })
                .collect::<Vec<_>>();
            list.sort();
            list
        }
    };
    let counter = InstructionCounter::open();
    let template = if evaluate {
        let started = Instant::now();
        let mut template = vakint::Vakint::new()?;
        if let Some(library) = &library {
            template.load_p_library_directory(library, false)?;
        }
        // warm-up (process-wide lazy statics) on a fresh engine, untimed
        eprintln!(
            "[idenso_prof] p library loaded in {:.3} s",
            started.elapsed().as_secs_f64()
        );
        Some(template)
    } else {
        None
    };
    let mut out = std::io::BufWriter::new(std::fs::File::create(&out_path)?);
    let metric_only = projection == Projection::Metric;
    let stems = stems.unwrap_or_else(|| graphs.iter().map(|g| format!("{g:06}")).collect());
    for stem in stems {
        // the graph number is the stem's leading digits; files and dumps are named by the stem
        let graph: usize = stem
            .split('.')
            .next()
            .and_then(|digits| digits.parse().ok())
            .ok_or_else(|| eyre::eyre!("a stem without a leading graph number: {stem}"))?;
        let file = directory.join(format!("diagram-{stem}.dot"));
        let dot = std::fs::read_to_string(&file)?;
        refpath::start_peak_rss();
        let started = Instant::now();
        let diagram = FeynmanDiagram::from_dot(Arc::clone(&model), &dot)?;
        let parse = started.elapsed().as_secs_f64();
        if refpath::scaleless_contact(&diagram) {
            writeln!(
                out,
                "{}",
                json!({"graph": graph, "stem": stem, "zero": "scaleless"})
            )?;
            continue;
        }
        let (branches, total) = refpath::colour_branches(&diagram)?;
        if branches.is_empty() {
            writeln!(
                out,
                "{}",
                json!({"graph": graph, "stem": stem, "zero": "colour", "branches_total": total})
            )?;
            continue;
        }
        let dimension = Atom::var(refpath::SYMBOLS.d);
        let kinematics = Kinematics::in_dimension(&dimension)?;
        let family = diagram.propagator_family(&kinematics)?;
        let [momentum] = family.external_momenta() else {
            eyre::bail!("one external momentum expected");
        };
        let momentum = momentum.clone();
        let q_squared = Kinematics::in_dimension(&dimension)?
            .with_momenta([momentum.clone()])?
            .scalar_product(&momentum, &momentum)?;

        // paired: R2's own contract (whole, timed) and the staged replay, alternating.
        let mut whole = Vec::new();
        let mut staged = Vec::new();
        let mut reference = None;
        let mut reference_terms: Option<usize> = None;
        let mut sizes = Sizes::default();
        let mut equal = true;
        let staged_reps = if no_staged { 0 } else { reps };
        for rep in 0..staged_reps {
            let run_whole =
                |whole: &mut Vec<(f64, f64)>| -> eyre::Result<Vec<(&'static str, Atom)>> {
                    let wall = Instant::now();
                    let cpu = cpu_seconds();
                    let value =
                        refpath::contract(&diagram, &branches, &momentum, &q_squared, projection)?;
                    whole.push((wall.elapsed().as_secs_f64(), cpu_seconds() - cpu));
                    Ok(value)
                };
            let mut clocks = Clocks::default();
            let mut rep_sizes = Sizes::default();
            let mut dump = (rep == 0 && dump_dir.is_some()).then(Vec::new);
            let mut run_staged = |clocks: &mut Clocks| {
                contract_staged(
                    &diagram,
                    &branches,
                    &momentum,
                    &q_squared,
                    metric_only,
                    clocks,
                    counter.as_ref(),
                    &mut rep_sizes,
                    dump.as_mut(),
                    None,
                )
            };
            let (value_whole, value_staged) = if rep % 2 == 0 {
                let a = run_whole(&mut whole)?;
                let b = run_staged(&mut clocks)?;
                (a, b)
            } else {
                let b = run_staged(&mut clocks)?;
                let a = run_whole(&mut whole)?;
                (a, b)
            };
            colour_staged(&diagram, &mut clocks, counter.as_ref())?;
            if let (Some(dir), Some(dump)) = (&dump_dir, &dump) {
                for (b, (input, output)) in dump.iter().enumerate() {
                    std::fs::write(dir.join(format!("{stem}.{b}.bin")), input)?;
                    std::fs::write(dir.join(format!("{stem}.{b}.out")), output)?;
                }
            }
            for ((_, a), (_, b)) in value_whole.iter().zip(&value_staged) {
                if a != b {
                    equal = false;
                }
            }
            if reference.is_none() {
                reference_terms =
                    Some(value_whole.iter().map(|(_, v)| terms_of(&v.expand())).sum());
                reference = Some(value_whole);
            }
            sizes = rep_sizes;
            staged.push(clocks);
        }
        let mut stage_json = serde_json::Map::new();
        for (index, name) in STAGES.iter().enumerate() {
            stage_json.insert(
                name.to_string(),
                json!({
                    "wall": median(staged.iter().map(|c| c.wall[index]).collect()),
                    "cpu": median(staged.iter().map(|c| c.cpu[index]).collect()),
                    "instructions": median(staged.iter().map(|c| c.instructions[index] as f64).collect()),
                }),
            );
        }
        // idenso alone: the edge-basis contraction (R3's) vs R2's routed contraction, paired
        let mut edge = serde_json::Value::Null;
        if edge_contract {
            let mut edge_runs = Vec::new();
            let mut edge_terms = 0usize;
            let (terms, _) = refpath::direct::colour_terms(&diagram, merge_colour)?;
            for _ in 0..reps {
                let wall = Instant::now();
                let cpu = cpu_seconds();
                let value = refpath::direct::contract_edges(
                    &diagram,
                    &terms,
                    &momentum,
                    &q_squared,
                    projection,
                    refpath::direct::DirectOptions {
                        dots: r3_dots,
                        closing: r3_closing,
                        frontier_order: r3_order,
                        routed: r3_routed,
                        compact: false,
                        prune_early: r3_prune,
                    },
                    None,
                )?;
                edge_runs.push((wall.elapsed().as_secs_f64(), cpu_seconds() - cpu));
                edge_terms = value.iter().map(|(_, v)| terms_of(&v.expand())).sum();
            }
            edge = json!({
                "wall": median(edge_runs.iter().map(|r| r.0).collect()),
                "cpu": median(edge_runs.iter().map(|r| r.1).collect()),
                "terms": edge_terms, "colour_terms": terms.len(), "dots": r3_dots, "routed": r3_routed, "merge_colour": merge_colour,
                "routed_terms": reference_terms.unwrap_or(0),
            });
        }
        // engine's colour stage (from_feynkit), per channel and merged
        let mut engine = serde_json::Map::new();
        if engine_colour {
            for (name, merge) in [("channels", false), ("merged", true)] {
                let mut walls = Vec::new();
                let mut values = 0;
                for _ in 0..reps {
                    let options = FromFeynkitOptions {
                        merge_colour: merge,
                        ..Default::default()
                    };
                    let t = Instant::now();
                    let result = FactorizedNumerator::from_feynkit(&diagram, &model, &options);
                    walls.push(t.elapsed().as_secs_f64());
                    values = result.map(|v| v.len()).unwrap_or(usize::MAX);
                }
                engine.insert(
                    name.into(),
                    json!({"wall": median(walls), "values": values}),
                );
            }
        }
        #[allow(unused_mut)]
        let mut phase_lines: Vec<String> = Vec::new();
        #[cfg(feature = "idenso-phase-timing")]
        if want_phases {
            let mut clocks = Clocks::default();
            let mut s = Sizes::default();
            contract_staged(
                &diagram,
                &branches,
                &momentum,
                &q_squared,
                metric_only,
                &mut clocks,
                counter.as_ref(),
                &mut s,
                None,
                Some(&mut phase_lines),
            )?;
        }
        let mut evaluation = serde_json::Value::Null;
        if let Some(template) = &template {
            let options = refpath::RunOptions {
                projection,
                extra_terms: 1,
                count_terms: false,
                contract_only: false,
            };
            let mut runs = Vec::new();
            let mut reference: Option<
                Vec<(
                    &'static str,
                    vakint::rustred_evaluation::p_gate::GateView,
                    Atom,
                )>,
            > = None;
            let mut direct_runs = Vec::new();
            let mut direct_gate = serde_json::Value::Null;
            let mut direct_terms = 0usize;
            let mut direct_counts = serde_json::Value::Null;
            let mut direct_colour = Vec::new();
            for rep in 0..reps {
                // paired, alternating: R2 then R3 on even repetitions, R3 then R2 on odd ones
                for turn in 0..2 {
                    let r2_turn = (turn == 0) == (rep % 2 == 0);
                    if r2_turn {
                        let engine = refpath::fresh_engine(template)?;
                        eprintln!("[idenso_prof] EVAL graph {stem} rep {rep}");
                        let wall = Instant::now();
                        let cpu = cpu_seconds();
                        let run = refpath::run_graph(&engine, &model, &dot, &options)?;
                        let (wall, cpu) = (wall.elapsed().as_secs_f64(), cpu_seconds() - cpu);
                        eprintln!("[idenso_prof] EVAL-END graph {stem} rep {rep}");
                        let keys: usize = run
                            .results
                            .iter()
                            .map(|r| refpath::gate_view(r).0.len())
                            .sum();
                        if reference.is_none() {
                            reference = Some(
                                run.results
                                    .iter()
                                    .map(|r| {
                                        (r.projection, refpath::gate_view(r).0, r.value.clone())
                                    })
                                    .collect(),
                            );
                        }
                        runs.push((run.seconds, wall, cpu, keys, run.contracted_terms));
                    } else if direct {
                        let engine = refpath::fresh_engine(template)?;
                        let settings = refpath::vakint_settings(diagram.loop_count() as i64 + 1);
                        let weight = feynkit_graph::expressions::evaluate_overall_factor(
                            diagram.overall_factor().as_view(),
                        ) * Atom::var(refpath::SYMBOLS.nf)
                            .pow(Atom::num(refpath::closed_quark_loops(&diagram)? as i64));
                        let mut seconds = refpath::direct::DirectSeconds::default();
                        let wall = Instant::now();
                        let cpu = cpu_seconds();
                        let colour_started = Instant::now();
                        let (terms, _) = refpath::direct::colour_terms(&diagram, merge_colour)?;
                        let colour_seconds = colour_started.elapsed().as_secs_f64();
                        let outcome = refpath::direct::lower_graph(
                            &engine,
                            &settings,
                            &diagram,
                            &terms,
                            &momentum,
                            &q_squared,
                            &weight,
                            projection,
                            true,
                            refpath::direct::DirectOptions {
                                dots: r3_dots,
                                closing: r3_closing,
                                frontier_order: r3_order,
                                routed: r3_routed,
                                compact: r3_compact,
                                prune_early: r3_prune,
                            },
                            &mut seconds,
                        )?;
                        direct_terms = terms.len();
                        direct_colour.push(colour_seconds);
                        let (wall, cpu) = (wall.elapsed().as_secs_f64(), cpu_seconds() - cpu);
                        let mut views = Vec::new();
                        let status = match &outcome {
                            refpath::direct::DirectOutcome::Zero(reason) => {
                                format!("zero: {reason}")
                            }
                            refpath::direct::DirectOutcome::Declined(reason) => {
                                format!("declined: {reason}")
                            }
                            refpath::direct::DirectOutcome::Lowered(projections) => {
                                for p in projections {
                                    let (view, _) =
                                        refpath::direct::gate_view(&engine, &settings, p)?;
                                    views.push((p.projection, view, p.value.clone()));
                                    let c = p.counts();
                                    direct_counts = json!({
                                        "numerator_terms": c.numerator_terms, "distinct_dots": c.distinct_dots,
                                        "canonical_terms": c.canonical_terms, "merged_terms": c.merged_terms,
                                        "distinct_monomials": c.distinct_monomials, "lowered": c.lowered,
                                        "folded": c.folded, "slot_products": c.slot_products, "pruned_early": c.pruned_early,
                                    });
                                }
                                "lowered".into()
                            }
                        };
                        direct_runs.push((seconds, wall, cpu, status, views));
                    }
                }
            }
            if direct {
                // gate R3 against R2's first repetition: views key by key, values exactly
                let reference = reference.clone().unwrap_or_default();
                let (_, _, _, status, views) = &direct_runs[0];
                let mut identical = 0;
                let mut differing = Vec::new();
                let mut values_equal = 0;
                for (name, view, value) in views {
                    match reference.iter().find(|(n, _, _)| n == name) {
                        Some((_, rview, rvalue)) => {
                            let difference =
                                vakint::rustred_evaluation::p_gate::compare_gate_views(rview, view);
                            if difference.is_empty() {
                                identical += 1;
                            } else {
                                differing.push(format!("{name}: {difference:?}"));
                            }
                            if value
                                .as_ref()
                                .is_some_and(|v| refpath::direct::values_equal(v, rvalue))
                            {
                                values_equal += 1;
                            }
                        }
                        None => differing.push(format!("{name}: not in R2")),
                    }
                }
                // and against lane refpath's STORED dump of the same graph and projection
                let mut stored = serde_json::Value::Null;
                if let Some(dir) = &ref_dir {
                    let mut stored_identical = 0;
                    let mut stored_missing = 0;
                    for (name, view, _) in views {
                        let file = dir.join("gate").join(format!("{stem}.{name}.gate"));
                        match std::fs::read_to_string(&file) {
                            Ok(text) => {
                                let (sview, _) =
                                    vakint::rustred_evaluation::p_gate::read_gate_view(&text)
                                        .map_err(|e| eyre::eyre!(e))?;
                                if vakint::rustred_evaluation::p_gate::compare_gate_views(
                                    &sview, view,
                                )
                                .is_empty()
                                {
                                    stored_identical += 1;
                                } else {
                                    differing.push(format!("{name}: differs from the stored dump"));
                                }
                            }
                            Err(_) => stored_missing += 1,
                        }
                    }
                    stored = json!({"identical": stored_identical, "missing": stored_missing});
                }
                let r2_nonzero = reference.iter().filter(|(_, v, _)| !v.is_empty()).count();
                let pass =
                    status == "lowered" && differing.is_empty() && identical == reference.len()
                        || (status.starts_with("zero") && r2_nonzero == 0);
                if !pass {
                    eprintln!(
                        "[idenso_prof] R3 GATE FAIL graph {stem}: {status} {differing:?} (R2 projections {})",
                        reference.len()
                    );
                }
                let d = |f: &dyn Fn(&refpath::direct::DirectSeconds) -> f64| {
                    median(direct_runs.iter().map(|r| f(&r.0)).collect())
                };
                direct_gate = json!({
                    "status": status, "pass": pass, "identical": identical, "values_equal": values_equal, "stored": stored,
                    "differing": differing,
                    "contract": d(&|s| s.contract), "matching": d(&|s| s.matching), "handover": d(&|s| s.handover),
                    "lower": d(&|s| s.lower), "reduce": d(&|s| s.reduce), "finalise": d(&|s| s.finalise),
                    "wall": median(direct_runs.iter().map(|r| r.1).collect()),
                    "cpu": median(direct_runs.iter().map(|r| r.2).collect()),
                    "colour": median(direct_colour.clone()), "terms": direct_terms, "merge_colour": merge_colour, "dots": r3_dots, "routed": r3_routed,
                    "compact": r3_compact, "prune_early": r3_prune, "routing": d(&|s| s.routing), "counts": direct_counts,
                    "routing_size": direct_runs.first().map(|r| r.0.routing_size),
                    "stored_routing_size": direct_runs.first().map(|r| r.0.stored_routing_size),
                    "trees": direct_runs.first().map(|r| r.0.trees),
                    "closing": r3_closing, "frontier_order": format!("{r3_order:?}"),
                });
            }
            let m = |f: &dyn Fn(&refpath::StageSeconds) -> f64| {
                median(runs.iter().map(|r| f(&r.0)).collect())
            };
            evaluation = json!({
                "parse": m(&|s| s.parse), "colour": m(&|s| s.colour), "contract": m(&|s| s.contract),
                "lower": m(&|s| s.lower), "reduce": m(&|s| s.reduce), "finalise": m(&|s| s.finalise),
                "ingress": m(&|s| s.ingress), "split": m(&|s| s.split), "lowering": m(&|s| s.lowering),
                "prepare": m(&|s| s.prepare), "library": m(&|s| s.library),
                "wall": median(runs.iter().map(|r| r.1).collect()),
                "cpu": median(runs.iter().map(|r| r.2).collect()),
                "keys_kept": runs.first().map(|r| r.3),
                "direct": direct_gate,
            });
        }
        let record = json!({
            "graph": graph,
            "stem": stem,
            "evaluation": evaluation,
            "branches": branches.len(),
            "branches_total": total,
            "parse": parse,
            "contract_whole_wall": median(whole.iter().map(|w| w.0).collect()),
            "contract_whole_cpu": median(whole.iter().map(|w| w.1).collect()),
            "stages": stage_json,
            "staged_equal_whole": equal,
            "cooked_bytes": sizes.cooked_bytes,
            "simplified_bytes": sizes.simplified_bytes,
            "scalar_terms": sizes.scalar_terms,
            "engine_colour": engine,
            "edge_contract": edge,
            "peak_rss": refpath::peak_rss_bytes(),
            "phases": phase_lines.iter().map(|l| serde_json::from_str::<serde_json::Value>(l).unwrap_or(json!(l))).collect::<Vec<_>>(),
            "instructions_available": counter.is_some(),
        });
        writeln!(out, "{record}")?;
        out.flush()?;
        if !equal {
            eprintln!("[idenso_prof] graph {graph}: staged replay differs from refpath::contract");
        }
    }
    Ok(())
}
