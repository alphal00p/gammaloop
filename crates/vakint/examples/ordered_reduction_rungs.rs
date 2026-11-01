//! Reduce native massless-quark rung fixtures with RustRed's ordered scheduler.
//! Build with `--release --features feynkit-ingress,experimental-rustred`;
//! run one weight arithmetic per fresh process for performance measurements, or
//! `--arithmetic both` for exact exact/modular parity. Reduction math and the
//! complete retained public call have separate timing scopes.

use ::vakint;
#[path = "../tests/fixtures/qcd_ladder/mod.rs"]
#[allow(dead_code)]
mod fixture;

use std::{io::Read, path::PathBuf, process::Command};
use vakint::{RustRedMonitorDisplay, RustRedMonitorOptions};

fn main() -> fixture::Result<()> {
    let mut arguments = std::env::args().skip(1);
    let mut work = None;
    let mut loop_count = 1_usize;
    let mut projection = fixture::Projection::ZeroMomentumMetric;
    let mut arithmetic = None;
    let mut artifacts = None;
    let mut prepare_only = false;
    let mut prepared_input = None;
    let mut numerical = true;
    let mut tui = false;
    let mut no_monitor = false;
    let mut snapshot = None;
    let mut build_profile = None;
    let mut compare_two_reports = None;
    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "--work-dir" => {
                work = Some(PathBuf::from(
                    arguments.next().ok_or("--work-dir needs a path")?,
                ))
            }
            "--loops" => {
                loop_count = arguments
                    .next()
                    .ok_or("--loops needs a positive integer")?
                    .parse()?
            }
            "--projection" => {
                projection = match arguments.next().as_deref() {
                    Some("zero-momentum") => fixture::Projection::ZeroMomentumMetric,
                    Some("quadratic") => fixture::Projection::QuadraticTransverse,
                    _ => return Err("--projection must be zero-momentum or quadratic".into()),
                }
            }
            "--arithmetic" => {
                arithmetic = Some(
                    arguments
                        .next()
                        .ok_or("--arithmetic needs exact, modular[:key=value,...], or both")?,
                )
            }
            "--artifacts" => {
                artifacts = Some(PathBuf::from(
                    arguments.next().ok_or("--artifacts needs a path")?,
                ))
            }
            "--prepare-only" => prepare_only = true,
            "--build-profile" => {
                let value = arguments
                    .next()
                    .ok_or("--build-profile needs release, ci-optim, or dev")?;
                if !["release", "ci-optim", "dev"].contains(&value.as_str()) {
                    return Err("unknown build profile declaration".into());
                }
                build_profile = Some(value);
            }
            "--compare-two-reports" => {
                compare_two_reports = Some((
                    PathBuf::from(
                        arguments
                            .next()
                            .ok_or("--compare-two-reports needs a reference report directory")?,
                    ),
                    PathBuf::from(
                        arguments
                            .next()
                            .ok_or("--compare-two-reports needs an ordered report directory")?,
                    ),
                ))
            }
            "--prepared-input" => {
                prepared_input = Some(PathBuf::from(
                    arguments.next().ok_or("--prepared-input needs a path")?,
                ))
            }
            "--no-numerical-reference" => numerical = false,
            "--monitor" => tui = true,
            "--no-monitor" => no_monitor = true,
            "--progress-snapshot" => {
                snapshot = Some(PathBuf::from(
                    arguments.next().ok_or("--progress-snapshot needs a path")?,
                ))
            }
            "--help" => {
                println!(
                    "ordered_reduction_rungs --work-dir PATH --build-profile release [--loops N] [--projection zero-momentum|quadratic] [--arithmetic exact|modular_sampled[:key=value,...]|modular_polynomial[:key=value,...]|both] [--artifacts PATH] [--prepared-input PATH] [--prepare-only] [--no-numerical-reference] [--no-monitor | --monitor] [--progress-snapshot PATH]\nFor a functional proof between saved reports: --work-dir OUT --compare-two-reports REFERENCE_DIRECTORY ORDERED_DIRECTORY, where the reference is another Ordered report or a preserved historical depth-first report; --no-numerical-reference explicitly compares retained expressions only without changing the saved reports.\nDefault: one loop, zero momentum; the arithmetic of VAKINT_ORDERED_WEIGHT_ARITHMETIC (Vakint's default when unset); tracked numerical reference through four loops. `both` runs exact weights and the modular arithmetic of VAKINT_ORDERED_WEIGHT_ARITHMETIC (modular_polynomial when it is unset or exact) in one process and checks that they retain the same expression. No prepared or previous-run cache is restored unless requested: Ordered step caches follow VAKINT_ORDERED_CACHE, VAKINT_ORDERED_CACHE_DIR and VAKINT_ORDERED_CACHE_SHIPPED (shipped caches are ignored unless it is on; a run that replays one records cache_policy.prepared_cache_loaded). Saved inputs under work-dir/input-LN-projection are authenticated and reused; --prepared-input selects a shared or legacy bundle. Reduction math is captured with the experimental-rustred feature; the complete retained public call is timed separately. Both timings exclude the subsequent numerical reference. Use --no-monitor for qualified performance comparisons: it disables both dashboard and snapshot callbacks and conflicts with --monitor and --progress-snapshot. Use a release build without profile overrides.\nEnvironment: VAKINT_N_CORES (worker cores; default min(10, available cores)) and VAKINT_ORDERED_WEIGHT_ARITHMETIC (exact | modular_sampled[:key=value,...] | modular_polynomial[:key=value,...]; default modular_polynomial; modular workers follow the cores unless threads is set); both are recorded under report.json execution and never change the exact results."
                );
                return Ok(());
            }
            _ => return Err(format!("unknown argument {argument}").into()),
        }
    }
    if no_monitor && (tui || snapshot.is_some()) {
        return Err("--no-monitor conflicts with --monitor and --progress-snapshot".into());
    }
    if loop_count == 0 {
        return Err("--loops must be positive".into());
    }
    let arms = match arithmetic.as_deref() {
        None => vec![(
            "ordered",
            fixture::comparison::arithmetic_from_environment()?,
        )],
        Some("both") => vec![
            ("exact", vakint::OrderedWeightArithmetic::exact()),
            (
                "modular",
                fixture::comparison::modular_arithmetic_from_environment()?,
            ),
        ],
        Some(text) => vec![("ordered", text.parse::<vakint::OrderedWeightArithmetic>()?)],
    };
    let work = work.ok_or("--work-dir is required; keep development data under workspace")?;
    std::fs::create_dir_all(&work)?;
    if let Some((reference, ordered)) = compare_two_reports {
        fixture::comparison::compare_two_reports(
            &reference,
            &ordered,
            &work.join("functional-pair.json"),
            numerical,
        )?;
        println!(
            "Authenticated functional two-report comparison passed; no performance qualification."
        );
        return Ok(());
    }
    if !prepare_only && build_profile.is_none() {
        return Err(
            "declare the actual Cargo profile with --build-profile release, ci-optim, or dev"
                .into(),
        );
    }
    let projection_name = match projection {
        fixture::Projection::ZeroMomentumMetric => "zero-momentum",
        fixture::Projection::QuadraticTransverse => "quadratic",
    };
    let input_directory = work.join(format!("input-L{loop_count}-{projection_name}"));
    let prepared = match prepared_input {
        Some(directory) => {
            fixture::LadderFixture::load_prepared_at_loop_count(&directory, projection, loop_count)?
        }
        None => fixture::comparison::prepared(&input_directory, projection, loop_count)?,
    };
    if prepare_only {
        println!(
            "Prepared {loop_count}L {projection_name}: {} contributions, {} stored scalar terms",
            prepared.contributions, prepared.scalar_terms
        );
        return Ok(());
    }
    let mut run_directories = Vec::new();
    // Development binaries contain gigabytes of debug information. Their
    // supervising launch receipt authenticates them; avoid repeating that
    // metadata work with unoptimized hashing. Release hashes stream once.
    let executable_sha256 = if cfg!(debug_assertions) {
        None
    } else {
        let mut file = std::fs::File::open(std::env::current_exe()?)?;
        let mut digest = sha2::Sha256::new();
        let mut buffer = [0_u8; 65_536];
        loop {
            let bytes = file.read(&mut buffer)?;
            if bytes == 0 {
                break;
            }
            digest.update(&buffer[..bytes]);
        }
        Some(format!("{:x}", digest.finalize()))
    };
    for (name, arithmetic) in &arms {
        let directory = work.join(format!("L{loop_count}-{projection_name}-{name}"));
        if directory.join("report.json").exists() {
            return Err(format!(
                "preserve existing {}; choose a new work directory",
                directory.display()
            )
            .into());
        }
        let options = fixture::comparison::RunOptions {
            arithmetic: *arithmetic,
            artifacts: artifacts.clone(),
            numerical_reference: numerical && loop_count <= 4,
            monitor: RustRedMonitorOptions {
                display: if !no_monitor && tui {
                    RustRedMonitorDisplay::Tui
                } else {
                    RustRedMonitorDisplay::Off
                },
                progress_snapshot: snapshot.as_ref().filter(|_| !no_monitor).map(|path| {
                    if arms.len() > 1 {
                        path.with_extension(format!("{name}.json"))
                    } else {
                        path.clone()
                    }
                }),
                ..Default::default()
            },
        };
        let mut run = fixture::comparison::run(&prepared, &directory, &options)?;
        let git = |args: &[&str]| {
            Command::new("git")
                .args(args)
                .output()
                .ok()
                .filter(|out| out.status.success())
                .map(|out| String::from_utf8_lossy(&out.stdout).trim().to_owned())
        };
        run.report["producer"] = serde_json::json!({"commit":git(&["rev-parse","HEAD"]),"debug_assertions":cfg!(debug_assertions),
            "build_profile":build_profile,
            "cargo_lock_sha256":format!("{:x}",sha2::Sha256::digest(std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"),"/../../Cargo.lock"))?)),
            "executable_sha256":executable_sha256,
            "executable_identity_scope":"release: streamed once; DEV: use the supervising original-process launch receipt",
            "native_runtime_arities":option_env!("RUSTRED_RUNTIME_ARITIES"),
            "allowed_cpus":std::fs::read_to_string("/proc/self/status").ok().and_then(|status| status.lines().find_map(|line| line.strip_prefix("Cpus_allowed_list:").map(|value|value.trim().to_owned()))),
            "profile_contract":"caller must build --release with no profile overrides; debug_assertions is not compiler-profile introspection",
            "comparison_process":if arms.len() > 1 {"two fresh Vakint instances in one process; use separate invocations for fresh-process timing"} else {"one engine in this process"}});
        std::fs::write(
            directory.join("report.json"),
            serde_json::to_vec_pretty(&run.report)?,
        )?;
        let execution = &run.report["execution"];
        println!(
            "{name}: {} weights on {} cores",
            execution["ordered_weight_arithmetic"]
                .as_str()
                .unwrap_or("unrecorded"),
            execution["n_cores"]
        );
        let public_seconds = run.report["retained_reduction_seconds"]
            .as_f64()
            .unwrap_or_default();
        if let Some(math_seconds) = run.report["reduction_math_seconds"].as_f64() {
            println!(
                "{name}: {math_seconds:.6}s reduction math; {public_seconds:.6}s complete retained public call"
            );
        } else {
            println!(
                "{name}: reduction math unavailable; {public_seconds:.6}s complete retained public call"
            );
        }
        run_directories.push(directory);
    }
    if let [exact, modular] = run_directories.as_slice() {
        fixture::comparison::compare_two_reports(
            exact,
            modular,
            &work.join(format!("L{loop_count}-{projection_name}-parity.json")),
            false,
        )?;
        println!("Exact and modular weights retained the same expression.");
    }
    Ok(())
}

use sha2::Digest;
