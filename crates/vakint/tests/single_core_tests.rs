//! `n_cores = 1` keeps every Vakint computation on the calling thread.
//!
//! The guarded test re-runs itself in a fresh process, where no thread pool
//! exists yet. Its test thread installs a seccomp filter that fails every
//! thread creation of that thread and of its descendants, then builds the
//! engine and evaluates Ordered RustRed reductions with each weight
//! arithmetic. Any spawned compute thread or thread pool fails the evaluation.
use symbolica::atom::Atom;
use vakint::{
    EvaluationMethod, EvaluationOrder, OrderedWeightArithmetic, RustRedEvaluationOptions,
    RustRedMasterMode, Vakint, VakintError, VakintSettings, vakint_parse,
};

/// Child-mode marker of the guarded test.
const CHILD: &str = "VAKINT_SINGLE_CORE_CHILD";
/// Printed by the child after the guarded evaluations completed.
const COMPLETED: &str = "single-core guarded evaluations completed";
/// Error of RustRed's polynomial modular arithmetic while it is a stub.
const POLYNOMIAL_STUB: &str = "modular_polynomial weight arithmetic is not implemented yet";

fn settings(n_cores: usize, arithmetic: OrderedWeightArithmetic) -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Retain,
            },
        )]),
        form_exe_path: "/single-core-tests-must-not-use-form".into(),
        use_dot_product_notation: true,
        n_cores,
        ordered_weight_arithmetic: arithmetic,
        ..Default::default()
    }
}

/// One Ordered batch per loop count through three loops, with numerators and
/// raised powers so that every request expands a nontrivial frontier.
fn lower_loop_input() -> Atom {
    vakint_parse!(
        "topo(I1L(muvsq,3))+dot(k(1),k(1))*topo(I2L(muvsq,2,1,1))\
         +topo(I3L(muvsq,2,1,1,1,1,1))+dot(k(1),k(2))*topo(I3L(muvsq,1,1,1,2,1,1))"
    )
    .unwrap()
}

fn four_loop_input() -> Atom {
    vakint_parse!("topo(I4L_H(muvsq,2,1,1,1,1,1,1,1,1))").unwrap()
}

/// Evaluate `input` with every arithmetic on `n_cores` and check that the
/// results agree with the exact one; the polynomial arithmetic may still be
/// RustRed's stub.
fn evaluate_with_every_arithmetic(vakint: &Vakint, input: &Atom, n_cores: usize) -> Atom {
    let exact = vakint
        .evaluate(
            &settings(n_cores, OrderedWeightArithmetic::exact()),
            input.as_view(),
        )
        .unwrap();
    assert!(!exact.is_zero());
    let sampled = vakint
        .evaluate(
            &settings(n_cores, OrderedWeightArithmetic::modular_sampled()),
            input.as_view(),
        )
        .unwrap();
    assert_eq!(sampled, exact, "modular_sampled on {n_cores} cores");
    match vakint.evaluate(
        &settings(n_cores, OrderedWeightArithmetic::modular_polynomial()),
        input.as_view(),
    ) {
        Ok(polynomial) => assert_eq!(polynomial, exact, "modular_polynomial on {n_cores} cores"),
        Err(error) => assert!(error.to_string().contains(POLYNOMIAL_STUB), "{error}"),
    }
    exact
}

/// Directory of a local external-data checkout with four-loop packs, if any.
fn external_data() -> Option<std::path::PathBuf> {
    std::env::var_os("VAKINT_TEST_EXTERNAL_DATA").map(Into::into)
}

fn load_four_loops(vakint: &mut Vakint, path: std::path::PathBuf, n_cores: usize) {
    vakint
        .load_official_external_data(
            4,
            &vakint::OfficialExternalDataOptions {
                local_artifacts_path: Some(path),
                load_prepared_caches: false,
                n_cores,
                ..Default::default()
            },
        )
        .unwrap();
}

#[cfg(all(
    target_os = "linux",
    any(target_arch = "x86_64", target_arch = "aarch64")
))]
mod thread_creation {
    /// `AUDIT_ARCH_*` of the compiled target (linux/audit.h).
    const AUDIT_ARCH: u32 = if cfg!(target_arch = "x86_64") {
        0xc000_003e
    } else {
        0xc000_00b7
    };
    const CLONE_THREAD: u32 = 0x0001_0000;

    fn statement(code: u32, k: u32) -> libc::sock_filter {
        libc::sock_filter {
            code: code as u16,
            jt: 0,
            jf: 0,
            k,
        }
    }

    fn jump(code: u32, k: u32, jt: u8, jf: u8) -> libc::sock_filter {
        libc::sock_filter {
            code: code as u16,
            jt,
            jf,
            k,
        }
    }

    /// Fail every thread creation of the calling thread and its descendants
    /// with `EPERM`; other system calls, including process creation, pass.
    /// `clone3` reports `ENOSYS` so that the C library falls back to `clone`,
    /// whose flags the filter can inspect.
    pub fn forbid() {
        use libc::{
            BPF_ABS, BPF_JEQ, BPF_JMP, BPF_JSET, BPF_K, BPF_LD, BPF_RET, BPF_W, ENOSYS, EPERM,
            SECCOMP_RET_ALLOW, SECCOMP_RET_ERRNO, SECCOMP_RET_KILL_PROCESS, SYS_clone, SYS_clone3,
        };
        // Offsets in `struct seccomp_data`: nr 0, arch 4, args[0] 16 (low
        // half on these little-endian targets).
        let mut filter = [
            statement(BPF_LD | BPF_W | BPF_ABS, 4),
            jump(BPF_JMP | BPF_JEQ | BPF_K, AUDIT_ARCH, 1, 0),
            statement(BPF_RET | BPF_K, SECCOMP_RET_KILL_PROCESS),
            statement(BPF_LD | BPF_W | BPF_ABS, 0),
            jump(BPF_JMP | BPF_JEQ | BPF_K, SYS_clone3 as u32, 0, 1),
            statement(BPF_RET | BPF_K, SECCOMP_RET_ERRNO | ENOSYS as u32),
            jump(BPF_JMP | BPF_JEQ | BPF_K, SYS_clone as u32, 1, 0),
            statement(BPF_RET | BPF_K, SECCOMP_RET_ALLOW),
            statement(BPF_LD | BPF_W | BPF_ABS, 16),
            jump(BPF_JMP | BPF_JSET | BPF_K, CLONE_THREAD, 0, 1),
            statement(BPF_RET | BPF_K, SECCOMP_RET_ERRNO | EPERM as u32),
            statement(BPF_RET | BPF_K, SECCOMP_RET_ALLOW),
        ];
        let program = libc::sock_fprog {
            len: filter.len() as u16,
            filter: filter.as_mut_ptr(),
        };
        // SAFETY: plain prctl calls; `program` and `filter` outlive them and
        // the kernel copies the filter.
        unsafe {
            assert_eq!(libc::prctl(libc::PR_SET_NO_NEW_PRIVS, 1, 0, 0, 0), 0);
            assert_eq!(
                libc::prctl(
                    libc::PR_SET_SECCOMP,
                    libc::SECCOMP_MODE_FILTER,
                    &program as *const libc::sock_fprog,
                ),
                0,
                "seccomp filter: {}",
                std::io::Error::last_os_error()
            );
        }
    }
}

#[cfg(all(
    target_os = "linux",
    any(target_arch = "x86_64", target_arch = "aarch64")
))]
#[test]
fn single_core_reductions_never_create_threads() {
    let name = "single_core_reductions_never_create_threads";
    if std::env::var_os(CHILD).is_none() {
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args(["--exact", name, "--nocapture", "--test-threads=1"])
            .env(CHILD, "1")
            .output()
            .unwrap();
        let stdout = String::from_utf8_lossy(&output.stdout);
        let stderr = String::from_utf8_lossy(&output.stderr);
        let completed = stdout.find(COMPLETED).filter(|_| output.status.success());
        let Some(start) = completed else {
            panic!(
                "guarded child failed ({}):\n{stdout}\n{stderr}",
                output.status
            );
        };
        println!("{}", stdout[start..].lines().next().unwrap_or_default());
        return;
    }
    // Symbolica checks a license key once per process, possibly verifying it
    // on a background thread: do it before thread creation is forbidden.
    // Without a key every check pings the license server on a new thread.
    if !symbolica::license::LicenseManager::is_licensed() {
        println!("{COMPLETED} (skipped: unlicensed Symbolica)");
        return;
    }
    let available = std::thread::available_parallelism().map_or(1, usize::from);
    thread_creation::forbid();
    assert!(
        std::thread::Builder::new().spawn(|| ()).is_err(),
        "the filter must forbid thread creation"
    );

    let mut vakint = Vakint::new().unwrap();
    let _ = evaluate_with_every_arithmetic(&vakint, &lower_loop_input(), 1);
    let mut evaluated_loops = "1-3";
    if let Some(path) = external_data() {
        load_four_loops(&mut vakint, path, 1);
        let _ = evaluate_with_every_arithmetic(&vakint, &four_loop_input(), 1);
        evaluated_loops = "1-4";
    }
    // Control: more cores reach RustRed's modular pool, which cannot start.
    if available > 1 {
        let error = vakint
            .evaluate(
                &settings(2, OrderedWeightArithmetic::modular_sampled()),
                lower_loop_input().as_view(),
            )
            .unwrap_err();
        assert!(error.to_string().contains("thread pool"), "{error}");
    }
    println!("{COMPLETED} ({evaluated_loops} loops)");
}

#[test]
fn results_do_not_depend_on_n_cores() {
    let available = std::thread::available_parallelism().map_or(1, usize::from);
    let mut vakint = Vakint::new().unwrap();
    let n_cores = vakint::default_n_cores().min(4).min(available);
    let input = lower_loop_input();
    assert_eq!(
        evaluate_with_every_arithmetic(&vakint, &input, 1),
        evaluate_with_every_arithmetic(&vakint, &input, n_cores)
    );
    if let Some(path) = external_data() {
        load_four_loops(&mut vakint, path, n_cores);
        let input = four_loop_input();
        assert_eq!(
            evaluate_with_every_arithmetic(&vakint, &input, 1),
            evaluate_with_every_arithmetic(&vakint, &input, n_cores)
        );
    }
}

#[test]
fn invalid_core_counts_and_arithmetics_are_rejected_before_any_work() {
    let vakint = Vakint::new().unwrap();
    let input = lower_loop_input();
    let error = vakint
        .evaluate(
            &settings(0, OrderedWeightArithmetic::exact()),
            input.as_view(),
        )
        .unwrap_err();
    assert!(
        matches!(error, VakintError::InvalidNCores { n_cores: 0, .. }),
        "{error}"
    );
    assert!(
        "modular_sampled:threads=0"
            .parse::<OrderedWeightArithmetic>()
            .is_err()
    );
    assert!("dfs".parse::<OrderedWeightArithmetic>().is_err());
}
