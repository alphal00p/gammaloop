//! Linux process observations shared by opt-in scalar benchmarks.
//!
//! CPU includes waited-for children (FORM); RSS is the parent snapshot only.

use std::time::Instant;
use symbolica::atom::Atom;
use vakint::{Vakint, VakintSettings};

pub(super) struct Observation {
    pub(super) value: Atom,
    wall_ns: u128,
    cpu_seconds: Option<f64>,
    rss_before_kib: Option<u64>,
    rss_after_kib: Option<u64>,
}

impl Observation {
    pub(super) fn measure(
        vakint: &Vakint,
        settings: &VakintSettings,
        scalar: &Atom,
        ticks: Option<f64>,
    ) -> Self {
        let rss_before_kib = resident_kib();
        let cpu_before = cpu_ticks();
        let started = Instant::now();
        let value = vakint
            .evaluate_integral(settings, scalar.as_view())
            .unwrap();
        let wall_ns = started.elapsed().as_nanos();
        let cpu_after = cpu_ticks();
        let rss_after_kib = resident_kib();
        let cpu_seconds =
            cpu_before
                .zip(cpu_after)
                .zip(ticks)
                .and_then(|((before, after), ticks)| {
                    after.checked_sub(before).map(|delta| delta as f64 / ticks)
                });
        Self {
            value,
            wall_ns,
            cpu_seconds,
            rss_before_kib,
            rss_after_kib,
        }
    }

    pub(super) fn report(&self, case: &str, phase: &str, lane: &str, repeat: usize) {
        println!(
            "PUBLIC_SCALAR_TIMING\t{case}\t{phase}\t{lane}\t{repeat}\t{}\t{:?}\t{:?}\t{:?}",
            self.wall_ns, self.cpu_seconds, self.rss_before_kib, self.rss_after_kib
        );
    }
}

pub fn cpu_ticks() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/stat").ok()?;
    let fields = status
        .rsplit_once(')')?
        .1
        .split_whitespace()
        .collect::<Vec<_>>();
    // Fields 14--17: user, system, waited-for child user and child system.
    [11, 12, 13, 14].into_iter().try_fold(0u64, |sum, field| {
        sum.checked_add(fields.get(field)?.parse::<u64>().ok()?)
    })
}

pub fn resident_kib() -> Option<u64> {
    std::fs::read_to_string("/proc/self/status")
        .ok()?
        .lines()
        .find_map(|line| line.strip_prefix("VmRSS:"))?
        .split_whitespace()
        .next()?
        .parse()
        .ok()
}

pub fn clock_ticks_per_second() -> Option<f64> {
    let output = std::process::Command::new("getconf")
        .arg("CLK_TCK")
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    String::from_utf8(output.stdout)
        .ok()?
        .trim()
        .parse::<f64>()
        .ok()
        .filter(|ticks| ticks.is_finite() && *ticks > 0.0)
}
