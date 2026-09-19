//! Linux process observations shared by opt-in scalar benchmarks.
//!
//! CPU includes waited-for children (FORM); RSS is the parent snapshot only.

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
