//! Opt-in reproducer for the R* hard-Taylor simplification
//! (`gammalooprs::uv::rstar::RStarIntegrand::simplify_view`) on a dumped
//! input: `IDENSO_SIMPLIFY_RECORD=<file with one canonical atom string>`
//! (written by GammaLoop's `GAMMALOOP_ALGEBRA_DUMP`). With
//! `IDENSO_SIMPLIFY_OUT=<file>` the canonical result is written for a byte
//! comparison of two builds. Prints the time and the process high-water mark.

use idenso::dirac::GammaSimplifySettings;
use idenso::tensor::{AlgebraSettings, SymbolicTensor};
use idenso::{CookMode, CookSettings};
use symbolica::atom::{Atom, AtomCore};

fn high_water_mark() -> String {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|status| {
            status
                .lines()
                .find_map(|line| Some(line.strip_prefix("VmHWM:")?.trim().to_string()))
        })
        .unwrap_or_default()
}

#[test]
fn simplify_dumped_hard_taylor_input() {
    let Some(path) = std::env::var_os("IDENSO_SIMPLIFY_RECORD") else {
        return;
    };
    idenso::representations::initialize();
    let _ = *idenso::epsilon::EPSILON_SYMBOL;
    let text = std::fs::read_to_string(&path).unwrap();
    let input = Atom::parse(
        text.trim(),
        "gammalooprs",
        symbolica::parser::ParseSettings::symbolica(),
    )
    .unwrap();
    eprintln!(
        "input {} bytes, high-water mark {}",
        input.as_view().get_byte_size(),
        high_water_mark()
    );
    let started = std::time::Instant::now();
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let simplified = SymbolicTensor::infer_indexed(cooking.try_cook(input.as_view()).unwrap())
        .unwrap()
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::canonical()),
            epsilon: true,
            ..Default::default()
        })
        .unwrap()
        .to_dots()
        .unwrap();
    let output = cooking.uncook(simplified.expression().as_view());
    eprintln!(
        "simplified to {} bytes in {:.2} s, high-water mark {}",
        output.as_view().get_byte_size(),
        started.elapsed().as_secs_f64(),
        high_water_mark()
    );
    if let Some(out) = std::env::var_os("IDENSO_SIMPLIFY_OUT") {
        std::fs::write(out, output.to_canonical_string()).unwrap();
    }
}
