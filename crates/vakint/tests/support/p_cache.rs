//! The fixed workload of the light p cache shipped in vakint-external-data
//! (M2-V6; M3-V10) for the p library program of scope `1..=scope`
//! ([`p_cache_workload`]), in its fixed order: the 45 one- to three-loop cells
//! of the benchmark grid (`tests/fixtures/benchmark_p_integrals.json`, table
//! order), the two- and three-loop physics cells (`benchmark_p_physics.json`),
//! the Forcer batteries (`forcer/p{1,2,3}l_battery.json`, 87 cases, including
//! the inputs carrying `Master(no)` and `Master(t105)`), the anchors
//! (`forcer/p_anchors.json`, 11) and, for a scope reaching four loops, the
//! four-loop anchors (`forcer/p4l_anchors.json`, 17: the four-loop masters and
//! three tops). The four-loop benchmark cells stay out: recording them takes
//! hours and gigabytes (M3-T19, plan M3 D-M3-8), far above the light budget.
//!
//! Include it with `#[path = "support/p_cache.rs"] mod p_cache;` (the example
//! `p_cache_maintenance` uses `#[path = "../tests/support/p_cache.rs"]`).

#![allow(dead_code)]

use std::path::PathBuf;

use sha2::{Digest, Sha256};
use symbolica::atom::{Atom, AtomCore};
use vakint::p_integral::forcer_fixtures::{BenchmarkInputs, ForcerFixture};
use vakint::vakint_parse;

/// One input of the workload.
pub struct PCacheInput {
    /// `<fixture stem>/<case>`.
    pub name: String,
    /// The fixture file, relative to `crates/vakint`.
    pub fixture: &'static str,
    pub case: String,
    pub loop_count: usize,
    pub input: Atom,
}

impl PCacheInput {
    /// The workload row of `ordered-cache.json` (the caller adds `inserts`
    /// and `seconds`).
    pub fn row(&self) -> serde_json::Value {
        serde_json::json!({
            "name": self.name,
            "class": "training",
            "fixture": self.fixture,
            "case": self.case,
            "loop_count": self.loop_count,
            "input_sha256": format!(
                "{:x}",
                Sha256::digest(self.input.to_canonical_string().as_bytes())
            ),
        })
    }
}

/// Number of inputs of [`p_cache_workload`] for a scope of three loops.
pub const P_CACHE_WORKLOAD_INPUTS: usize = 45 + 2 + 87 + 11;
/// The four-loop anchors a scope of four loops adds.
pub const P_CACHE_FOUR_LOOP_ANCHORS: usize = 17;

fn read(fixture: &str) -> Result<String, String> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(fixture);
    std::fs::read_to_string(&path).map_err(|error| format!("{}: {error}", path.display()))
}

fn input(fixture: &'static str, case: &str, loop_count: usize, input: Atom) -> PCacheInput {
    let stem = fixture.rsplit('/').next().unwrap_or(fixture);
    PCacheInput {
        name: format!("{}/{case}", stem.trim_end_matches(".json")),
        fixture,
        case: case.to_owned(),
        loop_count,
        input,
    }
}

/// The workload of the program of scope `1..=scope`, three or four loops
/// (module documentation).
pub fn p_cache_workload(scope: usize) -> Result<Vec<PCacheInput>, String> {
    let expected = match scope {
        3 => P_CACHE_WORKLOAD_INPUTS,
        4 => P_CACHE_WORKLOAD_INPUTS + P_CACHE_FOUR_LOOP_ANCHORS,
        _ => return Err(format!("no p cache workload for the scope 1-{scope}")),
    };
    let mut inputs = Vec::new();
    const GRID: &str = "tests/fixtures/benchmark_p_integrals.json";
    let grid = BenchmarkInputs::parse(&read(GRID)?).map_err(|error| error.to_string())?;
    for cell in grid.cells.iter().filter(|cell| cell.loops <= 3) {
        let atom = vakint_parse!(&cell.input).map_err(|error| format!("{}: {error}", cell.name))?;
        inputs.push(input(GRID, &cell.name, cell.loops, atom));
    }
    const PHYSICS: &str = "tests/fixtures/benchmark_p_physics.json";
    let physics: serde_json::Value =
        serde_json::from_str(&read(PHYSICS)?).map_err(|error| error.to_string())?;
    for cell in physics["cells"].as_array().ok_or("physics cells")? {
        let loops = cell["loops"].as_u64().ok_or("physics loops")? as usize;
        if !(2..=3).contains(&loops) {
            continue;
        }
        let name = cell["name"].as_str().ok_or("physics name")?;
        let text = cell["input"].as_str().ok_or("physics input")?;
        let atom = vakint_parse!(text).map_err(|error| format!("{name}: {error}"))?;
        inputs.push(input(PHYSICS, name, loops, atom));
    }
    let mut fixtures = vec![
        ("tests/fixtures/forcer/p1l_battery.json", Some(1)),
        ("tests/fixtures/forcer/p2l_battery.json", Some(2)),
        ("tests/fixtures/forcer/p3l_battery.json", Some(3)),
        ("tests/fixtures/forcer/p_anchors.json", None),
    ];
    if scope >= 4 {
        fixtures.push(("tests/fixtures/forcer/p4l_anchors.json", Some(4)));
    }
    for (fixture, loops) in fixtures {
        let parsed = ForcerFixture::parse(&read(fixture)?).map_err(|error| error.to_string())?;
        for case in &parsed.cases {
            let count = match loops {
                Some(loops) => loops,
                None => case.loops().map_err(|error| error.to_string())?,
            };
            let atom = case.expression().map_err(|error| error.to_string())?;
            inputs.push(input(fixture, &case.name, count, atom));
        }
    }
    if inputs.len() != expected {
        return Err(format!(
            "the p cache workload of the scope 1-{scope} has {} inputs, expected {expected}",
            inputs.len()
        ));
    }
    Ok(inputs)
}
