//! Frozen Forcer outputs: the fixture files `tests/fixtures/forcer/*.json`.
//!
//! A fixture freezes the exact Forcer output of canonical p-integral cases
//! (family, line powers, numerator in the canonical `k(j)` and the user head,
//! external momentum `σ·p(i)`), one raw `%E` text per Forcer local, together
//! with the Forcer commit and the FORM version that produced it. Tests compare
//! other backends with these values without running FORM; with FORM 5 the
//! writer regenerates a file (`VAKINT_FORCER_FIXTURES=write`).
//!
//! Files:
//! - `p{1,2,3,4}l_battery.json`: the M1 batteries (hand-picked cases);
//! - `p{1,2,3}l_grid.json`: [`grid_cases`], the benchmark cells
//!   (`tests/fixtures/benchmark_p_integrals.json`), every non-scaleless
//!   labelled sector of the registry families at rank 0 and 2, and dotted
//!   sectors;
//! - `p_anchors.json`: [`anchor_cases`], the normalization anchors of the M2
//!   plan (section 3.3), including the parents of the eps-line kite `t1star05`.
//!
//! This module is public for the integration tests only; it is not part of
//! vakint's API.

use std::sync::atomic::{AtomicUsize, Ordering};

use serde::{Deserialize, Serialize};
use symbolica::atom::{Atom, Symbol};

use super::{ForcerInput, combine_parts, forcer_input, run_forcer};
use crate::p_integral::{FORCER_UPSTREAM_COMMIT, P_FAMILIES, PFamilyData, ReducedPropagatorValue};
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{Vakint, VakintError, VakintSettings};

fn malformed(detail: impl std::fmt::Display) -> VakintError {
    VakintError::InvalidGenericExpression(format!("Forcer fixture: {detail}"))
}

/// One fixture file.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct ForcerFixture {
    pub forcer_commit: String,
    pub form_version: String,
    /// What the file holds and how it was made (newer files only).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub provenance: Option<String>,
    pub cases: Vec<ForcerCase>,
}

/// One canonical case and its frozen Forcer output.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct ForcerCase {
    pub name: String,
    /// A registry family name (`P1L`, `P2L`, `P3L_LA`, ...).
    pub family: String,
    /// The power of every canonical line (0 = pinched, negative = numerator).
    pub powers: Vec<i64>,
    /// A scalar polynomial in the canonical `k(j)` and the user head.
    pub numerator: String,
    /// The canonical `p(1)` written in the user head, `σ·p(i)`.
    pub external: String,
    /// The raw `%E` output of each Forcer local, in the order of the parts of
    /// the program; empty for a scaleless case.
    pub outputs: Vec<String>,
}

impl ForcerFixture {
    /// Parse a fixture file; its Forcer commit must be the vendored one.
    pub fn parse(json: &str) -> Result<Self, VakintError> {
        let fixture: ForcerFixture = serde_json::from_str(json).map_err(malformed)?;
        if fixture.forcer_commit != FORCER_UPSTREAM_COMMIT {
            return Err(malformed(format!(
                "Forcer commit {} instead of the vendored {FORCER_UPSTREAM_COMMIT}",
                fixture.forcer_commit
            )));
        }
        Ok(fixture)
    }

    /// The file text: pretty JSON and a final newline.
    pub fn to_json(&self) -> String {
        serde_json::to_string_pretty(self).expect("fixtures serialize") + "\n"
    }

    /// A fixture of `cases` without outputs, ready for [`Self::write_outputs`].
    pub fn unwritten(provenance: &str, cases: Vec<ForcerCase>) -> Self {
        ForcerFixture {
            forcer_commit: FORCER_UPSTREAM_COMMIT.into(),
            form_version: String::new(),
            provenance: Some(provenance.into()),
            cases,
        }
    }

    /// The cases with their outputs removed, to compare with a definition.
    pub fn definitions(&self) -> Vec<ForcerCase> {
        self.cases
            .iter()
            .map(|case| ForcerCase {
                outputs: Vec::new(),
                ..case.clone()
            })
            .collect()
    }

    /// Run Forcer (FORM 5 at `form`, one worker per FORM process) on every
    /// case with `jobs` FORM processes at a time and store the outputs and the
    /// FORM version.
    pub fn write_outputs(&mut self, form: &str, jobs: usize) -> Result<(), VakintError> {
        self.form_version =
            super::forcer_form_version(form).map_err(|reason| VakintError::ForcerUnavailable {
                form_exe_path: form.into(),
                reason,
            })?;
        let programs = self
            .cases
            .iter()
            .map(ForcerCase::input)
            .collect::<Result<Vec<_>, _>>()?;
        let settings = VakintSettings {
            form_exe_path: form.into(),
            n_cores: 1,
            ..VakintSettings::default()
        };
        let next = AtomicUsize::new(0);
        let mut outputs = (0..programs.len()).map(|_| None).collect::<Vec<_>>();
        std::thread::scope(|scope| {
            let workers = (0..jobs.max(1))
                .map(|_| {
                    scope.spawn(|| {
                        let mut done = Vec::new();
                        loop {
                            let index = next.fetch_add(1, Ordering::Relaxed);
                            let Some(program) = programs.get(index) else {
                                return done;
                            };
                            let result = match program {
                                None => Ok(Vec::new()),
                                Some(input) => {
                                    let names = input
                                        .parts
                                        .iter()
                                        .map(|part| part.name.as_str())
                                        .collect::<Vec<_>>();
                                    run_forcer(&settings, &input.program, &names)
                                }
                            };
                            done.push((index, result));
                        }
                    })
                })
                .collect::<Vec<_>>();
            for worker in workers {
                for (index, result) in worker.join().expect("Forcer fixture worker") {
                    outputs[index] = Some(result);
                }
            }
        });
        for (case, output) in self.cases.iter_mut().zip(outputs) {
            case.outputs = output
                .expect("every case ran")
                .map_err(|error| malformed(format!("{}: {error}", case.name)))?;
        }
        Ok(())
    }

    /// `blake3` over the name and the Forcer program of every case, in
    /// order (a scaleless case contributes `scaleless`): the guard that keeps
    /// the generated programs byte-identical.
    pub fn program_digest(&self) -> Result<String, VakintError> {
        let mut hasher = blake3::Hasher::new();
        for case in &self.cases {
            hasher.update(case.name.as_bytes());
            hasher.update(b"\0");
            match case.input()? {
                Some(input) => hasher.update(input.program.as_bytes()),
                None => hasher.update(b"scaleless"),
            };
            hasher.update(b"\n");
        }
        Ok(hasher.finalize().to_hex().to_string())
    }
}

fn parse(text: &str) -> Result<Atom, VakintError> {
    Vakint::initialize_vakint_symbols();
    vk_parse!(text).map_err(|error| malformed(format!("{text}: {error}")))
}

pub(crate) fn registry_family(name: &str) -> Result<&'static PFamilyData, VakintError> {
    P_FAMILIES
        .iter()
        .find(|family| family.name == name)
        .ok_or_else(|| malformed(format!("unknown family {name}")))
}

impl ForcerCase {
    /// A case without output.
    pub fn new(
        name: String,
        family: &str,
        powers: Vec<i64>,
        numerator: &str,
        external: &str,
    ) -> Self {
        ForcerCase {
            name,
            family: family.into(),
            powers,
            numerator: numerator.into(),
            external: external.into(),
            outputs: Vec::new(),
        }
    }

    /// The loop count of the family.
    pub fn loops(&self) -> Result<usize, VakintError> {
        Ok(registry_family(&self.family)?.loops)
    }

    /// The Forcer program and its parts, `None` for a scaleless case.
    pub(crate) fn input(&self) -> Result<Option<ForcerInput>, VakintError> {
        forcer_input(
            registry_family(&self.family)?,
            &self.powers,
            parse(&self.numerator)?.as_view(),
            &parse(&self.external)?,
        )
        .map_err(|error| malformed(format!("{}: {error}", self.name)))
    }

    /// The Forcer program of this case, `None` when it is scaleless.
    pub fn program(&self) -> Result<Option<String>, VakintError> {
        Ok(self.input()?.map(|input| input.program))
    }

    /// The frozen reduced values in `epsilon`: the backend contract of
    /// `docs/p_integrals.md` (section 5) at the frozen output.
    pub(crate) fn frozen(
        &self,
        epsilon: Symbol,
    ) -> Result<Vec<ReducedPropagatorValue>, VakintError> {
        match self.input()? {
            None if self.outputs.is_empty() => Ok(Vec::new()),
            None => Err(malformed(format!(
                "{}: outputs for a scaleless case",
                self.name
            ))),
            Some(input) if input.parts.len() != self.outputs.len() => Err(malformed(format!(
                "{}: {} outputs for {} parts",
                self.name,
                self.outputs.len(),
                input.parts.len()
            ))),
            Some(input) => combine_parts(&input.parts, &self.outputs, epsilon),
        }
    }

    /// The frozen `(scale_power, reduced)` pairs in the symbol `epsilon`,
    /// `reduced = (-1)^(Σa+r)·C_G` in the common p-master basis.
    pub fn frozen_values(&self, epsilon: &str) -> Result<Vec<(i64, Atom)>, VakintError> {
        Ok(self
            .frozen(vk_symbol!(epsilon))?
            .into_iter()
            .map(|value| (value.scale_power, value.reduced))
            .collect())
    }

    /// The frozen value finalized by the shared p-finalizer in `mode` with
    /// `settings` (normalization, Laurent depth, symbols), exactly as
    /// `Vakint::evaluate` with Forcer would return it.
    pub fn finalized(
        &self,
        settings: &VakintSettings,
        mode: crate::rustred_evaluation::RustRedMasterMode,
    ) -> Result<Atom, VakintError> {
        let values = self.frozen(vk_symbol!(&settings.epsilon_symbol))?;
        crate::p_integral::finalize::finalize(
            &values,
            self.loops()?,
            &parse(&self.external)?,
            settings,
            mode,
            "Forcer",
        )
    }

    /// The case as a vakint input, `(numerator)*topo(NAME(external, powers...))`.
    pub fn expression(&self) -> Result<Atom, VakintError> {
        let powers = self
            .powers
            .iter()
            .map(i64::to_string)
            .collect::<Vec<_>>()
            .join(",");
        parse(&format!(
            "({})*topo({}({},{powers}))",
            self.numerator, self.family, self.external
        ))
    }
}

/// One cell of `tests/fixtures/benchmark_p_integrals.json`.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct BenchmarkCell {
    pub name: String,
    pub loops: usize,
    pub dots: usize,
    pub rank: usize,
    pub numerator: String,
    pub tops: Vec<BenchmarkTop>,
    /// The frozen vakint input: `Σ_tops (numerator)*top`, long forms.
    pub input: String,
}

/// One top family of a benchmark cell and its line powers.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct BenchmarkTop {
    pub family: String,
    pub powers: Vec<i64>,
}

/// `tests/fixtures/benchmark_p_integrals.json`, written by
/// `tools/p_benchmark/generate_inputs.py`.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct BenchmarkInputs {
    pub schema: String,
    pub generator: String,
    pub grid: serde_json::Value,
    pub cells: Vec<BenchmarkCell>,
}

impl BenchmarkInputs {
    pub const SCHEMA: &'static str = "vakint.p-benchmark-inputs.v1";

    pub fn parse(json: &str) -> Result<Self, VakintError> {
        let inputs: BenchmarkInputs = serde_json::from_str(json).map_err(malformed)?;
        if inputs.schema != Self::SCHEMA {
            return Err(malformed(format!("benchmark schema {}", inputs.schema)));
        }
        Ok(inputs)
    }
}

/// Every labelled sector of `family` (line `i` at power 1 or 0) that is not
/// scaleless, as power vectors in increasing order of the line mask.
pub fn nonzero_sectors(family: &str) -> Result<Vec<Vec<i64>>, VakintError> {
    let data = registry_family(family)?;
    let lines = data.line_count();
    let mut sectors = Vec::new();
    for mask in 1u32..(1 << lines) {
        let powers = (0..lines)
            .map(|line| i64::from(mask >> line & 1))
            .collect::<Vec<_>>();
        if forcer_input(data, &powers, Atom::num(1).as_view(), &parse("p(1)")?)?.is_some() {
            sectors.push(powers);
        }
    }
    Ok(sectors)
}

fn power_label(powers: &[i64]) -> String {
    powers.iter().map(i64::to_string).collect()
}

/// The rank-2 scalar numerators of `loops` loops: `dot(k(i),p(1))*dot(k(j),p(1))`
/// and `dot(k(i),k(j))` for `i <= j`.
fn rank_two_numerators(loops: usize) -> Vec<String> {
    let mut numerators = Vec::new();
    for i in 1..=loops {
        for j in i..=loops {
            numerators.push(format!("dot(k({i}),p(1))*dot(k({j}),p(1))"));
            numerators.push(format!("dot(k({i}),k({j}))"));
        }
    }
    numerators
}

/// Number of dotted sectors per loop count in the grid fixtures (30 in all).
const DOTTED_SECTORS: [(usize, usize); 2] = [(2, 8), (3, 22)];

/// The cases of `p{loops}l_grid.json`, without outputs:
/// - `bench_<cell>_<family>`: every top of every benchmark cell of `loops`
///   loops (`benchmark` is the text of `benchmark_p_integrals.json`);
/// - `sector_<family>_<lines>_r0` and `_r2`: every non-scaleless labelled
///   sector of every family of `loops` loops (lines at power 1 or 0, the line
///   powers in the name), with numerator 1 and with one rank-2 numerator, the
///   `(7n+3) mod m`-th of the `m` rank-2 numerators for the `n`-th sector;
/// - `dotted_<family>_<lines>`: a spread sample of those sectors (2L: 8,
///   3L: 22, 30 in all) with every second positive line, starting with the
///   first or the second, at power 2.
pub fn grid_cases(loops: usize, benchmark: &str) -> Result<Vec<ForcerCase>, VakintError> {
    let mut cases = Vec::new();
    for cell in BenchmarkInputs::parse(benchmark)?
        .cells
        .iter()
        .filter(|cell| cell.loops == loops)
    {
        for top in &cell.tops {
            cases.push(ForcerCase::new(
                format!("bench_{}_{}", cell.name, top.family),
                &top.family,
                top.powers.clone(),
                &cell.numerator,
                "p(1)",
            ));
        }
    }
    let rank_two = rank_two_numerators(loops);
    let mut sectors = Vec::new();
    for family in P_FAMILIES.iter().filter(|family| family.loops == loops) {
        for powers in nonzero_sectors(family.name)? {
            sectors.push((family.name, powers));
        }
    }
    for (index, (family, powers)) in sectors.iter().enumerate() {
        let label = format!("sector_{family}_{}", power_label(powers));
        cases.push(ForcerCase::new(
            format!("{label}_r0"),
            family,
            powers.clone(),
            "1",
            "p(1)",
        ));
        cases.push(ForcerCase::new(
            format!("{label}_r2"),
            family,
            powers.clone(),
            &rank_two[(7 * index + 3) % rank_two.len()],
            "p(1)",
        ));
    }
    if let Some(&(_, sample)) = DOTTED_SECTORS.iter().find(|(l, _)| *l == loops) {
        for draw in 0..sample {
            let (family, powers) = &sectors[draw * sectors.len() / sample];
            let mut dotted = powers.clone();
            for (ordinal, power) in dotted.iter_mut().filter(|power| **power > 0).enumerate() {
                if (ordinal + draw) % 2 == 0 {
                    *power = 2;
                }
            }
            cases.push(ForcerCase::new(
                format!("dotted_{family}_{}", power_label(&dotted)),
                family,
                dotted,
                "1",
                "p(1)",
            ));
        }
    }
    Ok(cases)
}

/// The provenance line of the grid fixtures.
pub const GRID_PROVENANCE: &str = "M2-V2 grid: benchmark cells (tests/fixtures/benchmark_p_integrals.json), \
every non-scaleless labelled sector at rank 0 and 2, dotted sectors; \
fixtures::grid_cases in src/p_integral/forcer/fixtures.rs; VAKINT_FORCER_FIXTURES=write with FORM 5";

/// The provenance line of the anchors.
pub const ANCHOR_PROVENANCE: &str = "M2-V2 normalization anchors (M2 plan section 3.3); \
fixtures::anchor_cases in src/p_integral/forcer/fixtures.rs; VAKINT_FORCER_FIXTURES=write with FORM 5";

/// The normalization anchors of the M2 plan (section 3.3), without outputs:
/// - `p1l_11`: `P1L(1,1) = 1/ε` in the G-scheme;
/// - `p1l_kq_plus`, `p1l_kq_minus`: the bubble with `k·p` for `p(1) = ±p(1)`;
/// - `kite`: `P2L(1,1,1,1,1)`, whose finalized `ε^0` is `-6ζ3`;
/// - `kite_rank3_kp`, `kite_rank3_mixed`: rank-3 kites;
/// - `no`: `P3L_NO(1^8)` (the ISP at power 0), the master `no`;
/// - `t105_parent_*`: the plain three-loop kite with a doubled central line,
///   which integrates to `t1star05` with its ε-line at `n = 0` (the common
///   chart-), in the LA family (lines 1 and 3 are the central bubble) and in
///   the NO family, plus a dot on a bubble line (`n = 1`) and a rank-2
///   numerator.
pub fn anchor_cases() -> Vec<ForcerCase> {
    let kite_parent = vec![1, 1, 1, 0, 1, 1, 0, 1];
    let mut dotted_parent = kite_parent.clone();
    dotted_parent[0] = 2;
    vec![
        ForcerCase::new("p1l_11".into(), "P1L", vec![1, 1], "1", "p(1)"),
        ForcerCase::new(
            "p1l_kq_plus".into(),
            "P1L",
            vec![1, 1],
            "dot(k(1),p(1))",
            "p(1)",
        ),
        ForcerCase::new(
            "p1l_kq_minus".into(),
            "P1L",
            vec![1, 1],
            "dot(k(1),p(1))",
            "-p(1)",
        ),
        ForcerCase::new("kite".into(), "P2L", vec![1; 5], "1", "p(1)"),
        ForcerCase::new(
            "kite_rank3_kp".into(),
            "P2L",
            vec![1; 5],
            "dot(k(1),p(1))^2*dot(k(2),p(1))",
            "p(1)",
        ),
        ForcerCase::new(
            "kite_rank3_mixed".into(),
            "P2L",
            vec![1; 5],
            "dot(k(1),p(1))*dot(k(1),k(2))",
            "p(1)",
        ),
        ForcerCase::new("no".into(), "P3L_NO", vec![1; 8], "1", "p(1)"),
        ForcerCase::new(
            "t105_parent_la".into(),
            "P3L_LA",
            kite_parent.clone(),
            "1",
            "p(1)",
        ),
        ForcerCase::new(
            "t105_parent_no".into(),
            "P3L_NO",
            vec![1, 0, 1, 0, 1, 1, 1, 1],
            "1",
            "p(1)",
        ),
        ForcerCase::new(
            "t105_parent_la_dotted".into(),
            "P3L_LA",
            dotted_parent,
            "1",
            "p(1)",
        ),
        ForcerCase::new(
            "t105_parent_la_rank2".into(),
            "P3L_LA",
            kite_parent,
            "dot(k(1),p(1))*dot(k(2),p(1))",
            "p(1)",
        ),
    ]
}

/// The fixture at `path` after an optional regeneration: with `write`, the
/// cases `definitions` run through Forcer (FORM 5 at `form`, `jobs` FORM
/// processes at a time) and the file is rewritten; then the file is read and
/// its cases must be exactly `definitions` (outputs aside).
pub fn write_or_read(
    path: &std::path::Path,
    provenance: &str,
    definitions: Vec<ForcerCase>,
    write: Option<(&str, usize)>,
) -> Result<ForcerFixture, VakintError> {
    if let Some((form, jobs)) = write {
        let mut fixture = ForcerFixture::unwritten(provenance, definitions.clone());
        fixture.write_outputs(form, jobs)?;
        std::fs::write(path, fixture.to_json())?;
    }
    let text = std::fs::read_to_string(path)
        .map_err(|error| malformed(format!("{}: {error}", path.display())))?;
    let fixture = ForcerFixture::parse(&text)?;
    if fixture.definitions() != definitions {
        return Err(malformed(format!(
            "{}: the cases differ from their definition; regenerate with VAKINT_FORCER_FIXTURES=write",
            path.display()
        )));
    }
    if fixture.provenance.as_deref() != Some(provenance) {
        return Err(malformed(format!("{}: provenance differs", path.display())));
    }
    Ok(fixture)
}
