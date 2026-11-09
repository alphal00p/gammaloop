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
//!   `p4l_battery.json` extended by M3-V7 to all eleven four-loop families,
//!   with the depth of every case ([`CaseDepth`]), its value at
//!   [`I7_TERMS`] Laurent terms ([`ForcerCase::eleven_terms`]), M2's
//!   three-loop `Oep` boundary ([`ForcerFixture::boundary`]) and Forcer's
//!   whole-family blind spots ([`ForcerFixture::forcer_blind`]);
//! - `p{1,2,3}l_grid.json`: [`grid_cases`], the benchmark cells
//!   (`tests/fixtures/benchmark_p_integrals.json`), every non-scaleless
//!   labelled sector of the registry families at rank 0 and 2, and dotted
//!   sectors;
//! - `p_anchors.json`: [`anchor_cases`], the normalization anchors of the M2
//!   plan (section 3.3), including the parents of the eps-line kite `t1star05`;
//! - `p4l_anchors.json`: [`anchor_cases_4l`], the four-loop master anchors of
//!   the M3 plan (D-M3-5): every master integral of Forcer's four-loop basis
//!   at its Forcer definition (source line in `definition`), the ladder
//!   `la4_top`, and the `no3` and `bubu` tops;
//! - `p4l_grid.json`: [`grid_cases_4l`], the 165 four-loop benchmark tops and
//!   a seeded sample of labelled sectors, each with its depth.
//!
//! The depths and the finalized values are derived from the frozen outputs
//! by the shared p-finalizer ([`ForcerCase::finalized_at`], `Symbolic` mode,
//! [`finalize_settings`]); no FORM runs for them.
//!
//! This module is public for the integration tests only; it is not part of
//! vakint's API.

use std::sync::atomic::{AtomicUsize, Ordering};

use serde::{Deserialize, Serialize};
use symbolica::atom::{Atom, AtomCore, Symbol};

use super::{ForcerInput, combine_parts, forcer_input, run_forcer};
use crate::p_integral::{FORCER_UPSTREAM_COMMIT, P_FAMILIES, PFamilyData, ReducedPropagatorValue};
use crate::rustred_evaluation::RustRedMasterMode;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{LoopNormalizationFactor, Vakint, VakintError, VakintSettings};

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
    /// How the outputs were produced (files written since M3).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run: Option<ForcerRun>,
    pub cases: Vec<ForcerCase>,
    /// Cases of another loop count frozen with this file: in
    /// `p4l_battery.json`, M2's three-loop `Oep` boundary (the `P3L_NO` top,
    /// D-S2-2) with its value at [`I7_TERMS`] terms (I-7).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub boundary: Vec<ForcerCase>,
    /// The nodes Forcer cannot check whole-family (M3 plan section 2 item 4)
    /// and the cases of this file it does run for their masters.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub forcer_blind: Vec<ForcerBlind>,
}

/// A Forcer blind spot: a four-loop node whose whole-family comparison with
/// Forcer fails, and the per-integral cases Forcer runs for its master.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct ForcerBlind {
    /// The cascade node (`rules_study/targets.tsv`).
    pub node: String,
    /// Its Forcer master, `Master(master)`.
    pub master: String,
    /// What fails whole-family, and the source of that measurement.
    pub whole_family: String,
    /// Cases of the file whose Forcer output carries `Master(master)`.
    pub per_integral: Vec<String>,
}

/// Laurent terms of the I-7 values (M3, D-S2-2): `ε^-4..ε^6` at four loops,
/// `ε^-3..ε^7` at three.
pub const I7_TERMS: i64 = 11;

/// The most Laurent terms a depth is probed at: `ε^-4..ε^11` at four loops,
/// beyond every tabulated master order (`Masters.prc` stops at `ε^9`).
pub const DEPTH_CAP_TERMS: i64 = 16;

/// The settings of the frozen finalized values: the G-scheme, dot-product
/// notation, 40 digits, vakint's default symbols and `terms` Laurent terms.
pub fn finalize_settings(terms: i64) -> VakintSettings {
    VakintSettings {
        integral_normalization_factor: LoopNormalizationFactor::GScheme,
        number_of_terms_in_epsilon_expansion: terms,
        use_dot_product_notation: true,
        run_time_decimal_precision: 40,
        ..VakintSettings::default()
    }
}

/// A master whose unknown Laurent tail reaches a finalized value (one entry
/// of `VakintError::MissingMasterOrders`).
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct MissingOrder {
    /// The master as the `Oep` guard names it (`PMaster(name)`).
    pub master: String,
    /// Its first Laurent order that `Masters.prc` does not give (`ALARM`).
    pub first_missing_order: i64,
    /// The lowest order of the result that this tail affects.
    pub first_affected_result_order: i64,
}

/// A frozen value finalized at a Laurent depth in `Symbolic` mode
/// ([`ForcerCase::finalized_at`]).
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Finalized {
    /// The value's Laurent coefficients at unit scale
    /// ([`ForcerCase::unit_scale_coefficients`]).
    Symbolic(Vec<String>),
    /// The depth needs master orders that are unknown: an `Oep` case, an
    /// expected `VakintError::MissingMasterOrders`.
    MissingMasterOrders(Vec<MissingOrder>),
}

/// The deepest Laurent depth at which a frozen value finalizes (M3 plan G7:
/// "through the deepest order Forcer finalizes per case").
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct CaseDepth {
    /// `number_of_terms_in_epsilon_expansion` (from `ε^-L`), at most
    /// [`DEPTH_CAP_TERMS`].
    pub terms: i64,
    /// The highest order of the value, `ε^(terms-L-1)`.
    pub order: i64,
    /// The masters whose unknown orders stop one term more; empty when
    /// `terms` is the cap.
    pub limit: Vec<MissingOrder>,
    /// The `Symbolic` value at `terms` ([`finalize_settings`]): its Laurent
    /// coefficients `ε^-L..ε^order` at unit scale
    /// ([`ForcerCase::unit_scale_coefficients`]).
    pub symbolic: Vec<String>,
}

/// The run that wrote a fixture's outputs.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct ForcerRun {
    /// The FORM command and the Forcer call (`tform -w4`, default options).
    pub flags: String,
    /// FORM processes at a time.
    pub jobs: usize,
    /// Wall seconds of the whole regeneration, rounded.
    pub seconds: u64,
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
    /// For a master anchor: the Forcer master (`Master(name)`) the case is.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub master: Option<String>,
    /// For an anchor: its Forcer definition (source path:line, relative to
    /// the crate) and how the input realizes it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub definition: Option<String>,
    /// The deepest depth the frozen value finalizes at (M3 fixtures).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub depth: Option<CaseDepth>,
    /// The frozen value finalized at [`I7_TERMS`] terms (I-7): Forcer's
    /// `Symbolic` value, or the expected `MissingMasterOrders` of an `Oep`
    /// case.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub eleven_terms: Option<Finalized>,
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
            run: None,
            cases,
            boundary: Vec::new(),
            forcer_blind: Vec::new(),
        }
    }

    /// The cases with their outputs (and the values derived from them)
    /// removed, to compare with a definition.
    pub fn definitions(&self) -> Vec<ForcerCase> {
        self.cases.iter().map(ForcerCase::definition).collect()
    }

    /// The fixture as a definition, ready for [`Self::write_outputs_with`]:
    /// every case and boundary case without outputs, no FORM version, no
    /// run, the provenance `provenance`; the blind spots are kept.
    pub fn unwritten_like(&self, provenance: &str) -> Self {
        ForcerFixture {
            boundary: self.boundary.iter().map(ForcerCase::definition).collect(),
            forcer_blind: self.forcer_blind.clone(),
            ..Self::unwritten(provenance, self.definitions())
        }
    }

    /// Run Forcer (FORM 5 at `form`, one worker per FORM process) on every
    /// case with `jobs` FORM processes at a time and store the outputs and the
    /// FORM version.
    pub fn write_outputs(&mut self, form: &str, jobs: usize) -> Result<(), VakintError> {
        self.write_outputs_with(form, jobs, 1)
    }

    /// [`Self::write_outputs`] with `workers` TFORM workers per FORM process
    /// (`tform -w<workers>`; 1 is plain FORM); records the run.
    pub fn write_outputs_with(
        &mut self,
        form: &str,
        jobs: usize,
        workers: usize,
    ) -> Result<(), VakintError> {
        let started = std::time::Instant::now();
        self.form_version =
            super::forcer_form_version(form).map_err(|reason| VakintError::ForcerUnavailable {
                form_exe_path: form.into(),
                reason,
            })?;
        let programs = self
            .cases
            .iter()
            .chain(&self.boundary)
            .map(ForcerCase::input)
            .collect::<Result<Vec<_>, _>>()?;
        let settings = VakintSettings {
            form_exe_path: form.into(),
            n_cores: workers.max(1),
            ..VakintSettings::default()
        };
        let command = super::form_command(&settings);
        let program = std::path::Path::new(command.get_program())
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .into_owned();
        let flags = std::iter::once(program)
            .chain(
                command
                    .get_args()
                    .map(|arg| arg.to_string_lossy().into_owned()),
            )
            .collect::<Vec<_>>()
            .join(" ");
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
        for (case, output) in self.cases.iter_mut().chain(&mut self.boundary).zip(outputs) {
            case.outputs = output
                .expect("every case ran")
                .map_err(|error| malformed(format!("{}: {error}", case.name)))?;
        }
        self.run = Some(ForcerRun {
            flags: format!(
                "{flags} run_forcer.frm; #call Forcer() with default options (exact G-scheme reduction, no expansion)"
            ),
            jobs: jobs.max(1),
            seconds: started.elapsed().as_secs_f64().round() as u64,
        });
        Ok(())
    }

    /// Store the depth ([`ForcerCase::deepest`]) of every case and boundary
    /// case and, with `eleven`, its value at [`I7_TERMS`] terms, from the
    /// frozen outputs (no FORM), finalizing on `workers` threads.
    pub fn write_depths(&mut self, eleven: bool, workers: usize) -> Result<(), VakintError> {
        let next = AtomicUsize::new(0);
        let cases = self
            .cases
            .iter()
            .chain(&self.boundary)
            .cloned()
            .collect::<Vec<_>>();
        let derive = |case: &ForcerCase| -> Result<_, VakintError> {
            let depth = case.deepest()?;
            let at_eleven = if !eleven {
                None
            } else if depth.terms == I7_TERMS {
                Some(Finalized::Symbolic(depth.symbolic.clone()))
            } else {
                Some(case.finalized_at(I7_TERMS)?)
            };
            Ok((depth, at_eleven))
        };
        let mut derived = (0..cases.len()).map(|_| None).collect::<Vec<_>>();
        std::thread::scope(|scope| {
            let threads = (0..workers.max(1))
                .map(|_| {
                    scope.spawn(|| {
                        let mut done = Vec::new();
                        loop {
                            let index = next.fetch_add(1, Ordering::Relaxed);
                            let Some(case) = cases.get(index) else {
                                return done;
                            };
                            done.push((index, derive(case)));
                        }
                    })
                })
                .collect::<Vec<_>>();
            for thread in threads {
                for (index, result) in thread.join().expect("depth worker") {
                    derived[index] = Some(result);
                }
            }
        });
        for (case, result) in self.cases.iter_mut().chain(&mut self.boundary).zip(derived) {
            let (depth, at_eleven) = result
                .expect("every case finalized")
                .map_err(|error| malformed(format!("{}: {error}", case.name)))?;
            case.depth = Some(depth);
            case.eleven_terms = at_eleven;
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
            master: None,
            definition: None,
            depth: None,
            eleven_terms: None,
        }
    }

    /// The case without its outputs and the values derived from them.
    pub fn definition(&self) -> Self {
        ForcerCase {
            outputs: Vec::new(),
            depth: None,
            eleven_terms: None,
            ..self.clone()
        }
    }

    /// An anchor case: `self` labelled with the Forcer master it is (if any)
    /// and its definition.
    pub fn anchoring(mut self, master: Option<&str>, definition: &str) -> Self {
        self.master = master.map(Into::into);
        self.definition = Some(definition.into());
        self
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
        mode: RustRedMasterMode,
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

    /// The Laurent coefficients `ε^-L..ε^(terms-L-1)` of `value`, a value of
    /// this case finalized at `terms` terms with [`finalize_settings`], at
    /// unit scale: `Q² = μ² = 1` (`dot(p,p) = -1`, `mursq = 1`, so every
    /// logarithm of the scale vanishes), each in Symbolica's canonical text.
    /// The scale dependence is fixed by dimensional analysis; this is what the
    /// numerical comparisons evaluate (`tests/support/p_backends.rs`).
    pub fn unit_scale_coefficients(
        &self,
        value: &Atom,
        terms: i64,
    ) -> Result<Vec<String>, VakintError> {
        let settings = finalize_settings(terms);
        let loops = self.loops()? as i64;
        let (head, _) = crate::p_integral::external_head(&parse(&self.external)?)?;
        let e = Atom::var(vk_symbol!(&settings.epsilon_symbol));
        let shifted = (value
            .replace(Atom::var(vk_symbol!(&settings.mu_r_sq_symbol)).to_pattern())
            .with(Atom::num(1).to_pattern())
            .replace(crate::symbols::S.dot.call((&head, &head)).to_pattern())
            .with(Atom::num(-1).to_pattern())
            .replace(Symbol::LOG.call(Atom::num(1)).to_pattern())
            .with(Atom::Zero.to_pattern())
            * e.pow(Atom::num(loops)))
        .expand();
        let mut coefficients = vec![Atom::Zero; terms as usize];
        for (monomial, coefficient) in shifted.coefficient_list::<i8>(std::slice::from_ref(&e)) {
            let power = match monomial.as_view() {
                _ if monomial == Atom::num(1) => 0,
                _ if monomial == e => 1,
                symbolica::atom::AtomView::Pow(power) => {
                    crate::get_integer_from_atom(power.get_base_exp().1).unwrap_or(-1)
                }
                _ => -1,
            };
            let slot = usize::try_from(power)
                .ok()
                .and_then(|power| coefficients.get_mut(power))
                .ok_or_else(|| {
                    malformed(format!(
                        "{}: ε monomial {monomial} at unit scale",
                        self.name
                    ))
                })?;
            *slot = coefficient.expand();
        }
        if let Some(scale) = coefficients
            .iter()
            .find(|coefficient| coefficient.to_string().contains("log("))
        {
            return Err(malformed(format!(
                "{}: a scale logarithm at unit scale: {scale}",
                self.name
            )));
        }
        Ok(coefficients
            .iter()
            .map(AtomCore::to_canonical_string)
            .collect())
    }

    /// The frozen value finalized at `terms` Laurent terms in `Symbolic` mode
    /// with [`finalize_settings`], as unit-scale Laurent coefficients; a
    /// missing master order is a value too.
    pub fn finalized_at(&self, terms: i64) -> Result<Finalized, VakintError> {
        match self.finalized(&finalize_settings(terms), RustRedMasterMode::Symbolic) {
            Ok(value) => Ok(Finalized::Symbolic(
                self.unit_scale_coefficients(&value, terms)?,
            )),
            Err(VakintError::MissingMasterOrders(missing)) => Ok(Finalized::MissingMasterOrders(
                missing
                    .masters
                    .iter()
                    .map(|master| MissingOrder {
                        master: master.master.clone(),
                        first_missing_order: master.first_missing_order,
                        first_affected_result_order: master.first_affected_result_order,
                    })
                    .collect(),
            )),
            Err(error) => Err(error),
        }
    }

    /// The deepest depth at which the frozen value finalizes, at most
    /// [`DEPTH_CAP_TERMS`]: the cap, or one order below the lowest result
    /// order an unknown master tail affects at the cap.
    pub fn deepest(&self) -> Result<CaseDepth, VakintError> {
        let loops = self.loops()? as i64;
        let (terms, limit) = match self.finalized_at(DEPTH_CAP_TERMS)? {
            Finalized::Symbolic(symbolic) => {
                return Ok(CaseDepth {
                    terms: DEPTH_CAP_TERMS,
                    order: DEPTH_CAP_TERMS - loops - 1,
                    limit: Vec::new(),
                    symbolic,
                });
            }
            Finalized::MissingMasterOrders(missing) => {
                let first = missing
                    .iter()
                    .map(|master| master.first_affected_result_order)
                    .min()
                    .ok_or_else(|| malformed(format!("{}: no missing master", self.name)))?;
                let limit = missing
                    .into_iter()
                    .filter(|master| master.first_affected_result_order == first)
                    .collect();
                // Through ε^(first-1): terms - L - 1 = first - 1.
                (first + loops, limit)
            }
        };
        match self.finalized_at(terms)? {
            Finalized::Symbolic(symbolic) if terms > 0 => Ok(CaseDepth {
                terms,
                order: terms - loops - 1,
                limit,
                symbolic,
            }),
            other => Err(malformed(format!(
                "{}: does not finalize at {terms} terms: {other:?}",
                self.name
            ))),
        }
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

/// The labelled sectors (family and line powers) in which the fixed rank-2
/// numerator of [`grid_cases`] gives an exact 0, by cancelling a line into a
/// scaleless integral or by a symmetry of the sector, and how many numerators
/// further the first nonzero one is. A zero would carry no rank-2 information;
/// Forcer's silent zero for an unrecognized input is a different thing, refused
/// by `run_forcer`.
const RANK_TWO_STEPS: [(&str, usize); 29] = [
    ("P1L_11", 1),
    ("P2L_11001", 1),
    ("P2L_11011", 1),
    ("P2L_10111", 1),
    ("P3L_LA_11100110", 1),
    ("P3L_LA_10111110", 1),
    ("P3L_LA_10101101", 1),
    ("P3L_LA_11011101", 1),
    ("P3L_LA_11110011", 1),
    ("P3L_LA_01011111", 3),
    ("P3L_LA_10111111", 1),
    ("P3L_BE_11110010", 1),
    ("P3L_BE_11110110", 1),
    ("P3L_BE_10111110", 1),
    ("P3L_BE_11110001", 1),
    ("P3L_BE_11011001", 1),
    ("P3L_BE_11110101", 1),
    ("P3L_BE_01011101", 1),
    ("P3L_BE_11100011", 1),
    ("P3L_BE_01110111", 1),
    ("P3L_BE_11001111", 1),
    ("P3L_BE_11011111", 1),
    ("P3L_NO_11111010", 1),
    ("P3L_NO_10110110", 1),
    ("P3L_NO_11101001", 1),
    ("P3L_NO_11101101", 1),
    ("P3L_NO_10111101", 1),
    ("P3L_NO_00111011", 1),
    ("P3L_NO_00101111", 1),
];

/// Number of dotted sectors per loop count in the grid fixtures (30 in all).
const DOTTED_SECTORS: [(usize, usize); 2] = [(2, 8), (3, 22)];

/// The cases of `p{loops}l_grid.json`, without outputs:
/// - `bench_<cell>_<family>`: every top of every benchmark cell of `loops`
///   loops (`benchmark` is the text of `benchmark_p_integrals.json`);
/// - `sector_<family>_<lines>_r0` and `_r2`: every non-scaleless labelled
///   sector of every family of `loops` loops (lines at power 1 or 0, the line
///   powers in the name), with numerator 1 and with one rank-2 numerator, the
///   `(7n+3+s) mod m`-th of the `m` rank-2 numerators for the `n`-th sector,
///   `s` from [`RANK_TWO_STEPS`] (0 for most sectors), so that it is nonzero;
/// - `dotted_<family>_<lines>`: a spread sample of those sectors (2L: 8,
///   3L: 22, 30 in all) with every second positive line, starting with the
///   first or the second, at power 2.
pub fn grid_cases(loops: usize, benchmark: &str) -> Result<Vec<ForcerCase>, VakintError> {
    let mut cases = benchmark_tops(loops, benchmark)?;
    let rank_two = rank_two_numerators(loops);
    let mut sectors = Vec::new();
    for family in P_FAMILIES.iter().filter(|family| family.loops == loops) {
        for powers in nonzero_sectors(family.name)? {
            sectors.push((family.name, powers));
        }
    }
    for (index, (family, powers)) in sectors.iter().enumerate() {
        let sector = format!("{family}_{}", power_label(powers));
        let step = RANK_TWO_STEPS
            .iter()
            .find(|(listed, _)| *listed == sector)
            .map_or(0, |(_, step)| *step);
        let label = format!("sector_{sector}");
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
            &rank_two[(7 * index + 3 + step) % rank_two.len()],
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

/// `bench_<cell>_<family>`: every top of every benchmark cell of `loops`
/// loops (`benchmark` is the text of `benchmark_p_integrals.json`).
fn benchmark_tops(loops: usize, benchmark: &str) -> Result<Vec<ForcerCase>, VakintError> {
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
    Ok(cases)
}

/// The seed of the sampled sectors of `p4l_grid.json` ([`grid_cases_4l`]).
pub const GRID4_SEED: u64 = 20_261_109;

/// SplitMix64: the reproducible draw of [`grid_cases_4l`].
struct SplitMix64(u64);

impl SplitMix64 {
    fn below(&mut self, bound: usize) -> usize {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        ((z ^ (z >> 31)) % bound as u64) as usize
    }
}

/// The cases of `p4l_grid.json`, without outputs (M3-V7; not all 3,606
/// labelled four-loop sectors, which would take hours of FORM):
/// - `bench_<cell>_<family>`: the 165 tops of the 15 four-loop benchmark
///   cells;
/// - per four-loop family, three distinct non-scaleless labelled sectors
///   other than the top, drawn with [`GRID4_SEED`] (SplitMix64 seeded with
///   the seed plus the family's index among the four-loop families):
///   `sector_<family>_<lines>_r0` and `_r2` for the first two (numerator 1
///   and a drawn rank-2 numerator, shifted by [`RANK_TWO_STEPS_4L`] where it
///   would give an exact 0), and `dotted_<family>_<lines>` for the third,
///   every second positive line from the first at power 2.
pub fn grid_cases_4l(benchmark: &str) -> Result<Vec<ForcerCase>, VakintError> {
    let mut cases = benchmark_tops(4, benchmark)?;
    let rank_two = rank_two_numerators(4);
    for (index, family) in P_FAMILIES
        .iter()
        .filter(|family| family.loops == 4)
        .enumerate()
    {
        let sectors = nonzero_sectors(family.name)?
            .into_iter()
            .filter(|powers| powers.contains(&0))
            .collect::<Vec<_>>();
        let mut draw = SplitMix64(GRID4_SEED + index as u64);
        let mut picked = Vec::<usize>::new();
        while picked.len() < 3 {
            let sector = draw.below(sectors.len());
            if !picked.contains(&sector) {
                picked.push(sector);
            }
        }
        for &sector in &picked[..2] {
            let powers = &sectors[sector];
            let label = format!("sector_{}_{}", family.name, power_label(powers));
            let step = RANK_TWO_STEPS_4L
                .iter()
                .find(|(listed, _)| *listed == label)
                .map_or(0, |(_, step)| *step);
            let numerator = &rank_two[(draw.below(rank_two.len()) + step) % rank_two.len()];
            cases.push(ForcerCase::new(
                format!("{label}_r0"),
                family.name,
                powers.clone(),
                "1",
                "p(1)",
            ));
            cases.push(ForcerCase::new(
                format!("{label}_r2"),
                family.name,
                powers.clone(),
                numerator,
                "p(1)",
            ));
        }
        let mut dotted = sectors[picked[2]].clone();
        for (ordinal, power) in dotted.iter_mut().filter(|power| **power > 0).enumerate() {
            if ordinal % 2 == 0 {
                *power = 2;
            }
        }
        cases.push(ForcerCase::new(
            format!("dotted_{}_{}", family.name, power_label(&dotted)),
            family.name,
            dotted,
            "1",
            "p(1)",
        ));
    }
    Ok(cases)
}

/// The sampled sectors of [`grid_cases_4l`] whose drawn rank-2 numerator
/// gives an exact 0, and how many numerators further the first nonzero one
/// is.
const RANK_TWO_STEPS_4L: [(&str, usize); 2] = [
    ("sector_P4L_J11_11011000101", 1),
    ("sector_P4L_NO1_01111001101", 1),
];

/// The provenance line of `p4l_grid.json`.
pub const GRID4_PROVENANCE: &str = "M3-V7 four-loop grid: the 165 tops of the 15 four-loop benchmark cells \
(tests/fixtures/benchmark_p_integrals.json) and a seeded sample of labelled sectors (rank 0 and 2, dotted; \
fixtures::grid_cases_4l, GRID4_SEED) in src/p_integral/forcer/fixtures.rs; VAKINT_FORCER_FIXTURES=write with \
FORM 5 (run: flags, jobs, seconds); per case, depth: the deepest Laurent depth (at most DEPTH_CAP_TERMS = 16 \
terms) at which the frozen value finalizes in Symbolic mode (G-scheme, fixtures::finalize_settings), its order, \
the masters that stop one term more and the Symbolic value there (ForcerCase::deepest, no FORM); \
outputs are raw Forcer C_G before the Wick sign: reduced = (-1)^(sum a + r) C_G (ForcerCase::frozen_values)";

/// The provenance line of `p4l_battery.json`.
pub const BATTERY4_PROVENANCE: &str = "M1 four-loop battery (hand-picked cases) extended by M3-V7 to all \
eleven four-loop families (top, ISP at rank 2, rank 3 with p(1), a dot); boundary: M2's three-loop Oep \
boundary case (P3L_NO top, D-S2-2); forcer_blind: the nodes Forcer cannot check whole-family \
(rules_study/M2_M3_INPUTS.md section 4a) and the per-integral cases it runs; VAKINT_FORCER_FIXTURES=write \
with FORM 5 through tests/p_fixtures4_tests.rs (run: flags, jobs, seconds); per case, depth (as in \
p4l_grid.json, ForcerCase::deepest) and eleven_terms: the Symbolic value at 11 Laurent terms (I-7) or, for an \
Oep case, the expected MissingMasterOrders (ForcerCase::finalized_at, G-scheme, no FORM); \
outputs are raw Forcer C_G before the Wick sign: reduced = (-1)^(sum a + r) C_G (ForcerCase::frozen_values)";

/// The provenance line of the grid fixtures.
pub const GRID_PROVENANCE: &str = "M2-V2 grid: benchmark cells (tests/fixtures/benchmark_p_integrals.json), \
every non-scaleless labelled sector at rank 0 and 2, dotted sectors; \
fixtures::grid_cases in src/p_integral/forcer/fixtures.rs; VAKINT_FORCER_FIXTURES=write with FORM 5; \
outputs are raw Forcer C_G before the Wick sign: reduced = (-1)^(sum a + r) C_G (ForcerCase::frozen_values)";

/// The provenance line of the anchors.
pub const ANCHOR_PROVENANCE: &str = "M2-V2 normalization anchors (M2 plan section 3.3); \
fixtures::anchor_cases in src/p_integral/forcer/fixtures.rs; VAKINT_FORCER_FIXTURES=write with FORM 5; \
outputs are raw Forcer C_G before the Wick sign: reduced = (-1)^(sum a + r) C_G (ForcerCase::frozen_values)";

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

/// The provenance line of the four-loop anchors.
pub const ANCHOR4_PROVENANCE: &str = "M3-V7 four-loop normalization anchors (M3 plan D-M3-5, G5): \
every master of Forcer's four-loop basis at its Forcer definition (case field definition: source path:line \
relative to crates/vakint and the realization), la4_top (finalized eps^0 = -70 zeta(7) in the G-scheme), \
the no3 and bubu tops; fixtures::anchor_cases_4l in src/p_integral/forcer/fixtures.rs; \
VAKINT_FORCER_FIXTURES=write with FORM 5 (run: flags, jobs, seconds); \
outputs are raw Forcer C_G before the Wick sign: reduced = (-1)^(sum a + r) C_G (ForcerCase::frozen_values)";

/// Forcer's definition of each four-loop-relevant master that is not already
/// anchored by [`anchor_cases`] (`no`, `t105`), with the labelled sector that
/// realizes it: `(master, family, powers, definition)`. The sectors are the
/// first ones (registry order, then line mask) whose contracted graph is the
/// master's Forcer topology (`forcer_dag_68a30bf.tsv` notation), found by
/// graph isomorphism; an eps unit on a line (`*`) is realized at four loops
/// by a line in series with a one-loop bubble (power `1+ε`), and the two
/// units of `t155` by a bubble, a line and a bubble in series (power
/// `1+2ε`).
const FOUR_LOOP_MASTERS: [(&str, &str, [i64; 11], &str); 14] = [
    (
        "bebe",
        "P4L_J1",
        [1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1],
        "form_src/forcer/forcer/forcer-manual.h:9066: Z(1,1,1,1,1,1,1,1,0,0,0,0,0,0) = Master(bebe); \
         topology bebe 0-1 4-5 1-4 3-5 1-3 1-2 2-4 2-3 3-4 5-9",
    ),
    (
        "cross",
        "P4L_J2",
        [1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1],
        "form_src/forcer/forcer/forcer-manual.h:7383: Z(1,1,1,1,1,1,1,1,0,0,0,0,0,0) = Master(cross); \
         topology cross 0-1 3-5 1-3 4-5 1-4 2-5 1-2 2-3 2-4 5-9",
    ),
    (
        "haha",
        "P4L_HAHA",
        [1; 11],
        "form_src/forcer/forcer/forcer-manual.h:1886: Z(1,1,1,1,1,1,1,1,1,1,1,0,0,0) = Master(haha); \
         topology haha 0-3 5-8 4-5 3-4 2-3 1-2 1-8 5-6 1-6 4-7 2-7 6-7 8-9",
    ),
    (
        "lala",
        "P4L_LA4",
        [1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1],
        "form_src/forcer/forcer/forcer-manual.h:5008: Z(1,1,1,1,1,1,1,1,1,0,0,0,0,0) = Master(lala); \
         topology lala 0-3 5-8 4-5 3-4 2-3 1-2 1-8 2-4 1-4 1-5 8-9",
    ),
    (
        "no1",
        "P4L_NO1",
        [1; 11],
        "form_src/forcer/forcer/forcer-manual.h:2594, :2687: Z(1,1,1,1,1,1,1,1,1,1,1,0,0,0) = Master(no1); \
         topology no1 0-1 4-8 3-4 2-3 1-2 1-5 5-6 6-7 7-8 3-7 2-6 4-5 8-9",
    ),
    (
        "no2",
        "P4L_NO2",
        [1; 11],
        "form_src/forcer/forcer/forcer-rules.h:57: Zno2(1,1,1,1,1,1,1,1,1,1,1,0,0,0) = Master(no2) / Q2^3; \
         topology no2 0-1 5-8 4-5 2-4 1-2 1-6 3-6 3-7 7-8 5-6 3-4 2-7 8-9",
    ),
    (
        "no6",
        "P4L_NO1",
        [1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1],
        "form_src/forcer/forcer/forcer-manual.h:3851: Z(1,1,1,1,1,1,1,1,1,1,0,0,0,0) = Master(no6); \
         topology no6 0-1 7-8 6-7 2-6 1-2 1-7 3-7 3-4 4-8 4-6 2-3 8-9",
    ),
    (
        "nono",
        "P4L_J4",
        [1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1],
        "form_src/forcer/forcer/forcer-manual.h:6421: Z(1,1,1,1,1,1,1,1,1,0,0,0,0,0) = Master(nono); \
         topology nono 0-1 4-6 2-4 1-2 1-3 3-5 5-6 4-5 3-4 2-5 6-9",
    ),
    (
        "fastar2",
        "P4L_J2",
        [1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1],
        "form_src/forcer/forcer/forcer-rules.h:182 (FASTAR2SCHEME 1, the default): \
         Zfastar2(1,1,1,1,1,1,1,0,0) = Master(fastar2) / Q2 (scheme 0: forcer-manual.h:665); \
         topology fastar2 0-1 3-4 2-3* 1-2 1-5 4-5 3-5 2-5 4-9; the eps line 2-3* is a line and a one-loop bubble in series",
    ),
    (
        "nostar5",
        "P4L_J4",
        [1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1],
        "form_src/forcer/forcer/forcer-rules.h:201 (NOSTAR5SCHEME 1, the default): \
         Znostar5(1,1,1,1,1,1,1,1,0) = BasicNOIntegral5 / Q2^2 (scheme 0: :196, XXX: normalization?); \
         form_src/forcer/forcer.h:278: BasicNOIntegral5 = Master(nostar5); \
         topology nostar5 0-1 3-4 2-3 1-2 1-6 5-6* 4-5 3-6 2-5 4-9; the eps line 5-6* is a line and a one-loop bubble in series",
    ),
    (
        "nostar6",
        "P4L_J4",
        [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        "form_src/forcer/forcer/forcer-rules.h:221 (NOSTAR6SCHEME 1, the default): \
         Znostar6(1,1,1,1,1,1,1,1,0) = BasicNOIntegral6 / Q2^2 (scheme 0: :216, XXX: normalization?); \
         form_src/forcer/forcer.h:279: BasicNOIntegral6 = Master(nostar6); \
         topology nostar6 0-1 3-4 2-3 1-2 1-6 5-6 4-5* 3-6 2-5 4-9; the eps line 4-5* is a line and a one-loop bubble in series",
    ),
    (
        "t124",
        "P4L_J6",
        [0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1],
        "form_src/forcer/forcer/forcer-rules.h:252: Z(0,1,0,1,0,1,1,1,1,1) = BasicT1Integral24 / Q2 * Q2^5 \
         (XXX: normalization?); form_src/forcer/forcer.h:282: BasicT1Integral24 = Master(t124); \
         topology t1star24 0-1 2-3 1-2* 1-4 3-4* 2-4 3-9; each eps line is a line and a one-loop bubble in series",
    ),
    (
        "t145",
        "P4L_J1",
        [1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1],
        "form_src/forcer/forcer/forcer-rules.h:271: Z(0,0,0,1,1,1,1,1,1,1) = BasicT1Integral45 / Q2 * Q2^5 \
         (XXX: normalization?); form_src/forcer/forcer.h:281: BasicT1Integral45 = Master(t145); \
         topology t1star45 0-1 2-3 1-2 1-4 3-4* 2-4* 3-9; each eps line is a line and a one-loop bubble in series",
    ),
    (
        "t155",
        "P4L_J6",
        [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1],
        "form_src/forcer/forcer/forcer-rules.h:281: Z(0,0,0,0,2,1,1,1,1,1) = BasicT1Integral55 / Q2 * Q2^5 \
         (XXX: normalization?); form_src/forcer/forcer.h:280: BasicT1Integral55 = Master(t155); \
         topology t1star55 0-1 2-3 1-2 1-4 3-4 2-4** 3-9: line 2-4 at power 1+2eps (unit power, two eps units); \
         realized as a one-loop bubble, a line and a one-loop bubble in series",
    ),
];

/// The four-loop anchors of `p4l_anchors.json`, without outputs (M3 plan
/// D-M3-5): one case per [`FOUR_LOOP_MASTERS`] entry, named after its master,
/// with numerator 1 and `p(1)`; then
/// - `la4_top`: the four-loop ladder `P4L_LA4(1^11)`, whose finalized `ε^0`
///   is `-70ζ7` in the G-scheme (the case of `p4l_battery.json`);
/// - `no3_top`: `P4L_NO3(1^11)`, Forcer's `no3` with its `SchemeZno3`
///   override (the diamond lane's check);
/// - `bubu_top`: Forcer's knock-out topology `bubu` at unit powers
///   (forcer-rules.h:62-63), reduced onto `bebe` and `cross`.
pub fn anchor_cases_4l() -> Vec<ForcerCase> {
    let mut cases = FOUR_LOOP_MASTERS
        .iter()
        .map(|(master, family, powers, definition)| {
            ForcerCase::new((*master).into(), family, powers.to_vec(), "1", "p(1)")
                .anchoring(Some(master), definition)
        })
        .collect::<Vec<_>>();
    cases.push(
        ForcerCase::new("la4_top".into(), "P4L_LA4", vec![1; 11], "1", "p(1)").anchoring(
            None,
            "the four-loop ladder (topology la4 0-1 1-2 1-5 2-3 2-5 3-4 3-6 4-7 4-8 5-6 6-7 7-8 8-9): \
             C(8,4) zeta(7) = 70 zeta(7) in the G-scheme, finalized eps^0 = -70 zeta(7) \
             (src/p_integral/forcer/tests.rs, forcer_four_loop_battery)",
        ),
    );
    cases.push(
        ForcerCase::new("no3_top".into(), "P4L_NO3", vec![1; 11], "1", "p(1)").anchoring(
            None,
            "form_src/forcer/forcer/forcer-rules.h:7 (HASSCHEMEZno3, pathforcer.h:347) and SchemeZno3.prc: Forcer's no3 override; \
             topology no3 0-1 4-8 7-8 1-4 1-5 3-4 2-7 6-7 2-5 5-6 2-3 3-6 8-9 (diamond)",
        ),
    );
    cases.push(
        ForcerCase::new(
            "bubu_top".into(),
            "P4L_J4",
            vec![1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1],
            "1",
            "p(1)",
        )
        .anchoring(
            None,
            "form_src/forcer/forcer/forcer-rules.h:62-63: #case bubu, ReduceImpl2(bubu,Zbubu,14), no master of its own (knock-out); \
             topology bubu 0-4 4-5 2-4 3-4 5-6 2-5 3-5 1-6 1-2 1-3 6-9",
        ),
    );
    cases
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
    write_or_read_with(
        path,
        provenance,
        definitions,
        write.map(|(form, jobs)| (form, jobs, 1)),
    )
}

/// [`write_or_read`] with `(form, jobs, workers)`: each FORM process of a
/// regeneration runs `tform -w<workers>` (plain FORM for 1).
pub fn write_or_read_with(
    path: &std::path::Path,
    provenance: &str,
    definitions: Vec<ForcerCase>,
    write: Option<(&str, usize, usize)>,
) -> Result<ForcerFixture, VakintError> {
    if let Some((form, jobs, workers)) = write {
        let mut fixture = ForcerFixture::unwritten(provenance, definitions.clone());
        fixture.write_outputs_with(form, jobs, workers)?;
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
